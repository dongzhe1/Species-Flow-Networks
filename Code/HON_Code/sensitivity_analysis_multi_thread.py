"""
Sensitivity Analysis for Ballast Water Discharge Uncertainty
============================================================
Integrates pipelines from 02, 03, and 04 (first half) without
intermediate file I/O to maximize speed.

Optimizations vs. naive version:
  1. Fouling HON is computed ONCE (fouling does not use ballast discharge,
     so it is identical across all perturbation runs).
  2. Monte Carlo iterations run in parallel via multiprocessing.Pool,
     using all available CPU cores.

Perturbation applies ONLY to ballast water discharge values for non-US
ports, since that is the quantity estimated by the random forest model.

Both ballast and fouling results are reported separately.

Parameters (edit in CONFIG section below):
    SIGMA_VALUES     : perturbation magnitudes (e.g. [0.2, 0.5])
    N_SIMULATIONS    : number of Monte Carlo iterations
    N_WORKERS        : parallel processes (None = use all CPU cores)
    STAY_ATTENUATION : bool, corresponds to 'stayAttenuation' in original code
    TRADE_ADJUST     : bool, corresponds to 'tradeAdjust' in original code
    TOP_N            : number of top countries for stability reporting
"""

import csv
import math
import sys
import numpy as np
from collections import defaultdict
from decimal import Decimal, getcontext
from scipy.stats import spearmanr
import multiprocessing as mp
import gc

sys.setrecursionlimit(5000)
getcontext().prec = 50

# ============================================================
# CONFIG — edit these parameters
# ============================================================
SIGMA_VALUES     = [0.2]   # perturbation magnitudes to test
N_SIMULATIONS    = 8          # Monte Carlo iterations
N_WORKERS        = 8          # None = use all CPU cores (recommended)
STAY_ATTENUATION = False         # 'stayAttenuation' in original code
TRADE_ADJUST     = False         # 'tradeAdjust' in original code
TOP_N            = 10            # top-N countries for stability metric
SEED             = 42

RUN_BALLAST      = True          # run ballast Monte Carlo sensitivity analysis
RUN_FOULING      = False         # run fouling analysis (invariant to ballast
                                 # discharge perturbation — set True only if


# Fixed model parameters (do not change — must match original)
ORDER       = 16
R           = 'r0'
ENV         = 'env_'
ECO         = 'Eco_'
MIN_SUPPORT = 1e-20
B_EFFICACY  = {'r0': 0, 'r1': 0.7661, 'r2': 0.99150495, 'r3': 0.999371131}

# File paths — adjust to your environment
BASE                 = '../../data/'
clean_move           = BASE + 'moves/moves_cleaned_2018.txt'
port_data_file       = BASE + 'Places_allportdata_mergedSept2017.csv'
stay_file            = BASE + 'stay/stay.csv'
port_to_country_file = BASE + 'stay/portToCountry.csv'
export_file          = BASE + 'trade/Export_change_for_BW.csv'
import_file          = BASE + 'trade/Import_change_for_fouling.csv'

# ============================================================
# GLOBAL VARIABLES FOR WSL/LINUX MULTIPROCESSING
# ============================================================
global_probs_raw = {}
global_portToCountry = {}
global_exportDict = {}
global_importDict = {}

# ============================================================
# UTILITY FUNCTIONS  (from 02, logic unchanged)
# ============================================================

def get_port_data(fn, field, delim):
    ports = {}
    with open(fn) as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            ports[row[field]] = row
    return ports

def deg2rad(deg):
    return deg * (math.pi / 180)

def get_distance_km(ports, source, dest):
    lat1 = float(ports[source]['LATITUDE_DECIMAL'])
    lon1 = float(ports[source]['LONGITUDE_DECIMAL'])
    lat2 = float(ports[dest]['LATITUDE_DECIMAL'])
    lon2 = float(ports[dest]['LONGITUDE_DECIMAL'])
    R    = 6371
    dLat = deg2rad(lat2 - lat1)
    dLon = deg2rad(lon2 - lon1)
    a    = (math.sin(dLat/2)**2
            + math.cos(deg2rad(lat1)) * math.cos(deg2rad(lat2)) * math.sin(dLon/2)**2)
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

def filterbyeco(ports, source, dest):
    p_alien = 1
    s_eco   = ports[source]['MEOW_region']
    d_eco   = ports[dest]['MEOW_region']
    if s_eco != 'NA' and d_eco != 'NA':
        if s_eco == d_eco:
            p_alien = 0
        elif s_eco in ports[dest]['MEOW_Neighbors'].split('|'):
            p_alien = 0
    elif s_eco == 'NA' and d_eco != 'NA':
        s_eco = ports[source]['FEOW_region']
        if s_eco != 'NA' and (s_eco in ports[dest]['FEOW_Neighbors'].split('|')
                               or s_eco in ports[dest]['MEOW_Neighbors'].split('|')):
            p_alien = 0
    elif d_eco == 'NA' and s_eco != 'NA':
        d_eco = ports[dest]['FEOW_region']
        if d_eco != 'NA' and (s_eco in ports[dest]['FEOW_Neighbors'].split('|')
                               or s_eco in ports[dest]['MEOW_Neighbors'].split('|')):
            p_alien = 0
    return p_alien

def filterbyeco_same_only(ports, source, dest):
    try:
        s_eco = ports[source]['MEOW_region']
        d_eco = ports[dest]['MEOW_region']
        if s_eco != 'NA' and d_eco != 'NA':
            return 0 if s_eco == d_eco else 1
        return 1
    except:
        return -1

def establishment(ports, source, dest):
    try:
        sal_diff  = abs(float(ports[source]['Salinity'])  - float(ports[dest]['Salinity']))
        temp_diff = abs(float(ports[source]['YR_MEAN_T']) - float(ports[dest]['YR_MEAN_T']))
        return 0.00015 * math.exp(-0.5 * ((temp_diff/2)**2 + (sal_diff/10)**2))
    except:
        return -1

def fouling_risk(ports, source, stay_d, distance, trip_duration,
                 antifouling_p, portToCountry, countryStayRatio):
    if STAY_ATTENUATION and source in portToCountry:
        ratio  = countryStayRatio.get(portToCountry[source]['Country_name'], 1)
        stay_d = stay_d * ratio
    v = distance / trip_duration
    if abs(float(ports[source]['LATITUDE_DECIMAL'])) < 0.35:
        f = (0.000000129*stay_d**3 - 0.000083165*stay_d**2
             + 0.01495187*stay_d) * antifouling_p * math.exp(-0.008*v)
    else:
        f = (0.0000000014*stay_d**3 - 0.000016566*stay_d**2
             + 0.00519377*stay_d) * antifouling_p * math.exp(-0.008*v)
    return f

def adjust_trade_prob(source, dest, vessel_type, orig_prob, t,
                      portToCountry, exportDict, importDict):
    if orig_prob == 0 or not TRADE_ADJUST:
        return orig_prob
    trade_dict = importDict if t == 'fouling' else exportDict
    src_c = portToCountry[source]['Country_name']
    dst_c = portToCountry[dest]['Country_name']
    if src_c not in trade_dict:
        src_c = 'Rest of the World'
    elif dst_c not in trade_dict:
        dst_c = 'Rest of the World'
    try:
        ratio = Decimal(trade_dict[src_c][dst_c].get(vessel_type, 1))
    except:
        ratio = 1
    orig_prob = Decimal(orig_prob)
    return float(1 - (1 - orig_prob) ** ratio)

def is_valid_month(pair):
    first, second       = pair
    first_sail_year     = first['SAIL DATE'].split(' ')[0][0:4]
    second_arrival_year = second['ARRIVAL DATE'].split(' ')[0][0:4]
    if (first['VESSEL ID'] == second['VESSEL ID']
            and first_sail_year != '' and second_arrival_year != ''):
        if first['PLACE ID'] != second['PLACE ID']:
            valid_cross = (first_sail_year in
                           ('1997','1999','2002','2005','2008','2012','2015','2018')
                           and int(second_arrival_year) == int(first_sail_year) + 1)
            if first_sail_year == second_arrival_year or valid_cross:
                try:
                    trip_d = float(second['DURATION'])
                    stay_d = float(second['STAY DURATION'])
                    if trip_d > 0 and stay_d >= 0:
                        return True, trip_d, stay_d
                except:
                    return False, -1, -1
    return False, -1, -1

# ============================================================
# HON FUNCTIONS  (from 03, logic unchanged — no file I/O)
# ============================================================

def _kld(d1, d2):
    keys   = set(d1) | set(d2)
    result = 0.0
    for k in keys:
        p = d1.get(k, 1e-10)
        q = d2.get(k, 1e-10)
        if p > 0:
            result += p * math.log(p / q)
    return abs(result)


def build_hon_in_memory(trajectory_list, max_order=ORDER,
                        min_support=MIN_SUPPORT):
    import gc

    Distribution = defaultdict(lambda: defaultdict(list))
    for traj, prob in trajectory_list:
        target = traj[-1]
        source = tuple(traj[:-1])
        Distribution[source][target].append(prob)

    final_dist = defaultdict(dict)
    for source, targets in Distribution.items():
        for target, vals in targets.items():
            fp = float(np.mean(vals))
            if fp > min_support:
                final_dist[source][target] = fp

    del Distribution
    gc.collect()

    SourceToExtSource = defaultdict(lambda: defaultdict(set))
    for source in final_dist:
        if len(source) > 1:
            new_order = len(source)
            for start in range(1, len(source)):
                curr = source[start:]
                SourceToExtSource[curr][new_order].add(source)

    def kld_threshold(new_order, ext_source):
        return new_order / math.exp(1 + sum(final_dist[ext_source].values()))

    Rules = defaultdict(dict)

    def add_to_rules(source):
        if len(source) > 0:
            Rules[source] = final_dist[source]
            add_to_rules(source[:-1])

    def extend_rule(valid, curr, order):
        if order >= max_order:
            add_to_rules(valid)
            return
        distr = final_dist[valid]
        new_order = order + 1
        extended = SourceToExtSource[curr].get(new_order, [])
        if not extended:
            add_to_rules(valid)
            return
        for ext_source in extended:
            if _kld(final_dist[ext_source], distr) > kld_threshold(new_order,
                                                                   ext_source):
                extend_rule(ext_source, ext_source, new_order)
            else:
                extend_rule(valid, ext_source, new_order)

    for source in list(final_dist.keys()):
        if len(source) == 1:
            add_to_rules(source)
            extend_rule(source, source, 1)

    Graph = defaultdict(dict)
    for source in sorted(Rules, key=lambda x: len(x)):
        for target in Rules[source]:
            Graph[source][(target,)] = Rules[source][target]
            if len(source) > 1:
                prev_source = source[:-1]
                prev_target = (source[-1],)
                if prev_source in Graph and source not in Graph[prev_source]:
                    try:
                        Graph[prev_source][source] = Graph[prev_source][
                            prev_target]
                        del Graph[prev_source][prev_target]
                    except:
                        pass

    to_add, to_remove = [], []
    for source in Graph:
        for target in list(Graph[source]):
            if len(target) == 1:
                new_target = source + target
                while len(new_target) > 1:
                    if new_target in Graph:
                        to_add.append(
                            (source, new_target, Graph[source][target]))
                        to_remove.append((source, target))
                        break
                    new_target = new_target[1:]
    for (s, t, w) in to_add:
        Graph[s][t] = w
    for (s, t) in to_remove:
        Graph[s].pop(t, None)

    return Graph

# ============================================================
# RISK AGGREGATION  (from 04, logic unchanged)
# ============================================================

def compute_port_agg_risks(Graph):
    phys_net = defaultdict(list)
    for source_hon, targets in Graph.items():
        src_phys = source_hon[-1]
        for target_hon, weight in targets.items():
            dst_phys = target_hon[-1]
            if src_phys != dst_phys:
                phys_net[(src_phys, dst_phys)].append(weight)

    dest_risks = defaultdict(list)
    for (src, dst), weights in phys_net.items():
        dest_risks[dst].append(float(np.mean(weights)))

    return {port: float(np.mean(risks)) for port, risks in dest_risks.items()}

def agg_to_country(port_risks, portToCountry):
    country_risks = defaultdict(float)
    for port, risk in port_risks.items():
        if port in portToCountry:
            country_risks[portToCountry[port]['Country_name']] += risk
    return dict(country_risks)


def build_country_risks_from_items(items, portToCountry):
    if not items:
        return {}
    max_p = max(p for p, _ in items)
    if max_p == 0:
        return {}

    def traj_generator():
        for prob, subseqs in items:
            norm = prob / max_p
            if norm > 0:
                for subseq in subseqs:
                    yield (subseq, norm)

    graph = build_hon_in_memory(traj_generator())
    port_risks = compute_port_agg_risks(graph)
    return agg_to_country(port_risks, portToCountry)

# ============================================================
# PRE-COMPUTE FOULING ITEMS  (called once, result reused)
# ============================================================

def compute_fouling_items(probs_raw, ports, portToCountry,
                          countryStayRatio, exportDict, importDict):
    """
    Fouling does not use ballast discharge, so it never changes
    across perturbation runs. Compute once and reuse.
    Returns pre-normalised fouling country risks.
    """
    items = []
    for pair, move_list in probs_raw.items():
        source, dest = pair
        for move in move_list:
            if move['_p_estab'] == -1:
                continue
            p_f    = fouling_risk(ports, source,
                                  move['stay_duration'], move['distance'],
                                  move['trip_duration'], move['antifouling_p'],
                                  portToCountry, countryStayRatio)
            prob_f = move['_p_alien'] * p_f * move['_p_estab']
            prob_f = adjust_trade_prob(source, dest, move['vessel_type'],
                                       prob_f, 'fouling',
                                       portToCountry, exportDict, importDict)
            items.append((prob_f, move['subseq']))
    return build_country_risks_from_items(items, portToCountry)

# ============================================================
# SINGLE BALLAST SIMULATION  (called in parallel workers)
# ============================================================

def _run_one_ballast(args):
    """
    Worker function for multiprocessing.
    Each call runs one Monte Carlo iteration for ballast only.
    Returns dict {country: risk}.
    """
    sigma, seed = args
    rng   = np.random.RandomState(seed)
    items = []

    for pair, move_list in global_probs_raw.items():
        source, dest = pair
        for move in move_list:
            if move['_p_estab'] == -1:
                continue
            bw = move['ballast_discharge']
            if sigma > 0 and not move['is_us_source']:
                noise = rng.normal(1.0, sigma)
                noise = max(noise, 0.0)
                bw    = bw * noise
            bd      = 1 - math.exp(-3.22e-6 * bw)
            p_intro = (1 - B_EFFICACY[R]) * bd * math.exp(-0.02 * move['trip_duration'])
            prob_b  = move['_p_alien'] * p_intro * move['_p_estab']
            prob_b  = adjust_trade_prob(source, dest, move['vessel_type'],
                                        prob_b, 'ballast',
                                        global_portToCountry, global_exportDict, global_importDict)
            items.append((prob_b, move['subseq']))

    return build_country_risks_from_items(items, global_portToCountry)

# ============================================================
# STABILITY METRICS
# ============================================================

def make_baseline_info(country_risks):
    countries = sorted(country_risks.keys())
    n         = len(countries)
    values    = np.array([country_risks.get(c, 0.0) for c in countries])
    order     = values.argsort()[::-1]
    top_n     = set(np.array(countries)[order[:TOP_N]])
    topQ      = set(np.array(countries)[order[:n//4]])
    return countries, n, values, top_n, topQ

def compute_metrics(baseline_values, sim_values, all_countries,
                    top_n_baseline, topQ_baseline, n_countries):
    corr, _   = spearmanr(baseline_values, sim_values)
    sim_order  = sim_values.argsort()[::-1]
    countries  = np.array(all_countries)
    top_n_sim  = set(countries[sim_order[:TOP_N]])
    topQ_sim   = set(countries[sim_order[:n_countries//4]])
    top_n_pct  = len(top_n_baseline & top_n_sim) / TOP_N * 100
    topQ_pct   = len(topQ_baseline  & topQ_sim)  / len(topQ_baseline) * 100
    return corr, top_n_pct, topQ_pct

# ============================================================
# MAIN
# ============================================================

def main():
    global global_probs_raw, global_portToCountry, global_exportDict, global_importDict

    print("=" * 60)
    print("Sensitivity Analysis — Ballast Water Discharge Uncertainty")
    print("  Fouling computed once (no perturbation)")
    print("  Ballast Monte Carlo parallelised across CPU cores")
    print("=" * 60)
    n_workers = N_WORKERS or mp.cpu_count()
    print(f"sigma values:     {SIGMA_VALUES}")
    print(f"N simulations:    {N_SIMULATIONS}")
    print(f"CPU workers:      {n_workers}")
    print(f"stayAttenuation:  {STAY_ATTENUATION}")
    print(f"tradeAdjust:      {TRADE_ADJUST}")
    print()

    # ---- Load static data (once) ----
    print("Loading port data...")
    ports = get_port_data(port_data_file, 'ID', ',')

    countryStayRatio = {}
    with open(stay_file) as f:
        for row in csv.DictReader(f):
            countryStayRatio[row['Country']] = float(row['coefficient'])

    portToCountry = {}
    with open(port_to_country_file) as f:
        for row in csv.DictReader(f):
            portToCountry[row['Port_id']] = row

    us_ports = {pid for pid, info in portToCountry.items()
                if info.get('Country_name') == 'U.S.A.'}
    print(f"  US ports identified: {len(us_ports)}")

    # Use regular dicts (not defaultdict with lambda) so they can be
    # pickled and sent to multiprocessing worker processes on Windows.
    exportDict = {}
    importDict = {}
    with open(export_file) as f:
        for row in csv.DictReader(f):
            v = float(row['coefficient']) if row['coefficient'] != 'NA' else 1.0
            exportDict.setdefault(row['O_Countryname'], {}).setdefault(
                row['D_Countryname'], {})[row['VESSEL.TYPE']] = v
    with open(import_file) as f:
        for row in csv.DictReader(f):
            v = float(row['coefficient']) if row['coefficient'] != 'NA' else 1.0
            importDict.setdefault(row['O_Countryname'], {}).setdefault(
                row['D_Countryname'], {})[row['VESSEL.TYPE']] = v

    # ---- Load moves & build ProbsRaw (once) ----
    print("Loading moves...")
    moves       = []
    VALID_TYPES = {'Auto','Container','Bulk','Tanker','Chemical',
                   'Oil','General','Liquified-Gas','Refrigerated-Cargo'}
    with open(clean_move) as f:
        reader = csv.DictReader(f, delimiter='|')
        for row in reader:
            if (float(row['STAY DURATION']) < 20.0
                    and row['VESSEL TYPE'] in VALID_TYPES):
                parts = row['ROUT'].split('-')
                if len(parts) != 2:
                    continue
                orig, dest = parts
                if orig not in portToCountry or dest not in portToCountry:
                    continue
                if portToCountry[orig]['Country_id'] == portToCountry[dest]['Country_id']:
                    continue
                moves.append(row)
    print(f"  Moves loaded: {len(moves)}")

    print("Building ProbsRaw...")
    ship_id   = -1
    ship_traj = []
    probs_raw = {}

    for prev_idx in range(len(moves) - 1):
        prev_move = moves[prev_idx]
        next_move = moves[prev_idx + 1]
        source    = prev_move['PLACE ID']
        dest      = next_move['PLACE ID']
        pair      = (source, dest)
        is_valid, trip_dur, stay_dur = is_valid_month((prev_move, next_move))

        if is_valid:
            if prev_move['VESSEL ID'] != ship_id:
                ship_traj = [source, dest]
            else:
                ship_traj.append(dest)
            ship_id = prev_move['VESSEL ID']
            try:
                distance = get_distance_km(ports, source, dest)
                subseq   = [ship_traj[-lastn:]
                             for lastn in range(2, 1 + min(len(ship_traj), ORDER))]
                if pair not in probs_raw:
                    probs_raw[pair] = []
                probs_raw[pair].append({
                    'trip_duration':    float(trip_dur),
                    'stay_duration':    float(stay_dur),
                    'ballast_discharge':float(prev_move['BALLAST DISCHARGE']),
                    'antifouling_p':    float(prev_move['ANTIFOULING PROB']),
                    'distance':         float(distance),
                    'subseq':           subseq,
                    'vessel_type':      prev_move['VESSEL TYPE'],
                    'is_us_source':     source in us_ports,
                })
            except:
                pass
        else:
            ship_id = -1

    print(f"  ProbsRaw built: {len(probs_raw)} port pairs")

    # ---- Precompute static p_alien and p_estab ----
    print("Precomputing static probabilities (p_alien, p_estab)...")
    for pair, move_list in probs_raw.items():
        source, dest = pair
        p_alien = (1 if ECO == 'noEco_'
                   else filterbyeco_same_only(ports, source, dest) if ECO == 'sameEco_'
                   else filterbyeco(ports, source, dest))
        p_estab = 1 if ENV == 'noEnv_' else establishment(ports, source, dest)
        for move in move_list:
            move['_p_alien'] = p_alien
            move['_p_estab'] = p_estab
    print("  Done.")

    # Assign completed dictionaries to global variables for WSL fork sharing
    global_probs_raw     = probs_raw
    global_portToCountry = portToCountry
    global_exportDict    = exportDict
    global_importDict    = importDict

    # ---- Compute fouling ONCE (reused for all sigmas) ----
    if RUN_FOULING:
        print("\nComputing fouling risks (once, no perturbation)...")
        fouling_country_risks = compute_fouling_items(
            probs_raw, ports, portToCountry, countryStayRatio, exportDict, importDict)
        f_countries, f_n, f_values, f_top_n, f_topQ = make_baseline_info(fouling_country_risks)
        print(f"  Fouling — countries: {f_n}, top-{TOP_N}: {sorted(f_top_n)}")
    else:
        print("\nSkipping fouling (RUN_FOULING=False).")

    # ---- Baseline ballast run (sigma=0) ----
    if RUN_BALLAST:
        print("\nComputing baseline ballast (sigma=0)...")
        baseline_args = (0.0, SEED)
        base_b        = _run_one_ballast(baseline_args)
        b_countries, b_n, b_values, b_top_n, b_topQ = make_baseline_info(base_b)
        print(f"  Ballast — countries: {b_n}, top-{TOP_N}: {sorted(b_top_n)}")
    else:
        print("\nSkipping ballast (RUN_BALLAST=False).")

    # ---- Monte Carlo ----
    all_results = {}

    for sigma in SIGMA_VALUES:
        all_results[sigma] = {}

        # -- Ballast Monte Carlo (parallelised) --
        if RUN_BALLAST:
            print(f"\nMonte Carlo ballast: sigma={sigma}, N={N_SIMULATIONS}, "
                  f"workers={n_workers} ...")
            worker_args = [
                (sigma, SEED + i)
                for i in range(N_SIMULATIONS)
            ]

            gc.disable()
            gc.freeze()

            with mp.Pool(processes=n_workers) as pool:
                sim_results = pool.map(_run_one_ballast, worker_args)

            b_corr, b_top_n_pct, b_topQ_pct = [], [], []
            for sim_b in sim_results:
                sim_b_vals = np.array([sim_b.get(c, 0.0) for c in b_countries])
                c, tn, tq  = compute_metrics(b_values, sim_b_vals, b_countries,
                                             b_top_n, b_topQ, b_n)
                b_corr.append(c)
                b_top_n_pct.append(tn)
                b_topQ_pct.append(tq)

            all_results[sigma]['ballast'] = {
                'corr_mean':  np.mean(b_corr),
                'corr_std':   np.std(b_corr),
                'corr_min':   np.min(b_corr),
                'top_n_mean': np.mean(b_top_n_pct),
                'top_n_min':  np.min(b_top_n_pct),
                'topQ_mean':  np.mean(b_topQ_pct),
                'topQ_min':   np.min(b_topQ_pct),
            }

        # -- Fouling (computed once, identical across all runs) --
        if RUN_FOULING:
            sim_f_vals = f_values  # same as baseline — perturbation has no effect
            f_corr_val, f_tn, f_tq = compute_metrics(
                f_values, sim_f_vals, f_countries, f_top_n, f_topQ, f_n)
            all_results[sigma]['fouling'] = {
                'corr_mean':  f_corr_val,
                'corr_std':   0.0,
                'corr_min':   f_corr_val,
                'top_n_mean': f_tn,
                'top_n_min':  f_tn,
                'topQ_mean':  f_tq,
                'topQ_min':   f_tq,
                'note': ('Fouling risk does not depend on ballast discharge '
                         'estimates, so rankings are invariant to perturbation.'),
            }

    # ---- Print results ----
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    for sigma, res in all_results.items():
        print(f"\n  sigma = {sigma}  (+-{int(sigma*100)}% perturbation on non-US ballast discharge)")
        for label, key in [('Ballast', 'ballast'), ('Fouling', 'fouling')]:
            if key not in res:
                continue
            r = res[key]
            print(f"\n  [{label}]")
            if 'note' in r:
                print(f"    Note: {r['note']}")
            print(f"    Spearman:          {r['corr_mean']:.4f} +- {r['corr_std']:.4f}"
                  f"  (min = {r['corr_min']:.4f})")
            print(f"    Top-{TOP_N} stability: {r['top_n_mean']:.1f}% mean"
                  f"  (min = {r['top_n_min']:.1f}%)")
            print(f"    Top-Q stability:   {r['topQ_mean']:.1f}% mean"
                  f"  (min = {r['topQ_min']:.1f}%)")

if __name__ == '__main__':
    # Required for multiprocessing on Windows
    mp.freeze_support()
    main()