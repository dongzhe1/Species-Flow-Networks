"""
Sensitivity Analysis — Biofouling model parameters
====================================================
Mirrors the structure of bw_sensitivity_analysis.py but targets the
continuous parameters of the biofouling P(NIS) model rather than
ballast discharge volumes.

Two analyses are run:

1. MONTE CARLO PARAMETER PERTURBATION
   Each continuous parameter is independently perturbed each iteration:
       param_i ~ param_i_base * max(Normal(1.0, sigma), 0)
   Parameters perturbed:
       P(intro)  : lambda_, aT1, aT2, aT3, aE1, aE2, aE3
       P(estab)  : beta, delta_T, delta_S

2. LATITUDE THRESHOLD SCENARIO ANALYSIS
   Fixed parameters, tropical/temperate boundary set to:
       baseline  : 0.35 rad  (~20 deg)
       scenario A: 15 deg converted to radians
       scenario B: 25 deg converted to radians
   Each scenario is a single deterministic run (no Monte Carlo).

CONFIG section below controls all key settings.
"""

import csv
import math
import sys
import numpy as np
from collections import defaultdict
from decimal import Decimal, getcontext
from scipy.stats import spearmanr

sys.setrecursionlimit(5000)
getcontext().prec = 50

# ============================================================
# CONFIG — edit these parameters
# ============================================================
SIGMA_VALUES     = [0.2, 0.4]    # perturbation magnitudes to test
N_SIMULATIONS    = 100           # Monte Carlo iterations per sigma
STAY_ATTENUATION = False          # 'stayAttenuation' in original code
TRADE_ADJUST     = False          # 'tradeAdjust' in original code
TOP_N            = 10             # top-N countries for stability metric
SEED             = 42

# Latitude threshold scenarios (degrees -> radians)
# Baseline is 0.35 rad; scenarios A and B override this value
LAT_THRESHOLD_BASE = 0.35                      # ~20 deg, original value
LAT_THRESHOLD_15   = 15.0 * math.pi / 180.0   # 0.2618 rad
LAT_THRESHOLD_25   = 25.0 * math.pi / 180.0   # 0.4363 rad

# Fixed baseline model parameters
BASE_LAMBDA = 0.008          # speed decay rate in P(intro)
BASE_AT1    = 1.29e-7        # tropical cubic coefficient
BASE_AT2    = 8.316e-5       # tropical quadratic coefficient
BASE_AT3    = 0.01495187     # tropical linear coefficient
BASE_AE1    = 1.4e-9         # temperate cubic coefficient
BASE_AE2    = 1.6566e-5      # temperate quadratic coefficient
BASE_AE3    = 5.19377e-3     # temperate linear coefficient
BASE_BETA   = 0.00015        # establishment base probability
BASE_DELTA_T = 2.0           # temperature tolerance std (deg C)
BASE_DELTA_S = 10.0          # salinity tolerance std (ppt)

# Fixed structural settings
ORDER       = 16
R           = 'r0'
ENV         = 'env_'
ECO         = 'Eco_'
MIN_SUPPORT = 1e-20

# File paths — adjust to your environment
BASE_DIR             = '../../data/'
clean_move           = BASE_DIR + 'moves/moves_cleaned_2018.txt'
port_data_file       = BASE_DIR + 'Places_allportdata_mergedSept2017.csv'
stay_file            = BASE_DIR + 'stay/stay.csv'
port_to_country_file = BASE_DIR + 'stay/portToCountry.csv'
export_file          = BASE_DIR + 'trade/Export_change_for_BW.csv'
import_file          = BASE_DIR + 'trade/Import_change_for_fouling.csv'

# ============================================================
# UTILITY FUNCTIONS  (from 02_Prepare_HON_Input, logic unchanged)
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
    a    = (math.sin(dLat / 2) ** 2
            + math.cos(deg2rad(lat1)) * math.cos(deg2rad(lat2))
            * math.sin(dLon / 2) ** 2)
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

def establishment(ports, source, dest,
                  beta=BASE_BETA,
                  delta_T=BASE_DELTA_T,
                  delta_S=BASE_DELTA_S):
    """P(estab) with optionally perturbed parameters."""
    try:
        sal_diff  = abs(float(ports[source]['Salinity'])
                        - float(ports[dest]['Salinity']))
        temp_diff = abs(float(ports[source]['YR_MEAN_T'])
                        - float(ports[dest]['YR_MEAN_T']))
        return beta * math.exp(
            -0.5 * ((temp_diff / delta_T) ** 2 + (sal_diff / delta_S) ** 2)
        )
    except:
        return -1

def fouling_risk(ports, source, stay_d, distance, trip_duration,
                 antifouling_p, portToCountry, countryStayRatio,
                 lambda_=BASE_LAMBDA,
                 aT1=BASE_AT1, aT2=BASE_AT2, aT3=BASE_AT3,
                 aE1=BASE_AE1, aE2=BASE_AE2, aE3=BASE_AE3,
                 lat_threshold=LAT_THRESHOLD_BASE):
    """P(intro) with optionally perturbed parameters."""
    if STAY_ATTENUATION and source in portToCountry:
        ratio  = countryStayRatio.get(portToCountry[source]['Country_name'], 1)
        stay_d = stay_d * ratio
    v = distance / trip_duration
    if abs(float(ports[source]['LATITUDE_DECIMAL'])) < lat_threshold:
        # tropical
        f = ((aT1 * stay_d ** 3 - aT2 * stay_d ** 2 + aT3 * stay_d)
             * antifouling_p * math.exp(-lambda_ * v))
    else:
        # temperate
        f = ((aE1 * stay_d ** 3 - aE2 * stay_d ** 2 + aE3 * stay_d)
             * antifouling_p * math.exp(-lambda_ * v))
    return f

def adjust_trade_prob(source, dest, vessel_type, orig_prob,
                      portToCountry, exportDict, importDict):
    if orig_prob == 0 or not TRADE_ADJUST:
        return orig_prob
    src_c = portToCountry[source]['Country_name']
    dst_c = portToCountry[dest]['Country_name']
    if src_c not in importDict:
        src_c = 'Rest of the World'
    elif dst_c not in importDict:
        dst_c = 'Rest of the World'
    try:
        ratio = Decimal(importDict[src_c][dst_c].get(vessel_type, 1))
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
                           ('1997', '1999', '2002', '2005',
                            '2008', '2012', '2015', '2018')
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
# HON FUNCTIONS  (from 03_Build-SF-HON, logic unchanged)
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

def build_hon_in_memory(trajectory_list, max_order=ORDER, min_support=MIN_SUPPORT):
    # Build_Observations_Distributions
    Distribution = defaultdict(lambda: defaultdict(list))
    for traj, prob in trajectory_list:
        target = traj[-1]
        source = tuple(traj[:-1])
        Distribution[source][target].append(prob)

    # Aggregate_Probs
    final_dist = defaultdict(dict)
    for source, targets in Distribution.items():
        for target, vals in targets.items():
            fp = float(np.mean(vals))
            if fp > min_support:
                final_dist[source][target] = fp

    # SourceToExtSource cache
    SourceToExtSource = defaultdict(lambda: defaultdict(set))
    for source in final_dist:
        if len(source) > 1:
            new_order = len(source)
            for start in range(1, len(source)):
                curr = source[start:]
                SourceToExtSource[curr][new_order].add(source)

    def kld_threshold(new_order, ext_source):
        return new_order / math.exp(1 + sum(final_dist[ext_source].values()))

    # GenerateAllRules
    Rules = defaultdict(dict)

    def add_to_rules(source):
        if len(source) > 0:
            Rules[source] = final_dist[source]
            add_to_rules(source[:-1])

    def extend_rule(valid, curr, order):
        if order >= max_order:
            add_to_rules(valid)
            return
        distr     = final_dist[valid]
        new_order = order + 1
        extended  = SourceToExtSource[curr].get(new_order, [])
        if not extended:
            add_to_rules(valid)
            return
        for ext_source in extended:
            if _kld(final_dist[ext_source], distr) > kld_threshold(new_order, ext_source):
                extend_rule(ext_source, ext_source, new_order)
            else:
                extend_rule(valid, ext_source, new_order)

    for source in list(final_dist.keys()):
        if len(source) == 1:
            add_to_rules(source)
            extend_rule(source, source, 1)

    # BuildNetwork
    Graph = defaultdict(dict)
    for source in sorted(Rules, key=lambda x: len(x)):
        for target in Rules[source]:
            Graph[source][(target,)] = Rules[source][target]
            if len(source) > 1:
                prev_source = source[:-1]
                prev_target = (source[-1],)
                if prev_source in Graph and source not in Graph[prev_source]:
                    try:
                        Graph[prev_source][source] = Graph[prev_source][prev_target]
                        del Graph[prev_source][prev_target]
                    except:
                        pass

    # RewireTails
    to_add, to_remove = [], []
    for source in Graph:
        for target in list(Graph[source]):
            if len(target) == 1:
                new_target = source + target
                while len(new_target) > 1:
                    if new_target in Graph:
                        to_add.append((source, new_target, Graph[source][target]))
                        to_remove.append((source, target))
                        break
                    new_target = new_target[1:]
    for (s, t, w) in to_add:
        Graph[s][t] = w
    for (s, t) in to_remove:
        Graph[s].pop(t, None)

    return Graph

# ============================================================
# RISK AGGREGATION  (from 04_make_HONet_dict, logic unchanged)
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
    traj_list = []
    for prob, subseqs in items:
        norm = prob / max_p
        if norm > 0:
            for subseq in subseqs:
                traj_list.append((subseq, norm))
    graph      = build_hon_in_memory(traj_list)
    port_risks = compute_port_agg_risks(graph)
    return agg_to_country(port_risks, portToCountry)

# ============================================================
# PIPELINE
# ============================================================

def run_pipeline(probs_raw, ports, portToCountry, countryStayRatio,
                 exportDict, importDict,
                 rng=None, sigma=0.0,
                 lat_threshold=LAT_THRESHOLD_BASE):
    """
    One full fouling pipeline iteration.

    If sigma > 0 and rng is provided, each continuous parameter is
    independently perturbed:
        param_i = param_i_base * max(Normal(1.0, sigma), 0)

    If sigma == 0, baseline parameters are used (deterministic).
    lat_threshold controls the tropical/temperate boundary.

    Returns: country_risks dict
    """
    # Sample perturbed parameters (or use baseline)
    def _p(base):
        if sigma > 0 and rng is not None:
            return base * max(float(rng.normal(1.0, sigma)), 0.0)
        return base

    lambda_ = _p(BASE_LAMBDA)
    aT1     = _p(BASE_AT1)
    aT2     = _p(BASE_AT2)
    aT3     = _p(BASE_AT3)
    aE1     = _p(BASE_AE1)
    aE2     = _p(BASE_AE2)
    aE3     = _p(BASE_AE3)
    beta    = _p(BASE_BETA)
    delta_T = _p(BASE_DELTA_T)
    delta_S = _p(BASE_DELTA_S)

    fouling_items = []

    for pair, move_list in probs_raw.items():
        source, dest = pair

        p_alien = move_list[0]['_p_alien']  # precomputed, structural — not perturbed

        # P(estab) — re-computed with perturbed beta / delta_T / delta_S
        p_estab = establishment(ports, source, dest,
                                beta=beta,
                                delta_T=delta_T,
                                delta_S=delta_S)
        if p_estab == -1:
            continue

        for move in move_list:
            p_f = fouling_risk(
                ports, source,
                move['stay_duration'],
                move['distance'],
                move['trip_duration'],
                move['antifouling_p'],
                portToCountry, countryStayRatio,
                lambda_=lambda_,
                aT1=aT1, aT2=aT2, aT3=aT3,
                aE1=aE1, aE2=aE2, aE3=aE3,
                lat_threshold=lat_threshold,
            )
            prob_f = p_alien * p_f * p_estab
            prob_f = adjust_trade_prob(source, dest, move['vessel_type'],
                                       prob_f, portToCountry,
                                       exportDict, importDict)
            fouling_items.append((prob_f, move['subseq']))

    return build_country_risks_from_items(fouling_items, portToCountry)

# ============================================================
# STABILITY METRICS  (same as bw_sensitivity_analysis.py)
# ============================================================

def make_baseline_info(country_risks):
    countries = sorted(country_risks.keys())
    n         = len(countries)
    values    = np.array([country_risks.get(c, 0.0) for c in countries])
    order     = values.argsort()[::-1]
    top_n     = set(np.array(countries)[order[:TOP_N]])
    topQ      = set(np.array(countries)[order[:n // 4]])
    return countries, n, values, top_n, topQ

def compute_metrics(baseline_values, sim_values, all_countries,
                    top_n_baseline, topQ_baseline, n_countries):
    corr, _   = spearmanr(baseline_values, sim_values)
    sim_order  = sim_values.argsort()[::-1]
    countries  = np.array(all_countries)
    top_n_sim  = set(countries[sim_order[:TOP_N]])
    topQ_sim   = set(countries[sim_order[:n_countries // 4]])
    top_n_pct  = len(top_n_baseline & top_n_sim) / TOP_N * 100
    topQ_pct   = len(topQ_baseline  & topQ_sim)  / len(topQ_baseline) * 100
    return corr, top_n_pct, topQ_pct

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("Sensitivity Analysis — Biofouling model parameters")
    print("=" * 60)
    print(f"sigma values:     {SIGMA_VALUES}")
    print(f"N simulations:    {N_SIMULATIONS}")
    print(f"stayAttenuation:  {STAY_ATTENUATION}")
    print(f"tradeAdjust:      {TRADE_ADJUST}")
    print(f"Lat threshold (baseline): {LAT_THRESHOLD_BASE:.4f} rad"
          f"  (~{math.degrees(LAT_THRESHOLD_BASE):.1f} deg)")
    print(f"Lat threshold (15 deg):   {LAT_THRESHOLD_15:.4f} rad")
    print(f"Lat threshold (25 deg):   {LAT_THRESHOLD_25:.4f} rad")
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
    VALID_TYPES = {'Auto', 'Container', 'Bulk', 'Tanker', 'Chemical',
                   'Oil', 'General', 'Liquified-Gas', 'Refrigerated-Cargo'}
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
                    'trip_duration':  float(trip_dur),
                    'stay_duration':  float(stay_dur),
                    'antifouling_p':  float(prev_move['ANTIFOULING PROB']),
                    'distance':       float(distance),
                    'subseq':         subseq,
                    'vessel_type':    prev_move['VESSEL TYPE'],
                })
            except:
                pass
        else:
            ship_id = -1

    print(f"  ProbsRaw built: {len(probs_raw)} port pairs")

    # ---- Precompute p_alien (structural, never perturbed) ----
    print("Precomputing p_alien (structural)...")
    for pair, move_list in probs_raw.items():
        source, dest = pair
        p_alien = (1 if ECO == 'noEco_'
                   else filterbyeco_same_only(ports, source, dest) if ECO == 'sameEco_'
                   else filterbyeco(ports, source, dest))
        for move in move_list:
            move['_p_alien'] = p_alien
    print("  Done.")

    # ============================================================
    # PART 1 — MONTE CARLO PARAMETER PERTURBATION
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 1 — Monte Carlo parameter perturbation")
    print("=" * 60)

    print("\nComputing baseline (sigma=0, lat_threshold=baseline)...")
    base_risks = run_pipeline(
        probs_raw, ports, portToCountry, countryStayRatio,
        exportDict, importDict,
        rng=None, sigma=0.0,
        lat_threshold=LAT_THRESHOLD_BASE,
    )
    countries, n_c, base_values, base_top_n, base_topQ = make_baseline_info(base_risks)
    print(f"  Countries: {n_c}, top-{TOP_N}: {sorted(base_top_n)}")

    mc_results = {}
    for sigma in SIGMA_VALUES:
        print(f"\nMonte Carlo: sigma={sigma}, N={N_SIMULATIONS} ...")
        rng = np.random.RandomState(SEED)
        corr_list, top_n_list, topQ_list = [], [], []

        for i in range(N_SIMULATIONS):
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{N_SIMULATIONS}")
            sim_risks = run_pipeline(
                probs_raw, ports, portToCountry, countryStayRatio,
                exportDict, importDict,
                rng=rng, sigma=sigma,
                lat_threshold=LAT_THRESHOLD_BASE,
            )
            sim_vals = np.array([sim_risks.get(c, 0.0) for c in countries])
            corr, tn, tq = compute_metrics(base_values, sim_vals, countries,
                                           base_top_n, base_topQ, n_c)
            corr_list.append(corr)
            top_n_list.append(tn)
            topQ_list.append(tq)

        mc_results[sigma] = {
            'corr_mean':  np.mean(corr_list),
            'corr_std':   np.std(corr_list),
            'corr_min':   np.min(corr_list),
            'top_n_mean': np.mean(top_n_list),
            'top_n_min':  np.min(top_n_list),
            'topQ_mean':  np.mean(topQ_list),
            'topQ_min':   np.min(topQ_list),
        }

    # ============================================================
    # PART 2 — LATITUDE THRESHOLD SCENARIO ANALYSIS
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 2 — Latitude threshold scenario analysis")
    print("(fixed parameters, baseline sigma=0)")
    print("=" * 60)

    scenario_results = {}
    for label, threshold in [
        ('15 deg', LAT_THRESHOLD_15),
        ('25 deg', LAT_THRESHOLD_25),
    ]:
        print(f"\nScenario: lat_threshold = {label} ({threshold:.4f} rad)...")
        scen_risks = run_pipeline(
            probs_raw, ports, portToCountry, countryStayRatio,
            exportDict, importDict,
            rng=None, sigma=0.0,
            lat_threshold=threshold,
        )
        scen_vals = np.array([scen_risks.get(c, 0.0) for c in countries])
        corr, tn, tq = compute_metrics(base_values, scen_vals, countries,
                                       base_top_n, base_topQ, n_c)
        scen_order   = scen_vals.argsort()[::-1]
        scen_top_n   = list(np.array(countries)[scen_order[:TOP_N]])
        scenario_results[label] = {
            'corr':       corr,
            'top_n_pct':  tn,
            'topQ_pct':   tq,
            'top_n_list': scen_top_n,
        }

    # ============================================================
    # PRINT RESULTS
    # ============================================================
    print("\n" + "=" * 60)
    print("RESULTS — PART 1: Monte Carlo parameter perturbation")
    print("=" * 60)
    print(f"  Baseline top-{TOP_N}: {sorted(base_top_n)}\n")
    for sigma, r in mc_results.items():
        print(f"  sigma = {sigma}  (+-{int(sigma*100)}% independent perturbation"
              f" on each continuous parameter)")
        print(f"    Spearman:          {r['corr_mean']:.4f} +- {r['corr_std']:.4f}"
              f"  (min = {r['corr_min']:.4f})")
        print(f"    Top-{TOP_N} stability: {r['top_n_mean']:.1f}% mean"
              f"  (min = {r['top_n_min']:.1f}%)")
        print(f"    Top-Q stability:   {r['topQ_mean']:.1f}% mean"
              f"  (min = {r['topQ_min']:.1f}%)")
        print()

    print("=" * 60)
    print("RESULTS — PART 2: Latitude threshold scenarios")
    print("=" * 60)
    print(f"  Baseline (~20 deg) top-{TOP_N}: {sorted(base_top_n)}\n")
    for label, r in scenario_results.items():
        print(f"  Scenario: {label}")
        print(f"    Spearman vs baseline:  {r['corr']:.4f}")
        print(f"    Top-{TOP_N} overlap:       {r['top_n_pct']:.1f}%")
        print(f"    Top-Q overlap:         {r['topQ_pct']:.1f}%")
        print(f"    Top-{TOP_N} countries:     {r['top_n_list']}")
        print()


if __name__ == '__main__':
    main()