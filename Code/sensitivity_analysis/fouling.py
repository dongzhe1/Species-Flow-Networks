"""
Sensitivity Analysis — Biofouling model parameters
====================================================
Mirrors the structure of ballast_water.py but targets the continuous
parameters of the biofouling P(NIS) model.

Two analyses are run, each evaluated at three levels:
  A) Country-level  : 32 BRI countries only
  B) Pair-level (global) : all country pairs
  C) Pair-level (BRI)    : pairs where BOTH countries are BRI members

Analysis 1 — MONTE CARLO PARAMETER PERTURBATION
  Each continuous parameter is independently perturbed each iteration:
      param_i ~ param_i_base * max(Normal(1.0, sigma), 0)
  Parameters perturbed:
      P(intro)  : lambda_, aT1, aT2, aT3, aE1, aE2, aE3
      P(estab)  : beta, delta_T, delta_S

Analysis 2 — LATITUDE THRESHOLD SCENARIO ANALYSIS
  Fixed parameters, tropical/temperate boundary set to:
      baseline  : 0.35 rad  (~20 deg)
      scenario A: 15 deg  (narrower tropical zone)
      scenario B: 25 deg  (wider tropical zone)
  Each scenario is a single deterministic run (no Monte Carlo).

Metrics:
  Country-level (A)       : Spearman, Top-5, Top-25%, Category agreement rate
  Pair-level global (B)   : Spearman, Top-20, Top-25%
  Pair-level BRI (C)      : Spearman, Top-10, Top-25%, Category agreement rate

Risk categories (from paper):
  very high : p(NIS) >= 0.75
  high      : 0.50 <= p(NIS) < 0.75
  medium    : 0.25 <= p(NIS) < 0.50
  low       : p(NIS) < 0.25

Note: pipeline values are HON-normalized, not raw P(NIS). Category
agreement is computed on these normalized values using the same thresholds.
The metric may be less informative if most values fall below 0.25.
"""

import csv
import math
import sys
import warnings
import numpy as np
from collections import defaultdict
from decimal import Decimal, getcontext
from scipy.stats import spearmanr

sys.setrecursionlimit(5000)
getcontext().prec = 50

# ============================================================
# CONFIG — edit these parameters
# ============================================================
SIGMA_VALUES      = [0.1, 0.3]   # perturbation magnitudes to test
N_SIMULATIONS     = 100           # Monte Carlo iterations per sigma
STAY_ATTENUATION  = False         # 'stayAttenuation' in original code
TRADE_ADJUST      = False         # 'tradeAdjust' in original code
TOP_N_COUNTRY     = 5            # top-N for country-level (small sample of 32)
TOP_N_PAIR_GLOBAL = 20           # top-N for global country-pair stability
TOP_N_PAIR_BRI    = 10           # top-N for BRI country-pair stability
SEED              = 42

# 32 BRI Maritime Silk Road countries, identified by Country_id.
# Country names are resolved at runtime from portToCountry to avoid
# hard-coding strings that may differ from the actual data.
BRI_COUNTRY_IDS = {
    'CHN', 'IDN', 'PHL', 'VNM', 'ARE', 'KWT', 'OMN', 'IND', 'SGP', 'BGD',
    'LKA', 'KOR', 'MAR', 'EGY', 'PAK', 'THA', 'ITA', 'TUR', 'SAU', 'BHR',
    'QAT', 'GRC', 'KEN', 'NZL', 'MYS', 'RUS', 'PRT', 'MLT', 'IRN', 'MMR',
    'TUN', 'KHM',
}

# Latitude threshold scenarios (degrees -> radians)
LAT_THRESHOLD_BASE = 0.35                      # ~20 deg, original value
LAT_THRESHOLD_15   = 15.0 * math.pi / 180.0   # narrower tropical zone
LAT_THRESHOLD_25   = 25.0 * math.pi / 180.0   # wider tropical zone
RUN_LAT = False

# Fixed baseline model parameters
BASE_LAMBDA  = 0.008        # speed decay rate in P(intro)
BASE_AT1     = 1.29e-7      # tropical cubic coefficient
BASE_AT2     = 8.316e-5     # tropical quadratic coefficient
BASE_AT3     = 0.01495187   # tropical linear coefficient
BASE_AE1     = 1.4e-9       # temperate cubic coefficient
BASE_AE2     = 1.6566e-5    # temperate quadratic coefficient
BASE_AE3     = 5.19377e-3   # temperate linear coefficient
BASE_BETA    = 0.00015      # establishment base probability
BASE_DELTA_T = 2.0          # temperature tolerance std (deg C)
BASE_DELTA_S = 10.0         # salinity tolerance std (ppt)

# Fixed structural settings
ORDER       = 16
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
        f = ((aT1 * stay_d ** 3 - aT2 * stay_d ** 2 + aT3 * stay_d)
             * antifouling_p * math.exp(-lambda_ * v))
    else:
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
    """Destination port risk: mean of incoming HON edge weights."""
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

def compute_port_pair_risks(Graph):
    """Mean HON edge weight for each physical (src_port, dst_port) pair."""
    phys_net = defaultdict(list)
    for source_hon, targets in Graph.items():
        src_phys = source_hon[-1]
        for target_hon, weight in targets.items():
            dst_phys = target_hon[-1]
            if src_phys != dst_phys:
                phys_net[(src_phys, dst_phys)].append(weight)
    return {pair: float(np.mean(weights)) for pair, weights in phys_net.items()}

def agg_to_country(port_risks, portToCountry):
    country_risks = defaultdict(float)
    for port, risk in port_risks.items():
        if str(port) in portToCountry:
            country_risks[portToCountry[str(port)]['Country_name']] += risk
    return dict(country_risks)

def agg_to_country_pair(port_pair_risks, portToCountry):
    """Aggregate port-pair risks to (src_country, dst_country) pairs."""
    country_pair_risks = defaultdict(float)
    for (src_port, dst_port), risk in port_pair_risks.items():
        src_s = str(src_port)
        dst_s = str(dst_port)
        if src_s in portToCountry and dst_s in portToCountry:
            src_c = portToCountry[src_s]['Country_name']
            dst_c = portToCountry[dst_s]['Country_name']
            if src_c != dst_c:
                country_pair_risks[(src_c, dst_c)] += risk
    return dict(country_pair_risks)

def build_risks_from_items(items, portToCountry):
    """
    Build HON and return (country_risks, pair_risks).
    country_risks : {country: aggregated_risk}
    pair_risks    : {(src_country, dst_country): aggregated_risk}
    """
    if not items:
        return {}, {}
    max_p = max(p for p, _ in items)
    if max_p == 0:
        return {}, {}
    traj_list = []
    for prob, subseqs in items:
        norm = prob / max_p
        if norm > 0:
            for subseq in subseqs:
                traj_list.append((subseq, norm))
    graph           = build_hon_in_memory(traj_list)
    port_risks      = compute_port_agg_risks(graph)
    port_pair_risks = compute_port_pair_risks(graph)
    country_risks   = agg_to_country(port_risks, portToCountry)
    pair_risks      = agg_to_country_pair(port_pair_risks, portToCountry)
    return country_risks, pair_risks

# ============================================================
# PIPELINE
# ============================================================

def run_pipeline(probs_raw, ports, portToCountry, countryStayRatio,
                 exportDict, importDict,
                 rng=None, sigma=0.0,
                 lat_threshold=LAT_THRESHOLD_BASE):
    """
    One full fouling pipeline iteration.
    Returns (country_risks, pair_risks).
    sigma=0 -> baseline (deterministic). sigma>0 -> perturbed Monte Carlo.
    """
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
        p_alien = move_list[0]['_p_alien']
        p_estab = establishment(ports, source, dest,
                                beta=beta, delta_T=delta_T, delta_S=delta_S)
        if p_estab == -1:
            continue
        for move in move_list:
            p_f = fouling_risk(
                ports, source,
                move['stay_duration'], move['distance'], move['trip_duration'],
                move['antifouling_p'], portToCountry, countryStayRatio,
                lambda_=lambda_,
                aT1=aT1, aT2=aT2, aT3=aT3,
                aE1=aE1, aE2=aE2, aE3=aE3,
                lat_threshold=lat_threshold,
            )
            prob_f = p_alien * p_f * p_estab
            prob_f = adjust_trade_prob(source, dest, move['vessel_type'],
                                       prob_f, portToCountry, exportDict, importDict)
            fouling_items.append((prob_f, move['subseq']))

    return build_risks_from_items(fouling_items, portToCountry)

# ============================================================
# METRICS
# ============================================================

def spearman_safe(a, b):
    """Spearman correlation; returns nan without warning if input is constant."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        corr, _ = spearmanr(a, b)
    return float(corr)

def risk_category(value):
    """Paper risk classification (applied to normalized pipeline values)."""
    if value >= 0.75:   return 'very high'
    elif value >= 0.50: return 'high'
    elif value >= 0.25: return 'medium'
    else:               return 'low'

def category_agreement_rate(base_risks, sim_risks, entities):
    """
    Fraction (%) of entities whose risk category is unchanged after perturbation.
    entities: list of country names or (src, dst) tuples.
    """
    if not entities:
        return float('nan')
    n_agree = sum(
        1 for e in entities
        if risk_category(base_risks.get(e, 0.0)) == risk_category(sim_risks.get(e, 0.0))
    )
    return n_agree / len(entities) * 100

def make_baseline_info(risks, key_list, top_n):
    """
    key_list : sorted list of countries or pair-tuples to analyse.
    Returns  : (values_array, top_n_set, topQ_set)
    """
    values  = np.array([risks.get(k, 0.0) for k in key_list])
    n       = len(key_list)
    order   = values.argsort()[::-1]
    top_n_set = set(key_list[i] for i in order[:top_n])
    topQ_set  = set(key_list[i] for i in order[:max(1, n // 4)])
    return values, top_n_set, topQ_set

def compute_metrics(base_values, sim_risks, key_list,
                    base_top_n, base_topQ, top_n,
                    base_risks_for_cat=None, entities_for_cat=None):
    """
    General metric computation for both country-level and pair-level.
    If base_risks_for_cat and entities_for_cat are provided, also computes
    category agreement rate.
    Returns: (spearman, top_n_pct, topQ_pct, cat_agree_pct)
             cat_agree_pct is nan if not requested.
    """
    n        = len(key_list)
    sim_vals = np.array([sim_risks.get(k, 0.0) for k in key_list])
    corr     = spearman_safe(base_values, sim_vals)

    sim_order = sim_vals.argsort()[::-1]
    sim_top_n = set(key_list[i] for i in sim_order[:top_n])
    sim_topQ  = set(key_list[i] for i in sim_order[:max(1, n // 4)])

    top_n_pct = len(base_top_n & sim_top_n) / top_n * 100
    topQ_pct  = len(base_topQ  & sim_topQ)  / len(base_topQ) * 100

    cat_agree = float('nan')
    if base_risks_for_cat is not None and entities_for_cat:
        cat_agree = category_agreement_rate(base_risks_for_cat, sim_risks, entities_for_cat)

    return corr, top_n_pct, topQ_pct, cat_agree

# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 60)
    print("Sensitivity Analysis — Biofouling model parameters")
    print("=" * 60)
    print(f"sigma values:          {SIGMA_VALUES}")
    print(f"N simulations:         {N_SIMULATIONS}")
    print(f"stayAttenuation:       {STAY_ATTENUATION}")
    print(f"tradeAdjust:           {TRADE_ADJUST}")
    print(f"TOP_N_COUNTRY:         {TOP_N_COUNTRY}  (BRI countries only)")
    print(f"TOP_N_PAIR_GLOBAL:     {TOP_N_PAIR_GLOBAL}")
    print(f"TOP_N_PAIR_BRI:        {TOP_N_PAIR_BRI}")
    print(f"Lat threshold base:    {LAT_THRESHOLD_BASE:.4f} rad "
          f"(~{math.degrees(LAT_THRESHOLD_BASE):.1f} deg)")
    print(f"Lat threshold 15 deg:  {LAT_THRESHOLD_15:.4f} rad")
    print(f"Lat threshold 25 deg:  {LAT_THRESHOLD_25:.4f} rad")
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

    # Resolve BRI country names from IDs using the loaded portToCountry data
    _id_to_name = {}
    for row in portToCountry.values():
        _id_to_name[row['Country_id']] = row['Country_name']
    BRI_COUNTRY_NAMES = {_id_to_name[cid] for cid in BRI_COUNTRY_IDS
                         if cid in _id_to_name}
    missing_ids = BRI_COUNTRY_IDS - set(_id_to_name.keys())
    if missing_ids:
        print(f"  WARNING: BRI country IDs not found in portToCountry: {missing_ids}")
    print(f"  BRI country names resolved: {len(BRI_COUNTRY_NAMES)} / {len(BRI_COUNTRY_IDS)}")

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
                    'trip_duration': float(trip_dur),
                    'stay_duration': float(stay_dur),
                    'antifouling_p': float(prev_move['ANTIFOULING PROB']),
                    'distance':      float(distance),
                    'subseq':        subseq,
                    'vessel_type':   prev_move['VESSEL TYPE'],
                })
            except:
                pass
        else:
            ship_id = -1

    print(f"  ProbsRaw built: {len(probs_raw)} port pairs")

    # ---- Precompute p_alien (structural, never perturbed) ----
    print("Precomputing p_alien...")
    for pair, move_list in probs_raw.items():
        source, dest = pair
        p_alien = (1 if ECO == 'noEco_'
                   else filterbyeco_same_only(ports, source, dest) if ECO == 'sameEco_'
                   else filterbyeco(ports, source, dest))
        for move in move_list:
            move['_p_alien'] = p_alien
    print("  Done.")

    # ---- Baseline run ----
    print("\nComputing baseline (sigma=0, lat_threshold=baseline)...")
    base_country_risks, base_pair_risks = run_pipeline(
        probs_raw, ports, portToCountry, countryStayRatio,
        exportDict, importDict, rng=None, sigma=0.0,
        lat_threshold=LAT_THRESHOLD_BASE,
    )

    # --- Country-level setup (BRI only) ---
    bri_list = sorted([c for c in base_country_risks if c in BRI_COUNTRY_NAMES])
    print(f"  BRI countries found in results: {len(bri_list)} / {len(BRI_COUNTRY_NAMES)}")
    missing = BRI_COUNTRY_NAMES - set(bri_list)
    if missing:
        print(f"  WARNING: BRI countries not in results (check name mapping): {sorted(missing)}")
    base_vals_c, base_top5, base_topQ_c = make_baseline_info(
        base_country_risks, bri_list, TOP_N_COUNTRY)
    print(f"  Baseline top-{TOP_N_COUNTRY} (BRI): "
          f"{sorted(base_top5)}")

    # --- Pair-level setup (global) ---
    all_pairs = sorted(base_pair_risks.keys())
    base_vals_pg, base_top20, base_topQ_pg = make_baseline_info(
        base_pair_risks, all_pairs, TOP_N_PAIR_GLOBAL)
    print(f"  Global pairs: {len(all_pairs)}")

    # --- Pair-level setup (BRI only) ---
    bri_pairs = sorted([p for p in all_pairs
                        if p[0] in BRI_COUNTRY_NAMES and p[1] in BRI_COUNTRY_NAMES])
    base_vals_pb, base_top10_bri, base_topQ_pb = make_baseline_info(
        base_pair_risks, bri_pairs, TOP_N_PAIR_BRI)
    print(f"  BRI pairs: {len(bri_pairs)}")

    # ============================================================
    # PART 1 — MONTE CARLO PARAMETER PERTURBATION
    # ============================================================
    print("\n" + "=" * 60)
    print("PART 1 — Monte Carlo parameter perturbation")
    print("=" * 60)

    mc_results = {}
    for sigma in SIGMA_VALUES:
        print(f"\nMonte Carlo: sigma={sigma}, N={N_SIMULATIONS} ...")
        rng = np.random.RandomState(SEED)

        # collectors: country
        c_corr, c_top5, c_topQ, c_cat = [], [], [], []
        # collectors: global pair
        pg_corr, pg_top20, pg_topQ = [], [], []
        # collectors: BRI pair
        pb_corr, pb_top10, pb_topQ, pb_cat = [], [], [], []
        nan_count = 0

        for i in range(N_SIMULATIONS):
            if (i + 1) % 10 == 0:
                print(f"  {i+1}/{N_SIMULATIONS}")

            sim_c_risks, sim_p_risks = run_pipeline(
                probs_raw, ports, portToCountry, countryStayRatio,
                exportDict, importDict,
                rng=rng, sigma=sigma,
                lat_threshold=LAT_THRESHOLD_BASE,
            )

            # Country-level (BRI)
            corr, tn, tq, ca = compute_metrics(
                base_vals_c, sim_c_risks, bri_list,
                base_top5, base_topQ_c, TOP_N_COUNTRY,
                base_risks_for_cat=base_country_risks,
                entities_for_cat=bri_list,
            )
            if np.isnan(corr):
                nan_count += 1
            c_corr.append(corr); c_top5.append(tn)
            c_topQ.append(tq);   c_cat.append(ca)

            # Global pairs
            corr_pg, tn_pg, tq_pg, _ = compute_metrics(
                base_vals_pg, sim_p_risks, all_pairs,
                base_top20, base_topQ_pg, TOP_N_PAIR_GLOBAL,
            )
            pg_corr.append(corr_pg); pg_top20.append(tn_pg); pg_topQ.append(tq_pg)

            # BRI pairs
            corr_pb, tn_pb, tq_pb, ca_pb = compute_metrics(
                base_vals_pb, sim_p_risks, bri_pairs,
                base_top10_bri, base_topQ_pb, TOP_N_PAIR_BRI,
                base_risks_for_cat=base_pair_risks,
                entities_for_cat=bri_pairs,
            )
            pb_corr.append(corr_pb); pb_top10.append(tn_pb)
            pb_topQ.append(tq_pb);   pb_cat.append(ca_pb)

        mc_results[sigma] = {
            'nan_count': nan_count,
            'country': {
                'corr_mean':  np.nanmean(c_corr),
                'corr_std':   np.nanstd(c_corr),
                'corr_min':   np.nanmin(c_corr),
                'top_n_mean': np.mean(c_top5),
                'top_n_min':  np.min(c_top5),
                'topQ_mean':  np.mean(c_topQ),
                'topQ_min':   np.min(c_topQ),
                'cat_mean':   np.nanmean(c_cat),
                'cat_min':    np.nanmin(c_cat),
            },
            'pair_global': {
                'corr_mean':  np.nanmean(pg_corr),
                'corr_std':   np.nanstd(pg_corr),
                'corr_min':   np.nanmin(pg_corr),
                'top_n_mean': np.mean(pg_top20),
                'top_n_min':  np.min(pg_top20),
                'topQ_mean':  np.mean(pg_topQ),
                'topQ_min':   np.min(pg_topQ),
            },
            'pair_bri': {
                'corr_mean':  np.nanmean(pb_corr),
                'corr_std':   np.nanstd(pb_corr),
                'corr_min':   np.nanmin(pb_corr),
                'top_n_mean': np.mean(pb_top10),
                'top_n_min':  np.min(pb_top10),
                'topQ_mean':  np.mean(pb_topQ),
                'topQ_min':   np.min(pb_topQ),
                'cat_mean':   np.nanmean(pb_cat),
                'cat_min':    np.nanmin(pb_cat),
            },
        }

    # ============================================================
    # PART 2 — LATITUDE THRESHOLD SCENARIO ANALYSIS
    # ============================================================

    if RUN_LAT:
        print("\n" + "=" * 60)
        print("PART 2 — Latitude threshold scenario analysis")
        print("(fixed parameters, sigma=0)")
        print("=" * 60)

        scenario_results = {}
        for label, threshold in [('15 deg', LAT_THRESHOLD_15), ('25 deg', LAT_THRESHOLD_25)]:
            print(f"\nScenario: {label} ({threshold:.4f} rad)...")
            scen_c_risks, scen_p_risks = run_pipeline(
                probs_raw, ports, portToCountry, countryStayRatio,
                exportDict, importDict,
                rng=None, sigma=0.0, lat_threshold=threshold,
            )

            # Country-level (BRI)
            corr_c, tn_c, tq_c, ca_c = compute_metrics(
                base_vals_c, scen_c_risks, bri_list,
                base_top5, base_topQ_c, TOP_N_COUNTRY,
                base_risks_for_cat=base_country_risks,
                entities_for_cat=bri_list,
            )
            scen_c_order = np.array([scen_c_risks.get(c, 0.0) for c in bri_list]).argsort()[::-1]
            scen_top5_list = [bri_list[i] for i in scen_c_order[:TOP_N_COUNTRY]]

            # Global pairs
            corr_pg, tn_pg, tq_pg, _ = compute_metrics(
                base_vals_pg, scen_p_risks, all_pairs,
                base_top20, base_topQ_pg, TOP_N_PAIR_GLOBAL,
            )
            scen_pg_order = np.array([scen_p_risks.get(p, 0.0) for p in all_pairs]).argsort()[::-1]
            scen_top20_list = [all_pairs[i] for i in scen_pg_order[:TOP_N_PAIR_GLOBAL]]

            # BRI pairs
            corr_pb, tn_pb, tq_pb, ca_pb = compute_metrics(
                base_vals_pb, scen_p_risks, bri_pairs,
                base_top10_bri, base_topQ_pb, TOP_N_PAIR_BRI,
                base_risks_for_cat=base_pair_risks,
                entities_for_cat=bri_pairs,
            )
            scen_pb_order = np.array([scen_p_risks.get(p, 0.0) for p in bri_pairs]).argsort()[::-1]
            scen_top10_bri_list = [bri_pairs[i] for i in scen_pb_order[:TOP_N_PAIR_BRI]]

            scenario_results[label] = {
                'country':     {'corr': corr_c, 'top_n': tn_c, 'topQ': tq_c,
                                'cat': ca_c, 'top_n_list': scen_top5_list},
                'pair_global': {'corr': corr_pg, 'top_n': tn_pg, 'topQ': tq_pg,
                                'top_n_list': scen_top20_list},
                'pair_bri':    {'corr': corr_pb, 'top_n': tn_pb, 'topQ': tq_pb,
                                'cat': ca_pb, 'top_n_list': scen_top10_bri_list},
            }

    # ============================================================
    # PRINT RESULTS
    # ============================================================
    def _fmt(v, pct=False):
        if np.isnan(v):
            return '  nan  '
        return f"{v:.1f}%" if pct else f"{v:.4f}"

    print("\n" + "=" * 60)
    print("RESULTS — PART 1: Monte Carlo parameter perturbation")
    print("=" * 60)
    print(f"  Baseline top-{TOP_N_COUNTRY} BRI countries: {sorted(base_top5)}")
    print(f"  Baseline top-{TOP_N_PAIR_GLOBAL} global pairs: {sorted(base_top20)[:5]} ...")
    print(f"  Baseline top-{TOP_N_PAIR_BRI} BRI pairs:    {sorted(base_top10_bri)[:5]} ...\n")

    for sigma, res in mc_results.items():
        print(f"  sigma = {sigma}  (+-{int(sigma*100)}% independent perturbation)")
        if res['nan_count'] > 0:
            print(f"    [Note] {res['nan_count']} / {N_SIMULATIONS} runs produced NaN Spearman "
                  f"(collapsed probabilities); excluded from mean/std via nanmean.")

        r = res['country']
        print(f"\n  [A] Country-level — BRI ({len(bri_list)} countries)")
        print(f"    Spearman:            {_fmt(r['corr_mean'])} +- {_fmt(r['corr_std'])}"
              f"  (min = {_fmt(r['corr_min'])})")
        print(f"    Top-{TOP_N_COUNTRY} stability:    {_fmt(r['top_n_mean'], True)} mean"
              f"  (min = {_fmt(r['top_n_min'], True)})")
        print(f"    Top-Q stability:     {_fmt(r['topQ_mean'], True)} mean"
              f"  (min = {_fmt(r['topQ_min'], True)})")
        print(f"    Category agreement:  {_fmt(r['cat_mean'], True)} mean"
              f"  (min = {_fmt(r['cat_min'], True)})")

        r = res['pair_global']
        print(f"\n  [B] Pair-level — Global ({len(all_pairs)} pairs)")
        print(f"    Spearman:            {_fmt(r['corr_mean'])} +- {_fmt(r['corr_std'])}"
              f"  (min = {_fmt(r['corr_min'])})")
        print(f"    Top-{TOP_N_PAIR_GLOBAL} stability:   {_fmt(r['top_n_mean'], True)} mean"
              f"  (min = {_fmt(r['top_n_min'], True)})")
        print(f"    Top-Q stability:     {_fmt(r['topQ_mean'], True)} mean"
              f"  (min = {_fmt(r['topQ_min'], True)})")

        r = res['pair_bri']
        print(f"\n  [C] Pair-level — BRI ({len(bri_pairs)} pairs)")
        print(f"    Spearman:            {_fmt(r['corr_mean'])} +- {_fmt(r['corr_std'])}"
              f"  (min = {_fmt(r['corr_min'])})")
        print(f"    Top-{TOP_N_PAIR_BRI} stability:   {_fmt(r['top_n_mean'], True)} mean"
              f"  (min = {_fmt(r['top_n_min'], True)})")
        print(f"    Top-Q stability:     {_fmt(r['topQ_mean'], True)} mean"
              f"  (min = {_fmt(r['topQ_min'], True)})")
        print(f"    Category agreement:  {_fmt(r['cat_mean'], True)} mean"
              f"  (min = {_fmt(r['cat_min'], True)})")
        print()

    if RUN_LAT:
        print("=" * 60)
        print("RESULTS — PART 2: Latitude threshold scenarios")
        print("=" * 60)
        print(f"  Baseline (~20 deg) top-{TOP_N_COUNTRY} BRI: {sorted(base_top5)}\n")

        for label, res in scenario_results.items():
            print(f"  Scenario: {label}")

            r = res['country']
            print(f"  [A] Country-level — BRI")
            print(f"    Spearman:           {_fmt(r['corr'])}")
            print(f"    Top-{TOP_N_COUNTRY} overlap:     {_fmt(r['top_n'], True)}")
            print(f"    Top-Q overlap:      {_fmt(r['topQ'], True)}")
            print(f"    Category agreement: {_fmt(r['cat'], True)}")
            print(f"    Top-{TOP_N_COUNTRY} countries:   {r['top_n_list']}")

            r = res['pair_global']
            print(f"  [B] Pair-level — Global")
            print(f"    Spearman:           {_fmt(r['corr'])}")
            print(f"    Top-{TOP_N_PAIR_GLOBAL} overlap:    {_fmt(r['top_n'], True)}")
            print(f"    Top-Q overlap:      {_fmt(r['topQ'], True)}")
            print(f"    Top-{TOP_N_PAIR_GLOBAL} pairs:      {r['top_n_list'][:5]} ...")

            r = res['pair_bri']
            print(f"  [C] Pair-level — BRI")
            print(f"    Spearman:           {_fmt(r['corr'])}")
            print(f"    Top-{TOP_N_PAIR_BRI} overlap:    {_fmt(r['top_n'], True)}")
            print(f"    Top-Q overlap:      {_fmt(r['topQ'], True)}")
            print(f"    Category agreement: {_fmt(r['cat'], True)}")
            print(f"    Top-{TOP_N_PAIR_BRI} BRI pairs:  {r['top_n_list']}")
            print()


if __name__ == '__main__':
    main()
