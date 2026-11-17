#!/usr/bin/env python3
"""Species-flow higher-order network (SF-HON): ballast-water and biofouling introduction risk along ship trajectories.

Script version (2026-11-17) of the notebooks 02_Prepare_HON_Input, 03_Build-SF-HON, 04_make_HONet_dict and of
scripts/build_rules_utils.py, scripts/make_HONet_dic.py, in one file that needs only numpy. The notebooks were
removed on 2026-11-17 and are in tag v1.0.0. Step 0 (00_sort_clean_moves) is hon_moves.py. Not converted:
01_create_traces (its output is not read by later steps), 05_translate_clusters (GIS figures), and the trade and stay
adjustments of 02 (switched off there: tradeAdjust = stayAttenuation = False).

    python Code/HON_Code/hon_sf.py run --moves data/moves/moves_cleaned_2018.txt \
        --places data/Places_allportdata_mergedSept2017.csv --port-country data/stay/portToCountry.csv \
        --out <dir> [--legacy]
    (or the steps one by one: prepare, build, aggregate, with the same --out)

Steps:
  prepare    moves -> legs -> P(NIS) per leg for ballast water and biofouling, written as weighted sub-trajectories
             (orders 2..16) for the HON (02). A leg is two consecutive moves of a ship at different ports. Legs
             between countries, of the nine cargo types, after a stay under 20 days and with a trip of at least one
             hour get a risk record; the others keep their ports in the ship's trajectory.
  build      weighted sub-trajectories -> rules -> higher-order network (03, with the KLD test of the original code).
  aggregate  HON -> physical port network and port risks; link table; node dictionaries for clustering (04).

--legacy reproduces the notebooks as they were run for the 2018 outputs (02 with the moves filter: stays under 20
days, nine cargo types, no domestic legs), bugs included. Checked with compare_legacy.py: all 16 outputs have the same
lines as data/2018, values within 7e-13 (relative). Without --legacy the following are fixed (2026-11-17; found by
reading the code; the 2018 counts below are from the legacy run):
  1. prepare, ballast: P(intro) used the trip duration of the LAST leg in the file for every leg (a global variable
     left over from the loop); now each leg's own trip duration.
  2. prepare: the maximum used to scale P(NIS) to [0, 1] was taken with P(alien) = 0 for the same OR a neighbouring
     ecoregion, the written values with P(alien) = 0 for the same ecoregion only; 10 fouling values in 2018 were above
     1. Now one rule for both (--eco, default "eco": same or neighbouring ecoregion, as the "Eco_" setting).
  3. prepare, biofouling: the tropical test compared latitude in degrees with 0.35, a value in radians (~20 deg), so
     only ports within 0.35 deg of the equator counted as tropical. Now |lat| in radians < 0.35.
  4. prepare, biofouling: the stay entering the fouling model was the stay at the DESTINATION port (the move's own
     STAY DURATION), while 00 documents it as the stay at the source port. Now the source port's stay; the 20-day
     limit applies to it.
  5. prepare: dropping a move (domestic leg, stay >= 20 days) joined the ship's moves around it, creating legs
     never sailed (A->B->C became A->C, with B->C's duration and A-C's ports). Now a leg is two consecutive moves
     of the file (the leg's ROUT starts at the previous move's port). A domestic leg, or one after a stay of 20
     days or more, gets no risk record but its ports stay in the ship's trajectory (A, B, C), so the higher-order
     memory is kept. Legs whose risk cannot be computed (trip under one hour, which 00 rounds to 0 days; a bad
     stay) are treated the same way. The trajectory restarts only where the moves do not connect (new ship, a move
     missing, same port, or the year rule of 02).
  6. build: the first observation of a target under an existing source was stored twice, so it weighed double in
     the mean probability.
  7. aggregate: collapsing higher-order nodes summed probabilities (822 ballast and 82 fouling links above 1 in 2018,
     231 port risks outside [0, 1]); and 1 - prod(1 - r) lost every r below ~1e-16. Now P(at least one) computed as
     -expm1(sum(log1p(-r))), for the collapse and for port risks.
  8. aggregate: the ballast node dictionary for clustering was built from the FOULING network's nodes (a variable
     overwritten in 04), so ballast edges between other nodes were lost. Now each network's own nodes.
  9. data: nine Lloyd's places in Places_allportdata_mergedSept2017.csv have a sign or digit error in their
     coordinates (Panama Canal at 79.68 E, Vancouver Anch. at 123.3 E, ...; reference/place_coords_fix.csv, found
     2026-11-16). They enter leg distances (fouling speed term) and the link table. --coords-fix corrects them when
     the table still holds the recorded wrong values (not with --legacy).
Kept as in the original code (not bugs, but not the textbook HON): the KLD test uses p * exp(p / q) and the
threshold order / exp(1 + sum of the source's probabilities) (Saebi's code; the HON paper uses log and counts);
port risk "out" aggregates a port's outgoing links (04's agg_risk). Added: port risk "in" (incoming links).

Scale: 02 divides every P(NIS) by the largest in the file, so two years run separately are on different scales.
--scale-from <hon_summary.json ...> uses the largest maxima of the given runs instead (added 2026-11-17).

Outputs in --out: input_{ballast,fouling}.txt, rules_*.txt, honet_*.csv, net_*.txt (physical ports), links_*.csv,
port_risk_out_*.txt (04's format), port_risk.csv (out and in, both kinds), trans_honet_*.{net,dic}, hon_summary.json.
Moves and ports use Lloyd's place IDs (Places_allportdata_mergedSept2017.csv).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

TYPES9 = {"Auto", "Container", "Bulk", "Tanker", "Chemical", "Oil", "General", "Liquified-Gas", "Refrigerated-Cargo"}
EFFICACY = {"r0": 0, "r1": 0.7661, "r2": 0.99150495, "r3": 0.999371131}
CROSS_YEARS = ("1997", "1999", "2002", "2005", "2008", "2012", "2015", "2018")   # 02: legs may end in the next year
KINDS = ("ballast", "fouling")


def read_dict(fn, field, delim=","):
    with open(fn, newline="") as f:
        return {row[field]: row for row in csv.DictReader(f, delimiter=delim)}


def fix_coords(ports, path):
    """Correct places whose coordinates still equal the recorded wrong values (bug 9)."""
    done = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            q = ports.get(r["place_id"])
            if q is None:
                continue
            try:
                lat, lon = float(q["LATITUDE_DECIMAL"]), float(q["LONGITUDE_DECIMAL"])
            except ValueError:
                continue
            if math.isclose(lat, float(r["wrong_lat"]), abs_tol=1e-6) and \
                    math.isclose(lon, float(r["wrong_lon"]), abs_tol=1e-6):
                q["LATITUDE_DECIMAL"], q["LONGITUDE_DECIMAL"] = r["lat"], r["lon"]
                done.append(r["place_id"])
    return done


# ------------------------------------------------------------------ prepare (02) ---
class Model:
    """The per-leg risk terms of 02."""

    def __init__(self, ports, eco="eco", env=True, regime="r0", legacy=False):
        self.ports, self.eco, self.env, self.eff, self.legacy = ports, eco, env, EFFICACY[regime], legacy

    def distance(self, source, dest):
        p = self.ports
        lat1, lon1 = float(p[source]["LATITUDE_DECIMAL"]), float(p[source]["LONGITUDE_DECIMAL"])
        lat2, lon2 = float(p[dest]["LATITUDE_DECIMAL"]), float(p[dest]["LONGITUDE_DECIMAL"])
        d2r = lambda d: d * (math.pi / 180)
        dLat, dLon = d2r(lat2 - lat1), d2r(lon2 - lon1)
        a = math.sin(dLat / 2) * math.sin(dLat / 2) + math.cos(d2r(lat1)) * math.cos(d2r(lat2)) * math.sin(
            dLon / 2) * math.sin(dLon / 2)
        return 6371 * (2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)))

    def alien_eco(self, source, dest):
        """02 filterbyeco: 0 for the same or a neighbouring ecoregion (freshwater: FEOW), else 1."""
        p = self.ports
        p_alien = 1
        s_eco, d_eco = p[source]["MEOW_region"], p[dest]["MEOW_region"]
        if s_eco != "NA" and d_eco != "NA":
            if s_eco == d_eco or s_eco in p[dest]["MEOW_Neighbors"].split("|"):
                p_alien = 0
            else:
                p_alien = 1
        if s_eco == "NA" and d_eco != "NA":
            s_eco = p[source]["FEOW_region"]
            p_alien = 0 if (s_eco != "NA" and (s_eco in p[dest]["FEOW_Neighbors"].split("|")
                                                or s_eco in p[dest]["MEOW_Neighbors"].split("|"))) else 1
        if d_eco == "NA" and s_eco != "NA":
            d_eco = p[dest]["FEOW_region"]
            p_alien = 0 if (d_eco != "NA" and (s_eco in p[dest]["FEOW_Neighbors"].split("|")
                                                or s_eco in p[dest]["MEOW_Neighbors"].split("|"))) else 1
        return p_alien

    def alien_same(self, source, dest):
        """02 filterbyeco_same_only: 0 only for the same ecoregion."""
        try:
            s_eco, d_eco = self.ports[source]["MEOW_region"], self.ports[dest]["MEOW_region"]
            return 0 if (s_eco != "NA" and d_eco != "NA" and s_eco == d_eco) else 1
        except Exception:
            return -1

    def alien(self, source, dest, legacy_write=False):
        if self.eco == "none":
            return 1
        if self.eco == "same" or legacy_write:        # bug 2: 02 wrote with same-only whatever the setting
            return self.alien_same(source, dest)
        return self.alien_eco(source, dest)

    def establishment(self, source, dest):
        if not self.env:
            return 1
        try:
            p = self.ports
            sal_diff = abs(float(p[source]["Salinity"]) - float(p[dest]["Salinity"]))
            temp_diff = abs(float(p[source]["YR_MEAN_T"]) - float(p[dest]["YR_MEAN_T"]))
            return 0.00015 * math.exp(-0.5 * ((temp_diff / 2) ** 2 + (sal_diff / 10) ** 2))
        except Exception:
            return -1

    def intro_ballast(self, ballast, trip):
        ballast_discharge = (1 - math.exp(-3.22 * (10 ** -6) * ballast))
        return (1 - self.eff) * ballast_discharge * math.exp(-0.02 * trip)

    def fouling(self, source, duration, stay_d, distance, antifouling_p):
        v = distance / duration
        lat = float(self.ports[source]["LATITUDE_DECIMAL"])
        tropical = abs(lat) < 0.35 if self.legacy else abs(math.radians(lat)) < 0.35      # bug 3
        if tropical:
            return (0.000000129 * stay_d ** 3 - 0.000083165 * stay_d ** 2 + 0.01495187 * stay_d) * (
                antifouling_p) * math.exp(-0.008 * v)
        return (0.0000000014 * stay_d ** 3 - 0.000016566 * stay_d ** 2 + 0.00519377 * stay_d) * (
            antifouling_p) * math.exp(-0.008 * v)


def year_ok(first, second):
    fy = first["SAIL DATE"].split(" ")[0][0:4]
    sy = second["ARRIVAL DATE"].split(" ")[0][0:4]
    return fy != "" and sy != "" and (fy == sy or fy in CROSS_YEARS and int(sy) == int(fy) + 1)


def legs_legacy(rows, port_country):
    """02 cells 4-5: filter moves, then pair consecutive kept moves (bugs 4, 5)."""
    moves = []
    for row in rows:
        if float(row["STAY DURATION"]) < 20.0 and row["VESSEL TYPE"] in TYPES9:
            pp = row["ROUT"].split("-")
            if len(pp) != 2:
                continue
            if port_country[pp[0]]["Country_id"] == port_country[pp[1]]["Country_id"]:
                continue
            moves.append(row)
    out, last_trip = [], -1
    for k in range(len(moves) - 1):
        first, second = moves[k], moves[k + 1]
        ok = False
        if first["VESSEL ID"] == second["VESSEL ID"] and first["PLACE ID"] != second["PLACE ID"] \
                and year_ok(first, second):
            try:
                trip, stay = float(second["DURATION"]), float(second["STAY DURATION"])
                ok = trip > 0 and stay >= 0
            except ValueError:
                ok = False
        last_trip = trip if ok else -1
        out.append((first, second, ok, trip if ok else -1, stay if ok else -1))
    return out, {"moves_after_filter": len(moves)}, last_trip


def legs_fixed(rows, port_country, ports, stay_max=20.0):
    """A leg is two consecutive moves of a ship whose second move's ROUT starts at the first move's port (bug 5).
    Stay = the source port's stay (bug 4). Returns (first, second, connects, record, trip, stay): `connects` False
    restarts the trajectory; `record` False keeps the ports in the trajectory but writes no risk."""
    def country(pid):
        c = port_country.get(pid, {}).get("Country_id") or ports.get(pid, {}).get("COUNTRY_CODE")
        return c if c not in (None, "", "NA") else None
    out, why = [], defaultdict(int)
    for k in range(1, len(rows)):
        first, second = rows[k - 1], rows[k]
        pp = second["ROUT"].split("-")
        reason = None
        brk = skip = None
        trip = stay = -1
        if first["VESSEL ID"] != second["VESSEL ID"]:
            brk = "new_ship"
        elif len(pp) != 2 or pp[0] != first["PLACE ID"] or pp[1] != second["PLACE ID"]:
            brk = "not_consecutive"
        elif first["PLACE ID"] == second["PLACE ID"]:
            brk = "same_port"
        elif not year_ok(first, second):
            brk = "year"
        elif second["VESSEL TYPE"] not in TYPES9:
            skip = "type"
        elif country(pp[0]) is None or country(pp[1]) is None:
            skip = "no_country"
        elif country(pp[0]) == country(pp[1]):
            skip = "domestic"
        else:
            try:
                trip, stay = float(second["DURATION"]), float(first["STAY DURATION"])
                if not trip > 0:
                    skip = "trip_not_positive"          # under one hour (00 rounds down to whole hours)
                elif not 0 <= stay < stay_max:
                    skip = "source_stay"
            except ValueError:
                skip = "bad_duration"
        why[("break_" + brk) if brk else ("no_record_" + skip) if skip else "leg"] += 1
        out.append((first, second, brk is None, brk is None and skip is None, trip, stay))
    return out, dict(why)


def prepare(moves_path, ports, port_country, out: Path, order=16, eco="eco", env=True, regime="r0", legacy=False,
            scale=None):
    """scale: (max_ballast, max_fouling) to divide by instead of this file's maxima, so that runs on different
    years share one scale (the HON rule threshold and all risks depend on it); values may then exceed 1 only if a
    leg exceeds the given maxima, which the summary counts."""
    need = ("VESSEL ID", "PLACE ID", "ARRIVAL DATE", "SAIL DATE", "STAY DURATION", "ROUT", "DURATION", "VESSEL TYPE",
            "ANTIFOULING PROB", "BALLAST DISCHARGE")
    with open(moves_path, newline="") as f:
        rows = [{k: r[k] for k in need} for r in csv.DictReader(f, delimiter="|")]
    m = Model(ports, eco, env, regime, legacy)
    if legacy:
        pairs, info, last_trip = legs_legacy(rows, port_country)
        pairs = [(f, s, ok, ok, t, st) for f, s, ok, t, st in pairs]
    else:
        pairs, info = legs_fixed(rows, port_country, ports)
    # trajectories and per-leg records (02 cell 5)
    ship, traj, raw, err = -1, [], {}, 0
    for first, second, connects, record, trip, stay in pairs:
        if not connects:
            ship = -1
            continue
        source, dest = first["PLACE ID"], second["PLACE ID"]
        if first["VESSEL ID"] != ship:
            traj = [source, dest]
        else:
            traj.append(dest)
        ship = first["VESSEL ID"]
        if not record:
            continue
        try:
            dist = m.distance(source, dest)
            raw.setdefault((source, dest), []).append({
                "trip": trip, "stay": stay, "distance": float(dist), "traj": tuple(traj[-order:]),
                # 02 read the ship terms from the first move; they are per ship, so either move gives the same
                "ballast": float(first["BALLAST DISCHARGE"]), "antifouling": float(first["ANTIFOULING PROB"])})
        except Exception:
            err += 1
    # pass 1: maxima for scaling (02 cell 6)
    bw, fw = [], []
    for (s, d), recs in raw.items():
        p_alien, p_estab = m.alien(s, d), m.establishment(s, d)
        if p_estab == -1:
            continue
        for r in recs:
            bw.append(p_alien * m.intro_ballast(r["ballast"], last_trip if legacy else r["trip"]) * p_estab)
            fw.append(p_alien * m.fouling(s, r["trip"], r["stay"], r["distance"], r["antifouling"]) * p_estab)
    own = (max(bw), max(fw))
    max_b, max_f = scale if scale else own
    # pass 2: write
    out.mkdir(parents=True, exist_ok=True)
    n_written = {"ballast": 0, "fouling": 0}
    n_legs = {"ballast": 0, "fouling": 0}
    with open(out / "input_ballast.txt", "w") as fb, open(out / "input_fouling.txt", "w") as ff:
        for (s, d), recs in raw.items():
            p_alien, p_estab = m.alien(s, d, legacy_write=legacy), m.establishment(s, d)
            if p_estab == -1:
                continue
            for r in recs:
                subseq = [r["traj"][-n:] for n in range(2, len(r["traj"]) + 1)]
                pb = p_alien * m.intro_ballast(r["ballast"], last_trip if legacy else r["trip"]) * p_estab / max_b
                if pb > 0:
                    n_legs["ballast"] += 1
                    for q in subseq:
                        fb.write(" ".join(q) + " " + str(pb) + "\n")
                    n_written["ballast"] += len(subseq)
                pf = p_alien * m.fouling(s, r["trip"], r["stay"], r["distance"], r["antifouling"]) * p_estab / max_f
                if pf > 0:
                    n_legs["fouling"] += 1
                    for q in subseq:
                        ff.write(" ".join(q) + " " + str(pf) + "\n")
                    n_written["fouling"] += len(subseq)
    no_env = sum(len(v) for (s, d), v in raw.items() if m.establishment(s, d) == -1)
    return {"moves": len(rows), **info, "legs_valid": sum(1 for p in pairs if p[3]),
            "legs_no_port_data": err, "legs_with_record": sum(len(v) for v in raw.values()),
            "legs_no_env_data": no_env, "port_pairs": len(raw), "own_max_ballast": own[0], "own_max_fouling": own[1],
            "max_ballast": max_b, "max_fouling": max_f, "scaled_by": "given" if scale else "own maxima",
            "legs_above_scale": {"ballast": int(sum(v > max_b for v in bw)), "fouling": int(sum(v > max_f for v in fw))},
            "legs_written": n_legs, "lines_written": n_written, "order": order, "eco": eco, "env": env,
            "regime": regime}


# -------------------------------------------------------------------- build (03) ---
def kld(a, b):
    def prob(d, key):
        return d.get(str(key), 1e-6)
    divergence = 0
    for target in a:
        try:
            divergence += prob(a, target) * math.exp(prob(a, target) / prob(b, target))
        except OverflowError:
            divergence = float("inf")
    return divergence


def build_rules(path, max_order=16, min_support=1e-20, legacy=False):
    dist = defaultdict(dict)                                    # source -> target -> [p, ...]
    with open(path) as f:
        for line in f:
            fields = line.strip().split(" ")
            mv = fields[:-1]
            if len(mv) < 2:
                continue
            target, source, p = mv[-1], tuple(mv[:-1]), float(fields[-1])
            if legacy and source in dist:                       # bug 6
                if target not in dist[source]:
                    dist[source][str(target)] = [p]
                if target in dist[source]:
                    dist[source][target].append(p)
            else:
                dist[source].setdefault(target, []).append(p)
    final = defaultdict(dict)
    for source, value in dist.items():
        for target, val in value.items():
            fp = float(np.mean(val)) if len(val) > 1 else float(val[0])
            if fp > min_support:
                final[source][target] = fp
    del dist
    ext = {}                                                    # curr -> order -> sources ending in curr
    for source in list(final):
        if len(source) > 1:
            for start in range(1, len(source)):
                ext.setdefault(source[start:], {}).setdefault(len(source), set()).add(source)
    rules = {}

    def add(source):
        while len(source) > 0:
            rules[source] = final[source]
            source = source[:-1]

    def threshold(new_order, ext_source):
        return new_order / math.exp(1 + (sum(i for i in final[ext_source].values())))

    def extend(valid, curr, order):
        if order >= max_order:
            add(valid)
            return
        distr = final[valid]
        new_order = order + 1
        extended = ext.get(curr, {}).get(new_order, ())
        if len(extended) == 0:
            add(valid)
            return
        for ext_source in sorted(extended):                     # sorted: same rules, stable file order
            if kld(final[ext_source], distr) > threshold(new_order, ext_source):
                extend(ext_source, ext_source, new_order)
            else:
                extend(valid, ext_source, new_order)

    for source in list(final):
        if len(source) == 1:
            add(source)
            extend(source, source, 1)
    return rules


def build_network(rules):
    graph = defaultdict(dict)
    for source in sorted(rules, key=len):
        for target in rules[source]:
            graph[source][(target,)] = rules[source][target]
            if len(source) > 1:
                prev_source, prev_target = source[:-1], (source[-1],)
                if prev_source not in graph or source not in graph[prev_source]:
                    try:
                        graph[prev_source][source] = graph[prev_source][prev_target]
                        del graph[prev_source][prev_target]
                    except KeyError:
                        pass
    to_add, to_remove = [], []
    for source in graph:
        for target in graph[source]:
            if len(target) == 1:
                new_target = source + target
                while len(new_target) > 1:
                    if new_target in graph:
                        to_add.append((source, new_target, graph[source][target]))
                        to_remove.append((source, target))
                        break
                    new_target = new_target[1:]
    for source, target, w in to_add:
        graph[source][target] = w
    for source, target in to_remove:
        del graph[source][target]
    return graph


def seq_to_node(seq):
    node = seq[-1] + "|"
    seq = seq[:-1]
    while len(seq) > 0:
        node = node + seq[-1] + "."
        seq = seq[:-1]
    return node[:-1] if node[-1] == "." else node


def build(out: Path, max_order=16, min_support=1e-20, legacy=False):
    info = {}
    for kind in KINDS:
        rules = build_rules(out / f"input_{kind}.txt", max_order, min_support, legacy)
        net = build_network(rules)
        with open(out / f"rules_{kind}.txt", "w") as f:
            for s in rules:
                for t in rules[s]:
                    f.write(" ".join([" ".join(str(x) for x in s), "=>", t, str(rules[s][t])]) + "\n")
        with open(out / f"honet_{kind}.csv", "w") as f:
            for s in net:
                for t in net[s]:
                    f.write(",".join([seq_to_node(s), seq_to_node(t), str(net[s][t])]) + "\n")
        orders = defaultdict(int)
        for s in rules:
            orders[len(s)] += 1
        info[kind] = {"rules_sources": len(rules), "sources_by_order": dict(sorted(orders.items())),
                      "network_edges": sum(len(v) for v in net.values())}
    return {"max_order": max_order, "min_support": min_support, **info}


# ---------------------------------------------------------------- aggregate (04) ---
def union(rs):
    """P(at least one) = 1 - prod(1 - r), without losing r < 1e-16 (bug 7)."""
    r = np.asarray(rs, dtype=float)
    return float(-np.expm1(np.log1p(-np.clip(r, 0, 1)).sum()))


def read_honet(path):
    adj, nodes = defaultdict(dict), set()
    with open(path) as f:
        for line in f:
            a, b, w = line.split(",")[:3]
            a, b = a.strip(), b.strip()
            if a == b:
                continue
            w = float(w.strip())
            if w > 0:
                nodes.add(a)
                nodes.add(b)
                adj[a][b] = w
    phys = defaultdict(lambda: defaultdict(list))
    for a, tos in adj.items():
        for b, w in tos.items():
            phys[int(a.split("|")[0])][int(b.split("|")[0])].append(w)
    return phys, adj, nodes


def aggregate(out: Path, ports, legacy=False):
    nodes_of, adj_of, info = {}, {}, {}
    risks = {}
    for kind in KINDS:
        phys, adj, nodes = read_honet(out / f"honet_{kind}.csv")
        adj_of[kind], nodes_of[kind] = adj, nodes
        coll = {s: {d: (np.sum(v) if legacy else union(v)) for d, v in tos.items()} for s, tos in phys.items()}
        r_out = {}
        for s, tos in coll.items():
            if legacy:
                prod = 1
                for r in tos.values():
                    prod *= (1 - r)
                r_out[s] = 1 - prod
            else:
                r_out[s] = union(list(tos.values()))
        inc = defaultdict(list)
        for s, tos in coll.items():
            for d, r in tos.items():
                inc[d].append(r)
        r_in = {d: union(v) for d, v in inc.items()}
        risks[kind] = (r_out, r_in, coll)
        with open(out / f"net_{kind}.txt", "w") as f:
            for s, tos in coll.items():
                for d, r in tos.items():
                    f.write("{} {} {}\n".format(s, d, r))
        with open(out / f"port_risk_out_{kind}.txt", "w") as f:
            for p, r in sorted(r_out.items(), key=lambda x: x[1], reverse=True):
                f.write("{0},{1} \n".format(p, r))
        with open(out / f"links_{kind}.csv", "w", newline="") as f:
            w = csv.writer(f, lineterminator="\n") if not legacy else None
            head = ["source", "destination", "risk", "s_name", "d_name", "s_latt", "s_long", "d_latt", "d_long",
                    "s_country", "d_country"]
            f.write(",".join(head) + "\n")
            for s, tos in coll.items():
                for d, r in sorted(tos.items(), key=lambda x: x[1], reverse=True):
                    ps, pd_ = ports[str(s)], ports[str(d)]
                    row = [str(s), str(d), r, ps["NAME"], pd_["NAME"], ps["LATITUDE_DECIMAL"], ps["LONGITUDE_DECIMAL"],
                           pd_["LATITUDE_DECIMAL"], pd_["LONGITUDE_DECIMAL"], ps["COUNTRY_CODE"], pd_["COUNTRY_CODE"]]
                    if legacy:
                        f.write(",".join(map(str, row)) + "\n")
                    else:
                        w.writerow(row)
        vals = [r for tos in coll.values() for r in tos.values()]
        info[kind] = {"ports_out": len(r_out), "ports_in": len(r_in), "links": len(vals),
                      "links_above_1": int(sum(v > 1 for v in vals)),
                      "port_risk_out_outside_0_1": int(sum(not 0 <= v <= 1 for v in r_out.values())),
                      "port_risk_out_zero": int(sum(v == 0 for v in r_out.values()))}
    for kind in KINDS:                                          # bug 8: legacy used the fouling nodes for both
        nodes = nodes_of["fouling"] if legacy else nodes_of[kind]
        new = {old: i for i, old in enumerate(sorted(nodes))}
        with open(out / f"trans_honet_{kind}.dic", "w") as g:
            for old, i in new.items():
                g.write("{} {}\n".format(old, i))
        errors = 0
        with open(out / f"trans_honet_{kind}.net", "w") as h:
            for a, tos in adj_of[kind].items():
                for b, wt in tos.items():
                    if a in new and b in new:
                        h.write("{} {} {}\n".format(new[a], new[b], wt))
                    else:
                        errors += 1
        info[kind]["trans_edges_dropped"] = errors
    allp = sorted(set().union(*[set(risks[k][0]) | set(risks[k][1]) for k in KINDS]))
    with open(out / "port_risk.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["place_id", "name", "country", "lat", "lon", "ballast_out", "ballast_in", "fouling_out",
                    "fouling_in", "ballast_n_dest", "ballast_n_source", "fouling_n_dest", "fouling_n_source"])
        for p in allp:
            q = ports.get(str(p), {})
            rb, rf = risks["ballast"], risks["fouling"]
            w.writerow([p, q.get("NAME", ""), q.get("COUNTRY_CODE", ""), q.get("LATITUDE_DECIMAL", ""),
                        q.get("LONGITUDE_DECIMAL", ""), rb[0].get(p, 0.0), rb[1].get(p, 0.0), rf[0].get(p, 0.0),
                        rf[1].get(p, 0.0), len(rb[2].get(p, {})), sum(p in t for t in rb[2].values()),
                        len(rf[2].get(p, {})), sum(p in t for t in rf[2].values())])
    return info


# ------------------------------------------------------------------------- main ---
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("prepare", "build", "aggregate", "run"):
        s = sub.add_parser(name)
        s.add_argument("--out", type=Path, required=True)
        s.add_argument("--legacy", action="store_true", help="reproduce the notebooks, bugs included")
        if name in ("prepare", "aggregate", "run"):
            s.add_argument("--places", type=Path, required=True, help="Places_allportdata_mergedSept2017.csv")
            s.add_argument("--coords-fix", type=Path, default=Path(__file__).resolve().parent / "reference"
                           / "place_coords_fix.csv", help="place coordinate corrections (ignored with --legacy)")
        if name in ("prepare", "run"):
            s.add_argument("--moves", type=Path, required=True, help="moves_cleaned format ('|'-separated)")
            s.add_argument("--port-country", type=Path, required=True, help="portToCountry.csv")
            s.add_argument("--order", type=int, default=16)
            s.add_argument("--eco", choices=["eco", "same", "none"], default="eco")
            s.add_argument("--no-env", action="store_true")
            s.add_argument("--regime", choices=sorted(EFFICACY), default="r0", help="ballast treatment efficacy")
            s.add_argument("--scale-from", type=Path, nargs="*", default=None,
                           help="hon_summary.json of runs whose largest maxima set the scale (common scale)")
        if name in ("build", "run"):
            s.add_argument("--max-order", type=int, default=16)
            s.add_argument("--min-support", type=float, default=1e-20)
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    sf = a.out / "hon_summary.json"
    summary = json.loads(sf.read_text()) if sf.exists() else {}
    summary["legacy"] = a.legacy
    ports = read_dict(a.places, "ID") if hasattr(a, "places") else None
    if ports is not None and not a.legacy and a.coords_fix:
        summary["coords_fixed"] = fix_coords(ports, a.coords_fix)
    if a.cmd in ("prepare", "run"):
        summary["inputs"] = {"moves": str(a.moves), "places": str(a.places), "port_country": str(a.port_country)}
        scale = None
        if a.scale_from:
            ps = [json.loads(Path(f).read_text())["prepare"] for f in a.scale_from]
            scale = (max(q.get("own_max_ballast", q["max_ballast"]) for q in ps),
                     max(q.get("own_max_fouling", q["max_fouling"]) for q in ps))
            summary["scale_from"] = [str(f) for f in a.scale_from]
        summary["prepare"] = prepare(a.moves, ports, read_dict(a.port_country, "Port_id"), a.out, a.order, a.eco,
                                     not a.no_env, a.regime, a.legacy, scale)
        print("prepare", json.dumps(summary["prepare"]))
    if a.cmd in ("build", "run"):
        summary["build"] = build(a.out, a.max_order, a.min_support, a.legacy)
        print("build", json.dumps(summary["build"]))
    if a.cmd in ("aggregate", "run"):
        summary["aggregate"] = aggregate(a.out, ports, a.legacy)
        print("aggregate", json.dumps(summary["aggregate"]))
    sf.write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
