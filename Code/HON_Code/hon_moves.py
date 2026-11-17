#!/usr/bin/env python3
"""Step 0 of the SF-HON pipeline: raw Lloyd's moves -> moves_cleaned (one row per move, with the previous port, trip
and stay durations and the ship terms of the risk model), the input of hon_sf.py.

Script version (2026-11-17) of 00_sort_clean_moves.ipynb (its 2015-2018 branch); the notebook is in tag v1.0.0.

    python Code/HON_Code/hon_moves.py --raw moves_1819.lst --vessels <ship terms .csv> --out moves_cleaned_2018.txt
    python Code/HON_Code/hon_moves.py --raw moves_1819.lst --terms-from moves_cleaned_2018.txt --out <file>
    (--terms-from reuses the ship terms of an existing moves_cleaned file; --legacy keeps 00's behaviour)

Rules (as in 00):
  * moves without an arrival or a sail date are dropped; dates dd/mm/yyyy [HH:MM];
  * moves are grouped by ship and sorted by arrival; ships in the order of their ID as text;
  * a move is written if the previous move (in that order) is of the same ship, at another place, sailed no later
    than this arrival, and the ship has its terms (type with an antifouling probability, gross tonnage >= 0,
    ballast discharge);
  * STAY DURATION = this move's sail - arrival; DURATION = this arrival - previous sail; both as whole days plus whole
    hours / 24, rounded to 2 decimals (00's formula; a gap under one hour is 0); ROUT = previous place - this place;
  * the 14th raw column (POSSIBLE ATM RECORD) is dropped.
Ship terms (--vessels): CSV with VESSEL.ID, VESSEL.TYPE (one of the antifouling types below, or a Lloyd's type code),
GROSS, and the ballast discharge as `ballast_discharge`, or `z` and `release` (z * release, the 2015-2018 branch), or
`random_forest_br` (the 1997-2012 branch). --terms-from reads them from an existing moves_cleaned file instead.

Fixed (2026-11-17; --legacy keeps the old behaviour):
  1. 00 grouped moves by ship with `if p_vessel_id == VesselID: append`, so the first move of every ship (in file
     order) was never kept, and the next one, now without a previous move, was not written either: each ship lost a
     call and the leg into it. On the 2018-19 raw moves (moves_1819.lst): 63,370 ships lost their second move, and
     47,254 of these would have been rows (1.3% of the corrected file; 80% arriving in May 2018).
  2. 00 read dates through local time (time.mktime / fromtimestamp). It was run in a US time zone, so times between
     02:00 and 03:00 on 10 March 2019 (start of US daylight saving) were moved one hour later, changing those moves'
     stays and trips (644 rows of moves_cleaned_2018.txt). Dates are now read as written (no time zone).
--legacy reproduces moves_cleaned_2018.txt byte for byte from moves_1819.lst (ship terms taken from that file).
Not carried over: the 1997-2012 branch's own loop (its row index advanced only on written rows, so "the previous
move" drifted after a skipped row); its ship terms can be given through --vessels instead.
Note: STAY DURATION is the stay at THIS move's place (00's comment calls it the source port's stay); hon_sf.py takes
the source port's stay from the previous row.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

HEADER = ("VESSEL ID|MOVE REFNO|MOVE ID|PLACE ID|MOVETYPE|MOVETYPE QUALIFIER|ARRIVAL DATE|ARRIVAL DATE ESTIMATED|"
          "ARRIVAL DATE QUALIFIER|SAIL DATE|SAIL DATE ESTIMATED|SAIL DATE QUALIFIER|EDIT DATE|STAY DURATION|ROUT|"
          "DURATION|VESSEL TYPE|ANTIFOULING PROB|GWT|BALLAST DISCHARGE")
ANTIFOULING = {"Auto": 0.2, "Container": 0.19, "Tanker": 0.3, "Passenger": 0.31, "Bulk": 0.42, "General": 0.53,
               "Other": 0.63, "Chemical": 0.63, "Liquified-Gas": 0.63, "Oil": 0.63, "Refrigerated-Cargo": 0.63,
               "Fishing": 0.63, "Research": 0.63, "Yacht": 0.63}
TYPE_CODES = {"Auto": ["OLC", "PRR", "URC", "URR", "MVE"], "Container": ["UBC", "UCC", "UCR"],
              "Passenger": ["MPR", "OFY"], "Bulk": ["BBU", "BCE", "BOR", "BWC", "UBG"],
              "General": ["BCB", "GCT", "GGC", "GPC"],
              "Other": ["XTG", "OBA", "OTB", "XAH", "XTS", "OSY", "OPO", "DDR", "DTS", "XPT", "DSS", "DTD", "OCL",
                        "OHL", "DHD"],
              "Chemical": ["TAC", "TAS", "TBK", "TCH", "TPD"], "Liquified-Gas": ["LNG", "LNP", "LPG"],
              "Oil": ["TCR", "COO", "CBO", "TCO"], "Refrigerated-Cargo": ["GRF"],
              "Fishing": ["FFC", "FFF", "FFS", "FTR", "FWH"], "Research": ["RRE"], "Yacht": ["OYT"]}
# 00's code table, after its own overwrites (its Tanker codes end up as Chemical, Liquified-Gas or Oil; OYT as Yacht)
CODE_TO_TYPE = {c: t for t, cs in TYPE_CODES.items() for c in cs}
NAMES = {"General Cargo": "General", "RoRo": "Auto", "Reefer": "Refrigerated-Cargo", "Bulker": "Bulk",
         "Tanker": "Tanker"}                                  # 00's mapping of NBIC type names (1997-2012 branch)


def parse_date(s: str):
    s = s.strip()
    for fmt in ("%d/%m/%Y %H:%M", "%d/%m/%Y"):
        try:
            return datetime.strptime(s if fmt == "%d/%m/%Y %H:%M" else s.split(" ")[0], fmt)
        except ValueError:
            continue
    return None


def parse_date_legacy(s: str):
    """00: through local time (time.mktime, datetime.fromtimestamp); clean() sets a US time zone for --legacy."""
    s = s.strip()
    try:
        ts = time.strptime(s, "%d/%m/%Y %H:%M")
    except ValueError:
        try:
            ts = time.strptime(s.split(" ")[0], "%d/%m/%Y")
        except ValueError:
            return None
    return datetime.fromtimestamp(time.mktime(ts))


def days(td) -> float:
    return td.days + round(td.seconds // 3600 / 24, 2)


def read_terms(path: Path) -> dict:
    """Ship terms from a CSV (see the docstring): id -> (type, gross tonnage, ballast discharge)."""
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            t = r.get("VESSEL.TYPE", "").strip().strip('"')
            t = NAMES.get(t, CODE_TO_TYPE.get(t, t))
            if r.get("ballast_discharge") not in (None, "", "NA"):
                bd = r["ballast_discharge"]
            elif r.get("z") not in (None, "", "NA") and r.get("release") not in (None, "", "NA"):
                bd = float(r["z"]) * float(r["release"])
            elif r.get("random_forest_br") not in (None, "", "NA"):
                bd = r["random_forest_br"]
            else:
                continue
            if r.get("GROSS") in (None, "", "NA"):
                continue
            out[r["VESSEL.ID"].strip()] = (t, float(r["GROSS"]), float(bd))
    return out


def terms_from_cleaned(path: Path) -> dict:
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="|"):
            v = r["VESSEL ID"].strip()
            if v not in out:
                out[v] = (r["VESSEL TYPE"], float(r["GWT"]), float(r["BALLAST DISCHARGE"]))
    return out


def clean(raw: Path, terms: dict, out: Path, legacy=False) -> dict:
    n = defaultdict(int)
    parse = parse_date
    old_tz = os.environ.get("TZ")
    if legacy:                       # 00 was run in a US time zone (see bug 2); any US zone with DST gives the same
        os.environ["TZ"] = "America/New_York"
        time.tzset()
        parse = parse_date_legacy
    try:
        return _clean(raw, terms, out, legacy, parse, n)
    finally:
        if legacy:
            if old_tz is None:
                os.environ.pop("TZ", None)
            else:
                os.environ["TZ"] = old_tz
            time.tzset()


def _clean(raw, terms, out, legacy, parse, n):
    by_ship = defaultdict(list)
    seen = set()
    with open(raw, newline="") as f:
        next(f)                                                    # header
        for line in f:
            l = [x.strip() for x in line.rstrip("\n").split("|")]
            n["raw"] += 1
            if len(l) == 14:
                l = l[:-1]
            if l[6] == "" or l[9] == "":
                n["no_arrival_or_sail"] += 1
                continue
            a, s = parse(l[6]), parse(l[9])
            if a is None or s is None:
                n["bad_date"] += 1
                continue
            ship = l[0]
            # legacy bug 1: the first move of each ship in file order (00 sorted by ship id with a stable sort first)
            if legacy and ship not in seen:
                seen.add(ship)
                n["legacy_first_move_dropped"] += 1
                continue
            seen.add(ship)
            by_ship[ship].append((a, s, l))
    n["ships"] = len(by_ship)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as o:
        o.write(HEADER + "\n")
        for ship in sorted(by_ship):                                # ship ids as text, as in 00
            moves = sorted(by_ship[ship], key=lambda m: m[0])
            t = terms.get(ship)
            for k in range(1, len(moves)):
                (pa, ps, pl), (a, s, l) = moves[k - 1], moves[k]
                if t is None:
                    n["no_ship_terms"] += 1
                    continue
                if l[3] == pl[3]:
                    n["same_place"] += 1
                    continue
                vtype, gwt, bd = t
                if vtype not in ANTIFOULING:
                    n["type_without_antifouling"] += 1
                    continue
                if not (gwt >= 0 and a >= ps):
                    n["overlap_or_gwt"] += 1
                    continue
                row = list(l)
                row[6], row[9] = str(a), str(s)
                rout = f"{days(s - a)}|{pl[3]}-{l[3]}|{days(a - ps)}"
                o.write("|".join(row) + "|" + rout + "|" + f"{vtype}|{ANTIFOULING[vtype]}|{gwt}|{bd}" + "\n")
                n["written"] += 1
    return dict(n)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", type=Path, required=True, help="raw Lloyd's moves ('|'-separated, with header)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--vessels", type=Path, help="ship terms CSV")
    g.add_argument("--terms-from", type=Path, help="take the ship terms from an existing moves_cleaned file")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--legacy", action="store_true", help="keep 00's first-move bug")
    a = ap.parse_args(argv)
    terms = read_terms(a.vessels) if a.vessels else terms_from_cleaned(a.terms_from)
    s = clean(a.raw, terms, a.out, a.legacy)
    s.update({"raw": str(a.raw), "terms": str(a.vessels or a.terms_from), "ships_with_terms": len(terms),
              "legacy": a.legacy})
    Path(str(a.out) + ".summary.json").write_text(json.dumps(s, indent=1))
    print(json.dumps(s))


if __name__ == "__main__":
    main()
