"""hon_sf.py: one test per fixed bug, on synthetic ports and moves; the legacy path must keep the old behaviour.

    python -m pytest Code/HON_Code/tests
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hon_sf as hs  # noqa: E402

COLS = ["ID", "NAME", "ABBREV", "COUNTRY_CODE", "LATITUDE_DECIMAL", "LONGITUDE_DECIMAL", "LATDD", "LONGDD",
        "YR_MEAN_T", "Salinity", "MEOW_region", "MEOW_Neighbors", "FEOW_region", "FEOW_Neighbors"]


def port(pid, country, lat, lon, eco, t=20.0, s=35.0, nb="NA"):
    return {"ID": pid, "NAME": f"P{pid}", "ABBREV": country, "COUNTRY_CODE": country, "LATITUDE_DECIMAL": str(lat),
            "LONGITUDE_DECIMAL": str(lon), "LATDD": str(lat), "LONGDD": str(lon), "YR_MEAN_T": str(t),
            "Salinity": str(s), "MEOW_region": eco, "MEOW_Neighbors": nb, "FEOW_region": "NA", "FEOW_Neighbors": "NA"}


PORTS = {p["ID"]: p for p in [
    port("1", "CHN", 31.0, 121.0, "E1"), port("2", "CHN", 30.0, 122.0, "E1"), port("3", "USA", 34.0, -118.0, "E2"),
    port("4", "JPN", 35.0, 139.0, "E3", nb="E1"), port("5", "PER", -11.6, -77.3, "E4"),
    port("6", "ECU", -2.2, -79.9, "E5", t=26.0)]}
P2C = {k: {"Country_id": v["COUNTRY_CODE"]} for k, v in PORTS.items()}
HEAD = ["VESSEL ID", "PLACE ID", "ARRIVAL DATE", "SAIL DATE", "STAY DURATION", "ROUT", "DURATION", "VESSEL TYPE",
        "ANTIFOULING PROB", "GWT", "BALLAST DISCHARGE"]


def mv(ship, place, prev, arr, sail, stay, trip, typ="Container", bd="1000"):
    return {"VESSEL ID": ship, "PLACE ID": place, "ARRIVAL DATE": arr, "SAIL DATE": sail, "STAY DURATION": str(stay),
            "ROUT": f"{prev}-{place}", "DURATION": str(trip), "VESSEL TYPE": typ, "ANTIFOULING PROB": "0.19",
            "GWT": "5000", "BALLAST DISCHARGE": bd}


def write_moves(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEAD, delimiter="|")
        w.writeheader()
        w.writerows(rows)


def read_input(path):
    out = []
    for line in open(path):
        f = line.split()
        out.append((tuple(f[:-1]), float(f[-1])))
    return out


def test_tropical_latitude_in_radians():                       # bug 3
    m_fix, m_old = hs.Model(PORTS, legacy=False), hs.Model(PORTS, legacy=True)
    # port 6 at 2.2 deg S is tropical (< 0.35 rad = 20 deg); the old test (|lat| < 0.35 deg) says temperate
    trop = (0.000000129 * 5 ** 3 - 0.000083165 * 5 ** 2 + 0.01495187 * 5) * 0.19 * math.exp(-0.008 * 100)
    temp = (0.0000000014 * 5 ** 3 - 0.000016566 * 5 ** 2 + 0.00519377 * 5) * 0.19 * math.exp(-0.008 * 100)
    assert m_fix.fouling("6", 10, 5, 1000, 0.19) == pytest.approx(trop)
    assert m_old.fouling("6", 10, 5, 1000, 0.19) == pytest.approx(temp)
    assert m_fix.fouling("4", 10, 5, 1000, 0.19) == pytest.approx(temp)          # 35 N: temperate in both


def test_alien_rules():
    m = hs.Model(PORTS)
    assert m.alien_eco("1", "2") == 0 and m.alien_eco("1", "4") == 0 and m.alien_eco("1", "3") == 1
    assert m.alien_same("1", "4") == 1
    assert m.alien("1", "4") == 0 and m.alien("1", "4", legacy_write=True) == 1  # bug 2: legacy wrote same-only


def test_trip_stay_and_trajectory_break(tmp_path):              # bugs 1, 4, 5
    rows = [mv("S", "1", "5", "2018-01-01 00:00:00", "2018-01-03 00:00:00", 2, 1),
            mv("S", "2", "1", "2018-01-04 00:00:00", "2018-01-09 00:00:00", 5, 1),     # domestic 1 -> 2
            mv("S", "3", "2", "2018-01-20 00:00:00", "2018-01-23 00:00:00", 3, 11),    # 2 -> 3
            mv("S", "4", "3", "2018-02-02 00:00:00", "2018-02-06 00:00:00", 4, 10)]    # 3 -> 4
    write_moves(tmp_path / "m.txt", rows)
    fix, old = tmp_path / "fix", tmp_path / "old"
    s_fix = hs.prepare(tmp_path / "m.txt", PORTS, P2C, fix)
    s_old = hs.prepare(tmp_path / "m.txt", PORTS, P2C, old, legacy=True)
    b_fix, b_old = read_input(fix / "input_ballast.txt"), read_input(old / "input_ballast.txt")
    # fixed: legs 2->3 and 3->4 written; the domestic leg 1->2 has no record but 1 stays in the trajectory
    assert [q for q, _ in b_fix] == [("2", "3"), ("1", "2", "3"), ("3", "4"), ("2", "3", "4"), ("1", "2", "3", "4")]
    assert s_fix["legs_valid"] == 2 and s_fix["no_record_domestic"] == 1
    # legacy: the row arriving at 2 is dropped and 1 is joined to 3 (a leg never sailed)
    assert ("1", "3") in [q for q, _ in b_old]
    # bug 1: each leg its own trip (11 vs 10 days), so the two ballast values differ by exp(0.02)
    p23, p34 = b_fix[0][1], b_fix[2][1]
    assert p23 / p34 == pytest.approx(math.exp(-0.02 * 11) / math.exp(-0.02 * 10)
                                      * hs.Model(PORTS).establishment("2", "3") / hs.Model(PORTS).establishment("3", "4"))
    # bug 4: fouling of 3->4 uses the stay at 3 (3 days), not at 4 (4 days)
    m = hs.Model(PORTS)
    f = read_input(fix / "input_fouling.txt")
    ratio = f[2][1] / f[0][1]
    want = (m.fouling("3", 10, 3, m.distance("3", "4"), 0.19) * m.establishment("3", "4")) / (
        m.fouling("2", 11, 5, m.distance("2", "3"), 0.19) * m.establishment("2", "3"))
    assert ratio == pytest.approx(want)


def test_long_source_stay_no_record(tmp_path):
    rows = [mv("S", "3", "9", "2018-01-01 00:00:00", "2018-01-30 00:00:00", 29, 1),   # 29-day stay at 3
            mv("S", "4", "3", "2018-02-08 00:00:00", "2018-02-10 00:00:00", 2, 9),
            mv("S", "5", "4", "2018-03-08 00:00:00", "2018-03-10 00:00:00", 2, 26)]
    write_moves(tmp_path / "m.txt", rows)
    s = hs.prepare(tmp_path / "m.txt", PORTS, P2C, tmp_path / "o")
    assert s["no_record_source_stay"] == 1 and s["legs_valid"] == 1
    assert [q for q, _ in read_input(tmp_path / "o" / "input_ballast.txt")] == [("4", "5"), ("3", "4", "5")]


def test_scaled_to_one(tmp_path):                               # bug 2: max and written values use one rule
    rows = [mv("S", "3", "9", "2018-01-01 00:00:00", "2018-01-03 00:00:00", 2, 1),
            mv("S", "4", "3", "2018-01-12 00:00:00", "2018-01-14 00:00:00", 2, 9),
            mv("S", "1", "4", "2018-01-16 00:00:00", "2018-01-18 00:00:00", 2, 2)]      # 4 -> 1: neighbours
    write_moves(tmp_path / "m.txt", rows)
    hs.prepare(tmp_path / "m.txt", PORTS, P2C, tmp_path / "o")
    vals = [p for _, p in read_input(tmp_path / "o" / "input_fouling.txt")]
    assert max(vals) == pytest.approx(1.0) and all(0 < v <= 1 for v in vals)


def test_first_observation_counted_once(tmp_path):              # bug 6
    (tmp_path / "in.txt").write_text("1 2 0.2\n1 3 0.5\n1 3 0.1\n")
    fix = hs.build_rules(tmp_path / "in.txt")
    old = hs.build_rules(tmp_path / "in.txt", legacy=True)
    assert fix[("1",)]["3"] == pytest.approx(0.3)
    assert old[("1",)]["3"] == pytest.approx((0.5 + 0.5 + 0.1) / 3)


def test_higher_order_rule_and_network(tmp_path):
    # 1->2 is followed by 3 only after 4->1, by 5 only after 6->1: a second-order dependency
    lines = ["1 2 0.5", "2 3 0.5", "2 5 0.5", "4 1 0.5", "6 1 0.5"] + ["4 1 2 3 0.9"] * 3 + ["6 1 2 5 0.9"] * 3 + \
            ["1 2 3 0.9"] * 3 + ["1 2 5 0.9"] * 3
    (tmp_path / "in.txt").write_text("\n".join(lines) + "\n")
    rules = hs.build_rules(tmp_path / "in.txt")
    net = hs.build_network(rules)
    assert ("1",) in rules and ("2",) in rules
    nodes = {hs.seq_to_node(s) for s in net} | {hs.seq_to_node(t) for v in net.values() for t in v}
    assert "2|" in nodes and "1|" in nodes
    assert hs.seq_to_node(("4", "1", "2")) == "2|1.4"


def test_union_precision_and_range():                          # bug 7
    assert hs.union([1e-18, 1e-18]) == pytest.approx(2e-18)
    assert 1 - (1 - 1e-18) * (1 - 1e-18) == 0.0                # the old formula
    assert hs.union([0.7, 0.8]) == pytest.approx(1 - 0.3 * 0.2)
    assert hs.union([0.9, 0.9, 0.9]) <= 1


def test_aggregate_own_nodes_and_in_out(tmp_path):            # bug 8
    (tmp_path / "honet_ballast.csv").write_text("1|,3|,0.5\n3|,4|,0.25\n5|,4|,0.5\n")
    (tmp_path / "honet_fouling.csv").write_text("1|,3|,0.5\n")
    info = hs.aggregate(tmp_path, PORTS)
    assert info["ballast"]["trans_edges_dropped"] == 0
    pr = pd.read_csv(tmp_path / "port_risk.csv").set_index("place_id")
    assert pr.loc[4, "ballast_in"] == pytest.approx(1 - 0.75 * 0.5)
    assert pr.loc[3, "ballast_out"] == pytest.approx(0.25) and pr.loc[3, "ballast_in"] == pytest.approx(0.5)
    old = hs.aggregate(tmp_path, PORTS, legacy=True)
    assert old["ballast"]["trans_edges_dropped"] == 2         # the fouling nodes miss 4 and 5


def test_fix_coords(tmp_path):
    ports = {"5569": port("5569", "PAN", 9.08, 79.68, "E9")}
    (tmp_path / "fix.csv").write_text("place_id,place_name,wrong_lat,wrong_lon,lat,lon,reason\n"
                                      "5569,Panama Canal,9.08,79.68,9.08,-79.68,sign\n")
    assert hs.fix_coords(ports, tmp_path / "fix.csv") == ["5569"]
    assert ports["5569"]["LONGITUDE_DECIMAL"] == "-79.68"
    assert hs.fix_coords(ports, tmp_path / "fix.csv") == []      # only while the wrong value is there


def test_common_scale(tmp_path):
    rows = [mv("S", "3", "9", "2018-01-01 00:00:00", "2018-01-03 00:00:00", 2, 1),
            mv("S", "4", "3", "2018-01-12 00:00:00", "2018-01-14 00:00:00", 2, 9)]
    write_moves(tmp_path / "m.txt", rows)
    own = hs.prepare(tmp_path / "m.txt", PORTS, P2C, tmp_path / "a")
    half = hs.prepare(tmp_path / "m.txt", PORTS, P2C, tmp_path / "b",
                      scale=(own["own_max_ballast"] * 2, own["own_max_fouling"] * 2))
    assert read_input(tmp_path / "a" / "input_ballast.txt")[0][1] == pytest.approx(1.0)
    assert read_input(tmp_path / "b" / "input_ballast.txt")[0][1] == pytest.approx(0.5)
    assert half["legs_above_scale"] == {"ballast": 0, "fouling": 0}


def test_short_hop_keeps_trajectory(tmp_path):
    # 3 -> 4 within the hour (DURATION 0) gets no record, but 3 stays in the trajectory of 4 -> 5
    rows = [mv("S", "3", "9", "2018-01-01 00:00:00", "2018-01-03 00:00:00", 2, 1),
            mv("S", "4", "3", "2018-01-03 00:30:00", "2018-01-05 00:00:00", 2, 0),
            mv("S", "5", "4", "2018-01-30 00:00:00", "2018-02-01 00:00:00", 2, 25)]
    write_moves(tmp_path / "m.txt", rows)
    s = hs.prepare(tmp_path / "m.txt", PORTS, P2C, tmp_path / "o")
    assert s["no_record_trip_not_positive"] == 1
    assert [q for q, _ in read_input(tmp_path / "o" / "input_ballast.txt")] == [("4", "5"), ("3", "4", "5")]
