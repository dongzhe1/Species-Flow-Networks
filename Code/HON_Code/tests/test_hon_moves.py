"""hon_moves.py (step 0 of SF-HON, from 00_sort_clean_moves): the two fixed bugs and the cleaning rules.

    python -m pytest Code/HON_Code/tests
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import hon_moves as hm  # noqa: E402

HEAD = ("VESSEL ID|MOVE REFNO|MOVE ID|PLACE ID|MOVE TYPE|MOVE QUALIFIER|ARRIVAL DATE|ARRIVAL ESTIMATED|"
        "ARRIVAL QUALIFIER|SAILING DATE|SAILING ESTIMATED|SAILING QUALIFIER|RECORD LAST UPDATED|POSSIBLE ATM RECORD")


def raw_row(ship, ref, place, arr, sail, atm="Y"):
    return f"{ship}|{ref}|1|{place}|||{arr}|||{sail}|||01/06/2018|{atm}"


def run(tmp_path, rows, terms, legacy=False):
    (tmp_path / "raw.lst").write_text(HEAD + "\n" + "\n".join(rows) + "\n")
    out = tmp_path / ("legacy.txt" if legacy else "fixed.txt")
    n = hm.clean(tmp_path / "raw.lst", terms, out, legacy)
    with open(out, newline="") as f:
        return list(csv.DictReader(f, delimiter="|")), n


TERMS = {"7": ("Container", 5000.0, 100.0), "8": ("Bulk", 9000.0, 2000.0)}
ROWS = [raw_row("7", "a1", "10", "01/05/2018 00:00", "02/05/2018 06:00"),
        raw_row("7", "a2", "20", "05/05/2018 09:30", "06/05/2018 12:00"),
        raw_row("7", "a3", "20", "07/05/2018 00:00", "07/05/2018 08:00"),     # same place: not a row
        raw_row("7", "a4", "30", "10/05/2018 08:00", "11/05/2018 08:00"),
        raw_row("8", "b1", "10", "03/05/2018", "04/05/2018"),               # date only
        raw_row("8", "b2", "40", "06/05/2018 00:00", "07/05/2018 00:00"),
        raw_row("9", "c1", "10", "03/05/2018 00:00", "04/05/2018 00:00"),     # no ship terms
        raw_row("9", "c2", "40", "06/05/2018 00:00", "07/05/2018 00:00")]


def test_first_move_kept(tmp_path):                            # bug 1
    fixed, n = run(tmp_path, ROWS, TERMS)
    old, n_old = run(tmp_path, ROWS, TERMS, legacy=True)
    assert [r["MOVE REFNO"] for r in fixed] == ["a2", "a4", "b2"]
    assert [r["MOVE REFNO"] for r in old] == ["a4"]              # a1, b1 dropped, so a2 and b2 have no previous move
    assert n["no_ship_terms"] == 1 and n["same_place"] == 1 and n_old["legacy_first_move_dropped"] == 3


def test_fields_and_durations(tmp_path):
    fixed, _ = run(tmp_path, ROWS, TERMS)
    a2 = fixed[0]
    assert a2["ROUT"] == "10-20" and a2["ARRIVAL DATE"] == "2018-05-05 09:30:00"
    assert a2["DURATION"] == "3.12"                              # 3 d 3 h 30 min -> 3 + round(3/24, 2)
    assert a2["STAY DURATION"] == "1.08"                         # 1 d 2 h 30 min -> 1 + round(2/24, 2)
    assert (a2["VESSEL TYPE"], a2["ANTIFOULING PROB"], a2["GWT"], a2["BALLAST DISCHARGE"]) == \
        ("Container", "0.19", "5000.0", "100.0")
    a4 = fixed[1]
    assert a4["ROUT"] == "20-30" and a4["DURATION"] == "3.0"     # from a3's sail (same place as a2)
    assert len(a4) == 20 and fixed[2]["ARRIVAL DATE"] == "2018-05-06 00:00:00"


def test_overlap_not_written(tmp_path):
    rows = [raw_row("7", "a1", "10", "01/05/2018 00:00", "05/05/2018 00:00"),
            raw_row("7", "a2", "20", "04/05/2018 00:00", "06/05/2018 00:00")]   # arrives before a1 sails
    fixed, n = run(tmp_path, rows, TERMS)
    assert fixed == [] and n["overlap_or_gwt"] == 1


def test_dst_gap(tmp_path):                                    # bug 2
    rows = [raw_row("7", "a0", "40", "01/03/2019 00:00", "02/03/2019 00:00"),
            raw_row("7", "a1", "10", "08/03/2019 00:00", "10/03/2019 02:30"),   # in the US daylight-saving gap
            raw_row("7", "a2", "20", "11/03/2019 02:30", "12/03/2019 00:00")]
    fixed, _ = run(tmp_path, rows, TERMS)
    old, _ = run(tmp_path, rows, TERMS, legacy=True)
    assert [r["MOVE REFNO"] for r in fixed] == ["a1", "a2"] and fixed[1]["DURATION"] == "1.0"
    assert [r["MOVE REFNO"] for r in old] == ["a2"] and old[0]["DURATION"] == "0.96"   # read as 03:30 -> 23 h
    assert str(hm.parse_date("10/03/2019 02:30")) == "2019-03-10 02:30:00"


def test_read_terms(tmp_path):
    (tmp_path / "v.csv").write_text("VESSEL.ID,VESSEL.TYPE,GROSS,z,release,random_forest_br\n"
                                    "7,UCC,5000,10,2,\n8,Bulker,9000,NA,NA,300\n9,General Cargo,,1,1,\n")
    t = hm.read_terms(tmp_path / "v.csv")
    assert t == {"7": ("Container", 5000.0, 20.0), "8": ("Bulk", 9000.0, 300.0)}
