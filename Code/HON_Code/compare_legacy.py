"""Check hon_sf.py --legacy against the notebook outputs of 2018 (data/2018).

    python Code/HON_Code/compare_legacy.py <dir of hon_sf.py run --legacy> data/2018


Lines are compared after sorting (the notebooks iterate sets, so their line order is arbitrary). Where a file differs,
lines are matched on their non-numeric key and the numbers compared with relative tolerance 1e-9 (sums taken in a
different order differ in the last bits)."""
import math
import sys
from pathlib import Path

NEW = Path(sys.argv[1])
OLD = Path(sys.argv[2])
PAIRS = {"input_ballast.txt": "l_InputForHON_Ballast_env_Eco_16.csv",
         "input_fouling.txt": "l_InputForHON_Fouling_env_Eco_16.csv",
         "rules_ballast.txt": "Rules_Ballast_env_Eco_2018_16.csv",
         "rules_fouling.txt": "Rules_Fouling_env_Eco_2018_16.csv",
         "honet_ballast.csv": "HONet_Ballast_env_Eco_2018_16.csv",
         "honet_fouling.csv": "HONet_Fouling_env_Eco_2018_16.csv",
         "net_ballast.txt": "hon_ballast_net_env_Eco_2018.txt",
         "net_fouling.txt": "hon_fouling_net_env_Eco_2018.txt",
         "port_risk_out_ballast.txt": "hon_ballast_agg_risk_env_Eco_2018.txt",
         "port_risk_out_fouling.txt": "hon_fouling_agg_risk_env_Eco_2018.txt",
         "links_ballast.csv": "HONlinks_Ballast_env_Eco_2018_16.csv",
         "links_fouling.csv": "HONlinks_Fouling_env_Eco_2018_16.csv",
         "trans_honet_ballast.dic": "trans_HONet_Ballast_env_Eco_2018_16.dic",
         "trans_honet_fouling.dic": "trans_HONet_Fouling_env_Eco_2018_16.dic",
         "trans_honet_ballast.net": "trans_HONet_Ballast_env_Eco_2018_16.net",
         "trans_honet_fouling.net": "trans_HONet_Fouling_env_Eco_2018_16.net"}


def split(line):
    toks = line.replace(",", " ").split()
    key, nums = [], []
    for t in toks:
        try:
            v = float(t)
            nums.append(v) if ("." in t or "e" in t) else key.append(t)
        except ValueError:
            key.append(t)
    return tuple(key), nums


for new, old in PAIRS.items():
    a = sorted(l.rstrip("\n").rstrip() for l in open(NEW / new))
    b = sorted(l.rstrip("\n").rstrip() for l in open(OLD / old))
    if a == b:
        print(f"{new:28s} identical ({len(a):,} lines)")
        continue
    ka, kb = {}, {}
    for l in a:
        k, v = split(l)
        ka.setdefault(k, []).append(v)
    for l in b:
        k, v = split(l)
        kb.setdefault(k, []).append(v)
    only_a, only_b = set(ka) - set(kb), set(kb) - set(ka)
    worst, nbad = 0.0, 0
    for k in set(ka) & set(kb):
        for x, y in zip(sorted(ka[k]), sorted(kb[k])):
            for p, q in zip(x, y):
                r = abs(p - q) / max(abs(p), abs(q), 1e-300)
                worst = max(worst, r)
                nbad += r > 1e-9
    print(f"{new:28s} lines {len(a):,} vs {len(b):,}; keys only new {len(only_a)}, only old {len(only_b)}; "
          f"max rel diff {worst:.2e}; values off by >1e-9: {nbad}")
