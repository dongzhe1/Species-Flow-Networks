# NIS Spread Risk on Ship Trajectories (SF-HON)

This repository computes ballast-water and biofouling-mediated non-indigenous species (NIS) spread risk on a
species-flow higher-order network (SF-HON) built from ship movements. It is used in the paper *"Impact of trade
facilitation on biofouling-mediated species spread risk"* and builds on the SF-HON code of Saebi et al.
([msaebi1993/Species-Flow-Networks](https://github.com/msaebi1993/Species-Flow-Networks)).

Since 2026-11-17 the pipeline is two Python scripts with fixes to the original code (section 5). The Jupyter
notebooks and their helper modules (`scripts/`) were removed; they are in tag `v1.0.0`.

---

## 1. Repository structure

```
Code/
├── HON_Code/
│   ├── hon_moves.py                # Step 0: raw Lloyd's moves -> moves_cleaned
│   ├── hon_sf.py                   # Steps 2-4: risk per leg -> SF-HON -> port network and port risks
│   ├── compare_legacy.py           # Compares hon_sf.py --legacy with the notebook outputs (data/2018)
│   ├── p_intro.py                  # Older standalone step 2 for ballast treatment scenarios (not updated, section 7)
│   ├── reference/
│   │   └── place_coords_fix.csv    # Known coordinate errors in the place table
│   └── tests/                      # One test per fixed problem
└── sensitivity_analysis/           # Not updated (section 7)
    ├── ballast_water.py
    └── fouling.py

data/
└── Places_allportdata_mergedSept2017.csv   # Lloyd's places with temperature, salinity, ecoregions
```

Other inputs are not in the repository: the raw Lloyd's moves (e.g. `moves_1819.lst`), ship terms (type, gross
tonnage, ballast discharge per Lloyd's vessel ID) and `data/stay/portToCountry.csv` (place ID to country).

---

## 2. Pipeline

```
raw Lloyd's moves + ship terms
        │  hon_moves.py                              (step 0)
        ▼
moves_cleaned_<year>.txt: one row per move, with previous port, trip and stay, ship terms
        │  hon_sf.py prepare                         (step 2)
        ▼
legs -> P(NIS) per leg for ballast water and biofouling -> weighted sub-trajectories (orders 2-16)
        │  hon_sf.py build                           (step 3)
        ▼
rules and higher-order network (SF-HON)
        │  hon_sf.py aggregate                       (step 4)
        ▼
physical port network, port risks (outgoing and incoming), link table, files for clustering
```

```
python Code/HON_Code/hon_moves.py --raw moves_1819.lst --vessels <ship terms .csv> --out data/moves/moves_cleaned_2018.txt
python Code/HON_Code/hon_sf.py run --moves data/moves/moves_cleaned_2018.txt \
    --places data/Places_allportdata_mergedSept2017.csv --port-country data/stay/portToCountry.csv --out <dir>
```

- `hon_moves.py`:
  - The ship terms file is a CSV with `VESSEL.ID`, `VESSEL.TYPE` (a type name or a Lloyd's type code), `GROSS`, and
    the ballast discharge as `ballast_discharge`, as `z` and `release` (their product), or as `random_forest_br`.
  - `--terms-from <moves_cleaned file>` reuses the ship terms of an existing cleaned file instead.
- `hon_sf.py`:
  - The steps can be run one by one (`prepare`, `build`, `aggregate`, with the same `--out`).
  - Options: `--eco`, `--no-env`, `--regime`, `--order`, `--max-order`, `--min-support`, `--coords-fix`,
    `--scale-from` and `--legacy`.
- Both scripts have `--legacy`, which reproduces the notebooks with their problems.

Outputs of `hon_sf.py`:

| File | Content |
|------|---------|
| `input_{ballast,fouling}.txt` | Weighted sub-trajectories (HON input) |
| `rules_*.txt` | HON rules |
| `honet_*.csv` | Higher-order network |
| `net_*.txt` | Physical port network |
| `links_*.csv` | Link table with port names, coordinates and countries |
| `port_risk_out_*.txt` | Outgoing port risk, in the format of the old step 4 |
| `port_risk.csv` | Outgoing and incoming port risk, ballast and fouling |
| `trans_honet_*.{net,dic}` | Node dictionary and network for clustering |
| `hon_summary.json` | Counts and settings |

Not converted: step 1 (per-ship traces; not read by later steps), step 5 (cluster translation for GIS figures) and
the trade and stay adjustments of step 2 (switched off in the notebooks).

---

## 3. Model

```
P(NIS) = P(alien) × P(intro) × P(estab), per leg (two consecutive moves of a ship at ports in different countries)
```

- **P(alien)**:
  - `--eco eco` (default): 0 when the source and destination are in the same or neighbouring marine ecoregions
    (MEOW; FEOW for freshwater ports), else 1.
  - `--eco same`: 0 only for the same ecoregion.
  - `--eco none`: 1 everywhere.
- **P(intro), ballast water**: (1 − efficacy) × (1 − exp(−3.22×10⁻⁶ × ballast discharge)) × exp(−0.02 × trip days).
  The efficacy is set by `--regime` (r0–r3: 0, 0.766, 0.992, 0.999).
- **P(intro), biofouling**: A(m) × poly(stay) × exp(−0.008 × v).
  - A(m): the share of ships of the type without working antifouling.
  - stay: days at the source port.
  - v: leg distance / trip days (km/day).
  - poly: a cubic in stay. A source port is tropical when |latitude| < 0.35 rad (≈ 20°).

  | Coefficient | Tropical | Temperate |
  |---|---|---|
  | Cubic | 1.29×10⁻⁷ | 1.4×10⁻⁹ |
  | Quadratic | −8.3165×10⁻⁵ | −1.6566×10⁻⁵ |
  | Linear | 1.495187×10⁻² | 5.19377×10⁻³ |
- **P(estab)**: 1.5×10⁻⁴ × exp(−0.5 × [(ΔT/2 °C)² + (ΔS/10 ppt)²]), from the yearly mean temperature and the
  salinity of the two ports (`--no-env`: 1).

Which legs get a risk record, and how values are scaled:
- Legs get a risk record for nine cargo types (Auto, Container, Bulk, Tanker, Chemical, Oil, General, Liquified-Gas,
  Refrigerated-Cargo), between countries, after a stay under 20 days, and with a trip of at least one hour.
- Other legs keep their ports in the ship's trajectory but get no record.
- Values are divided by the largest in the run. Use `--scale-from` to share one scale across runs (section 8).

HON (step 3):
- Sub-trajectories of length 2 to 16 carry their leg's P(NIS); a source's distribution is the mean per target, and
  MinSupport is 1e-20.
- A higher-order rule is kept when KLD(extended, current) > order / exp(1 + Σp of the extended source), with
  KLD = Σ p · exp(p/q). This is Saebi's form; the HON paper uses log and counts.

Aggregation (step 4):
- Higher-order nodes are collapsed to physical ports as P(at least one) = 1 − ∏(1 − r).
- Port risk is the same combination over a port's outgoing links (`out`) or incoming links (`in`).

---

## 4. Check against the notebooks

- `hon_moves.py --legacy` rebuilds `moves_cleaned_2018.txt` from `moves_1819.lst` byte for byte (ship terms taken
  from that file).
- `hon_sf.py --legacy` reproduces all 16 outputs of the notebooks in `data/2018` (local outputs, not in the
  repository): same lines, values within 7×10⁻¹³ (relative). Check with:

  ```
  python Code/HON_Code/compare_legacy.py <dir of hon_sf.py run --legacy> data/2018
  ```

---

## 5. Problems fixed

Numbers are from the 2018 legacy runs.

| # | Step | Problem in the notebooks | Fix | 2018 |
|---|------|--------------------------|-----|------|
| 1 | 2 | The tropical test compared latitude in degrees with 0.35, a value in radians: only ports within 0.35° of the equator used the tropical polynomial | Latitude in radians < 0.35 | — |
| 2 | 2 | Ballast P(intro) used the trip duration of the last leg in the file for every leg (a variable left over from a loop) | Each leg's own trip | — |
| 3 | 2 | Dropping a move (domestic leg, stay ≥ 20 days) joined the moves around it: A→B→C became A→C, with B→C's duration | A leg is two consecutive moves; legs without a record keep their ports in the trajectory | — |
| 4 | 4 | Collapsing higher-order nodes summed their risks and divided by a fixed 36.51 (later switched off), so values could exceed 1 before 1 − ∏(1 − r) was applied | P(at least one) = 1 − ∏(1 − r) (≈ the sum for small values, never above 1) | largest sums 45.04 (ballast), 8.93 (fouling); with the division off, 822 ballast and 82 fouling links > 1 and 231 port risks outside [0, 1] |
| 5 | 2 | Biofouling used the stay at the destination port; step 0 documents it as the stay at the source port | Source port's stay; the 20-day limit applies to it | — |
| 6 | 0 | Grouping moves by ship skipped each ship's first move, so the next one, without a previous move, was not written either | All moves kept | 63,370 ships lost their second move; 47,254 rows (1.3%) missing, 80% in May 2018 |
| 7 | 2 | The maximum used for scaling took P(alien) = 0 for the same or neighbouring ecoregions; the written values for the same ecoregion only | One rule (`--eco`) for both | 10 fouling values > 1 |
| 8 | 3 | The first observation of a new target under an existing source was stored twice and weighed double in the mean | Stored once | — |
| 9 | 0 | Dates went through local time; the 2018 run was in a US time zone, so times between 02:00 and 03:00 on 10 March 2019 moved one hour later | Dates read as written | 644 rows |
| 10 | 4 | 1 − ∏(1 − r) lost every r below about 10⁻¹⁶ | −expm1(Σ log1p(−r)) | — |
| 11 | 4 | The ballast node dictionary for clustering was built from the fouling network's nodes | Each network's own nodes | 41,486 ballast edges not written |
| 12 | data | Nine places in `Places_allportdata_mergedSept2017.csv` have a sign or digit error (Panama Canal at 79.68 E, Vancouver Anch. at 123.3 E, ...) | Corrected from `reference/place_coords_fix.csv` while the table holds the wrong values | — |

---

## 6. Effect on 2018

Port risks of the legacy and fixed runs of `hon_sf.py` on the original `moves_cleaned_2018.txt` (before fixes 6 and
9), Spearman correlation, ballast/fouling:

| P(alien) rule | Outgoing | Incoming |
|---|---|---|
| `--eco eco` | 0.77/0.63 | 0.83/0.67 |
| `--eco same` | 0.87/0.83 | 0.96/0.88 |

The largest ports' risks are close to 1 in both runs, so their order at the top is not informative.

---

## 7. Code not updated

- **`p_intro.py`**: an older standalone version of step 2 for ballast treatment scenarios (s1–s8). It still has
  problems 1, 2, 5 and 7.
- **`sensitivity_analysis/fouling.py` and `ballast_water.py`**: they still compare latitude in degrees with 0.35 rad,
  use the destination stay, and drop moves before pairing (problems 1, 5 and 3). Their results (README of
  `v1.0.0`, section 5) were computed with these problems. The latitude scenarios there (15°, 20° and 25°) moved the
  tropical boundary only within about 0.5° of the equator.

---

## 8. Comparing years

Step 2 divides every P(NIS) by the largest in the run, so years run separately are on different scales. The HON rule
threshold also depends on that scale. `--scale-from <hon_summary.json ...>` uses the largest maxima of the given runs,
so that several years share one scale.

---

## 9. Tests

```
python -m pytest Code/HON_Code/tests
```
