# NIS Biofouling Spread Risk — Pipeline & Sensitivity Analysis

This repository contains the code for computing biofouling-mediated
non-indigenous species (NIS) spread risk among Maritime Silk Road countries,
as described in the paper *"Impact of trade facilitation on biofouling-mediated
species spread risk"*.

---

## 1. Repository Structure

```
HON_Code/
│
├── scripts/
│   ├── build_rules_utils.py        # HON rule-building utilities (KLD, sequence I/O)
│   ├── make_HONet_dic.py           # Physical network construction and port risk aggregation
│   ├── translate_clusters.py       # HON cluster translation for visualization
│   └── utils.py                    # Port data loader
│
├── 00_sort_clean_moves.ipynb       # Step 0: sort and clean raw vessel movement records
├── 01_create_traces.ipynb          # Step 1: construct per-ship movement traces
├── 02_Prepare_HON_Input.ipynb      # Step 2: compute P(NIS), write HON input sequences
├── 03_Build-SF-HON.ipynb           # Step 3: build SF-HON from weighted sequences
├── 04_make_HONet_dict.ipynb        # Step 4: aggregate risks to country level
├── 05_translate_clusters-HON.ipynb # Step 5: cluster translation for GIS visualization
├── p_intro.py                      # Standalone P(intro) calculation utilities
├── hon_sf.py                       # Steps 2-4 as one script (no notebooks), with bug fixes (section 6)
├── compare_legacy.py               # Checks hon_sf.py --legacy against the notebook outputs
├── reference/place_coords_fix.csv  # Known coordinate errors in the place table
└── tests/test_hon_sf.py            # Tests, one per fixed bug

data/
└── 2008/  (and other years)

sensitivity_analysis/
├── ballast_water.py                # Sensitivity: ballast discharge input data
└── fouling.py                      # Sensitivity: biofouling model parameters
```

---

## 2. Pipeline Overview

```
Raw vessel movement records (Lloyd's List Intelligence)
        │
        ▼
00_sort_clean_moves.ipynb
   Sort and clean raw movement records → moves_cleaned_2018.txt
        │
        ▼
01_create_traces.ipynb
   Construct continuous per-ship movement traces
        │
        ▼
02_Prepare_HON_Input.ipynb
   Compute P(NIS) = P(alien) × P(intro) × P(estab) per vessel movement
   → write weighted trajectory sequences as HON input
        │
        ▼
03_Build-SF-HON.ipynb
   Build Species Flow Higher-Order Network (SF-HON)
   from weighted sequences using KLD-based rule extension
        │
        ▼
04_make_HONet_dict.ipynb
   Collapse HON nodes to physical ports (sum + normalize)
   → aggregate port-level risks → country-level risks
        │
        ▼
05_translate_clusters-HON.ipynb
   Translate HON clusters back to physical ports for GIS visualization
        │
        ▼
   Country-level NIS spread risk (biofouling / ballast water)
```

---

## 3. Step-by-Step Description

### Step 0 — `00_sort_clean_moves.ipynb`

Sorts raw Lloyd's List Intelligence vessel movement records chronologically
and removes invalid entries. Filters to vessel types relevant to biofouling
and ballast water risk, and removes domestic voyages (same-country port pairs).
Output: `moves_cleaned_2018.txt`.

---

### Step 1 — `01_create_traces.ipynb`

Reconstructs continuous movement traces from the cleaned records by linking
consecutive voyages of the same vessel. This trace information is used in
Step 2 to extract sub-trajectory sequences for HON input.

---

### Step 2 — `02_Prepare_HON_Input.ipynb`

Reads cleaned movement records and computes, for each consecutive port-pair
visited by the same ship, the probability of NIS spread:

```
P(NIS) = P(alien) × P(intro) × P(estab)
```

For each valid movement, all sub-trajectories of length 2 up to `order` are
extracted and written as weighted sequence lines to the HON input file.
Each line encodes a sub-trajectory and its associated risk weight.

**Key outputs:**
- `l_InputForHON_Fouling_env_Eco_16.csv` — weighted trajectories for biofouling
- `l_InputForHON_Ballast_env_Eco_16.csv` — weighted trajectories for ballast water

---

### Step 3 — `03_Build-SF-HON.ipynb`

Builds the Species Flow Higher-Order Network (SF-HON) from the weighted
trajectory sequences produced in Step 2.

The algorithm detects higher-order dependencies in vessel movement patterns
using KL-divergence (KLD): for each source node, it tests whether extending
the rule to a longer historical context is statistically justified. If the
KLD between the extended and current distribution exceeds the threshold
`NewOrder / exp(1 + sum_of_weights)`, the higher-order rule is adopted.
This captures the fact that NIS introduction risk from port i depends not
only on the immediately preceding port but on the full recent trajectory.

Multiple probabilities mapped to the same source–target pair are averaged
(`Aggragate_Probs`). Rules below `MinSupport` are discarded.

**Key outputs:**
- `HONet_Fouling_env_Eco_2018_16.csv` — SF-HON edge list for biofouling
- `HONet_Ballast_env_Eco_2018_16.csv` — SF-HON edge list for ballast water

---

### Step 4 — `04_make_HONet_dict.ipynb`

Collapses the HON (whose nodes encode trajectory contexts) back to the
physical port network. For each physical port pair, all HON edge weights
mapping to that pair are summed (`list_sum_dic`) and then normalized by the
global maximum value. Port-level risk is aggregated using `risk_aggregator`
(complement of the product of non-occurrence probabilities over all incoming
edges). Country-level risk is obtained by summing port risks within each
country.

**Key outputs:**
- `hon_fouling_net_env_Eco_2018.txt` — normalized physical adjacency matrix
- `hon_fouling_agg_risk_env_Eco_2018.txt` — per-port aggregated risk
- `HONlinks_Fouling_env_Eco_2018_16.csv` — edge list for GIS visualization

---

### Step 5 — `05_translate_clusters-HON.ipynb`

Translates the InfoMap clustering output (`.tree`, `.clu`) from
integer-indexed HON nodes back to physical port IDs and names,
for use in GIS visualization.

---

## 4. Key Parameters

### 4.1 Structural / pipeline parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `order` (MaxOrder) | 16 | Maximum trajectory length when building HON rules. Controls how far back in a ship's history the model looks when assessing NIS spread risk. |
| `MinSupport` | 1e-20 | Minimum probability threshold; source–target pairs below this value are discarded from the HON to reduce noise. |
| `r` | `r0` | Ballast water treatment efficacy scenario. `r0` = no treatment (efficacy = 0); `r1`–`r3` represent increasing treatment levels (0.766, 0.992, 0.999). |
| `eco` | `Eco_` | Ecoregion filter mode for P(alien). `Eco_` uses both MEOW and FEOW ecoregions; `sameEco_` uses same-region only; `noEco_` disables the filter (P(alien) = 1 for all pairs). |
| `env` | `env_` | Environmental similarity filter for P(estab). `env_` computes P(estab) from temperature and salinity differences; `noEnv_` sets P(estab) = 1 for all pairs. |

---

### 4.2 P(alien) — probability of being alien

Species are considered non-alien (P(alien) = 0) if the source and destination
ports belong to the same or neighboring marine ecoregion (MEOW), or to
connected freshwater ecoregions (FEOW) for freshwater ports.
Otherwise P(alien) = 1. This component has no continuous free parameters
and is not subject to perturbation in the sensitivity analysis.

---

### 4.3 P(intro) — biofouling introduction probability

Computed by `fouling()` in `02_Prepare_HON_Input.ipynb`:

```
P(intro) = A(m) × poly(stay_d) × exp(−λ × v)
```

where `v = distance / trip_duration` is average sailing speed (km/day),
`stay_d` is port residence time (days), and `poly(stay_d)` is a cubic
polynomial fitted separately for tropical and temperate source ports.

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `λ` (lambda) | 0.008 | Exponential decay rate of biofouling species survival with increasing ship speed. A higher value means organisms die off faster as speed increases. |
| `αT1` | 1.29 × 10⁻⁷ | Cubic coefficient of the **tropical** port-stay polynomial. Controls biofouling accumulation rate at very long stay durations in tropical ports. |
| `αT2` | 8.316 × 10⁻⁵ | Quadratic coefficient of the tropical polynomial. |
| `αT3` | 1.495 × 10⁻² | Linear coefficient of the tropical polynomial. Dominates at short stay durations. |
| `αE1` | 1.40 × 10⁻⁹ | Cubic coefficient of the **temperate** port-stay polynomial. |
| `αE2` | 1.657 × 10⁻⁵ | Quadratic coefficient of the temperate polynomial. |
| `αE3` | 5.194 × 10⁻³ | Linear coefficient of the temperate polynomial. |
| Lat threshold | 0.35 rad (≈ 20°) | Boundary between tropical and temperate port classification. Ports with absolute latitude < 0.35 rad use the tropical polynomial; all others use the temperate polynomial. |
| `A(m)` (`antifouling_p`) | ship-type dependent | Proportion of ships of a given type **without** a functional antifouling system. Read directly from the vessel movement records; not a fixed model parameter and therefore not perturbed in the sensitivity analysis. |

---

### 4.4 P(estab) — establishment probability

Computed by `establishment()` in `02_Prepare_HON_Input.ipynb`:

```
P(estab) = β × exp(−0.5 × [(ΔT / δT)² + (ΔS / δS)²])
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `β` (beta) | 1.5 × 10⁻⁴ | Baseline establishment probability when source and destination environments are identical (ΔT = 0, ΔS = 0). |
| `δT` (delta_T) | 2 °C | Temperature tolerance: the Gaussian kernel standard deviation for temperature difference. At ΔT = δT the establishment probability falls to ~61% of its peak. |
| `δS` (delta_S) | 10 ppt | Salinity tolerance: the Gaussian kernel standard deviation for salinity difference. |

---

## 5. Sensitivity Analyses

Two independent sensitivity analysis scripts are provided under
`sensitivity_analysis/`, both mirroring the full pipeline structure
(data loading → ProbsRaw construction → HON build → risk aggregation).

---

### 5.1 `ballast_water.py` — Ballast water input data sensitivity

**Motivation:** Ballast water discharge volumes outside the United States are
model-predicted values rather than direct measurements. In contrast, US ports
provide true measured discharge volumes. A reviewer raised concerns about
whether prediction uncertainty in these non-US input values materially affects
the country-level risk rankings. This analysis directly quantifies that
sensitivity.

**What is perturbed:** The ballast water discharge volume of each individual
vessel movement originating from a non-US port is independently multiplied
by a random noise factor:

```
bw_i = bw_i_baseline × max(Normal(1.0, σ), 0)
```

US-source port movements are held fixed at their measured values.
Each vessel movement receives an independent noise draw (per-record
perturbation), reflecting that prediction errors are independent across vessels.

**Metrics reported** (over N Monte Carlo iterations per σ value):
- Spearman rank correlation between baseline and perturbed country risk rankings
- Top-10 country stability: % of baseline top-10 countries retained
- Top-Q country stability: % of baseline top-25% countries retained

**Results (N = 100):**

| σ | Spearman (mean ± std) | Spearman min | Top-10 mean | Top-10 min | Top-Q mean | Top-Q min |
|---|----------------------|--------------|-------------|------------|------------|-----------|
| 0.2 (±20%) | 0.9986 ± 0.0003 | 0.9974 | 99.4% | 90.0% | 97.6% | 93.3% |
| 0.4 (±40%) | 0.9954 ± 0.0016 | 0.9849 | 96.3% | 80.0% | 95.5% | 91.1% |

Country risk rankings are highly stable under both perturbation levels.
At ±40%, Spearman correlation remains above 0.98 in all runs and the
top-25% country set is retained at 95.5% on average, indicating that
prediction uncertainty in non-US ballast discharge volumes does not
materially affect the conclusions.

---

### 5.2 `fouling.py` — Biofouling model parameter sensitivity

**Motivation:** The continuous parameters of the P(NIS) biofouling model
(polynomial coefficients, decay rate, establishment parameters) were
calibrated from external empirical literature. Their uncertainty could affect
the risk rankings. This analysis quantifies robustness to that parameter
uncertainty at two levels: 32 BRI countries and country pairs.

The script runs two analyses:

#### Part 1 — Monte Carlo parameter perturbation

In each iteration, all continuous parameters are independently perturbed:

```
param_i = param_i_baseline × max(Normal(1.0, σ), 0)
```

Each parameter receives an independent noise draw **once per iteration**,
shared across all vessel movements in that run. This reflects the
interpretation that within one simulation the true parameter value is fixed
but unknown; uncertainty is across simulations.

Parameters perturbed, with their baseline values:

**P(intro) parameters:**

| Parameter | Baseline value | Meaning |
|-----------|---------------|---------|
| `λ` | 0.008 | Speed-dependent survival decay rate |
| `αT1` | 1.29 × 10⁻⁷ | Tropical polynomial cubic coefficient |
| `αT2` | 8.316 × 10⁻⁵ | Tropical polynomial quadratic coefficient |
| `αT3` | 1.495 × 10⁻² | Tropical polynomial linear coefficient |
| `αE1` | 1.40 × 10⁻⁹ | Temperate polynomial cubic coefficient |
| `αE2` | 1.657 × 10⁻⁵ | Temperate polynomial quadratic coefficient |
| `αE3` | 5.194 × 10⁻³ | Temperate polynomial linear coefficient |

**P(estab) parameters:**

| Parameter | Baseline value | Meaning |
|-----------|---------------|---------|
| `β` | 1.5 × 10⁻⁴ | Baseline establishment probability |
| `δT` | 2 °C | Temperature tolerance standard deviation |
| `δS` | 10 ppt | Salinity tolerance standard deviation |

P(alien) is not perturbed (binary structural classification, no continuous
free parameters). The latitude threshold (0.35 rad) is not perturbed in the
Monte Carlo — see Part 2.

**Metrics reported** (over N Monte Carlo iterations per σ value), evaluated
at three levels:

[A] **Country-level — 32 BRI countries:**
- Spearman rank correlation
- Top-5 stability: % of baseline top-5 BRI countries retained
- Top-Q stability: % of baseline top-25% BRI countries retained
- Category agreement rate: % of countries retaining the same risk category
  (very high / high / medium / low) as baseline ¹

[B] **Pair-level — all global country pairs:**
- Spearman rank correlation
- Top-20 stability
- Top-Q stability

[C] **Pair-level — 32 BRI country pairs only** (both source and destination
belong to the 32 BRI countries):
- Spearman rank correlation
- Top-10 stability
- Top-Q stability
- Category agreement rate ¹

¹ Category agreement rate is computed on HON-normalized output values using
the paper's classification thresholds (0.75 / 0.50 / 0.25). Since the
HON-normalized values differ from raw P(NIS), this metric is indicative
rather than strict and should be interpreted with caution.

**Results (N = 100):**

Baseline top-5 BRI countries: India, Indonesia, Iran, Italy, Philippines

**[A] Country-level — BRI (32 countries)**

| σ | Spearman mean ± std | Spearman min | Top-5 mean | Top-5 min | Top-Q mean | Top-Q min | Cat. agree mean | Cat. agree min |
|---|---------------------|--------------|------------|-----------|------------|-----------|----------------|----------------|
| 0.1 (±10%) | 0.9956 ± 0.0029 | 0.9879 | 98.2% | 80.0% | 97.1% | 87.5% | 97.9% | 84.4% |
| 0.2 (±20%) | 0.9901 ± 0.0068 | 0.9615 | 91.4% | 80.0% | 94.4% | 75.0% | 94.2% | 68.8% |
| 0.3 (±30%) | 0.9826 ± 0.0172 | 0.8592 | 87.4% | 60.0% | 91.8% | 75.0% | 89.5% | 56.2% |
| 0.4 (±40%) ² | 0.9657 ± 0.0757 | 0.2658 | 83.0% | 0.0% | 87.1% | 12.5% | 86.5% | 56.2% |

**[B] Pair-level — Global (4,492 pairs)**

| σ | Spearman mean ± std | Spearman min | Top-20 mean | Top-20 min | Top-Q mean | Top-Q min |
|---|---------------------|--------------|-------------|------------|------------|-----------|
| 0.1 (±10%) | 0.9983 ± 0.0021 | 0.9878 | 94.2% | 90.0% | 97.6% | 91.9% |
| 0.2 (±20%) | 0.9930 ± 0.0085 | 0.9433 | 91.3% | 70.0% | 95.3% | 83.5% |
| 0.3 (±30%) | 0.9831 ± 0.0237 | 0.8174 | 88.0% | 60.0% | 92.9% | 71.2% |
| 0.4 (±40%) ² | 0.9605 ± 0.0989 | 0.0545 | 82.0% | 0.0% | 89.3% | 28.7% |

**[C] Pair-level — BRI (471 pairs)**

| σ | Spearman mean ± std | Spearman min | Top-10 mean | Top-10 min | Top-Q mean | Top-Q min | Cat. agree mean | Cat. agree min |
|---|---------------------|--------------|-------------|------------|------------|-----------|----------------|----------------|
| 0.1 (±10%) | 0.9989 ± 0.0013 | 0.9918 | 95.7% | 80.0% | 98.3% | 93.2% | 97.9% | 92.8% |
| 0.2 (±20%) | 0.9955 ± 0.0058 | 0.9590 | 92.2% | 80.0% | 96.0% | 88.0% | 95.9% | 86.6% |
| 0.3 (±30%) | 0.9893 ± 0.0152 | 0.8835 | 87.9% | 70.0% | 93.8% | 78.6% | 94.1% | 79.2% |
| 0.4 (±40%) ² | 0.9725 ± 0.0868 | 0.1463 | 81.0% | 10.0% | 90.2% | 25.6% | 92.7% | 79.6% |

² At σ = 0.4, it may produce a NaN Spearman correlation because
extreme parameter draws caused all biofouling probabilities to collapse to
zero, yielding a degenerate constant risk vector. This run is excluded from
mean and std via nanmean. The σ = 0.4 scenario represents a physically
unrealistic degree of parameter uncertainty and is reported for completeness.

At σ = 0.1 (±10%), all metrics are highly stable across all three levels.
At σ = 0.2 (±20%), Spearman remains above 0.96 and Top-Q above 75% in all
runs. Results degrade gradually with increasing σ, but rankings remain
broadly consistent even at σ = 0.3 (±30%), confirming that the conclusions
are robust to plausible parameter uncertainty.

#### Part 2 — Latitude threshold scenario analysis

The tropical/temperate boundary is a structural categorical assumption rather
than a continuously uncertain quantity, so it is assessed through deterministic
scenario runs rather than Monte Carlo perturbation. All model parameters are
held at their baseline values. Results are evaluated at all three levels (A/B/C).

| Scenario | Threshold | Rationale |
|----------|-----------|-----------| 
| Baseline | 0.35 rad ≈ 20° | Original model value |
| Scenario A | 15° (0.262 rad) | Narrower tropical zone |
| Scenario B | 25° (0.436 rad) | Wider tropical zone |

**Results:**

Baseline top-5 BRI countries (~20°): India, Indonesia, Iran, Italy, Philippines

**[A] Country-level — BRI**

| Scenario | Spearman | Top-5 overlap | Top-Q overlap | Cat. agreement |
|----------|----------|--------------|---------------|----------------|
| 15° (narrower tropical) | 1.0000 | 100.0% | 100.0% | 100.0% |
| 25° (wider tropical) | 1.0000 | 100.0% | 100.0% | 100.0% |

**[B] Pair-level — Global**

| Scenario | Spearman | Top-20 overlap | Top-Q overlap |
|----------|----------|----------------|---------------|
| 15° (narrower tropical) | 1.0000 | 100.0% | 99.8% |
| 25° (wider tropical) | 1.0000 | 95.0% | 100.0% |

**[C] Pair-level — BRI**

| Scenario | Spearman | Top-10 overlap | Top-Q overlap | Cat. agreement |
|----------|----------|----------------|---------------|----------------|
| 15° (narrower tropical) | 1.0000 | 100.0% | 100.0% | 100.0% |
| 25° (wider tropical) | 1.0000 | 100.0% | 100.0% | 100.0% |

The latitude threshold has negligible effect on risk rankings at all levels.
Country-level and BRI pair-level results are identical across both scenarios.
At the global pair level, the 25° scenario shows a 5% reduction in Top-20
overlap, driven by a small number of border-zone ports switching classification,
but Spearman remains perfect and Top-Q overlap is 100%. This confirms that
the tropical/temperate classification boundary is not a sensitive assumption
for the study's conclusions.

---

## 6. Script version: `hon_sf.py` (2026-11-17)

`Code/HON_Code/hon_sf.py` runs steps 2–4 (`02_Prepare_HON_Input`, `03_Build-SF-HON`, `04_make_HONet_dict`, with
`scripts/build_rules_utils.py` and `scripts/make_HONet_dic.py`) as one script. It needs only numpy. The notebooks
are unchanged.

```
python Code/HON_Code/hon_sf.py run --moves data/moves/moves_cleaned_2018.txt \
    --places data/Places_allportdata_mergedSept2017.csv --port-country data/stay/portToCountry.csv \
    --out <dir>
```

The steps can also be run one by one (`prepare`, `build`, `aggregate`, same `--out`). Options: `--eco {eco,same,none}`
(P(alien): same or neighbouring ecoregion, same ecoregion only, none; default `eco`, the `Eco_` setting of
section 4.1), `--no-env`, `--regime {r0,r1,r2,r3}`, `--order`, `--max-order`, `--min-support`, `--coords-fix`,
`--scale-from` (see 6.4) and `--legacy`.

Outputs: `input_{ballast,fouling}.txt` (weighted sub-trajectories), `rules_*.txt`, `honet_*.csv`, `net_*.txt`
(physical port network), `links_*.csv`, `port_risk_out_*.txt` (the agg_risk files of step 4), `port_risk.csv`
(outgoing and incoming port risk for both vectors), `trans_honet_*.{net,dic}` (for clustering) and
`hon_summary.json` (counts and settings).

Not converted: step 0 (raw moves to `moves_cleaned`), step 1 (its output is not read by later steps), step 5 (GIS
figures), and the trade and stay adjustments of step 2 (switched off there).

### 6.1 Check against the notebooks

`--legacy` reproduces the notebooks as they were run for the 2018 outputs (step 2 with the moves filter: stays under
20 days, nine cargo types, no domestic legs), bugs included. On `moves_cleaned_2018.txt`, all 16 outputs have the
same lines as the notebook outputs in `data/2018`, with values within 7e-13 (relative):

```
python Code/HON_Code/hon_sf.py run --legacy ... --out <dir>
python Code/HON_Code/compare_legacy.py <dir> data/2018
```

### 6.2 Fixed without `--legacy`

| # | Step | Problem in the notebooks | Fix | 2018 (legacy run) |
|---|------|--------------------------|-----|-------------------|
| 1 | 2 | Tropical test compares latitude in degrees with 0.35, a value in radians (≈ 20°, section 4.3): only ports within 0.35° of the equator use the tropical polynomial | Latitude in radians < 0.35 | — |
| 2 | 2 | Ballast P(intro) uses the trip duration of the last leg in the file for every leg (a variable left over from a loop) | Each leg's own trip duration | — |
| 3 | 2 | Dropping a move (domestic leg, stay ≥ 20 days) joins the moves around it: A→B→C becomes A→C, with B→C's duration | A leg is two consecutive moves. Domestic legs, legs after a stay of 20 days or more and legs under one hour get no risk record, but their ports stay in the ship's trajectory | — |
| 4 | 4 | Collapsing higher-order nodes to ports sums probabilities; the sum can exceed 1 | P(at least one) = 1 − ∏(1 − r) | 822 ballast and 82 fouling links > 1; 231 port risks outside [0, 1] |
| 5 | 2 | Biofouling uses the stay at the destination port; step 0 documents it as the stay at the source port | Source port's stay; the 20-day limit applies to it | — |
| 6 | 2 | The maximum used to scale P(NIS) takes P(alien) = 0 for the same or a neighbouring ecoregion; the written values for the same ecoregion only | One rule (`--eco`) for both | 10 fouling values > 1 |
| 7 | 3 | The first observation of a new target under an existing source is stored twice and weighs double in the mean | Stored once | — |
| 8 | 4 | 1 − ∏(1 − r) loses every r below ~1e-16 | Computed as −expm1(Σ log1p(−r)) | — |
| 9 | 4 | The ballast node dictionary for clustering is built from the fouling network's nodes | Each network's own nodes | 41,486 ballast edges not written |
| 10 | data | Nine places in `Places_allportdata_mergedSept2017.csv` have a sign or digit error (Panama Canal at 79.68 E, Vancouver Anch. at 123.3 E, ...); they enter leg distances and the link table | Corrected from `reference/place_coords_fix.csv` while the table still holds the wrong values | — |

Effect on 2018 (Spearman correlation of port risks, legacy vs fixed, ballast/fouling): with `--eco eco`, outgoing
0.77/0.63 and incoming 0.83/0.67; with `--eco same`, outgoing 0.87/0.83 and incoming 0.96/0.88. The large ports'
risks are close to 1 in both, so their order at the top is not informative.

### 6.3 Not changed

- The KLD test uses p · exp(p/q) and the threshold order / exp(1 + Σp), as in the original code (the HON paper uses
  log and counts). Kept as a method choice.
- Step 0 (not converted) drops the first move of every ship when grouping moves by ship, and in its 1997–2012 branch
  the row index advances only on valid rows. Files made by step 0 (`moves_cleaned_*`) carry the first problem.
- `sensitivity_analysis/fouling.py` and `ballast_water.py` compare latitude in degrees with 0.35 rad, take the
  destination stay, and drop moves before pairing (problems 1, 5 and 3 above). The latitude scenarios of section 5.2
  (Part 2) therefore moved the tropical boundary only within about 0.5° of the equator.

### 6.4 Comparing years

Step 2 divides every P(NIS) by the largest in the file, so years run separately are on different scales, and the
HON rule threshold depends on that scale. `--scale-from <hon_summary.json ...>` uses the largest maxima of the given
runs, so that several years share one scale.

### 6.5 Tests

```
python -m pytest Code/HON_Code/tests
```
