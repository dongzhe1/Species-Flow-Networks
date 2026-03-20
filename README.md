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
└── p_intro.py                      # Standalone P(intro) calculation utilities

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

² At σ = 0.4, 1 out of 100 runs produced a NaN Spearman correlation because
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
