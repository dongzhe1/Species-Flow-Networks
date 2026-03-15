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
the country-level risk rankings. This analysis quantifies robustness to
that parameter uncertainty.

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

**Metrics reported:** Same as `ballast_water.py`:
- Spearman rank correlation
- Top-10 stability
- Top-Q stability

**Results (N = 100):**

Baseline top-10 countries: Faroe Islands, Germany, Greece, India, Indonesia,
Iran, Italy, Malaysia, Philippines, U.S.A.

| σ | Spearman (mean ± std) | Spearman min | Top-10 mean | Top-10 min | Top-Q mean | Top-Q min |
|---|----------------------|--------------|-------------|------------|------------|-----------|
| 0.2 (±20%) | 0.9931 ± 0.0076 | 0.9458 | 91.5% | 70.0% | 94.6% | 86.7% |
| 0.4 (±40%) | — ¹ | — | 82.5% | 10.0% | 88.6% | 33.3% |

¹ Spearman correlation is undefined (NaN) at σ = 0.4 because extreme parameter
draws in some iterations cause all biofouling probabilities to collapse to zero,
producing a degenerate constant risk vector for which rank correlation is
undefined. Top-N and Top-Q overlap metrics are unaffected and remain computable.
At σ = 0.2, results are stable: Spearman exceeds 0.94 in all runs and the
top-25% country set is retained at 94.6% on average. The σ = 0.4 scenario
represents a physically unrealistic degree of parameter uncertainty (each
parameter independently varying by up to ±40%), and is reported for
completeness.

#### Part 2 — Latitude threshold scenario analysis

The tropical/temperate boundary is a structural categorical assumption rather
than a continuously uncertain quantity, so it is assessed through deterministic
scenario runs rather than Monte Carlo perturbation. All model parameters are
held at their baseline values.

| Scenario | Threshold | Rationale |
|----------|-----------|-----------|
| Baseline | 0.35 rad ≈ 20° | Original model value |
| Scenario A | 15° (0.262 rad) | Wider tropical zone |
| Scenario B | 25° (0.436 rad) | Narrower tropical zone |

Each scenario produces a single deterministic country risk ranking, compared
to baseline using Spearman correlation and Top-10 / Top-Q overlap.

**Results:**

Baseline top-10 (~20°): Faroe Islands, Germany, Greece, India, Indonesia,
Iran, Italy, Malaysia, Philippines, U.S.A.

| Scenario | Spearman | Top-10 overlap | Top-Q overlap |
|----------|----------|---------------|---------------|
| 15° (wider tropical zone) | 1.0000 | 100.0% | 100.0% |
| 25° (narrower tropical zone) | 1.0000 | 100.0% | 100.0% |

The tropical/temperate boundary has no effect on the country-level risk
rankings: both alternative thresholds produce identical top-10 and top-Q
country sets, with perfect Spearman correlation against the baseline.
This confirms that the categorical classification of ports as tropical or
temperate is not a sensitive assumption for the conclusions of the study.
