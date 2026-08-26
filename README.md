# baria_muscle

Analysis code for the **BARIA gut–muscle project**: associations between the
faecal gut microbiome (shotgun metagenomics) and skeletal muscle mass following
bariatric surgery, with **fat-free mass index (FFMI)** as the primary outcome.

The differentiating angle relative to the wider bariatric–microbiome literature
is the **gut–muscle axis**: skeletal muscle as a tissue-level link between
post-surgical microbiome shifts and metabolic outcomes.

- **Cohort:** BARIA longitudinal bariatric-surgery cohort (Amsterdam UMC)
- **Microbiome:** shotgun metagenomics — MetaPhlAn taxonomic profiles and HUMAnN functional-pathway profiles
- **Clinical:** clinical, anthropometric and BIA-derived measures
- **Language:** R

The code currently covers data cleaning, a baseline descriptive table, the
baseline microbiome descriptive analyses (composition and diversity), a
longitudinal mixed-model layer (`4b`) relating baseline species abundance to
FFMI trajectory after surgery, and a functional-pathway layer (`5a`, HUMAnN)
that runs both cross-sectional and longitudinal models.

---

## Repository structure

```
baria_muscle/
└── scripts/
    ├── 0a_datacleaning.R                 # Clinical + microbiome cleaning, derived variables, cohort definition
    ├── 1a_tableone_v0.R                   # Table 1, stratified by baseline FFMI status
    ├── 2a_mb_alphadiversity_v0.R         # Alpha diversity (Shannon, Simpson, observed richness)
    ├── 2b_mb_betadiversity_v0.R          # Beta diversity (Bray–Curtis PCoA)
    ├── 2c_mb_composition_v0.R            # Compositional bar plots (top 20 species, by baseline FFMI group)
    ├── 4b_lmm_mb_ffmi_trajectories.R     # Mixed models: baseline species abundance × FFMI trajectory
    └── 5a_humann_pathways.R              # Functional pathways (HUMAnN): cleaning + cross-sectional + longitudinal models
```

Only `scripts/` is tracked. At run time the scripts read from and write to a
`data/` tree, and the plotting scripts write figures to `graphs/` and Table 1
to `tables/` (paths are relative to the repo root, created alongside the repo):

```
data/raw_data/         # inputs (not tracked — see Data)
data/processed_data/   # written by 0a
graphs/                # figures, in per-analysis subfolders (alphadiversity, betadiversity, composition, lmm_species_ffmi)
tables/                # Table 1 (t1_v0_mb_matrix.csv)
```

### Order of execution

`0a` runs first — it writes the processed datasets every other script reads.
`1a` and the `2*` descriptive scripts are independent of one another. `4b` reads
the processed long dataset and the cleaned microbiome object from `0a`. `5a`
reads the raw HUMAnN profiles together with the processed data from `0a`, and is
otherwise independent of the taxonomic scripts.

> **Path note:** `0a` currently writes datestamped files (`260824_…`), while the
> downstream scripts still read an earlier datestamp (`260818_…`). Align the
> read paths (or the write datestamp) before running the pipeline end-to-end.

---

## Data

**Inputs** (placed in `data/raw_data/`, not tracked in the repository):

| File | Description |
|------|-------------|
| `BARIA.clinical.<date>.RDS` | Clinical, anthropometric and BIA data |
| `ps.BARIA.metaphlan.<...>.RDS` | MetaPhlAn shotgun phyloseq object (per-sample relative abundances) |
| `BARIA.humann4.profiles.<...>.RDS` | HUMAnN 4.0 community functional-pathway profiles (per-sample pathway abundances, CPM) |

**Outputs** (written by `0a` to `data/processed_data/`):

| File | Description |
|------|-------------|
| `<date>_BARIA_muscle_long.{RDS,csv}` | Long clinical dataset (one row per participant × visit) |
| `<date>_BARIA_muscle_wide.{RDS,csv}` | Wide clinical dataset (one row per participant) |
| `<date>_BARIA_mb_clean.RDS` | MetaPhlAn phyloseq restricted to the analysis cohort, with FFMI metadata joined |

Visit coding: `v0` = baseline, `v2`–`v5` = yearly follow-ups, `v6` = 5 years.

---

## Scripts

### `0a_datacleaning.R`

Cleans the raw clinical export and the MetaPhlAn phyloseq object, derives the
analysis variables (see **Formulas** and **Derived classifications** below),
defines the analysis cohort (see below), joins FFMI grouping metadata into the
phyloseq `sample_data`, and writes the processed long / wide / microbiome
datasets read by every other script.

### `1a_tableone_v0.R`

Builds Table 1 (baseline characteristics) for participants with available
shotgun data, stratified by baseline FFMI status (`low_ffmi_v0`), using the
`tableone` package; writes `tables/t1_v0_mb_matrix.csv`.

### `2a`–`2c` — baseline microbiome descriptives

Composition, alpha and beta diversity use the **full shotgun table without
prefiltering**; per-sample relative abundances sum to ~100%. All three read the
cleaned phyloseq object from `0a` (which already carries the FFMI grouping in
its metadata) and share a custom `theme_minimal_custom` / `theme_minimal_composition`.

- **`2a_mb_alphadiversity_v0.R`** — Shannon, Simpson and observed richness by
  baseline FFMI group (box/violin panels).
- **`2b_mb_betadiversity_v0.R`** — Bray–Curtis distances with PCoA.
- **`2c_mb_composition_v0.R`** — top-20-species composition bar plots, stratified
  by baseline FFMI group.

### `4b_lmm_mb_ffmi_trajectories.R` — longitudinal species → FFMI trajectory

- **Species set.** Poorly annotated taxa (uncharacterised `GGB` genome bins) are
  dropped first. A species is then kept if detected in **≥ 50%** of baseline
  samples **and** above a mean relative-abundance floor (`mean_abundance_v0 >= 0.1`).
  Prevalence and abundance thresholds are set blind to the FFMI results
  (independent filtering) and locked before any inferential step. Abundances
  enter on a `log10` scale with a species-specific half-minimum pseudocount,
  `log10(abundance + min_abundance/2)`.
- **Model.** Baseline abundance against FFMI trajectory, fitted per species with
  `lmerTest` as
  `ffmi ~ log10_baseline_abundance * visit + age_centered_v0 + sex + perc_change_weight_kg + (1 | id)`.
  The term of interest is the `log10_baseline_abundance × visit` interaction —
  whether baseline abundance predicts the *change* in FFMI. Run separately for
  each post-surgical follow-up (`v4`, `v5`), with p-values FDR-adjusted
  (Benjamini–Hochberg, ungrouped) within each interval.
- **Plots.** For a significant species, observed-trajectory and %-FFMI-change
  plots comparing high vs low baseline-abundance tertile groups
  (`graphs/lmm_species_ffmi/`).

Covariate choice follows DAG-style reasoning (blocking back-door paths from
microbe → FFMI without conditioning on mediators or colliders; treating insulin
resistance as a co-outcome; avoiding BMI as an over-adjustment since FFMI already
height-normalises). No standalone DAG script is included.

### `5a_humann_pathways.R` — functional pathways (HUMAnN)

Community-level functional-pathway profiles from **HUMAnN 4.0** (MetaCyc pathways),
the encoded-function layer sitting on top of the taxonomic analyses.

- **Cleaning / filtering.** Profiles are reshaped to long form; sample names are
  harmonised to the MetaPhlAn convention and restricted to the cleaned cohort.
  `UNMAPPED` and `UNINTEGRATED` are removed, and each pathway string is split into
  identifier and name. Undetected pathways are set to zero. A pathway is retained
  if it exceeds **5 CPM** in **≥ 50%** of baseline samples; abundances are
  `log10`-transformed with a `+1` pseudocount.
- **Cross-sectional models.** Baseline pathway abundance against baseline FFMI,
  across three nested adjustment sets:
  - `model1`: `ffmi_v0 ~ log10_pathway_abundance + age_v0 + sex`
  - `model2`: `+ fmi_v0`
  - `model3`: `+ t2d_v0 + dm_meds_v0 + statins_v0`

  with BH FDR adjustment across pathways.
- **Longitudinal model.** A pathway × visit LMM mirroring `4b`:
  `ffmi ~ log10_baseline_pathway_abundance * visit + age_centered_v0 + sex + perc_change_weight_kg + (1 | id)`,
  run per follow-up (`v4`, `v5`), BH-FDR-adjusted on the interaction term.

---

## Analysis cohort definition

A participant enters the final analysis cohort if **all three** hold at baseline
(`v0`):

1. **Valid baseline BIA** — non-missing FFM (kg) and FM (kg), the two BIA
   consistency checks passed, and a valid resistance reading (see
   *Body-composition quality control*).
2. **No antibiotics at baseline** (`abx_v0 == "no"`), to avoid antibiotic
   confounding of the microbiome.
3. **Available baseline shotgun sample** — a run-1 MetaPhlAn sample at
   `Time_Point == "V-1"`. For participants sequenced twice, the single run or the
   first of the duplicated runs is kept.

The cohort is the intersection of these three ID sets.

---

## Body-composition quality control

Body composition is measured by multi-frequency **bioelectrical impedance
analysis (BIA)**. Raw resistance, reactance, impedance and phase angle are
recorded at 5, 50, 100 and 200 kHz; the **50 kHz resistance** is the value
carried into the fat-free-mass / skeletal-muscle-mass equations, as is standard
for whole-body prediction.

### Plausibility filtering

**Internal-consistency and device-range checks (in `0a`)**, applied per baseline
record:

- **Completeness** — both FFM (kg) and FM (kg) must be present.
- **Percentage closure** — `|FFM% + FM% − 100| ≤ 5` percentage points.
- **Mass closure** — `|FFM_kg + FM_kg − weight_kg| ≤ 5` kg.
- **Device resistance range** — the 50 kHz resistance is accepted only within
  **110–1000 Ω**. Values outside this window, or missing, are treated as invalid
  and the muscle-mass estimate is not computed.

Records failing any check at baseline are excluded from the analysis cohort. A
small number of individually reviewed records are additionally hard-excluded
(their BIA marked invalid) in `0a`, flagged in-line for supervisor review.

### Permitted manual corrections

The correcting principle separates **data-quality errors** (impossible or
internally inconsistent values from data entry) from **true biological outliers**
(plausible but extreme values). Only the former are altered, and only when the
value is physiologically implausible, the error mechanism is identifiable, and
the correct value can be recovered (from internal consistency or the source
record):

- **Unit / decimal-scale fixes** — a field off by a consistent factor of 10 or
  100 is rescaled so it agrees with the participant's other measures.
- **Source-verified replacements** — a small number of clearly erroneous entries
  are replaced with the value read back from the source BIA record.

Plausible-but-extreme values are retained, not edited. All manual corrections are
hard-coded per participant-visit in `0a` and flagged in-line. No participant
identifiers are reproduced in this document.

---

## Formulas

All formulas applied in `0a_datacleaning.R`. Superscripts refer to **References**.

### Body composition

- **Fat-free mass index** — `FFMI = FFM_kg / height_m²` <sup>[1]</sup>
- **Fat mass index** — `FMI = FM_kg / height_m²` <sup>[2]</sup>
- **Skeletal muscle mass (Janssen BIA equation)** <sup>[3]</sup>

  ```
  SMM_kg = (height_cm² / R50 × 0.401)
         + (age × −0.071)
         + 5.102
         + (sex: male = 3.825, female = 0)
  ```

  where `R50` is the 50 kHz BIA resistance (Ω). Computed only where the resistance
  reading is valid.
- **Weight-normalised muscle mass** — `SMM_by_weight = SMM_kg / weight_kg` <sup>[4]</sup>

### Glycaemia and insulin

- **HbA1c unit conversion** (rounded IFCC–NGSP master equation) <sup>[6]</sup>

  ```
  HbA1c_%        = 0.0915 × HbA1c_mmol/mol + 2.15
  HbA1c_mmol/mol = 10.93 × HbA1c_% − 23.5
  ```

  In cleaning, a value below 15 is assumed already in %, otherwise it is converted
  from mmol/mol; missing mmol/mol values are back-calculated from %.
- **Insulin unit conversion** — `insulin_µU/mL = insulin_pmol/L / 6.945`
- **HOMA-IR** — `insulin_µU/mL × glucose_mmol/L / 22.5` <sup>[5]</sup>
- **HOMA-B (%)** — `20 × insulin_µU/mL / (glucose_mmol/L − 3.5)` <sup>[5]</sup>

### Time and trajectory

- **Follow-up time** — `n_years_from_v0 = (date − date_baseline) / 365.25`
- **Age at visit** — `age = age_v0 + n_years_from_v0`
- **Percentage change** (weight, FFM, FM) — `(X_vt − X_v0) / X_v0 × 100`
- **FFMI change** — `delta_ffmi = FFMI_vt − FFMI_v0`
- **FFMI percentage change** — `(FFMI_vt − FFMI_v0) / FFMI_v0 × 100`

---

## Derived classifications, definitions & cut-offs

### Missing-value and placeholder codes

- Numeric sentinels **-99, -98, -97** are recoded to `NA`.
- Placeholder dates **`01-01-2999`, `01-01-2997`, `01-01-2995`** are recoded to `NA`.

### Categorical recoding

| Variable | Coding |
|----------|--------|
| `sex` | 1 → male, 2 → female |
| `t2d_v0` | 1 → yes, 2 → no |
| `aht` (hypertension) | 1 → yes, 2 → no |
| `medication_binary_v0` | 1 → yes, 2 → no |
| `sport_v0` | 1 → yes, 2 → no |
| `sg_type` (surgery) | 1 → RYGB, 2 → omega-loop, 3 → sleeve gastrectomy |

### Glycaemic status (lab-based)

Applied using harmonised HbA1c (%) and fasting glucose. Thresholds follow the ADA
Standards of Care <sup>[7]</sup>.

- **T2D (`t2d_labs`)** — HbA1c **≥ 6.5%** **or** fasting glucose **≥ 7.0 mmol/L**.
- **Prediabetes (`prediab_labs`)**, among non-T2D — HbA1c **5.7–6.4%** **or**
  fasting glucose **5.6–6.9 mmol/L**.
- **Incident prediabetes / T2D** at follow-up (`new_*`) is defined only among
  participants normoglycaemic at baseline on reported T2D, lab-based T2D and
  lab-based prediabetes.

### Low-muscle-mass cut-offs

Muscle-mass status is defined by the **lowest sex- and visit-specific tertile**
(the ⅓ quantile within each `sex × visit` stratum):

- `low_ffmi` — `FFMI ≤ sex/visit ⅓ quantile`
- `low_smm` — `SMM_kg ≤ sex/visit ⅓ quantile`
- `low_smm_by_weight` — `SMM_by_weight ≤ sex/visit ⅓ quantile`

### FFMI trajectory groups

Percentage change in FFMI from baseline (at `v4` and `v5`) is grouped within
**sex**: the lowest tertile (most negative change) is labelled **"high"** (high
FFMI loss); the remainder **"low/modest"**.

### Medication categories

Free-text baseline medication is lower-cased, typo-corrected, and matched with
word-boundary regular expressions to derive binary (`yes`/`no`) indicators per
drug and per class. Term lists (Dutch and English generics and brand names) are
defined inline in `0a`; the categories are:

| Category | Sub-classes captured |
|----------|----------------------|
| **Diabetes** | metformin, sulfonylureas, DPP-4i, GLP-1RA, SGLT2i, TZD, insulin |
| **Antihypertensive** | ACE-i, ARB, CCB, β-blockers, central α2-agonists, diuretics, combinations |
| **Lipid-lowering** | statins, ezetimibe, PCSK9i, inclisiran, bempedoic acid, bile-acid sequestrants, fibrates, omega-3, niacin |
| **Thyroid** | levothyroxine |
| **Psychiatric** (weight-relevant) | SSRI, TCA, SNRI, NDRI, atypical antipsychotics, mood stabilisers, ADHD stimulants, hypnotics |
| **PPI** | omeprazole, pantoprazole, esomeprazole |
| **Antibiotics** | penicillins, tetracyclines, macrolides, fluoroquinolones, lincosamides, nitroimidazoles, urinary-tract agents, sulfonamides |

Two carry into the modelling: **antibiotics** are a cohort exclusion (above), and
**statins** enter the `5a` cross-sectional extended adjustment set (`model3`).

---

## Dependencies

R, with (across scripts): `tidyverse`, `phyloseq`, `vegan`, `tableone`,
`lmerTest`, `broom.mixed`, `ggpubr`, `ggrepel`, `ggthemes`, `ggsci`, `patchwork`,
`grid`, and the palette package `MetBrewer`.

---

## References

1. Kawakami R, Tanisawa K, Ito T, et al. Fat-Free Mass Index as a surrogate marker
   of appendicular skeletal muscle mass index for low muscle mass screening in
   sarcopenia. *J Am Med Dir Assoc.* 2022;23(12):1955–1961.e3.
   doi:10.1016/j.jamda.2022.08.016.
2. VanItallie TB, Yang MU, Heymsfield SB, Funk RC, Boileau RA. Height-normalized
   indices of the body's fat-free mass and fat mass. *Am J Clin Nutr.*
   1990;52(6):953–959. doi:10.1093/ajcn/52.6.953.
3. Janssen I, Heymsfield SB, Baumgartner RN, Ross R. Estimation of skeletal muscle
   mass by bioelectrical impedance analysis. *J Appl Physiol.* 2000;89(2):465–471.
   doi:10.1152/jappl.2000.89.2.465.
4. Donini LM, Busetto L, Bischoff SC, et al. Definition and diagnostic criteria for
   sarcopenic obesity: ESPEN and EASO consensus statement. *Obes Facts.*
   2022;15(3):321–335. doi:10.1159/000521241.
5. Matthews DR, Hosker JP, Rudenski AS, et al. Homeostasis model assessment:
   insulin resistance and beta-cell function from fasting plasma glucose and
   insulin concentrations in man. *Diabetologia.* 1985;28(7):412–419.
   doi:10.1007/BF00280883.
6. Hoelzel W, Weykamp C, Jeppsson JO, et al. IFCC Reference System for measurement
   of hemoglobin A1c in human blood. *Clin Chem.* 2004;50(1):166–174.
   doi:10.1373/clinchem.2003.024802.
7. American Diabetes Association. Standards of Care in Diabetes — 2026.
   *Diabetes Care.* 2026;49(Suppl 1).
