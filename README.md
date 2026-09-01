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

The code currently covers data cleaning and cohort definition, split into a
clinical layer (`0a`) and an omics layer — microbiome, HUMAnN pathways,
metabolomics (`0b`); Table 1 (`1a`); the baseline microbiome descriptive
analyses — composition and diversity (`2a`–`2c`); a longitudinal mixed-model
layer relating baseline species abundance to FFMI trajectory after surgery
(`4b`); a functional-pathway layer that filters and models the HUMAnN profiles
against FFMI (`5a`); a gradient-boosting (XGBeast/XGBoost) pipeline that ranks
species and pathways by how well they predict FFMI and its post-surgical
change, with LM/LMM forest plots confirming the top features (`6a`–`6d`); and
correlation analyses relating the top predictive species to the serum
metabolome and to baseline dietary macronutrient intake (`7a`–`7b`).

---

## Repository structure

```
baria_muscle/
├── scripts/
│   ├── 0a_datacleaning_clinical.R          # Clinical cleaning, derived variables, cohort definition
│   ├── 0b_datacleaning_omics.R             # Microbiome, HUMAnN pathway & metabolomics cleaning
│   ├── 1a_tableone_v0.R                    # Table 1, stratified by baseline FFMI status
│   ├── 2a_mb_alphadiversity_v0.R           # Alpha diversity (Shannon, Simpson, observed richness)
│   ├── 2b_mb_betadiversity_v0.R            # Beta diversity (Bray–Curtis PCoA)
│   ├── 2c_mb_composition_v0.R              # Compositional bar plots (top 20 species, by baseline FFMI group)
│   ├── 4b_lmm_mb_ffmi_trajectories.R       # Mixed models: baseline species abundance × FFMI trajectory
│   ├── 5a_humann_pathways.R                # HUMAnN pathways: filtering + FFMI models
│   ├── 6a_ml_prep.R                        # Build XGBeast input tables per outcome/subgroup
│   ├── 6b_ml_run.sh                        # SLURM job: fit XGBeast models for every outcome/subgroup
│   ├── 6c_ml_process.R                     # Feature importance, LM/LMM confirmation, forest plots
│   ├── 6d_ml_explained_variance_violin.R   # Violin plots of explained variance across CV iterations
│   ├── 7a_mb_metabolome_correlations.R     # Top-15 predictive species × serum metabolome correlations
│   ├── 7b_mb_diet_correlations.R           # Top-15 predictive species × dietary macronutrient correlations
│   └── assets/
│       ├── functions.R                     # Shared helper functions (plotting, XGBeast I/O, forest models)
│       ├── XGBeast_new.py                  # Repeated-CV XGBoost + feature-importance framework
│       ├── param_grid.json                 # XGBoost hyperparameter grid used by 6b
│       └── grid_tuning.json                # Grid-search tuning config
├── pixi.toml                               # Environment manifest
├── pixi.lock                               # Locked environment
├── .gitignore
└── README.md
```

Scripts are numbered and run in that order. Gaps in the numbering (there is no
`3*` and no `4a`) correspond to superseded scripts kept locally in a
**git-ignored** `scripts/archive/` and are not tracked in this repository.

The `data/` and `results/` directories, and the pixi environment
(`.pixi/`), are also git-ignored: they are created and populated locally when the
scripts run (see **Data** below) and are not part of the repository.

### Order of execution

`0a` runs first (clinical cleaning) and writes the processed clinical datasets;
`0b` runs next and reads `0a`'s output alongside the raw microbiome, HUMAnN and
metabolomics files to write the processed omics datasets. `1a` and the `2*`
descriptive scripts are independent of one another and read the processed data
written by `0a`/`0b`. `4b` and `5a` fit their models independently of each
other and of the `2*` scripts.

`6a` (build ML input) → `6b` (fit XGBeast models, typically on a SLURM cluster)
→ `6c` (feature importance + confirmatory LM/LMM forest plots) → `6d`
(explained-variance violin plots) must run in that order, since each step
reads the previous step's output under `results/mlmodels/`. `7a` and `7b` are
independent of one another but both depend on `6b`'s output for the
`delta_ffmi_v4`/`all` model, from which they take the top 15 predictive
species.

---

## Data

**Inputs** (placed in `data/raw_data/`, not tracked in the repository):

| File | Description |
|------|-------------|
| `BARIA.clinical.2024-12-09.723.2043.RDS` | Clinical, anthropometric and BIA data |
| `ps.BARIA.metaphlan.706.2548.RDS` | MetaPhlAn shotgun phyloseq object (per-sample relative abundances, 0–100%) |
| `BARIA.humann4.profiles.2026.581.910.RDS` | HUMAnN 4.0 community functional-pathway profiles (per-sample pathway abundances, CPM) |
| `BARIA.metabolon.1158.V12.BatchNorm.RDS` | Metabolon serum metabolomics, batch-normalised phyloseq object |
| `BARIA.metabolon.1158.V12.Peak.Area.RDS` | Metabolon serum metabolomics, raw peak-area phyloseq object (used only to flag imputed/undetected values) |
| `251002_BARIA_macronutrients.RDS` | Baseline dietary macronutrient intake, per diary tool (`V1`) |

**Outputs written by `0a`** to `data/processed_data/` (clinical):

| File | Description |
|------|-------------|
| `BARIA_muscle_long.{RDS,csv}` | Long clinical dataset (one row per participant × visit) |
| `BARIA_muscle_wide.{RDS,csv}` | Wide clinical dataset (one row per participant) |

**Outputs written by `0b`** to `data/processed_data/` (omics; `0b` reads `0a`'s
`BARIA_muscle_wide.RDS` to restrict samples to the analysis cohort):

| File | Description |
|------|-------------|
| `BARIA_mb_clean.RDS` | MetaPhlAn phyloseq restricted to the analysis cohort, with FFMI metadata joined; Eukaryota and unannotated `GGB` genera dropped |
| `BARIA_mb_clean_unfiltered.RDS` | Same as above but before dropping Eukaryota/`GGB` genera — used for diversity |
| `BARIA_mb_baseline.RDS` | `BARIA_mb_clean.RDS` restricted to `v0` samples |
| `BARIA_humann_pathways_long.RDS` | HUMAnN pathway table reshaped to long form, restricted to faecal samples in the analysis cohort, `UNMAPPED`/`UNINTEGRATED` removed |
| `BARIA_metabolon_clean.RDS` | Baseline serum metabolomics phyloseq object: unidentified/sparse (>10% imputed) compounds dropped, log10-transformed and z-scored |

Visit coding: `v0` = baseline (raw `V-1` / MetaPhlAn `Time_Point == "V-1"`),
`v2`–`v5` = yearly follow-ups, `v6` = 5 years, `v7` = 10 years.

---

## Analysis cohort definition

A participant enters the final analysis cohort if **all three** hold at baseline
(`v0`):

1. **Valid baseline BIA** — non-missing FFM (kg) and FM (kg), the two BIA
   consistency checks passed, and a valid resistance reading (see
   *Body-composition quality control*).
2. **No antibiotics at baseline** (`abx == "no"`), to avoid antibiotic
   confounding of the microbiome.
3. **Available baseline shotgun sample** — a run-1 MetaPhlAn sample at
   `Time_Point == "V-1"`. For participants sequenced twice, the single run or the
   first of the duplicated runs is kept.

The cohort is the intersection of these three ID sets.

---

## Body-composition quality control

Body composition is measured by multi-frequency **bioelectrical impedance
analysis (BIA)** (Maltron BioScan 920). The **50 kHz resistance** is the value
carried into the fat-free-mass / skeletal-muscle-mass equations, as is standard
for whole-body prediction.

### Internal-consistency and device-range checks (in `0a`)

Applied per baseline record:

- **Completeness** — both FFM (kg) and FM (kg) must be present.
- **Percentage closure** — `|FFM% + FM% − 100| ≤ 5` percentage points. The two
  compartments are reported as complementary percentages and should sum to ~100%;
  a larger discrepancy indicates an inconsistent or mis-entered reading.
- **Mass closure** — `|FFM_kg + FM_kg − weight_kg| ≤ 5` kg. The two compartment
  masses should reconstruct the measured body weight to within a small tolerance.
- **Device resistance range** — the 50 kHz resistance is accepted only within
  **110–1000 Ω** (the plausible operating range for the device). Values outside
  this window, or missing, are treated as invalid and the muscle-mass estimate is
  not computed.

Records failing any check at baseline are excluded from the analysis cohort
(together with the antibiotic and shotgun-availability criteria above). A small
number of individually reviewed records are additionally hard-excluded (their BIA
marked invalid) in `0a`, flagged in-line for supervisor review.

### Permitted manual corrections

The correcting principle is to separate **data-quality errors** (impossible or
internally inconsistent values arising from data entry) from **true biological
outliers** (plausible but extreme values). Only the former are altered.

A value is corrected only when all of the following hold: it is physiologically
implausible, the error mechanism is identifiable, and the correct value can be
recovered — either from internal consistency with the participant's own
body-composition fields (percentage and mass closure above) or by checking the
source record. Two correction types are used:

- **Unit / decimal-scale fixes** — an isolated field off by a consistent factor of
  10 or 100 (a misplaced decimal or wrong unit scale) is rescaled by that factor
  so that it agrees with the participant's other measures and restores percentage
  and mass closure.
- **Source-verified replacements** — a small number of clearly erroneous entries
  are replaced with the correct value read back from the source BIA record.

Values that are plausible but extreme (true statistical outliers) are **not**
edited; they are retained and handled analytically rather than corrected. All
manual corrections are hard-coded per participant-visit in `0a` and flagged
in-line for supervisor review. No participant identifiers are reproduced in this
document.

---

## Formulas

All formulas applied in `0a_datacleaning_clinical.R`. Superscripts refer to the
**References** section.

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
- **Weight-normalised muscle mass** — `SMM_by_weight = SMM_kg / weight_kg`
  (the ESPEN/EASO-favoured BIA index for low muscle mass in obesity) <sup>[4]</sup>

### Glycaemia and insulin

- **HbA1c unit conversion** (rounded IFCC–NGSP master equation) <sup>[6]</sup>

  ```
  HbA1c_%        = 0.0915 × HbA1c_mmol/mol + 2.15
  HbA1c_mmol/mol = 10.93 × HbA1c_% − 23.5
  ```

  In cleaning, a value below 15 is assumed already to be in %, otherwise it is
  converted from mmol/mol; missing mmol/mol values are back-calculated from %.
- **HOMA-IR** — `(insulin_pmol/L / 6.945) × glucose_mmol/L / 22.5` <sup>[5]</sup>
- **HOMA-B (%)** — `20 × (insulin_pmol/L / 6.945) / (glucose_mmol/L − 3.5)` <sup>[5]</sup>

  (insulin is converted from pmol/L to µU/mL, 1 µU/mL = 6.945 pmol/L)

### Time and trajectory

- **Follow-up time** — `n_years_from_v0 = (date − date_baseline) / 365.25`
- **Age at visit** — `age = age_v0 + n_years_from_v0`
- **Percentage change** (weight, FFM, FM) — `(X_vt − X_v0) / X_v0 × 100`
- **FFMI change** — `delta_ffmi = FFMI_vt − FFMI_v0`
- **FFMI percentage change** — `(FFMI_vt − FFMI_v0) / FFMI_v0 × 100`

---

## Derived classifications, definitions & cut-offs

### Missing-value and placeholder codes

- Numeric sentinel values **-99, -98, -97** are recoded to `NA`.
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

Applied using harmonised HbA1c (%) and fasting glucose (see Formulas). Thresholds
follow the ADA Standards of Care <sup>[7]</sup>.

- **T2D (`t2d_labs`)** — HbA1c **≥ 6.5%** **or** fasting glucose **≥ 7.0 mmol/L**.
- **Prediabetes (`prediab_labs`)**, among non-T2D — HbA1c **5.7–6.4%**
  **or** fasting glucose **5.6–6.9 mmol/L**.
- **Incident prediabetes / T2D** at follow-up (`new_*`) is defined only among
  participants normoglycaemic at baseline on all of: reported T2D, lab-based T2D,
  and lab-based prediabetes.

### Low-muscle-mass cut-offs

Muscle-mass status is defined by the **lowest sex- and visit-specific tertile**
(the ⅓ quantile within each `sex × visit` stratum), consistent with the common
practice in sarcopenic-obesity research of using the lowest sex-stratified tertile:

- `low_ffmi` — `FFMI ≤ sex/visit ⅓ quantile`
- `low_smm` — `SMM_kg ≤ sex/visit ⅓ quantile`
- `low_smm_by_weight` — `SMM_by_weight ≤ sex/visit ⅓ quantile`

### FFMI trajectory groups

Percentage change in FFMI from baseline (at `v4` and `v5`) is grouped within
**sex**: the lowest tertile (most negative change) is labelled **"high"**
(high FFMI loss); the remainder **"modest/low"**.

### Medication categories

Free-text baseline medication is lower-cased, typo-corrected, and matched with
word-boundary regular expressions to derive binary (`yes`/`no`) indicators per
drug and per class. The operational term lists (Dutch and English generics and
brand names) are defined inline in `0a`; the categories are:

| Category | Sub-classes captured |
|----------|----------------------|
| **Diabetes** | metformin, sulfonylureas, DPP-4i, GLP-1RA, SGLT2i, TZD, insulin |
| **Antihypertensive** | ACE-i, ARB, CCB, β-blockers, central α2-agonists, diuretics, combinations |
| **Lipid-lowering** | statins, ezetimibe, PCSK9i, inclisiran, bempedoic acid, bile-acid sequestrants, fibrates, omega-3, niacin |
| **Thyroid** | levothyroxine |
| **Psychiatric** (weight-relevant) | SSRI, TCA, SNRI, NDRI, atypical antipsychotics, mood stabilisers, ADHD stimulants, hypnotics |
| **PPI** | omeprazole, pantoprazole, esomeprazole (microbiome-relevant) |
| **Antibiotics** | penicillins, tetracyclines, macrolides, fluoroquinolones, lincosamides, nitroimidazoles, urinary-tract agents, sulfonamides |

Two of these carry into the modelling: **antibiotics** are a cohort exclusion
(above), and diabetes medication and **statins** enter the `5a` adjustment set.

---

## Microbiome analysis

### Descriptive analyses (`2a`–`2c`)

- Composition, alpha and beta diversity use the **full shotgun table without
  prefiltering**; per-sample relative abundances sum to ~100%.
- Alpha diversity: Shannon, Simpson and observed richness. Beta diversity:
  Bray–Curtis distances with PCoA. Composition: the top 20 species, stratified by
  baseline FFMI group.

### Longitudinal association model (`4b`)

- **Species set.** Poorly annotated taxa (uncharacterised `GGB` genome bins) are
  dropped first. A species is then kept if detected in **≥ 50%** of baseline
  samples **and** with **mean relative abundance ≥ 0.1** (on the 0–100% scale).
  Prevalence and abundance thresholds are set blind to the FFMI results
  (independent filtering) and locked before any inferential step. Abundances enter
  on a `log10` scale with a species-specific half-minimum pseudocount,
  `log10(abundance + min_abundance/2)`.
- **Mixed models.** Baseline abundance against FFMI trajectory, fitted per species
  with `lmerTest` (ML) as
  `ffmi ~ log10_baseline_abundance * visit + age_centred + sex + %weight_change +
  (1 | id)` (random intercept per participant). The term of interest is the
  `log10_baseline_abundance × visit` interaction — whether baseline abundance
  predicts the *change* in FFMI. The model is run separately for each post-surgical
  follow-up (`v4`, `v5`; 1-year and 2-year), with p-values FDR-adjusted
  (Benjamini–Hochberg) across species — ungrouped — within each interval.
- **Plots.** An observed-FFMI-trajectory plot across baseline / 1 year / 2 years;
  and, for any species significant at 1 year (`v4`), a plot of %FFMI change by
  baseline-abundance tertile group (low vs high). Figures are written to
  `results/graphs/lmm_species_ffmi/`.

Covariate choice adjusts for age, sex and %weight change; BMI is deliberately not
added, since FFMI already height-normalises and adjusting for it would be an
over-adjustment.

### Functional pathways (`5a`)

Community-level functional-pathway profiles from **HUMAnN 4.0** (MetaCyc pathways).
Reshaping to long form, sample-name harmonisation to the MetaPhlAn convention,
restriction to **faecal** samples in the cleaned cohort, removal of the
`UNMAPPED`/`UNINTEGRATED` categories, splitting each pathway string into its
identifier and name, and zero-filling undetected pathways all happen in `0b`
(see **Data** above); `5a` itself starts from that cleaned long table.

- A pathway is retained if it exceeds **5 CPM** in **≥ 20%** of baseline samples,
  and abundances are `log10`-transformed with a species-specific half-minimum
  pseudocount, `log10(abundance + min_abundance/2)`.
- **Baseline pathway ~ FFMI linear models**, fitted per pathway across three nested
  adjustment sets (M1: age + sex; M2: + FMI; M3: + T2D + diabetes medication +
  statins), with Benjamini–Hochberg FDR across pathways.
- **Baseline pathway abundance × FFMI trajectory LMMs** at `v4` and `v5`, mirroring
  the `4b` model structure, again with BH-FDR across pathways, plus a check of the
  nominal cross-interval (`v4` vs `v5`) overlap.

No pathways survive FDR correction in either the baseline LMs or the trajectory
LMMs.

---

## Machine-learning pipeline (`6a`–`6d`)

A gradient-boosting layer that complements the linear/mixed models above by
letting all baseline features compete jointly (rather than one species/pathway
at a time) and ranking them by predictive contribution, using **XGBeast**
(`scripts/assets/XGBeast_new.py`) — a repeated shuffle-split **XGBoost**
regression framework with grid-search hyperparameter tuning that reports
explained variance and permutation feature importance across CV iterations.

- **`6a_ml_prep.R`** builds one input table per outcome × subgroup
  (`all`/`male`/`female`) under `results/mlmodels/<outcome>/<subgroup>/`.
  Predictors are baseline (`v0`) species (prevalence/abundance-filtered as in
  `4b`) or, for the `*_path` outcomes, baseline HUMAnN pathways
  (filtered analogously). Outcomes cover baseline FFMI/FMI, baseline FFMI
  adjusted for FMI (residualised), the same two at `v4` (cross-sectional,
  predicted from **`v4` species**, not baseline), baseline FFMI/FMI restricted
  to the subjects with `v4` microbiome data (`*_matched`, to isolate subject-
  subset effects from timepoint effects), and 1-year change in FFMI
  (`delta_ffmi_v4`, `perc_change_ffmi_v4`) with and without FMI adjustment.
- **`6b_ml_run.sh`** is a SLURM batch script that calls `XGBeast_new.py` once
  per outcome × subgroup (17 outcomes × 3 subgroups) with 200 iterations, an
  8-thread grid search (`param_grid.json`) and a fixed random seed.
- **`6c_ml_process.R`** reads each model's XGBeast output and, per
  outcome/subgroup: plots the explained-variance distribution across CV
  iterations, extracts the top-15 features by importance, and **confirms them**
  with a classical model on the same predictors — an LM for cross-sectional
  outcomes or an LMM (mirroring `4b`'s baseline-abundance × visit structure)
  for change outcomes — reporting a forest plot of per-feature effect estimates
  with BH-FDR across the top-15 features.
- **`6d_ml_explained_variance_violin.R`** plots the explained-variance
  distribution across CV iterations for the six headline FFMI/delta-FFMI
  models (± FMI adjustment), for the species- and pathway-abundance model sets
  separately, to `results/graphs/ml_explained_variance/`.

---

## Correlation analyses (`7a`–`7b`)

Both scripts take the **top-15 species by feature importance** from the
`delta_ffmi_v4`/`all` XGBeast model (`6b`'s output) — i.e. the baseline species
most predictive of 1-year FFMI change — and correlate their `log10`-transformed
baseline relative abundance against an external baseline dataset, using
Spearman correlations with Benjamini–Hochberg FDR across all bug × feature
pairs.

- **`7a_mb_metabolome_correlations.R`** — correlates the top-15 species against
  every baseline serum metabolite in `BARIA_metabolon_clean.RDS`. Writes the
  full correlation table, a scatter plot per species (one panel per
  FDR-significant metabolite), and `ComplexHeatmap` bugs × metabolites heatmaps
  (all FDR-significant metabolites, and a top-25-by-hit-count compact version),
  annotated by Metabolon super-pathway, to `results/graphs/mb_metabolome_correlations/`.
- **`7b_mb_diet_correlations.R`** — correlates the top-15 species against
  baseline energy-normalised macronutrient intake (carbohydrate, protein, fat,
  saturated fat, fibre, alcohol, plus total energy) from
  `251002_BARIA_macronutrients.RDS`, run separately per dietary-diary tool.
  Species are additionally required to be present in **≥ 30%** of participants
  within each tool's subset before testing. Because little survives FDR
  correction at this sample size, significance here is reported **nominal**
  (uncorrected `p < 0.05`), not FDR-corrected. Outputs (correlation table,
  scatter plots, heatmap) are written per tool to
  `results/graphs/mb_diet_correlations/<tool>/`.

---

## Dependencies

The environment is managed with **pixi** (`pixi.toml`, `pixi.lock`), and defines
pixi tasks for every script plus grouped tasks (`datacleaning`, `mb`, `ml`,
`mb-correlations`) that chain the steps within each stage — run e.g.
`pixi run ml` or `pixi run datacleaning-clinical`; see `pixi.toml` for the full
list. R packages loaded across the scripts include: `tidyverse`, `phyloseq`,
`vegan`, `tableone`, `ggpubr`, `patchwork`, `ggthemes`, `ggsci`, `ggrepel`,
`grid`, `MetBrewer`, `lmerTest`, `broom` / `broom.mixed`, and (`7a`/`7b`)
`ComplexHeatmap` / `circlize`. The ML pipeline (`6b`) additionally uses Python:
`xgboost`, `scikit-learn`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `shap`,
`tqdm`.

> **Note:** `ComplexHeatmap` and `circlize` (used by `7a`/`7b`) are not
> currently pinned in `pixi.toml` — install them separately (e.g. via
> Bioconductor) until the manifest is updated.

Expected working-directory layout (script paths are relative to the repo root):

```
data/raw_data/          # inputs (see Data)
data/processed_data/    # created by 0a (clinical) and 0b (omics)
results/graphs/         # figures, in per-analysis subfolders
results/tables/         # Table 1 output
results/mlmodels/       # ML input data + XGBeast output, per outcome/subgroup
```

Each script creates its own output folders (via `dir.create(..., recursive = TRUE)`)
before writing to them, so `results/` does not need to exist beforehand.

---

## References

1. Kawakami R, Tanisawa K, Ito T, et al. Fat-Free Mass Index as a surrogate marker
   of appendicular skeletal muscle mass index for low muscle mass screening in
   sarcopenia. *J Am Med Dir Assoc.* 2022;23(12):1955–1961.e3.
   doi:10.1016/j.jamda.2022.08.016.
2. VanItallie TB, Yang MU, Heymsfield SB, Funk RC, Boileau RA. Height-normalized
   indices of the body's fat-free mass and fat mass: potentially useful indicators
   of nutritional status. *Am J Clin Nutr.* 1990;52(6):953–959.
   doi:10.1093/ajcn/52.6.953.
3. Janssen I, Heymsfield SB, Baumgartner RN, Ross R. Estimation of skeletal muscle
   mass by bioelectrical impedance analysis. *J Appl Physiol.* 2000;89(2):465–471.
   doi:10.1152/jappl.2000.89.2.465.
4. Donini LM, Busetto L, Bischoff SC, et al. Definition and diagnostic criteria for
   sarcopenic obesity: ESPEN and EASO consensus statement. *Obes Facts.*
   2022;15(3):321–335. doi:10.1159/000521241.
5. Matthews DR, Hosker JP, Rudenski AS, Naylor BA, Treacher DF, Turner RC.
   Homeostasis model assessment: insulin resistance and beta-cell function from
   fasting plasma glucose and insulin concentrations in man. *Diabetologia.*
   1985;28(7):412–419. doi:10.1007/BF00280883.
6. Hoelzel W, Weykamp C, Jeppsson JO, et al. IFCC Reference System for measurement
   of hemoglobin A1c in human blood. *Clin Chem.* 2004;50(1):166–174.
   doi:10.1373/clinchem.2003.024802.
7. American Diabetes Association. Standards of Care in Diabetes — 2026.
   *Diabetes Care.* 2026;49(Suppl 1). (Diagnostic thresholds for T2D and prediabetes.)
