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
- **Other omics:** bulk RNA-seq (liver, jejunum, visceral/subcutaneous fat); serum metabolomics; dietary macronutrient intake
- **Language:** R

The code covers data cleaning and cohort definition, split into a clinical
layer (`0a`), a microbiome/HUMAnN/metabolomics layer (`0b`), and an RNA-seq
layer (`0c`); Table 1, stratified by 1-year %FFMI change (`1a`); baseline
microbiome descriptive analyses — alpha diversity, beta diversity and
composition, stratified the same way (`2a`–`2c`); a gradient-boosting
(XGBeast/XGBoost) pipeline that ranks baseline species by how well they
predict FFMI and its post-surgical change, with each model's top-15 features
confirmed by LM/LMM forest plots (`3a`–`3c`); correlation analyses relating
the species most predictive of 1-year FFMI change to HUMAnN pathways, the
serum metabolome, and baseline dietary macronutrient intake (`4a`–`4c`); and
liver and visceral adipose tissue RNA-seq (`5a`–`5b`), still under active
development.

---

## Repository structure

```
baria_muscle/
├── scripts/
│   ├── 0a_datacleaning_clinical.R          # Clinical cleaning, derived variables, cohort definition
│   ├── 0b_datacleaning_omics.R             # Microbiome, HUMAnN pathway & metabolomics cleaning
│   ├── 0c_datacleaning_rnaseq.R            # Bulk RNA-seq cleaning, split by tissue
│   ├── 1a_tableone_ffmi_1y.R               # Table 1, stratified by 1-year %FFMI change group
│   ├── 2a_mb_alphadiversity_ffmi_1y.R      # Alpha diversity (Shannon, Simpson, observed richness)
│   ├── 2b_mb_betadiversity_ffmi_1y.R       # Beta diversity (Bray–Curtis PCoA)
│   ├── 2c_mb_composition_ffmi_1y.R         # Compositional bar plots (top 20 species)
│   ├── 3a_ml_prep.R                        # Build XGBeast input tables per outcome/subgroup
│   ├── 3b_ml_run.sh                        # SLURM job: fit XGBeast models for every outcome/subgroup
│   ├── 3c_ml_process.R                     # Feature importance, LM/LMM confirmation + forest plots, explained-variance violin plots
│   ├── 4a_humann_pathways.R                # Top-15 species x HUMAnN pathway correlations
│   ├── 4b_mb_metabolome_correlations.R     # Top-15 species x serum metabolome correlations
│   ├── 4c_mb_diet_correlations.R           # Top-15 species x dietary macronutrient correlations
│   ├── 5a_liver_rnaseq.R                   # Liver RNA-seq vs. top-15 species / FFMI change — in progress
│   ├── 5b_vadipose_rnaseq.R                # Visceral adipose RNA-seq vs. top-15 species / FFMI change — in progress
│   └── assets/
│       ├── functions.R                     # Shared helper functions (plotting, XGBeast I/O, forest models)
│       ├── XGBeast_new.py                  # Repeated-CV XGBoost + feature-importance framework
│       └── param_grid.json                 # XGBoost hyperparameter grid used by 3b
├── pixi.toml                               # Environment manifest
├── pixi.lock                               # Locked environment
├── .gitignore
└── README.md
```

Scripts are numbered and run in that order. Earlier versions of several
scripts — baseline-FFMI-status Table 1/diversity/composition scripts, and a
standalone species/pathway-vs-FFMI mixed-model screen that predates the
ML-first pipeline — are kept locally in a **git-ignored** `scripts/archive/`
and are not tracked in this repository.

The `data/` and `results/` directories, and the pixi environment
(`.pixi/`), are also git-ignored: they are created and populated locally when the
scripts run (see **Data** below) and are not part of the repository.

### Order of execution

`0a` runs first (clinical cleaning) and writes the processed clinical
datasets; `0b` runs next and reads `0a`'s output alongside the raw
microbiome, HUMAnN and metabolomics files to write the processed omics
datasets; `0c` also reads `0a`'s output (to restrict to the analysis cohort)
alongside the raw RNA-seq counts, runs a per-tissue sample QC (library size and
PCA) and writes one processed dataset per tissue (see **Tissue RNA-seq**).
`1a` and the `2*` descriptive scripts are independent of one another and read
the processed data written by `0a`/`0b`.

`3a` (build ML input) → `3b` (fit XGBeast models, typically on a SLURM
cluster) → `3c` (feature importance, LM/LMM confirmation of each model's
top-15 features, and explained-variance violin plots) must run in that
order, since each step reads the previous step's output under
`results/mlmodels/`. `4a`, `4b` and `4c` are independent of one another but
all three depend on `3b`'s output for the `perc_change_ffmi_v4`/`all` model,
from which they take the top-15 species most predictive of 1-year FFMI
change. `5a` and `5b` additionally depend on `0c` (liver and visceral adipose RNA-seq,
respectively) and on `3c`'s forest-plot output for that same model, to colour
species by the direction of their association with FFMI change.

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
| `RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS` | Bulk RNA-seq (kallisto) count table, one column per sample across liver, subcutaneous fat, visceral fat and jejunum |

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

**Outputs written by `0c`** to `data/processed_data/` (RNA-seq; `0c` reads
`0a`'s `BARIA_muscle_wide.RDS` to restrict samples to the analysis cohort,
keeps baseline (`V1`) samples only, and drops duplicated `REFER` runs and
non-numeric/comma-decimal artefacts in the raw counts):

| File | Description |
|------|-------------|
| `BARIA_Liver_RNAseq.RDS` | Liver RNA-seq raw counts, **all** samples (one row per participant, `ENSG*` gene columns) plus a `qc_exclude` flag (TRUE = excluded at sample QC) |
| `BARIA_Jejunum_RNAseq.RDS` | Jejunum RNA-seq raw counts, same structure |
| `BARIA_vFat_RNAseq.RDS` | Visceral fat RNA-seq raw counts, same structure |
| `BARIA_subFat_RNAseq.RDS` | Subcutaneous fat RNA-seq raw counts, same structure |
| `BARIA_Liver_RNAseq_vst.RDS` | Liver: DESeq2 size factor + VST-normalised expression for QC-retained samples only; size factors, expression filter (≥10 normalised reads in ≥50% of participants) and VST computed on the retained samples — used by `5a` |
| `BARIA_Jejunum_RNAseq_vst.RDS` | Jejunum: same processing |
| `BARIA_vFat_RNAseq_vst.RDS` | Visceral fat: same processing — used by `5b` |
| `BARIA_subFat_RNAseq_vst.RDS` | Subcutaneous fat: same processing |

`0c` also writes its sample-QC outputs: `results/graphs/RNAseq/qc_librarysize_<tissue>.pdf`,
`results/graphs/RNAseq/qc_pca_<tissue>.pdf` and `results/tables/rnaseq_qc_samples.csv`
(see **Tissue RNA-seq**).

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

### Insulin sensitivity, hepatic insulin resistance and insulin clearance

Fasting values are the `0`-min mixed-meal test (MMT) samples. Unit conversions:
glucose mmol/L × 18.016 → mg/dL; insulin pmol/L / 6.945 → µU/mL; C-peptide
nmol/L × 1000 → pmol/L; triglycerides mmol/L × 88.57 → mg/dL. MMT-derived indices
were originally validated on the OGTT and are applied here to the MMT; each is
computed only where all required timepoints are available.

- **Triglyceride–glucose index (TyG)** — `ln(TG_mg/dL × glucose_mg/dL / 2)` <sup>[14]</sup>
- **Matsuda index** (whole-body insulin sensitivity; 0/30/60/90/120 min) <sup>[16]</sup>

  ```
  Matsuda = 10000 / √(G0 × I0 × mean(G0–120) × mean(I0–120))      [G mg/dL, I µU/mL]
  ```
- **Hepatic insulin resistance index** (0/10/20/30 min, trapezoidal AUC) <sup>[15]</sup>

  ```
  hepatic_ir_index = AUC0–30 glucose [mg/dL·min] × AUC0–30 insulin [µU/mL·min] / 10^6
  ```
- **Hepatic insulin clearance, fasting** — `C-peptide_pmol/L / insulin_pmol/L` (molar ratio) <sup>[17, 18]</sup>
- **Hepatic insulin clearance, MMT** — ratio of trapezoidal AUC0–120 C-peptide to
  AUC0–120 insulin (0/10/20/30/60/90/120 min) <sup>[18]</sup>

### Liver non-invasive tests (NITs)

- **FIB-4** — `age × AST / (platelets × √ALT)` <sup>[8]</sup>
- **NAFLD fibrosis score (NFS)** <sup>[11]</sup>

  ```
  NFS = −1.675 + 0.037 × age + 0.094 × BMI + 1.13 × (IFG or diabetes)
        + 0.99 × (AST/ALT) − 0.013 × platelets − 0.66 × albumin_g/dL
  ```

  IFG = fasting glucose ≥ 6.1 mmol/L; diabetes = `t2d_any`; albumin g/L / 10 → g/dL.
- **Fatty liver index (FLI)** <sup>[12]</sup>

  ```
  FLI = 100 × e^y / (1 + e^y),
  y   = 0.953 × ln(TG_mg/dL) + 0.139 × BMI + 0.718 × ln(GGT) + 0.053 × waist_cm − 15.745
  ```
- **Hepatic steatosis index (HSI)** — `8 × (ALT/AST) + BMI (+2 if female; +2 if diabetes)` <sup>[13]</sup>

Units: AST, ALT, GGT in U/L; platelets ×10⁹/L; age in years at the visit.

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
- Platelet counts at `v4`/`v5` are stored under two complementary (non-overlapping)
  field spellings in the raw export (`trombocytes` / `thrombocytes`); sentinel codes are
  set to `NA` first, then the two are merged with `coalesce()`.

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
- **T2D, any (`t2d_any`)** — used in the liver NITs. At baseline: reported T2D
  (`t2d_v0`) **or** lab-based T2D. At follow-up: lab-based T2D **or** use of
  glucose-lowering medication (the baseline report is not carried forward, so that
  post-surgical remission is reflected).

### Liver NIT categories

- **FIB-4** (`fib4_cat`) — < 1.30 low, 1.30–2.67 indeterminate, > 2.67 high risk of
  advanced fibrosis <sup>[9]</sup>. The age-adjusted low cut-off of 2.0 for ages ≥ 65
  <sup>[10]</sup> is not applied.
- **NFS** (`nfs_cat`) — < −1.455 low, −1.455 to 0.676 indeterminate, > 0.676 high <sup>[11]</sup>.
- FLI (< 30 rules out, ≥ 60 rules in steatosis) <sup>[12]</sup> and HSI (< 30 rules out,
  > 36 rules in) <sup>[13]</sup> are kept continuous.

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
(high FFMI loss); the remainder **"modest/low"**. This is the grouping
(`perc_change_ffmi_v4_group`) used to stratify Table 1 (`1a`) and the
microbiome descriptive analyses (`2a`–`2c`).

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
(above), and diabetes medication and **statins** enter the `1a` Table 1
variable set.

---

## Table 1 (1a)

Baseline and 1-year characteristics, stratified by the **1-year %FFMI change
group** (`perc_change_ffmi_v4_group`, see *FFMI trajectory groups* above):
demographics, body composition, glycaemic/lipid/renal labs, and
lipid-lowering/antidiabetic medication use. Continuous variables are reported
as mean (SD) or median [IQR] (`tableone`; non-normal variables tested with
Wilcoxon), categorical variables with χ². A combined baseline + 1-year table
additionally runs paired baseline→1-year tests (paired t-test / Wilcoxon /
McNemar) within each FFMI-change group, alongside the between-group tests at
each visit — four comparisons per row, flagged with superscript letters and
significance stars.

Output as CSV (`results/tables/t1_ffmi_loss_*.csv`) and as a styled `gt` HTML
table (`results/tables/t1_ffmi_loss_v0_v4_gt.html`).

---

## Microbiome analysis

### Descriptive analyses (2a–2c)

- Composition, alpha and beta diversity use the **full shotgun table without
  prefiltering**; per-sample relative abundances sum to ~100%.
- Alpha diversity: Shannon, Simpson and observed richness. Beta diversity:
  Bray–Curtis distances with PCoA. Composition: the top 20 species.
- All three are stratified by the **1-year %FFMI change group** (high vs.
  modest/low loss — the same grouping as `1a`'s Table 1), rather than by
  baseline FFMI status.

The baseline-FFMI-status versions of these three scripts, and a standalone
mixed-model screen relating baseline species abundance to FFMI trajectory
across the whole filtered species set (LINDA/LM/LMM), have been superseded by
the ML-first workflow below (`3a`–`3c`) and are kept only in
`scripts/archive/`. Confirmatory LM/LMM testing of baseline-abundance effects
on FFMI now runs on each XGBeast model's top-15 features (`3c`), rather than
as an independent full-species screen.

---

## Machine-learning pipeline (3a–3c)

A gradient-boosting layer that lets all baseline species compete jointly
(rather than one at a time) and ranks them by predictive contribution, using
**XGBeast** (`scripts/assets/XGBeast_new.py`) — a repeated shuffle-split
**XGBoost** regression framework with grid-search hyperparameter tuning that
reports explained variance and permutation feature importance across CV
iterations. Confirmatory classical modelling of each model's top-ranked
features is folded into the same pipeline (`3c`), replacing the standalone
mixed-model species screen described above.

- **`3a_ml_prep.R`** builds one input table per outcome × subgroup
  (`all`/`male`/`female`) under `results/mlmodels/<outcome>/<subgroup>/`.
  Predictors are baseline (`v0`) species, kept if present at >0.05% relative
  abundance in ≥20% of samples (or the equivalent `v4` species table for the
  two `v4` outcomes). The 11 outcomes are: baseline FFMI/FMI, baseline FFMI
  adjusted for FMI (residualised), the same two at `v4` (cross-sectional,
  predicted from **`v4` species**, not baseline), baseline FFMI/FMI
  restricted to the subjects with `v4` microbiome data (`*_matched`, to
  isolate subject-subset effects from timepoint effects), and 1-year change
  in FFMI (`delta_ffmi_v4`, `perc_change_ffmi_v4`) with and without FMI
  adjustment.
- **`3b_ml_run.sh`** is a SLURM batch script that calls `XGBeast_new.py` once
  per outcome × subgroup (11 outcomes × 3 subgroups) with a grid search
  (`param_grid.json`) and a fixed random seed.
- **`3c_ml_process.R`** reads each model's XGBeast output and, per
  outcome/subgroup: plots the explained-variance distribution across CV
  iterations, extracts the top-15 features by importance, plots per-feature
  abundance-vs-outcome correlation panels as a sanity check, and **confirms
  them** with a classical model on the same predictors — an LM for
  cross-sectional outcomes (age/sex-adjusted, plus FMI where applicable) or
  an LMM (baseline-abundance × visit interaction, adjusted for age, sex,
  %weight change, and FMI for `delta_ffmi_v4`) for change outcomes —
  reporting a forest plot of per-feature effect estimates with BH-FDR across
  the top-15 features. It finishes by plotting the explained-variance
  distribution across CV iterations for the six headline FFMI/delta-FFMI
  models (± FMI adjustment) to `results/graphs/ml_explained_variance/`.

---

## Correlation analyses (4a–4c)

All three scripts take the **top-15 species by feature importance** from the
`perc_change_ffmi_v4`/`all` XGBeast model (`3b`'s output) — i.e. the baseline
species most predictive of 1-year FFMI change — and correlate their
`log10`-transformed baseline relative abundance against an external baseline
dataset, using Spearman correlations. Species labels in every heatmap are
coloured by the direction of that species' effect on FFMI change (`3c`'s
forest-plot estimate for this model).

- **`4a_humann_pathways.R`** — correlates the top-15 species against baseline
  HUMAnN pathways retained at >5 CPM in ≥30% of baseline samples, with
  Benjamini–Hochberg FDR applied per species across pathways. Writes a
  species × pathway `ComplexHeatmap` (all FDR-significant pairs, plus a
  compact top-15-pathway version) to `results/graphs/HUMAnN/`, and saves the
  set of implicated pathways to `data/processed_data/HUMAnN_selected_pathways.RDS`.
- **`4b_mb_metabolome_correlations.R`** — correlates the top-15 species
  against every baseline serum metabolite in `BARIA_metabolon_clean.RDS`,
  with Benjamini–Hochberg FDR across all species × metabolite pairs. Writes
  the full correlation table, a scatter plot per species (one panel per
  FDR-significant metabolite), and `ComplexHeatmap` species × metabolite
  heatmaps — all FDR-significant metabolites, a top-25-by-hit-count compact
  version, and a further compact version restricted to species with ≥10
  significant hits among those top metabolites — annotated by Metabolon
  super-pathway, to `results/graphs/mb_metabolome_correlations/`.
- **`4c_mb_diet_correlations.R`** — correlates the top-15 species against
  baseline energy-normalised macronutrient intake (carbohydrate, protein,
  fat, saturated fat, fibre, alcohol, plus total energy) from
  `251002_BARIA_macronutrients.RDS`, run separately per dietary-diary tool
  and once pooled ("overall"). Species are additionally required to be
  present in **≥ 30%** of participants within each group's subset before
  testing. Because little survives FDR correction at this sample size,
  significance here is reported **nominal** (uncorrected `p < 0.05`), not
  FDR-corrected. Outputs (correlation table, scatter plots, heatmap) are
  written per group to `results/graphs/mb_diet_correlations/<group>/`.

---

## Tissue RNA-seq (0c, 5a–5b) — in progress

`0c_datacleaning_rnaseq.R` cleans the raw kallisto count table into one
VST-normalised dataset per tissue (DESeq2 size factors, an expression filter
of ≥10 normalised reads in ≥50% of participants, variance-stabilising
transformation; see **Data**), after a per-tissue sample QC:

- **QC metrics.** For each tissue, library size (total raw reads per sample)
  and a principal component analysis on the VST of all samples, restricted to
  the 500 most variable genes, centred and not scaled — the convention of
  DESeq2's `plotPCA()` (`ntop = 500`) <sup>[19, 20]</sup>. The % variance shown
  on the PCA axes refers to these 500 genes.
- **Exclusion.** Samples forming a distinct cluster outside the main cloud on
  the QC PCA, together with small library size, were excluded on visual
  inspection, following the approach used for the same BARIA tissue RNA-seq by
  Zwartjes et al. <sup>[21]</sup> (see also <sup>[22]</sup>). This applied to
  liver only: 20 samples formed a separate low-library-size cluster (the same
  participants' visceral adipose samples did not separate, indicating a
  liver-specific sample/library issue). Jejunum, visceral and subcutaneous fat
  showed no distinct cluster and had no exclusions. Retained samples: liver
  340, jejunum 212, visceral fat 227, subcutaneous fat 220.
- **Outputs.** QC plots per tissue (`qc_librarysize_<tissue>.pdf`,
  `qc_pca_<tissue>.pdf`, excluded samples labelled) and one QC table across
  tissues (`rnaseq_qc_samples.csv`: library size, size factor, PC1/PC2, % variance,
  exclusion flag). Excluded samples are removed before the size factors,
  expression filter and VST are recomputed, so they do not influence the
  normalisation of retained samples; the raw-count files keep all samples with
  a `qc_exclude` flag.

`5a_liver_rnaseq.R` and `5b_vadipose_rnaseq.R`
each start from one tissue's VST-normalised expression and the same top-15
species (`perc_change_ffmi_v4`/`all` model) used in the correlation analyses
above, and run the same two-part analysis: an untargeted genome-wide Spearman
screen (species × all expressed genes, BH-FDR per species, with the top 25
annotated genes by |rho| among FDR-significant hits taken forward to a
`ComplexHeatmap`) and a targeted screen restricted to a curated, tissue-specific
gene panel (grouped into biological categories, each tested the same way, with
individual correlation scatter plots for FDR-significant pairs). The liver
panel covers hepatokines, fatty-acid oxidation, mitochondrial OXPHOS,
oxidative stress, BCAA turnover, gluconeogenesis/glycogen metabolism,
glycolysis, de-novo lipogenesis/VLDL, insulin signaling and glucagon/incretin
signaling. The visceral adipose panel covers adipokines, adipogenesis,
lipolysis/lipid storage, lipogenesis/fatty-acid uptake, mitochondrial OXPHOS,
insulin signaling, macrophage/inflammation, fibrosis/ECM, oxidative stress,
BCAA turnover and glycolysis — browning/thermogenesis markers were
deliberately left out, since omental/mesenteric visceral fat has minimal
beiging potential in adult humans (browning localises mainly to subcutaneous
and supraclavicular/cervical/perirenal depots). Both scripts guard against
the case where no gene (untargeted) or no targeted gene reaches FDR < 0.05 for
any species: rather than erroring, they skip the corresponding heatmap(s) and
print the closest-to-significant associations by raw p-value instead.

---

## Dependencies

The environment is managed with **pixi** (`pixi.toml`, `pixi.lock`), and defines
pixi tasks for every script plus grouped tasks (`datacleaning`, `mb`, `ml`,
`mb-correlations`, `all`) that chain the steps within each stage — run e.g.
`pixi run ml` or `pixi run datacleaning-clinical`; see `pixi.toml` for the
full list. `liver-rnaseq` (`5a`) is defined but excluded from
`all`, since that script is still in progress; `5b_vadipose_rnaseq.R` has no
pixi task yet and is run directly (`Rscript scripts/5b_vadipose_rnaseq.R`). R
packages loaded across the scripts include: `tidyverse`, `phyloseq`, `vegan`,
`tableone`, `gt`, `ggpubr`, `patchwork`, `ggthemes`, `ggsci`, `ggrepel`,
`ape`, `grid`, `MetBrewer`, `lmerTest`, `broom` / `broom.mixed`,
`ComplexHeatmap` / `circlize`, and (`5a`, `5b`) `annotables` and `DESeq2`
(the latter via `0c`; `0c` also uses `matrixStats`, installed as a DESeq2
dependency, and `ggrepel`). The ML pipeline (`3b`) additionally uses Python:
`xgboost`, `scikit-learn`, `numpy`, `pandas`, `matplotlib`, `seaborn`, `shap`,
`tqdm`.

> **Note:** `ggrepel` and `annotables` are not currently pinned in
> `pixi.toml` — install them separately until the manifest is updated.
> `ComplexHeatmap`, `circlize`, `gt`, `ape`, `ggsci` and `DESeq2` *are*
> already pinned, unlike in an earlier version of this README.

Expected working-directory layout (script paths are relative to the repo root):

```
data/raw_data/          # inputs (see Data)
data/processed_data/    # created by 0a (clinical), 0b (omics) and 0c (RNA-seq)
results/graphs/          # figures, in per-analysis subfolders
results/tables/          # Table 1 output; RNA-seq sample-QC table (0c)
results/mlmodels/        # ML input data + XGBeast output, per outcome/subgroup
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
8. Sterling RK, Lissen E, Clumeck N, et al. Development of a simple noninvasive index
   to predict significant fibrosis in patients with HIV/HCV coinfection. *Hepatology.*
   2006;43(6):1317–1325. doi:10.1002/hep.21178.
9. Shah AG, Lydecker A, Murray K, et al. Comparison of noninvasive markers of fibrosis
   in patients with nonalcoholic fatty liver disease. *Clin Gastroenterol Hepatol.*
   2009;7(10):1104–1112. doi:10.1016/j.cgh.2009.05.033.
10. McPherson S, Hardy T, Dufour JF, et al. Age as a confounding factor for the
    accurate non-invasive diagnosis of advanced NAFLD fibrosis. *Am J Gastroenterol.*
    2017;112(5):740–751. doi:10.1038/ajg.2016.453.
11. Angulo P, Hui JM, Marchesini G, et al. The NAFLD fibrosis score: a noninvasive
    system that identifies liver fibrosis in patients with NAFLD. *Hepatology.*
    2007;45(4):846–854. doi:10.1002/hep.21496.
12. Bedogni G, Bellentani S, Miglioli L, et al. The Fatty Liver Index: a simple and
    accurate predictor of hepatic steatosis in the general population. *BMC
    Gastroenterol.* 2006;6:33. doi:10.1186/1471-230X-6-33.
13. Lee JH, Kim D, Kim HJ, et al. Hepatic steatosis index: a simple screening tool
    reflecting nonalcoholic fatty liver disease. *Dig Liver Dis.* 2010;42(7):503–508.
    doi:10.1016/j.dld.2009.08.002.
14. Simental-Mendía LE, Rodríguez-Morán M, Guerrero-Romero F. The product of fasting
    glucose and triglycerides as surrogate for identifying insulin resistance in
    apparently healthy subjects. *Metab Syndr Relat Disord.* 2008;6(4):299–304.
    doi:10.1089/met.2008.0034.
15. Abdul-Ghani MA, Matsuda M, Balas B, DeFronzo RA. Muscle and liver insulin
    resistance indexes derived from the oral glucose tolerance test. *Diabetes Care.*
    2007;30(1):89–94. doi:10.2337/dc06-1519.
16. Matsuda M, DeFronzo RA. Insulin sensitivity indices obtained from oral glucose
    tolerance testing: comparison with the euglycemic insulin clamp. *Diabetes Care.*
    1999;22(9):1462–1470. doi:10.2337/diacare.22.9.1462.
17. Matsubayashi Y, Yoshida A, Suganami H, et al. Role of fatty liver in the
    association between obesity and reduced hepatic insulin clearance. *Diabetes
    Metab.* 2018;44(2):135–142. doi:10.1016/j.diabet.2017.12.003.
18. Tricò D, Galderisi A, Mari A, et al. Intrahepatic fat, irrespective of ethnicity,
    is associated with reduced endogenous insulin clearance and hepatic insulin
    resistance in obese youths. *Diabetes Obes Metab.* 2020;22(9):1628–1638.
    doi:10.1111/dom.14076.
19. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion
    for RNA-seq data with DESeq2. *Genome Biol.* 2014;15(12):550.
    doi:10.1186/s13059-014-0550-8.
20. Love MI, Anders S, Huber W. Analyzing RNA-seq data with DESeq2 (package
    vignette, section "Principal component plot of the samples").
    https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html;
    `plotPCA()` source: https://github.com/thelovelab/DESeq2/blob/devel/R/plots.R.
21. Zwartjes MSZ, de Jonge PA, van de Laar AW, et al. Adipose tissue inflammation,
    oxidative stress, and altered adipogenesis are associated with dyslipidemia in
    obesity: a multiomics profiling study. *J Am Heart Assoc.* 2026;15(11):e047397
    (Supplemental Methods). doi:10.1161/JAHA.125.047397.
22. Meijnikman AS, Davids M, Herrema H, et al. Microbiome-derived ethanol in
    nonalcoholic fatty liver disease. *Nat Med.* 2022;28(10):2100–2106.
    doi:10.1038/s41591-022-02016-6.
