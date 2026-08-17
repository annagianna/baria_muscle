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

The code covers data cleaning and the baseline (cross-sectional) analyses —
microbiome composition and diversity, and per-species associations with FFMI —
together with a longitudinal mixed-model layer relating baseline species abundance
to FFMI trajectory after surgery. A functional-pathway layer (HUMAnN) is being
added on top of the taxonomic analyses, starting with pathway cleaning and
prevalence filtering.

---

## Repository structure

```
baria_muscle/
├── scripts/
│   ├── 0a_datacleaning.R              # Clinical + microbiome cleaning, derived variables, cohort definition
│   ├── 1a_tableone.R                  # Table 1, stratified by baseline FFMI status
│   ├── 2a_mb_alphadiversity.R         # Alpha diversity (Shannon, Simpson, observed richness)
│   ├── 2b_mb_betadiversity.R          # Beta diversity (Bray–Curtis PCoA + PERMANOVA / dispersion)
│   ├── 2c_mb_composition.R            # Compositional bar plots (top 20 species, by baseline FFMI group)
│   ├── 4a_lm_mb_baseline_ffmi.R       # Cross-sectional linear models: species abundance ~ baseline FFMI
│   ├── 4b_lmm_mb_ffmi_trajectories.R  # Mixed models: baseline abundance × FFMI trajectory + overlap with 4a
│   ├── 4c_lm_lmm_ffmi_dags.R          # DAGs + adjustment sets for the models in 4a / 4b
│   └── 5a_humann_pathways.R           # Functional pathways (HUMAnN): cleaning + prevalence filtering
├── data/
│   ├── raw_data/                      # Inputs (not tracked — see Data)
│   └── processed_data/                # Cleaned datasets written by 0a
├── tables/                            # Table 1, saved model-result objects, overlap tables
├── graphs/                            # Figures, in per-analysis subfolders
│   ├── alphadiversity/                #   2a
│   ├── betadiversity/                 #   2b
│   ├── composition/                   #   2c
│   ├── microbe_associations/          #   4a combined forest plot
│   └── lmm_species_ffmi/              #   4b trajectory / overlap plots
└── README.md
```

### Order of execution

`0a` runs first — it writes the processed datasets every other script reads.
`1a` and the `2*` descriptive scripts are independent of one another. The `4*`
association scripts depend on the processed data from `0a`; `4c` defines the DAGs
and adjustment sets that the models in `4a` / `4b` implement. `4a` runs before
`4b`: it saves the baseline significant-species object that `4b` reads to
cross-reference the cross-sectional and longitudinal hits. `5a` reads the raw
HUMAnN profiles together with the processed data from `0a`, and is otherwise
independent of the taxonomic scripts.

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

Visit coding: `v0` = baseline (raw `V-1` / MetaPhlAn `Time_Point == "V-1"`),
`v2`–`v5` = yearly follow-ups, `v6` = 5 years, `v7` = 10 years.

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
analysis (BIA)** (Maltron BioScan 920). Raw resistance, reactance, impedance and
phase angle are recorded at 5, 50, 100 and 200 kHz; the **50 kHz resistance** is
the value carried into the fat-free-mass / skeletal-muscle-mass equations, as is
standard for whole-body prediction.

### Plausibility filtering

Two layers of checks are applied before a baseline BIA record is accepted into the
cohort:

**1. Entry-range checks (upstream, at data capture).**
The data-capture system enforces plausible entry ranges on the primary anthropometry that feeds
the body-composition indices — weight is bounded to 50–250 kg and height to
100–300 cm — so grossly out-of-range height/weight are caught at source.

**2. Internal-consistency and device-range checks (in `0a`).**
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
source record. Two correction types are permitted:

- **Unit / decimal-scale fixes** — an isolated field off by a consistent factor of
  10 or 100 (a misplaced decimal or wrong unit scale) is rescaled by that factor
  so that it agrees with the participant's other measures and restores percentage
  and mass closure.
- **Source-verified replacements** — a small number of clearly erroneous entries
  are replaced with the correct value read back from the source BIA record.

Values that are plausible but extreme (true statistical outliers) are **not**
edited; they are retained and handled analytically (e.g. sensitivity analysis)
rather than corrected. All manual corrections are hard-coded per participant-visit
in `0a` and flagged in-line for supervisor review. No participant identifiers are
reproduced in this document.

---

## Formulas

All formulas applied in `0a_datacleaning.R`. Superscripts refer to the
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
- **Insulin unit conversion** — `insulin_µU/mL = insulin_pmol/L / 6.945`
  (standard SI conversion, 1 µU/mL = 6.945 pmol/L)
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

- Numeric sentinel values **-99, -98, -97** are recoded to `NA` (they denote
  missingness for different reasons in the BARIA export).
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

- **T2D (`t2d_labs`)** — HbA1c **≥ 6.5%** (≥ 48 mmol/mol) **or** fasting glucose
  **≥ 7.0 mmol/L**.
- **Prediabetes (`prediab_labs`)**, among non-T2D — HbA1c **5.7–6.4%**
  (39–47 mmol/mol) **or** fasting glucose **5.6–6.9 mmol/L** (IFG).
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
(high FFMI loss); the remainder **"low/modest"**.

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
(above), and **statins** enter the extended sensitivity adjustment set.

---

## Microbiome analysis

### Descriptive analyses (`2a`–`2c`)

- Composition, alpha and beta diversity use the **full shotgun table without
  prefiltering**; per-sample relative abundances sum to ~100%.
- Alpha diversity: Shannon, Simpson and observed richness. Beta diversity:
  Bray–Curtis distances with PCoA. Composition: the top 20 species, stratified by
  baseline FFMI group.

### Association models (`4a`, `4b`)

- **Species set (for the models only).** Poorly annotated taxa (uncharacterised
  `GGB` genome bins) are dropped first. A species is then kept if detected in
  **≥ 20%** of baseline samples **and** its mean relative abundance is
  **≥ 0.01%**. Abundances enter the models on a `log10` scale with a
  species-specific half-minimum pseudocount, `log10(abundance + min_abundance/2)`.
  p-values are FDR-adjusted (Benjamini–Hochberg) across species within each model.
- **Cross-sectional linear models (`4a`).** Fitted per species as
  `FFMI_v0 ~ log10_abundance + covariates`, in three nested specifications that
  implement the three DAGs in `4c`:
  - *Model 1* — `+ age + sex`
  - *Model 2* (adiposity sensitivity) — `+ age + sex + FMI`
  - *Model 3* (extended sensitivity) — `+ age + sex + FMI + T2D + antidiabetic
    medication + statins`

  The three models are drawn together in a combined forest plot
  (`graphs/microbe_associations/`), and the Model-1 significant-species object is
  saved to `tables/` for use by `4b`.
- **Longitudinal mixed models (`4b`).** Baseline abundance against FFMI
  trajectory, fitted per species with `lmerTest` (ML) as
  `FFMI ~ log10_baseline_abundance * visit + age_centred + sex + %weight_change +
  (1 | id)` (random intercept per participant). The term of interest is the
  `log10_baseline_abundance × visit` interaction — whether baseline abundance
  predicts the *change* in FFMI. The model is run separately for each
  post-surgical follow-up (`v4`, `v5`; labelled 1-year and 2-year in the figures),
  and species are flagged as significant at each interval, at **both** intervals
  (overlap), or at one only.
- **Cross-sectional × longitudinal overlap (`4b`).** The longitudinal overlap
  hits are then cross-referenced against the baseline hits from `4a`, giving a
  combined table of baseline, 1-year and 2-year effect sizes
  (`tables/lm_lmm_ffmi_overlap_table.csv`). For those species, `4b` also draws
  observed-trajectory and %-FFMI-change plots comparing high vs low
  baseline-abundance tertile groups (`graphs/lmm_species_ffmi/`).

### DAG-based adjustment (`4c`)

Adjustment sets are **DAG-derived**, blocking back-door paths from microbe → FFMI
without conditioning on mediators or colliders. Insulin resistance is treated as a
co-outcome rather than a confounder; BMI is avoided as an over-adjustment, since
FFMI already height-normalises. Model 2 adds adiposity (FMI); Model 3 additionally
conditions on T2D, antidiabetic medication and statins.

### Functional pathways (`5a`)

Community-level functional-pathway profiles from **HUMAnN 4.0** (MetaCyc pathways,
HUMAnN's native output), the encoded-function layer that sits on top of the
taxonomic analyses. This script currently covers cleaning and filtering; modelling
against FFMI is not yet implemented.

- Profiles are reshaped to long form; sample names are harmonised to the MetaPhlAn
  convention and restricted to **faecal** samples present in the cleaned cohort.
- The `UNMAPPED` and `UNINTEGRATED` categories are removed, and each pathway string
  is split into its identifier and name.
- Undetected pathways are set to zero. A pathway is retained if it exceeds
  **5 CPM** in **≥ 50%** of baseline samples, and abundances are `log10`-transformed
  with a `+1` pseudocount.

---

## Dependencies

R, with (across scripts): `tidyverse`, `phyloseq`, `vegan`, `ggpubr`, `tableone`,
`lmerTest`, `broom` / `broom.mixed`, `dagitty`, `ggdag`, `ggrepel`, `ggthemes`,
`patchwork`, `grid`, and the palette packages `MetBrewer` / `ggsci`.

Expected working-directory layout (script paths are relative to the repo root):

```
data/raw_data/         # inputs (see Data)
data/processed_data/   # created by 0a
```

---

## References

1. Kawakami R, Tanisawa K, Ito T, et al. Fat-Free Mass Index as a surrogate marker
   of appendicular skeletal muscle mass index for low muscle mass screening in
   sarcopenia. *J Am Med Dir Assoc.* 2022;23(12):1955–1961.e3.
   doi:10.1016/j.jamda.2022.08.016. (Defines FFMI = FFM/height²; basis for using
   BIA FFMI as the muscle-mass proxy.)
2. VanItallie TB, Yang MU, Heymsfield SB, Funk RC, Boileau RA. Height-normalized
   indices of the body's fat-free mass and fat mass: potentially useful indicators
   of nutritional status. *Am J Clin Nutr.* 1990;52(6):953–959.
   doi:10.1093/ajcn/52.6.953. (Origin of the height-normalized FFM and FM indices.)
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
   of hemoglobin A1c in human blood and the national standardization schemes in the
   United States, Japan, and Sweden: a method-comparison study. *Clin Chem.*
   2004;50(1):166–174. doi:10.1373/clinchem.2003.024802. (Basis of the IFCC–NGSP
   master equation, `NGSP% = 0.09148 × IFCC + 2.152`.)
7. American Diabetes Association. Standards of Care in Diabetes — 2026.
   *Diabetes Care.* 2026;49(Suppl 1). (Diagnostic thresholds for T2D and prediabetes.)
