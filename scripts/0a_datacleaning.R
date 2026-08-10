# Baria project: Muscle mass trajectories and gut microbiota following bariatric surgery - Data cleaning
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)
library(phyloseq)

# Open Data and see properties
baria_clinical_data_raw <- readRDS("./data/raw_data/BARIA.clinical.2024-12-09.723.2043.RDS")
baria_mb <- readRDS("data/raw_data/ps.BARIA.metaphlan.706.2548.RDS")

## Cut-offs used in this script (qc)

## Formulas used in this script
# Hba1c(%) = (0,0915 * HbA1c (mmol/mol) + 2,15 (from diabetesfonds.nl)
# Hba1c(mmol/mol) = (10,93 x Hba1c (%)) - 23,5

# DOI: 10.2337/diacare.21.12.2191 for HOMA calculations
# HOMA-IR = (Fasting insulin (µU/mL) × Fasting glucose (mmol/L)) / 22.5 (for glucose in mmol/L)
# HOMA-B = (20 × Fasting insulin (µU/mL)) / (Fasting glucose (mmol/L) − 3.5)
# Insulin was converted from pmol/l to μU/ml

# Skeletal muscle mass (SMM) by Janssen: SMM [kg] = (height^2 [cm] / BIA-resistance [Ohms] X 0.401) + (gender x 3.825) + (age [years] x - 0.071)] + 5.102 (men = 1; women = 0)

# Longitudinal vars
long_vars <- c("date",
  # Body composition
  "bmi", "weight_kg", "wc_cm", "fm_kg", "fm_percent", "ffm_kg", "ffm_percent", "bia_resistance_50khz", "upperleg_cm",

  # Lab/glycemia parameters
  "fasting_glucose_mmoll_mmt", "fasting_insulin_pmoll_mmt", "fasting_glucagon_ngl_mmt", "fasting_cpeptide_nmoll_mmt",
  "hba1c", "hba1c_mmolmol", "crp_mgl",

  # Nexfin
  "nexfin_hr", "nexfin_dpdt", "nexfin_sv", "nexfin_svi", "nexfin_co", "nexfin_ci", "nexfin_svr", "nexfin_svri"
)

long_vars_pattern <- str_c("^(", str_c(long_vars, collapse = "|"), ")_(v\\d+)$")

# Clinical data
baria_muscle_vars <- baria_clinical_data_raw |>
  select(
    id = Subject_ID,

    # Visits
    date_v0 = date,
    date_v2 = V2_date,
    date_v3 = V3_date,
    date_v4 = V4_date,
    date_v5 = V5_date,
    date_v6 = V6_date_1,
    date_v7 = V7_date,
    sg_type = type_surgery,

    # Baseline (v0)
    age_v0 = Age,
    sex,
    t2d_v0 = dm,
    bmi_v0 = bmi,
    weight_kg_v0 = weight,
    height_cm = height,
    wc_cm_v0 = taille,
    fm_kg_v0 = tbf,
    fm_percent_v0 = tbf_percent,
    ffm_kg_v0 = ffm,
    ffm_percent_v0 = ffm_percent,
    bia_resistance_50khz_v0 = rawdata_50Khz_Resistance,
    upperleg_cm_v0 = upperleg,
    aht = hypertension,
    medication_binary_v0 = meds,
    medication_list_v0 = medication_v0_freetext,
    sport_v0 = scre_sport,
    systolic_bp_mmhg_v0 = systolic_pressure_v0,
    diastolic_bp_mmhg_v0 = diastolic_pressure_v0,
    glucose_mmoll_v0 = glucose,
    fasting_glucose_mmoll_mmt_v0 = min0gluc,
    fasting_insulin_pmoll_mmt_v0 = min0insulin,
    fasting_glucagon_ngl_mmt_v0 = min0glucagon,
    fasting_cpeptide_nmoll_mmt_v0 = min0cpept,
    hba1c_v0 = hba1c,
    hba1c_mmolmol_v0 = hba1c__IFCC_mmolmol,
    gammagt_ul_v0 = ggt,
    alat_ul_v0 = alat,
    asat_ul_v0 = asat,
    triglycerides_mmoll_v0 = triglycerides,
    crp_mgl_v0 = crp,
    tsh_miul_v0 = tsh,
    ft4_pmoll_v0 = ft4,
    nexfin_hr_v0 = nexfin_hr_v0,
    nexfin_dpdt_v0 = nexfin_dpdt_v0,
    nexfin_sv_v0 = nexfin_sv_v0,
    nexfin_svi_v0 = nexfin_svi_v0,
    nexfin_co_v0 = nexfin_CO_v0,
    nexfin_ci_v0 = nexfin_CI_v0,
    nexfin_svr_v0 = nexfin_svr_v0,
    nexfin_svri_v0 = nexfin_svri_v0,

    # Follow-Up Data
    # v2 = 6 weeks
    bmi_v2 = V2_bmi,
    weight_kg_v2 = V2_weight,
    wc_cm_v2 = V2_taille,
    upperleg_cm_v2 = V2_upperleg,
    hba1c_v2 = V2_hba1c,
    hba1c_mmolmol_v2 = V2_hba1c__IFCC_mmolmol,
    crp_mgl_v2 = V2_crp,
      
    # v3 = 6 months
    bmi_v3 = V3_bmi,
    weight_kg_v3 = V3_weight,
    wc_cm_v3 = V3_taille,
    upperleg_cm_v3 = V3_upperleg,
    hba1c_v3 = V3_hba1c,
    hba1c_mmolmol_v3 = V3_hba1c__IFCC_mmolmol,
    crp_mgl_v3 = V3_crp,

    # v4 = 1 year
    bmi_v4 = V4_bmi,
    weight_kg_v4 = V4_weight,
    wc_cm_v4 = V4_taille,
    fm_kg_v4 = V4_tbf,
    fm_percent_v4 = V4_tbf_percent,
    ffm_kg_v4 = V4_ffm,
    ffm_percent_v4 = V4_ffm_percent,
    bia_resistance_50khz_v4 = V4_rawdata_50Khz_Resistance,
    upperleg_cm_v4 = V4_upperleg,
    fasting_glucose_mmoll_mmt_v4 = V4_min0gluc,
    fasting_insulin_pmoll_mmt_v4 = V4_min0insulin,
    fasting_glucagon_ngl_mmt_v4 = V4_min0glucagon,
    fasting_cpeptide_nmoll_mmt_v4 = V4_min0cpept,
    hba1c_v4 = V4_hba1c,
    hba1c_mmolmol_v4 = V4_hba1c__IFCC_mmolmol,
    crp_mgl_v4 = V4_crp,
    nexfin_hr_v4 = nexfin_hr_v4,
    nexfin_dpdt_v4 = nexfin_dpdt_v4,
    nexfin_sv_v4 = nexfin_sv_v4,
    nexfin_svi_v4 = nexfin_svi_v4,
    nexfin_co_v4 = nexfin_co_v4,
    nexfin_ci_v4 = nexfin_ci_v4,
    nexfin_svr_v4 = nexfin_svr_v4,
    nexfin_svri_v4 = nexfin_svri_v4,

    # v5 = 2 years
    bmi_v5 = V5_bmi,
    weight_kg_v5 = V5_weight,
    wc_cm_v5 = V5_taille,
    fm_kg_v5 = V5_tbf,
    fm_percent_v5 = V5_tbf_percent,
    ffm_kg_v5 = V5_ffm,
    ffm_percent_v5 = V5_ffm_percent,
    bia_resistance_50khz_v5 = V5_rawdata_50Khz_Resistance,
    upperleg_cm_v5 = V5_upperleg,
    fasting_glucose_mmoll_mmt_v5 = V5_min0gluc,
    fasting_insulin_pmoll_mmt_v5 = V5_min0insulin,
    fasting_glucagon_ngl_mmt_v5 = V5_min0glucagon,
    fasting_cpeptide_nmoll_mmt_v5 = V5_min0cpept,
    hba1c_v5 = V5_hba1c,
    hba1c_mmolmol_v5 = V5_hba1c__IFCC_mmolmol,
    crp_mgl_v5 = V5_crp,
    nexfin_hr_v5 = nexfin_hr_v5,
    nexfin_dpdt_v5 = nexfin_dpdt_v5,
    nexfin_sv_v5 = nexfin_sv_v5,
    nexfin_svi_v5 = nexfin_svi_v5,
    nexfin_co_v5 = nexfin_co_v5,
    nexfin_ci_v5 = nexfin_ci_v5,
    nexfin_svr_v5 = nexfin_svr_v5,
    nexfin_svri_v5 = nexfin_svri_v5,

    # v6 = 5 years
    bmi_v6 = V6_bmi_1,
    weight_kg_v6 = V6_weight_1,
    wc_cm_v6 = V6_taille_1,
    fm_kg_v6 = V6_tbf_1,
    fm_percent_v6 = V6_tbf_percent_1,
    ffm_kg_v6 = V6_ffm_1,
    ffm_percent_v6 = V6_ffm_percent_1,
    bia_resistance_50khz_v6 = V6_rawdata_1_50Khz_Resistance,
    upperleg_cm_v6 = V6_upperleg_1,
    hba1c_v6 = V6_hba1c_1,
    hba1c_mmolmol_v6 = V6_hba1c_IFCC_mmolmol,
    crp_mgl_v6 = V6_crp_1,
    nexfin_hr_v6 = nexfin_hr_v6_1,
    nexfin_dpdt_v6 = nexfin_dpdt_v6_1,
    nexfin_sv_v6 = nexfin_sv_v6_1,
    nexfin_svi_v6 = nexfin_svi_v6_1,
    nexfin_co_v6 = nexfin_co_v6_1,
    nexfin_ci_v6 = nexfin_ci_v6_1,
    nexfin_svr_v6 = nexfin_svr_v6_1,
    nexfin_svri_v6 = nexfin_svri_v6_1,

    # v7 = 10 years
    bmi_v7 = V7_BMI,
    weight_kg_v7 = V7_weight,
    wc_cm_v7 = V7_taille,
    upperleg_cm_v7 = V7_upperleg,
    hba1c_v7 = V7_hba1c,
    hba1c_mmolmol_v7 = V7_hba1c_IFCC_mmolmol,
    crp_mgl_v7 = V7_CRP
  ) |>
   mutate(
    across(everything(), ~ replace(.x, .x %in% c(-99, -98, -97), NA)), # these values are NA for different reasons in the Baria dataset
    across(.cols = starts_with("date"), ~ if_else(.x %in% c("01-01-2999", "01-01-2997", "01-01-2995"), NA_character_, .x)),
    sex = case_when(sex == "1" ~ "male", sex == "2" ~ "female"),
    t2d_v0 = case_when(t2d_v0 == "1" ~ "yes", t2d_v0 == "2" ~ "no"),
    aht = case_when(aht == "1" ~ "yes", aht == "2" ~ "no"),
    medication_binary_v0 = case_when(medication_binary_v0 == "1" ~ "yes", medication_binary_v0 == "2" ~ "no"),
    sport_v0 = case_when(sport_v0 == "1" ~ "yes", sport_v0 == "2" ~ "no"),
    sg_type = case_when(sg_type == "1" ~ "rygb", sg_type == "2" ~ "omegaloop", sg_type == "3" ~ "sg")
  ) |>
  pivot_longer(
    cols = matches(long_vars_pattern),
    names_to = c(".value", "visit"), 
    names_pattern = "(.+)_(v\\d+)$"
  ) |> 
  mutate(date = dmy(date)) |> 
  group_by(id) |> 
  mutate(
    date_baseline = date[visit == "v0"],
    n_years_from_v0 = as.numeric(date - date_baseline) / 365.25,
    age = age_v0 + n_years_from_v0
  )  |> 
  ungroup() |> 
  mutate(
    hba1c_percent = if_else(hba1c < 15, hba1c, hba1c * 0.0915 + 2.15),
    hba1c_mmolmol = if_else(is.na(hba1c_mmolmol) == FALSE, hba1c_mmolmol, 10.93 * hba1c_percent - 23.5),

    # HOMA-IR & HOMA-2B (insulin unit conversion from pmol/l to uU/ml)
    homa_ir = (fasting_insulin_pmoll_mmt / 6.945) * fasting_glucose_mmoll_mmt / 22.5,
    homa_b = (20 * (fasting_insulin_pmoll_mmt / 6.945)) / (fasting_glucose_mmoll_mmt - 3.5),

    # T2D prevalence at follow-up: i. A1C ≥ 6.5% (≥ 48 mmol/mol) OR ii. FPG ≥ 126 mg/dL (≥ 7.0 mmol/L) (2h OGTT or random plasma glucose not measured)
    t2d_labs = case_when(
      is.na(hba1c_percent) & is.na(fasting_glucose_mmoll_mmt) ~ NA_character_,
      hba1c_percent >= 6.5 | fasting_glucose_mmoll_mmt >= 7.0 ~ "yes",
      TRUE ~ "no"
    ),

    # Prediabetes (ADA SOC 2026): i. A1C 5.7–6.4% (39–47 mmol/mol) OR ii. FPG 100 mg/dL (5.6 mmol/L) to 125 mg/dL (6.9 mmol/L) (IFG) OR 2h gluc (OGTT not available)
    prediab_labs = case_when(
      t2d_labs == "yes" ~ "no",
      is.na(hba1c_percent) & is.na(fasting_glucose_mmoll_mmt) ~ NA_character_,
      (hba1c_percent >= 5.7 & hba1c_percent <= 6.4) | (fasting_glucose_mmoll_mmt >= 5.6 & fasting_glucose_mmoll_mmt <= 6.9) ~ "yes",
      TRUE ~ "no"
    ),

    # de novo occurence(1 & 2y post-surgery, if NGT at baseline/previous visits) 
    # This I need to fix
    # denovo_prediab = case_when(
     # is.na(t2d_v0) | (is.na(prediab_labs) & visit == "v0") | (is.na(prediab_labs) & visit == "v0") ~ NA_character_,
      #(t2d_v0 == "no" & prediab_v0 == "no" & prediab_v4 == "yes") ~ "yes",
      #TRUE ~ "no")
  ) |>
  select(-date_baseline, -age_v0) |> 
  ungroup() |> 
  mutate(# manual corrections for incorrect BIA data entries
    id = as.character(id),
    fm_kg = case_when(
      id == "212" & visit == "v0" ~ fm_kg / 100,
      id == "600" & visit == "v0" ~ fm_kg * 10,
      id == "30004" & visit == "v4" ~ fm_kg / 100,
      id == "567" & visit == "v4" ~ fm_kg / 10,
      id == "335" & visit == "v4" ~ fm_kg * 10,
      id == "612" & visit == "v0" ~ 52.0,
      id == "23" & visit == "v5" ~ 13.8,
      TRUE ~ fm_kg
    ),
    ffm_kg = case_when(
      id == "238" & visit == "v0" ~ ffm_kg * 10,
      id == "30016" & visit == "v0" ~ ffm_kg * 10,
      id == "50" & visit == "v5" ~ 79.16,
      id == "484" & visit == "v4" ~ 51.3,
      id == "520" & visit == "v4" ~ 57.5,
      id == "451" & visit == "v4" ~ 39.7,
      TRUE ~ ffm_kg
    ),
    fm_percent = case_when(
      id == "551" & visit == "v0" ~ fm_percent / 100,
      id == "572" & visit == "v4" ~ fm_percent / 100,
      id == "120" & visit == "v5" ~ fm_percent / 100,
      TRUE ~ fm_percent
    ),
    ffm_percent = case_when(
      id == "590" & visit == "v0" ~ ffm_percent * 10,
      id == "95"  & visit == "v5" ~ ffm_percent / 100,
      id == "130" & visit == "v5" ~ ffm_percent / 100,
      id == "496" & visit == "v5" ~ ffm_percent / 100,
      id == "364" & visit == "v5" ~ ffm_percent / 10,
      id == "217" & visit == "v5" ~ 67.0,
      id == "257" & visit == "v0" ~ 54.4,
      TRUE ~ ffm_percent
    ),
    # BIA quality & consistency checks
    bia_perc_diff = abs((ffm_percent + fm_percent) - 100),
    bia_kg_diff = abs((ffm_kg + fm_kg) - weight_kg),
    bia_resistance_valid = case_when(
    is.na(bia_resistance_50khz) ~ NA,
    between(bia_resistance_50khz, 110, 1000) ~ TRUE, # Device's resistance range (Maltron BioScan 920)
    TRUE ~ FALSE
  )
)

# Define cohort: valid baseline BIA & available baseline shotgun data & not taking antibiotics
# 1. Valid baseline BIA
bia_v0_ids <- baria_muscle_vars |>
  filter(
    visit == "v0",
    !is.na(ffm_kg),
    !is.na(fm_kg),
    bia_perc_diff <= 5,
    bia_kg_diff <= 5
  ) |>
  pull(id) |> 
  unique()

# 2. Available baseline microbiome/shotgun data
mb_v0_ids <- sample_data(baria_mb) |> 
  data.frame() |> 
  filter(Time_Point == "V-1") |> 
  pull(Subject_ID) |> 
  as.character() |>
  unique()

bia_mb_ids <- intersect(bia_v0_ids, mb_v0_ids)

#### Stopped here #####
baria_muscle_final_cohort <- baria_muscle_vars |> 
  filter(id %in% bia_mb_ids) |> 
  mutate(
    ffmi = ffm_kg / ((height_cm / 100)^2),
    smm_kg = ((height_cm^2) / bia_resistance_50khz * 0.401) + (age * -0.071) + 5.102 + if_else(sex == "male", 3.825, 0),
    smm_by_weight = smm_kg / weight_kg
  ) |> 
  group_by(sex, visit) |> # calculate sex-specific cut-offs
  mutate(
    # Calculate tertiles for FFMI and SMM (raw and indexed)
    ffmi_tertile = quantile(ffmi, probs = 1/3, na.rm = TRUE),
    smm_kg_tertile = quantile(smm_kg, probs = 1/3, na.rm = TRUE),
    smm_by_weight_tertile = quantile(smm_by_weight, probs = 1/3, na.rm = TRUE),
    low_ffmi = if_else(ffmi <= ffmi_tertile, "yes", "no"),
    low_smm = if_else(smm_kg <= smm_kg_tertile, "yes", "no"),
    low_smm_by_weight = if_else(smm_by_weight <= smm_by_weight_tertile, "yes", "no")
  ) |> 
  ungroup() 
  #pivot_wider(# empty for now) |> 
  mutate(
    # %BW change from baseline to 1, 2 and 5 years
    perc_weight_change_v4 = (weight_kg_v4 - weight_kg_v0) / weight_kg_v0 * 100,
    perc_weight_change_v5 = (weight_kg_v5 - weight_kg_v0) / weight_kg_v0 * 100,
    perc_weight_change_v6 = (weight_kg_v6 - weight_kg_v0) / weight_kg_v0 * 100,

    # %FFM change from baseline to 1, 2 and 5 years
    perc_ffm_change_v4 = (ffm_kg_v4 - ffm_kg_v0) / ffm_kg_v0 * 100,
    perc_ffm_change_v5 = (ffm_kg_v5 - ffm_kg_v0) / ffm_kg_v0 * 100,
    perc_ffm_change_v6 = (ffm_kg_v6 - ffm_kg_v0) / ffm_kg_v0 * 100,

    # ΔFFMI and %FFMI change from baseline to 1, 2 and 5 years
    delta_ffmi_v4 = ffmi_v4 - ffmi_v0, # 1y
    delta_ffmi_v5 = ffmi_v5 - ffmi_v0, # 2y
    delta_ffmi_v6 = ffmi_v6 - ffmi_v0, # 5y

    perc_ffmi_change_v4 = (ffmi_v4 - ffmi_v0) / ffmi_v0 * 100,
    perc_ffmi_change_v5 = (ffmi_v5 - ffmi_v0) / ffmi_v0 * 100,
    perc_ffmi_change_v6 = (ffmi_v6 - ffmi_v0) / ffmi_v0 * 100,

    # %FM change from baseline to 1, 2 and 5 years
    perc_fm_change_v4 = (fm_kg_v4 - fm_kg_v0) / fm_kg_v0 * 100,
    perc_fm_change_v5 = (fm_kg_v5 - fm_kg_v0) / fm_kg_v0 * 100,
    perc_fm_change_v6 = (fm_kg_v6 - fm_kg_v0) / fm_kg_v0 * 100,

    # Calculate ΔASM and %ASM change from baseline to 1, 2 and 5 years
    delta_asm_v4 = asm_kg_v4 - asm_kg_v0, # 1y
    delta_asm_v5 = asm_kg_v5 - asm_kg_v0, # 2y
    delta_asm_v6 = asm_kg_v6 - asm_kg_v0, # 5y
    
    perc_asm_change_v4 = (asm_kg_v4 - asm_kg_v0) / asm_kg_v0 * 100, 
    perc_asm_change_v5 = (asm_kg_v5 - asm_kg_v0) / asm_kg_v0 * 100,
    perc_asm_change_v6 = (asm_kg_v6 - asm_kg_v0) / asm_kg_v0 * 100) |> 
  group_by(sex) |> 
  mutate(
    # Calculate tertiles for %FFMI change
    ffmi_change_v4_tertile = quantile(perc_ffmi_change_v4, probs = 1/3, na.rm = TRUE),
    ffmi_change_v5_tertile = quantile(perc_ffmi_change_v5, probs = 1/3, na.rm = TRUE),
    ffmi_change_v6_tertile = quantile(perc_ffmi_change_v6, probs = 1/3, na.rm = TRUE),

    # Grouping based on %FFMI change tertile
    ffmi_change_v4_group = if_else(perc_ffmi_change_v4 <= ffmi_change_v4_tertile, "high", "low/modest"),
    ffmi_change_v5_group = if_else(perc_ffmi_change_v5 <= ffmi_change_v5_tertile, "high", "low/modest"),
    ffmi_change_v6_group = if_else(perc_ffmi_change_v6 <= ffmi_change_v6_tertile, "high", "low/modest"), 
   
    # Calculate tertiles for %ASM change
    asm_change_v4_tertile = quantile(perc_asm_change_v4, probs = 1/3, na.rm = TRUE),
    asm_change_v5_tertile = quantile(perc_asm_change_v5, probs = 1/3, na.rm = TRUE),
    asm_change_v6_tertile = quantile(perc_asm_change_v6, probs = 1/3, na.rm = TRUE),

    # pool low/modest together
    asm_change_v4_group = if_else(perc_asm_change_v4 <= asm_change_v4_tertile, "high", "low/modest"), 
    asm_change_v5_group = if_else(perc_asm_change_v5 <= asm_change_v5_tertile, "high", "low/modest"),
    asm_change_v6_group = if_else(perc_asm_change_v6 <= asm_change_v6_tertile, "high", "low/modest")
  ) |> 
  ungroup() |> 
  mutate(across(where(is.character) & !id, as.factor))

## Medication
baria_muscle_clinical_with_medication_notypos <- baria_muscle_clinical |>
  mutate(medication_list_v0 = str_to_lower(medication_list_v0)) |> 
  mutate(
    medication_list_v0 = if_else(
      medication_list_v0 %in% c("geen", "geen medicatie"),
      "none",
      medication_list_v0 # replaced "geen" and "geen medicatie" with english "none"
      ), # fix typos I noticed in the medication list
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bamitriptiline\\b", "amitriptyline"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bduloxatine\\b", "duloxetine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bariprazol\\b", "aripiprazol"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bnifidipine\\b", "nifedipine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bevothyroxine\\b", "levothyroxine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bglicliazide\\b", "gliclazide"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bglicazide\\b", "gliclazide"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bamlopdipine\\b", "amlodipine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bzitalopram\\b", "citalopram"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bperinopril\\b", "perindopril"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bvenlaflaxine\\b", "venlafaxine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bglimerpiride\\b", "glimepiride"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bhydrchloorthiazide\\b", "hct"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bhydrochloorthiazide\\b", "hct"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bpnatoprazol\\b", "pantoprazole"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bomeprazol\\b", "omeprazole"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\besomeprazol\\b", "esomeprazole"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bparvastatine\\b", "pravastatine"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bcandarsetan\\b", "candesartan"),
    medication_list_v0 = str_replace_all(medication_list_v0,"\\bezatimibe\\b", "ezetimibe")
  )

# Define medication categories
# Diabetes medication
dm_meds <- list(
  # Metformin
  metformin = c("metformine", "glucophage"),
  
  # Sulfonylureas
  sus = c("gliclazide", "diamicron", "glibenclamide", "daonil", "glimepiride", "amaryl", 
  "glipizide", "minodiab", "tolbutamide"),
  
  # DPP-4 inhibitors
  dpp4is = c("sitagliptine", "januvia", "vildagliptine", "galvus", "saxagliptine", 
  "onglyza", "linagliptine", "trajenta", "alogliptine", "vipidia"),
  
  # GLP-1 receptor agonists
  glp1ras = c("liraglutide", "victoza", "semaglutide", "ozempic", "rybelsus",
  "exenatide", "byetta", "bydureon", "dulaglutide", "trulicity", "lyxumia", "lixisenatide"),
  
  # SGLT2-inhibitors
  sglt2is = c("dapagliflozine", "forxiga", "empagliflozine", "jardiance",
  "canagliflozine", "invokana", "ertugliflozine", "steglatro", "steeglatro"),

  # thiazolidinediones
  tzds = c("pioglitazon", "pioglitazone"),
  
  #insulins (short-, medium- and long-acting)
  insulin = c("insuline", "lantus", "levemir", "novorapid", "apidra", "toujeo", "tresiba", "degludec",
  "humalog", "novomix", "fiasp", "actrapid", "isofaan", "insulatard", "glargine", "aspart")
)

# create regex patterns for each class of diabetes medication
dm_meds_patterns <- lapply(dm_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"), ")\\b")
})
dm_meds_any_pattern <- str_c("\\b(", str_c((unlist(dm_meds)), collapse = "|"), ")\\b")

# Antihypertensive medication
aht_meds <- list(
  # ACE inhibitors
  ace_inhibitors = c("perindopril", "candesartan", "enalapril", "ramipril", "lisinopril", "captopril", "fosinopril",
  "quinapril"),

  # angiotensin receptor blockers
  arbs = c("losartan", "valsartan", "olmesartan", "irbesartan", "telmisartan", "eprosartan"),

  # calcium channel blockers
  ccbs = c("amlodipine", "lercandipine", "lerdip", "barnidipine", "cyress", "nifedipine", "adalat", "verapamil",
  "felodipine", "nicardipine", "isradipine", "diltiazem"),

  # bblockers
  bblockers = c("metoprolol", "metoprololtart", "metoprololsucc", "selokeen", "carvedilol", "labetalol", "bisoprolol",
  "atenolol", "propanonol", "nebivolol", "solatol"),

  # central a2 agonists
  a2_agonists = c("clonidine", "moxonidine", "methyldopa"),     

  # diuretics
  diuretics = c("furosemide", "hct", "hydrochloorthiazide", "spironolacton", "triamteren", "amiloride",
  "dytenzide", "bumetadine", "chloortalidon", "indapamide", "eplerenone"),

  # combination/others
  combi_aht_meds = c("losartan/hydrochloorthiazide", "lodoz", "preterax", "moduretic")
)

aht_meds_patterns <- lapply(aht_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"),")\\b")
})
aht_meds_any_pattern <- str_c("\\b(", str_c(unlist(aht_meds), collapse = "|"),")\\b")

# Lipid lowering
lipidlowering_meds <- list(
  # Statins (HMG-CoA reductase inhibitors)
  statins = c(
  "simvastatine", "atorvastatine", "pravastatine", "rosuvastatine", "fluvastatine", 
  "crestor", "lipitor", "zocor", "pravachol", "lescol"),

  # Cholesterol absorption inhibitors
  ezetimibe = c("ezetimib", "ezetrol"),

  # PCSK9 inhibitors
  pcsk9is = c("alirocumab", "praluent", "evolocumab", "repatha"),

  # siRNA against PCSK9
  inclisiran = c("inclisiran", "leqvio"),

  # ATP citrate lyase inhibitors
  bempedoic_acid = c("bempedoïnezuur", "nexletol", "nustendi"),

  # Bile acid sequestrants
  bile_acid_sequestrants = c("colestyramine", "questran"),

  # Fibrates (PPAR-α agonists)
  fibrates = c("bezalip", "fibraat", "fenofibraat", "lipanthyl", "tricor", "gemfibrozil", "lopid"),

  # Omega-3 fatty acids
  omega3 = c("omega-3-vetzuren", "omacor", "epanova"),

  # Nicotinic acid
  niacin = c("nicotinezuur", "niaspan")
)

# create regex patterns for each class of lipid-lowering medication
lipidlowering_meds_patterns <- lapply(lipidlowering_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"), ")\\b")
})
lipidlowering_meds_any_pattern <- str_c("\\b(", str_c((unlist(lipidlowering_meds)), collapse = "|"), ")\\b")

# Thyroid medication
thyroid_meds <- c("levothyroxine", "thyrax", "euthyrox")
thyroid_meds_pattern <- str_c("\\b(", str_c(thyroid_meds, collapse = "|"),")\\b")

# Psychiatric medication (that may cause weight loss/gain)
psychiatric_meds <- list(
  # SSRIS
  ssris = c("escitalopram", "citalopram", "cipramil", "sertraline", "paroxetine", "fluoxetine", "fluvoxamine"),

  # tricyclic antidepressants
  tcas = c("amitriptyline", "nortriptyline", "clomipramine"),

  # SNRIs
  snris = c("venlafaxine", "efexor", "duloxetine"),

  # NDRIs
  ndris = c("bupropion"),

  # atypical antipsychotics
  antipsychotics = c("aripiprazol", "quetiapine", "olanzapine"),

  # mood stabilizers (anticonvulsants)
  moods = c("lamotrigine", "pregabalin"),

  # ADHD medication/stimulants
  adhd_meds = c("methylfenidaat", "concerta", "elvanse"),

  # sleep medication / hypnotics
  hypnotics = c("zopiclone", "zolpidem")
)

psychiatric_meds_patterns <- lapply(psychiatric_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"),")\\b")
})
psychiatric_meds_any_pattern <- str_c("\\b(", str_c(unlist(psychiatric_meds), collapse = "|"), ")\\b")

# PPIs (effects on gut microbiota)
ppi <- c("omeprazole", "pantoprazole", "esomeprazole")
ppi_pattern <- str_c("\\b(", str_c(ppi, collapse = "|"), ")\\b")

# Antibiotics
antibiotics <- c(
  #Penicillins
  "amoxicilline", "amoxicillin", "augmentin", "co-amoxiclav", "flucloxacilline", "flucloxacillin", "penicilline", "penicillin",

  # Tetracyclines
  "doxycycline", "minocycline",

  # Macrolides
  "azitromycine", "azithromycin", "claritromycine", "clarithromycin", "erytromycine", "erythromycin",

  # Fluoroquinolones
  "ciprofloxacine", "ciprofloxacin", "levofloxacine", "levofloxacin", "moxifloxacine", "moxifloxacin",

  # Lincosamides
  "clindamycine", "clindamycin",

  # Nitroimidazoles
  "metronidazol", "metronidazole",

  # Urinary tract antibiotics
  "nitrofurantoine", "nitrofurantoin", "fosfomycine", "fosfomycin",

  # Sulfonamides
  "co-trimoxazol", "cotrimoxazol", "trimethoprim", "sulfamethoxazol", "sulfamethoxazole"
)
antibiotics_pattern <- str_c("\\b(", str_c(antibiotics, collapse = "|"), ")\\b")

# Use regex patterns to create clean medication columns for all medication categories and subcategories
baria_muscle_clean <- baria_muscle_clinical_with_medication_notypos |> 
  mutate(
    # Diabetes medication classes
    metformin_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$metformin), "yes", "no"),
    sus_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$sus), "yes", "no"),
    dpp4is_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$dpp4is), "yes", "no"),
    glp1ras_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$glp1ras), "yes", "no"),
    sglt2is_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$sglt2is), "yes", "no"),
    tzds_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$tzds), "yes", "no"),
    insulin_v0 = if_else(str_detect(medication_list_v0, dm_meds_patterns$insulin), "yes", "no"),
    
    # Use of any diabetes medication
    dm_meds_v0 = if_else(str_detect(medication_list_v0, dm_meds_any_pattern), "yes", "no"),

    # Antihypertensive medication classes
    ace_inhibitors_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$ace_inhibitors), "yes", "no"),
    arbs_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$arbs), "yes", "no"),
    ccbs_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$ccbs), "yes", "no"),
    bblockers_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$bblockers), "yes", "no"),
    a2_agonists_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$a2_agonists), "yes", "no"),
    diuretics_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$diuretics), "yes", "no"),
    combi_aht_meds_v0 = if_else(str_detect(medication_list_v0, aht_meds_patterns$combi_aht_meds), "yes", "no"),

    # general antihypertensive medication
    aht_meds_v0 = if_else(str_detect(medication_list_v0, aht_meds_any_pattern), "yes", "no"),

    # Statins
    statins_v0 = if_else(str_detect(medication_list_v0, lipidlowering_meds_patterns$statins), "yes", "no"),

    # Use of any lipid-lowering medication
    lipidlowering_meds_v0 = if_else(str_detect(medication_list_v0, lipidlowering_meds_any_pattern), "yes", "no"),

    # thyroid medication
    thyroid_meds_v0 = if_else(str_detect(medication_list_v0, thyroid_meds_pattern), "yes", "no"),

    # psychiatric medication
    ssris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$ssris), "yes", "no"),
    tcas_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$tcas), "yes", "no"),
    snris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$snris), "yes", "no"),
    ndris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$ndris), "yes", "no"),
    antipsychotics_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$antipsychotics), "yes", "no"),
    moods_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$moods), "yes", "no"),
    adhd_meds_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$adhd_meds), "yes", "no"),
    hypnotics_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_patterns$hypnotics), "yes", "no"),

    psychiatric_meds_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_any_pattern), "yes", "no"),

    # PPIs
    ppi_v0 = if_else(str_detect(medication_list_v0, ppi_pattern), "yes", "no"),

    # Antibiotics
    abx_v0 = if_else(str_detect(medication_list_v0, antibiotics_pattern), "yes", "no")
  ) |> 
  # 3. filter for final analysis cohort: antibiotics exclusion
  filter(abx_v0 == "no") |> 
  select(-medication_list_v0) |> 
  arrange(date_v0)


# then save as both RDS and csv files
write.csv(baria_muscle_clean, "data/processed_data/260810_BARIA_muscle_clinical.csv")
saveRDS(baria_muscle_clean, "data/processed_data/260810_BARIA_muscle_clinical.RDS")