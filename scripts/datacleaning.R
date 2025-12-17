# Data cleaning for BARIA Project on muscle mass & weight trajectories
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(lubridate)
library(ggpubr)
library(purrr)

# Open Data and see properties

baria_clinical_data_raw <- readRDS("./data/BARIA.clinical.2024-12-09.723.2043.RDS")
View(baria_clinical_data_raw)

# roughly check for incorrect data entries of BIA data (repeat for v0, v4-v6)
baria_clinical_data_raw |> 
  select(Subject_ID, tbf, ffm, weight) |> 
  filter(!is.na(tbf)) |> 
  mutate(
    sum = tbf + ffm,
    nomatch = if_else((abs((tbf + ffm) - weight) > 1), "nomatch", "match")
  ) |> 
  filter(nomatch == "nomatch") |> 
  View()

# clinical data
# Need: clinical data, anthropometry, body composition, medication, basic lab and cardiometabolic risk factors, diabetes, medication
baria_muscle_clinical <- baria_clinical_data_raw |>
  select(id = Subject_ID,

    # visits
    date_v0 = date,
    date_v2 = V2_date,
    date_v3 = V3_date,
    date_v4 = V4_date,
    date_v5 = V5_date,
    date_v6 = V6_date_1,
    date_v7 = V7_date,

    # Baseline (v0)
    age_v0 = Age,
    sex,
    contains("race"),
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
    across(
      everything(),
      ~ replace(.x, .x %in% c(-99, -98, -99), NA)), # these values are NA for different reasons in the Baria dataset
    across(
     .cols = starts_with("date"),
      ~ na_if(., "01-01-2999")),
    across(
      .cols = starts_with("date"),
      ~ na_if(., "01-01-2997"))
    ) |> 
  mutate( # time diff from baseline visit in weeks
    across(date_v0:date_v7, dmy),
    v0_to_v2_weeks = difftime(date_v2, date_v0, units = "weeks"),
    v0_to_v3_weeks = difftime(date_v3, date_v0, units = "weeks"),
    v0_to_v4_weeks = difftime(date_v4, date_v0, units = "weeks"),
    v0_to_v5_weeks = difftime(date_v5, date_v0, units = "weeks"),
    v0_to_v6_weeks = difftime(date_v6, date_v0, units = "weeks"),
    v0_to_v7_weeks = difftime(date_v7, date_v0, units = "weeks")
  ) |> 
  mutate( # calculate age at each visit to be used in the SMM and ASM formulas
    # first convert difftime in weeks to numeric in years
    v0_to_v2_numeric_y = as.numeric(v0_to_v2_weeks)/(365.25 / 7),
    v0_to_v3_numeric_y = as.numeric(v0_to_v3_weeks)/(365.25 / 7),
    v0_to_v4_numeric_y = as.numeric(v0_to_v4_weeks)/(365.25 / 7),
    v0_to_v5_numeric_y = as.numeric(v0_to_v5_weeks)/(365.25 / 7),
    v0_to_v6_numeric_y = as.numeric(v0_to_v6_weeks)/(365.25 / 7),
    v0_to_v7_numeric_y = as.numeric(v0_to_v7_weeks)/(365.25 / 7)
  ) |>
  filter(!is.na(age_v0)) |> 
  mutate(# calculate the age at every visit
    age_v2 = age_v0 + v0_to_v2_numeric_y,
    age_v3 = age_v0 + v0_to_v3_numeric_y,
    age_v4 = age_v0 + v0_to_v4_numeric_y,
    age_v5 = age_v0 + v0_to_v5_numeric_y,
    age_v6 = age_v0 + v0_to_v6_numeric_y,
    age_v7 = age_v0 + v0_to_v7_numeric_y
  ) |> 
  mutate(
    sex = case_when(sex == "1" ~ "male", sex == "2" ~ "female"),
    t2d_v0 = case_when(t2d_v0 == "1" ~ "yes", t2d_v0 == "2" ~ "no"),
    aht = case_when(aht == "1" ~ "yes", aht == "2" ~ "no"),
    medication_binary_v0 = case_when(medication_binary_v0 == "1" ~ "yes", medication_binary_v0 == "2" ~ "no"),
    sport_v0 = case_when(sport_v0 == "1" ~ "yes", sport_v0 == "2" ~ "no")
  ) |> 
  mutate(
    # fix Hba1c units
    hba1c_percent_v0 = if_else(
      hba1c_v0 < 15,
      hba1c_v0,
      hba1c_v0 * 0.0915 + 2.15
    ),
    hba1c_percent_v2 = if_else(
      hba1c_v2 < 15,
      hba1c_v2,
      hba1c_v2 * 0.0915 + 2.15
    ),
    hba1c_percent_v3 = if_else(
      hba1c_v3 < 15,
      hba1c_v3,
      hba1c_v3 * 0.0915 + 2.15
    ),
    hba1c_percent_v4 = if_else(
      hba1c_v4 < 15,
      hba1c_v4,
      hba1c_v4 * 0.0915 + 2.15
    ),
    hba1c_percent_v5 = if_else(
      hba1c_v5 < 15,
      hba1c_v5,
      hba1c_v5 * 0.0915 + 2.15
    ),
    hba1c_percent_v6 = if_else(
      hba1c_v6 < 15,
      hba1c_v6,
      hba1c_v6 * 0.0915 + 2.15
    ),
    hba1c_percent_v7 = if_else(
      hba1c_v7 < 15,
      hba1c_v7,
      hba1c_v7 * 0.0915 + 2.15
    ),
    # hba1c in mmol/l
    hba1c_mmolmol_v0 = if_else(
      is.na(hba1c_mmolmol_v0) == FALSE,
      hba1c_mmolmol_v0,
      10.93 * hba1c_percent_v0 - 23,5
    ),
    hba1c_mmolmol_v2 = if_else(
      is.na(hba1c_mmolmol_v2) == FALSE,
      hba1c_mmolmol_v2,
      10.93 * hba1c_percent_v2 - 23,5
    ),
    hba1c_mmolmol_v3 = if_else(
      is.na(hba1c_mmolmol_v3) == FALSE,
      hba1c_mmolmol_v3,
      10.93 * hba1c_percent_v3 - 23,5
    ),
    hba1c_mmolmol_v4 = if_else(
      is.na(hba1c_mmolmol_v4) == FALSE,
      hba1c_mmolmol_v4,
      10.93 * hba1c_percent_v4 - 23,5
    ),
    hba1c_mmolmol_v5 = if_else(
      is.na(hba1c_mmolmol_v5) == FALSE,
      hba1c_mmolmol_v5,
      10.93 * hba1c_percent_v5 - 23,5
    ),
    hba1c_mmolmol_v6 = if_else(
      is.na(hba1c_mmolmol_v6) == FALSE,
      hba1c_mmolmol_v6,
      10.93 * hba1c_percent_v6 - 23,5
    ),
    hba1c_mmolmol_v7 = if_else(
      is.na(hba1c_mmolmol_v7) == FALSE,
      hba1c_mmolmol_v7,
      10.93 * hba1c_percent_v7 - 23,5
    ),

    # HOMA-IR & HOMA-2B (only v0, v4-v6)
    # insulin unit conversion from pmol/l to uU/ml
    homa_ir_v0 = (fasting_insulin_pmoll_mmt_v0 / 6.945) * fasting_glucose_mmoll_mmt_v0 / 22.5,
    homa_b_v0 = (20 * (fasting_insulin_pmoll_mmt_v0 / 6.945)) / (fasting_glucose_mmoll_mmt_v0 - 3.5),
    homa_ir_v4 = (fasting_insulin_pmoll_mmt_v4 / 6.945) * fasting_glucose_mmoll_mmt_v4 / 22.5,
    homa_b_v4 = (20 * (fasting_insulin_pmoll_mmt_v4 / 6.945)) / (fasting_glucose_mmoll_mmt_v4 - 3.5),
    homa_ir_v5 = (fasting_insulin_pmoll_mmt_v5 / 6.945) * fasting_glucose_mmoll_mmt_v5 / 22.5,
    homa_b_v5 = (20 * (fasting_insulin_pmoll_mmt_v5 / 6.945)) / (fasting_glucose_mmoll_mmt_v5 - 3.5)
  ) |> 
  # body composition
  # fix incorrect data entries for bia data
  mutate(
    id = as.character(id),
    fm_kg_v0 = if_else(id == "212", fm_kg_v0/100, fm_kg_v0),
    fm_kg_v0 = if_else(id == "600", fm_kg_v0*10, fm_kg_v0),
    ffm_kg_v0 = if_else(id == "238", ffm_kg_v0*10, ffm_kg_v0),
    fm_kg_v4 = if_else(id == "30004", fm_kg_v4/100, fm_kg_v4),
    fm_kg_v4 = if_else(id == "335", fm_kg_v4*10, fm_kg_v4),
    ffm_kg_v5 = if_else(id == "50", 79.16, ffm_kg_v5)
  ) |> 
  filter( # filtering for device's (Maltron BioScan 920) resolution range excluding high and low outliers
    # qc baseline BIA results
    bia_resistance_50khz_v0 >= 110,
    bia_resistance_50khz_v0 <= 1000,
    !is.na(ffm_percent_v0), !is.na(fm_percent_v0), # baseline bia is needed

    # check if percentages add up
    # baseline
    abs(ffm_percent_v0 + fm_percent_v0 - 100) < 0.1,
    # follow-up: keep if follow-up bia is missing
    # if follow-up bia available, both ffm and fm must be available and add up to 100
    (is.na(ffm_percent_v4) & is.na(fm_percent_v4)) |
    (!is.na(ffm_percent_v4) & !is.na(fm_percent_v4) & abs(ffm_percent_v4 + fm_percent_v4 - 100) < 0.1),
    (is.na(ffm_percent_v5) & is.na(fm_percent_v5)) |
    (!is.na(ffm_percent_v5) & !is.na(fm_percent_v5) & abs(ffm_percent_v5 + fm_percent_v5 - 100) < 0.1),
    (is.na(ffm_percent_v6) & is.na(fm_percent_v6)) |
    (!is.na(ffm_percent_v6) & !is.na(fm_percent_v6) & abs(ffm_percent_v6 + fm_percent_v6 - 100) < 0.1),

    # qc for v4, v5, v6 if not NA (50 khz is used in smm calculation later on)
    is.na(bia_resistance_50khz_v4) | (bia_resistance_50khz_v4 >= 110 & bia_resistance_50khz_v4 <= 1000),
    is.na(bia_resistance_50khz_v5)| (bia_resistance_50khz_v5 >= 110 & bia_resistance_50khz_v5 <= 1000),
    is.na(bia_resistance_50khz_v6) | (bia_resistance_50khz_v6 >= 110 & bia_resistance_50khz_v6 <= 1000)
  )|> 
  mutate(
    # ffm and fm indices v0, v4, v5, v6
    ffmi_v0 = ffm_kg_v0 / ((height_cm / 100)^2),
    fmi_v0 = fm_kg_v0 / ((height_cm / 100)^2),
    ffmi_v4 = ffm_kg_v4 / ((height_cm / 100)^2),
    fmi_v4 = fm_kg_v4 / ((height_cm / 100)^2),
    ffmi_v5 = ffm_kg_v5 / ((height_cm / 100)^2),
    fmi_v5 = fm_kg_v5 / ((height_cm / 100)^2),
    ffmi_v6 = ffm_kg_v6 / ((height_cm / 100)^2),
    fmi_v6 = fm_kg_v6 / ((height_cm / 100)^2),

    # calculate muscle mass indices
    # SMM by Janssen: SMM [kg] = (height^2 [cm] / BIA-resistance [Ohms] X 0.401) + (gender x 3.825) + (age [years] x - 0.071)] + 5.102 (men = 1; women = 0)
    smm_kg_v0 = if_else(
      sex == 1, # male
      0.401 * ((height_cm^2) / bia_resistance_50khz_v0) + (1 * 3.825) + (age_v0 * -0.071) + 5.102,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v0) + (age_v0 * -0.071) + 5.102
    ),
    smm_kg_v4 = if_else(
      sex == 1,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v4) + (1 * 3.825) + (age_v4 * -0.071) + 5.102,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v4) + (age_v4 * -0.071) + 5.102
    ),
    smm_kg_v5 = if_else(
      sex == 1,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v5) + (1 * 3.825) + (age_v5 * -0.071) + 5.102,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v5) + (age_v5 * -0.071) + 5.102
    ),
    smm_kg_v6 = if_else(
      sex == 1,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v6) + (1 * 3.825) + (age_v6 * -0.071) + 5.102,
      0.401 * ((height_cm^2) / bia_resistance_50khz_v6) + (age_v6 * -0.071) + 5.102
    )
  ) |> 
  mutate(
    # (sex: women = 0, men = 1; race: White or Hispanic = 0, Black = 1.9, Asian = −1.6) 
    # race in a format to be used for Janssen formula (White/Hispanic vs. Black vs. Asian)
    race = case_when(
      `race#Kaukasisch` == "1" ~ "white/hispanic",
      `race#Mediterraans` == "1" ~ "white/hispanic",
      `race#Midden-Aziatisch` == "1" ~ "asian",
      `race#Negroïde` == "1" ~ "black",
      `race#Noord-Afrikaans` == "1" ~ "white/hispanic", # for the purposes of the jannsen formula
      `race#Oost-Aziatisch` == "1" ~ "asian",
      `race#Overig` == "1" ~ "white/hispanic", # to be able to use the formula
      `race#Slavisch` == "1" ~ "white/hispanic",
      `race#West-Aziatisch` == "1" ~ "asian",
      `race#Zuid-Amerikaans` == "1" ~ "white/hispanic"
    ),
    # Janssen formula: ASM = (0.244 × weight [kg]) + (7.8 × height [m]) + (6.6 × sex) – (0.098 × age [years]) + (race – 3.3); 
    # (sex: women = 0, men = 1; race: White or Hispanic = 0, Black = 1.9, Asian = −1.6)
    race_num = as.numeric(case_when(race == "white/hispanic" ~ "0", race == "black" ~ "1.9", race == "asian" ~ "-1.6"))
  ) |> 
  select(-contains("race#")) |>
  mutate(
     asm_kg_v0 = if_else(
      sex == 1, # male
      (0.244 * weight_kg_v0) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v0) + (race_num - 3.3), # height must be in m not cm
      (0.244 * weight_kg_v0) + (7.8 * (height_cm/100)) - (0.098 * age_v0) + (race_num - 3.3)
     ),
     asm_kg_v2 = if_else(
      sex == 1,
      (0.244 * weight_kg_v2) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v2) + (race_num - 3.3),
      (0.244 * weight_kg_v2) + (7.8 * (height_cm/100)) - (0.098 * age_v2) + (race_num - 3.3)
     ),
     asm_kg_v3 = if_else(
      sex == 1,
      (0.244 * weight_kg_v3) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v3) + (race_num - 3.3),
      (0.244 * weight_kg_v3) + (7.8 * (height_cm/100)) - (0.098 * age_v3) + (race_num - 3.3)
     ),
     asm_kg_v4 = if_else(
      sex == 1,
      (0.244 * weight_kg_v4) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v4) + (race_num - 3.3),
      (0.244 * weight_kg_v4) + (7.8 * (height_cm/100)) - (0.098 * age_v4) + (race_num - 3.3)
     ),
     asm_kg_v5 = if_else(
      sex == 1,
      (0.244 * weight_kg_v5) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v5) + (race_num - 3.3),
      (0.244 * weight_kg_v5) + (7.8 * (height_cm/100)) - (0.098 * age_v5) + (race_num - 3.3)
     ),
     asm_kg_v6 = if_else(
      sex == 1,
      (0.244 * weight_kg_v6) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v6) + (race_num - 3.3),
      (0.244 * weight_kg_v6) + (7.8 * (height_cm/100)) - (0.098 * age_v6) + (race_num - 3.3)
     ),
     asm_kg_v7 = if_else(
      sex == 1,
      (0.244 * weight_kg_v7) + (7.8 * (height_cm/100)) + 6.6 - (0.098 * age_v7) + (race_num - 3.3),
      (0.244 * weight_kg_v7) + (7.8 * (height_cm/100)) - (0.098 * age_v7) + (race_num - 3.3)
    ) ,
    # smm/weight (recommended for BIA by ESPEN/EASO consensus statement)
    smm_by_weight_v0 = smm_kg_v0 / weight_kg_v0,
    smm_by_weight_v4 = smm_kg_v4 / weight_kg_v4,
    smm_by_weight_v5 = smm_kg_v5 / weight_kg_v5,
    smm_by_weight_v6 = smm_kg_v6 / weight_kg_v6,
  ) |> 
  group_by(sex) |> # calculate cut-offs for female and male participants
  mutate(
    # Calculate the quintiles
    ffm_kg_v0_quintile = quantile(ffm_kg_v0, probs = 0.20, na.rm = TRUE), 
    smm_kg_v0_quintile = quantile(smm_kg_v0, probs = 0.20, na.rm = TRUE), 
    asm_kg_v0_quintile = quantile(asm_kg_v0, probs = 0.20, na.rm = TRUE),
    smm_by_weight_v0_quintile = quantile(smm_by_weight_v0, probs = 0.20, na.rm = TRUE),

    low_ffm_v0 = if_else(ffm_kg_v0 <= ffm_kg_v0_quintile, 1, 0),
    low_smm_v0 = if_else(smm_kg_v0 <= smm_kg_v0_quintile, 1, 0),
    low_asm_v0 = if_else(asm_kg_v0 <= asm_kg_v0_quintile, 1, 0),
    low_smm_by_weight_v0 = if_else(smm_by_weight_v0 <= smm_by_weight_v0_quintile, 1, 0),

    low_ffm_v0 = case_when(low_ffm_v0 == 1 ~ "yes", low_ffm_v0 == 0 ~ "no"),
    low_smm_v0 = case_when(low_smm_v0 == 1 ~ "yes", low_smm_v0 == 0 ~ "no"),
    low_asm_v0 = case_when(low_asm_v0 == 1 ~ "yes", low_asm_v0 == 0 ~ "no"), 
    low_smm_by_weight_v0 = case_when(low_smm_by_weight_v0 == 1 ~ "yes", low_smm_by_weight_v0 == 0 ~ "no"),

    # 2 year follow-up (v5) only asm, smm and smm/weight
    smm_kg_v5_quintile = quantile(smm_kg_v5, probs = 0.20, na.rm = TRUE), 
    asm_kg_v5_quintile = quantile(asm_kg_v5, probs = 0.20, na.rm = TRUE),
    smm_by_weight_v5_quintile = quantile(smm_by_weight_v5, probs = 0.20, na.rm = TRUE),

    low_smm_v5 = if_else(smm_kg_v5 <= smm_kg_v5_quintile, 1, 0),
    low_asm_v5 = if_else(asm_kg_v5 <= asm_kg_v5_quintile, 1, 0),
    low_smm_by_weight_v5 = if_else(smm_by_weight_v5 <= smm_by_weight_v5_quintile, 1, 0),

    low_smm_v5 = case_when(low_smm_v5 == 1 ~ "yes", low_smm_v5 == 0 ~ "no"),
    low_asm_v5 = case_when(low_asm_v5 == 1 ~ "yes", low_asm_v5 == 0 ~ "no"), 
    low_smm_by_weight_v5 = case_when(low_smm_by_weight_v5 == 1 ~ "yes", low_smm_by_weight_v5 == 0 ~ "no"),

    # 5 year follow-up (v6) only asm, smm and smm/weight
    smm_kg_v6_quintile = quantile(smm_kg_v6, probs = 0.20, na.rm = TRUE), 
    asm_kg_v6_quintile = quantile(asm_kg_v6, probs = 0.20, na.rm = TRUE),
    smm_by_weight_v6_quintile = quantile(smm_by_weight_v6, probs = 0.20, na.rm = TRUE),

    low_smm_v6 = if_else(smm_kg_v6 <= smm_kg_v6_quintile, 1, 0),
    low_asm_v6 = if_else(asm_kg_v6 <= asm_kg_v6_quintile, 1, 0),
    low_smm_by_weight_v6 = if_else(smm_by_weight_v6 <= smm_by_weight_v6_quintile, 1, 0),

    low_smm_v6 = case_when(low_smm_v6 == 1 ~ "yes", low_smm_v6 == 0 ~ "no"),
    low_asm_v6 = case_when(low_asm_v6 == 1 ~ "yes", low_asm_v6 == 0 ~ "no"), 
    low_smm_by_weight_v6 = case_when(low_smm_by_weight_v6 == 1 ~ "yes", low_smm_by_weight_v6 == 0 ~ "no")
  ) |> 
  ungroup() |> 
  relocate(race, .after = sex) |> 
  relocate(c(hba1c_percent_v0, hba1c_mmolmol_v0, homa_ir_v0, homa_b_v0), .after = hba1c_mmolmol_v0) |> # Hba1c 
  relocate(hba1c_percent_v2, hba1c_mmolmol_v2, .after = hba1c_mmolmol_v2) |>
  relocate(hba1c_percent_v3, hba1c_mmolmol_v3, .after = hba1c_mmolmol_v3) |>
  relocate(c(hba1c_percent_v4, hba1c_mmolmol_v4, homa_ir_v4, homa_b_v4), .after = hba1c_mmolmol_v4) |>
  relocate(c(hba1c_percent_v5, hba1c_mmolmol_v5, homa_ir_v5, homa_b_v5), .after = hba1c_mmolmol_v5) |>
  relocate(hba1c_percent_v6, hba1c_mmolmol_v6, .after = hba1c_mmolmol_v6) |>
  relocate(hba1c_percent_v7, hba1c_mmolmol_v7, .after = hba1c_mmolmol_v7) |>
  relocate(c(ffmi_v0, fmi_v0), .after = ffm_percent_v0) |> # FFMI & FMI
  relocate(c(ffmi_v4, fmi_v4), .after = ffm_percent_v4) |>
  relocate(c(ffmi_v5, fmi_v5), .after = ffm_percent_v5) |>
  relocate(c(ffmi_v6, fmi_v6), .after = ffm_percent_v6) |>
  relocate(
    c(v0_to_v2_weeks, v0_to_v3_weeks, v0_to_v4_weeks, v0_to_v5_weeks, v0_to_v6_weeks, v0_to_v7_weeks),
    .after = date_v7
  ) |>
  relocate(age_v2:age_v7, .after = age_v0) |>
  relocate(v0_to_v2_numeric_y:v0_to_v7_numeric_y, .after = v0_to_v7_weeks) |> 
  relocate(c(smm_kg_v0, asm_kg_v0, smm_by_weight_v0), .after = upperleg_cm_v0) |> # BIA
  relocate(asm_kg_v2, .after = upperleg_cm_v2) |> 
  relocate(asm_kg_v3, .after = upperleg_cm_v3) |> 
  relocate(c(smm_kg_v4, asm_kg_v4, smm_by_weight_v4), .after = upperleg_cm_v4) |> # BIA
  relocate(c(smm_kg_v5, asm_kg_v5, smm_by_weight_v5), .after = upperleg_cm_v5) |> # BIA
  relocate(c(smm_kg_v6, asm_kg_v6, smm_by_weight_v6), .after = upperleg_cm_v6) |> # BIA
  relocate(asm_kg_v7, .after = upperleg_cm_v7) |>
  relocate(ffm_kg_v0_quintile:smm_by_weight_v0_quintile, .after = smm_kg_v0) |> 
  relocate(low_ffm_v0:low_smm_by_weight_v0, .after = smm_by_weight_v0_quintile) |> 
  mutate(across(where(is.character), as.factor)) |>
  print()

nrow(baria_muscle_clinical)
View(baria_muscle_clinical)

# to explore the range of hba1c in the population (solve mixed variable containing different units)
ggplot(baria_clinical_data_raw, aes(hba1c)) +
  geom_histogram(binwidth = 1.2) # no values between ~15-25

## Formulas used:
# Hba1c(%) = (0,0915 * HbA1c (mmol/mol) + 2,15 (from diabetesfonds.nl)
# Hba1c(mmol/mol) = (10,93 x Hba1c (%)) - 23,5

# DOI: 10.2337/diacare.21.12.2191 for HOMA calculations
# HOMA-IR = (Fasting insulin (µU/mL) × Fasting glucose (mmol/L)) / 22.5 (for glucose in mmol/L)
# HOMA-B = (20 × Fasting insulin (µU/mL)) / (Fasting glucose (mmol/L) − 3.5)
# Insulin was converted from pmol/l to μU/ml

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
dm_meds_patterns

# Antihypertensive medication

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
dm_meds_pattern <- lapply(dm_meds, function(x) {
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

aht_meds_pattern <- lapply(aht_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"),")\\b")
})
aht_meds_any_pattern <- str_c("\\b(", str_c(unlist(aht_meds), collapse = "|"),")\\b")

# Lipid lowering

lipidlowering_meds <- c("simvastatine", "atorvastatine", "pravastatine", "rosuvastatine", "bezalip", "fibraat",
                        "crestor", "lipitor", "zocor", "pravachol", "fluvastatine", "lescol", "ezetimib", "ezetrol",
                        "alirocumab", "praluent", "evolocumab", "repatha", "inclisiran", "leqvio", "bempedoïnezuur", "nexletol", 
                        "nustendi", "colestyramine", "questran","fenofibraat", "lipanthyl", "tricor",
                        "gemfibrozil", "lopid","omega-3-vetzuren", "omacor", "epanova","nicotinezuur", "niaspan" 
)

str_c(lipidlowering_meds, collapse = "|") # inner
lipidlowering_meds_pattern <- str_c("\\b(", str_c(lipidlowering_meds, collapse = "|"),")\\b")

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

psychiatric_meds_pattern <- lapply(psychiatric_meds, function(x) {
  str_c("\\b(", str_c(x, collapse = "|"),")\\b")
})
psychiatric_meds_any_pattern <- str_c("\\b(", str_c(unlist(psychiatric_meds), collapse = "|"), ")\\b")

# PPIs (effects on gut microbiota)
ppi <- c("omeprazole", "pantoprazole", "esomeprazole")
ppi_pattern <- str_c("\\b(", str_c(ppi, collapse = "|"), ")\\b")

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
    ace_inhibitors_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$ace_inhibitors), "yes", "no"),
    arbs_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$arbs), "yes", "no"),
    ccbs_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$ccbs), "yes", "no"),
    bblockers_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$bblockers), "yes", "no"),
    a2_agonists_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$a2_agonists), "yes", "no"),
    diuretics_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$diuretics), "yes", "no"),
    combi_aht_meds_v0 = if_else(str_detect(medication_list_v0, aht_meds_pattern$combi_aht_meds), "yes", "no"),

    # general antihypertensive medication
    aht_meds_v0 = if_else(str_detect(medication_list_v0, aht_meds_any_pattern), "yes", "no"),

    # lipid lowering
    lipidlowering_meds_v0 = if_else(str_detect(medication_list_v0, lipidlowering_meds_pattern), "yes", "no"),

    # thyroid medication
    thyroid_meds_v0 = if_else(str_detect(medication_list_v0, thyroid_meds_pattern), "yes", "no"),

    # psychiartric medication
    ssris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$ssris), "yes", "no"),
    tcas_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$tcas), "yes", "no"),
    snris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$snris), "yes", "no"),
    ndris_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$ndris), "yes", "no"),
    antipsychotics_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$antipsychotics), "yes", "no"),
    moods_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$moods), "yes", "no"),
    adhd_meds = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$adhd_meds), "yes", "no"),
    hypnotics_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_pattern$hypnotics), "yes", "no"),

    psychiatric_meds_v0 = if_else(str_detect(medication_list_v0, psychiatric_meds_any_pattern), "yes", "no"),

    # PPIs
    ppi_v0 = if_else(str_detect(medication_list_v0, ppi_pattern), "yes", "no")
  ) |> 
  select(-medication_list_v0) |> 
  arrange(date_v0) |> 
  print()

View(baria_muscle_clean)
nrow(baria_muscle_clean)

# then save as both RDS and csv files
write.csv(baria_muscle_clean, "data/251217_BARIA_muscle_clinical.csv")
saveRDS(baria_muscle_clean, "data/251217_BARIA_muscle_clinical.RDS")
