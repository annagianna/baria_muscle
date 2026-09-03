## Table 1 stratified by baseline FFMI status
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")

#### Table 1 grouped by FFMI status at baseline (participants with available shotgun data) ####
# Choose vars
t1_v0_mb_vars <- baria_muscle_wide |> 
  select(
    age_v0, sex, bmi_v0, wc_cm_v0, fm_kg_v0, ffm_kg_v0, ffmi_v0, low_ffmi_v0, prediab_labs_v0, t2d_v0, 
    hba1c_mmolmol_v0, glucose_mmoll_mmt_0_v0, insulin_pmoll_mmt_0_v0, homa_ir_v0, homa_b_v0
  ) |> 
  filter(!is.na(low_ffmi_v0)) |> 
  mutate(
    low_ffmi_v0 = factor(
      low_ffmi_v0,
      levels = c("yes", "no"),
      labels = c("Low FFMI", "Moderate/high FFMI")
    )
  ) |> 
  rename(
    `Age (years)` = age_v0,
    `Sex` = sex,
    `BMI (kg/m²)` = bmi_v0,
    `Waist circumference (cm)` = wc_cm_v0,
    `Fat mass (kg)` = fm_kg_v0,
    `Fat-free mass (kg)` = ffm_kg_v0,
    `Fat-free mass index (kg/m²)` = ffmi_v0,
    `FFMI group` = low_ffmi_v0,
    `Prediabetes` = prediab_labs_v0,
    `T2D` = t2d_v0,
    `HbA1c (mmol/mol)` = hba1c_mmolmol_v0,
    `Fasting glucose (mmol/L)` = glucose_mmoll_mmt_0_v0,
    `Fasting insulin (pmol/L)` = insulin_pmoll_mmt_0_v0,
    `HOMA-IR` = homa_ir_v0,
    `HOMA-B` = homa_b_v0
  )

# Create table
t1_v0_mb <- CreateTableOne(
  vars = names(t1_v0_mb_vars)[names(t1_v0_mb_vars) != "FFMI group"],
  strata = "FFMI group",
  data = t1_v0_mb_vars, 
  factorVars = c("Sex", "Prediabetes", "T2D"),
  test = TRUE
)

# Print & format table
t1_v0_mb_matrix <- print(
  t1_v0_mb,
  nonnormal = c(
    "Fasting glucose (mmol/L)", 
    "Fasting insulin (pmol/L)",
    "HbA1c (mmol/mol)",
    "HOMA-IR",
    "HOMA-B"
  ),
  showAllLevels = FALSE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1,
  catDigits = 1,
  pDigits = 3
  )

# Write table
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(t1_v0_mb_matrix, file = "results/tables/t1_v0_mb_matrix.csv", row.names = TRUE)
