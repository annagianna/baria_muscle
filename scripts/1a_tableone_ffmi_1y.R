## Table 1 stratified by 1-year FFMI loss group
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")

# Baseline vars
t1_ffmi_loss_vars_long <- baria_muscle_wide |>
  select(
    id, age_v0, sex, perc_change_ffmi_v4_group,
    matches(
      "^(bmi|wc_cm|fm_kg|ffm_kg|ffmi|prediab_labs|t2d_labs|hba1c_mmolmol|glucose_mmoll_mmt_0|insulin_pmoll_mmt_0|homa_ir|homa_b|total_cholesterol_mmoll|ldl_cholesterol_mmoll|hdl_cholesterol_mmoll|triglycerides_mmoll|lipidlowering_meds|dm_meds)_v(0|4)$"
    )
  ) |>
  pivot_longer(
    cols = -c(id, age_v0, sex, perc_change_ffmi_v4_group),
    names_pattern = "(.*)_(v0|v4)$",
    names_to = c(".value", "visit")
  ) |> 
  filter(!is.na(perc_change_ffmi_v4_group)) |>
  mutate(
    visit = factor(visit, levels = c("v0", "v4"), labels = c("Baseline", "1 year")),
    perc_change_ffmi_v4_group = factor(perc_change_ffmi_v4_group, levels = c("high", "modest/low"), labels = c("High FFMI loss", "Modest/low FFMI loss"))
  ) |> 
  rename(
    `Age (years)` = age_v0,
    `Sex` = sex,
    `BMI (kg/m²)` = bmi,
    `Waist circumference (cm)` = wc_cm,
    `Fat mass (kg)` = fm_kg,
    `Fat-free mass (kg)` = ffm_kg,
    `Fat-free mass index (kg/m²)` = ffmi,
    `Prediabetes` = prediab_labs,
    `T2D` = t2d_labs,
    `HbA1c (mmol/mol)` = hba1c_mmolmol,
    `Fasting glucose (mmol/L)` = glucose_mmoll_mmt_0,
    `Fasting insulin (pmol/L)` = insulin_pmoll_mmt_0,
    `HOMA-IR` = homa_ir,
    `HOMA-B` = homa_b,
    `Total cholesterol (mmol/L)` = total_cholesterol_mmoll,
    `LDL cholesterol (mmol/L)` = ldl_cholesterol_mmoll,
    `HDL cholesterol (mmol/L)` = hdl_cholesterol_mmoll,
    `Triglycerides (mmol/L)` = triglycerides_mmoll,
    `Lipid-lowering medication` = lipidlowering_meds,
    `Antidiabetic medication` = dm_meds,
    `1-year FFMI loss group` = perc_change_ffmi_v4_group
  )

# Table 1 vars (both v0 and v4)
t1_vars <- c(
  "Age (years)", "Sex", "BMI (kg/m²)", "Waist circumference (cm)",
  "Fat mass (kg)", "Fat-free mass (kg)", "Fat-free mass index (kg/m²)",
  "Prediabetes", "T2D", "HbA1c (mmol/mol)", "Fasting glucose (mmol/L)", "Fasting insulin (pmol/L)", "HOMA-IR", "HOMA-B",
  "Total cholesterol (mmol/L)", "LDL cholesterol (mmol/L)", "HDL cholesterol (mmol/L)", "Triglycerides (mmol/L)",
  "Lipid-lowering medication", "Antidiabetic medication"
)

## Baseline ##
# Baseline table stratified by %FFMI change group at 1y
t1_ffmi_loss_v0 <- CreateTableOne(
  vars = t1_vars,
  strata = "1-year FFMI loss group",
  data = t1_ffmi_loss_vars_long |>
    filter(visit == "Baseline"),
  factorVars = c("Sex", "Prediabetes", "T2D", "Lipid-lowering medication", "Antidiabetic medication"),
  test = TRUE
)
# Print baseline table
t1_ffmi_loss_v0_matrix <- print(
  t1_ffmi_loss_v0,
  nonnormal = c("Fasting glucose (mmol/L)", "Fasting insulin (pmol/L)", "HbA1c (mmol/mol)", "HOMA-IR", "HOMA-B", "Triglycerides (mmol/L)"),
  showAllLevels = FALSE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1,
  catDigits = 1,
  pDigits = 3
)

## 1y ##
# 1y table stratified by %FFMI change group at 1y
t1_ffmi_loss_v4 <- CreateTableOne(
  vars = t1_vars,
  strata = "1-year FFMI loss group",
  data = t1_ffmi_loss_vars_long |>
    filter(visit == "1 year"),
  factorVars = c("Sex", "Prediabetes", "T2D", "Lipid-lowering medication", "Antidiabetic medication"),
  test = TRUE
)

# Print 1y table
t1_ffmi_loss_v4_matrix <- print(
  t1_ffmi_loss_v4,
  nonnormal = c("Fasting glucose (mmol/L)", "Fasting insulin (pmol/L)", "HbA1c (mmol/mol)", "HOMA-IR", "HOMA-B", "Triglycerides (mmol/L)"),
  showAllLevels = FALSE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1,
  catDigits = 1,
  pDigits = 3
)

# Save tables
write.csv(t1_ffmi_loss_v0_matrix, "results/tables/t1_ffmi_loss_v0_matrix.csv", row.names = TRUE)
write.csv(t1_ffmi_loss_v4_matrix, "results/tables/t1_ffmi_loss_v4_matrix.csv", row.names = TRUE)