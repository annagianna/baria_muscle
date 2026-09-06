## Table 1: Baseline characteristics stratified by 1-year FFMI loss group
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")

# Choose baseline variables
t1_ffmi_loss_vars_v0 <- baria_muscle_wide |>
  select(
    age_v0, sex, bmi_v0,
    wc_cm_v0, fm_kg_v0, ffm_kg_v0,ffmi_v0,perc_change_ffmi_v4_group,
    prediab_labs_v0, t2d_v0, hba1c_mmolmol_v0, glucose_mmoll_mmt_0_v0, insulin_pmoll_mmt_0_v0, homa_ir_v0, homa_b_v0,
    total_cholesterol_mmoll_v0, ldl_cholesterol_mmoll_v0, hdl_cholesterol_mmoll_v0, triglycerides_mmoll_v0, lipidlowering_meds_v0, dm_meds_v0
  ) |>
  filter(!is.na(perc_change_ffmi_v4_group)) |>
  mutate(perc_change_ffmi_v4_group = factor(perc_change_ffmi_v4_group, levels = c("high", "modest/low"), labels = c("High FFMI loss", "Modest/low FFMI loss"))) |>
  rename(
    `Age (years)` = age_v0,
    `Sex` = sex,
    `BMI (kg/m²)` = bmi_v0,
    `Waist circumference (cm)` = wc_cm_v0,
    `Fat mass (kg)` = fm_kg_v0,
    `Fat-free mass (kg)` = ffm_kg_v0,
    `Fat-free mass index (kg/m²)` = ffmi_v0,
    `1-year FFMI loss group` = perc_change_ffmi_v4_group,
    `Prediabetes` = prediab_labs_v0,
    `T2D` = t2d_v0,
    `HbA1c (mmol/mol)` = hba1c_mmolmol_v0,
    `Fasting glucose (mmol/L)` = glucose_mmoll_mmt_0_v0,
    `Fasting insulin (pmol/L)` = insulin_pmoll_mmt_0_v0,
    `HOMA-IR` = homa_ir_v0,
    `HOMA-B` = homa_b_v0,
    `Total cholesterol (mmol/L)` = total_cholesterol_mmoll_v0,
    `LDL cholesterol (mmol/L)` = ldl_cholesterol_mmoll_v0,
    `HDL cholesterol (mmol/L)` = hdl_cholesterol_mmoll_v0,
    `Triglycerides (mmol/L)` = triglycerides_mmoll_v0,
    `Lipid-lowering medication` = lipidlowering_meds_v0,
    `Antidiabetic medication` = dm_meds_v0
  )

# Create Table 1
t1_ffmi_loss_1y_v0 <- CreateTableOne(
  vars = names(t1_ffmi_loss_vars_v0)[
    names(t1_ffmi_loss_vars_v0) != "1-year FFMI loss group"
  ],
  strata = "1-year FFMI loss group",
  data = t1_ffmi_loss_vars_v0,
  factorVars = c("Sex", "Prediabetes", "T2D", "Lipid-lowering medication", "Antidiabetic medication"),
  test = TRUE
)

# Print and format table
t1_ffmi_loss_1y_v0_matrix <- print(
  t1_ffmi_loss_1y_v0,
  nonnormal = c("Fasting glucose (mmol/L)", "Fasting insulin (pmol/L)", "HbA1c (mmol/mol)", "HOMA-IR", "HOMA-B", "Triglycerides (mmol/L)"),
  showAllLevels = FALSE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1,
  catDigits = 1,
  pDigits = 3
)

# Inspect
t1_ffmi_loss_1y_v0_matrix

# Write table
write.csv(t1_ffmi_loss_1y_v0_matrix, file = "results/tables/t1_ffmi_loss_1y_v0_matrix.csv", row.names = TRUE)
