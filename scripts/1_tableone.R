## Table 1 for BARIA muscle project
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(phyloseq)
library(tableone)

# Data
baria_muscle_ab <- readRDS("data/20260722_BARIA_muscle_clinical.RDS")
baria_mb <- readRDS("data/ps.BARIA.metaphlan.706.2548.RDS")

# Remove participants taking antibiotics
baria_muscle <- baria_muscle_ab |> 
  filter(abx_v0 == "no")

#### Table 1 grouped by FFMI status at baseline (clinical cohort) ####
# Choose vars
t1_v0_vars <- baria_muscle |> 
  select(age_v0, sex, bmi_v0, wc_cm_v0, fm_kg_v0, ffm_kg_v0, ffmi_v0, low_ffmi_v0, prediab_v0, t2d_v0, 
    hba1c_percent_v0, fasting_glucose_mmoll_mmt_v0, fasting_insulin_pmoll_mmt_v0, homa_ir_v0, homa_b_v0) |> 
  filter(!is.na(low_ffmi_v0)) |> 
  mutate(
    across(where(is.character), as.factor),
    low_ffmi_v0 = factor(
      low_ffmi_v0,
      levels = c("yes", "no"),
      labels = c("Low FFMI", "High/moderate FFMI")
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
    `Prediabetes` = prediab_v0,
    `T2D` = t2d_v0,
    `HbA1c%` = hba1c_percent_v0,
    `Fasting glucose (mmol/L)` = fasting_glucose_mmoll_mmt_v0,
    `Fasting insulin (pmol/L)` = fasting_insulin_pmoll_mmt_v0,
    `HOMA-IR` = homa_ir_v0,
    `HOMA-B` = homa_b_v0
  )

# Check distribution of continuous vars
t1_v0_vars |>
  select(`FFMI group`, where(is.numeric)) |>
  pivot_longer(
    cols = -`FFMI group`,
    names_to = "variable",
    values_to = "value"
  ) |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ variable, scales = "free")

# Create table
t1_v0 <- CreateTableOne(
  vars = names(t1_v0_vars)[names(t1_v0_vars) != "FFMI group"],
  strata = "FFMI group",
  data = t1_v0_vars, 
  factorVars = c("Sex", "Prediabetes", "T2D"),
  test = TRUE
)

# Print & format table
t1_v0_matrix <- print(
  t1_v0,
  nonnormal = c(
    "Fasting glucose (mmol/L)", 
    "Fasting insulin (pmol/L)",
    "HbA1c%",
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
write.csv(t1_v0_matrix, file = "tables/t1_v0_matrix.csv", row.names = TRUE)

#### Table 1 grouped by FFMI status at baseline (participants with available shotgun data) ####
# Filter patient IDs out of phyloseq object
shotgun_ids <- sample_data(baria_mb)$Subject_ID

# Filter metadata only to include ids that also appear in the phyloseq object
baria_muscle_mb <- baria_muscle |> 
  filter(id %in% shotgun_ids)

# Checks
nrow(baria_muscle)
nrow(baria_muscle_mb)

n_distinct(baria_muscle$id)
n_distinct(baria_muscle_mb$id)

# Choose vars
t1_v0_mb_vars <- baria_muscle_mb |> 
  select(age_v0, sex, bmi_v0, wc_cm_v0, fm_kg_v0, ffm_kg_v0, ffmi_v0, low_ffmi_v0, prediab_v0, t2d_v0, 
    hba1c_percent_v0, fasting_glucose_mmoll_mmt_v0, fasting_insulin_pmoll_mmt_v0, homa_ir_v0, homa_b_v0) |> 
  filter(!is.na(low_ffmi_v0)) |> 
  mutate(
    across(where(is.character), as.factor),
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
    `Prediabetes` = prediab_v0,
    `T2D` = t2d_v0,
    `HbA1c%` = hba1c_percent_v0,
    `Fasting glucose (mmol/L)` = fasting_glucose_mmoll_mmt_v0,
    `Fasting insulin (pmol/L)` = fasting_insulin_pmoll_mmt_v0,
    `HOMA-IR` = homa_ir_v0,
    `HOMA-B` = homa_b_v0
  )

# Check distribution of continuous vars
t1_v0_mb_vars |>
  select(`FFMI group`, where(is.numeric)) |>
  pivot_longer(
    cols = -`FFMI group`,
    names_to = "variable",
    values_to = "value"
  ) |>
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ variable, scales = "free")

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
    "HbA1c%",
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
write.csv(t1_v0_mb_matrix, file = "tables/t1_v0_mb_matrix.csv", row.names = TRUE)
