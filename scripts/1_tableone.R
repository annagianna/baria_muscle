## Table 1 for BARIA muscle project
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)

# Data
baria_muscle <- readRDS("data/20260624_BARIA_muscle_clinical.RDS")

# Table 1 grouped by FFMI status at baseline
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
  )

# Check distribution of continuous vars
t1_v0_vars |>
  select(low_ffmi_v0, where(is.numeric)) |>
  pivot_longer(
    cols = -low_ffmi_v0,
    names_to = "variable",
    values_to = "value"
  ) |>
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~ variable, scales = "free")

# Create table
t1_v0 <- CreateTableOne(
  vars = names(t1_v0_vars)[names(t1_v0_vars) != "low_ffmi_v0"],
  strata = "low_ffmi_v0",
  data = t1_v0_vars, 
  factorVars = c("sex", "prediab_v0", "t2d_v0"),
  test = TRUE
)

# Print & format table
t1_v0_matrix <- print(
  t1_v0,
  nonnormal = c(
    "fasting_glucose_mmoll_mmt_v0", 
    "fasting_insulin_pmoll_mmt_v0",
    "hba1c_percent_v0",
    "homa_ir_v0",
    "homa_b_v0"
    ),
    showAllLevels = TRUE,
    quote = FALSE,
    noSpaces = TRUE,
    printToggle = FALSE,
    contDigits = 1,
    catDigits = 1,
    pDigitsw = 3
)

write.csv(t1_v0_matrix, file = "tables/t1_v0_matrix.csv", row.names = TRUE)
