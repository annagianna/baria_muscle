## Table 1 for BARIA muscle project
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)

# Data
baria_muscle <- readRDS("data/20260624_BARIA_muscle_clinical.RDS")

# Table 1 grouped by FFMI status at baseline


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

t1_v0 <- CreateTableOne(
  vars = names(t1d_v0_vars)[[names(t1_v0_vars) != "low_ffmi_v0"]],
  strata = "low_ffmi_v0",
  data = t1_v0_vars, 
  factorVars = c("sex", "prediab_v0", "t2d_v0"),
  test = TRUE
)




# Export


print(t1_v0)

t1_v0_matrix <- print()

write.csv(t1_v0_matrix, file = "tables/t1_v0.csv")