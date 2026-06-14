## Table 1 for BARIA muscle project
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(gtsummary)
library(gt)

# Data
baria_muscle <- readRDS("data/20260613_BARIA_muscle_clinical.RDS")

### Baseline low muscle mass / FFMI ###
# Table 1 grouped by FFMI status at baseline
t1_v0 <- baria_muscle |> 
  filter(!is.na(low_ffmi_v0)) |> 
  mutate(
    low_ffmi_v0 = fct_relevel(low_ffmi_v0, "yes"),
    t2d_v0 = fct_relevel(t2d_v0, "no"),
    prediab_v0 = fct_relevel(prediab_v0, "no"),
  ) |>
  tbl_summary(
    by = low_ffmi_v0,
    include = c("age_v0", "sex", "bmi_v0", "wc_cm_v0", "fm_kg_v0", "ffm_kg_v0", "ffmi_v0", "t2d_v0", "prediab_v0", "hba1c_percent_v0", "fasting_glucose_mmoll_mmt_v0", "fasting_insulin_pmoll_mmt_v0", "homa_ir_v0", "homa_b_v0"),
    type  = list(sex ~ "dichotomous"),
    value = list(sex ~ "female"),
    label = list(
      "age_v0" = "Age (years)",
      "sex" = "Sex, female n(%)",
      "bmi_v0" = "BMI (kg/m²)",
      "wc_cm_v0" = "WC (cm)",
      "fm_kg_v0" = "FM (kg)",
      "ffm_kg_v0" = "FFM (kg)",
      "ffmi_v0" = "FFMI",
      "t2d_v0" = "T2D (n (%))",
      "prediab_v0" = "Prediabetes (n (%))",
      "hba1c_percent_v0" = "HbA1c (%)",
      "fasting_glucose_mmoll_mmt_v0" = "FPG (mmol/L)",
      "fasting_insulin_pmoll_mmt_v0" = "Fasting insulin (pmol/L)",
      "homa_ir_v0" = "HOMA-IR",
      "homa_b_v0" = "HOMA-B"
    ),
    missing = c("no")
  ) |> 
  add_p() |> 
  bold_p() |> 
  modify_header(
    label = "**Variable**",
    stat_1 = "**Low FFMI**, \nn = {n}",
    stat_2 = "**High/moderate FFMI**, \nn = {n}",
    p.value = "**P**"
  )

# write table
t1_v0 |> 
  as_gt() |> 
  gtsave("tables/t1_v0.html")

t1_v0 |>
  as_tibble() |>
  write_tsv("tables/t1_v0.tsv")

### Low muscle mass / FFMI at year 1 post-operatively ###