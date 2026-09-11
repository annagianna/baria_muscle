## Table 1 stratified by 1-year %FFMI change group
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)
library(gt)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")

# Baseline & 1y vars
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
    perc_change_ffmi_v4_group = factor(perc_change_ffmi_v4_group, levels = c("high", "modest/low"), labels = c("High %FFMI change", "Modest/low %FFMI change")),
    ffmi_group_visit = factor(
      interaction(perc_change_ffmi_v4_group, visit, sep = " - "),
      levels = c("High %FFMI change - Baseline", "Modest/low %FFMI change - Baseline", "High %FFMI change - 1 year", "Modest/low %FFMI change - 1 year")
    )
  ) |> 
  rename(
    `Age at baseline (years)` = age_v0,
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
    `1-year %FFMI change group` = perc_change_ffmi_v4_group
  )

# Table 1 vars (both v0 and v4)
t1_vars <- c(
  "Age at baseline (years)", "Sex", 
  "BMI (kg/m²)", "Waist circumference (cm)", "Fat mass (kg)", "Fat-free mass (kg)", "Fat-free mass index (kg/m²)",
  "Prediabetes", "T2D", "HbA1c (mmol/mol)", "Fasting glucose (mmol/L)", "Fasting insulin (pmol/L)", "HOMA-IR", "HOMA-B",
  "Total cholesterol (mmol/L)", "LDL cholesterol (mmol/L)", "HDL cholesterol (mmol/L)", "Triglycerides (mmol/L)",
  "Lipid-lowering medication", "Antidiabetic medication"
)

# Variable types
nonnormal_vars <- c(
  "HbA1c (mmol/mol)", "Fasting glucose (mmol/L)", "Fasting insulin (pmol/L)", "HOMA-IR", "HOMA-B",
  "Total cholesterol (mmol/L)", "LDL cholesterol (mmol/L)", "HDL cholesterol (mmol/L)", "Triglycerides (mmol/L)"
)
categorical_vars <- c("Sex", "Prediabetes", "T2D", "Lipid-lowering medication", "Antidiabetic medication")
normal_vars <- setdiff(t1_vars, c(nonnormal_vars, categorical_vars))

## Baseline ##
# Baseline table stratified by %FFMI change group at 1y
t1_ffmi_loss_v0 <- CreateTableOne(
  vars = t1_vars,
  strata = "1-year %FFMI change group",
  data = t1_ffmi_loss_vars_long |>
    filter(visit == "Baseline"),
  factorVars = categorical_vars,
  test = TRUE
)
# Print baseline table
t1_ffmi_loss_v0_matrix <- print(
  t1_ffmi_loss_v0,
  nonnormal = nonnormal_vars,
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
  strata = "1-year %FFMI change group",
  data = t1_ffmi_loss_vars_long |>
    filter(visit == "1 year"),
  factorVars = categorical_vars,
  test = TRUE
)

# Print 1y table
t1_ffmi_loss_v4_matrix <- print(
  t1_ffmi_loss_v4,
  nonnormal = nonnormal_vars,
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

#### Combined baseline & 1y table ####
t1_ffmi_loss_v0_v4 <- CreateTableOne(
  vars = t1_vars,
  strata = "ffmi_group_visit",
  data = t1_ffmi_loss_vars_long,
  factorVars = categorical_vars,
  test = FALSE
)

# Print
t1_ffmi_loss_v0_v4_matrix <- print(
  t1_ffmi_loss_v0_v4,
  nonnormal = nonnormal_vars,
  showAllLevels = FALSE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1,
  catDigits = 1
)
write.csv(t1_ffmi_loss_v0_v4_matrix, "results/tables/t1_ffmi_loss_v0_v4_matrix.csv", row.names = TRUE)

# Wide data for paired comparisons (v0 vs. v4)
t1_ffmi_loss_vars_wide <- t1_ffmi_loss_vars_long |> 
  select(id, visit, `1-year %FFMI change group`, all_of(t1_vars)) |> 
  pivot_wider(
    names_from = visit,
    values_from = all_of(t1_vars),
    names_glue = "{.value}_{visit}"
  )

### Tests ###
# High vs. modest/low FFMI groups within one visit
run_test_between <- function(var, visit_name, var_type) {

  visit_data <- t1_ffmi_loss_vars_long |>
    filter(visit == visit_name)

  if (var_type == "normal") {
    test <- t.test(visit_data[[var]] ~ visit_data[["1-year %FFMI change group"]])
  } else if (var_type == "nonnormal") {
    test <- wilcox.test(visit_data[[var]] ~ visit_data[["1-year %FFMI change group"]], exact = FALSE)
  } else if (var_type == "categorical") {
    test <- chisq.test(table(visit_data[[var]], visit_data[["1-year %FFMI change group"]]))
  }
  test$p.value
}

# Paired comparisons baseline -> 1 year
run_test_paired <- function(var, ffmi_group, var_type) {

  group_data <- t1_ffmi_loss_vars_wide |>
    filter(`1-year %FFMI change group` == ffmi_group)

  var_v0 <- paste0(var, "_Baseline")
  var_v4 <- paste0(var, "_1 year")

  if (var_type == "normal") {
    test <- t.test(group_data[[var_v0]], group_data[[var_v4]], paired = TRUE)
  } else if (var_type == "nonnormal") {
    test <- wilcox.test(group_data[[var_v0]], group_data[[var_v4]], paired = TRUE, exact = FALSE)
  } else if (var_type == "categorical") {
    test <- mcnemar.test(table(group_data[[var_v0]], group_data[[var_v4]]))
  }
  test$p.value
}

# Function to run comparisons across var types
run_all_tests <- function(vars, var_type) {

  tibble(
    var = vars,
    p_v0_high_low = map_dbl(vars, \(var) run_test_between(var, "Baseline", var_type)),
    p_v4_high_low = map_dbl(vars, \(var) run_test_between(var, "1 year", var_type)),
    p_high_v0_v4 = map_dbl(
      vars, \(var) {if (var %in% c("Age at baseline (years)", "Sex")) {NA_real_} else {run_test_paired(var, "High %FFMI change", var_type)}}
    ),
    p_low_v0_v4 = map_dbl(
      vars, \(var) {if (var %in% c("Age at baseline (years)", "Sex")) {NA_real_} else {run_test_paired(var, "Modest/low %FFMI change", var_type)}}
    )
  )
}

## Run tests ##
p_stars <- function(p) {case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")}

test_results <- bind_rows(
  run_all_tests(normal_vars, "normal"),
  run_all_tests(nonnormal_vars, "nonnormal"),
  run_all_tests(categorical_vars, "categorical")
) |>
  mutate(
    sig_v0_high_low = if_else(!is.na(p_v0_high_low) & p_v0_high_low < 0.05, paste0("ᵃ", p_stars(p_v0_high_low)), ""),
    sig_high_v0_v4 = if_else(!is.na(p_high_v0_v4) & p_high_v0_v4 < 0.05, paste0("ᵇ", p_stars(p_high_v0_v4)), ""),
    sig_low_v0_v4 = if_else(!is.na(p_low_v0_v4) & p_low_v0_v4 < 0.05, paste0("ᶜ", p_stars(p_low_v0_v4)), ""),
    sig_v4_high_low = if_else(!is.na(p_v4_high_low) & p_v4_high_low < 0.05, paste0("ᵈ", p_stars(p_v4_high_low)), "")
  )

# Format final table data
t1_ffmi_loss_v0_v4_final <- t1_ffmi_loss_v0_v4_matrix |> 
  as.data.frame() |> 
  rownames_to_column(var = "Variable") |> 
  mutate(
    var = Variable |>
      str_remove(" \\(mean \\(SD\\)\\)$") |>
      str_remove(" \\(median \\[IQR\\]\\)$") |>
      str_remove(" =.*$")
  ) |>
  rename(
    v0_high = `High %FFMI change - Baseline`,
    v0_modest_low = `Modest/low %FFMI change - Baseline`,
    v4_high = `High %FFMI change - 1 year`,
    v4_modest_low = `Modest/low %FFMI change - 1 year`
  ) |> 
  left_join(
    test_results |> 
      select(var, sig_v0_high_low, sig_high_v0_v4, sig_low_v0_v4, sig_v4_high_low),
    by = "var"
  ) |> 
  mutate(across(starts_with("sig_"), ~ replace_na(.x, ""))) |> 
  mutate(
    v0_modest_low = paste0(v0_modest_low, sig_v0_high_low),
    v4_high = paste0(v4_high, sig_high_v0_v4),
    v4_modest_low = paste0(v4_modest_low, sig_low_v0_v4, sig_v4_high_low)
  ) |>
  select(Variable, v0_high, v0_modest_low, v4_high, v4_modest_low)

# Save csv
write.csv(t1_ffmi_loss_v0_v4_final, "results/tables/t1_ffmi_loss_v0_v4_final.csv", row.names = FALSE)

# Format as gt
t1_ffmi_loss_v0_v4_gt <- t1_ffmi_loss_v0_v4_final |> 
  gt(rowname_col = "Variable") |> 
  tab_spanner(label = "%FFMI change", columns = everything(), level = 1) |> 
  tab_spanner(label = "Baseline", columns = contains("v0"), level = 2) |> 
  tab_spanner(label = "1 year", columns = contains("v4"), level = 2) |> 
  cols_label(
    v0_high = "High",
    v0_modest_low = "Modest/low",
    v4_high = "High",
    v4_modest_low = "Modest/low"
  ) |> 
  tab_source_note(
    source_note = "Groups are based on 1-year %FFMI change. ᵃ between-group comparison at baseline; ᵇ within high group from baseline to 1 year; ᶜ within modest/low group from baseline to 1 year; ᵈ between-group comparison at 1 year. * p < 0.05; ** p < 0.01; *** p < 0.001."
  )

# Save final gt table
Sys.unsetenv("CHROMOTE_CHROME")
gtsave(t1_ffmi_loss_v0_v4_gt, "results/tables/t1_ffmi_loss_v0_v4_gt.html")