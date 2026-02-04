## Table 1 for BARIA muscle project
## Anna Giannakogeorgou

# load packages
library(tableone)
library(tidyverse)

# Data
baria_muscle_clean <- readRDS("data/260203_BARIA_muscle_clinical.RDS")
?CreateTableOne

levels(baria_muscle_clean$sex)
levels(baria_muscle_clean$t2d_v0)
?factor

# Table 1 grouped by %ASM change at 1y group
t1_1y <- baria_muscle_clean |> 
  filter(!is.na(asm_change_v4_group)) |> 
  mutate(
    asm_change_v4_group = fct_relevel(asm_change_v4_group, "high"),
    sex = fct_relevel(sex, "male"),
    t2d_v0 = fct_relevel(t2d_v0, "no")
  ) |>
  CreateTableOne(
    vars = c("age_v0", "sex", "bmi_v0", "wc_cm_v0", "fm_kg_v0", "ffm_kg_v0", "asm_kg_v0", "t2d_v0", "hba1c_percent_v0", "fasting_glucose_mmoll_mmt_v0", "homa_ir_v0"),
    factorVars = c("sex", "t2d_v0"),
    strata = "asm_change_v4_group",
    test = TRUE,
    addOverall = TRUE
  ) |> 
  print(
    nonnormal = c("bmi_v0", "wc_cm_v0", "fm_kg_v0", "ffm_kg_v0", "asm_kg_v0", "t2d_v0", "hba1c_percent_v0", "fasting_glucose_mmoll_mmt_v0", "homa_ir_v0"),
    showAllLevels = FALSE,
    noSpaces = TRUE,
    pDigits = 3,
    contDigits = 1,
    #missing = TRUE,
    printToggle = FALSE
  ) |>
  as.data.frame() |> 
  select(-Overall, -test) |> 
  rename(
    `High %ASM loss` = high,
    `Low/modest %ASM loss` = `low/modest`,
    `P value` = p
  )

# replace rownames for table export
rownames(t1_1y) <- c("n", "Age (years)", "Sex", "BMI (kg/m²)", "WC (cm)", "FM (kg)", "FFM (kg)", "ASM (kg)", "T2D", "HbA1c (%)", "FPG (mmol/L)", "HOMA-IR")

# write table
write.table(t1_1y, "tables/t1_1y.tsv", sep = "\t")