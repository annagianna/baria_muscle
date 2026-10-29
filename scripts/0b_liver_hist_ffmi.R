## Liver histology table stratified by 1-year %FFMI change group
## Anna Giannakogeorgou

# Packages
library(tidyverse)
library(tableone)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

#### Baseline liver histology by 1-year %FFMI change group ####
histo_labels <- c(
  `Liver histology category` = "masld_category_histo_v0",
  `Steatosis grade` = "steatosis_grade_histo_v0",
  `Ballooning grade` = "ballooning_grade_histo_v0",
  `Lobular inflammation grade` = "inflammation_grade_histo_v0",
  `Fibrosis stage` = "fibrosis_stage_histo_v0",
  `Steatosis (HistoIndex)` = "steatosis_continuous_value_histo_v0",
  `Ballooning (HistoIndex)` = "ballooning_continuous_value_histo_v0",
  `Lobular inflammation (HistoIndex)` = "inflammation_continuous_value_histo_v0",
  `Fibrosis (HistoIndex)` = "fibrosis_continuous_value_histo_v0"
)

# Participants with both a 1-year %FFMI change group and baseline histology
histo_data <- baria_muscle_wide |>
  filter(!is.na(perc_change_ffmi_v4_group), !is.na(masld_category_histo_v0)) |>
  mutate(
    `1-year %FFMI change group` = factor(
      perc_change_ffmi_v4_group,
      levels = c("high", "modest/low"),
      labels = c("High %FFMI change", "Modest/low %FFMI change") # (same labels as Table 1)
    )
  ) |>
  select(`1-year %FFMI change group`, all_of(histo_labels))
count(histo_data, `1-year %FFMI change group`) # n per FFMI group

# Variables
histo_vars <- names(histo_labels)
histo_categorical <- c(
  "Liver histology category", "Steatosis grade", "Ballooning grade",
  "Lobular inflammation grade", "Fibrosis stage"
)
histo_continuous <- setdiff(histo_vars, histo_categorical) # HistoIndex continuous values: median [IQR], Wilcoxon

##  Table
t_histo_ffmi <- CreateTableOne(
  vars = histo_vars,
  strata = "1-year %FFMI change group",
  data = histo_data,
  factorVars = histo_categorical,
  test = TRUE
)

t_histo_ffmi_matrix <- print(
  t_histo_ffmi,
  nonnormal = histo_continuous,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE,
  printToggle = FALSE,
  contDigits = 1,
  catDigits = 1,
  pDigits = 3
)

write.csv(t_histo_ffmi_matrix, "results/tables/ts_histology_by_ffmi_group.csv", row.names = TRUE)
