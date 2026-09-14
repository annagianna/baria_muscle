# RNA seq data cleaning
# Anna Giannakogeorgou

# Packages
library(tidyverse)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
rnaseq <- readRDS("data/raw_data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")

## Cleaning ##
baria_muscle_rnaseq <- rnaseq |> 
  select(matches("\\.(Liver|subFat|vFat|Jejunum)\\.V1$")) |> 
  select(-matches("^REFER")) |> # exclude duplicated "REFER samples"
  t() |> 
  as.data.frame() |> 
  rownames_to_column("Sample") |> 
  mutate(
    id = as.numeric(str_extract(Sample, "^\\d+")),
    tissue = str_extract(Sample, "(Liver|subFat|vFat|Jejunum)")
  ) |>
  filter(id %in% baria_muscle_wide$id) |> 
  mutate(across(starts_with("ENSG"), ~ as.numeric(str_replace(as.character(.x), ",", ".")))) # fix inconsistent decimals
  
# Save clean RNAseq data for each tissue
save_rnaseq_tissue <- function(rnaseq_data, tissue_name) {
  rnaseq_data |> 
    filter(tissue == tissue_name) |> 
    select(id, everything(), -Sample, -tissue) |> 
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq.RDS"))
}

walk(c("Liver", "Jejunum", "vFat", "subFat"), ~ save_rnaseq_tissue(baria_muscle_rnaseq, .x))

# Save clean RNAseq data for all tissues
saveRDS(baria_muscle_rnaseq, file = "data/processed_data/BARIA_muscle_RNAseq_clean.RDS")
