# RNA seq data cleaning
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(DESeq2)

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
  
## Save RNA-seq data per tissue
# i. Raw counts
# ii. DESeq2-normalized (VST)
save_rnaseq_tissue <- function(rnaseq_data, tissue_name, min_count = 10, min_prop = 0.5) {

  # Select one tissue; one sample per id 
  tissue_data <- rnaseq_data |> 
    filter(tissue == tissue_name)
  stopifnot(anyDuplicated(tissue_data$id) == 0)

  ## i. Save raw count data
  tissue_data |> 
    select(id, everything(), -Sample, -tissue) |> 
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq.RDS"))

  ## ii. Create DESeq2-style count matrix (genes x samples)
  count_matrix <- tissue_data |> 
    select(starts_with("ENSG")) |> 
    as.matrix() |> 
    t()
  colnames(count_matrix) <- tissue_data$id

  # Create DESeq2 object: rounded kallisto raw order counts & size factors (seq depths & composition)
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_matrix),
    colData = data.frame(id = colnames(count_matrix)),
    design = ~ 1
  )
  dds <- estimateSizeFactors(dds)

  # Gene expression filter:
  # >= 10 normalized reads (min_count)
  # >= 50% of patients (min_prop)
  genes_keep <- rowMeans(counts(dds, normalized = TRUE) >= min_count) >= min_prop
  dds <- dds[genes_keep, ]

  # Variance-stabilizing transformation
  vst_matrix <- assay(vst(dds, blind = TRUE)) # pull out the plain matrix

  # Save VST-normalized matrix
  vst_matrix |> 
    t() |> # id x gene
    as.data.frame() |> 
    rownames_to_column(var = "id") |> 
    mutate(
      id = as.numeric(id),
      size_factor = sizeFactors(dds)
    ) |> 
    relocate(id, size_factor) |> 
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq_vst.RDS"))

}

walk(c("Liver", "Jejunum", "vFat", "subFat"), ~ save_rnaseq_tissue(baria_muscle_rnaseq, .x))

# Save clean RNAseq data for all tissues
saveRDS(baria_muscle_rnaseq, file = "data/processed_data/BARIA_muscle_RNAseq_clean.RDS")