# RNA seq data cleaning
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(DESeq2)
library(MetBrewer)

source("scripts/assets/functions.R")

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
rnaseq <- readRDS("data/raw_data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

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

# Expression filter
min_count <- 10 # >= 10 normalised reads per gene
min_prop <- 0.5 # in >= 50% of samples

# Helper: size factors, expression filter, VST normalization (from DESEq2)
run_vst <- function(counts) {
  dds <- DESeqDataSetFromMatrix(round(counts), colData = data.frame(id = colnames(counts)), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  genes_keep <- rowMeans(counts(dds, normalized = TRUE) >= min_count) >= min_prop
  dds <- dds[genes_keep, ]
  list(dds = dds, vst = assay(vst(dds, blind = TRUE)))
}

# Helper: count matrix (genes x samples) for each tissue separately
make_count_matrix <- function(rnaseq_data, tissue_name) {

  tissue_data <- rnaseq_data |>
    filter(tissue == tissue_name)
  stopifnot(anyDuplicated(tissue_data$id) == 0)   # one sample per patient

  count_matrix <- tissue_data |>
    select(starts_with("ENSG")) |>
    as.matrix() |>
    t()                                          # transpose: genes become rows, samples columns
  colnames(count_matrix) <- tissue_data$id      # column names = patient ids

  count_matrix

}

# QC per tissue: library size + PCA on VST of all samples
qc_rnaseq_tissue <- function(rnaseq_data, tissue_name) {

  count_matrix <- make_count_matrix(rnaseq_data, tissue_name)
  
  # VST on all samples (no exclusions yet)
  qc <- run_vst(count_matrix)

  # PCA: samples as rows
  pca <- prcomp(t(qc$vst))
  var_expl <- round(100 * summary(pca)$importance[2, 1:2], 1)   # % variance explained by PC1 and PC2

  # QC table: one row per sample
  qc_tbl <- tibble(
    tissue = tissue_name,
    id = as.numeric(colnames(count_matrix)),
    library_size = colSums(count_matrix), # total reads per sample (raw counts)
    size_factor = sizeFactors(qc$dds),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  ## QC Plots ##
  # Plot 1: library size per sample (sorted)
  plot_library <- qc_tbl |>
    mutate(id = fct_reorder(as.character(id), library_size)) |>
    ggplot(aes(id, library_size / 1e6)) +
    geom_col() +
    labs(
      title = paste0(tissue_name, ": library size per sample"),
      x = "Sample (sorted)", y = "Library size (million reads)"
    ) +
    theme_minimal_custom() +
    theme(axis.text.x = element_blank())

  # Plot 2: PCA
  plot_pca <- ggplot(qc_tbl, aes(PC1, PC2, colour = library_size / 1e6)) +
    geom_point(size = 2, alpha = 0.8) +
    ggrepel::geom_text_repel(
      aes(label = id),
      size = 1.8,
      colour = "grey30",
      max.overlaps = Inf,       # label every point, even where they crowd
      segment.size = 0.2,       # thin line from label to its point
      min.segment.length = 0    # always draw that line, so each label is traceable to its dot
    ) +
    scale_colour_gradient(low = renoir_15[3], high = renoir_15[11]) +
    labs(
      title = paste0(tissue_name, ": PCA on VST (all samples)"),
      x = paste0("PC1 (", var_expl[1], "%)"),
      y = paste0("PC2 (", var_expl[2], "%)"),
      colour = "Library size\n(million reads)"
    ) +
    theme_minimal_custom()

  ggsave(paste0("results/graphs/RNAseq/qc_librarysize_", tissue_name, ".pdf"), plot_library, width = 8, height = 4)
  ggsave(paste0("results/graphs/RNAseq/qc_pca_", tissue_name, ".pdf"), plot_pca, width = 9, height = 7)
  qc_tbl

}

# Run QC for all tissues
tissues <- c("Liver", "Jejunum", "vFat", "subFat")
qc_all <- map(tissues, ~ qc_rnaseq_tissue(baria_muscle_rnaseq, .x)) |>
  set_names(tissues)

## Outlier exclusions
# Liver: distinct cluster on PC1 (small library size; visible cluster in plot vs. main cloud on inspection)
liver_exclude <- c(361, 362, 377, 382, 383, 386, 388, 390, 395, 399, 401, 402, 405, 407, 408, 409, 413, 415, 421, 424)

# Inspect library sizes + PC coordinates
qc_all$Liver |> 
  filter(id %in% liver_exclude) |> 
  arrange(desc(library_size))   

qc_exclude <- list(Liver = liver_exclude, Jejunum = numeric(0), vFat = numeric(0), subFat = numeric(0))

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