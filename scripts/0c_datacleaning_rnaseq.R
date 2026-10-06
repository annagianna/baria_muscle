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

  # PCA on the 500 most variable genes
  # centred, not scaled (as in DESeq2::plotPCA, ntop = 500)
  top_genes <- order(matrixStats::rowVars(qc$vst), decreasing = TRUE)[1:500]
  pca <- prcomp(t(qc$vst[top_genes, ]))
  var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)[1:2] # %variance explained by PC1 and PC2

  # QC table: one row per sample
  qc_tbl <- tibble(
    tissue = tissue_name,
    id = as.numeric(colnames(count_matrix)),
    library_size = colSums(count_matrix), # total reads per sample (raw counts)
    size_factor = sizeFactors(qc$dds),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_var = var_expl[1], # same value in every row; used for the axis labels
    PC2_var  = var_expl[2]
  )
}

## QC Plots ##
# QC plots per tissue: 1. library size, 2. PCA;
plot_qc_tissue <- function(qc_tbl, label_ids = NULL) { # label_ids = ids to label (NULL = label all)

  tissue_name <- unique(qc_tbl$tissue) # gets tissue name out of the table

  # Plot 1: library size per sample (sorted)
  plot_library <- qc_tbl |>
    mutate(id = fct_reorder(as.character(id), library_size)) |>
    ggplot(aes(x = id, y = library_size / 1e6)) +
    geom_col(fill = renoir_15[2]) +
    labs(
      title = paste0(tissue_name, ": library size per sample"),
      x = "Sample (sorted)",
      y = "Library size (million reads)"
    ) +
    theme_minimal_custom() +
    theme(axis.text.x = element_blank())

  # Which points get an ID label
  label_data <- if (is.null(label_ids)) qc_tbl else filter(qc_tbl, id %in% label_ids) # label only outlier cluster IDs with ggrepel

  # Plot 2: PCA, coloured by library size
  plot_pca <- ggplot(qc_tbl, aes(x = PC1, y = PC2, colour = library_size / 1e6)) +
    geom_point(size = 2, alpha = 0.8) +
    ggrepel::geom_label_repel(
      data = label_data,
      aes(label = id),
      size = 2, colour = "grey20", fill = "white", label.size = 0.2,
      max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0
    ) +
    scale_colour_gradient(low = renoir_15[1], high = renoir_15[12]) +
    labs(
      title = paste0(tissue_name, ": PCA on VST (top 500 variable genes, all samples)"),
      x = paste0("PC1 (", unique(qc_tbl$PC1_var), "%)"),
      y = paste0("PC2 (", unique(qc_tbl$PC2_var), "%)"),
      colour = "Library size\n(million reads)"
    ) +
    theme_minimal_custom()

  ggsave(paste0("results/graphs/RNAseq/qc_librarysize_", tissue_name, ".pdf"), plot_library, width = 8, height = 4)
  ggsave(paste0("results/graphs/RNAseq/qc_pca_", tissue_name, ".pdf"), plot_pca, width = 9, height = 7)
}

# Run QC for all tissues (computes VST + PCA once per tissue)
tissues <- c("Liver", "Jejunum", "vFat", "subFat")
qc_all <- map(tissues, ~ qc_rnaseq_tissue(baria_muscle_rnaseq, .x)) |>
  set_names(tissues)

## Outlier exclusions (IDs in gitignored file, not in public repo)
# Liver
qc_reasons <- read_csv("data/processed_data/rnaseq_qc_exclusions.csv", col_types = "cdc") # tissue, id, qc_reason

# i. low depth: samples forming a distinct cluster outside the main cloud on the QC PCA, all with small library size
liver_exclude_lowdepth <- qc_reasons |> 
  filter(tissue == "Liver", qc_reason == "low_depth") |> 
  pull(id)

# ii. adipose tissue profile: clear adipose-tissue expression profile (ADIPOQ, LEP, PLIN1, FABP4, CIDEA; CPM from raw counts)
liver_exclude_adipose <- qc_reasons |> 
  filter(tissue == "Liver", qc_reason == "adipose_profile") |> 
  pull(id)

liver_exclude <- c(liver_exclude_lowdepth, liver_exclude_adipose)

stopifnot(
  length(liver_exclude_lowdepth) == 20,
  length(liver_exclude_adipose) == 2,
  !anyDuplicated(liver_exclude)# no sample listed under both reasons
)

qc_exclude <- map(set_names(tissues), ~ qc_reasons$id[qc_reasons$tissue == .x])

# QC plots: excluded samples labelled
walk(tissues, ~ plot_qc_tissue(qc_all[[.x]], label_ids = qc_exclude[[.x]]))

# QC table: library size, size factor, PC coordinates, exclusion flag per sample and qc reason
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
bind_rows(qc_all) |>
  left_join(qc_reasons, by = c("tissue", "id")) |>   # qc_reason: NA = passed QC
  mutate(qc_excluded = !is.na(qc_reason)) |>
  write_csv("results/tables/rnaseq_qc_samples.csv")

## Save RNA-seq data per tissue
# i. Raw counts: all samples + qc_exclude flag (5a: excluded in main analysis, included in sensitivity analysis)
# ii. VST: excluded samples removed; size factors, expression filter and VST recomputed on the remaining samples
save_rnaseq_tissue <- function(rnaseq_data, tissue_name, exclude_ids, reasons) {

  count_matrix <- make_count_matrix(rnaseq_data, tissue_name)

  # i. Raw counts (samples as rows) with QC flag
  count_matrix |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id") |>
    mutate(id = as.numeric(id)) |>
    left_join(reasons |> filter(tissue == tissue_name) |> select(id, qc_reason), by = "id") |>
    mutate(qc_exclude = id %in% exclude_ids) |>
    relocate(id, qc_exclude, qc_reason) |>
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq.RDS"))

  # ii. VST without excluded samples
  keep_samples <- !(colnames(count_matrix) %in% as.character(exclude_ids))
  clean <- run_vst(count_matrix[, keep_samples, drop = FALSE])

  clean$vst |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id") |>
    mutate(
      id = as.numeric(id),
      size_factor = sizeFactors(clean$dds)
    ) |>
    relocate(id, size_factor) |>
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq_vst.RDS"))

  message(tissue_name, ": ", sum(!keep_samples), " excluded, ", sum(keep_samples), " retained")
}

walk(tissues, ~ save_rnaseq_tissue(baria_muscle_rnaseq, .x, qc_exclude[[.x]], qc_reasons))