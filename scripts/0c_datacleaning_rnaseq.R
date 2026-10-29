# RNA seq data cleaning
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(DESeq2)
library(MetBrewer)
library(annotables)
library(readxl)
library(limma)
source("scripts/assets/functions.R")

# Bioconductor packages (DESeq2/S4Vectors/matrixStats) mask some dplyr verbs; always use dplyr's
select <- dplyr::select
rename <- dplyr::rename
count  <- dplyr::count
filter <- dplyr::filter

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
rnaseq <- readRDS("data/raw_data/RNASeq.Counttable.kallisto.39546.1453.2025-01.29.RDS")
omics_metadata <- "data/raw_data/BARIA_omics_metadata_V19.xlsx"
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

## Sample metadata excel table
rnaseq_meta <- read_excel(omics_metadata, sheet = "RNAseq") |>
  filter(Tissue == "Liver", Time_Point == "V-1", Study == "BARIA") |> # Liver (v1) only for now
  mutate(
    label = Sample_id |> # label = original sample name as in the kallisto count table
      str_remove("^BARIA\\.RNAseq\\.") |> 
      str_remove("\\.V-1\\.NA$"),
    id_corrected = as.numeric(coalesce(na_if(as.character(`SNP check`), "Failed"), as.character(Subject_ID))), # SNP-verified ID where available; "Failed" = identity not verifiable, label ID kept
    batch = Batch,
    exclude = Exclude %in% "T",
    comment = Comment,
    .keep = "none" # keep only the new columns
  )
stopifnot(!anyDuplicated(rnaseq_meta$label)) # one metadata row per sample

## QC exclusions (own QC, in addition to the omics metadata exclusions)
# Low-quality liver samples: distinct cluster on the QC PCA (PC1) with low library size, all within sequencing batch 2024
# (batch 2024: n = 135; these 20 differ from the other 115 samples in the same batch, so not a batch effect)
low_depth_ids <- c(361, 362, 377, 382, 383, 386, 388, 390, 395, 399, 401, 402, 405, 407, 408, 409, 413, 415, 421, 424)

## Cleaning ##
baria_muscle_rnaseq <- rnaseq |> 
  select(matches("\\.(Liver|subFat|vFat|Jejunum)\\.V1$")) |> 
  select(-matches("^REFER")) |> # exclude duplicated "REFER samples"
  t() |> 
  as.data.frame() |> 
  rownames_to_column("Sample") |> 
  mutate(
    label = str_remove(Sample, "\\.V1$"),
    id = as.numeric(str_extract(Sample, "^\\d+")),
    tissue = str_extract(Sample, "(Liver|subFat|vFat|Jejunum)")
  ) |>
  left_join(rnaseq_meta, by = "label") |>
  mutate(
    # identity corrections from xlsx SNP check col / subject IDs; samples flagged for exclusion keep their label id
    id = if_else(!is.na(id_corrected) & !exclude, id_corrected, id),
    exclude = coalesce(exclude, FALSE),
    exclude = exclude | (tissue == "Liver" & batch %in% "RNA_Seq_2022_Fat"),
    comment = if_else(tissue == "Liver" & batch %in% "RNA_Seq_2022_Fat", "Sequencing batch with n = 2 liver samples; batch effect not estimable", comment),
    exclude = exclude | (tissue == "Liver" & id %in% low_depth_ids),
    comment = if_else(tissue == "Liver" & id %in% low_depth_ids, "Low sequencing depth; distinct cluster on QC PCA (within batch 2024)", comment)
  ) |>
  filter(id %in% baria_muscle_wide$id) |>
  mutate(across(starts_with("ENSG"), ~ as.numeric(str_replace(as.character(.x), ",", ".")))) # fix inconsistent decimals

# Checks (liver)
liver_meta_check <- baria_muscle_rnaseq |> filter(tissue == "Liver")
count(liver_meta_check, batch)                                   # every liver sample needs a batch
stopifnot(!any(is.na(liver_meta_check$batch)))
liver_meta_check |>
  filter(id != as.numeric(str_extract(label, "^\\d+")) | exclude) |>
  select(label, id, batch, exclude, comment)                    # corrected and excluded samples

# Sample table: one row per sample with batch and exclusion flag
sample_meta <- baria_muscle_rnaseq |>
  transmute(
    tissue, id,
    batch = factor(batch),
    qc_exclude = exclude,
    qc_reason = if_else(exclude, comment, NA_character_)
  )

# Expression filter
min_count <- 10 # >= 10 normalised reads per gene
min_prop <- 0.5 # in >= 50% of samples

# Remove unnamed genes
named_genes <- grch38 |> 
  filter(!is.na(symbol), str_trim(symbol) != "") |> 
  pull(ensgene) |> 
  unique()

ensg_cols <- baria_muscle_rnaseq |> 
  select(starts_with("ENSG")) |> 
  names()

unnamed_genes <- ensg_cols[!str_remove(ensg_cols, "\\.\\d+") %in% named_genes]

baria_muscle_rnaseq_named <- baria_muscle_rnaseq |> 
  select(-all_of(unnamed_genes))

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
    t() # transpose: genes x samples
  colnames(count_matrix) <- tissue_data$id

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
  ) |>
    left_join(sample_meta, by = c("tissue", "id"))
}

## QC Plots ##
# QC plots per tissue: 1. library size, 2. PCA coloured by sequencing batch; excluded samples labelled
plot_qc_tissue <- function(qc_tbl) {

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

  # Plot 2: PCA, coloured by sequencing batch
  plot_pca <- ggplot(qc_tbl, aes(x = PC1, y = PC2, colour = batch)) +
    geom_point(size = 2, alpha = 0.8) +
    ggrepel::geom_label_repel(
      data = filter(qc_tbl, qc_exclude),
      aes(label = id),
      size = 2, colour = "grey20", fill = "white", label.size = 0.2,
      max.overlaps = Inf, segment.size = 0.2, min.segment.length = 0
    ) +
    scale_colour_manual(values = renoir_15[c(2, 6, 11, 14)], na.value = "grey60") +
    labs(
      title = paste0(tissue_name, ": PCA on VST (top 500 variable genes, all samples)"),
      x = paste0("PC1 (", unique(qc_tbl$PC1_var), "%)"),
      y = paste0("PC2 (", unique(qc_tbl$PC2_var), "%)"),
      colour = "Sequencing batch"
    ) +
    theme_minimal_custom()

  ggsave(paste0("results/graphs/RNAseq/qc_librarysize_", tissue_name, ".pdf"), plot_library, width = 8, height = 4)
  ggsave(paste0("results/graphs/RNAseq/qc_pca_", tissue_name, ".pdf"), plot_pca, width = 9, height = 7)
}

# Run QC for all tissues (computes VST + PCA once per tissue)
tissues <- c("Liver", "Jejunum", "vFat", "subFat")
qc_all <- map(tissues, ~ qc_rnaseq_tissue(baria_muscle_rnaseq_named, .x)) |>
  set_names(tissues)

# QC check: low-depth samples vs the rest of their batch (library size, PC1)
qc_all$Liver |>
  mutate(low_depth = id %in% low_depth_ids) |>
  group_by(batch, low_depth) |>
  summarise(n = n(), median_lib_M = median(library_size) / 1e6, median_PC1 = median(PC1), .groups = "drop")
stopifnot(sum(qc_all$Liver$id %in% low_depth_ids) == 20)

# QC plots and QC table (library size, size factor, PC coordinates, batch, exclusion flag and reason)
walk(qc_all, plot_qc_tissue)

bind_rows(qc_all) |>
  write_csv("results/tables/rnaseq_qc_samples.csv")

## Save RNA-seq data per tissue
# i. Raw counts: all samples + batch + qc_exclude flag (DESeq2 in 5a/5b: excluded samples dropped, batch in the design)
# ii. VST: excluded samples removed; size factors, expression filter and VST recomputed on the remaining samples
# iii. Batch-corrected VST (limma::removeBatchEffect) for correlation analyses; only for tissues with > 1 known batch
save_rnaseq_tissue <- function(rnaseq_data, tissue_name) {

  count_matrix <- make_count_matrix(rnaseq_data, tissue_name)
  meta <- sample_meta |>
    filter(tissue == tissue_name) |>
    dplyr::slice(match(as.numeric(colnames(count_matrix)), id)) # same order as the count matrix columns
  stopifnot(identical(meta$id, as.numeric(colnames(count_matrix))))

  # i. Raw counts (samples as rows) with batch and QC flag
  count_matrix |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id") |>
    mutate(id = as.numeric(id)) |>
    left_join(meta |> select(id, batch, qc_exclude, qc_reason), by = "id") |>
    relocate(id, batch, qc_exclude, qc_reason) |>
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq.RDS"))

  # ii. VST without excluded samples
  keep_samples <- !meta$qc_exclude
  meta_keep <- meta[keep_samples, ]
  clean <- run_vst(count_matrix[, keep_samples, drop = FALSE])

  clean$vst |>
    t() |>
    as.data.frame() |>
    rownames_to_column("id") |>
    mutate(
      id = as.numeric(id),
      size_factor = sizeFactors(clean$dds),
      batch = meta_keep$batch
    ) |>
    relocate(id, size_factor, batch) |>
    saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq_vst.RDS"))

  # iii. Batch-corrected VST for correlation analyses
  if (!anyNA(meta_keep$batch) && n_distinct(meta_keep$batch) > 1) {
    removeBatchEffect(clean$vst, batch = droplevels(meta_keep$batch)) |>
      t() |>
      as.data.frame() |>
      rownames_to_column("id") |>
      mutate(
        id = as.numeric(id),
        size_factor = sizeFactors(clean$dds),
        batch = meta_keep$batch
      ) |>
      relocate(id, size_factor, batch) |>
      saveRDS(paste0("data/processed_data/BARIA_", tissue_name, "_RNAseq_vst_batchcorrected.RDS"))
  }

  message(tissue_name, ": ", sum(!keep_samples), " excluded, ", sum(keep_samples), " retained")
}

walk(tissues, ~ save_rnaseq_tissue(baria_muscle_rnaseq_named, .x))
