# Jejunum RNA-seq - DE high vs. modest/low %FFMI loss at 1y
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(MetBrewer)
library(annotables)
library(ggpubr)
library(broom)
library(DESeq2)
library(apeglm)
library(ggrepel)
source("scripts/assets/functions.R")
dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

# Bioconductor packages (DESeq2/S4Vectors/matrixStats) mask some dplyr verbs; always use dplyr's
select <- dplyr::select
rename <- dplyr::rename
count  <- dplyr::count
filter <- dplyr::filter

tissue <- "Jejunum"
prefix <- "jejunum"
tissue_label <- "Jejunum"

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)
renoir_cols_20 <- met.brewer("Renoir", n = 20)
de_cols <- c( # same high/low colors used in the alpha, beta diveristy scripts
  "Higher in high loss" = renoir_cols_20[18],
  "Higher in modest/low loss"  = renoir_cols_20[5],
  "not significant" = "grey80"
)

# Significance asterisks based on p_fdr
stars <- function(p_fdr) ifelse(p_fdr < 0.001, "***", ifelse(p_fdr < 0.01, "**", ifelse(p_fdr < 0.05, "*", "")))

# Long cor_block output -> wide matrix (rows = x, cols = y)
to_matrix <- function(df, value) {
  df |>
    dplyr::select(x, y, all_of(value)) |>
    pivot_wider(names_from = y, values_from = all_of(value)) |>
    column_to_rownames("x") |>
    as.matrix()
}

# Helper function: Spearman correlations for every x–y pair, BH-correction within the block
cor_block <- function(df, x_vars, y_vars) {
  expand_grid(x = x_vars, y = y_vars) |>
    mutate(
      n = map2_int(x, y, ~ sum(complete.cases(df[[.x]], df[[.y]]))),
      test = map2(x, y, ~ cor.test(df[[.x]], df[[.y]], method = "spearman", exact = FALSE)),
      rho = map_dbl(test, "estimate"),
      p.value = map_dbl(test, "p.value")
    ) |>
    dplyr::select(-test) |>
    mutate(p_fdr = p.adjust(p.value, method = "BH")) |>
    arrange(p.value)
}

# Data
tissue_rnaseq <- readRDS(paste0("data/processed_data/BARIA_", tissue, "_RNAseq.RDS")) # raw counts -> for DESeq2
vst_file <- paste0("data/processed_data/BARIA_", tissue, "_RNAseq_vst_batchcorrected.RDS") # 0c only writes this for tissues with > 1 batch
if (!file.exists(vst_file)) vst_file <- paste0("data/processed_data/BARIA_", tissue, "_RNAseq_vst.RDS")
tissue_rnaseq_vst <- readRDS(vst_file) # (batch-corrected) VST -> for correlations
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_mb_v0 <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
metab <- readRDS("data/processed_data/BARIA_metabolon_clean.RDS")

# Prep data
# Top-15 species (dichotomous ML model %FFMI change at 1y)
top15_species_labels <- get_top_species()
top15_species <- top15_species_labels$species

# Baseline abundance of top-15 species
top15_species_v0 <- t(as(otu_table(baria_mb_v0), "matrix")) |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  mutate(id = Sample |> 
    str_remove("^BARIA_") |> 
    str_remove("_v0$") |> 
    as.numeric()
) |>
  dplyr::select(id, all_of(top15_species))

# Log10-transform abundance & species-specific pseudocount
top15_species_log10 <- top15_species_v0 |>
  mutate(across(all_of(top15_species), ~ {
    pseudo <- min(.x[.x > 0]) / 2
    log10(.x + pseudo)
  }))

# Map gene ids to their names
gene_ids <- tibble(
  ensembl_gene_id = tissue_rnaseq_vst |> 
    dplyr::select(starts_with("ENSG")) |> 
    colnames()
) |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$"))

# Gene annotations
gene_annotations <- grch38 |>
  filter(ensgene %in% gene_ids$ensgene) |>
  dplyr::select(ensgene, symbol, biotype, description) |>
  distinct()

# Display labels for heatmaps: "SYMBOL – short description" / only symbol / only description
gene_label <- function(symbols, type = c("both", "symbol", "description"), max_char = 45) {

  type <- match.arg(type)

  desc <- gene_annotations$description[match(symbols, gene_annotations$symbol)] |>
    str_remove("\\s*\\[Source.*\\]$") |> # grch38 descriptions can carry a "[Source:...]" suffix
    str_trunc(max_char) |>
    coalesce(symbols) # fall back to symbol if no description

  switch(type,
         both = paste0(symbols, " - ", desc),
         symbol = symbols,
         description = desc)
}

# Species label colours by direction of association with high %FFMI loss at 1y (logistic OR)
species_direction_colors <- c(`Higher in modest/low loss` = renoir_15[15], `Higher in high loss` = renoir_15[9])
species_label_colors <- setNames(species_direction_colors[top15_species_labels$estimate_direction], top15_species_labels$species_label)

#### DE: high vs modest/low %FFMI loss (1y) ####
# Prep data
coldata_grp_pre <- baria_muscle_wide |>
  mutate(id = as.numeric(as.character(id))) |>
  inner_join(tissue_rnaseq |> filter(!qc_exclude) |> select(id, batch), by = "id") |>
  filter(!is.na(perc_change_ffmi_v4_group))

# Complete cases for age, FMI and %TWL: main and sensitivity models use the same samples
coldata_grp <- coldata_grp_pre |>
  filter(!is.na(age_v0), !is.na(fmi_v0), !is.na(perc_change_weight_kg_v4)) |>
  mutate(
    ffmi_group = factor(if_else(perc_change_ffmi_v4_group == "high", "high", "modest_low"), levels = c("modest_low", "high")), # Positive LFC = higher in high loss
    age_v0_z = as.numeric(scale(age_v0)),
    fmi_v0_z = as.numeric(scale(fmi_v0)),
    twl_z = as.numeric(scale(perc_change_weight_kg_v4)),
    batch = droplevels(batch) # drop batch levels without samples in the analyzed population
  ) |>
  arrange(id) |>
  select(id, sex, age_v0_z, fmi_v0_z, batch, ffmi_group, twl_z)
message("Complete-case filter drops ", nrow(coldata_grp_pre) - nrow(coldata_grp), " of ", nrow(coldata_grp_pre), " participants")
count(coldata_grp, batch, ffmi_group)
has_batch <- n_distinct(coldata_grp$batch) > 1 # batch only in the designs if > 1 batch in the analysed population

count_data_grp <- tissue_rnaseq |>
  filter(!qc_exclude, id %in% coldata_grp$id) |>
  arrange(id) |>
  column_to_rownames(var = "id") |>
  select(starts_with("ENSG")) |>
  as.matrix() |> 
  t() |> 
  round()
stopifnot(identical(colnames(count_data_grp), as.character(coldata_grp$id)))

design_grp <- if (has_batch) ~ sex + age_v0_z + batch + twl_z + ffmi_group else ~ sex + age_v0_z + twl_z + ffmi_group
dds_grp <- DESeqDataSetFromMatrix(count_data_grp, coldata_grp, design = design_grp)
dds_grp <- estimateSizeFactors(dds_grp)
dds_grp <- dds_grp[rowMeans(counts(dds_grp, normalized = TRUE) >= 10) >= 0.5, ]
nrow(dds_grp) # n of genes tested
dds_grp <- DESeq(dds_grp)
resultsNames(dds_grp)

res_grp <- results(dds_grp, name = "ffmi_group_high_vs_modest_low", alpha = 0.05)
res_grp_shr <- lfcShrink(dds_grp, coef = "ffmi_group_high_vs_modest_low", type = "apeglm", res = res_grp) # apeglm shrunken results; padj unchanged

# DE results (BH across all tested genes)
de_grp <- as.data.frame(res_grp) |>
  rownames_to_column("ensembl_gene_id") |>
  select(ensembl_gene_id, baseMean, log2FC_mle = log2FoldChange, pvalue, padj) |>
  left_join(
    as.data.frame(res_grp_shr) |> 
      rownames_to_column("ensembl_gene_id") |>
      select(ensembl_gene_id, log2FC_shr = log2FoldChange, lfcSE_shr = lfcSE),
    by = "ensembl_gene_id") |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
  left_join(gene_annotations |> distinct(ensgene, .keep_all = TRUE) |> select(ensgene, symbol, description), by = "ensgene") |>
  mutate(symbol = if_else(is.na(symbol) | symbol == "", ensgene, symbol)) |> # unannotated genes keep their Ensembl id
  arrange(pvalue)

de_sig <- de_grp |> 
  filter(padj < 0.05)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write_csv(de_grp, paste0("results/tables/", prefix, "_de_grp_all_genes.csv"))

de_sig |> 
  select(symbol, description, log2FC_mle, log2FC_shr, padj)

# MA plot
ma_df <- bind_rows(
  plotMA(res_grp, returnData = TRUE) |> mutate(type = "Before shrinkage (MLE)"),
  plotMA(res_grp_shr, returnData = TRUE) |> mutate(type = "After shrinkage (apeglm)")
) |>
  mutate(type = factor(type, levels = c("Before shrinkage (MLE)", "After shrinkage (apeglm)")))

ma_plot_grp <- ggplot(ma_df, aes(mean, lfc, colour = isDE)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_x_log10() +
  scale_colour_manual(
    values = c(`FALSE` = "grey75", `TRUE` = renoir_15[11]),
    labels = c("not significant", "q < 0.05"),
    name = NULL
  ) +
  facet_wrap(~ type) +
  labs(x = "Mean of normalised counts", y = "log2FC high vs modest/low loss") +
  theme_minimal_custom()
ggsave(paste0("results/graphs/RNAseq/", prefix, "_de_grp_MA_before_after_shrinkage.pdf"), ma_plot_grp, width = 10, height = 5)

## Forest plot (group DE)
if (nrow(de_sig) > 0) {
  # Data
  forest_de_grp_df <- de_sig |> 
    mutate(
      direction = if_else(log2FC_shr > 0, "Higher in high loss", "Higher in modest/low loss"),
      lower = log2FC_shr - 1.96 * lfcSE_shr, # approx. 95% interval
      upper = log2FC_shr + 1.96 * lfcSE_shr,
      label = fct_reorder(gene_label(symbol), log2FC_shr)
    )

  # Plot
  forest_de_grp_plot <- forest_de_grp_df |> 
    ggplot(aes(x = log2FC_shr, y = label, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = de_cols, name = NULL) +
    labs(
      title = "Differences in gene expression",
      x = "log2FC high vs modest/low loss (approx. 95% CI)",
      y = NULL
    ) +
    theme_minimal_custom()
  ggsave(paste0("results/graphs/RNAseq/", prefix, "_forest_de_grp_plot.pdf"), forest_de_grp_plot, width = 9, height = 1.5 + 0.4 * nrow(forest_de_grp_df))
} else {
  message("No significant genes in the main DE: forest plot skipped")
}

## Volcano plot (group DE)
# Data
volcano_df <- de_grp |>
  filter(!is.na(padj)) |>
  mutate(
    direction = case_when(
      padj < 0.05 & log2FC_shr > 0 ~ "Higher in high loss",
      padj < 0.05 & log2FC_shr < 0 ~ "Higher in modest/low loss",
      TRUE ~ "not significant"
    )
  )

# Plot
volcano_plot <- volcano_df |>
  ggplot(aes(x = log2FC_shr, y = -log10(pvalue), colour = direction)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = volcano_df |> slice_min(padj, n = 15, with_ties = FALSE),
    aes(label = symbol),
    size = 2.5,
    max.overlaps = Inf,
    show.legend = FALSE, 
    seed = 13
  ) +
  scale_colour_manual(values = de_cols, name = NULL) +
  labs(
    title = "Differences in gene expression",
    x = "log2FC high vs modest/low loss (apeglm)",
    y = "-log10(p)"
  ) +
  theme_minimal_custom()
ggsave(paste0("results/graphs/RNAseq/", prefix, "_de_grp_volcano.pdf"), volcano_plot, width = 7, height = 6)

# Sensitivity model: main model additionally adj. for baseline FMI
dds_grp_fmi <- DESeqDataSetFromMatrix(
  count_data_grp,
  coldata_grp,
  design = if (has_batch) ~ sex + age_v0_z + batch + twl_z + fmi_v0_z + ffmi_group else ~ sex + age_v0_z + twl_z + fmi_v0_z + ffmi_group
)
dds_grp_fmi <- estimateSizeFactors(dds_grp_fmi)
dds_grp_fmi <- dds_grp_fmi[rowMeans(counts(dds_grp_fmi, normalized = TRUE) >= 10) >= 0.5, ]
dds_grp_fmi <- DESeq(dds_grp_fmi)
res_grp_fmi <- results(dds_grp_fmi, name = "ffmi_group_high_vs_modest_low", alpha = 0.05)

de_grp_fmi <- as.data.frame(res_grp_fmi) |>
  rownames_to_column("ensembl_gene_id") |>
  select(
    ensembl_gene_id,
    log2FC_mle_fmi = log2FoldChange,
    pvalue_fmi = pvalue,
    padj_fmi = padj
  )

# Compare main vs. FMI-adj. DE results
de_fmi_compare <- de_sig |>
  select(ensembl_gene_id, symbol, log2FC_mle, padj) |>
  left_join(de_grp_fmi, by = "ensembl_gene_id") |>
  mutate(
    same_direction = sign(log2FC_mle) == sign(log2FC_mle_fmi),
    still_sig = padj_fmi < 0.05
  )
de_fmi_compare

## Batch diagnostics: genes differing by sequencing batch (LRT: full model vs reduced model without batch)
if (has_batch) {
  dds_lrt <- DESeq(dds_grp, test = "LRT", reduced = ~ sex + age_v0_z + twl_z + ffmi_group)
  res_lrt <- results(dds_lrt, alpha = 0.05)
  c(genes_tested = nrow(res_lrt), genes_differing_by_batch = sum(res_lrt$padj < 0.05, na.rm = TRUE))
} else {
  message("Single sequencing batch for ", tissue, ": batch LRT skipped")
}

#### Tissue RNA-seq x top15 species ####
# Top15 species x tissue RNA seq
tissue_mb <- tissue_rnaseq_vst |>
  mutate(id = as.numeric(as.character(id))) |>
  inner_join(top15_species_log10, by = "id") |>
  semi_join(
    baria_muscle_wide |>
      filter(!is.na(perc_change_ffmi_v4_group)) |>
      mutate(id = as.numeric(as.character(id))),
    by = "id"
  )

# QC: sequencing depth (DESeq2 size factor) vs. species abundance (no expected association)
size_factor_qc <- map_dfr(top15_species,
  ~ broom::tidy(cor.test(tissue_mb$size_factor, tissue_mb[[.x]], method = "spearman", exact = FALSE)
) |> mutate(species = .x)) |>
  dplyr::select(species, rho = estimate, p.value)
size_factor_qc

# Spearman correlations: top-15 species x DE genes (BH within the block)
de_genes_vst <- intersect(de_sig$ensembl_gene_id, names(tissue_mb)) # DE genes present in the VST file
cor_de_species <- cor_block(tissue_mb, top15_species, de_genes_vst) |>
  rename(species = x, ensembl_gene_id = y) |>
  left_join(
    de_sig |> 
      select(ensembl_gene_id, symbol), 
    by = "ensembl_gene_id"
  )

#### I. DE genes x top 15 species ####
# FDR-significant species x DE gene associations
sig_de_species <- cor_de_species |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

if (nrow(sig_de_species) > 0) {
  ## Heatmap: DE genes (>= 1 FDR-signif. species association) x all top-15 species
  genes_hm <- unique(sig_de_species$ensembl_gene_id)

  rho_de_species <- cor_de_species |>
    filter(ensembl_gene_id %in% genes_hm) |>
    select(ensembl_gene_id, species, rho) |>
    pivot_wider(names_from = species, values_from = rho) |>
    column_to_rownames("ensembl_gene_id") |>
    as.matrix()
  rho_de_species <- rho_de_species[, top15_species, drop = FALSE]

  p_fdr_de_species <- cor_de_species |>
    filter(ensembl_gene_id %in% genes_hm) |>
    select(ensembl_gene_id, species, p_fdr) |>
    pivot_wider(names_from = species, values_from = p_fdr) |>
    column_to_rownames("ensembl_gene_id") |>
    as.matrix()
  p_fdr_de_species <- p_fdr_de_species[rownames(rho_de_species), colnames(rho_de_species), drop = FALSE]

  # Labels
  hm_symbols <- de_sig$symbol[match(rownames(rho_de_species), de_sig$ensembl_gene_id)]
  hm_species_labels <- top15_species_labels$species_label[match(colnames(rho_de_species), top15_species_labels$species)]

  # DE direction annotation
  hm_de_direction <- if_else(
    de_sig$log2FC_shr[match(rownames(rho_de_species), de_sig$ensembl_gene_id)] > 0,
    "Higher in high loss", "Higher in modest/low loss"
  )
  ha_de_direction <- rowAnnotation(
    `DE direction` = hm_de_direction,
    col = list(`DE direction` = de_cols),
    show_annotation_name = FALSE
  )

  rho_max <- max(abs(rho_de_species))

  ht_de_species <- Heatmap(
    rho_de_species,
    name = "Spearman\nrho",
    col = colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[3], "white", renoir_15[11])),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_labels = gene_label(hm_symbols),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_labels = hm_species_labels,
    column_names_rot = 45,
    column_names_gp = gpar(fontface = "italic", col = species_label_colors[hm_species_labels]),
    right_annotation = ha_de_direction,
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(stars(p_fdr_de_species[i, j]), x, y, gp = gpar(fontsize = 8))
    }
  )

  pdf(paste0("results/graphs/RNAseq/", prefix, "_de_genes_species_heatmap.pdf"),
      width = 5 + 0.35 * ncol(rho_de_species),
      height = 2.5 + 0.2 * nrow(rho_de_species))
  draw(ht_de_species)
  dev.off()

  # Individual correlation plots (supplementary)
  cor_plots <- pmap(sig_de_species, function(species, ensembl_gene_id, n, rho, p.value, p_fdr, symbol) {

      plot_data <- tissue_mb |>
        dplyr::select(
          abundance = all_of(species),
          expression = all_of(ensembl_gene_id)
        )

      species_name <- top15_species_labels$species_label[match(species, top15_species_labels$species)]

      ggplot(plot_data, aes(x = abundance, y = expression)) +
        geom_point(alpha = 0.7, size = 1.8) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7, color = renoir_15[3]) +
        labs(
          title = gene_label(symbol),
          subtitle = paste0("Spearman rho = ", round(rho, 2), " | q = ", signif(p_fdr, 2)),
          x = paste0(species_name, "\nlog10 relative abundance"),
          y = paste0(tissue_label, " gene expression (VST)")
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 8),
          axis.title.x = element_text(face = "italic"),
          panel.grid.minor = element_blank()
        )
    }
  )

  # Arrange individual cor plots
  n_cols_cor_plots <- 5
  n_rows_cor_plots <- ceiling(length(cor_plots) / n_cols_cor_plots) # grid grows with the number of significant pairs
  cor_plots_arranged <- ggarrange(plotlist = cor_plots, ncol = n_cols_cor_plots, nrow = n_rows_cor_plots)
  ggsave(
    paste0("results/graphs/RNAseq/", prefix, "_de_genes_species_significant_correlations.pdf"),
     cor_plots_arranged,
     width = 12,
     height = 2.5 * n_rows_cor_plots, # keeps each individual cor plot the same size
     limitsize = FALSE
    )
} else {
  message("No FDR-significant DE gene x top-15 species pairs: heatmap and individual plots skipped")
}

#### Metabolomics x DE genes ####
# Prep metab data
metab_df <- as(otu_table(metab), "matrix") |> # samples x metabolites
  as_tibble(rownames = "sample") |>
  mutate(id = sample |> str_remove("^BARIA_") |> str_remove("_v0$") |> as.numeric()) |>
  dplyr::select(-sample)
metab_names <- setdiff(names(metab_df), "id")

if (nrow(sig_de_species) > 0) {
  # DE genes with >=1 FDR-significant species association (from I)
  genes_I <- sig_de_species |> 
    distinct(ensembl_gene_id, symbol)

  tissue_metab_df <- tissue_mb |>
    dplyr::select(id, all_of(genes_I$ensembl_gene_id)) |>
    rename_with(~ genes_I$symbol[match(.x, genes_I$ensembl_gene_id)], starts_with("ENSG")) |>
    inner_join(metab_df, by = "id")
  nrow(tissue_metab_df)  # overlap n

  # Correlation tissue genes x metabolites
  cor_genes_metab <- cor_block(tissue_metab_df, genes_I$symbol, metab_names)
  cor_genes_metab |> 
    filter(p_fdr < 0.05) |> 
    dplyr::count(x, sort = TRUE) # hits per gene

  cor_genes_metab |> 
    filter(p_fdr < 0.05) |> 
    slice_head(n = 30)

  ## Plots ##
  # Super-pathway strip: same as in 4b_mb_metabolome_correlations.R
  pathway_colors <- c(
    "Amino Acid"                         = "#3C9189",
    "Carbohydrate"                       = "#765135",
    "Cofactors and Vitamins"             = "#385C7E",
    "Energy"                             = "#3E6330",
    "Lipid"                              = "#CCA38D",
    "Nucleotide"                         = "#8EAFCF",
    "Peptide"                            = "#96B485",
    "Xenobiotics"                        = "#814956",
    "Partially Characterized Molecules"  = "#6A4F7E"
  )
  metab_tax <- as(tax_table(metab), "matrix")

  # Genes x top metabolites heatmap (20 metabolites with the most FDR-significant associations)
  top_metabs <- cor_genes_metab |>
    group_by(y) |>
    summarize(n_sig = sum(p_fdr < 0.05), max_abs_rho = max(abs(rho)), .groups = "drop") |>
    filter(n_sig > 0) |>
    arrange(desc(n_sig), desc(max_abs_rho)) |>
    slice_head(n = 20) |>
    pull(y)

  metab_sub <- cor_genes_metab |> 
    filter(y %in% top_metabs) |> 
    dplyr::rename(x = y, y = x)  # metabolites as rows
  rho_metab <- to_matrix(metab_sub, "rho")
  fdr_metab <- to_matrix(metab_sub, "p_fdr")[rownames(rho_metab), colnames(rho_metab)]

  # Drop gene columns with no hit among those 20
  genes_keep_metab <- colnames(fdr_metab)[colSums(fdr_metab < 0.05) > 0]
  rho_metab <- rho_metab[, genes_keep_metab, drop = FALSE]
  fdr_metab <- fdr_metab[, genes_keep_metab, drop = FALSE]

  ht_metab <- Heatmap(
    rho_metab,
    name = "Spearman\nrho",
    col = colorRamp2(c(-0.4, 0, 0.4), c(renoir_15[3], "white", renoir_15[11])),
    left_annotation = rowAnnotation(
      `Super pathway` = metab_tax[rownames(rho_metab), "SUPER_PATHWAY"],
      col = list(`Super pathway` = pathway_colors),
      show_annotation_name = TRUE,
      annotation_name_gp = gpar(fontsize = 8)
    ),
    row_labels = str_trunc(rownames(rho_metab), 45),
    row_names_side = "left",
    row_dend_side = "right",
    column_labels = gene_label(colnames(rho_metab)),
    column_names_max_height = max_text_width(gene_label(colnames(rho_metab))),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) grid.text(stars(fdr_metab[i, j]), x, y, gp = gpar(fontsize = 8))
  )

  pdf(paste0("results/graphs/RNAseq/", prefix, "_genes_metabolites_heatmap.pdf"), width = 9, height = 12)
  draw(ht_metab)
  dev.off()

  ## Compact genes × metabolites heatmap (for slides)
  n_metabs_compact <- 10 # top metabolites by number of FDR-signif. gene correlations

  metabs_compact <- cor_genes_metab |>
    group_by(y) |>
    summarize(n_sig = sum(p_fdr < 0.05), max_abs_rho = max(abs(rho)), .groups = "drop") |>
    filter(n_sig > 0) |>
    arrange(desc(n_sig), desc(max_abs_rho)) |>
    slice_head(n = n_metabs_compact) |>
    pull(y)

  metab_sub_compact <- cor_genes_metab |> filter(y %in% metabs_compact) |> dplyr::rename(x = y, y = x)
  rho_metab_compact <- to_matrix(metab_sub_compact, "rho")
  fdr_metab_compact <- to_matrix(metab_sub_compact, "p_fdr")[rownames(rho_metab_compact), colnames(rho_metab_compact)]

  n_genes_compact <- 12 # genes with the most FDR-signif. correlations among the shown metabolites
  genes_keep_metab_compact <- names(sort(colSums(fdr_metab_compact < 0.05), decreasing = TRUE))[1:min(n_genes_compact, ncol(fdr_metab_compact))]
  rho_metab_compact <- rho_metab_compact[, genes_keep_metab_compact, drop = FALSE]
  fdr_metab_compact <- fdr_metab_compact[, genes_keep_metab_compact, drop = FALSE]

  dim(rho_metab_compact) # check before plotting

  ht_metab_compact <- Heatmap(
    rho_metab_compact,
    name = "Spearman\nrho",
    col = colorRamp2(c(-0.4, 0, 0.4), c(renoir_15[3], "white", renoir_15[11])),
    left_annotation = rowAnnotation(
      `Super pathway` = metab_tax[rownames(rho_metab_compact), "SUPER_PATHWAY"],
      col = list(`Super pathway` = pathway_colors),
      annotation_name_gp = gpar(fontsize = 8)
    ),
    row_labels = str_trunc(rownames(rho_metab_compact), 40),
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    show_row_dend = FALSE, show_column_dend = FALSE,
    column_labels = colnames(rho_metab_compact), # symbols only
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) grid.text(stars(fdr_metab_compact[i, j]), x, y, gp = gpar(fontsize = 8))
  )

  pdf(paste0("results/graphs/RNAseq/", prefix, "_genes_metabolites_heatmap_compact.pdf"),
      width = 5 + 0.3 * ncol(rho_metab_compact), height = 2.5 + 0.3 * nrow(rho_metab_compact))
  draw(ht_metab_compact)
  dev.off()
} else {
  message("No DE genes with species associations: genes x metabolites skipped")
}