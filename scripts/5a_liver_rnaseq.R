# Liver RNA-seq
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
source("scripts/assets/functions.R")
dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

# Bioconductor packages (DESeq2/S4Vectors/matrixStats) mask some dplyr verbs; always use dplyr's
select <- dplyr::select
rename <- dplyr::rename
count  <- dplyr::count
filter <- dplyr::filter

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
liver_rnaseq <- readRDS("data/processed_data/BARIA_Liver_RNAseq.RDS") # raw counts for DESeq2
liver_rnaseq_vst <- readRDS("data/processed_data/BARIA_Liver_RNAseq_vst.RDS") # VST normalized
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
  ensembl_gene_id = liver_rnaseq_vst |> 
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
         both = paste0(symbols, " – ", desc),
         symbol = symbols,
         description = desc)
}

# Diverging colour scale: negative = Renoir[3], positive = Renoir[11]
rho_max <- max(abs(rho_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[3], "white", renoir_15[11]))

# Species label colours by direction of association with high %FFMI loss at 1y (logistic OR)
species_direction_colors <- c(`Higher in modest/low loss` = renoir_15[15], `Higher in high loss` = renoir_15[9])
species_label_colors <- setNames(species_direction_colors[top15_species_labels$estimate_direction], top15_species_labels$species_label)

#### DE: high vs modest/low %FFMI loss (1y) ####
# MA plot
ma_df <- bind_rows(
  plotMA(res_cont, returnData = TRUE) |> mutate(type = "Before shrinkage (MLE)"),
  plotMA(res_cont_shr, returnData = TRUE) |> mutate(type = "After shrinkage (apeglm)")
) |>
  mutate(type = factor(type, levels = c("Before shrinkage (MLE)", "After shrinkage (apeglm)")))

ma_plot_cont <- ggplot(ma_df, aes(mean, lfc, colour = isDE)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_x_log10() +
  scale_colour_manual(
    values = c(`FALSE` = "grey75", `TRUE` = renoir_15[11]),
    labels = c("not significant", "FDR < 0.05 (all genes)"),
    name = NULL
  ) +
  facet_wrap(~ type) +
  labs(x = "Mean of normalised counts", y = "log2FC per 6% FFMI loss") +
  theme_minimal_custom()
ggsave("results/graphs/RNAseq/liver_de_cont_MA_before_after_shrinkage.pdf", ma_plot_cont, width = 10, height = 5)

## Forest plot (continuous DE)
# Data
forest_de_cont_df <- de_cont_targeted |> 
  filter(padj_targeted < 0.05) |> 
  mutate(
    direction = if_else(log2FC_shr > 0, "higher with more FFMI loss", "lower with more FFMI loss"),
    lower = log2FC_shr - 1.96 * lfcSE_shr, # approx. 95% interval
    upper = log2FC_shr + 1.96 * lfcSE_shr,
    label = fct_reorder(gene_label(symbol), log2FC_shr)
  )

# Plot
forest_de_cont_plot <- forest_de_cont_df |> 
  ggplot(aes(x = log2FC_shr, y = label, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = de_cols, name = NULL) +
  labs(
    title = "Differences in gene expression",
    x = paste0("log2FC per ", ffmi_unit, "% FFMI loss (approx. 95% CI)"),
    y = NULL
  ) +
  theme_minimal_custom()
ggsave("results/graphs/RNAseq/forest_de_cont_plot.pdf", forest_de_cont_plot, width = 9, height = 1.5 + 0.4 * nrow(forest_de_cont_df))

# Sensitivity: DE adj. for %TWL
dds_cont_twl <- DESeqDataSetFromMatrix(
  count_data_cont,
  coldata_cont, 
  design = ~ sex + age_v0_z + fmi_v0_z + twl_z + ffmi_loss_u
)
dds_cont_twl <- estimateSizeFactors(dds_cont_twl)
dds_cont_twl <- dds_cont_twl[rowMeans(counts(dds_cont_twl, normalized = TRUE) >= 10) >= 0.5, ]
dds_cont_twl <- DESeq(dds_cont_twl)
res_cont_twl <- results(dds_cont_twl, name = "ffmi_loss_u")

de_cont_targeted_twl <- as.data.frame(res_cont_twl) |>
  rownames_to_column("ensembl_gene_id") |>
  select(
    ensembl_gene_id,
    log2FC_mle_twl = log2FoldChange,
    pvalue_twl = pvalue
  ) |>
  inner_join(target_genes |> 
    distinct(ensembl_gene_id, symbol), 
  by = "ensembl_gene_id") |>
  mutate(padj_targeted_twl = p.adjust(pvalue_twl, method = "BH"))

# Compare unadj. vs. %TWL-adj. DE results
de_twl_compare <- de_signif |>
  select(ensembl_gene_id, symbol, log2FC_mle, padj_targeted) |>
  left_join(de_cont_targeted_twl |> 
    select(ensembl_gene_id, log2FC_mle_twl, padj_targeted_twl), 
  by = "ensembl_gene_id") |>
  mutate(
    same_direction = sign(log2FC_mle) == sign(log2FC_mle_twl),
    still_sig = padj_targeted_twl < 0.05
  )
de_twl_compare # all TRUE

#### Liver RNA-seq x top15 species ####
# Top15 species x liver RNA seq
liver_mb <- liver_rnaseq_vst |>
  mutate(id = as.numeric(as.character(id))) |> 
  inner_join(top15_species_log10, by = "id")

# QC: sequencing depth (DESeq2 size factor) vs. species abundance (no expected association)
size_factor_qc <- map_dfr(top15_species,
  ~ broom::tidy(cor.test(liver_mb$size_factor, liver_mb[[.x]], method = "spearman", exact = FALSE)
) |> mutate(species = .x)) |>
  dplyr::select(species, rho = estimate, p.value)
size_factor_qc

#### Individual correlation plots ####
# FDR-significant targeted species x gene associations
sig_targeted <- cor_targeted |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

cor_plots <- pmap(
  sig_targeted,

  function(species, ensembl_gene_id, rho, p.value, p_fdr, symbol, category, description) {

    plot_data <- liver_mb |>
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
        subtitle = paste0(
          "Spearman rho = ", round(rho, 2),
          " | q = ", signif(p_fdr, 2)
        ),
        x = paste0(species_name, "\nlog10 relative abundance"),
        y = "Liver gene expression (VST)"
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
  "results/graphs/RNAseq/liver_targeted_significant_correlations.pdf",
   cor_plots_arranged,
   width = 12,
   height = 2.5 * n_rows_cor_plots, # keeps each individual cor plot the same size
   limitsize = FALSE
  )

#### Metabolomics ####
# Prep metab data
metab_df <- as(otu_table(metab), "matrix") |> # samples x metabolites
  as_tibble(rownames = "sample") |>
  mutate(id = sample |> str_remove("^BARIA_") |> str_remove("_v0$") |> as.numeric()) |>
  dplyr::select(-sample)
metab_names <- setdiff(names(metab_df), "id")

# Correlation liver genes x metabolites
cor_genes_metab <- cor_block(liver_metab_df, genes_keep_targeted, metab_names)
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

# Figure 3: liver genes × top metabolites heatmap (shows 20 metabolites with the strongest signif associations)
top_metabs <- cor_genes_metab |>
  group_by(y) |>
  summarize(n_sig = sum(p_fdr < 0.05), max_abs_rho = max(abs(rho)), .groups = "drop") |>
  filter(n_sig > 0) |>
  arrange(desc(n_sig), desc(max_abs_rho)) |>
  slice_head(n = 20) |>
  pull(y)

metab_sub <- cor_genes_metab |> filter(y %in% top_metabs) |> dplyr::rename(x = y, y = x)  # metabolites as rows
rho_metab <- to_matrix(metab_sub, "rho")
fdr_metab <- to_matrix(metab_sub, "p_fdr")[rownames(rho_metab), colnames(rho_metab)]

# Drop gene columns with no hit among those 20
genes_keep_metab <- colnames(fdr_metab)[colSums(fdr_metab < 0.05) > 0]
rho_metab <- rho_metab[, genes_keep_metab, drop = FALSE]
fdr_metab <- fdr_metab[, genes_keep_metab, drop = FALSE]
metab_gene_categories <- target_genes$category[match(colnames(rho_metab), target_genes$symbol)]

ht_metab <- Heatmap(
  rho_metab,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.4, 0, 0.4), c(renoir_15[3], "white", renoir_15[11])),
  top_annotation = HeatmapAnnotation(Category = metab_gene_categories, col = list(Category = category_cols), show_annotation_name = FALSE),
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

pdf("results/graphs/RNAseq/liver_genes_metabolites_heatmap.pdf", width = 9, height = 12)
draw(ht_metab)
dev.off()

## Compact genes × metabolites heatmap (for slides)
n_metabs_compact <- 10 # top metabolites by number of FDR-signif. gene correlations
min_hits_gene_metab <- 2 # genes with >= 2 FDR-signif. correlations among the shown metabolites

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
metab_gene_categories_compact <- target_genes$category[match(colnames(rho_metab_compact), target_genes$symbol)]

dim(rho_metab_compact) # check before plotting

ht_metab_compact <- Heatmap(
  rho_metab_compact,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.4, 0, 0.4), c(renoir_15[3], "white", renoir_15[11])),
  top_annotation = HeatmapAnnotation(Category = metab_gene_categories_compact, col = list(Category = category_cols), show_annotation_name = FALSE),
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

pdf("results/graphs/RNAseq/liver_genes_metabolites_heatmap_compact.pdf",
    width = 5 + 0.3 * ncol(rho_metab_compact), height = 2.5 + 0.3 * nrow(rho_metab_compact))
draw(ht_metab_compact)
dev.off()
