# Visceral adipose tissue RNA-seq
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
source("scripts/assets/functions.R")
dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

# Data
vfat_rnaseq <- readRDS("data/processed_data/BARIA_vFat_RNAseq_vst.RDS") # VST normalized
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_mb_v0 <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
forest_perc_change_ffmi_v4 <- read.csv("results/mlmodels/perc_change_ffmi_v4/all/forest_results_top15.csv")

# Prep data
## Top 15 species by ML feature importance for %FFMI change at 1y
fi <- get_feature_importance("results/mlmodels/perc_change_ffmi_v4/all", "perc_change_ffmi_v4_all", "reg")
top15_species <- top_features(fi, n = 15)$FeatName
top15_species_labels <- tibble(species = top15_species, species_label = species_label(top15_species)) |>
  left_join(
    forest_perc_change_ffmi_v4 |>
      dplyr::select(species, estimate, p_fdr) |>
      mutate(estimate_direction = if_else(estimate > 0, "positive", "negative")),
    by = "species"
  )

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
  ensembl_gene_id = vfat_rnaseq |>
    dplyr::select(starts_with("ENSG")) |>
    colnames()
) |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$"))

# Gene annotations
gene_annotations <- grch38 |>
  filter(ensgene %in% gene_ids$ensgene) |>
  dplyr::select(ensgene, symbol, biotype) |>
  distinct()

### Visceral adipose tissue RNA-seq x top15 species ###
# Top15 species x vFat RNA seq
vfat_mb <- vfat_rnaseq |>
  mutate(id = as.numeric(id)) |>
  inner_join(top15_species_log10, by = "id")

# QC: sequencing depth (DESeq2 size factor) vs. species abundance (no expected association)
size_factor_qc <- map_dfr(top15_species,
  ~ broom::tidy(cor.test(vfat_mb$size_factor, vfat_mb[[.x]], method = "spearman", exact = FALSE)
) |> mutate(species = .x)) |>
  dplyr::select(species, rho = estimate, p.value)
size_factor_qc

# Genes in the normalised file (already expression-filtered in cleaning script 0c)
vfat_genes <- vfat_rnaseq |>
  dplyr::select(starts_with("ENSG")) |>
  colnames()

# Species label colours by direction of association with 1-year %FFMI change
species_direction_colors <- c(positive = renoir_15[15], negative = renoir_15[9])
species_label_colors <- setNames(species_direction_colors[top15_species_labels$estimate_direction], top15_species_labels$species_label)

#### Untargeted ####
### Spearman correlations top15 species ###
# Align species order
species_mat <- vfat_mb |>
  dplyr::select(all_of(top15_species)) |>
  as.matrix()

gene_mat <- vfat_mb |>
  dplyr::select(all_of(vfat_genes)) |>
  as.matrix()

cor_results <- expand_grid(
  species = top15_species,
  ensembl_gene_id = colnames(gene_mat)
) |>
  mutate(
    test = map2(species, ensembl_gene_id, ~ cor.test(
      species_mat[, .x],
      gene_mat[, .y],
      method = "spearman",
      exact = FALSE
    )),
    rho = map_dbl(test, "estimate"),
    p.value = map_dbl(test, "p.value")
  ) |>
  dplyr::select(-test) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH"))

# Top 25 (annotated) associated genes by |rho| (among FDR-signif associations)
genes_heatmap <- cor_results |>
  group_by(ensembl_gene_id) |>
  summarize(
    min_fdr = min(p_fdr),
    n_sig = sum(p_fdr < 0.05),
    max_abs_rho = max(abs(rho)),
    .groups = "drop"
  ) |>
  filter(min_fdr < 0.05) |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
  left_join(
    gene_annotations |> distinct(ensgene, .keep_all = TRUE) |> dplyr::select(ensgene, symbol),
    by = "ensgene"
  ) |>
  filter(!is.na(symbol), symbol != "") |> # annotated genes only, selected after FDR correction
  arrange(desc(max_abs_rho)) |>
  slice_head(n = 25) |>
  mutate(gene_label = symbol)

if (nrow(genes_heatmap) == 0) {

  message(
    "No gene reached FDR < 0.05 for any species in the untargeted visceral adipose scan; ",
    "skipping untargeted heatmaps. Closest associations (by raw p-value):"
  )
  cor_results |>
    arrange(p.value) |>
    head(20) |>
    print(n = Inf)

} else {

  # Heatmap data
  heatmap_data <- cor_results |>
    filter(ensembl_gene_id %in% genes_heatmap$ensembl_gene_id)

  rho_mat <- heatmap_data |>
    dplyr::select(species, ensembl_gene_id, rho) |>
    pivot_wider(names_from = ensembl_gene_id, values_from = rho) |>
    column_to_rownames("species") |>
    as.matrix()

  fdr_mat <- heatmap_data |>
    dplyr::select(species, ensembl_gene_id, p_fdr) |>
    pivot_wider(names_from = ensembl_gene_id, values_from = p_fdr) |>
    column_to_rownames("species") |>
    as.matrix()

  # Replace gene labels
  colnames(rho_mat) <- genes_heatmap$gene_label[match(colnames(rho_mat), genes_heatmap$ensembl_gene_id)]
  colnames(fdr_mat) <- colnames(rho_mat)

  # Add species labels
  rownames(rho_mat) <- top15_species_labels$species_label[match(rownames(rho_mat), top15_species_labels$species)]
  rownames(fdr_mat) <- rownames(rho_mat)

  # Diverging colour scale: negative = Renoir[3], positive = Renoir[11]
  rho_max <- max(abs(rho_mat), na.rm = TRUE)
  col_fun <- colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[3], "white", renoir_15[11]))

  ## Complex Heatmap ##
  rho_mat_t <- t(rho_mat)
  fdr_mat_t <- t(fdr_mat)

  ht <- Heatmap(
    rho_mat_t,
    name = "Spearman\nrho",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_dend_side = "right",
    row_names_side = "left",
    column_names_rot = 45,
    column_names_gp = gpar(
      fontsize = 10,
      fontface = "italic",
      col = species_label_colors[colnames(rho_mat_t)]
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (fdr_mat_t[i, j] < 0.001) {
        grid.text("***", x, y)
      } else if (fdr_mat_t[i, j] < 0.01) {
        grid.text("**", x, y)
      } else if (fdr_mat_t[i, j] < 0.05) {
        grid.text("*", x, y)
      }
    }
  )

  pdf("results/graphs/RNAseq/vfat_top15_species_genes_heatmap.pdf", width = 10, height = 8)
  draw(ht)
  dev.off()

  # More compact heatmap: only species with >=1 signif association among the selected genes
  species_keep <- rownames(fdr_mat)[apply(fdr_mat < 0.05, 1, any)]
  rho_mat_compact <- rho_mat[species_keep, , drop = FALSE]
  fdr_mat_compact <- fdr_mat[species_keep, , drop = FALSE]
  rho_mat_compact <- t(rho_mat_compact)
  fdr_mat_compact <- t(fdr_mat_compact)

  # Compact heatmap: only species with >=1 significant association
  ht_compact <- Heatmap(
    rho_mat_compact,
    name = "Spearman\nrho",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_dend_side = "right",
    row_names_side = "left",
    column_names_rot = 45,
    column_names_gp = gpar(
      fontsize = 10,
      fontface = "italic",
      col = species_label_colors[colnames(rho_mat_compact)]
    ),
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (fdr_mat_compact[i, j] < 0.001) {
        grid.text("***", x, y)
      } else if (fdr_mat_compact[i, j] < 0.01) {
        grid.text("**", x, y)
      } else if (fdr_mat_compact[i, j] < 0.05) {
        grid.text("*", x, y)
      }
    }
  )

  # Save
  pdf("results/graphs/RNAseq/vfat_species_genes_heatmap_compact.pdf", width = 9, height = 8)
  draw(ht_compact)
  dev.off()

}

#### Targeted genes #####
# Create list of selected genes
core_genes <- list(

  adipokines = c(
    "ADIPOQ","LEP","RETN","NAMPT","RBP4","LCN2","RARRES2","SERPINA12",
    "ITLN1","SFRP5","ANGPTL4","ANGPTL8","WISP1","GDF15","FGF21","FABP4","CFD"),

  adipogenesis = c(
    "PPARG","CEBPA","CEBPB","CEBPD","KLF15","SREBF1","DLK1","WNT10B",
    "PLIN1","CIDEC","LPL"),

  lipolysis_lipid_storage = c(
    "PLIN2","PLIN4","PLIN5","LIPE","PNPLA2","MGLL","ABHD5","PDE3B",
    "PRKAR2B","AQP7"),

  lipogenesis_fa_uptake = c(
    "FASN","ACACA","ACACB","SCD","ELOVL6","MLXIPL","CD36",
    "DGAT1","DGAT2","GPAM","AGPAT2"),

  mito_oxphos = c(
    "PPARGC1A","PPARGC1B","NRF1","GABPA","ESRRA","ESRRG","TFAM","TFB2M","POLG",
    "MFN1","MFN2","OPA1","DNM1L","SDHA","SDHB","SDHC","SDHD","CYCS"),

  insulin_signaling = c(
    "INSR","IRS1","IRS2","PIK3R1","PIK3CA","PDPK1","AKT1","AKT2",
    "GSK3B","MTOR","RPTOR","RPS6KB1","EIF4EBP1","TSC1","TSC2",
    "PRKAA1","PRKAA2","STK11","PTEN","PTPN1","SLC2A4","TBC1D4","SORBS1"),

  inflammation_macrophage = c(
    "TNF","IL6","IL1B","IL10","CCL2","CD68","ITGAX","ADGRE1","NLRP3",
    "TLR4","SAA1","SAA2","MMP9","LBP","NFKB1","CRP"),

  fibrosis_ecm = c(
    "COL1A1","COL1A2","COL3A1","COL6A1","COL6A2","COL6A3","TGFB1",
    "MMP2","TIMP1","SPARC","FN1"),

  oxidative_stress = c(
    "NFE2L2","KEAP1","NQO1","HMOX1","GCLC","GCLM","GSR","SOD1","SOD2",
    "SOD3","CAT","GPX1","GPX3","GPX4","TXN","TXNRD1","PRDX1","PRDX3",
    "SESN2","FOXO3"),

  aa_bcaa_turnover = c(
    "BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","BCKDK","PPM1K",
    "GLUD1","GLUL","GOT1","GOT2"),

  glycolysis = c(
    "HK2","GPI","PFKL","ALDOA","GAPDH","PGK1","PGAM1","ENO1","PKM",
    "PDHA1","PDHB","LDHA","LDHB","PDK4","SLC2A1")
)
stopifnot(!anyDuplicated(unlist(core_genes))) # each gene in exactly one category

target_genes <- enframe(core_genes, name = "category", value = "symbol") |>
  unnest(symbol) |>
  left_join(
    grch38 |>
      dplyr::select(symbol, ensgene, description) |>
      distinct(),
    by = "symbol"
  ) |>
  inner_join(gene_ids, by = "ensgene")

# Spearman correlations: top-15 species x targeted genes
target_gene_mat <- vfat_mb |>
  dplyr::select(all_of(target_genes$ensembl_gene_id)) |>
  as.matrix()

cor_targeted <- expand_grid(
  species = top15_species,
  ensembl_gene_id = colnames(target_gene_mat)
) |>
  mutate(
    test = map2(species, ensembl_gene_id, ~ cor.test(
      species_mat[, .x],
      target_gene_mat[, .y],
      method = "spearman",
      exact = FALSE
    )),
    rho = map_dbl(test, "estimate"),
    p.value = map_dbl(test, "p.value")
  ) |>
  dplyr::select(-test) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH")) |> # global FDR p for all selected genes
  left_join(
    target_genes |>
      dplyr::select(ensembl_gene_id, symbol, category, description),
    by = "ensembl_gene_id"
  )

# Genes with >=1 FDR-significant association
target_genes_heatmap <- cor_targeted |>
  filter(p_fdr < 0.05) |>
  distinct(ensembl_gene_id, symbol, description, category)

if (nrow(target_genes_heatmap) == 0) {

  message(
    "No targeted gene-species associations reached FDR < 0.05 in visceral adipose tissue; ",
    "skipping targeted heatmaps/plots. Closest associations (by raw p-value):"
  )
  cor_targeted |>
    arrange(p.value) |>
    dplyr::select(category, symbol, species, rho, p.value, p_fdr, description) |>
    head(20) |>
    print(n = Inf)

} else {

  # Build matrices
  rho_targeted <- cor_targeted |>
    filter(ensembl_gene_id %in% target_genes_heatmap$ensembl_gene_id) |>
    dplyr::select(species, symbol, rho) |>
    pivot_wider(names_from = species, values_from = rho) |>
    column_to_rownames("symbol") |>
    as.matrix()

  fdr_targeted <- cor_targeted |>
    filter(ensembl_gene_id %in% target_genes_heatmap$ensembl_gene_id) |>
    dplyr::select(species, symbol, p_fdr) |>
    pivot_wider(names_from = species, values_from = p_fdr) |>
    column_to_rownames("symbol") |>
    as.matrix()

  # Gene categories
  gene_categories <- target_genes_heatmap$category[match(rownames(rho_targeted), target_genes_heatmap$symbol)]
  category_cols <- setNames(met.brewer("Renoir", n = length(core_genes)), names(core_genes))
  category_anno <- rowAnnotation(Category = gene_categories, col = list(Category = category_cols), show_annotation_name = FALSE)

  # Species labels
  species_labels <- top15_species_labels$species_label[match(colnames(rho_targeted), top15_species_labels$species)]

  # Significance asterisks
  stars_targeted <- ifelse(
    fdr_targeted < 0.001, "***",
    ifelse(fdr_targeted < 0.01, "**",
           ifelse(fdr_targeted < 0.05, "*", ""))
  )

  # Full targeted heatmap
  ht_targeted <- Heatmap(
    rho_targeted,
    name = "Spearman\nrho",
    col = colorRamp2(c(-0.3, 0, 0.3), c(renoir_15[3], "white", renoir_15[11])),
    left_annotation = category_anno,
    row_split = gene_categories,
    row_title = NULL,
    row_names_side = "left",
    row_dend_side = "right",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_labels = species_labels,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9, fontface = "italic", col = species_label_colors[species_labels]),
    row_names_gp = gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, width, height, fill) {grid.text(stars_targeted[i, j], x, y, gp = gpar(fontsize = 10))}
  )

  # Save
  pdf("results/graphs/RNAseq/vfat_top15_species_targeted_genes_heatmap.pdf", width = 13, height = 8)
  draw(ht_targeted, newpage = FALSE)
  dev.off()

  #### Individual correlation plots ####
  # FDR-significant targeted species x gene associations
  sig_targeted <- cor_targeted |>
    filter(p_fdr < 0.05) |>
    arrange(p_fdr)

  cor_plots <- pmap(
    sig_targeted,
    function(species, ensembl_gene_id, rho, p.value, p_fdr, symbol, category, description) {

      plot_data <- vfat_mb |>
        dplyr::select(
          abundance = all_of(species),
          expression = all_of(ensembl_gene_id)
        )

      species_name <- top15_species_labels$species_label[match(species, top15_species_labels$species)]

      ggplot(plot_data, aes(x = abundance, y = expression)) +
        geom_point(alpha = 0.7, size = 1.8) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7, color = renoir_15[3]) +
        labs(
          title = description,
          subtitle = paste0(
            "Spearman rho = ", round(rho, 2),
            " | FDR = ", signif(p_fdr, 2)
          ),
          x = paste0(species_name, "\nlog10 relative abundance"),
          y = "Visceral adipose gene expression (VST)"
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
  cor_plots_arranged <- ggarrange(plotlist = cor_plots, ncol = 5, nrow = ceiling(length(cor_plots) / 5))
  ggsave("results/graphs/RNAseq/vfat_targeted_significant_correlations.pdf", cor_plots_arranged, width = 12, height = 10)

  ## Compact targeted heatmap ##
  # Species with >=1 FDR-significant association x visceral adipose gene (targeted)
  species_keep_targeted <- colnames(fdr_targeted)[
    apply(fdr_targeted < 0.05, 2, any)
  ]

  rho_targeted_compact <- rho_targeted[, species_keep_targeted, drop = FALSE]
  fdr_targeted_compact <- fdr_targeted[, species_keep_targeted, drop = FALSE]

  # Species labels
  species_labels_compact <- top15_species_labels$species_label[match(colnames(rho_targeted_compact), top15_species_labels$species)]

  # Significance asterisks
  stars_targeted_compact <- ifelse(
    fdr_targeted_compact < 0.001, "***",
    ifelse(fdr_targeted_compact < 0.01, "**",
           ifelse(fdr_targeted_compact < 0.05, "*", ""))
  )

  # Heatmap
  ht_targeted_compact <- Heatmap(
    rho_targeted_compact,
    name = "Spearman\nrho",
    col = colorRamp2(
      c(-0.3, 0, 0.3),
      c(renoir_15[3], "white", renoir_15[11])
    ),
    left_annotation = category_anno,
    row_split = gene_categories,
    row_title = NULL,
    row_names_side = "left",
    row_dend_side = "right",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_labels = species_labels_compact,
    column_names_rot = 45,
    column_names_gp = gpar(
      fontsize = 10,
      fontface = "italic",
      col = species_label_colors[species_labels_compact]
    ),
    row_names_gp = gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        stars_targeted_compact[i, j],
        x, y,
        gp = gpar(fontsize = 10)
      )
    }
  )

  # Save
  pdf("results/graphs/RNAseq/vfat_species_targeted_genes_heatmap_compact.pdf", width = 7, height = 8)
  draw(ht_targeted_compact, newpage = FALSE)
  dev.off()

  sig_targeted |>
    left_join(top15_species_labels |> dplyr::select(species, species_label, estimate_direction), by = "species") |>
    dplyr::select(category, symbol, species_label, rho, p_fdr, estimate_direction, description) |>
    arrange(species_label, category, p_fdr) |>
    print(n = Inf)

}


nrow(vfat_mb)
ncol(gene_mat)
