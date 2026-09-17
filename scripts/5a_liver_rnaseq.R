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
source("scripts/assets/functions.R")

dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

# Data
liver_rnaseq <- readRDS("data/processed_data/BARIA_Liver_RNAseq.RDS")
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_mb_v0 <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
forest_perc_change_ffmi_v4 <- read.csv("results/mlmodels/perc_change_ffmi_v4/all/forest_results_top15.csv")

# Prep data
## Top15 species
# Top 15 species by ML feature importance for %FFMI change at 1y
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
  ensembl_gene_id = liver_rnaseq |> 
    dplyr::select(starts_with("ENSG")) |> 
    colnames()
) |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$"))

# Gene annotations
gene_annotations <- grch38 |>
  filter(ensgene %in% gene_ids$ensgene) |>
  dplyr::select(ensgene, symbol, biotype) |>
  distinct()

### Liver RNA-seq x top15 species ###
# Top15 species x liver RNA seq
liver_mb <- liver_rnaseq |>
  mutate(id = as.numeric(id)) |> 
  inner_join(top15_species_log10, by = "id")

liver_expr <- liver_rnaseq |>
  dplyr::select(starts_with("ENSG")) |>
  as.matrix()

# Filter genes with non-zero expression in >= 50% of participants
gene_prevalence <- apply(liver_expr, 2, \(x) mean(x > 0))
genes_keep <- names(gene_prevalence[gene_prevalence >= 0.5])

#### Untargeted ####
### Spearman correlations top15 species ###
# Align species order
species_mat <- liver_mb |>
  dplyr::select(all_of(top15_species)) |>
  as.matrix()

gene_mat <- liver_mb |>
  dplyr::select(all_of(genes_keep)) |>
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

gene_rank <- cor_results |>
  filter(p_fdr < 0.05) |> 
  group_by(ensembl_gene_id) |>
  summarise(
    n_sig = n(),
    max_abs_rho = max(abs(rho)),
    .groups = "drop"
  ) |>
  arrange(desc(n_sig), desc(max_abs_rho))

# Select genes with the greatest overlap across top-15 species:
# 5 genes were FDR-significant for 5 species and 15 genes for 4 species,
#  giving a natural set of 20 genes associated with >=4/15 species
# Add annotations
genes_heatmap <- gene_rank |>
  filter(n_sig >= 4) |>
  mutate(ensembl_gene_id_clean = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
  left_join(
    gene_annotations |> 
      dplyr::select(ensgene, symbol),
    by = c("ensembl_gene_id_clean" = "ensgene")
  ) |>
  filter(!is.na(symbol), symbol != "") |> # Filter out unmapped/ENSG genes (but only after the correlations + FDR correction)
  mutate(gene_label = symbol)

# Heatmao data
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


# Species label colours by direction of association with 1-year %FFMI change
species_direction_colors <- c(positive = renoir_15[15],negative = renoir_15[9])
species_label_colors <- setNames(species_direction_colors[top15_species_labels$estimate_direction], top15_species_labels$species_label)

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

pdf("results/graphs/RNAseq/liver_top15_species_genes_heatmap.pdf", width = 10, height = 8)
draw(ht)
dev.off()

# More compact heatmap: only species with ≥1 significant association among these 20 genes
# Keep species with >=1 FDR-significant association among selected genes
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
draw(ht_compact)

# Save
pdf("results/graphs/RNAseq/liver_species_genes_heatmap_compact.pdf",width = 9,height = 8)
draw(ht_compact)
dev.off()

#### Targeted genes #####
# Create list of selected genes
core_genes <- list(

  hepatokines = c(
    "FGF21","SELENOP","LECT2","ANGPTL3","ANGPTL4","ANGPTL8","AHSG","FETUB",
    "FGL1","RBP4","IGF1","IGFBP1","IGFBP2","IGFBP3","INHBA","FST","FSTL3",
    "GDF15","LEAP2","GPLD1","ENHO","TSKU","SHBG","SMOC1","APOA5"),

  beta_oxidation = c(
    "CPT1A","CPT1B","CPT2","CRAT","CROT","SLC25A20","SLC22A5",
    "ACADVL","ACADL","ACADM","ACADS","ACADSB","ACAD9","ACAD8","ACAD11",
    "HADHA","HADHB","HADH","ECHS1","EHHADH","ACAA1","ACAA2","DECR1",
    "ECI1","ECI2","ETFA","ETFB","ETFDH","MLYCD","ACOX1","ACOX2",
    "HSD17B4","SCP2","PPARA"),

  mito_oxphos = c(
    "PPARGC1A","PPARGC1B","PPARD","ESRRA","ESRRG","NRF1","GABPA",
    "TFAM","TFB2M","POLG","MFN1","MFN2","OPA1","DNM1L",
    "SDHA","SDHB","SDHC","SDHD","CYCS"),

  oxidative_stress = c(
    "NFE2L2","KEAP1","NQO1","HMOX1","GCLC","GCLM","GSR","GSTP1","GSTA1",
    "GSTM1","SLC7A11","TXN","TXNRD1","PRDX1","PRDX3","PRDX5","PRDX6",
    "SOD1","SOD2","SOD3","CAT","GPX1","GPX3","GPX4","SRXN1","G6PD",
    "SESN2","FOXO3"),

  aa_bcaa_turnover = c(
    "BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","BCKDK","PPM1K",
    "HIBCH","HIBADH","IVD","MCCC1","MCCC2","AUH","HMGCL","PCCA","PCCB",
    "MMUT","ALDH6A1",
    "CPS1","OTC","ASS1","ASL","ARG1","NAGS","SLC25A15","SLC25A13",
    "GLUD1","GLS","GLS2","GLUL","GOT1","GOT2","GPT","GPT2"),

  gluconeogenesis = c(
    "PCK1","PCK2","G6PC1","SLC37A4","FBP1","FBP2","PC","MDH1","MDH2",
    "PDK1","PDK2","PDK4","FOXO1"),

  glycogen_metabolism = c(
    "GYS2","GYG1","GBE1","UGP2","PGM1","PYGL","AGL",
    "PHKA2","PHKB","PHKG2","PPP1R3B","PPP1R3C"),

  glycolysis_fructose = c(
    "GCK","HK1","HK2","GPI","PFKL","ALDOA","ALDOB","TPI1","GAPDH",
    "PGK1","PGAM1","ENO1","PKLR","PKM","PDHA1","PDHB","DLAT",
    "LDHA","LDHB","KHK","TKFC","SORD","AKR1B1","SLC2A2","SLC2A5"),

  insulin_signaling = c(
    "INSR","IGF1R","IRS1","IRS2","PIK3R1","PIK3CA","PDPK1","AKT1","AKT2",
    "GSK3B","MTOR","RPTOR","RPS6KB1","EIF4EBP1","TSC1","TSC2",
    "PRKAA1","PRKAA2","STK11","PTEN","PTPN1","MLXIPL","SREBF1")
)

target_genes <- enframe(core_genes, name = "category", value = "symbol") |> 
  unnest(symbol) |> 
  left_join(
    grch38 |> 
      dplyr::select(symbol, ensgene, description) |> 
      distinct(),
    by = "symbol"
  ) |> 
  inner_join(gene_ids, by = "ensgene") |> 
  filter(ensembl_gene_id %in% genes_keep) # adds prevalence filter (non-zero expression in >= 50% of participants)

# Spearman correlations: top-15 species x targeted genes
target_gene_mat <- liver_mb |>
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

# Build matrices
rho_targeted <- cor_targeted |>
  filter(ensembl_gene_id %in% target_genes_heatmap$ensembl_gene_id) |>
  dplyr::select(species, description, rho) |>
  pivot_wider(names_from = species, values_from = rho) |>
  column_to_rownames("description") |>
  as.matrix()

fdr_targeted <- cor_targeted |>
  filter(ensembl_gene_id %in% target_genes_heatmap$ensembl_gene_id) |>
  dplyr::select(species, description, p_fdr) |>
  pivot_wider(names_from = species, values_from = p_fdr) |>
  column_to_rownames("description") |>
  as.matrix()

# Gene categories
gene_categories <- target_genes_heatmap$category[
  match(rownames(rho_targeted), target_genes_heatmap$description)
]

category_cols <- setNames(
  renoir_15[c(1, 5, 7, 10, 13, 2, 6, 12, 14)],
  names(core_genes)
)

category_anno <- rowAnnotation(
  Category = gene_categories,
  col = list(Category = category_cols),
  show_annotation_name = FALSE
)

# Species labels
species_labels <- top15_species_labels$species_label[
  match(colnames(rho_targeted), top15_species_labels$species)
]

# Significance asterisks
stars_targeted <- ifelse(
  fdr_targeted < 0.001, "***",
  ifelse(fdr_targeted < 0.01, "**",
         ifelse(fdr_targeted < 0.05, "*", ""))
)

# Full targeted heatmap
# Full targeted heatmap
ht_targeted <- Heatmap(
  rho_targeted,
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
  column_labels = species_labels,
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = 9,
    fontface = "italic",
    col = species_label_colors[species_labels]
  ),
  row_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(stars_targeted[i, j], x, y, gp = gpar(fontsize = 10))
  }
)

# Save
pdf(
  "results/graphs/RNAseq/liver_top15_species_targeted_genes_heatmap.pdf",
  width = 13,
  height = 8
)
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

    plot_data <- liver_mb |>
      dplyr::select(
        abundance = all_of(species),
        expression = all_of(ensembl_gene_id)
      )

    species_name <- top15_species_labels$species_label[
      match(species, top15_species_labels$species)
    ]

    ggplot(plot_data, aes(x = abundance, y = expression)) +
      geom_point(alpha = 0.7, size = 1.8) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = renoir_15[3]) +
      labs(
        title = description,
        subtitle = paste0(
          "Spearman rho = ", round(rho, 2),
          " | FDR = ", signif(p_fdr, 2)
        ),
        x = paste0(species_name, "\nlog10 relative abundance"),
        y = "Liver gene expression"
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
cor_plots_arranged <- ggarrange(plotlist = cor_plots, ncol = 3, nrow = 4)
ggsave("results/graphs/RNAseq/liver_targeted_significant_correlations.pdf", cor_plots_arranged,  width = 12, height = 12)

## Compact targeted heatmap ##
# Species with >=1 FDR-significant association x liver gene (targeted)
species_keep_targeted <- colnames(fdr_targeted)[
  apply(fdr_targeted < 0.05, 2, any)
]

rho_targeted_compact <- rho_targeted[, species_keep_targeted, drop = FALSE]
fdr_targeted_compact <- fdr_targeted[, species_keep_targeted, drop = FALSE]

# Species labels
species_labels_compact <- top15_species_labels$species_label[
  match(colnames(rho_targeted_compact), top15_species_labels$species)
]

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
pdf(
  "results/graphs/RNAseq/liver_species_targeted_genes_heatmap_compact.pdf",
  width = 11,
  height = 8
)
draw(ht_targeted_compact, newpage = FALSE)
dev.off()