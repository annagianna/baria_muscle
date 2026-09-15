# Liver RNA-seq analysis
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(MetBrewer)
library(biomaRt)

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
ensembl_ids <- liver_rnaseq |> 
  dplyr::select(starts_with("ENSG")) |> 
  colnames() |> 
  str_remove("\\.\\d+$")

# Connect to Ensembl BioMart database
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "description", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

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
      dplyr::select(ensembl_gene_id, hgnc_symbol, description),
    by = c("ensembl_gene_id_clean" = "ensembl_gene_id")
  ) |> 
   mutate(gene_label = if_else(
    is.na(hgnc_symbol) | hgnc_symbol == "",
    ensembl_gene_id_clean,
    hgnc_symbol
  ))

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