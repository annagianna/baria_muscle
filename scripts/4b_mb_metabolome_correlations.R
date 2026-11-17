# Correlations: top-15 ML species and top-15 ML metabolites (1-year %FFMI loss group)
# Barbara Verhaar, Anna Giannakogeorgou

# Packages
library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(MetBrewer)
library(grid)
source("scripts/assets/functions.R")

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)

# Data
baria_mb_baseline <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
metab <- readRDS("data/processed_data/BARIA_metabolon_clean.RDS") # baseline metabolomics

out_dir <- "results/graphs/mb_metabolome_correlations/ml_selected"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#### ML-selected features from dichotomous models ####
# Top-15 species
top15_species_labels <- get_top_species()
top15_species <- top15_species_labels$species

# Top-15 metabolites
top15_metabolites_info <- get_top_metabolites()
top15_metabolites <- top15_metabolites_info$metabolite

#### Baseline microbiome ####
species_mat <- t(as(baria_mb_baseline@otu_table, "matrix"))
stopifnot(all(top15_species %in% colnames(species_mat)))

species_wide <- as.data.frame(species_mat[, top15_species, drop = FALSE]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v0$")) |>
  select(-sampleid)

# log10-transform (half-minimum pseudocount)
species_log10 <- species_wide |> select(id, all_of(top15_species))
species_log10[top15_species] <- map(species_log10[top15_species], ~ {
  pseudo <- min(.x[.x > 0]) / 2
  log10(.x + pseudo)
})

#### Baseline metabolomics ####
metab_mat <- as(otu_table(metab), "matrix")
metab_meta <- as(sample_data(metab), "data.frame")
rownames(metab_mat) <- metab_meta$id

# Check and retain ML-selected metabolites
missing_metabs <- setdiff(top15_metabolites, colnames(metab_mat))
if (length(missing_metabs) > 0) {
  print(missing_metabs)
  stop("ML metabolites missing from Metabolon data")
}
metab_mat <- metab_mat[, top15_metabolites, drop = FALSE]

#### Match participants ####
shared_ids <- intersect(species_log10$id, rownames(metab_mat))
cat(length(shared_ids), "participants with baseline microbiome and metabolome data\n")

species_log10_mat <- species_log10 |>
  filter(id %in% shared_ids) |>
  column_to_rownames("id")
species_log10_mat <- as.matrix(species_log10_mat[shared_ids, , drop = FALSE])
metab_mat <- metab_mat[shared_ids, , drop = FALSE]

stopifnot(
  identical(rownames(species_log10_mat), rownames(metab_mat)),
  !anyNA(species_log10_mat),
  !anyNA(metab_mat)
)

#### Species x metabolite correlations ####
cor_results <- expand_grid(species = top15_species, metabolite = top15_metabolites) |>
  mutate(
    test = map2(species, metabolite, ~ cor.test(species_log10_mat[, .x], metab_mat[, .y], method = "spearman", exact = FALSE)),
    rho = map_dbl(test, "estimate"),
    p.value = map_dbl(test, "p.value")
  ) |>
  select(-test) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH"))

write.csv(cor_results, file.path(out_dir, "ml_species_metabolite_correlations.csv"), row.names = FALSE)

sig <- cor_results |> filter(p_fdr < 0.05)
cat(nrow(sig), "/", nrow(cor_results), "species-metabolite correlations FDR < 0.05\n")

# Summary of significant associations
sig_summary <- sig |>
  group_by(species) |>
  summarise(n_sig = n(), n_positive = sum(rho > 0), n_negative = sum(rho < 0), .groups = "drop") |>
  arrange(desc(n_sig))

#### Scatter plots: significant associations ####
species_labels <- setNames(top15_species_labels$species_label, top15_species_labels$species)
format_fdr <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))

for (species_id in unique(sig$species)) {
  sig_metabs <- sig |> filter(species == species_id) |> arrange(p_fdr) |> pull(metabolite)
  ncol_used <- min(5, length(sig_metabs))
  nrow_used <- ceiling(length(sig_metabs) / ncol_used)

  df <- tibble(id = rownames(species_log10_mat), abundance = species_log10_mat[, species_id]) |>
    left_join(
      metab_mat[, sig_metabs, drop = FALSE] |>
        as.data.frame() |>
        rownames_to_column("id") |>
        pivot_longer(-id, names_to = "metabolite", values_to = "level"),
      by = "id"
    ) |>
    left_join(sig |> filter(species == species_id) |> select(metabolite, rho, p_fdr), by = "metabolite") |>
    mutate(facet_label = str_wrap(metabolite_label(metabolite), 22))

  ann <- df |>
    distinct(facet_label, rho, p_fdr) |>
    mutate(label = sprintf("rho = %.2f\nq %s", rho,
                           if_else(p_fdr < 0.001, format_fdr(p_fdr), paste0("= ", format_fdr(p_fdr)))))

  p <- ggplot(df, aes(x = abundance, y = level)) +
    geom_point(alpha = 0.5, size = 1, color = renoir_15[14]) +
    geom_smooth(method = "lm", se = FALSE, color = renoir_15[6], linewidth = 0.6) +
    geom_text(data = ann, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
              inherit.aes = FALSE, size = 2.5) +
    facet_wrap(~facet_label, scales = "free", ncol = ncol_used) +
    labs(title = species_labels[[species_id]], x = "log10 relative abundance (baseline)",
         y = "Metabolite level (z-scored log10)") +
    theme_minimal_custom() +
    theme(strip.text = element_text(size = 8), plot.title = element_text(size = 10))

  fname <- str_replace_all(species_labels[[species_id]], " ", "_")
  ggsave(file.path(out_dir, sprintf("correlation_scatter_%s.pdf", fname)), p,
         width = 2.6 * ncol_used + 0.5, height = 2.6 * nrow_used + 0.6)
}

#### ComplexHeatmap: all top-15 species x top-15 metabolites ####
rho_max <- max(abs(cor_results$rho), na.rm = TRUE)
col_fun <- colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[14], "white", renoir_15[6]))

# Species and metabolite label colours by logistic regression association with FFMI loss
direction_colors <- c(
  "Higher in modest/low loss" = renoir_15[15],
  "Higher in high loss" = renoir_15[9],
  "Not significant" = "grey50"
)

# Species directions
species_directions <- top15_species_labels |>
  mutate(direction = if_else(!is.na(p_fdr) & p_fdr < 0.05, estimate_direction, "Not significant"))

species_label_colors <- setNames(
  direction_colors[species_directions$direction],
  species_directions$species_label
)

# Metabolite directions
forest_metab <- read.csv("results/mlmodels_metab/perc_change_ffmi_v4_group/all/forest_results_top15.csv")

metabolite_directions <- forest_metab |>
  mutate(
    direction = case_when(
      is.na(p_fdr) | p_fdr >= 0.05 ~ "Not significant",
      estimate > 1 ~ "Higher in high loss",
      estimate < 1 ~ "Higher in modest/low loss",
      TRUE ~ "Not significant"
    ),
    label = metabolite_label(metabolite)
  )

metabolite_label_colors <- setNames(
  direction_colors[metabolite_directions$direction],
  metabolite_directions$label
)

# Metabolite super-pathway colours for strip
pathway_colors <- c(
  "Amino Acid"                        = "#3C9189",
  "Carbohydrate"                      = "#765135",
  "Cofactors and Vitamins"            = "#385C7E",
  "Energy"                            = "#3E6330",
  "Lipid"                             = "#CCA38D",
  "Nucleotide"                        = "#8EAFCF",
  "Peptide"                           = "#96B485",
  "Xenobiotics"                       = "#814956",
  "Partially Characterized Molecules" = "#6A4F7E"
)

tax <- as(tax_table(metab), "matrix")

# Matrices: metabolites in rows, species in columns
rho_mat <- cor_results |>
  select(species, metabolite, rho) |>
  pivot_wider(names_from = species, values_from = rho) |>
  column_to_rownames("metabolite") |>
  as.matrix()

fdr_mat <- cor_results |>
  select(species, metabolite, p_fdr) |>
  pivot_wider(names_from = species, values_from = p_fdr) |>
  column_to_rownames("metabolite") |>
  as.matrix()

rho_mat <- rho_mat[top15_metabolites, top15_species, drop = FALSE]
fdr_mat <- fdr_mat[top15_metabolites, top15_species, drop = FALSE]

# Super-pathway annotation
pathway_anno <- as.character(tax[top15_metabolites, "SUPER_PATHWAY"])
pathway_anno[is.na(pathway_anno)] <- "Unknown"
if ("Unknown" %in% pathway_anno) pathway_colors <- c(pathway_colors, "Unknown" = "#BBBBBB")

side_anno <- rowAnnotation(
  `Super pathway` = pathway_anno,
  col = list(`Super pathway` = pathway_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "bottom",
  annotation_name_rot = 45
)

# Readable labels
rownames(rho_mat) <- metabolite_label(rownames(rho_mat))
colnames(rho_mat) <- species_labels[colnames(rho_mat)]
rownames(fdr_mat) <- rownames(rho_mat)
colnames(fdr_mat) <- colnames(rho_mat)

# Colour both species and metabolite names by FFMI-loss association
species_gp <- gpar(fontsize = 9, fontface = "italic", col = species_label_colors[colnames(rho_mat)])
metabolite_gp <- gpar(fontsize = 9, col = metabolite_label_colors[rownames(rho_mat)])

stopifnot(
  !anyNA(species_label_colors[colnames(rho_mat)]),
  !anyNA(metabolite_label_colors[rownames(rho_mat)])
)

lgd_sig <- Legend(
  labels = c("q < 0.05", "q < 0.01", "q < 0.001"),
  title = "Significance",
  type = "points",
  pch = c("*", "**", "***"),
  legend_gp = gpar(col = "black")
)

lgd_direction <- Legend(
  labels = names(direction_colors),
  title = "Association with 1-year FFMI loss",
  legend_gp = gpar(fill = direction_colors)
)

ht <- Heatmap(
  rho_mat,
  name = "Spearman rho",
  col = col_fun,
  left_annotation = side_anno,
  rect_gp = gpar(col = "white", lwd = 0.7),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- fdr_mat[i, j]
    stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    if (stars != "") grid.text(stars, x, y, gp = gpar(fontsize = 10, col = "black"))
  },
  column_names_gp = species_gp,
  column_names_rot = 45,
  row_names_gp = metabolite_gp,
  row_names_side = "left",
  row_names_max_width = unit(14, "cm"),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_title = "ML-selected species x metabolites"
)

pdf(file.path(out_dir, "ml_species_metabolite_heatmap.pdf"), width = 14, height = 10)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lgd_sig, lgd_direction),
  padding = unit(c(20, 10, 4, 4), "mm")
)

dev.off()

#### Compact ComplexHeatmap: species and metabolites with >=1 significant correlation ####
sig_pairs <- cor_results |> filter(p_fdr < 0.05)

species_keep <- top15_species[top15_species %in% sig_pairs$species]
metabs_keep <- top15_metabolites[top15_metabolites %in% sig_pairs$metabolite]

if (length(species_keep) > 0 && length(metabs_keep) > 0) {
  rho_compact <- rho_mat[metabolite_label(metabs_keep), species_labels[species_keep], drop = FALSE]
  fdr_compact <- fdr_mat[metabolite_label(metabs_keep), species_labels[species_keep], drop = FALSE]

  pathway_compact <- as.character(tax[metabs_keep, "SUPER_PATHWAY"])
  pathway_compact[is.na(pathway_compact)] <- "Unknown"

  side_anno_compact <- rowAnnotation(
    `Super pathway` = pathway_compact,
    col = list(`Super pathway` = pathway_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8),
    annotation_name_side = "bottom",
    annotation_name_rot = 45
  )

  ht_compact <- Heatmap(
    rho_compact,
    name = "Spearman rho",
    col = col_fun,
    left_annotation = side_anno_compact,
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- fdr_compact[i, j]
      stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      if (stars != "") grid.text(stars, x, y, gp = gpar(fontsize = 10, col = "black"))
    },
    column_names_gp = gpar(fontsize = 9, fontface = "italic", col = species_label_colors[colnames(rho_compact)]),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 9, col = metabolite_label_colors[rownames(rho_compact)]),
    row_names_side = "left",
    row_names_max_width = unit(14, "cm"),
    cluster_rows = nrow(rho_compact) > 1,
    cluster_columns = ncol(rho_compact) > 1,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    column_title = "ML-selected species and metabolites with >=1 significant correlation"
  )

  pdf(file.path(out_dir, "ml_species_metabolite_heatmap_compact.pdf"), width = 14, height = 10)
  draw(
    ht_compact,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(lgd_sig, lgd_direction),
    padding = unit(c(20, 10, 4, 4), "mm")
  )
  dev.off()
}

cat("Done. Outputs in", out_dir, "\n")
