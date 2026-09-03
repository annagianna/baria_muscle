# Functional pathways (HUMAnN) and top 15 species
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(tidyverse)
library(grid)
library(phyloseq)
library(MetBrewer)
library(ComplexHeatmap)
library(circlize)

dir.create("results/graphs/HUMAnN", recursive = TRUE, showWarnings = FALSE)

renoir_15 <- met.brewer("Renoir", n = 15)

source("scripts/assets/functions.R")

# Data
humann_long <- readRDS("data/processed_data/BARIA_humann_pathways_long.RDS")
baria_mb <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
baria_muscle_long <- readRDS("data/processed_data/BARIA_muscle_long.RDS")

# MetaCyc/HUMAnN pathway names embed literal HTML entities (e.g. "&beta;")
# instead of the Greek letters themselves - decode once so every downstream
# plot gets clean text
html_entities <- c("&alpha;" = "α", "&beta;" = "β", "&gamma;" = "γ", "&delta;" = "δ")
humann_long <- humann_long |>
  mutate(pathway_name = str_replace_all(pathway_name, html_entities))

# Top-15 species by ML feature importance for deltaFFMI
fi <- get_feature_importance("results/mlmodels/perc_change_ffmi_v4/all", "perc_change_ffmi_v4_all", "reg")
feats <- top_features(fi, n = 15)$FeatName

species_lookup <- tibble(
  species = feats,
  species_label = species_label(feats)
)

# Multi-visit abundance for the top 15 species (baria_mb spans v0/v4/v5,
# unlike the baseline-only object used in 4b/4c)
top_species_long <- as(otu_table(baria_mb), "matrix") |>
  as.data.frame() |>
  rownames_to_column(var = "species") |>
  inner_join(species_lookup, by = "species") |>
  pivot_longer(
    cols = -c(species, species_label),
    names_to = "Sample",
    values_to = "species_abundance"
  ) |>
  left_join(
    as(sample_data(baria_mb), "data.frame") |>
      rownames_to_column("Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  )

# Filter HUMAnN pathways based on baseline prevalence (>5 CPM in at least 30% of samples)
humann_keep_v0 <- humann_long |> 
  filter(visit == "v0") |> 
  group_by(pathway_id, pathway_name) |> 
  summarize(prevalence = mean(pathway_abundance > 5), .groups = "drop") |> 
  filter(prevalence >= 0.30)

# Keep filtered baseline pathways across all visits (v0, v4, v5) & join top-15 ML species
humann_top_species <- humann_long |>
  filter(pathway_id %in% humann_keep_v0$pathway_id) |>
  inner_join(top_species_long, by = c("Sample", "id", "visit"), relationship = "many-to-many")

#### Correlations ####
#### Species-pathway Spearman correlations per visit ####
humann_species_cor <- humann_top_species |>
  group_by(visit, species, species_label, pathway_id, pathway_name) |>
  summarise(
    n = sum(complete.cases(species_abundance, pathway_abundance)),
    test = list(cor.test(species_abundance, pathway_abundance, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) |>
  mutate(
    rho = map_dbl(test, ~ unname(.x$estimate)),
    p = map_dbl(test, ~ .x$p.value)
  ) |>
  select(-test) |>
  group_by(visit, species) |>
  mutate(padj = p.adjust(p, method = "BH")) |>
  ungroup()

# Strongest significant baseline pathway associations per species
humann_species_cor_top10 <- humann_species_cor |>
  filter(visit == "v0",padj < 0.05) |>
  mutate(direction = if_else(rho > 0, "positive", "negative")) |>
  group_by(species_label, direction) |> 
  arrange(desc(abs(rho)), .by_group = TRUE) |> 
  slice_head(n = 10) |> 
  mutate(rank_direction = row_number()) |> 
  ungroup() |> 
  group_by(species_label) |> 
  mutate(# Select up to 10 strongest signif baseline pathway associations per species (up to 5 from each direction where available)
    priority = if_else(rank_direction <= 5, 1L, 2L)
  ) |>
  arrange(priority, desc(abs(rho)), .by_group = TRUE) |> 
  slice_head(n = 10) |> 
  ungroup() |> 
  select(species_label, pathway_id, pathway_name, rho, padj) |> 
  print(n = 40)

# Save selected pathways for downstream analyses
humann_pathways_top <- humann_species_cor_top10 |>
  distinct(pathway_id, pathway_name)
saveRDS(humann_pathways_top, "data/processed_data/HUMAnN_selected_pathways.RDS")

# Pull pathway IDs
pathways_top <- humann_pathways_top |>
  pull(pathway_id)

### Plot ###
# Plot data
humann_species_plot_data <- humann_species_cor |>
  filter(visit == "v0", pathway_id %in% pathways_top) |>
  mutate(
    signif = padj < 0.05,
    species_label = factor(species_label, levels = species_lookup$species_label)
  )

rho_max <- max(abs(humann_species_plot_data$rho), na.rm = TRUE)

## ComplexHeatmap: species x pathway
col_fun_species <- colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[14], "white", renoir_15[6]))

# Tiered significance legend (drawn once, reused by both orientations and by
# the pathway x clinical heatmap below) - stars belong in a proper legend,
# not spelled out in the title.
lgd_sig_pathway <- Legend(
  labels = c("FDR < 0.05", "FDR < 0.01", "FDR < 0.001"),
  title = "Significance",
  type = "points",
  pch = c("*", "**", "***"),
  legend_gp = gpar(col = "black"),
  background = "white"
)

rho_mat <- humann_species_plot_data |>
  select(species_label, pathway_name, rho) |>
  pivot_wider(names_from = pathway_name, values_from = rho) |>
  column_to_rownames("species_label") |>
  as.matrix()

padj_mat <- humann_species_plot_data |>
  select(species_label, pathway_name, padj) |>
  pivot_wider(names_from = pathway_name, values_from = padj) |>
  column_to_rownames("species_label") |>
  as.matrix()
padj_mat <- padj_mat[rownames(rho_mat), colnames(rho_mat)]

# Species labels are always italic, pathway labels are not - built as rows =
# species, columns = pathway, then optionally transposed for the wide/pptx
# orientation, swapping which axis carries the italic font.
plot_species_pathway_heatmap <- function(rho_mat, padj_mat, file, transpose, width, height) {
  if (transpose) {
    rho_mat <- t(rho_mat)
    padj_mat <- t(padj_mat)
  }
  species_gp <- gpar(fontsize = 8, fontface = "italic")
  pathway_gp <- gpar(fontsize = 7)

  ht <- Heatmap(
    rho_mat,
    name = "Spearman\nrho",
    col = col_fun_species,
    rect_gp = gpar(col = "white", lwd = 0.3),
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- padj_mat[i, j]
      stars <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      if (stars != "") grid::grid.text(stars, x, y, gp = grid::gpar(fontsize = 9, col = "black"))
    },
    row_names_gp = if (transpose) pathway_gp else species_gp,
    column_names_gp = if (transpose) species_gp else pathway_gp,
    column_names_rot = 45,
    row_names_side = "left",
    row_names_max_width = unit(10, "cm"),
    row_dend_side = "right",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_dend = FALSE
  )

  # Bottom padding clears the 45-degree column labels' left-hanging overhang;
  # left padding clears the row names, now on the left (wider when those rows
  # are the long HUMAnN pathway names, i.e. the transposed orientation).
  left_padding <- if (transpose) 14 else 6

  cairo_pdf(file, width = width, height = height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       annotation_legend_list = list(lgd_sig_pathway),
       padding = unit(c(20, left_padding, 4, 4), "mm"))
  dev.off()
}

# Species on x axis (pathway rows) - wide, for the long pathway names now on the left
plot_species_pathway_heatmap(
  rho_mat, padj_mat, "results/graphs/HUMAnN/HUMAnN_species_pathway_heatmap_x.pdf",
  transpose = TRUE, width = 14, height = 17
)

# Species on y axis (for pptx)
plot_species_pathway_heatmap(
  rho_mat, padj_mat, "results/graphs/HUMAnN/HUMAnN_species_pathway_heatmap_y.pdf",
  transpose = FALSE, width = 20, height = 8
)

# Rank a pathway by how many sig corr
rank_pathways <- function(df) {
  df |>
    group_by(pathway_id, pathway_name) |>
    summarise(n_sig = sum(signif), max_abs_rho = max(abs(rho)), .groups = "drop") |>
    arrange(desc(n_sig), desc(max_abs_rho))
}

# Compact heatmap: the 15 pathways implicated for the most species
pathway_rank <- rank_pathways(humann_species_plot_data)
pathways_compact <- pathway_rank |> slice_head(n = 15) |> pull(pathway_name)

plot_species_pathway_heatmap(
  rho_mat[, pathways_compact, drop = FALSE], padj_mat[, pathways_compact, drop = FALSE],
  "results/graphs/HUMAnN/HUMAnN_species_pathway_heatmap_top15.pdf",
  transpose = TRUE, width = 12, height = 6
)

#### Correlations pathways-clinical features ####
# Clinical vars of interest
clinical_vars <- c(
  "ffmi", "fmi", "bmi", "wc_cm", "hba1c_mmolmol", "glucose_mmoll_mmt_0", "homa_ir", "homa_b",
  "systolic_bp_mmhg", "diastolic_bp_mmhg",
  "total_cholesterol_mmoll", "ldl_cholesterol_mmoll", "hdl_cholesterol_mmoll", "triglycerides_mmoll",
  "crp_mgl", "creatinine_umoll", "egfr_mlmin", "gammagt_ul", "asat_ul", "alat_ul"
)

# Join HUMAnN data with clinical vars of interest
humann_clinical <- humann_long |>
  filter(
    visit == "v0",
    pathway_id %in% pathways_top
  ) |>
  inner_join(
    baria_muscle_long |>
      filter(visit == "v0") |>
      select(id, visit, all_of(clinical_vars)),
    by = c("id", "visit")
  ) |>
  pivot_longer(
    cols = all_of(clinical_vars),
    names_to = "clinical_feature",
    values_to = "clinical_value"
  )

# Spearman correlations between pathways and clinical features
humann_clinical_cor <- humann_clinical |>
  group_by(clinical_feature, pathway_id, pathway_name) |>
  summarise(
    n = sum(complete.cases(pathway_abundance, clinical_value)),
    test = list(
      cor.test(pathway_abundance, clinical_value, method = "spearman", exact = FALSE)),
    .groups = "drop"
  ) |>
  mutate(
    rho = map_dbl(test, ~ unname(.x$estimate)),
    p = map_dbl(test, ~ .x$p.value)
  ) |>
  select(-test) |>
  ungroup() |> 
  mutate(padj = p.adjust(p, method = "BH"))


### Plot ###
# Plot data
humann_clinical_plot_data <- humann_clinical_cor |>
  mutate(
    signif = padj < 0.05,
    clinical_feature = dplyr::recode(
      clinical_feature,
      ffmi = "FFMI",
      fmi = "FMI",
      bmi = "BMI",
      wc_cm = "Waist circumference",
      hba1c_mmolmol = "HbA1c",
      glucose_mmoll_mmt_0 = "Fasting glucose",
      homa_ir = "HOMA-IR",
      homa_b = "HOMA-B",
      systolic_bp_mmhg = "Systolic BP",
      diastolic_bp_mmhg = "Diastolic BP",
      total_cholesterol_mmoll = "Total cholesterol",
      ldl_cholesterol_mmoll = "LDL cholesterol",
      hdl_cholesterol_mmoll = "HDL cholesterol",
      triglycerides_mmoll = "Triglycerides",
      crp_mgl = "CRP",
      creatinine_umoll = "Creatinine",
      egfr_mlmin = "eGFR",
      gammagt_ul = "Gamma-GT",
      asat_ul = "ASAT",
      alat_ul = "ALAT"
    )
  )

# Symmetric colour scale
rho_max_clinical <- max(abs(humann_clinical_plot_data$rho), na.rm = TRUE)

col_fun_clinical <- colorRamp2(c(-rho_max_clinical, 0, rho_max_clinical), c(renoir_15[2], "white", renoir_15[6]))

## ComplexHeatmap: pathway x clinical feature
rho_mat_clinical <- humann_clinical_plot_data |>
  select(clinical_feature, pathway_name, rho) |>
  pivot_wider(names_from = clinical_feature, values_from = rho) |>
  column_to_rownames("pathway_name") |>
  as.matrix()

padj_mat_clinical <- humann_clinical_plot_data |>
  select(clinical_feature, pathway_name, padj) |>
  pivot_wider(names_from = clinical_feature, values_from = padj) |>
  column_to_rownames("pathway_name") |>
  as.matrix()
padj_mat_clinical <- padj_mat_clinical[rownames(rho_mat_clinical), colnames(rho_mat_clinical)]

plot_clinical_pathway_heatmap <- function(rho_mat, padj_mat, file, width, height) {
  ht <- Heatmap(
    rho_mat,
    name = "Spearman\nrho",
    col = col_fun_clinical,
    rect_gp = gpar(col = "white", lwd = 0.3),
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- padj_mat[i, j]
      stars <- if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      if (stars != "") grid::grid.text(stars, x, y, gp = grid::gpar(fontsize = 9, col = "black"))
    },
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    row_names_side = "left",
    row_names_max_width = unit(10, "cm"),
    row_dend_side = "right",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_dend = FALSE
  )

  cairo_pdf(file, width = width, height = height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       annotation_legend_list = list(lgd_sig_pathway),
       padding = unit(c(20, 10, 4, 4), "mm"))
  dev.off()
}

plot_clinical_pathway_heatmap(
  rho_mat_clinical, padj_mat_clinical,
  "results/graphs/HUMAnN/HUMAnN_clinical_pathway_heatmap.pdf",
  width = 15, height = 10
)

# Compact heatmap: the 15 pathways implicated for the most clinical features
pathway_rank_clinical <- rank_pathways(humann_clinical_plot_data)
pathways_compact_clinical <- pathway_rank_clinical |> slice_head(n = 15) |> pull(pathway_name)

plot_clinical_pathway_heatmap(
  rho_mat_clinical[pathways_compact_clinical, , drop = FALSE],
  padj_mat_clinical[pathways_compact_clinical, , drop = FALSE],
  "results/graphs/HUMAnN/HUMAnN_clinical_pathway_heatmap_top15.pdf",
  width = 10, height = 6
)
