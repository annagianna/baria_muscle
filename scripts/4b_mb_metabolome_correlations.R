# Correlations: top-15 species (deltaFFMI) and metabolites
# Barbara Verhaar

# Packages
library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(MetBrewer)

source("scripts/assets/functions.R")
out_dir <- "results/graphs/mb_metabolome_correlations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Top-15 species by ML feature importance for deltaFFMI
fi <- get_feature_importance("results/mlmodels/perc_change_ffmi_v4/all", "perc_change_ffmi_v4_all", "reg")
top15_species <- top_features(fi, n = 15)$FeatName

# Baseline (v0) abundance for the top 15
baria_mb_baseline <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
species_mat <- t(as(baria_mb_baseline@otu_table, "matrix"))
stopifnot(all(top15_species %in% colnames(species_mat)))
species_wide <- as.data.frame(species_mat[, top15_species, drop = FALSE]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v0$")) |>
  select(-sampleid)

# log10-transform (pseudocount)
species_log10 <- species_wide |> select(id, all_of(top15_species))
species_log10[top15_species] <- map(species_log10[top15_species], ~ {
  pseudo <- min(.x[.x > 0]) / 2
  log10(.x + pseudo)
})

#### Metabolomics ####
metab <- readRDS("data/processed_data/BARIA_metabolon_clean.RDS")
metab_mat <- as(otu_table(metab), "matrix")
metab_meta <- as(sample_data(metab), "data.frame")
rownames(metab_mat) <- metab_meta$id

#### Join on id, correlate every species x metabolite pair ####
shared_ids <- intersect(species_log10$id, rownames(metab_mat))
cat(length(shared_ids), "participants with both baseline microbiome and metabolome data\n")

species_log10_mat <- species_log10 |>
  filter(id %in% shared_ids) |>
  column_to_rownames("id")
species_log10_mat <- as.matrix(species_log10_mat[shared_ids, , drop = FALSE])
metab_mat <- metab_mat[shared_ids, , drop = FALSE]
stopifnot(identical(rownames(species_log10_mat), rownames(metab_mat)), !anyNA(species_log10_mat), !anyNA(metab_mat))

cor_results <- expand_grid(species = top15_species, metabolite = colnames(metab_mat)) |>
  mutate(
    test = map2(species, metabolite, ~ cor.test(species_log10_mat[, .x], metab_mat[, .y], method = "spearman", exact = FALSE)),
    rho = map_dbl(test, "estimate"),
    p.value = map_dbl(test, "p.value")
  ) |>
  select(-test) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH"))

write.csv(cor_results, file.path(out_dir, "microbe_metabolite_correlations.csv"), row.names = FALSE)

sig <- cor_results |> filter(p_fdr < 0.05)
cat(nrow(sig), "/", nrow(cor_results), "microbe-metabolite pairs FDR < 0.05\n")

# FDR values below 0.001 round to "0.000" at 3 decimals - report a floor instead
format_fdr <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))

# Sphingolipid names shortening
metabolite_label <- function(x) {
  str_replace(
    x,
    "(\\(d\\d+:\\d+/\\d+:\\d+(?:\\([^()]*\\))?)(?:,\\s*d\\d+:\\d+/\\d+:\\d+(?:\\([^()]*\\))?)+(\\))",
    "\\1, ...\\2"
  )
}

#### Scatter plots: every FDR-significant metabolite, one panel per species ####
species_labels <- species_label(top15_species)
names(species_labels) <- top15_species

for (species_id in unique(sig$species)) {
  # Cap panels per page at the 30 strongest hits - a species with hundreds of
  # FDR-significant metabolites would otherwise blow the PDF past its 200in
  # page-size limit (nrow_used unbounded x panel_size) and fail at dev.off()
  sig_metabs <- sig |> filter(species == species_id) |> arrange(p_fdr) |> slice_head(n = 30) |> pull(metabolite)
  if (length(sig_metabs) == 0) next

  # Layout scales with how many panels there actually are, so e.g. a single
  # significant metabolite doesn't get stretched across a fixed-width page
  panel_size <- 2.6
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
    mutate(facet_label = fct_reorder(str_wrap(metabolite_label(metabolite), 22), -rho))

  # rho/FDR text drawn inside each panel (not the facet strip) - one label
  # per facet, positioned in the top-left corner of its own (free) scale
  ann <- df |> distinct(facet_label, rho, p_fdr) |>
    mutate(label = sprintf("rho = %.2f\nFDR %s", rho,
                            if_else(p_fdr < 0.001, format_fdr(p_fdr), paste0("= ", format_fdr(p_fdr)))))

  p <- ggplot(df, aes(x = abundance, y = level)) +
    geom_point(alpha = 0.5, size = 1, color = "#2a78d6") +
    geom_smooth(method = "lm", se = FALSE, color = "#e34948", linewidth = 0.6) +
    geom_text(
      data = ann, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
      inherit.aes = FALSE, size = 2.5
    ) +
    facet_wrap(~facet_label, scales = "free", ncol = ncol_used) +
    labs(
      title = species_labels[[species_id]],
      x = "log10 relative abundance (baseline)", y = "Metabolite level (z-scored log10)"
    ) +
    theme_Publication() +
    theme(strip.text = element_text(size = 8), plot.title = element_text(size = 10))

  fname <- str_replace_all(species_labels[[species_id]], " ", "_")
  ggsave(
    file.path(out_dir, sprintf("correlation_scatter_%s.pdf", fname)),
    p, width = panel_size * ncol_used + 0.5, height = panel_size * nrow_used + 0.6, limitsize = FALSE
  )
}

#### ComplexHeatmap: species x metabolites ####
# Diverging colour scale matching the species x pathway heatmap in
# 4a_humann_pathways.R (MetBrewer Renoir[14]/[6]); negative = Renoir[14],
# positive = Renoir[6], scale bounded by the actual max |rho|
renoir_15 <- met.brewer("Renoir", n = 15)
rho_max <- max(abs(cor_results$rho), na.rm = TRUE)
col_fun <- colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[14], "white", renoir_15[6]))
# Super pathway strip: a muted, earthy 9-colour set (moderate chroma ~35,
# lightness alternating so adjacent picks aren't just hue-apart) chosen to
# stay clear of both diverging-scale endpoints (Lab distance >=23 from each)
# while maximizing separation between categories themselves (min pairwise Lab
# distance ~23) - straight Renoir subsampling put 6 of 9 categories within a
# 60-degree hue wedge (pink/red/orange/gold/mustard) and was hard to tell
# apart; this keeps the same muted/painterly character without that crowding
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
tax <- as(tax_table(metab), "matrix")

# Significance legend (drawn once, reused by every heatmap) - stars belong in
# a proper legend, not spelled out in the title.
lgd_sig <- Legend(
  labels = c("FDR < 0.05", "FDR < 0.01", "FDR < 0.001"),
  title = "Significance",
  type = "points",
  pch = c("*", "**", "***"),
  legend_gp = gpar(col = "black"),
  background = "white"
)

# Build & save one species x `metabs` correlation heatmap
plot_corr_heatmap <- function(metabs, file, column_title, col_width, pdf_width_min,
                               column_fontsize = 8, cell_fontsize = 10,
                               transpose = FALSE, row_height = 0.35) {
  rho_mat <- cor_results |>
    filter(metabolite %in% metabs) |>
    select(species, metabolite, rho) |>
    pivot_wider(names_from = metabolite, values_from = rho) |>
    column_to_rownames("species") |>
    as.matrix()
  rho_mat <- rho_mat[top15_species, metabs, drop = FALSE]

  fdr_mat <- cor_results |>
    filter(metabolite %in% metabs) |>
    select(species, metabolite, p_fdr) |>
    pivot_wider(names_from = metabolite, values_from = p_fdr) |>
    column_to_rownames("species") |>
    as.matrix()
  fdr_mat <- fdr_mat[top15_species, metabs, drop = FALSE]

  species_active <- rownames(fdr_mat)[apply(fdr_mat < 0.05, 1, any)]
  cat(file, ":", length(species_active), "/", length(top15_species), "species have >=1 FDR-significant correlation\n")
  rho_mat <- rho_mat[species_active, , drop = FALSE]
  fdr_mat <- fdr_mat[species_active, , drop = FALSE]

  rownames(rho_mat) <- species_labels[rownames(rho_mat)]
  colnames(rho_mat) <- metabolite_label(colnames(rho_mat))
  rownames(fdr_mat) <- rownames(rho_mat)
  colnames(fdr_mat) <- colnames(rho_mat)

  # Super pathway strip
  pathway_anno <- tax[metabs, "SUPER_PATHWAY"]

  if (transpose) {
    rho_mat <- t(rho_mat)
    fdr_mat <- t(fdr_mat)
    side_anno <- rowAnnotation(
      `Super pathway` = pathway_anno,
      col = list(`Super pathway` = pathway_colors),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 8),
      annotation_name_side = "bottom",
      annotation_name_rot = 45
    )
  } else {
    side_anno <- HeatmapAnnotation(
      `Super pathway` = pathway_anno,
      col = list(`Super pathway` = pathway_colors),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 8)
    )
  }

  ht <- Heatmap(
    rho_mat,
    name = "Spearman rho",
    col = col_fun,
    top_annotation = if (!transpose) side_anno else NULL,
    left_annotation = if (transpose) side_anno else NULL,
    rect_gp = gpar(col = "white", lwd = 0.7),
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- fdr_mat[i, j]
      stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      if (stars != "") grid::grid.text(stars, x, y, gp = grid::gpar(fontsize = cell_fontsize, col = "black"))
    },
    column_names_gp = grid::gpar(fontsize = column_fontsize),
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 10),
    row_names_side = "left",
    row_names_max_width = unit(14, "cm"),
    row_dend_side = "right",
    column_title = column_title,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_dend = FALSE
  )

  # Transposed heatmaps need extra width for the (now much longer) metabolite
  # row names plus the super-pathway strip alongside them.
  width <- max(pdf_width_min, 3 + ncol(rho_mat) * col_width) + if (transpose) 4.5 else 0
  height <- if (transpose) max(6, 3 + nrow(rho_mat) * row_height) else 9
  left_padding <- if (transpose) 10 else 8

  # PDF page dimensions are capped at 200in (14400pt) by the format itself -
  # a heatmap with many rows/columns can otherwise silently blow past that
  # and fail at dev.off() with a generic "write failed"
  cat(file, ": requested pdf size", round(width, 1), "x", round(height, 1), "in (", nrow(rho_mat), "rows x", ncol(rho_mat), "cols )\n")
  if (width > 200 || height > 200) {
    warning(file, ": requested size exceeds the 200in PDF page limit (", round(width, 1), "x", round(height, 1), "in) - clamping")
  }
  width <- min(width, 200)
  height <- min(height, 200)

  quartz(type = "pdf", file = file, width = width, height = height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       annotation_legend_list = list(lgd_sig),
       padding = unit(c(20, left_padding, 4, 4), "mm"))
  dev.off()
}

sig_metabs_all <- unique(sig$metabolite)
cat(length(sig_metabs_all), "distinct metabolites significant for >=1 species\n")

if (length(sig_metabs_all) > 0) { # if anything is sig

  # Full heatmap: every metabolite significant for >=1 species
  plot_corr_heatmap(
    sig_metabs_all, file.path(out_dir, "microbe_metabolite_heatmap_vertical.pdf"),
    column_title = "All significant metabolites",
    col_width = 0.45, pdf_width_min = 3, column_fontsize = 10,
    transpose = TRUE, row_height = 0.18
  )

  # Compact version: the 20 metabolites correlated with the most species
  metab_rank <- sig |>
    group_by(metabolite) |>
    summarise(n_sig = n(), max_abs_rho = max(abs(rho)), .groups = "drop") |>
    arrange(desc(n_sig), desc(max_abs_rho))
  sig_metabs_small <- metab_rank |> slice_head(n = 20) |> pull(metabolite)

  plot_corr_heatmap(
    sig_metabs_small, file.path(out_dir, "microbe_metabolite_heatmap_top20_vertical.pdf"),
    column_title = sprintf("Top %d metabolites", length(sig_metabs_small)),
    col_width = 0.45, pdf_width_min = 3, column_fontsize = 10,
    transpose = TRUE, row_height = 0.35
  )
}

cat("Done. Outputs in", out_dir, "\n")
