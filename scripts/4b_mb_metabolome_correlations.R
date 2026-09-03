# Correlations: top-15 species (deltaFFMI) and metabolites
# Barbara Verhaar

library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)

source("scripts/assets/functions.R")

out_dir <- "results/graphs/mb_metabolome_correlations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#### Microbes ####
fi <- get_feature_importance("results/mlmodels/delta_ffmi_v4/all", "delta_ffmi_v4_all", "reg")
feats <- top_features(fi, n = 15)$FeatName

# Baseline (v0) abundance for the top 15
mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb_mat <- t(as(mb@otu_table, "matrix"))
stopifnot(all(feats %in% colnames(mb_mat)))
species_wide <- as.data.frame(mb_mat[, feats, drop = FALSE]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v0$")) |>
  select(-sampleid)

# log10-transform (pseudocount)
bugs <- species_wide |> select(id, all_of(feats))
bugs[feats] <- map(bugs[feats], ~ {
  pseudo <- min(.x[.x > 0]) / 2
  log10(.x + pseudo)
})

#### Metabolomics ###
metab <- readRDS("data/processed_data/BARIA_metabolon_clean.RDS")
metab_mat <- as(otu_table(metab), "matrix")
metab_meta <- as(sample_data(metab), "data.frame")
rownames(metab_mat) <- metab_meta$id

#### Join on id, correlate every bug x metabolite pair ####
shared_ids <- intersect(bugs$id, rownames(metab_mat))
cat(length(shared_ids), "participants with both baseline microbiome and metabolome data\n")

bugs_mat <- bugs |>
  filter(id %in% shared_ids) |>
  column_to_rownames("id")
bugs_mat <- as.matrix(bugs_mat[shared_ids, , drop = FALSE])
metab_mat <- metab_mat[shared_ids, , drop = FALSE]
stopifnot(identical(rownames(bugs_mat), rownames(metab_mat)), !anyNA(bugs_mat), !anyNA(metab_mat))

cor_results <- expand_grid(bug = feats, metabolite = colnames(metab_mat)) |>
  mutate(
    test = map2(bug, metabolite, ~ cor.test(bugs_mat[, .x], metab_mat[, .y], method = "spearman", exact = FALSE)),
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

#### Scatter plots: every FDR-significant metabolite, one panel per bug ####
bug_labels <- species_label(feats)
names(bug_labels) <- feats

for (feat in unique(sig$bug)) {
  sig_metabs <- sig |> filter(bug == feat) |> arrange(p_fdr) |> pull(metabolite)
  if (length(sig_metabs) == 0) next

  # Layout scales with how many panels there actually are, so e.g. a single
  # significant metabolite doesn't get stretched across a fixed-width page
  panel_size <- 2.6
  ncol_used <- min(5, length(sig_metabs))
  nrow_used <- ceiling(length(sig_metabs) / ncol_used)

  df <- tibble(id = rownames(bugs_mat), abundance = bugs_mat[, feat]) |>
    left_join(
      metab_mat[, sig_metabs, drop = FALSE] |>
        as.data.frame() |>
        rownames_to_column("id") |>
        pivot_longer(-id, names_to = "metabolite", values_to = "level"),
      by = "id"
    ) |>
    left_join(sig |> filter(bug == feat) |> select(metabolite, rho, p_fdr), by = "metabolite") |>
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
      title = bug_labels[[feat]],
      x = "log10 relative abundance (baseline)", y = "Metabolite level (z-scored log10)"
    ) +
    theme_Publication() +
    theme(strip.text = element_text(size = 8), plot.title = element_text(size = 10))

  fname <- str_replace_all(bug_labels[[feat]], " ", "_")
  ggsave(
    file.path(out_dir, sprintf("correlation_scatter_%s.pdf", fname)),
    p, width = panel_size * ncol_used + 0.5, height = panel_size * nrow_used + 0.6, limitsize = FALSE
  )
}

#### ComplexHeatmap: bugs x metabolites ####
# Diverging olive<->mauve pair, gray neutral midpoint - matches the
# species x pathway heatmap colours in 4a_humann_pathways.R (MetBrewer
# Renoir[14]/[6]); olive = negative, mauve = positive correlation
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#939336", "#f0efec", "#AE7B9E"))
pathway_colors <- c(
  "Amino Acid"                         = "#2a78d6",
  "Carbohydrate"                       = "#eb6834",
  "Cofactors and Vitamins"             = "#1baf7a",
  "Energy"                             = "#eda100",
  "Lipid"                              = "#e87ba4",
  "Nucleotide"                         = "#008300",
  "Peptide"                            = "#4a3aa7",
  "Xenobiotics"                        = "#e34948",
  "Partially Characterized Molecules"  = "#898781"
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

# Build & save one bugs x `metabs` correlation heatmap
plot_corr_heatmap <- function(metabs, file, column_title, col_width, pdf_width_min,
                               column_fontsize = 8, cell_fontsize = 10,
                               transpose = FALSE, row_height = 0.35) {
  rho_mat <- cor_results |>
    filter(metabolite %in% metabs) |>
    select(bug, metabolite, rho) |>
    pivot_wider(names_from = metabolite, values_from = rho) |>
    column_to_rownames("bug") |>
    as.matrix()
  rho_mat <- rho_mat[feats, metabs, drop = FALSE]

  fdr_mat <- cor_results |>
    filter(metabolite %in% metabs) |>
    select(bug, metabolite, p_fdr) |>
    pivot_wider(names_from = metabolite, values_from = p_fdr) |>
    column_to_rownames("bug") |>
    as.matrix()
  fdr_mat <- fdr_mat[feats, metabs, drop = FALSE]

  bugs_active <- rownames(fdr_mat)[apply(fdr_mat < 0.05, 1, any)]
  cat(file, ":", length(bugs_active), "/", length(feats), "bugs have >=1 FDR-significant correlation\n")
  rho_mat <- rho_mat[bugs_active, , drop = FALSE]
  fdr_mat <- fdr_mat[bugs_active, , drop = FALSE]

  rownames(rho_mat) <- bug_labels[rownames(rho_mat)]
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

  pdf(file, width = width, height = height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       annotation_legend_list = list(lgd_sig),
       padding = unit(c(20, left_padding, 4, 4), "mm"))
  dev.off()
}

sig_metabs_all <- unique(sig$metabolite)
cat(length(sig_metabs_all), "distinct metabolites significant for >=1 bug\n")

if (length(sig_metabs_all) > 0) { # if anything is sig

  # Full heatmap: every metabolite significant for >=1 bug
  plot_corr_heatmap(
    sig_metabs_all, file.path(out_dir, "microbe_metabolite_heatmap_vertical.pdf"),
    column_title = "All significant metabolites",
    col_width = 0.45, pdf_width_min = 3, column_fontsize = 10,
    transpose = TRUE, row_height = 0.18
  )

  # Compact version: the 20 metabolites correlated with the most bugs
  metab_rank <- sig |>
    group_by(metabolite) |>
    summarise(n_sig = n(), max_abs_rho = max(abs(rho)), .groups = "drop") |>
    arrange(desc(n_sig), desc(max_abs_rho))
  sig_metabs_small <- metab_rank |> slice_head(n = 25) |> pull(metabolite)

  plot_corr_heatmap(
    sig_metabs_small, file.path(out_dir, "microbe_metabolite_heatmap_top20_vertical.pdf"),
    column_title = sprintf("Top %d metabolites", length(sig_metabs_small)),
    col_width = 0.45, pdf_width_min = 3, column_fontsize = 10,
    transpose = TRUE, row_height = 0.35
  )
}

cat("Done. Outputs in", out_dir, "\n")
