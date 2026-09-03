# Correlations: top-15 species and macronutrients
# Barbara Verhaar

library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
source("scripts/assets/functions.R")

base_out_dir <- "results/graphs/mb_diet_correlations"

#### Microbes top 15 ####
fi <- get_feature_importance("results/mlmodels/perc_change_ffmi_v4/all", "perc_change_ffmi_v4_all", "reg")
feats <- top_features(fi, n = 15)$FeatName
mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb_mat <- t(as(mb@otu_table, "matrix"))
stopifnot(all(feats %in% colnames(mb_mat)))
species_wide <- as.data.frame(mb_mat[, feats, drop = FALSE]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v0$")) |>
  select(-sampleid)

# log10-transform (pseudocount = min(nonzero)/2 per species
bugs <- species_wide |> select(id, all_of(feats))
bugs[feats] <- map(bugs[feats], ~ {
  pseudo <- min(.x[.x > 0]) / 2
  log10(.x + pseudo)
})

#### Diet: baseline macronutrient intake energy-normalized ####
kcal_norm_vars <- c("Carbs_g", "Protein_g", "Fat_g", "SatFat_g", "Fibers_g", "Alcohol_g")
diet_vars <- c("TotalCal_kcal", str_c(kcal_norm_vars, "_per1000kcal"))
diet_all <- readRDS("data/raw_data/251002_BARIA_macronutrients.RDS") |>
  ungroup() |>
  filter(Visit == "V1") |>
  mutate(
    id = str_remove(ID, "^BARIA_"),
    across(all_of(kcal_norm_vars), ~ .x / TotalCal_kcal * 1000, .names = "{.col}_per1000kcal")
  ) |>
  select(id, diary_tool, all_of(diet_vars))

# Format p value because nothing left after FDR
format_pval <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))
bug_labels <- species_label(feats) # as a make unique
names(bug_labels) <- feats

# Color pal - matches the species x pathway heatmap in 4a_humann_pathways.R
# (MetBrewer Renoir[14]/[6]); olive = negative, mauve = positive correlation
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("#939336", "#f0efec", "#AE7B9E"))
lgd_sig <- Legend(
  labels = c("p < 0.05", "p < 0.01", "p < 0.001"),
  title = "Significance (nominal)",
  type = "points",
  pch = c("*", "**", "***"),
  legend_gp = gpar(col = "black"),
  background = "white"
)

# One pooled analysis across both diary tools (macronutrient variables are
# defined identically regardless of tool, and no participant used both), plus
# one analysis per tool
diet_groups <- c(levels(diet_all$diary_tool), "overall")

for (group in diet_groups) {
  label <- if (group == "overall") "all tools" else group
  cat("\n====", label, "====\n")
  out_dir <- file.path(base_out_dir, group)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  diet_mat <- diet_all |>
    {\(d) if (group == "overall") d else filter(d, diary_tool == group)}() |>
    select(-diary_tool) |>
    column_to_rownames("id") |>
    as.matrix()

  #### Join on id, correlate every bug x nutrient pair ####
  shared_ids <- intersect(bugs$id, rownames(diet_mat))
  cat(length(shared_ids), "participants with both baseline microbiome and", label, "dietary data\n")

  ## Filter because correlations are meaningless if not prevalent enough
  raw_bugs_tool <- species_wide |>
    filter(id %in% shared_ids) |>
    column_to_rownames("id")
  raw_bugs_tool <- as.matrix(raw_bugs_tool[shared_ids, , drop = FALSE])
  prevalent_feats <- colnames(raw_bugs_tool)[colMeans(raw_bugs_tool > 0.05) >= 0.3]
  cat(length(prevalent_feats), "/", length(feats), "bugs present in >=30% of", label, "participants\n")


  bugs_mat <- bugs |>
    filter(id %in% shared_ids) |>
    column_to_rownames("id")
  bugs_mat <- as.matrix(bugs_mat[shared_ids, prevalent_feats, drop = FALSE])
  diet_mat_tool <- diet_mat[shared_ids, , drop = FALSE]

  cor_results <- expand_grid(bug = prevalent_feats, nutrient = colnames(diet_mat_tool)) |>
    mutate(
      test = map2(bug, nutrient, ~ cor.test(bugs_mat[, .x], diet_mat_tool[, .y], method = "spearman", exact = FALSE)),
      rho = map_dbl(test, "estimate"),
      p.value = map_dbl(test, "p.value")
    ) |>
    select(-test) |>
    mutate(p_fdr = p.adjust(p.value, method = "BH"))

  write.csv(cor_results, file.path(out_dir, "microbe_diet_correlations.csv"), row.names = FALSE)

  sig <- cor_results |> filter(p.value < 0.05)
  cat(nrow(sig), "/", nrow(cor_results), "microbe-macronutrient pairs p < 0.05 (nominal)\n")

  #### Scatter plots: every nominally significant nutrient, one panel per bug ####
  for (feat in unique(sig$bug)) {
    sig_nutrients <- sig |> filter(bug == feat) |> arrange(p.value) |> pull(nutrient)
    if (length(sig_nutrients) == 0) next

    # Layout scales with how many panels there actually are, so e.g. a single
    # significant nutrient doesn't get stretched across a fixed-width page
    panel_size <- 2.6
    ncol_used <- min(5, length(sig_nutrients))
    nrow_used <- ceiling(length(sig_nutrients) / ncol_used)

    df <- tibble(id = rownames(bugs_mat), abundance = bugs_mat[, feat]) |>
      left_join(
        diet_mat_tool[, sig_nutrients, drop = FALSE] |>
          as.data.frame() |>
          rownames_to_column("id") |>
          pivot_longer(-id, names_to = "nutrient", values_to = "intake"),
        by = "id"
      ) |>
      left_join(sig |> filter(bug == feat) |> select(nutrient, rho, p.value), by = "nutrient") |>
      mutate(facet_label = fct_reorder(str_wrap(nutrient, 22), -rho))

    # rho/p text drawn inside each panel
    ann <- df |> distinct(facet_label, rho, p.value) |>
      mutate(label = sprintf("rho = %.2f\np %s", rho,
                              if_else(p.value < 0.001, format_pval(p.value), paste0("= ", format_pval(p.value)))))

    p <- ggplot(df, aes(x = abundance, y = intake)) +
      geom_point(alpha = 0.5, size = 1, color = "#2a78d6") +
      geom_smooth(method = "lm", se = FALSE, color = "#e34948", linewidth = 0.6) +
      geom_text(
        data = ann, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
        inherit.aes = FALSE, size = 2.5
      ) +
      facet_wrap(~facet_label, scales = "free", ncol = ncol_used) +
      labs(
        title = bug_labels[[feat]],
        x = "log10 relative abundance (baseline)", y = "Daily intake (baseline diary)"
      ) +
      theme_Publication() +
      theme(strip.text = element_text(size = 8), plot.title = element_text(size = 10))

    fname <- str_replace_all(bug_labels[[feat]], " ", "_")
    ggsave(
      file.path(out_dir, sprintf("correlation_scatter_%s.pdf", fname)),
      p, width = panel_size * ncol_used + 0.5, height = panel_size * nrow_used + 0.6, limitsize = FALSE
    )
  }

  #### ComplexHeatmap: nutrients (rows) x bugs (columns), vertical ####
  rho_mat <- cor_results |>
    select(bug, nutrient, rho) |>
    pivot_wider(names_from = nutrient, values_from = rho) |>
    column_to_rownames("bug") |>
    as.matrix()
  rho_mat <- rho_mat[prevalent_feats, diet_vars, drop = FALSE]

  p_mat <- cor_results |>
    select(bug, nutrient, p.value) |>
    pivot_wider(names_from = nutrient, values_from = p.value) |>
    column_to_rownames("bug") |>
    as.matrix()
  p_mat <- p_mat[prevalent_feats, diet_vars, drop = FALSE]

  # Bugs with no nominally significant nutrient correlation are dropped from
  # the row (pre-transpose)/column (post-transpose) axis
  bugs_active <- rownames(p_mat)[apply(p_mat < 0.05, 1, any)]
  cat(length(bugs_active), "/", length(prevalent_feats), "prevalent bugs have >=1 nominally significant nutrient correlation\n")
  rho_mat <- rho_mat[bugs_active, , drop = FALSE]
  p_mat <- p_mat[bugs_active, , drop = FALSE]

  if (length(bugs_active) > 0) {
    rownames(rho_mat) <- bug_labels[rownames(rho_mat)]
    rownames(p_mat) <- rownames(rho_mat)

    # Vertical: nutrients as rows, bugs as columns
    rho_mat <- t(rho_mat)
    p_mat <- t(p_mat)

    ht <- Heatmap(
      rho_mat,
      name = "Spearman rho",
      col = col_fun,
      rect_gp = gpar(col = "white", lwd = 0.7),
      cell_fun = function(j, i, x, y, width, height, fill) {
        p <- p_mat[i, j]
        stars <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
        if (stars != "") grid::grid.text(stars, x, y, gp = grid::gpar(fontsize = 10, col = "black"))
      },
      column_names_gp = grid::gpar(fontsize = 9),
      column_names_rot = 45,
      row_names_gp = grid::gpar(fontsize = 10),
      row_names_side = "left",
      row_dend_side = "right",
      column_title = sprintf("Baseline dietary macronutrients (%s)", label),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_dend = FALSE
    )

    cairo_pdf(file.path(out_dir, "microbe_diet_heatmap_vertical.pdf"),
        width = max(5, 3 + ncol(rho_mat) * 0.55), height = max(4, 3 + nrow(rho_mat) * 0.4))
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         annotation_legend_list = list(lgd_sig), merge_legends = TRUE,
         padding = unit(c(20, 8, 4, 4), "mm"))
    dev.off()
  } else {
    cat("No bug has a nominally significant nutrient correlation - skipping heatmap\n")
  }
}

cat("\nDone. Outputs in", base_out_dir, "\n")
