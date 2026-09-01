# Violin plot: distribution exp var across ML models
# Barbara Verhaar

library(tidyverse)
source("scripts/assets/functions.R")

out_dir <- "results/graphs/ml_explained_variance"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_defs_base <- tribble(
  ~outcome,                       ~label,                          ~group,        ~adj,
  "ffmi",                         "FFMI (v0)",                     "FFMI",        FALSE,
  "ffmi_adj_fmi",                 "FFMI (v0)\nadj. FMI",           "FFMI",        TRUE,
  "delta_ffmi_v4",                "Delta FFMI (1y)",               "Delta FFMI",  FALSE,
  "delta_ffmi_v4_adj_fmi",        "Delta FFMI (1y)\nadj. FMI",     "Delta FFMI",  TRUE,
  "perc_change_ffmi_v4",          "%change FFMI (1y)",             "%Delta FFMI", FALSE,
  "perc_change_ffmi_v4_adj_fmi",  "%change FFMI (1y)\nadj. FMI",   "%Delta FFMI", TRUE,
)

# Pair each outcome with its FMI-adjusted counterpart via a light/dark shade
# of the same hue, so the six violins read as three related pairs rather than
# six arbitrary colors.
fill_colors <- c(
  "FFMI.FALSE"        = "#8ec6f2",
  "FFMI.TRUE"          = "#2a78d6",
  "Delta FFMI.FALSE"   = "#a8dfc4",
  "Delta FFMI.TRUE"    = "#1baf7a",
  "%Delta FFMI.FALSE"  = "#f5c396",
  "%Delta FFMI.TRUE"   = "#eb6834"
)

plot_ev_violin <- function(plot_title, out_file) {
  df <- pmap_dfr(model_defs_base, function(outcome, label, group, adj) {
    it <- get_iterations_reg(file.path("results/mlmodels", outcome, "all"), paste0(outcome, "_all"))
    if (is.null(it)) {
      cat("  no XGBeast output found for", outcome, "(all), skipping\n")
      return(NULL)
    }
    tibble(label = label, group = group, adj = adj, explained_variance = it$Explained.Variance * 100)
  })
  if (nrow(df) == 0) {
    cat("  no models found, skipping plot\n")
    return(invisible(NULL))
  }

  df <- df |>
    mutate(
      # reversed so, after coord_flip(), FFMI (v0) reads at the top
      label = factor(label, levels = rev(model_defs_base$label)),
      fill_key = paste(group, adj, sep = ".")
    )

  medians <- df |>
    summarise(med = median(explained_variance), max_val = max(explained_variance), .by = label) |>
    mutate(label_text = sprintf("%.1f%%", med))

  p <- ggplot(df, aes(x = label, y = explained_variance, fill = fill_key)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_violin(trim = FALSE, alpha = 0.9, color = "grey30", linewidth = 0.3) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey20", fill = "white") +
    geom_text(
      data = medians, aes(x = label, y = max_val + 2, label = label_text),
      inherit.aes = FALSE, position = position_nudge(x = 0.32),
      hjust = 0, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    coord_flip() +
    labs(title = plot_title, x = "", y = "Explained variance (%)") +
    theme_Publication() +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(out_dir, out_file), p, width = 9, height = 6)
  cat("Saved", file.path(out_dir, out_file), "\n")
}

plot_ev_violin(
  "Explained variance across CV iterations - FFMI models, species (all subjects)",
  "ffmi_models_explained_variance_violin.pdf"
)

cat("Done. Output in", out_dir, "\n")
