# Process XGBeast models and make plots for each model
# Barbara Verhaar

# Packages
library(tidyverse)
library(phyloseq)
library(broom)
library(MetBrewer)

# Theme/colors
renoir_15 <- met.brewer("Renoir", n = 15)

source("scripts/assets/functions.R")

# Read aggregated classification metrics
get_metrics_class <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "class")
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "aggregated_metrics_classification.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Read per-iteration classification metrics
get_iterations_class <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "class")
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "model_results_per_iteration.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Read the input data for the ML models
read_input_data <- function(model_path) {
  d <- file.path(model_path, "input_data")
  X <- as.matrix(read.delim(file.path(d, "X_data.txt"), header = FALSE, sep = "\t"))
  feat_ids <- readLines(file.path(d, "feat_ids.txt"))
  subject_ids <- readLines(file.path(d, "subject_ids.txt"))
  colnames(X) <- feat_ids
  rownames(X) <- subject_ids
  y_file <- if (file.exists(file.path(d, "y_reg.txt"))) "y_reg.txt" else "y_binary.txt"
  y <- scan(file.path(d, y_file), quiet = TRUE)
  list(X = X, y = y, subject_ids = subject_ids)
}

# Logistic regression per feature (for binary ML models): odds ratio of the
# positive class per 1 log10 unit of species abundance
run_glm_forest <- function(data, feats, outcome, covariates = character(0)) {
  f <- reformulate(c("log10_abundance", covariates), response = outcome)
  purrr::map_dfr(feats, function(feat) {
    x <- data[[feat]]
    pseudo <- min(x[x > 0], na.rm = TRUE) / 2
    data$log10_abundance <- log10(x + pseudo)
    fit <- glm(f, data = data, family = binomial)
    broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) |>
      filter(term == "log10_abundance") |>
      mutate(species = feat)
  }) |>
    mutate(p_fdr = p.adjust(p.value, method = "BH"))
}

# Forest plot of a run_glm_forest() result: one row per feature, odds ratio
# + 95% CI (log scale), ordered by effect size.
plot_forest <- function(forest, title = "") {
  forest |>
    mutate(
      label = species_label(species),
      label = forcats::fct_reorder(label, estimate, .desc = TRUE)
    ) |>
    ggplot(aes(x = estimate, y = label, color = p_fdr < 0.05)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(size = 2) +
    scale_x_log10() +
    scale_color_manual(values = c(`TRUE` = renoir_15[6], `FALSE` = "grey40"), guide = "none") +
    labs(title = title, x = "Odds ratio per log10 abundance (95% CI)", y = NULL) +
    theme_Publication()
}

# Tidy top-n feature importance bar plot from top_features().
plot_feature_importance <- function(fi_top, title = "") {
  fi_top |>
    mutate(
      label = species_label(FeatName),
      label = forcats::fct_reorder(label, RelFeatImp)
    ) |>
    ggplot(aes(x = RelFeatImp, y = label)) +
    geom_col(fill = renoir_15[6], width = 0.7) +
    labs(title = title, x = "Relative importance", y = NULL) +
    theme_Publication()
}

# Abundance of the top features per outcome class, to check direction
plot_feature_by_group <- function(X, y, feats, class_labels, title = "") {
  labels <- species_label(feats)
  names(labels) <- feats

  df <- as.data.frame(X[, feats, drop = FALSE])
  colnames(df) <- feats
  df$.y <- factor(y, levels = c(0, 1), labels = class_labels)
  df <- df |>
    tidyr::pivot_longer(-.y, names_to = "species", values_to = "abundance") |>
    group_by(species) |>
    mutate(diff_median = median(abundance[.y == class_labels[2]]) - median(abundance[.y == class_labels[1]])) |>
    ungroup() |>
    mutate(
      facet_label = str_wrap(labels[species], 22),
      facet_label = forcats::fct_reorder(facet_label, -diff_median)
    )

  ggplot(df, aes(x = .y, y = abundance, fill = .y)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.6, color = "grey30") +
    scale_fill_manual(values = c(renoir_15[14], renoir_15[6]), guide = "none") +
    facet_wrap(~facet_label, scales = "free_y", ncol = 5) +
    labs(title = title, x = "", y = "Relative abundance") +
    theme_Publication() +
    theme(strip.text = element_text(size = 8), axis.text.x = element_text(size = 7))
}

groups <- c("all", "male", "female")
class_labels <- c("Modest/low %FFMI loss", "High %FFMI loss") # y = 0, y = 1

model_defs <- tribble(
  ~outcome,                             ~label,                                     ~is_adj,
  "perc_change_ffmi_v4_group",          "High vs. modest/low %FFMI loss (1y)",            FALSE,
  "perc_change_ffmi_v4_group_adj_fmi",  "High vs. modest/low %FFMI loss (1y, adj. FMI)",  TRUE,
)

### Data shared by the logistic models ###

# Covariates per subject; the outcome itself is taken from each model's
# y_binary.txt so the glm uses exactly the grouping XGBeast was trained on
meta_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS") |>
  select(id, sex, age_v0, fmi_v0, perc_change_weight_kg_v4) |>
  mutate(id = as.character(id))

# Baseline (v0) species abundance (same filtering as 3a_ml_prep.R), wide with `id`
mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb_mat <- t(as(mb@otu_table, "matrix"))
tk <- apply(mb_mat, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
species_wide <- as.data.frame(mb_mat[, tk]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v0$")) |>
  select(-sampleid)

#### Process each outcome ####
for (i in seq_len(nrow(model_defs))) {
  outcome <- model_defs$outcome[i]
  label   <- model_defs$label[i]
  is_adj  <- model_defs$is_adj[i]
  base <- file.path("results/mlmodels", outcome)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  cat("\n====", label, "====\n")

  ## Aggregate metrics + AUC summary across subgroups
  metrics <- map_dfr(groups, function(g) {
    m <- get_metrics_class(file.path(base, g), paste0(outcome, "_", g))
    if (is.null(m)) return(NULL)
    m |> mutate(subgroup = g)
  })
  if (nrow(metrics) > 0) {
    write.csv(metrics, file.path(base, paste0(outcome, "_metrics.csv")), row.names = FALSE)
    pl <- ggplot(metrics, aes(x = subgroup, y = Median.AUC)) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
      geom_point(size = 2) +
      geom_segment(aes(xend = subgroup, yend = 0.5)) +
      coord_flip() +
      labs(title = label, x = "", y = "Median ROC AUC") +
      theme_Publication()
    ggsave(file.path(base, paste0(outcome, "_auc.pdf")), pl, width = 4, height = 3)
    print(metrics)
  }

  ## Per-subgroup: feature importance, abundance per group, forest plot
  for (g in groups) {
    model_path <- file.path(base, g)
    fi <- get_feature_importance(model_path, paste0(outcome, "_", g), mode = "class")
    if (is.null(fi)) {
      cat(sprintf("  %-12s no XGBeast output found, skipping\n", g))
      next
    }
    fi_top <- top_features(fi, n = 15)
    feats <- fi_top$FeatName
    subtitle <- paste0(label, " (", g, ")")
    input <- read_input_data(model_path)

    p_imp <- plot_feature_importance(fi_top, title = subtitle)
    ggsave(file.path(model_path, "feature_importance_top15.pdf"), p_imp, width = 7, height = 5)

    p_box <- plot_feature_by_group(input$X, input$y, feats, class_labels, title = subtitle)
    ggsave(file.path(model_path, "feature_by_group_top15.pdf"), p_box, width = 12, height = 9)

    # Logistic regression on the exact subjects/outcome XGBeast trained on
    covariates <- c("age_v0", "perc_change_weight_kg_v4")
    if (g == "all") covariates <- c(covariates, "sex")
    if (is_adj) covariates <- c(covariates, "fmi_v0")
    model_data <- tibble(
      id = input$subject_ids |> str_remove("^BARIA_") |> str_remove("_v0$"),
      high_loss = input$y
    ) |>
      inner_join(species_wide, by = "id") |>
      inner_join(meta_wide, by = "id")
    forest <- run_glm_forest(model_data, feats, "high_loss", covariates)
    write.csv(forest, file.path(model_path, "forest_results_top15.csv"), row.names = FALSE)

    p_forest <- plot_forest(forest, title = subtitle)
    ggsave(file.path(model_path, "forest_plot_top15.pdf"), p_forest, width = 7, height = 6)

    cat(sprintf("  %-12s %d/%d features FDR < 0.05 (logistic)\n", g,
                sum(forest$p_fdr < 0.05), nrow(forest)))
  }
}

### AUC violin ###
violin_out_dir <- "results/graphs/ml_auc"
dir.create(violin_out_dir, recursive = TRUE, showWarnings = FALSE)

violin_defs <- tribble(
  ~ outcome,                            ~ violin_label,
  "perc_change_ffmi_v4_group",          "High %FFMI loss (1y)",
  "perc_change_ffmi_v4_group_adj_fmi",  "High %FFMI loss (1y)\nadj. FMI",
) |>
  left_join(model_defs |> select(outcome, is_adj), by = "outcome")

violin_fill_colors <- c(`FALSE` = renoir_15[8], `TRUE` = renoir_15[9])

auc_df <- pmap_dfr(violin_defs, function(outcome, violin_label, is_adj) {
  map_dfr(groups, function(g) {
    it <- get_iterations_class(file.path("results/mlmodels", outcome, g), paste0(outcome, "_", g))
    if (is.null(it)) {
      cat("  no XGBeast output found for", outcome, "(", g, "), skipping violin\n")
      return(NULL)
    }
    tibble(label = violin_label, subgroup = g, adj = is_adj, auc = it$ROC_AUC_scores)
  })
})

if (nrow(auc_df) == 0) {
  cat("  no models found, skipping AUC violin plot\n")
} else {
  auc_df <- auc_df |>
    mutate(
      # reversed so, after coord_flip(), the unadjusted model reads at the top
      label = factor(label, levels = rev(violin_defs$violin_label)),
      subgroup = factor(subgroup, levels = groups)
    )

  auc_medians <- auc_df |>
    summarise(med = median(auc), max_val = max(auc), .by = c(label, subgroup)) |>
    mutate(label_text = sprintf("%.2f", med))

  p_violin <- ggplot(auc_df, aes(x = label, y = auc, fill = adj)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_violin(trim = FALSE, alpha = 0.9, color = "grey30", linewidth = 0.3) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey20", fill = "white") +
    geom_text(
      data = auc_medians, aes(x = label, y = max_val + 0.03, label = label_text),
      inherit.aes = FALSE, position = position_nudge(x = 0.32),
      hjust = 0, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = violin_fill_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(~subgroup, ncol = 1) +
    coord_flip() +
    labs(
      title = "ROC AUC across CV iterations - %FFMI loss group, species",
      x = "", y = "ROC AUC"
    ) +
    theme_Publication() +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(violin_out_dir, "ffmi_group_models_auc_violin.pdf"), p_violin, width = 9, height = 8)
  cat("Saved", file.path(violin_out_dir, "ffmi_group_models_auc_violin.pdf"), "\n")
}
