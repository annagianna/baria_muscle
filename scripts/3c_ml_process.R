# Process XGBeast models and make plots for each model
# Barbara Verhaar

library(tidyverse)
library(phyloseq)
library(lmerTest)
library(broom.mixed)
library(broom)
source("scripts/assets/functions.R")

# Read aggregated regression metrics
get_metrics_reg <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "reg")
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "aggregated_metrics_regression.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Read per-iteration regression metrics
get_iterations_reg <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "reg")
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

# Cross-sectional linear model (for cross sectional ML models)
run_lm_forest <- function(data, feats, outcome, covariates = character(0)) {
  f <- reformulate(c("log10_abundance", covariates), response = outcome)
  purrr::map_dfr(feats, function(feat) {
    x <- data[[feat]]
    pseudo <- min(x[x > 0], na.rm = TRUE) / 2
    data$log10_abundance <- log10(x + pseudo)
    fit <- lm(f, data = data)
    broom::tidy(fit, conf.int = TRUE) |>
      filter(term == "log10_abundance") |>
      mutate(species = feat)
  }) |>
    mutate(p_fdr = p.adjust(p.value, method = "BH"))
}

# LMM model (for delta ML models)
run_lmm_forest <- function(species_wide, long_data, feats, outcome, follow_up, covariates = character(0)) {
  f <- reformulate(c("log10_abundance * visit", covariates, "(1 | id)"), response = outcome)
  model_data <- species_wide |>
    select(id, all_of(feats)) |>
    tidyr::pivot_longer(all_of(feats), names_to = "species", values_to = "baseline_abundance") |>
    inner_join(long_data, by = "id", relationship = "many-to-many") |>
    filter(visit %in% c("v0", follow_up)) |>
    droplevels() |>
    group_by(species) |>
    mutate(
      pseudo = min(baseline_abundance[baseline_abundance > 0], na.rm = TRUE) / 2,
      log10_abundance = log10(baseline_abundance + pseudo)
    ) |>
    ungroup()

  model_data |>
    group_by(species) |>
    tidyr::nest() |>
    mutate(
      model = purrr::map(data, ~ lmerTest::lmer(f, data = .x, REML = FALSE)),
      results = purrr::map(model, ~ broom.mixed::tidy(.x, effects = "fixed", conf.int = TRUE))
    ) |>
    select(species, results) |>
    tidyr::unnest(results) |>
    filter(str_detect(term, "log10_abundance:visit")) |>
    ungroup() |>
    mutate(p_fdr = p.adjust(p.value, method = "BH"))
}

# Forest plot of a run_lm_forest()/run_lmm_forest() result: one row per
# feature, point estimate + 95% CI, ordered by effect size.
plot_forest <- function(forest, title = "") {
  forest |>
    mutate(
      label = species_label(species),
      label = forcats::fct_reorder(label, estimate, .desc = TRUE)
    ) |>
    ggplot(aes(x = estimate, y = label, color = p_fdr < 0.05)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(size = 2) +
    scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), guide = "none") +
    labs(title = title, x = "Effect estimate (95% CI)", y = NULL) +
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
    geom_col(fill = "#2c7fb8", width = 0.7) +
    labs(title = title, x = "Relative importance", y = NULL) +
    theme_Publication()
}

# Correlation plots to check dis
plot_feature_correlations <- function(X, y, feats, title = "") {
  labels <- species_label(feats)
  names(labels) <- feats

  df <- as.data.frame(X[, feats, drop = FALSE])
  colnames(df) <- feats
  df$.y <- y
  df <- df |>
    tidyr::pivot_longer(-.y, names_to = "species", values_to = "abundance") |>
    group_by(species) |>
    mutate(rho = cor(abundance, .y, method = "spearman", use = "complete.obs")) |>
    ungroup() |>
    mutate(
      facet_label = str_wrap(labels[species], 22),
      facet_label = forcats::fct_reorder(facet_label, -rho)
    )

  # rho text drawn inside each panel (not the facet strip) - one label per facet
  ann <- df |> distinct(facet_label, rho) |>
    mutate(label = sprintf("rho = %.2f", rho))

  ggplot(df, aes(x = abundance, y = .y)) +
    geom_point(alpha = 0.5, size = 1, color = "#2c7fb8") +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.6) +
    geom_text(
      data = ann, aes(label = label), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
      inherit.aes = FALSE, size = 2.5
    ) +
    facet_wrap(~facet_label, scales = "free_x", ncol = 5) +
    labs(title = title, x = "Relative abundance", y = "Outcome") +
    theme_Publication() +
    theme(strip.text = element_text(size = 8))
}

groups <- c("all", "male", "female")

model_defs <- tribble(
  ~outcome,                       ~label,                             ~is_change, ~response_var, ~timepoint, ~is_adj,
  "ffmi",                         "FFMI (v0)",                        FALSE,      "ffmi",         "v0",       FALSE,
  "ffmi_adj_fmi",                 "FFMI (v0, adj. FMI)",              FALSE,      "ffmi",         "v0",       TRUE,
  "fmi_v0",                       "FMI (v0)",                         FALSE,      "fmi",          "v0",       FALSE,
  "ffmi_v0_matched",              "FFMI (v0, matched to v4 mb subset)", FALSE,    "ffmi",         "v0",       FALSE,
  "fmi_v0_matched",               "FMI (v0, matched to v4 mb subset)",  FALSE,    "fmi",          "v0",       FALSE,
  "ffmi_v4",                      "FFMI (v4)",                        FALSE,      "ffmi",         "v4",       FALSE,
  "fmi_v4",                       "FMI (v4)",                         FALSE,      "fmi",          "v4",       FALSE,
  "delta_ffmi_v4",                "Delta FFMI (1 year change)",               TRUE,       "ffmi",         NA,         FALSE,
  "delta_ffmi_v4_adj_fmi",        "Delta FFMI (1 year change, adj. FMI)",     TRUE,       "ffmi",         NA,         TRUE,
  "perc_change_ffmi_v4",          "%change FFMI (1 year change)",             TRUE,       "ffmi",         NA,         FALSE,
  "perc_change_ffmi_v4_adj_fmi",  "%change FFMI (1 year change, adj. FMI)",   TRUE,       "ffmi",         NA,         TRUE,
)

### Data shared by theLM/LMM models ###

meta_full <- readRDS("data/processed_data/BARIA_muscle_long.RDS")
baseline_fmi <- meta_full |> filter(visit == "v0") |> select(id, fmi_v0 = fmi)

# Cross-sectional outcomes at v0 (ffmi, ffmi_adj_fmi, fmi_v0) or v4 (ffmi_v4,
# fmi_v4); lm() drops NA rows per-formula, so no need to pre-filter on
# whichever of ffmi/fmi happens to be this outcome's response_var.
meta_cs_by_timepoint_group <- list(
  v0 = list(
    all    = meta_full |> filter(visit == "v0"),
    male   = meta_full |> filter(visit == "v0", sex == "male"),
    female = meta_full |> filter(visit == "v0", sex == "female")
  ),
  v4 = list(
    all    = meta_full |> filter(visit == "v4"),
    male   = meta_full |> filter(visit == "v4", sex == "male"),
    female = meta_full |> filter(visit == "v4", sex == "female")
  )
)
age_var_by_timepoint <- c(v0 = "age_v0", v4 = "age")

# Change outcomes: delta_ffmi_v4, perc_change_ffmi_v4 (+ adj_fmi variants).
# Long format (v0 + v4) so the LMM can test a baseline-abundance x visit
# interaction, as in 4b_lmm_mb_ffmi_trajectories.R.
meta_long <- meta_full |>
  filter(visit %in% c("v0", "v4"), !is.na(ffmi)) |>
  mutate(
    visit = factor(visit, levels = c("v0", "v4")),
    age_centered_v0 = age_v0 - mean(age_v0, na.rm = TRUE)
  ) |>
  left_join(baseline_fmi, by = "id")
meta_long_by_group <- list(
  all    = meta_long,
  male   = meta_long |> filter(sex == "male"),
  female = meta_long |> filter(sex == "female")
)

# Baseline (v0) species abundance (same filtering as 6a_ml_prep.R), wide with `id`
mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb_mat <- t(as(mb@otu_table, "matrix"))
tk <- apply(mb_mat, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
species_wide <- as.data.frame(mb_mat[, tk]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v0$")) |>
  select(-sampleid)

# v4 species abundance (same filtering as 6a_ml_prep.R's mb2_v4), wide with
# `id` - needed because ffmi_v4/fmi_v4's top features come from the v4
# species table, not the baseline one
mb_v4 <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
mb_v4 <- phyloseq::prune_samples(phyloseq::sample_data(mb_v4)$visit == "v4", mb_v4)
mb_v4_mat <- t(as(mb_v4@otu_table, "matrix"))
tk_v4 <- apply(mb_v4_mat, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
species_wide_v4 <- as.data.frame(mb_v4_mat[, tk_v4]) |>
  rownames_to_column("sampleid") |>
  mutate(id = sampleid |> str_remove("^BARIA_") |> str_remove("_v4$")) |>
  select(-sampleid)

species_wide_by_timepoint <- list(v0 = species_wide, v4 = species_wide_v4)

#### Process each outcome ####
for (i in seq_len(nrow(model_defs))) {
  outcome      <- model_defs$outcome[i]
  label        <- model_defs$label[i]
  is_change    <- model_defs$is_change[i]
  response_var <- model_defs$response_var[i]
  timepoint    <- model_defs$timepoint[i]
  is_adj       <- model_defs$is_adj[i]
  base <- file.path("results/mlmodels", outcome)
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  cat("\n====", label, "====\n")

  ## Aggregate metrics + explained variance summary across subgroups
  metrics <- map_dfr(groups, function(g) {
    m <- get_metrics_reg(file.path(base, g), paste0(outcome, "_", g))
    if (is.null(m)) return(NULL)
    m |> mutate(subgroup = g)
  })
  if (nrow(metrics) > 0) {
    write.csv(metrics, file.path(base, paste0(outcome, "_metrics.csv")), row.names = FALSE)
    pl <- ggplot(metrics, aes(x = subgroup, y = Median.Explained.Variance * 100)) +
      geom_point(size = 2) +
      geom_segment(aes(xend = subgroup, yend = 0)) +
      coord_flip() +
      labs(title = label, x = "", y = "Explained variance (%)") +
      theme_Publication()
    ggsave(file.path(base, paste0(outcome, "_explained_variance.pdf")), pl, width = 4, height = 3)
    print(metrics)
  }

  ## Per-subgroup: feature importance, correlation sanity checks, forest plot
  for (g in groups) {
    model_path <- file.path(base, g)
    fi <- get_feature_importance(model_path, paste0(outcome, "_", g))
    if (is.null(fi)) {
      cat(sprintf("  %-12s no XGBeast output found, skipping\n", g))
      next
    }
    fi_top <- top_features(fi, n = 15)
    feats <- fi_top$FeatName
    subtitle <- paste0(label, " (", g, ")")
    input <- read_input_data(model_path)

    dir.create(model_path, recursive = TRUE, showWarnings = FALSE)
    p_imp <- plot_feature_importance(fi_top, title = subtitle)
    ggsave(file.path(model_path, "feature_importance_top15.pdf"), p_imp, width = 7, height = 5)

    p_cor <- plot_feature_correlations(input$X, input$y, feats, title = subtitle)
    ggsave(file.path(model_path, "feature_correlations_top15.pdf"), p_cor, width = 12, height = 9)

    covariates <- if (g == "all") "sex" else character(0)
    forest_title <- subtitle
    if (is_change) {
      covariates <- c("age_centered_v0", "perc_change_weight_kg", covariates)
      if (is_adj) covariates <- c(covariates, "fmi_v0")
      # delta_ffmi_v4: the LMM this confirms whether baseline microbiome
      # predicts the FFMI trajectory, adjusted for concurrent FMI
      if (outcome == "delta_ffmi_v4") {
        covariates <- c(covariates, "fmi")
        forest_title <- paste0("FFMI over 1 year (baseline→v4) ~ baseline microbiome",
                                ", LMM, adj. FMI (", g, ")")
      }
      forest <- run_lmm_forest(species_wide, meta_long_by_group[[g]], feats, "ffmi", "v4", covariates)
    } else {
      covariates <- c(age_var_by_timepoint[[timepoint]], covariates)
      if (is_adj) covariates <- c(covariates, "fmi")
      # Restrict to the exact subjects XGBeast trained on for this model (for matched subj model)
      ids_used <- input$subject_ids |> str_remove("^BARIA_") |> str_remove(paste0("_", timepoint, "$"))
      forest <- run_lm_forest(
        species_wide_by_timepoint[[timepoint]] |>
          filter(id %in% ids_used) |>
          inner_join(meta_cs_by_timepoint_group[[timepoint]][[g]], by = "id"),
        feats, response_var, covariates
      )
    }
    p_forest <- plot_forest(forest, title = forest_title)
    ggsave(file.path(model_path, "forest_plot_top15.pdf"), p_forest, width = 7, height = 6)

    cat(sprintf("  %-12s %d/%d features FDR < 0.05 (%s)\n", g,
                sum(forest$p_fdr < 0.05), nrow(forest),
                if (is_change) "LMM" else "LM"))
  }
}

### Explained variance violin ###
violin_out_dir <- "results/graphs/ml_explained_variance"
dir.create(violin_out_dir, recursive = TRUE, showWarnings = FALSE)

violin_defs <- tribble(
  ~outcome,                       ~violin_label,                   ~violin_group,
  "ffmi",                         "FFMI (v0)",                     "FFMI",
  "ffmi_adj_fmi",                 "FFMI (v0)\nadj. FMI",           "FFMI",
  "delta_ffmi_v4",                "Delta FFMI (1y)",               "Delta FFMI",
  "delta_ffmi_v4_adj_fmi",        "Delta FFMI (1y)\nadj. FMI",     "Delta FFMI",
  "perc_change_ffmi_v4",          "%change FFMI (1y)",             "%Delta FFMI",
  "perc_change_ffmi_v4_adj_fmi",  "%change FFMI (1y)\nadj. FMI",   "%Delta FFMI",
) |>
  left_join(model_defs |> select(outcome, is_adj), by = "outcome")

# color palette
violin_fill_colors <- c(
  "FFMI.FALSE"        = "#8ec6f2",
  "FFMI.TRUE"          = "#2a78d6",
  "Delta FFMI.FALSE"   = "#a8dfc4",
  "Delta FFMI.TRUE"    = "#1baf7a",
  "%Delta FFMI.FALSE"  = "#f5c396",
  "%Delta FFMI.TRUE"   = "#eb6834"
)

ev_df <- pmap_dfr(violin_defs, function(outcome, violin_label, violin_group, is_adj) {
  it <- get_iterations_reg(file.path("results/mlmodels", outcome, "all"), paste0(outcome, "_all"))
  if (is.null(it)) {
    cat("  no XGBeast output found for", outcome, "(all), skipping violin\n")
    return(NULL)
  }
  tibble(label = violin_label, group = violin_group, adj = is_adj, explained_variance = it$Explained.Variance * 100)
})

if (nrow(ev_df) == 0) {
  cat("  no models found, skipping explained variance violin plot\n")
} else {
  ev_df <- ev_df |>
    mutate(
      # reversed so, after coord_flip(), FFMI (v0) reads at the top
      label = factor(label, levels = rev(violin_defs$violin_label)),
      fill_key = paste(group, adj, sep = ".")
    )

  ev_medians <- ev_df |>
    summarise(med = median(explained_variance), max_val = max(explained_variance), .by = label) |>
    mutate(label_text = sprintf("%.1f%%", med))

  p_violin <- ggplot(ev_df, aes(x = label, y = explained_variance, fill = fill_key)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_violin(trim = FALSE, alpha = 0.9, color = "grey30", linewidth = 0.3) +
    geom_boxplot(width = 0.08, outlier.shape = NA, color = "grey20", fill = "white") +
    geom_text(
      data = ev_medians, aes(x = label, y = max_val + 2, label = label_text),
      inherit.aes = FALSE, position = position_nudge(x = 0.32),
      hjust = 0, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = violin_fill_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    coord_flip() +
    labs(
      title = "Explained variance across CV iterations - FFMI models, species (all subjects)",
      x = "", y = "Explained variance (%)"
    ) +
    theme_Publication() +
    theme(axis.text.y = element_text(size = 9))

  ggsave(file.path(violin_out_dir, "ffmi_models_explained_variance_violin.pdf"), p_violin, width = 9, height = 6)
  cat("Saved", file.path(violin_out_dir, "ffmi_models_explained_variance_violin.pdf"), "\n")
}
