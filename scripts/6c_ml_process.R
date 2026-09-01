# Process XGBeast models and make plots for each model
# Barbara Verhaar

library(tidyverse)
library(phyloseq)
library(lmerTest)
library(broom.mixed)
library(broom)
source("scripts/assets/functions.R")

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

#### Data shared by the confirmatory LM/LMM models (not the ML input data
#### itself, which is read per-model via read_input_data()) ####

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
