# Shared helpers for the 4_ml_crossectional pipeline (batch2 microbiome ~
# aging-variable models): prep (write_y/write_data/build_input_data),
# model-output discovery, and result plotting.
#
# Sourced by every script in 1_prep/ and 3_process/.
#
# Barbara Verhaar

library(dplyr)

# ── prep helpers ────────────────────────────────────────────────────────────

write_y <- function(x, name_y, data_path) {
  if (!name_y %in% c("y_binary.txt", "y_reg.txt")) {
    stop('name_y must be "y_binary.txt" or "y_reg.txt"')
  }
  if (any(is.na(x))) stop("There are missing values in the outcome data!")
  data_path <- file.path(data_path, "input_data")
  dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = file.path(data_path, name_y),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

write_data <- function(x, data_path) {
  x <- as.matrix(x)
  if (any(is.na(x))) stop("There are missing values in the input data!")
  data_path <- file.path(data_path, "input_data")
  dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
  write.table(x, file.path(data_path, "X_data.txt"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  write.table(colnames(x), file.path(data_path, "feat_ids.txt"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  write.table(rownames(x), file.path(data_path, "subject_ids.txt"),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}


# Build one XGBeast input_data folder: intersect a metadata subgroup with the
# species table on sampleID_mb, write X_data/feat_ids/subject_ids + y_reg.txt
# (mode = "reg") or y_binary.txt coded 1 = pos_class, 0 = other (mode = "class").
build_input_data <- function(meta_sub, mb2, var, out_path,
                              mode = c("reg", "class"), pos_class = NULL) {
  mode <- match.arg(mode)
  X <- mb2[rownames(mb2) %in% meta_sub$sampleid, , drop = FALSE]
  meta_sub <- meta_sub[match(rownames(X), meta_sub$sampleid), ]
  stopifnot(all(meta_sub$sampleid == rownames(X)))

  if (mode == "reg") {
    y <- as.data.frame(as.numeric(meta_sub[[var]]))
    y_name <- "y_reg.txt"
  } else {
    y <- as.data.frame(as.integer(meta_sub[[var]] == pos_class))
    y_name <- "y_binary.txt"
  }

  write_data(X, out_path)
  write_y(y, name_y = y_name, out_path)
  cat(sprintf("  %-55s %4d samples, %4d taxa\n", out_path, nrow(X), ncol(X)))
  invisible(nrow(X))
}

# Regress `outcome` on `covariate` and add the residuals as column `new_var`,
# e.g. to build an outcome "corrected for" FMI (muscle mass independent of
# adiposity) before handing it to build_input_data().
add_residual_var <- function(data, outcome, covariate, new_var) {
  fit <- lm(reformulate(covariate, outcome), data = data)
  data[[new_var]] <- as.numeric(resid(fit))
  data
}

# ── process helpers ───────────────────────────────────────────────────────────

theme_Publication <- function(base_size = 12, base_family = "sans") {
  library(grid); library(ggthemes)
  theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
      panel.background = element_rect(colour = NA, fill = NA),
      plot.background = element_rect(colour = NA, fill = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(0.8)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.line.y = element_line(colour = "black"),
      axis.title.x = element_text(vjust = -0.2),
      axis.line.x = element_line(colour = "black"),
      axis.ticks.x = element_line(),
      axis.ticks.y = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "right",
      legend.key.size = unit(0.2, "cm"),
      legend.spacing = unit(0, "cm"),
      plot.margin = unit(c(5, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(face = "italic", size = rel(0.6))
    )
}

# Find the (non-PERMUTED) XGBeast output folder for a given subgroup path +
# model name, e.g. find_output_folder("results/ml_crossectional/phenoage/all", "phenoage_all", "reg")
find_output_folder <- function(base_path, name, mode = c("reg", "class")) {
  mode <- match.arg(mode)
  if (!dir.exists(base_path)) return(NA_character_)
  prefix <- paste0("output_XGB_", mode, "_", name)
  li <- list.files(base_path)
  hit <- li[startsWith(li, prefix) & !grepl("PERMUTED", li)]
  if (length(hit) == 0) return(NA_character_)
  file.path(base_path, hit[1])
}

# Read aggregated regression metrics ("Median R2", "Median Explained
# Variance", "Median RMSE", "Median MAE") for one subgroup; NULL if the model
# hasn't been run yet.
get_metrics_reg <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "reg")
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "aggregated_metrics_regression.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Read per-iteration regression metrics ("R2", "Explained Variance", "RMSE",
# "MAE"), one row per CV iteration, for one subgroup; NULL if the model
# hasn't been run yet.
get_iterations_reg <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "reg")
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "model_results_per_iteration.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Read aggregated classification metrics ("Median AUC", "Median Accuracy",
# "Median Precision", "Median Recall", "Median F1-score", "Median Avg.
# Precision Score") for one subgroup; NULL if the model hasn't been run yet.
get_metrics_class <- function(base_path, name) {
  folder <- find_output_folder(base_path, name, "class")
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "aggregated_metrics_classification.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Read XGBeast's per-feature relative importance ("FeatName", "RelFeatImp",
# 0-100, includes random_variable1/2) for one subgroup; NULL if the model
# hasn't been run yet.
get_feature_importance <- function(base_path, name, mode = c("reg", "class")) {
  mode <- match.arg(mode)
  folder <- find_output_folder(base_path, name, mode)
  if (is.na(folder)) return(NULL)
  f <- file.path(folder, "feature_importance.txt")
  if (!file.exists(f)) return(NULL)
  read.delim(f)
}

# Top n real (non-random-variable) features by relative importance.
top_features <- function(feature_importance, n = 15) {
  feature_importance |>
    filter(!FeatName %in% c("random_variable1", "random_variable2")) |>
    arrange(desc(RelFeatImp)) |>
    head(n)
}

# Read the exact X/y XGBeast trained on for one subgroup (input_data/
# X_data.txt, feat_ids.txt, subject_ids.txt, y_reg.txt or y_binary.txt).
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

# Short, human-readable species label from a "Genus_species_SGB####[_group]"
# name. Keeps the SGB id (e.g. "Faecalibacterium prausnitzii (SGB15318)")
# since several SGBs can share the same genus_species name; dropping it
# would silently collapse distinct features onto the same plot row/facet.
species_label <- function(x) {
  sgb <- str_extract(x, "SGB\\d+(_group)?$")
  base <- str_remove(x, "_SGB\\d+(_group)?$") |> str_replace_all("_", " ")
  paste0(base, " (", sgb, ")")
}

# Cross-sectional confirmatory model: for each feature, lm(outcome ~
# log10(abundance) + covariates), returning the log10_abundance term as one
# forest-plot row per feature (species, estimate, conf.low/high, p.value,
# p_fdr). `data` must have one row per subject with columns `outcome`,
# `covariates`, and every name in `feats`.
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

# Change-over-time confirmatory model: for each feature, lmer(outcome ~
# log10(baseline_abundance) * visit + covariates + (1 | id)), returning the
# log10_abundance:visit interaction term as one forest-plot row per feature.
# `species_wide` has one row per subject (column `id`) with baseline
# abundance columns named in `feats`; `long_data` is long-format (v0 +
# follow_up) with columns `id`, `visit`, `outcome`, `covariates`.
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
      label = forcats::fct_reorder(label, estimate)
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

# Sanity-check scatter grid: the outcome actually used by the model (from
# read_input_data()$y) against each top feature's abundance (from
# read_input_data()$X), one panel per feature, annotated with Spearman rho.
plot_feature_correlations <- function(X, y, feats, title = "") {
  df <- as.data.frame(X[, feats, drop = FALSE])
  colnames(df) <- feats
  df$.y <- y
  df |>
    tidyr::pivot_longer(-.y, names_to = "species", values_to = "abundance") |>
    group_by(species) |>
    mutate(rho = cor(abundance, .y, method = "spearman", use = "complete.obs")) |>
    ungroup() |>
    mutate(
      facet_label = sprintf("%s\n(rho = %.2f)", str_wrap(species_label(species), 22), rho),
      facet_label = forcats::fct_reorder(facet_label, -rho)
    ) |>
    ggplot(aes(x = abundance, y = .y)) +
    geom_point(alpha = 0.5, size = 1, color = "#2c7fb8") +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.6) +
    facet_wrap(~facet_label, scales = "free_x", ncol = 5) +
    labs(title = title, x = "Relative abundance (as used in model)", y = "Outcome (as used in model)") +
    theme_Publication() +
    theme(strip.text = element_text(size = 8))
}
