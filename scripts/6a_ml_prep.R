# Prep XGBeast input data: microbiome species ~ FFMI
# Barbara Verhaar

library(stringr)
library(tidyr)
library(tibble)
library(dplyr)

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

meta <- readRDS("data/processed_data/BARIA_muscle_long.RDS")
meta$sampleid <- str_c("BARIA_", meta$id, "_", meta$visit)
dim(meta)
names(meta)
cat("FFMI: ", sum(!is.na(meta[["ffmi"]])), "/", nrow(meta), "non-missing\n")

meta <- meta |> filter(visit == "v0")
baseline_fmi <- meta |> select(id, fmi_v0 = fmi)

meta_var <- meta |> filter(!is.na(ffmi))
meta_var_male <- meta_var |> filter(sex == "male")
meta_var_female <- meta_var |> filter(sex == "female")

mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb <- t(as(mb@otu_table, "matrix"))
tk <- apply(mb, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
mb2 <- as.data.frame(mb[, tk])
dim(mb2) # 147 species left, 365 subjects
mb2[1:5,1:5]

# v4 species table, for outcomes measured at v4 that should be predicted
# from the microbiome at the same visit rather than the baseline microbiome
mb_v4 <- readRDS("data/processed_data/BARIA_mb_clean.RDS")
mb_v4 <- phyloseq::prune_samples(phyloseq::sample_data(mb_v4)$visit == "v4", mb_v4)
mb_v4 <- t(as(mb_v4@otu_table, "matrix"))
tk_v4 <- apply(mb_v4, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
mb2_v4 <- as.data.frame(mb_v4[, tk_v4])
dim(mb2_v4)
mb2_v4[1:5,1:5]

base <- file.path("results/mlmodels/ffmi")
build_input_data(meta_var,        mb2, "ffmi", file.path(base, "all"),    mode = "reg")
build_input_data(meta_var_male,   mb2, "ffmi", file.path(base, "male"),   mode = "reg")
build_input_data(meta_var_female, mb2, "ffmi", file.path(base, "female"), mode = "reg")

#### FFMI corrected for FMI ####
# Outcome: FFMI residualized on baseline FMI (fat mass index), i.e. muscle
# mass independent of adiposity; predictors: baseline (v0) species (mb2)
meta_var_fmi <- meta_var |> filter(!is.na(fmi)) |>
  add_residual_var("ffmi", "fmi", "ffmi_adj_fmi")
meta_var_fmi_male <- meta_var_male |> filter(!is.na(fmi)) |>
  add_residual_var("ffmi", "fmi", "ffmi_adj_fmi")
meta_var_fmi_female <- meta_var_female |> filter(!is.na(fmi)) |>
  add_residual_var("ffmi", "fmi", "ffmi_adj_fmi")

base_adj <- file.path("results/mlmodels/ffmi_adj_fmi")
build_input_data(meta_var_fmi,        mb2, "ffmi_adj_fmi", file.path(base_adj, "all"),    mode = "reg")
build_input_data(meta_var_fmi_male,   mb2, "ffmi_adj_fmi", file.path(base_adj, "male"),   mode = "reg")
build_input_data(meta_var_fmi_female, mb2, "ffmi_adj_fmi", file.path(base_adj, "female"), mode = "reg")

#### FMI at v0 ####
# Outcome: baseline FMI (fat mass index); predictors: baseline (v0) species
# (mb2). Alongside ffmi/ffmi_v4/fmi_v4, lets us compare how well the baseline
# microbiome predicts muscle vs. fat mass, before vs. after bariatric surgery.
meta_fmi_var <- meta |> filter(!is.na(fmi))
meta_fmi_var_male <- meta_fmi_var |> filter(sex == "male")
meta_fmi_var_female <- meta_fmi_var |> filter(sex == "female")

base_fmi_v0 <- file.path("results/mlmodels/fmi_v0")
build_input_data(meta_fmi_var,        mb2, "fmi", file.path(base_fmi_v0, "all"),    mode = "reg")
build_input_data(meta_fmi_var_male,   mb2, "fmi", file.path(base_fmi_v0, "male"),   mode = "reg")
build_input_data(meta_fmi_var_female, mb2, "fmi", file.path(base_fmi_v0, "female"), mode = "reg")

#### Cross-sectional at v4 (post-surgery) ####
# FFMI/FMI at v4 are predicted from v4 species (mb2_v4); delta/percent-change
# outcomes below are predicted from baseline (v0) species (mb2), since those
# are genuinely prospective (baseline microbiome -> future change).
meta_delta <- readRDS("data/processed_data/BARIA_muscle_long.RDS") |> filter(visit == "v4")
meta_delta$sampleid <- str_c("BARIA_", meta_delta$id, "_v0") # match against baseline mb sample ids
meta_delta$sampleid_v4 <- str_c("BARIA_", meta_delta$id, "_v4") # match against v4 mb sample ids
meta_delta <- meta_delta |> left_join(baseline_fmi, by = "id")

# FFMI at v4
cat("ffmi_v4: ", sum(!is.na(meta_delta[["ffmi"]])), "/", nrow(meta_delta), "non-missing\n")
meta_v4_ffmi_var <- meta_delta |> filter(!is.na(ffmi)) |> mutate(sampleid = sampleid_v4)
meta_v4_ffmi_var_male <- meta_v4_ffmi_var |> filter(sex == "male")
meta_v4_ffmi_var_female <- meta_v4_ffmi_var |> filter(sex == "female")

base_ffmi_v4 <- file.path("results/mlmodels/ffmi_v4")
build_input_data(meta_v4_ffmi_var,        mb2_v4, "ffmi", file.path(base_ffmi_v4, "all"),    mode = "reg")
build_input_data(meta_v4_ffmi_var_male,   mb2_v4, "ffmi", file.path(base_ffmi_v4, "male"),   mode = "reg")
build_input_data(meta_v4_ffmi_var_female, mb2_v4, "ffmi", file.path(base_ffmi_v4, "female"), mode = "reg")

# FMI at v4
cat("fmi_v4: ", sum(!is.na(meta_delta[["fmi"]])), "/", nrow(meta_delta), "non-missing\n")
meta_v4_fmi_var <- meta_delta |> filter(!is.na(fmi)) |> mutate(sampleid = sampleid_v4)
meta_v4_fmi_var_male <- meta_v4_fmi_var |> filter(sex == "male")
meta_v4_fmi_var_female <- meta_v4_fmi_var |> filter(sex == "female")

base_fmi_v4 <- file.path("results/mlmodels/fmi_v4")
build_input_data(meta_v4_fmi_var,        mb2_v4, "fmi", file.path(base_fmi_v4, "all"),    mode = "reg")
build_input_data(meta_v4_fmi_var_male,   mb2_v4, "fmi", file.path(base_fmi_v4, "male"),   mode = "reg")
build_input_data(meta_v4_fmi_var_female, mb2_v4, "fmi", file.path(base_fmi_v4, "female"), mode = "reg")

#### FFMI/FMI at v0, restricted to the v4-microbiome subset (matched) ####
# Same outcome & v0 species (mb2) as the "ffmi"/"fmi_v0" models above, but
# limited to the subjects who also have v4 microbiome data. Comparing these
# against ffmi_v4/fmi_v4 isolates whether those results differ from the
# full-cohort v0 models because of the subject subset or because of the
# predictor/outcome timepoint.
ffmi_v4_ids <- meta_v4_ffmi_var$id[meta_v4_ffmi_var$sampleid %in% rownames(mb2_v4)]
fmi_v4_ids <- meta_v4_fmi_var$id[meta_v4_fmi_var$sampleid %in% rownames(mb2_v4)]

meta_var_v0matched <- meta_var |> filter(id %in% ffmi_v4_ids)
meta_var_v0matched_male <- meta_var_v0matched |> filter(sex == "male")
meta_var_v0matched_female <- meta_var_v0matched |> filter(sex == "female")

base_ffmi_v0matched <- file.path("results/mlmodels/ffmi_v0_matched")
build_input_data(meta_var_v0matched,        mb2, "ffmi", file.path(base_ffmi_v0matched, "all"),    mode = "reg")
build_input_data(meta_var_v0matched_male,   mb2, "ffmi", file.path(base_ffmi_v0matched, "male"),   mode = "reg")
build_input_data(meta_var_v0matched_female, mb2, "ffmi", file.path(base_ffmi_v0matched, "female"), mode = "reg")

meta_fmi_var_v0matched <- meta_fmi_var |> filter(id %in% fmi_v4_ids)
meta_fmi_var_v0matched_male <- meta_fmi_var_v0matched |> filter(sex == "male")
meta_fmi_var_v0matched_female <- meta_fmi_var_v0matched |> filter(sex == "female")

base_fmi_v0matched <- file.path("results/mlmodels/fmi_v0_matched")
build_input_data(meta_fmi_var_v0matched,        mb2, "fmi", file.path(base_fmi_v0matched, "all"),    mode = "reg")
build_input_data(meta_fmi_var_v0matched_male,   mb2, "fmi", file.path(base_fmi_v0matched, "male"),   mode = "reg")
build_input_data(meta_fmi_var_v0matched_female, mb2, "fmi", file.path(base_fmi_v0matched, "female"), mode = "reg")

#### Delta FFMI v0 -> v4 ####
# Outcome: change in FFMI from baseline to v4; predictors: baseline (v0) species (mb2)
cat("delta_ffmi_v4: ", sum(!is.na(meta_delta[["delta_ffmi"]])), "/", nrow(meta_delta), "non-missing\n")

meta_delta_var <- meta_delta |> filter(!is.na(delta_ffmi))
meta_delta_var_male <- meta_delta_var |> filter(sex == "male")
meta_delta_var_female <- meta_delta_var |> filter(sex == "female")

base_delta <- file.path("results/mlmodels/delta_ffmi_v4")
build_input_data(meta_delta_var,        mb2, "delta_ffmi", file.path(base_delta, "all"),    mode = "reg")
build_input_data(meta_delta_var_male,   mb2, "delta_ffmi", file.path(base_delta, "male"),   mode = "reg")
build_input_data(meta_delta_var_female, mb2, "delta_ffmi", file.path(base_delta, "female"), mode = "reg")

#### Delta FFMI corrected for FMI ####
# Outcome: delta FFMI residualized on baseline FMI; predictors: baseline (v0) species (mb2)
meta_delta_var_fmi <- meta_delta_var |> filter(!is.na(fmi_v0)) |>
  add_residual_var("delta_ffmi", "fmi_v0", "delta_ffmi_adj_fmi")
meta_delta_var_fmi_male <- meta_delta_var_male |> filter(!is.na(fmi_v0)) |>
  add_residual_var("delta_ffmi", "fmi_v0", "delta_ffmi_adj_fmi")
meta_delta_var_fmi_female <- meta_delta_var_female |> filter(!is.na(fmi_v0)) |>
  add_residual_var("delta_ffmi", "fmi_v0", "delta_ffmi_adj_fmi")

base_delta_adj <- file.path("results/mlmodels/delta_ffmi_v4_adj_fmi")
build_input_data(meta_delta_var_fmi,        mb2, "delta_ffmi_adj_fmi", file.path(base_delta_adj, "all"),    mode = "reg")
build_input_data(meta_delta_var_fmi_male,   mb2, "delta_ffmi_adj_fmi", file.path(base_delta_adj, "male"),   mode = "reg")
build_input_data(meta_delta_var_fmi_female, mb2, "delta_ffmi_adj_fmi", file.path(base_delta_adj, "female"), mode = "reg")

#### Percentual change FFMI v0 -> v4 ####
# Outcome: percentual change in FFMI from baseline to v4; predictors: baseline (v0) species (mb2)
cat("perc_change_ffmi_v4: ", sum(!is.na(meta_delta[["perc_change_ffmi"]])), "/", nrow(meta_delta), "non-missing\n")

meta_pchange_var <- meta_delta |> filter(!is.na(perc_change_ffmi))
meta_pchange_var_male <- meta_pchange_var |> filter(sex == "male")
meta_pchange_var_female <- meta_pchange_var |> filter(sex == "female")

base_pchange <- file.path("results/mlmodels/perc_change_ffmi_v4")
build_input_data(meta_pchange_var,        mb2, "perc_change_ffmi", file.path(base_pchange, "all"),    mode = "reg")
build_input_data(meta_pchange_var_male,   mb2, "perc_change_ffmi", file.path(base_pchange, "male"),   mode = "reg")
build_input_data(meta_pchange_var_female, mb2, "perc_change_ffmi", file.path(base_pchange, "female"), mode = "reg")

#### Percentual change FFMI corrected for FMI ####
# Outcome: percentual change FFMI residualized on baseline FMI; predictors: baseline (v0) species (mb2)
meta_pchange_var_fmi <- meta_pchange_var |> filter(!is.na(fmi_v0)) |>
  add_residual_var("perc_change_ffmi", "fmi_v0", "perc_change_ffmi_adj_fmi")
meta_pchange_var_fmi_male <- meta_pchange_var_male |> filter(!is.na(fmi_v0)) |>
  add_residual_var("perc_change_ffmi", "fmi_v0", "perc_change_ffmi_adj_fmi")
meta_pchange_var_fmi_female <- meta_pchange_var_female |> filter(!is.na(fmi_v0)) |>
  add_residual_var("perc_change_ffmi", "fmi_v0", "perc_change_ffmi_adj_fmi")

base_pchange_adj <- file.path("results/mlmodels/perc_change_ffmi_v4_adj_fmi")
build_input_data(meta_pchange_var_fmi,        mb2, "perc_change_ffmi_adj_fmi", file.path(base_pchange_adj, "all"),    mode = "reg")
build_input_data(meta_pchange_var_fmi_male,   mb2, "perc_change_ffmi_adj_fmi", file.path(base_pchange_adj, "male"),   mode = "reg")
build_input_data(meta_pchange_var_fmi_female, mb2, "perc_change_ffmi_adj_fmi", file.path(base_pchange_adj, "female"), mode = "reg")
