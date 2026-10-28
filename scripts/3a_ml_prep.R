# Prep XGBeast input data: microbiome species ~ 1-year %FFMI loss group
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

# Add `new_var` = sex-specific lowest tertile of `var` ("high") vs. the rest
# ("modest/low"), same definition as perc_change_ffmi_v4_group in 0a.
add_tertile_group <- function(data, var, new_var) {
  data |>
    group_by(sex) |>
    mutate(!!new_var := if_else(
      .data[[var]] <= quantile(.data[[var]], probs = 1/3, na.rm = TRUE),
      "high", "modest/low"
    )) |>
    ungroup()
}

meta <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
meta$sampleid <- str_c("BARIA_", meta$id, "_v0") # match against baseline mb sample ids
dim(meta)

mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb <- t(as(mb@otu_table, "matrix"))
tk <- apply(mb, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
mb2 <- as.data.frame(mb[, tk])
dim(mb2) # 147 species left, 365 subjects
mb2[1:5,1:5]

#### 1-year %FFMI loss group (high vs. modest/low) ####
# Outcome: perc_change_ffmi_v4_group from 0a (lowest sex-specific tertile of
# %FFMI change v0 -> v4 = "high" loss, cut-off from the full clinical cohort),
# coded 1 = high loss; predictors: baseline (v0) species (mb2)
cat("perc_change_ffmi_v4_group: ", sum(!is.na(meta[["perc_change_ffmi_v4_group"]])), "/", nrow(meta), "non-missing\n")

meta_group_var <- meta |> filter(!is.na(perc_change_ffmi_v4_group))
meta_group_var_male <- meta_group_var |> filter(sex == "male")
meta_group_var_female <- meta_group_var |> filter(sex == "female")

base_group <- file.path("results/mlmodels/perc_change_ffmi_v4_group")
build_input_data(meta_group_var,        mb2, "perc_change_ffmi_v4_group", file.path(base_group, "all"),    mode = "class", pos_class = "high")
build_input_data(meta_group_var_male,   mb2, "perc_change_ffmi_v4_group", file.path(base_group, "male"),   mode = "class", pos_class = "high")
build_input_data(meta_group_var_female, mb2, "perc_change_ffmi_v4_group", file.path(base_group, "female"), mode = "class", pos_class = "high")

#### 1-year %FFMI loss group, corrected for baseline FMI ####
# %FFMI change residualized on baseline FMI within each sex, then dichotomized
# at the lowest sex-specific tertile of the residual (full clinical cohort, as
# above); one definition per subject, shared by the all/male/female models
meta_group_var_fmi <- meta |>
  filter(!is.na(perc_change_ffmi_v4), !is.na(fmi_v0)) |>
  group_split(sex) |>
  lapply(add_residual_var, "perc_change_ffmi_v4", "fmi_v0", "perc_change_ffmi_v4_adj_fmi") |>
  bind_rows() |>
  add_tertile_group("perc_change_ffmi_v4_adj_fmi", "perc_change_ffmi_v4_group_adj_fmi")
meta_group_var_fmi_male <- meta_group_var_fmi |> filter(sex == "male")
meta_group_var_fmi_female <- meta_group_var_fmi |> filter(sex == "female")

# agreement with the unadjusted grouping
print(table(unadjusted = meta_group_var_fmi$perc_change_ffmi_v4_group,
            adj_fmi = meta_group_var_fmi$perc_change_ffmi_v4_group_adj_fmi))

base_group_adj <- file.path("results/mlmodels/perc_change_ffmi_v4_group_adj_fmi")
build_input_data(meta_group_var_fmi,        mb2, "perc_change_ffmi_v4_group_adj_fmi", file.path(base_group_adj, "all"),    mode = "class", pos_class = "high")
build_input_data(meta_group_var_fmi_male,   mb2, "perc_change_ffmi_v4_group_adj_fmi", file.path(base_group_adj, "male"),   mode = "class", pos_class = "high")
build_input_data(meta_group_var_fmi_female, mb2, "perc_change_ffmi_v4_group_adj_fmi", file.path(base_group_adj, "female"), mode = "class", pos_class = "high")
