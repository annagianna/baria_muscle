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
