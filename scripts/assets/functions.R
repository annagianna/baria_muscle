# Shared helpers genuinely reused across pipelines: XGBeast output-folder
# discovery, the shared ggplot theme, and species-name/plotting helpers used
# by both 6c_ml_process.R and 7a/7b_mb_*_correlations.R. Anything used by
# only one script lives in that script instead.
#
# Barbara Verhaar

library(dplyr)

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
# model name, e.g. find_output_folder("results/ml_crossectional/phenoage/all", "phenoage_all", "reg").
# If the model has been rerun, multiple timestamped folders can exist for the
# same name; since the timestamp suffix (YYYY_MM_DD__HH-MM-SS) sorts
# lexicographically in chronological order, take the last one so a rerun
# always wins over a stale one.
find_output_folder <- function(base_path, name, mode = c("reg", "class")) {
  mode <- match.arg(mode)
  if (!dir.exists(base_path)) return(NA_character_)
  prefix <- paste0("output_XGB_", mode, "_", name)
  li <- list.files(base_path)
  hit <- sort(li[startsWith(li, prefix) & !grepl("PERMUTED", li)])
  if (length(hit) == 0) return(NA_character_)
  file.path(base_path, hit[length(hit)])
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

# Short, human-readable species label from a "Genus_species_SGB####[_group]"
# name, e.g. "Faecalibacterium prausnitzii". The SGB id is meaningless to
# most readers so it's dropped; make.unique() disambiguates the rare case
# where several SGBs share the same genus_species name, instead of silently
# collapsing distinct features onto the same plot row/facet.
species_label <- function(x) {
  base <- str_remove(x, "_SGB\\d+(_group)?$") |> str_replace_all("_", " ")
  make.unique(base, sep = " ")
}
