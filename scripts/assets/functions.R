
## Shared helpers reused across pipelines:
# XGBeast output-folder discovery, shared ggplot themes,
# top species/metabolite loaders, and species/metabolite-name helpers
# used by 3c_ml_process.R, 4a–4c and 5a/5b.
# Anything used by only one script lives in that script instead.

# Barbara Verhaar, Anna Giannakogeorgou

#### Packages ####
library(dplyr)

#### Themes ####
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

theme_minimal_custom <- function(base_size = 14, base_family = "sans") {

  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(0.8), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(0.8)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(colour = "black"),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major = element_line(
        colour = "#dddddd",
        linewidth = 0.4,
        linetype = "22"
      ),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "bottom",
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}

theme_minimal_composition <- function(base_size = 14, base_family = "sans") {

  theme_minimal(base_size = base_size, base_family = base_family) +

    theme(
      plot.title = element_text(face = "bold", size = rel(1), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      panel.grid.major.y = element_line(colour = "#dddddd", linewidth = 0.4, linetype = "22"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.text = element_text(size = rel(0.7)),
      plot.margin = unit(c(10, 5, 5, 5), "mm")
    )
}

#### XGBeast ####
# Find the latest non-permuted XGBeast output folder
find_output_folder <- function(base_path, name, mode = c("reg", "class")) {
  mode <- match.arg(mode)
  if (!dir.exists(base_path)) return(NA_character_)

  prefix <- paste0("output_XGB_", mode, "_", name)
  li <- list.files(base_path)
  hit <- sort(li[startsWith(li, prefix) & !grepl("PERMUTED", li)])

  if (length(hit) == 0) return(NA_character_)
  file.path(base_path, hit[length(hit)])
}

# Read feature importance from XGBeast output
get_feature_importance <- function(base_path, name, mode = c("reg", "class")) {
  mode <- match.arg(mode)
  folder <- find_output_folder(base_path, name, mode)
  if (is.na(folder)) return(NULL)

  f <- file.path(folder, "feature_importance.txt")
  if (!file.exists(f)) return(NULL)

  read.delim(f)
}

# Top n real features by relative importance
top_features <- function(feature_importance, n = 15) {
  feature_importance |>
    filter(!FeatName %in% c("random_variable1", "random_variable2")) |>
    arrange(desc(RelFeatImp)) |>
    head(n)
}

#### ML-selected features ####

# Top-15 species from dichotomous ML model (high vs. modest/low %FFMI loss at 1y)
# Direction from logistic forest: OR > 1 = higher in high loss
get_top_species <- function(path = "results/mlmodels/perc_change_ffmi_v4_group/all") {

  fi <- get_feature_importance(path, "perc_change_ffmi_v4_group_all", "class")
  forest <- read.csv(file.path(path, "forest_results_top15.csv"))

  tibble(species = top_features(fi, n = 15)$FeatName) |>
    mutate(species_label = species_label(species)) |>
    left_join(
      forest |>
        select(species, estimate, p_fdr),
      by = "species"
    ) |>
    mutate(estimate_direction = if_else(estimate > 1, "Higher in high loss", "Higher in modest/low loss"))
}

# Top metabolites from dichotomous ML model (high vs. modest/low %FFMI loss at 1y)
# No effect direction yet: metabolite ML-processing has not been run
get_top_metabolites <- function(path = "results/mlmodels_metab/perc_change_ffmi_v4_group/all", n = 15) {

  fi <- get_feature_importance(path, "perc_change_ffmi_v4_group_all", "class")

  tibble(metabolite = top_features(fi, n = n)$FeatName) |>
    left_join(
      fi |>
        select(FeatName, RelFeatImp),
      by = c("metabolite" = "FeatName")
    )
}

#### Feature labels ####
# Short, readable species label from a "Genus_species_SGB####[_group]" name.
# make.unique() disambiguates species with identical genus/species labels.
species_label <- function(x) {
  base <- str_remove(x, "_SGB\\d+(_group)?$") |> str_replace_all("_", " ")
  make.unique(base, sep = " ")
}

# Shorten long sphingolipid names listing multiple isomers
metabolite_label <- function(x) {
  str_replace(
    x,
    "(\\(d\\d+:\\d+/\\d+:\\d+(?:\\([^()]*\\))?)(?:,\\s*d\\d+:\\d+/\\d+:\\d+(?:\\([^()]*\\))?)+(\\))",
    "\\1, ...\\2"
  )
}

#### Gene annotation ####
# Human gene annotation from org.Hs.eg.db
# Called with :: so AnnotationDbi::select() -> does not mask dplyr::select().
get_gene_annotations <- function() {
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  AnnotationDbi::select(
    orgdb,
    keys = AnnotationDbi::keys(orgdb, keytype = "ENSEMBL"),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENENAME", "GENETYPE")
  ) |>
    dplyr::as_tibble() |>
    dplyr::rename(ensgene = ENSEMBL, symbol = SYMBOL, description = GENENAME, biotype = GENETYPE) |>
    dplyr::distinct()
}