# Liver RNA-seq GSVA - %FFMI loss groups at 1y
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(DESeq2)
library(GSVA)
library(MetBrewer)
source("scripts/assets/functions.R")

# Output directories
plot_dir <- "results/graphs/RNAseq/liver/gsva"
table_dir <- "results/tables/RNAseq/liver/gsva"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mito_category_dir, recursive = TRUE, showWarnings = FALSE)