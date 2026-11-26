# Liver RNA-seq GSVA - %FFMI loss groups at 1y
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(DESeq2)
library(GSVA)
library(MetBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggtext)
library(broom)
source("scripts/assets/functions.R")

# Output directories
plot_dir <- "results/graphs/RNAseq/liver/gsva"
table_dir <- "results/tables/RNAseq/liver/gsva"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Theme 
renoir_15 <- met.brewer("Renoir", n = 15)
group_cols <- c("modest_low" = renoir_15[6], "high" = renoir_15[11])
database_cols <- c("Hallmark" = renoir_15[2], "Reactome" = renoir_15[4], "MitoPathways" = renoir_15[7])

# Data
#### GSEA-derived pathway selection ####
gsva_input <- readRDS("results/tables/liver_gsva_input.rds")

