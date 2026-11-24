# Liver RNA-seq MitoCarta genes DE- %FFMI loss groups at 1y (exploratory)
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(readxl)
library(MetBrewer)
source("scripts/assets/functions.R")
dir.create("results/graphs/RNAseq/liver/mito", recursive = TRUE, showWarnings = FALSE)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)
renoir_cols_20 <- met.brewer("Renoir", n = 20)
de_cols <- c(
  "Higher in high loss" = renoir_cols_20[18],
  "Higher in modest/low loss" = renoir_cols_20[5],
  "not significant" = "grey80"
)

# Data
de_grp <- read_csv("results/tables/liver_de_grp_all_genes.csv")
de_grp_fmi <- read_csv("results/tables/liver_de_grp_fmi_all_genes.csv")
mito_file <- "data/raw_data/MitoCarta/Human.MitoCarta3.0.xls"

# Check that same genes are present in both DE models
stopifnot(setequal(de_grp$ensembl_gene_id, de_grp_fmi$ensembl_gene_id))
stopifnot(!anyDuplicated(de_grp$ensembl_gene_id))
stopifnot(!anyDuplicated(de_grp_fmi$ensembl_gene_id))

#### MitoCarta3.0 gene inventory ####
excel_sheets(mito_file)
mito_raw <- read_excel(
  mito_file,
  sheet = "A Human MitoCarta3.0",
  na = c("", "NA")
)

# MitoCarta gene symbols
mito_genes <- mito_raw |>
  pull(Symbol) |>
  na.omit() |>
  unique()
stopifnot(length(mito_genes) == 1136)

#### Targeted DE: MitoCarta3.0 genes ####
# Main model: sex + age + batch + %TWL
de_grp_mito <- de_grp |>
  filter(symbol %in% mito_genes) |>
  mutate(
    padj_mito = p.adjust(pvalue, method = "BH"),
    perc_diff = (2^log2FC_mle - 1) * 100,
    direction = case_when(
      is.na(pvalue) ~ NA_character_,
      padj_mito < 0.05 & log2FC_mle > 0 ~ "Higher in high loss",
      padj_mito < 0.05 & log2FC_mle < 0 ~ "Higher in modest/low loss",
      TRUE ~ "not significant"
    )
  ) |>
  arrange(pvalue)

de_sig_mito <- de_grp_mito |>
  filter(padj_mito < 0.05)

# Sensitivity model: additional adjustment for baseline FMI
de_grp_fmi_mito <- de_grp_fmi |>
  inner_join(
    de_grp_mito |>
      select(ensembl_gene_id, symbol, description),
    by = "ensembl_gene_id"
  ) |>
  mutate(
    padj_mito_fmi = p.adjust(pvalue_fmi, method = "BH"),
    perc_diff_fmi = (2^log2FC_mle_fmi - 1) * 100,
    direction = case_when(
      is.na(pvalue_fmi) ~ NA_character_,
      padj_mito_fmi < 0.05 & log2FC_mle_fmi > 0 ~ "Higher in high loss",
      padj_mito_fmi < 0.05 & log2FC_mle_fmi < 0 ~ "Higher in modest/low loss",
      TRUE ~ "not significant"
    )
  ) |>
  arrange(pvalue_fmi)

de_sig_fmi_mito <- de_grp_fmi_mito |>
  filter(padj_mito_fmi < 0.05)

#### Comparison: main vs. FMI-adjusted DE ####
de_mito_compare <- de_grp_mito |>
  select(
    ensembl_gene_id, symbol, description,
    baseMean, log2FC_mle, log2FC_shr, perc_diff,
    pvalue, padj, padj_mito
  ) |>
  left_join(
    de_grp_fmi_mito |>
      select(
        ensembl_gene_id,
        log2FC_mle_fmi, perc_diff_fmi,
        pvalue_fmi, padj_fmi, padj_mito_fmi
      ),
    by = "ensembl_gene_id"
  ) |>
  mutate(
    same_direction = sign(log2FC_mle) == sign(log2FC_mle_fmi),
    sig_main = padj_mito < 0.05,
    sig_fmi = padj_mito_fmi < 0.05,
    significance = case_when(
      is.na(sig_main) | is.na(sig_fmi) ~ "Missing test result",
      sig_main & sig_fmi ~ "Both models",
      sig_main ~ "Main model only",
      sig_fmi ~ "FMI-adjusted only",
      TRUE ~ "Neither model"
    )
  ) |>
  arrange(padj_mito)

#### Checks ####
# MitoCarta coverage
length(mito_genes)
n_distinct(de_grp_mito$symbol)

# Duplicate symbols
de_grp_mito |> count(symbol) |> filter(n > 1)

# Significant genes
table(de_grp_mito$direction, useNA = "ifany")
table(de_grp_fmi_mito$direction, useNA = "ifany")
table(de_mito_compare$significance)

de_sig_mito
de_sig_fmi_mito