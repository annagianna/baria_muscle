# Liver RNA-seq GSEA - %FFMI loss groups at 1y
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(fgsea)
library(msigdbr)
library(MetBrewer)
library(readxl)
source("scripts/assets/functions.R")
dir.create("results/graphs/RNAseq/liver/gsea", recursive = TRUE, showWarnings = FALSE)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)
renoir_cols_20 <- met.brewer("Renoir", n = 20)
gsea_cols <- c( "Higher in high loss" = renoir_cols_20[18], "Higher in modest/low loss" = renoir_cols_20[5])

# Data
de_grp <- read_csv("results/tables/liver_de_grp_all_genes.csv")
de_grp_fmi <- read_csv("results/tables/liver_de_grp_fmi_all_genes.csv")
wald_grp <- read_csv("results/tables/liver_de_grp_wald_stats.csv")
wald_grp_fmi <- read_csv("results/tables/liver_de_grp_fmi_wald_stats.csv")
mito_file <- "data/raw_data/MitoCarta/Human.MitoPathways3.0.gmx"
mitocarta_file <- "data/raw_data/MitoCarta/Human.MitoCarta3.0.xls"

# Same genes in both DE models
stopifnot(setequal(wald_grp$ensembl_gene_id, wald_grp_fmi$ensembl_gene_id))

#### Ranked gene statistics ####
# pos = higher expression in high %FFMI loss at 1y; neg = higher expression in modest/low loss
make_ranks <- function(df) {

  df <- df |>
    mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
    filter(is.finite(stat)) |>
    arrange(desc(abs(stat))) |>
    distinct(ensgene, .keep_all = TRUE)

  ranks <- df$stat
  names(ranks) <- df$ensgene

  sort(ranks, decreasing = TRUE)
}

ranks_grp <- make_ranks(wald_grp)
ranks_grp_fmi <- make_ranks(wald_grp_fmi)

stopifnot(!anyDuplicated(names(ranks_grp)), !anyNA(ranks_grp))
stopifnot(!anyDuplicated(names(ranks_grp_fmi)), !anyNA(ranks_grp_fmi))

#### Hallmark and Reactome gene sets ####
hallmark <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = "H")
reactome <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME" )

make_gene_sets <- function(df) {

  df <- df |>
    filter(!is.na(ensembl_gene), ensembl_gene != "") |>
    distinct(gs_name, ensembl_gene)

  split(df$ensembl_gene, df$gs_name)
}

hallmark_sets <- make_gene_sets(hallmark)
reactome_sets <- make_gene_sets(reactome)

#### MitoPathways3.0 gene sets ####
mito_gmx <- read.delim(mito_file, header = FALSE, fill = TRUE, check.names = FALSE, quote = "")
mito_symbols <- setNames(
  lapply(seq_len(ncol(mito_gmx)), function(i) {
    genes <- as.character(mito_gmx[-c(1, 2), i])
    unique(genes[!is.na(genes) & genes != ""])
  }),
  as.character(mito_gmx[1, ])
)

# Map gene symbols to Ensembl IDs
symbol_map <- get_gene_annotations() |>
  select(symbol, ensgene) |>
  filter(!is.na(symbol), symbol != "") |>
  distinct() |>
  group_by(symbol) |>
  filter(n() == 1) |>
  ungroup()

mito_sets <- lapply(mito_symbols, function(genes) {
  ids <- symbol_map$ensgene[match(genes, symbol_map$symbol)]
  unique(ids[!is.na(ids)])
})

# Combine gene sets
gene_sets <- list(
  Hallmark = hallmark_sets,
  Reactome = reactome_sets,
  MitoCarta = mito_sets
)

# All available pathways
all_pathways <- imap_dfr(
  gene_sets,
  ~ tibble(collection = .y, pathway = names(.x))
) |>
  arrange(collection, pathway)

#### Preranked GSEA ####
# BH-FDR within each collection, SEPARATELY for each model
run_gsea <- function(sets, collection, ranks) {

  set.seed(13)

  fgsea::fgseaMultilevel(
    pathways = sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    BPPARAM = BiocParallel::SerialParam()
  ) |>
    as_tibble() |>
    mutate(
      collection = collection,
      direction = if_else(
        NES > 0,
        "Higher in high loss",
        "Higher in modest/low loss"
      ),
      .before = 1
    ) |>
    arrange(pval)
}

#### Hallmark GSEA ####
# Model 1: sex + age + batch + %TWL
gsea_hallmark <- run_gsea(hallmark_sets, "Hallmark", ranks_grp)

# Model 2: adj. baseline FMI
gsea_hallmark_fmi <- run_gsea(hallmark_sets, "Hallmark", ranks_grp_fmi)

#### Reactome GSEA ####
gsea_reactome <- run_gsea(reactome_sets, "Reactome", ranks_grp)
gsea_reactome_fmi <- run_gsea(reactome_sets, "Reactome", ranks_grp_fmi)

#### MitoPathways GSEA ####
gsea_mito <- run_gsea(mito_sets, "MitoCarta", ranks_grp)
gsea_mito_fmi <- run_gsea(mito_sets, "MitoCarta", ranks_grp_fmi)

#### Combine GSEA results ####
gsea_grp <- bind_rows(gsea_hallmark, gsea_reactome, gsea_mito)
gsea_grp_fmi <- bind_rows(gsea_hallmark_fmi, gsea_reactome_fmi, gsea_mito_fmi)

# Significant pathways per model
gsea_grp_sig <- gsea_grp |>
  filter(padj < 0.05)
gsea_grp_fmi_sig <- gsea_grp_fmi |>
  filter(padj < 0.05)

#### Models 1 & 2 (FMI-adj.) comparison ####
compare_gsea <- function(gsea_main, gsea_fmi) {

  gsea_main |>
    select(collection, pathway, size, NES, pval, padj, leadingEdge) |>
    full_join(
      gsea_fmi |>
        select(
          collection, pathway,
          size_fmi = size,
          NES_fmi = NES,
          pval_fmi = pval,
          padj_fmi = padj,
          leadingEdge_fmi = leadingEdge
        ),
      by = c("collection", "pathway")
    ) |>
    mutate(
      same_direction = sign(NES) == sign(NES_fmi),
      sig_main = padj < 0.05,
      sig_fmi = padj_fmi < 0.05,
      significance = case_when(
        is.na(sig_main) | is.na(sig_fmi) ~ "Missing test result",
        sig_main & sig_fmi ~ "Both models",
        sig_main ~ "Main model only",
        sig_fmi ~ "FMI-adjusted only",
        TRUE ~ "Neither model"
      )
    ) |>
    arrange(collection, padj)
}

gsea_compare <- compare_gsea(gsea_grp, gsea_grp_fmi)

# Significant in both models with concordant direction
gsea_robust <- gsea_compare |>
  filter(sig_main, sig_fmi, same_direction)

# Comparison per collection
hallmark_compare <- gsea_compare |> filter(collection == "Hallmark")
reactome_compare <- gsea_compare |> filter(collection == "Reactome")
mito_compare <- gsea_compare |> filter(collection == "MitoCarta")

# Significant and directionally concordant per collection
hallmark_robust <- gsea_robust |> filter(collection == "Hallmark")
reactome_robust <- gsea_robust |> filter(collection == "Reactome")
mito_robust <- gsea_robust |> filter(collection == "MitoCarta")

#### MitoCarta3.0 functional annotations ####
mitocarta <- read_excel(mitocarta_file, sheet = "A Human MitoCarta3.0", col_types = "text", na = c("", "NA"))

# MitoCarta gene inventory
mito_genes <- mitocarta |>
  pull(Symbol) |>
  na.omit() |>
  unique()
stopifnot(length(mito_genes) == 1136)

# Gene-level functional annotations
mito_annotations <- mitocarta |>
  transmute(
    symbol = Symbol,
    mito_localization = MitoCarta3.0_SubMitoLocalization,
    mito_function = MitoCarta3.0_MitoPathways
  ) |>
  filter(!is.na(symbol), symbol != "") |>
  distinct()
stopifnot(!anyDuplicated(mito_annotations$symbol))

#### MitoPathways3.0 hierarchy ####
mito_hierarchy <- read_excel(mitocarta_file, sheet = "C MitoPathways", col_types = "text", na = c("", "NA")) |>
  select(
    pathway = MitoPathway,
    hierarchy = `MitoPathways Hierarchy`,
    genes = Genes
  ) |>
  filter(!is.na(pathway), !is.na(hierarchy)) |>
  separate_wider_delim(
    hierarchy,
    delim = " > ",
    names = c("level1", "level2", "level3"),
    too_few = "align_start",
    cols_remove = FALSE
  ) |>
  mutate(
    pathway_gmx = pathway |>
      str_remove_all(",") |>
      str_squish() |>
      str_replace_all(" ", "_")
  )

# All available mitochondrial pathways
all_mito_pathways <- all_pathways |> filter(collection == "MitoCarta")

# Check correspondence between official hierarchy and GMX
stopifnot(all(all_mito_pathways$pathway %in% mito_hierarchy$pathway_gmx))
stopifnot(!anyDuplicated(mito_hierarchy$pathway_gmx))

#### Leading-edge (LE) genes ####
# Original gene-level DESeq2 results
gene_de <- de_grp |>
  select(
    ensembl_gene_id, symbol, description, baseMean, log2FC_mle, 
    log2FC_shr, pvalue, padj_gene = padj
  ) |>
  left_join(
    de_grp_fmi |>
      select(ensembl_gene_id, log2FC_mle_fmi, pvalue_fmi, padj_gene_fmi = padj_fmi),
    by = "ensembl_gene_id"
  ) |>
  mutate(
    ensgene = str_remove(ensembl_gene_id, "\\.\\d+$"),
    perc_diff = (2^log2FC_mle - 1) * 100,
    perc_diff_fmi = (2^log2FC_mle_fmi - 1) * 100,
    mitocarta_gene = symbol %in% mito_genes
  ) |>
  left_join(mito_annotations, by = "symbol")

# Match gene-level results to original ranked statistics
make_LE_annotation <- function(wald) {

  wald |>
    mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
    filter(is.finite(stat)) |>
    arrange(desc(abs(stat))) |>
    distinct(ensgene, .keep_all = TRUE) |>
    select(ensembl_gene_id, stat) |>
    left_join(gene_de, by = "ensembl_gene_id")
}

gene_de_grp <- make_LE_annotation(wald_grp)
gene_de_fmi <- make_LE_annotation(wald_grp_fmi)

# Extract LE genes from significant pathways
extract_LE <- function(gsea, annotation) {

  gsea |>
    filter(padj < 0.05) |>
    select(collection, pathway, size, NES, padj_pathway = padj, leadingEdge) |>
    unnest_longer(leadingEdge, values_to = "ensgene") |>
    left_join(annotation, by = "ensgene") |>
    arrange(collection, pathway, desc(abs(stat)))
}

LE_grp <- extract_LE(gsea_grp, gene_de_grp)
LE_grp_fmi <- extract_LE(gsea_grp_fmi, gene_de_fmi)

#### LE overlap: main vs. FMI-adjusted model ####
# Jaccard similarity
jaccard <- function(a, b) {

  if (length(union(a, b)) == 0) return(NA_real_)
  
  length(intersect(a, b)) / length(union(a, b))
}

LE_overlap <- gsea_robust |>
  mutate(
    n_LE = lengths(leadingEdge),
    n_LE_fmi = lengths(leadingEdge_fmi),
    shared_genes = map2(leadingEdge, leadingEdge_fmi, intersect),
    n_shared = lengths(shared_genes),
    jaccard = map2_dbl(leadingEdge, leadingEdge_fmi, jaccard)
  ) |>
  select(
    collection, pathway,
    NES, NES_fmi, padj, padj_fmi,
    n_LE, n_LE_fmi, n_shared,
    jaccard, shared_genes
  )

# Shared LE genes between models
LE_shared <- LE_overlap |>
  unnest_longer(shared_genes, values_to = "ensgene") |>
  left_join(gene_de_grp, by = "ensgene") |>
  arrange(collection, pathway, desc(abs(stat)))

#### OXPHOS and lipid metabolism/FAO focus (following GSEA) ####
# based on the official MitoPathways3.0 hierarchy following
mito_candidates <- mito_hierarchy |>
  filter(level1 == "OXPHOS" | level2 == "Lipid metabolism") |>
  transmute(
    collection = "MitoCarta",
    pathway = pathway_gmx,
    level1, level2, level3,
    process = case_when(
      level1 == "OXPHOS" ~ "OXPHOS",
      level2 == "Lipid metabolism" ~ "Lipid metabolism"
    )
  ) |>
  distinct()

# Number of ranked genes represented in each pathway
mito_pathway_sizes <- tibble(
  pathway = names(mito_sets),
  n_genes_main = map_int(
    mito_sets,
    ~ length(intersect(.x, names(ranks_grp)))
  ),
  n_genes_fmi = map_int(
    mito_sets,
    ~ length(intersect(.x, names(ranks_grp_fmi)))
  )
)

# Eligibility and significance of candidate pathways
mito_candidates_check <- mito_candidates |>
  left_join(
    mito_compare |>
      select(collection, pathway, size, size_fmi, sig_main, sig_fmi, same_direction),
    by = c("collection", "pathway")
  ) |>
  left_join(mito_pathway_sizes, by = "pathway") |>
  mutate(
    tested_main = !is.na(size),
    tested_fmi = !is.na(size_fmi)
  )

# Select pathways significant in both models with concordant direction
mito_focus <- mito_candidates |>
  semi_join(gsea_robust, by = c("collection", "pathway"))

#### GSEA and LE results: mitochondrial focus ####
gsea_mito_focus <- gsea_robust |>
  inner_join(mito_focus, by = c("collection", "pathway"))

# Main model LE genes
LE_mito_focus <- LE_grp |>
  inner_join(mito_focus, by = c("collection", "pathway"))

# FMI-adjusted LE genes
LE_mito_focus_fmi <- LE_grp_fmi |>
  inner_join(mito_focus, by = c("collection", "pathway"))

# Shared LE genes between models
LE_mito_focus_shared <- LE_shared |>
  inner_join(mito_focus, by = c("collection", "pathway"))

# LE overlap between models
LE_mito_focus_overlap <- LE_overlap |>
  inner_join(mito_focus, by = c("collection", "pathway"))

#### CHECKS ####

# MSigDB version
unique(hallmark$db_version)
unique(reactome$db_version)

# Number of ranked genes
length(ranks_grp)
length(ranks_grp_fmi)
length(intersect(names(ranks_grp), names(ranks_grp_fmi)))

# Available gene sets
all_pathways |>
  count(collection)

# Number of significant pathways per collection and model
gsea_compare |>
  count(collection, significance) |>
  pivot_wider(
    names_from = significance,
    values_from = n,
    values_fill = 0
  ) |>
  print(width = Inf)

# NES concordance between models
gsea_compare |>
  group_by(collection) |>
  summarise(
    rho = cor(
      NES, NES_fmi,
      method = "spearman",
      use = "complete.obs"
    ),
    .groups = "drop"
  )

#### GSEA results checks ####

# Inspect pathway names, direction and significance
show_gsea <- function(df) {

  df |>
    filter(sig_main | sig_fmi) |>
    transmute(
      pathway,
      direction_main = case_when(
        NES > 0 ~ "+",
        NES < 0 ~ "-",
        TRUE ~ NA_character_
      ),
      significant_main = sig_main,
      direction_fmi = case_when(
        NES_fmi > 0 ~ "+",
        NES_fmi < 0 ~ "-",
        TRUE ~ NA_character_
      ),
      significant_fmi = sig_fmi
    ) |>
    print(n = Inf, width = Inf)
}

show_gsea(hallmark_compare)
show_gsea(reactome_compare)
show_gsea(mito_compare)

#### MitoPathways hierarchy checks ####

# Unmatched GMX pathway names
all_mito_pathways |>
  anti_join(
    mito_hierarchy,
    by = c("pathway" = "pathway_gmx")
  )

# Duplicate pathway names
mito_hierarchy |>
  count(pathway_gmx) |>
  filter(n > 1)

# All top-level mitochondrial processes
mito_hierarchy |>
  distinct(level1)

# OXPHOS and lipid metabolism candidates
mito_candidates |>
  count(process)

# Candidate pathway eligibility and significance
mito_candidates_check |>
  select(
    pathway, process,
    tested_main, sig_main,
    tested_fmi, sig_fmi
  ) |>
  print(n = Inf, width = Inf)

# Pathways absent from one or both GSEA result tables
mito_candidates_check |>
  filter(!tested_main | !tested_fmi) |>
  select(
    pathway,
    n_genes_main, n_genes_fmi,
    tested_main, tested_fmi
  ) |>
  print(n = Inf, width = Inf)

# Selected mitochondrial pathways
mito_focus |>
  print(n = Inf, width = Inf)

# Number of selected pathways
mito_focus |>
  count(process)

#### LE results checks ####

# Annotation completeness
sum(is.na(LE_grp$symbol))
sum(is.na(LE_grp_fmi$symbol))

sum(is.na(LE_mito_focus$symbol))
sum(is.na(LE_mito_focus_fmi$symbol))

# MitoCarta functional annotation
LE_mito_focus |>
  distinct(ensgene, symbol, mitocarta_gene, mito_function) |>
  summarise(
    n_genes = n(),
    n_mitocarta = sum(mitocarta_gene, na.rm = TRUE),
    n_annotated = sum(!is.na(mito_function))
  )

# Jaccard: same pathway across models
LE_mito_focus_overlap |>
  select(
    pathway, process,
    n_LE, n_LE_fmi, n_shared, jaccard
  ) |>
  arrange(process, desc(jaccard)) |>
  print(n = Inf, width = Inf)

# Jaccard summary per biological process
LE_mito_focus_overlap |>
  group_by(process) |>
  summarise(
    n_pathways = n(),
    median_jaccard = median(jaccard),
    min_jaccard = min(jaccard),
    max_jaccard = max(jaccard),
    .groups = "drop"
  )

warnings()