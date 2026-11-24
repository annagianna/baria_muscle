
# Liver RNA-seq GSEA and ML-selected metabolite correlations
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(fgsea)
library(msigdbr)
library(phyloseq)
library(MetBrewer)
source("scripts/assets/functions.R")

# Output
table_dir <- "results/tables/RNAseq/gsea"
graph_dir <- "results/graphs/RNAseq/gsea"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)
gsea_cols <- c(
  "Higher in high loss" = renoir_15[11],
  "Higher in modest/low loss" = renoir_15[3]
)

# Data
de_grp <- read_csv("results/tables/liver_de_grp_wald_stats.csv")
liver_rnaseq_vst <- readRDS("data/processed_data/BARIA_Liver_RNAseq_vst_batchcorrected.RDS")
metab <- readRDS("data/processed_data/BARIA_metabolon_clean.RDS")
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
mito_file <- "data/raw_data/MitoCarta/Human.MitoPathways3.0.gmx"

#### Ranked gene statistics ####
# Positive = higher expression in high FFMI loss
# Negative = higher expression in modest/low FFMI loss

make_ranks <- function(df) {
  df <- df |>
    mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
    filter(!is.na(stat), is.finite(stat)) |>
    arrange(desc(abs(stat))) |>
    distinct(ensgene, .keep_all = TRUE)

  ranks <- df$stat
  names(ranks) <- df$ensgene
  sort(ranks, decreasing = TRUE)
}

ranks_grp <- make_ranks(de_grp)
stopifnot(!anyDuplicated(names(ranks_grp)), !anyNA(ranks_grp))

#### Hallmark and Reactome gene sets ####

hallmark <- msigdbr(
  db_species = "HS", species = "Homo sapiens",
  collection = "H"
)

reactome <- msigdbr(
  db_species = "HS", species = "Homo sapiens",
  collection = "C2", subcollection = "CP:REACTOME"
)

make_gene_sets <- function(df) {
  df <- df |>
    filter(!is.na(ensembl_gene), ensembl_gene != "") |>
    distinct(gs_name, ensembl_gene)

  split(df$ensembl_gene, df$gs_name)
}

hallmark_sets <- make_gene_sets(hallmark)
reactome_sets <- make_gene_sets(reactome)

#### MitoCarta 3.0 gene sets ####

dir.create(dirname(mito_file), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(mito_file)) {
  download.file(
    "https://personal.broadinstitute.org/scalvo/MitoCarta_Download/Human.MitoPathways3.0.gmx",
    mito_file, mode = "wb"
  )
}

mito_gmx <- read.delim(
  mito_file, header = FALSE, fill = TRUE,
  check.names = FALSE, quote = ""
)

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
  unique(ids[!is.na(ids) & ids %in% names(ranks_grp)])
})

gene_sets <- list(
  Hallmark = hallmark_sets,
  Reactome = reactome_sets,
  MitoCarta = mito_sets
)

#### Preranked GSEA ####
# Main model only; BH-FDR within each collection

run_gsea <- function(sets, collection) {
  set.seed(13)

  fgsea::fgseaMultilevel(
    pathways = sets,
    stats = ranks_grp,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    BPPARAM = BiocParallel::SerialParam()
  ) |>
    as_tibble() |>
    mutate(
      collection = collection,
      direction = if_else(NES > 0, "Up", "Down"),
      .before = 1
    ) |>
    arrange(pval)
}

gsea_results <- imap_dfr(gene_sets, ~ run_gsea(.x, .y))

# Save complete GSEA results
gsea_results |>
  mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";"))) |>
  write_csv(file.path(table_dir, "liver_gsea_main_all_collections.csv"))

# Significant pathways
selected_pathways <- gsea_results |>
  filter(padj < 0.05) |>
  select(collection, pathway, size, NES, padj) |>
  arrange(collection, padj)

stopifnot(nrow(selected_pathways) > 0)

write_csv(
  selected_pathways,
  file.path(table_dir, "liver_gsea_main_significant_pathways.csv")
)

#### Plot: GSEA enrichment by collection ####
# All significant Hallmark and MitoCarta pathways
# Top 10 per enrichment direction for Reactome

gsea_plot_df <- selected_pathways |>
  mutate(
    direction = if_else(
      NES > 0, "Higher in high loss", "Higher in modest/low loss"
    ),
    pathway_label = pathway |>
      str_remove("^HALLMARK_|^REACTOME_") |>
      str_replace_all("_", " ") |>
      str_to_sentence()
  ) |>
  group_by(collection, direction) |>
  arrange(padj, .by_group = TRUE) |>
  filter(collection != "Reactome" | row_number() <= 10) |>
  ungroup()

for (coll in unique(gsea_plot_df$collection)) {

  plot_df <- gsea_plot_df |>
    filter(collection == coll) |>
    mutate(pathway_label = fct_reorder(pathway_label, NES))

  p <- ggplot(plot_df, aes(x = NES, y = pathway_label, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    scale_fill_manual(values = gsea_cols, name = NULL) +
    labs(
      title = paste("Hepatic", coll, "pathway enrichment"),
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    theme_minimal_custom() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9, colour = "black"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )

  ggsave(
    file.path(graph_dir, paste0("liver_", tolower(coll), "_gsea_enrichment.pdf")),
    p,
    width = 11,
    height = max(5, 0.32 * nrow(plot_df) + 2)
  )
}

#### Participant-level hepatic pathway scores ####
# Mean z-standardized expression across genes within each pathway

expr_ids <- names(liver_rnaseq_vst) |>
  str_remove("\\.\\d+$")

pathway_scores <- pmap_dfr(
  selected_pathways |> select(collection, pathway),
  function(collection, pathway) {

    genes <- names(liver_rnaseq_vst)[
      expr_ids %in% gene_sets[[collection]][[pathway]]
    ]

    expr <- as.matrix(liver_rnaseq_vst[, genes, drop = FALSE])
    expr <- expr[, apply(expr, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

    if (ncol(expr) == 0) stop("No expressed genes for pathway: ", pathway)

    tibble(
      id = as.numeric(as.character(liver_rnaseq_vst$id)),
      collection = collection,
      pathway = pathway,
      n_genes = ncol(expr),
      score = rowMeans(scale(expr), na.rm = TRUE)
    )
  }
)

#### Baseline top-15 ML metabolites ####

top15_metabolites <- get_top_metabolites()$metabolite

metab_mat <- as(otu_table(metab), "matrix")
if (taxa_are_rows(metab)) metab_mat <- t(metab_mat)

stopifnot(all(top15_metabolites %in% colnames(metab_mat)))

metab_df <- as.data.frame(metab_mat[, top15_metabolites, drop = FALSE]) |>
  rownames_to_column("sample") |>
  mutate(id = sample |> str_remove("^BARIA_") |> str_remove("_v0$") |> as.numeric()) |>
  select(-sample)

stopifnot(!anyDuplicated(metab_df$id))

#### Match liver RNA-seq and metabolomics ####

pathway_metab <- pathway_scores |>
  inner_join(metab_df, by = "id") |>
  semi_join(
    baria_muscle_wide |>
      filter(!is.na(perc_change_ffmi_v4_group)) |>
      mutate(id = as.numeric(as.character(id))) |>
      distinct(id),
    by = "id"
  )

#### Pathway x ML metabolite correlations ####
# BH-FDR across all Hallmark, Reactome and MitoCarta pairs

cor_pathway_metab <- pathway_metab |>
  group_by(collection, pathway) |>
  group_modify(~ map_dfr(top15_metabolites, function(m) {

    valid <- complete.cases(.x$score, .x[[m]])

    if (sum(valid) < 3 || n_distinct(.x$score[valid]) < 2 ||
        n_distinct(.x[[m]][valid]) < 2) {
      return(tibble(
        metabolite = m, n = sum(valid),
        rho = NA_real_, p.value = NA_real_
      ))
    }

    test <- cor.test(
      .x$score[valid], .x[[m]][valid],
      method = "spearman", exact = FALSE
    )

    tibble(
      metabolite = m,
      n = sum(valid),
      rho = unname(test$estimate),
      p.value = test$p.value
    )
  })) |>
  ungroup() |>
  mutate(p_fdr = p.adjust(p.value, method = "BH")) |>
  arrange(p_fdr)

sig_pathway_metab <- cor_pathway_metab |> filter(p_fdr < 0.05)

write_csv(
  cor_pathway_metab,
  file.path(table_dir, "liver_gsea_ml_metabolite_correlations.csv")
)

#### Plot: significant pathway x ML metabolite correlations ####
# Cleveland dot plots; colours identify metabolites

metab_plot_colors <- setNames(
  renoir_15[seq_along(top15_metabolites)],
  metabolite_label(top15_metabolites)
)

for (coll in c("Hallmark", "MitoCarta", "Reactome")) {

  plot_df <- sig_pathway_metab |>
    filter(collection == coll) |>
    mutate(
      pathway_label = pathway |>
        str_remove("^HALLMARK_|^REACTOME_") |>
        str_replace_all("_", " ") |>
        str_to_sentence(),
      metabolite_name = metabolite_label(metabolite)
    )

  if (nrow(plot_df) == 0) next

  if (coll == "Reactome") {

    # Divide pathways into two similarly sized panels
    pathway_panels <- plot_df |>
      group_by(pathway_label) |>
      summarise(mean_rho = mean(rho), .groups = "drop") |>
      arrange(mean_rho) |>
      mutate(
        panel = if_else(row_number() <= ceiling(n() / 2),
                        "Panel 1", "Panel 2"),
        plot_label = paste0(pathway_label, "___", panel)
      )

    plot_df <- plot_df |>
      left_join(
        pathway_panels |> select(pathway_label, panel, plot_label),
        by = "pathway_label"
      ) |>
      mutate(plot_label = factor(
        plot_label, levels = pathway_panels$plot_label
      ))

  } else {
    plot_df <- plot_df |>
      mutate(plot_label = fct_reorder(pathway_label, rho))
  }

  p <- ggplot(plot_df, aes(x = rho, y = plot_label)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_point(aes(colour = metabolite_name), size = 3) +
    scale_colour_manual(values = metab_plot_colors, name = "Metabolite") +
    labs(
      title = paste("Hepatic", coll, "pathways x ML metabolites"),
      x = "Spearman correlation coefficient (rho)",
      y = NULL
    ) +
    theme_minimal_custom() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9, colour = "black"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )

  if (coll == "Reactome") {
    p <- p +
      facet_wrap(~panel, ncol = 2, scales = "free_y") +
      scale_y_discrete(
        labels = function(x) str_wrap(
          str_remove(x, "___Panel [12]$"), width = 30
        )
      ) +
      theme(strip.text = element_blank())
  }

  ggsave(
    file.path(graph_dir, paste0("liver_", tolower(coll), "_ml_metabolite_dotplot.pdf")),
    p,
    width = if (coll == "Reactome") 15 else if (coll == "MitoCarta") 7 else 6.5,
    height = if (coll == "Reactome") 7.5 else max(4, 0.4 * n_distinct(plot_df$pathway) + 2)
  )
}

#### Results checks ####

# GSEA pathways
gsea_results |>
  group_by(collection) |>
  summarise(
    n_tested = n(),
    n_significant = sum(padj < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) |>
  print()

# Participant overlap
pathway_metab |>
  summarise(
    n_participants = n_distinct(id),
    n_pathways = n_distinct(pathway),
    n_metabolites = length(top15_metabolites)
  ) |>
  print()

# Significant pathway-metabolite correlations
cor_pathway_metab |>
  group_by(collection) |>
  summarise(
    n_tests = n(),
    n_significant = sum(p_fdr < 0.05, na.rm = TRUE),
    n_pathways = n_distinct(pathway[p_fdr < 0.05 & !is.na(p_fdr)]),
    n_metabolites = n_distinct(metabolite[p_fdr < 0.05 & !is.na(p_fdr)]),
    .groups = "drop"
  ) |>
  print()