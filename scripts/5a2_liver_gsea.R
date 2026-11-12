
# Liver RNA-seq GSEA - %FFMI loss groups at 1y
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(fgsea)
library(msigdbr)
library(MetBrewer)
source("scripts/assets/functions.R")

# Output
dir.create("results/tables/RNAseq/gsea", recursive = TRUE, showWarnings = FALSE)
dir.create("results/graphs/RNAseq/gsea", recursive = TRUE, showWarnings = FALSE)

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)
gsea_cols <- c("Up" = renoir_15[11], "Down" = renoir_15[3])

#### Data: DESeq2 Wald statistics ####
# Original Wald statistics from 5a1_liver_deseq2.R
de_grp <- read_csv("results/tables/liver_de_grp_wald_stats.csv")
de_grp_fmi <- read_csv("results/tables/liver_de_grp_fmi_wald_stats.csv")

stopifnot(setequal(de_grp$ensembl_gene_id, de_grp_fmi$ensembl_gene_id))

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
ranks_fmi <- make_ranks(de_grp_fmi)

# Checks
length(ranks_grp)
length(ranks_fmi)
summary(ranks_grp)
summary(ranks_fmi)
stopifnot(!anyDuplicated(names(ranks_grp)))
stopifnot(!anyNA(ranks_grp))

#### MSigDB gene sets ####
hallmark <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = "H")
reactome <- msigdbr(db_species = "HS", species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")

# Database version
unique(hallmark$db_version)
unique(reactome$db_version)

# Convert gene sets into lists of Ensembl IDs
make_gene_sets <- function(df) {

  df <- df |>
    dplyr::filter(!is.na(ensembl_gene), ensembl_gene != "") |>
    distinct(gs_name, ensembl_gene)

  split(df$ensembl_gene, df$gs_name)
}

hallmark_sets <- make_gene_sets(hallmark)
reactome_sets <- make_gene_sets(reactome)

# Gene-set coverage
tibble(
  collection = c("Hallmark", "Reactome"),
  available_sets = c(length(hallmark_sets), length(reactome_sets)),
  median_genes = c(
    median(lengths(hallmark_sets)),
    median(lengths(reactome_sets))
  )
)

#### GSEA: Hallmark and Reactome ####
# Preranked GSEA; BH correction within each collection
run_gsea <- function(gene_sets, ranks, collection, model) {

  set.seed(13)

  fgsea::fgseaMultilevel(
    pathways = gene_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    BPPARAM = BiocParallel::SerialParam()
  ) |>
    as_tibble() |>
    mutate(
      collection = collection,
      model = model,
      direction = if_else(NES > 0, "Up", "Down"),
      .before = 1
    ) |>
    arrange(pval)
}

# Main model
gsea_hallmark <- run_gsea(hallmark_sets, ranks_grp, "Hallmark", "Main")
gsea_reactome <- run_gsea(reactome_sets, ranks_grp, "Reactome", "Main")

# FMI-adjusted sensitivity model
gsea_hallmark_fmi <- run_gsea(hallmark_sets, ranks_fmi, "Hallmark", "FMI-adjusted")
gsea_reactome_fmi <- run_gsea(reactome_sets, ranks_fmi, "Reactome", "FMI-adjusted")

#### Results ####
# Summary
gsea_summary <- bind_rows(
  gsea_hallmark,
  gsea_reactome,
  gsea_hallmark_fmi,
  gsea_reactome_fmi
) |>
  group_by(model, collection) |>
  summarise(
    pathways_tested = n(),
    significant_fdr = sum(padj < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# Significant pathways
gsea_hallmark |> dplyr::filter(padj < 0.05)
gsea_reactome |> dplyr::filter(padj < 0.05)

# Top-ranked pathways
gsea_hallmark |>
  dplyr::select(pathway, size, ES, NES, pval, padj) |>
  slice_head(n = 15)

gsea_reactome |>
  dplyr::select(pathway, size, ES, NES, pval, padj) |>
  slice_head(n = 15)

#### Mitochondrial pathways ####
# Explore mito pathways already tested in Hallmark and Reactome
mito_pattern <- "OXIDATIVE_PHOSPHORYLATION|RESPIRAT|ELECTRON_TRANSPORT|TCA_CYCLE|MITOCHONDR|ATP_SYNTH"

gsea_mito <- bind_rows(gsea_hallmark, gsea_reactome) |>
  dplyr::filter(str_detect(pathway, mito_pattern)) |>
  arrange(pval)

gsea_mito |>
  dplyr::select(collection, pathway, size, NES, pval, padj)

#### Save results ####
# Convert leading-edge gene lists into character strings for CSV export
save_gsea <- function(df, filename) {

  df |>
    mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";"))) |>
    write_csv(filename)
}

save_gsea(gsea_hallmark, "results/tables/RNAseq/gsea/liver_gsea_hallmark.csv")
save_gsea(gsea_reactome, "results/tables/RNAseq/gsea/liver_gsea_reactome.csv")
save_gsea(gsea_hallmark_fmi, "results/tables/RNAseq/gsea/liver_gsea_hallmark_fmi.csv")
save_gsea(gsea_reactome_fmi, "results/tables/RNAseq/gsea/liver_gsea_reactome_fmi.csv")
save_gsea(gsea_mito, "results/tables/RNAseq/gsea/liver_gsea_mito.csv")
write_csv(gsea_summary, "results/tables/RNAseq/gsea/liver_gsea_summary.csv")

#### Plots ####
# FDR-significant Hallmark pathways
hallmark_focus <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_FATTY_ACID_METABOLISM")

hallmark_plot_df <- gsea_hallmark |>
  filter(padj < 0.05) |>
  mutate(
    pathway_label = str_remove(pathway, "^HALLMARK_") |>
      str_replace_all("_", " "),
    pathway_label = fct_reorder(pathway_label, NES)
  )

p_hallmark <- ggplot(hallmark_plot_df, aes(x = NES, y = pathway_label, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(
    values = gsea_cols,
    labels = c("Up" = "Higher in high loss", "Down" = "Higher in modest/low loss"),
    name = NULL
  ) +
  scale_y_discrete(
    labels = parse(text = ifelse(
      levels(hallmark_plot_df$pathway_label) %in% c("OXIDATIVE PHOSPHORYLATION", "FATTY ACID METABOLISM"),
      paste0("bold('", levels(hallmark_plot_df$pathway_label), "')"),
      paste0("'", levels(hallmark_plot_df$pathway_label), "'")
    ))
  ) +
  labs(
    title = "Hepatic Hallmark pathway enrichment",
    x = "Normalized enrichment score (NES)",
    y = NULL
  ) +
  theme_minimal_custom() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10, colour = "black"),
    legend.position = "bottom"
  )
ggsave("results/graphs/RNAseq/gsea/liver_gsea_hallmark_significant.pdf", p_hallmark, width = 9, height = 6)

#### GSEA enrichment plots (OXPHOS & FAO) ####
plot_gsea_enrichment <- function(pathway, title, gene_sets, ranks, results) {

  res <- results |>
    filter(.data$pathway == !!pathway) |>
    slice_head(n = 1)

  stopifnot(nrow(res) == 1)

  pd <- fgsea::plotEnrichmentData(gene_sets[[pathway]], ranks)

  line_col <- if (res$NES < 0) renoir_15[3] else renoir_15[11]

  # Separate strip for pathway gene positions
  y_range <- max(pd$spreadES, 0.1)
  tick_bottom <- min(0, pd$negES) - 0.14 * y_range
  tick_top <- tick_bottom + 0.12 * y_range
  strip_line <- tick_top + 0.02 * y_range

  # Statistics inside the plot
  stat_label <- sprintf(
    "NES = %.2f\nFDR = %.3g\nLeading edge = %d/%d genes",
    res$NES, res$padj, length(res$leadingEdge[[1]]), res$size
  )

  ggplot(pd$curve, aes(x = rank, y = ES)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_line(colour = line_col, linewidth = 1.1) +

    # Pathway gene ticks
    geom_segment(
      data = pd$ticks,
      aes(x = rank, xend = rank),
      y = tick_bottom, yend = tick_top,
      linewidth = 0.45, colour = "grey30",
      inherit.aes = FALSE
    ) +

    # Horizontal line separating gene positions from curve
    geom_segment(
      aes(
        x = 1, xend = length(ranks),
        y = strip_line, yend = strip_line
      ),
      inherit.aes = FALSE,
      linewidth = 0.45, colour = "grey40"
    ) +
    
    # Bottom border of gene-position strip
    geom_hline(yintercept = tick_bottom, linewidth = 0.45, colour = "grey40") +

    # NES and FDR above the enrichment curve
    annotate(
      "label",
      x = 0.045 * length(ranks),
      y = max(0, pd$posES) + 0.38 * y_range,
      label = stat_label,
      hjust = 0, vjust = 1,
      size = 3.4,
      label.size = 0,
      fill = "white",
      colour = "grey20"
    ) +

    scale_x_continuous(
      breaks = c(1, length(ranks)),
      labels = c("Higher in high loss", "Higher in modest/low loss")
    ) +
    coord_cartesian(ylim = c(tick_bottom - 0.015 * y_range, max(0, pd$posES) + 0.42 * y_range)) +
    labs(
      title = title,
      x = "Genes ranked by DESeq2 Wald statistic",
      y = "Running enrichment score"
    ) +
    theme_minimal_custom() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

## i. OXPHOS ##
p_oxphos <- plot_gsea_enrichment("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Hepatic OXPHOS", hallmark_sets, ranks_grp, gsea_hallmark)
ggsave("results/graphs/RNAseq/gsea/liver_gsea_oxphos.pdf", p_oxphos, width = 8.5, height = 5)

## ii. Fatty-acid metabolism ##
p_fa <- plot_gsea_enrichment("HALLMARK_FATTY_ACID_METABOLISM", "Hepatic fatty-acid metabolism", hallmark_sets, ranks_grp, gsea_hallmark)
ggsave("results/graphs/RNAseq/gsea/liver_gsea_fa_metabolism.pdf", p_fa, width = 8.5, height = 5)