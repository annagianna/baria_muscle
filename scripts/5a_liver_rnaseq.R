# Liver RNA-seq
# Anna Giannakogeorgou

# Packages
library(tidyverse)
library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(MetBrewer)
library(annotables)
library(ggpubr)
library(broom)
library(DESeq2)
library(apeglm)
library(ppcor)
source("scripts/assets/functions.R")
dir.create("results/graphs/RNAseq", recursive = TRUE, showWarnings = FALSE)

# Bioconductor packages (DESeq2/S4Vectors/matrixStats) mask some dplyr verbs; always use dplyr's
select <- dplyr::select
rename <- dplyr::rename
count  <- dplyr::count
filter <- dplyr::filter

# Theme
renoir_15 <- met.brewer("Renoir", n = 15)
renoir_cols_20 <- met.brewer("Renoir", n = 20)
de_cols <- c( # same high/low colors used in the alpha, beta diveristy scripts
  "higher in high FFMI loss" = renoir_cols_20[18],
  "lower in high FFMI loss"  = renoir_cols_20[5],
  "not significant" = "grey80"
)

# Significance asterisks based on p_fdr
stars <- function(p_fdr) ifelse(p_fdr < 0.001, "***", ifelse(p_fdr < 0.01, "**", ifelse(p_fdr < 0.05, "*", "")))

# Data
liver_rnaseq <- readRDS("data/processed_data/BARIA_Liver_RNAseq.RDS") # raw counts for DESeq2
liver_rnaseq_vst <- readRDS("data/processed_data/BARIA_Liver_RNAseq_vst.RDS") # VST normalized
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_mb_v0 <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
forest_perc_change_ffmi_v4 <- read.csv("results/mlmodels/perc_change_ffmi_v4/all/forest_results_top15.csv")
metab <- readRDS("data/processed_data/BARIA_metabolon_clean.RDS")

# Prep data
## Top15 species
# Top 15 species by ML feature importance for %FFMI change at 1y
fi <- get_feature_importance("results/mlmodels/perc_change_ffmi_v4/all", "perc_change_ffmi_v4_all", "reg")
top15_species <- top_features(fi, n = 15)$FeatName
top15_species_labels <- tibble(species = top15_species, species_label = species_label(top15_species)) |> 
  left_join(
    forest_perc_change_ffmi_v4 |> 
      dplyr::select(species, estimate, p_fdr) |> 
      mutate(estimate_direction = if_else(estimate > 0, "positive", "negative")),
    by = "species"
  )

# Baseline abundance of top-15 species
top15_species_v0 <- t(as(otu_table(baria_mb_v0), "matrix")) |>
  as.data.frame() |>
  rownames_to_column("Sample") |>
  mutate(id = Sample |> 
    str_remove("^BARIA_") |> 
    str_remove("_v0$") |> 
    as.numeric()
) |>
  dplyr::select(id, all_of(top15_species))

# Log10-transform abundance & species-specific pseudocount
top15_species_log10 <- top15_species_v0 |>
  mutate(across(all_of(top15_species), ~ {
    pseudo <- min(.x[.x > 0]) / 2
    log10(.x + pseudo)
  }))

# Map gene ids to their names
gene_ids <- tibble(
  ensembl_gene_id = liver_rnaseq_vst |> 
    dplyr::select(starts_with("ENSG")) |> 
    colnames()
) |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$"))

# Gene annotations
gene_annotations <- grch38 |>
  filter(ensgene %in% gene_ids$ensgene) |>
  dplyr::select(ensgene, symbol, biotype, description) |>
  distinct()

# Display labels for heatmaps: "SYMBOL – short description" / only symbol / only description
gene_label <- function(symbols, type = c("both", "symbol", "description"), max_char = 45) {

  type <- match.arg(type)

  desc <- gene_annotations$description[match(symbols, gene_annotations$symbol)] |>
    str_remove("\\s*\\[Source.*\\]$") |> # grch38 descriptions can carry a "[Source:...]" suffix
    str_trunc(max_char) |>
    coalesce(symbols) # fall back to symbol if no description

  switch(type,
         both = paste0(symbols, " – ", desc),
         symbol = symbols,
         description = desc)
}

#### Liver RNA-seq x top15 species ####
# Top15 species x liver RNA seq
liver_mb <- liver_rnaseq_vst |>
  mutate(id = as.numeric(as.character(id))) |> 
  inner_join(top15_species_log10, by = "id")

# QC: sequencing depth (DESeq2 size factor) vs. species abundance (no expected association)
size_factor_qc <- map_dfr(top15_species,
  ~ broom::tidy(cor.test(liver_mb$size_factor, liver_mb[[.x]], method = "spearman", exact = FALSE)
) |> mutate(species = .x)) |>
  dplyr::select(species, rho = estimate, p.value)
size_factor_qc

# Genes in the (VST-)normalized file (already expression-filtered in cleaning script 0c)
liver_genes <- liver_rnaseq_vst |> 
  dplyr::select(starts_with("ENSG")) |>
  colnames()

### Untargeted ###
# Spearman correlations top15 species
# Align species order
species_mat <- liver_mb |>
  dplyr::select(all_of(top15_species)) |>
  as.matrix()

gene_mat <- liver_mb |>
  dplyr::select(all_of(liver_genes)) |>
  as.matrix()

cor_results <- expand_grid(
  species = top15_species,
  ensembl_gene_id = colnames(gene_mat)
) |>
  mutate(
    test = map2(species, ensembl_gene_id, ~ cor.test(
      species_mat[, .x],
      gene_mat[, .y],
      method = "spearman",
      exact = FALSE
    )),
    rho = map_dbl(test, "estimate"),
    p.value = map_dbl(test, "p.value")
  ) |>
  dplyr::select(-test) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH"))

# Top 25 (annotated) associated genes by |rho| (among FDR-signif associations)
genes_heatmap <- cor_results |>
  group_by(ensembl_gene_id) |>
  summarize(
    min_fdr = min(p_fdr),
    n_sig = sum(p_fdr < 0.05),
    max_abs_rho = max(abs(rho)),
    .groups = "drop"
  ) |>
  filter(min_fdr < 0.05) |>
  mutate(ensgene = str_remove(ensembl_gene_id, "\\.\\d+$")) |>
  left_join(
    gene_annotations |> distinct(ensgene, .keep_all = TRUE) |> dplyr::select(ensgene, symbol, description),
    by = "ensgene"
  ) |>
  filter(!is.na(symbol), symbol != "") |> # annotated genes only, selected after FDR correction
  arrange(desc(max_abs_rho)) |>
  slice_head(n = 25) |>
  mutate(gene_label = symbol)

# Heatmap data
heatmap_data <- cor_results |>
  filter(ensembl_gene_id %in% genes_heatmap$ensembl_gene_id)

rho_mat <- heatmap_data |>
  dplyr::select(species, ensembl_gene_id, rho) |>
  pivot_wider(names_from = ensembl_gene_id, values_from = rho) |>
  column_to_rownames("species") |>
  as.matrix()

fdr_mat <- heatmap_data |> 
  dplyr::select(species, ensembl_gene_id, p_fdr) |> 
  pivot_wider(names_from = ensembl_gene_id, values_from = p_fdr) |>
  column_to_rownames("species") |> 
  as.matrix()

# Replace gene labels
colnames(rho_mat) <- genes_heatmap$gene_label[match(colnames(rho_mat), genes_heatmap$ensembl_gene_id)]
colnames(fdr_mat) <- colnames(rho_mat)

# Add species labels
rownames(rho_mat) <- top15_species_labels$species_label[match(rownames(rho_mat), top15_species_labels$species)]
rownames(fdr_mat) <- rownames(rho_mat)

# Diverging colour scale: negative = Renoir[3], positive = Renoir[11]
rho_max <- max(abs(rho_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-rho_max, 0, rho_max), c(renoir_15[3], "white", renoir_15[11]))


# Species label colours by direction of association with 1-year %FFMI change
species_direction_colors <- c(positive = renoir_15[15],negative = renoir_15[9])
species_label_colors <- setNames(species_direction_colors[top15_species_labels$estimate_direction], top15_species_labels$species_label)

## Complex Heatmap ##
rho_mat_t <- t(rho_mat)
fdr_mat_t <- t(fdr_mat)

ht <- Heatmap(
  rho_mat_t,
  name = "Spearman\nrho",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_dend_side = "right",
  row_names_side = "left",
  row_labels = gene_label(rownames(rho_mat_t)),
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = max_text_width(gene_label(rownames(rho_mat_t)), gp = gpar(fontsize = 9)),
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "italic",
    col = species_label_colors[colnames(rho_mat_t)]
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {grid.text(stars(fdr_mat_t[i, j]), x, y)}
)

pdf("results/graphs/RNAseq/liver_top15_species_genes_heatmap.pdf", width = 10, height = 8)
draw(ht)
dev.off()

# More compact heatmap: only species with >=1 signif association among the selected genes
species_keep <- rownames(fdr_mat)[apply(fdr_mat < 0.05, 1, any)]
rho_mat_compact <- rho_mat[species_keep, , drop = FALSE]
fdr_mat_compact <- fdr_mat[species_keep, , drop = FALSE]
rho_mat_compact <- t(rho_mat_compact)
fdr_mat_compact <- t(fdr_mat_compact)

# Compact heatmap: only species with >=1 significant association
ht_compact <- Heatmap(
  rho_mat_compact,
  name = "Spearman\nrho",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_dend_side = "right",
  row_names_side = "left",
  row_labels = gene_label(rownames(rho_mat_compact)),
  row_names_gp = gpar(fontsize = 9),
  row_names_max_width = max_text_width(gene_label(rownames(rho_mat_compact)), gp = gpar(fontsize = 9)),
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "italic",
    col = species_label_colors[colnames(rho_mat_compact)]
  ),
  rect_gp = gpar(col = "white", lwd = 0.7),
  cell_fun = function(j, i, x, y, width, height, fill) {grid.text(stars(fdr_mat_compact[i, j]), x, y)}
)

# Save
pdf("results/graphs/RNAseq/liver_species_genes_heatmap_compact.pdf", width = 9, height = 8)
draw(ht_compact)
dev.off()

#### Targeted genes #####
# Create list of selected genes
core_genes <- list(

  hepatokines = c(
    "FGF21","SELENOP","LECT2","ANGPTL3","ANGPTL4","ANGPTL8","AHSG","FETUB",
    "FGL1","RBP4","IGF1","IGFBP1","IGFBP2","IGFBP3","INHBA","INHBE","FST","FSTL3",
    "GDF15","LEAP2","GPLD1","ENHO","TSKU","SHBG","SMOC1","APOA5"),

  # Fatty acid oxidation & ketogenesis
  fao = c(
    "CPT1A","CPT1B","CPT2","CRAT","CROT","SLC25A20","SLC22A5",
    "ACADVL","ACADL","ACADM","ACADS","ACADSB","ACAD9","ACAD8","ACAD11",
    "HADHA","HADHB","HADH","ECHS1","EHHADH","ACAA1","ACAA2","DECR1",
    "ECI1","ECI2","ETFA","ETFB","ETFDH","MLYCD","ACOX1","ACOX2",
    "HSD17B4","SCP2","PPARA","HMGCS2","BDH1"),

  # Mitochondrial oxitative phosphorylation
  mito_oxphos = c(
    "PPARGC1A","PPARGC1B","PPARD","ESRRA","ESRRG","NRF1","GABPA", "TFAM","TFB2M",
    "POLG","MFN1","MFN2","OPA1","DNM1L","SDHA","SDHB","SDHC","SDHD","CYCS"),

  oxidative_stress = c(
    "NFE2L2","KEAP1","NQO1","HMOX1","GCLC","GCLM","GSR","GSTP1","GSTA1",
    "GSTM1","SLC7A11","TXN","TXNRD1","PRDX1","PRDX3","PRDX5","PRDX6",
    "SOD1","SOD2","SOD3","CAT","GPX1","GPX3","GPX4","SRXN1","G6PD",
    "SESN2","FOXO3"),

  aa_bcaa_turnover = c(
    "BCAT1","BCAT2","BCKDHA","BCKDHB","DBT","DLD","BCKDK","PPM1K", "HIBCH","HIBADH",
    "IVD","MCCC1","MCCC2","AUH","HMGCL","PCCA","PCCB", "MMUT","ALDH6A1",
    "CPS1","OTC","ASS1","ASL","ARG1","NAGS","SLC25A15","SLC25A13",
    "GLUD1","GLS","GLS2","GLUL","GOT1","GOT2","GPT","GPT2"),

  gluconeogenesis = c(
    "PCK1","PCK2","G6PC1","SLC37A4","FBP1","FBP2","PC","MDH1","MDH2",
    "PDK1","PDK2","PDK4","FOXO1","HNF4A","CREB3L3","NR3C1"),

  glycogen_synthesis = c("GYS2","GYG1","GBE1","UGP2","PGM1","PPP1R3B","PPP1R3C","PPP1R3G"),

  glycogenolysis = c("PYGL","AGL","PHKA2","PHKB","PHKG2"),

  glycolysis = c(
    "GCK","GCKR","HK1","HK2","GPI","PFKL","ALDOA","TPI1","GAPDH", "PGK1","PGAM1",
    "ENO1","PKLR","PKM","PDHA1","PDHB","DLAT", "LDHA","LDHB","SLC2A2"),

  fructose_metabolism = c("KHK","ALDOB","TKFC","SORD","AKR1B1","SLC2A5"),

  # De-novo lipogenesis
  dnl = c(
    "SREBF1","MLXIPL","NR1H3","INSIG1","SCAP","THRSP","MID1IP1",
    "ACLY","ACSS2","ACACA","ACACB","FASN","ELOVL6","SCD","ME1"),

  # Triglyceride synthesis/VLDL
  tg_vldl = c("GPAM","AGPAT2","LPIN1","DGAT1","DGAT2","MTTP","APOB","PNPLA3","TM6SF2"),

  ins_signaling = c(
    "INSR","IGF1R","IRS1","IRS2","PIK3R1","PIK3CA","PDPK1","AKT1","AKT2",
    "GSK3B","MTOR","RPTOR","RPS6KB1","EIF4EBP1","TSC1","TSC2",
    "PRKAA1","PRKAA2","STK11","PTEN","PTPN1","PRKCE","SOCS3","TRIB3"),

  ins_clearance = c("CEACAM1","IDE"),

  # Glucagon/incretin receptors + cAMP/PKA/CREB + glucagon-driven amino-acid uptake
  gluc_incr_signaling = c( "GCGR","GIPR","GLP1R","GLP2R","DPP4", "GNAS","PRKACA","CREB1","CRTC2", "SLC7A2","SLC38A4","SLC38A5")
)
stopifnot(!anyDuplicated(unlist(core_genes))) # each gene in exactly one category

target_genes <- enframe(core_genes, name = "category", value = "symbol") |> 
  unnest(symbol) |> 
  left_join(gene_annotations, by = "symbol") |> 
  inner_join(gene_ids, by = "ensgene")

# Spearman correlations: top-15 species x targeted genes
target_gene_mat <- liver_mb |>
  dplyr::select(all_of(target_genes$ensembl_gene_id)) |>
  as.matrix()

cor_targeted <- expand_grid(
  species = top15_species,
  ensembl_gene_id = colnames(target_gene_mat)
) |>
  mutate(
    test = map2(species, ensembl_gene_id, ~ cor.test(
      species_mat[, .x],
      target_gene_mat[, .y],
      method = "spearman",
      exact = FALSE
    )),
    rho = map_dbl(test, "estimate"),
    p.value = map_dbl(test, "p.value")
  ) |>
  dplyr::select(-test) |>
  mutate(p_fdr = p.adjust(p.value, method = "BH")) |> # global FDR p for all selected genes
  left_join(
    target_genes |>
      dplyr::select(ensembl_gene_id, symbol, category, description),
    by = "ensembl_gene_id"
  )

# Genes with >=1 FDR-significant association
target_genes_heatmap <- cor_targeted |>
  filter(p_fdr < 0.05) |>
  distinct(ensembl_gene_id, symbol, description, category)

# Build matrices
rho_targeted <- cor_targeted |>
  filter(ensembl_gene_id %in% target_genes_heatmap$ensembl_gene_id) |>
  dplyr::select(species, symbol, rho) |>
  pivot_wider(names_from = species, values_from = rho) |>
  column_to_rownames("symbol") |>
  as.matrix()

fdr_targeted <- cor_targeted |>
  filter(ensembl_gene_id %in% target_genes_heatmap$ensembl_gene_id) |>
  dplyr::select(species, symbol, p_fdr) |>
  pivot_wider(names_from = species, values_from = p_fdr) |>
  column_to_rownames("symbol") |>
  as.matrix()

# Gene categories
gene_categories <- target_genes_heatmap$category[match(rownames(rho_targeted), target_genes_heatmap$symbol)]

# Category colours: Renoir palette, alternating between its two halves so neighbouring categories contrast
k <- length(core_genes)
half <- ceiling(k / 2)
idx <- as.vector(rbind(1:half, half + 1:half)) # 1, 9, 2, 10, 3, 11, ...
idx <- idx[idx <= k]
category_cols <- setNames(met.brewer("Renoir", n = k)[idx], names(core_genes))
category_anno <- rowAnnotation(Category = gene_categories, col = list(Category = category_cols), show_annotation_name = FALSE)

# Species labels
species_labels <- top15_species_labels$species_label[match(colnames(rho_targeted), top15_species_labels$species)]

# Full targeted heatmap
ht_targeted <- Heatmap(
  rho_targeted,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.3, 0, 0.3), c(renoir_15[3], "white", renoir_15[11])),
  left_annotation = category_anno,
  row_split = gene_categories,
  row_title = NULL,
  row_names_side = "left",
  row_dend_side = "right",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_labels = gene_label(rownames(rho_targeted)),
  row_names_max_width = max_text_width(gene_label(rownames(rho_targeted)), gp = gpar(fontsize = 8)),
  column_labels = species_labels,
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9, fontface = "italic", col = species_label_colors[species_labels]),
  row_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {grid.text(stars(fdr_targeted[i, j]), x, y, gp = gpar(fontsize = 10))}
)

# Save
pdf(
  "results/graphs/RNAseq/liver_top15_species_targeted_genes_heatmap.pdf",
    width = 15, height = 2 + 0.18 * nrow(rho_targeted)
  )
draw(ht_targeted, newpage = FALSE)
dev.off()

#### Individual correlation plots ####
# FDR-significant targeted species x gene associations
sig_targeted <- cor_targeted |>
  filter(p_fdr < 0.05) |>
  arrange(p_fdr)

cor_plots <- pmap(
  sig_targeted,

  function(species, ensembl_gene_id, rho, p.value, p_fdr, symbol, category, description) {

    plot_data <- liver_mb |>
      dplyr::select(
        abundance = all_of(species),
        expression = all_of(ensembl_gene_id)
      )

    species_name <- top15_species_labels$species_label[match(species, top15_species_labels$species)]

    ggplot(plot_data, aes(x = abundance, y = expression)) +
      geom_point(alpha = 0.7, size = 1.8) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7, color = renoir_15[3]) +
      labs(
        title = gene_label(symbol),
        subtitle = paste0(
          "Spearman rho = ", round(rho, 2),
          " | FDR = ", signif(p_fdr, 2)
        ),
        x = paste0(species_name, "\nlog10 relative abundance"),
        y = "Liver gene expression (VST)"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8),
        axis.title.x = element_text(face = "italic"),
        panel.grid.minor = element_blank()
      )
  }
)

# Arrange individual cor plots
n_cols_cor_plots <- 5
n_rows_cor_plots <- ceiling(length(cor_plots) / n_cols_cor_plots) # grid grows with the number of significant pairs
cor_plots_arranged <- ggarrange(plotlist = cor_plots, ncol = n_cols_cor_plots, nrow = n_rows_cor_plots)
ggsave(
  "results/graphs/RNAseq/liver_targeted_significant_correlations.pdf",
   cor_plots_arranged,
   width = 12,
   height = 2.5 * n_rows_cor_plots, # keeps each individual cor plot the same size
   limitsize = FALSE
  )

## Compact targeted heatmap ##
# Display threshold for slides: genes and species with >= min_hits FDR-significant associations
min_hits <- 1
sig_targeted_mat <- fdr_targeted < 0.05
genes_keep_targeted   <- rownames(fdr_targeted)[rowSums(sig_targeted_mat) >= min_hits]
species_keep_targeted <- colnames(fdr_targeted)[colSums(sig_targeted_mat) >= min_hits]

c(genes = length(genes_keep_targeted), species = length(species_keep_targeted))  # check before plotting

rho_targeted_compact <- rho_targeted[genes_keep_targeted, species_keep_targeted, drop = FALSE]
fdr_targeted_compact <- fdr_targeted[genes_keep_targeted, species_keep_targeted, drop = FALSE]

# Category annotation for the kept genes only
gene_categories_compact <- gene_categories[match(genes_keep_targeted, rownames(rho_targeted))]
category_anno_compact <- rowAnnotation(Category = gene_categories_compact, col = list(Category = category_cols), show_annotation_name = FALSE
)

# Species labels
species_labels_compact <- top15_species_labels$species_label[match(colnames(rho_targeted_compact), top15_species_labels$species)]

# Heatmap
ht_targeted_compact <- Heatmap(
  rho_targeted_compact,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.3, 0, 0.3), c(renoir_15[3], "white", renoir_15[11])),
  left_annotation = category_anno_compact,
  row_split = gene_categories_compact,
  row_title = NULL,
  row_names_side = "left",
  row_dend_side = "right",
  row_labels = gene_label(rownames(rho_targeted_compact)),
  row_names_max_width = max_text_width(gene_label(rownames(rho_targeted_compact)), gp = gpar(fontsize = 8)),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_labels = species_labels_compact,
  column_names_rot = 45,
  column_names_gp = gpar(
    fontsize = 10,
    fontface = "italic",
    col = species_label_colors[species_labels_compact]
  ),
  row_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      stars(fdr_targeted_compact[i, j]),
      x, y,
      gp = gpar(fontsize = 10)
    )
  }
)

# Save
# Size grows with the number of genes (rows) and species (columns)
pdf("results/graphs/RNAseq/liver_species_targeted_genes_heatmap_compact.pdf", width = 6 + 0.45 * ncol(rho_targeted_compact), height = 2 + 0.18 * nrow(rho_targeted_compact))
draw(ht_targeted_compact, newpage = FALSE)
dev.off()

# Helper function: Spearman correlations for every x–y pair, BH-correction within the block
cor_block <- function(df, x_vars, y_vars) {
  expand_grid(x = x_vars, y = y_vars) |>
    mutate(
      n = map2_int(x, y, ~ sum(complete.cases(df[[.x]], df[[.y]]))),
      test = map2(x, y, ~ cor.test(df[[.x]], df[[.y]], method = "spearman", exact = FALSE)),
      rho = map_dbl(test, "estimate"),
      p.value = map_dbl(test, "p.value")
    ) |>
    dplyr::select(-test) |>
    mutate(p_fdr = p.adjust(p.value, method = "BH")) |>
    arrange(p.value)
}

liver_df <- liver_mb |>
  dplyr::select(id, all_of(sig_genes$ensembl_gene_id)) |>
  rename_with(~ sig_genes$symbol[match(.x, sig_genes$ensembl_gene_id)], starts_with("ENSG")) |>
  left_join(clin, by = "id")

#### Liver RNA seq x %FFMI change 1y ####


#### Metabolomics ####
# Prep metab data
metab_df <- as(otu_table(metab), "matrix") |> # samples x metabolites
  as_tibble(rownames = "sample") |>
  mutate(id = sample |> str_remove("^BARIA_") |> str_remove("_v0$") |> as.numeric()) |>
  dplyr::select(-sample)
metab_names <- setdiff(names(metab_df), "id")
liver_metab_df <- liver_df |> inner_join(metab_df, by = "id")
nrow(liver_metab_df)  # overlap n

# Correlation liver genes x metabolites
cor_genes_metab <- cor_block(liver_metab_df, sig_genes$symbol, metab_names)
cor_genes_metab |> 
  filter(p_fdr < 0.05) |> 
  dplyr::count(x, sort = TRUE) # hits per gene

cor_genes_metab |> 
  filter(p_fdr < 0.05) |> 
  slice_head(n = 30)

# Correlation metabolites × outcome (%FFMI change at 1y) -> triangulation
metab_hits <- cor_genes_metab |> filter(p_fdr < 0.05) |> distinct(y) |> pull(y)
cor_metab_ffmi <- cor_block(liver_metab_df, metab_hits, "perc_change_ffmi_v4")
cor_metab_ffmi |> slice_head(n = 20)

## Plots ##
# Helper functions
# Long cor_block output -> wide matrix (rows = x, cols = y)
to_matrix <- function(df, value) {
  df |>
    dplyr::select(x, y, all_of(value)) |>
    pivot_wider(names_from = y, values_from = all_of(value)) |>
    column_to_rownames("x") |>
    as.matrix()
}

# Super-pathway strip: same as in 4b_mb_metabolome_correlations.R
pathway_colors <- c(
  "Amino Acid"                         = "#3C9189",
  "Carbohydrate"                       = "#765135",
  "Cofactors and Vitamins"             = "#385C7E",
  "Energy"                             = "#3E6330",
  "Lipid"                              = "#CCA38D",
  "Nucleotide"                         = "#8EAFCF",
  "Peptide"                            = "#96B485",
  "Xenobiotics"                        = "#814956",
  "Partially Characterized Molecules"  = "#6A4F7E"
)
metab_tax <- as(tax_table(metab), "matrix")

# Figure 3: liver genes × top metabolites heatmap (shows 30 metabolites with the strongest signif associations)
top_metabs <- cor_genes_metab |>
  group_by(y) |>
  summarize(min_fdr = min(p_fdr), max_abs_rho = max(abs(rho)), .groups = "drop") |>
  filter(min_fdr < 0.05) |>
  slice_max(max_abs_rho, n = 30) |>
  pull(y)

metab_sub <- cor_genes_metab |> filter(y %in% top_metabs) |> dplyr::rename(x = y, y = x)  # metabolites as rows
rho_metab <- to_matrix(metab_sub, "rho")
fdr_metab <- to_matrix(metab_sub, "p_fdr")[rownames(rho_metab), colnames(rho_metab)]
metab_gene_categories <- target_genes$category[match(colnames(rho_metab), target_genes$symbol)]

ht_metab <- Heatmap(
  rho_metab,
  name = "Spearman\nrho",
  col = colorRamp2(c(-0.4, 0, 0.4), c(renoir_15[3], "white", renoir_15[11])),
  top_annotation = HeatmapAnnotation(Category = metab_gene_categories, col = list(Category = category_cols), show_annotation_name = FALSE),
  left_annotation = rowAnnotation(
    `Super pathway` = metab_tax[rownames(rho_metab), "SUPER_PATHWAY"],
    col = list(`Super pathway` = pathway_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 8)
  ),
  row_labels = str_trunc(rownames(rho_metab), 45),
  row_names_side = "left",
  row_dend_side = "right",
  column_labels = gene_label(colnames(rho_metab)),
  column_names_max_height = max_text_width(gene_label(colnames(rho_metab))),
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  rect_gp = gpar(col = "white", lwd = 0.7),
  cell_fun = function(j, i, x, y, width, height, fill) grid.text(stars(fdr_metab[i, j]), x, y, gp = gpar(fontsize = 8))
)

pdf("results/graphs/RNAseq/liver_genes_metabolites_heatmap.pdf", width = 9, height = 12)
draw(ht_metab)
dev.off()

#### DE: %FFMI-change groups (1y) ####
coldata <- baria_muscle_wide |>
  mutate(id = as.numeric(as.character(id))) |>
  filter(
    id %in% liver_rnaseq$id[!liver_rnaseq$qc_exclude],
    !is.na(perc_change_ffmi_v4_group), !is.na(age_v0), !is.na(fmi_v0) # DESeq2 can't handle NAs
  ) |>
  mutate(
    ffmi_group_1y = factor(
      if_else(perc_change_ffmi_v4_group == "high", "high_loss", "moderate_low_loss"),
      levels = c("moderate_low_loss", "high_loss") # moderate_low_loss is the reference
    ),
    age_v0_z = as.numeric(scale(age_v0)),  # centre on the final sample set
    fmi_v0_z = as.numeric(scale(fmi_v0))
  ) |>
  arrange(id) |>  # same order as count_data
  select(id, sex, age_v0_z, fmi_v0_z, ffmi_group_1y)

count_data <- liver_rnaseq |>
  filter(!qc_exclude, id %in% coldata$id) |> # QC-excluded samples out (main analysis)
  arrange(id) |>
  column_to_rownames(var = "id") |>
  dplyr::select(starts_with("ENSG")) |> # drops qc_exclude, keeps genes only
  as.matrix() |>
  t() |>
  round()
stopifnot(identical(colnames(count_data), as.character(coldata$id)))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design =  ~ sex + age_v0_z + fmi_v0_z + ffmi_group_1y
)

# Pre-filter: >= 10 normalised counts in >= 50% of samples (same filter as in 0c_datacleaning_rnaseq.R)
dds <- estimateSizeFactors(dds)
keep <- rowMeans(counts(dds, normalized = TRUE) >= 10) >= 0.5
dds <- dds[keep, ]
nrow(dds)  # n of genes tested

# DE analysis
dds <- DESeq(dds)
resultsNames(dds)  # exact coefficient name for the group effect



pdf("results/graphs/RNAseq/liver_de_genes_top15_species_heatmap.pdf", width = 11, height = 2.5 + 0.35 * nrow(rho_de_sp))
draw(ht_de_species)
dev.off()
