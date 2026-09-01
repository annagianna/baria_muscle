# Baria project: Muscle mass trajectories and gut microbiota following bariatric surgery - Omics data cleaning
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Libraries
library(tidyverse)
library(phyloseq)

# Open cleaned clinical data (see 0a_datacleaning.R) & raw microbiome data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_mb <- readRDS("data/raw_data/ps.BARIA.metaphlan.706.2548.RDS")
bia_ids <- baria_muscle_wide |> pull(id) |> as.character() |> unique()

# Keep single runs and first run of duplicated samples (after comparing 1st and 2nd runs)
run1_mb <- prune_samples(
  (sample_data(baria_mb)$Extra_data == "NA" | sample_data(baria_mb)$Extra_data == "rep1") & 
    !str_detect(sample_data(baria_mb)$Time_Point, "^V\\d+re$"),
  baria_mb
)

# Baseline microbiome sample ids (all subjects with valid baseline mb, before cohort restriction)
mb_v0_ids <- sample_data(run1_mb) |>
  data.frame() |>
  filter(Time_Point == "V-1") |>
  pull(Subject_ID) |>
  as.character() |>
  unique()

# Clean microbiome metadata
sample_data(run1_mb)$visit <- case_when(
  sample_data(run1_mb)$Time_Point == "V-1" ~ "v0",
  sample_data(run1_mb)$Time_Point == "V4" ~ "v4",
  sample_data(run1_mb)$Time_Point == "V5" ~ "v5",
  TRUE ~ NA_character_
)

# Clean sample IDs in compositional table
sample_data(run1_mb)$id <- as.character(sample_data(run1_mb)$Subject_ID)
sample_names(run1_mb) <- str_c("BARIA_", str_remove(str_remove(str_remove(sample_names(run1_mb), "BARIA.Metagenome."), ".Fecal"), ".NA"))
sample_names(run1_mb) <- str_replace(str_replace(str_replace(sample_names(run1_mb), "V-1", "v0"), "V4", "v4"), "V5", "v5")
sample_names(run1_mb) <- str_replace(sample_names(run1_mb), "\\.", "_")
sample_names(run1_mb) <- str_replace(sample_names(run1_mb), ".rep\\d", "") # to remove .rep1

all(sample_names(run1_mb) == str_c("BARIA_", sample_data(run1_mb)$id, "_", sample_data(run1_mb)$visit))
# so we couldve also done sample_names(run1_mb) <- str_c("BARIA_", sample_data(run1_mb)$id, "_", sample_data(run1_mb)$visit)
# but that feels less safe

# Restrict to final Baria muscle cohort
baria_mb_clean <- prune_samples(sample_data(run1_mb)$id %in% bia_ids, run1_mb)

# Add relevant metadata to mb
baria_mb_metadata <- as(sample_data(baria_mb_clean), "data.frame") |>
  rownames_to_column(var = "Sample") |>
  select(-any_of(c("sex", "ffmi_v0", "low_ffmi_v0", "perc_change_ffmi_v4_group", "perc_change_ffmi_v5_group"))) |>
  left_join(
    baria_muscle_wide |>
      select(id, sex, ffmi_v0, low_ffmi_v0, perc_change_ffmi_v4_group, perc_change_ffmi_v5_group),
    by = "id"
  ) |>
  column_to_rownames("Sample")

sample_data(baria_mb_clean) <- sample_data(baria_mb_metadata)

# Add a clean species label to the tax table
tax_table(baria_mb_clean) <- cbind(
  tax_table(baria_mb_clean),
  Tax = tax_table(baria_mb_clean)[, "Species"] |>
    str_remove("^s__") |>
    str_replace_all("_", " ")
)

# Shorten taxa names from the full taxonomy string to "species_SGB####"
species_label <- taxa_names(baria_mb_clean) |>
  str_extract("s__[^|]+") |>
  str_remove("^s__")
sgb <- taxa_names(baria_mb_clean) |>
  str_extract("t__.*$") |>
  str_remove("^t__")
taxa_names(baria_mb_clean) <- if_else(
  str_ends(species_label, sgb), species_label, str_c(species_label, "_", sgb)
)

# Save unfiltered (incl. Eukaryota & GGB-labelled genera) for diversity
baria_mb_unfiltered <- baria_mb_clean

# Drop eukaryotes
baria_mb_clean <- subset_taxa(baria_mb_clean, Kingdom != "k__Eukaryota")
# Drop GGB genera (without proper taxonomy)
baria_mb_clean <- subset_taxa(baria_mb_clean, !str_detect(Genus, "^g__GGB"))

baria_mb_baseline <- prune_samples(baria_mb_clean@sam_data$visit == "v0", baria_mb_clean)
baria_mb_baseline

#### Pathway data cleaning ####
humann <- readRDS("data/raw_data/BARIA.humann4.profiles.2026.581.910.RDS")

humann_long <- humann |>
  rownames_to_column(var = "pathway") |>
  pivot_longer(
    cols = -pathway,
    names_to = "Sample",
    values_to = "pathway_abundance"
  ) |>
  mutate(
    # Raw column names look like "BARIA.Metagenome.<id>.Fecal.V.1.NA_Abundance";
    # rebuild them into the same "BARIA_<id>_<visit>" ids used for the mb
    # sample_names (see the run1_mb renaming above)
    fecal_sample = str_detect(Sample, "\\.Fecal\\."),
    Sample = Sample |>
      str_remove("_Abundance$") |>
      str_replace("V\\.1", "V-1") |>
      str_remove("BARIA.Metagenome.") |>
      str_remove(".Fecal") |>
      str_remove(".NA"),
    Sample = str_c("BARIA_", Sample),
    Sample = Sample |>
      str_replace("V-1", "v0") |>
      str_replace("V4", "v4") |>
      str_replace("V5", "v5") |>
      str_replace("\\.", "_") |>
      str_replace(".rep\\d", "")
  ) |>
  filter(
    fecal_sample,
    Sample %in% sample_names(baria_mb_clean),
    !pathway %in% c("UNMAPPED", "UNINTEGRATED")
  ) |>
  select(-fecal_sample) |>
  mutate(
    pathway_id = str_extract(pathway, "[A-Z0-9-]+(?=:)"),
    pathway_name = str_extract(pathway, "(?<=:).*")
      |> trimws()
  ) |>
  left_join(
    as(sample_data(baria_mb_clean), "data.frame") |>
      rownames_to_column(var = "Sample") |>
      select(Sample, id, visit),
    by = "Sample"
  ) |>
  mutate(pathway_abundance = replace_na(pathway_abundance, 0)) # Treat undetected pathways as zero abundance

#### Metabolomics cleaning ####
## These are phyloseq objects
metab_batchnorm <- readRDS("data/raw_data/BARIA.metabolon.1158.V12.BatchNorm.RDS")
metab_peakarea <- readRDS("data/raw_data/BARIA.metabolon.1158.V12.Peak.Area.RDS")

## Sample IDs
metab_meta <- as(sample_data(metab_batchnorm), "data.frame") |> rownames_to_column("Sample")

## Visit nrs: filter for baseline
metab_baseline <- metab_meta |>
  filter(
    Time_Point == "V-1",
    Type == "Fasted",
    is.na(Exclude) | Exclude != "T",
    is.na(BOX_NUMBER) | !str_detect(BOX_NUMBER, "5 jaar")
  )

## ~36 subjects still have two baseline runs, filter out
metab_baseline <- metab_baseline |>
  group_by(Subject_ID) |>
  filter(!(n() > 1 & Batch == "Batch3")) |>
  ungroup() |>
  mutate(id = as.character(Subject_ID))

# Selection subgroup subjects FFMI and microbiome
metab_ffmi_ids <- baria_muscle_wide |> filter(!is.na(ffmi_v0)) |> pull(id) |> as.character()
metab_target_ids <- intersect(mb_v0_ids, metab_ffmi_ids)

metab_baseline <- metab_baseline |> filter(id %in% metab_target_ids)

metab_batchnorm_v0 <- prune_samples(sample_names(metab_batchnorm) %in% metab_baseline$Sample, metab_batchnorm)
metab_peakarea_v0 <- prune_samples(sample_names(metab_peakarea) %in% metab_baseline$Sample, metab_peakarea)

## Make matrix from phyloseq tables
bn_mat <- as(otu_table(metab_batchnorm_v0), "matrix")
pa_mat <- as(otu_table(metab_peakarea_v0), "matrix")

## Peak is 0 when not detected: % zero == % imputed per metabolite
pct_imputed <- colMeans(pa_mat == 0) * 100
metab_var <- apply(bn_mat, 2, var)

## Drop unidentified compounds
metab_chem_name <- as(tax_table(metab_batchnorm_v0), "matrix")[, "CHEMICAL_NAME"]
metab_named <- names(metab_chem_name)[!str_detect(metab_chem_name, "^X\\s*-\\s*\\d+$")]

## Keep metabolites that are less than 10% imputed
metab_keep <- names(pct_imputed)[pct_imputed <= 10 & metab_var > 0 & names(pct_imputed) %in% metab_named]
bn_keep <- bn_mat[, metab_keep, drop = FALSE]

## Impute + log10-transform + zero mean/unit variance scale
metab_scaled <- apply(bn_keep, 2, function(x) {as.numeric(scale(log10(x + min(x[x > 0]))))})
rownames(metab_scaled) <- rownames(bn_keep)

# Fix rownames and colnames of metabolomics data
metab_tax <- as(tax_table(metab_batchnorm_v0), "matrix")[metab_keep, , drop = FALSE]
metab_tax <- cbind(metab_tax, pct_imputed = round(pct_imputed[metab_keep], 1))

chem_name <- metab_tax[, "CHEMICAL_NAME"]
colnames(metab_scaled) <- chem_name
rownames(metab_tax) <- chem_name

metab_sample_data <- metab_baseline |>
  select(Sample, id, Batch, Subject_ID) |>
  mutate(visit = "v0") |>
  arrange(match(Sample, rownames(metab_scaled))) |>
  column_to_rownames("Sample")
rownames(metab_scaled) <- str_c("BARIA_", metab_sample_data$id, "_v0")
rownames(metab_sample_data) <- rownames(metab_scaled)

baria_metab_clean <- phyloseq( # questionable whether this is needed
  otu_table(metab_scaled, taxa_are_rows = FALSE),
  sample_data(metab_sample_data),
  tax_table(metab_tax)
)
baria_metab_clean

# Save microbiome, pathway & metabolomics data as RDS files
dir.create("data/processed_data", recursive = TRUE, showWarnings = FALSE)

# Save mb data
saveRDS(baria_mb_clean, "data/processed_data/BARIA_mb_clean.RDS")
saveRDS(baria_mb_unfiltered, "data/processed_data/BARIA_mb_clean_unfiltered.RDS")
saveRDS(baria_mb_baseline, "data/processed_data/BARIA_mb_baseline.RDS")

# Save pathway data
saveRDS(humann_long, "data/processed_data/BARIA_humann_pathways_long.RDS")

# Save metabolomics data
saveRDS(baria_metab_clean, "data/processed_data/BARIA_metabolon_clean.RDS")
