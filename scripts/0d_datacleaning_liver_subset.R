# Baria project: Liver subcohort (participants with available BIA, baseline shotgun, liver histology and RNA-seq)
# Anna Giannakogeorgou

# Libraries
library(tidyverse)

# Data
baria_muscle_wide <- readRDS("data/processed_data/BARIA_muscle_wide.RDS")
baria_muscle_long <- readRDS("data/processed_data/BARIA_muscle_long.RDS")
liver_rnaseq <- readRDS("data/processed_data/BARIA_Liver_RNAseq.RDS")
liver_rnaseq_vst <- readRDS("data/processed_data/BARIA_Liver_RNAseq_vst.RDS")
liver_rnaseq_vst_bc <- readRDS("data/processed_data/BARIA_Liver_RNAseq_vst_batchcorrected.RDS")
baria_mb_v0 <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

#### IDs per data layer ####
# 1. Valid baseline BIA + baseline shotgun
bia_mb_ids <- baria_muscle_wide |> 
  pull(id) |> 
  as.character() |> 
  as.numeric()

# Check: every cohort ID has a baseline shotgun sample
mb_v0_ids <- as(phyloseq::sample_data(baria_mb_v0), "data.frame") |> 
  pull(id) |> 
  as.numeric()
stopifnot(all(bia_mb_ids %in% mb_v0_ids))

# 2. Baseline liver histology available
histo_ids <- baria_muscle_wide |>
  filter(!is.na(masld_category_histo_v0)) |>
  pull(id) |> 
  as.character() |> 
  as.numeric()

# 3. Baseline liver RNA-seq passing QC (from 0c RNAseq cleaning script)
liver_rnaseq_ids <- liver_rnaseq |>
  filter(!qc_exclude) |>
  pull(id)

# Subcohort: all three
liver_subcohort_ids <- reduce(list(bia_mb_ids, histo_ids, liver_rnaseq_ids), intersect)

#### Flowchart of participants (for methods) ####
liver_flow <- tibble(
  step = c(
    "Valid baseline BIA + baseline shotgun (analysis cohort)",
    "+ baseline liver histology",
    "+ baseline liver RNA-seq (passing QC)",
    "+ 1-year %FFMI change group available"
  ),
  n = c(
    length(bia_mb_ids),
    length(intersect(bia_mb_ids, histo_ids)),
    length(liver_subcohort_ids),
    baria_muscle_wide |>
      filter(id %in% liver_subcohort_ids, !is.na(perc_change_ffmi_v4)) |>
      nrow()
  )
)
liver_flow

# Who drops out where (histology without RNA-seq, and vice versa)
tibble(
  histo_no_rnaseq = length(setdiff(intersect(bia_mb_ids, histo_ids), liver_rnaseq_ids)),
  rnaseq_no_histo = length(setdiff(intersect(bia_mb_ids, liver_rnaseq_ids), histo_ids))
)

#### Wide data ####
# Sex-specific %FFMI change tertiles recomputed WITHIN the subcohort
# full-cohort cut-offs and groups from 0a kept as *_fullcohort for comparison
baria_muscle_wide_liver <- baria_muscle_wide |>
  filter(id %in% liver_subcohort_ids) |>
  rename_with(~ paste0(.x, "_fullcohort"), starts_with("perc_change_ffmi_v") & ends_with(c("_tertile", "_group"))) |>
  group_by(sex) |>
  mutate(
    across(
      starts_with("perc_change_ffmi_v") & !ends_with(c("_tertile", "_group", "_fullcohort")),
      ~ quantile(.x, probs = 1/3, na.rm = TRUE),
      .names = "{.col}_tertile"
    ),
    across(
      starts_with("perc_change_ffmi_v") & !ends_with(c("_tertile", "_group", "_fullcohort")),
      ~ if_else(.x <= quantile(.x, probs = 1/3, na.rm = TRUE), "high", "modest/low"),
      .names = "{.col}_group"
    )
  ) |>
  ungroup() |>
  mutate(across(ends_with("_group") & !ends_with("_fullcohort"), as.factor))

## Checks
# Cut-offs per sex: full cohort vs. subcohort
baria_muscle_wide_liver |> distinct(sex, perc_change_ffmi_v4_tertile_fullcohort, perc_change_ffmi_v4_tertile)

# Participants who switch group with the subcohort cut-offs
baria_muscle_wide_liver |>
  count(perc_change_ffmi_v4_group_fullcohort, perc_change_ffmi_v4_group)

# Group sizes and histology by 1-year %FFMI change group
baria_muscle_wide_liver |> count(sex, perc_change_ffmi_v4_group)
baria_muscle_wide_liver |> count(perc_change_ffmi_v4_group, masld_category_histo_v0)

#### Long data ####
baria_muscle_long_liver <- baria_muscle_long |>
  filter(id %in% liver_subcohort_ids)

#### Liver RNA-seq: subcohort ####
# Rows subset only; VST and batch correction as computed in 0c on all QC-passing liver samples
liver_rnaseq_sub <- liver_rnaseq |> filter(id %in% liver_subcohort_ids)
liver_rnaseq_vst_sub <- liver_rnaseq_vst |> filter(id %in% liver_subcohort_ids)
liver_rnaseq_vst_bc_sub <- liver_rnaseq_vst_bc |> filter(id %in% liver_subcohort_ids)

# Checks
stopifnot(
  !anyDuplicated(baria_muscle_wide_liver$id),
  nrow(baria_muscle_wide_liver) == length(liver_subcohort_ids),
  setequal(liver_rnaseq_sub |> filter(!qc_exclude) |> pull(id), liver_subcohort_ids),
  setequal(liver_rnaseq_vst_sub$id, liver_subcohort_ids)
)

#### Save ####
write_csv(liver_flow, "results/tables/liver_subcohort_flow.csv")
write_csv(tibble(id = sort(liver_subcohort_ids)), "data/processed_data/BARIA_liver_subcohort_ids.csv")

# Metadata & liver histology
saveRDS(baria_muscle_wide_liver, "data/processed_data/BARIA_muscle_wide_livercohort.RDS")
saveRDS(baria_muscle_long_liver, "data/processed_data/BARIA_muscle_long_livercohort.RDS")

# Liver RNAseq data
saveRDS(liver_rnaseq_sub, "data/processed_data/BARIA_Liver_RNAseq_livercohort.RDS")
saveRDS(liver_rnaseq_vst_sub, "data/processed_data/BARIA_Liver_RNAseq_vst_livercohort.RDS")
saveRDS(liver_rnaseq_vst_bc_sub, "data/processed_data/BARIA_Liver_RNAseq_vst_batchcorrected_livercohort.RDS")