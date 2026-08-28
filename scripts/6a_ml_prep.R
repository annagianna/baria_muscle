# Prep XGBeast input data: microbiome species ~ FFMI
# Barbara Verhaar

library(stringr)
source("scripts/assets/functions.R")

meta <- readRDS("data/processed_data/BARIA_muscle_long.RDS")
meta$sampleid <- str_c("BARIA_", meta$id, "_", meta$visit)
dim(meta)
names(meta)
cat("FFMI: ", sum(!is.na(meta[["ffmi"]])), "/", nrow(meta), "non-missing\n")

meta <- meta |> filter(visit == "v0")
baseline_fmi <- meta |> select(id, fmi_v0 = fmi)

meta_var <- meta |> filter(!is.na(ffmi))
meta_var_male <- meta_var |> filter(sex == "male")
meta_var_female <- meta_var |> filter(sex == "female")

mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb <- t(as(mb@otu_table, "matrix"))
tk <- apply(mb, 2, function(x) sum(x > 0.05) > (0.2 * length(x)))
mb2 <- as.data.frame(mb[, tk])
dim(mb2) # 147 species left, 365 subjects
mb2[1:5,1:5]

base <- file.path("results/mlmodels/ffmi")
build_input_data(meta_var,        mb2, "ffmi", file.path(base, "all"),    mode = "reg")
build_input_data(meta_var_male,   mb2, "ffmi", file.path(base, "male"),   mode = "reg")
build_input_data(meta_var_female, mb2, "ffmi", file.path(base, "female"), mode = "reg")

#### FFMI corrected for FMI ####
# Outcome: FFMI residualized on baseline FMI (fat mass index), i.e. muscle
# mass independent of adiposity; predictors: baseline (v0) species (mb2)
meta_var_fmi <- meta_var |> filter(!is.na(fmi)) |>
  add_residual_var("ffmi", "fmi", "ffmi_adj_fmi")
meta_var_fmi_male <- meta_var_male |> filter(!is.na(fmi)) |>
  add_residual_var("ffmi", "fmi", "ffmi_adj_fmi")
meta_var_fmi_female <- meta_var_female |> filter(!is.na(fmi)) |>
  add_residual_var("ffmi", "fmi", "ffmi_adj_fmi")

base_adj <- file.path("results/mlmodels/ffmi_adj_fmi")
build_input_data(meta_var_fmi,        mb2, "ffmi_adj_fmi", file.path(base_adj, "all"),    mode = "reg")
build_input_data(meta_var_fmi_male,   mb2, "ffmi_adj_fmi", file.path(base_adj, "male"),   mode = "reg")
build_input_data(meta_var_fmi_female, mb2, "ffmi_adj_fmi", file.path(base_adj, "female"), mode = "reg")

#### Delta FFMI v0 -> v4 ####
# Outcome: change in FFMI from baseline to v4; predictors: baseline (v0) species (mb2)
meta_delta <- readRDS("data/processed_data/BARIA_muscle_long.RDS") |> filter(visit == "v4")
meta_delta$sampleid <- str_c("BARIA_", meta_delta$id, "_v0") # match against baseline mb sample ids
meta_delta <- meta_delta |> left_join(baseline_fmi, by = "id")
cat("delta_ffmi_v4: ", sum(!is.na(meta_delta[["delta_ffmi"]])), "/", nrow(meta_delta), "non-missing\n")

meta_delta_var <- meta_delta |> filter(!is.na(delta_ffmi))
meta_delta_var_male <- meta_delta_var |> filter(sex == "male")
meta_delta_var_female <- meta_delta_var |> filter(sex == "female")

base_delta <- file.path("results/mlmodels/delta_ffmi_v4")
build_input_data(meta_delta_var,        mb2, "delta_ffmi", file.path(base_delta, "all"),    mode = "reg")
build_input_data(meta_delta_var_male,   mb2, "delta_ffmi", file.path(base_delta, "male"),   mode = "reg")
build_input_data(meta_delta_var_female, mb2, "delta_ffmi", file.path(base_delta, "female"), mode = "reg")

#### Delta FFMI corrected for FMI ####
# Outcome: delta FFMI residualized on baseline FMI; predictors: baseline (v0) species (mb2)
meta_delta_var_fmi <- meta_delta_var |> filter(!is.na(fmi_v0)) |>
  add_residual_var("delta_ffmi", "fmi_v0", "delta_ffmi_adj_fmi")
meta_delta_var_fmi_male <- meta_delta_var_male |> filter(!is.na(fmi_v0)) |>
  add_residual_var("delta_ffmi", "fmi_v0", "delta_ffmi_adj_fmi")
meta_delta_var_fmi_female <- meta_delta_var_female |> filter(!is.na(fmi_v0)) |>
  add_residual_var("delta_ffmi", "fmi_v0", "delta_ffmi_adj_fmi")

base_delta_adj <- file.path("results/mlmodels/delta_ffmi_v4_adj_fmi")
build_input_data(meta_delta_var_fmi,        mb2, "delta_ffmi_adj_fmi", file.path(base_delta_adj, "all"),    mode = "reg")
build_input_data(meta_delta_var_fmi_male,   mb2, "delta_ffmi_adj_fmi", file.path(base_delta_adj, "male"),   mode = "reg")
build_input_data(meta_delta_var_fmi_female, mb2, "delta_ffmi_adj_fmi", file.path(base_delta_adj, "female"), mode = "reg")

#### Percentual change FFMI v0 -> v4 ####
# Outcome: percentual change in FFMI from baseline to v4; predictors: baseline (v0) species (mb2)
cat("perc_change_ffmi_v4: ", sum(!is.na(meta_delta[["perc_change_ffmi"]])), "/", nrow(meta_delta), "non-missing\n")

meta_pchange_var <- meta_delta |> filter(!is.na(perc_change_ffmi))
meta_pchange_var_male <- meta_pchange_var |> filter(sex == "male")
meta_pchange_var_female <- meta_pchange_var |> filter(sex == "female")

base_pchange <- file.path("results/mlmodels/perc_change_ffmi_v4")
build_input_data(meta_pchange_var,        mb2, "perc_change_ffmi", file.path(base_pchange, "all"),    mode = "reg")
build_input_data(meta_pchange_var_male,   mb2, "perc_change_ffmi", file.path(base_pchange, "male"),   mode = "reg")
build_input_data(meta_pchange_var_female, mb2, "perc_change_ffmi", file.path(base_pchange, "female"), mode = "reg")

#### Percentual change FFMI corrected for FMI ####
# Outcome: percentual change FFMI residualized on baseline FMI; predictors: baseline (v0) species (mb2)
meta_pchange_var_fmi <- meta_pchange_var |> filter(!is.na(fmi_v0)) |>
  add_residual_var("perc_change_ffmi", "fmi_v0", "perc_change_ffmi_adj_fmi")
meta_pchange_var_fmi_male <- meta_pchange_var_male |> filter(!is.na(fmi_v0)) |>
  add_residual_var("perc_change_ffmi", "fmi_v0", "perc_change_ffmi_adj_fmi")
meta_pchange_var_fmi_female <- meta_pchange_var_female |> filter(!is.na(fmi_v0)) |>
  add_residual_var("perc_change_ffmi", "fmi_v0", "perc_change_ffmi_adj_fmi")

base_pchange_adj <- file.path("results/mlmodels/perc_change_ffmi_v4_adj_fmi")
build_input_data(meta_pchange_var_fmi,        mb2, "perc_change_ffmi_adj_fmi", file.path(base_pchange_adj, "all"),    mode = "reg")
build_input_data(meta_pchange_var_fmi_male,   mb2, "perc_change_ffmi_adj_fmi", file.path(base_pchange_adj, "male"),   mode = "reg")
build_input_data(meta_pchange_var_fmi_female, mb2, "perc_change_ffmi_adj_fmi", file.path(base_pchange_adj, "female"), mode = "reg")
