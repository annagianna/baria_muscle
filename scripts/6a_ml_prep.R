# Prep XGBeast input data: microbiome species ~ FFMI
# Barbara Verhaar

source("scripts/assets/functions.R")

meta <- readRDS("data/processed_data/BARIA_muscle_long.RDS")
meta$sampleid <- str_c("BARIA_", meta$id, "_", meta$visit)
dim(meta)
names(meta)
cat("FFMI: ", sum(!is.na(meta[["ffmi"]])), "/", nrow(meta), "non-missing\n")

meta <- meta |> filter(visit == "v0")
meta_var <- meta |> filter(!is.na(ffmi))
meta_var_male <- meta_var |> filter(sex == "male")
meta_var_female <- meta_var |> filter(sex == "female")

mb <- readRDS("data/processed_data/BARIA_mb_baseline.RDS")
mb <- t(as(mb@otu_table, "matrix"))
tk <- apply(mb, 2, function(x) sum(x > 0.05) > (0.3 * length(x)))
mb2 <- as.data.frame(mb[, tk])
dim(mb2)

base <- file.path("results/mlmodels/ffmi")
build_input_data(meta_var,        mb2, "ffmi", file.path(base, "all"),    mode = "reg")
build_input_data(meta_var_male,   mb2, "ffmi", file.path(base, "male"),   mode = "reg")
build_input_data(meta_var_female, mb2, "ffmi", file.path(base, "female"), mode = "reg")
