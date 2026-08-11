# DAGs for associations between species abundance and FFMI
# Anna Giannakogeorgou, a.gianna@amsterdamumc.nl

# Packages
library(dagitty)
library(ggdag)

### DAGs ###
# Model 1: Primary DAG
dag_ffmi_v0_1 <- dagitty::dagitty("
dag {
microbiome [exposure]
FFMI [outcome]

age -> microbiome
age -> FFMI

sex -> microbiome
sex -> FFMI

microbiome -> metabolites
metabolites -> FFMI
}")

dagitty::adjustmentSets(
  dag_ffmi_v0_1,
  exposure = "microbiome",
  outcome = "FFMI",
  type = "minimal",
  effect = "total"
)

dagitty::paths(
  dag_ffmi_v0_1,
  from = "microbiome",
  to = "FFMI"
)

# Model 2: Sensitivity DAG (with adjustment for adiposity/FMI)
dag_ffmi_v0_2s <- dagitty::dagitty("
dag {
microbiome [exposure]
FFMI [outcome]

age -> microbiome
age -> FFMI

sex -> microbiome
sex -> FFMI

adiposity -> microbiome
adiposity -> FFMI

microbiome -> metabolites
metabolites -> FFMI
}")

dagitty::adjustmentSets(
  dag_ffmi_v0_2s,
  exposure = "microbiome",
  outcome = "FFMI",
  type = "minimal",
  effect = "total"
)

# Model 3: Extensive sensitivity DAG
# Additional adjustment for T2D, antidiabetic medication and statins
dag_ffmi_v0_3s <- dagitty::dagitty("
dag {
species_abundance [exposure]
FFMI [outcome]

age -> species_abundance
age -> FFMI
age -> T2D
age -> statins

sex -> species_abundance
sex -> FFMI
sex -> T2D

adiposity -> species_abundance
adiposity -> FFMI
adiposity -> T2D
adiposity -> antidiabetic_drugs
adiposity -> statins

T2D -> species_abundance
T2D -> FFMI
T2D -> antidiabetic_drugs
T2D -> statins

antidiabetic_drugs -> species_abundance
antidiabetic_drugs -> FFMI

statins -> species_abundance
statins -> FFMI

species_abundance -> FFMI
}
")

dagitty::adjustmentSets(
  dag_ffmi_v0_3s,
  exposure = "species_abundance",
  outcome = "FFMI",
  type = "minimal",
  effect = "total"
)
