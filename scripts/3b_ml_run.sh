#!/bin/bash
#SBATCH --job-name=bariaml
#SBATCH --partition=genoa
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --output=bariaml_%j.out
#SBATCH --error=bariaml_%j.err

BASE=$PWD
cd "$BASE"
export PATH="$HOME/.pixi/bin:$PATH"   # SLURM batch shells don't source ~/.bashrc
SCRIPT=$BASE/scripts/assets/XGBeast_new.py
PARAM=$BASE/scripts/assets/param_grid.json

# Predictor set (results tree from 3a_ml_prep.R):
#   sbatch scripts/3b_ml_run.sh          -> species     (results/mlmodels)
#   sbatch scripts/3b_ml_run.sh metab    -> metabolites (results/mlmodels_metab)
DATASET=${1:-species}
case "$DATASET" in
    species) ROOT=$BASE/results/mlmodels ;;
    metab)   ROOT=$BASE/results/mlmodels_metab ;;
    *) echo "Unknown dataset '$DATASET' (use: species | metab)" >&2; exit 1 ;;
esac

# Every outcome produced by 3a_ml_prep.R, each split into all/male/female
OUTCOMES=(
    perc_change_ffmi_v4_group
    perc_change_ffmi_v4_group_adj_fmi
)
SUBGROUPS=(all male female)

for outcome in "${OUTCOMES[@]}"; do
    for group in "${SUBGROUPS[@]}"; do
        pixi run python $SCRIPT \
            -name "${outcome}_${group}" \
            -path "$ROOT/${outcome}/${group}" \
            -x class \
            -n 200 \
            -t 6 \
            -rand_seed 1234 \
            -param $BASE/scripts/assets/param_grid.json
    done
done
