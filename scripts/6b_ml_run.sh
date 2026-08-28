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

# Every outcome produced by 6a_ml_prep.R, each split into all/male/female
OUTCOMES=(
    ffmi
    ffmi_adj_fmi
    fmi_v0
    ffmi_v4
    fmi_v4
    delta_ffmi_v4
    delta_ffmi_v4_adj_fmi
    perc_change_ffmi_v4
    perc_change_ffmi_v4_adj_fmi
)
SUBGROUPS=(all male female)

for outcome in "${OUTCOMES[@]}"; do
    for group in "${SUBGROUPS[@]}"; do
        pixi run python $SCRIPT \
            -name "${outcome}_${group}" \
            -path "$BASE/results/mlmodels/${outcome}/${group}" \
            -x reg \
            -n 200 \
            -t 8 \
            -rand_seed 1234 \
            -param $PARAM
    done
done
