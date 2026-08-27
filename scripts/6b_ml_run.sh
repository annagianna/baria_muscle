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

pixi run python $SCRIPT \
    -name ffmi_all \
    -path $BASE/results/mlmodels/ffmi/all \
    -x reg \
    -n 200 \
    -t 8 \
    -rand_seed 1234 \
    -param $PARAM

pixi run python $SCRIPT \
    -name ffmi_male \
    -path $BASE/results/mlmodels/ffmi/male \
    -x reg \
    -n 200 \
    -t 8 \
    -rand_seed 1234 \
    -param $PARAM

pixi run python $SCRIPT \
    -name ffmi_female \
    -path $BASE/results/mlmodels/ffmi/female \
    -x reg \
    -n 200 \
    -t 8 \
    -rand_seed 1234 \
    -param $PARAM
