#!/usr/bin/env bash
#SBATCH --job-name=c5g7-cs-fullcurve
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=node3
#SBATCH --gres=gpu:6
#SBATCH --cpus-per-task=12
#SBATCH --hint=nomultithread
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err
#SBATCH --export=ALL

# Dedicated Slurm entry point for the additional ARO-to-ARI curve points.
# The common runner performs compilation, GPU-ID verification, task dispatch,
# monitoring, and postprocessing. Optional NAME=value overrides passed after
# this .sbatch filename are forwarded unchanged.

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"

export CONFIG_FILE="scripts/cs_prediction_validation_more_points_slurm.env"
exec bash scripts/run_cs_prediction_validation_slurm.sbatch "$@"
