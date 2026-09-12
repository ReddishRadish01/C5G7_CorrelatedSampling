#!/usr/bin/env bash
#SBATCH --job-name=c5g7-bisect-resume
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err
#SBATCH --export=ALL

set -euo pipefail

# Example:
# sbatch scripts/run_critical_rod_binary_search_resume_slurm.sbatch \
#   RESUME_DIR=results/critical_rod_binary_search_slurm_4606_20260910_150206
for setting in "$@"; do
  if [[ "$setting" =~ ^[A-Z_][A-Z0-9_]*=.*$ ]]; then
    export "$setting"
  else
    echo "Invalid setting '$setting'; expected NAME=value." >&2
    exit 2
  fi
done

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$PROJECT_DIR"

if [[ ! -f Makefile || ! -f scripts/run_critical_rod_binary_search_resume.sh ]]; then
  echo "C5G7 project or resume script not found at PROJECT_DIR=$PROJECT_DIR" >&2
  exit 2
fi
if [[ -z "${RESUME_DIR:-}" ]]; then
  echo "RESUME_DIR is required; point it to the existing eval_01...eval_13 result directory." >&2
  exit 2
fi

if [[ -n "${CUDA_MODULE:-}" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if ! type module >/dev/null 2>&1; then
    echo "CUDA_MODULE=$CUDA_MODULE was requested, but the module command is unavailable." >&2
    exit 1
  fi
  module load "$CUDA_MODULE"
fi

export CUDA_ARCH="${CUDA_ARCH:-native}"
export RUN_TAG="${RUN_TAG:-slurm_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)}"
export LOG_DIR="${LOG_DIR:-results/critical_rod_binary_search_resume_${RUN_TAG}}"

echo "SLURM job: ${SLURM_JOB_ID:-not-set}"
echo "Node list: ${SLURM_JOB_NODELIST:-not-set}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-not-set}"
echo "Project: $PROJECT_DIR"
echo "Resume source: $RESUME_DIR"

if command -v nvidia-smi >/dev/null 2>&1; then nvidia-smi -L; fi

bash scripts/run_critical_rod_binary_search_resume.sh
