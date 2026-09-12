#!/usr/bin/env bash
#SBATCH --job-name=c5g7-newton-cs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=results/slurm-%x-%j.out
#SBATCH --error=results/slurm-%x-%j.err
#SBATCH --export=ALL

set -euo pipefail

# Simulation/Newton defaults are collected in:
#   scripts/critical_rod_newton_cs.env
# Edit that file for persistent defaults, or pass NAME=value arguments below.
# Slurm resources (GPU/CPU/memory/time) are the #SBATCH lines above.
# Example:
# sbatch scripts/run_critical_rod_newton_cs_slurm.sbatch \
#   INITIAL_DEPTH_CM=64.24 PROBE_STEP_CM=0.25 PARTICLES=2000000
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

if [[ ! -f Makefile || ! -f scripts/run_critical_rod_newton_cs.sh ]]; then
  echo "C5G7 project not found at PROJECT_DIR=$PROJECT_DIR" >&2
  echo "Submit from the project root or pass PROJECT_DIR=/absolute/path." >&2
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

# Load the same defaults that the Newton driver will use. Explicit NAME=value
# arguments were exported above, so they continue to take precedence.
CONFIG_FILE="${CONFIG_FILE:-scripts/critical_rod_newton_cs.env}"
if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

export CUDA_ARCH="${CUDA_ARCH:-native}"
export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_MODULE_LOADING="${CUDA_MODULE_LOADING:-LAZY}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export RUN_TAG="${RUN_TAG:-slurm_${SLURM_JOB_ID:-manual}_$(date +%Y%m%d_%H%M%S)}"
export LOG_DIR="${LOG_DIR:-results/critical_rod_newton_cs_${RUN_TAG}}"

mkdir -p results

echo "SLURM job: ${SLURM_JOB_ID:-not-set}"
echo "Node list: ${SLURM_JOB_NODELIST:-not-set}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-not-set}"
echo "Project: $PROJECT_DIR"
echo "CUDA_ARCH: $CUDA_ARCH"
echo "Particles/cycle: ${PARTICLES:-unset}"
echo "Cycles: ${INACTIVE:-unset} inactive + ${ACTIVE:-unset} active"
echo "Threads/block: ${THREADS_PER_BLOCK:-unset}"
echo "Fast math / device LTO: ${FAST_MATH:-unset} / ${DEVICE_LTO:-unset}"
echo "GPU binding: ${GPU_BIND:-closest}"
echo "High clock request: ${REQUEST_HIGH_CLOCKS:-1}"

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi -L
  nvidia-smi --query-gpu=index,name,pstate,clocks.current.sm,clocks.current.memory,power.limit,pcie.link.gen.current,pcie.link.width.current \
    --format=csv || true
fi

# srun creates the actual job step and keeps the process inside Slurm's GPU
# allocation/binding environment.
srun_args=(
  srun
  --ntasks=1
  "--cpus-per-task=${SLURM_CPUS_PER_TASK:-1}"
  "--gpu-bind=${GPU_BIND:-closest}"
)

if [[ "${REQUEST_HIGH_CLOCKS:-1}" == "1" ]]; then
  # These are Slurm requests, not direct nvidia-smi mutations. Unsupported
  # clock controls are reported by Slurm and the hardware/admin limits win.
  srun_args+=(--gpu-freq=high,memory=high,verbose --cpu-freq=high)
elif [[ "${REQUEST_HIGH_CLOCKS:-1}" != "0" ]]; then
  echo "REQUEST_HIGH_CLOCKS must be 0 or 1." >&2
  exit 2
fi

"${srun_args[@]}" bash scripts/run_critical_rod_newton_cs.sh
