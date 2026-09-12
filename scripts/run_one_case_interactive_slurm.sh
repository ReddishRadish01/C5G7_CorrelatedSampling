#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
SCRIPT_PATH="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"

if [[ "${1:-}" != "--inside-slurm-step" ]]; then
  if [[ ! -t 0 ]]; then
    echo "This launcher needs an interactive terminal." >&2
    echo "Run it directly from the server shell; do not submit it with sbatch." >&2
    exit 2
  fi

  if ! command -v srun >/dev/null 2>&1; then
    echo "srun was not found. Run this script on a Slurm login server." >&2
    exit 1
  fi

  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    exec srun --pty bash "$SCRIPT_PATH" --inside-slurm-step
  fi

  srun_args=(
    srun
    "--job-name=${SLURM_JOB_NAME:-c5g7-one-case}"
    "--gres=${SLURM_GRES:-gpu:1}"
    "--cpus-per-task=${SLURM_CPUS_PER_TASK:-4}"
    "--mem=${SLURM_MEM:-16G}"
    "--time=${SLURM_TIME:-24:00:00}"
  )

  if [[ -n "${SLURM_PARTITION:-}" ]]; then
    srun_args+=("--partition=$SLURM_PARTITION")
  fi
  if [[ -n "${SLURM_ACCOUNT:-}" ]]; then
    srun_args+=("--account=$SLURM_ACCOUNT")
  fi
  if [[ -n "${SLURM_QOS:-}" ]]; then
    srun_args+=("--qos=$SLURM_QOS")
  fi
  if [[ -n "${SLURM_CONSTRAINT:-}" ]]; then
    srun_args+=("--constraint=$SLURM_CONSTRAINT")
  fi

  echo "Requesting an interactive Slurm GPU session..."
  exec "${srun_args[@]}" --pty bash "$SCRIPT_PATH" --inside-slurm-step
fi

shift
cd "$PROJECT_DIR"

if [[ ! -f Makefile || ! -f main.cu ]]; then
  echo "C5G7 project not found at PROJECT_DIR=$PROJECT_DIR" >&2
  exit 2
fi

if [[ -n "${CUDA_MODULE:-}" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if ! type module >/dev/null 2>&1; then
    echo "CUDA_MODULE=$CUDA_MODULE was requested, but module is unavailable." >&2
    exit 1
  fi
  module load "$CUDA_MODULE"
fi

echo "Slurm job: ${SLURM_JOB_ID:-not-set}"
echo "Node: ${SLURMD_NODENAME:-${SLURM_JOB_NODELIST:-not-set}}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-not-set}"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi -L
fi

echo "Building the current source..."
make \
  CUDA_ARCH="${CUDA_ARCH:-native}" \
  THREADS_PER_BLOCK="${THREADS_PER_BLOCK:-32}" \
  DEVICE_LTO="${DEVICE_LTO:-1}" \
  FAST_MATH="${FAST_MATH:-0}"

echo
echo "Available runCase values:"
echo "  unrodded, rodded-a, rodded-b, reference2003, critical, aro, test"

while true; do
  read -r -p "runCase [reference2003]: " run_case
  run_case="${run_case:-reference2003}"
  case_key="${run_case,,}"
  case_key="${case_key//[-_]/}"

  case "$case_key" in
    unrodded)
      run_case="unrodded"
      break
      ;;
    roddeda)
      run_case="rodded-a"
      break
      ;;
    roddedb)
      run_case="rodded-b"
      break
      ;;
    reference|reference2003)
      run_case="reference2003"
      break
      ;;
    critical|criticalsearchbinary)
      run_case="critical"
      break
      ;;
    aro)
      run_case="aro"
      break
      ;;
    test)
      run_case="test"
      break
      ;;
    *)
      echo "Unknown runCase '$run_case'. Please choose one of the values above." >&2
      ;;
  esac
done

run_args=(--case "$run_case")
if [[ "$run_case" == "critical" ]]; then
  read -r -p "Rod insertion depth from top in cm [0.0]: " rod_depth
  run_args+=(--rod-depth "${rod_depth:-0.0}")
fi

echo
printf 'Running: ./build/c5g7'
printf ' %q' "${run_args[@]}"
printf '\n\n'
exec ./build/c5g7 "${run_args[@]}"
