#!/usr/bin/env bash
# Submit the additional-depth CS validation to all six GPUs on node3.
# Run from the C5G7_CorrelatedSampling project root:
#   bash scripts/submit_cs_prediction_validation_more_points_node3.sh
# Optional overrides can be appended, for example:
#   bash scripts/submit_cs_prediction_validation_more_points_node3.sh NO_BUILD=1

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"
cd "$PROJECT_DIR"

exec sbatch scripts/run_cs_prediction_validation_slurm.sbatch \
  CONFIG_FILE=scripts/cs_prediction_validation_more_points_slurm.env \
  "$@"
