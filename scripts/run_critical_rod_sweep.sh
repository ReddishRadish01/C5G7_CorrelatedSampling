#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

EXE="${EXE:-build/c5g7}"
START_CM="${START_CM:-28.0}"
END_CM="${END_CM:-32.0}"
CELL_HEIGHT_CM="${CELL_HEIGHT_CM:-0.21}"
PARTICLES="${PARTICLES:-500000}"
INACTIVE="${INACTIVE:-300}"
ACTIVE="${ACTIVE:-200}"
LOG_DIR="${LOG_DIR:-results/critical_rod_sweep}"
SUMMARY="${SUMMARY:-$LOG_DIR/summary.csv}"
CASE_ARGS="${CASE_ARGS:---critical-search-binary}"

mkdir -p "$LOG_DIR"

make "$EXE"

start_idx="$(awk -v start="$START_CM" -v h="$CELL_HEIGHT_CM" 'BEGIN { print int(start / h + 0.999999999) }')"
end_idx="$(awk -v end="$END_CM" -v h="$CELL_HEIGHT_CM" 'BEGIN { print int(end / h + 0.000000001) }')"

if (( start_idx > end_idx )); then
  echo "No quantized rod depths found in [$START_CM, $END_CM] with cell height $CELL_HEIGHT_CM cm." >&2
  exit 1
fi

echo "depth_cm,actual_depth_cm,k_mean,k_stderr_mean,status,log_file" > "$SUMMARY"
echo "Running CriticalSearch_Binary rod-depth sweep"
echo "  actual depth range: $(awk -v i="$start_idx" -v h="$CELL_HEIGHT_CM" 'BEGIN { printf "%.6f", i*h }') to $(awk -v i="$end_idx" -v h="$CELL_HEIGHT_CM" 'BEGIN { printf "%.6f", i*h }') cm"
echo "  cell height: $CELL_HEIGHT_CM cm"
echo "  particles/inactive/active: $PARTICLES/$INACTIVE/$ACTIVE"
echo "  logs: $LOG_DIR"

for (( idx = start_idx; idx <= end_idx; ++idx )); do
  depth="$(awk -v i="$idx" -v h="$CELL_HEIGHT_CM" 'BEGIN { printf "%.6f", i*h }')"
  depth_tag="${depth//./p}"
  log_file="$LOG_DIR/rod_${depth_tag}cm.log"

  echo "==> rod depth ${depth} cm"
  set +e
  "$EXE" \
    $CASE_ARGS \
    --rod-depth="$depth" \
    --particles="$PARTICLES" \
    --inactive="$INACTIVE" \
    --active="$ACTIVE" \
    > "$log_file" 2>&1
  status=$?
  set -e

  actual="$(awk -F': ' '/Actual (modeled|quantized) rod depth:/ { gsub(/ cm/, "", $2); value=$2 } END { print value }' "$log_file")"
  k_mean="$(awk -F': ' '/^k_mean:/ { value=$2 } END { print value }' "$log_file")"
  k_stderr="$(awk -F': ' '/^k_stderr\(mean\):/ { value=$2 } END { print value }' "$log_file")"

  if (( status == 0 )); then
    echo "$depth,$actual,$k_mean,$k_stderr,OK,$log_file" >> "$SUMMARY"
  else
    echo "    failed with exit code $status; see $log_file" >&2
    echo "$depth,$actual,$k_mean,$k_stderr,FAILED:$status,$log_file" >> "$SUMMARY"
  fi
done

echo "Done. Summary written to $SUMMARY"
