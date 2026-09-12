#!/usr/bin/env bash
set -uo pipefail

# Continue an existing CriticalSearch_Binary run without repeating its ARO check
# or already-sampled bracket endpoints.  The next evaluation number, current
# bracket, target, seed base, and statistical settings are read from the prior
# result directory.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

RESUME_DIR="${RESUME_DIR:-${1:-}}"
if (( $# > 1 )); then
  echo "Usage: RESUME_DIR=results/<old-run> bash $0" >&2
  echo "   or: bash $0 results/<old-run>" >&2
  exit 2
fi
if [[ -z "$RESUME_DIR" ]]; then
  echo "RESUME_DIR is required." >&2
  echo "Example: RESUME_DIR=results/critical_rod_binary_search_slurm_4606_20260910_150206 bash $0" >&2
  exit 2
fi

if [[ "$RESUME_DIR" = /* ]]; then
  SOURCE_DIR="$RESUME_DIR"
else
  SOURCE_DIR="$ROOT_DIR/$RESUME_DIR"
fi
SOURCE_SUMMARY="$SOURCE_DIR/summary.csv"
SOURCE_REPLICATES="$SOURCE_DIR/replicates.csv"
SOURCE_RESULT="$SOURCE_DIR/result.txt"

for source_file in "$SOURCE_SUMMARY" "$SOURCE_REPLICATES" "$SOURCE_RESULT"; do
  if [[ ! -f "$source_file" ]]; then
    echo "Required resume file not found: $source_file" >&2
    exit 2
  fi
done

result_value() {
  local key="$1"
  awk -v wanted="$key" '{
    for (i=1; i<=NF; ++i) {
      split($i, pair, "=")
      if (pair[1] == wanted) {
        sub(/^[^=]*=/, "", $i)
        print $i
        exit
      }
    }
  }' "$SOURCE_RESULT"
}

is_number() {
  awk -v value="$1" 'BEGIN {
    exit !(value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/)
  }'
}

is_nonnegative_integer() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

is_positive_integer() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

# Find the last successful old log.  This lets the continuation inherit the
# particle/cycle settings actually used, even though the old result.txt did not
# record them explicitly.
last_old_log_entry="$(awk -F, 'NR>1 && $8 == "OK" { value=$9 } END { print value }' "$SOURCE_REPLICATES")"
LAST_OLD_LOG=""
if [[ -n "$last_old_log_entry" ]]; then
  if [[ "$last_old_log_entry" = /* && -f "$last_old_log_entry" ]]; then
    LAST_OLD_LOG="$last_old_log_entry"
  elif [[ -f "$ROOT_DIR/$last_old_log_entry" ]]; then
    LAST_OLD_LOG="$ROOT_DIR/$last_old_log_entry"
  elif [[ -f "$SOURCE_DIR/$(basename "$last_old_log_entry")" ]]; then
    LAST_OLD_LOG="$SOURCE_DIR/$(basename "$last_old_log_entry")"
  fi
fi

detected_particles=""
detected_active=""
detected_inactive=""
if [[ -n "$LAST_OLD_LOG" ]]; then
  detected_particles="$(awk -F'initial Weight: ' '/initial Weight:/ {
    split($2, value, ","); printf "%.0f", value[1]; exit
  }' "$LAST_OLD_LOG")"
  detected_active="$(awk -F': ' '/^Active cycles:/ { value=$2 } END {
    gsub(/[^0-9].*/, "", value); print value
  }' "$LAST_OLD_LOG")"
  detected_total_cycles="$(awk '/^Cycle [0-9]+,/ {
    value=$2; gsub(/,/, "", value)
  } END { print value }' "$LAST_OLD_LOG")"
  if is_nonnegative_integer "${detected_total_cycles:-x}" && is_nonnegative_integer "${detected_active:-x}" && \
     (( detected_total_cycles >= detected_active )); then
    detected_inactive=$((detected_total_cycles - detected_active))
  fi
fi

source_target="$(result_value target_k)"
source_replicates="$(result_value replicates_per_depth)"
source_base_seed="$(result_value base_seed)"
source_confidence_sigma="$(result_value confidence_sigma)"
source_stop_on_ci="$(result_value stop_on_confidence_interval)"

EXE="${EXE:-build/c5g7}"
NVCC="${NVCC:-nvcc}"
TARGET_K="${TARGET_K:-${source_target:-1.0}}"
DEPTH_TOL_CM="${DEPTH_TOL_CM:-0.01}"
K_TOL_PCM="${K_TOL_PCM:-0.0}"
REPLICATES_PER_DEPTH="${REPLICATES_PER_DEPTH:-${source_replicates:-1}}"
MAX_REPLICATE_FAILURES="${MAX_REPLICATE_FAILURES:-2}"
BASE_SEED="${BASE_SEED:-${source_base_seed:-92235922383}}"
CONFIDENCE_SIGMA="${CONFIDENCE_SIGMA:-${source_confidence_sigma:-2.0}}"
STOP_ON_CONFIDENCE_INTERVAL="${STOP_ON_CONFIDENCE_INTERVAL:-${source_stop_on_ci:-0}}"
SCATTER_WARNING_RATIO="${SCATTER_WARNING_RATIO:-2.0}"
MAX_ADDITIONAL_ITERATIONS="${MAX_ADDITIONAL_ITERATIONS:-20}"
PARTICLES="${PARTICLES:-${detected_particles:-10000000}}"
INACTIVE="${INACTIVE:-${detected_inactive:-150}}"
ACTIVE="${ACTIVE:-${detected_active:-200}}"
CUDA_ARCH="${CUDA_ARCH:-native}"
THREADS_PER_BLOCK="${THREADS_PER_BLOCK:-32}"
DEVICE_LTO="${DEVICE_LTO:-1}"
FAST_MATH="${FAST_MATH:-0}"
NO_BUILD="${NO_BUILD:-0}"
CASE_ARGS="${CASE_ARGS:---critical-search-binary}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-results/critical_rod_binary_search_resume_${RUN_TAG}}"
SUMMARY="${SUMMARY:-$LOG_DIR/summary.csv}"
REPLICATE_SUMMARY="${REPLICATE_SUMMARY:-$LOG_DIR/replicates.csv}"
RESULT="${RESULT:-$LOG_DIR/result.txt}"

for value in "$TARGET_K" "$DEPTH_TOL_CM" "$K_TOL_PCM" "$CONFIDENCE_SIGMA" "$SCATTER_WARNING_RATIO"; do
  if ! is_number "$value"; then
    echo "Invalid numeric setting: $value" >&2
    exit 2
  fi
done
if ! awk -v value="$DEPTH_TOL_CM" 'BEGIN { exit !(value > 0) }'; then
  echo "DEPTH_TOL_CM must be greater than zero." >&2
  exit 2
fi
if ! awk -v value="$K_TOL_PCM" 'BEGIN { exit !(value >= 0) }'; then
  echo "K_TOL_PCM cannot be negative." >&2
  exit 2
fi
if ! is_positive_integer "$REPLICATES_PER_DEPTH" || ! is_nonnegative_integer "$MAX_REPLICATE_FAILURES" || \
   ! is_positive_integer "$MAX_ADDITIONAL_ITERATIONS" || ! is_positive_integer "$PARTICLES" || \
   ! is_nonnegative_integer "$INACTIVE" || ! is_positive_integer "$ACTIVE" || (( ACTIVE < 2 )); then
  echo "Invalid replicate/iteration/particle/cycle setting." >&2
  exit 2
fi
if ! [[ "$THREADS_PER_BLOCK" =~ ^(32|64|128|256|512|1024)$ ]]; then
  echo "THREADS_PER_BLOCK must be one of 32, 64, 128, 256, 512, or 1024." >&2
  exit 2
fi
if ! [[ "$NO_BUILD" =~ ^[01]$ && "$DEVICE_LTO" =~ ^[01]$ && "$FAST_MATH" =~ ^[01]$ && "$STOP_ON_CONFIDENCE_INTERVAL" =~ ^[01]$ ]]; then
  echo "NO_BUILD, DEVICE_LTO, FAST_MATH, and STOP_ON_CONFIDENCE_INTERVAL must be 0 or 1." >&2
  exit 2
fi

if (( NO_BUILD == 0 )); then
  if ! command -v "$NVCC" >/dev/null 2>&1; then
    echo "nvcc not found: $NVCC" >&2
    echo "Load the cluster CUDA module or set NVCC to the full nvcc path." >&2
    exit 1
  fi
  if ! [[ "$CUDA_ARCH" =~ ^(native|sm_[0-9]+[a-z]?|compute_[0-9]+[a-z]?)$ ]]; then
    echo "CUDA_ARCH must be native, sm_NN, or compute_NN; got: $CUDA_ARCH" >&2
    exit 2
  fi
  echo "Building $EXE"
  if ! make TARGET="$EXE" NVCC="$NVCC" CUDA_ARCH="$CUDA_ARCH" \
      THREADS_PER_BLOCK="$THREADS_PER_BLOCK" DEVICE_LTO="$DEVICE_LTO" \
      FAST_MATH="$FAST_MATH" "$EXE"; then
    echo "Build failed." >&2
    exit 1
  fi
fi
if [[ ! -x "$EXE" ]]; then
  echo "Executable not found or not executable: $EXE" >&2
  exit 1
fi

mkdir -p "$LOG_DIR"
source_real="$(cd "$SOURCE_DIR" && pwd -P)"
output_real="$(cd "$LOG_DIR" && pwd -P)"
if [[ "$source_real" == "$output_real" ]]; then
  echo "LOG_DIR must differ from RESUME_DIR; refusing to overwrite the old run." >&2
  exit 2
fi
cp "$SOURCE_SUMMARY" "$SUMMARY"
cp "$SOURCE_REPLICATES" "$REPLICATE_SUMMARY"

read -r -a case_args <<< "$CASE_ARGS"

evaluation="$(awk -F, 'NR>1 && $1 ~ /^[0-9]+$/ && $1>maximum { maximum=$1 } END { print maximum+0 }' "$SOURCE_SUMMARY")"
starting_evaluation="$evaluation"
low_depth="$(result_value final_lower_depth_cm)"
low_k="$(result_value final_lower_k)"
high_depth="$(result_value final_upper_depth_cm)"
high_k="$(result_value final_upper_k)"

for value in "$low_depth" "$low_k" "$high_depth" "$high_k"; do
  if ! is_number "$value"; then
    echo "Could not read a valid final bracket from $SOURCE_RESULT" >&2
    exit 2
  fi
done
if ! awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { exit !(lo <= hi) }'; then
  echo "Invalid resumed bracket: $low_depth to $high_depth" >&2
  exit 2
fi

best_line="$(awk -F, -v target="$TARGET_K" 'NR>1 && $10 == "OK" {
  d=$4-target; if (d<0) d=-d
  if (!found || d<best) {
    found=1; best=d
    line=$2 " " $3 " " $4 " " $5 " " $6 " " $7 " " $9
  }
} END { print line }' "$SOURCE_SUMMARY")"
read -r best_depth best_actual best_k best_stderr best_within_stderr best_between_stderr best_scatter_ratio <<< "$best_line"
best_delta="$(awk -v k="$best_k" -v target="$TARGET_K" 'BEGIN { d=k-target; if (d<0) d=-d; printf "%.15g", d }')"

SINGLE_VALID=0
SINGLE_DEPTH=""
SINGLE_ACTUAL=""
SINGLE_K=""
SINGLE_STDERR=""
SINGLE_LOG=""
SINGLE_SEED=""

RUN_VALID=0
RUN_DEPTH=""
RUN_ACTUAL=""
RUN_K=""
RUN_STDERR=""
RUN_WITHIN_STDERR=""
RUN_BETWEEN_STDERR=""
RUN_SCATTER_RATIO=""
RUN_REPLICATE_COUNT=0
RUN_LOGS=""
RUN_SEEDS=""

format_depth() { awk -v value="$1" 'BEGIN { printf "%.8f", value }'; }
format_k() { awk -v value="$1" 'BEGIN { printf "%.7f", value }'; }
format_general() { awk -v value="$1" 'BEGIN { printf "%.12g", value }'; }

independent_seed() {
  local depth_key
  depth_key="$(awk -v depth="$1" 'BEGIN { printf "%.0f", depth * 1000000.0 }')"
  echo $(( (BASE_SEED % 281474976710656 + (depth_key * 1000003) % 281474976710656 + ($2 * 104729) % 281474976710656) % 281474976710656 ))
}

invoke_single_rod_case() {
  local requested="$1" eval_no="$2" attempt="$3" seed="$4"
  local depth depth_tag status status_text
  depth="$(format_depth "$requested")"
  depth_tag="${depth//./p}"
  SINGLE_LOG="$LOG_DIR/eval_$(printf '%02d' "$eval_no")_attempt_$(printf '%02d' "$attempt")_seed_${seed}_rod_${depth_tag}cm.log"

  echo "    replicate attempt $attempt, seed $seed"
  "$EXE" "${case_args[@]}" --rod-depth="$depth" --particles="$PARTICLES" \
    --inactive="$INACTIVE" --active="$ACTIVE" --seed="$seed" > "$SINGLE_LOG" 2>&1
  status=$?

  SINGLE_DEPTH="$depth"
  SINGLE_ACTUAL="$(awk -F': ' '/Actual (modeled|quantized) rod depth:/ { gsub(/ cm.*/, "", $2); value=$2 } END { print value }' "$SINGLE_LOG")"
  SINGLE_K="$(awk -F': ' '/^k_mean:/ { value=$2 } END { print value }' "$SINGLE_LOG")"
  SINGLE_STDERR="$(awk -F': ' '/^k_stderr\(mean\):/ { value=$2 } END { print value }' "$SINGLE_LOG")"
  SINGLE_SEED="$seed"

  if (( status == 0 )) && is_number "$SINGLE_K"; then
    SINGLE_VALID=1
    status_text="OK"
    echo "        k = $SINGLE_K, stderr = ${SINGLE_STDERR:-unknown}"
  else
    SINGLE_VALID=0
    status_text="FAILED:$status"
    echo "WARNING: depth $depth seed $seed failed with exit code $status; keeping $SINGLE_LOG." >&2
  fi
  echo "$eval_no,$attempt,$seed,$depth,$SINGLE_ACTUAL,$SINGLE_K,$SINGLE_STDERR,$status_text,$SINGLE_LOG" >> "$REPLICATE_SUMMARY"
}

update_best_run() {
  local delta
  delta="$(awk -v k="$RUN_K" -v target="$TARGET_K" 'BEGIN { d=k-target; if (d<0) d=-d; printf "%.15g", d }')"
  if [[ -z "$best_delta" ]] || awk -v candidate="$delta" -v best="$best_delta" 'BEGIN { exit !(candidate < best) }'; then
    best_depth="$RUN_DEPTH"
    best_actual="$RUN_ACTUAL"
    best_k="$RUN_K"
    best_stderr="$RUN_STDERR"
    best_within_stderr="$RUN_WITHIN_STDERR"
    best_between_stderr="$RUN_BETWEEN_STDERR"
    best_scatter_ratio="$RUN_SCATTER_RATIO"
    best_delta="$delta"
  fi
}

run_case() {
  local requested="$1"
  local depth attempt max_attempts seed count stats
  local -a k_values=() stderr_values=() actual_values=() seeds=() logs=()

  evaluation=$((evaluation + 1))
  depth="$(format_depth "$requested")"
  echo "==> evaluation $evaluation: rod depth $depth cm, replicates $REPLICATES_PER_DEPTH"

  max_attempts=$((REPLICATES_PER_DEPTH + MAX_REPLICATE_FAILURES))
  for (( attempt=1; attempt<=max_attempts && ${#k_values[@]}<REPLICATES_PER_DEPTH; ++attempt )); do
    seed="$(independent_seed "$depth" "$attempt")"
    if [[ "$seed" == "0" ]]; then seed=1; fi
    invoke_single_rod_case "$depth" "$evaluation" "$attempt" "$seed"
    if (( SINGLE_VALID == 1 )); then
      k_values+=("$SINGLE_K")
      stderr_values+=("${SINGLE_STDERR:-nan}")
      actual_values+=("$SINGLE_ACTUAL")
      seeds+=("$SINGLE_SEED")
      logs+=("$SINGLE_LOG")
    fi
  done

  count=${#k_values[@]}
  RUN_DEPTH="$depth"
  RUN_REPLICATE_COUNT=$count
  RUN_SEEDS="$(IFS=';'; echo "${seeds[*]-}")"
  RUN_LOGS="$(IFS=';'; echo "${logs[*]-}")"
  if (( count < REPLICATES_PER_DEPTH )); then
    RUN_VALID=0
    RUN_ACTUAL=""
    RUN_K="nan"
    RUN_STDERR="nan"
    RUN_WITHIN_STDERR="nan"
    RUN_BETWEEN_STDERR="nan"
    RUN_SCATTER_RATIO="nan"
    echo "$evaluation,$depth,,NaN,NaN,NaN,NaN,$count,NaN,FAILED,$RUN_SEEDS" >> "$SUMMARY"
    return 1
  fi

  stats="$(awk -v ks="${k_values[*]}" -v ses="${stderr_values[*]}" 'BEGIN {
    n=split(ks, kval, " "); split(ses, serr, " "); mean=0.0
    for (i=1; i<=n; ++i) mean+=kval[i]; mean/=n
    ss=0.0; se2=0.0; finite=0
    for (i=1; i<=n; ++i) {
      d=kval[i]-mean; ss+=d*d
      if (serr[i] ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
        se2+=serr[i]*serr[i]; finite++
      }
    }
    if (n == 1) {
      within=(finite>0) ? sqrt(se2) : -1.0
      printf "%.17g %.17g %.17g -1 -1", mean, within, within; exit
    }
    sample_sd=sqrt(ss/(n-1)); between=sample_sd/sqrt(n)
    if (finite>0) {
      within=sqrt(se2)/n; rms=sqrt(se2/finite)
      scatter=(rms>0.0) ? sample_sd/rms : 0.0
      conservative=(within>between) ? within : between
    } else { within=-1.0; scatter=-1.0; conservative=between }
    printf "%.17g %.17g %.17g %.17g %.17g", mean, conservative, within, between, scatter
  }')"
  read -r RUN_K RUN_STDERR RUN_WITHIN_STDERR RUN_BETWEEN_STDERR RUN_SCATTER_RATIO <<< "$stats"
  RUN_ACTUAL="${actual_values[0]}"
  RUN_VALID=1
  echo "$evaluation,$depth,$RUN_ACTUAL,$(format_k "$RUN_K"),$(format_general "$RUN_STDERR"),$(format_general "$RUN_WITHIN_STDERR"),$(format_general "$RUN_BETWEEN_STDERR"),$count,$(awk -v value="$RUN_SCATTER_RATIO" 'BEGIN { printf "%.3f", value }'),OK,$RUN_SEEDS" >> "$SUMMARY"
  echo "    k = $(format_k "$RUN_K") +/- $(format_general "$RUN_STDERR")"
  if awk -v ratio="$RUN_SCATTER_RATIO" -v warning="$SCATTER_WARNING_RATIO" 'BEGIN { exit !(ratio >= 0 && ratio > warning) }'; then
    echo "WARNING: replicate scatter is high: $(format_general "$RUN_SCATTER_RATIO")" >&2
  fi
  update_best_run
  return 0
}

brackets_target() {
  awk -v first="$1" -v second="$2" -v target="$TARGET_K" 'BEGIN {
    a=first-target; b=second-target
    exit !((a<=0 && b>=0) || (a>=0 && b<=0))
  }'
}

if ! brackets_target "$low_k" "$high_k"; then
  echo "The saved interval does not bracket target k=$TARGET_K." >&2
  echo "lower: $low_depth cm, k=$low_k; upper: $high_depth cm, k=$high_k" >&2
  exit 4
fi

echo "CriticalSearch_Binary continuation"
echo "  source: $SOURCE_DIR"
echo "  existing evaluations: 1-$starting_evaluation; next evaluation: $((starting_evaluation + 1))"
echo "  resumed bracket: $low_depth cm (k=$(format_k "$low_k")) to $high_depth cm (k=$(format_k "$high_k"))"
echo "  target k: $TARGET_K; depth tolerance: $DEPTH_TOL_CM cm; k tolerance: $K_TOL_PCM pcm"
echo "  particles/inactive/active: $PARTICLES/$INACTIVE/$ACTIVE"
echo "  base seed: $BASE_SEED; confidence sigma: $CONFIDENCE_SIGMA; confidence early stop: $STOP_ON_CONFIDENCE_INTERVAL"
echo "  output: $LOG_DIR"

stop_reason="maximum additional iterations reached"
search_status="COMPLETE"
additional_iterations=0
for (( iteration=1; iteration<=MAX_ADDITIONAL_ITERATIONS; ++iteration )); do
  width="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.15g", hi-lo }')"
  if awk -v width="$width" -v tolerance="$DEPTH_TOL_CM" 'BEGIN { exit !(width <= tolerance) }'; then
    stop_reason="depth bracket tolerance reached"
    break
  fi
  additional_iterations=$iteration

  midpoint="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.12g", (lo+hi)/2.0 }')"
  valid=0
  if run_case "$midpoint"; then
    valid=1
  else
    left_candidate="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.12g", lo+0.45*(hi-lo) }')"
    right_candidate="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.12g", lo+0.55*(hi-lo) }')"
    echo "WARNING: midpoint failed; trying $left_candidate cm and then $right_candidate cm." >&2
    if run_case "$left_candidate" || run_case "$right_candidate"; then valid=1; fi
  fi
  if (( valid == 0 )); then
    stop_reason="three candidate simulations failed in the current bracket"
    search_status="PARTIAL"
    break
  fi

  mid_depth="$RUN_DEPTH"
  mid_k="$RUN_K"
  if awk -v k="$mid_k" -v target="$TARGET_K" -v sigma="$CONFIDENCE_SIGMA" -v se="$RUN_STDERR" 'BEGIN {
    d=k-target; if (d<0) d=-d; exit !(d <= sigma*se)
  }'; then
    if (( STOP_ON_CONFIDENCE_INTERVAL == 1 )); then
      stop_reason="target k lies inside the averaged confidence interval"
      low_depth="$mid_depth"; high_depth="$mid_depth"; low_k="$mid_k"; high_k="$mid_k"
      break
    fi
    echo "    target lies inside +/- $CONFIDENCE_SIGMA stderr, but confidence-interval early stop is disabled; continuing."
  fi

  delta_pcm="$(awk -v k="$mid_k" -v target="$TARGET_K" 'BEGIN { d=k-target; if (d<0) d=-d; printf "%.8f", d*100000.0 }')"
  if awk -v delta="$delta_pcm" -v tolerance="$K_TOL_PCM" 'BEGIN { exit !(tolerance>0 && delta<=tolerance) }'; then
    stop_reason="k tolerance reached"
    low_depth="$mid_depth"; high_depth="$mid_depth"; low_k="$mid_k"; high_k="$mid_k"
    break
  fi

  if brackets_target "$low_k" "$mid_k"; then
    high_depth="$mid_depth"; high_k="$mid_k"
  else
    low_depth="$mid_depth"; low_k="$mid_k"
  fi
  echo "    bracket is now $low_depth cm (k=$(format_k "$low_k")) to $high_depth cm (k=$(format_k "$high_k"))"
done

# The original loop only checked bracket width at the beginning of an iteration.
# Check once more after the last new point so an adequate final bracket is
# reported as tolerance convergence rather than as a safety-cap termination.
final_width="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.15g", hi-lo }')"
if [[ "$search_status" == "COMPLETE" ]] && \
   awk -v width="$final_width" -v tolerance="$DEPTH_TOL_CM" 'BEGIN { exit !(width <= tolerance) }'; then
  stop_reason="depth bracket tolerance reached"
fi

if awk -v low="$low_k" -v high="$high_k" 'BEGIN { exit !(low != high) }'; then
  estimated_depth="$(awk -v lo="$low_depth" -v hi="$high_depth" -v kl="$low_k" -v kh="$high_k" -v target="$TARGET_K" 'BEGIN {
    estimate=lo+(target-kl)*(hi-lo)/(kh-kl)
    if (estimate<lo) estimate=lo; if (estimate>hi) estimate=hi
    printf "%.8f", estimate
  }')"
else
  estimated_depth="$low_depth"
fi

aro_check_status="$(result_value aro_check_status)"
aro_check_depth="$(result_value aro_check_depth_cm)"
aro_expected_k="$(result_value aro_expected_k)"
aro_observed_k="$(result_value aro_observed_k)"
aro_stderr="$(result_value aro_k_stderr)"
aro_difference="$(result_value aro_difference_pcm)"
aro_allowed="$(result_value aro_allowed_pcm)"

{
  echo "STATUS=$search_status"
  echo "stop_reason=$stop_reason"
  echo "resumed_from=$RESUME_DIR"
  echo "starting_evaluation=$starting_evaluation"
  echo "last_evaluation=$evaluation"
  echo "additional_bisection_iterations=$additional_iterations"
  echo "target_k=$(format_k "$TARGET_K")"
  echo "depth_tolerance_cm=$DEPTH_TOL_CM"
  echo "k_tolerance_pcm=$K_TOL_PCM"
  echo "particles=$PARTICLES"
  echo "inactive_cycles=$INACTIVE"
  echo "active_cycles=$ACTIVE"
  echo "replicates_per_depth=$REPLICATES_PER_DEPTH"
  echo "base_seed=$BASE_SEED"
  echo "confidence_sigma=$CONFIDENCE_SIGMA"
  echo "stop_on_confidence_interval=$STOP_ON_CONFIDENCE_INTERVAL"
  echo "aro_check_status=$aro_check_status"
  echo "aro_check_depth_cm=$aro_check_depth"
  echo "aro_expected_k=$aro_expected_k"
  echo "aro_observed_k=$aro_observed_k"
  echo "aro_k_stderr=$aro_stderr"
  echo "aro_difference_pcm=$aro_difference"
  echo "aro_allowed_pcm=$aro_allowed"
  echo "final_lower_depth_cm=$low_depth final_lower_k=$(format_k "$low_k")"
  echo "final_upper_depth_cm=$high_depth final_upper_k=$(format_k "$high_k")"
  echo "final_bracket_width_cm=$final_width"
  echo "estimated_critical_depth_cm=$estimated_depth"
  echo "estimated_fuel_insertion_below_top_moderator_cm=$(awk -v d="$estimated_depth" 'BEGIN { printf "%.8f", d-21.42 }')"
  echo "best_sampled_depth_cm=$best_depth"
  echo "best_sampled_actual_depth_cm=$best_actual"
  echo "best_sampled_k=$(format_k "$best_k")"
  echo "best_sampled_k_stderr=$(format_general "$best_stderr")"
  echo "best_sampled_within_stderr=$(format_general "$best_within_stderr")"
  echo "best_sampled_between_stderr=$(format_general "$best_between_stderr")"
  echo "best_sampled_scatter_ratio=$(awk -v value="$best_scatter_ratio" 'BEGIN { printf "%.3f", value }')"
  echo "summary=$SUMMARY"
  echo "replicates=$REPLICATE_SUMMARY"
} | tee "$RESULT"
