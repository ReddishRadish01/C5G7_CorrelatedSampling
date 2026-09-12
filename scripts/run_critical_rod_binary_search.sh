#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

EXE="${EXE:-build/c5g7}"
NVCC="${NVCC:-nvcc}"
LOWER_CM="${LOWER_CM:-82.0}"
UPPER_CM="${UPPER_CM:-84.0}"
TARGET_K="${TARGET_K:-1.0}"
DEPTH_TOL_CM="${DEPTH_TOL_CM:-0.01}"
K_TOL_PCM="${K_TOL_PCM:-0.0}"
ARO_CHECK_DEPTH_CM="${ARO_CHECK_DEPTH_CM:-21.42}"
ARO_EXPECTED_K="${ARO_EXPECTED_K:-1.14390}"
ARO_MIN_TOLERANCE_PCM="${ARO_MIN_TOLERANCE_PCM:-20.0}"
REPLICATES_PER_DEPTH="${REPLICATES_PER_DEPTH:-1}"
MAX_REPLICATE_FAILURES="${MAX_REPLICATE_FAILURES:-2}"
BASE_SEED="${BASE_SEED:-92235922383}"
CONFIDENCE_SIGMA="${CONFIDENCE_SIGMA:-2.0}"
STOP_ON_CONFIDENCE_INTERVAL="${STOP_ON_CONFIDENCE_INTERVAL:-0}"
SCATTER_WARNING_RATIO="${SCATTER_WARNING_RATIO:-2.0}"
MAX_ITER="${MAX_ITER:-10}"
PARTICLES="${PARTICLES:-20000000}"
INACTIVE="${INACTIVE:-150}"
ACTIVE="${ACTIVE:-200}"
CUDA_ARCH="${CUDA_ARCH:-native}"
THREADS_PER_BLOCK="${THREADS_PER_BLOCK:-32}"
DEVICE_LTO="${DEVICE_LTO:-1}"
FAST_MATH="${FAST_MATH:-0}"
NO_BUILD="${NO_BUILD:-0}"
SKIP_ARO_CHECK="${SKIP_ARO_CHECK:-0}"
CASE_ARGS="${CASE_ARGS:---critical-search-binary}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-results/critical_rod_binary_search_${RUN_TAG}}"
SUMMARY="${SUMMARY:-$LOG_DIR/summary.csv}"
REPLICATE_SUMMARY="${REPLICATE_SUMMARY:-$LOG_DIR/replicates.csv}"
RESULT="${RESULT:-$LOG_DIR/result.txt}"

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

for value in \
  "$LOWER_CM" "$UPPER_CM" "$TARGET_K" "$DEPTH_TOL_CM" "$K_TOL_PCM" \
  "$ARO_CHECK_DEPTH_CM" "$ARO_EXPECTED_K" "$ARO_MIN_TOLERANCE_PCM" \
  "$CONFIDENCE_SIGMA" "$SCATTER_WARNING_RATIO"; do
  if ! is_number "$value"; then
    echo "Invalid numeric setting: $value" >&2
    exit 2
  fi
done

if ! awk -v lo="$LOWER_CM" -v hi="$UPPER_CM" 'BEGIN { exit !(lo < hi) }'; then
  echo "LOWER_CM must be less than UPPER_CM." >&2
  exit 2
fi
if ! awk -v value="$DEPTH_TOL_CM" 'BEGIN { exit !(value > 0) }'; then
  echo "DEPTH_TOL_CM must be greater than zero." >&2
  exit 2
fi
if ! awk -v value="$K_TOL_PCM" 'BEGIN { exit !(value >= 0) }'; then
  echo "K_TOL_PCM cannot be negative." >&2
  exit 2
fi
if ! awk -v d="$ARO_CHECK_DEPTH_CM" -v k="$ARO_EXPECTED_K" -v t="$ARO_MIN_TOLERANCE_PCM" \
  'BEGIN { exit !(d >= 0 && k > 0 && t >= 0) }'; then
  echo "ARO check depth/tolerance cannot be negative, and ARO_EXPECTED_K must be positive." >&2
  exit 2
fi
if ! is_positive_integer "$REPLICATES_PER_DEPTH"; then
  echo "REPLICATES_PER_DEPTH must be a positive integer." >&2
  exit 2
fi
if ! is_nonnegative_integer "$MAX_REPLICATE_FAILURES"; then
  echo "MAX_REPLICATE_FAILURES must be a nonnegative integer." >&2
  exit 2
fi
if ! is_positive_integer "$MAX_ITER" || ! is_positive_integer "$PARTICLES" || \
  ! is_nonnegative_integer "$INACTIVE" || ! is_positive_integer "$ACTIVE" || (( ACTIVE < 2 )); then
  echo "MAX_ITER/PARTICLES must be positive, INACTIVE nonnegative, and ACTIVE at least 2." >&2
  exit 2
fi
if ! [[ "$THREADS_PER_BLOCK" =~ ^(32|64|128|256|512|1024)$ ]]; then
  echo "THREADS_PER_BLOCK must be one of 32, 64, 128, 256, 512, or 1024." >&2
  exit 2
fi
if ! [[ "$NO_BUILD" =~ ^[01]$ && "$SKIP_ARO_CHECK" =~ ^[01]$ && "$DEVICE_LTO" =~ ^[01]$ && "$FAST_MATH" =~ ^[01]$ && "$STOP_ON_CONFIDENCE_INTERVAL" =~ ^[01]$ ]]; then
  echo "NO_BUILD, SKIP_ARO_CHECK, DEVICE_LTO, FAST_MATH, and STOP_ON_CONFIDENCE_INTERVAL must be 0 or 1." >&2
  exit 2
fi
if ! awk -v c="$CONFIDENCE_SIGMA" -v s="$SCATTER_WARNING_RATIO" 'BEGIN { exit !(c > 0 && s > 0) }'; then
  echo "CONFIDENCE_SIGMA and SCATTER_WARNING_RATIO must be greater than zero." >&2
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
  echo "  nvcc: $(command -v "$NVCC")"
  echo "  CUDA_ARCH: $CUDA_ARCH"
  "$NVCC" --version | tail -n 1
  if ! make \
    TARGET="$EXE" \
    NVCC="$NVCC" \
    CUDA_ARCH="$CUDA_ARCH" \
    THREADS_PER_BLOCK="$THREADS_PER_BLOCK" \
    DEVICE_LTO="$DEVICE_LTO" \
    FAST_MATH="$FAST_MATH" \
    "$EXE"; then
    echo "Build failed." >&2
    exit 1
  fi
fi

if [[ ! -x "$EXE" ]]; then
  echo "Executable not found or not executable: $EXE" >&2
  echo "Run with NO_BUILD=0, or set EXE to an existing Linux binary." >&2
  exit 1
fi

mkdir -p "$LOG_DIR"
echo "evaluation,requested_depth_cm,actual_depth_cm,k_mean,k_stderr_conservative,k_stderr_within,k_stderr_between,replicates,scatter_ratio,status,seeds" > "$SUMMARY"
echo "evaluation,attempt,seed,requested_depth_cm,actual_depth_cm,k_mean,k_stderr_mean,status,log_file" > "$REPLICATE_SUMMARY"

read -r -a case_args <<< "$CASE_ARGS"

evaluation=0
best_depth=""
best_actual=""
best_k=""
best_stderr=""
best_within_stderr=""
best_between_stderr=""
best_scatter_ratio=""
best_delta=""

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

format_depth() {
  awk -v value="$1" 'BEGIN { printf "%.8f", value }'
}

format_k() {
  awk -v value="$1" 'BEGIN { printf "%.7f", value }'
}

format_general() {
  awk -v value="$1" 'BEGIN { printf "%.12g", value }'
}

independent_seed() {
  local depth_key seed
  depth_key="$(awk -v depth="$1" 'BEGIN { printf "%.0f", depth * 1000000.0 }')"
  seed=$(( (BASE_SEED % 281474976710656 + (depth_key * 1000003) % 281474976710656 + ($2 * 104729) % 281474976710656) % 281474976710656 ))
  if (( seed % 2 == 0 )); then seed=$((seed + 1)); fi
  echo "$seed"
}

invoke_single_rod_case() {
  local requested="$1"
  local eval_no="$2"
  local attempt="$3"
  local seed="$4"
  local depth depth_tag status status_text

  depth="$(format_depth "$requested")"
  depth_tag="${depth//./p}"
  SINGLE_LOG="$LOG_DIR/eval_$(printf '%02d' "$eval_no")_attempt_$(printf '%02d' "$attempt")_seed_${seed}_rod_${depth_tag}cm.log"

  echo "    replicate attempt $attempt, seed $seed"
  "$EXE" \
    "${case_args[@]}" \
    --rod-depth="$depth" \
    --particles="$PARTICLES" \
    --inactive="$INACTIVE" \
    --active="$ACTIVE" \
    --seed="$seed" \
    > "$SINGLE_LOG" 2>&1
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
  local -a k_values=()
  local -a stderr_values=()
  local -a actual_values=()
  local -a seeds=()
  local -a logs=()

  evaluation=$((evaluation + 1))
  depth="$(format_depth "$requested")"
  if (( REPLICATES_PER_DEPTH == 1 )); then
    echo "==> evaluation $evaluation: rod depth $depth cm, running once"
  else
    echo "==> evaluation $evaluation: rod depth $depth cm, averaging $REPLICATES_PER_DEPTH independent runs"
  fi

  max_attempts=$((REPLICATES_PER_DEPTH + MAX_REPLICATE_FAILURES))
  for (( attempt=1; attempt<=max_attempts && ${#k_values[@]}<REPLICATES_PER_DEPTH; ++attempt )); do
    seed="$(independent_seed "$depth" "$attempt")"
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
    n=split(ks, kval, " ");
    split(ses, serr, " ");
    mean=0.0;
    for (i=1; i<=n; ++i) mean+=kval[i];
    mean/=n;
    ss=0.0;
    se2=0.0;
    finite=0;
    for (i=1; i<=n; ++i) {
      d=kval[i]-mean;
      ss+=d*d;
      if (serr[i] ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/) {
        se2+=serr[i]*serr[i];
        finite++;
      }
    }
    if (n == 1) {
      within=(finite>0) ? sqrt(se2) : -1.0;
      conservative=within;
      between=-1.0;
      scatter=-1.0;
      printf "%.17g %.17g %.17g %.17g %.17g", mean, conservative, within, between, scatter;
      exit;
    }
    sample_sd=sqrt(ss/(n-1));
    between=sample_sd/sqrt(n);
    if (finite>0) {
      within=sqrt(se2)/n;
      rms=sqrt(se2/finite);
      scatter=(rms>0.0) ? sample_sd/rms : 0.0;
      conservative=(within>between) ? within : between;
    } else {
      within=-1.0;
      scatter=-1.0;
      conservative=between;
    }
    printf "%.17g %.17g %.17g %.17g %.17g", mean, conservative, within, between, scatter;
  }')"
  read -r RUN_K RUN_STDERR RUN_WITHIN_STDERR RUN_BETWEEN_STDERR RUN_SCATTER_RATIO <<< "$stats"
  RUN_ACTUAL="${actual_values[0]}"
  RUN_VALID=1

  echo "$evaluation,$depth,$RUN_ACTUAL,$(format_k "$RUN_K"),$(format_general "$RUN_STDERR"),$(format_general "$RUN_WITHIN_STDERR"),$(format_general "$RUN_BETWEEN_STDERR"),$count,$(awk -v value="$RUN_SCATTER_RATIO" 'BEGIN { printf "%.3f", value }'),OK,$RUN_SEEDS" >> "$SUMMARY"

  if (( count == 1 )); then
    echo "    k = $(format_k "$RUN_K") +/- $(format_general "$RUN_STDERR")"
    echo "    within-run SE = $(format_general "$RUN_WITHIN_STDERR"); between-run SE unavailable (single run)"
  else
    echo "    averaged k = $(format_k "$RUN_K") +/- $(format_general "$RUN_STDERR")"
    echo "    within-run SE = $(format_general "$RUN_WITHIN_STDERR"); between-run SE = $(format_general "$RUN_BETWEEN_STDERR"); scatter ratio = $(awk -v value="$RUN_SCATTER_RATIO" 'BEGIN { printf "%.2f", value }')"
  fi
  if awk -v ratio="$RUN_SCATTER_RATIO" -v warning="$SCATTER_WARNING_RATIO" 'BEGIN { exit !(ratio >= 0 && ratio > warning) }'; then
    echo "WARNING: replicate scatter is $(awk -v value="$RUN_SCATTER_RATIO" 'BEGIN { printf "%.2f", value }')x the runs' reported stderr. This may indicate GPU nondeterminism or underestimated uncertainty." >&2
  fi

  update_best_run
  return 0
}

run_endpoint() {
  local requested="$1"
  local side="$2"
  local fallback

  if run_case "$requested"; then
    return 0
  fi

  if [[ "$side" == "lower" ]]; then
    fallback="$(awk -v lo="$LOWER_CM" -v hi="$UPPER_CM" 'BEGIN { printf "%.12g", lo + 0.05*(hi-lo) }')"
  else
    fallback="$(awk -v lo="$LOWER_CM" -v hi="$UPPER_CM" 'BEGIN { printf "%.12g", hi - 0.05*(hi-lo) }')"
  fi
  echo "WARNING: endpoint failed; trying nearby depth $fallback cm." >&2
  run_case "$fallback"
}

brackets_target() {
  awk -v first="$1" -v second="$2" -v target="$TARGET_K" 'BEGIN {
    a=first-target;
    b=second-target;
    exit !((a<=0 && b>=0) || (a>=0 && b<=0));
  }'
}

echo "CriticalSearch_Binary rod-depth search (Linux/WSL direct NVCC)"
echo "  exe: $EXE"
echo "  initial bracket: $LOWER_CM to $UPPER_CM cm"
echo "  fuel insertion below top moderator: $(awk -v d="$LOWER_CM" 'BEGIN {printf "%.2f", d-21.42}') to $(awk -v d="$UPPER_CM" 'BEGIN {printf "%.2f", d-21.42}') cm"
echo "  target k: $TARGET_K"
echo "  depth tolerance: $DEPTH_TOL_CM cm; k tolerance: $K_TOL_PCM pcm"
echo "  particles/inactive/active: $PARTICLES/$INACTIVE/$ACTIVE"
echo "  independent replicates/depth: $REPLICATES_PER_DEPTH; base seed: $BASE_SEED"
echo "  statistical decision interval: +/- $CONFIDENCE_SIGMA aggregate stderr"
echo "  stop when target lies inside confidence interval: $STOP_ON_CONFIDENCE_INTERVAL"
echo "  output: $LOG_DIR"

aro_check_status="SKIPPED"
aro_check_k=""
aro_check_stderr=""
aro_check_delta_pcm=""
aro_check_limit_pcm=""

if (( SKIP_ARO_CHECK == 0 )); then
  echo
  echo "Running ARO equivalence check at $ARO_CHECK_DEPTH_CM cm..."
  if ! run_case "$ARO_CHECK_DEPTH_CM"; then
    {
      echo "STATUS=ARO_CHECK_FAILED"
      echo "reason=simulation failed"
      echo "aro_check_depth_cm=$ARO_CHECK_DEPTH_CM"
      echo "aro_expected_k=$ARO_EXPECTED_K"
      echo "logs=$RUN_LOGS"
    } | tee "$RESULT"
    exit 5
  fi

  aro_check_k="$(format_k "$RUN_K")"
  aro_check_stderr="$(format_general "$RUN_STDERR")"
  aro_check_delta_pcm="$(awk -v k="$RUN_K" -v expected="$ARO_EXPECTED_K" 'BEGIN {
    d=(k-1.0)/k-(expected-1.0)/expected;
    if (d<0) d=-d;
    printf "%.3f", d*100000.0;
  }')"
  aro_sigma_pcm="$(awk -v se="$RUN_STDERR" -v k="$RUN_K" 'BEGIN { printf "%.12g", se/(k*k)*100000.0 }')"
  aro_check_limit_pcm="$(awk -v minimum="$ARO_MIN_TOLERANCE_PCM" -v sigma="$aro_sigma_pcm" 'BEGIN {
    limit=3.0*sigma;
    if (limit<minimum) limit=minimum;
    printf "%.3f", limit;
  }')"

  echo "    expected k = $ARO_EXPECTED_K"
  echo "    observed k = $aro_check_k +/- $aro_check_stderr"
  echo "    difference = $aro_check_delta_pcm pcm; allowed = $aro_check_limit_pcm pcm"

  if awk -v delta="$aro_check_delta_pcm" -v limit="$aro_check_limit_pcm" 'BEGIN { exit !(delta > limit) }'; then
    {
      echo "STATUS=ARO_CHECK_FAILED"
      echo "reason=multiplication factor mismatch"
      echo "aro_check_depth_cm=$ARO_CHECK_DEPTH_CM"
      echo "aro_expected_k=$ARO_EXPECTED_K"
      echo "aro_observed_k=$aro_check_k"
      echo "aro_k_stderr=$aro_check_stderr"
      echo "aro_difference_pcm=$aro_check_delta_pcm"
      echo "aro_allowed_pcm=$aro_check_limit_pcm"
      echo "logs=$RUN_LOGS"
    } | tee "$RESULT"
    exit 5
  fi

  aro_check_status="PASSED"
  echo "    ARO equivalence check passed."
  echo
fi

if ! run_endpoint "$LOWER_CM" lower; then
  printf 'STATUS=FAILED\nCould not obtain a valid lower-endpoint result.\n' | tee "$RESULT" >&2
  exit 3
fi
low_depth="$RUN_DEPTH"
low_k="$RUN_K"

if ! run_endpoint "$UPPER_CM" upper; then
  printf 'STATUS=FAILED\nCould not obtain a valid upper-endpoint result.\n' | tee "$RESULT" >&2
  exit 3
fi
high_depth="$RUN_DEPTH"
high_k="$RUN_K"

if ! brackets_target "$low_k" "$high_k"; then
  {
    echo "STATUS=NOT_BRACKETED"
    echo "The initial interval does not bracket target k=$TARGET_K."
    echo "lower_depth_cm=$low_depth lower_k=$(format_k "$low_k")"
    echo "upper_depth_cm=$high_depth upper_k=$(format_k "$high_k")"
    echo "Expand LOWER_CM/UPPER_CM and run again."
  } | tee "$RESULT" >&2
  exit 4
fi

stop_reason="maximum iterations reached"
search_status="COMPLETE"
for (( iteration=1; iteration<=MAX_ITER; ++iteration )); do
  width="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.15g", hi-lo }')"
  if awk -v width="$width" -v tolerance="$DEPTH_TOL_CM" 'BEGIN { exit !(width <= tolerance) }'; then
    stop_reason="depth bracket tolerance reached"
    break
  fi

  midpoint="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.12g", (lo+hi)/2.0 }')"
  valid=0
  if run_case "$midpoint"; then
    valid=1
  else
    left_candidate="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.12g", lo+0.45*(hi-lo) }')"
    right_candidate="$(awk -v lo="$low_depth" -v hi="$high_depth" 'BEGIN { printf "%.12g", lo+0.55*(hi-lo) }')"
    echo "WARNING: midpoint failed; trying $left_candidate cm and then $right_candidate cm." >&2
    if run_case "$left_candidate" || run_case "$right_candidate"; then
      valid=1
    fi
  fi

  if (( valid == 0 )); then
    stop_reason="three candidate simulations failed in the current bracket"
    search_status="PARTIAL"
    break
  fi

  mid_depth="$RUN_DEPTH"
  mid_k="$RUN_K"
  if awk -v k="$mid_k" -v target="$TARGET_K" -v sigma="$CONFIDENCE_SIGMA" -v se="$RUN_STDERR" 'BEGIN {
    d=k-target;
    if (d<0) d=-d;
    exit !(d <= sigma*se);
  }'; then
    if (( STOP_ON_CONFIDENCE_INTERVAL == 1 )); then
      stop_reason="target k lies inside the averaged confidence interval"
      low_depth="$mid_depth"
      high_depth="$mid_depth"
      low_k="$mid_k"
      high_k="$mid_k"
      break
    fi
    echo "    target lies inside +/- $CONFIDENCE_SIGMA stderr, but confidence-interval early stop is disabled; continuing bisection using the averaged k."
  fi

  delta_pcm="$(awk -v k="$mid_k" -v target="$TARGET_K" 'BEGIN { d=k-target; if (d<0) d=-d; printf "%.8f", d*100000.0 }')"
  if awk -v delta="$delta_pcm" -v tolerance="$K_TOL_PCM" 'BEGIN { exit !(tolerance>0 && delta<=tolerance) }'; then
    stop_reason="k tolerance reached"
    low_depth="$mid_depth"
    high_depth="$mid_depth"
    low_k="$mid_k"
    high_k="$mid_k"
    break
  fi

  if brackets_target "$low_k" "$mid_k"; then
    high_depth="$mid_depth"
    high_k="$mid_k"
  else
    low_depth="$mid_depth"
    low_k="$mid_k"
  fi

  echo "    bracket is now $low_depth cm (k=$(format_k "$low_k")) to $high_depth cm (k=$(format_k "$high_k"))"
done

if awk -v low="$low_k" -v high="$high_k" 'BEGIN { exit !(low != high) }'; then
  estimated_depth="$(awk -v lo="$low_depth" -v hi="$high_depth" -v kl="$low_k" -v kh="$high_k" -v target="$TARGET_K" 'BEGIN {
    estimate=lo+(target-kl)*(hi-lo)/(kh-kl);
    if (estimate<lo) estimate=lo;
    if (estimate>hi) estimate=hi;
    printf "%.8f", estimate;
  }')"
else
  estimated_depth="$low_depth"
fi

{
  echo "STATUS=$search_status"
  echo "stop_reason=$stop_reason"
  echo "target_k=$(format_k "$TARGET_K")"
  echo "replicates_per_depth=$REPLICATES_PER_DEPTH"
  echo "base_seed=$BASE_SEED"
  echo "confidence_sigma=$CONFIDENCE_SIGMA"
  echo "stop_on_confidence_interval=$STOP_ON_CONFIDENCE_INTERVAL"
  echo "aro_check_status=$aro_check_status"
  echo "aro_check_depth_cm=$ARO_CHECK_DEPTH_CM"
  echo "aro_expected_k=$(format_k "$ARO_EXPECTED_K")"
  echo "aro_observed_k=$aro_check_k"
  echo "aro_k_stderr=$aro_check_stderr"
  echo "aro_difference_pcm=$aro_check_delta_pcm"
  echo "aro_allowed_pcm=$aro_check_limit_pcm"
  echo "final_lower_depth_cm=$low_depth final_lower_k=$(format_k "$low_k")"
  echo "final_upper_depth_cm=$high_depth final_upper_k=$(format_k "$high_k")"
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
