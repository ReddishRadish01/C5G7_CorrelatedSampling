#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

# Keep routinely edited physics/statistics settings in one small file. Values
# already present in the environment (including sbatch NAME=value arguments)
# are preserved by the default assignments in that file.
CONFIG_FILE="${CONFIG_FILE:-scripts/critical_rod_newton_cs.env}"
if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
fi

# Stochastic, safeguarded Newton iteration using one-sided correlated sampling.
# All settings can be overridden as environment variables or as NAME=value
# arguments to the accompanying Slurm batch script.
EXE="${EXE:-build/c5g7}"
NVCC="${NVCC:-nvcc}"
INITIAL_DEPTH_CM="${INITIAL_DEPTH_CM:-64.24}"
LOWER_CM="${LOWER_CM:-21.42}"
UPPER_CM="${UPPER_CM:-128.52}"
TARGET_K="${TARGET_K:-1.0}"
PROBE_STEP_CM="${PROBE_STEP_CM:-0.25}"
MIN_PROBE_STEP_CM="${MIN_PROBE_STEP_CM:-0.03125}"
MAX_PROBE_RETRIES="${MAX_PROBE_RETRIES:-3}"
DAMPING="${DAMPING:-1.0}"
MAX_NEWTON_STEP_CM="${MAX_NEWTON_STEP_CM:-30.0}"
MIN_ABS_SLOPE="${MIN_ABS_SLOPE:-1.0e-5}"
MIN_SLOPE_SIGMA="${MIN_SLOPE_SIGMA:-3.0}"
EXPECTED_SLOPE_SIGN="${EXPECTED_SLOPE_SIGN:--1}"
DEPTH_TOL_CM="${DEPTH_TOL_CM:-0.01}"
K_TOL_PCM="${K_TOL_PCM:-10.0}"
MAX_K_STDERR_PCM="${MAX_K_STDERR_PCM:-10.0}"
CONFIDENCE_SIGMA="${CONFIDENCE_SIGMA:-2.0}"
SOURCE_TAIL_ABS_PCM="${SOURCE_TAIL_ABS_PCM:-1.0}"
SOURCE_TAIL_REL_FRACTION="${SOURCE_TAIL_REL_FRACTION:-0.02}"
MAX_SLOPE_RELATIVE_DIFFERENCE="${MAX_SLOPE_RELATIVE_DIFFERENCE:-0.05}"
MAX_ITER="${MAX_ITER:-8}"
PARTICLES="${PARTICLES:-2000000}"
INACTIVE="${INACTIVE:-100}"
ACTIVE="${ACTIVE:-200}"
CS_GENERATIONS="${CS_GENERATIONS:-10}"
BASE_SEED="${BASE_SEED:-92235922383}"
REQUIRE_SOURCE_CONVERGENCE="${REQUIRE_SOURCE_CONVERGENCE:-1}"
FINAL_DIRECT_CHECK="${FINAL_DIRECT_CHECK:-1}"
CUDA_ARCH="${CUDA_ARCH:-native}"
THREADS_PER_BLOCK="${THREADS_PER_BLOCK:-64}"
DEVICE_LTO="${DEVICE_LTO:-1}"
FAST_MATH="${FAST_MATH:-1}"
NO_BUILD="${NO_BUILD:-0}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
LOG_DIR="${LOG_DIR:-results/critical_rod_newton_cs_${RUN_TAG}}"
SUMMARY="${SUMMARY:-$LOG_DIR/newton_history.csv}"
RESULT="${RESULT:-$LOG_DIR/result.txt}"

is_number() {
  awk -v value="$1" 'BEGIN {
    exit !(value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/)
  }'
}

is_positive_integer() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

is_nonnegative_integer() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

odd_48bit_seed() {
  local value=$(( $1 % 281474976710656 ))
  if (( value % 2 == 0 )); then
    value=$((value + 1))
  fi
  printf '%s\n' "$value"
}

abs_value() {
  awk -v value="$1" 'BEGIN { if (value < 0) value = -value; printf "%.12g", value }'
}

format_depth() {
  awk -v value="$1" 'BEGIN { printf "%.8f", value }'
}

clamp() {
  awk -v value="$1" -v lo="$2" -v hi="$3" 'BEGIN {
    if (value < lo) value = lo;
    if (value > hi) value = hi;
    printf "%.12g", value;
  }'
}

for value in \
  "$INITIAL_DEPTH_CM" "$LOWER_CM" "$UPPER_CM" "$TARGET_K" \
  "$PROBE_STEP_CM" "$MIN_PROBE_STEP_CM" "$DAMPING" \
  "$MAX_NEWTON_STEP_CM" "$MIN_ABS_SLOPE" "$MIN_SLOPE_SIGMA" "$DEPTH_TOL_CM" \
  "$K_TOL_PCM" "$MAX_K_STDERR_PCM" "$CONFIDENCE_SIGMA" "$SOURCE_TAIL_ABS_PCM" \
  "$SOURCE_TAIL_REL_FRACTION" "$MAX_SLOPE_RELATIVE_DIFFERENCE"; do
  if ! is_number "$value"; then
    echo "Invalid numeric setting: $value" >&2
    exit 2
  fi
done

if ! awk -v lo="$LOWER_CM" -v hi="$UPPER_CM" 'BEGIN { exit !(lo < hi) }'; then
  echo "LOWER_CM must be less than UPPER_CM." >&2
  exit 2
fi
if ! awk -v h="$PROBE_STEP_CM" -v hmin="$MIN_PROBE_STEP_CM" \
  'BEGIN { exit !(h > 0 && hmin > 0 && hmin <= h) }'; then
  echo "Require 0 < MIN_PROBE_STEP_CM <= PROBE_STEP_CM." >&2
  exit 2
fi
if ! awk -v d="$DAMPING" -v m="$MAX_NEWTON_STEP_CM" -v s="$MIN_ABS_SLOPE" -v ss="$MIN_SLOPE_SIGMA" \
  -v dt="$DEPTH_TOL_CM" -v kt="$K_TOL_PCM" -v se="$MAX_K_STDERR_PCM" -v c="$CONFIDENCE_SIGMA" \
  -v ta="$SOURCE_TAIL_ABS_PCM" -v tr="$SOURCE_TAIL_REL_FRACTION" -v sr="$MAX_SLOPE_RELATIVE_DIFFERENCE" \
  'BEGIN { exit !(d > 0 && d <= 1 && m > 0 && s > 0 && ss > 0 && dt > 0 && kt >= 0 && se >= 0 && c > 0 && ta >= 0 && tr >= 0 && sr >= 0 && (ta > 0 || tr > 0)) }'; then
  echo "Invalid damping, step, slope, tolerance, or confidence setting." >&2
  exit 2
fi
if ! [[ "$EXPECTED_SLOPE_SIGN" =~ ^(-1|0|1)$ ]]; then
  echo "EXPECTED_SLOPE_SIGN must be -1, 0, or 1." >&2
  exit 2
fi
if ! is_positive_integer "$MAX_ITER" || ! is_positive_integer "$PARTICLES" || \
  ! is_nonnegative_integer "$INACTIVE" || ! is_positive_integer "$ACTIVE" || (( ACTIVE < 2 )) || \
  ! is_positive_integer "$CS_GENERATIONS" || ! is_nonnegative_integer "$MAX_PROBE_RETRIES"; then
  echo "Invalid iteration, population, cycle, propagation, or retry count." >&2
  exit 2
fi
if ! is_nonnegative_integer "$BASE_SEED"; then
  echo "BASE_SEED must be a nonnegative integer." >&2
  exit 2
fi
if (( BASE_SEED > 281474976710655 )); then
  echo "BASE_SEED must fit in the RNG's 48-bit state (0..281474976710655)." >&2
  exit 2
fi
if ! [[ "$NO_BUILD" =~ ^[01]$ && "$DEVICE_LTO" =~ ^[01]$ && "$FAST_MATH" =~ ^[01]$ && \
  "$REQUIRE_SOURCE_CONVERGENCE" =~ ^[01]$ && "$FINAL_DIRECT_CHECK" =~ ^[01]$ ]]; then
  echo "Boolean settings must be 0 or 1." >&2
  exit 2
fi
if ! [[ "$THREADS_PER_BLOCK" =~ ^(32|64|128|256|512|1024)$ ]]; then
  echo "THREADS_PER_BLOCK must be one of 32, 64, 128, 256, 512, or 1024." >&2
  exit 2
fi
if (( CS_GENERATIONS > 10 )); then
  echo "CS_GENERATIONS cannot exceed this build's compiled maximum of 10." >&2
  exit 2
fi

if (( NO_BUILD == 0 )); then
  if ! command -v "$NVCC" >/dev/null 2>&1; then
    echo "nvcc not found: $NVCC. Load a CUDA module or set NVCC." >&2
    exit 1
  fi
  if ! [[ "$CUDA_ARCH" =~ ^(native|sm_[0-9]+[a-z]?|compute_[0-9]+[a-z]?)$ ]]; then
    echo "CUDA_ARCH must be native, sm_NN, or compute_NN; got $CUDA_ARCH." >&2
    exit 2
  fi
  echo "Building $EXE (CUDA_ARCH=$CUDA_ARCH)"
  make \
    TARGET="$EXE" \
    NVCC="$NVCC" \
    CUDA_ARCH="$CUDA_ARCH" \
    THREADS_PER_BLOCK="$THREADS_PER_BLOCK" \
    DEVICE_LTO="$DEVICE_LTO" \
    FAST_MATH="$FAST_MATH" \
    "$EXE"
fi

if [[ ! -x "$EXE" ]]; then
  echo "Executable not found or not executable: $EXE" >&2
  exit 1
fi

mkdir -p "$LOG_DIR"
echo "iteration,seed,requested_depth_cm,reference_depth_cm,perturbed_depth_cm,probe_step_cm,k_current,k_stderr,k_reference,k_perturbed,delta_k,delta_k_stderr,slope_per_cm,slope_stderr_per_cm,residual_pcm,last_source_increment,last_source_increment_stderr,source_tail_pcm,source_tail_upper_pcm,source_tail_limit_pcm,source_tail_fraction,probe_retry,probe_assessment,slope_consistency_difference,slope_consistency_limit,raw_newton_depth_cm,next_depth_cm,applied_step_cm,status,log_file" > "$SUMMARY"

depth="$(clamp "$INITIAL_DEPTH_CM" "$LOWER_CM" "$UPPER_CM")"
converged=0
last_step=""
last_k=""
last_k_stderr=""
last_slope=""
last_slope_stderr=""
last_log=""

echo "Correlated-sampling Newton critical-depth search"
echo "  initial/bounds: $depth cm / [$LOWER_CM, $UPPER_CM] cm"
echo "  target/tolerance: k=$TARGET_K / $K_TOL_PCM delta-k pcm"
echo "  population/cycles: $PARTICLES / $INACTIVE inactive + $ACTIVE active"
echo "  CS propagation generations: $CS_GENERATIONS"
echo "  base seed: $BASE_SEED (every generated run seed is forced odd)"
echo "  source-tail budget: max($SOURCE_TAIL_ABS_PCM pcm, $(awk -v x="$SOURCE_TAIL_REL_FRACTION" 'BEGIN {printf "%.6g", 100*x}')% of |delta k|)"
echo "  consecutive-probe slope tolerance: $(awk -v x="$MAX_SLOPE_RELATIVE_DIFFERENCE" 'BEGIN {printf "%.6g", 100*x}')% or $CONFIDENCE_SIGMA-sigma"
echo "  output: $LOG_DIR"

for ((iteration = 1; iteration <= MAX_ITER; ++iteration)); do
  probe_step="$PROBE_STEP_CM"
  retry=0
  previous_probe_slope=""
  previous_probe_slope_stderr=""
  slope_consistency_difference=""
  slope_consistency_limit=""
  slope_consistent=1
  probe_assessment="not-assessed"

  while true; do
    # The implementation supports only an increasing-depth CS perturbation.
    # Use a backward stencil only at the upper physical boundary.
    if awk -v d="$depth" -v h="$probe_step" -v hi="$UPPER_CM" 'BEGIN { exit !(d + h <= hi) }'; then
      reference_depth="$depth"
      perturbed_depth="$(awk -v d="$depth" -v h="$probe_step" 'BEGIN { printf "%.12g", d + h }')"
      current_side="reference"
    else
      reference_depth="$(awk -v d="$depth" -v h="$probe_step" -v lo="$LOWER_CM" \
        'BEGIN { value=d-h; if (value < lo) value=lo; printf "%.12g", value }')"
      perturbed_depth="$depth"
      current_side="perturbed"
    fi

    actual_probe_step="$(awk -v a="$reference_depth" -v b="$perturbed_depth" 'BEGIN { printf "%.12g", b-a }')"
    if ! awk -v h="$actual_probe_step" 'BEGIN { exit !(h > 0) }'; then
      echo "Cannot form a positive CS stencil at depth $depth cm." >&2
      exit 3
    fi

    seed="$(odd_48bit_seed "$((BASE_SEED + iteration * 104729 + retry * 1009))")"
    log_file="$LOG_DIR/iter_$(printf '%02d' "$iteration")_retry_${retry}_d$(format_depth "$depth").log"

    echo
    echo "Iteration $iteration, retry $retry: d=$depth cm, CS [$reference_depth -> $perturbed_depth] cm, seed=$seed"
    set +e
    "$EXE" \
      --critical-search-binary \
      --rod-depth="$reference_depth" \
      --cs-depth="$perturbed_depth" \
      --cs-generations="$CS_GENERATIONS" \
      --particles="$PARTICLES" \
      --inactive="$INACTIVE" \
      --active="$ACTIVE" \
      --seed="$seed" 2>&1 | tee "$log_file"
    run_status=${PIPESTATUS[0]}
    set -e
    if (( run_status != 0 )); then
      echo "Simulation failed with exit code $run_status; see $log_file" >&2
      exit 3
    fi

    run_seconds="$(awk '/^Total Time:/ {value=$3; sub(/seconds.*/, "", value); print value; exit}' "$log_file")"
    run_seconds_per_cycle="$(awk '/^Total Time:/ {print $5; exit}' "$log_file")"
    if is_number "$run_seconds" && is_number "$run_seconds_per_cycle" && \
      awk -v value="$run_seconds_per_cycle" 'BEGIN { exit !(value > 0) }'; then
      run_million_particles_per_second="$(awk -v particles="$PARTICLES" -v seconds="$run_seconds_per_cycle" \
        'BEGIN { printf "%.6f", particles/(seconds*1.0e6) }')"
      echo "  performance: $run_seconds s total, $run_seconds_per_cycle s/cycle, $run_million_particles_per_second million particles/s"
    fi

    k_reference="$(awk -F': ' '/^Paired reference k:/ {value=$2} END {print value}' "$log_file")"
    k_reference_stderr="$(awk -F': ' '/^k_stderr\(mean\):/ {value=$2} END {print value}' "$log_file")"
    k_perturbed="$(awk -F': ' '/^Perturbed k, NPS \+ fission-source propagation/ {value=$2} END {print value}' "$log_file")"
    read -r delta_k delta_k_stderr < <(awk '/^Total delta k:/ {
      for (i=1; i<=NF; ++i) if ($i == "+/-") {mean=$(i-1); stderr=$(i+1)}
    } END {print mean, stderr}' "$log_file")
    read -r source_increment source_increment_stderr < <(awk '/^Last propagation increment,/ {
      for (i=1; i<=NF; ++i) if ($i == "+/-") {mean=$(i-1); stderr=$(i+1)}
    } END {print mean, stderr}' "$log_file")
    diagnostic_count="$(awk -F': ' '/^Support\/non-finite diagnostic count:/ {value=$2} END {print value}' "$log_file")"

    for parsed in "$k_reference" "$k_reference_stderr" "$k_perturbed" "$delta_k" \
      "$delta_k_stderr" "$source_increment" "$source_increment_stderr"; do
      if ! is_number "$parsed"; then
        echo "Could not parse a finite CS result from $log_file" >&2
        exit 3
      fi
    done
    if ! [[ "$diagnostic_count" =~ ^[0-9]+$ ]] || (( diagnostic_count != 0 )); then
      echo "CS support/non-finite diagnostic count is ${diagnostic_count:-missing}; refusing this estimate." >&2
      exit 3
    fi

    slope="$(awk -v dk="$delta_k" -v h="$actual_probe_step" 'BEGIN { printf "%.12g", dk/h }')"
    slope_stderr="$(awk -v se="$delta_k_stderr" -v h="$actual_probe_step" 'BEGIN { printf "%.12g", se/h }')"
    read -r source_tail_pcm source_tail_upper_pcm source_tail_limit_pcm source_tail_fraction \
      source_tail_significant source_tail_practical < <(awk \
      -v inc="$source_increment" -v incse="$source_increment_stderr" -v dk="$delta_k" \
      -v sigma="$CONFIDENCE_SIGMA" -v absolute="$SOURCE_TAIL_ABS_PCM" -v relative="$SOURCE_TAIL_REL_FRACTION" \
      'BEGIN {
        if (inc < 0) inc=-inc;
        if (dk < 0) dk=-dk;
        tail=inc*1.0e5;
        upper=(inc+sigma*incse)*1.0e5;
        relative_limit=relative*dk*1.0e5;
        limit=(absolute > relative_limit ? absolute : relative_limit);
        fraction=(dk > 0 ? inc/dk : 1.0e300);
        significant=(inc > sigma*incse ? 1 : 0);
        practical=(upper <= limit ? 1 : 0);
        printf "%.12g %.12g %.12g %.12g %d %d\n", tail, upper, limit, fraction, significant, practical;
      }')

    slope_consistent=1
    if (( retry > 0 )); then
      read -r slope_consistency_difference slope_consistency_limit slope_consistent < <(awk \
        -v current="$slope" -v previous="$previous_probe_slope" \
        -v current_se="$slope_stderr" -v previous_se="$previous_probe_slope_stderr" \
        -v sigma="$CONFIDENCE_SIGMA" -v relative="$MAX_SLOPE_RELATIVE_DIFFERENCE" \
        'BEGIN {
          difference=current-previous; if (difference < 0) difference=-difference;
          scale=current; if (scale < 0) scale=-scale;
          previous_scale=previous; if (previous_scale < 0) previous_scale=-previous_scale;
          if (previous_scale > scale) scale=previous_scale;
          statistical=sigma*sqrt(current_se*current_se+previous_se*previous_se);
          practical=relative*scale;
          limit=(statistical > practical ? statistical : practical);
          printf "%.12g %.12g %d\n", difference, limit, (difference <= limit ? 1 : 0);
        }')
    fi

    echo "  source tail: $source_tail_pcm pcm; confidence upper=$source_tail_upper_pcm pcm; allowed=$source_tail_limit_pcm pcm"
    if (( retry > 0 )); then
      echo "  consecutive-probe slope difference=$slope_consistency_difference; allowed=$slope_consistency_limit per cm"
    fi

    needs_smaller_probe=0
    validation_reason=""
    if (( REQUIRE_SOURCE_CONVERGENCE == 1 )); then
      if (( source_tail_practical == 0 )); then
        needs_smaller_probe=1
        validation_reason="the source-tail confidence bound exceeds its practical error budget"
      elif (( retry == 0 && source_tail_significant == 1 )); then
        needs_smaller_probe=1
        validation_reason="the source tail is statistically resolved, so a half-step slope check is required"
      elif (( retry > 0 && slope_consistent == 0 )); then
        needs_smaller_probe=1
        validation_reason="the consecutive probe slopes are inconsistent"
      fi
    fi

    if (( needs_smaller_probe == 1 )); then
      next_probe="$(awk -v h="$probe_step" 'BEGIN { printf "%.12g", h/2 }')"
      if (( retry < MAX_PROBE_RETRIES )) && \
        awk -v candidate="$next_probe" -v minimum="$MIN_PROBE_STEP_CM" 'BEGIN { exit !(candidate >= minimum) }'; then
        echo "$validation_reason; validating with h=$next_probe cm."
        previous_probe_slope="$slope"
        previous_probe_slope_stderr="$slope_stderr"
        probe_step="$next_probe"
        retry=$((retry + 1))
        continue
      fi
      {
        echo "status=CS_PROBE_VALIDATION_FAILED"
        echo "reason=$validation_reason"
        echo "depth_cm=$depth"
        echo "probe_step_cm=$actual_probe_step"
        echo "source_tail_upper_pcm=$source_tail_upper_pcm"
        echo "source_tail_limit_pcm=$source_tail_limit_pcm"
        echo "slope_consistency_difference=${slope_consistency_difference:-unavailable}"
        echo "slope_consistency_limit=${slope_consistency_limit:-unavailable}"
        echo "last_log=$log_file"
      } > "$RESULT"
      echo "CS probe validation failed: $validation_reason. See $log_file" >&2
      exit 3
    fi

    if (( REQUIRE_SOURCE_CONVERGENCE == 0 )); then
      probe_assessment="checks-disabled"
    elif (( retry > 0 )); then
      probe_assessment="practical-tail-and-slope-consistent"
    else
      probe_assessment="practical-tail-not-statistically-resolved"
    fi
    break
  done

  slope="$(awk -v dk="$delta_k" -v h="$actual_probe_step" 'BEGIN { printf "%.12g", dk/h }')"
  slope_stderr="$(awk -v se="$delta_k_stderr" -v h="$actual_probe_step" 'BEGIN { printf "%.12g", se/h }')"

  if [[ "$current_side" == "reference" ]]; then
    k_current="$k_reference"
    k_current_stderr="$k_reference_stderr"
  else
    k_current="$k_perturbed"
    # Conservative because the perturbed-k covariance is not printed separately.
    k_current_stderr="$(awk -v a="$k_reference_stderr" -v b="$delta_k_stderr" \
      'BEGIN { printf "%.12g", sqrt(a*a+b*b) }')"
  fi

  residual="$(awk -v k="$k_current" -v target="$TARGET_K" 'BEGIN { printf "%.12g", k-target }')"
  residual_pcm="$(awk -v value="$residual" 'BEGIN { printf "%.6f", value*1.0e5 }')"
  abs_residual_pcm="$(abs_value "$residual_pcm")"
  k_stderr_pcm="$(awk -v value="$k_current_stderr" 'BEGIN { printf "%.6f", value*1.0e5 }')"
  abs_slope="$(abs_value "$slope")"

  echo "  k(d)=$k_current +/- $k_current_stderr; residual=$residual_pcm pcm"
  echo "  dk/dd=$slope +/- $slope_stderr per cm"

  stderr_ok=1
  if awk -v limit="$MAX_K_STDERR_PCM" 'BEGIN { exit !(limit > 0) }' && \
    ! awk -v se="$k_stderr_pcm" -v limit="$MAX_K_STDERR_PCM" 'BEGIN { exit !(se <= limit) }'; then
    stderr_ok=0
  fi
  if awk -v error="$abs_residual_pcm" -v tolerance="$K_TOL_PCM" 'BEGIN { exit !(error <= tolerance) }' && \
    (( stderr_ok == 1 )); then
    converged=1
    last_step="0"
    last_k="$k_current"
    last_k_stderr="$k_current_stderr"
    last_slope="$slope"
    last_slope_stderr="$slope_stderr"
    last_log="$log_file"
    echo "$iteration,$seed,$depth,$reference_depth,$perturbed_depth,$actual_probe_step,$k_current,$k_current_stderr,$k_reference,$k_perturbed,$delta_k,$delta_k_stderr,$slope,$slope_stderr,$residual_pcm,$source_increment,$source_increment_stderr,$source_tail_pcm,$source_tail_upper_pcm,$source_tail_limit_pcm,$source_tail_fraction,$retry,$probe_assessment,$slope_consistency_difference,$slope_consistency_limit,,,$last_step,converged,$log_file" >> "$SUMMARY"
    echo "Converged at d=$depth cm."
    break
  fi

  if ! awk -v slope="$abs_slope" -v minimum="$MIN_ABS_SLOPE" 'BEGIN { exit !(slope >= minimum) }'; then
    echo "|dk/dd|=$abs_slope is below MIN_ABS_SLOPE=$MIN_ABS_SLOPE; Newton update is unsafe." >&2
    exit 3
  fi
  if ! awk -v slope="$abs_slope" -v se="$slope_stderr" -v sigma="$MIN_SLOPE_SIGMA" \
    'BEGIN { exit !(se <= 0 || slope >= sigma*se) }'; then
    echo "The CS slope is not significant at MIN_SLOPE_SIGMA=$MIN_SLOPE_SIGMA; increase PROBE_STEP_CM, PARTICLES, or ACTIVE." >&2
    exit 3
  fi
  if (( EXPECTED_SLOPE_SIGN != 0 )) && ! awk -v slope="$slope" -v sign="$EXPECTED_SLOPE_SIGN" \
    'BEGIN { exit !((sign < 0 && slope < 0) || (sign > 0 && slope > 0)) }'; then
    echo "Slope $slope has the wrong sign (EXPECTED_SLOPE_SIGN=$EXPECTED_SLOPE_SIGN)." >&2
    exit 3
  fi

  raw_newton_depth="$(awk -v d="$depth" -v r="$residual" -v slope="$slope" \
    'BEGIN { printf "%.12g", d-r/slope }')"
  raw_step="$(awk -v candidate="$raw_newton_depth" -v d="$depth" 'BEGIN { printf "%.12g", candidate-d }')"
  limited_step="$(awk -v step="$raw_step" -v maximum="$MAX_NEWTON_STEP_CM" 'BEGIN {
    if (step > maximum) step=maximum;
    if (step < -maximum) step=-maximum;
    printf "%.12g", step;
  }')"
  applied_step="$(awk -v step="$limited_step" -v damping="$DAMPING" 'BEGIN { printf "%.12g", damping*step }')"
  next_depth="$(awk -v d="$depth" -v step="$applied_step" 'BEGIN { printf "%.12g", d+step }')"
  next_depth="$(clamp "$next_depth" "$LOWER_CM" "$UPPER_CM")"
  applied_step="$(awk -v candidate="$next_depth" -v d="$depth" 'BEGIN { printf "%.12g", candidate-d }')"
  abs_applied_step="$(abs_value "$applied_step")"

  echo "$iteration,$seed,$depth,$reference_depth,$perturbed_depth,$actual_probe_step,$k_current,$k_current_stderr,$k_reference,$k_perturbed,$delta_k,$delta_k_stderr,$slope,$slope_stderr,$residual_pcm,$source_increment,$source_increment_stderr,$source_tail_pcm,$source_tail_upper_pcm,$source_tail_limit_pcm,$source_tail_fraction,$retry,$probe_assessment,$slope_consistency_difference,$slope_consistency_limit,$raw_newton_depth,$next_depth,$applied_step,continue,$log_file" >> "$SUMMARY"
  echo "  raw Newton depth=$raw_newton_depth cm; accepted next depth=$next_depth cm"

  if awk -v step="$abs_applied_step" -v tolerance="$DEPTH_TOL_CM" 'BEGIN { exit !(step <= tolerance) }'; then
    echo "Newton step is below DEPTH_TOL_CM, but the k/statistical tolerance is not met." >&2
    last_step="$applied_step"
    last_k="$k_current"
    last_k_stderr="$k_current_stderr"
    last_slope="$slope"
    last_slope_stderr="$slope_stderr"
    last_log="$log_file"
    break
  fi

  last_step="$applied_step"
  last_k="$k_current"
  last_k_stderr="$k_current_stderr"
  last_slope="$slope"
  last_slope_stderr="$slope_stderr"
  last_log="$log_file"
  depth="$next_depth"
done

if (( converged == 0 )); then
  {
    echo "status=NOT_CONVERGED"
    echo "depth_cm=$depth"
    echo "last_k=${last_k:-unavailable}"
    echo "last_k_stderr=${last_k_stderr:-unavailable}"
    echo "last_slope_per_cm=${last_slope:-unavailable}"
    echo "source_tail_upper_pcm=${source_tail_upper_pcm:-unavailable}"
    echo "source_tail_limit_pcm=${source_tail_limit_pcm:-unavailable}"
    echo "probe_assessment=${probe_assessment:-unavailable}"
    echo "history=$SUMMARY"
    echo "last_log=${last_log:-unavailable}"
  } > "$RESULT"
  echo "Newton search did not meet the requested tolerances. Result: $RESULT" >&2
  exit 4
fi

verification_status="SKIPPED"
verification_k=""
verification_stderr=""
verification_error_pcm=""
verification_log=""

if (( FINAL_DIRECT_CHECK == 1 )); then
  verification_seed="$(odd_48bit_seed "$((BASE_SEED + (MAX_ITER + 1) * 104729))")"
  verification_log="$LOG_DIR/final_direct_d$(format_depth "$depth").log"
  echo
  echo "Independent direct verification at d=$depth cm, seed=$verification_seed"
  set +e
  "$EXE" \
    --critical-search-binary \
    --rod-depth="$depth" \
    --particles="$PARTICLES" \
    --inactive="$INACTIVE" \
    --active="$ACTIVE" \
    --seed="$verification_seed" 2>&1 | tee "$verification_log"
  run_status=${PIPESTATUS[0]}
  set -e
  if (( run_status != 0 )); then
    echo "Final direct verification failed; see $verification_log" >&2
    exit 5
  fi

  verification_k="$(awk -F': ' '/^k_mean:/ {value=$2} END {print value}' "$verification_log")"
  verification_stderr="$(awk -F': ' '/^k_stderr\(mean\):/ {value=$2} END {print value}' "$verification_log")"
  if ! is_number "$verification_k" || ! is_number "$verification_stderr"; then
    echo "Could not parse final direct verification output." >&2
    exit 5
  fi
  verification_error_pcm="$(awk -v k="$verification_k" -v target="$TARGET_K" \
    'BEGIN { value=(k-target)*1.0e5; if (value < 0) value=-value; printf "%.6f", value }')"
  verification_limit_pcm="$(awk -v tolerance="$K_TOL_PCM" -v sigma="$CONFIDENCE_SIGMA" -v se="$verification_stderr" \
    'BEGIN { printf "%.6f", tolerance+sigma*se*1.0e5 }')"
  if awk -v error="$verification_error_pcm" -v limit="$verification_limit_pcm" 'BEGIN { exit !(error <= limit) }'; then
    verification_status="PASSED"
  else
    verification_status="FAILED"
  fi
  echo "  direct k=$verification_k +/- $verification_stderr"
  echo "  |k-target|=$verification_error_pcm pcm; acceptance limit=$verification_limit_pcm pcm"
fi

depth_uncertainty_cm="$(awk -v se="$last_k_stderr" -v slope="$last_slope" 'BEGIN {
  if (slope < 0) slope=-slope;
  if (slope > 0) printf "%.12g", se/slope; else print "nan";
}')"

{
  echo "status=CONVERGED"
  echo "critical_depth_cm=$depth"
  echo "cs_k=$last_k"
  echo "cs_k_stderr=$last_k_stderr"
  echo "cs_slope_per_cm=$last_slope"
  echo "cs_slope_stderr_per_cm=$last_slope_stderr"
  echo "source_tail_upper_pcm=$source_tail_upper_pcm"
  echo "source_tail_limit_pcm=$source_tail_limit_pcm"
  echo "probe_assessment=$probe_assessment"
  echo "approx_depth_stderr_cm=$depth_uncertainty_cm"
  echo "target_k=$TARGET_K"
  echo "verification_status=$verification_status"
  echo "verification_k=${verification_k:-not_run}"
  echo "verification_k_stderr=${verification_stderr:-not_run}"
  echo "verification_abs_error_pcm=${verification_error_pcm:-not_run}"
  echo "history=$SUMMARY"
  echo "last_cs_log=$last_log"
  echo "verification_log=${verification_log:-not_run}"
} > "$RESULT"

echo
echo "Critical rod depth: $depth cm"
echo "Approximate one-sigma depth uncertainty: $depth_uncertainty_cm cm"
echo "Verification: $verification_status"
echo "Result: $RESULT"

if [[ "$verification_status" == "FAILED" ]]; then
  echo "The CS Newton estimate converged, but the independent direct run did not verify it statistically." >&2
  exit 5
fi
