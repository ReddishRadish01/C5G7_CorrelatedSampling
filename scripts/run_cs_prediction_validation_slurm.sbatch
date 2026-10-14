#!/usr/bin/env bash
#SBATCH --job-name=c5g7-cs-valid
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

set -euo pipefail

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
PROJECT_DIR="$(pwd -P)"

CONFIG_FILE="${CONFIG_FILE:-scripts/cs_prediction_validation_slurm.env}"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Missing configuration file: $CONFIG_FILE" >&2
  exit 2
fi
# Values already exported through NAME=value arguments keep precedence because
# the configuration file supplies defaults with ${NAME:=value}.
# shellcheck disable=SC1090
source "$CONFIG_FILE"

is_positive_integer() { [[ "$1" =~ ^[1-9][0-9]*$ ]]; }
is_nonnegative_integer() { [[ "$1" =~ ^[0-9]+$ ]]; }
is_finite_number() {
  [[ "$1" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]
}
make_odd_seed() {
  awk -v raw="$1" 'BEGIN {
    modulus = 281474976710656;
    value = raw % modulus;
    if (value < 0) value += modulus;
    if ((value % 2) == 0) value += 1;
    printf "%.0f", value;
  }'
}
depth_tag() {
  awk -v d="$1" 'BEGIN {
    text = sprintf("%012.6f", d + 0);
    gsub(/-/, "m", text); gsub(/\./, "p", text);
    print text;
  }'
}
normalize_numeric_gpu_list() {
  local raw="${1//[[:space:]]/}" token first last value
  local -a tokens=() expanded=() sorted=()
  IFS=',' read -r -a tokens <<< "$raw"
  for token in "${tokens[@]}"; do
    if [[ "$token" =~ ^[0-9]+$ ]]; then
      expanded+=("$((10#$token))")
    elif [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      first=$((10#${BASH_REMATCH[1]}))
      last=$((10#${BASH_REMATCH[2]}))
      (( first <= last )) || return 1
      for ((value = first; value <= last; ++value)); do expanded+=("$value"); done
    else
      return 1
    fi
  done
  (( ${#expanded[@]} > 0 )) || return 1
  mapfile -t sorted < <(printf '%s\n' "${expanded[@]}" | sort -n -u)
  local joined=""
  for value in "${sorted[@]}"; do joined+="${joined:+,}$value"; done
  printf '%s' "$joined"
}

if [[ ! -f Makefile || ! -f main.cu || ! -f scripts/analyze_cs_prediction_validation_slurm.py ]]; then
  echo "The current directory is not a complete C5G7 project: $PROJECT_DIR" >&2
  exit 2
fi
if ! is_positive_integer "$PARTICLES" || ! is_nonnegative_integer "$INACTIVE" ||
   ! is_positive_integer "$ACTIVE" || ! is_positive_integer "$CS_GENERATIONS" ||
   ! is_positive_integer "$GPU_WORKERS" || ! is_positive_integer "$MONITOR_INTERVAL_SEC"; then
  echo "Invalid particle, cycle, CS level, worker, or monitor setting." >&2
  exit 2
fi
if (( GPU_WORKERS != 6 )); then
  echo "GPU_WORKERS must remain exactly 6 for the enforced node3 full-GPU allocation." >&2
  exit 2
fi
if (( CS_GENERATIONS > 10 )); then
  echo "CS_GENERATIONS cannot exceed the compiled maximum of 10." >&2
  exit 2
fi
if (( INACTIVE + ACTIVE <= CS_GENERATIONS )); then
  echo "Not enough cycles for the requested source-propagation level." >&2
  exit 2
fi
if ! is_finite_number "$PROBE_STEP_CM" ||
   ! awk -v x="$PROBE_STEP_CM" 'BEGIN { exit !(x > 0) }'; then
  echo "PROBE_STEP_CM must be a positive finite number." >&2
  exit 2
fi
if ! is_finite_number "$CONFIDENCE_SIGMA" ||
   ! is_finite_number "$SOURCE_TAIL_ABS_PCM" ||
   ! is_finite_number "$SOURCE_TAIL_REL_FRACTION" ||
   ! awk -v c="$CONFIDENCE_SIGMA" -v a="$SOURCE_TAIL_ABS_PCM" -v r="$SOURCE_TAIL_REL_FRACTION" \
      'BEGIN { exit !(c > 0 && a >= 0 && r >= 0) }'; then
  echo "Invalid confidence or source-tail setting." >&2
  exit 2
fi
if [[ "$EXPECTED_DELTA_K_SIGN" != "-1" && "$EXPECTED_DELTA_K_SIGN" != "1" ]]; then
  echo "EXPECTED_DELTA_K_SIGN must be -1 or 1." >&2
  exit 2
fi

read -r -a REFERENCE_DEPTHS <<< "$REFERENCE_DEPTHS_CM"
if (( ${#REFERENCE_DEPTHS[@]} == 0 )); then
  echo "REFERENCE_DEPTHS_CM is empty." >&2
  exit 2
fi
for depth in "${REFERENCE_DEPTHS[@]}"; do
  if ! is_finite_number "$depth"; then
    echo "Invalid reference depth: $depth" >&2
    exit 2
  fi
  if ! awk -v d="$depth" -v h="$PROBE_STEP_CM" -v upper="$MODELED_UPPER_DEPTH_CM" \
      'BEGIN { exit !(d >= 0 && d + h <= upper) }'; then
    echo "Depth pair $depth + $PROBE_STEP_CM cm is outside the modeled interval." >&2
    exit 2
  fi
done

# Fail closed before compilation or simulation. A generic `--gres=gpu:6`
# request does not guarantee physical indices, so the job proceeds only when
# Slurm assigned exactly all six node3 GPUs: 0,1,2,3,4,5.
REQUIRED_PHYSICAL_GPUS="0,1,2,3,4,5"
if [[ "${SLURM_JOB_NODELIST:-}" != "node3" ]]; then
  echo "SAFETY ABORT: expected node3, got '${SLURM_JOB_NODELIST:-unset}'." >&2
  echo "No compilation or simulation process was started." >&2
  exit 3
fi
if [[ -z "${SLURM_JOB_GPUS:-}" ]]; then
  echo "SAFETY ABORT: SLURM_JOB_GPUS is unavailable, so physical GPU IDs cannot be verified." >&2
  echo "No compilation or simulation process was started; the full six-GPU allocation is required." >&2
  exit 3
fi
if ! ALLOCATED_PHYSICAL_GPUS="$(normalize_numeric_gpu_list "$SLURM_JOB_GPUS")"; then
  echo "SAFETY ABORT: cannot interpret SLURM_JOB_GPUS='$SLURM_JOB_GPUS' as physical GPU indices." >&2
  echo "No compilation or simulation process was started." >&2
  exit 3
fi
if [[ "$ALLOCATED_PHYSICAL_GPUS" != "$REQUIRED_PHYSICAL_GPUS" ]]; then
  echo "SAFETY ABORT: Slurm allocated physical GPUs [$ALLOCATED_PHYSICAL_GPUS]," >&2
  echo "but this job is allowed to use only [$REQUIRED_PHYSICAL_GPUS]." >&2
  echo "No compilation or simulation process was started; resubmit after all node3 GPUs are available." >&2
  exit 3
fi

if [[ -n "${GCC_MODULE:-}" || -n "${CUDA_MODULE:-}" ]]; then
  if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if ! type module >/dev/null 2>&1; then
    echo "A compiler/CUDA module was requested, but the module command is unavailable." >&2
    exit 1
  fi
fi
if [[ -n "${GCC_MODULE:-}" ]]; then
  module load "$GCC_MODULE"
fi
if [[ -n "${CUDA_MODULE:-}" ]]; then
  module load "$CUDA_MODULE"
fi

export CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
export CUDA_MODULE_LOADING="${CUDA_MODULE_LOADING:-LAZY}"
export OMP_NUM_THREADS="${OMP_THREADS_PER_GPU:-2}"

mkdir -p results
if [[ -z "${RUN_TAG:-}" ]]; then
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    RUN_TAG="slurm_${SLURM_JOB_ID}"
  else
    RUN_TAG="manual_$(date +%Y%m%d_%H%M%S)"
  fi
fi
RUN_DIR="${RUN_DIR:-results/cs_prediction_validation_${RUN_TAG}}"
mkdir -p "$RUN_DIR/logs" "$RUN_DIR/work" "$RUN_DIR/status"

echo "Building validation executable..."
if [[ "$NO_BUILD" == "1" ]]; then
  if [[ ! -x "$BUILD_TARGET" ]]; then
    echo "NO_BUILD=1, but executable is missing: $BUILD_TARGET" >&2
    exit 2
  fi
elif [[ "$NO_BUILD" == "0" ]]; then
  if ! command -v "$NVCC_BIN" >/dev/null 2>&1; then
    echo "CUDA compiler not found: $NVCC_BIN" >&2
    echo "Load the server CUDA module or set NVCC_BIN=/absolute/path/to/nvcc." >&2
    exit 1
  fi
  if [[ -z "$HOST_COMPILER" ]]; then
    HOST_COMPILER="$(command -v g++ || true)"
  fi
  if [[ -z "$HOST_COMPILER" || ! -x "$HOST_COMPILER" ]]; then
    echo "A working host C++ compiler was not found." >&2
    echo "Load a GCC module (GCC_MODULE=...) or set HOST_COMPILER=/path/to/g++." >&2
    echo "If build/c5g7 is already current, resubmit with NO_BUILD=1 BUILD_TARGET=build/c5g7." >&2
    exit 1
  fi
  cc1plus_path="$("$HOST_COMPILER" -print-prog-name=cc1plus 2>/dev/null || true)"
  if [[ -z "$cc1plus_path" ]] ||
     { [[ ! -x "$cc1plus_path" ]] && ! command -v "$cc1plus_path" >/dev/null 2>&1; }; then
    echo "Host compiler is incomplete: $HOST_COMPILER cannot locate cc1plus." >&2
    echo "This is a GCC/G++ installation or module issue, not a CUDA source error." >&2
    echo "Use 'module avail gcc' and resubmit with GCC_MODULE=<available-module>." >&2
    echo "Alternatively reuse a current binary with NO_BUILD=1 BUILD_TARGET=build/c5g7." >&2
    exit 1
  fi
  nvcc_make_command="$NVCC_BIN -ccbin $HOST_COMPILER"
  make_args=(
    -j "$BUILD_JOBS"
    "NVCC=$nvcc_make_command"
    "TARGET=$BUILD_TARGET"
    "CUDA_ARCH=$CUDA_ARCH"
    "THREADS_PER_BLOCK=$THREADS_PER_BLOCK"
    "DEVICE_LTO=$DEVICE_LTO"
    "FAST_MATH=$FAST_MATH"
  )
  if [[ "$FORCE_BUILD" == "1" ]]; then
    make_args+=(-B)
  elif [[ "$FORCE_BUILD" != "0" ]]; then
    echo "FORCE_BUILD must be 0 or 1." >&2
    exit 2
  fi
  make "${make_args[@]}"
else
  echo "NO_BUILD must be 0 or 1." >&2
  exit 2
fi

EXE="$(readlink -f "$BUILD_TARGET")"
RUN_DIR_ABS="$(readlink -f "$RUN_DIR")"
MANIFEST="$RUN_DIR_ABS/tasks.tsv"
SETTINGS="$RUN_DIR_ABS/run_settings.txt"

{
  echo "project_dir=$PROJECT_DIR"
  echo "slurm_job_id=${SLURM_JOB_ID:-manual}"
  echo "node=${SLURM_JOB_NODELIST:-unknown}"
  echo "reference_depths_cm=$REFERENCE_DEPTHS_CM"
  echo "probe_step_cm=$PROBE_STEP_CM"
  echo "particles=$PARTICLES"
  echo "inactive=$INACTIVE"
  echo "active=$ACTIVE"
  echo "cs_generations=$CS_GENERATIONS"
  echo "base_seed=$BASE_SEED"
  echo "confidence_sigma=$CONFIDENCE_SIGMA"
  echo "expected_delta_k_sign=$EXPECTED_DELTA_K_SIGN"
  echo "source_tail_abs_pcm=$SOURCE_TAIL_ABS_PCM"
  echo "source_tail_relative_fraction=$SOURCE_TAIL_REL_FRACTION"
  echo "executable=$EXE"
  echo "cuda_arch=$CUDA_ARCH"
  echo "threads_per_block=$THREADS_PER_BLOCK"
  echo "device_lto=$DEVICE_LTO"
  echo "fast_math=$FAST_MATH"
  echo "nvcc_bin=$NVCC_BIN"
  echo "host_compiler=${HOST_COMPILER:-not-used}"
  echo "gcc_module=${GCC_MODULE:-not-loaded}"
  echo "cuda_module=${CUDA_MODULE:-not-loaded}"
} > "$SETTINGS"

printf 'task_id\tcase_index\trole\treference_depth_cm\ttarget_depth_cm\tseed\tlog_path\twork_dir\n' > "$MANIFEST"
task_id=0

# Role-major ordering spreads the expensive CS cases over the five fixed GPU
# workers before assigning the direct calculations.
for role in cs direct_ref direct_target; do
  for ((case_index = 0; case_index < ${#REFERENCE_DEPTHS[@]}; ++case_index)); do
    reference="${REFERENCE_DEPTHS[$case_index]}"
    target="$(awk -v d="$reference" -v h="$PROBE_STEP_CM" 'BEGIN {printf "%.12g", d+h}')"
    case_tag="$(depth_tag "$reference")"
    case_seed_base="$(awk -v b="$BASE_SEED" -v i="$case_index" 'BEGIN {printf "%.0f", b + 6*i}')"
    case "$role" in
      direct_ref) seed_offset=0 ;;
      direct_target) seed_offset=2 ;;
      cs) seed_offset=4 ;;
    esac
    raw_seed="$(awk -v b="$case_seed_base" -v o="$seed_offset" 'BEGIN {printf "%.0f", b+o}')"
    seed="$(make_odd_seed "$raw_seed")"
    log_path="$RUN_DIR_ABS/logs/case_$(printf '%02d' "$case_index")_${case_tag}_${role}.log"
    work_dir="$RUN_DIR_ABS/work/case_$(printf '%02d' "$case_index")_${case_tag}_${role}"
    printf '%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$task_id" "$case_index" "$role" "$reference" "$target" "$seed" "$log_path" "$work_dir" >> "$MANIFEST"
    task_id=$((task_id + 1))
  done
done

# Pin by immutable NVIDIA UUID rather than CUDA ordinal. CUDA ordinals can be
# renumbered by Slurm/cgroups; UUIDs still identify physical GPUs 0-5 exactly.
if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "SAFETY ABORT: nvidia-smi is unavailable, so GPU UUIDs cannot be verified." >&2
  echo "No simulation process was started." >&2
  exit 3
fi
declare -A gpu_uuid_by_physical_index=()
while IFS=',' read -r physical_index gpu_uuid; do
  physical_index="${physical_index//[[:space:]]/}"
  gpu_uuid="${gpu_uuid//[[:space:]]/}"
  [[ "$physical_index" =~ ^[0-9]+$ && "$gpu_uuid" == GPU-* ]] || continue
  gpu_uuid_by_physical_index[$physical_index]="$gpu_uuid"
done < <(nvidia-smi --query-gpu=index,uuid --format=csv,noheader,nounits)

GPU_TOKENS=()
for physical_index in 0 1 2 3 4 5; do
  gpu_uuid="${gpu_uuid_by_physical_index[$physical_index]:-}"
  if [[ -z "$gpu_uuid" ]]; then
    echo "SAFETY ABORT: could not resolve the UUID of physical GPU $physical_index." >&2
    echo "No simulation process was started." >&2
    exit 3
  fi
  GPU_TOKENS+=("$gpu_uuid")
done
declare -A seen_gpu_tokens=()
for gpu_token in "${GPU_TOKENS[@]}"; do
  if [[ -n "${seen_gpu_tokens[$gpu_token]:-}" ]]; then
    echo "SAFETY ABORT: duplicate CUDA worker token '$gpu_token'. No simulation was started." >&2
    exit 3
  fi
  seen_gpu_tokens[$gpu_token]=1
done

TASK_LINES=()
while IFS= read -r line; do TASK_LINES+=("$line"); done < <(tail -n +2 "$MANIFEST")
if (( ${#TASK_LINES[@]} < GPU_WORKERS )); then
  GPU_WORKERS=${#TASK_LINES[@]}
fi

prepare_work_dir() {
  local work_dir="$1"
  mkdir -p "$work_dir/results"
  local path base
  for path in "$PROJECT_DIR"/C5txt "$PROJECT_DIR"/C5Unity "$PROJECT_DIR"/Geometry_*; do
    [[ -e "$path" ]] || continue
    base="$(basename "$path")"
    [[ -e "$work_dir/$base" || -L "$work_dir/$base" ]] || ln -s "$path" "$work_dir/$base"
  done
}

run_task() {
  local line="$1" gpu_token="$2"
  local task_id case_index role reference target seed log_path work_dir
  IFS=$'\t' read -r task_id case_index role reference target seed log_path work_dir <<< "$line"
  prepare_work_dir "$work_dir"

  local -a args=(
    --critical-search-binary
    "--rod-depth=$reference"
    "--particles=$PARTICLES"
    "--inactive=$INACTIVE"
    "--active=$ACTIVE"
    "--seed=$seed"
  )
  if [[ "$role" == "direct_target" ]]; then
    args[1]="--rod-depth=$target"
  elif [[ "$role" == "cs" ]]; then
    args+=("--cs-depth=$target" "--cs-generations=$CS_GENERATIONS")
  fi

  echo "$(date '+%F %T') START task=$task_id case=$case_index role=$role d=$reference target=$target seed=$seed gpu=$gpu_token"
  printf 'RUNNING\n' > "$RUN_DIR_ABS/status/task_$(printf '%03d' "$task_id").status"
  (
    cd "$work_dir"
    export CUDA_VISIBLE_DEVICES="$gpu_token"
    "$EXE" "${args[@]}"
  ) > "$log_path" 2>&1
  local rc=$?
  if (( rc == 0 )); then
    printf 'DONE\n' > "$RUN_DIR_ABS/status/task_$(printf '%03d' "$task_id").status"
    echo "$(date '+%F %T') DONE  task=$task_id case=$case_index role=$role gpu=$gpu_token"
  else
    printf 'FAILED rc=%d\n' "$rc" > "$RUN_DIR_ABS/status/task_$(printf '%03d' "$task_id").status"
    echo "$(date '+%F %T') FAIL  task=$task_id case=$case_index role=$role rc=$rc log=$log_path" >&2
    return "$rc"
  fi
}

gpu_worker() {
  local worker_index="$1" gpu_token="$2" task_index
  for ((task_index = worker_index; task_index < ${#TASK_LINES[@]}; task_index += GPU_WORKERS)); do
    run_task "${TASK_LINES[$task_index]}" "$gpu_token" || return $?
  done
}

monitor_progress() {
  while true; do
    sleep "$MONITOR_INTERVAL_SEC"
    echo "---- progress $(date '+%F %T') ----"
    local line task_id case_index role reference target seed log_path work_dir last_cycle state
    for line in "${TASK_LINES[@]}"; do
      IFS=$'\t' read -r task_id case_index role reference target seed log_path work_dir <<< "$line"
      state="PENDING"
      [[ -f "$RUN_DIR_ABS/status/task_$(printf '%03d' "$task_id").status" ]] &&
        state="$(tr -d '\r\n' < "$RUN_DIR_ABS/status/task_$(printf '%03d' "$task_id").status")"
      last_cycle="-"
      if [[ -f "$log_path" ]]; then
        last_cycle="$(awk '/^Cycle [0-9]+,/ {value=$2; gsub(/,/, "", value)} END {if (value != "") print value; else print "-"}' "$log_path")"
      fi
      printf 'task=%02d case=%02d %-13s cycle=%s/%d %s\n' \
        "$task_id" "$case_index" "$role" "$last_cycle" "$((INACTIVE + ACTIVE))" "$state"
    done
  done
}

echo "CS validation job ${SLURM_JOB_ID:-manual}"
echo "  node: ${SLURM_JOB_NODELIST:-unknown}"
echo "  verified full node3 GPU allocation: $ALLOCATED_PHYSICAL_GPUS"
echo "  pinned GPU UUIDs for physical GPUs 0-5: ${GPU_TOKENS[*]}"
echo "  workers: $GPU_WORKERS"
echo "  reference depths: $REFERENCE_DEPTHS_CM cm"
echo "  step: $PROBE_STEP_CM cm"
echo "  histories/cycle: $PARTICLES"
echo "  cycles: $INACTIVE inactive + $ACTIVE active"
echo "  output: $RUN_DIR_ABS"
if command -v nvidia-smi >/dev/null 2>&1; then nvidia-smi -L; fi

monitor_progress &
monitor_pid=$!
worker_pids=()
for ((worker = 0; worker < GPU_WORKERS; ++worker)); do
  gpu_worker "$worker" "${GPU_TOKENS[$worker]}" &
  worker_pids+=("$!")
done

worker_failure=0
for pid in "${worker_pids[@]}"; do
  if ! wait "$pid"; then worker_failure=1; fi
done
kill "$monitor_pid" 2>/dev/null || true
wait "$monitor_pid" 2>/dev/null || true

if (( worker_failure != 0 )); then
  echo "At least one GPU worker failed. Check $RUN_DIR_ABS/status and logs." >&2
  exit 1
fi

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Python interpreter not found: $PYTHON_BIN" >&2
  exit 1
fi
"$PYTHON_BIN" scripts/analyze_cs_prediction_validation_slurm.py \
  --run-dir "$RUN_DIR_ABS" \
  --inactive "$INACTIVE" \
  --cs-generations "$CS_GENERATIONS" \
  --confidence-sigma "$CONFIDENCE_SIGMA" \
  --expected-sign "$EXPECTED_DELTA_K_SIGN" \
  --source-tail-abs-pcm "$SOURCE_TAIL_ABS_PCM" \
  --source-tail-relative-fraction "$SOURCE_TAIL_REL_FRACTION"

echo "Validation complete."
echo "  CSV:     $RUN_DIR_ABS/validation.csv"
echo "  summary: $RUN_DIR_ABS/summary.md"
echo "  plot:    $RUN_DIR_ABS/cs_validation.svg"
