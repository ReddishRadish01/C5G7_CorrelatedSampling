#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 )); then
  echo "Usage: $0 ROD_DEPTH_CM [extra program args...]" >&2
  echo "Example: $0 34.65 --particles=100000 --inactive=100 --active=100" >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

EXE="${EXE:-build/c5g7}"
CASE_ARGS="${CASE_ARGS:---critical-search-binary}"
DEPTH="$1"
shift

make "$EXE"

"$EXE" $CASE_ARGS --rod-depth="$DEPTH" "$@"
