#!/usr/bin/env bash
# Run all tools on the example data and compare the results with the
# reference outputs in tests/reference.
#
# Usage:
#   tests/check_examples.sh BIN_DIR [THREADS]
#
# Numeric values may differ within a small tolerance (e.g. between
# compilers); anything else counts as a regression.

set -euo pipefail

if [[ $# -lt 1 ]]; then
    sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi

HERE=$(cd "$(dirname "$0")" && pwd)
OUT=$(mktemp -d "${TMPDIR:-/tmp}/bamm-check.XXXXXX")
trap 'rm -rf "$OUT"' EXIT

"$HERE/run_examples.sh" "$1" "$OUT" "${2:-2}"
python3 "$HERE/compare_outputs.py" "$HERE/reference" "$OUT"
