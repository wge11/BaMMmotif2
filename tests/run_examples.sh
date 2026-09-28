#!/usr/bin/env bash
# Run every BaMM!motif tool on the bundled example data.
#
# Usage:
#   tests/run_examples.sh BIN_DIR OUT_DIR [THREADS]
#
# Each test case writes into its own sub-directory of OUT_DIR, so two runs
# (e.g. before and after a code change) can be compared with
#   tests/compare_outputs.py OUT_DIR_A OUT_DIR_B
#
# THREADS defaults to 1, which makes all runs fully deterministic.

set -euo pipefail

if [[ $# -lt 2 ]]; then
    sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi

BIN=$(cd "$1" && pwd)
OUT=$2
THREADS=${3:-1}
EXAMPLE=$(cd "$(dirname "$0")/../example" && pwd)

SEQ="$EXAMPLE/JunD.fasta"
PWM="$EXAMPLE/PWM_peng10.meme"

mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)

FAILURES=0

run() {
    local name=$1; shift
    echo "== $name"
    mkdir -p "$OUT/$name"
    if ! "$@" > "$OUT/$name/stdout.log" 2>&1; then
        echo "   FAILED (see $OUT/$name/stdout.log)"
        FAILURES=$((FAILURES + 1))
    fi
}

run em_pwm        "$BIN/BaMMmotif" "$OUT/em_pwm" "$SEQ" --PWMFile "$PWM" --maxPWM 2 \
                  --EM --FDR --scoreSeqset --threads "$THREADS"

run cgs_pwm       "$BIN/BaMMmotif" "$OUT/cgs_pwm" "$SEQ" --PWMFile "$PWM" --maxPWM 1 \
                  --CGS --FDR --threads "$THREADS"

run em_sites_ss   "$BIN/BaMMmotif" "$OUT/em_sites_ss" "$SEQ" --ss -k 1 \
                  --bindingSiteFile "$EXAMPLE/bindingsites.block" --EM --scoreSeqset \
                  --threads "$THREADS"

run advance_em    "$BIN/BaMMmotif" "$OUT/advance_em" "$SEQ" --PWMFile "$PWM" --maxPWM 1 \
                  --EM --advanceEM --threads "$THREADS"

run score_bamm    "$BIN/BaMMmotif" "$OUT/score_bamm" "$SEQ" \
                  --BaMMFile "$EXAMPLE/JunD_motif_1.ihbcp" --bgModelFile "$EXAMPLE/JunD.hbcp" \
                  --scoreSeqset --threads "$THREADS"

run scan          "$BIN/BaMMScan" "$OUT/scan" "$SEQ" --PWMFile "$PWM" --maxPWM 2

run fdr           "$BIN/FDR" "$OUT/fdr" "$SEQ" --PWMFile "$PWM" --maxPWM 1 --EM \
                  --threads "$THREADS"

run extract_probs "$BIN/extractProbs" "$OUT/extract_probs" \
                  "$EXAMPLE/JunD_motif_1.ihbcp" "$EXAMPLE/JunD.hbcp"

if [[ $FAILURES -gt 0 ]]; then
    echo "$FAILURES example run(s) failed: $OUT"
    exit 1
fi
echo "All example runs finished: $OUT"
