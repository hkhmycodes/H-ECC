#!/usr/bin/env bash
# EXP-4 driver (rev 2).
#
# Invokes exp4_openmp ONE thread count per process, with a cooldown
# between invocations. This is the protocol that produced the tightest
# per-pass spreads (5-19%) on the test laptop, versus 30-100% for the
# "all T inside one process" variant.
#
# Usage:
#   ./run_exp4.sh [outdir]
#
# Override the grid with environment variables:
#   NS="1024 4096" SEEDS="0 1" THREADS="1 2 4" COOLDOWN=30 ./run_exp4.sh
#
# Default grid: N in {1024 4096}, T in {1 2 4 8}, seed 0, cooldown 60 s.
#
# Stops at the first nonzero exit code of exp4_openmp:
#   0 ok
#   1 usage
#   2 pre-timing validation failed
#   3 outputs differ across thread counts
#   4 pK/Jacobian oracle mismatch
#   5 OpenMP team size != requested T

set -u

OUT="${1:-results}"
NS="${NS:-1024 4096}"
THREADS="${THREADS:-1 2 4 8}"
SEEDS="${SEEDS:-0}"
COOLDOWN="${COOLDOWN:-60}"
SMALL_NS="${SMALL_NS:-}"          # set SMALL_NS="16 64 256" to include crossover

BIN="./exp4_openmp"
[ -x "$BIN" ] || BIN="./exp4_openmp.exe"
if [ ! -x "$BIN" ]; then
  echo "exp4_openmp not found; compile first:"
  echo "  g++ -std=c++17 -O3 -fopenmp -o exp4_openmp exp4_openmp.cpp -lgmpxx -lgmp"
  exit 1
fi

mkdir -p "$OUT"

# ------------------------------------------------------------------
# One invocation, one thread count, one JSON + one log.
# ------------------------------------------------------------------
run_one () {
  local N="$1" seed="$2" T="$3" tag="$4"
  local base="$OUT/exp4_${tag}N${N}_T${T}_s${seed}"
  printf '[%s] N=%-6s T=%-2s seed=%s -> %s.json\n' \
         "$(date +%H:%M:%S)" "$N" "$T" "$seed" "$base"

  "$BIN" "$N" "$seed" "$T" > "$base.json" 2> "$base.log"
  local rc=$?
  if [ $rc -ne 0 ]; then
    echo "  FAILED (exit $rc). See $base.log"
    exit $rc
  fi
  # Echo the summary line from stderr for the console record.
  grep -E '^\[exp4\] (pk|jac)' "$base.log" | sed 's/^/  /'
}

# ------------------------------------------------------------------
# Full grid: for each seed, for each N, for each T, run once.
# Cooldown between every invocation to let the CPU cool.
# ------------------------------------------------------------------
for seed in $SEEDS; do
  for N in $NS; do
    for T in $THREADS; do
      run_one "$N" "$seed" "$T" ""
      if [ "$COOLDOWN" -gt 0 ]; then
        sleep "$COOLDOWN"
      fi
    done
  done
done

# Optional small-N crossover runs at seed 0.
if [ -n "$SMALL_NS" ]; then
  for N in $SMALL_NS; do
    for T in $THREADS; do
      run_one "$N" 0 "$T" "small_"
      if [ "$COOLDOWN" -gt 0 ]; then
        sleep "$COOLDOWN"
      fi
    done
  done
fi

echo
echo "all runs finished OK -> $OUT"
echo "next: python3 check.py $OUT/exp4_*.json"