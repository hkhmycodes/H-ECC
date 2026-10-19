#!/usr/bin/env bash
# run_check2.sh -- MSYS2 orchestration for the pK affine/cached cross-check
set -euo pipefail
cd "$(dirname "$0")"

# -----------------------------------------------------------------------------
# 0. Environment guard
# -----------------------------------------------------------------------------
[[ -n "${MSYSTEM:-}" ]] || {
  echo "error: not in an MSYS2 shell.  Open 'MSYS2 UCRT64' from Start." >&2
  exit 1
}
command -v gcc >/dev/null || {
  echo "error: gcc not found.  Run:  pacman -S mingw-w64-ucrt-x86_64-gcc" >&2
  exit 1
}
echo "[env]  MSYSTEM=$MSYSTEM"
echo "[env]  $(gcc --version | head -n1)"
echo "[env]  $(nproc 2>/dev/null || echo '?') threads"

# -----------------------------------------------------------------------------
# 1. Inputs
# -----------------------------------------------------------------------------
for f in bench.c check2.c r2.h aff.h; do
  [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
done

# -----------------------------------------------------------------------------
# 2. Build
#     -O3 -march=native mirrors the paper's setup (gcc -O3 -march=native)
#     -fno-strict-aliasing : the bit-twiddling in fe* overlaps unions
#     -fno-plt             : keeps __uint128_t mul helpers in-register
# -----------------------------------------------------------------------------
CFLAGS=(-O3 -march=native -mtune=native -std=gnu11
        -fno-strict-aliasing -fno-plt
        -Wall -Wextra -Wno-unused-parameter -Wno-unused-function)

echo "[build] check2.exe"
gcc "${CFLAGS[@]}" -o check2.exe check2.c

# only build bench.exe if bench.c exposes its own main()
if grep -qE '^[[:space:]]*int[[:space:]]+main[[:space:]]*\(' bench.c; then
  echo "[build] bench.exe"
  gcc "${CFLAGS[@]}" -o bench.exe bench.c
fi

# -----------------------------------------------------------------------------
# 3. Run the cross-check
# -----------------------------------------------------------------------------
echo "[run]   check2.exe"
out=$(./check2.exe)
printf '%s\n' "$out"

u0=$(sed -n 's/^u0=//p' <<<"$out")
u1=$(sed -n 's/^u1=//p' <<<"$out")
u2=$(sed -n 's/^u2=//p' <<<"$out")
m0=$(sed -n 's/^m0=//p' <<<"$out")
m1=$(sed -n 's/^m1=//p' <<<"$out")
m2=$(sed -n 's/^m2=//p' <<<"$out")

if [[ "$u0" == "$m0" && "$u1" == "$m1" && "$u2" == "$m2" ]]; then
  echo "[OK]    affine path  ==  cached path  (all three coords)"
else
  echo "[FAIL]  paths diverge:" >&2
  diff <(printf 'u0=%s\nu1=%s\nu2=%s\n' "$u0" "$u1" "$u2") \
       <(printf 'm0=%s\nm1=%s\nm2=%s\n' "$m0" "$m1" "$m2") >&2 || true
  exit 2
fi

# -----------------------------------------------------------------------------
# 4. Optional timed benchmark (bench.exe, if built)
# -----------------------------------------------------------------------------
if [[ -x bench.exe ]]; then
  echo
  echo "[run]   bench.exe  (timed benchmark)"
  ./bench.exe
fi