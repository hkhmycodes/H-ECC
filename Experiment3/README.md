# EXP-3 — CPU arithmetic cost of pK vs Weierstrass–Jacobian

Two-stage prototype:

- **Stage A (validation):** every benchmarked input is verified against the
  Experiment-1 oracle in **native Python arithmetic**, before any timing.
- **Stage B (timing):** the same formulas are timed with either the pure
  Python field backend or the gmpy2/GMP backend, selected with `--field`.

## What is compared

On the 256-bit test vector of the manuscript's Section `sec:testvector`:

- **pK current** — the formulas of Sections `sec:add`–`sec:double`, exactly
  as validated in Experiment 1.
- **Weierstrass–Jacobian generic** — `add-2007-bl` and `dbl-2007-bl`
  (general `a`) on \(E_W : Y^2 = X^3 - \tfrac{1}{36}X\).
- **Weierstrass–Jacobian mixed** — `madd-2007-bl`, used where one operand
  has \(Z = 1\).

## Methodology

1. **Oracle-only input generation.** Every test point is obtained as
   \(P = \Phi_D^{-1}([k]\,\Phi_D(G))\), with the scalar multiplication
   performed by the Weierstrass backend, never by the pK side.

2. **Exceptional-pair rejection.** Each addition pair is rejected if
   \(P = Q\) or \(P = -Q\); rejections are counted by category.

3. **Four-way scalar validation.** For every scalar \(k\) used:
   \[
   W_{\rm oracle} = \text{oracle.scalar\_mul}(k, \Phi_D(G)),
   \]
   and the following three values must all equal \(W_{\rm oracle}\):
   \[
   \Phi_D(\mathrm{pK\_scalar}(k, G)),\qquad
   \text{WS-RL}(k, \Phi_D(G)),\qquad
   \text{WS-mixed}(k, \Phi_D(G)).
   \]
   The oracle used here is the split-cubic scalar multiplication from
   `pk_oracle.py`, not this package's own `ws_scalar_mul`.

4. **Separation of validation from timing.** Validation runs in native
   Python arithmetic; only after it passes are the input points converted
   to the active field backend (py or gmp) for the timing pass.

5. **Warm-up + repeated timing.** One untimed pass over each argument list,
   then median of `--repeats` passes (default 5) and `--scalar-repeats`
   (default 5) for scalar multiplication.

6. **Failure on insufficient samples.** If fewer than `--N` usable pairs or
   `--N-scalar` scalars are obtained, the run raises.

7. **\(T_M\) microbenchmark.** Median time for one field multiplication,
   reported separately and used to normalise primitive times.

## What is *not* claimed

- Stage B with `--field py` measures **Python** timings, and with
  `--field gmp` measures **GMP-backed Python** timings — i.e. GMP field
  arithmetic with Python dispatch on top of each call.
- `T_M` normalisation is meaningful within one backend and one implementation
  but is **not** portable across machines or implementations.
- In the timing path, `cmul` and small-integer multiplications execute as
  ordinary modular multiplications, so the timed costs do not reflect a
  hardware model with \(C < M\). Counted operations and wall-clock times
  are two independent metrics.
- The generic pK-vs-WS comparison uses the same right-to-left binary
  double-and-add strategy. The mixed-Jacobian baseline is a **separate
  left-to-right traversal** with `madd-2007-bl`, and is reported as a
  practical Weierstrass baseline, not as a formula-cost control.
- The C++/GMP port for a publication-quality absolute CPU table is deferred
  until the gmpy2 run shows whether the pK/WS ratios are backend-stable.

## Usage

```bash
cd Experiment3

# Pure Python (development)
python run_benchmark.py --field py \
    --seeds 0,1,2,3,4 --N 1000 --N-scalar 50 \
    --repeats 5 --scalar-repeats 5

# GMP-backed Python (prototype for the paper's CPU table)
pip install gmpy2
python run_benchmark.py --field gmp \
    --seeds 0,1,2,3,4 --N 1000 --N-scalar 50 \
    --repeats 5 --scalar-repeats 5
```

Outputs:

- `results/report.json` — per-seed reports plus an aggregate
- `results/report.txt`  — human-readable aggregate

## Dependencies

- Python 3.10+ (required).
- `gmpy2` for `--field gmp` (optional).
- Requires `Experiment1/` to be a sibling directory.