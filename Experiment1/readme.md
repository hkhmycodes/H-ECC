## Methodology (matching the manuscript)

- **Uniform sampling.** Every random trial samples a scalar uniformly from
  `[0, q)` — the same convention as the manuscript.
- **Successful-trials loop.** Each step runs until `N` *successful* trials have
  been recorded, with a hard cap of `10N` attempts. A step reports PASS only
  if `successful >= N` **and** `failures == 0`; a step that reaches the cap
  with fewer successes is FAIL even without any explicit mismatch.
- **Per-step counters.** Every step reports `n_requested`, `n_attempts`,
  `n_trials` (successes), `n_skipped`, `n_failed`.
- **Test points generated on the oracle.** Every random `P` is
  `Φ_D⁻¹([k]·Φ_D(G))`, never produced by the pK implementation under test.
- **A, B, C constructions in Step 12.** For each base-locus point `X ∈ {A,B,C}`,
  the harness constructs `N` pairs `(P,Q)` whose secant line meets `E` at `X`,
  verifies `P,Q ∈ E`, `X,P,Q` collinear, `P,Q ≠ X`, and records whether
  `curve.add(P,Q)` raised the expected undefined-map condition.
- **D\* handled explicitly.** `pk_core.neg(D*) = D*` (2-torsion), so
  `pk_neg` is defined at `D*` and is tested deterministically in Step 7.

## Suggested commands

Smoke test:

    python run_experiment.py --N-add 10 --N-double 10 --N-neg 10 \
                             --N-assoc 5 --N-step-45 10 --N-step-12 5 \
                             --seed 0

Full run:

    python run_experiment.py --N-add 500 --N-double 500 --N-neg 500 \
                             --N-assoc 200 --N-step-45 500 --N-step-12 50 \
                             --seed 0

Robustness:

    python run_experiment.py ... --seed 1
    python run_experiment.py ... --seed 2

## What success looks like

For every step with random trials:

    n_requested == n_trials   (successful)
    n_failed   == 0
    n_skipped  is small and, ideally, 0

The A/B/C construction table should show

    A: n_undefined_as_expected == n_constructions
    B: n_undefined_as_expected == n_constructions
    C: n_undefined_as_expected == n_constructions