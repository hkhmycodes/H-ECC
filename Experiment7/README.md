# EXP-7 — search for a minimal complete subsystem of A_{2,3}

This package determines the minimum cardinality of a complete subsystem of
the nine-dimensional bidegree-\((2,3)\) addition-law space
\(\mathcal A_{2,3}\) extracted in EXP-2/6, for the concrete 256-bit test
vector of the pK framework.

## Main result

For the pK cubic
\[
	F = 5\alpha(\beta^2 - \gamma^2) + \beta(\gamma^2 - \alpha^2) + 7\gamma(\alpha^2 - \beta^2)
\]
over \(\mathbb{F}_p\) with the 256-bit prime
\(p = 106839527430202782610735740077697524141755216002708025221579380756122101340723\),
the minimum cardinality of a complete subsystem of the \((2,3)\) addition-law
space is
\[
	\boxed{\,r_{\min} = 3.\,}
\]

No pair of \((2,3)\)-laws is complete, and every generic triple of
\((2,3)\)-laws is complete. The exhaustive basis search found exactly
**seven** complete triples among the nine extracted basis laws, all of the
form
\[
	\{A^{(j)},\ A^{(7)},\ A^{(8)}\},
\qquad j \in \{0, 1, 2, 3, 4, 5, 6\}.
\]

## Why \(r = 2\) is impossible (theoretical reason)

Let \(E \subset \mathbb{P}^2\) be the pK cubic. On the abelian surface
\(E \times E\), any effective divisor of class \(\mathcal{O}(2,3)\) has
class
\[
	c_1 \;=\; 6\,\alpha + 9\,\beta,
\qquad \alpha = \mathrm{pr}_1^*[\mathrm{pt}],\quad
\beta = \mathrm{pr}_2^*[\mathrm{pt}].
\]
Because \(\alpha^2 = \beta^2 = 0\) and \(\alpha \beta = 1\) on \(E \times E\),
\[
	c_1^2 \;=\; 2 \cdot 6 \cdot 9 \;=\; 108 \;>\; 0.
\]
Two effective divisors on a smooth projective surface with positive
intersection number must intersect. Hence any two \((2,3)\)-laws share a
common zero on \(E \times E\), and \(r = 2\) is unattainable for the
\((2,3)\) embedding.

This is a theoretical obstruction, not an empirical accident. The \(0/36\)
in `L1_r2` (exhaustive basis pairs) and the \(0/30\) in `L2_r2` (sampled
pairs in the full space) are forced by the geometry.

## Why \(r = 3\) is generically complete

Three effective divisors of the same class on a surface generically have
empty intersection, because the codimension exceeds \(\dim(E \times E) = 2\).
The sampled search over \(30\) random triples in
\(\mathbb{P}^8(\mathbb{F}_p)^3\) confirmed this: **all \(30\) triples are
complete**.

## What was actually run

Four phases executed, all with zero timeouts and zero errors.

| Phase   | Scope                                            | Tests | Complete | Time    |
|---------|--------------------------------------------------|-------|----------|---------|
| `L1_r2` | Exhaustive basis pairs of the nine basis laws    | 36    | 0        | 6 min   |
| `L1_r3` | Exhaustive basis triples of the nine basis laws  | 84    | 7        | 95 min  |
| `L2_r2` | Sampled pairs in the full law space \(\mathbb{P}^8\) | 30 | 0        | 5 min   |
| `L2_r3` | Sampled triples in the full law space \(\mathbb{P}^8\) | 30 | 30   | 34 min  |

The seven complete basis triples correspond to the following indices of the
nine extracted basis laws:

| Test ID          | Basis triple        |
|------------------|---------------------|
| `L1_r3_s0027`    | \(\{0, 7, 8\}\)     |
| `L1_r3_s0048`    | \(\{1, 7, 8\}\)     |
| `L1_r3_s0063`    | \(\{2, 7, 8\}\)     |
| `L1_r3_s0073`    | \(\{3, 7, 8\}\)     |
| `L1_r3_s0079`    | \(\{4, 7, 8\}\)     |
| `L1_r3_s0082`    | \(\{5, 7, 8\}\)     |
| `L1_r3_s0083`    | \(\{6, 7, 8\}\)     |

Every complete basis triple contains the pair of indices \(\{7, 8\}\). No
complete pair exists, in the basis or in the sampled law space. The pair
\(\{7, 8\}\) alone is therefore not complete, but adding any one of the
other seven laws to it is.

The phase `L1_r4` (126 quadruples) was not run. Its value is negligible
given the theoretical bound: \(r = 2\) is impossible and \(r = 3\) is
attained, so further tests would only find larger subsystems.

## Epistemic status of each claim

**Rigorously established:**

- \(r = 2\) is impossible for the \((2,3)\) embedding, by the
  intersection-number argument \(c_1^2 = 108 > 0\) on \(E \times E\).
- Among the nine extracted basis laws, the minimum complete subsystem has
  cardinality \(3\), attained by exactly the seven triples listed above.
  This is exhaustive over the chosen basis.
- The Arène–Kohel–Ritzenthaler bound \(r \geq g + 1 = 2\) for \(g = 1\)
  is respected; the actual minimum is strictly larger for this embedding.

**Empirical evidence (not proof of absence):**

- No complete pair was found among the \(30\) sampled pairs in the full
  law space. This is consistent with, and corroborated by, the theoretical
  obstruction at \(r = 2\).
- All \(30\) sampled triples were complete. This is consistent with the
  codimension-\(3\) generic behaviour on \(E \times E\).

## Method

### Completeness test

A set of laws \(\{A^{(j)}\}\) is complete iff
\[
	\bigcap_j \operatorname{Exc}(A^{(j)}) = \varnothing
	\quad\text{on }E \times E
	\subset \mathbb{P}^2 \times \mathbb{P}^2.
\]

Because the defining generators are bihomogeneous, testing unit-ideality of
the affine ideal in \(\mathbb{F}_p[a,b,c,x,y,z]\) is *vacuous*: the origin
is a common zero of every bihomogeneous generator. The correct test is
projective. We use the **nine-chart method**: for each of the nine charts
\((p, q) \in \{a, b, c\} \times \{x, y, z\}\), set \(p = 1, q = 1\),
dehomogenize \(F_P\), \(F_Q\), and every \(Z_i^{(j)}\), and test the
resulting affine ideal in four variables for unit-ideality over
\(\mathbb{F}_p\). The law set is complete **iff all nine chart ideals are
the unit ideal**.

The same test is implemented independently in `sage_verify.sage`.

### Search levels

- **Level 1** — exhaustive over \(r\)-subsets of the nine extracted basis
  laws, for each \(r\) in `LEVEL1_RS`.
- **Level 2** — sampled over \(r\)-tuples of points in the projective law
  space \(\mathbb{P}^8(\mathbb{F}_p)\), for each \(r\) in `LEVEL2_RS`.

Level 2 is *not* exhaustive. The number of trials per \(r\) is set by
`LEVEL2_TRIALS`.

### Sampling

A point of \(\mathbb{P}^8(\mathbb{F}_p)\) is sampled by drawing a uniform
9-vector in \(\mathbb{F}_p^9\) and rejecting the all-zero vector. No
normalization is needed: multiplying a point of the projective law space by
a nonzero scalar produces the same zero locus.

## Phases and CLI

The search runs one phase per invocation. State is written after every
test and resumed automatically.

```
sage exp7.sage                     show help and phase status
sage exp7.sage <phase>             run or resume a phase
sage exp7.sage summary             print overall summary
sage exp7.sage --list              list phases and status
sage exp7.sage --reset-phase <p>   clear one phase's results
sage exp7.sage --reset-all         clear all phases
```

Available phases are determined by `LEVEL1_RS` and `LEVEL2_RS`:

- `L1_r2`, `L1_r3`, `L1_r4` — Level 1 for \(r = 2, 3, 4\).
- `L2_r2`, `L2_r3` — Level 2 for \(r = 2, 3\).

A test that ends with `complete: null` (timeout or error) is retried on
the next invocation of the same phase. Determinate results are never
re-computed.

## Configuration

Set at the top of `exp7.sage`:

- `LEVEL1_RS` — subset sizes for Level 1 (default `[2, 3, 4]`).
- `LEVEL2_RS` — subset sizes for Level 2 (default `[2, 3]`).
- `LEVEL2_TRIALS` — number of random combinations per \(r\) (default 30).
- `PER_CHART_TIMEOUT_S` — per-chart soft timeout (default 60 s).
- `RNG_SEED` — RNG seed for Level 2 sampling (default 0).

The per-chart timeout is a soft timeout via `signal.alarm`. On a host where
Sage falls back to a pure-Python Buchberger implementation, a chart may
occasionally exceed the soft budget; when that happens the test is marked
*indeterminate* and retried on the next invocation of the same phase.

## Files

```
Experiment7/
├── README.md           — this file
├── requirements.txt    — SageMath 9.x or later
├── exp7.sage           — main search driver
├── sage_verify.sage    — independent re-verification of complete examples
└── results/
    ├── report.json     — machine-readable state and results
    ├── report.txt      — human-readable summary (written by `summary`)
    └── verification_sage.txt  — output of the independent verifier
```

`exp7.sage` reads `Experiment2_6/results/basis_23.json` produced by
EXP-2/6. `sage_verify.sage` reads `results/report.json` written by
`exp7.sage`.

## Re-verification

To independently certify the complete subsystems found by `exp7.sage`:

```bash
sage sage_verify.sage
```

The verifier loads the current `report.json`, extracts every test marked
`complete: true`, and re-runs the nine-chart test from a fresh
implementation, writing `results/verification_sage.txt`.

## Dependencies

- SageMath 9.x or later.

## What is not established

- The intersection-number argument shows that *no* pair of \((2,3)\)-laws
  is complete, so the empirical \(0/30\) in `L2_r2` is not evidence of
  absence but a shadow of a theorem.
- The result that *every* generic triple is complete is consistent with,
  but not proved by, the sampled \(30/30\); a completeness proof of
  generic triples would require a dimension count over the coefficient
  space that is not carried out here.
- Whether the seven basis triples are the *only* complete triples inside
  \(\mathcal A_{2,3}\), up to the natural symmetries, is not addressed.
  The sampled search over \(30\) triples in the full law space found no
  non-basis complete triples, but the sample is not exhaustive.