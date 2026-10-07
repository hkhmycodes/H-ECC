# Experiment 8 — Complete (2,3) Addition-Law Benchmark

Benchmark of the seven complete (2,3) addition-law triples identified by
Experiment 7, against the current composite pK formulas and against
Weierstrass--Jacobian baselines on the birationally equivalent curve.

This directory accompanies Section 10 (Experimental program, Stage II) of
the manuscript "A Projective Framework for Inversion-Free Elliptic Curve
Cryptography on Isogonal Pivotal Cubics".


## 1. Scientific question

Experiment 7 established that the pK (2,3)-addition-law space

  - is exactly nine-dimensional (dim A_{2,3} = 9),
  - requires THREE individual laws for a geometrically complete system
    (r_min = 3), by an intersection-number obstruction
    (c_1(M_{2,3})^2 = 18) combined with Lange-Ruppert base-point-freeness,
  - and has SEVEN complete basis triples among the C(9,3) = 84 subsets of
    the extracted nullspace basis.

Experiment 8 answers the follow-up question:

  Does the theoretically correct complete (2,3) representation actually
  reduce the arithmetic cost of pK addition and doubling, compared with
  the present composite pK formula P (+) Q = R* ?

Both outcomes are scientifically acceptable:

  - If faster, the paper's emphasis shifts to the (2,3) formulation as the
    recommended pK arithmetic backend.
  - If slower, the paper retains a complete correctness, extraction,
    minimality, and cost study of the present pK representation; the result
    motivates the small-coefficient (LLL-reduced) basis search as future
    work.


## 2. Inputs and outputs

### Inputs (must exist before running)

  ../Experiment2_6/results/basis_23.json   produced by EXP-2 / EXP-6
      Nine basis laws as 162-vectors; Q_P and C_Q exponent lists.

  ../Experiment7/results/report.json       produced by EXP-7
      Completion flags for every L1_r3 subset; the seven complete == true
      entries define the triples to benchmark.

basis_23.json layout (assumed by the loader):

  {
    "Q_P": [[2,0,0], [1,1,0], [1,0,1], [0,2,0], [0,1,1], [0,0,2]],
    "C_Q": [[0,0,3], [0,1,2], [0,2,1], [0,3,0],
            [1,0,2], [1,1,1], [1,2,0], [2,0,1], [3,0,0]],
    "laws": [ {"coeffs": [162 ints]}, ... ]
  }

The coefficient layout of each 162-vector is

  idx = coord_idx * 54 + mu * 9 + nu

where coord_idx in {0,1,2} selects the projective coordinate, mu indexes
Q_P, and nu indexes C_Q. This layout is asserted at load time against the
hard-coded exponent lists.

### Outputs (written to results/)

  exp8_verify.json     produced by verify
  exp8_profile.json    produced by profile
  exp8_analysis.json   produced by analyze
  exp8_bench.json      produced by bench

All outputs are written atomically (temp file + os.replace), so a killed
run leaves the previous file intact.


## 3. Files in this directory

  Experiment8/
    README.md          this file
    exp8.py            the experiment script (pure Python 3)
    results/           created automatically on first write
      exp8_verify.json
      exp8_profile.json
      exp8_analysis.json
      exp8_bench.json

There is no .sage file: the whole experiment is runnable with plain
Python 3, because Sage's preparser would turn every integer literal into
a Sage Integer and make the arithmetic incomparable with the pure-Python
numbers of EXP-3.


## 4. Requirements

  - Python 3.8+ (the version is recorded in exp8_bench.json)
  - No mandatory third-party packages
  - Optional: gmpy2 for the second arithmetic backend. If gmpy2 is not
    installed, the gmpy2 backend is silently skipped and the int backend
    results are still complete.

SageMath is NOT required to run EXP-8. The only Sage-side work referenced
by the manuscript (LLL reduction of the nullspace basis; symbolic
reduction of the diagonal modulo F(P)) is left as future work and is not
invoked here.


## 5. How to run

Run the phases in this order. The bench phase refuses to run unless verify
has produced a clean exp8_verify.json.

  cd ~/Documents/HGECC/Experiment8

  python3 exp8.py verify      correctness on on-curve and inverse pairs
  python3 exp8.py profile     coefficient profile of the nine basis laws
  python3 exp8.py analyze     static CSE-aware cost, both evaluation orders
  python3 exp8.py bench       wall-clock timing, int and (optional) gmpy2
  python3 exp8.py summary     aggregate report

Expected runtimes (single machine, pure-Python backend):

  verify    under 1 minute
  profile   under 10 seconds
  analyze   under 10 seconds (three codegens per triple, no algebra)
  bench     several minutes, dominated by 5 seeds x 1000 inputs x 5 passes
            x (2 baselines + 2 backends + 7 triples)
  summary   under 1 second

The benchmark does NOT use any symbolic algebra. It is a straight-line
evaluation timing, not a Groebner computation.


## 6. What each phase does

### verify - correctness gate

This is the only phase that can reject a triple. It must exit with code 0
before bench will run. It performs four independent checks, in order:

1. Generator and pivot sanity. Confirms F(G) = 0, that the image of the
   pK generator lies on E_W : y^2 = x^3 - x/36, and that the pivot image
   W_D lies on E_W. This is the earliest place the file can abort, and
   the message names the failing object.

2. Oracle versus composite. Runs 30 addition comparisons and 15 doubling
   comparisons of the independent E_W oracle oracle_sum_pk against
   composite_pk_add / composite_pk_double. EXP-1 already established that
   the composite formula agrees with its own Phi_D oracle, so any mismatch
   here points at the E_W map code, not at the triples.

3. Triple-level checks. For each of the seven complete triples:

     - Addition on random pairs: at least MIN_ADD_OK = 30 successful
       comparisons against oracle_sum_pk.
     - Addition on inverse pairs: at least MIN_INV_OK = 10 comparisons
       against D = (5:1:7). Inverse pairs are constructed by
       phi_0(Q) = 2*W_D - phi_0(P).
     - Doubling: at least MIN_DBL_OK = 15 comparisons against
       oracle_sum_pk(P, P) on the diagonal codegen.

   Every check verifies ALL THREE laws of the triple (not only the one
   returned by ct_select_nonzero), and rejects the zero vector as a
   spurious match. Sample points come from two sources: multiples of the
   generator image on E_W, and random_curve_point from the full E(F_p)
   via the a = 1 quadratic chart.

4. Exit code. Returns 0 only if failures == 0 and every threshold is met.
   Any failure increments failures and is reported on stdout.

### profile - coefficient profile

Reports, for each of the nine basis laws and each coordinate:

  - number of nonzero coefficients,
  - number of +1 and -1 entries (these are free in codegen),
  - number of "small" entries (1 < |c| < 2^16, weighted 0.10),
  - number of "full" entries (|c| >= 2^16, weighted 1.00),
  - min/median/max bit length of the signed representative.

The totals line is the decisive output:

  Totals:  +1 <N>   -1 <N>   small <N>   full <N>

  - If full dominates, the current arbitrary nullspace basis is dense in
    coefficients, and the small-coefficient LLL basis search is the
    natural next step.
  - If small or +/-1 dominate, the current basis is already sparse and
    the static cost model will reflect that.

### analyze - static CSE-aware cost

Computes three static costs per triple:

  - Order A (product-first): shares all distinct P_mu * Q_nu products
    across the nine coordinates, then accumulates coefficient multiplies
    and additions per coordinate.
  - Order B (coefficient-first): for each coordinate and each mu group,
    builds W = sum_nu c_{mu,nu} Q_nu, then sums Pm[mu] * W over mu. This
    avoids materialising the shared products.
  - Diagonal (Q = P): derives the degree-5 collapse from the codegen
    itself, so the static count matches what bench will time. The raw
    monomials are not reduced modulo F(P); the reduction is a Sage-side
    step listed as future work.

Each static cost is scored by a configurable heuristic (HEUR_WEIGHTS at
the top of the file). The defaults are:

  {"M_general": 1.00, "C_full": 1.00, "C_small": 0.10,
   "S": 0.80, "A": 0.05}

### bench - wall-clock timing

Runs, for each of five seeds (SEEDS = [0,1,2,3,4]):

  - N_WARMUP = 100 untimed calls to warm the caches,
  - N_PASSES = 5 timed passes over N_TRIALS = 1000 inputs,
  - the reported figure is the median across the five seeds of the median
    across the five passes.

For each triple the phase measures:

  - t_add: full addition - evaluate all three laws and select a nonzero
    output via ct_select_nonzero (the cryptographically meaningful
    number),
  - t_dbl: full doubling - same but on the diagonal codegen,
  - t_sel: selection alone, timed on precomputed triple outputs, so the
    completeness-handling cost is separable from the polynomial cost.

Baselines measured in the same loop, on the same random inputs:

  - Composite pK addition and doubling (Section 4 of the manuscript),
  - Jacobian addition (add-2007-bl) and doubling (dbl-2007-bl, general a).

The int backend runs on Python int. The optional gmpy2 backend runs on
gmpy2.mpz, with the composite and Jacobian baselines rewritten directly
on mpz (no conversion inside the timed loop), and all coefficient
literals pre-bound as mpz in the exec namespace.

### summary - aggregate report

Prints the static rankings (A, B, diagonal) and the wall-clock rankings,
and reports the composite and Jacobian baselines side by side.


## 7. Timing methodology (matching EXP-3)

The benchmark deliberately reuses EXP-3's methodology so the two are
directly comparable:

  - Same 256-bit prime p and same q (from the reproducibility vector of
    Section 6.5 of the manuscript).
  - Same composite pK baselines and same Jacobian baselines on the
    birationally equivalent E_W : y^2 = x^3 - x/36.
  - Same warm-up / timed-pass / median structure.
  - Five independent seeds for the scalar loop, matching EXP-3.
  - Garbage collection is disabled around every timed pass
    (gc.disable() / gc.enable()); the flag is recorded in
    exp8_bench.json under gc_disabled_during_timing.

Known, documented biases. Two small timing biases push the numbers
AGAINST the (2,3) side:

  1. The triple path wraps each call in an extra Python function and then
     calls ct_select_nonzero. The baselines are called directly.
  2. ct_select_nonzero performs % P_MOD with a Python int on every
     selected coordinate, even in the gmpy2 backend.

Both are small and are stated on stdout at the start of bench. They do
not change the qualitative conclusion.


## 8. Interpreting the results

The headline ratios the paper uses are

  T_pK_current / T_pK_(2,3)       improvement of the new representation
  T_pK_(2,3)   / T_Jacobian       absolute standing

Read them together with the coefficient profile:

  Profile                    Expected outcome              Next step
  -------------------------------------------------------------------
  Mostly +1 / -1 / small     (2,3) competitive             Enable fixed-
                             with composite                base variant
  Mostly full 256-bit        (2,3) loses to composite      Run LLL
                                                           basis search
  Mixed                      Depends on n_products         Compare both
                             and C_full                    orders and
                                                           both backends

Regardless of the outcome, the paper's structural results
(dim A_{2,3} = 9, r_min = 3, the seven explicit complete triples) are
unaffected. Only the arithmetic-cost comparison changes.


## 9. Deliberately not implemented

The following are out of scope for EXP-8's first run, and are labelled as
future work in the manuscript:

  - Fixed-base mixed variant. For scalar multiplication with a fixed
    affine base G (z = 1), the cubic monomials in Q and the products
    P_mu * Q_nu become precomputed constants, changing the cost profile
    substantially. This is the realistic scalar-multiplication scenario
    and is intended as a follow-up to EXP-8.

  - LLL / small-coefficient basis. Replacing the arbitrary nullspace
    basis with an LLL-reduced basis over Q (denominators cleared),
    re-certifying the few chosen triples with the EXP-7 Groebner test.
    Recommended only if the profile phase shows predominantly full-size
    coefficients.

  - Reduction of the degree-5 diagonal modulo F(P). The static diagonal
    cost in analyze is the raw degree-5 collapse; the reduction to a
    canonical 15-monomial form requires polynomial division in Sage.

  - Scalar multiplication benchmarking. Deferred until the addition and
    doubling numbers are stable and understood.

  - Version A (single-law benchmark). Skipped; the cost of completeness
    shows up directly in the triple timings and in t_sel.


## 10. Reproducibility checklist

A publication-quality EXP-8 run should record:

  [ ] python3 --version
  [ ] CPU model and frequency (lscpu / /proc/cpuinfo)
  [ ] Whether gmpy2 is installed; if so, its version
  [ ] Whether the CPU was pinned during timing (e.g. taskset)
  [ ] SHA-256 of exp8.py and of basis_23.json
  [ ] SHA-256 of report.json (the EXP-7 frozen snapshot)
  [ ] The exp8_verify.json file (must show "failures": 0)
  [ ] The exp8_profile.json totals line
  [ ] All four exp8_*.json files, unedited

The EXP-7 snapshot should be frozen before EXP-8 is run; EXP-8 reads but
never writes to the EXP-7 directory.


## 11. Failure modes and what they mean

  Q_P order mismatch / C_Q order mismatch at load
      basis_23.json stores a different monomial order than the hard-coded
      layout in exp8.py; adjust the asserts to the stored order.

  generator image fails the Weierstrass equation
      phi0_W or _SHIFT sign is off; check against the worked example
      D -> (1/3, -1/6) in the manuscript.

  oracle disagrees with composite add
      Bug in oracle_sum_pk or _W_D_on_E_W, not in the triples; the
      triples have not yet been tested at that point.

  add MISMATCH vs oracle for one triple only
      Decoding bug in decode_coordinate or a genuinely incomplete triple
      (should not happen for the seven EXP-7 triples).

  add MISMATCH vs oracle for all triples
      Layout bug: the 162-vector is being read with the wrong mu/nu
      strides.

  triple output is the zero vector
      Codegen bug: the parenthesised % m was lost, or a coefficient was
      emitted with the wrong sign.

  bench exits with "verify has not passed"
      exp8_verify.json is missing or has "failures" != 0; run verify
      first and read the failure log.

  gmpy2 not installed
      Expected; the int backend still runs. Install gmpy2 if the second
      backend is wanted.

  Verbatim "verbose 0 (3722: multi_polynomial_ideal.py...)" lines
      Not from EXP-8; EXP-8 does not invoke any Groebner computation.


## 12. Relation to the manuscript

  Manuscript location          What EXP-8 contributes
  ---------------------------------------------------------------------
  Section 10, EXP-7            Confirms the frozen minimality result
                               r_min = 3 and the seven complete triples
  Section 10, EXP-8            Static CSE-aware cost, wall-clock timing,
  (planned -> executed)        coefficient profile
  Section 10, tab:variants     Supplies the measured numbers for the
                               pK-vs-Jacobian row
  Section 10, Exceptional      Confirms that the complete triple handles
  states                       the base-locus cases the composite cannot
  Section 12, Conclusion       Provides the arithmetic-cost comparison
                               that decides whether the (2,3) route is
                               the recommended pK backend

If profile shows predominantly full-size coefficients, the honest
conclusion to write into the manuscript is:

  The (2,3) addition-law space is structurally correct, complete, and
  minimal at r_min = 3, but on the arbitrary nullspace basis extracted
  by EXP-6 its arithmetic cost is dominated by full-width constant
  multiplications. Closing the gap to the composite pK formula requires
  a small-coefficient basis (LLL over Q, re-certified by EXP-7's
  Groebner test), which is left as future work.

If profile shows mostly +1 / -1 / small coefficients, the honest
conclusion is the opposite, and the paper's emphasis shifts to the (2,3)
realisation.


## 13. License and provenance

Part of the manuscript's reproducibility package. The experiment is
designed to be re-run end to end on a fresh checkout given the EXP-2/6
and EXP-7 result files. No global state is mutated outside results/.