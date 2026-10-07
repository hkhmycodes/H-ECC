# EXP-2 / EXP-6 — Extraction and symbolic certification of the
# nine-dimensional (2,3)/(3,2) addition-law basis

This package produces two artifacts from a single computation:

- **EXP-6** — the extraction of the nine-dimensional \((2,3)\) addition-law
  space for the concrete 256-bit test vector of the manuscript.
- **EXP-2** — the independent symbolic certification of the extracted laws
  in SageMath.

The extraction is a custom Python pipeline that assembles a 756 × 162
matrix over \(\mathbb{F}_p\), computes its nullspace by exact RREF, and
outputs the nine basis laws. The certification is a separate SageMath
implementation that loads the coefficient vectors and verifies each law by
exact ideal membership modulo \((F(P), F(Q))\).

Extraction and certification do **not** share the reduction engine: the
Python extractor uses a hand-written memoized substitution, while the Sage
verifier uses Sage's polynomial ideal membership. Agreement of the two is
the certificate.

## Test vector

    p = 106839527430202782610735740077697524141755216002708025221579380756122101340723
    (u : v : w) = (5 : 1 : 7)

The curve is

    F(a, b, c) = 5a(b² − c²) + b(c² − a²) + 7c(a² − b²) = 0.

## Method

1. **Sanity checks before any RREF.** The script verifies, in order:
   - `F(P)` reduces to zero under the substitution
         a²b → 5ab² − 5ac² + bc² + 7a²c − 7b²c,
   - `F(Q)` reduces to zero under the same substitution on `(x, y, z)`,
   - each `S_i` of the composite group sum has bidegree exactly (4, 4)
     and is nonzero,
   - `F(S_0, S_1, S_2)` reduces to zero modulo `(F(P), F(Q))`,
   - the assembled matrix has shape exactly 756 × 162.

2. **Composite group sum.** Compute
       A = P·∇F(Q)      (bidegree (1, 2))
       B = Q·∇F(P)      (bidegree (2, 1))
       R = A·P − B·Q    (bidegree (2, 2))
       S = (R₁R₂, R₂R₀, R₀R₁)   (bidegree (4, 4))

3. **Candidate space.** Each coordinate Z_i is a linear combination of
   the 6 · 9 = 54 basis monomials μ·ν with μ ∈ Q_P (six quadratic
   monomials in (a,b,c)) and ν ∈ C_Q (nine reduced cubic monomials in
   (x,y,z)). This gives 3 · 54 = 162 unknowns.

4. **Linear system.** For each μν, reduce S_i · μν modulo `(F(P), F(Q))`
   to a 378-vector and fill matrix rows:
       G_1 = S_1 Z_0 − S_0 Z_1  (rows 0..377)
       G_2 = S_2 Z_0 − S_0 Z_2  (rows 378..755)

5. **Nullspace.** RREF over \(\mathbb{F}_p\); the free columns give the
   kernel basis.

6. **Certification.** Each law is certified by exact reduction modulo
   `(F(P), F(Q))`. Any monomial of unexpected degree raises an error
   instead of being silently dropped.

7. **(3,2) basis.** `A_{3,2}(P, Q) = A_{2,3}(Q, P)` re-expressed in the
   bidegree-(3,2) monomial basis, by explicit coefficient permutation.

## Usage

    cd Experiment2_6
    python exp6.py
    sage sage_verify.sage

The Python script produces the extraction; the Sage script produces the
independent certification. Both write into `results/`.

## Expected result

    matrix shape : 756 × 162
    rank         : 153
    nullity      : 9
    all certified: True

The nullity equals `dim A_{2,3} = 9` as predicted by the exact-dimension
proposition of the manuscript.

## Outputs

    results/
    ├── parameters.json         — p, monomial bases, matrix shape, rank, nullity
    ├── basis_23.json           — nine (2,3) laws as coefficient vectors
    ├── basis_32.json           — nine (3,2) laws as coefficient vectors
    ├── matrix_rank.txt         — shape, rank, nullity
    ├── verification_sage.txt   — output of the independent Sage verifier
    └── certificates/
        └── law_00.txt … law_08.txt

## Dependencies

- Python 3.10+ (no external packages for the extractor).
- SageMath (for the independent verifier).