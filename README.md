# H-ECC: Pivotal Isogonal Cubics as Elliptic Curves

Code, data and certificates accompanying the manuscript

**"Pivotal Isogonal Cubics as Elliptic Curves: Group Law, Addition-Law Spaces, and the Cost of the Native Plane Embedding"**, Hussein Khayou (2026).

A pivotal isogonal (pK) cubic with pivot `D` is a plane elliptic curve whose group law is built from the collinearity of `D`, a point `P` and its isogonal conjugate `P*`. The article studies arithmetic on this native plane model: its group-law structure, its addition-law spaces, and what the native embedding costs compared with Jacobian coordinates.

The result is structural, not a speed-up. In a common C backend the best pK formulas we found are about 2x slower than Jacobian coordinates. The pK model is not a new hardness assumption: discrete-logarithm security is inherited from the Weierstrass model through the oracle isomorphism of the article.

## Repository layout

```
H-ECC/
├── pk_opt/               # Article: optimized inversion-free formulas, operation counts,
│                         #   I-centred law, fixed-point translation test, law-space row analysis,
│                         #   C benchmark pK vs. Jacobian
│                         # Cross-platform reproduction of the benchmark table
│                         #   (bench.c, r2.h, aff.h, check2.c, pk_opt_formulas.py, run_check2.sh)
├── Experiment1/          # EXP-1: correctness of pK primitives against the Weierstrass oracle
├── Experiment2_6/        # EXP-2 / EXP-6: extraction and certification of the (2,3)/(3,2) basis
├── Experiment3/          # EXP-3: earlier CPU cost study (unoptimized composite pK vs. Jacobian)
├── Experiment4/          # EXP-4: earlier CPU parallel-scaling study (OpenMP)
├── Experiment7/          # EXP-7: minimal complete subcollections of A_{2,3}
├── Experiment8/          # EXP-8: static cost analysis of the (2,3) triple (+ earlier Python timing)
├── ecpp_certificates.gp  # ECPP primality certificates for p and q (PARI/GP)
├── LICENSE               # MIT License
└── README.md
```

### Old experiments

The `Experiment*` folders come from an earlier manuscript on the same model ("A Projective Framework for Inversion-Free Elliptic Curve Cryptography on Isogonal Pivotal Cubics"). They keep their original labels, which is why the numbering is not consecutive (EXP-2 and EXP-6 are reported together). Each folder has its own `README.md`.

Not every number in those folders appears in the article, and none of them should be compared with the article's C benchmark:

| Folder | Role for the article |
| --- | --- |
| `pk_opt/` | **Main source** of the optimized formulas, operation counts and the pK vs. Jacobian benchmark. |
| `Experiment1/` | Supports the correctness checks of the primitives against the Weierstrass oracle. |
| `Experiment2_6/` | Supports the (2,3)/(3,2) law-space computation (dimension 9) and its certification. |
| `Experiment7/` | Supports the completeness subsets and `r_min = 3`. |
| `Experiment8/` | Static cost analysis supports the cost of the (2,3) triple. Its Python wall-clock timings are earlier results, not used in the article. |
| `Experiment3/`, `Experiment4/` | Earlier unoptimized timing and OpenMP-scaling studies. Kept for the record, not used in the article. |

## Test vector

All experiments and the article use the same 256-bit **arithmetic-only** reproducibility vector.

- **Prime field:** `p = 106839527430202782610735740077697524141755216002708025221579380756122101340723`
- **Prime subgroup order:** `q = 26709881857550695652683935019424381035438804000677006305394845189030525335181`
- **Relation:** `p = 4q - 1`, so `#E(F_p) = p + 1 = 4q`
- **Pivot:** `D = (5 : 1 : 7)`
- **pK cubic:** `5α(β² − γ²) + β(γ² − α²) + 7γ(α² − β²) = 0`
- **Weierstrass model:** `E_W : Y² = X³ − X/36`, `j = 1728`
- **Generator** of the order-`q` subgroup in the `D`-centred group (`[q]G = D`):
  `G = (1 : 43885198696659819331868805857976446729580581201517994933321256281291666494433 : 52096118237682591038514212999531378342298973921021687132305677862559506645034)`

> **Security notice.** The vector has `j = 1728` and `p ≡ 3 (mod 4)`, so the curve is **supersingular** with embedding degree 2. The MOV and Frey–Rück reductions map the discrete logarithm into `F_{p²}` (512 bits), so the 254-bit `q` gives nowhere near 127-bit security. The vector is for testing formulas only and must **not** be used as a security parameter set. A security-grade instance needs an ordinary pK curve with a prime-order subgroup of cofactor divisible by 4; none is constructed here.

## Main results of the article

- **Group law.** With identity `D`, `P1 ⊕ P2 = Third(P1,P2)*` and `ι(P) = D* ⊖ P`. At `I = (1:1:1)` the isogonal conjugation `ι` is the group inverse, and the `I`-centred sum is `Third(Third(P,Q), I)`.
- **Torsion.** Every pK curve has full rational 2-torsion, so `4 | #E(F_p)` and every odd prime-order subgroup has cofactor divisible by 4.
- **Symmetry obstruction.** The plane line bundle is never symmetric for the `D`-centred law. Hence there are no bidegree-(2,2) addition laws, `dim A_{2,3} = 9`, no pair of (2,3) laws is complete, and three suffice generically. For the test vector `r_min = 3`, attained by seven explicit rational triples. No degree-4 doubling exists.
- **Formulas.** Inversion-free doubling in `15M + 5S + 5C` and mixed addition in `16M + 3S + 4C` (cached gradient of the addend), for the pivot `(5:1:7)`.
- **Cost.** Against Jacobian coordinates on the isomorphic curve `Y² = X³ − X`, in the same C backend, on an Intel Xeon workstation and a Windows MINGW64 host: doubling about 2.05–2.13x slower, mixed affine addition about 1.78–1.92x slower, scalar multiplication (model: 254 doublings + 51 table additions) about 2.0–2.1x slower.

## Reproducing the results

Requirements: Python 3, SageMath, PARI/GP, and a C compiler (GCC).

```
# EXP-1: correctness against the Weierstrass oracle
cd Experiment1
python run_experiment.py --N-add 500 --N-double 500 --N-neg 500 \
    --N-assoc 200 --N-step-45 500 --N-step-12 50 --seed 0

# EXP-2 / EXP-6: extraction and certification of the (2,3)/(3,2) basis
cd Experiment2_6
python exp6.py
sage sage_verify.sage

# EXP-7: completeness / minimality search
cd Experiment7
sage exp7.sage

# EXP-8: static cost analysis
cd Experiment8
python exp8.py
```

**Optimized formulas and C benchmark** (`pk_opt/`). The scripts for the optimized formulas and operation counts, the `I`-centred law, the bilinear translation test, the law-space row analysis and the C benchmark are in `pk_opt/`. The cross-platform reproduction of the benchmark table, including the exact cross-check of the affine and cached addition paths on the test vector, is in `pk_opt/reproduce/` (see `run_check2.sh`). Flags used in the article: `gcc -O3 -march=native` on the Xeon host; on MINGW64 (UCRT64, GCC 16.2.0) the same plus `-std=gnu11 -fno-strict-aliasing -fno-plt`.

**Primality certificates.**

```
gp -q ecpp_certificates.gp
```

This prints `1 1 1 1` when all certificates verify. The certificates were generated with PARI/GP 2.15.4.

## License

Released under the MIT License. See `LICENSE`.

## Citation

```
@misc{Khayou2026PivotalIsogonalCubics,
  author = {Khayou, Hussein},
  title  = {Pivotal Isogonal Cubics as Elliptic Curves:
            Group Law, Addition-Law Spaces, and the
            Cost of the Native Plane Embedding},
  year   = {2026},
  note   = {Manuscript; accompanying code and reproducibility artifacts},
  url    = {https://github.com/hkhmycodes/H-ECC/}
}
```
