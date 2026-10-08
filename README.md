# H-ECC: Projective Framework for Inversion-Free ECC on Isogonal Pivotal Cubics

This repository contains the code, data, and certificates accompanying the paper **"A Projective Framework for Inversion-Free Elliptic Curve Cryptography on Isogonal Pivotal Cubics."**

The pK representation is an alternative arithmetic representation of an elliptic-curve group. It is **not** a new hardness assumption. The repository is intended to make every numerical and algebraic claim in the paper reproducible.

## Repository layout

```
H-ECC/
├── Experiment1/          # EXP-1: Algebraic correctness against the Weierstrass oracle
├── Experiment2_6/        # EXP-2 / EXP-6: Extraction and certification of the (2,3)/(3,2) basis
├── Experiment3/          # EXP-3: CPU arithmetic cost (unoptimized composite pK vs. Jacobian)
├── Experiment4/          # EXP-4: CPU parallel scaling (OpenMP)
├── Experiment7/          # EXP-7: Minimality search for a complete subcollection of A_{2,3}
├── Experiment8/          # EXP-8: Static cost analysis and wall-clock benchmark of the (2,3) triple
├── ecpp_certificates.gp  # ECPP primality certificates for p and q (PARI/GP)
├── LICENSE               # MIT License
└── README.md
```

Each experiment folder contains its own `README.md` with methodology, usage, and expected output.

## Test vector

All experiments use the same 256-bit arithmetic-only reproducibility vector.

- **Prime field:**
  `p = 106839527430202782610735740077697524141755216002708025221579380756122101340723`
- **Prime subgroup order:**
  `q = 26709881857550695652683935019424381035438804000677006305394845189030525335181`
- **Relation:** `p = 4q - 1`
- **pK cubic:**
  `5α(β² − γ²) + β(γ² − α²) + 7γ(α² − β²) = 0`
- **Pivot:** `D = (5 : 1 : 7)`
- **Birationally equivalent Weierstrass curve:**
  `E_W : Y² = X³ − X/36`

> **Security notice.** The vector has `j = 1728` and `p ≡ 3 (mod 4)`. The curve is **supersingular**, with embedding degree 2. It is supplied solely as a reproducibility vector for the algebraic formulas and the addition-law extraction. It must **not** be used as a security parameter set.

## Experiments at a glance

| Folder | Experiment | Purpose |
|---|---|---|
| `Experiment1/` | EXP-1 | Validate pK primitives against an independent Weierstrass-model oracle. |
| `Experiment2_6/` | EXP-2 / EXP-6 | Extract the nine-dimensional (2,3) addition-law basis; certify each law independently in SageMath. |
| `Experiment3/` | EXP-3 | Measure the cost of the composite pK formulas against generic Jacobian baselines. |
| `Experiment4/` | EXP-4 | Measure OpenMP parallel scaling of the pK representation against Jacobian. |
| `Experiment7/` | EXP-7 | Determine the minimal cardinality of a complete subsystem of A_{2,3}. |
| `Experiment8/` | EXP-8 | Static CSE-aware cost model and wall-clock benchmark of the complete (2,3) triple. |

Experiment labels follow the identifiers in the paper and are not consecutive; EXP-2 and EXP-6 are reported together.

## Reproducing the results

The experiments use Python 3, SageMath, and PARI/GP. Each folder contains its own instructions, but the typical workflow is:

```bash
# EXP-1: correctness
cd Experiment1
python run_experiment.py --N-add 500 --N-double 500 --N-neg 500 \
    --N-assoc 200 --N-step-45 500 --N-step-12 50 --seed 0

# EXP-2 / EXP-6: extraction and certification
cd Experiment2_6
python exp6.py
sage sage_verify.sage

# EXP-7: minimality search
cd Experiment7
sage exp7.sage

# EXP-8: static and wall-clock benchmark
cd Experiment8
python exp8.py
```

The ECPP certificates can be verified independently:

```bash
gp -q ecpp_certificates.gp
```

This prints `1 1 1 1` when all certificates verify.

## Key results

- **Composite pK formulas** are `2.7×–3.8×` slower than Weierstrass–Jacobian in the unoptimized EXP-3 implementation.
- After CSE optimization (EXP-8), composite pK is about `1.4×` slower on addition and `1.9×` on doubling than Jacobian in pure Python.
- The complete **(2,3) triple** is about `9.5×` slower than the optimized composite pK on addition and `5.2×` on doubling; against Jacobian it is about `13.5×` and `10.1×` slower, respectively.
- **CPU parallelism does not close the gap** (EXP-4): the relative throughput `R ≈ 0.42–0.46` is constant across `T = 1, 2, 4, 8`.
- **No pair of (2,3) laws is complete** (EXP-7); the minimal complete subcollection has `r_min = 3`.

## Primality certificates

`ecpp_certificates.gp` contains ECPP primality certificates for `p` and `q` generated with PARI/GP 2.15.4. It is included at the repository root for convenience.

## License

This project is released under the MIT License. See `LICENSE` for details.

## Citation

If you use this code or data in your work, please cite the accompanying paper:

```bibtex
@article{Khayou2026HECC,
  title   = {A Projective Framework for Inversion-Free Elliptic Curve
             Cryptography on Isogonal Pivotal Cubics},
  author  = {Khayou, Hussein},
  journal = {Designs, Codes and Cryptography},
  year    = {2026},
  note    = {Code and data available at
             \url{https://github.com/hkhmycodes/H-ECC/}}
}
```