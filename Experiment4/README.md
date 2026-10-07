# Experiment 4 — CPU parallel scaling of the pK representation

Benchmark of the composite pK formulas of Sections 4.1--4.2 of the
manuscript against the Weierstrass--Jacobian implementation of EXP-3,
on the same 256-bit reproducibility vector of Section 6.5, under
OpenMP parallelisation.

This directory accompanies **Section 10** (Experimental program,
Stage I) of the manuscript *A Projective Framework for Inversion-Free
Elliptic Curve Cryptography on Isogonal Pivotal Cubics*.

---

## 1. Scientific question

The composite pK formulas of Sections 4.1--4.2 are approximately
2.7x--3.1x slower than Weierstrass--Jacobian on addition and 3.0x slower
on scalar multiplication (EXP-3). The pK representation is proposed as
a parallelizable alternative arithmetic substrate. EXP-4 asks:

> Does CPU parallelism close the pK/Jacobian arithmetic gap?

The headline quantity is the relative throughput

    R(N, T) = throughput_pK(N, T) / throughput_jac(N, T)

measured on the same points and the same scalars in both representations.
If R is constant in T, parallelism scales both representations
identically and does not recover the arithmetic gap. If R rises with T,
parallelism partially closes the gap. If R falls with T, parallelism
widens it.

---

## 2. Inputs and outputs

### Inputs

None external. The test vector is hard-coded:

- p = 106839527430202782610735740077697524141755216002708025221579380756122101340723
- q = 26709881857550695652683935019424381035438804000677006305394845189030525335181
- pK pivot D = (5 : 1 : 7), generator G of order q in the D-centered group
- E_W : Y^2 = X^3 - X/36

### Outputs (written to results/)

For each (N, T, seed) configuration:

    results/exp4_N<N>_T<T>_s<seed>.json     JSON record (stdout)
    results/exp4_N<N>_T<T>_s<seed>.log      stderr progress log

    results/env.txt                         machine + toolchain record

### JSON schema

```json
{
  "version": "rev2",
  "n": 4096,
  "seed": 0,
  "repeats": 5,
  "workload_generation_s": 13.2,
  "validation_subset": 20,
  "validation_passed": true,
  "oracle_checked": 4096,
  "oracle_skipped": 0,
  "oracle_mismatches": 0,
  "oracle_passed": true,
  "outputs_identical_across_threads": true,
  "team_sizes_as_requested": true,
  "env": { ... },
  "threads": [1, 2, 4, 8],
  "rows": [
    { "kind": "pk",  "T": 1, "team": 1,
      "T_median_s": 14.4,
      "throughput_ops": 284.3,
      "speedup": 1.0, "efficiency": 1.0,
      "hash": "86714717b8632d8b",
      "identical_to_T1": true,
      "samples_s": [ ... ] },
    ...
  ]
}