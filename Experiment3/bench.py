"""Timing and operation-count harness for EXP-3.

Methodology:
  1. Test points generated via the Experiment-1 oracle.
  2. Exceptional addition pairs (P = Q, P = -Q) rejected and counted.
  3. Four-way scalar cross-check against the oracle, performed in native
     Python arithmetic, before any timing.
  4. Conversion to the active field backend (py or gmp) happens only after
     validation passes.
  5. Warm-up pass, then median of `repeats` passes.
  6. Failure if insufficient inputs.
  7. T_M microbenchmark.

Representation conventions (important):
  - The Experiment-1 oracle works with **affine** Weierstrass points:
        W = (X, Y)       or None for the identity
  - The Weierstrass arithmetic in ws_arith.py works with **Jacobian**
    points:
        P = (X, Y, Z)    with identity (0, 1, 0)
  - Mixed-addition routines take one Jacobian operand and one **affine**
    2-tuple operand.
This module keeps both conventions explicit and converts at the boundary.
"""

from __future__ import annotations
import json
import os
import random
import sys
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / "Experiment1"))

# Backend switch must run before importing field / pk_arith / ws_arith.
if os.environ.get("EXP3_FIELD", "py") == "gmp":
    import field_gmp
    sys.modules["field"] = field_gmp

from pk_core import PKCurve                    # noqa: E402
from pk_oracle import Oracle, MapUndefined     # noqa: E402

import field as f                              # noqa: E402
import pk_arith as pa                          # noqa: E402
import ws_arith as wa                          # noqa: E402


P_PRIME = 106839527430202782610735740077697524141755216002708025221579380756122101340723
Q_ORDER = 26709881857550695652683935019424381035438804000677006305394845189030525335181
U, V, W = 5, 1, 7
G_PK = (1,
        43885198696659819331868805857976446729580581201517994933321256281291666494433,
        52096118237682591038514212999531378342298973921021687132305677862559506645034)


def setup():
    p, q = P_PRIME, Q_ORDER
    curve = PKCurve(U, V, W, p)
    oracle = Oracle(curve)
    uvw = (curve.u, curve.v, curve.w)
    Ds = curve.Ds
    gDs = pa.pk_grad(Ds, uvw, p)
    a_curve = (-pow(36, p - 2, p)) % p
    return curve, oracle, uvw, Ds, gDs, a_curve, p, q


def _to_field_point(P: Tuple) -> Tuple:
    return tuple(f.from_int(x) for x in P)


def _affine_to_jac(W_aff: Optional[Tuple[int, int]], p: int) -> Tuple:
    """Convert an affine Weierstrass point (or None) to Jacobian."""
    if W_aff is None:
        return wa.IDENTITY
    return (W_aff[0] % p, W_aff[1] % p, 1)


# ---------------------------------------------------------------------------
# Oracle-generated test points and pairs (native Python arithmetic)
# ---------------------------------------------------------------------------

def _gen_one_point(rng, q, W_G_aff, oracle, p):
    """Sample k, compute [k] W_G_aff on the oracle, and return
    (P_pk, W_jac) where W_jac is Jacobian (X, Y, 1)."""
    k = rng.randrange(1, q)
    W_k = oracle.split.scalar_mul(k, W_G_aff)   # affine or None
    if W_k is None:
        return None
    try:
        P_pk = oracle.from_w(W_k)
    except MapUndefined:
        return None
    if P_pk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        return None
    return (P_pk, (W_k[0] % p, W_k[1] % p, 1))


def gen_pairs(N, rng, q, G_pk, oracle, p):
    """Returns (pairs, counts) where each pair is
    (P_pk, Q_pk, W_P_jac, W_Q_jac) in native Python arithmetic."""
    W_G_aff = oracle.to_w(G_pk)
    pairs: List[Tuple] = []
    counts = {
        "attempts": 0,
        "rejected_map": 0,
        "rejected_base_locus": 0,
        "rejected_equal": 0,
        "rejected_inverse": 0,
    }
    max_attempts = 100 * (2 * N)
    while len(pairs) < N and counts["attempts"] < max_attempts:
        counts["attempts"] += 1
        first = _gen_one_point(rng, q, W_G_aff, oracle, p)
        if first is None:
            counts["rejected_map"] += 1
            continue
        second = _gen_one_point(rng, q, W_G_aff, oracle, p)
        if second is None:
            counts["rejected_map"] += 1
            continue

        P_pk, W_P_jac = first
        Q_pk, W_Q_jac = second

        if W_P_jac == W_Q_jac:
            counts["rejected_equal"] += 1
            continue
        if (W_P_jac[0] == W_Q_jac[0]
                and ((W_P_jac[1] + W_Q_jac[1]) % p == 0)):
            counts["rejected_inverse"] += 1
            continue

        pairs.append((P_pk, Q_pk, W_P_jac, W_Q_jac))

    if len(pairs) < N:
        raise RuntimeError(
            f"gen_pairs: only {len(pairs)}/{N} usable pairs "
            f"after {counts['attempts']} attempts ({counts})"
        )
    return pairs, counts


def gen_scalars(N_scalar, rng, q, G_pk, oracle, p):
    """Returns (scalars, W_G_jac).

    W_G_jac is the Jacobian form of Phi_D(G), used by the WS scalar routines.
    The affine form can always be recovered via wa.ws_to_affine when needed.
    """
    W_G_aff = oracle.to_w(G_pk)
    W_G_jac = _affine_to_jac(W_G_aff, p)
    scalars = [rng.randrange(1, q) for _ in range(N_scalar)]
    return scalars, W_G_jac


# ---------------------------------------------------------------------------
# Four-way correctness checks (native Python arithmetic)
# ---------------------------------------------------------------------------

def crosscheck_pairs(pairs, oracle, uvw, Ds, gDs, a_curve, p):
    for idx, (P_pk, Q_pk, W_P_jac, W_Q_jac) in enumerate(pairs):
        # pK side
        R_pk = pa.pk_add(P_pk, Q_pk, uvw, p)
        D_pk = pa.pk_double(P_pk, uvw, p)
        N_pk = pa.pk_neg(P_pk, uvw, Ds, gDs, p)

        # Weierstrass side (Jacobian input, affine output)
        R_ws = wa.ws_add(W_P_jac, W_Q_jac, p)
        D_ws = wa.ws_double(W_P_jac, a_curve, p)
        N_ws = wa.ws_neg(W_P_jac, p)

        if oracle.to_w(R_pk) != wa.ws_to_affine(R_ws, p):
            raise RuntimeError(f"crosscheck add mismatch at idx={idx}")
        if oracle.to_w(D_pk) != wa.ws_to_affine(D_ws, p):
            raise RuntimeError(f"crosscheck double mismatch at idx={idx}")
        if oracle.to_w(N_pk) != wa.ws_to_affine(N_ws, p):
            raise RuntimeError(f"crosscheck neg mismatch at idx={idx}")


def crosscheck_scalars(scalars, G_pk, W_G_jac, oracle, uvw, Ds, a_curve, p):
    """Four-way scalar check:

      W_oracle = oracle.split.scalar_mul(k, Phi_D(G))     (independent)
      W_pk     = Phi_D(pk_scalar(k, G_pk))                (pK side)
      W_rtl    = ws_scalar_mul(k, Phi_D(G))               (WS right-to-left)
      W_mix    = ws_scalar_mul_mixed(k, Phi_D(G))         (WS mixed)

    All four must coincide.  Note: the oracle consumes the **affine** form
    of Phi_D(G); the WS routines consume Jacobian and affine respectively.
    """
    W_G_affine = wa.ws_to_affine(W_G_jac, p)
    for idx, k in enumerate(scalars):
        W_oracle = oracle.split.scalar_mul(k, W_G_affine)
        if W_oracle is None:
            raise RuntimeError(f"oracle scalar produced identity at idx={idx}")
        W_oracle_aff = (W_oracle[0] % p, W_oracle[1] % p)

        S_pk = pa.pk_scalar_mul(k, G_pk, uvw, p, Ds)
        W_pk = oracle.to_w(S_pk)

        W_rtl = wa.ws_to_affine(
            wa.ws_scalar_mul(k, W_G_jac, a_curve, p), p)
        W_mix = wa.ws_to_affine(
            wa.ws_scalar_mul_mixed(k, W_G_affine, a_curve, p), p)

        if W_pk != W_oracle_aff:
            raise RuntimeError(
                f"scalar check (pK) failed at idx={idx}, k={k}")
        if W_rtl != W_oracle_aff:
            raise RuntimeError(
                f"scalar check (WS-RL) failed at idx={idx}, k={k}")
        if W_mix != W_oracle_aff:
            raise RuntimeError(
                f"scalar check (WS-mixed) failed at idx={idx}, k={k}")


# ---------------------------------------------------------------------------
# Timing and counting
# ---------------------------------------------------------------------------

def time_calls(fn: Callable, args_list: List[Tuple], repeats: int) -> float:
    for args in args_list:
        fn(*args)
    samples = []
    for _ in range(repeats):
        t0 = time.perf_counter_ns()
        for args in args_list:
            fn(*args)
        t1 = time.perf_counter_ns()
        samples.append((t1 - t0) / len(args_list))
    samples.sort()
    return samples[len(samples) // 2]


def time_field_mul(N: int, rng: random.Random, p: int,
                   repeats: int) -> float:
    pairs = [(rng.randrange(1, int(p)), rng.randrange(1, int(p)))
             for _ in range(N)]
    return time_calls(f.mul, [(a, b, p) for a, b in pairs], repeats)


def count_ops(fn: Callable, args: Tuple) -> Dict[str, int]:
    f.start_tracking()
    fn(*args)
    return f.stop_tracking()


# ---------------------------------------------------------------------------
# Single-seed run
# ---------------------------------------------------------------------------

def run_one_seed(seed: int, N: int, N_scalar: int,
                 repeats: int, scalar_repeats: int,
                 verbose: bool = True) -> Dict[str, Any]:
    curve, oracle, uvw, Ds, gDs, a_curve, p, q = setup()
    rng = random.Random(seed)
    backend = os.environ.get("EXP3_FIELD", "py")

    if verbose:
        print(f"[seed {seed}] generating {N} pairs, {N_scalar} scalars "
              f"via oracle (backend={backend}) ...")

    # --- native-arithmetic generation & validation ---------------------
    pairs_py, pair_counts = gen_pairs(N, rng, q, G_PK, oracle, p)
    scalars_py, W_G_jac_native = gen_scalars(
        N_scalar, rng, q, G_PK, oracle, p)

    if verbose:
        print(f"[seed {seed}] pair counts: {pair_counts}")
        print(f"[seed {seed}] four-way cross-check (native Python) ...")

    crosscheck_pairs(pairs_py, oracle, uvw, Ds, gDs, a_curve, p)
    crosscheck_scalars(scalars_py, G_PK, W_G_jac_native,
                       oracle, uvw, Ds, a_curve, p)

    if verbose:
        print(f"[seed {seed}] validation passed")

    # --- convert to active field backend for the timing pass -----------
    P_pk = [_to_field_point(pr[0]) for pr in pairs_py]
    Q_pk = [_to_field_point(pr[1]) for pr in pairs_py]
    W_P_jac = [_to_field_point(pr[2]) for pr in pairs_py]
    W_Q_jac = [_to_field_point(pr[3]) for pr in pairs_py]

    uvw_f = tuple(f.from_int(x) for x in uvw)
    Ds_f = _to_field_point(Ds)
    gDs_f = _to_field_point(gDs)
    G_pk_f = _to_field_point(G_PK)
    W_G_jac_f = _to_field_point(W_G_jac_native)
    a_curve_f = f.from_int(a_curve)
    p_f = f.from_int(p)

    W_G_affine_f = wa.ws_to_affine(W_G_jac_f, p_f)
    scalars = list(scalars_py)

    # --- argument lists -------------------------------------------------
    add_pk_args = [(P_pk[i], Q_pk[i], uvw_f, p_f) for i in range(N)]

    # ws_add expects Jacobian on both sides
    add_ws_args = [(W_P_jac[i], W_Q_jac[i], p_f) for i in range(N)]

    # ws_madd expects (Jacobian, affine 2-tuple)
    add_mix_args = [(W_P_jac[i], (W_Q_jac[i][0], W_Q_jac[i][1]),
                     a_curve_f, p_f)
                    for i in range(N)]

    dbl_pk_args = [(P_pk[i], uvw_f, p_f) for i in range(N)]
    dbl_ws_args = [(W_P_jac[i], a_curve_f, p_f) for i in range(N)]

    neg_pk_args = [(P_pk[i], uvw_f, Ds_f, gDs_f, p_f) for i in range(N)]
    neg_ws_args = [(W_P_jac[i], p_f) for i in range(N)]

    sm_pk_args = [(scalars[i], G_pk_f, uvw_f, p_f, Ds_f)
                  for i in range(N_scalar)]

    # ws_scalar_mul takes a Jacobian base point
    sm_ws_args = [(scalars[i], W_G_jac_f, a_curve_f, p_f)
                  for i in range(N_scalar)]

    # ws_scalar_mul_mixed takes an affine base point
    sm_mix_args = [(scalars[i], W_G_affine_f, a_curve_f, p_f)
                   for i in range(N_scalar)]

    # --- timing ---------------------------------------------------------
    if verbose:
        print(f"[seed {seed}] T_M microbenchmark ...")
    t_M = time_field_mul(N, rng, p_f, repeats)

    if verbose:
        print(f"[seed {seed}] timing ...")
    t_pk_add  = time_calls(pa.pk_add,   add_pk_args, repeats)
    t_ws_add  = time_calls(wa.ws_add,   add_ws_args, repeats)
    t_mix_add = time_calls(wa.ws_madd,  add_mix_args, repeats)

    t_pk_dbl  = time_calls(pa.pk_double, dbl_pk_args, repeats)
    t_ws_dbl  = time_calls(wa.ws_double, dbl_ws_args, repeats)

    t_pk_neg  = time_calls(pa.pk_neg,   neg_pk_args, repeats)
    t_ws_neg  = time_calls(wa.ws_neg,   neg_ws_args, repeats)

    t_pk_sm   = time_calls(pa.pk_scalar_mul,       sm_pk_args,  scalar_repeats)
    t_ws_sm   = time_calls(wa.ws_scalar_mul,       sm_ws_args,  scalar_repeats)
    t_mix_sm  = time_calls(wa.ws_scalar_mul_mixed, sm_mix_args, scalar_repeats)

    # --- op counts ------------------------------------------------------
    c_pk_add  = count_ops(pa.pk_add,   (P_pk[0], Q_pk[0], uvw_f, p_f))
    c_ws_add  = count_ops(wa.ws_add,   (W_P_jac[0], W_Q_jac[0], p_f))
    c_mix_add = count_ops(wa.ws_madd,  (W_P_jac[0],
                                        (W_Q_jac[0][0], W_Q_jac[0][1]),
                                        a_curve_f, p_f))
    c_pk_dbl  = count_ops(pa.pk_double, (P_pk[0], uvw_f, p_f))
    c_ws_dbl  = count_ops(wa.ws_double, (W_P_jac[0], a_curve_f, p_f))
    c_pk_neg  = count_ops(pa.pk_neg,   (P_pk[0], uvw_f, Ds_f, gDs_f, p_f))
    c_ws_neg  = count_ops(wa.ws_neg,   (W_P_jac[0], p_f))

    return {
        "experiment": "EXP-3",
        "field_backend": backend,
        "seed": seed,
        "N": N, "N_scalar": N_scalar,
        "repeats": repeats, "scalar_repeats": scalar_repeats,
        "pair_counts": pair_counts,
        "params": {"p": int(p), "q": int(q),
                   "u": int(U), "v": int(V), "w": int(W),
                   "a_curve": int(a_curve)},
        "T_M_ns": t_M,
        "timing_ns": {
            "pk_add": t_pk_add, "ws_add": t_ws_add, "mix_add": t_mix_add,
            "pk_double": t_pk_dbl, "ws_double": t_ws_dbl,
            "pk_neg": t_pk_neg, "ws_neg": t_ws_neg,
            "pk_scalar": t_pk_sm, "ws_scalar": t_ws_sm,
            "mix_scalar": t_mix_sm,
        },
        "time_over_TM": {
            k: v / t_M for k, v in {
                "pk_add": t_pk_add, "ws_add": t_ws_add, "mix_add": t_mix_add,
                "pk_double": t_pk_dbl, "ws_double": t_ws_dbl,
                "pk_neg": t_pk_neg, "ws_neg": t_ws_neg,
                "pk_scalar": t_pk_sm, "ws_scalar": t_ws_sm,
                "mix_scalar": t_mix_sm,
            }.items()
        },
        "ratios_pk_over_ws": {
            "add":  t_pk_add  / t_ws_add,
            "add_vs_mixed": t_pk_add / t_mix_add,
            "double": t_pk_dbl / t_ws_dbl,
            "neg":  t_pk_neg  / t_ws_neg,
            "scalar": t_pk_sm / t_ws_sm,
            "scalar_vs_mixed": t_pk_sm / t_mix_sm,
        },
        "op_counts": {
            "pk_add": c_pk_add, "ws_add": c_ws_add, "mix_add": c_mix_add,
            "pk_double": c_pk_dbl, "ws_double": c_ws_dbl,
            "pk_neg": c_pk_neg, "ws_neg": c_ws_neg,
        },
    }


# ---------------------------------------------------------------------------
# Multi-seed aggregation
# ---------------------------------------------------------------------------

def run(seeds: List[int], N: int, N_scalar: int,
        repeats: int, scalar_repeats: int,
        outdir: str = "results") -> Dict[str, Any]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    per_seed = [run_one_seed(s, N, N_scalar, repeats, scalar_repeats)
                for s in seeds]

    def _median(xs):
        s = sorted(xs)
        return s[len(s) // 2]

    timing_keys = ["pk_add", "ws_add", "mix_add", "pk_double", "ws_double",
                   "pk_neg", "ws_neg", "pk_scalar", "ws_scalar", "mix_scalar"]
    ratio_keys = ["add", "add_vs_mixed", "double", "neg",
                  "scalar", "scalar_vs_mixed"]

    aggregate = {
        "experiment": "EXP-3",
        "scope": "GMP-backed Python prototype; not the C++/GMP benchmark",
        "field_backend": os.environ.get("EXP3_FIELD", "py"),
        "seeds": seeds,
        "N": N, "N_scalar": N_scalar,
        "repeats": repeats, "scalar_repeats": scalar_repeats,
        "T_M_ns_median": _median([r["T_M_ns"] for r in per_seed]),
        "timing_ns_median": {
            k: _median([r["timing_ns"][k] for r in per_seed])
            for k in timing_keys
        },
        "time_over_TM_median": {
            k: _median([r["time_over_TM"][k] for r in per_seed])
            for k in timing_keys
        },
        "ratios_pk_over_ws_median": {
            k: _median([r["ratios_pk_over_ws"][k] for r in per_seed])
            for k in ratio_keys
        },
        "op_counts": per_seed[0]["op_counts"],
        "pair_counts_per_seed": [r["pair_counts"] for r in per_seed],
    }

    report = {"per_seed": per_seed, "aggregate": aggregate}
    (out / "report.json").write_text(json.dumps(report, indent=2))
    lines = _format_aggregate(aggregate)
    (out / "report.txt").write_text("\n".join(lines))
    print()
    print("\n".join(lines))
    return report


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt_time(ns: float) -> str:
    if ns < 1_000:
        return f"{ns:.1f} ns"
    if ns < 1_000_000:
        return f"{ns / 1_000:.2f} us"
    return f"{ns / 1_000_000:.2f} ms"


def _fmt_ops(c: Dict[str, int]) -> str:
    return (f"{c['M']:3d}M + {c['S']:3d}S + "
            f"{c['C']:3d}C + {c['A']:4d}A")


def _format_aggregate(a: Dict[str, Any]) -> List[str]:
    t = a["timing_ns_median"]
    r = a["ratios_pk_over_ws_median"]
    over = a["time_over_TM_median"]
    oc = a["op_counts"]
    TM = a["T_M_ns_median"]
    lines = []
    lines.append("=" * 78)
    lines.append("EXP-3  pK current vs Weierstrass-Jacobian")
    lines.append(f"       backend = {a['field_backend']}  "
                 f"seeds = {a['seeds']}")
    lines.append("=" * 78)
    lines.append(f"N               : {a['N']}")
    lines.append(f"N_scalar        : {a['N_scalar']}")
    lines.append(f"repeats         : {a['repeats']}   "
                 f"scalar_repeats: {a['scalar_repeats']}")
    lines.append(f"T_M (field mul) : {_fmt_time(TM)}")
    lines.append("")
    lines.append(f"{'Primitive':<14} {'pK':>14} {'WS gen':>14} "
                 f"{'WS mixed':>14} {'pK/WS':>8} {'pK/mix':>8}")
    lines.append("-" * 78)
    lines.append(f"{'addition':<14} {_fmt_time(t['pk_add']):>14} "
                 f"{_fmt_time(t['ws_add']):>14} {_fmt_time(t['mix_add']):>14} "
                 f"{r['add']:>7.2f}x {r['add_vs_mixed']:>7.2f}x")
    lines.append(f"{'doubling':<14} {_fmt_time(t['pk_double']):>14} "
                 f"{_fmt_time(t['ws_double']):>14} {'---':>14} "
                 f"{r['double']:>7.2f}x {'---':>8}")
    lines.append(f"{'negation':<14} {_fmt_time(t['pk_neg']):>14} "
                 f"{_fmt_time(t['ws_neg']):>14} {'---':>14} "
                 f"{r['neg']:>7.2f}x {'---':>8}")
    lines.append(f"{'scalar-mul':<14} {_fmt_time(t['pk_scalar']):>14} "
                 f"{_fmt_time(t['ws_scalar']):>14} "
                 f"{_fmt_time(t['mix_scalar']):>14} "
                 f"{r['scalar']:>7.2f}x {r['scalar_vs_mixed']:>7.2f}x")
    lines.append("")
    lines.append("Timing normalized by measured field-multiplication baseline:")
    for label in ["pk_add", "ws_add", "mix_add",
                  "pk_double", "ws_double",
                  "pk_neg", "ws_neg",
                  "pk_scalar", "ws_scalar", "mix_scalar"]:
        lines.append(f"  {label:<10}: {over[label]:>7.2f} T_M")
    lines.append("")
    lines.append("Operation counts (per call):")
    for label, key in [("pK add", "pk_add"), ("WS add", "ws_add"),
                       ("WS mixed", "mix_add"),
                       ("pK double", "pk_double"), ("WS double", "ws_double"),
                       ("pK neg", "pk_neg"), ("WS neg", "ws_neg")]:
        lines.append(f"  {label:<10}: {_fmt_ops(oc[key])}")
    lines.append("")
    lines.append("Note: operation counts and wall-clock times are two")
    lines.append("independent metrics. The generic pK-vs-WS comparison uses")
    lines.append("the same right-to-left binary double-and-add; the mixed")
    lines.append("scalar-mul baseline is a separate left-to-right traversal")
    lines.append("and is reported as a practical baseline, not as a formula-")
    lines.append("cost control.")
    return lines