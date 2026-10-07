"""All test steps for Experiment 1.

Every random test:
  - samples scalars uniformly from [0, q);
  - loops until `N` *successful* trials are obtained (hard cap 10*N attempts);
  - PASS requires >= N successes and zero failures (no vacuous PASS).

Reported counters per step:
    n_requested, n_attempts, n_trials, n_skipped, n_failed.
"""

from __future__ import annotations
import random
from typing import Any, Dict, List

from pk_core import (
    PKCurve, proj_eq, proj_normalize, is_zero_point, are_collinear, is_prime,
)
from pk_oracle import Oracle, MapUndefined, phi_0, phi_0_inv


def _record(failures, **kw):
    failures.append(kw)


def _w_eq(P, Q, p: int) -> bool:
    if P is None and Q is None:
        return True
    if P is None or Q is None:
        return False
    return (P[0] - Q[0]) % p == 0 and (P[1] - Q[1]) % p == 0


def _finalize(res: Dict[str, Any], N: int) -> Dict[str, Any]:
    """PASS iff n_trials >= N and n_failed == 0."""
    res["n_requested"] = N
    ok = (res.get("n_failed", 0) == 0) and (res.get("n_trials", 0) >= N)
    res["passed"] = ok
    if not ok:
        if res.get("n_failed", 0) > 0:
            res["reason"] = f"{res['n_failed']} failures"
        else:
            res["reason"] = (f"only {res.get('n_trials', 0)} of {N} "
                             f"successful trials "
                             f"(attempts={res.get('n_attempts', 0)})")
    return res


def _new_res(name: str) -> Dict[str, Any]:
    return {"name": name, "n_requested": 0, "n_attempts": 0,
            "n_trials": 0, "n_skipped": 0, "n_failed": 0}


# ---------------------------------------------------------------------------
# Step 0 — parameters, primality, generator order
# ---------------------------------------------------------------------------

def step_0_parameters(curve: PKCurve, q: int, D, G) -> Dict[str, Any]:
    res: Dict[str, Any] = {"name": "step_0_parameters", "checks": []}
    p = curve.p

    res["checks"].append(("p is prime", is_prime(p)))
    res["checks"].append(("q is prime", is_prime(q)))
    res["checks"].append(("p = 4q - 1", p == 4 * q - 1))
    res["checks"].append(("F(D) = 0", curve.F(D) == 0))
    res["checks"].append(("F(G) = 0", curve.F(G) == 0))
    res["checks"].append(("D != G", not proj_eq(D, G, p)))

    oracle = Oracle(curve)
    try:
        W_G = oracle.to_w(G)
        W_qG = oracle.split.scalar_mul(q, W_G)
        G_back = oracle.from_w(W_qG)
        res["checks"].append(("[q]G = D via oracle", proj_eq(G_back, D, p)))
    except (MapUndefined, Exception) as e:
        res["checks"].append((f"oracle scalar mul failed: {e}", False))

    res["passed"] = all(ok for _, ok in res["checks"])
    return res


# ---------------------------------------------------------------------------
# Step 4 — birational maps (both directions, explicit)
# ---------------------------------------------------------------------------

def step_4_birational_maps(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = _new_res("step_4_birational_maps")
    oracle = Oracle(curve)

    W_D = phi_0(curve.D, curve, oracle.split)
    if W_D is None or not oracle.split.is_on_curve(W_D):
        _record(failures, step="4", kind="phi_0(D) not on split cubic")
        res["n_failed"] += 1
        return _finalize(res, N)

    W_G = oracle.to_w(G)
    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        k = rng.randrange(0, q)
        W = oracle.split.scalar_mul(k, W_G)
        if W is None:
            res["n_skipped"] += 1
            continue

        try:
            P = phi_0_inv(W, curve, oracle.split)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not curve.is_on_curve(P):
            _record(failures, step="4", kind="phi_0_inv off pK curve",
                    k=k, W=W, P=P, F=curve.F(P))
            res["n_failed"] += 1
            continue

        try:
            W0 = phi_0(P, curve, oracle.split)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if W0 is None or not oracle.split.is_on_curve(W0):
            _record(failures, step="4", kind="phi_0 off split cubic",
                    k=k, W=W, P=P, W0=W0)
            res["n_failed"] += 1
            continue

        try:
            P_back = phi_0_inv(W0, curve, oracle.split)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not proj_eq(P_back, P, curve.p):
            _record(failures, step="4", kind="phi_0_inv o phi_0 != id",
                    k=k, P=P, P_back=P_back)
            res["n_failed"] += 1
            continue

        if not _w_eq(W0, W, curve.p):
            _record(failures, step="4", kind="phi_0 o phi_0_inv != id",
                    k=k, W=W, W0=W0)
            res["n_failed"] += 1
            continue

        try:
            W_Phi = oracle.to_w(P)
            P_Phi = oracle.from_w(W_Phi)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not proj_eq(P_Phi, P, curve.p):
            _record(failures, step="4", kind="Phi_D round trip failed",
                    k=k, P=P, P_Phi=P_Phi)
            res["n_failed"] += 1
            continue

        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 5 — generator order via oracle
# ---------------------------------------------------------------------------

def step_5_generator(curve, q, D, G, N, rng, failures) -> Dict[str, Any]:
    res = _new_res("step_5_generator")
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)

    for k in [0, 1, 2, q - 1]:
        res["n_attempts"] += 1
        try:
            P = oracle.from_w(oracle.split.scalar_mul(k, W_G))
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not curve.is_on_curve(P):
            _record(failures, step="5", kind="off-curve", k=k, P=P)
            res["n_failed"] += 1

    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        k = rng.randrange(0, q)
        try:
            P = oracle.from_w(oracle.split.scalar_mul(k, W_G))
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not curve.is_on_curve(P):
            _record(failures, step="5", kind="off-curve", k=k, P=P)
            res["n_failed"] += 1
            continue
        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 6 — isogonal map
# ---------------------------------------------------------------------------

def step_6_isogonal(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = _new_res("step_6_isogonal")
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)
    Ds = curve.Ds

    for Ti in [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1)]:
        if not proj_eq(curve.isogonal(Ti), Ti, curve.p):
            _record(failures, step="6", kind="Ti* != Ti", Ti=Ti)
            res["n_failed"] += 1

    if not proj_eq(curve.isogonal(curve.D), Ds, curve.p):
        _record(failures, step="6", kind="iota(D) != D*",
                iota_D=curve.isogonal(curve.D), Ds=Ds)
        res["n_failed"] += 1

    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        k = rng.randrange(0, q)
        try:
            P = oracle.from_w(oracle.split.scalar_mul(k, W_G))
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if curve.is_base_locus(P):
            res["n_skipped"] += 1
            continue

        Pstar = curve.isogonal(P)
        if not curve.is_on_curve(Pstar):
            _record(failures, step="6", kind="iota off-curve",
                    P=P, Pstar=Pstar)
            res["n_failed"] += 1
            continue

        P2 = curve.isogonal(Pstar)
        if not proj_eq(P2, P, curve.p):
            _record(failures, step="6", kind="iota^2 != id", P=P, P2=P2)
            res["n_failed"] += 1
            continue

        try:
            W_Ps = oracle.to_w(Pstar)
            W_expected = oracle.split.add(oracle.to_w(Ds),
                                          oracle.split.neg(oracle.to_w(P)))
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not _w_eq(W_Ps, W_expected, curve.p):
            _record(failures, step="6", kind="iota mismatch in oracle",
                    P=P, Pstar=Pstar, W_Ps=W_Ps, W_expected=W_expected)
            res["n_failed"] += 1
            continue
        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 7 — negation
# ---------------------------------------------------------------------------

def step_7_negation(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = _new_res("step_7_negation")
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)
    D = curve.D
    Ds = curve.Ds

    # Deterministic: -D* == D* (2-torsion)
    if not proj_eq(curve.neg(Ds), Ds, curve.p):
        _record(failures, step="7", kind="-D* != D*",
                neg_Ds=curve.neg(Ds), Ds=Ds)
        res["n_failed"] += 1

    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        k = rng.randrange(0, q)
        try:
            P = oracle.from_w(oracle.split.scalar_mul(k, W_G))
        except MapUndefined:
            res["n_skipped"] += 1
            continue

        try:
            N_pk = curve.neg(P)
        except (ValueError, ZeroDivisionError):
            res["n_skipped"] += 1
            continue
        if not curve.is_on_curve(N_pk):
            _record(failures, step="7", kind="neg off-curve",
                    P=P, N_pk=N_pk, F=curve.F(N_pk))
            res["n_failed"] += 1
            continue

        try:
            N_or = oracle.oracle_neg(P)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not proj_eq(N_pk, N_or, curve.p):
            _record(failures, step="7", kind="neg mismatch",
                    k=k, P=P, N_pk=N_pk, N_or=N_or)
            res["n_failed"] += 1
            continue

        try:
            S = curve.add(P, N_pk)
        except (ValueError, ZeroDivisionError):
            res["n_skipped"] += 1
            continue
        if not proj_eq(S, D, curve.p):
            _record(failures, step="7", kind="P + (-P) != D",
                    P=P, N_pk=N_pk, S=S)
            res["n_failed"] += 1
            continue

        try:
            NN = curve.neg(N_pk)
        except (ValueError, ZeroDivisionError):
            res["n_skipped"] += 1
            continue
        if not proj_eq(NN, P, curve.p):
            _record(failures, step="7", kind="-(-P) != P", P=P, NN=NN)
            res["n_failed"] += 1
            continue
        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 8 — addition
# ---------------------------------------------------------------------------

def step_8_addition(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = _new_res("step_8_addition")
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)

    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        i = rng.randrange(0, q)
        j = rng.randrange(0, q)
        if i == j:
            res["n_skipped"] += 1
            continue
        try:
            W_P = oracle.split.scalar_mul(i, W_G)
            W_Q = oracle.split.scalar_mul(j, W_G)
            P = oracle.from_w(W_P)
            Q = oracle.from_w(W_Q)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if proj_eq(P, Q, curve.p):
            res["n_skipped"] += 1
            continue

        try:
            R_pk = curve.add(P, Q)
        except (ValueError, ZeroDivisionError):
            res["n_skipped"] += 1
            continue
        if not curve.is_on_curve(R_pk):
            _record(failures, step="8", kind="add off-curve",
                    i=i, j=j, P=P, Q=Q, R_pk=R_pk, F=curve.F(R_pk))
            res["n_failed"] += 1
            continue

        try:
            R_or = oracle.oracle_add(P, Q)
        except MapUndefined:
            res["n_skipped"] += 1
            continue

        if not proj_eq(R_pk, R_or, curve.p):
            _record(failures, step="8", kind="add mismatch",
                    i=i, j=j, P=P, Q=Q, R_pk=R_pk, R_or=R_or,
                    P_norm=proj_normalize(P, curve.p),
                    Q_norm=proj_normalize(Q, curve.p),
                    R_pk_norm=proj_normalize(R_pk, curve.p),
                    R_or_norm=proj_normalize(R_or, curve.p))
            res["n_failed"] += 1
            continue
        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 9 — doubling
# ---------------------------------------------------------------------------

def step_9_doubling(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = _new_res("step_9_doubling")
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)

    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        i = rng.randrange(0, q)
        try:
            P = oracle.from_w(oracle.split.scalar_mul(i, W_G))
        except MapUndefined:
            res["n_skipped"] += 1
            continue

        try:
            D_pk = curve.double(P)
        except (ValueError, ZeroDivisionError):
            res["n_skipped"] += 1
            continue
        if not curve.is_on_curve(D_pk):
            _record(failures, step="9", kind="double off-curve",
                    i=i, P=P, D_pk=D_pk)
            res["n_failed"] += 1
            continue

        try:
            D_or = oracle.oracle_double(P)
        except MapUndefined:
            res["n_skipped"] += 1
            continue
        if not proj_eq(D_pk, D_or, curve.p):
            _record(failures, step="9", kind="double mismatch",
                    i=i, P=P, D_pk=D_pk, D_or=D_or)
            res["n_failed"] += 1
            continue

        try:
            D_add = curve.add(P, P)
        except (ValueError, ZeroDivisionError):
            D_add = None
        if D_add is not None and not proj_eq(D_pk, D_add, curve.p):
            _record(failures, step="9", kind="double != add(P, P)",
                    i=i, P=P, D_pk=D_pk, D_add=D_add)
            res["n_failed"] += 1
            continue
        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 10 — group identities (internal consistency)
# ---------------------------------------------------------------------------

def step_10_group_identities(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = _new_res("step_10_group_identities")
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)
    D = curve.D

    max_attempts = 10 * N
    while res["n_trials"] < N and res["n_attempts"] < max_attempts:
        res["n_attempts"] += 1
        i = rng.randrange(0, q)
        j = rng.randrange(0, q)
        k = rng.randrange(0, q)
        try:
            P = oracle.from_w(oracle.split.scalar_mul(i, W_G))
            Q = oracle.from_w(oracle.split.scalar_mul(j, W_G))
            R = oracle.from_w(oracle.split.scalar_mul(k, W_G))
        except MapUndefined:
            res["n_skipped"] += 1
            continue

        try:
            # Identity: D + P = P and P + D = P
            if not proj_eq(curve.add(D, P), P, curve.p):
                _record(failures, step="10", kind="D + P != P", P=P)
                res["n_failed"] += 1
                continue
            if not proj_eq(curve.add(P, D), P, curve.p):
                _record(failures, step="10", kind="P + D != P", P=P)
                res["n_failed"] += 1
                continue

            # Inverse: P + (-P) = D
            Np = curve.neg(P)
            if not proj_eq(curve.add(P, Np), D, curve.p):
                _record(failures, step="10", kind="P + (-P) != D", P=P)
                res["n_failed"] += 1
                continue

            # Commutativity
            L = curve.add(P, Q)
            Rr = curve.add(Q, P)
            if not proj_eq(L, Rr, curve.p):
                _record(failures, step="10", kind="commutativity",
                        P=P, Q=Q, L=L, R=Rr)
                res["n_failed"] += 1
                continue

            # Associativity
            S1 = curve.add(curve.add(P, Q), R)
            S2 = curve.add(P, curve.add(Q, R))
            if not proj_eq(S1, S2, curve.p):
                _record(failures, step="10", kind="associativity",
                        P=P, Q=Q, R=R, S1=S1, S2=S2)
                res["n_failed"] += 1
                continue

            # Doubling consistency: 2P = P + P
            if not proj_eq(curve.double(P), curve.add(P, P), curve.p):
                _record(failures, step="10", kind="double != add", P=P)
                res["n_failed"] += 1
                continue
        except (ValueError, ZeroDivisionError):
            res["n_skipped"] += 1
            continue

        res["n_trials"] += 1

    return _finalize(res, N)


# ---------------------------------------------------------------------------
# Step 11 — exceptional points with expected statuses
# ---------------------------------------------------------------------------

def step_11_exceptional(curve, failures) -> Dict[str, Any]:
    res = {"name": "step_11_exceptional", "checks": [], "n_failed": 0}

    A, B, C = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    D = curve.D
    Ds = curve.Ds
    Ts = [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1)]

    expected_base_locus = {
        "A": True, "B": True, "C": True,
        "D": False, "D*": False,
        "T0": False, "T1": False, "T2": False, "T3": False,
    }
    named = {
        "A": A, "B": B, "C": C, "D": D, "D*": Ds,
        "T0": Ts[0], "T1": Ts[1], "T2": Ts[2], "T3": Ts[3],
    }

    for name, Pt in named.items():
        on = curve.is_on_curve(Pt)
        res["checks"].append((f"F({name}) = 0", on, True))
        if not on:
            _record(failures, step="11", kind=f"{name} off-curve",
                    Pt=Pt, F=curve.F(Pt))
            res["n_failed"] += 1

        bl = curve.is_base_locus(Pt)
        exp = expected_base_locus[name]
        res["checks"].append((f"{name} in base locus", bl, exp))
        if bl != exp:
            _record(failures, step="11",
                    kind=f"{name} base locus mismatch",
                    expected=exp, got=bl, Pt=Pt)
            res["n_failed"] += 1

    for i, Ti in enumerate(Ts):
        if not proj_eq(curve.isogonal(Ti), Ti, curve.p):
            _record(failures, step="11", kind=f"T{i}* != T{i}", Ti=Ti)
            res["n_failed"] += 1

    res["passed"] = (res["n_failed"] == 0)
    return res


# ---------------------------------------------------------------------------
# Step 12 — deliberately construct A, B, C exceptional inputs
# ---------------------------------------------------------------------------

def _third_intersection(curve, X, Y):
    Ap = curve.dot(X, curve.grad_F(Y))
    Bp = curve.dot(Y, curve.grad_F(X))
    p = curve.p
    return tuple((Ap * X[i] - Bp * Y[i]) % p for i in range(3))


def _abc_case(curve, oracle, W_G, rng, q, X, X_name, N):
    p = curve.p
    counters = {"n_constructions": 0,
                "n_undefined_as_expected": 0,
                "n_unexpectedly_defined": 0,
                "n_failed_construction": 0}
    failures = []

    max_attempts = 20 * N
    attempts = 0
    while (counters["n_constructions"] < N and attempts < max_attempts):
        attempts += 1
        k = rng.randrange(0, q)
        try:
            P = oracle.from_w(oracle.split.scalar_mul(k, W_G))
        except MapUndefined:
            continue
        if curve.is_base_locus(P):
            continue
        if proj_eq(P, X, p):
            continue

        Q = _third_intersection(curve, X, P)
        if is_zero_point(Q, p):
            continue
        if proj_eq(Q, P, p) or proj_eq(Q, X, p):
            continue

        if not curve.is_on_curve(Q):
            failures.append({"kind": f"{X_name} construction: Q off-curve",
                             "P": P, "Q": Q, "F_Q": curve.F(Q)})
            counters["n_failed_construction"] += 1
            continue

        if not are_collinear(X, P, Q, p):
            failures.append({"kind": f"{X_name} construction: not collinear",
                             "X": X, "P": P, "Q": Q})
            counters["n_failed_construction"] += 1
            continue

        R = _third_intersection(curve, P, Q)
        if not proj_eq(R, X, p):
            failures.append({"kind": f"{X_name} third intersection mismatch",
                             "P": P, "Q": Q, "R": R, "X": X})
            counters["n_failed_construction"] += 1
            continue

        try:
            S = curve.add(P, Q)
            counters["n_unexpectedly_defined"] += 1
            failures.append({"kind": f"{X_name} formula unexpectedly returned",
                             "X": X, "P": P, "Q": Q, "S": S})
        except (ValueError, ZeroDivisionError):
            counters["n_undefined_as_expected"] += 1

        counters["n_constructions"] += 1

    counters["failures"] = failures
    counters["X_name"] = X_name
    return counters


def step_12_abc_hunt(curve, G, N, rng, failures, q) -> Dict[str, Any]:
    res = {"name": "step_12_abc_hunt", "n_trials": 0, "n_failed": 0,
           "per_point": {}}
    oracle = Oracle(curve)
    W_G = oracle.to_w(G)

    A = (1, 0, 0)
    B = (0, 1, 0)
    C = (0, 0, 1)

    for X, name in [(A, "A"), (B, "B"), (C, "C")]:
        info = _abc_case(curve, oracle, W_G, rng, q, X, name, N)
        res["per_point"][name] = {
            "n_constructions": info["n_constructions"],
            "n_undefined_as_expected": info["n_undefined_as_expected"],
            "n_unexpectedly_defined": info["n_unexpectedly_defined"],
            "n_failed_construction": info["n_failed_construction"],
        }
        res["n_trials"] += info["n_constructions"]
        res["n_failed"] += info["n_failed_construction"]
        for f in info["failures"]:
            _record(failures, step="12", **f)

    enough = all(res["per_point"][name]["n_constructions"] >= N
                 for name in ("A", "B", "C"))
    ok = enough and (res["n_failed"] == 0)
    res["passed"] = ok
    res["n_requested_per_point"] = N
    if not enough:
        res["reason"] = "could not construct enough exceptional inputs"
    return res


# ---------------------------------------------------------------------------
# Step 13 — special group identities
# ---------------------------------------------------------------------------

def step_13_special_identities(curve, failures) -> Dict[str, Any]:
    res = {"name": "step_13_special_identities", "checks": [], "n_failed": 0}
    Ts = [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1)]
    Ds = curve.Ds

    for i, Ti in enumerate(Ts):
        fixed = proj_eq(curve.isogonal(Ti), Ti, curve.p)
        res["checks"].append((f"iota(T{i}) = T{i}", fixed, True))
        if not fixed:
            res["n_failed"] += 1

    for i, Ti in enumerate(Ts):
        try:
            Ti2 = curve.double(Ti)
        except (ValueError, ZeroDivisionError) as e:
            _record(failures, step="13", kind=f"double(T{i}) undefined",
                    Ti=Ti, err=str(e))
            res["n_failed"] += 1
            continue
        ok = proj_eq(Ti2, Ds, curve.p)
        res["checks"].append((f"2 T{i} = D* (pK)", ok, True))
        if not ok:
            _record(failures, step="13", kind=f"2 T{i} != D*",
                    Ti=Ti, Ti2=Ti2, Ds=Ds)
            res["n_failed"] += 1

    res["passed"] = (res["n_failed"] == 0)
    return res