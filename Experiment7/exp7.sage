"""
EXP-7 — incremental projective completeness search for A_{2,3}.

Runs one phase per invocation.  State is written after every test and
resumed automatically on the next invocation of the same phase, including
retries of tests that were previously indeterminate (timeout or error).

Phases:
    L1_r{r}   — Level 1: exhaustive over r-subsets of the nine basis laws.
    L2_r{r}   — Level 2: sampled r-tuples of points in P^8(F_p).

Usage:
    sage exp7.sage                     show help and phase status
    sage exp7.sage <phase>             run or resume a phase
    sage exp7.sage summary             print overall summary
    sage exp7.sage --list              list phases and status
    sage exp7.sage --reset-phase <p>   clear one phase's results
    sage exp7.sage --reset-all         clear all phases
"""

import gc
import json
import os
import random
import signal
import sys
import tempfile
import time
from itertools import combinations
from pathlib import Path

from sage.all import GF, PolynomialRing, matrix


# ====================================================================
# Configuration
# ====================================================================
HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
RESULTS_DIR.mkdir(exist_ok=True)
REPORT_JSON = RESULTS_DIR / "report.json"
REPORT_TXT = RESULTS_DIR / "report.txt"

P_MOD = int("106839527430202782610735740077697524141755216002708025221579380756122101340723")
U = int(5)
V = int(1)
W = int(7)

LEVEL1_RS           = [2, 3, 4]
LEVEL2_RS           = [2, 3]
LEVEL2_TRIALS       = 30
PER_CHART_TIMEOUT_S = 60
RNG_SEED            = 0

SCHEMA = "exp7_v3"    # bumped: pending predicate now retries indeterminates


# ====================================================================
# JSON helpers
# ====================================================================
def _json_safe(o):
    if isinstance(o, (int, float, str, bool)) or o is None:
        return o
    try:
        return int(o)
    except (TypeError, ValueError):
        return str(o)


def _atomic_write_json(path, obj):
    text = json.dumps(obj, indent=2, default=_json_safe)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp_report_", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
        os.replace(tmp, str(path))
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


# ====================================================================
# Phase table
# ====================================================================
def all_phase_names():
    names = []
    for r in LEVEL1_RS:
        names.append(f"L1_r{r}")
    for r in LEVEL2_RS:
        names.append(f"L2_r{r}")
    return names


def phase_spec(name):
    if name.startswith("L1_r"):
        r = int(name[4:])
        if r in LEVEL1_RS:
            return ("L1", r)
    if name.startswith("L2_r"):
        r = int(name[4:])
        if r in LEVEL2_RS:
            return ("L2", r)
    return None


# ====================================================================
# Basis loading
# ====================================================================
def find_basis_file():
    candidates = [
        HERE / ".." / "Experiment2_6" / "results" / "basis_23.json",
        HERE / ".." / ".." / "Experiment2_6" / "results" / "basis_23.json",
        HERE / "basis_23.json",
    ]
    for c in candidates:
        c = c.resolve()
        if c.exists():
            return c
    raise FileNotFoundError(
        "Could not locate Experiment2_6/results/basis_23.json; "
        "run EXP-2/6 first."
    )


basis_path = find_basis_file()
print(f"[setup] basis file: {basis_path}")
with open(basis_path) as f:
    basis_data = json.load(f)

Q_P = basis_data["Q_P"]
C_Q = basis_data["C_Q"]
LAWS_COEFFS = [[int(x) for x in law["coeffs"]] for law in basis_data["laws"]]
assert len(LAWS_COEFFS) == 9, f"expected 9 laws, got {len(LAWS_COEFFS)}"
for k, c in enumerate(LAWS_COEFFS):
    assert len(c) == 162, f"law {k} has {len(c)} coeffs, expected 162"

F_p = GF(P_MOD)
R = PolynomialRing(F_p, names=("a", "b", "c", "x", "y", "z"))
a, b, c, x, y, z = R.gens()

rank_coeffs = int(matrix(F_p, LAWS_COEFFS).rank())
assert rank_coeffs == 9, (
    f"basis_23.json has rank {rank_coeffs} over F_p, expected 9"
)
print(f"[setup] basis rank: {rank_coeffs}")

F_P = U * a * (b ** 2 - c ** 2) + V * b * (c ** 2 - a ** 2) + W * c * (a ** 2 - b ** 2)
F_Q = U * x * (y ** 2 - z ** 2) + V * y * (z ** 2 - x ** 2) + W * z * (x ** 2 - y ** 2)


def monomial_poly(exponents, gens):
    result = R(1)
    for i, e in enumerate(exponents):
        if e:
            result *= gens[i] ** e
    return result


Q_P_polys = [monomial_poly(m, (a, b, c)) for m in Q_P]
C_Q_polys = [monomial_poly(m, (x, y, z)) for m in C_Q]


def law_polys_from_coeffs(coeffs):
    out = []
    for i in range(3):
        Zi = R(0)
        base = i * 54
        for mu_idx in range(6):
            qm = Q_P_polys[mu_idx]
            sub = base + mu_idx * 9
            for nu_idx in range(9):
                cc = int(coeffs[sub + nu_idx]) % P_MOD
                if cc:
                    Zi += cc * qm * C_Q_polys[nu_idx]
        out.append(Zi)
    return out


print("[setup] reconstructing the nine law polynomial triples ...")
LAWS = [law_polys_from_coeffs(coeffs) for coeffs in LAWS_COEFFS]


# ====================================================================
# Nine-chart data
# ====================================================================
def build_charts():
    names6 = ["a", "b", "c", "x", "y", "z"]
    charts = {}
    for p_idx in (0, 1, 2):
        for q_idx in (3, 4, 5):
            kept = [names6[i] for i in range(6) if i not in (p_idx, q_idx)]
            R4 = PolynomialRing(F_p, names=kept)
            images = []
            for i in range(6):
                if i in (p_idx, q_idx):
                    images.append(R4(1))
                else:
                    images.append(R4.gen(kept.index(names6[i])))
            phi = R.hom(images, R4)
            charts[(p_idx, q_idx)] = (R4, phi, phi(F_P), phi(F_Q))
    return charts


CHART_KEYS = [(p, q) for p in (0, 1, 2) for q in (3, 4, 5)]
CHARTS = build_charts()


# ====================================================================
# Per-chart soft timeout
# ====================================================================
class ChartTimeout(Exception):
    pass


def _chart_timeout_handler(signum, frame):
    raise ChartTimeout()


def chart_is_unit_ideal(key, restrictions_for_key, timeout_s):
    R4, _, FP_c, FQ_c = CHARTS[key]
    gens = [FP_c, FQ_c] + restrictions_for_key
    I = R4.ideal(gens)
    old = signal.signal(signal.SIGALRM, _chart_timeout_handler)
    signal.alarm(int(timeout_s))
    try:
        gb = I.groebner_basis()
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
        unit = any(g.is_constant() for g in gb) if gb else False
        del gb, I, gens
        return (unit, None)
    except ChartTimeout:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
        del I, gens
        return (None, "timeout")
    except Exception as exc:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
        del I, gens
        return (None, f"error: {exc}")


def law_polys_to_chart_restrictions(law_polys_6):
    restrictions = {}
    for key in CHART_KEYS:
        R4, phi, _, _ = CHARTS[key]
        restr = []
        for zs in law_polys_6:
            for zz in zs:
                restr.append(phi(zz))
        restrictions[key] = restr
    return restrictions


def is_complete(law_polys_6, per_chart_timeout):
    restrictions = law_polys_to_chart_restrictions(law_polys_6)
    for key in CHART_KEYS:
        unit, note = chart_is_unit_ideal(
            key, restrictions[key], per_chart_timeout)
        if unit is None:
            del restrictions
            return (None, note)
        if unit is False:
            del restrictions
            return (False, None)
    del restrictions
    return (True, None)


# ====================================================================
# Test builders
# ====================================================================
def random_point_P8(rng, p_mod):
    while True:
        v = [int(rng.randrange(p_mod)) for _ in range(9)]
        if any(x != 0 for x in v):
            return v


def build_level1_tests(r):
    out = []
    for idx, subset in enumerate(combinations(range(9), r)):
        out.append({
            "id": f"L1_r{r}_s{idx:04d}",
            "kind": "indices",
            "indices": [int(j) for j in subset],
            "r": int(r),
        })
    return out


def build_level2_tests(r, trials):
    # Deterministic per r, so the same lambdas are regenerated on resume.
    rng = random.Random(RNG_SEED + 1000 * r)
    out = []
    for trial in range(trials):
        lambdas = [random_point_P8(rng, P_MOD) for _ in range(r)]
        out.append({
            "id": f"L2_r{r}_t{trial:04d}",
            "kind": "combination",
            "lambdas": lambdas,
            "r": int(r),
        })
    return out


def tests_for_phase(name):
    spec = phase_spec(name)
    if spec is None:
        return None
    kind, r = spec
    if kind == "L1":
        return build_level1_tests(r)
    else:
        return build_level2_tests(r, LEVEL2_TRIALS)


# ====================================================================
# Test evaluation
# ====================================================================
def evaluate_test(t):
    kind = t["kind"]

    if kind == "indices":
        law_polys_6 = [LAWS[j] for j in t["indices"]]

    elif kind == "combination":
        # Direct polynomial-level combination of the nine basis laws.
        # Mathematically equivalent to combining the 162-coefficient
        # vectors, but avoids rebuilding a 54-monomial polynomial
        # for every sampled law.
        law_polys_6 = []
        for lam in t["lambdas"]:
            combined = []
            for j in range(3):
                Z = R(0)
                for k in range(9):
                    lk = int(lam[k]) % P_MOD
                    if lk:
                        Z += lk * LAWS[k][j]
                combined.append(Z)
            law_polys_6.append(combined)

    else:
        return {"id": t["id"], "complete": None,
                "error": f"unknown kind {kind}"}

    verdict, note = is_complete(law_polys_6, PER_CHART_TIMEOUT_S)
    result = {"id": t["id"], "complete": verdict}
    if note is not None:
        result["error"] = note
    del law_polys_6
    return result


# ====================================================================
# State
# ====================================================================
def empty_state():
    return {
        "schema": SCHEMA,
        "p_mod": int(P_MOD),
        "U": int(U), "V": int(V), "W": int(W),
        "rng_seed": int(RNG_SEED),
        "level1_rs": [int(r) for r in LEVEL1_RS],
        "level2_rs": [int(r) for r in LEVEL2_RS],
        "level2_trials": int(LEVEL2_TRIALS),
        "per_chart_timeout_s": int(PER_CHART_TIMEOUT_S),
        "phases": {name: {"results_by_id": {}}
                   for name in all_phase_names()},
        "last_updated_s": 0.0,
    }


def load_state():
    if REPORT_JSON.exists():
        try:
            s = json.loads(REPORT_JSON.read_text())
        except Exception as exc:
            print(f"[state] could not read {REPORT_JSON}: {exc}")
            return empty_state()
        if s.get("schema") != SCHEMA:
            print(f"[state] schema mismatch; starting fresh")
            return empty_state()
        if s.get("p_mod") != int(P_MOD) or \
           s.get("U") != int(U) or \
           s.get("V") != int(V) or \
           s.get("W") != int(W) or \
           s.get("rng_seed") != int(RNG_SEED) or \
           s.get("level1_rs") != [int(r) for r in LEVEL1_RS] or \
           s.get("level2_rs") != [int(r) for r in LEVEL2_RS] or \
           s.get("level2_trials") != int(LEVEL2_TRIALS):
            print("[state] configuration differs from stored state; "
                  "starting fresh")
            return empty_state()
        for name in all_phase_names():
            s.setdefault("phases", {}).setdefault(
                name, {"results_by_id": {}})
        return s
    return empty_state()


def save_state(state, t0):
    state["last_updated_s"] = float(time.time() - t0)
    _atomic_write_json(REPORT_JSON, state)


# ====================================================================
# Phase status / summary
# ====================================================================
def phase_status(state, name):
    spec = phase_spec(name)
    if spec is None:
        return None
    tests = tests_for_phase(name)
    if tests is None:
        return None
    n_total = len(tests)
    results = state["phases"].get(name, {}).get("results_by_id", {})
    n_stored = len(results)
    n_determinate = sum(
        1 for t in tests
        if results.get(t["id"], {}).get("complete") in (True, False)
    )
    n_indeterminate = n_stored - n_determinate
    if n_determinate == n_total:
        status = "done"
    elif n_stored == 0:
        status = "pending"
    else:
        status = "partial"
    return {
        "name": name,
        "status": status,
        "n_total": n_total,
        "n_stored": n_stored,
        "n_determinate": n_determinate,
        "n_indeterminate": n_indeterminate,
        "n_complete": sum(
            1 for t in tests
            if results.get(t["id"], {}).get("complete") is True
        ),
        "n_incomplete": sum(
            1 for t in tests
            if results.get(t["id"], {}).get("complete") is False
        ),
    }


def print_help_and_status(state):
    print("=" * 72)
    print("EXP-7  —  incremental completeness search for A_{2,3}")
    print("=" * 72)
    print(f"p_mod: {P_MOD.bit_length()} bits")
    print(f"per_chart_timeout: {PER_CHART_TIMEOUT_S} s")
    print("")
    print("Available phases:")
    for name in all_phase_names():
        s = phase_status(state, name)
        print(f"  {name:<10}  [{s['status']:<7}]  "
              f"{s['n_determinate']}/{s['n_total']} determinate, "
              f"{s['n_indeterminate']} indeterminate, "
              f"{s['n_complete']} complete")
    print("")
    print("Usage:")
    print("  sage exp7.sage <phase>            run/resume phase")
    print("  sage exp7.sage summary            print overall summary")
    print("  sage exp7.sage --list             list phases")
    print("  sage exp7.sage --reset-phase <p>  clear one phase")
    print("  sage exp7.sage --reset-all        clear all phases")


def print_summary(state):
    print("=" * 72)
    print("EXP-7  —  overall summary")
    print("=" * 72)
    for name in all_phase_names():
        s = phase_status(state, name)
        print(f"  {name:<10}  [{s['status']:<7}]  "
              f"{s['n_determinate']}/{s['n_total']} determinate, "
              f"{s['n_indeterminate']} indeterminate, "
              f"{s['n_complete']} complete")
    r1_complete = []
    for r in LEVEL1_RS:
        s = phase_status(state, f"L1_r{r}")
        if s["n_complete"] > 0:
            r1_complete.append(r)
    if r1_complete:
        print(f"\n  smallest r with complete L1 subset: {min(r1_complete)}")
    r2_complete = []
    for r in LEVEL2_RS:
        s = phase_status(state, f"L2_r{r}")
        if s["n_complete"] > 0:
            r2_complete.append(r)
    if r2_complete:
        print(f"  smallest r with complete L2 sample: {min(r2_complete)}")
    print("")
    print("Minimality note: any complete pair is globally minimal "
          "(Arene-Kohel-Ritzenthaler, g=1).")
    lines = ["EXP-7 overall summary"]
    for name in all_phase_names():
        s = phase_status(state, name)
        lines.append(f"  {name}: status={s['status']} "
                     f"complete={s['n_complete']} "
                     f"determinate={s['n_determinate']}/{s['n_total']} "
                     f"indeterminate={s['n_indeterminate']}")
    REPORT_TXT.write_text("\n".join(lines) + "\n")


# ====================================================================
# Phase runner
# ====================================================================
def run_phase(name, state):
    spec = phase_spec(name)
    if spec is None:
        print(f"[error] unknown phase: {name}")
        return 1
    tests = tests_for_phase(name)
    if tests is None:
        print(f"[error] could not build tests for {name}")
        return 1

    results = state["phases"][name]["results_by_id"]

    # Pending = missing OR previously indeterminate (complete not in
    # (True, False)).  This makes timeouts and errors retryable on the
    # next invocation of the same phase.
    pending = [
        t for t in tests
        if results.get(t["id"], {}).get("complete") not in (True, False)
    ]

    print(f"[phase {name}] {len(pending)} pending/indeterminate "
          f"of {len(tests)}")
    if not pending:
        print(f"[phase {name}] already complete")
        return 0

    t0 = time.time()
    for cnt, t in enumerate(pending):
        try:
            r = evaluate_test(t)
        except Exception as exc:
            r = {"id": t["id"], "complete": None, "error": str(exc)}
        results[t["id"]] = r
        del t
        gc.collect()
        # Save after every test (bounded loss on interruption).
        save_state(state, t0)
        if (cnt + 1) % 5 == 0 or (cnt + 1) == len(pending):
            dt = time.time() - t0
            rate = (cnt + 1) / dt if dt > 0 else 0.0
            print(f"  [{name}] {cnt + 1}/{len(pending)} done in "
                  f"{dt:.1f}s ({rate:.2f} tests/s)", flush=True)

    s = phase_status(state, name)
    print(f"[phase {name}] status={s['status']} "
          f"complete={s['n_complete']} "
          f"determinate={s['n_determinate']}/{s['n_total']} "
          f"indeterminate={s['n_indeterminate']}")
    return 0


# ====================================================================
# Main
# ====================================================================
def main():
    args = sys.argv[1:]
    state = load_state()

    if not args:
        print_help_and_status(state)
        return 0

    if args[0] == "--list":
        print_help_and_status(state)
        return 0

    if args[0] == "--reset-all":
        state = empty_state()
        save_state(state, time.time())
        print("[reset] all phases cleared")
        return 0

    if args[0] == "--reset-phase":
        if len(args) < 2:
            print("usage: sage exp7.sage --reset-phase <name>")
            return 1
        name = args[1]
        if name not in state["phases"]:
            print(f"[error] unknown phase: {name}")
            return 1
        state["phases"][name]["results_by_id"] = {}
        save_state(state, time.time())
        print(f"[reset] {name} cleared")
        return 0

    if args[0] == "summary":
        print_summary(state)
        return 0

    name = args[0]
    if name not in state["phases"]:
        print(f"[error] unknown phase: {name}")
        print_help_and_status(state)
        return 1
    t0 = time.time()
    rc = run_phase(name, state)
    save_state(state, t0)
    return rc


if __name__ == "__main__":
    sys.exit(main())