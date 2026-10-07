#!/usr/bin/env python3
"""Orchestrator for Experiment 1 — pK arithmetic correctness."""

from __future__ import annotations
import argparse
import json
import platform
import random
import sys
import traceback
from pathlib import Path

from pk_core import PKCurve, is_prime, check_curve_params, proj_eq
from pk_oracle import Oracle, MapUndefined
import pk_tests as T


P_PRIME = 106839527430202782610735740077697524141755216002708025221579380756122101340723
Q_ORDER = 26709881857550695652683935019424381035438804000677006305394845189030525335181
U, V, W = 5, 1, 7
G_ALPHA = 1
G_BETA  = 43885198696659819331868805857976446729580581201517994933321256281291666494433
G_GAMMA = 52096118237682591038514212999531378342298973921021687132305677862559506645034


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N-add",     type=int, default=500)
    ap.add_argument("--N-double",  type=int, default=500)
    ap.add_argument("--N-neg",     type=int, default=500)
    ap.add_argument("--N-assoc",   type=int, default=200)
    ap.add_argument("--N-step-45", type=int, default=500)
    ap.add_argument("--N-step-12", type=int, default=50)
    ap.add_argument("--seed",      type=int, default=0)
    ap.add_argument("--outdir",    type=str, default="results")
    return ap.parse_args()


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    curve = PKCurve(U, V, W, P_PRIME)
    G = (G_ALPHA, G_BETA, G_GAMMA)

    try:
        check_curve_params(U, V, W, P_PRIME, Q_ORDER)
    except ValueError as e:
        print(f"[fatal] parameter check failed: {e}")
        return 2

    print(f"[params] p = {P_PRIME}")
    print(f"[params] q = {Q_ORDER}")
    print(f"[params] p - 4q + 1 = {P_PRIME - 4 * Q_ORDER + 1}")
    print(f"[params] F(G) = {curve.F(G)}  (must be 0)")

    failures: list = []
    results: list = []

    def run(name, fn, *a, **kw):
        print(f"\n===== {name} =====")
        try:
            r = fn(*a, **kw)
        except Exception as e:
            traceback.print_exc()
            r = {"name": name, "passed": False,
                 "exception": f"{type(e).__name__}: {e}"}
        results.append(r)
        status = "PASS" if r.get("passed") else "FAIL"
        summary = {k: v for k, v in r.items() if k not in ("checks",)}
        print(f"[{status}] {name}  {summary}")
        if "checks" in r:
            for entry in r["checks"]:
                if len(entry) == 2:
                    label, ok = entry
                    print(f"    - {label}: {'OK' if ok else 'FAIL'}")
                else:
                    label, got, exp = entry
                    ok = (got == exp)
                    print(f"    - {label}: {'OK' if ok else 'FAIL'}"
                          f" (got {got}, expected {exp})")
        return r

    run("step_0_parameters", T.step_0_parameters, curve, Q_ORDER, curve.D, G)

    run("step_4_birational_maps", T.step_4_birational_maps,
        curve, G, args.N_step_45, rng, failures, Q_ORDER)

    run("step_5_generator", T.step_5_generator,
        curve, Q_ORDER, curve.D, G, args.N_step_45, rng, failures)

    run("step_6_isogonal", T.step_6_isogonal,
        curve, G, args.N_step_45, rng, failures, Q_ORDER)

    run("step_7_negation", T.step_7_negation,
        curve, G, args.N_neg, rng, failures, Q_ORDER)

    run("step_8_addition", T.step_8_addition,
        curve, G, args.N_add, rng, failures, Q_ORDER)

    run("step_9_doubling", T.step_9_doubling,
        curve, G, args.N_double, rng, failures, Q_ORDER)

    run("step_10_group_identities", T.step_10_group_identities,
        curve, G, args.N_assoc, rng, failures, Q_ORDER)

    run("step_11_exceptional", T.step_11_exceptional,
        curve, failures)

    run("step_12_abc_hunt", T.step_12_abc_hunt,
        curve, G, args.N_step_12, rng, failures, Q_ORDER)

    run("step_13_special_identities", T.step_13_special_identities,
        curve, failures)

    # ----------------------------------------------------------------
    # Report
    # ----------------------------------------------------------------
    report = {
        "experiment": "EXP-1",
        "seed": args.seed,
        "python_version": sys.version,
        "platform": platform.platform(),
        "parameters": {
            "p": P_PRIME,
            "q": Q_ORDER,
            "u": U, "v": V, "w": W,
            "p_is_prime": is_prime(P_PRIME),
            "q_is_prime": is_prime(Q_ORDER),
            "p_eq_4q_minus_1": (P_PRIME == 4 * Q_ORDER - 1),
            "F_of_G": curve.F(G),
        },
        "requested_trials": {
            "N_add": args.N_add,
            "N_double": args.N_double,
            "N_neg": args.N_neg,
            "N_assoc": args.N_assoc,
            "N_step_45": args.N_step_45,
            "N_step_12": args.N_step_12,
        },
        "results": results,
        "n_failures": len(failures),
    }
    (outdir / "report.json").write_text(
        json.dumps(report, indent=2, default=str))

    with open(outdir / "failures.jsonl", "w") as f:
        for rec in failures:
            f.write(json.dumps(rec, default=str) + "\n")

    lines = ["Experiment 1 - pK arithmetic correctness",
             f"seed = {args.seed}", ""]
    overall = True
    for r in results:
        status = "PASS" if r.get("passed") else "FAIL"
        if not r.get("passed"):
            overall = False
        extra = ""
        if "n_requested" in r and "n_trials" in r:
            extra = (f" requested={r['n_requested']}"
                     f" attempts={r.get('n_attempts', '-')}"
                     f" successful={r['n_trials']}"
                     f" skipped={r.get('n_skipped', '-')}"
                     f" failed={r.get('n_failed', '-')}")
        lines.append(f"  [{status}] {r['name']}{extra}")
    lines.append("")
    lines.append(f"Total failures recorded: {len(failures)}")
    lines.append(f"Overall: {'PASS' if overall else 'FAIL'}")
    (outdir / "report.txt").write_text("\n".join(lines))

    print("\n===== SUMMARY =====")
    print("\n".join(lines))
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())