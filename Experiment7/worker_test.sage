"""
EXP-7 subprocess worker.

Reads a JSON task from stdin and writes a JSON result to stdout.
For each test the worker runs the nine-chart projective completeness
test on P^2 x P^2. The test is complete iff all nine chart ideals are
the unit ideal.
"""

import json
import signal
import sys

from sage.all import GF, PolynomialRing


class ChartTimeout(Exception):
    pass


def _chart_timeout_handler(signum, frame):
    raise ChartTimeout()


def build_ring_and_basis(basis_path, p_mod, U, V, W):
    with open(basis_path) as f:
        data = json.load(f)
    Q_P = data["Q_P"]
    C_Q = data["C_Q"]
    laws_coeffs = [law["coeffs"] for law in data["laws"]]

    F_p = GF(p_mod)
    R = PolynomialRing(F_p, names=("a", "b", "c", "x", "y", "z"))
    a, b, c, x, y, z = R.gens()
    F_P = U * a * (b ** 2 - c ** 2) + V * b * (c ** 2 - a ** 2) + W * c * (a ** 2 - b ** 2)
    F_Q = U * x * (y ** 2 - z ** 2) + V * y * (z ** 2 - x ** 2) + W * z * (x ** 2 - y ** 2)
    return R, F_P, F_Q, Q_P, C_Q, laws_coeffs


def monomial_poly(R, exponents, gens):
    result = R(1)
    for i, e in enumerate(exponents):
        if e:
            result *= gens[i] ** e
    return result


def law_polys_from_coeffs(R, Q_P_polys, C_Q_polys, coeffs, p_mod):
    out = []
    for i in range(3):
        Zi = R(0)
        base = i * 54
        for mu_idx in range(6):
            qm = Q_P_polys[mu_idx]
            sub = base + mu_idx * 9
            for nu_idx in range(9):
                cc = int(coeffs[sub + nu_idx]) % p_mod
                if cc:
                    Zi += cc * qm * C_Q_polys[nu_idx]
        out.append(Zi)
    return out


def make_chart_data(R, F_p, F_P, F_Q):
    """Return charts[(p_idx, q_idx)] = (R4, phi, phi(F_P), phi(F_Q))."""
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


def coeffs_from_test(task, laws_coeffs, p_mod):
    kind = task["kind"]
    if kind == "indices":
        return [laws_coeffs[i] for i in task["indices"]]
    if kind == "combination":
        out = []
        for lam in task["lambdas"]:
            new_coeffs = [0] * 162
            for k in range(9):
                lk = int(lam[k]) % p_mod
                if lk == 0:
                    continue
                old = laws_coeffs[k]
                for i in range(162):
                    new_coeffs[i] = (new_coeffs[i] + lk * int(old[i])) % p_mod
            out.append(new_coeffs)
        return out
    raise ValueError(f"unknown test kind: {kind}")


def test_is_complete(coeffs_list, R, F_p, Q_P_polys, C_Q_polys,
                     charts, per_chart_timeout):
    """Nine-chart projective completeness test.

    Returns (verdict, error), verdict in {True, False, None}:
      True  : all nine charts yield the unit ideal  -> complete
      False : at least one chart is not the unit ideal
      None  : a chart timed out or errored
    """
    law_polys_6 = [
        law_polys_from_coeffs(R, Q_P_polys, C_Q_polys, c, F_p.characteristic())
        for c in coeffs_list
    ]
    for (p_idx, q_idx), (R4, phi, FP_c, FQ_c) in charts.items():
        law_polys_4v = []
        for zs in law_polys_6:
            for z in zs:
                law_polys_4v.append(phi(z))
        gens = [FP_c, FQ_c] + law_polys_4v
        I = R4.ideal(gens)

        old = signal.signal(signal.SIGALRM, _chart_timeout_handler)
        signal.alarm(per_chart_timeout)
        try:
            gb = I.groebner_basis()
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
            is_unit = any(g.is_constant() for g in gb) if gb else False
        except ChartTimeout:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
            return (None, f"chart ({p_idx},{q_idx}) timeout")
        except Exception as exc:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old)
            return (None, f"chart ({p_idx},{q_idx}) error: {exc}")
        if not is_unit:
            return (False, None)
    return (True, None)


def main():
    task = json.load(sys.stdin)
    basis_path = task["basis_path"]
    p_mod = int(task["p_mod"])
    U = int(task["U"]); V = int(task["V"]); W = int(task["W"])
    per_chart_timeout = int(task.get("per_chart_timeout_s", 60))
    tests = task["tests"]

    R, F_P, F_Q, Q_P, C_Q, laws_coeffs = build_ring_and_basis(
        basis_path, p_mod, U, V, W)
    F_p = R.base_ring()
    gens = R.gens()
    Q_P_polys = [monomial_poly(R, m, gens[:3]) for m in Q_P]
    C_Q_polys = [monomial_poly(R, m, gens[3:]) for m in C_Q]
    charts = make_chart_data(R, F_p, F_P, F_Q)

    results = []
    for t in tests:
        tid = t["id"]
        try:
            coeffs_list = coeffs_from_test(t, laws_coeffs, p_mod)
            verdict, err = test_is_complete(
                coeffs_list, R, F_p, Q_P_polys, C_Q_polys,
                charts, per_chart_timeout)
            results.append({"id": tid, "complete": verdict, "error": err})
        except Exception as exc:
            results.append({"id": tid, "complete": None, "error": str(exc)})

    json.dump({"results": results}, sys.stdout)


if __name__ == "__main__":
    main()