"""
EXP-2 / EXP-6 independent SageMath certification.

Loads results/basis_23.json and verifies that each extracted law satisfies

    S_1 * Z_0 - S_0 * Z_1  ∈  (F(P), F(Q))
    S_2 * Z_0 - S_0 * Z_2  ∈  (F(P), F(Q))

by exact ideal membership. The reduction engine is Sage's own, which is
independent of the memoized substitution used by the Python extractor.

Usage:
    sage sage_verify.sage

Output:
    results/verification_sage.txt
"""

import json
import sys
from pathlib import Path

from sage.all import GF, PolynomialRing

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"

P_MOD = 106839527430202782610735740077697524141755216002708025221579380756122101340723
U, V, W = 5, 1, 7

R = PolynomialRing(GF(P_MOD), names=("a", "b", "c", "x", "y", "z"))
a, b, c, x, y, z = R.gens()

F_P = U * a * (b**2 - c**2) + V * b * (c**2 - a**2) + W * c * (a**2 - b**2)
F_Q = U * x * (y**2 - z**2) + V * y * (z**2 - x**2) + W * z * (x**2 - y**2)
I = R.ideal([F_P, F_Q])

A = a * F_Q.derivative(x) + b * F_Q.derivative(y) + c * F_Q.derivative(z)
B = x * F_P.derivative(a) + y * F_P.derivative(b) + z * F_P.derivative(c)
R_0 = A * a - B * x
R_1 = A * b - B * y
R_2 = A * c - B * z
S_0 = R_1 * R_2
S_1 = R_2 * R_0
S_2 = R_0 * R_1

Q_P = [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)]
C_Q = [(0, 0, 3), (0, 1, 2), (0, 2, 1), (0, 3, 0),
       (1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (3, 0, 0)]


def mono_to_poly(m, vars):
    result = R(1)
    for i, e in enumerate(m):
        if e:
            result *= vars[i] ** e
    return result


Q_P_polys = [mono_to_poly(m, (a, b, c)) for m in Q_P]
C_Q_polys = [mono_to_poly(m, (x, y, z)) for m in C_Q]


def main():
    with open(RESULTS / "basis_23.json") as f:
        data = json.load(f)

    lines = ["Sage independent verification of EXP-2/6 (2,3) basis",
             "=" * 60]
    all_ok = True
    for law in data["laws"]:
        j = law["index"]
        coeffs = law["coeffs"]
        Z = [R(0), R(0), R(0)]
        for i in range(3):
            for mu_idx in range(6):
                for nu_idx in range(9):
                    cc = int(coeffs[i * 54 + mu_idx * 9 + nu_idx]) % P_MOD
                    if cc:
                        Z[i] += cc * Q_P_polys[mu_idx] * C_Q_polys[nu_idx]

        G1 = S_1 * Z[0] - S_0 * Z[1]
        G2 = S_2 * Z[0] - S_0 * Z[2]
        ok1 = (G1 in I)
        ok2 = (G2 in I)
        ok = ok1 and ok2
        all_ok = all_ok and ok

        line = (f"law {j:02d}:  G1 in I = {ok1},  G2 in I = {ok2},  "
                f"overall = {ok}")
        print(line)
        lines.append(line)

    lines.append("=" * 60)
    lines.append(f"all certified: {all_ok}")
    (RESULTS / "verification_sage.txt").write_text("\n".join(lines) + "\n")

    print()
    print("\n".join(lines))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())