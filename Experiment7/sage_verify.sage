"""
Independent re-verification of EXP-7 reported complete examples.

Loads results/report.json, checks the reported basis and its rank,
extracts every reported complete example, and re-runs the nine-chart
projective completeness test from a fresh implementation.

Writes results/verification_sage.txt.
"""

import json
import signal
import sys
from pathlib import Path

from sage.all import GF, PolynomialRing, matrix


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


class _Timeout(Exception):
    pass


def _handler(signum, frame):
    raise _Timeout()


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
    raise FileNotFoundError("Could not locate Experiment2_6/results/basis_23.json")


class ProjectiveTester:
    def __init__(self, basis_path, p_mod, U, V, W, per_chart_timeout=120):
        with open(basis_path) as f:
            data = json.load(f)
        self.Q_P = data["Q_P"]
        self.C_Q = data["C_Q"]
        self.laws_coeffs = [[int(x) for x in law["coeffs"]]
                            for law in data["laws"]]

        # ---- basis / rank assertions -------------------------------
        assert len(self.laws_coeffs) == 9, (
            f"expected 9 basis laws, got {len(self.laws_coeffs)}"
        )
        for k, c in enumerate(self.laws_coeffs):
            assert len(c) == 162, (
                f"law {k} has {len(c)} coefficients, expected 162"
            )
        self.p_mod = int(p_mod)
        F_p = GF(self.p_mod)
        M = matrix(F_p, self.laws_coeffs)
        r = int(M.rank())
        assert r == 9, (
            f"basis_23.json has rank {r} over F_p, expected 9"
        )
        print(f"[verifier] basis shape: {M.nrows()} x {M.ncols()}, rank {r}")

        self.F_p = F_p
        self.R = PolynomialRing(F_p, names=("a", "b", "c", "x", "y", "z"))
        gens = self.R.gens()
        a, b, c, x, y, z = gens
        self.F_P = U * a * (b ** 2 - c ** 2) + V * b * (c ** 2 - a ** 2) + W * c * (a ** 2 - b ** 2)
        self.F_Q = U * x * (y ** 2 - z ** 2) + V * y * (z ** 2 - x ** 2) + W * z * (x ** 2 - y ** 2)
        self.Q_P_polys = [self._monomial(m, gens[:3]) for m in self.Q_P]
        self.C_Q_polys = [self._monomial(m, gens[3:]) for m in self.C_Q]
        self.per_chart_timeout = int(per_chart_timeout)

    def _monomial(self, exponents, gens):
        result = self.R(1)
        for i, e in enumerate(exponents):
            if e:
                result *= gens[i] ** e
        return result

    def _law_polys(self, coeffs):
        out = []
        for i in range(3):
            Zi = self.R(0)
            base = i * 54
            for mu_idx in range(6):
                qm = self.Q_P_polys[mu_idx]
                sub = base + mu_idx * 9
                for nu_idx in range(9):
                    cc = int(coeffs[sub + nu_idx]) % self.p_mod
                    if cc:
                        Zi += cc * qm * self.C_Q_polys[nu_idx]
            out.append(Zi)
        return out

    def _dehomogenize(self, poly, p_idx, q_idx):
        names6 = ["a", "b", "c", "x", "y", "z"]
        kept = [names6[i] for i in range(6) if i not in (p_idx, q_idx)]
        R4 = PolynomialRing(self.F_p, names=kept)
        images = []
        for i in range(6):
            if i in (p_idx, q_idx):
                images.append(R4(1))
            else:
                images.append(R4.gen(kept.index(names6[i])))
        phi = self.R.hom(images, R4)
        return R4, phi(poly)

    def is_complete(self, coeffs_list):
        law_polys = [self._law_polys(c) for c in coeffs_list]
        for p_idx in (0, 1, 2):
            for q_idx in (3, 4, 5):
                R4, FP_c = self._dehomogenize(self.F_P, p_idx, q_idx)
                _, FQ_c = self._dehomogenize(self.F_Q, p_idx, q_idx)
                gens = [FP_c, FQ_c]
                for zs in law_polys:
                    for z in zs:
                        _, z4 = self._dehomogenize(z, p_idx, q_idx)
                        gens.append(z4)
                I = R4.ideal(gens)
                old = signal.signal(signal.SIGALRM, _handler)
                signal.alarm(self.per_chart_timeout)
                try:
                    gb = I.groebner_basis()
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, old)
                except _Timeout:
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, old)
                    return None
                if not (any(g.is_constant() for g in gb) if gb else False):
                    return False
        return True


def main():
    report_path = RESULTS / "report.json"
    if not report_path.exists():
        print(f"[error] {report_path} not found; run exp7.sage first.")
        return 1
    report = json.loads(report_path.read_text())

    basis_path = find_basis_file()
    tester = ProjectiveTester(
        basis_path, report["p_mod"], report["U"], report["V"], report["W"],
        per_chart_timeout=120)

    lines = ["Independent projective re-verification of EXP-7",
             "=" * 60]
    all_ok = True

    for level_key, label in (("level1", "L1"), ("level2", "L2")):
        block = report.get(level_key, {})
        for t in block.get("tests_with_results", []):
            if t.get("complete") is not True:
                continue
            tid = t["id"]
            if t["kind"] == "indices":
                coeffs_list = [tester.laws_coeffs[i] for i in t["indices"]]
            else:
                coeffs_list = []
                for lam in t["lambdas"]:
                    new_coeffs = [0] * 162
                    for k in range(9):
                        lk = int(lam[k]) % report["p_mod"]
                        if lk == 0:
                            continue
                        old = tester.laws_coeffs[k]
                        for i in range(162):
                            new_coeffs[i] = (
                                new_coeffs[i] + lk * old[i]
                            ) % report["p_mod"]
                    coeffs_list.append(new_coeffs)
            ok = tester.is_complete(coeffs_list)
            all_ok = all_ok and (ok is True)
            status = "OK" if ok is True else ("TIMEOUT" if ok is None else "FAIL")
            line = f"  {label} {tid}: {status}"
            print(line)
            lines.append(line)

    lines.append("=" * 60)
    lines.append(f"all verified: {all_ok}")
    (RESULTS / "verification_sage.txt").write_text("\n".join(lines) + "\n")
    print()
    print("\n".join(lines))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())