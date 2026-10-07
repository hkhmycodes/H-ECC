#!/usr/bin/env python3
"""
EXP-8 (revised) -- static CSE-aware cost analysis, correctness verification,
and wall-clock benchmark of the seven complete (2,3) addition-law triples
identified by EXP-7.

Pure Python (no Sage preparser).

This revision applies:
  1. inverse_pair now computes phi_0(Q) = 2*W_D - phi_0(P), not W_D - phi_0(P).
  2. Inverse-pair checking has its own counter and threshold (MIN_INV_OK).
  3. verify pre-checks that the E_W oracle agrees with the composite
     pK formula on random pairs, before testing any triple.
  4. bench refuses to run unless verify has exited with zero failures.
  5. Comment about "order-q subgroup" corrected; formula left as-is.
  6. bench output notes the timing bias against the (2,3) side.
"""

import gc
import json
import os
import random
import re
import sys
import tempfile
import time
from itertools import combinations
from pathlib import Path
from statistics import median


# =====================================================================
# Configuration
# =====================================================================
P_MOD = int(
    "106839527430202782610735740077697524141755216002708025221579380756122101340723"
)
Q_ORDER = int(
    "26709881857550695652683935019424381035438804000677006305394845189030525335181"
)
U, V, W = 5, 1, 7

INV36 = pow(36, -1, P_MOD)
A_W = (-INV36) % P_MOD
B_W = 0

N_TRIALS = 1000
N_WARMUP = 100
N_PASSES = 5
SEEDS    = [0, 1, 2, 3, 4]

HEUR_WEIGHTS = {
    "M_general": 1.00,
    "C_full":    1.00,
    "C_small":   0.10,
    "S":         0.80,
    "A":         0.05,
}

SMALL_CONST_BITS = 16

# Verify thresholds (minimum successful comparisons per triple)
MIN_ADD_OK = 30
MIN_DBL_OK = 15
MIN_INV_OK = 10       # FIX 2: inverse-pair threshold


# =====================================================================
# JSON / path helpers
# =====================================================================
def atomic_write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp_", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=2, default=str)
        os.replace(tmp, str(path))
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def find_file(candidates, name):
    for c in candidates:
        c = c.resolve()
        if c.exists():
            return c
    raise FileNotFoundError(
        f"Could not locate {name}; tried {[str(c) for c in candidates]}"
    )


HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
RESULTS_DIR.mkdir(exist_ok=True)

ANALYSIS_JSON = RESULTS_DIR / "exp8_analysis.json"
BENCH_JSON    = RESULTS_DIR / "exp8_bench.json"
VERIFY_JSON   = RESULTS_DIR / "exp8_verify.json"
PROFILE_JSON  = RESULTS_DIR / "exp8_profile.json"


# =====================================================================
# Load basis and EXP-7 report
# =====================================================================
print(f"[setup] python {sys.version.split()[0]}, "
      f"p = {P_MOD.bit_length()} bits")

_basis_path = find_file([
    HERE / ".." / "Experiment2_6" / "results" / "basis_23.json",
    HERE / ".." / ".." / "Experiment2_6" / "results" / "basis_23.json",
    HERE / "basis_23.json",
], "Experiment2_6/results/basis_23.json")
print(f"[setup] basis:  {_basis_path}")

with open(_basis_path) as _fh:
    _basis_data = json.load(_fh)

Q_P = [tuple(int(x) for x in e) for e in _basis_data["Q_P"]]
C_Q = [tuple(int(x) for x in e) for e in _basis_data["C_Q"]]
LAWS_COEFFS = [[int(x) for x in law["coeffs"]] for law in _basis_data["laws"]]

_QR_EXPECTED = [(2, 0, 0), (1, 1, 0), (1, 0, 1),
                (0, 2, 0), (0, 1, 1), (0, 0, 2)]
_CQ_EXPECTED = [(0, 0, 3), (0, 1, 2), (0, 2, 1), (0, 3, 0),
                (1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (3, 0, 0)]
assert Q_P == _QR_EXPECTED, f"Q_P order mismatch: {Q_P}"
assert C_Q == _CQ_EXPECTED, f"C_Q order mismatch: {C_Q}"

assert len(LAWS_COEFFS) == 9, f"expected 9 laws, got {len(LAWS_COEFFS)}"
for k, c in enumerate(LAWS_COEFFS):
    assert len(c) == 162, f"law {k}: {len(c)} coeffs, expected 162"

_report_path = find_file([
    HERE / ".." / "Experiment7" / "results" / "report.json",
    HERE / ".." / ".." / "Experiment7" / "results" / "report.json",
    HERE / "report.json",
], "Experiment7/results/report.json")
print(f"[setup] report: {_report_path}")

with open(_report_path) as _fh:
    _report = json.load(_fh)


# =====================================================================
# Decode each law once
# =====================================================================
def decode_coordinate(coeffs, coord_idx):
    base = coord_idx * 54
    out = []
    for mu in range(6):
        sub = base + mu * 9
        for nu in range(9):
            c = int(coeffs[sub + nu]) % P_MOD
            if c:
                out.append((mu, nu, c))
    return out


def decode_law(coeffs):
    return [decode_coordinate(coeffs, k) for k in range(3)]


DECODED_LAWS = [decode_law(c) for c in LAWS_COEFFS]


# =====================================================================
# Recover the seven complete triples from EXP-7
# =====================================================================
RE_L1R3 = re.compile(r"^L1_r3_s(\d+)$")
_triples_all = list(combinations(range(9), 3))
complete_triples = []
for entry_id, entry in _report["phases"]["L1_r3"]["results_by_id"].items():
    if entry.get("complete") is True:
        m = RE_L1R3.match(entry_id)
        assert m, f"unexpected test id: {entry_id}"
        idx = int(m.group(1))
        if "indices" in entry:
            assert list(entry["indices"]) == list(_triples_all[idx]), \
                f"index mismatch for {entry_id}"
        complete_triples.append({
            "exp7_id": entry_id,
            "combinatorial_index": idx,
            "indices": list(_triples_all[idx]),
        })
complete_triples.sort(key=lambda t: t["combinatorial_index"])
assert len(complete_triples) == 7, \
    f"expected 7 complete triples, got {len(complete_triples)}"

print(f"[setup] complete triples from EXP-7: {len(complete_triples)}")
for t in complete_triples:
    print(f"  {t['exp7_id']}  indices={t['indices']}")


# =====================================================================
# Signed representative / coefficient classification
# =====================================================================
def signed_rep(c, p=P_MOD):
    c %= p
    return c - p if c > p // 2 else c


def classify_const(c):
    cs = signed_rep(c)
    if cs == 0:
        return "zero"
    if abs(cs) == 1:
        return "one"
    if abs(cs) < (1 << SMALL_CONST_BITS):
        return "small"
    return "full"


# =====================================================================
# Operation-count references
# =====================================================================
COMPOSITE_ADD = {"M": 21, "S": 6, "C": 12, "A": 31}
COMPOSITE_DBL = {"M": 24, "S": 9, "C": 21, "A": 34}
JACOBIAN_ADD  = {"M": 11, "S": 5, "C": 0, "A": 13}
JACOBIAN_DBL  = {"M": 1,  "S": 8, "C": 1, "A": 14}


# =====================================================================
# Weierstrass affine oracle (E_W)
# =====================================================================
def w_neg(P):
    if P is None:
        return None
    x, y = P
    return (x, (-y) % P_MOD)


def w_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % P_MOD == 0:
        return None
    if x1 == x2 and y1 == y2:
        if y1 == 0:
            return None
        lam = (3 * x1 * x1 + A_W) % P_MOD
        lam = lam * pow(2 * y1 % P_MOD, -1, P_MOD) % P_MOD
    else:
        lam = (y2 - y1) % P_MOD
        lam = lam * pow((x2 - x1) % P_MOD, -1, P_MOD) % P_MOD
    x3 = (lam * lam - x1 - x2) % P_MOD
    y3 = (lam * (x1 - x3) - y1) % P_MOD
    return (x3, y3)


def w_scalar_mul(k, P):
    R = None
    while k > 0:
        if k & 1:
            R = w_add(R, P)
        P = w_add(P, P)
        k >>= 1
    return R


# ---- pK <-> E_W (phi_0 composed with the shift sigma) ----
_INVW  = pow(W, -1, P_MOD)
_U_N   = U * _INVW % P_MOD
_V_N   = V * _INVW % P_MOD
_A_MAP = (_U_N + 1) * pow((_U_N + _V_N) % P_MOD, -1, P_MOD) % P_MOD
_D_MAP = (_U_N - 1) * pow((_U_N + 1) % P_MOD, -1, P_MOD) % P_MOD
_LAM   = (_V_N * _V_N - _U_N * _U_N) % P_MOD \
         * pow((1 - _U_N * _U_N) % P_MOD, -1, P_MOD) % P_MOD
_SHIFT = _D_MAP * (1 + _LAM) % P_MOD * pow(3, -1, P_MOD) % P_MOD


def phi0_W(x_pk, y_pk):
    if (x_pk + 1) % P_MOD == 0 or (y_pk + 1) % P_MOD == 0:
        return None
    r = (x_pk - 1) * pow((x_pk + 1) % P_MOD, -1, P_MOD) % P_MOD
    s = (y_pk - 1) * pow((y_pk + 1) % P_MOD, -1, P_MOD) % P_MOD
    if s == 0:
        return None
    t = r * pow(s, -1, P_MOD) % P_MOD
    T = t * pow(_A_MAP, -1, P_MOD) % P_MOD
    X = _D_MAP * T % P_MOD
    Y = _D_MAP * s % P_MOD * T % P_MOD * ((T - 1) % P_MOD) % P_MOD
    XW = (X - _SHIFT) % P_MOD
    return (XW, Y)


def phi0_W_inv(XW, YW):
    X = (XW + _SHIFT) % P_MOD
    if X == 0:
        return None
    T = X * pow(_D_MAP, -1, P_MOD) % P_MOD
    denom = _D_MAP * T % P_MOD * ((T - 1) % P_MOD) % P_MOD
    if denom == 0:
        return None
    s = YW * pow(denom, -1, P_MOD) % P_MOD
    r = _A_MAP * T % P_MOD * s % P_MOD
    if (1 - r) % P_MOD == 0 or (1 - s) % P_MOD == 0:
        return None
    x = (1 + r) * pow((1 - r) % P_MOD, -1, P_MOD) % P_MOD
    y = (1 + s) * pow((1 - s) % P_MOD, -1, P_MOD) % P_MOD
    return (x, y)


G_PK = (
    1,
    43885198696659819331868805857976446729580581201517994933321256281291666494433,
    52096118237682591038514212999531378342298973921021687132305677862559506645034,
)


def G_pk_to_W():
    alpha, beta, gamma = G_PK
    inv_a = pow(alpha, -1, P_MOD)
    x = beta * inv_a % P_MOD
    y = gamma * inv_a % P_MOD
    return phi0_W(x, y)


def _W_D_on_E_W():
    xD = V * pow(U, -1, P_MOD) % P_MOD
    yD = W * pow(U, -1, P_MOD) % P_MOD
    return phi0_W(xD, yD)


W_D_W = _W_D_on_E_W()


# =====================================================================
# On-curve checks and pK helpers
# =====================================================================
def on_curve_pk(P):
    a, b, c = P
    return (U * a * (b * b - c * c)
            + V * b * (c * c - a * a)
            + W * c * (a * a - b * b)) % P_MOD == 0


def pk_to_affine(P):
    a, b, c = P
    if a == 0:
        return None
    ia = pow(a, -1, P_MOD)
    return (b * ia % P_MOD, c * ia % P_MOD)


def composite_pk_add(P, Q):
    a, b, c = P
    x, y, z = Q
    Fa = (U * (b * b - c * c) + 2 * a * (W * c - V * b)) % P_MOD
    Fb = (V * (c * c - a * a) + 2 * b * (U * a - W * c)) % P_MOD
    Fc = (W * (a * a - b * b) + 2 * c * (V * b - U * a)) % P_MOD
    Ga = (U * (y * y - z * z) + 2 * x * (W * z - V * y)) % P_MOD
    Gb = (V * (z * z - x * x) + 2 * y * (U * x - W * z)) % P_MOD
    Gc = (W * (x * x - y * y) + 2 * z * (V * y - U * x)) % P_MOD
    A_ = (a * Ga + b * Gb + c * Gc) % P_MOD
    B_ = (x * Fa + y * Fb + z * Fc) % P_MOD
    R0 = (A_ * a - B_ * x) % P_MOD
    R1 = (A_ * b - B_ * y) % P_MOD
    R2 = (A_ * c - B_ * z) % P_MOD
    return ((R1 * R2) % P_MOD, (R2 * R0) % P_MOD, (R0 * R1) % P_MOD)


def composite_pk_double(P):
    a, b, c = P
    Fa = (U * (b * b - c * c) + 2 * a * (W * c - V * b)) % P_MOD
    Fb = (V * (c * c - a * a) + 2 * b * (U * a - W * c)) % P_MOD
    Fc = (W * (a * a - b * b) + 2 * c * (V * b - U * a)) % P_MOD
    B_ = Fb
    G_ = Fc
    Ft = (V * G_ * B_ * B_ + W * B_ * G_ * G_) % P_MOD
    H  = ((U * a - W * c) * G_ * G_
          - 2 * (V * c - W * b) * B_ * G_
          + (V * b - U * a) * B_ * B_) % P_MOD
    R0 = (-Ft * a) % P_MOD
    R1 = (-Ft * b + H * G_) % P_MOD
    R2 = (-Ft * c - H * B_) % P_MOD
    return ((R1 * R2) % P_MOD, (R2 * R0) % P_MOD, (R0 * R1) % P_MOD)


def projectively_equal(P, Q):
    x1, y1, z1 = P
    x2, y2, z2 = Q
    return ((y1 * z2 - z1 * y2) % P_MOD == 0 and
            (z1 * x2 - x1 * z2) % P_MOD == 0 and
            (x1 * y2 - y1 * x2) % P_MOD == 0)


def _is_nonzero_triple(T):
    return any(c % P_MOD for c in T)


def oracle_sum_pk(P, Q):
    """Independent oracle: phi_0(S) = phi_0(P) + phi_0(Q) - W_D on E_W.

    This equals phi_D(P) + phi_D(Q) + W_D in the D-centered model.
    """
    aff_P = pk_to_affine(P)
    aff_Q = pk_to_affine(Q)
    if aff_P is None or aff_Q is None:
        return None
    PW = phi0_W(*aff_P)
    QW = phi0_W(*aff_Q)
    if PW is None or QW is None or W_D_W is None:
        return None
    SW = w_add(w_add(PW, QW), w_neg(W_D_W))
    if SW is None:
        return None
    back = phi0_W_inv(*SW)
    if back is None:
        return None
    xS, yS = back
    P_oracle = (1, xS, yS)
    assert on_curve_pk(P_oracle), \
        f"oracle output not on curve: {P_oracle}"
    return P_oracle


# --- Point sampling for verify ---
def random_curve_point(rng):
    """Uniformly-flavoured sample of E(F_p) with the a = 1 chart.

    Curve: (b-5)c^2 + 7(1-b^2)c + (5b^2-b) = 0.
    """
    for _ in range(50):
        b = rng.randrange(P_MOD)
        A2 = (b - 5) % P_MOD
        B1 = 7 * (1 - b * b) % P_MOD
        C0 = (5 * b * b - b) % P_MOD
        if A2 == 0:
            continue
        disc = (B1 * B1 - 4 * A2 * C0) % P_MOD
        s = pow(disc, (P_MOD + 1) // 4, P_MOD)
        if s * s % P_MOD != disc:
            continue
        c = (-B1 + (s if rng.randrange(2) else (-s) % P_MOD)) \
            * pow(2 * A2, -1, P_MOD) % P_MOD
        P = (1, b, c)
        if on_curve_pk(P):
            return P
    return None


def _collect_on_curve_points(rng, target=40, max_attempts=2000):
    """Sample on-curve pK points.

    Note: the first half is derived from multiples of the pK generator G
    on E_W, i.e. phi_0(G_W) = [k] phi_0(G), so those points lie in the
    image of the order-q subgroup of E_W.  The second half is drawn from
    the full E(F_p) via the (a=1) quadratic chart.  Both sets are valid
    on-curve pK points; only the first is confined to the order-q
    subgroup.  (The reviewer's earlier "order-q subgroup sampling"
    comment refers to a D-centered subgroup statement; corrected here.)
    """
    G_W = G_pk_to_W()
    assert G_W is not None, "generator image unavailable"
    pts = []
    attempts = 0
    while len(pts) < target // 2 and attempts < max_attempts:
        attempts += 1
        k = rng.randrange(1, Q_ORDER)
        PW = w_scalar_mul(k, G_W)
        if PW is None:
            continue
        aff = phi0_W_inv(*PW)
        if aff is None:
            continue
        x, y = aff
        if x == 0 or y == 0:
            continue
        P_pk = (1, x, y)
        if on_curve_pk(P_pk):
            pts.append(P_pk)
    while len(pts) < target and attempts < max_attempts:
        attempts += 1
        P = random_curve_point(rng)
        if P is not None:
            pts.append(P)
    return pts


def inverse_pair(P):
    """Q such that P + Q = D in the D-centered group.

    D-centered group law: phi_D(P) = phi_0(P) - W_D.
    P + Q = D  <=>  (phi_0(P) - W_D) + (phi_0(Q) - W_D) = 0
                <=>  phi_0(Q) = 2*W_D - phi_0(P).

    FIX 1: the previous version computed phi_0(Q) = W_D - phi_0(P),
    which gives the 2-torsion point T_1 = (1:-1:1), not D.
    """
    aff = pk_to_affine(P)
    if aff is None:
        return None
    PW = phi0_W(*aff)
    if PW is None:
        return None
    QW = w_add(w_add(W_D_W, W_D_W), w_neg(PW))    # 2*W_D - phi_0(P)
    if QW is None:
        return None
    back = phi0_W_inv(*QW)
    if back is None:
        return None
    xQ, yQ = back
    Q_pk = (1, xQ, yQ)
    if on_curve_pk(Q_pk):
        return Q_pk
    return None


# =====================================================================
# Codegen: addition
# =====================================================================
def _collect_constants(triple_indices):
    cs = set()
    for li in triple_indices:
        for k in range(3):
            for _, _, c in DECODED_LAWS[li][k]:
                cs.add(signed_rep(c))
    cs.discard(0)
    cs.discard(1)
    cs.discard(-1)
    return cs


def codegen_triple(triple_indices):
    """Straight-line source + the set of literal constants used.

    Coefficient literals are emitted as signed representatives.
    Caller binds each constant name to an int (or mpz) in the exec
    namespace.
    """
    consts = _collect_constants(triple_indices)
    cnames = {}
    for i, c in enumerate(sorted(consts)):
        cnames[c] = f"C_{i}"

    L = ["def _f(P, Q, m=m):",
         "    a, b, c = P",
         "    x, y, z = Q",
         "    a2 = a*a % m; b2 = b*b % m; c2 = c*c % m",
         "    ab = a*b % m; ac = a*c % m; bc = b*c % m",
         "    Pm0, Pm1, Pm2, Pm3, Pm4, Pm5 = a2, ab, ac, b2, bc, c2",
         "    x2 = x*x % m; y2 = y*y % m; z2 = z*z % m; xy = x*y % m",
         "    Qm0 = z2*z % m; Qm1 = y*z2 % m; Qm2 = y2*z % m; Qm3 = y2*y % m",
         "    Qm4 = x*z2 % m; Qm5 = xy*z % m; Qm6 = x*y2 % m; Qm7 = x2*z % m;"
         " Qm8 = x2*x % m",
         "    Pm = (Pm0, Pm1, Pm2, Pm3, Pm4, Pm5)",
         "    Qm = (Qm0, Qm1, Qm2, Qm3, Qm4, Qm5, Qm6, Qm7, Qm8)"]

    all_pairs = set()
    for li in triple_indices:
        for k in range(3):
            for mu, nu, _ in DECODED_LAWS[li][k]:
                all_pairs.add((mu, nu))
    for (mu, nu) in sorted(all_pairs):
        L.append(f"    t_{mu}_{nu} = Pm[{mu}] * Qm[{nu}] % m")

    for li in triple_indices:
        for k in range(3):
            entries = DECODED_LAWS[li][k]
            terms = []
            for mu, nu, c in entries:
                cs = signed_rep(c)
                tn = f"t_{mu}_{nu}"
                if cs == 1:
                    terms.append(tn)
                elif cs == -1:
                    terms.append(f"(-{tn})")
                else:
                    terms.append(f"({cnames[cs]}*{tn})")
            if not terms:
                L.append(f"    Z{li}_{k} = 0")
            else:
                L.append(f"    Z{li}_{k} = (" + " + ".join(terms) + ") % m")

    L.append("    return (")
    for li in triple_indices:
        L.append(f"        (Z{li}_0, Z{li}_1, Z{li}_2),")
    L.append("    )")
    return "\n".join(L), {c: cnames[c] for c in cnames}


def compile_triple_exec(triple_indices, backend="int"):
    src, cnames = codegen_triple(triple_indices)
    if backend == "int":
        ns = {"m": P_MOD}
        for c, name in cnames.items():
            ns[name] = c
    else:
        import gmpy2
        ns = {"m": gmpy2.mpz(P_MOD)}
        for c, name in cnames.items():
            ns[name] = gmpy2.mpz(c)
    exec(src, ns)
    return ns["_f"], src


# =====================================================================
# Codegen: diagonal Q = P
# =====================================================================
def codegen_diag(triple_indices):
    per_coord = []
    for li in triple_indices:
        for k in range(3):
            d = {}
            for mu, nu, c in DECODED_LAWS[li][k]:
                e = tuple(Q_P[mu][i] + C_Q[nu][i] for i in range(3))
                d[e] = (d.get(e, 0) + c) % P_MOD
            per_coord.append({e: v for e, v in d.items() if v})

    exps = sorted({e for d in per_coord for e in d})
    max_exp = {"a": 0, "b": 0, "c": 0}
    for e in exps:
        for i, v in enumerate("abc"):
            if e[i] > max_exp[v]:
                max_exp[v] = e[i]

    consts = set()
    for d in per_coord:
        for c in d.values():
            cs = signed_rep(c)
            if cs not in (0, 1, -1):
                consts.add(cs)
    cnames = {c: f"D_{i}" for i, c in enumerate(sorted(consts))}

    L = ["def _g(P, m=m):", "    a, b, c = P"]
    nS = nM = 0
    for v in "abc":
        me = max_exp[v]
        if me >= 2:
            L.append(f"    {v}2={v}*{v}%m"); nS += 1
        if me >= 3:
            L.append(f"    {v}3={v}2*{v}%m"); nM += 1
        if me >= 4:
            L.append(f"    {v}4={v}2*{v}2%m"); nS += 1
        if me >= 5:
            L.append(f"    {v}5={v}4*{v}%m"); nM += 1

    seen = set()
    mono = {}
    for e in exps:
        fs = [f"{v}{i}" if i > 1 else v for v, i in zip("abc", e) if i]
        if not fs:
            continue
        cur = fs[0]
        for f in fs[1:]:
            key = f"{cur}_{f}"
            if key not in seen:
                L.append(f"    {key}={cur}*{f}%m")
                seen.add(key)
                nM += 1
            cur = key
        mono[e] = cur

    for i, d in enumerate(per_coord):
        terms = []
        for e, c in d.items():
            if e not in mono:
                continue
            cs = signed_rep(c)
            if cs == 1:
                terms.append(mono[e])
            elif cs == -1:
                terms.append(f"(-{mono[e]})")
            else:
                terms.append(f"({cnames[cs]}*{mono[e]})")
        if not terms:
            L.append(f"    Z{i} = 0")
        else:
            L.append(f"    Z{i} = (" + " + ".join(terms) + ") % m")

    L.append("    return ((Z0, Z1, Z2), (Z3, Z4, Z5), (Z6, Z7, Z8))")

    C_one = C_small = C_full = A = 0
    for d in per_coord:
        if not d:
            continue
        A += len(d) - 1
        for c in d.values():
            cls = classify_const(c)
            if cls == "one":
                C_one += 1
            elif cls == "small":
                C_small += 1
            elif cls == "full":
                C_full += 1

    return "\n".join(L), cnames, {
        "S": nS, "M_general": nM,
        "C_one": C_one, "C_small": C_small, "C_full": C_full,
        "A": A, "distinct_monomials": len(exps),
    }


def compile_diag_exec(triple_indices, backend="int"):
    src, cnames, counts = codegen_diag(triple_indices)
    if backend == "int":
        ns = {"m": P_MOD}
        for c, name in cnames.items():
            ns[name] = c
    else:
        import gmpy2
        ns = {"m": gmpy2.mpz(P_MOD)}
        for c, name in cnames.items():
            ns[name] = gmpy2.mpz(c)
    exec(src, ns)
    return ns["_g"], counts


# =====================================================================
# Static cost models
# =====================================================================
def p_monomial_cost():
    return {"S": 3, "M_general": 3}


def q_monomial_cost():
    return {"S": 3, "M_general": 10}


def static_order_A(triple_indices):
    coord_entries = []
    for li in triple_indices:
        for k in range(3):
            coord_entries.append(DECODED_LAWS[li][k])

    union_products = set()
    for entries in coord_entries:
        for mu, nu, _ in entries:
            union_products.add((mu, nu))

    C_one = C_small = C_full = A = 0
    for entries in coord_entries:
        if not entries:
            continue
        A += len(entries) - 1
        for _, _, c in entries:
            cls = classify_const(c)
            if cls == "one":
                C_one += 1
            elif cls == "small":
                C_small += 1
            elif cls == "full":
                C_full += 1

    pm = p_monomial_cost()
    qm = q_monomial_cost()
    return {
        "order": "A",
        "M_general": len(union_products) + pm["M_general"] + qm["M_general"],
        "S": pm["S"] + qm["S"],
        "C_one": C_one, "C_small": C_small, "C_full": C_full,
        "A": A,
        "distinct_products": len(union_products),
    }


def static_order_B(triple_indices):
    coord_entries = []
    for li in triple_indices:
        for k in range(3):
            coord_entries.append(DECODED_LAWS[li][k])

    M_general = 0
    C_one = C_small = C_full = A = 0
    for entries in coord_entries:
        if not entries:
            continue
        by_mu = {}
        for mu, nu, c in entries:
            by_mu.setdefault(mu, []).append((nu, c))
        n_mu_groups = len(by_mu)

        for _, nu_list in by_mu.items():
            A += len(nu_list) - 1
            for _, c in nu_list:
                cls = classify_const(c)
                if cls == "one":
                    C_one += 1
                elif cls == "small":
                    C_small += 1
                elif cls == "full":
                    C_full += 1
        M_general += n_mu_groups
        A += n_mu_groups - 1

    pm = p_monomial_cost()
    qm = q_monomial_cost()
    M_general += pm["M_general"] + qm["M_general"]
    S = pm["S"] + qm["S"]

    return {
        "order": "B", "M_general": M_general, "S": S,
        "C_one": C_one, "C_small": C_small, "C_full": C_full, "A": A,
    }


def heuristic(costs):
    w = HEUR_WEIGHTS
    return (w["M_general"] * costs["M_general"]
            + w["C_full"]  * costs["C_full"]
            + w["C_small"] * costs["C_small"]
            + w["S"]       * costs["S"]
            + w["A"]       * costs["A"])


# =====================================================================
# PHASE  profile
# =====================================================================
def phase_profile():
    print("=" * 72)
    print("EXP-8  —  coefficient profile")
    print("=" * 72)

    rows = []
    for li in range(9):
        for k in range(3):
            entries = DECODED_LAWS[li][k]
            signed = [signed_rep(c) for _, _, c in entries]
            nz = len(entries)
            n_pos1 = sum(1 for c in signed if c == 1)
            n_neg1 = sum(1 for c in signed if c == -1)
            n_small = sum(1 for c in signed
                          if 1 < abs(c) < (1 << SMALL_CONST_BITS))
            n_full = sum(1 for c in signed
                         if abs(c) >= (1 << SMALL_CONST_BITS))
            bitlens = sorted(abs(c).bit_length() for c in signed)
            rows.append({
                "law": li, "coord": k, "nz": nz,
                "n_pos1": n_pos1, "n_neg1": n_neg1,
                "n_small": n_small, "n_full": n_full,
                "bitlen_min": bitlens[0] if bitlens else 0,
                "bitlen_med": bitlens[len(bitlens) // 2] if bitlens else 0,
                "bitlen_max": bitlens[-1] if bitlens else 0,
            })

    print(f"{'law':<4}{'coord':<6}{'nz':<5}{'+1':<4}{'-1':<4}"
          f"{'small':<7}{'full':<6}{'bits(min/med/max)':<20}")
    print("-" * 72)
    for r in rows:
        print(f"{r['law']:<4}{r['coord']:<6}{r['nz']:<5}"
              f"{r['n_pos1']:<4}{r['n_neg1']:<4}"
              f"{r['n_small']:<7}{r['n_full']:<6}"
              f"{r['bitlen_min']:>3}/{r['bitlen_med']:>3}/{r['bitlen_max']:<4}")

    tot_pos1  = sum(r["n_pos1"] for r in rows)
    tot_neg1  = sum(r["n_neg1"] for r in rows)
    tot_small = sum(r["n_small"] for r in rows)
    tot_full  = sum(r["n_full"] for r in rows)
    print(f"\nTotals:  +1 {tot_pos1}   -1 {tot_neg1}   "
          f"small {tot_small}   full {tot_full}")
    if tot_full > 0:
        print("Note: this basis is dense in the coefficient-cost sense; "
              "an LLL-reduced basis over Q is the natural next step.")

    atomic_write_json(PROFILE_JSON, {
        "p_mod_bits": P_MOD.bit_length(),
        "small_const_bits": SMALL_CONST_BITS,
        "rows": rows,
        "totals": {"n_pos1": tot_pos1, "n_neg1": tot_neg1,
                   "n_small": tot_small, "n_full": tot_full},
    })
    print(f"[profile] wrote {PROFILE_JSON}")
    return 0


# =====================================================================
# PHASE  analyze
# =====================================================================
def phase_analyze():
    print("=" * 72)
    print("EXP-8.1  —  static CSE-aware cost analysis")
    print("=" * 72)

    print(f"\nComposite pK baseline (Table 3):")
    print(f"  addition : {COMPOSITE_ADD['M']}M + {COMPOSITE_ADD['S']}S "
          f"+ {COMPOSITE_ADD['C']}C + {COMPOSITE_ADD['A']}A")
    print(f"  doubling : {COMPOSITE_DBL['M']}M + {COMPOSITE_DBL['S']}S "
          f"+ {COMPOSITE_DBL['C']}C + {COMPOSITE_DBL['A']}A")
    print(f"Jacobian reference (EXP-3):")
    print(f"  addition : {JACOBIAN_ADD['M']}M + {JACOBIAN_ADD['S']}S "
          f"+ {JACOBIAN_ADD['A']}A")
    print(f"  doubling : {JACOBIAN_DBL['M']}M + {JACOBIAN_DBL['S']}S "
          f"+ {JACOBIAN_DBL['C']}C + {JACOBIAN_DBL['A']}A")

    out = {"triples": [], "ranking_A": [], "ranking_B": [], "ranking_diag": []}

    print(f"\nStatic cost — order A (product-first, shared Pm*Qm):")
    header = (f"{'#':<3}{'exp7_id':<16}{'indices':<14}"
              f"{'Mgen':<7}{'S':<5}{'C+/-1':<7}{'Csm':<6}{'Cfull':<7}"
              f"{'A':<6}{'heur':<8}")
    print(header); print("-" * 72)
    for i, t in enumerate(complete_triples):
        cA = static_order_A(t["indices"])
        cA["heuristic"] = heuristic(cA)
        cA["exp7_id"] = t["exp7_id"]
        cA["indices"] = t["indices"]
        out["triples"].append({"order_A": cA})
        print(f"{i:<3}{t['exp7_id']:<16}{str(t['indices']):<14}"
              f"{cA['M_general']:<7}{cA['S']:<5}{cA['C_one']:<7}"
              f"{cA['C_small']:<6}{cA['C_full']:<7}{cA['A']:<6}"
              f"{cA['heuristic']:<8.1f}")
    ranked_A = sorted(out["triples"], key=lambda x: x["order_A"]["heuristic"])
    out["ranking_A"] = [x["order_A"]["exp7_id"] for x in ranked_A]

    print(f"\nStatic cost — order B (coefficient-first, no shared products):")
    print(header); print("-" * 72)
    for i, t in enumerate(complete_triples):
        cB = static_order_B(t["indices"])
        cB["heuristic"] = heuristic(cB)
        cB["exp7_id"] = t["exp7_id"]
        cB["indices"] = t["indices"]
        out["triples"][i]["order_B"] = cB
        print(f"{i:<3}{t['exp7_id']:<16}{str(t['indices']):<14}"
              f"{cB['M_general']:<7}{cB['S']:<5}{cB['C_one']:<7}"
              f"{cB['C_small']:<6}{cB['C_full']:<7}{cB['A']:<6}"
              f"{cB['heuristic']:<8.1f}")
    ranked_B = sorted(out["triples"], key=lambda x: x["order_B"]["heuristic"])
    out["ranking_B"] = [x["order_B"]["exp7_id"] for x in ranked_B]

    print(f"\nStatic cost — diagonal Q = P (degree-5 collapse, codegen-derived):")
    print(f"{'#':<3}{'exp7_id':<16}{'indices':<14}"
          f"{'Mgen':<7}{'S':<5}{'Cfull':<7}{'A':<6}{'#mon':<7}{'heur':<8}")
    print("-" * 72)
    for i, t in enumerate(complete_triples):
        _, counts = compile_diag_exec(t["indices"])
        counts["heuristic"] = heuristic(counts)
        counts["exp7_id"] = t["exp7_id"]
        counts["indices"] = t["indices"]
        out["triples"][i]["diag"] = counts
        print(f"{i:<3}{t['exp7_id']:<16}{str(t['indices']):<14}"
              f"{counts['M_general']:<7}{counts['S']:<5}"
              f"{counts['C_full']:<7}{counts['A']:<6}"
              f"{counts['distinct_monomials']:<7}"
              f"{counts['heuristic']:<8.1f}")
    ranked_D = sorted(out["triples"], key=lambda x: x["diag"]["heuristic"])
    out["ranking_diag"] = [x["diag"]["exp7_id"] for x in ranked_D]

    out["weights"] = HEUR_WEIGHTS
    out["composite_add"] = COMPOSITE_ADD
    out["composite_dbl"] = COMPOSITE_DBL
    out["jacobian_add"]  = JACOBIAN_ADD
    out["jacobian_dbl"]  = JACOBIAN_DBL

    atomic_write_json(ANALYSIS_JSON, out)
    print(f"\n[analyze] wrote {ANALYSIS_JSON}")
    return 0


# =====================================================================
# PHASE  verify
# =====================================================================
def check_all_laws(outs, want):
    n = 0
    for T in outs:
        if _is_nonzero_triple(T):
            if not projectively_equal(T, want):
                return False, n
            n += 1
    return (n > 0), n


def phase_verify():
    print("=" * 72)
    print("EXP-8  —  correctness verification")
    print("=" * 72)

    assert on_curve_pk(G_PK), "generator does not satisfy the pK equation"
    print("[verify] generator satisfies F(G) = 0")

    G_W = G_pk_to_W()
    assert G_W is not None, "could not map generator to E_W"
    xW, yW = G_W
    assert (yW * yW - (xW * xW * xW + A_W * xW + B_W)) % P_MOD == 0, \
        "generator image fails the Weierstrass equation"
    print(f"[verify] generator image on E_W: x = {xW}")
    print(f"[verify] W_D on E_W: {W_D_W}")
    assert W_D_W is not None, "W_D is None"
    xD, yD = W_D_W
    assert (yD * yD - (xD * xD * xD + A_W * xD + B_W)) % P_MOD == 0, \
        "W_D is not on E_W"

    rng = random.Random(0)
    pts = _collect_on_curve_points(rng, target=40)
    assert len(pts) >= 20, \
        f"only {len(pts)} on-curve points generated; aborting"
    print(f"[verify] generated {len(pts)} on-curve pK points")

    # FIX 3: pre-check the E_W oracle against the composite pK formula.
    # EXP-1 already established that the composite formula agrees with
    # its own Phi_D oracle; a failure here means the E_W map in this
    # file is wrong, not the triples.
    agree_add = 0
    attempts = 0
    while agree_add < 30 and attempts < 300:
        attempts += 1
        P = pts[rng.randrange(len(pts))]
        Q = pts[rng.randrange(len(pts))]
        if P == Q:
            continue
        want = oracle_sum_pk(P, Q)
        comp = composite_pk_add(P, Q)
        if want is not None and _is_nonzero_triple(comp):
            assert projectively_equal(want, comp), \
                f"oracle disagrees with composite add on P={P}, Q={Q}"
            agree_add += 1
    assert agree_add >= 30, \
        f"only {agree_add} oracle/composite addition comparisons"

    agree_dbl = 0
    attempts = 0
    while agree_dbl < 15 and attempts < 300:
        attempts += 1
        P = pts[rng.randrange(len(pts))]
        want = oracle_sum_pk(P, P)
        comp = composite_pk_double(P)
        if want is not None and _is_nonzero_triple(comp):
            assert projectively_equal(want, comp), \
                f"oracle disagrees with composite double on P={P}"
            agree_dbl += 1
    assert agree_dbl >= 15, \
        f"only {agree_dbl} oracle/composite doubling comparisons"
    print(f"[verify] oracle/composite agreement: "
          f"add {agree_add}, dbl {agree_dbl}")

    # Inverse pairs: P + Q = D
    inv_pairs = []
    for _ in range(40):
        P = pts[rng.randrange(len(pts))]
        Q = inverse_pair(P)
        if Q is not None:
            inv_pairs.append((P, Q))
    assert len(inv_pairs) >= MIN_INV_OK, \
        f"only {len(inv_pairs)} inverse pairs constructed"
    print(f"[verify] generated {len(inv_pairs)} inverse pairs (P+Q=D)")

    failures = 0
    for t in complete_triples:
        exp7_id = t["exp7_id"]
        indices = t["indices"]
        add_fn, _ = compile_triple_exec(indices)
        dbl_fn, _ = compile_diag_exec(indices)

        add_ok = add_n = 0
        dbl_ok = dbl_n = 0
        inv_ok = 0

        # --- Addition on random pairs ---
        for _ in range(80):
            if add_ok >= MIN_ADD_OK:
                break
            P = pts[rng.randrange(len(pts))]
            Q = pts[rng.randrange(len(pts))]
            if P == Q:
                continue
            try:
                outs = add_fn(P, Q)
            except Exception as exc:
                print(f"  [verify] {exp7_id} add: exception {exc}")
                failures += 1
                break
            oracle = oracle_sum_pk(P, Q)
            if oracle is None:
                continue
            ok, n = check_all_laws(outs, oracle)
            if not ok:
                print(f"  [verify] {exp7_id}: add MISMATCH vs oracle")
                print(f"    P = {P}")
                print(f"    Q = {Q}")
                print(f"    oracle = {oracle}")
                for j, T in enumerate(outs):
                    print(f"    law{j} = {T}")
                failures += 1
                break
            add_n += n
            add_ok += 1

        # --- Addition on inverse pairs --- FIX 2: own counter
        D_pk = (U, V, W)
        for (P, Q) in inv_pairs:
            if inv_ok >= MIN_INV_OK:
                break
            try:
                outs = add_fn(P, Q)
            except Exception as exc:
                print(f"  [verify] {exp7_id} inv-add: exception {exc}")
                failures += 1
                break
            ok, n = check_all_laws(outs, D_pk)
            if not ok:
                print(f"  [verify] {exp7_id}: inverse-pair MISMATCH")
                print(f"    P = {P}")
                print(f"    Q = {Q}")
                for j, T in enumerate(outs):
                    print(f"    law{j} = {T}")
                failures += 1
                break
            inv_ok += 1

        # --- Doubling ---
        for _ in range(60):
            if dbl_ok >= MIN_DBL_OK:
                break
            P = pts[rng.randrange(len(pts))]
            try:
                outs = dbl_fn(P)
            except Exception as exc:
                print(f"  [verify] {exp7_id} dbl: exception {exc}")
                failures += 1
                break
            oracle = oracle_sum_pk(P, P)
            if oracle is None:
                continue
            ok, n = check_all_laws(outs, oracle)
            if not ok:
                print(f"  [verify] {exp7_id}: dbl MISMATCH vs oracle")
                print(f"    P = {P}")
                print(f"    oracle = {oracle}")
                for j, T in enumerate(outs):
                    print(f"    law{j} = {T}")
                failures += 1
                break
            dbl_n += n
            dbl_ok += 1

        status = "OK"
        if (add_ok < MIN_ADD_OK or dbl_ok < MIN_DBL_OK
                or inv_ok < MIN_INV_OK):
            status = "FAIL"
            failures += 1
        print(f"[verify] {exp7_id}  "
              f"add={add_ok}/{MIN_ADD_OK} (laws {add_n})  "
              f"dbl={dbl_ok}/{MIN_DBL_OK} (laws {dbl_n})  "
              f"inv={inv_ok}/{MIN_INV_OK}  {status}")

    atomic_write_json(VERIFY_JSON, {
        "n_on_curve_points": len(pts),
        "n_inverse_pairs": len(inv_pairs),
        "oracle_composite_agreement": {
            "add": agree_add, "dbl": agree_dbl,
        },
        "triples_verified": [t["exp7_id"] for t in complete_triples],
        "failures": failures,
    })
    print(f"\n[verify] wrote {VERIFY_JSON}")
    return 0 if failures == 0 else 1


# =====================================================================
# ct_select_nonzero
# =====================================================================
def ct_select_nonzero(results):
    X = Y = Z = 0
    done = 0
    for (Xi, Yi, Zi) in results:
        nz = 0 if (Xi == 0 and Yi == 0 and Zi == 0) else 1
        take = nz if not done else 0
        X = (X + take * Xi) % P_MOD
        Y = (Y + take * Yi) % P_MOD
        Z = (Z + take * Zi) % P_MOD
        done = done or take
    return (X, Y, Z)


# =====================================================================
# Jacobian baselines
# =====================================================================
def jac_add(P, Q):
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q
    m = P_MOD
    Z1Z1 = Z1 * Z1 % m
    Z2Z2 = Z2 * Z2 % m
    U1 = X1 * Z2Z2 % m
    U2 = X2 * Z1Z1 % m
    S1 = Y1 * Z2 % m * Z2Z2 % m
    S2 = Y2 * Z1 % m * Z1Z1 % m
    H = (U2 - U1) % m
    I = 4 * H * H % m
    J = H * I % m
    r = 2 * (S2 - S1) % m
    V = U1 * I % m
    X3 = (r * r - J - 2 * V) % m
    Y3 = (r * (V - X3) - 2 * S1 * J) % m
    Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H % m
    return X3, Y3, Z3


def jac_dbl(P):
    X1, Y1, Z1 = P
    m = P_MOD
    XX = X1 * X1 % m
    YY = Y1 * Y1 % m
    YYYY = YY * YY % m
    ZZ = Z1 * Z1 % m
    S = 2 * ((X1 + YY) ** 2 - XX - YYYY) % m
    M = (3 * XX + A_W * ZZ * ZZ) % m
    T = (M * M - 2 * S) % m
    return T, (M * (S - T) - 8 * YYYY) % m, ((Y1 + Z1) ** 2 - YY - ZZ) % m


# =====================================================================
# Timing
# =====================================================================
def time_op(fn, inputs, warmup, passes):
    for inp in inputs[:warmup]:
        fn(*inp)
    timed = inputs[warmup:]
    samples = []
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for _ in range(passes):
            t0 = time.perf_counter()
            for inp in timed:
                fn(*inp)
            t1 = time.perf_counter()
            samples.append((t1 - t0) / len(timed))
    finally:
        if was_enabled:
            gc.enable()
    return median(samples)


# =====================================================================
# PHASE  bench
# =====================================================================
def _check_verify_passed():
    """FIX 4: refuse to benchmark without a passing verify."""
    if not VERIFY_JSON.exists():
        print(f"[bench] {VERIFY_JSON} not found; "
              f"run `python3 exp8.py verify` first")
        return False
    with open(VERIFY_JSON) as fh:
        v = json.load(fh)
    if v.get("failures", 1) != 0:
        print(f"[bench] verify reported {v.get('failures')} failures; "
              f"refusing to benchmark")
        return False
    print(f"[bench] verify passed: {VERIFY_JSON}")
    return True


def _bench_backend_int():
    print(f"\n[bench:int] {N_TRIALS} timed inputs x {N_PASSES} passes, "
          f"seeds = {SEEDS}")

    compiled_add = {t["exp7_id"]: compile_triple_exec(t["indices"])[0]
                    for t in complete_triples}
    compiled_dbl = {t["exp7_id"]: compile_diag_exec(t["indices"])[0]
                    for t in complete_triples}

    baseline_add = []
    baseline_dbl = []
    baseline_jadd = []
    baseline_jdbl = []
    for seed in SEEDS:
        rng = random.Random(seed)
        N = N_TRIALS + N_WARMUP
        inputs_add, inputs_dbl = [], []
        for _ in range(N):
            P = (rng.randrange(P_MOD), rng.randrange(P_MOD),
                 rng.randrange(P_MOD))
            Q = (rng.randrange(P_MOD), rng.randrange(P_MOD),
                 rng.randrange(P_MOD))
            inputs_add.append((P, Q))
            inputs_dbl.append((P,))
        baseline_add.append(
            time_op(composite_pk_add,    inputs_add, N_WARMUP, N_PASSES))
        baseline_dbl.append(
            time_op(composite_pk_double, inputs_dbl, N_WARMUP, N_PASSES))
        baseline_jadd.append(
            time_op(jac_add,             inputs_add, N_WARMUP, N_PASSES))
        baseline_jdbl.append(
            time_op(jac_dbl,             inputs_dbl, N_WARMUP, N_PASSES))

    t_base_add  = median(baseline_add)
    t_base_dbl  = median(baseline_dbl)
    t_jac_add   = median(baseline_jadd)
    t_jac_dbl   = median(baseline_jdbl)
    print(f"[bench:int] composite pK  add {t_base_add * 1e6:.2f} us, "
          f"dbl {t_base_dbl * 1e6:.2f} us")
    print(f"[bench:int] Jacobian      add {t_jac_add * 1e6:.2f} us, "
          f"dbl {t_jac_dbl * 1e6:.2f} us")

    results = []
    for t in complete_triples:
        eid = t["exp7_id"]
        add_fn = compiled_add[eid]
        dbl_fn = compiled_dbl[eid]

        per_seed_add = []
        per_seed_dbl = []
        per_seed_sel = []
        for seed in SEEDS:
            rng = random.Random(seed)
            N = N_TRIALS + N_WARMUP
            inputs_add, inputs_dbl = [], []
            for _ in range(N):
                P = (rng.randrange(P_MOD), rng.randrange(P_MOD),
                     rng.randrange(P_MOD))
                Q = (rng.randrange(P_MOD), rng.randrange(P_MOD),
                     rng.randrange(P_MOD))
                inputs_add.append((P, Q))
                inputs_dbl.append((P,))

            def add_full(P, Q, f=add_fn):
                return ct_select_nonzero(f(P, Q))

            def dbl_full(P, f=dbl_fn):
                return ct_select_nonzero(f(P))

            t_add = time_op(add_full, inputs_add, N_WARMUP, N_PASSES)
            t_dbl = time_op(dbl_full, inputs_dbl, N_WARMUP, N_PASSES)

            precomp = [(add_fn(P, Q),) for (P, Q) in inputs_add]
            t_sel = time_op(ct_select_nonzero, precomp, N_WARMUP, N_PASSES)

            per_seed_add.append(t_add)
            per_seed_dbl.append(t_dbl)
            per_seed_sel.append(t_sel)

        results.append({
            "exp7_id": eid,
            "indices": t["indices"],
            "t_add_us_median": median(per_seed_add) * 1e6,
            "t_dbl_us_median": median(per_seed_dbl) * 1e6,
            "t_sel_us_median": median(per_seed_sel) * 1e6,
            "per_seed_add_us": [x * 1e6 for x in per_seed_add],
            "per_seed_dbl_us": [x * 1e6 for x in per_seed_dbl],
        })

    results.sort(key=lambda r: r["t_add_us_median"])
    return results, {
        "composite_add_us": t_base_add * 1e6,
        "composite_dbl_us": t_base_dbl * 1e6,
        "jacobian_add_us":  t_jac_add * 1e6,
        "jacobian_dbl_us":  t_jac_dbl * 1e6,
    }


def _bench_backend_gmp():
    try:
        import gmpy2
    except ImportError:
        return None, None

    m_g = gmpy2.mpz(P_MOD)
    AW_g = gmpy2.mpz(A_W)
    print(f"\n[bench:gmpy2] {N_TRIALS} timed inputs x {N_PASSES} passes, "
          f"seeds = {SEEDS}")

    compiled_add = {}
    compiled_dbl = {}
    for t in complete_triples:
        eid = t["exp7_id"]
        f, _ = compile_triple_exec(t["indices"], backend="gmpy2")
        compiled_add[eid] = f
        g, _ = compile_diag_exec(t["indices"], backend="gmpy2")
        compiled_dbl[eid] = g

    def composite_add_g(P, Q):
        a, b, c = P
        x, y, z = Q
        two = gmpy2.mpz(2)
        Fa = (5 * (b * b - c * c) + two * a * (7 * c - b)) % m_g
        Fb = ((c * c - a * a) + two * b * (5 * a - 7 * c)) % m_g
        Fc = (7 * (a * a - b * b) + two * c * (b - 5 * a)) % m_g
        Ga = (5 * (y * y - z * z) + two * x * (7 * z - y)) % m_g
        Gb = ((z * z - x * x) + two * y * (5 * x - 7 * z)) % m_g
        Gc = (7 * (x * x - y * y) + two * z * (y - 5 * x)) % m_g
        A_ = (a * Ga + b * Gb + c * Gc) % m_g
        B_ = (x * Fa + y * Fb + z * Fc) % m_g
        R0 = (A_ * a - B_ * x) % m_g
        R1 = (A_ * b - B_ * y) % m_g
        R2 = (A_ * c - B_ * z) % m_g
        return ((R1 * R2) % m_g, (R2 * R0) % m_g, (R0 * R1) % m_g)

    def composite_dbl_g(P):
        a, b, c = P
        two = gmpy2.mpz(2)
        Fa = (5 * (b * b - c * c) + two * a * (7 * c - b)) % m_g
        Fb = ((c * c - a * a) + two * b * (5 * a - 7 * c)) % m_g
        Fc = (7 * (a * a - b * b) + two * c * (b - 5 * a)) % m_g
        B_ = Fb
        G_ = Fc
        Ft = (G_ * B_ * B_ + 7 * B_ * G_ * G_) % m_g
        H = ((5 * a - 7 * c) * G_ * G_
             - two * (c - 7 * b) * B_ * G_
             + (b - 5 * a) * B_ * B_) % m_g
        R0 = (-Ft * a) % m_g
        R1 = (-Ft * b + H * G_) % m_g
        R2 = (-Ft * c - H * B_) % m_g
        return ((R1 * R2) % m_g, (R2 * R0) % m_g, (R0 * R1) % m_g)

    def jac_add_g(P, Q):
        X1, Y1, Z1 = P
        X2, Y2, Z2 = Q
        four = gmpy2.mpz(4)
        two = gmpy2.mpz(2)
        Z1Z1 = Z1 * Z1 % m_g
        Z2Z2 = Z2 * Z2 % m_g
        U1 = X1 * Z2Z2 % m_g
        U2 = X2 * Z1Z1 % m_g
        S1 = Y1 * Z2 % m_g * Z2Z2 % m_g
        S2 = Y2 * Z1 % m_g * Z1Z1 % m_g
        H = (U2 - U1) % m_g
        I = four * H * H % m_g
        J = H * I % m_g
        r = two * (S2 - S1) % m_g
        V = U1 * I % m_g
        X3 = (r * r - J - two * V) % m_g
        Y3 = (r * (V - X3) - two * S1 * J) % m_g
        Z3 = ((Z1 + Z2) ** 2 - Z1Z1 - Z2Z2) * H % m_g
        return X3, Y3, Z3

    def jac_dbl_g(P):
        X1, Y1, Z1 = P
        two = gmpy2.mpz(2)
        three = gmpy2.mpz(3)
        eight = gmpy2.mpz(8)
        XX = X1 * X1 % m_g
        YY = Y1 * Y1 % m_g
        YYYY = YY * YY % m_g
        ZZ = Z1 * Z1 % m_g
        S = two * ((X1 + YY) ** 2 - XX - YYYY) % m_g
        M = (three * XX + AW_g * ZZ * ZZ) % m_g
        T = (M * M - two * S) % m_g
        return T, (M * (S - T) - eight * YYYY) % m_g, \
               ((Y1 + Z1) ** 2 - YY - ZZ) % m_g

    baseline_add = []
    baseline_dbl = []
    baseline_jadd = []
    baseline_jdbl = []
    all_inputs_per_seed = []
    for seed in SEEDS:
        rng = random.Random(seed)
        N = N_TRIALS + N_WARMUP
        inputs_add, inputs_dbl = [], []
        for _ in range(N):
            P = tuple(gmpy2.mpz(rng.randrange(P_MOD)) for _ in range(3))
            Q = tuple(gmpy2.mpz(rng.randrange(P_MOD)) for _ in range(3))
            inputs_add.append((P, Q))
            inputs_dbl.append((P,))
        baseline_add.append(
            time_op(composite_add_g, inputs_add, N_WARMUP, N_PASSES))
        baseline_dbl.append(
            time_op(composite_dbl_g, inputs_dbl, N_WARMUP, N_PASSES))
        baseline_jadd.append(
            time_op(jac_add_g, inputs_add, N_WARMUP, N_PASSES))
        baseline_jdbl.append(
            time_op(jac_dbl_g, inputs_dbl, N_WARMUP, N_PASSES))
        all_inputs_per_seed.append((inputs_add, inputs_dbl))

    t_base_add = median(baseline_add)
    t_base_dbl = median(baseline_dbl)
    t_jac_add  = median(baseline_jadd)
    t_jac_dbl  = median(baseline_jdbl)
    print(f"[bench:gmpy2] composite pK  add {t_base_add * 1e6:.2f} us, "
          f"dbl {t_base_dbl * 1e6:.2f} us")
    print(f"[bench:gmpy2] Jacobian      add {t_jac_add * 1e6:.2f} us, "
          f"dbl {t_jac_dbl * 1e6:.2f} us")

    results = []
    for t in complete_triples:
        eid = t["exp7_id"]
        add_fn = compiled_add[eid]
        dbl_fn = compiled_dbl[eid]

        per_seed_add = []
        per_seed_dbl = []
        per_seed_sel = []
        for seed, (inputs_add, inputs_dbl) in zip(SEEDS, all_inputs_per_seed):
            def add_full(P, Q, f=add_fn):
                return ct_select_nonzero(f(P, Q))
            def dbl_full(P, f=dbl_fn):
                return ct_select_nonzero(f(P))
            per_seed_add.append(time_op(add_full, inputs_add, N_WARMUP, N_PASSES))
            per_seed_dbl.append(time_op(dbl_full, inputs_dbl, N_WARMUP, N_PASSES))

            precomp = [(add_fn(P, Q),) for (P, Q) in inputs_add]
            per_seed_sel.append(
                time_op(ct_select_nonzero, precomp, N_WARMUP, N_PASSES))

        results.append({
            "exp7_id": eid,
            "indices": t["indices"],
            "t_add_us_median": median(per_seed_add) * 1e6,
            "t_dbl_us_median": median(per_seed_dbl) * 1e6,
            "t_sel_us_median": median(per_seed_sel) * 1e6,
        })

    results.sort(key=lambda r: r["t_add_us_median"])
    return results, {
        "composite_add_us": t_base_add * 1e6,
        "composite_dbl_us": t_base_dbl * 1e6,
        "jacobian_add_us":  t_jac_add * 1e6,
        "jacobian_dbl_us":  t_jac_dbl * 1e6,
    }


def phase_bench():
    print("=" * 72)
    print("EXP-8.2 / EXP-8.3  —  wall-clock benchmark of complete triples")
    print("=" * 72)

    # FIX 4: hard gate on verify having passed
    if not _check_verify_passed():
        return 1

    # FIX 6: note the timing bias
    print("[bench] note: triple timings include an extra Python function")
    print("[bench]       wrapper and a ct_select_nonzero call, while the")
    print("[bench]       composite/Jacobian baselines are called directly.")
    print("[bench]       Both biases are small and against the (2,3) side.")

    results_int, base_int = _bench_backend_int()

    print(f"\n[int backend] addition ranking (best first):")
    print(f"{'#':<3}{'exp7_id':<16}{'indices':<14}"
          f"{'add(us)':<11}{'dbl(us)':<11}{'sel(us)':<11}"
          f"{'x/comp':<9}{'x/jac':<9}")
    print("-" * 80)
    for i, r in enumerate(results_int, 1):
        print(f"{i:<3}{r['exp7_id']:<16}{str(r['indices']):<14}"
              f"{r['t_add_us_median']:<11.2f}"
              f"{r['t_dbl_us_median']:<11.2f}"
              f"{r['t_sel_us_median']:<11.2f}"
              f"{r['t_add_us_median']/base_int['composite_add_us']:<9.3f}"
              f"{r['t_add_us_median']/base_int['jacobian_add_us']:<9.3f}")

    results_gmp, base_gmp = _bench_backend_gmp()
    if results_gmp is not None:
        print(f"\n[gmpy2 backend] addition ranking (best first):")
        print(f"{'#':<3}{'exp7_id':<16}{'indices':<14}"
              f"{'add(us)':<11}{'dbl(us)':<11}{'sel(us)':<11}"
              f"{'x/comp':<9}{'x/jac':<9}")
        print("-" * 80)
        for i, r in enumerate(results_gmp, 1):
            print(f"{i:<3}{r['exp7_id']:<16}{str(r['indices']):<14}"
                  f"{r['t_add_us_median']:<11.2f}"
                  f"{r['t_dbl_us_median']:<11.2f}"
                  f"{r['t_sel_us_median']:<11.2f}"
                  f"{r['t_add_us_median']/base_gmp['composite_add_us']:<9.3f}"
                  f"{r['t_add_us_median']/base_gmp['jacobian_add_us']:<9.3f}")

    atomic_write_json(BENCH_JSON, {
        "n_trials": N_TRIALS, "n_warmup": N_WARMUP,
        "n_passes": N_PASSES, "seeds": SEEDS,
        "python_version": sys.version.split()[0],
        "gc_disabled_during_timing": True,
        "baseline_int": base_int,
        "baseline_gmp": base_gmp,
        "results_int": results_int,
        "results_gmp": results_gmp,
    })
    print(f"\n[bench] wrote {BENCH_JSON}")
    return 0


# =====================================================================
# PHASE  summary
# =====================================================================
def phase_summary():
    print("=" * 72)
    print("EXP-8  —  summary")
    print("=" * 72)

    need = [ANALYSIS_JSON, BENCH_JSON]
    missing = [p for p in need if not p.exists()]
    if missing:
        print(f"[summary] missing: {missing}; run `analyze` and `bench`")
        return 1

    with open(ANALYSIS_JSON) as fh:
        A = json.load(fh)
    with open(BENCH_JSON) as fh:
        B = json.load(fh)

    print(f"\nStatic ranking (order A):")
    for i, eid in enumerate(A["ranking_A"], 1):
        t = next(x["order_A"] for x in A["triples"]
                 if x["order_A"]["exp7_id"] == eid)
        print(f"  {i}. {eid}  {t['M_general']}M + {t['C_full']}Cfull + "
              f"{t['C_small']}Csm + {t['S']}S + {t['A']}A")

    print(f"\nStatic ranking (order B):")
    for i, eid in enumerate(A["ranking_B"], 1):
        t = next(x["order_B"] for x in A["triples"]
                 if x["order_B"]["exp7_id"] == eid)
        print(f"  {i}. {eid}  {t['M_general']}M + {t['C_full']}Cfull + "
              f"{t['C_small']}Csm + {t['S']}S + {t['A']}A")

    print(f"\nStatic ranking (diagonal):")
    for i, eid in enumerate(A["ranking_diag"], 1):
        t = next(x["diag"] for x in A["triples"]
                 if x["diag"]["exp7_id"] == eid)
        print(f"  {i}. {eid}  {t['M_general']}M + {t['C_full']}Cfull + "
              f"{t['S']}S + {t['A']}A  (#mon={t['distinct_monomials']})")

    print(f"\nWall-clock ranking (int):")
    for i, r in enumerate(B["results_int"], 1):
        print(f"  {i}. {r['exp7_id']}  "
              f"add={r['t_add_us_median']:.2f}  "
              f"dbl={r['t_dbl_us_median']:.2f}  "
              f"sel={r['t_sel_us_median']:.2f} us  "
              f"(x comp add "
              f"{r['t_add_us_median']/B['baseline_int']['composite_add_us']:.3f}, "
              f"x jac add "
              f"{r['t_add_us_median']/B['baseline_int']['jacobian_add_us']:.3f})")

    if B.get("results_gmp"):
        best_g = B["results_gmp"][0]
        print(f"\nWall-clock best (gmpy2): {best_g['exp7_id']}  "
              f"add={best_g['t_add_us_median']:.2f} us  "
              f"(x comp "
              f"{best_g['t_add_us_median']/B['baseline_gmp']['composite_add_us']:.3f})")

    print(f"\nComposite pK baselines (int): "
          f"add {B['baseline_int']['composite_add_us']:.2f} us, "
          f"dbl {B['baseline_int']['composite_dbl_us']:.2f} us")
    print(f"Jacobian baselines    (int): "
          f"add {B['baseline_int']['jacobian_add_us']:.2f} us, "
          f"dbl {B['baseline_int']['jacobian_dbl_us']:.2f} us")
    return 0


# =====================================================================
# Main
# =====================================================================
def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("Phases: verify, profile, analyze, bench, summary")
        return 0
    phase = sys.argv[1]
    if phase == "verify":  return phase_verify()
    if phase == "profile": return phase_profile()
    if phase == "analyze": return phase_analyze()
    if phase == "bench":   return phase_bench()
    if phase == "summary": return phase_summary()
    print(f"Unknown phase: {phase}")
    return 1


if __name__ == "__main__":
    sys.exit(main())