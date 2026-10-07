#!/usr/bin/env python3
"""EXP-2 / EXP-6: extraction of the nine-dimensional (2,3)/(3,2)
addition-law basis for the 256-bit test vector of the manuscript.

Runs the full extraction pipeline:
    1. Sanity checks (F(P), F(Q) reduce to zero; S has bidegree (4,4);
       F(S) reduces to zero modulo (F(P), F(Q)); matrix shape).
    2. Assemble M_{2,3}  (756 x 162 over F_p).
    3. RREF; nullspace.
    4. Certify each law by exact polynomial reduction.
    5. Transpose to the (3,2) basis by explicit coefficient permutation.

Independent Sage certification is provided separately in sage_verify.sage.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.setrecursionlimit(100000)

P_MOD = 106839527430202782610735740077697524141755216002708025221579380756122101340723
U, V, W = 5, 1, 7
OUTDIR = Path(__file__).parent / "results"

Mono3 = Tuple[int, int, int]
Mono6 = Tuple[int, int, int, int, int, int]

# Substitution used for reduction modulo F = 0, on either (a,b,c) or (x,y,z):
#     first-coord^2 * second-coord  =  5 * first * second^2
#                                    - 5 * first * third^2
#                                    + second * third^2
#                                    + 7 * first^2 * third
#                                    - 7 * second^2 * third
_REL = [(1, 2, 0, 5), (1, 0, 2, -5), (0, 1, 2, 1),
        (2, 0, 1, 7), (0, 2, 1, -7)]


# ---------------------------------------------------------------------------
# Monomial reduction with memoization
# ---------------------------------------------------------------------------

_P_cache: Dict[Mono3, Dict[Mono3, int]] = {}


def reduce_P(m: Mono3) -> Dict[Mono3, int]:
    """Reduce a monomial (i,j,k) in (a,b,c) modulo a^2 b."""
    if m in _P_cache:
        return _P_cache[m]
    i, j, k = m
    if i < 2 or j < 1:
        r = {m: 1}
    else:
        r: Dict[Mono3, int] = {}
        for da, db, dc, c in _REL:
            nm = (i - 2 + da, j - 1 + db, k + dc)
            for mm, cc in reduce_P(nm).items():
                r[mm] = (r.get(mm, 0) + c * cc) % P_MOD
        r = {mm: cc for mm, cc in r.items() if cc}
    _P_cache[m] = r
    return r


_Q_cache: Dict[Mono3, Dict[Mono3, int]] = {}


def reduce_Q(m: Mono3) -> Dict[Mono3, int]:
    """Reduce a monomial (i,j,k) in (x,y,z) modulo x^2 y."""
    if m in _Q_cache:
        return _Q_cache[m]
    i, j, k = m
    if i < 2 or j < 1:
        r = {m: 1}
    else:
        r: Dict[Mono3, int] = {}
        for da, db, dc, c in _REL:
            nm = (i - 2 + da, j - 1 + db, k + dc)
            for mm, cc in reduce_Q(nm).items():
                r[mm] = (r.get(mm, 0) + c * cc) % P_MOD
        r = {mm: cc for mm, cc in r.items() if cc}
    _Q_cache[m] = r
    return r


# ---------------------------------------------------------------------------
# Polynomial arithmetic in 6 variables
# ---------------------------------------------------------------------------

def p_add(a, b):
    r = dict(a)
    for m, c in b.items():
        v = (r.get(m, 0) + c) % P_MOD
        if v:
            r[m] = v
        elif m in r:
            del r[m]
    return r


def p_sub(a, b):
    r = dict(a)
    for m, c in b.items():
        v = (r.get(m, 0) - c) % P_MOD
        if v:
            r[m] = v
        elif m in r:
            del r[m]
    return r


def p_scale(a, s):
    s %= P_MOD
    if s == 0:
        return {}
    return {m: (c * s) % P_MOD for m, c in a.items() if (c * s) % P_MOD}


def p_mul(a, b):
    r = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = tuple(x + y for x, y in zip(ma, mb))
            r[m] = (r.get(m, 0) + ca * cb) % P_MOD
    return {m: c for m, c in r.items() if c}


def var(i):
    m = [0] * 6
    m[i] = 1
    return {tuple(m): 1}


# ---------------------------------------------------------------------------
# Curve, gradient, composite S
# ---------------------------------------------------------------------------

def F_P():
    a, b, c = var(0), var(1), var(2)
    a2, b2, c2 = p_mul(a, a), p_mul(b, b), p_mul(c, c)
    return p_add(
        p_add(p_scale(p_mul(a, p_sub(b2, c2)), U),
              p_scale(p_mul(b, p_sub(c2, a2)), V)),
        p_scale(p_mul(c, p_sub(a2, b2)), W))


def F_Q():
    x, y, z = var(3), var(4), var(5)
    x2, y2, z2 = p_mul(x, x), p_mul(y, y), p_mul(z, z)
    return p_add(
        p_add(p_scale(p_mul(x, p_sub(y2, z2)), U),
              p_scale(p_mul(y, p_sub(z2, x2)), V)),
        p_scale(p_mul(z, p_sub(x2, y2)), W))


def grad_F_P():
    a, b, c = var(0), var(1), var(2)
    a2, b2, c2 = p_mul(a, a), p_mul(b, b), p_mul(c, c)
    ab, ac, bc = p_mul(a, b), p_mul(a, c), p_mul(b, c)
    Fa = p_add(p_scale(p_sub(b2, c2), U),
               p_scale(p_mul(a, p_sub(p_scale(c, W), p_scale(b, V))), 2))
    Fb = p_add(p_add(p_scale(ab, 2 * U), p_scale(p_sub(c2, a2), V)),
               p_scale(bc, -2 * W))
    Fc = p_add(p_add(p_scale(ac, -2 * U), p_scale(bc, 2 * V)),
               p_scale(p_sub(a2, b2), W))
    return Fa, Fb, Fc


def grad_F_Q():
    x, y, z = var(3), var(4), var(5)
    x2, y2, z2 = p_mul(x, x), p_mul(y, y), p_mul(z, z)
    xy, xz, yz = p_mul(x, y), p_mul(x, z), p_mul(y, z)
    Fx = p_add(p_scale(p_sub(y2, z2), U),
               p_scale(p_mul(x, p_sub(p_scale(z, W), p_scale(y, V))), 2))
    Fy = p_add(p_add(p_scale(xy, 2 * U), p_scale(p_sub(z2, x2), V)),
               p_scale(yz, -2 * W))
    Fz = p_add(p_add(p_scale(xz, -2 * U), p_scale(yz, 2 * V)),
               p_scale(p_sub(x2, y2), W))
    return Fx, Fy, Fz


def compute_S():
    a, b, c = var(0), var(1), var(2)
    x, y, z = var(3), var(4), var(5)
    Fx, Fy, Fz = grad_F_Q()
    Fa, Fb, Fc = grad_F_P()
    A = p_add(p_add(p_mul(a, Fx), p_mul(b, Fy)), p_mul(c, Fz))
    B = p_add(p_add(p_mul(x, Fa), p_mul(y, Fb)), p_mul(z, Fc))
    R0 = p_sub(p_mul(A, a), p_mul(B, x))
    R1 = p_sub(p_mul(A, b), p_mul(B, y))
    R2 = p_sub(p_mul(A, c), p_mul(B, z))
    return p_mul(R1, R2), p_mul(R2, R0), p_mul(R0, R1)


# ---------------------------------------------------------------------------
# Bidegree and nonzero tests
# ---------------------------------------------------------------------------

def pdeg(poly):
    dp = dq = 0
    for m in poly:
        dp = max(dp, sum(m[:3]))
        dq = max(dq, sum(m[3:]))
    return dp, dq


def is_nonzero(poly):
    return any(c % P_MOD != 0 for c in poly.values())


# ---------------------------------------------------------------------------
# Reduction of any 6-var polynomial modulo (F(P), F(Q))
# ---------------------------------------------------------------------------

def reduce_poly6(poly):
    result: Dict[Mono6, int] = {}
    for m, c in poly.items():
        for mp, cp in reduce_P(m[:3]).items():
            for mq, cq in reduce_Q(m[3:]).items():
                mm = mp + mq
                result[mm] = (result.get(mm, 0) + c * cp * cq) % P_MOD
    return {m: c for m, c in result.items() if c}


def is_zero_mod_ideal(poly):
    return not reduce_poly6(poly)


# ---------------------------------------------------------------------------
# Reduced bases for degree 6 in P and degree 7 in Q
# ---------------------------------------------------------------------------

def mons3(d):
    return [(i, j, d - i - j) for i in range(d + 1) for j in range(d - i + 1)]


RED_P_6 = sorted(m for m in mons3(6) if not (m[0] >= 2 and m[1] >= 1))
RED_Q_7 = sorted(m for m in mons3(7) if not (m[0] >= 2 and m[1] >= 1))

assert len(RED_P_6) == 18, len(RED_P_6)
assert len(RED_Q_7) == 21, len(RED_Q_7)

IDX_P_6 = {m: i for i, m in enumerate(RED_P_6)}
IDX_Q_7 = {m: i for i, m in enumerate(RED_Q_7)}

R_P_6 = {m: {IDX_P_6[mm]: c for mm, c in reduce_P(m).items()}
         for m in mons3(6)}
R_Q_7 = {m: {IDX_Q_7[mm]: c for mm, c in reduce_Q(m).items()}
         for m in mons3(7)}


# ---------------------------------------------------------------------------
# Candidate-space bases
# ---------------------------------------------------------------------------

# Six quadratic monomials in the first point
Q_P = [(2, 0, 0), (1, 1, 0), (1, 0, 1),
       (0, 2, 0), (0, 1, 1), (0, 0, 2)]

# Nine reduced cubic monomials in the second point
C_Q = [(0, 0, 3), (0, 1, 2), (0, 2, 1), (0, 3, 0),
       (1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1), (3, 0, 0)]


def reduced_vector(S_i, mu, nu):
    """Reduced 378-vector of S_i · μ · ν modulo (F(P), F(Q))."""
    n_Q = len(RED_Q_7)  # 21
    vec = [0] * (len(RED_P_6) * n_Q)
    for mono6, s in S_i.items():
        i, j, k = mono6[:3]
        l, m, n = mono6[3:]
        Pa = (i + mu[0], j + mu[1], k + mu[2])
        Qa = (l + nu[0], m + nu[1], n + nu[2])
        if sum(Pa) != 6 or sum(Qa) != 7:
            raise ValueError(
                f"Unexpected monomial degree in S·μ·ν: "
                f"sum(Pa)={sum(Pa)}, sum(Qa)={sum(Qa)}, "
                f"Pa={Pa}, Qa={Qa}")
        for ip, cp in R_P_6[Pa].items():
            base = ip * n_Q
            for iq, cq in R_Q_7[Qa].items():
                vec[base + iq] = (vec[base + iq] + s * cp * cq) % P_MOD
    return vec


# ---------------------------------------------------------------------------
# Matrix assembly
# ---------------------------------------------------------------------------

def build_M(S_list):
    S0, S1, S2 = S_list
    n_red = len(RED_P_6) * len(RED_Q_7)  # 378
    M = [[0] * 162 for _ in range(2 * n_red)]
    for mu_idx, mu in enumerate(Q_P):
        for nu_idx, nu in enumerate(C_Q):
            base = mu_idx * 9 + nu_idx
            v0 = reduced_vector(S0, mu, nu)
            v1 = reduced_vector(S1, mu, nu)
            v2 = reduced_vector(S2, mu, nu)
            # column (0, μν): +S_1 into G_1 rows; +S_2 into G_2 rows
            col = 0 * 54 + base
            for r in range(n_red):
                if v1[r]:
                    M[r][col] = (M[r][col] + v1[r]) % P_MOD
                if v2[r]:
                    M[n_red + r][col] = (M[n_red + r][col] + v2[r]) % P_MOD
            # column (1, μν): -S_0 into G_1 rows
            col = 1 * 54 + base
            for r in range(n_red):
                if v0[r]:
                    M[r][col] = (M[r][col] - v0[r]) % P_MOD
            # column (2, μν): -S_0 into G_2 rows
            col = 2 * 54 + base
            for r in range(n_red):
                if v0[r]:
                    M[n_red + r][col] = (M[n_red + r][col] - v0[r]) % P_MOD
    return M


# ---------------------------------------------------------------------------
# RREF and nullspace over F_p
# ---------------------------------------------------------------------------

def rref(M, p):
    M = [list(r) for r in M]
    m, n = len(M), len(M[0])
    row = 0
    pivots = []
    for col in range(n):
        if row >= m:
            break
        pr = next((r for r in range(row, m) if M[r][col]), -1)
        if pr == -1:
            continue
        M[row], M[pr] = M[pr], M[row]
        inv = pow(M[row][col], -1, p)
        M[row] = [(x * inv) % p for x in M[row]]
        for r in range(m):
            if r != row and M[r][col]:
                f = M[r][col]
                M[r] = [(M[r][c] - f * M[row][c]) % p for c in range(n)]
        pivots.append((row, col))
        row += 1
    return M, pivots, row


def nullspace(M_rref, pivots, n, p):
    pivot_cols = {c for _, c in pivots}
    basis = []
    for fc in [c for c in range(n) if c not in pivot_cols]:
        v = [0] * n
        v[fc] = 1
        for pr, pc in pivots:
            v[pc] = (-M_rref[pr][fc]) % p
        basis.append(v)
    return basis


# ---------------------------------------------------------------------------
# Certification
# ---------------------------------------------------------------------------

def certify(coeffs, S_list):
    """Verify G_1 = S_1 Z_0 - S_0 Z_1 ≡ 0 and G_2 = S_2 Z_0 - S_0 Z_2 ≡ 0
    modulo (F(P), F(Q)).  Raises on any unexpected monomial degree."""
    S0, S1, S2 = S_list
    Z = [{}, {}, {}]
    for i in range(3):
        for mu_idx, mu in enumerate(Q_P):
            for nu_idx, nu in enumerate(C_Q):
                c = coeffs[i * 54 + mu_idx * 9 + nu_idx] % P_MOD
                if c:
                    m = (mu[0], mu[1], mu[2], nu[0], nu[1], nu[2])
                    Z[i][m] = c
    G1 = p_sub(p_mul(S1, Z[0]), p_mul(S0, Z[1]))
    G2 = p_sub(p_mul(S2, Z[0]), p_mul(S0, Z[2]))
    for G in (G1, G2):
        for m in G:
            if sum(m[:3]) != 6 or sum(m[3:]) != 7:
                raise ValueError(
                    f"Unexpected monomial degree in certificate: "
                    f"P-deg={sum(m[:3])}, Q-deg={sum(m[3:])}, mono={m}")
    return is_zero_mod_ideal(G1) and is_zero_mod_ideal(G2)


# ---------------------------------------------------------------------------
# (3,2) transposition
# ---------------------------------------------------------------------------

def transpose_23_to_32(c):
    """A_{3,2}(P,Q) = A_{2,3}(Q,P) written in the (3,2) monomial basis.

    (2,3) layout: col = i*54 + mu_idx*9 + nu_idx  (mu ∈ Q_P, nu ∈ C_Q)
    (3,2) layout: col = i*54 + nu_idx*6 + mu_idx  (nu ∈ C_P, mu ∈ Q_Q)

    Under (a,b,c) ↔ (x,y,z), the monomial index sets Q_P and Q_Q are the
    same, and so are C_Q and C_P. The transposition is therefore a pure
    permutation of coefficients.
    """
    out = [0] * 162
    for i in range(3):
        for mu_idx in range(6):
            for nu_idx in range(9):
                out[i * 54 + nu_idx * 6 + mu_idx] = \
                    c[i * 54 + mu_idx * 9 + nu_idx]
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"[EXP-2/6] p = {P_MOD} ({P_MOD.bit_length()} bits)")
    print(f"[EXP-2/6] (u:v:w) = ({U}:{V}:{W})")
    print()

    # --- sanity checks -------------------------------------------------
    print("[check] F(P) reduces to zero ...")
    assert is_zero_mod_ideal(F_P()), "F(P) does not reduce to zero"

    print("[check] F(Q) reduces to zero ...")
    assert is_zero_mod_ideal(F_Q()), "F(Q) does not reduce to zero"

    print("[check] computing S ...")
    S_list = compute_S()
    S0, S1, S2 = S_list
    for name, S in zip(("S_0", "S_1", "S_2"), S_list):
        dp, dq = pdeg(S)
        print(f"        bideg({name}) = ({dp}, {dq})")
        assert (dp, dq) == (4, 4), f"{name} has wrong bidegree"
        assert is_nonzero(S), f"{name} is zero"

    print("[check] F(S) reduces to zero modulo (F(P), F(Q)) ...")
    FS = p_add(
        p_add(p_scale(p_mul(S0, p_sub(p_mul(S1, S1), p_mul(S2, S2))), U),
              p_scale(p_mul(S1, p_sub(p_mul(S2, S2), p_mul(S0, S0))), V)),
        p_scale(p_mul(S2, p_sub(p_mul(S0, S0), p_mul(S1, S1))), W))
    dp, dq = pdeg(FS)
    print(f"        bideg(F(S)) = ({dp}, {dq})")
    assert is_zero_mod_ideal(FS), "F(S) does not reduce to zero"

    # --- matrix --------------------------------------------------------
    print("[EXP-2/6] assembling M_{2,3} ...")
    M = build_M(S_list)
    shape = (len(M), len(M[0]))
    print(f"        matrix shape: {shape[0]} x {shape[1]}")
    assert shape == (756, 162), f"unexpected matrix shape {shape}"

    print("[EXP-2/6] RREF ...")
    M_rref, pivots, rank = rref(M, P_MOD)
    nullity = 162 - rank
    print(f"        rank = {rank}, nullity = {nullity}")

    print("[EXP-2/6] nullspace ...")
    basis = nullspace(M_rref, pivots, 162, P_MOD)
    print(f"        basis size = {len(basis)}")

    print("[EXP-2/6] certifying each (2,3) law ...")
    certs = [certify(c, S_list) for c in basis]
    for j, ok in enumerate(certs):
        print(f"        law {j:02d}: {'OK' if ok else 'FAIL'}")

    print("[EXP-2/6] transposing to (3,2) ...")
    basis_32 = [transpose_23_to_32(c) for c in basis]

    # --- outputs -------------------------------------------------------
    OUTDIR.mkdir(exist_ok=True)

    parameters = {
        "p": P_MOD, "U": U, "V": V, "W": W,
        "Q_P": [list(m) for m in Q_P],
        "C_Q": [list(m) for m in C_Q],
        "RED_P_6": [list(m) for m in RED_P_6],
        "RED_Q_7": [list(m) for m in RED_Q_7],
        "matrix_shape": list(shape),
        "matrix_rank": rank,
        "nullity": nullity,
    }
    (OUTDIR / "parameters.json").write_text(json.dumps(parameters, indent=2))

    basis_23_json = {
        "law_count": len(basis),
        "column_layout": "col = i*54 + mu_idx*9 + nu_idx  (mu ∈ Q_P, nu ∈ C_Q)",
        "Q_P": [list(m) for m in Q_P],
        "C_Q": [list(m) for m in C_Q],
        "laws": [{"index": j, "coeffs": [int(x) for x in c]}
                 for j, c in enumerate(basis)],
    }
    (OUTDIR / "basis_23.json").write_text(json.dumps(basis_23_json, indent=2))

    basis_32_json = {
        "law_count": len(basis_32),
        "column_layout": "col = i*54 + nu_idx*6 + mu_idx  (nu ∈ C_P, mu ∈ Q_Q)",
        "C_P": [[m[2], m[1], m[0]] for m in C_Q],
        "Q_Q": [[m[2], m[1], m[0]] for m in Q_P],
        "note": "A_{3,2}(P,Q) = A_{2,3}(Q,P) expressed in the (3,2) basis",
        "laws": [{"index": j, "coeffs": [int(x) for x in c]}
                 for j, c in enumerate(basis_32)],
    }
    (OUTDIR / "basis_32.json").write_text(json.dumps(basis_32_json, indent=2))

    (OUTDIR / "matrix_rank.txt").write_text(
        f"matrix shape: {shape[0]} x {shape[1]}\n"
        f"rank        : {rank}\n"
        f"nullity     : {nullity}\n")

    cert_dir = OUTDIR / "certificates"
    cert_dir.mkdir(exist_ok=True)
    for j, ok in enumerate(certs):
        (cert_dir / f"law_{j:02d}.txt").write_text(
            f"law {j:02d}\n"
            f"certified (S_1 Z_0 - S_0 Z_1 = 0, S_2 Z_0 - S_0 Z_2 = 0 "
            f"mod (F(P), F(Q))): {ok}\n")

    success = (len(basis) == 9 and all(certs))
    print()
    print("=" * 60)
    print(f"matrix shape    : {shape[0]} x {shape[1]}")
    print(f"rank            : {rank}")
    print(f"nullity         : {nullity}")
    print(f"laws found      : {len(basis)}")
    print(f"all certified   : {all(certs)}")
    print(f"output          : {OUTDIR}/")
    print("=" * 60)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())