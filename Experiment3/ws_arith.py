"""Weierstrass–Jacobian arithmetic on E_W : Y² = X³ + aX + b.

Backends:
  - ws_add              : generic Jacobian addition (add-2007-bl)
  - ws_madd             : mixed Jacobian-affine addition (madd-2007-bl)
                          with internal exceptional handling (P = Q → double,
                          P = -Q → identity)
  - ws_double           : Jacobian doubling, general a (dbl-2007-bl)
  - ws_scalar_mul       : right-to-left binary double-and-add, generic add
  - ws_scalar_mul_mixed : left-to-right binary double-and-add, mixed add
"""

from __future__ import annotations
from typing import Optional, Tuple

import field as f


Point = Tuple[int, int, int]
IDENTITY: Point = (0, 1, 0)


def ws_from_affine(x: int, y: int, p: int) -> Point:
    return (x % p, y % p, 1)


def ws_to_affine(P: Point, p: int) -> Optional[Tuple[int, int]]:
    X, Y, Z = P
    if Z % p == 0:
        return None
    Zi = f.inv(Z, p)
    Zi2 = f.sqr(Zi, p)
    x = f.mul(X, Zi2, p)
    y = f.mul(Y, f.mul(Zi2, Zi, p), p)
    return (x, y)


def ws_neg(P: Point, p: int) -> Point:
    X, Y, Z = P
    return (X, (-Y) % p, Z)


def ws_double(P: Point, a: int, p: int) -> Point:
    X1, Y1, Z1 = P
    if Z1 % p == 0 or Y1 % p == 0:
        return IDENTITY

    XX = f.sqr(X1, p)
    YY = f.sqr(Y1, p)
    YYYY = f.sqr(YY, p)
    ZZ = f.sqr(Z1, p)

    t = f.add(X1, YY, p)
    t = f.sqr(t, p)
    t = f.sub(f.sub(t, XX, p), YYYY, p)
    S = f.smul(t, 2, p)

    M = f.add(f.smul(XX, 3, p), f.cmul(f.sqr(ZZ, p), a, p), p)

    T = f.sub(f.sqr(M, p), f.smul(S, 2, p), p)
    X3 = T

    t = f.sub(S, T, p)
    t = f.mul(M, t, p)
    Y3 = f.sub(t, f.smul(YYYY, 8, p), p)

    t = f.add(Y1, Z1, p)
    t = f.sqr(t, p)
    t = f.sub(f.sub(t, YY, p), ZZ, p)
    Z3 = t

    return (X3, Y3, Z3)


def ws_add(P: Point, Q: Point, p: int) -> Point:
    """Generic Jacobian addition.  Assumes P ≠ Q, P ≠ −Q, neither identity."""
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q

    Z1Z1 = f.sqr(Z1, p)
    Z2Z2 = f.sqr(Z2, p)
    U1 = f.mul(X1, Z2Z2, p)
    U2 = f.mul(X2, Z1Z1, p)
    S1 = f.mul(f.mul(Y1, Z2, p), Z2Z2, p)
    S2 = f.mul(f.mul(Y2, Z1, p), Z1Z1, p)

    H = f.sub(U2, U1, p)
    I = f.sqr(f.smul(H, 2, p), p)
    J = f.mul(H, I, p)
    r = f.smul(f.sub(S2, S1, p), 2, p)
    V = f.mul(U1, I, p)

    X3 = f.sub(f.sub(f.sqr(r, p), J, p), f.smul(V, 2, p), p)

    t = f.sub(V, X3, p)
    t = f.mul(r, t, p)
    t2 = f.mul(S1, J, p)
    Y3 = f.sub(t, f.smul(t2, 2, p), p)

    t = f.add(Z1, Z2, p)
    t = f.sqr(t, p)
    t = f.sub(f.sub(t, Z1Z1, p), Z2Z2, p)
    Z3 = f.mul(t, H, p)

    return (X3, Y3, Z3)


def ws_madd(P: Point, Q_affine: Tuple[int, int],
            a: int, p: int) -> Point:
    """Mixed Jacobian-affine addition (madd-2007-bl) with internal
    exceptional handling.

    Cases:
      - P = Q  (H = 0, S2 = Y1)  → doubling
      - P = -Q (H = 0, S2 ≠ Y1)  → identity
      - generic (H ≠ 0)          → normal formula

    The generic path incurs only the zero-test on H, which is already
    computed for the standard formula; the exceptional branch is entered
    only when H = 0.
    """
    X1, Y1, Z1 = P
    x2, y2 = Q_affine

    Z1Z1 = f.sqr(Z1, p)
    U2 = f.mul(x2, Z1Z1, p)
    S2 = f.mul(y2, f.mul(Z1, Z1Z1, p), p)

    H = f.sub(U2, X1, p)
    if H == 0:
        if S2 == Y1:
            return ws_double(P, a, p)
        return IDENTITY

    HH = f.sqr(H, p)
    I = f.smul(HH, 4, p)
    J = f.mul(H, I, p)
    r = f.smul(f.sub(S2, Y1, p), 2, p)
    V = f.mul(X1, I, p)

    X3 = f.sub(f.sub(f.sqr(r, p), J, p), f.smul(V, 2, p), p)

    t = f.sub(V, X3, p)
    t = f.mul(r, t, p)
    t2 = f.mul(Y1, J, p)
    Y3 = f.sub(t, f.smul(t2, 2, p), p)

    Z3 = f.mul(f.smul(Z1, 2, p), H, p)

    return (X3, Y3, Z3)


def ws_scalar_mul(k: int, P: Point, a: int, p: int) -> Point:
    """Right-to-left binary double-and-add using generic ws_add."""
    if k == 0:
        return IDENTITY
    R: Optional[Point] = None
    Q = P
    while k:
        if k & 1:
            R = Q if R is None else ws_add(R, Q, p)
        Q = ws_double(Q, a, p)
        k >>= 1
    return R  # type: ignore


def ws_scalar_mul_mixed(k: int, P_affine: Tuple[int, int],
                        a: int, p: int) -> Point:
    """Left-to-right binary double-and-add using mixed addition where the
    accumulator is combined with the affine base point.

    Exceptional cases are handled inside ws_madd.  The only remaining check
    on this side is the identity test R[2] == 0 on the accumulator, which
    is a modular reduction and a comparison, not a field multiplication.
    """
    if k == 0:
        return IDENTITY
    R: Point = IDENTITY
    for bit in bin(k)[2:]:
        R = ws_double(R, a, p)
        if bit == '1':
            if R[2] % p == 0:
                R = ws_from_affine(P_affine[0], P_affine[1], p)
            else:
                R = ws_madd(R, P_affine, a, p)
    return R