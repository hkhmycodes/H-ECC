"""pK arithmetic on E:  uα(β²−γ²) + vβ(γ²−α²) + wγ(α²−β²) = 0.

Same mathematics as `pk_core.py` in Experiment 1, written against the
tracked `field` primitives.  No CSE optimisation is applied on purpose:
this file is the Stage-I formula baseline that EXP-3 measures.
"""

from __future__ import annotations
from typing import Optional, Tuple

import field as f


Point = Tuple[int, int, int]
UVW = Tuple[int, int, int]


def pk_F(P: Point, uvw: UVW, p: int):
    u, v, w = uvw
    a, b, c = P
    a2 = f.sqr(a, p)
    b2 = f.sqr(b, p)
    c2 = f.sqr(c, p)
    t1 = f.mul(f.cmul(a, u, p), f.sub(b2, c2, p), p)
    t2 = f.mul(f.cmul(b, v, p), f.sub(c2, a2, p), p)
    t3 = f.mul(f.cmul(c, w, p), f.sub(a2, b2, p), p)
    return f.add(f.add(t1, t2, p), t3, p)


def pk_grad(P: Point, uvw: UVW, p: int) -> Point:
    u, v, w = uvw
    a, b, c = P
    a2 = f.sqr(a, p)
    b2 = f.sqr(b, p)
    c2 = f.sqr(c, p)
    ua = f.cmul(a, u, p)
    vb = f.cmul(b, v, p)
    wc = f.cmul(c, w, p)

    Fa = f.cmul(f.sub(b2, c2, p), u, p)
    t = f.mul(a, f.sub(wc, vb, p), p)
    Fa = f.add(Fa, f.smul(t, 2, p), p)

    Fb = f.cmul(f.sub(c2, a2, p), v, p)
    t = f.mul(b, f.sub(ua, wc, p), p)
    Fb = f.add(Fb, f.smul(t, 2, p), p)

    Fc = f.cmul(f.sub(a2, b2, p), w, p)
    t = f.mul(c, f.sub(vb, ua, p), p)
    Fc = f.add(Fc, f.smul(t, 2, p), p)

    return (Fa, Fb, Fc)


def pk_dot(P: Point, Q: Point, p: int):
    return f.add(f.add(f.mul(P[0], Q[0], p),
                       f.mul(P[1], Q[1], p), p),
                 f.mul(P[2], Q[2], p), p)


def pk_isogonal(P: Point, p: int) -> Point:
    a, b, c = P
    return (f.mul(b, c, p), f.mul(c, a, p), f.mul(a, b, p))


def pk_add(P: Point, Q: Point, uvw: UVW, p: int) -> Point:
    gQ = pk_grad(Q, uvw, p)
    gP = pk_grad(P, uvw, p)
    A = pk_dot(P, gQ, p)
    B = pk_dot(Q, gP, p)
    R = (
        f.sub(f.mul(A, P[0], p), f.mul(B, Q[0], p), p),
        f.sub(f.mul(A, P[1], p), f.mul(B, Q[1], p), p),
        f.sub(f.mul(A, P[2], p), f.mul(B, Q[2], p), p),
    )
    return pk_isogonal(R, p)


def pk_double(P: Point, uvw: UVW, p: int) -> Point:
    u, v, w = uvw
    a, b, c = P
    Fa, Fb, Fc = pk_grad(P, uvw, p)

    if a % p != 0:
        t = (0, Fc, (-Fb) % p)
    elif b % p != 0:
        t = (Fc, 0, (-Fa) % p)
    else:
        t = ((-Fb) % p, Fa, 0)
    t0, t1, t2 = t

    H00 = f.sub(f.cmul(c, w, p), f.cmul(b, v, p), p)
    H11 = f.sub(f.cmul(a, u, p), f.cmul(c, w, p), p)
    H22 = f.sub(f.cmul(b, v, p), f.cmul(a, u, p), p)
    H01 = f.sub(f.cmul(b, u, p), f.cmul(a, v, p), p)
    H02 = f.sub(f.cmul(a, w, p), f.cmul(c, u, p), p)
    H12 = f.sub(f.cmul(c, v, p), f.cmul(b, w, p), p)

    H = f.add(
        f.add(f.mul(H00, f.sqr(t0, p), p),
              f.mul(H11, f.sqr(t1, p), p), p),
        f.mul(H22, f.sqr(t2, p), p), p)
    H = f.add(H, f.smul(f.mul(H01, f.mul(t0, t1, p), p), 2, p), p)
    H = f.add(H, f.smul(f.mul(H02, f.mul(t0, t2, p), p), 2, p), p)
    H = f.add(H, f.smul(f.mul(H12, f.mul(t1, t2, p), p), 2, p), p)

    Ft = pk_F(t, uvw, p)
    mFt = (-Ft) % p
    R = (
        f.add(f.mul(mFt, a, p), f.mul(H, t0, p), p),
        f.add(f.mul(mFt, b, p), f.mul(H, t1, p), p),
        f.add(f.mul(mFt, c, p), f.mul(H, t2, p), p),
    )
    return pk_isogonal(R, p)


def pk_neg(P: Point, uvw: UVW, Ds: Point, gDs: Point, p: int) -> Point:
    if P == Ds:
        return Ds
    Ap = pk_dot(P, gDs, p)
    Bp = pk_dot(Ds, pk_grad(P, uvw, p), p)
    return (
        f.sub(f.mul(Ap, P[0], p), f.mul(Bp, Ds[0], p), p),
        f.sub(f.mul(Ap, P[1], p), f.mul(Bp, Ds[1], p), p),
        f.sub(f.mul(Ap, P[2], p), f.mul(Bp, Ds[2], p), p),
    )


def pk_scalar_mul(k: int, P: Point, uvw: UVW, p: int, D: Point) -> Point:
    if k == 0:
        return D
    R: Optional[Point] = None
    Q = P
    while k:
        if k & 1:
            R = Q if R is None else pk_add(R, Q, uvw, p)
        Q = pk_double(Q, uvw, p)
        k >>= 1
    return R  # type: ignore