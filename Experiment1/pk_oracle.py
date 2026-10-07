"""Weierstrass (split cubic) oracle and the birational map Phi_D.

Split cubic:  Y^2 = X(X - d)(X - d*lam)
Points are affine pairs (X, Y) or None for the identity O.
Birational map:  Phi_D = (translation by -W_D) o phi_0,
where W_D = phi_0(D).
"""

from __future__ import annotations
from typing import Optional, Tuple
from pk_core import modinv, PKCurve, proj_eq


WPoint = Optional[Tuple[int, int]]


class SplitCubic:
    def __init__(self, d: int, lam: int, p: int):
        self.d = d % p
        self.lam = lam % p
        self.p = p
        self.a2 = (-self.d * (1 + self.lam)) % p
        self.a4 = (self.d * self.d * self.lam) % p

    def is_on_curve(self, P: WPoint) -> bool:
        if P is None:
            return True
        X, Y = P
        p = self.p
        return (Y * Y - (X * X * X + self.a2 * X * X + self.a4 * X)) % p == 0

    def neg(self, P: WPoint) -> WPoint:
        if P is None:
            return None
        X, Y = P
        return (X, (-Y) % self.p)

    def add(self, P: WPoint, Q: WPoint) -> WPoint:
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        X1, Y1 = P
        X2, Y2 = Q
        if X1 == X2:
            if (Y1 + Y2) % p == 0:
                return None
            return self.double(P)
        s = ((Y2 - Y1) * modinv(X2 - X1, p)) % p
        X3 = (s * s - self.a2 - X1 - X2) % p
        Y3 = (-Y1 + s * (X1 - X3)) % p
        return (X3, Y3)

    def double(self, P: WPoint) -> WPoint:
        p = self.p
        if P is None:
            return None
        X, Y = P
        if Y % p == 0:
            return None
        s = ((3 * X * X + 2 * self.a2 * X + self.a4) * modinv(2 * Y, p)) % p
        X3 = (s * s - self.a2 - 2 * X) % p
        Y3 = (-Y + s * (X - X3)) % p
        return (X3, Y3)

    def scalar_mul(self, k: int, P: WPoint) -> WPoint:
        if k < 0:
            return self.scalar_mul(-k, self.neg(P))
        R: WPoint = None
        Q: WPoint = P
        while k > 0:
            if k & 1:
                R = self.add(R, Q)
            Q = self.double(Q)
            k >>= 1
        return R


def curve_params(curve: PKCurve) -> Tuple[int, int, int]:
    p = curve.p
    u, v = curve.u, curve.v
    d = ((u - 1) * modinv(u + 1, p)) % p
    lam = ((v * v - u * u) * modinv(1 - u * u, p)) % p
    a = ((u + 1) * modinv(u + v, p)) % p
    return d, lam, a


class MapUndefined(Exception):
    pass


def phi_0(P, curve: PKCurve, split: SplitCubic) -> WPoint:
    p = curve.p
    a, b, c = P
    if a % p == 0:
        raise MapUndefined("phi_0: alpha = 0")
    ai = modinv(a, p)
    x = (b * ai) % p
    y = (c * ai) % p

    d, lam, aa = curve_params(curve)

    if (x + 1) % p == 0:
        return None
    if (y + 1) % p == 0:
        return (d, 0)

    r = ((x - 1) * modinv(x + 1, p)) % p
    s = ((y - 1) * modinv(y + 1, p)) % p

    if s % p == 0:
        return ((d * lam) % p, 0)

    t = (r * modinv(s, p)) % p
    T = (t * modinv(aa, p)) % p
    X = (d * T) % p
    Y = (d * s * T * ((T - 1) % p)) % p
    return (X, Y)


def phi_0_inv(W: WPoint, curve: PKCurve, split: SplitCubic):
    p = curve.p
    d, lam, aa = curve_params(curve)
    dd = (d * lam) % p

    if W is None:
        return (1, -1, 1)

    X, Y = W
    if Y % p == 0:
        if X % p == 0:
            return (-1, -1, 1)
        if (X - d) % p == 0:
            return (-1, 1, 1)
        if (X - dd) % p == 0:
            return (1, 1, 1)
        raise MapUndefined("phi_0_inv: Y = 0 off 2-torsion")

    if X % p == 0 or (X - d) % p == 0:
        raise MapUndefined("phi_0_inv: X at pole")

    s = (d * Y * modinv(X * (X - d), p)) % p
    r = (aa * Y * modinv(X - d, p)) % p

    if (1 - r) % p == 0 or (1 - s) % p == 0:
        raise MapUndefined("phi_0_inv: r or s at pole")

    x = ((1 + r) * modinv(1 - r, p)) % p
    y = ((1 + s) * modinv(1 - s, p)) % p
    return (1, x, y)


class Oracle:
    def __init__(self, curve: PKCurve):
        d, lam, _ = curve_params(curve)
        self.curve = curve
        self.split = SplitCubic(d, lam, curve.p)
        self.W_D = phi_0(curve.D, curve, self.split)
        if self.W_D is None or not self.split.is_on_curve(self.W_D):
            raise RuntimeError("phi_0(D) must be finite and on the split cubic")

    def to_w(self, P) -> WPoint:
        W0 = phi_0(P, self.curve, self.split)
        return self.split.add(W0, self.split.neg(self.W_D))

    def from_w(self, W: WPoint):
        Ws = self.split.add(W, self.W_D)
        return phi_0_inv(Ws, self.curve, self.split)

    def oracle_add(self, P, Q):
        return self.from_w(self.split.add(self.to_w(P), self.to_w(Q)))

    def oracle_double(self, P):
        return self.from_w(self.split.double(self.to_w(P)))

    def oracle_neg(self, P):
        return self.from_w(self.split.neg(self.to_w(P)))