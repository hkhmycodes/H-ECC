"""Core arithmetic for the pK cubic over F_p."""

from __future__ import annotations
from typing import Optional, Tuple


Point = Tuple[int, int, int]


# ---------------------------------------------------------------------------
# Modular helpers
# ---------------------------------------------------------------------------

def modinv(a: int, p: int) -> int:
    a %= p
    if a == 0:
        raise ZeroDivisionError("modular inverse of zero")
    return pow(a, p - 2, p)


def is_prime(n: int, witnesses=None) -> bool:
    """Deterministic Miller-Rabin for the sizes we use."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False

    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    if witnesses is None:
        witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
                     41, 43, 47, 53, 59, 61, 67, 71]

    for a in witnesses:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def is_zero_point(P: Point, p: int) -> bool:
    return all(c % p == 0 for c in P)


def proj_eq(P: Point, Q: Point, p: int) -> bool:
    """P == Q in P^2(F_p), via vanishing cross product.

    Returns False if either input is the zero vector (not a valid projective
    point).
    """
    if is_zero_point(P, p) or is_zero_point(Q, p):
        return False
    x0, x1, x2 = P
    y0, y1, y2 = Q
    return (
        (x1 * y2 - x2 * y1) % p == 0
        and (x2 * y0 - x0 * y2) % p == 0
        and (x0 * y1 - x1 * y0) % p == 0
    )


def proj_normalize(P: Point, p: int) -> Point:
    for c in P:
        if c % p != 0:
            ci = modinv(c, p)
            return tuple((x * ci) % p for x in P)  # type: ignore
    return P


def are_collinear(P: Point, Q: Point, R: Point, p: int) -> bool:
    """Determinant test for collinearity of three projective points."""
    return (
        (P[0] * (Q[1] * R[2] - Q[2] * R[1])
         - P[1] * (Q[0] * R[2] - Q[2] * R[0])
         + P[2] * (Q[0] * R[1] - Q[1] * R[0])) % p
    ) == 0


# ---------------------------------------------------------------------------
# Curve
# ---------------------------------------------------------------------------

def check_curve_params(u: int, v: int, w: int, p: int, q: int) -> None:
    """Raise ValueError if the parameter set violates a required assumption."""
    if not is_prime(p):
        raise ValueError("p is not prime")
    if not is_prime(q):
        raise ValueError("q is not prime")
    w %= p
    if w == 0:
        raise ValueError("w == 0 mod p")
    wi = modinv(w, p)
    u, v = (u * wi) % p, (v * wi) % p
    if (u + 1) % p == 0:
        raise ValueError("u + 1 == 0 (d undefined)")
    if (1 - u * u) % p == 0:
        raise ValueError("1 - u^2 == 0 (lambda undefined)")
    if (u + v) % p == 0:
        raise ValueError("u + v == 0 (a undefined)")
    if (u * u - v * v) % p == 0 or (v * v - 1) % p == 0 or (1 - u * u) % p == 0:
        raise ValueError("curve is singular")


class PKCurve:
    """The pK cubic E, normalized so that w = 1."""

    def __init__(self, u: int, v: int, w: int, p: int):
        w %= p
        if w == 0:
            raise ValueError("w must be nonzero")
        wi = modinv(w, p)
        self.u = (u * wi) % p
        self.v = (v * wi) % p
        self.w = 1
        self.p = p
        self.D = (self.u, self.v, 1)
        self.Ds = (self.v, self.u, (self.u * self.v) % p)

    def F(self, P: Point) -> int:
        u, v, w, p = self.u, self.v, self.w, self.p
        a, b, c = P
        return (u * a * (b * b - c * c)
                + v * b * (c * c - a * a)
                + w * c * (a * a - b * b)) % p

    def grad_F(self, P: Point) -> Point:
        u, v, w, p = self.u, self.v, self.w, self.p
        a, b, c = P
        Fa = (u * (b * b - c * c) + 2 * a * (w * c - v * b)) % p
        Fb = (v * (c * c - a * a) + 2 * b * (u * a - w * c)) % p
        Fc = (w * (a * a - b * b) + 2 * c * (v * b - u * a)) % p
        return (Fa, Fb, Fc)

    def dot(self, P: Point, Q: Point) -> int:
        return (P[0] * Q[0] + P[1] * Q[1] + P[2] * Q[2]) % self.p

    def is_on_curve(self, P: Point) -> bool:
        return self.F(P) == 0

    def is_base_locus(self, P: Point) -> bool:
        """True iff P is (projectively) one of A=(1:0:0), B=(0:1:0), C=(0:0:1)."""
        p = self.p
        if is_zero_point(P, p):
            return False
        return sum(1 for c in P if c % p == 0) >= 2

    def isogonal(self, P: Point) -> Point:
        a, b, c = P
        p = self.p
        return ((b * c) % p, (c * a) % p, (a * b) % p)

    def neg(self, P: Point) -> Point:
        """-P.

        Generic case: -P = A_P P - B_P D* with
            A_P = P . grad F(D*),  B_P = D* . grad F(P).
        Exceptional case: D* is 2-torsion, so -D* = D*.
        """
        if proj_eq(P, self.Ds, self.p):
            return self.Ds
        Ap = self.dot(P, self.grad_F(self.Ds))
        Bp = self.dot(self.Ds, self.grad_F(P))
        p = self.p
        return tuple((Ap * P[i] - Bp * self.Ds[i]) % p for i in range(3))  # type: ignore

    def third_intersection(self, P: Point, Q: Point) -> Point:
        """Third intersection of line PQ with E (secant formula).

        Does NOT apply the isogonal map.  Returns the point R such that
        P, Q, R are collinear.  Does not check for base locus.
        """
        Ap = self.dot(P, self.grad_F(Q))
        Bp = self.dot(Q, self.grad_F(P))
        p = self.p
        return tuple((Ap * P[i] - Bp * Q[i]) % p for i in range(3))  # type: ignore

    def add(self, P: Point, Q: Point) -> Point:
        if proj_eq(P, Q, self.p):
            return self.double(P)
        R = self.third_intersection(P, Q)
        if is_zero_point(R, self.p):
            raise ValueError("secant third point is zero vector")
        if self.is_base_locus(R):
            raise ValueError("secant third point lands in base locus")
        return self.isogonal(R)

    def double(self, P: Point) -> Point:
        p = self.p
        a, b, c = P
        if is_zero_point(P, p):
            raise ValueError("zero point")
        Fa, Fb, Fc = self.grad_F(P)
        if a % p != 0:
            t = (0, Fc, (-Fb) % p)
        elif b % p != 0:
            t = (Fc, 0, (-Fa) % p)
        elif c % p != 0:
            t = ((-Fb) % p, Fa, 0)
        else:
            raise ValueError("zero point")

        H = self._quad_Hcal(P, t)
        Ft = self.F(t)
        R = tuple(((-Ft) * P[i] + H * t[i]) % p for i in range(3))
        if is_zero_point(R, p):
            raise ValueError("tangent third point is zero vector")
        if self.is_base_locus(R):
            raise ValueError("tangent third point lands in base locus")
        return self.isogonal(R)

    def _quad_Hcal(self, P: Point, t: Point) -> int:
        u, v, w, p = self.u, self.v, self.w, self.p
        a, b, c = P
        t0, t1, t2 = t
        H00 = (w * c - v * b) % p
        H11 = (u * a - w * c) % p
        H22 = (v * b - u * a) % p
        H01 = (u * b - v * a) % p
        H02 = (w * a - u * c) % p
        H12 = (v * c - w * b) % p
        return (H00 * t0 * t0 + H11 * t1 * t1 + H22 * t2 * t2
                + 2 * H01 * t0 * t1
                + 2 * H02 * t0 * t2
                + 2 * H12 * t1 * t2) % p