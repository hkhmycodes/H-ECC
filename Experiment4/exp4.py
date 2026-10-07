#!/usr/bin/env python3
"""
EXP-4 smoke-test implementation.

Runs independent batches of k_i * P_i on the composite pK and on the
birationally equivalent Jacobian curve, sequentially and with T worker
processes.  Uses the SAME underlying points via Phi_D, matching the
C++/OpenMP version.

For the paper, compile and run exp4_openmp.cpp.  This script is a
cross-check: it validates the workload generation, the sequential/parallel
consistency, and the JSON schema on Windows without needing MSYS2.

Usage:
    python3 exp4.py --kind both --out results/exp4_smoke.json
"""

import argparse
import json
import multiprocessing as mp
import os
import random
import sys
import time
from pathlib import Path
from statistics import median


# ----------------------------------------------------------------
# Constants (test vector of Section 6.5)
# ----------------------------------------------------------------
P_MOD = int(
    "106839527430202782610735740077697524141755216002708025221579380756122101340723"
)
Q_ORDER = int(
    "26709881857550695652683935019424381035438804000677006305394845189030525335181"
)
G_B = int(
    "43885198696659819331868805857976446729580581201517994933321256281291666494433"
)
G_C = int(
    "52096118237682591038514212999531378342298973921021687132305677862559506645034"
)

INV36 = pow(36, -1, P_MOD)
A_W = (-INV36) % P_MOD
B_W = 0

INV2 = pow(2, -1, P_MOD)
INV3 = pow(3, -1, P_MOD)
INV6 = pow(6, -1, P_MOD)


# ----------------------------------------------------------------
# Arithmetic
# ----------------------------------------------------------------
def pk_add(P, Q):
    a, b, c = P
    x, y, z = Q
    Fa = (5 * (b * b - c * c) + 2 * a * (7 * c - b)) % P_MOD
    Fb = (    (c * c - a * a) + 2 * b * (5 * a - 7 * c)) % P_MOD
    Fc = (7 * (a * a - b * b) + 2 * c * (b - 5 * a)) % P_MOD
    Ga = (5 * (y * y - z * z) + 2 * x * (7 * z - y)) % P_MOD
    Gb = (    (z * z - x * x) + 2 * y * (5 * x - 7 * z)) % P_MOD
    Gc = (7 * (x * x - y * y) + 2 * z * (y - 5 * x)) % P_MOD
    A_ = (a * Ga + b * Gb + c * Gc) % P_MOD
    B_ = (x * Fa + y * Fb + z * Fc) % P_MOD
    R0 = (A_ * a - B_ * x) % P_MOD
    R1 = (A_ * b - B_ * y) % P_MOD
    R2 = (A_ * c - B_ * z) % P_MOD
    return ((R1 * R2) % P_MOD, (R2 * R0) % P_MOD, (R0 * R1) % P_MOD)


def pk_double(P):
    a, b, c = P
    Fa = (5 * (b * b - c * c) + 2 * a * (7 * c - b)) % P_MOD
    Fb = (    (c * c - a * a) + 2 * b * (5 * a - 7 * c)) % P_MOD
    Fc = (7 * (a * a - b * b) + 2 * c * (b - 5 * a)) % P_MOD
    B_ = Fb
    G_ = Fc
    Ft = (G_ * B_ * B_ + 7 * B_ * G_ * G_) % P_MOD
    H  = ((5 * a - 7 * c) * G_ * G_
          - 2 * (c - 7 * b) * B_ * G_
          + (b - 5 * a) * B_ * B_) % P_MOD
    R0 = (-Ft * a) % P_MOD
    R1 = (-Ft * b + H * G_) % P_MOD
    R2 = (-Ft * c - H * B_) % P_MOD
    return ((R1 * R2) % P_MOD, (R2 * R0) % P_MOD, (R0 * R1) % P_MOD)


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
    VV = U1 * I % m
    X3 = (r * r - J - 2 * VV) % m
    Y3 = (r * (VV - X3) - 2 * S1 * J) % m
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


def pk_scalar_mul(k, P):
    D = (5, 1, 7)
    if k == 0:
        return D
    R = D
    Q = P
    while k > 0:
        if k & 1:
            R = pk_add(R, Q)
        Q = pk_double(Q)
        k >>= 1
    return R


def jac_scalar_mul(k, P):
    if k == 0:
        return None
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R = Q if R is None else jac_add(R, Q)
        Q = jac_dbl(Q)
        k >>= 1
    return R


# ----------------------------------------------------------------
# Phi_D: pK -> E_W (affine), for workload generation only
# ----------------------------------------------------------------
def w_add_affine(x1, y1, x2, y2):
    if (x1 - x2) % P_MOD == 0 and (y1 + y2) % P_MOD == 0:
        return None
    if x1 == x2 and y1 == y2:
        num = (3 * x1 * x1 + A_W) % P_MOD
        den = (2 * y1) % P_MOD
        lam = num * pow(den, -1, P_MOD) % P_MOD
    else:
        num = (y2 - y1) % P_MOD
        den = (x2 - x1) % P_MOD
        lam = num * pow(den, -1, P_MOD) % P_MOD
    x3 = (lam * lam - x1 - x2) % P_MOD
    y3 = (lam * (x1 - x3) - y1) % P_MOD
    return (x3, y3)


def phi_D_to_W(p):
    a, b, c = p
    if a == 0:
        return None
    ai = pow(a, -1, P_MOD)
    x = b * ai % P_MOD
    y = c * ai % P_MOD
    if (x + 1) % P_MOD == 0 or (y + 1) % P_MOD == 0:
        return None
    r = (x - 1) * pow((x + 1) % P_MOD, -1, P_MOD) % P_MOD
    s = (y - 1) * pow((y + 1) % P_MOD, -1, P_MOD) % P_MOD
    if s == 0:
        return None
    t = r * pow(s, -1, P_MOD) % P_MOD
    T = t * INV2 % P_MOD
    XW = (-T * INV6) % P_MOD
    Tm1 = (T - 1) % P_MOD
    Y = (-s * T % P_MOD * Tm1 % P_MOD * INV6) % P_MOD
    # subtract W_D = (1/3, -1/6)
    WDx = INV3
    WDy = (-INV6) % P_MOD
    nWDy = (-WDy) % P_MOD
    res = w_add_affine(XW, Y, WDx, nWDy)
    if res is None:
        return None
    return (res[0], res[1], 1)


# ----------------------------------------------------------------
# Workload generation (identical structure to the C++ version)
# ----------------------------------------------------------------
def make_workload(N, seed):
    rng = random.Random(seed)
    G = (1, G_B, G_C)
    scalars, pk_pts, jac_pts = [], [], []
    retries = 0
    while len(pk_pts) < N:
        s = rng.randrange(1, Q_ORDER)
        Pi = pk_scalar_mul(s, G)
        Wi = phi_D_to_W(Pi)
        if Wi is None:
            retries += 1
            continue
        pk_pts.append(Pi)
        jac_pts.append(Wi)
        scalars.append(rng.randrange(1, Q_ORDER))
    return scalars, pk_pts, jac_pts, retries


# ----------------------------------------------------------------
# Workers and timing (checksum INSIDE the timed region for both seq and par)
# ----------------------------------------------------------------
def _checksum_pk(outs):
    h = 0
    for pt in outs:
        a, b, _ = pt
        h ^= (a ^ (b << 1))
    return h % P_MOD


def _checksum_jac(outs):
    h = 0
    for pt in outs:
        if pt is None:
            continue
        X, Y, _ = pt
        h ^= (X ^ (Y << 1))
    return h % P_MOD


def worker_pk(args):
    scalars, points = args
    outs = [pk_scalar_mul(k, P) for k, P in zip(scalars, points)]
    return _checksum_pk(outs)


def worker_jac(args):
    scalars, points = args
    outs = [jac_scalar_mul(k, P) for k, P in zip(scalars, points)]
    return _checksum_jac(outs)


def split_chunks(items, n):
    base, extra = divmod(len(items), n)
    start = 0
    for i in range(n):
        size = base + (1 if i < extra else 0)
        yield items[start:start + size]
        start += size


def time_sequential(kind, scalars, points, repeats):
    samples = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        if kind == "pk":
            outs = [pk_scalar_mul(k, P) for k, P in zip(scalars, points)]
            _ = _checksum_pk(outs)
        else:
            outs = [jac_scalar_mul(k, P) for k, P in zip(scalars, points)]
            _ = _checksum_jac(outs)
        t1 = time.perf_counter()
        samples.append(t1 - t0)
    samples.sort()
    return samples[len(samples) // 2]


def time_parallel(kind, scalars, points, T, repeats):
    runner = worker_pk if kind == "pk" else worker_jac
    chunks = list(zip(split_chunks(scalars, T), split_chunks(points, T)))
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=T) as pool:
        pool.map(runner, [([], []) for _ in range(T)])   # warm-up
        samples = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            _ = pool.map(runner, chunks)
            t1 = time.perf_counter()
            samples.append(t1 - t0)
    samples.sort()
    return samples[len(samples) // 2]


# ----------------------------------------------------------------
# Driver
# ----------------------------------------------------------------
def main():
    cpu_logical = os.cpu_count() or 1
    cpu_physical = max(1, cpu_logical // 2)     # rough heuristic; refine by hand

    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["pk", "jac", "both"], default="both")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--threads", type=int, nargs="+",
                    default=[1, 2, 4, min(8, cpu_logical)])
    ap.add_argument("--batch-sizes", type=int, nargs="+",
                    default=[1, 16, 64, 256, 1024, 4096])
    ap.add_argument("--out", default="results/exp4_smoke.json")
    args = ap.parse_args()

    print(f"[setup] python {sys.version.split()[0]}")
    print(f"[setup] logical cpus = {cpu_logical}, "
          f"physical estimate = {cpu_physical}")
    print(f"[setup] threads: {args.threads}")
    print(f"[setup] batch sizes: {args.batch_sizes}")

    payload = {
        "python_version": sys.version.split()[0],
        "cpu_logical": cpu_logical,
        "cpu_physical_estimate": cpu_physical,
        "seed": args.seed,
        "repeats": args.repeats,
        "threads": args.threads,
        "batch_sizes": args.batch_sizes,
        "rows": [],
    }

    kinds = ["pk", "jac"] if args.kind == "both" else [args.kind]
    for kind in kinds:
        for N in args.batch_sizes:
            scalars, pk_pts, jac_pts, retries = make_workload(N, args.seed)
            if retries:
                print(f"[{kind}] N={N}: {retries} retries in workload gen")
            points = pk_pts if kind == "pk" else jac_pts

            t_seq = time_sequential(kind, scalars, points, args.repeats)
            thr_seq = N / t_seq
            print(f"[{kind}] N={N:>6}  T= 1  "
                  f"T_med={t_seq*1e3:>9.2f} ms  thr={thr_seq:>9.1f} op/s")

            for T in args.threads:
                if T <= 1 or T > N:
                    continue
                t_par = time_parallel(kind, scalars, points, T, args.repeats)
                speedup = t_seq / t_par
                eff = speedup / T
                thr_par = N / t_par
                print(f"[{kind}] N={N:>6}  T={T:>2}  "
                      f"T_med={t_par*1e3:>9.2f} ms  "
                      f"speedup={speedup:>5.2f}x  "
                      f"eff={eff*100:>5.1f}%  "
                      f"thr={thr_par:>9.1f} op/s")
                payload["rows"].append({
                    "kind": kind, "N": N, "T": T,
                    "T_seq_s": t_seq, "T_par_s": t_par,
                    "speedup": speedup, "efficiency": eff,
                    "thr_seq_ops": thr_seq, "thr_par_ops": thr_par,
                })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\n[exp4] wrote {out}")


if __name__ == "__main__":
    main()