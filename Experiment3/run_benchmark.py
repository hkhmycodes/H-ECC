#!/usr/bin/env python3
"""Orchestrator for EXP-3 — CPU arithmetic cost benchmark."""

from __future__ import annotations
import argparse
import os
import sys


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", choices=["py", "gmp"], default="py",
                    help="field backend: py (pure Python) or gmp (gmpy2)")
    ap.add_argument("--seeds", default="0",
                    help="comma-separated list of seeds, e.g. 0,1,2,3,4")
    ap.add_argument("--N",              type=int, default=1000)
    ap.add_argument("--N-scalar",       type=int, default=50)
    ap.add_argument("--repeats",        type=int, default=5)
    ap.add_argument("--scalar-repeats", type=int, default=5)
    ap.add_argument("--outdir",         type=str, default="results")
    return ap.parse_args()


def main():
    args = parse_args()
    os.environ["EXP3_FIELD"] = args.field

    # Import bench *after* setting EXP3_FIELD, since bench.py reads it at
    # import time to select the field backend.
    import bench

    seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
    bench.run(seeds=seeds,
              N=args.N, N_scalar=args.N_scalar,
              repeats=args.repeats, scalar_repeats=args.scalar_repeats,
              outdir=args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())