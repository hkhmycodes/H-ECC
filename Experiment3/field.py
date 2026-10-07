"""F_p arithmetic primitives with optional operation counting.

Two backends are supported:
  - field.py       : pure Python ints (default)
  - field_gmp.py   : gmpy2.mpz (selected by setting EXP3_FIELD=gmp before
                     importing field / pk_arith / ws_arith)

Both backends expose the same function names.
"""

from __future__ import annotations
from typing import Dict


_stats: Dict[str, int] = {'M': 0, 'S': 0, 'A': 0, 'C': 0, 'I': 0}
_tracking = False


def from_int(x: int) -> int:
    return int(x)


def start_tracking() -> None:
    global _tracking
    _tracking = True
    for k in _stats:
        _stats[k] = 0


def stop_tracking() -> Dict[str, int]:
    global _tracking
    _tracking = False
    return dict(_stats)


def mul(a: int, b: int, p: int) -> int:
    if _tracking:
        _stats['M'] += 1
    return (a * b) % p


def sqr(a: int, p: int) -> int:
    if _tracking:
        _stats['S'] += 1
    return (a * a) % p


def add(a: int, b: int, p: int) -> int:
    if _tracking:
        _stats['A'] += 1
    return (a + b) % p


def sub(a: int, b: int, p: int) -> int:
    if _tracking:
        _stats['A'] += 1
    return (a - b) % p


def cmul(a: int, k: int, p: int) -> int:
    if _tracking:
        _stats['C'] += 1
    return (a * k) % p


def smul(a: int, k: int, p: int) -> int:
    if _tracking:
        _stats['A'] += 1
    return (a * k) % p


def inv(a: int, p: int) -> int:
    if _tracking:
        _stats['I'] += 1
    return pow(a, p - 2, p)