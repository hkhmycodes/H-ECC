"""F_p arithmetic using GMP via gmpy2.

Same interface as field.py.  Activated by setting EXP3_FIELD=gmp
before importing field / pk_arith / ws_arith.
"""

from __future__ import annotations
from typing import Dict

import gmpy2
from gmpy2 import mpz


_stats: Dict[str, int] = {'M': 0, 'S': 0, 'A': 0, 'C': 0, 'I': 0}
_tracking = False


def from_int(x) -> mpz:
    return mpz(x)


def start_tracking() -> None:
    global _tracking
    _tracking = True
    for k in _stats:
        _stats[k] = 0


def stop_tracking() -> Dict[str, int]:
    global _tracking
    _tracking = False
    return dict(_stats)


def mul(a, b, p):
    if _tracking:
        _stats['M'] += 1
    return (a * b) % p


def sqr(a, p):
    if _tracking:
        _stats['S'] += 1
    return (a * a) % p


def add(a, b, p):
    if _tracking:
        _stats['A'] += 1
    return (a + b) % p


def sub(a, b, p):
    if _tracking:
        _stats['A'] += 1
    return (a - b) % p


def cmul(a, k, p):
    if _tracking:
        _stats['C'] += 1
    return (a * k) % p


def smul(a, k, p):
    if _tracking:
        _stats['A'] += 1
    return (a * k) % p


def inv(a, p):
    if _tracking:
        _stats['I'] += 1
    return gmpy2.invert(a, p)