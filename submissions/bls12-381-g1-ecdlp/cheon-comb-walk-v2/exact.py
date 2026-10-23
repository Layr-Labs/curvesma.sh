"""Exact accounting helpers for curvesma.sh claims (stdlib only).

Work is expressed in target-group-operation equivalents ("units"): one reference
operation is 1 unit and every other scalar AArch64 instruction is 1/C units.
"""
import math
from fractions import Fraction


def iroot_ceil(n, k):
    """Smallest integer x with x^k >= n (integer Newton iteration)."""
    x = 1 << (n.bit_length() // k + 1)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            break
        x = y
    while x ** k < n:
        x += 1
    return x


def pow2_ceil(exponent, denominator=100):
    """Smallest integer >= 2^(exponent), exponent given as a decimal (e.g. 96.75).

    Uses 2^(a/d) >= v  <=>  2^a >= v^d, exact on integers.
    """
    e = Fraction(exponent).limit_denominator(denominator)
    return iroot_ceil(2 ** e.numerator, e.denominator) if e.denominator > 1 else 2 ** e.numerator


def log2(v):
    v = Fraction(v)
    return math.log2(v.numerator) - math.log2(v.denominator)


def log2_up(v, places=6):
    """log2 of an exact positive rational, rounded up at the given decimal place.

    math.log2 on Python integers is accurate to ~1e-15 relative; the 1e-4 ulp
    guard keeps the rounding upward.
    """
    approx = log2(v)
    scale = 10 ** places
    q = math.ceil(approx * scale + 1e-4)
    assert q / scale >= approx
    return Fraction(q, scale)


def units(instructions, C):
    """Convert a scalar instruction count to units at reference price C."""
    return Fraction(instructions, C)


class Ledger:
    """Named cost components in units; reports exact totals and rounded-up log2."""

    def __init__(self, C):
        self.C = C
        self.components = {}
        self.preprocessing = set()

    def add(self, name, value_units, preprocessing=False):
        self.components[name] = Fraction(value_units)
        if preprocessing:
            self.preprocessing.add(name)
        return self

    def add_instructions(self, name, instructions, preprocessing=False):
        return self.add(name, units(instructions, self.C), preprocessing)

    def total(self):
        return sum(self.components.values())

    def preprocessing_total(self):
        return sum((self.components[k] for k in self.preprocessing), Fraction(0))

    def report(self):
        total = self.total()
        pre = self.preprocessing_total()
        score = log2_up(total)
        return {
            "components_log2_units": {k: round(log2(v), 6) for k, v in self.components.items()},
            "time_log2": str(score),
            "time_log2_decimal": float(score),
            "time_log2_float": log2(total),
            "score_rational": {"numerator": str(score.numerator), "denominator": str(score.denominator)},
            "preprocessing_log2": str(log2_up(pre)) if pre > 0 else "0",
        }


def birthday_success(T, N, eps=Fraction(0), delta=Fraction(0)):
    """Random-mapping collision bound: 1 - exp(-T^2/(2N) * (1 - eps)) - delta.

    T: total fresh steps, N: number of classes, eps: relative exponent loss
    (walk penalty, fruitless steps), delta: additive failure terms.
    """
    tau = Fraction(T * T, 2 * N)
    p = 1 - math.exp(-float(tau * (1 - Fraction(eps))))
    return {"tau": float(tau), "eps": float(eps), "model": p, "lower_bound": p - float(delta)}


def tau_for_success(p):
    """tau = T^2/(2N) at which 1 - exp(-tau) = p."""
    return -math.log(1 - p)
