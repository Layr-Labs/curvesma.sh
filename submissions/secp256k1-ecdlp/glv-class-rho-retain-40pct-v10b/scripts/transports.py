"""Exact per-scale delta_transport values (Fractions) for the v10 envelope.

The v10b ledger uses transport = 64*DP^2/N with fixed production DP=2^64
and class count N=(r-1)/6. The same selected-pair transport formula at actual
toy spacing s=2^d is F(64*(2^d)^2,N_toy). This is NOT an eviction bound.
The default large toy subgroup order is Q=523969; each scale supplies N.
"""
from fractions import Fraction as F

Q = 523969
N_toy = (Q - 1) // 6

def delta_transport(d, classes=N_toy):
    """Exact same production formula, parameterized by actual toy class count."""
    assert isinstance(d,int) and d >= 0 and classes > 0
    return F(64 * (2 ** d) ** 2, classes)

if __name__ == "__main__":
    for d in (2, 3, 4, 5):
        v = delta_transport(d)
        print(d, v.numerator / v.denominator)
