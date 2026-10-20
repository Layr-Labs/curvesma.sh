"""Organizer checks of analytic package data; never executes candidate content."""
from bisect import bisect_left
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("secp256k1-ecdlp", "bls12-381-invert-g1")


def finite_bsgs(r, generator, target, identity, combine, inverse):
    """Small organizer fixture for the published decomposition, not a production run."""
    m = 1 << (r.bit_length() // 2)
    h = (r + m - 1) // m
    table, value = [], identity
    for j in range(m):
        table.append((value, j))
        value = combine(value, generator)
    table.sort()
    step, gamma = inverse(value), target
    for i in range(h):
        index = bisect_left(table, (gamma, -1))
        if index < m and table[index][0] == gamma:
            return (i*m + table[index][1]) % r, i+1, m, h
        gamma = combine(gamma, step)
    raise AssertionError("Complete decomposition failed")


class BaselinePackages(unittest.TestCase):
    def test_exhaustive_finite_groups_including_last_partial_block(self):
        for r in (3, 5, 7, 11, 23, 101):
            for x in range(r):
                with self.subTest(r=r, x=x):
                    recovered, steps, m, h = finite_bsgs(
                        r, 1, x, 0, lambda a,b: (a+b) % r, lambda a: -a % r)
                    self.assertEqual(recovered, x)
                    self.assertLessEqual(steps, h)
                    self.assertGreaterEqual(m*h, r)
                    self.assertLess((h-1)*m, r)

    def test_multiplicative_reduction_to_point_scalar(self):
        # Generator 2 has order 11 modulo 23, mirroring the GT reduction.
        for x in range(11):
            recovered, _, _, _ = finite_bsgs(
                11, 2, pow(2, x, 23), 1,
                lambda a,b: a*b % 23, lambda a: pow(a, 10, 23))
            self.assertEqual(recovered, x)
            self.assertEqual(pow(2, recovered, 23), pow(2, x, 23))

if __name__ == "__main__":
    unittest.main()
