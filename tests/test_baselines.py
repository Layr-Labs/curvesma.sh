"""Organizer checks of analytic package data; never executes candidate content."""
from bisect import bisect_left
from decimal import Decimal, localcontext
from fractions import Fraction
import re
import json
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

    def test_exact_accounting_and_logarithmic_bounds(self):
        from test_profiles import PROFILES
        profiles = PROFILES
        for target in profiles["targets"]:
            name = target["id"]
            with self.subTest(name=name):
                directory = ROOT / "tracks" / name / "submission"
                manifest = json.loads((directory / "submission.json").read_text())
                declared = manifest['claims'][0]
                value = declared['score']
                score = Fraction(int(value['numerator']), int(value['denominator']))
                statement = declared['statement']
                claim = {'claim': {'time_log2': str(Decimal(score.numerator)/score.denominator),
                    'memory_log2_bytes': re.search(r'Peak memory is at most 2\^([0-9.]+) bytes', statement)[1],
                    'preprocessing_log2': re.search(r'preprocessing is at most 2\^([0-9.]+) target', statement)[1]}}

                ledger = json.loads((ROOT / "tests/fixtures/baseline-accounting-v1.json").read_text())[name]
                cost = json.loads((ROOT / "tracks" / name / "rules.md").read_text().split("```json\n")[1].split("```", 1)[0])["accounting"]
                C = cost["reference_operation_costs"][name]
                work, parameters = ledger["work"], ledger["parameters"]
                total = int(work["total"])
                self.assertEqual(total, sum(int(work[k]) for k in
                    ("setup", "baby_table", "merge_sort", "giant_search", "verification_and_output")))
                self.assertEqual(int(work["preprocessing"]), sum(int(work[k]) for k in
                    ("setup", "baby_table", "merge_sort")))
                self.assertIn(str(total), (directory / "proof.md").read_text())
                self.assertNotIn("resource_ledger", claim)
                m, h = int(parameters["baby_steps"]), int(parameters["giant_steps"])
                r = int(target["subgroup_order"], 16)
                self.assertGreaterEqual(m*h, r)
                self.assertLess((h-1)*m, r)
                self.assertEqual(m, 1 << parameters["split_bits"])
                memory = int(ledger["peak_memory_bytes"])
                self.assertEqual(memory, 2*m*parameters["record_words"]*32+(1<<20))
                self.assertLess(memory//32, 1<<256)
                for key, bound in (("time_log2", Decimal(total)/C), ("memory_log2_bytes", memory),
                                   ("preprocessing_log2", Decimal(work["preprocessing"])/C)):
                    with localcontext() as ctx:
                        ctx.prec = 100
                        true_log = Decimal(bound).ln()/Decimal(2).ln()
                        reported = Decimal(str(claim["claim"][key]))
                        self.assertGreaterEqual(reported, true_log)
                        self.assertLess(reported-true_log, Decimal('0.0000011'))



if __name__ == "__main__":
    unittest.main()
