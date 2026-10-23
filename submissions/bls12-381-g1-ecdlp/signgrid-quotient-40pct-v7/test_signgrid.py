"""Regression gates for deterministic orientation-aware quotient recovery."""
import importlib.util
import unittest
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


class SignGridTests(unittest.TestCase):
    def test_exhaustive_real_ec_sign_coverage(self):
        spec = importlib.util.spec_from_file_location('toy', HERE / 'toy.py')
        self.assertIsNotNone(spec)
        self.assertTrue((HERE / 'toy.py').exists(), 'orientation-aware real-EC toy is missing')
        toy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(toy)
        curve = toy.Curve()
        P, r, g, d, a, b = curve.P, curve.r, 11, 64, 6, 7
        q, zeta = (r - 1) // d, pow(g, d, r)
        hits = 0
        for k in range(q):
            A = curve.mul(pow(zeta, k, r), P)
            result = toy.scalar_grid(curve, A, zeta, q, a, b)
            self.assertEqual(result['actions'], a + b)
            self.assertEqual(result['duplicate_keys'], 0)
            self.assertEqual(result['k2'] is not None, k % (q // 2) < a * b)
            if result['k2'] is not None:
                self.assertEqual(result['k2'], k)
                hits += 1
        self.assertEqual(hits, 84)


    def test_supplied_power_lift(self):
        import toy
        self.assertTrue(callable(getattr(toy, 'dense_lift', None)), 'supplied-power lift is missing')
        curve, g, d = toy.Curve(), 11, 64
        q = (curve.r - 1) // d
        for k1 in range(d):
            for k2 in [0, 41, 96, 137]:
                s = pow(g, q * k1 + k2, curve.r)
                powers = [curve.mul(pow(s, e, curve.r), curve.P) for e in range(1, d + 1)]
                found, bits = toy.dense_lift(curve, powers, k2, g, d)
                self.assertEqual(found, s)
                self.assertEqual(len(bits), 6)
        C = curve.mul(11, curve.P)
        self.assertNotEqual(C, curve.neg(C))
        self.assertEqual(curve.mul(11 + curve.r, curve.P), C)


    def test_multi_configuration_evidence(self):
        import toy
        self.assertTrue(callable(getattr(toy, 'run_matrix', None)), 'multi-configuration runner is missing')
        evidence = toy.run_matrix(5, 20261023, full=False)
        self.assertGreaterEqual(len(evidence['configurations']), 7)
        for row in evidence['configurations']:
            self.assertEqual(row['random_trials'], 5)
            self.assertEqual(row['wrong_outputs'], 0)
            self.assertEqual(row['coverage_prediction_errors'], 0)
            self.assertEqual(row['exhaustive_k2_hits'], 2 * row['a'] * row['b'])
            self.assertEqual(row['exhaustive_k2_trials'], row['q'])
            self.assertEqual(sum(row['k2_histogram']), 5)
            self.assertEqual(sum(row['hit_k2_histogram']), row['random_hits'])


if __name__ == '__main__':
    unittest.main()
