"""Regression tests for real-EC orbit metadata and exact exponent recovery."""
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent


def load_toy(test):
    path = HERE / 'toy.py'
    test.assertTrue(path.is_file(), 'six-class mechanism is not implemented')
    spec = importlib.util.spec_from_file_location('sixclass_toy', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SixClassTests(unittest.TestCase):
    def test_nonzero_orbits_have_six_distinct_points_and_exact_digits(self):
        toy = load_toy(self)
        curve = toy.Curve()
        points = [None]
        for _ in range(1, curve.r):
            points.append(curve.add(points[-1], curve.P))
        self.assertIsNone(curve.add(points[-1], curve.P))
        self.assertEqual(len(set(points)), curve.r)
        self.assertEqual(curve.mul(curve.lam, curve.P), curve.tau(curve.P))
        for k, point in enumerate(points[1:], 1):
            orbit = curve.orbit(point)
            self.assertEqual(len(set(orbit)), 6)
            self.assertEqual(points[curve.lam * k % curve.r], curve.tau(point))
            keys = []
            for transformed in orbit:
                rec = curve.canonical_record(transformed, 1, 123)
                self.assertEqual(len(rec), 128)
                self.assertEqual(int.from_bytes(rec[104:120], 'big'), 123)
                digit = int.from_bytes(rec[120:], 'big') >> 1
                self.assertIn(digit, range(6))
                x, y = int.from_bytes(rec[:48], 'big'), int.from_bytes(rec[48:96], 'big')
                self.assertEqual(curve.orbit(transformed)[digit], (x, y))
                keys.append(rec[:104])
            self.assertEqual(len(set(keys)), 1)
        self.assertEqual(curve.canonical_record(None, 0, 0)[96:104], (1).to_bytes(8, 'big'))

    def test_grid_recovers_exact_exponent_on_each_of_six_lifts(self):
        toy = load_toy(self)
        self.assertTrue(hasattr(toy, 'scalar_grid'), 'six-class grid is not implemented')
        curve, g, d, a, b = toy.Curve(), 11, 64, 3, 5
        q = (curve.r - 1) // d
        zeta, Q = pow(g, d, curve.r), q // 6
        hits = 0
        for k2 in range(q):
            target = curve.mul(pow(zeta, k2, curve.r), curve.P)
            result = toy.scalar_grid(curve, target, zeta, q, a, b)
            self.assertEqual(result['actions'], a + b)
            self.assertEqual(result['duplicate_keys'], 0)
            self.assertEqual(result['k2'] is not None, k2 % Q < a * b)
            if result['k2'] is not None:
                self.assertEqual(result['k2'], k2)
                hits += 1
        self.assertEqual(hits, 6 * a * b)
        with self.assertRaises(AssertionError):
            toy.scalar_grid(curve, curve.P, zeta, q, 4, 8)

    def test_literal_dense_powers_recover_all_root_bits_after_sixclass_match(self):
        toy = load_toy(self)
        self.assertTrue(hasattr(toy, 'dense_lift'), 'dense supplied-power lifting is absent')
        curve, g, d, k2 = toy.Curve(), 11, 64, 190
        q = (curve.r - 1) // d
        for k1 in range(d):
            secret = pow(g, q * k1 + k2, curve.r)
            powers = [curve.mul(pow(secret, e, curve.r), curve.P) for e in range(1, d + 1)]
            recovered, bits = toy.dense_lift(curve, powers, k2, g, d)
            self.assertEqual(recovered, secret)
            self.assertEqual(sum(bit << t for t, bit in enumerate(bits)), k1)
        with self.assertRaises(ValueError):
            toy.dense_lift(curve, powers, k2 - 1, g, d)

    def test_captured_matrix_checks_all_digits_orbits_and_coverage(self):
        toy = load_toy(self)
        self.assertTrue(hasattr(toy, 'run_matrix'), 'captured six-class matrix is absent')
        result = toy.run_matrix(trials=5, seed=71, full=False)
        audit = result['orbit_audit']
        self.assertEqual(audit['nonzero_points'], 12288)
        self.assertEqual(audit['orbit_size_exceptions'], 0)
        self.assertEqual(audit['digit_pair_tests'], 2048 * 36)
        self.assertEqual(audit['wrong_digit_relations'], 0)
        self.assertEqual(audit['fixed_nonidentity_counts'], [0] * 5)
        self.assertGreaterEqual(len(result['configurations']), 7)
        for row in result['configurations']:
            self.assertEqual(row['exhaustive_k2_hits'], 6 * row['a'] * row['b'])
            self.assertEqual(row['wrong_outputs'], 0)
            self.assertEqual(row['coverage_prediction_errors'], 0)
            self.assertEqual(row['random_trials'], 5)


if __name__ == '__main__':
    unittest.main()
