"""Independent small exact tests; run in a credential-free Python container."""
import importlib.util
from pathlib import Path
import unittest


class DenseTracer(unittest.TestCase):
    def load(self):
        path = Path(__file__).with_name('toy.py')
        if not path.exists():
            self.fail('dense supplied-power lift is not implemented')
        spec = importlib.util.spec_from_file_location('dense_toy', path)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_dense_lift_recovers_all_binary_root_choices(self):
        mod = self.load()
        r, g, d = 97, 5, 16
        q = (r - 1) // d
        mul = lambda k, point: k * point % r
        neg = lambda point: -point % r
        for k1 in range(d):
            for k2 in range(q):
                secret = pow(g, q * k1 + k2, r)
                powers = [pow(secret, i, r) for i in range(1, d + 1)]
                recovered, bits = mod.dense_lift(1, secret, powers, k2, g, d, r, mul, neg)
                self.assertEqual(recovered, secret)
                self.assertEqual(bits, [(k1 >> t) & 1 for t in range(4)])


    def test_scalar_grid_has_exact_partial_interval_coverage(self):
        mod = self.load()
        self.assertTrue(hasattr(mod, 'scalar_grid'), 'scalar-action grid is not implemented')
        r, g, d, m = 97, 5, 16, 2
        q, b = (r - 1) // d, pow(g, d, r)
        for k2 in range(q):
            k = mod.scalar_grid(1, pow(b, k2, r), b, r, m,
                                lambda c, point: c * point % r,
                                lambda point: point.to_bytes(1, 'big'))
            self.assertEqual(k, k2 if k2 < m * m else None)


    def test_real_curve_pipeline_preserves_exact_coverage(self):
        mod = self.load()
        self.assertTrue(hasattr(mod, 'run_toy'), 'real-curve dense pipeline is not implemented')
        report = mod.run_toy(40, 20261023)
        self.assertEqual(report['trials'], 40)
        self.assertEqual(report['wrong_outputs'], 0)
        self.assertEqual(report['coverage_prediction_errors'], 0)
        self.assertEqual(sum(report['k2_histogram']), 40)
        self.assertTrue(report['full_xy_negation_distinct'])


if __name__ == '__main__':
    unittest.main()
