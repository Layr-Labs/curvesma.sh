"""Exact ledger tracer, independent of the certificate generator."""
import importlib.util
from pathlib import Path
import unittest


class CertificateTracer(unittest.TestCase):
    def test_exact_budget_and_coverage_meet_rank_gate(self):
        path = Path(__file__).with_name('verify.py')
        if not path.exists():
            self.fail('v6 certificate generator is not implemented')
        spec = importlib.util.spec_from_file_location('v6_certificate', path)
        assert spec is not None and spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        c = mod.certificate()
        self.assertTrue(all(c['checks'].values()))
        self.assertEqual(c['E'], 79650816)
        self.assertEqual(c['claim'], [669, 5])
        R, C = c['R'], c['C']
        self.assertEqual(R, sum(c['rows'].values()))
        self.assertLess(R ** 5, C ** 5 * 2 ** 669)
        self.assertLess((6 * R) ** 10, (5 * C) ** 10 * 2 ** 1337)
        m, q = c['m'], c['q']
        self.assertLess(m * m, q)
        self.assertGreaterEqual(200 * m * m, 81 * q)
        self.assertEqual(c['grid_records'], 2 * m)


if __name__ == '__main__':
    unittest.main()
