"""Exact certificate qualification, including negative pricing and claim gates."""
import importlib.util
from fractions import Fraction as F
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent


class CertificateTests(unittest.TestCase):
    def test_sixclass_price_and_adjacent_hundredth_gate(self):
        path = HERE / 'verify.py'
        self.assertTrue(path.is_file(), 'six-class certificate is not implemented')
        spec = importlib.util.spec_from_file_location('sixclass_verify', path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        cert = module.certificate()
        self.assertTrue(all(cert['checks'].values()))
        self.assertEqual(cert['claim'], {'numerator': '2649', 'denominator': '20'})
        self.assertEqual(cert['orbit_extra_per_record'], 8192)
        self.assertEqual(cert['R'], 11500574304558579830758679329268362223869952)
        self.assertEqual(cert['policy_floor_interval'], ['132.441864108206', '132.441864108207'])
        with self.assertRaises(AssertionError):
            module.certificate(extra_per_record=256)
        with self.assertRaises(AssertionError):
            module.certificate(claim=F(3311, 25))  # 132.44 fails mandatory stress headroom.

    def test_package_binding_rejects_absent_documentation(self):
        path = HERE / 'verify.py'
        spec = importlib.util.spec_from_file_location('sixclass_verify', path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(hasattr(module, 'package_checks'), 'final package binding is absent')
        # The full package is checked only after its manifest and documentation exist.
        if not (HERE / 'submission.json').exists():
            with self.assertRaises(FileNotFoundError):
                module.package_checks(module.certificate())
        else:
            cert = module.certificate()
            module.package_checks(cert)
            self.assertTrue(all(cert['checks'].values()))


if __name__ == '__main__':
    unittest.main()
