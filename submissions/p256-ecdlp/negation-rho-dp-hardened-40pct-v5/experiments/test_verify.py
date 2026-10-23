#!/usr/bin/env python3
"""Regression checks for the inherited coefficient and certificate defects."""
import json
from pathlib import Path
import runpy
import unittest

HERE = Path(__file__).resolve().parent

class Coefficients(unittest.TestCase):
    def test_full_width_start_has_bounded_five_limb_update(self):
        module = runpy.run_path(str(HERE / 'verify.py'))
        update = module.get('update_coeff')
        self.assertTrue(callable(update), 'full-width five-limb update is missing')
        n = int('ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551', 16)
        for a in (0, 1, (1 << 255), n - 1):
            for c in (0, 1, n - 1):
                for sign in (1, -1):
                    got, temporary, subtractions = update(a, c, sign, n)
                    self.assertEqual(got, (sign * (2 * a + c)) % n)
                    self.assertLess(temporary, 1 << 258)
                    self.assertLessEqual(subtractions, 2)
                    self.assertLess(got, n)

    def test_shipped_certificate_has_no_failed_numeric_bounds(self):
        cert = json.loads((HERE.parent / 'certificate.json').read_text())
        self.assertTrue(cert.get('all_checks_pass'), 'inherited certificate contains FAIL bounds')
        self.assertTrue(cert.get('checks'))
        self.assertTrue(all(v == 'PASS' for v in cert['checks'].values()))

if __name__ == '__main__':
    unittest.main()
