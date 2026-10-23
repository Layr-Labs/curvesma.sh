"""Exact production ledger and fail-closed certificate regressions."""
import unittest
from fractions import Fraction as F
from pathlib import Path


class CertificateTests(unittest.TestCase):
    def test_exact_sign_quotient_budget(self):
        self.assertTrue((Path(__file__).parent / 'verify.py').exists(), 'exact verifier is missing')
        import verify
        cert = verify.certificate()
        self.assertEqual(cert['R'], 19918011671547009368992591426495508435762944)
        self.assertEqual(cert['C'], 2000)
        self.assertEqual(cert['baby'], 98748750602204457565883422358598783)
        self.assertEqual(cert['giant'], 98748750602204457565883422358598784)
        self.assertEqual(cert['claim'], {'numerator': '3331', 'denominator': '25'})
        self.assertTrue(all(cert['checks'].values()))
        self.assertEqual(sum(cert['rows'].values()), cert['R'])
        stress = F(6 * cert['R'], 5 * cert['C'])
        self.assertTrue(verify.under_power(stress, F(13314, 100)))
        self.assertFalse(verify.under_power(stress, F(13313, 100)))
        lo, hi = verify.log2_bounds(F(cert['R'], cert['C']))
        self.assertLess(lo, hi)
        self.assertLess(hi - lo, F(1, 10 ** 25))
        with self.assertRaises(AssertionError):
            verify.certificate(extra_per_record=0)
        with self.assertRaises(AssertionError):
            verify.certificate(claim=F(13323, 100))


if __name__ == '__main__':
    unittest.main()
