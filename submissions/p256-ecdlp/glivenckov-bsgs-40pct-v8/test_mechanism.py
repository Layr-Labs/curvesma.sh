"""Regression guard for the changed baby labels and exact coverage."""
import unittest
from mechanism import Field, affine_scalar, recover

class OddLiftCoverage(unittest.TestCase):
    def test_scalar_80_in_full_odd_lift_grid(self):
        # 8 odd babies and ceil(107/16)=7 giants must cover ALL 107 scalars.
        p, r, gen = 97, 107, (0, 43)
        got, count = recover(affine_scalar(80, gen, p), gen, p, r, 8, 7, Field(p))
        self.assertEqual(got, 80)
        self.assertLessEqual(count, 2)

if __name__ == '__main__':
    unittest.main()
