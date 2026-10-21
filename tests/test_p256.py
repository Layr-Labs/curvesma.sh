"""Independent P-256 group arithmetic; no candidate content is executed."""
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PROFILE = json.loads((ROOT / 'tracks/p256-ecdlp/target-profile.json').read_text())
P = int(PROFILE['field_modulus'], 16)
A = int(PROFILE['coefficient_a'], 16)
B = int(PROFILE['coefficient_b'], 16)
R = int(PROFILE['subgroup_order'], 16)
G = tuple(int(PROFILE['generator_P'][c], 16) for c in ('x', 'y'))


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u and (y + v) % P == 0:
        return None
    slope = ((3*x*x + A) * pow(2*y, -1, P) if left == right
             else (v-y) * pow(u-x, -1, P)) % P
    z = (slope*slope-x-u) % P
    return z, (slope*(x-z)-y) % P


def multiply(point, scalar):
    result = None
    while scalar:
        if scalar & 1:
            result = add(result, point)
        point = add(point, point)
        scalar >>= 1
    return result


class P256Profile(unittest.TestCase):
    def test_sec2_parameters_and_group_order(self):
        self.assertEqual(P, 2**256 - 2**224 + 2**192 + 2**96 - 1)
        self.assertEqual(A, P-3)
        self.assertEqual(PROFILE['cofactor'], 1)
        self.assertEqual((G[1]**2 - G[0]**3 - A*G[0] - B) % P, 0)
        self.assertIsNone(multiply(G, R))
        self.assertEqual(multiply(G, R-1), (G[0], (-G[1]) % P))
        for key in ('field_modulus', 'subgroup_order', 'coefficient_a', 'coefficient_b'):
            self.assertEqual(len(PROFILE[key]), 64)

    def test_doubling_vector_and_scalar_edges(self):
        doubled = (
            int('7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978', 16),
            int('07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1', 16),
        )
        self.assertEqual(multiply(G, 2), doubled)
        self.assertIsNone(multiply(G, 0))
        self.assertEqual(multiply(G, 1), G)
        self.assertNotEqual(multiply(G, 7), multiply(G, 8))
        self.assertEqual(set(PROFILE['input']), {'Q'})
        self.assertEqual(set(PROFILE['output']), {'x'})
        self.assertNotIn('auxiliary_inputs', PROFILE)
        self.assertNotEqual((G[1]**2 - G[0]**3 - 7) % P, 0)


if __name__ == '__main__':
    unittest.main()
