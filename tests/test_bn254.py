"""Organizer checks of BN254 GT parameters, encoding and reduction equations."""
import hashlib
import json
from pathlib import Path
import unittest

import py_ecc
from py_ecc.optimized_bn128 import FQ12, G1, G2, curve_order, field_modulus, multiply, pairing

ROOT = Path(__file__).resolve().parents[1]
NAME = 'bn254-gt-dlp'
PROFILE = json.loads((ROOT / 'tracks' / NAME / 'target-profile.json').read_text())


def decode_gt(value):
    if not isinstance(value, list) or len(value) != 12:
        raise ValueError('coefficient count')
    for c in value:
        if not isinstance(c, str) or len(c) != 64 or any(v not in '0123456789abcdef' for v in c):
            raise ValueError('noncanonical coefficient')
        if int(c, 16) >= field_modulus:
            raise ValueError('coefficient out of range')
    z = FQ12([int(c, 16) for c in value])
    if z == FQ12.zero() or z ** curve_order != FQ12.one():
        raise ValueError('outside GT')
    return z


class BN254TargetGroup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g = pairing(G2, G1)

    def test_pinned_sources_parameters_and_generator(self):
        self.assertEqual(int(PROFILE['field_modulus'], 16), field_modulus)
        self.assertEqual(int(PROFILE['subgroup_order'], 16), curve_order)
        u = int(PROFILE['parameter_x'], 16)
        self.assertEqual(36*u**4+36*u**3+24*u**2+6*u+1, field_modulus)
        self.assertEqual(36*u**4+36*u**3+18*u**2+6*u+1, curve_order)
        for path, expected in PROFILE['pairing']['reference']['source_sha256'].items():
            source = Path(py_ecc.__file__).parent / path.removeprefix('py_ecc/')
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), expected)
        self.assertEqual(decode_gt(PROFILE['generator_g']['coefficients']), self.g)
        self.assertNotEqual(self.g, FQ12.one())
        self.assertEqual(self.g ** curve_order, FQ12.one())
        u_field = FQ12([-9, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
        self.assertEqual(u_field*u_field, -FQ12.one())

    def test_dlp_equation_and_identity(self):
        self.assertEqual(set(PROFILE['input']), {'z'})
        self.assertEqual(set(PROFILE['output']), {'x'})
        self.assertNotIn('auxiliary_inputs', PROFILE)
        for x in (0, 1, 7, curve_order - 1):
            z = self.g ** x
            encoded = [format(int(c), '064x') for c in z.coeffs]
            self.assertEqual(decode_gt(encoded), z)
            self.assertNotEqual(self.g ** ((x+1) % curve_order), z)
            self.assertEqual(z == FQ12.one(), x == 0)

    def test_reject_noncanonical_and_non_subgroup_inputs(self):
        valid = PROFILE['generator_g']['coefficients']
        cases = [valid[:-1], valid+[valid[0]], ['0']+valid[1:],
                 [format(field_modulus, '064x')]+valid[1:],
                 ['A'*64]+valid[1:], ['0'*64]*12,
                 [format(2, '064x')]+['0'*64]*11]
        for value in cases:
            with self.subTest(value=value), self.assertRaises(ValueError):
                decode_gt(value)

    def test_source_group_and_fixed_argument_reductions(self):
        a, b = 7, 11
        self.assertEqual(pairing(G2, multiply(G1, a)), self.g ** a)
        # Fixed-generator logs log_g(z)=a and log_g(e(P,Q0))=b
        # recover the pairing preimage scalar a/b modulo r.
        scalar = a * pow(b, -1, curve_order) % curve_order
        self.assertEqual(pairing(multiply(G2, scalar), multiply(G1, b)), self.g ** a)
        self.assertEqual(pairing(multiply(G2, b), multiply(G1, scalar)), self.g ** a)

    def test_single_track_and_borrowed_normalization(self):
        manifest = json.loads((ROOT / 'benchmark.json').read_text())
        self.assertEqual([t['name'] for t in manifest['tracks'] if t['name'].startswith('bn254-')], [NAME])
        costs = json.loads((ROOT / 'tracks/classical-group-work-v2.json').read_text())
        self.assertEqual(costs['reference_operation_costs'][NAME], 20000)
        self.assertIn('Borrowed', costs['reference_operations'][NAME])
