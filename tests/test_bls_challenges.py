"""Checks of the published equations against the pinned arithmetic reference."""
import hashlib
import json
from pathlib import Path
import unittest

import py_ecc
from py_ecc.optimized_bls12_381 import (
    FQ, FQ2, FQ12, G1, G2, Z1, Z2, b, curve_order, eq,
    field_modulus, final_exponentiate, is_inf, is_on_curve,
    multiply, normalize, pairing,
)

ROOT = Path(__file__).resolve().parents[1]
PAIR = 'bls12-381-invert-g1'
DLOG = 'bls12-381-g1-ecdlp'
MILLER = 'bls12-381-miller-inversion'


def profile(name):
    return json.loads((ROOT / 'tracks' / name / 'target-profile.json').read_text())


def miller(point):
    if is_inf(point):
        return FQ12.one()
    x, y = normalize(point)
    return pairing(G2, (x, y, FQ.one()), final_exponentiate=False)


class BLSChallenges(unittest.TestCase):
    def test_reference_sources_and_curve_parameters(self):
        for name in (PAIR, DLOG, MILLER):
            p = profile(name)
            self.assertEqual(int(p['field_modulus'], 16), field_modulus)
            self.assertEqual(int(p['subgroup_order'], 16), curve_order)
            self.assertEqual(tuple(int(p['generator_P0'][c], 16) for c in ('x', 'y')),
                             tuple(int(c) for c in normalize(G1)))
            self.assertEqual(tuple(tuple(int(v, 16) for v in p['generator_Q0'][c])
                                   for c in ('x', 'y')),
                             tuple(tuple(c.coeffs) for c in normalize(G2)))
            for path, digest in p['pairing']['reference']['source_sha256'].items():
                source = Path(py_ecc.__file__).parent / path.removeprefix('py_ecc/')
                self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), digest)

    def test_auxiliary_powers_and_cheon_precondition(self):
        s, n = 11, 4
        powers = [multiply(G1, pow(s, i, curve_order)) for i in range(1, n + 1)]
        self.assertTrue(eq(powers[0], multiply(G1, s)))
        for left, right in zip(powers, powers[1:]):
            self.assertTrue(eq(multiply(left, s), right))
        for name in (PAIR, DLOG, MILLER):
            p = profile(name)
            self.assertEqual(p['auxiliary_inputs']['scoring_n'], 1 << 20)
            self.assertEqual(set(p['input']) - {'n', 'A'},
                             {'T'} if name == DLOG else {'z'} if name == PAIR else {'y'})
        self.assertEqual((curve_order - 1) % (1 << 20), 0)
        self.assertEqual(set(profile(DLOG)['output']), {'s'})
        self.assertFalse(eq(multiply(G1, s + 1), powers[0]))

    def test_general_pairing_accepts_a_nonfixed_g2_argument(self):
        a, q = 7, 3
        z = pairing(G2, multiply(G1, a))
        P = multiply(G1, a * pow(q, -1, curve_order) % curve_order)
        Q = multiply(G2, q)
        self.assertFalse(eq(Q, G2))
        self.assertEqual(pairing(Q, P), z)
        self.assertNotEqual(pairing(Q, multiply(P, 2)), z)
        self.assertEqual(set(profile(PAIR)['output']), {'P', 'Q'})
        self.assertEqual(pairing(G2, Z1), FQ12.one())
        self.assertEqual(pairing(Z2, G1), FQ12.one())

    def test_miller_requires_raw_equality_not_only_final_exponentiation(self):
        P = multiply(G1, 7)
        y = miller(P)
        self.assertNotEqual(y, FQ12.zero())
        self.assertNotEqual(y ** curve_order, FQ12.one())
        self.assertEqual(final_exponentiate(y), pairing(G2, P))
        self.assertNotEqual(miller(multiply(P, 2)), y)
        # Base-field factors disappear under final exponentiation. Such a
        # changed value must not pass an exact raw-Miller equality check.
        changed = y * FQ12([2] + [0] * 11)
        self.assertNotEqual(changed, y)
        self.assertEqual(final_exponentiate(changed), final_exponentiate(y))
        self.assertEqual(miller(Z1), FQ12.one())
        self.assertEqual(set(profile(MILLER)['output']), {'P'})

    def test_subgroup_requirement_is_material(self):
        point = (FQ(0), FQ(2), FQ.one())
        self.assertTrue(is_on_curve(point, b))
        self.assertFalse(is_inf(multiply(point, curve_order)))
        self.assertTrue(is_inf(multiply(G1, curve_order)))
        self.assertTrue(is_inf(multiply(G2, curve_order)))

    def test_yukon_tracks_keep_the_existing_interface(self):
        manifest = json.loads((ROOT / 'benchmark.json').read_text())
        self.assertEqual({t['name'] for t in manifest['tracks']},
                         {'secp256k1-ecdlp', PAIR, DLOG, MILLER})
        for track in manifest['tracks']:
            name = track['name']
            c = json.loads((ROOT / 'tracks' / name / 'verifier.json').read_text())
            self.assertEqual(track['objectives'], [{'name': 'time_log2', 'direction': '-'}])
            self.assertEqual(c['objectives'][0]['name'], 'time_log2')
            self.assertEqual(c['objectives'][0]['score_source'], {'claim': 0})
            self.assertEqual(track['paretoPolicies'], c['policies'])
            self.assertEqual(c['acceptance_policy'], 'exploratory-v1')
            self.assertEqual(track['promotionMode'], 'manual')
            self.assertEqual(track['setupCommand'][:2], ['llm-verifier', 'intake'])
            self.assertEqual(track['benchmarkCommand'][:2], ['llm-verifier', 'finalize'])
            for path in [c['rules_path'], *c['rules_files']]:
                self.assertTrue((ROOT / path).is_file(), path)
            selector = c['submission_selector']
            selected = (ROOT / selector['path']).read_text().strip()
            self.assertTrue(selected == selector['directory'] or
                            selected.startswith(selector['directory'] + '/'))
            self.assertTrue((ROOT / selected / 'submission.json').is_file())
            workflow = (ROOT / '.github/workflows' / track['runner']['workflow']).read_text()
            self.assertIn('verify.yml@64b3a2b01f21f1281383c673503964d2fc4c5c47', workflow)
            self.assertIn(f'config_path: tracks/{name}/verifier.json', workflow)
            self.assertIn(f'score_path: {track["scorePath"]}', workflow)
            self.assertIn('yukon_evaluation: ${{ inputs.yukon_evaluation }}', workflow)


if __name__ == '__main__':
    unittest.main()
