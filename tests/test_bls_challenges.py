"""Checks of the published equations against the pinned arithmetic reference."""
import hashlib
import json
from pathlib import Path
import unittest

import py_ecc
from py_ecc.optimized_bls12_381 import (
    FQ, FQ2, FQ12, G1, G2, Z1, Z2, b, b2, curve_order, eq,
    field_modulus, is_inf, is_on_curve,
    multiply, normalize, pairing,
)

ROOT = Path(__file__).resolve().parents[1]
DLOG = 'bls12-381-g1-ecdlp'
FIXED = 'bls12-381-fixed-g1-invert-g2'
REVERSE = 'bls12-381-fixed-g2-invert-g1'
PLAIN = 'bls12-381-g1-ecdlp-plain'


def profile(name):
    return json.loads((ROOT / 'tracks' / name / 'target-profile.json').read_text())


class BLSChallenges(unittest.TestCase):
    def test_reference_sources_and_curve_parameters(self):
        for name in (DLOG, FIXED, REVERSE, PLAIN):
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
        for name in (DLOG,):
            p = profile(name)
            self.assertEqual(p['auxiliary_inputs']['scoring_n'], 1 << 20)
            self.assertEqual(set(p['input']) - {'n', 'A'},
                             {'T'})
        self.assertEqual((curve_order - 1) % (1 << 20), 0)
        self.assertEqual(set(profile(DLOG)['output']), {'s'})
        self.assertFalse(eq(multiply(G1, s + 1), powers[0]))

    def test_fixed_g1_inversion_requires_the_supplied_point(self):
        p = profile(FIXED)
        self.assertEqual(set(p['input']), {'P', 'z'})
        self.assertEqual(set(p['output']), {'Q'})
        self.assertNotIn('auxiliary_inputs', p)
        # Distinct supplied P values, including the generator; a=0 covers z=1.
        for b_scalar, a in ((1, 0), (3, 7), (11, curve_order - 1)):
            P = multiply(G1, b_scalar)
            Q = multiply(G2, a)
            z = pairing(G2, P) ** a
            self.assertFalse(is_inf(P))
            self.assertTrue(is_on_curve(Q, b2))
            self.assertTrue(is_inf(multiply(Q, curve_order)))
            self.assertEqual(pairing(Q, P), z)
            self.assertEqual(is_inf(Q), z == FQ12.one())
        # A pair acceptable to the any-pair track cannot replace the input P.
        P = multiply(G1, 3)
        Q = multiply(G2, 7)
        z = pairing(Q, P)
        replacement_P = multiply(P, 2)
        replacement_Q = multiply(Q, pow(2, -1, curve_order))
        self.assertEqual(pairing(replacement_Q, replacement_P), z)
        self.assertNotEqual(pairing(replacement_Q, P), z)
        self.assertNotEqual(pairing(Z2, P), z)

    def test_fixed_g2_inversion_requires_the_supplied_point(self):
        p = profile(REVERSE)
        self.assertEqual(set(p['input']), {'Q', 'z'})
        self.assertEqual(set(p['output']), {'P'})
        self.assertNotIn('auxiliary_inputs', p)
        for q_scalar, a in ((1, 0), (3, 7), (11, curve_order - 1)):
            Q = multiply(G2, q_scalar)
            P = multiply(G1, a)
            z = pairing(Q, G1) ** a
            self.assertFalse(is_inf(Q))
            self.assertTrue(is_on_curve(P, b))
            self.assertTrue(is_inf(multiply(P, curve_order)))
            self.assertEqual(pairing(Q, P), z)
            self.assertEqual(is_inf(P), z == FQ12.one())
        P, Q = multiply(G1, 7), multiply(G2, 3)
        z = pairing(Q, P)
        replacement_Q = multiply(Q, 2)
        replacement_P = multiply(P, pow(2, -1, curve_order))
        self.assertEqual(pairing(replacement_Q, replacement_P), z)
        self.assertNotEqual(pairing(Q, replacement_P), z)
        self.assertNotEqual(pairing(Q, Z1), z)

    def test_plain_g1_discrete_log_and_pairing_transfer(self):
        p = profile(PLAIN)
        self.assertEqual(set(p['input']), {'T'})
        self.assertEqual(set(p['output']), {'x'})
        self.assertNotIn('auxiliary_inputs', p)
        generator = pairing(G2, G1)
        for x in (0, 1, 7, curve_order - 1):
            T = multiply(G1, x)
            self.assertTrue(is_on_curve(T, b))
            self.assertTrue(is_inf(multiply(T, curve_order)))
            self.assertEqual(pairing(G2, T), generator ** x)
            self.assertFalse(eq(T, multiply(G1, (x + 1) % curve_order)))
            self.assertEqual(is_inf(T), x == 0)

    def test_subgroup_requirement_is_material(self):
        point = (FQ(0), FQ(2), FQ.one())
        self.assertTrue(is_on_curve(point, b))
        self.assertFalse(is_inf(multiply(point, curve_order)))
        self.assertTrue(is_inf(multiply(G1, curve_order)))
        self.assertTrue(is_inf(multiply(G2, curve_order)))

    def test_all_tracks_budget_initial_and_final_review_retries(self):
        manifest = json.loads((ROOT / 'benchmark.json').read_text())
        for track in manifest['tracks']:
            with self.subTest(track=track['name']):
                config = json.loads((ROOT / 'tracks' / track['name'] / 'verifier.json').read_text())
                provider, limits = config['provider'], config['limits']
                self.assertEqual(config['manifest_preprocessor'], 'curvesmash-v1')
                self.assertEqual(provider['response_format'], 'json_object')
                self.assertEqual(provider['max_attempts'], 3)
                self.assertEqual(provider['timeout_seconds'], 300)
                self.assertEqual(provider['max_completion_tokens'], 32768)
                self.assertEqual(limits['max_followups'], 1)
                self.assertEqual(limits['max_initial_runs'], 1)
                self.assertEqual(limits['total_provider_requests'], 2 * provider['max_attempts'])
                self.assertEqual(limits['total_review_seconds'], 2000)
                self.assertGreater(limits['total_review_seconds'],
                                   limits['total_provider_requests'] * provider['timeout_seconds'])

    def test_yukon_tracks_keep_the_existing_interface(self):
        manifest = json.loads((ROOT / 'benchmark.json').read_text())
        self.assertEqual({t['name'] for t in manifest['tracks']},
                         {'secp256k1-ecdlp', 'p256-ecdlp', DLOG, FIXED, REVERSE, PLAIN})
        for retired in ('bls12-381-invert-g1', 'bls12-381-miller-inversion'):
            self.assertFalse((ROOT / 'tracks' / retired).exists())
            self.assertFalse((ROOT / '.github/workflows' / (retired + '.yml')).exists())
            self.assertNotIn(retired, (ROOT / 'scripts/run_experiments.sh').read_text())
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
            self.assertIn(name, (ROOT / 'scripts/run_experiments.sh').read_text())
            workflow = (ROOT / '.github/workflows' / track['runner']['workflow']).read_text()
            self.assertIn('verify.yml@f7b7a9972f5dda458bb34a02d5b04abad8f2121a', workflow)
            self.assertIn('verifier_revision: f7b7a9972f5dda458bb34a02d5b04abad8f2121a', workflow)
            self.assertIn(f'config_path: tracks/{name}/verifier.json', workflow)
            self.assertIn(f'score_path: {track["scorePath"]}', workflow)
            self.assertIn('yukon_evaluation: ${{ inputs.yukon_evaluation }}', workflow)


if __name__ == '__main__':
    unittest.main()
