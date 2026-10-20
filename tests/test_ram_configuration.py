"""Protected RAM scoring inputs and submission-selector consistency."""
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = 'tracks/docs/COST_MODEL.md'
POLICY = 'bernstein-lange-ram'


def read_json(path):
    return json.loads((ROOT / path).read_text())


class RAMConfiguration(unittest.TestCase):
    def test_all_tracks_bind_the_same_protected_model_and_objective(self):
        manifest = read_json('benchmark.json')
        self.assertEqual({t['name'] for t in manifest['tracks']},
                         {'secp256k1-ecdlp', 'bls12-381-invert-g1'})
        for track in manifest['tracks']:
            name = track['name']
            with self.subTest(track=name):
                config = read_json(f'tracks/{name}/verifier.json')
                self.assertEqual(track['objectives'], [{'name': 'ram_log2', 'direction': '-'}])
                self.assertEqual(len(config['objectives']), 1)
                objective = config['objectives'][0]
                self.assertEqual(objective['name'], 'ram_log2')
                self.assertEqual(objective['direction'], '-')
                self.assertEqual(objective['score_source'], {'claim': 0})
                self.assertFalse(config['allow_unscored'])
                self.assertEqual(config['policies'], track['paretoPolicies'])
                self.assertEqual(config['policies']['scoring'], {'version': POLICY})
                self.assertEqual(config['supported_policies']['scoring'][0]['selector'],
                                 {'version': POLICY})
                self.assertEqual(config['rules_files'],
                                 [f'tracks/{name}/target-profile.json', MODEL_PATH])
                for path in [config['rules_path'], *config['rules_files']]:
                    self.assertTrue((ROOT / path).is_file(), path)
                self.assertNotIn(MODEL_PATH, track['editablePaths'])
                rules = (ROOT / config['rules_path']).read_text()
                self.assertIn('../docs/COST_MODEL.md', rules)
                self.assertIn('39%', rules)

    def test_model_pins_the_paper_and_machine_assumptions(self):
        model = (ROOT / MODEL_PATH).read_text()
        for term in ('nonuniform-20130305.pdf', 'Appendix A.1', '§5.3',
                     '128-bit pointers', '256-bit vector operations',
                     '256-bit floating-point multiplications', 'W = T + L',
                     'ram_log2 = log2(W)', '39%', 'preprocessing',
                     'description', 'memory', 'heuristic'):
            self.assertIn(term, model)

    def test_drafts_have_no_numeric_attack_claim(self):
        for track in read_json('benchmark.json')['tracks']:
            name = track['name']
            with self.subTest(track=name):
                config = read_json(f'tracks/{name}/verifier.json')
                selector = config['submission_selector']
                candidate = (ROOT / selector['path']).read_text().strip()
                self.assertTrue(candidate.startswith(selector['directory'] + '/'))
                manifest = read_json(candidate + '/submission.json')
                self.assertEqual(manifest['version'], 1)
                self.assertTrue(manifest['claims'])
                self.assertTrue(all(isinstance(c, str) for c in manifest['claims']))
                self.assertEqual(manifest['experiments'], [])
                self.assertIn('no score', (ROOT / candidate / 'README.md').read_text())

    def test_shared_harness_and_exploratory_review_are_preserved(self):
        for track in read_json('benchmark.json')['tracks']:
            name = track['name']
            with self.subTest(track=name):
                config = read_json(f'tracks/{name}/verifier.json')
                self.assertEqual(config['acceptance_policy'], 'exploratory-v1')
                self.assertEqual(config['limits']['max_followups'], 1)
                workflow = (ROOT / '.github/workflows' / track['runner']['workflow']).read_text()
                self.assertIn('Layr-Labs/llm-verifier/.github/workflows/verify.yml@', workflow)
                self.assertIn(f'config_path: tracks/{name}/verifier.json', workflow)
                self.assertIn(f'score_path: {track["scorePath"]}', workflow)
                self.assertEqual(track['promotionMode'], 'manual')
                self.assertEqual(track['setupCommand'][0:2], ['llm-verifier', 'intake'])
                self.assertEqual(track['benchmarkCommand'][0:2], ['llm-verifier', 'finalize'])


if __name__ == '__main__':
    unittest.main()
