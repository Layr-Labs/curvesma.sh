#!/usr/bin/env python3
"""Tests for this package's own standard-library diagnostic, not production code."""
import json
from pathlib import Path
import subprocess
import sys
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


class DiagnosticTests(unittest.TestCase):
    def test_aligned_recovery_exhaustive_fixture(self):
        result = subprocess.run([sys.executable, '-I', str(HERE / 'toy.py'), '--alignment'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = json.loads(result.stdout)
        self.assertGreater(out['usable_pairs'], 0)
        self.assertEqual(out['aligned_correct'], out['usable_pairs'])
        self.assertGreater(out['unaligned_wrong'], 0)
        self.assertEqual(out['normalized_correct'], out['usable_pairs'])


    def test_floyd_on_every_four_vertex_functional_graph(self):
        import importlib.util
        import itertools
        spec = importlib.util.spec_from_file_location('candidate_toy', HERE / 'toy.py')
        assert spec is not None and spec.loader is not None
        toy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(toy)
        self.assertTrue(callable(getattr(toy, 'floyd_iterations', None)), 'Floyd diagnostic is missing')
        for mapping in itertools.product(range(4), repeat=4):
            for start in range(4):
                seen = {}
                y = start
                while y not in seen:
                    seen[y] = len(seen)
                    y = mapping[y]
                rho = len(seen)
                self.assertLessEqual(toy.floyd_iterations(mapping, start), rho)


    def test_shared_map_pipeline_counts_all_three_evaluations(self):
        result = subprocess.run([sys.executable, '-I', str(HERE / 'mechanism_toy.py'), '--profile', 'smoke', '--seed', '83'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = json.loads(result.stdout)
        self.assertEqual(out['R'], 1024)
        self.assertEqual(out['densities_bits'], [8, 10, 12])
        for row in out['pipeline']:
            self.assertEqual(row['evaluations'], 3 * row['tortoise_iterations'])
            self.assertEqual(row['dense_multiplication_analogues'], 2 * row['evaluations'])
        self.assertEqual(out['incorrect_recoveries'], 0)
        self.assertEqual(out['equivariance_failures'], 0)


    def test_certificate_recomputes_exact_bounds(self):
        result = subprocess.run([sys.executable, '-I', str(HERE / 'verify_exact.py')], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        out = json.loads(result.stdout)
        self.assertGreaterEqual(len(out['checks']), 60)
        self.assertTrue(all(out['checks'].values()))
        self.assertEqual(out['score'], {'numerator': '6443', 'denominator': '50'})


if __name__ == '__main__':
    unittest.main()
