#!/usr/bin/env python3
"""Regression tests for the actual-field minimum-coordinate diagnostic."""
import importlib.util
from pathlib import Path
import random
import sys
import unittest
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent

class MechanismTests(unittest.TestCase):
    def test_actual_order12_field_and_four_product_minimum(self):
        path = HERE / 'mechanism_toy.py'
        self.assertTrue(path.exists(), 'actual-field mechanism diagnostic is missing')
        spec = importlib.util.spec_from_file_location('mechanism', path)
        assert spec is not None and spec.loader is not None
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        field = m.Tower(m.PARAMETERS[0])
        rng = random.Random(83)
        for _ in range(20):
            e = rng.randrange(1, field.r)
            y = field.element(e)
            self.assertEqual(field.power(y, field.r), field.one)
            self.assertEqual(field.frobenius(y, 1), field.power(y, field.p))
            key, k = field.canonical(y)
            slow = [(field.frobenius(y, j)[2] * field.mont % field.p, j) for j in range(12)]
            self.assertEqual(key, min(z[0] for z in slow))
            for j in range(12):
                ky, kj = field.canonical(field.frobenius(y, j))
                self.assertEqual(ky, key)
                self.assertEqual(field.lp[kj] * field.lp[j] * e % field.r, field.lp[k] * e % field.r)
            self.assertEqual(field.mark(key, 8), ((key >> 10) & 255) == 0)
            self.assertEqual(field.mark(key, 12), ((key >> 10) & 4095) == 0)

    def test_field_states_match_tracked_labels_on_every_scale(self):
        spec = importlib.util.spec_from_file_location('mechanism', HERE / 'mechanism_toy.py')
        assert spec is not None and spec.loader is not None
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        rng = random.Random(97)
        for parameter in m.PARAMETERS:
            f = m.Tower(parameter)
            walk = m.SharedMap(f, rng)
            state = walk.fresh()
            for _ in range(40):
                y,a,b,K,k = state
                self.assertEqual(y, f.element((a+walk.x*b) % f.r))
                c,A,B = walk.normalized(state)
                self.assertEqual(f.frobenius(y,k), f.element(c))
                self.assertEqual(c, (A+walk.x*B) % f.r)
                state = walk.step(state)

    def test_stationary_shared_map_long_cycles_and_sparse_mark_counts(self):
        spec = importlib.util.spec_from_file_location('mechanism', HERE / 'mechanism_toy.py')
        assert spec is not None and spec.loader is not None
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        self.assertTrue(callable(getattr(m, 'experiment', None)), 'cycle/DP experiment missing')
        out = m.experiment(seed=83, profile='smoke', wall_seconds=30)
        self.assertEqual(out['partition'], 'minimum real w coefficient over all 12 actual Frobenius conjugates; low 10 bits')
        self.assertEqual(out['equivariance_failures'], 0)
        self.assertEqual(out['incorrect_recoveries'], 0)
        self.assertGreaterEqual(out['cycles_observed'], 32)
        self.assertEqual(out['densities_bits'], [8, 10, 12])
        for row in out['pipeline']:
            self.assertEqual(row['evaluations'], 3*row['tortoise_iterations'])
            self.assertEqual(row['dense_multiplication_analogues'], 2*row['evaluations'])
            self.assertGreater(row['tortoise_iterations'], 81)
        self.assertTrue(out['completed'])

    def test_cdf_grid_uses_exact_integer_comparisons(self):
        path = HERE / 'analyze_evidence.py'
        self.assertTrue(path.exists(), 'exact empirical-CDF analyzer missing')
        spec = importlib.util.spec_from_file_location('evidence_analysis', path)
        assert spec is not None and spec.loader is not None
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        self.assertEqual(m.grid_count([{'rho':1},{'rho':2},{'rho':3}], 4, 1, 1, 'rho'), 2)
        self.assertEqual(m.grid_count([{'rho':1},{'rho':2},{'rho':3}], 4, 1, 2, 'rho'), 1)
        import json
        captured = HERE / 'mechanism-output.json'
        if captured.exists():
            summary = m.summarize(json.loads(captured.read_text()))
            self.assertEqual(summary['independent_table_first_cycles'], 240)
            self.assertEqual(len(summary['empirical_CDF']), 15)

if __name__ == '__main__':
    unittest.main()
