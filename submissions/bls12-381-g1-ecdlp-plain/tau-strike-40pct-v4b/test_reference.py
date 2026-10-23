"""Reference regressions: same affine mapping, differently priced arithmetic."""
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent

class JacobianTests(unittest.TestCase):
    def load(self):
        path = HERE/'reference.py'
        self.assertTrue(path.exists(), 'new Jacobian transition is not implemented')
        spec = importlib.util.spec_from_file_location('strike_reference', path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_projective_action_matches_affine_for_every_toy_scalar(self):
        m = self.load()
        self.assertTrue(hasattr(m, 'projective_action'), 'optimized action is missing')
        c = m.toy_context(1543, 523, 3)
        for k in range(c['r']):
            for a in (0, 1, 7, c['r']-1):
                Q = m.mul(a, c['G'], c['p'])
                actual, counts = m.projective_action(Q, k, c['G'], c['p'], census=True)
                expected = m.add(Q, m.mul(k, c['G'], c['p']), c['p'])
                self.assertEqual(actual, expected, (k,a))
                self.assertLessEqual(counts['doublings'], c['r'].bit_length()+1)
                self.assertLessEqual(counts['mixed_calls'], c['r'].bit_length()+1)
        self.assertEqual(m.projective_action(None, 0, c['G'], c['p']), None)

    def test_solver_preserves_affine_mapping_with_identical_coins(self):
        m = self.load()
        self.assertTrue(hasattr(m, 'solve'), 'bounded optimized solver is missing')
        c = m.toy_context(1543, 523, 3)
        seen = []
        for i in range(60):
            x = (i*89) % c['r']
            T = m.mul(x,c['G'],c['p'])
            seed = i.to_bytes(32,'big')
            new = m.solve(T,c,seed,check=True,observe=lambda h,u:seen.append((h,u)))
            old = m.solve(T,c,seed,check=True,affine_steps=True)
            self.assertEqual(new,old)
            self.assertLessEqual(new['evaluations'],6*m.budget(c['r']))
            if new['x'] is not None:
                self.assertEqual(new['x'],x)
        self.assertTrue(seen)
        self.assertEqual(m.wire_output(m.solve(None,c,bytes(32))['x'],c['r']),{'x':'0'*64})

if __name__ == '__main__':
    unittest.main()
