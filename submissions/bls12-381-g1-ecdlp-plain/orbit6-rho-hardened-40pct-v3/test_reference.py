"""Local regression tests; no production-sized walk is run."""
import importlib.util
from pathlib import Path
import unittest

HERE = Path(__file__).resolve().parent

class ReferenceTests(unittest.TestCase):
    def load(self):
        path = HERE / 'reference.py'
        self.assertTrue(path.exists(), 'reference arithmetic has not been implemented')
        spec = importlib.util.spec_from_file_location('orbit6_reference', path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module

    def test_affine_inverse_and_identity(self):
        m = self.load()
        p = 1543
        P = (0, 2)
        self.assertEqual(m.add(P, (0, p-2), p), None)
        self.assertEqual(m.add(None, P, p), P)
        self.assertEqual(m.mul(3, P, p), None)

    def test_signed_class_canonicalizer_is_complete(self):
        m = self.load()
        self.assertTrue(hasattr(m, 'toy_context'), 'class canonicalizer has not been implemented')
        c = m.toy_context(1543, 523, 3)
        P = c['G']
        key = m.canonical(P, c)[0]
        for j in range(3):
            for sign in (1, -1):
                Q = (P[0]*pow(c['beta'], j, c['p']) % c['p'], sign*P[1] % c['p'])
                K, R, u = m.canonical(Q, c)
                self.assertEqual(K, key)
                self.assertEqual(R, m.mul(u, Q, c['p']))
                self.assertEqual(R, m.canonical(P, c)[1])
        self.assertNotEqual(m.canonical(None, c)[0], key)

    def test_floyd_pipeline_recovers_only_verified_scalars(self):
        m = self.load()
        self.assertTrue(hasattr(m, 'solve'), 'bounded Floyd solver has not been implemented')
        c = m.toy_context(1543, 523, 3)
        recovered = 0
        for seed in range(48):
            x = seed % c['r']
            T = m.mul(x, c['G'], c['p'])
            result = m.solve(T, c, seed.to_bytes(32, 'big'), check=True)
            self.assertLessEqual(result['evaluations'], 6*m.budget(c['r']))
            if result['x'] is not None:
                self.assertEqual(result['x'], x)
                recovered += 1
        self.assertGreater(recovered, 1)
        self.assertEqual(m.solve(None, c, bytes(32))['x'], 0)

    def test_wire_output_has_exact_scalar_encoding(self):
        m = self.load()
        self.assertTrue(hasattr(m, 'wire_output'), 'wire scalar encoder has not been implemented')
        self.assertEqual(m.wire_output(0, 523), {'x': '0'*64})
        self.assertEqual(m.wire_output(522, 523), {'x': '0'*61+'20a'})
        with self.assertRaises(ValueError):
            m.wire_output(523, 523)

if __name__ == '__main__':
    unittest.main()
