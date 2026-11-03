"""E[N] bound for the two-log descent node counter (DESCENT.md, proof.md s7):"""
def expected_nodes(Pm, dres):
 bins = dres['bins']; DET = dres['DET']; C = {}
 def cmax(side, n):
  return max([C[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
 def miss_c():
  return 1.0   # a missing-log child is one prime special-q node without...
 for b in bins:
  for side in ('f', 'g'):
   d = DET[(b, side)]
   sub = d['kA'] * (cmax('f', d['nA']) + cmax('g', d['nA'])) + d['kB'] * (cmax('f', d['n']) + cmax('g', d['n'])) + d['miss_children'] * miss_c()
   C[(b, side)] = 1.0 + d['retry'] * sub
 bt = dres['boot']
 sub = bt['kA'] * cmax('g', bt['nA']) + cmax('g', float(bt['n1'])) + bt['miss'] * miss_c()
 return 2 * bt['retry'] * sub, C
