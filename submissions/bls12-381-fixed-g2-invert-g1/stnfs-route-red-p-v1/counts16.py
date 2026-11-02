"""expected number of special-q nodes in the two-log descent (DESCENT.md, proof.md section 7), by the same
finite-tree recurrence that descent16.evaluate uses for work, with every node's own cost replaced by 1:
C(node) = 1 + 1/(1 - psi) x (sum over its child slots of the largest C among the bins the child can fall in,
plus the missing-log children at the largest C of a node with q <= 2^(THIN_MAX + 5)), and for one log
retry x (kA C_max(g, nA) + C_max(g, n1) + miss C_miss) below the boot.  Under the same H5 premises that give the
expected two-log work, twice the boot value is an upper bound on E[N], the expected number of special-q nodes
(ordinary and missing-log) run for both logarithms."""
def expected_nodes(Pm, dres):
  bins = dres['bins']; DET = dres['DET']; C = {}
  def cmax(side, n):
    return max([C[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
  def miss_c():
    return max([C[k] for k in C if k[0][0] <= Pm.THIN_MAX] + [0.0])
  for b in bins:
    for side in ('f', 'g'):
      d = DET[(b, side)]
      sub = d['kA'] * (cmax('f', d['nA']) + cmax('g', d['nA'])) + d['kB'] * (cmax('f', d['n']) + cmax('g', d['n'])) + d['miss_children'] * miss_c()
      C[(b, side)] = 1.0 + d['retry'] * sub
  bt = dres['boot']
  sub = bt['kA'] * cmax('g', bt['nA']) + cmax('g', float(bt['n1'])) + bt['miss'] * miss_c()
  return 2 * bt['retry'] * sub, C
