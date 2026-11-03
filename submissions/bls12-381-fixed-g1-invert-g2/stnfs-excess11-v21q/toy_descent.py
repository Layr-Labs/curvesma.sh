"""Toy special-q descent: measured dependence of recursive failures (DESCENT.md, H5c).

A scaled-down two-sided descent with the production structure:
- two linear "sides" N_s(a, b) = |a + b m_s| (s = 0, 1) with fixed m_s;
- a factor base of primes <= B; a main phase collects B-smooth relations from a random
  region and marks as missing every prime with <= 2 occurrences (the production rule);
- a node is a prime q on one side; its relations are the points (a, b) of the q-lattice
  {a + b m_s = 0 mod q} in a radius-R ball whose norms (q removed) have every prime < q / 2^DROP;
  children = primes > B, plus missing factor-base primes (descended in their own lattice);
- a node succeeds iff one of its relations has all children succeeding (recursively);
- a boot relation is a fresh random integer of BOOT_BITS bits that is q_max-smooth; its
  children are its primes > B (and missing ones).
Each table (main phase) is drawn independently. For each table we evaluate every node's relations
(all of them, not stopping at the first success) and record, per node, the relation failure
indicators; per boot relation, its failure indicator.
Measured:
  * node ratio: sum over nodes of 1[node fails] / sum over nodes of psi_b^M, psi_b the bin's
    marginal relation failure rate and M the node's relation count (the model's phi = E[psi^M]);
  * sibling pair ratio: P(F_i and F_j) / P(F_i) P(F_j) over sibling relation pairs, per bin pooled;
  * boot ratio: P(k boot relations all fail) / p^k for k = 2, 3, 4.
Ratios are pooled over tables; 95% intervals by bootstrap over tables.
Stdin JSON {"tables": n, "seed": s} (defaults 6, 1); stdout one JSON line."""
import json, math, random, sys, time
from functools import lru_cache

def is_prime(n):
	if n < 2: return False
	for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
		if n % p == 0: return n == p
	d, s = n - 1, 0
	while d % 2 == 0: d //= 2; s += 1
	for a in (2, 3, 5, 7, 11, 13, 17):
		x = pow(a, d, n)
		if x in (1, n - 1): continue
		for _ in range(s - 1):
			x = x * x % n
			if x == n - 1: break
		else:
			return False
	return True

def rho(n):
	if n % 2 == 0: return 2
	for c in range(1, 50):
		x = y = 2; d = 1
		f = lambda v: (v * v + c) % n
		while d == 1:
			x = f(x); y = f(f(y)); d = math.gcd(abs(x - y), n)
		if d != n: return d
	raise ValueError(n)

SMALL = [p for p in range(2, 1 << 10) if is_prime(p)]
def factor(n, bound):
	"""prime factors of n if every prime is < bound, else None."""
	out = []
	for p in SMALL:
		if p * p > n: break
		while n % p == 0: out.append(p); n //= p
	if n == 1: return out
	stack = [n]
	while stack:
		m = stack.pop()
		if m == 1: continue
		if is_prime(m):
			if m >= bound: return None
			out.append(m); continue
		d = rho(m); stack += [d, m // d]
	return out

class Toy:
	def __init__(self, P, table_seed):
		self.P = P
		rng = random.Random(table_seed)
		B = P['B']; self.B = B
		self.m = P['m']
		# main phase: random region points, B-smooth on both sides
		occ = [dict(), dict()]
		got = 0; tries = 0
		while got < P['main_rels']:
			tries += 1
			a = rng.randint(-P['main_A'], P['main_A']); b = rng.randint(1, P['main_A'])
			if math.gcd(a, b) != 1: continue
			fs = [factor(abs(a + b * m), B + 1) for m in self.m]
			if fs[0] is None or fs[1] is None or 0 in (a + b * self.m[0], a + b * self.m[1]): continue
			got += 1
			for s in (0, 1):
				for p in set(fs[s]): occ[s][p] = occ[s].get(p, 0) + 1
		self.missing = [set(p for p in SMALL + [p for p in range(1 << 10, B + 1) if is_prime(p)] if p <= B and occ[s].get(p, 0) <= 2) for s in (0, 1)]
		self.cache = {}

	def lattice(self, q, s):
		# basis of {(a, b): a + b m = 0 mod q}: (q, 0), (-m mod q, 1); Gauss-reduce
		u = (q, 0); v = ((-self.m[s]) % q, 1)
		def n2(w): return w[0] * w[0] + w[1] * w[1]
		if n2(u) < n2(v): u, v = v, u
		while True:
			if n2(u) < n2(v): u, v = v, u
			k = round((u[0] * v[0] + u[1] * v[1]) / n2(v)) if n2(v) else 0
			w = (u[0] - k * v[0], u[1] - k * v[1])
			if n2(w) >= n2(v): return v, w
			u = w

	def relations(self, q, s):
		"""relations of node (q, s): list of child lists."""
		R = self.P['R']; drop = 2 ** self.P['DROP']
		v, w = self.lattice(q, s)
		bound = max(self.B + 1, int(q / drop))
		rels = []; seen = set()
		for i in range(-R, R + 1):
			for j in range(0, R + 1):
				if i * i + j * j > R * R or (j == 0 and i <= 0): continue
				a = i * v[0] + j * w[0]; b = i * v[1] + j * w[1]
				if b < 0: a, b = -a, -b
				if b == 0 or math.gcd(a, b) != 1 or (a, b) in seen: continue
				seen.add((a, b))
				n0 = abs(a + b * self.m[0]); n1 = abs(a + b * self.m[1])
				if n0 == 0 or n1 == 0: continue
				ns = [n0, n1]; ns[s] //= q
				fs = [factor(n, bound) for n in ns]
				if fs[0] is None or fs[1] is None: continue
				kids = set()
				for t in (0, 1):
					for p in set(fs[t]):
						if p > self.B or p in self.missing[t]: kids.add((p, t))
				kids.discard((q, s))
				rels.append(sorted(kids))
		return rels

	def node(self, q, s, depth=0):
		key = (q, s)
		if key in self.cache: return self.cache[key]['ok']
		self.cache[key] = dict(ok=False, rels=[], fails=[])   # cycle guard: in-progress counts as failure
		if depth > 40: return False
		rels = self.relations(q, s)
		fails = []
		for kids in rels:
			ok = True
			for (p, t) in kids:
				if not self.node(p, t, depth + 1): ok = False
			fails.append(not ok)
		res = dict(ok=any(not f for f in fails), rels=len(rels), fails=fails)
		self.cache[key] = res
		return res['ok']

def run(P, tables, seed, workers=1):
	t0 = time.time()
	per_table = []
	jobs = [(P, seed, ti) for ti in range(tables)]
	if workers > 1:
		from multiprocessing import Pool
		with Pool(workers) as pool: per_table = pool.map(one_table, jobs)
	else:
		per_table = [one_table(j) for j in jobs]
	return per_table, time.time() - t0
def one_table(job):
	P, seed, ti = job
	if True:
		T = Toy(P, 1000003 * seed + ti)
		rng = random.Random(7 * seed + 31 * ti + 5)
		boots = []
		nb = 0
		while len(boots) < P['boot_rels'] and nb < 200000:
			nb += 1
			z = rng.getrandbits(P['BOOT_BITS']) | (1 << (P['BOOT_BITS'] - 1))
			fz = factor(z, P['qmax'])
			if fz is None: continue
			kids = sorted(set((p, 0) for p in fz if p > T.B or p in T.missing[0]))
			boots.append(not all(T.node(p, t) for p, t in kids))
		nodes = [(k, v) for k, v in T.cache.items() if isinstance(v['rels'], int)]
		return dict(boots=boots, nodes=[(k[0].bit_length(), v['fails']) for k, v in nodes])

def stats(per_table, idx):
	# pooled marginals per bit size
	cnt = {}; fl = {}
	for i in idx:
		for bits, fails in per_table[i]['nodes']:
			cnt[bits] = cnt.get(bits, 0) + len(fails); fl[bits] = fl.get(bits, 0) + sum(fails)
	psi = {b: fl[b] / cnt[b] for b in cnt if cnt[b]}
	obs = pred = 0.0; pj = pp = 0.0; tobs = tmod = 0.0; pbin = {}
	for i in idx:
		for bits, fails in per_table[i]['nodes']:
			M = len(fails)
			if not M: continue
			S = M - sum(fails); ps = psi[bits]
			obs += S == 0; pred += ps ** M
			x = pbin.setdefault(bits, [0.0, 0.0]); x[0] += S == 0; x[1] += ps ** M
			if ps < 1: tobs += (M + 1) / (S + 1) if S else M; tmod += 1 / (1 - ps)
			for a in range(M):
				for b in range(a + 1, M):
					pj += fails[a] and fails[b]; pp += psi[bits] ** 2
	# chain-rule conditional ratios: P(F_k | F_1..F_{k-1}) / psi, siblings in node order, k = 2..4
	num = {k: 0.0 for k in (2, 3, 4)}; den = {k: 0.0 for k in (2, 3, 4)}
	for i in idx:
		for bits, fails in per_table[i]['nodes']:
			for k in (2, 3, 4):
				if len(fails) >= k and all(fails[:k - 1]):
					num[k] += fails[k - 1]; den[k] += psi[bits]
	chain = {k: (num[k] / den[k] if den[k] else float('nan')) for k in (2, 3, 4)}
	boots = [x for i in idx for x in per_table[i]['boots']]
	p = sum(boots) / len(boots)
	bk = {}
	for k in (2, 3, 4):
		allk = 0; n = 0
		for i in idx:
			bt = per_table[i]['boots']
			for j in range(0, len(bt) - k + 1, k):
				n += 1; allk += all(bt[j:j + k])
		bk[k] = (allk / n) / (p ** k) if n and p > 0 else float('nan')
	ps_t = [sum(per_table[i]['boots']) / len(per_table[i]['boots']) for i in idx]
	pb_ = sum(sum(per_table[i]['boots']) for i in idx) / sum(len(per_table[i]['boots']) for i in idx)
	excess = (sum((x - pb_) ** 2 for x in ps_t) / (len(ps_t) - 1) - sum(pb_ * (1 - pb_) / len(per_table[i]['boots']) for i in idx) / len(idx)) if len(idx) > 1 else float('nan')
	return dict(node_ratio=obs / pred if pred else float('nan'), node_ratio_by_bits={b: round(v[0] / v[1], 4) for b, v in sorted(pbin.items()) if v[1] > 0},
	            tries_ratio=tobs / tmod if tmod else float('nan'), boot_excess_var=excess, pair_ratio=pj / pp if pp else float('nan'), chain_ratio=chain,
	            boot_p=p, boot_ratio=bk, boot_table_max_ratio=max(ps_t) / pb_ if pb_ > 0 else float('nan'), nodes=sum(len(per_table[i]['nodes']) for i in idx))

def boot_ci(per_table, key, reps=400, seed=1):
	rng = random.Random(seed); n = len(per_table); vals = []
	for _ in range(reps):
		idx = [rng.randrange(n) for _ in range(n)]
		s = stats(per_table, idx)
		v = s[key] if not isinstance(key, tuple) else s[key[0]][key[1]]
		if v == v: vals.append(v)
	vals.sort()
	return (vals[int(0.025 * len(vals))], vals[int(0.975 * len(vals)) - 1]) if vals else (float('nan'),) * 2

DEFAULT = dict(B=2 ** 8, m=(40961, 65537 * 3 + 2), main_rels=600, main_A=3000, R=6, DROP=1, qmax=2 ** 15,
               BOOT_BITS=40, boot_rels=40)

if __name__ == "__main__":
	req = json.load(sys.stdin) if not sys.stdin.isatty() else {}
	P = dict(DEFAULT); P.update(req.get("params", {}))
	per_table, dt = run(P, int(req.get("tables", 6)), int(req.get("seed", 1)), int(req.get("workers", 1)))
	idx = list(range(len(per_table))); s = stats(per_table, idx)
	out = dict(params={k: (list(v) if isinstance(v, tuple) else v) for k, v in P.items()}, tables=len(per_table), seconds=round(dt, 1), **{k: (v if not isinstance(v, float) else round(v, 4)) for k, v in s.items()})
	allb = {}
	for t in per_table:
		for bits, fails in t['nodes']:
			c = allb.setdefault(bits, [0, 0, 0]); c[0] += 1; c[1] += len(fails); c[2] += sum(fails)
	out['bins'] = {b: (c[0], round(c[1] / c[0], 2), round(c[2] / max(c[1], 1), 3)) for b, c in sorted(allb.items())}
	out['node_ratio_ci95'] = boot_ci(per_table, 'node_ratio'); out['tries_ratio_ci95'] = boot_ci(per_table, 'tries_ratio'); out['boot_excess_var_ci95'] = boot_ci(per_table, 'boot_excess_var'); out['pair_ratio_ci95'] = boot_ci(per_table, 'pair_ratio')
	out['boot_ratio_ci95'] = {k: boot_ci(per_table, ('boot_ratio', k)) for k in (2, 3, 4)}
	out['chain_ratio_ci95'] = {k: boot_ci(per_table, ('chain_ratio', k)) for k in (2, 3, 4)}
	print(json.dumps(out))
