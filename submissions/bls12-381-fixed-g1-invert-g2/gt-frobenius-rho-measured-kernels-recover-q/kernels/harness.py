"""The four shipped kernels (mul.s, sqr.s from gen.py; lazy.s, root.s from gen_lazy.py) generated in memory,
compiled for emu.py and run on BLS12-381 values held in experiment.py's profile tower (tuples of six Fp2 w^m
coefficients); experiment.py parts K and C drive them through Kernels."""
import emu, gen, gen_lazy

M64 = (1 << 64) - 1
p, MR384 = gen_lazy.P, gen_lazy.RM
RINV = pow(MR384, -1, p)

class Kernels:
	"""Generate mul.s, sqr.s (gen.py) and lazy.s, root.s (gen_lazy.py) in memory, compile them for emu.py
    and run them on BLS12-381 values held in the profile tower (tuples of six Fp2 w^m coefficients)."""

	def __init__(self):
		self.emu, self.gen, self.L = emu, gen, gen_lazy
		a_mul, self.lay_mul, _ = gen.mul_kernel()
		a_sqr, self.lay_sqr, _ = gen.sqr_kernel()
		a_lazy, _ = gen_lazy.lazy_kernel()
		a_root, _ = gen_lazy.root_kernel()
		self.prog, self.static = {}, {}
		for name, a in (("mul", a_mul), ("sqr", a_sqr), ("lazy", a_lazy), ("root", a_root)):
			self.prog[name] = emu.compile_prog(emu.parse(a.text()))
			self.static[name] = sum(a.count().values())
		base = gen_lazy.pool_words()
		self.pool = {}
		for name, lay in (("mul", self.lay_mul), ("sqr", self.lay_sqr), ("lazy", None), ("root", None)):
			w = dict(base)
			if lay is not None:
				for k, Kc in lay["K"].items():
					for j in range(12):
						w[gen.POOL["K"] + 96 * int(k) + 8 * j] = (int(Kc) >> (64 * j)) & M64
			self.pool[name] = w
		self.executed = {k: set() for k in self.prog}

	OUT, XB, EB, PB, SP = 0x10000, 0x20000, 0x30000, 0x40000, 0x80000

	@staticmethod
	def flat(y):
		"""profile tower element -> blst flat coordinates (stored Montgomery form)."""
		v = [0] * 12
		for m in range(6):
			for c in range(2):
				v[6 * (m % 2) + 2 * (m // 2) + c] = y[m][c] * MR384 % p
		return v

	@staticmethod
	def unflat(v):
		return tuple((v[6 * (m % 2) + 2 * (m // 2)] * RINV % p, v[6 * (m % 2) + 2 * (m // 2) + 1] * RINV % p)
					for m in range(6))

	def _run(self, name, mem):
		m = self.emu.Machine()
		m.mem.update({self.PB + k: v for k, v in self.pool[name].items()})
		m.mem.update(mem)
		m.r[0], m.r[1], m.r[2], m.r[3], m.r[31] = self.OUT, self.XB, self.EB, self.PB, self.SP
		n = self.emu.run(self.prog[name], m)
		self.executed[name].add(n)
		return m

	@staticmethod
	def _put(mem, base, vals):
		for i, v in enumerate(vals):
			for j in range(6):
				mem[base + 48 * i + 8 * j] = (v >> (64 * j)) & M64

	def _get(self, m, off, n=1):
		return [sum(m.mem[self.OUT + off + 48 * i + 8 * j] << (64 * j) for j in range(6)) for i in range(n)]

	def mul_entry(self, mval):
		st = self.flat(mval)
		env = {"m%d" % i: st[i] for i in range(12)}
		return [sum(c * env[n] for n, c in L.items()) % p for L in self.lay_mul["entry"]]

	def mul(self, x, entry):
		mem = {}
		self._put(mem, self.XB, self.flat(x))
		self._put(mem, self.EB, entry)
		return self.unflat(self._get(self._run("mul", mem), 0, 12))

	def sqr(self, x):
		mem = {}
		self._put(mem, self.XB, self.flat(x))
		return self.unflat(self._get(self._run("sqr", mem), 0, 12))

	def gamma(self, Pi):
		"""trie-node entry: gamma_0 = Pi_1, gamma_1 = Pi_0, gamma_k = xi Pi_{7-k}, stored; gen_lazy.node_vector."""
		xi = lambda c: ((c[0] - c[1]) % p, (c[0] + c[1]) % p)
		g = [Pi[1], Pi[0]] + [xi(Pi[7 - k]) for k in range(2, 6)]
		return self.L.node_vector([(re * MR384 % p, im * MR384 % p) for re, im in g])

	def base_data(self, x):
		"""per-base data root.s writes: six pre-sums mod p, Xa, Xd (gen_lazy.block_vectors)."""
		return self.L.block_vectors(self.flat(x))

	def key(self, x, gamma=None):
		"""(K, kappa, a, b, per-base data written) from root.s (gamma None) or lazy.s (node vector gamma)."""
		mem = {}
		self._put(mem, self.XB, self.flat(x))
		if gamma is None:
			m = self._run("root", mem)
		else:
			self._put(mem, self.XB + self.L.XS, self.base_data(x))
			self._put(mem, self.EB, gamma)
			m = self._run("lazy", mem)
		o = self.L.OUT
		K = self._get(m, o["K"])[0]
		kap = m.mem[self.OUT + o["kap"]]
		ps = [sum(m.mem[self.XB + self.L.XS + 48 * k + 8 * j] << (64 * j) for j in range(6)) for k in range(10)]
		return K, kap, self._get(m, o["a"])[0], self._get(m, o["b"])[0], ps
