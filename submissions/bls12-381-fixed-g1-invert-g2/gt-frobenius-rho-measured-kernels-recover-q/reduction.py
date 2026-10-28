"""Part F: the instance reduction and the output on BLS12-381, run by experiment.py (proof.md sections 1 and 5).

g = e(P,Q0) for the supplied P (the pinned order pairing(Q0,P)), the operation-by-operation count of that
pairing, Q = [x]Q0 in G2 with the hard maximum of its counted ladder, the verification e(P,Q) = z, and the
identity target z = 1. experiment.py passes p, r and its Tower class in; nothing else is shared.
Every counted operation (an Fp, Fp2 or Fp12 multiplication, addition, subtraction, scaling, comparison or
twist) is charged one full unit, a dense Fp12 multiplication, whatever the field: an overcount."""

P0_XY = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
		0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
Q0_XY = ((0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,
		 0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e),
		(0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,
		 0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be))
ATE_PBE = [0] * 16 + [1] + [0] * 31 + [1] + [0] * 8 + [1, 0, 0, 1, 0, 1, 1]   # the pinned pseudo_binary_encoding

class Fld:
	"""One field of the pinned library's point and line formulas; every call is one counted operation."""

	def __init__(self, kind):
		self.kind, self.n = kind, 0
		if kind == 1:
			self.zero, self.one = 0, 1
		elif kind == 2:
			self.zero, self.one = (0, 0), (1, 0)
		else:
			self.zero, self.one = (0,) * 12, (1,) + (0,) * 11

	def mul(self, a, b):
		self.n += 1
		if self.kind == 1:
			return a * b % p
		if self.kind == 2:
			return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
		acc = [0] * 23
		for i in range(12):
			if a[i]:
				for j in range(12):
					acc[i + j] += a[i] * b[j]
		for k in range(22, 11, -1):            # w^12 = 2 w^6 - 2 (the profile modulus)
			c = acc[k]
			acc[k - 6] += 2 * c
			acc[k - 12] -= 2 * c
		return tuple(v % p for v in acc[:12])

	def lin(self, a, b, sgn):
		self.n += 1
		if self.kind == 1:
			return (a + sgn * b) % p
		return tuple((u + sgn * v) % p for u, v in zip(a, b))

	def add(self, a, b):
		return self.lin(a, b, 1)

	def sub(self, a, b):
		return self.lin(a, b, -1)

	def smul(self, k, a):
		self.n += 1
		if self.kind == 1:
			return k * a % p
		return tuple(k * v % p for v in a)

	def eq(self, a, b):
		self.n += 1
		return a == b

def pt_double(F, pt):                          # pinned optimized_curve.double, operation for operation
	x, y, z = pt
	W = F.mul(F.smul(3, x), x)
	S = F.mul(y, z)
	B = F.mul(F.mul(x, y), S)
	H = F.sub(F.mul(W, W), F.smul(8, B))
	S2 = F.mul(S, S)
	nx = F.mul(F.smul(2, H), S)
	ny = F.sub(F.mul(W, F.sub(F.smul(4, B), H)), F.mul(F.mul(F.smul(8, y), y), S2))
	nz = F.mul(F.smul(8, S), S2)
	return (nx, ny, nz)

def pt_add(F, p1, p2):                         # pinned optimized_curve.add, operation for operation
	if F.eq(p1[2], F.zero) or F.eq(p2[2], F.zero):
		return p1 if p2[2] == F.zero else p2
	x1, y1, z1 = p1
	x2, y2, z2 = p2
	U1, U2, V1, V2 = F.mul(y2, z1), F.mul(y1, z2), F.mul(x2, z1), F.mul(x1, z2)
	if F.eq(V1, V2) and F.eq(U1, U2):
		return pt_double(F, p1)
	elif F.eq(V1, V2):
		return (F.one, F.one, F.zero)
	U, V = F.sub(U1, U2), F.sub(V1, V2)
	V2s = F.mul(V, V)
	V2sV2 = F.mul(V2s, V2)
	V3 = F.mul(V, V2s)
	W = F.mul(z1, z2)
	A = F.sub(F.sub(F.mul(F.mul(U, U), W), V3), F.smul(2, V2sV2))
	return (F.mul(V, A), F.sub(F.mul(U, F.sub(V2sV2, A)), F.mul(V3, U2)), F.mul(V3, W))

def twist12(F, pt):                            # pinned optimized_curve.twist (coefficient placement)
	F.n += 1
	(x0, x1), (y0, y1), (z0, z1) = pt
	nx = [0] * 12; nx[1], nx[7] = (x0 - x1) % p, x1
	ny = [0] * 12; ny[0], ny[6] = (y0 - y1) % p, y1
	nz = [0] * 12; nz[3], nz[9] = (z0 - z1) % p, z1
	return (tuple(nx), tuple(ny), tuple(nz))

def linefunc(F, P1, P2, T):                    # pinned optimized_pairing.linefunc
	x1, y1, z1 = P1
	x2, y2, z2 = P2
	xt, yt, zt = T
	mn = F.sub(F.mul(y2, z1), F.mul(y1, z2))
	md = F.sub(F.mul(x2, z1), F.mul(x1, z2))
	if not F.eq(md, F.zero):
		pass
	elif F.eq(mn, F.zero):
		mn = F.mul(F.smul(3, x1), x1)
		md = F.mul(F.smul(2, y1), z1)
	else:
		return F.sub(F.mul(xt, z1), F.mul(x1, zt)), F.mul(z1, zt)
	return (F.sub(F.mul(mn, F.sub(F.mul(xt, z1), F.mul(x1, zt))), F.mul(md, F.sub(F.mul(yt, z1), F.mul(y1, zt)))),
			F.mul(F.mul(md, zt), z1))

def to_tower(c):                               # profile basis -> blst tower (u = w^6 - 1): a ring isomorphism
	return tuple(((c[m] + c[m + 6]) % p, c[m + 6]) for m in range(6))

def from_tower(y):
	return tuple((a - b) % p for a, b in y) + tuple(b for _, b in y)

def tpow(Tw, x, e, cnt):                       # square-and-multiply as Tower.pow, one unit per product
	out = Tw.one
	while e:
		if e & 1:
			out = Tw.mul(out, x); cnt[0] += 1
		x = Tw.mul(x, x); cnt[0] += 1
		e >>= 1
	return out

def pairing_pinned(Tw, Qpt, Ppt):
	"""e(P,Q) as the pinned library defines it: pairing(Q,P) = miller_loop(Q,P) then f**((p^12-1)/r).
    Qpt is projective over Fp2, Ppt projective over Fp; None stands for infinity. Returns (value in the
    tower basis, counted units)."""
	if Qpt is None or Ppt is None:
		return Tw.one, 0
	F2, F12 = Fld(2), Fld(12)
	F12.n += 1                                 # cast_point_to_fq12
	cP = tuple((c,) + (0,) * 11 for c in Ppt)
	tR = tQ = twist12(F12, Qpt)
	R = Qpt
	fn, fd = F12.one, F12.one
	for v in ATE_PBE[62::-1]:
		n_, d_ = linefunc(F12, tR, tR, cP)
		fn = F12.mul(F12.mul(fn, fn), n_)
		fd = F12.mul(F12.mul(fd, fd), d_)
		R = pt_double(F2, R)
		tR = twist12(F12, R)
		if v == 1:
			n_, d_ = linefunc(F12, tR, tQ, cP)
			fn = F12.mul(fn, n_)
			fd = F12.mul(fd, d_)
			R = pt_add(F2, R, Qpt)
			tR = twist12(F12, R)
	cnt = [0]
	a, d = to_tower(fn), to_tower(fd)          # the change of basis, charged as two operations
	cnt[0] += 2
	conj = Tw.one                              # f_den^-1 = (prod of its 11 other conjugates) / norm, norm in Fp
	for k in range(1, 12):
		conj = Tw.mul(conj, Tw.frobk(d, k)); cnt[0] += 2
	nrm = Tw.mul(conj, d); cnt[0] += 1
	if any(nrm[m] != (0, 0) for m in range(1, 6)) or nrm[0][1] != 0:
		raise RuntimeError('norm of f_den is not in Fp')
	e, inv = p - 2, 1                          # Fermat inversion in Fp, one unit per Fp product
	base = nrm[0][0]
	while e:
		if e & 1:
			inv = inv * base % p; cnt[0] += 1
		base = base * base % p; cnt[0] += 1
		e >>= 1
	dinv = tuple(((u * inv) % p, (v * inv) % p) for u, v in conj); cnt[0] += 1
	f = Tw.mul(a, dinv); cnt[0] += 1
	g = tpow(Tw, f, (p ** 12 - 1) // r, cnt)
	return g, F2.n + F12.n + cnt[0]

def g2_mul_counted(x):
	"""Q = [x]Q0 by left-to-right double-and-add with the pinned projective formulas over Fp2, then the
    affine encoding (one Fp2 inversion). Returns (affine Q or None for infinity, ladder units,
    encoding units, doublings, additions)."""
	F = Fld(2)
	Q0p = (Q0_XY[0], Q0_XY[1], (1, 0))
	if x == 0:
		return None, 0, 0, 0, 0
	Rp, nd, na = Q0p, 0, 0
	for bit in bin(x)[3:]:
		Rp = pt_double(F, Rp); nd += 1
		if bit == "1":
			Rp = pt_add(F, Rp, Q0p); na += 1
	ladder = F.n
	zx = Rp[2]
	if F.eq(zx, F.zero):
		return None, ladder, F.n - ladder, nd, na
	nrm = F.add(F.mul((zx[0], 0), (zx[0], 0)), F.mul((zx[1], 0), (zx[1], 0)))[0]
	e, inv, base = p - 2, 1, nrm
	while e:
		if e & 1:
			inv = inv * base % p; F.n += 1
		base = base * base % p; F.n += 1
		e >>= 1
	zi = F.mul((zx[0], (-zx[1]) % p), (inv, 0))
	return (F.mul(Rp[0], zi), F.mul(Rp[1], zi)), ladder, F.n - ladder, nd, na

def on_g2_curve(Qa):
	F = Fld(2)
	x, y = Qa
	return F.mul(y, y) == F.add(F.mul(F.mul(x, x), x), (4, 4))

def g2_branch_maxima():
	"""Operation count of each branch of the pinned double/add on Fp2 points, so the hard maximum per
    ladder step is the largest branch, not the branch the sample happened to take."""
	F = Fld(2)
	Q0p = (Q0_XY[0], Q0_XY[1], (1, 0))
	negQ0 = (Q0_XY[0], (-Q0_XY[1][0] % p, -Q0_XY[1][1] % p), (1, 0))
	D = pt_double(F, Q0p)
	out = {}
	for name, a, b in (("generic", D, Q0p), ("equal_points", Q0p, Q0p), ("opposite_points", Q0p, negQ0),
					 ("infinity_operand", (F.one, F.one, F.zero), Q0p)):
		n0 = F.n; pt_add(F, a, b); out["add_" + name] = F.n - n0
	n0 = F.n; pt_double(F, D); out["double"] = F.n - n0
	return out

def run(rng):
	"""Instance reduction and output on BLS12-381: g = e(P,Q0) per instance (pinned order pairing(Q0,P)),
    Q = [x]Q0 in G2 with its counted cost, verification e(P,Q) = z, and the identity target z = 1."""
	Tw = Tower(p)
	Fp1 = Fld(1)
	b = rng.randrange(1, r)
	Pp, acc = None, (P0_XY[0], P0_XY[1], 1)
	for bit in bin(b)[:1:-1]:                  # the challenge's own sampling of P = [b]P0 (not attack work)
		if bit == "1":
			Pp = acc if Pp is None else pt_add(Fp1, Pp, acc)
		acc = pt_double(Fp1, acc)
	zi = pow(Pp[2], p - 2, p)
	Pa = (Pp[0] * zi % p, Pp[1] * zi % p)
	P_ok = (Pa[1] * Pa[1] - Pa[0] ** 3 - 4) % p == 0
	g, units_g = pairing_pinned(Tw, (Q0_XY[0], Q0_XY[1], (1, 0)), (Pa[0], Pa[1], 1))
	g_ok = tpow(Tw, g, r, [0]) == Tw.one and g != Tw.one
	x = rng.randrange(1, r)                    # stands in for the logarithm the walk returns
	z = tpow(Tw, g, x, [0])
	Qa, lad_q, enc_q, nd, na = g2_mul_counted(x)
	Qr = g2_mul_counted(r)[0]                  # acceptance-side check [r]Q0 = O (not attack work)
	Qrx = g2_mul_counted(r * x)[0] if Qa is not None else None
	v, units_v = pairing_pinned(Tw, (Qa[0], Qa[1], (1, 0)), (Pa[0], Pa[1], 1))
	one_v, _ = pairing_pinned(Tw, None, (Pa[0], Pa[1], 1))
	br = g2_branch_maxima()
	per_add = max(v_ for k_, v_ in br.items() if k_.startswith("add_"))
	steps = (r - 1).bit_length() - 1           # x <= r-1 < 2^255: at most 254 doublings and 254 additions
	q_max = steps * (br["double"] + per_add) + enc_q
	return {
		"P_on_G1_curve_and_nonidentity": P_ok and Pp[2] != 0,
		"g_in_GT_and_not_1": g_ok,
		"pairing_units_counted_g": units_g,
		"pairing_units_counted_verification": units_v,
		"pairing_counts_equal": units_g == units_v,
		"final_exponent_bits": ((p ** 12 - 1) // r).bit_length(),
		"Q_on_G2_curve": on_g2_curve(Qa),
		"r_times_Q0_is_infinity": Qr is None,
		"r_times_Q_is_infinity": Qrx is None,
		"verification_e(P,Q)_equals_z": v == z,
		"g2_ladder_doublings_additions_this_x": [nd, na],
		"g2_ladder_units_this_x": lad_q,
		"g2_ops_per_double": br["double"],
		"g2_ops_per_add_by_branch": {k_[4:]: v_ for k_, v_ in br.items() if k_.startswith("add_")},
		"g2_affine_encoding_units": enc_q,
		"g2_ladder_steps_max": steps,
		"g2_output_units_hard_max_255_bit_scalar": q_max,
		"identity_target": {"z_equals_1_test": "before any walk", "x": 0, "Q": "infinity",
							"e(P,O)_equals_1": one_v == Tw.one},
	}
