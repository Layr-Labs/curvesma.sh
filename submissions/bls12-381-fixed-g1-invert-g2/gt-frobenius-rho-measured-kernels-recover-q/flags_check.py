"""Instruction-level check of the multiword compares and increments of proof.md 0-i and 0-iv.

A small AArch64 interpreter (64-bit X registers, NZCV) runs each listed sequence exactly as proof.md
writes it and asserts that it implements its stated predicate, on edge cases and on random cases.
Flags follow the architecture's AddWithCarry pseudocode:
  ADDS/ADCS: AddWithCarry(a, b, carry_in), carry_in 0 for ADDS, C for ADCS (ADC: same sum, no flags);
  SUBS/SBCS/CMP: AddWithCarry(a, NOT b, carry_in), carry_in 1 for SUBS/CMP, C for SBCS (SBC: same, no flags);
  N = bit 63, Z = (64-bit result == 0), C = carry out of bit 63, V = signed overflow.
Conditions: hs = C, lo = !C, hi = C && !Z, ls = !(C && !Z), eq = Z, ne = !Z.
After a borrow chain C is the full-width "no borrow" result, while Z reflects only the last word, so only
a carry-only condition is exact there; the old b.hi meter compare is run too and shown to fail.
The heapsort level's two 10-word record compares (proof.md 0-ii) run on 80-byte records in memory,
and the key combination and key minimum of section 5 (KEY_SEG1, KEY_SEG2, KEY_MIN) on frame, pool and
walk-state memory.
Run by experiment.py (output key nzcv_flags_check); standalone: python3 -I flags_check.py.
"""
import json, random

M64 = (1 << 64) - 1
COND = {"hs": lambda f: f["C"], "lo": lambda f: not f["C"], "hi": lambda f: f["C"] and not f["Z"],
		"ls": lambda f: not (f["C"] and not f["Z"]), "eq": lambda f: f["Z"], "ne": lambda f: not f["Z"]}

def awc(a, b, cin):
	"""AddWithCarry on 64-bit words: (result, flags)."""
	u = a + b + cin
	r = u & M64
	sa, sb, sr = a >> 63, b >> 63, r >> 63
	return r, {"N": bool(sr), "Z": r == 0, "C": u > M64, "V": sa == sb and sr != sa}

def run_seq(seq, regs, mem=None):
	"""Run a ';'-separated listing; return (registers, whether the final conditional branch was taken).
    mem maps byte addresses to 64-bit words (ldr, ldp read it; stp writes it); adr reads a label from regs."""
	regs, flags, taken = dict(regs), {"N": False, "Z": False, "C": False, "V": False}, None

	def val(t):
		t = t.strip()
		if t == "xzr":
			return 0
		if t.startswith("#"):
			return int(t[1:], 0) & M64
		return regs[t]

	def put(d, v):
		if d != "xzr":
			regs[d] = v & M64

	for ins in [s.strip() for s in seq.split(";") if s.strip()]:
		op, _, rest = ins.partition(" ")
		if op in ("ldr", "ldp", "stp"):
			rs, _, m = rest.partition("[")
			base, _, off = m.strip("]").partition(",")
			ad, rs = (val(base) + (val(off) if off else 0)) & M64, [x.strip() for x in rs.split(",") if x.strip()]
			for k, x in enumerate(rs):
				if op == "stp":
					mem[ad + 8 * k] = val(x)
				else:
					put(x, mem[ad + 8 * k])
			continue
		a = [x.strip() for x in rest.split(",")]
		if op.startswith("b."):
			assert taken is None, "one conditional branch, at the end"
			taken = bool(COND[op[2:]](flags))
		elif op == "mov":
			put(a[0], val(a[1]))
		elif op == "and":
			put(a[0], val(a[1]) & val(a[2]))
		elif op == "add":
			put(a[0], val(a[1]) + val(a[2]))
		elif op == "adc":
			put(a[0], val(a[1]) + val(a[2]) + int(flags["C"]))
		elif op == "sbc":
			put(a[0], val(a[1]) + (~val(a[2]) & M64) + int(flags["C"]))
		elif op == "adr":
			put(a[0], regs[a[1]])
		elif op in ("adds", "adcs"):
			r, flags = awc(val(a[1]), val(a[2]), 0 if op == "adds" else int(flags["C"]))
			put(a[0], r)
		elif op in ("subs", "sbcs"):
			r, flags = awc(val(a[1]), ~val(a[2]) & M64, 1 if op == "subs" else int(flags["C"]))
			put(a[0], r)
		elif op == "csel":
			put(a[0], val(a[1]) if COND[a[3]](flags) else val(a[2]))
		elif op == "cmp":
			_, flags = awc(val(a[0]), ~val(a[1]) & M64, 1)
		else:
			raise ValueError("unsupported instruction " + ins)
	return regs, taken

def words(v, n):
	return [(v >> (64 * i)) & M64 for i in range(n)]

def join(ws):
	return sum(w << (64 * i) for i, w in enumerate(ws))

METER_CMP = "cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7; b.hs end_walk_phase"
OLD_METER_CMP = "cmp x0, x5; sbcs xzr, x1, x6; sbcs xzr, x2, x7; b.hi end_walk_phase"
METER_ADD = "adds x0, x0, x4; adcs x1, x1, xzr; adc x2, x2, xzr"
STUCK = "adds x8, x8, #1; adc x9, x9, xzr; mov x16, #5120; cmp x9, x16; b.hs restart"
REPORT_CMP = "cmp x10, x11; b.hs"
CHUNK = "subs x12, x12, #1; b.eq switch"
STUCK_CAP = 5120 << 64
# Sampler compare (proof.md 0-iv): accept iff X < r (X in x3..x0, r in x7..x4, top word first);
# after the chain C = 1 iff X >= r, so b.lo branches exactly when X < r.
R_MOD = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
SAMPLE_CMP = "cmp x0, x4; sbcs xzr, x1, x5; sbcs xzr, x2, x6; sbcs xzr, x3, x7; b.lo accept"
# Heapsort level (proof.md 0-ii): ten words w0..w9 at offsets 0..72, w0 most significant (F, alpha,
# beta, high word first); x0, x1 the children, x13 the sifted temporary.
CHILD_CMP = "; ".join(["ldr x%d, [x0, #%d]" % (2 + k, 8 * k) for k in range(10)]
					  + ["ldr x%d, [x1, #%d]" % (19 + k, 8 * k) for k in range(10)] + ["cmp x11, x28"]
					  + ["sbcs xzr, x%d, x%d" % (2 + k, 19 + k) for k in range(8, -1, -1)] + ["csel x12, x1, x0, lo"]
					  + ["csel x%d, x%d, x%d, lo" % (2 + k, 19 + k, 2 + k) for k in range(10)])
PARENT_CMP = "; ".join(["ldr x%d, [x13, #%d]" % (19 + k, 8 * k) for k in range(10)] + ["cmp x28, x11"]
					   + ["sbcs xzr, x%d, x%d" % (19 + k, 2 + k) for k in range(8, -1, -1)] + ["b.hs sift_done"])

# Key combination of a multiply step (proof.md 5): x17 = sp + 304, key slot k (1-11) at x17 + 48k - 288;
# a = Re(c_1) at byte 96 of the walk state, whose address is at [sp, #8]; segment 1 after the third Fp
# multiplication (t1, t2 in slots 8, 9), segment 2 after the fourth (B, za, zB in 1, 2, 3); p_pool: the 1p entry.
X = lambda *n: ["x%d" % i for i in n]
P6, A6, T2, U2 = X(0, 1, 2, 3, 4, 5), X(6, 7, 8, 9, 10, 11), X(12, 13) * 3, X(6, 7) * 3
O = lambda k, i: "[x17, #%d]" % (48 * k - 288 + 8 * i)
SUB, ADD = "subs sbcs sbcs sbcs sbcs sbc", "adds adcs adcs adcs adcs adc"

def C(d, m, s, f=SUB, ld=None, st=None):
	"""6-limb carry chain d = m op s; ld, st: (registers, slot) loaded before even, stored after odd limbs."""
	o = []
	for i, g in enumerate(f.split()):
		o += ["ldp %s, %s, %s" % (ld[0][i], ld[0][i + 1], O(ld[1], i))] if ld and i % 2 == 0 else []
		o += ["%s %s, %s, %s" % (g, d[i], m[i], s[i])]
		o += ["stp %s, %s, %s" % (st[0][i - 1], st[0][i], O(st[1], i - 1))] if st and i % 2 else []
	return o
LD = lambda r, k, op="ldp": ["%s %s, %s, %s" % (op, r[i], r[i + 1], O(k, i)) for i in (0, 2, 4)]
HEAD = ["add x17, sp, #304", "adr x16, p_pool"] + ["ldp x%d, x%d, [x16, #%d]" % (i, i + 1, 8 * i) for i in (0, 2, 4)]

def SIDE(z, u, nx, nz, nu):
	"""x in x6..x11: -x, u = z - x (masked add of p), -u, -z into slots nx, u, nu, nz."""
	m = C(A6, T2, A6, SUB[:-3] + "sbcs", ld=(T2, z)) + ["sbc x12, xzr, xzr"]
	for i, g in enumerate(ADD.split()):
		m += ["and x13, x%d, x12" % i, "%s %s, %s, x13" % (g, A6[i], A6[i])]
	return C(T2, P6, A6, st=(T2, nx)) + m + LD(A6, u, "stp") + C(T2, P6, A6, st=(T2, nu)) + C(U2, P6, U2, ld=(U2, z), st=(U2, nz))
KEY_SEG1 = (HEAD + LD(A6, 8) + C(A6, A6, T2, ADD, ld=(T2, 9)) + C(P6, A6, P6, SUB[:-3] + "sbcs")
			+ ["csel %s, x%d, %s, hs" % (r, i, r) for i, r in enumerate(A6)] + LD(A6, 1, "stp"))
KEY_SEG2 = (HEAD + ["ldr x16, [sp, #8]"] + ["ldp x%d, x%d, [x16, #%d]" % (6 + i, 7 + i, 96 + 8 * i) for i in (0, 2, 4)]
			+ SIDE(2, 4, 6, 8, 10) + LD(A6, 1) + SIDE(3, 5, 7, 9, 11))
KEY_COMB_INSTR = len(KEY_SEG1) + len(KEY_SEG2)
# Full-width minimum of the 12 keys (proof.md 5), after segment 2: key 0 = a, key k in slot k; x0-x5 hold the
# minimum, x6-x11 the candidate, x12 its index. subs/sbcs xzr, cand - min: C = 0 iff cand < min on all six limbs.
KEY_MIN = ["add x17, sp, #304", "ldr x16, [sp, #8]"] + ["ldp x%d, x%d, [x16, #%d]" % (i, i + 1, 96 + 8 * i) for i in (0, 2, 4)]
for k in range(1, 12):
	KEY_MIN += (LD(A6, k) + C(["xzr"] * 6, A6, P6, SUB[:-3] + "sbcs") + ["csel %s, %s, %s, lo" % (m, c, m) for c, m in zip(A6, P6)]
				+ ["mov x13, #%d" % k, "csel x12, x13, %s, lo" % ("x12" if k > 1 else "xzr")])
KEY_MIN_INSTR = len(KEY_MIN)

def edge3():
	"""192-bit (m, l) pairs: equal top words with larger, equal and smaller low 128 bits; a borrow through
    every word; all words equal; all-ones; zero; one-word differences at each position."""
	T = 1 << 128
	out = []
	for top in (0, 1, 5, M64 - 1, M64):
		for lo_m, lo_l in ((7, 3), (3, 7), (5, 5), (T - 1, 0), (0, T - 1), (1 << 64, (1 << 64) - 1)):
			out.append((top * T + lo_m, top * T + lo_l))
	full = (1 << 192) - 1
	out += [(0, 0), (full, full), (full, 0), (0, full), (1 << 128, (1 << 128) - 1), ((1 << 128) - 1, 1 << 128),
			(0, 1), (1, 0), (1 << 64, 1), (1, 1 << 64)]
	for i in range(3):
		for d in (-1, 1):
			b = 1 << (64 * i)
			out += [((1 << 191) + d * b, 1 << 191), (1 << 191, (1 << 191) + d * b)]
	return out

def run(n_random=4000, seed=20261026):
	rng = random.Random(seed)
	res = {}

	def rnd192():
		k = rng.randrange(6)
		if k == 0:
			return rng.getrandbits(192)
		base = rng.getrandbits(192)
		return base, base ^ (rng.getrandbits(64 * rng.randrange(1, 3)) if k < 3 else rng.randrange(2))

	def regs3(m, l):
		r = {"x%d" % i: w for i, w in zip((0, 1, 2), words(m, 3))}
		r.update({"x%d" % i: w for i, w in zip((5, 6, 7), words(l, 3))})
		return r

	# 1. meter compare: branch to end_walk_phase exactly when meter >= L (192-bit unsigned)
	cases = edge3()
	for _ in range(n_random):
		x = rnd192()
		cases.append(x if isinstance(x, tuple) else (x, rng.getrandbits(192)))
	ok = sum(run_seq(METER_CMP, regs3(m, l))[1] == (m >= l) for m, l in cases)
	old_bad = [(m, l) for m, l in cases if run_seq(OLD_METER_CMP, regs3(m, l))[1] != (m > l)]
	verdict_case = ((3 << 128) + 9, (3 << 128) + 4)  # equal top words, larger low 128 bits
	vc_old = run_seq(OLD_METER_CMP, regs3(*verdict_case))[1]
	vc_new = run_seq(METER_CMP, regs3(*verdict_case))[1]
	assert ok == len(cases) and vc_new is True and vc_old is False and old_bad
	res["meter_compare"] = {"listing": METER_CMP, "predicate": "branch iff meter >= L (192-bit unsigned)",
							"cases": len(cases), "pass": True,
							"old_b.hi_listing": OLD_METER_CMP, "old_b.hi_claimed_predicate": "meter > L",
							"old_b.hi_wrong_cases": len(old_bad),
							"verdict_case_equal_top_words_larger_low_128_bits":
								{"old_b.hi_branches": vc_old, "b.hs_branches": vc_new}}

	# 2. meter add: x2:x1:x0 += x4 (charge < 2^25), exact while the sum stays below 2^192
	n = 0
	for m, _ in cases:
		for c in (0, 1, 22402, 10086, 4000, 600, (1 << 25) - 1):
			if m + c < 1 << 192:
				r = regs3(m, 0)
				r["x4"] = c
				out, _ = run_seq(METER_ADD, r)
				assert join([out["x0"], out["x1"], out["x2"]]) == m + c
				n += 1
	res["meter_add"] = {"listing": METER_ADD, "predicate": "meter' = meter + charge (192-bit, no wrap below 2^192)",
						"cases": n, "pass": True}

	# 3. stuck-walk counter: x9:x8 += 1, then branch to restart iff the new counter >= 5120 * 2^64
	sc = [0, 1, M64, M64 + 1, STUCK_CAP - 2, STUCK_CAP - 1, STUCK_CAP, STUCK_CAP + M64, (5119 << 64) + M64,
		  (5121 << 64), (1 << 127), (1 << 128) - 2]
	sc += [rng.getrandbits(rng.choice((13, 64, 77, 127))) for _ in range(n_random)]
	for v in sc:
		r, t = run_seq(STUCK, {"x8": v & M64, "x9": v >> 64, "x16": 0})
		assert join([r["x8"], r["x9"]]) == v + 1 and t == (v + 1 >= STUCK_CAP)
	res["stuck_walk_counter"] = {"listing": STUCK, "predicate": "counter' = counter + 1 (128-bit); branch iff "
								 "counter' >= 5120 * 2^64 = 20/theta", "cases": len(sc), "pass": True}

	# 4. one-word counters (report against D = 2^55; restart and pair counters use the same cmp/b.hs form)
	rc = [(0, 1 << 55), ((1 << 55) - 1, 1 << 55), (1 << 55, 1 << 55), ((1 << 55) + 1, 1 << 55), (M64, 0), (0, M64)]
	rc += [(rng.getrandbits(56), rng.getrandbits(56)) for _ in range(n_random)]
	for a, b in rc:
		assert run_seq(REPORT_CMP + " cap", {"x10": a, "x11": b})[1] == (a >= b)
	res["one_word_counter_compare"] = {"listing": REPORT_CMP, "predicate": "branch iff count >= cap (64-bit unsigned)",
									   "cases": len(rc), "pass": True}

	# 5. chunk counter: x12 -= 1, branch iff the result is zero (single-word subs, so Z is exact)
	cc = [1, 2, 0, M64, 1 << 63] + [rng.getrandbits(64) for _ in range(n_random)]
	for v in cc:
		r, t = run_seq(CHUNK, {"x12": v})
		assert r["x12"] == (v - 1) & M64 and t == (v == 1)
	res["chunk_counter"] = {"listing": CHUNK, "predicate": "x12' = x12 - 1; branch iff x12' == 0 (one word)",
							"cases": len(cc), "pass": True}
	# 5b. sampler compare: accept iff the four-word X < r
	def regs4(X, M):
		d = {"x%d" % i: w for i, w in zip((0, 1, 2, 3), words(X, 4))}
		d.update({"x%d" % i: w for i, w in zip((4, 5, 6, 7), words(M, 4))})
		return d
	sm = [(0, R_MOD), (R_MOD - 1, R_MOD), (R_MOD, R_MOD), (2**255 - 1, R_MOD), (0, 0), (1, R_MOD), (R_MOD, 1)]
	for i in range(4):
		b = 1 << (64 * i)
		sm += [(R_MOD - b, R_MOD), (R_MOD, R_MOD - b)]
	for _ in range(n_random):
		sm += [(rng.getrandbits(255), R_MOD), (rng.getrandbits(255), rng.getrandbits(255))]
	assert all(run_seq(SAMPLE_CMP, regs4(X, M))[1] == (X < M) for X, M in sm)
	res["sampler_compare"] = {"listing": SAMPLE_CMP, "predicate": "branch to accept iff X < r (255-bit, four words)",
							  "cases": len(sm), "pass": True}
	# 5c. branch-field compare (proof.md section 3): square iff (T mod 64) < 37.
	BRANCH_CMP = "and x0, x0, #63; cmp x0, #37; b.lo square"
	bc = list(range(128)) + [(1 << 381) - 1, (1 << 381) - 2, (1 << 381) - 64, 1 << 300]
	bc += [rng.getrandbits(381) for _ in range(n_random)]
	assert all(run_seq(BRANCH_CMP, {"x0": v})[1] == ((v & 63) < 37) for v in bc)
	res["branch_field_compare"] = {"listing": BRANCH_CMP, "predicate": "branch to square iff (T mod 64) < 37",
								   "cases": len(bc), "pass": True}
	# 6. heapsort level compares on 80-byte records (640-bit unsigned, w0 most significant)
	W = 1 << 640
	rec = [(0, 0), (W - 1, W - 1), (5, 4), (4, 5), (W - 1, W - 2), (W - 2, W - 1)]
	for hi in (0, 1 << 62, M64):
		a = (hi << 576) + (5 << 512)  # equal top words; w1 differs by 1 and w2..w9 are all-zero vs all-ones,
		rec += [(a, a - 1), (a - 1, a)]  # so the borrow runs through every word
	for i in range(10):
		b = 1 << (64 * i)
		rec += [(W // 2 + b, W // 2), (W // 2, W // 2 + b), (W // 2 - b, W // 2)]
	for _ in range(n_random):
		x = rng.getrandbits(640)
		k = rng.randrange(4)
		rec.append((x, rng.getrandbits(640) if k == 0 else x ^ rng.getrandbits(64 * rng.randrange(1, 10)) if k == 1
					else x ^ rng.randrange(2) if k == 2 else x + rng.choice((1, -1)) * (rng.getrandbits(64) + 1)))
	rec = [(a % W, b % W) for a, b in rec]
	nc = 0
	for left, right in rec:
		mem = {}
		for base, v in ((0x1000, left), (0x1050, right), (0x2000, right)):
			for k in range(10):
				mem[base + 8 * k] = (v >> (64 * (9 - k))) & M64
		r, t = run_seq(CHILD_CMP, {"x0": 0x1000, "x1": 0x1050}, mem)
		big = right if right > left else left
		assert t is None and r["x12"] == (0x1050 if right > left else 0x1000)
		assert sum(r["x%d" % (2 + k)] << (64 * (9 - k)) for k in range(10)) == big
		for tmp in (left, right):  # parent compare: temporary (at x13) against the larger child in x2..x11
			for k in range(10):
				mem[0x3000 + 8 * k] = (tmp >> (64 * (9 - k))) & M64
			_, t = run_seq(PARENT_CMP, dict(r, x13=0x3000), mem)
			assert t == (tmp >= big)
			nc += 1
	res["heapsort_child_compare"] = {"listing": CHILD_CMP, "instructions": CHILD_CMP.count(";") + 1,
									 "predicate": "x12 = right if right > left else left (640-bit unsigned, "
									 "a tie picks left); x2..x11 = that child's words", "cases": len(rec), "pass": True}
	res["heapsort_parent_compare"] = {"listing": PARENT_CMP, "instructions": PARENT_CMP.count(";") + 1,
									  "predicate": "branch to sift_done (no swap) iff temporary >= larger child; "
									  "fall through to the swap iff child > temporary", "cases": nc, "pass": True}
	assert CHILD_CMP.count(";") + 1 == 41 and PARENT_CMP.count(";") + 1 == 21
	# 7. key combination, both segments with the fourth multiplication's effect between them (x0-x17
	# overwritten, zB stored); slots against values mod p on stored limbs (x -> x 2^384 mod p is linear).
	PP = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
	ZT = next(h for h in (pow(g, PP // 6, PP) for g in range(2, 50)) if (h * h - h + 1) % PP == 0)
	SP, WS, POOL = 1 << 16, 1 << 18, 1 << 19
	for s in KEY_SEG1 + KEY_SEG2 + KEY_MIN:  # immediates: ldp/stp imm7 x 8, ldr imm12 x 8, add imm12, mov
		n = int(s.rpartition("#")[2].rstrip("]")) if "#" in s else 0
		assert (n % 8 == 0 or s[:3] == "mov") and (-512 <= n <= 504 if s[:3] in ("ldp", "stp") else 0 <= n < 4096)
	S = lambda k: SP + 16 + 48 * k
	get = lambda m, k: join([m[S(k) + 8 * i] for i in range(6)])
	E = [0, 1, 2, PP - 1, PP - 2, PP // 2, PP // 2 + 1, 1 << 64, (1 << 320) - 1, PP - (1 << 64)]
	cs = [(x, y, v, rng.choice(E + [x])) for x in E for y in E for v in (0, PP - 1, (PP - y) % PP)]
	cs += [tuple(rng.randrange(PP) for _ in "abcd") for _ in range(300)]
	live, nt = X(14, 15, *range(18, 31)) + ["sp"], 0
	for tied in (1, 0):
		for a, t1, t2, w in cs:
			m = {SP + 8 * i: rng.getrandbits(64) for i in range(80)}
			m.update({WS + 8 * i: rng.getrandbits(64) for i in range(256)})
			put = lambda ad, v: m.update({ad + 8 * i: x for i, x in enumerate(words(v, 6))})
			put(WS + 96, a), put(S(8), t1), put(S(9), t2), put(POOL, PP), put(S(2), ZT * a % PP if tied else w)
			m[SP + 8], r, ch = WS, {"x%d" % i: rng.getrandbits(64) for i in range(31)}, []
			r.update(sp=SP, p_pool=POOL)
			for seg in (KEY_SEG1, KEY_SEG2):
				if ch:
					r.update({"x%d" % i: rng.getrandbits(64) for i in range(18)})
					put(S(3), ZT * get(m, 1) % PP if tied else rng.choice((0, get(m, 1), w)))
				m0, r0 = dict(m), dict(r)
				r = run_seq("; ".join(seg), r, m)[0]
				assert all(r[k] == r0[k] for k in live)
				ch.append({k for k in m if m[k] != m0[k]})
			assert get(m, 1) == (t1 + t2) % PP and ch[0] <= {S(1) + 8 * i for i in range(6)}
			assert ch[1] <= {S(k) + 8 * i for k in range(4, 12) for i in range(6)}
			for x, z, k in ((a, get(m, 2), 0), (get(m, 1), get(m, 3), 1)):
				odd = 0 in (x, z, (z - x) % PP)
				for j, v in ((4, z - x), (6, -x), (8, -z), (10, x - z)):
					assert get(m, k + j) in (v % PP, PP if odd and v % PP == 0 else -1)
				assert not (tied and odd) or x == z == 0  # keys x and z both 0: a charged tie
				nt += tied and odd
	res["key_combination"] = {"segment_1": KEY_SEG1, "segment_2": KEY_SEG2, "instructions": KEY_COMB_INSTR,
							  "cases": 2 * len(cs), "tie_cases": nt, "pass": True}
	assert KEY_COMB_INSTR == 154
	# 8. key minimum: full width, the lowest index on an exact tie, reads only a, [sp, #8] and slots 1-11, no store
	src = {WS + 96 + 8 * i for i in range(6)} | {S(k) + 8 * i for k in range(1, 12) for i in range(6)} | {SP + 8}

	class Rd(dict):
		def __getitem__(s, k):
			s.bad += [] if k in src else [k]
			return dict.__getitem__(s, k)
	km = [list(range(12)), list(range(11, -1, -1)), [PP - 1] * 12, [0] * 12, [5] * 6 + [4] * 6, [PP - 1] * 11 + [0]]
	for d in range(7):  # keys equal above limb d, a near pair min +- 2^(64e) (e <= d), ties; d = 6: random
		for _ in range(30 if d < 6 else 150):
			b = rng.randrange(PP // 2) >> (64 * d + 64) << (64 * d + 64)
			ks = [b + rng.getrandbits(64 * d + 64) if d < 5 else rng.randrange(PP) for _ in range(12)]
			v = min(ks) + rng.choice((1, -1)) * (1 << 64 * rng.randrange(min(d, 5) + 1))
			ks[rng.randrange(12)] = v if 0 <= v < PP else ks[0]
			ks[rng.randrange(12)] = rng.choice(ks)
			km.append(ks)

	def kmin(L, ks):
		m = Rd({SP + 8: WS})
		m.bad = []
		for k, v in enumerate(ks):
			m.update({(S(k) if k else WS + 96) + 8 * i: w for i, w in enumerate(words(v, 6))})
		r0 = {"x%d" % i: rng.getrandbits(64) for i in range(31)}
		r0.update(sp=SP)
		m0, r = dict(m), run_seq("; ".join(L), r0, m)[0]
		assert dict(m) == m0 and not m.bad and all(r[k] == r0[k] for k in r0 if k not in X(*range(14), 16, 17))
		return r["x12"] == ks.index(min(ks)) and join([r[x] for x in P6]) == min(ks)
	mut = {"top_limb_only": [s.replace("sbcs", "subs") for s in KEY_MIN if s[:4] not in ("subs", "sbcs") or "x11, x5" in s],
		   "z_flag": [s.replace(", lo", ", ls") for s in KEY_MIN], "dropped_csel": KEY_MIN[:100] + KEY_MIN[101:]}
	assert all(kmin(KEY_MIN, ks) for ks in km) and all(not all(kmin(L, ks) for ks in km) for L in mut.values())
	res["key_minimum"] = {"listing": KEY_MIN, "instructions": KEY_MIN_INSTR, "cases": len(km), "pass": True,
						  "failing_mutants": sorted(mut)}
	assert KEY_MIN_INSTR == 192
	res["all_pass"] = True
	return res

if __name__ == "__main__":
	print(json.dumps(run(), indent=1))
