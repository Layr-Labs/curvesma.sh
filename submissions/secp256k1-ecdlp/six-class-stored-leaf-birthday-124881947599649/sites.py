"""Site checks for experiment.py (proof.md 3.4): each site type and fixup of leaf.s alone on built inputs, and F's
field routines on weak and edge inputs (helpers: E)."""
import math
import os

SP = 0x7000000      # stack pointer of a site test; its frame's fixup word is at SP - 16

SITE_KERNELS = """
\t.macro ENTER
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstp xzr,x30,[sp,#-16]!
\t.endm
\t.macro LEAVE
\tldp xzr,x30,[sp],#16
\tret
\t.endm
t_shx:
\tENTER
\tSHX x5,x25,x28,x26,x27
\tLEAVE
t_shxt:
\tENTER
\tSHXT x25,x28,x26,x27
\tLEAVE
t_fdx:
\tENTER
\tFDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
\tLEAVE
t_cpf:
\tENTER
\tAD4S x1,128,144,x9,x10,x13,x14,x15,x16,x12,x18
\tFDX x13,x14,x12,x18,x13,x14,x19,x12,x15,x16
\tLEAVE
t_squ:
\tENTER
\tSQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
\tSQF x21,x22,x23,x24,x5,x6,x7,x8,x9
\tLEAVE
t_sqf:
\tENTER
\tSQF x21,x22,x23,x24,x5,x6,x7,x8,x9
\tLEAVE
t_sqv:
\tENTER
\tmov x30,x1
\tSQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
\tSQF x25,x26,x0,x2,x5,x6,x7,x8,x9
\tLEAVE
t_fs:
\tENTER
\tFS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
\tLEAVE
t_faddw:
\tENTER
\tFADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
\tLEAVE
"""
# a test's count without a fixup, beyond ENTER 3 and LEAVE 2 (AD4S 8; SQ13 41, SQF 20; t_sqv's MOV of x30), and what
# a fixup adds (proof.md 3.4)
COMMON = dict(t_shx=4, t_shxt=3, t_fdx=6, t_cpf=8 + 6, t_squ=41 + 20, t_sqv=1 + 41 + 20, t_sqf=20, t_fs=7, t_faddw=7)
EXTRA = dict(t_shx=3, t_shxt=3, t_fdx=3, t_cpf=3, t_squ=6, t_sqv=6, t_sqf=6, t_fs=4, t_faddw=4)


def run_site(E, prog, labels, entry, regs, mem=None):
    """Run one site test with a sentinel x30. returns: (registers, count, fixup word)."""
    mem = {} if mem is None else mem
    out = {}
    n = E.a64.a64_run(prog, labels, entry, dict(regs, sp=SP, x30=E.SENTINEL), mem, out_regs=out)
    return out, n, mem[SP - 16]


def site_checks(E, rng, fprog, flabels):
    """Each site type of 3.4 deviates exactly when its predicate holds and matches its model; a fixup runs
    exactly on its carry, adds C0 and its instructions, its output in 3.4's set; FS, FADDW, fe_mul_f and fe_sqr_f
    are exact. returns: dict."""
    p, C0, M64, put, get = E.p, E.C0, E.M64, E.put, E.get
    reg_int, as_regs = E.reg_int, E.as_regs
    prog, labels = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read() + SITE_KERNELS)
    res, fix = {}, {}

    def site(entry, regs, mem=None):
        """Run a site test. returns: (registers, count, fixup word)."""
        return run_site(E, prog, labels, entry, regs, mem)

    def tally(name, cases):
        """Record (predicate, deviates, matches model) triples. name: site; cases: list. returns: None."""
        res[name] = dict(cases=len(cases), predicate_holds=sum(c[0] for c in cases),
                         deviates_iff_predicate=all(c[0] == c[1] for c in cases), model=all(c[2] for c in cases))

    def fixed(entry, f, n, word, in_set):
        """Record one fixup observation. entry: test; f: the fixup predicate; n: count; word: fixup word; in_set: the
        output (or row) in 3.4's word set when f. returns: None."""
        d = fix.setdefault(entry, dict(cases=0, fixups=0, word_ok=True, count_ok=True, word_set_ok=True))
        d["cases"] += 1
        d["fixups"] += f
        d["word_ok"] &= word == C0 * f
        d["count_ok"] &= n == 5 + COMMON[entry] + EXTRA[entry] * f
        d["word_set_ok"] &= in_set or not f

    def fold(v):
        """One fold of a 5- or 6-word value: (v mod 2^256) + floor(v / 2^256) C0, unreduced. returns: integer."""
        return (v % 2 ** 256) + (v >> 256) * C0

    def low_carry(v, k):
        """Whether the fold of v carries out of its low k words. returns: bool."""
        return (v % 2 ** (64 * k)) + (v >> 256) * C0 >= 2 ** (64 * k)

    cases = []
    for kind in (0, 1, 2, 3):       # built to carry through, just below the carry, random, only the first ADDS carries
        for _ in range(48):
            b3 = rng.randrange(2 ** 40, 2 ** 64)
            hi = (b3 * C0) >> 64
            low = (2 ** 192 - hi + rng.randrange(hi), 2 ** 192 - hi - 1 - rng.randrange(2 ** 20), rng.getrandbits(192),
                   rng.getrandbits(126) << 64 | (2 ** 64 - hi + rng.randrange(hi)))[kind]
            B = low + (b3 << 192)
            o, n, word = site("t_shx", as_regs(B, ("x25", "x26", "x27", "x28")))
            o = reg_int(o, ("x5", "x25", "x26", "x27"))
            want = B * 2 ** 64 - b3 * p
            cases.append((low + hi >= 2 ** 192, o != want, o == want % 2 ** 256))
            fixed("t_shx", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
            o, n, word = site("t_shxt", as_regs(B, ("x25", "x26", "x27", "x28")))   # SHXT: words 1..3 of SHX's
            o = reg_int(o, ("x25", "x26", "x27"))
            cases.append((low + hi >= 2 ** 192, o != want >> 64, o == want % 2 ** 256 >> 64))
            fixed("t_shxt", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
    tally("row_shift", cases)
    for name, entry, top in (("fold", "t_fdx", 2 ** 66), ("fold_of_a_constant_plus_tabled_rows", "t_cpf", 2 ** 66 + 1)):
        cases = []
        for kind in (0, 1, 2, 3):   # built to carry through, just below, random, carrying only into limb 2
            for _ in range(48):
                u = rng.randrange(1 + (entry == "t_cpf"), top + 1)
                m, mem = u * C0, {}
                low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256),
                       rng.getrandbits(127) << 128 | (2 ** 128 - m + rng.randrange(m)))[kind]
                Z = low + (u << 256)
                if entry == "t_fdx":
                    regs = as_regs(low, ("x10", "x11", "x12", "x13"))
                    regs.update(x7=u & M64, x24=u >> 64)
                    outs = ("x10", "x11", "x12", "x13")
                else:       # AD4S adds a constant below p, from memory, to the product's six words, as for E, F
                    K = rng.randrange(p)
                    regs = as_regs(Z - K, ("x13", "x14", "x15", "x16", "x12", "x18"))
                    regs["x1"] = 0x10000
                    put(mem, 0x10000 + 128, K)
                    outs = ("x13", "x14", "x15", "x16")
                o, n, word = site(entry, regs, mem)
                o = reg_int(o, outs)
                t = fold(Z)
                cases.append((t >= 2 ** 256, o != t, o == t % 2 ** 256))
                fixed(entry, low_carry(Z, 2), n, word, t % 2 ** 128 < u * C0 and t >= 2 ** 128)
        tally(name, cases)
    for name, entry, a_regs, o_regs in (
            ("square_fold_u", "t_squ", ("x17", "x18", "x19", "x20"), ("x21", "x22", "x23", "x24")),
            ("square_fold_v", "t_sqv", ("x27", "x28", "x29", "x1"), ("x25", "x26", "x0", "x2"))):
        cases = []
        for kind in (0, 1, 2, 3):   # a^2 just below k 2^448 (step 1 drops), 4 (k-1) C0 2^192 below it, residue in
            for _ in range(32):     # (2^65, 2^80) (step 2's fixup), random
                k = rng.randrange(2, 2 ** 64)
                if kind < 2:
                    a = math.isqrt(k * 2 ** 448 - 1 - kind * 4 * (k - 1) * C0 * 2 ** 192)
                else:
                    while True:
                        T = (rng.randrange(2 ** 65, 2 ** 80), rng.randrange(p))[kind - 2]
                        a = E.sqrt_mod_p(T)
                        if a is not None:
                            break
                    a += rng.choice((0, p)) if a < 2 ** 256 - p else 0
                o, n, word = site(entry, as_regs(a, a_regs))
                o = reg_int(o, o_regs)
                t, lost, f = E.sqf(a * a)
                cases.append((lost, (o - a * a) % p != 0, o == t))
                fixed(entry, f, n, word, o < C0 * 2 ** 192)
        tally(name, cases)
    cases = []
    for kind in (0, 1, 2, 3):       # SQF alone on built words: step 1 drops, just below, step 2 carries, random
        for _ in range(48):
            c7 = rng.randrange(1, 2 ** 64)
            m = c7 * C0
            if kind == 2:               # c7 = 0, words 0..3 within h C0 of 2^256, h < 2^192
                h = rng.getrandbits(192) | 1
                S = (h << 256) + 2 ** 256 - h * C0 + rng.randrange(h * C0)
            else:
                mid = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), 0,
                       rng.getrandbits(256))[kind]
                S = (c7 << 448) + (mid << 192) + rng.getrandbits(192)
            o, n, word = site("t_sqf", as_regs(S, ("x21", "x22", "x23", "x24", "x5", "x6", "x7", "x8")))
            o = reg_int(o, ("x21", "x22", "x23", "x24"))
            t, lost, f = E.sqf(S)
            cases.append((lost, (o - S) % p != 0, o == t))
            fixed("t_sqf", f, n, word, o < C0 * 2 ** 192)
    tally("square_fold_built", cases)
    edge = [0, 1, p - 1, C0, 2 ** 255]
    pairs = [(a, b) for a in edge for b in edge] + [(rng.randrange(p), rng.randrange(p)) for _ in range(48)]
    for _ in range(24):             # canonical inputs that take FS's borrow fixup and FADDW's carry fixup
        a = rng.randrange(p >> 1)
        pairs.append((a, a + (rng.randrange(1, 2 ** 190) << 64) - rng.randrange(C0)))
        b = rng.randrange(p >> 1, p)
        pairs.append((2 ** 256 - b + (rng.randrange(1, 2 ** 120) << 64) + 2 ** 64 - 1 - rng.randrange(C0), b))
    fs = fa = True
    for a, b in pairs:
        if a >= p or b >= p:
            continue
        regs = as_regs(a, ("x5", "x6", "x7", "x8"))
        regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
        o, n, word = site("t_fs", regs)
        o = reg_int(o, ("x9", "x10", "x11", "x12"))
        fs &= o == (a - b) % p
        fixed("t_fs", a < b and (a - b) % 2 ** 64 < C0, n, word, o % 2 ** 64 >= 2 ** 64 - C0)
        regs = as_regs(a, ("x13", "x14", "x15", "x16"))
        regs.update(as_regs(b, ("x17", "x18", "x19", "x20")))
        o, n, word = site("t_faddw", regs)
        o = reg_int(o, ("x13", "x14", "x15", "x16"))
        fa &= o % p == (a + b) % p and o < 2 ** 256
        fixed("t_faddw", a + b >= 2 ** 256 and (a + b) % 2 ** 64 + C0 >= 2 ** 64, n, word, o % 2 ** 64 < C0)
    res["fs_canonical_exact"] = dict(cases=len(pairs), exact=fs)
    res["faddw_canonical_exact"] = dict(cases=len(pairs), exact=fa)
    res["fixups"] = fix
    res["fixups_ok"] = all(d["word_ok"] and d["count_ok"] and d["word_set_ok"] and d["fixups"] > 0 for d in fix.values())
    weak = [0, 1, C0, p - 1, p, p + 1, 2 ** 256 - 1] + [rng.getrandbits(256) for _ in range(12)]
    mul_in = [(a, b) for a in weak[:7] for b in weak[:7]] + [(a, rng.getrandbits(256)) for a in weak[7:]]
    for _ in range(12):             # products with a small residue, which take the final carry of the reduction
        a = rng.randrange(1, p)
        mul_in.append((a, rng.randrange(1, 2 ** 64) * pow(a, -1, p) % p))
    sqr_in = weak + [p - rng.randrange(1, 2 ** 100) for _ in range(12)] + [2 ** 256 - rng.randrange(1, 2 ** 40)
                                                                         for _ in range(6)]
    while len(sqr_in) < len(weak) + 42:     # square roots of small residues reach the final +C0 of the reduction
        s_ = rng.randrange(1, 2 ** 64)
        if pow(s_, (p - 1) // 2, p) == 1:
            sqr_in.append(pow(s_, (p + 1) // 4, p))
    ok_m = ok_s = True
    for a, b in mul_in:
        mem = {}
        put(mem, 0x1000, a)
        put(mem, 0x1100, b)
        E.a64.a64_run(fprog, flabels, "fe_mul_f", dict(x0=0x2000, x1=0x1000, x2=0x1100, sp=0x100000, x30=E.SENTINEL), mem)
        ok_m &= get(mem, 0x2000) % p == a * b % p
    for a in sqr_in:
        mem = {}
        put(mem, 0x1000, a)
        E.a64.a64_run(fprog, flabels, "fe_sqr_f", dict(x0=0x2000, x1=0x1000, sp=0x100000, x30=E.SENTINEL), mem)
        ok_s &= get(mem, 0x2000) % p == a * a % p

    def final_carry(v):
        """Whether reducing the 512-bit v takes the routines' final +C0 (the second fold reaches 2^256). returns: bool."""
        t = (v % 2 ** 256) + (v >> 256) * C0
        return (t % 2 ** 256) + (t >> 256) * C0 >= 2 ** 256
    res["fe_mul_f_exact"] = dict(cases=len(mul_in), exact=ok_m, final_carry_cases=sum(final_carry(a * b) for a, b in mul_in))
    res["fe_sqr_f_exact"] = dict(cases=len(sqr_in), exact=ok_s, final_carry_cases=sum(final_carry(a * a) for a in sqr_in))
    res["all_ok"] = res["fixups_ok"] and all(v.get("deviates_iff_predicate", True) and v.get("model", True) and v.get("exact", True)
                                    and v.get("predicate_holds", 1) > 0 for k, v in res.items()
                                    if isinstance(v, dict) and k != "fixups")
    return res
