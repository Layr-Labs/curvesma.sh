"""Site checks for experiment.py (proof.md 3.3 and 3.4): each arithmetic site type, rare-carry fixup and key path of
kernels/leaf.s runs alone on built inputs against exact arithmetic and its bit-level model, and F's field routines run
on weak and edge inputs. The shared helpers and constants come from experiment.py, passed in as the module E."""
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
t_fdx:
\tENTER
\tFDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
\tLEAVE
t_squ:
\tENTER
\tSQ5 x14,x15,x16,x17,x5,x6,x7,x8,x9,x0,x18,x29,x30,x19,x20,x21,x22,x23,x24
\tF5 x8,x5,x6,x9,x0,x19,x20,x18,x29
\tLEAVE
t_sqv:
\tENTER
\tSQ5 x10,x11,x12,x13,x5,x21,x7,x8,x22,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
\tF5 x8,x5,x21,x22,x14,x10,x11,x15,x16
\tLEAVE
t_sqc:
\tENTER
\tSQ5 x10,x11,x12,x13,x5,x6,x7,x8,x9,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
\tldp x15,x16,[x2,#192]
\tldp x17,x18,[x2,#208]
\tAC4 x8,x5,x6,x9,x14,x15,x16,x17,x18
\tF5 x8,x5,x6,x9,x14,x10,x11,x15,x16
\tLEAVE
t_sqd:
\tENTER
\tSQ5 x19,x20,x21,x22,x5,x12,x7,x8,x13,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
\tldp x15,x16,[x2,#224]
\tldp x17,x18,[x2,#240]
\tAC4 x8,x5,x12,x13,x14,x15,x16,x17,x18
\tF5 x8,x5,x12,x13,x14,x19,x20,x15,x16
\tLEAVE
t_sqt:
\tENTER
\tSQ5 x10,x11,x6,x9,x5,x12,x7,x8,x13,x14,x15,x16,x17,x23,x24,x0,x18,x29,x30
\tF5 x8,x5,x12,x13,x14,x10,x11,x15,x16
\tLEAVE
t_ltf:
\tENTER
\tSB6 x14,x15,x16,x17,x5,x29,x6,x7,x10,x11
\tFDX x14,x15,x5,x29,x14,x15,x18,x19,x16,x17
\tLEAVE
t_kj1:
\tstp xzr,x30,[sp,#-16]!
\tmov x16,#1
\tKEY
\tJ1 7
e7:
\tmov x16,#0
kd7:
\tLEAVE
t_kj2:
\tstp xzr,x30,[sp,#-16]!
\tmov x16,#1
\tKEY
\tJ2 8
e8:
\tmov x16,#0
kd8:
\tLEAVE
t_fs:
\tENTER
\tFS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
\tLEAVE
t_faddw:
\tENTER
\tFADDW x23,x24,x25,x26,x5,x6,x7,x8,x14,x15,x16,x17,x9
\tLEAVE
"""
# instructions of a site test without a fixup: ENTER 3 and LEAVE 2 around the macro's common path (SQ5 58, AC4 5 and
# its two LDPs where present; SB6 6), and the instructions a fixup adds (proof.md 3.4)
COMMON = dict(t_shx=4, t_fdx=6, t_squ=58 + 5, t_sqv=58 + 5, t_sqc=58 + 2 + 5 + 5, t_sqd=58 + 2 + 5 + 5, t_sqt=58 + 5,
              t_ltf=6 + 6, t_fs=7, t_faddw=7)
EXTRA = dict(t_shx=3, t_fdx=3, t_squ=3, t_sqv=3, t_sqc=3, t_sqd=3, t_sqt=3, t_ltf=3, t_fs=4, t_faddw=4)


def run_site(E, prog, labels, entry, regs, mem=None):
    """Run one site-test entry with a sentinel x30. E: experiment module; prog, labels; entry; regs; mem or None.
    returns: (final registers, executed instructions, the frame's fixup word)."""
    mem = {} if mem is None else mem
    out = {}
    n = E.a64.a64_run(prog, labels, entry, dict(regs, sp=SP, x30=E.SENTINEL), mem, out_regs=out)
    return out, n, mem[SP - 16]


def site_checks(E, rng, fprog, flabels):
    """Each site type of proof.md 3.4 and F's field routines against exact arithmetic: inputs built to meet a site's
    predicate must deviate, boundary and random ones must not, and each must match the bit-level model; a fixup runs
    exactly when its carry or borrow occurs, adds C0 to the fixup word and its instructions to the count, and its output
    (or the row it shifts) lies in the word set of 3.4; each key path matches its bit-exact model, adds its length to the
    meter and drops a leaf exactly when its window holds; FS,
    FADDW, fe_mul_f and fe_sqr_f are exact. E: experiment module; rng: Random; fprog, flabels: loaded secp_fe.s.
    returns: dict.
    """
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
    tally("row_shift", cases)
    cases = []
    for kind in (0, 1, 2, 3):       # the fold: built to carry through, just below, random, carrying only into limb 2
        for _ in range(48):
            u = rng.randrange(1, 2 ** 66)
            m = u * C0
            low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256),
                   rng.getrandbits(127) << 128 | (2 ** 128 - m + rng.randrange(m)))[kind]
            Z = low + (u << 256)
            regs = as_regs(low, ("x10", "x11", "x12", "x13"))
            regs.update(x7=u & M64, x24=u >> 64)
            o, n, word = site("t_fdx", regs)
            o = reg_int(o, ("x10", "x11", "x12", "x13"))
            t = fold(Z)
            cases.append((t >= 2 ** 256, o != t, o == t % 2 ** 256))
            fixed("t_fdx", low_carry(Z, 2), n, word, t % 2 ** 128 < u * C0 and t >= 2 ** 128)
    tally("fold", cases)
    cases = []
    for kind in (0, 1, 2, 3):       # the fused fold: the same four kinds
        for _ in range(48):
            u = rng.randrange(1, 2 ** 66 + 2)
            m = u * C0
            low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256),
                   rng.getrandbits(127) << 128 | (2 ** 128 - m + rng.randrange(m)))[kind]
            sub = rng.getrandbits(256)
            Z = low + (u << 256) + sub - 2 ** 256      # the accumulator before SB6, which subtracts sub and adds 2^256
            regs = as_regs(Z, ("x14", "x15", "x16", "x17", "x5", "x29"))
            regs.update(as_regs(sub, ("x6", "x7", "x10", "x11")))
            o, n, word = site("t_ltf", regs)
            o = reg_int(o, ("x14", "x15", "x16", "x17"))
            t = fold(low + (u << 256))
            cases.append((t >= 2 ** 256, o != t, o == t % 2 ** 256))
            fixed("t_ltf", low_carry(low + (u << 256), 2), n, word, t % 2 ** 128 < u * C0 and t >= 2 ** 128)
    tally("fused_fold", cases)
    Hb = 0x100000
    for name, entry, a_regs, o_regs, off in (
            ("square_fold_u", "t_squ", ("x14", "x15", "x16", "x17"), ("x19", "x20", "x6", "x9"), None),
            ("square_fold_v", "t_sqv", ("x10", "x11", "x12", "x13"), ("x10", "x11", "x21", "x22"), None),
            ("square_fold_q_t", "t_sqt", ("x10", "x11", "x6", "x9"), ("x10", "x11", "x12", "x13"), None),
            ("square_plus_constant_q_0", "t_sqc", ("x10", "x11", "x12", "x13"), ("x10", "x11", "x6", "x9"), 192),
            ("square_plus_constant_q_1", "t_sqd", ("x19", "x20", "x21", "x22"), ("x19", "x20", "x12", "x13"), 224)):
        cases = []
        for kind in (0, 1, 2):      # true value below 2^64, in (2^65, 2^80), random
            for _ in range(32):
                while True:
                    T = (rng.randrange(2 ** 64), rng.randrange(2 ** 65, 2 ** 80), rng.randrange(p))[kind]
                    if off is None:
                        a = E.sqrt_mod_p(T)
                        if a is None:
                            continue
                        a += rng.choice((0, p)) if a < 2 ** 256 - p else 0
                        c = 0
                    else:
                        a = rng.getrandbits(256)
                        c = (T - a * a) % p
                    break
                mem = {}
                put(mem, Hb + (off or 192), c)
                regs = as_regs(a, a_regs)
                regs.update(x2=Hb)
                o, n, word = site(entry, regs, mem)
                o = reg_int(o, o_regs)
                V = fold(a * a) + c
                t = fold(V)
                cases.append((t >= p, o != T, o == t % 2 ** 256))
                fixed(entry, low_carry(V, 2), n, word, t % 2 ** 128 < 2 ** 65)
        tally(name, cases)
    cases = []
    for kind in (0, 1, 2, 3):       # FS(a, b), b a fold output: b - a > p, b >= p without it, b random, borrow fixup
        for _ in range(32):
            if kind < 3:
                b = (rng.randrange(p, 2 ** 256), rng.randrange(p, 2 ** 256), rng.randrange(p))[kind]
                a = (rng.randrange(b - p) if b > p else 0, rng.randrange(b - p, p), rng.randrange(p))[kind]
            else:
                a = rng.randrange(p >> 1)
                b = a + (rng.randrange(1, 2 ** 190) << 64) - rng.randrange(C0)
            regs = as_regs(a, ("x5", "x6", "x7", "x8"))
            regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
            o, n, word = site("t_fs", regs)
            o = reg_int(o, ("x9", "x10", "x11", "x12"))
            d = (a - b) % 2 ** 256
            mdl = (d - C0) % 2 ** 256 if a < b else d
            cases.append((b - a > p, o % p != (a - b) % p, o == mdl))
            fixed("t_fs", a < b and d % 2 ** 64 < C0, n, word, mdl % 2 ** 64 >= 2 ** 64 - C0)
    tally("fs_on_fold_output", cases)
    kp = dict(cases=0, model=True, dropped_iff_window=True, counts_ok=True, meter_ok=True, class_key_ok=True,
              window_drops=[0, 0], exact_index_members=0)
    regs0 = dict(x7=E.GL[1], x8=E.GL[2], x9=E.GL[3], x10=E.GL[4], x11=E.GL[5], x14=E.GG[0], x15=E.GG[1], x18=0)
    for kind in range(3):           # random x, members of random classes, members of window-edge classes
        for _ in range(40):
            v = None
            if kind == 0:
                x = rng.randrange(1, p)
            else:
                v = E.edge_class(rng) if kind == 2 else E.exact_key(rng.randrange(1, p))[0]
                x = v * pow(E.BETA, rng.choice((1, 2)), p) % p      # a member of index 2 or 1
            for jj, entry, length in ((1, "t_kj1", 43), (2, "t_kj2", 47)):
                o, n, _ = site(entry, dict(regs0, **as_regs(x, ("x20", "x21", "x22", "x23"))))
                dropped = o["x16"] == 1
                a = reg_int(o, ("x21", "x22", "x23"))
                mdl = E.key_path(x, jj)
                kp["cases"] += 1
                kp["window_drops"][jj - 1] += dropped
                kp["dropped_iff_window"] &= dropped == (mdl is None)
                kp["model"] &= dropped or a == mdl
                kp["meter_ok"] &= o["x18"] == length
                kp["counts_ok"] &= n == 2 + 16 + length + (0 if dropped else 1) + 2 - (1 if dropped else 0) \
                    + (1 if jj == 2 else 0)
                if v is not None and E.exact_key(x)[1] == jj and E.rule_j(x) == jj:
                    kp["exact_index_members"] += 1
                    kp["class_key_ok"] &= dropped or a >> 16 == v >> 80
    res["key_paths"] = kp
    edge = [0, 1, p - 1, C0, 2 ** 255]
    pairs = [(a, b) for a in edge for b in edge] + [(rng.randrange(p), rng.randrange(p)) for _ in range(48)]
    for _ in range(24):             # canonical inputs that take FS's borrow fixup and FADDW's carry fixup
        a = rng.randrange(p >> 1)
        pairs.append((a, a + (rng.randrange(1, 2 ** 190) << 64) - rng.randrange(C0)))
        b = rng.randrange(p >> 1, p)
        pairs.append((2 ** 256 - b + (rng.randrange(1, 2 ** 120) << 64) + 2 ** 64 - 1 - rng.randrange(C0), b))
    weak_b = [p, p + 1, 2 ** 256 - 1, 2 ** 256 - C0] + [rng.randrange(p - C0, p) + C0 for _ in range(12)]
    fs = fa = fw = True
    for a, b in pairs:
        if a >= p:
            continue
        regs = as_regs(a, ("x5", "x6", "x7", "x8"))
        regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
        o, n, word = site("t_fs", regs)
        o = reg_int(o, ("x9", "x10", "x11", "x12"))
        fs &= o == (a - b) % p
        if b < p:
            fixed("t_fs", a < b and (a - b) % 2 ** 64 < C0, n, word, o % 2 ** 64 >= 2 ** 64 - C0)
        regs = as_regs(a, ("x5", "x6", "x7", "x8"))
        regs.update(as_regs(b, ("x14", "x15", "x16", "x17")))
        o, n, word = site("t_faddw", regs)
        o = reg_int(o, ("x23", "x24", "x25", "x26"))
        fa &= o % p == (a + b) % p
        if b < p:
            fixed("t_faddw", a + b >= 2 ** 256 and (a + b) % 2 ** 64 + C0 >= 2 ** 64, n, word, o % 2 ** 64 < C0)
    for a in edge + [rng.randrange(p) for _ in range(8)]:
        for b in weak_b:
            regs = as_regs(a, ("x5", "x6", "x7", "x8"))
            regs.update(as_regs(b, ("x14", "x15", "x16", "x17")))
            o, n, word = site("t_faddw", regs)
            fw &= reg_int(o, ("x23", "x24", "x25", "x26")) % p == (a + b) % p
            fixed("t_faddw", a + b >= 2 ** 256 and (a + b) % 2 ** 64 + C0 >= 2 ** 64, n, word,
                  reg_int(o, ("x23", "x24", "x25", "x26")) % 2 ** 64 < C0)
    res["fs_canonical_exact"] = dict(cases=len(pairs), exact=fs)
    res["faddw_exact"] = dict(cases=len(pairs), exact=fa)
    res["faddw_canonical_plus_weak_exact"] = dict(cases=len(weak_b) * (len(edge) + 8), exact=fw)
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
    kp = res["key_paths"]
    res["all_ok"] = res["fixups_ok"] and kp["dropped_iff_window"] and kp["counts_ok"] and kp["meter_ok"] \
        and kp["class_key_ok"] and min(kp["window_drops"]) > 0 and kp["exact_index_members"] > 0 and all(v.get("deviates_iff_predicate", True) and v.get("model", True) and v.get("exact", True)
                                    and v.get("predicate_holds", 1) > 0 for k, v in res.items()
                                    if isinstance(v, dict) and k != "fixups")
    return res
