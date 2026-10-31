"""Instruction-level check of the multiword compares and increments of proof.md 0-i and 0-iv.

A small AArch64 interpreter (64-bit X registers, NZCV) runs each listed sequence exactly as proof.md
writes it and asserts that it implements its stated predicate, on edge cases and on random cases.
Flags follow the architecture's AddWithCarry pseudocode:
  ADDS/ADCS: AddWithCarry(a, b, carry_in), carry_in 0 for ADDS, C for ADCS (ADC: same sum, no flags);
  SUBS/SBCS/CMP: AddWithCarry(a, NOT b, carry_in), carry_in 1 for SUBS/CMP, C for SBCS;
  N = bit 63, Z = (64-bit result == 0), C = carry out of bit 63, V = signed overflow.
Conditions: hs = C, lo = !C, hi = C && !Z, ls = !(C && !Z), eq = Z, ne = !Z.
After a borrow chain C is the full-width "no borrow" result, while Z reflects only the last word, so only
a carry-only condition is exact there; the old b.hi meter compare is run too and shown to fail.
The heapsort level's two 10-word record compares (proof.md 0-ii) run on 80-byte records in memory.
The two-word doubling countdown and interval doubling of proof-impl.md 0-i run across the 2^64 boundary,
and the one-word form (lsl / subs, b.eq) is run to show it wrapping at 2^63 -> 0.
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
    mem maps byte addresses to 64-bit words (ldr only)."""
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
        if op == "ldr":
            d, _, m = rest.partition(",")
            base, _, off = m.strip().strip("[]").partition(",")
            put(d.strip(), mem[val(base) + (val(off) if off else 0)])
            continue
        a = [x.strip() for x in rest.split(",")]
        if op.startswith("b.") or op == "cbz":
            assert taken is None, "one conditional branch, at the end"
            taken = val(a[0]) == 0 if op == "cbz" else bool(COND[op[2:]](flags))
        elif op == "mov":
            put(a[0], val(a[1]))
        elif op == "and":
            put(a[0], val(a[1]) & val(a[2]))
        elif op == "add":
            put(a[0], val(a[1]) + val(a[2]))
        elif op == "adc":
            put(a[0], val(a[1]) + val(a[2]) + int(flags["C"]))
        elif op == "sbc":
            put(a[0], val(a[1]) - val(a[2]) - 1 + int(flags["C"]))
        elif op == "orr":
            put(a[0], val(a[1]) | val(a[2]))
        elif op == "lsl":
            put(a[0], val(a[1]) << int(a[2][1:], 0))
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
STUCK = "adds x8, x8, #1; adc x9, x9, xzr; mov x16, #0x1E000000; cmp x9, x16; b.hs restart"
REPORT_CMP = "cmp x10, x11; b.hs"
CHUNK = "subs x12, x12, #1; b.eq switch"
STUCK_HI = 120 << 22          # 120 / theta = 120 * 2^86 = STUCK_HI * 2^64; movz #0x1E00, lsl 16
STUCK_CAP = STUCK_HI << 64

R_MOD = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
SAMPLE_CMP = "cmp x0, x4; sbcs xzr, x1, x5; sbcs xzr, x2, x6; sbcs xzr, x3, x7; b.lo accept"

CHILD_CMP = "; ".join(["ldr x%d, [x0, #%d]" % (2 + k, 8 * k) for k in range(10)]
                      + ["ldr x%d, [x1, #%d]" % (19 + k, 8 * k) for k in range(10)] + ["cmp x11, x28"]
                      + ["sbcs xzr, x%d, x%d" % (2 + k, 19 + k) for k in range(8, -1, -1)] + ["csel x12, x1, x0, lo"]
                      + ["csel x%d, x%d, x%d, lo" % (2 + k, 19 + k, 2 + k) for k in range(10)])
PARENT_CMP = "; ".join(["ldr x%d, [x13, #%d]" % (19 + k, 8 * k) for k in range(10)] + ["cmp x28, x11"]
                       + ["sbcs xzr, x%d, x%d" % (19 + k, 2 + k) for k in range(8, -1, -1)] + ["b.hs sift_done"])

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

    cases = edge3()
    for _ in range(n_random):
        x = rnd192()
        cases.append(x if isinstance(x, tuple) else (x, rng.getrandbits(192)))
    ok = sum(run_seq(METER_CMP, regs3(m, l))[1] == (m >= l) for m, l in cases)
    old_bad = [(m, l) for m, l in cases if run_seq(OLD_METER_CMP, regs3(m, l))[1] != (m > l)]
    verdict_case = ((3 << 128) + 9, (3 << 128) + 4)
    vc_old = run_seq(OLD_METER_CMP, regs3(*verdict_case))[1]
    vc_new = run_seq(METER_CMP, regs3(*verdict_case))[1]
    assert ok == len(cases) and vc_new is True and vc_old is False and old_bad
    res["meter_compare"] = {"listing": METER_CMP, "predicate": "branch iff meter >= L (192-bit unsigned)",
                            "cases": len(cases), "pass": True,
                            "old_b.hi_listing": OLD_METER_CMP, "old_b.hi_claimed_predicate": "meter > L",
                            "old_b.hi_wrong_cases": len(old_bad),
                            "verdict_case_equal_top_words_larger_low_128_bits":
                                {"old_b.hi_branches": vc_old, "b.hs_branches": vc_new}}

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

    sc = [0, 1, M64, M64 + 1, STUCK_CAP - 2, STUCK_CAP - 1, STUCK_CAP, STUCK_CAP + M64, ((STUCK_HI - 1) << 64) + M64,
          ((STUCK_HI + 1) << 64), (1 << 127), (1 << 128) - 2]
    sc += [rng.getrandbits(rng.choice((13, 64, 77, 127))) for _ in range(n_random)]
    for v in sc:
        r, t = run_seq(STUCK, {"x8": v & M64, "x9": v >> 64, "x16": 0})
        assert join([r["x8"], r["x9"]]) == v + 1 and t == (v + 1 >= STUCK_CAP)
    res["stuck_walk_counter"] = {"listing": STUCK, "predicate": "counter' = counter + 1 (128-bit); branch iff "
                                 "counter' >= 0x1E000000 * 2^64 = 120/theta", "cases": len(sc), "pass": True}

    rc = [(0, 1 << 55), ((1 << 55) - 1, 1 << 55), (1 << 55, 1 << 55), ((1 << 55) + 1, 1 << 55), (M64, 0), (0, M64)]
    rc += [(rng.getrandbits(56), rng.getrandbits(56)) for _ in range(n_random)]
    for a, b in rc:
        assert run_seq(REPORT_CMP + " cap", {"x10": a, "x11": b})[1] == (a >= b)
    res["one_word_counter_compare"] = {"listing": REPORT_CMP, "predicate": "branch iff count >= cap (64-bit unsigned)",
                                       "cases": len(rc), "pass": True}

    cc = [1, 2, 0, M64, 1 << 63] + [rng.getrandbits(64) for _ in range(n_random)]
    for v in cc:
        r, t = run_seq(CHUNK, {"x12": v})
        assert r["x12"] == (v - 1) & M64 and t == (v == 1)
    res["chunk_counter"] = {"listing": CHUNK, "predicate": "x12' = x12 - 1; branch iff x12' == 0 (one word)",
                            "cases": len(cc), "pass": True}

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

    W = 1 << 640
    rec = [(0, 0), (W - 1, W - 1), (5, 4), (4, 5), (W - 1, W - 2), (W - 2, W - 1)]
    for hi in (0, 1 << 62, M64):
        a = (hi << 576) + (5 << 512)
        rec += [(a, a - 1), (a - 1, a)]
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
        for tmp in (left, right):
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
    doubling_checks(res, rng)
    res["all_pass"] = True
    return res

COUNTDOWN2 = "subs x20, x20, #1; sbc x21, x21, xzr; orr x22, x20, x21; cbz x22, store"
DOUBLE2 = "mov x20, x23; mov x21, x24; adds x23, x23, x23; adc x24, x24, x24"
COUNTDOWN1, DOUBLE1 = "subs x12, x12, #1; b.eq store", "lsl x13, x13, #1; mov x12, x13"

def doubling_checks(res, rng):
    """proof-impl.md 0-i: countdown (x20 low, x21 high) decrements as one 128-bit value, branches iff it reaches 0; the
    interval (x23, x24) is copied into the countdown and then doubled as a 128-bit value 2^63 .. 2^94; from the 0-iv start
    values (both 1) the stores fall at segment positions 2^j - 1. One-word form wraps."""
    vals = [1, 2, (1 << 64) - 1, 1 << 64, (1 << 64) + 1, 1 << 65, (1 << 92) - 1, 1 << 92, (1 << 93) - 1]
    vals += [rng.getrandbits(93) | 1 << rng.randrange(93) for _ in range(3000)]
    for v in vals:
        r, t = run_seq(COUNTDOWN2, {"x20": v & M64, "x21": v >> 64})
        assert r["x20"] + (r["x21"] << 64) == v - 1 and t == (v == 1), v
    x23, x24 = 1 << 63, 0
    for j in range(64, 95):
        r, t = run_seq(DOUBLE2, {"x23": x23, "x24": x24, "x20": 0, "x21": 0})
        assert t is None and r["x20"] + (r["x21"] << 64) == x23 + (x24 << 64) == 1 << (j - 1), j
        x23, x24 = r["x23"], r["x24"]
        assert x23 + (x24 << 64) == 1 << j, j
    st = {"x20": 1, "x21": 0, "x23": 1, "x24": 0}
    stores = []
    for pos in range(128):
        r, t = run_seq(COUNTDOWN2, st)
        st.update({"x20": r["x20"], "x21": r["x21"]})
        if t:
            stores.append(pos)
            r, _ = run_seq(DOUBLE2, st)
            st.update({k: r[k] for k in ("x20", "x21", "x23", "x24")})
    assert stores == [(1 << j) - 1 for j in range(8)], stores
    st = {"x20": 3, "x21": 0, "x23": 0, "x24": 1}
    fired = []
    for pos in range(3):
        r, t = run_seq(COUNTDOWN2, st)
        st.update({"x20": r["x20"], "x21": r["x21"]})
        fired.append(t)
    r, _ = run_seq(DOUBLE2, st)
    assert fired == [False, False, True] and r["x20"] + (r["x21"] << 64) == 1 << 64 and r["x23"] + (r["x24"] << 64) == 1 << 65
    r1, _ = run_seq(DOUBLE1, {"x13": 1 << 63, "x12": 0})
    _, t1 = run_seq(COUNTDOWN1, {"x12": r1["x12"]})
    assert r1["x13"] == 0 and r1["x12"] == 0 and t1 is False
    _, t2 = run_seq(COUNTDOWN1, {"x12": 1})
    assert t2 is True
    assert COUNTDOWN2.count(";") + 1 == 4 and DOUBLE2.count(";") + 1 == 4
    res["doubling_countdown_two_word"] = {"listing": COUNTDOWN2, "instructions": 4, "cases": len(vals), "pass": True,
                                          "predicate": "128-bit decrement; branch iff the countdown reaches 0"}
    res["doubling_interval_two_word"] = {"listing": DOUBLE2, "instructions": 4, "cases": 31, "pass": True,
                                         "predicate": "128-bit interval copied into the countdown, then doubled, 2^63 .. 2^94"}
    res["doubling_schedule"] = {"start": "countdown 1, interval 1 (0-iv)", "stores_at": stores, "pass": True,
                                "predicate": "stores at segment positions 2^j - 1; a countdown of 3 with interval 2^64 fires on the third step and reloads 2^64"}
    res["doubling_one_word_wraps"] = {"listing": DOUBLE1 + "; " + COUNTDOWN1, "lsl_from_2^63": r1["x13"],
                                      "countdown_fires": t1, "pass": True, "predicate": "one-word form wraps (not used)"}

if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
