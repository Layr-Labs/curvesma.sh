"""Checks of the shared table banks, run by experiment.py in its own namespace (after its pipeline, proof_txt and certificate exist):
producer macros against the frame builders, the carry counters of P1_3WIN, the counts for 1 to 3 Y points and the traffic table of proof.md Section 3.4."""
# 4a. The Y producer macros (secp_fe.s WORDS: OFFS3_86 for E and G, OFF3FF for F and fwd) store exactly the pointer and
#     tag words that the kernels consume above: executed on 200 random variables and compared with yptr_words.
YFRAME_DOC = [96, 104, 112, 120, 272, 280, 288, 296, 304, 312, 368, 408, 416, 424, 432, 440, 464, 472, 480, 488, 496, 504]
XFRAME_DOC = {0: [320, 328, 336, 344, 352, 360, 376], 1: [40, 48, 320, 328, 336, 344, 352, 360], 2: [320, 328, 336, 344, 352, 360]}
for j_ in range(6): assert sorted(a_ for a_ in xframe_words(j_, 0)) == XFRAME_DOC[0 if j_ == 0 else 1 if j_ & 1 else 2], j_
_SECP = (K / "secp_fe.s").read_text()
_DRV = '''
.globl t_e
t_e:
OFFS3_86 x10,x11,x12,x13,280,288,296,272,304,312,x14,x15,x16,0,16,32
ret
.globl t_f
t_f:
OFF3FF
ret
.globl t_g
t_g:
OFFS3_86 x14,x15,x16,x17,472,480,488,464,496,504,x10,x18,x19,48,64,80
ret
'''
code_drv, lbl_drv = a64.a64_load(_SECP + _DRV)
class _ZM(dict):
    def __missing__(self, a): return 0
prod_rng = random.Random(SEED + 7); prod_counts = {}
for _ in range(200):
    bb_ = prod_rng.randrange(p); b2_ = 3 * bb_ * bb_ % p; b3_ = pow(bb_, 3, p)
    e_ = Y_VIEW + 512 * prod_rng.randrange(1 << 20)
    want = yptr_words(bb_, b2_, b3_ & ((1 << 223) - 1), e_)
    zm = _ZM()
    for ent, v, lim in (("t_e", bb_, "x10"), ("t_f", b3_, "x10"), ("t_g", b2_, "x14")):
        regs = {"x8": e_, "x30": SNT}
        for i in range(4): regs["x%d" % (int(lim[1:]) + i)] = (v >> (64 * i)) & M64
        prod_counts[ent] = a64.a64_run(code_drv, lbl_drv, ent, regs, zm, out_regs=dict(regs))
    got = {a_: zm[a_] for a_ in want}
    assert got == want and set(zm) == set(want), [hex(a_ - e_) for a_ in set(want) | set(zm) if dict.get(zm, a_) != want.get(a_)]   # exactly these words, no other store
    assert sorted(a_ - e_ for a_ in want) == YFRAME_DOC, "Y-frame words differ from the table of proof.md Section 3.4"
assert prod_counts["t_e"] == prod_counts["t_g"]
res_producer = dict(vectors=200, insns_e=prod_counts["t_e"] - 1, insns_f=prod_counts["t_f"] - 1, insns_g=prod_counts["t_g"] - 1)
assert res_producer['insns_e'] + res_producer['insns_f'] + res_producer['insns_g'] == cert.YR_WORDS_NEW   # the ret of each driver is not counted

# 4c. The two carry folds of P1_3WIN (Lemma V's c3 and bc): the shipped macro is run up to its two madds; the carry counters x0 (the s chain, with the
#     addend 1/D) and x2 (the t chain) reach exactly 3 and 2 when every row is p - 1 and 1/D = 2^256 - 1, and stay within 3 and 2 on random inputs.
_LT = (K / "leaf.s").read_text()
code_p1, lbl_p1 = a64.a64_load(_LT[:_LT.index(".globl fwd")] + "\n.globl t_p1\nt_p1:\nP1_3WIN\nret\n")
_mad = [i for i, ins in enumerate(code_p1) if ins in (("madd", ["x30", "x2", "x4", "x26"]), ("madd", ["x2", "x0", "x4", "x8"]))]
assert len(_mad) == 2
_stop = next(i for i in range(_mad[-1], len(code_p1)) if code_p1[i][0] == "ld64b")
class _CM(GMem):
    seen = None
    def instruction(self, pc, J, o, R, ev):
        if pc == _mad[0]: self.seen["x2"] = max(self.seen.get("x2", 0), R["x2"])
        if pc == _mad[1]: self.seen["x0"] = max(self.seen.get("x0", 0), R["x0"])
def p1_counters(z):
    mem = _CM(); mem.seen = {}; fr = Y_VIEW + (1 << 20)
    for a_, v_ in xframe_words(0, fr).items(): mem[a_] = v_
    regs = {"x%d" % (8 + i): (z >> (64 * i)) & M64 for i in range(4)}
    regs.update(x4=C0, x29=fr + 320, sp=Y_VIEW + (1 << 22), x30=SNT)
    G_LD64_P1 = [(lo, hi) for _, lo, hi in cert.G_MEM_REGIONS] + [(cert.T_AP_LO, cert.T_AP_HI)]
    a64.a64_run(code_p1, lbl_p1, "t_p1", regs, mem, ls64=(G_LD64_P1, []), stop_pc=_stop)
    return mem.seen["x0"], mem.seen["x2"]
carry_rng = random.Random(SEED + 11)
WORST[0] = True; carry_worst = p1_counters((1 << 256) - 1); WORST[0] = False
carry_random = [max(c_) for c_ in zip(*[p1_counters(carry_rng.getrandbits(256)) for _ in range(100)])]
assert carry_worst == (3, 2) and carry_random[0] <= 3 and carry_random[1] <= 2, (carry_worst, carry_random)
res_carry = dict(worst_rows=list(carry_worst), random_max=carry_random, bounds=[3, 2])

# 4d. Both outcomes of the two folds of P1_3WIN on exact integers (proof.md Section 5). The model returns the P1 rows that CRAFT names; the s chain sums the input z = 1/D and the
#     three variant-A rows, the t chain the three variant-B rows. With S the exact sum, x0 = S >> 256 the carry count and R = (S mod 2^256) + x0 C0 = S - x0 p, the final madd drops
#     the carry out of the low limb: the kernel must return R or R - 2^64, a drop needs s = S mod p with s mod 2^64 in [0, (x0 + 1) C0), and R = s + p needs s < (x0 + 1) C0.
_GLD = [(lo, hi) for _, lo, hi in cert.G_MEM_REGIONS] + [(cert.T_AP_LO, cert.T_AP_HI)]
def p1_fold(z, A, B):
    CRAFT[0] = {k_: (A[k_], B[k_]) for k_ in range(3)}
    try:
        mem = GMem(); fr = Y_VIEW + (1 << 20); out = {}
        for a_, v_ in xframe_words(0, fr).items(): mem[a_] = v_
        regs = {"x%d" % (8 + i): (z >> (64 * i)) & M64 for i in range(4)}
        regs.update(x4=C0, x29=fr + 320, sp=Y_VIEW + (1 << 22), x30=SNT)
        a64.a64_run(code_p1, lbl_p1, "t_p1", regs, mem, ls64=(_GLD, []), stop_pc=_stop, out_regs=out)
    finally: CRAFT[0] = None
    L_ = lambda *ns: sum(out[n_] << (64 * i) for i, n_ in enumerate(ns))
    return L_("x2", "x9", "x10", "x11"), L_("x30", "x27", "x28", "x29"), out["x0"]          # s chain, t chain, s chain's counter (the t chain's was folded into x30)
def fold_model(S):
    x0, Lw = S >> 256, S & ((1 << 256) - 1)
    R = Lw + x0 * C0; drop = ((Lw & M64) + x0 * C0) >> 64; s = S % p
    assert R in (s, s + p) and (R == s or s < (x0 + 1) * C0)           # R = S - x0 p lies in [x0 C0, p + (x0 + 1) C0)
    assert not drop or (s & M64) < (x0 + 1) * C0                       # a drop needs s mod 2^64 in the window
    assert drop or s < (x0 + 1) * C0 or R == s                         # outside the window, without a drop, the result is the canonical residue
    return R - (drop << 64), drop
fold_rng = random.Random(SEED + 13); HI = ((1 << 192) - 2) << 64
fold_cases = [(fold_rng.getrandbits(256), [fold_rng.randrange(p) for _ in range(3)], [fold_rng.randrange(p) for _ in range(3)]) for _ in range(12)]
for tS, tT in ((M64, M64), (M64 + 1 - 3 * C0, M64 + 1 - 2 * C0), (M64 - 3 * C0, M64 - 2 * C0), (M64, M64 - 2 * C0), (M64 - 3 * C0, M64), (0, 0)):   # last low limb of the sum: inside, at the edge of, just outside the window
    b_ = fold_rng.getrandbits(64); z_ = (fold_rng.getrandbits(255) | 1 << 255) & ~M64 | ((tS - 2 * (p - 1) - HI - b_) & M64)
    A_ = [p - 1, p - 1, HI + b_]; bT = (tT - 2 * (p - 1) - HI) & M64
    fold_cases.append((z_, A_, [p - 1, p - 1, HI + bT]))
fold_stat = dict(cases=0, drops=0, no_drops=0, max_x0=0, max_x2=0)
for z_, A_, B_ in fold_cases:
    rs_, rt_, x0_ = p1_fold(z_, A_, B_)
    ms_, ds_ = fold_model(z_ + sum(A_)); mt_, dt_ = fold_model(sum(B_))
    assert (rs_, rt_) == (ms_, mt_), "the P1_3WIN fold differs from the exact integers"
    assert x0_ == (z_ + sum(A_)) >> 256, "the carry counter differs from floor(S / 2^256)"
    fold_stat["cases"] += 2; fold_stat["drops"] += ds_ + dt_; fold_stat["no_drops"] += 2 - ds_ - dt_
    fold_stat["max_x0"] = max(fold_stat["max_x0"], x0_); fold_stat["max_x2"] = max(fold_stat["max_x2"], sum(B_) >> 256)
assert fold_stat["drops"] >= 4 and fold_stat["no_drops"] >= 4 and fold_stat["max_x0"] == 3 and fold_stat["max_x2"] == 2, fold_stat
res_fold = fold_stat

# 4e. The writer of the new X-entry words (Row 7: 12 an entry on top of 8cc64cb7's 18): eighteen consecutive 512-byte entries, j = 0 .. 5 three times.
_ENT = ".globl t_entry\nt_entry:\nmov x4,#0x4000000000000000\nmov x5,#0x3d1\nmovk x5,#1,lsl #32\n1:\n"
for j_ in range(6):
    for k_ in range(3): _ENT += "mov x%d,#0x%x\n" % (1 + k_, (6 + 3 * j_ + k_) << 36)
    _ENT += "stp x1,x2,[x0,#320]\nstp x3,x4,[x0,#336]\nstp x4,x4,[x0,#352]\n"
    if j_ == 0: _ENT += "str x5,[x0,#376]\n"
    if j_ & 1: _ENT += "mov x7,#0x%x\nmov x8,#0x%x\nstp x7,x8,[x0,#40]\n" % (25 << 36, 26 << 36)
    _ENT += "add x0,x0,#512\n"
_ENT += "subs x6,x6,#1\nb.ne 1b\nret\n"
code_ent, lbl_ent = a64.a64_load(_ENT)
mem_ent = _ZM(); ent0 = Y_VIEW + (1 << 24); ent_iters = 3
cnt_ent = a64.a64_run(code_ent, lbl_ent, "t_entry", dict(x0=ent0, x6=ent_iters, x30=SNT), mem_ent)
want_ent = {}
for l_ in range(6 * ent_iters): want_ent.update(xframe_words(l_ % 6, ent0 + 512 * l_))
assert dict(mem_ent) == want_ent, "t_entry stores other words than the X frame words"
ent_loop = (cnt_ent - 4) // ent_iters                        # three set-up instructions and the ret are not in the loop
assert (cnt_ent - 4) % ent_iters == 0 and ent_loop == 54 and ent_loop <= 12 * 6, cnt_ent
assert "54 instructions for six entries" in proof_txt
res_entry = dict(entries=6 * ent_iters, loop_insns_per_six_entries=ent_loop, charged_per_entry=cert.ENTRY_X - 18)

# 4f. Random words in the unused half of every BT row: the 144-leaf labels still match (BT rows: the builder writes 32 of their 64 bytes).
POISON[0] = True
try: r_poison = run_leaf_pipeline(nys=1)
finally: POISON[0] = False
assert r_poison["matched"] == 36, r_poison["matched"]
res_poison = r_poison["matched"]


for n_ in (1, 2, 3):                              # certificate.py's FWD_REM and BWD_REM are these measured counts
    r_ = run_leaf_pipeline(nys=n_)
    assert r_["c_fwd"] == cert.fwd_count(6 * n_) and r_["c_bwd"] == cert.bwd_count(6 * n_) and r_["matched"] == 36 * n_, n_
_ld = res_ok["mem"].lds                           # every table sum reads three windows: BLK3F 3 x 2 a incidence, P1_3WIN 3, JK 3 x 2 lines x 2 passes
assert (_ld["BT"], _ld["P1"], _ld["E"], _ld["G"], _ld["F"], _ld["FWD"]) == (24 * 2 * 3, 24 * 3, 4 * 2 * 3 * 2, 4 * 2 * 3 * 2, 4 * 2 * 2 * 2, 4 * 2 * 2), _ld
assert res_ok["mem"].nloads == 24 * cert.T_LOADS_INC, (res_ok["mem"].nloads, cert.T_LOADS_INC)

_sto, _lds = res_ok["mem"].sto, res_ok["mem"].lds
for needle in (f"| `BT` | {_sto['BT']} | {_lds['BT']} |", f"| `P1` | {_sto['P1']} | {_lds['P1']} |", f"| `E`, `G`, `F` (in `JK`) | {_sto['E'] + _sto['G'] + _sto['F']} | {_lds['E'] + _lds['G'] + _lds['F']} |", f"| `fwd` | {_sto['FWD']} | {_lds['FWD']} |", f"| all | {sum(_sto.values())} | {sum(_lds.values())} |"):
    assert needle in proof_txt, f"Missing {needle} in proof.md"
