"""Section 9 of experiment.py, kept in its own file because experiment.py travels as an archive member capped at
65,536 bytes. experiment.py executes this file in its own namespace after section 7, so it uses that section's
kernels, models and fixtures (code_rt, Relay, RMem, PMem, g_run, h_fixture, the ycmd_gen fixture).
"""
# 9. Checks of the reply handle end to end, the Phase-0 pre-ticket filter, the top index of each bit-field
#    insert, CN's boundary, and a static check that each map region an encoding carries fits that encoding's range.
RT_SRC = (K / "route.s").read_text()
def rt_part(a, b):
    i = RT_SRC.index(a); return RT_SRC[i:RT_SRC.index(b, i)]
# A frame-logging generic_send that uses only x16-x18 (which GENERIC_COMMIT already clobbers): each send appends its
# 80-byte frame (route lo/hi, command, pad, words 0..5) at the cursor kept at SEND_LOG.
SEND_LOG = 0x10000
LOG_STUB = ("\n.text\ngeneric_send:\nmov x18,#0x10000\nldr x18,[x18]\n"
            + "".join("ldp x16,x17,[x0,#%d]\nstp x16,x17,[x18,#%d]\n" % (16 * i, 16 * i) for i in range(5))
            + "add x18,x18,#80\nmov x16,#0x10000\nstr x18,[x16]\nret\n")
HP_MACROS = (rt_part(".macro GENERIC_COMMIT", ".text") + rt_part(".macro REQPUT", ".macro LINE64END")
             + rt_part(".macro LINE64END", ".macro OR18CBZ"))
def hp_prog(role, body): return a64.a64_load(materialize(HP_MACROS + "\n.text\n" + body + LOG_STUB, role=role))
code_hH, lbl_hH = hp_prog(2, rt_part(".globl reqh", ".globl fwdq\n"))
code_hR, lbl_hR = hp_prog(1, rt_part(".globl dlvq", ".globl reqh") + rt_part(".globl fwdq\n", ".globl sendr")
                          + rt_part(".globl fwdr\n", "// Retained historical"))
code_hG, lbl_hG = hp_prog(0, rt_part(".globl dlvq", ".globl reqh") + rt_part(".globl sendr", ".globl fwdr\n"))
class M0(dict):
    def __missing__(self, a): return 0
def run_log(code, lbl, entry, regs, mem):
    mem[SEND_LOG] = SEND_LOG + 64
    c = a64.a64_run(code, lbl, entry, dict(regs, x30=SNT, sp=0x50000), mem)
    n = (mem[SEND_LOG] - SEND_LOG - 64) // 80
    return c, [[mem[SEND_LOG + 64 + 80 * f + 8 * w] for w in range(10)] for f in range(n)]
def region(role, name): return next((lo, hi) for nm, lo, hi in cert.ROLE_MAPS[role] if nm == name)
RCB_LO, RCB_END = region("R", "DPW_RCB"); RCB0 = RCB_LO + 65536
assert RCB0 % 64 == 0 and (RCB_END - RCB0) % 64 == 0 and (RCB_END - RCB0) // 64 == cert.RELAY_CAP and RCB_END <= 1 << 64
def rcb_entry_ok(e): return RCB0 <= e < RCB_END and e % 64 == 0
# The relay's command-1 cursor: generic_handlers_init stores x3 + 8 at x3, so with x3 = RCB0 - 8 the cursor word sits
# just below the first entry and holds RCB0; its descriptor carries the 64-B stride and the exclusive limit.
ghm = M0()
a64.a64_run(code_rt, lbl_rt, "generic_handlers_init", {"x0": 0x90000, "x1": 0x1000, "x2": 0x2000, "x3": RCB0 - 8, "x4": RCB_END,
            "x5": 0x91000, "x30": SNT}, ghm)
assert ghm[RCB0 - 8] == RCB0 and all(ghm[0x90000 + 64 * i + 8] == 64 and ghm[0x90000 + 64 * i + 16] == RCB0 - 8
                                     and ghm[0x90000 + 64 * i + 24] == RCB_END for i in range(4))
HCC_LO, HCC_HI = region("H", "H_HCC_CAND")
HCC_REPLY0 = HCC_LO + 0x105000 + 16 * cert.CAND_CAP          # header, then the candidate list, then the reply slots
assert 0x105000 == (1 << 20) + 8 * 2048 + 4096 and HCC_REPLY0 + 3072 * cert.CAND_CAP == HCC_HI

# 9a. Handle round trips. fwdq1 (ordinary) and fwdq (Phase 0) put the 64-B-aligned request-entry address itself into
#     the request; fwdr strips the 6-bit message index, reads the entry's home word and reply slot and adds 32 a
#     message. At DPW_RCB's first and last entries and at 2^56 - 64, 2^56 and 2^56 + 64.
HOME_W, SLOT_T = (1 << 45) + 8 * 12345, (1 << 62) + 0x100000
def handle_of(entry, fn):
    m = M0(); m[entry + 8], m[entry + 16] = (5 << 3) | (7 << 13), 1 << 9; m[entry + 24], m[entry + 32] = HOME_W, SLOT_T
    _, fr = run_log(code_hR, lbl_hR, fn, {"x0": 0x30000, "x1": entry, "x6": entry + 64, "x7": 42}, m)
    assert len(fr) == 1 and fr[0][2] == 1
    return fr[0][7], m
def fwdr_frame(m, word0, k):
    rep = 0x40000; m[rep] = word0
    for i in range(1, 5): m[rep + 8 * i] = 0x1000 + 16 * k + i
    m[BASES[1] + 24] = 777
    _, fr = run_log(code_hR, lbl_hR, "fwdr", {"x0": 0x30000, "x1": rep, "x6": rep + 64}, m)
    assert len(fr) == 1
    return fr[0]
HANDLE_AT = [RCB0, RCB_END - 64, (1 << 56) - 64, 1 << 56, (1 << 56) + 64]
for e in HANDLE_AT:
    assert rcb_entry_ok(e), hex(e)
    for fn in ("fwdq1", "fwdq"):
        h, m = handle_of(e, fn)
        assert h == e, (fn, hex(e), hex(h))
        for k in (0, 47):
            f = fwdr_frame(m, h + k, k)
            assert f[2] == 0 and f[4] == HOME_W and f[5] == SLOT_T + 32 * k and f[6:10] == [0x1000 + 16 * k + i for i in (1, 2, 3, 4)]
# Negative control: dlvq with the cursor at DPW_RCB's end clamps (no entry written, cursor left at the limit).
cm_ = M0(); cm_[RCB0 - 8] = RCB_END; fb_ = 0x32000
for w in range(6): cm_[fb_ + 8 * w] = 0x7700 + w
a64.a64_run(code_hR, lbl_hR, "dlvq", {"x0": fb_, "x9": 64, "x10": RCB0 - 8, "x11": RCB_END, "x29": 0x71000, "x28": 1, "x30": SNT}, cm_)
nc_rcb_end_clamped = cm_[RCB0 - 8] == RCB_END and all(cm_[RCB_END + 8 * w] == 0 for w in range(5))
doc_retired_shift8_2p56 = (((1 << 56) << 8) & M64) >> 8 != 1 << 56                # documentation: the old entry << 8 handle

# 9b. One request's 48 reply messages end to end: reqh sends the candidate's two requests (reply slots S and
#     S + 1536); the relay's dlvq queues the first in DPW_RCB; fwdq1 forwards it with its handle; the server's sendr
#     sends 48 messages (six 256-B replies); dlvq queues each; fwdr forwards them, and message k lands at S + 32k.
#     Negative control: a 49th message lands at S + 1536, the next request's slot.
hm = M0(); LST, RU, RV, S0 = 0x60000, 0x61000, 0x62000, HCC_REPLY0
hm[LST], hm[LST + 8] = RU, RV
hm[RU + 16], hm[RU + 24], hm[RV + 16], hm[RV + 24] = (5 << 3) | (7 << 13), 1 << 9, (6 << 3) | (8 << 13), 2 << 9
hm[BASES[2] + 24], hm[BASES[2] + 32] = 3 << 42, 0x12
_, rq = run_log(code_hH, lbl_hH, "reqh", {"x0": 0x30000, "x12": LST, "x9": 1, "x13": S0, "x10": HOME_W}, hm)
assert len(rq) == 2 and [f[2] for f in rq] == [1, 1] and [f[9] for f in rq] == [S0, S0 + 1536] and all(f[8] == HOME_W for f in rq)
assert [f[6:8] for f in rq] == [[hm[RU + 16], hm[RU + 24]], [hm[RV + 16], hm[RV + 24]]]
E0 = RCB0; CUR, CNT = RCB0 - 8, 0x71000
def rcb_chain(n_msgs):
    rm = M0(); rm[CUR] = ghm[RCB0 - 8]
    def put_q(frame):
        fb = 0x32000
        for w in range(6): rm[fb + 8 * w] = frame[4 + w]
        a64.a64_run(code_hR, lbl_hR, "dlvq", {"x0": fb, "x9": 64, "x10": CUR, "x11": RCB_END, "x29": CNT, "x28": 1, "x30": SNT}, rm)
    put_q(rq[0])
    assert rm[CUR] == E0 + 64 and [rm[E0 + 8 * w] for w in range(5)] == rq[0][5:10]
    _, fq = run_log(code_hR, lbl_hR, "fwdq1", {"x0": 0x30000, "x1": E0, "x6": E0 + 64, "x7": 42 << 3}, rm)
    assert len(fq) == 1 and fq[0][7] == E0
    gq = M0(); PAY = 0x80000; GA = 0xd800_0000_0000_0000; gq[GA - 8] = GA        # the G's request arena (callbacks.md)
    for w in range(6): gq[0x32000 + 8 * w] = fq[0][4 + w]
    a64.a64_run(code_hG, lbl_hG, "dlvq", {"x0": 0x32000, "x9": 64, "x10": GA - 8, "x11": GA + (64 << 20), "x29": 0x71000, "x28": 1, "x30": SNT}, gq)
    assert gq[GA + 16] == E0 and gq[GA + 24] == 42 << 3                          # handle and relay word, as the server reads them
    for k in range(n_msgs):
        for i in range(4): gq[PAY + 32 * k + 8 * i] = (k << 8) | i
    c_s, rp = run_log(code_hG, lbl_hG, "sendr", {"x0": 0x30000, "x1": PAY, "x7": gq[GA + 24], "x8": gq[GA + 16], "x9": n_msgs}, gq)
    assert len(rp) == n_msgs and [f[5] for f in rp] == [E0 + k for k in range(n_msgs)]
    for f in rp: put_q(f)
    rm[BASES[1] + 24] = 777
    c_r, rr = run_log(code_hR, lbl_hR, "fwdr", {"x0": 0x30000, "x1": E0 + 64, "x6": E0 + 64 + 64 * n_msgs}, rm)
    assert len(rr) == n_msgs
    for k, f in enumerate(rr):
        assert f[2] == 0 and f[4] == HOME_W and f[6:10] == [(k << 8) | i for i in range(4)], k
    return [f[5] for f in rr], c_s, c_r
dst48, c_sendr48, c_fwdr48 = rcb_chain(48)
assert dst48 == [S0 + 32 * k for k in range(48)] and all(S0 <= d < S0 + 1536 for d in dst48)
dst49, _, _ = rcb_chain(49)
nc_reply_49_outside_slot = dst49[48] == S0 + 1536 and not (S0 <= dst49[48] < S0 + 1536)

# 9c. pfilt, the pre-ticket close-pair filter, on accept, reject and swap paths, against psel's per-report test (a
#     model, and psel itself on the oriented report); then ptick over 64 candidates. A(U) is message 0 of the
#     candidate's first REQPUT slot, A(V) of the next (+1536); identities come from the candidate's record pointers.
def pf_ref(Aa, Ab, ia, ib):
    lo, hi = min(Aa, Ab), max(Aa, Ab); return int(lo != 0 and 0 < hi - lo < 2**63 and ia != ib)
def pf_setup(m, slot, lst, recs, Aa, Ab, ia, ib):
    put(m, slot, Aa); put(m, slot + 1536, Ab); m[lst], m[lst + 8] = recs
    m[recs[0] + 16], m[recs[0] + 24] = ia; m[recs[1] + 16], m[recs[1] + 24] = ib
pf_cases = []
for _ in range(300):
    Aa = rng.randrange(1, p); kind = rng.randrange(7)
    if kind == 0: Ab = Aa
    elif kind == 1: Ab = (Aa + rng.randrange(1, 2**63)) % p
    elif kind == 2: Ab = (Aa - rng.randrange(1, 2**63)) % p
    elif kind == 3: Ab = (Aa + 2**63 + rng.randrange(0, 2**20)) % p
    elif kind == 4: Ab = Aa ^ (1 << rng.randrange(64, 255))
    elif kind == 5: Aa, Ab = 0, rng.randrange(1, 2**63)
    else: Ab = rng.randrange(0, p)
    ia = (rng.getrandbits(64), rng.getrandbits(64)); ib = ia if rng.random() < 0.1 else (rng.getrandbits(64), rng.getrandbits(64))
    pf_cases.append((Aa, Ab, ia, ib))
pf_cases += [(5, 5 + 2**63 - 1, (1, 2), (3, 4)), (5 + 2**63 - 1, 5, (1, 2), (3, 4)), (5, 5 + 2**63, (1, 2), (3, 4)),
             (5 + 2**63, 5, (1, 2), (3, 4)), (2**64 - 1, 2**64, (1, 2), (3, 4)), (2**64, 2**64 - 1, (1, 2), (3, 4)),
             (7, 8, (9, 9), (9, 9)), (0, 0, (1, 2), (3, 4)), (2**192, 2**192 + 1, (1, 2), (1, 3))]
pf_steps, pf_paths = [], {}
for Aa, Ab, ia, ib in pf_cases:
    m = M0(); pf_setup(m, 0x100000, 0x200000, (0x300000, 0x301000), Aa, Ab, ia, ib); o = {}
    st = a64.a64_run(code_rt, lbl_rt, "pfilt", {"x0": 0x100000, "x18": 0x200000, "x30": SNT}, m, out_regs=o)
    want = pf_ref(Aa, Ab, ia, ib)
    assert o["x0"] == want, (hex(Aa), hex(Ab))
    pf_steps.append(st); pf_paths.setdefault(("pass" if want else "reject") + ("_swap" if Ab < Aa else ""), set()).add(st)
    pm = M0(); rb, ob = 0x200000, 0x300000; lo, hi = min(Aa, Ab), max(Aa, Ab)
    put(pm, rb, lo); put(pm, rb + 32, hi); pm[rb + 64], pm[rb + 72] = ia; pm[rb + 80], pm[rb + 88] = ib
    a64.a64_run(code_rt, lbl_rt, "psel", {"x0": rb, "x1": 1, "x2": ob, "x30": SNT, "sp": 0x400000}, pm)
    took = int(any(pm[ob + 8 * w] for w in range(4))) if lo else 0
    assert took == want, (hex(Aa), hex(Ab))
pf_max = max(pf_steps)
assert pf_max == 36 and max(pf_paths.get("pass", {0})) <= pf_max
mt = M0(); NPT = 64; npass = 0; pf_sum = 0
for c, (Aa, Ab, ia, ib) in enumerate(pf_cases[:NPT]):
    pf_setup(mt, 0x1000000 + 3072 * c, 0x2000000 + 16 * c, (0x3000000 + 64 * c, 0x3100000 + 64 * c), Aa, Ab, ia, ib)
    npass += pf_ref(Aa, Ab, ia, ib); pf_sum += pf_steps[c]
o = {}
c_ptick = a64.a64_run(code_rt, lbl_rt, "ptick", {"x1": 0x1000000, "x18": 0x2000000, "x2": NPT, "x30": SNT}, mt, out_regs=o)
ptick_overhead = c_ptick - pf_sum
assert o["x3"] == npass and ptick_overhead == 7 * NPT + 4 and cert.PF_CAND == pf_max + 7 + 4 + 2 * cert.SQRC

# 9d. The top index of each bit-field insert. H2B at d2 = 2^42 - 1 (LC2 2^46 - 8, RC2 2^47 - 8, HREADY offset
#     2^45 - 8, window 2^63 + 2^62 - 2^20); xemit and an isolated EM site at d1 = 2^40 - 1, each with a block end
#     (LC1 2^43 - 8, READY_R 2^44 - 8, slot 2^56 + 2^58 - 2^18); fin1 on that slot; HCHGI at HH index 2^50 - 1;
#     relay_slot_bind at member 1023.
D2T, D1T = (1 << 42) - 1, (1 << 40) - 1
mt2 = RMem(); rlt = Relay(mt2); mt2[RSLOT0 + 8] = 32; mt2[RSLOT0 + 16] = 1
for w, v in enumerate((D2T, 0x11, 0x22, 0x33)): mt2[RSLOT0 + 32 + 8 * w] = v
mt2[LC2B + 8 * D2T] = H2UB + (D2T << 20) + 524288
rlt.receive(capture=False)
assert LC2B + 8 * D2T == (1 << 46) - 8 and RC2B + 8 * D2T == (1 << 47) - 8 and 8 * D2T == (1 << 45) - 8
assert H2UB + (D2T << 20) == (1 << 63) + (1 << 62) - (1 << 20) and rlt.hpush == [(8 * D2T, 524288 | 3)]
assert mt2[LC2B + 8 * D2T] == H2UB + (D2T << 20) + 32 and mt2[RC2B + 8 * D2T] == 524288
assert [mt2[H2UB + (D2T << 20) + 524288 + 8 * w] for w in range(4)] == [D2T, 0x11, 0x22, 0x33]
SLOT_TOP = H1_BASE + (D1T << 18)
assert 8 * D1T == (1 << 43) - 8 and SLOT_TOP == (1 << 56) + (1 << 58) - (1 << 18) and READY + 8 * D1T == (1 << 44) - 8
gt = PMem({}); rt_ = {}; FR0 = 40 << 40; A_TOP = (D1T << 103) | (((1 << 49) - 1) << 207)
TAG_TOP = ((1 << 14) - 1) << 49                                                   # xemit's x6: pool << 49, pool < 2^14
for i in range(4): gt[FR0 + 376 + 8 * i] = (A_TOP >> (64 * i)) & M64
gt[8 * D1T] = SLOT_TOP + 32 + 32 * (K1 - 1)
_, pu_t = g_run(gt, rt_, "xemit", {"x1": FR0, "x3": 1, "x4": 32, "x6": TAG_TOP, "x8": 0x7fffffffff8, "x2": 0x5a5a, "x26": 0})
assert pu_t == [(8 * D1T, RTOKEN)] and gt[8 * D1T] == SLOT_TOP + 32 and gt.relays[D1T][RSLOT0 + 8] == 32 * K1
assert gt.relays[D1T][RSLOT0 + 131072 + 8] == TAG_TOP | ((1 << 49) - 1)           # record word 1 = u = pool << 49 | A[255:207]
LF_SRC = (K / "leaf.s").read_text()
LF_MACROS = "\n".join(mm.group(0) for mm in re.finditer(r"(?ms)^\.macro .*?^\.endm[ \t]*$", LF_SRC))
code_em, lbl_em = a64.a64_load(LF_MACROS + "\n.text\n.globl emtest\nemtest:\nand x0,x1,x7,lsr #21\nEM x5,x9,x10,x11,t0\nret\n")
ge = PMem({}); re_ = {}; ge[8 * D1T] = SLOT_TOP + 32 + 32 * (K1 - 1)
_, pu_e = g_run(ge, re_, "emtest", {"x1": 0x7fffffffff8, "x7": D1T << 24, "x4": 32, "x9": 0x99, "x10": 0xAA, "x11": 0xBB, "x2": 0x77},
                code=(code_em, lbl_em))
assert pu_e == [(8 * D1T, RTOKEN)] and ge[8 * D1T] == SLOT_TOP + 32 and ge.relays[D1T][RSLOT0 + 8] == 32 * K1
gf_ = PMem({}); rf_ = {}; gf_[8 * D1T] = SLOT_TOP + 32 + 32 * 3
for i in range(3):
    for w, v in enumerate(rec(i, 0x44)): gf_[SLOT_TOP + 32 + 32 * i + 8 * w] = v
_, pu_f = g_run(gf_, rf_, "fin1", {"x1": 1, "x7": 8 * D1T})
assert pu_f == [(8 * D1T, RTOKEN)] and gf_[8 * D1T] == SLOT_TOP + 32
assert sorted((d2, w0, w1) for d2, w0, w1, w2, w3 in rf_[D1T].delivered()) == sorted((r_[0], r_[0], r_[1]) for r_ in (rec(i, 0x44) for i in range(3)))
hmt, rht = h_fixture(); hmt[r2_base + 8] = ((1 << 50) - 1) | (1 << 55)
a64.a64_run(code_rt, lbl_rt, "dlvh", rht, hmt, ls64=H_LS64, limit=10**6)
assert hh_base + 8 * ((1 << 50) - 1) == (0xffa0 << 48) - 8 and hmt[(0xffa0 << 48) - 8] == r2_base
msb = RMem(); msb[RTBL0 + (1023 << 15) + 16384] = 0xABC; rsb = {"x4": 1023, "x30": SNT}
a64.a64_run(code_rt, lbl_rt, "relay_slot_bind", rsb, msb, out_regs=rsb)
assert rsb["sp"] == RTBL0 + (1023 << 15) and rsb["x0"] == 0xABC

# 9e. CN at its boundary: k = Q1GD writes CN's last entry (7b above, under the map's CN range); k = Q1GD + 1 passes
#     ycmd_gen's k - 1 < 2^45 guard, and its LS64 store falls outside CN, which the map's range rejects.
CN_LO, CN_HI = region("G", "CN")
assert (CN_LO, CN_HI) == (84 << 40, (84 << 40) + 64 * cert.Q1GD) and CN_LS64[1] == [(CN_LO, CN_HI)]
cmem[ycmd_ctrl + 24] = 1; cmem[ycmd_in + 8] = cert.Q1GD + 1
nc_cn_q1gd_plus_1 = False
try:
    a64.a64_run(code_ys, lbl_ys, "ycmd_gen", {"sp": (34 << 40) + 0x7000000, "x30": SNT}, cmem, ls64=CN_LS64)
except AssertionError as exc:
    nc_cn_q1gd_plus_1 = "LS64 outside a permitted region" in str(exc)

# 9f. Each map region an encoding carries fits that encoding's range (a cap change that breaks one fails here).
ENCODINGS = [
    ("G", "LC1", 0, 1 << 43, "EM/xemit/FLUSH7/fin1: LC1 = 8 d1 (mask 0x7fffffffff8)"),
    ("G", "H1_CONTIG", 1 << 56, (1 << 56) + (1 << 58), "fin1 and ~0x3ffff, EM/xemit tbnz #17, endpoint d1 = (a - 2^56) >> 18"),
    ("G", "CN", 84 << 40, (84 << 40) + (64 << 45), "ycmd_gen 84 TiB + 64 (k - 1), k - 1 < 2^45"),
    ("R", "LC2", 1 << 45, 1 << 46, "H2B bfi x10,d2,#3,#42 on 2^45"),
    ("R", "RC2", 3 << 45, 1 << 47, "hfshared LC2 | 2^46"),
    ("R", "R_SLOT_TABLES", 0xC000_0200_0000, 0xC000_0200_0000 + (1024 << 15), "relay_slot_bind bfi x6,x4,#15,#10"),
    ("R", "R_SLOTS_1024x256K", 196 << 40, (196 << 40) + (1024 << 18), "member slot 196 TiB + (m << 18)"),
    ("R", "H2U_CONTIG", 1 << 63, (1 << 63) + (1 << 62), "H2U endpoint d2 = (a - 2^63) >> 20, tbnz #19"),
    ("R", "DPW_RCB", 0, 1 << 64, "reply handle: 64-B-aligned entry address | 6-bit message index"),
    ("H", "HH_2p50", 0xff80 << 48, 0xffa0 << 48, "HCHGI bfi x10,u,#3,#50"),
    ("B", "B_DATA", 0, 1 << 62, "yproduce slot << 35, slot < 2^27"),
    ("B", "B_META", 1 << 62, (1 << 62) + (16 << 27), "2^62 + 16 slot"),
]
for role_, nm_, lo_, hi_, _ in ENCODINGS:
    a_, b_ = region(role_, nm_); assert lo_ <= a_ and b_ <= hi_, (nm_, hex(a_), hex(b_))
assert 6 * 8 * 32 == 1536 and 6 * 8 <= 64     # six 256-B replies = 48 messages, inside the 6-bit index
assert TAG_TOP < 1 << 63 and TAG_TOP & ((1 << 49) - 1) == 0 and (1 << 14) == cert.N_POOLS   # pool tag above A[255:207]

# 9g. dlvh's segmented path (final partial blocks and blocks below the R2 frontier) at K2B = 16384: aligned and odd
#     starts, lengths up to 16383, and full blocks with the frontier inside. Every record is inserted, and the
#     non-HCHGI count stays within DLVHS.
def dlvh_seg(start, length, frontier, nrec=16400):
    hm2 = ZMem()
    hm2[(1 << 55) - 160], hm2[(1 << 55) - 152] = HCC_LIST, HCC_LIST + (1 << 30)
    hm2[(1 << 55) - 144], hm2[(1 << 55) - 136] = nx_base, r2_base + (1 << 20)
    hm2[h_win - 24], hm2[h_win - 16], hm2[h_win - 8] = start, length, 0x7ffd_0000_0000_00
    hm2[0x7ffd_0000_0000_00] = frontier
    for i in range(nrec):
        hm2[r2_base + 32 * i], hm2[r2_base + 32 * i + 8] = 0x900000 + i, (i + 1) | (1 << 55)
        hm2[r2_base + 32 * i + 16], hm2[r2_base + 32 * i + 24] = 0xA00000 + i, 0xB00000 + i
    rh2 = {"x30": SNT}
    a64.a64_run(code_rt, lbl_rt, "dlvh0", rh2, hm2, ls64=H_LS64, out_regs=rh2)
    rh2["x29"], rh2["x30"] = 0x4000_0000_0000_0000, SNT
    c_ = a64.a64_run(code_rt, lbl_rt, "dlvh", rh2, hm2, ls64=H_LS64, limit=10**6)
    first, last = (max(start, frontier) - r2_base) // 32, (start + length - r2_base) // 32
    assert [i for i in range(nrec) if hm2[hh_base + 8 * (i + 1)] != 0] == list(range(first, last)), (start, length, frontier)
    return c_ - cert.HCHG_NEW * (last - first)
seg_shapes = ([(r2_base, 32 * n, r2_base) for n in (1, 2, 3, 255, 511, 16382, 16383)]
              + [(r2_base + 32, 32 * n, r2_base) for n in (1, 2, 3, 16382, 16383)]
              + [(r2_base, 32 * K2B, r2_base + 32 * k) for k in (1, 2, 511, 16383)] + [(r2_base, 32 * 16383, r2_base + 32)])
dlvhs_max = max(dlvh_seg(*s_) for s_ in seg_shapes)
assert dlvhs_max <= cert.DLVHS

# 9h. The relay's last expected final in a group (h2s0's finals counter at 1), then h2sd: against the same block as
#     an ordinary final, it adds exactly H2S_LAST (mov x6; sub x14) + H2S_DONE (h2sd), and h2sd marks the relay done.
def h2s_fin(n, finals):
    m = RMem(); rl = Relay(m, finals=finals); m[RSLOT0 + 8] = 32 * n; m[RSLOT0 + 16] = 1
    for i in range(n):
        for w, v in enumerate(rec(i, 0x31)): m[RSLOT0 + 32 + 32 * i + 8 * w] = v
    c_ = rl.receive() - 15
    assert sorted(rl.delivered()) == sorted((rec(i, 0x31)[0],) + rec(i, 0x31) for i in range(n))
    return c_, m[(1 << 48) - 4096]
last_final = {}
for n in (0, 1, 2, 3, 1023, 4094, 4095):
    c_mid, d_mid = h2s_fin(n, 100); c_last, d_last = h2s_fin(n, 1)
    assert d_mid == 0 and d_last == 1 and c_last - c_mid == cert.H2S_LAST + cert.H2S_DONE, (n, c_mid, c_last)
    last_final[n] = c_last
assert max(c_ - 6 * n for n, c_ in last_final.items()) == cert.H2S_PFIX + cert.H2S_LAST + cert.H2S_DONE

# 9i. init3 and tinit, the relay and home table initializers: INIT3_BKT and TINIT_ENTRY a word, plus the return.
im_ = M0(); c_init3 = a64.a64_run(code_rt, lbl_rt, "init3", {"x0": 5, "x1": 0x40000, "x2": 0x100, "x3": 64, "x30": SNT}, im_)
assert c_init3 == 5 * cert.INIT3_BKT + 1 and [im_[0x40000 + 8 * i] for i in range(5)] == [0x100 + 64 * i for i in range(5)]
tm_ = M0(); c_tinit = a64.a64_run(code_rt, lbl_rt, "tinit", {"x0": 5, "x1": 0x40000, "x2": 0x77, "x30": SNT}, tm_)
assert c_tinit == 5 * cert.TINIT_ENTRY + 1 and [tm_[0x40000 + 8 * i] for i in range(5)] == [0x77] * 5

