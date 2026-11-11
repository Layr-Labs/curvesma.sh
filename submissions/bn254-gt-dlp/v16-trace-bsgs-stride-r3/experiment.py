"""Supporting checks for the BN254 GT trace-recurrence BSGS claim (random giant stride, recheck against
all 12 Frobenius conjugates, table-lookup trace step, executed message transport). Standard library only.

Modes (input JSON on stdin):
  kernels    execute every priced program at production shapes against Python and big integers (field
             kernels, table step, table build, emit, walk loop, send, relay and home handlers with their END
             paths, flushes, reports, waits, radix passes, zero and prefix loops, merge) and print the counts.
  transport  a multi-node simulation: 2 walkers walk real BN254 traces with the table step, 4 relays and 16
             homes run the executed send and handler listings under two schedules; a planted target is found
             by a home and reported over the return route to walker 0. Every executed count is compared with
             the ledger's cost functions.
  algebra    constants from the shipped protected profile, orders, trace formula, disjoint Frobenius images
             of the baby range, the exact coverage inequality for m and its failure for m - 1.
  ledger     recompute every ledger integer and the score; compare every figure stated in submission.json,
             proof.md and machine.md; a mutated step cost must be detected. Exit status 1 on any mismatch.
  toy        the whole algorithm (random stride, hashed routing, caps, sort, merge, 12-conjugate recheck) in a
             BN-family toy field; solved must equal exact ground truth when no list overflows.
  all        everything above.
Every mode exits with status 1 when one of its checks fails. Toy and simulation runs are evidence at small
scale only. Nothing here runs a production attack.
"""
import json
import math
import random
import re
import sys
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
if (HERE / "files" / "progs.py").exists():          # verifier sandbox: the snapshot's files are under files/
    HERE = HERE / "files"
sys.path.insert(0, str(HERE))
from progs import *                                  # noqa: E402,F401,F403

C_EQ = 20000
LAM = P % R
PROD = dict(K=256, relays=1 << 36, homes_per=1 << 36, S1=28, A=36, B=28, C=56, F=28)
SIMR = dict(K=8, relays=4, homes_per=4, S1=62, A=63, B=62, C=63, F=62)
ST, RING, BLK, CURB, SBLKA = 1 << 33, 1 << 26, 1 << 50, 1 << 40, 1 << 45
PORTM, PORTD = ST + 2048, ST + 2064
BABY, GIANT, HCUR, SCRATCH, SORTB, CNTB = 1 << 60, (1 << 60) + (1 << 59), ST + 3072, ST + 3200, 1 << 61, 1 << 36
EXP_PRODUCTS = 508                                   # products in one exponentiation below 2^254 (2 x 254)
CHARGE = {"coin_word": 16, "delta_attempts": 128, "delta_compare": 64, "salt_words": 2,
          "setup_products": 4096, "setup_extra": 1 << 24, "table_init": 1 << 10, "fold_table": 1 << 10,
          "exp_products": EXP_PRODUCTS, "recheck_cap": 1 << 10, "recheck_products": 3 * EXP_PRODUCTS + 13, "recheck_extra": 1 << 12, "pair_expand": 64,
          "home_driver": 1 << 12, "node_fixed": 1 << 20, "output": 1 << 16, "wait_passes": 3}


def build(rt):
    pr = Program()
    gen_fp_mul(pr, 'fp_mul', 0); gen_fp_add(pr, 'fp_add', 0); gen_fp_sub(pr, 'fp_sub', 0)
    t12 = gen_f12_mul(pr)
    for ph in range(12):
        pr.label('tstep%d' % ph); gen_tabstep(pr, ph); pr.emit('ret')
    gen_tabbuild(pr)
    gen_emit(pr, rt['S1'], rt['A'], rt['B'], rt['C'])
    gen_transport(pr, rt['K'], rt['F'])
    gen_walk(pr)
    for word, shift in ((0, 0), (0, 32), (1, 0), (1, 32)):
        gen_radix(pr, word, shift)
    gen_table_ops(pr); gen_merge(pr)
    return pr, t12


def setv(m_, addr, x):
    for k, w in enumerate(limbs(x)):
        m_.st(addr + 8 * k, w)


def getv(m_, addr):
    return sum(m_.ld(addr + 8 * k) << (64 * k) for k in range(4))


def consts(m_):
    setv(m_, CONST, P)
    m_.st(CONST + 32, (-pow(P, -1, 1 << 64)) & M64)
    m_.r['x28'] = CONST


def walker_state(pr_, mem, rt, listbit, s1, s2, z, delta, nrelays, coeffs):
    M_ = Proc(pr_, mem, coeffs, 1 << 41)
    for k, v in enumerate([s1, s2, 0xbf58476d1ce4e5b9, 0x94d049bb133111eb, 0, 0, BLK, stride(rt['K']), CURB,
                           listbit, 40, rt['K'], PORTM, nrelays]):
        M_.st(ST + 8 * k, v)
    for k, c_ in enumerate(flat(z)):
        setv(M_, ST + PAR + 32 * k, c_)
    setv(M_, ST + DLT, delta)
    M_.st(ST + REP + NEXP, nrelays); M_.st(ST + REP + SBLK, SBLKA)
    consts(M_)
    M_.r.update({'x29': ST, 'x27': ST, 'x19': RING, 'x26': TB, 'x24': TB + 480, 'x23': 960, 'x25': FOLD})
    return M_


def relay_state(m_, rt, nh):
    for k, v in enumerate([BLK, stride(rt['K']), CURB, 40, rt['K'], 0, 0, 0, PORTD]):
        m_[ST + 8 * k] = v
    m_[ST + REP + NH], m_[ST + REP + NEXP], m_[ST + REP + RBLK], m_[ST + REP + MPORT] = nh, nh, SBLKA, PORTM


def home_state(m_, cap):
    for k, v in enumerate([BABY, GIANT - BABY, HCUR, 40, cap, 0, SCRATCH, 0, PORTD]):
        m_[ST + 8 * k] = v
    m_[ST + REP + RBLK], m_[ST + REP + MPORT] = SBLKA, PORTM


POISON = 0x0bad0000
WAKE = WakeMem.lo = ST + REP + READY                 # each node's wake word, alone in its 2^11-byte block


def dproc(pr_, mem, sp0):
    """A delivery processor: x27 (state block) and x1 (window base, used by h_rep and h_setup, which never
    write it) are set once, when the processor is created."""
    d = Proc(pr_, mem, None, sp0)
    d.r['x27'], d.r['x1'] = ST, WIN
    return d


def deliver_to(dst, label, smem, base, size):
    """The receiver's delivery processor runs its handler with its window on the sender's block. Nothing is
    reloaded between messages: x27 and x1 must still hold their creation values, and x1 is poisoned before
    the record handlers, which derive the window base themselves."""
    assert not getattr(dst, 'busy', False), 'a delivery processor runs one handler at a time'
    assert dst.r['x27'] == ST, 'state register changed between messages'
    if label in ('h_relay', 'h_home'):
        dst.r['x1'] = POISON
    else:
        assert dst.r['x1'] == WIN, 'window register changed between messages'
    dst.busy = True
    dst.win = (smem, base, size)
    n = dst.call(label)
    dst.win, dst.busy = None, False
    if label not in ('h_relay', 'h_home'):
        assert dst.r['x1'] == WIN
    return n


def run_merge(mg, A2, na, B2, nb):
    """Drives the merge listing over two sorted lists; equal-key runs are expanded by the pair routine.
    Returns (executed count, advancing steps, calls, equal-key exits, end extra, pairs)."""
    mg.r.update({'x1': A2, 'x2': B2, 'x3': A2 + 40 * na, 'x4': B2 + 40 * nb})
    total = steps = calls = eqx = 0
    pairs = []
    while mg.r['x1'] < A2 + 40 * na and mg.r['x2'] < B2 + 40 * nb:
        x1_, x2_ = mg.r['x1'], mg.r['x2']
        total += mg.call('merge'); calls += 1
        steps += (mg.r['x1'] - x1_ + mg.r['x2'] - x2_) // 40
        if mg.r['x9']:
            eqx += 1
            kk = (mg.ld(mg.r['x1']), mg.ld(mg.r['x1'] + 8))
            a0, b0 = mg.r['x1'], mg.r['x2']
            while mg.r['x1'] < A2 + 40 * na and (mg.ld(mg.r['x1']), mg.ld(mg.r['x1'] + 8)) == kk: mg.r['x1'] += 40
            while mg.r['x2'] < B2 + 40 * nb and (mg.ld(mg.r['x2']), mg.ld(mg.r['x2'] + 8)) == kk: mg.r['x2'] += 40
            pairs += [(a_, b_) for a_ in range(a0, mg.r['x1'], 40) for b_ in range(b0, mg.r['x2'], 40)]
    return total, steps, calls, eqx, pairs


# ---------------------------------------------------------------- kernels at production shapes

def kernels(seed):
    rng = random.Random(seed)
    TW = Tower(P, 9)
    g = TW.to_tower(G_PINNED)
    pr, t12 = build(PROD)
    rt = PROD
    sb = stride(rt['K'])
    out, bad = {}, 0
    RM = 1 << 256
    mont, unmont = (lambda x: x * RM % P), (lambda x: x * pow(RM, -1, P) % P)
    mc = Proc(pr); consts(mc)
    A_, B_, O_ = 1 << 20, (1 << 20) + 64, (1 << 20) + 128
    edge = [(2 * P - 1, 2 * P - 1), (2 * P - 1, 1), (0, 2 * P - 1), (P, P), (P - 1, 2 * P - 1)]
    for it in range(120):
        x, y = edge[it] if it < len(edge) else (rng.randrange(2 * P), rng.randrange(2 * P))
        setv(mc, A_, x); setv(mc, B_, y); mc.r['x0'], mc.r['x1'], mc.r['x2'] = O_, A_, B_
        out['fp_mul'] = mc.call('fp_mul'); bad += getv(mc, O_) != x * y * pow(RM, -1, P) % P
    YA, YB, YO = 1 << 24, (1 << 24) + 4096, (1 << 24) + 8192
    for _ in range(3):
        xa, xb = [rng.randrange(P) for _ in range(12)], [rng.randrange(P) for _ in range(12)]
        a, b = TW.to_tower(xa), TW.to_tower(xb)
        for k, v in enumerate(mem_order(a)): setv(mc, YA + 32 * k, mont(v))
        for k, v in enumerate(mem_order(b)): setv(mc, YB + 32 * k, mont(v))
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = YO, YA, YB
        out['f12_mul'] = mc.call('f12_mul')
        bad += TW.from_tower(from_mem([unmont(getv(mc, YO + 32 * k)) for k in range(12)])) != TW.schoolbook(xa, xb)
    # table step on the production baby walk (y = g, start 1) and giant walk (y = D = g^-delta, start z)
    delta = rng.randrange(1, R)
    D = TW.pow(g, (-delta) % R)
    z = TW.pow(g, rng.randrange(1, R))
    step_bad = 0
    for y0, yy in ((TW.one, g), (z, D)):
        c = charpoly(TW, yy, P)
        bad += not (c[0] == 1 and c[12] == 1 and all(c[i] == c[12 - i] for i in range(13)))
        tr, cur = [], y0
        for _ in range(12 + 26):
            tr.append(12 * cur[0][0][0] % P); cur = TW.mul(cur, yy)
        tm = TabMachine(pr, c); consts(tm)
        for j in range(12): setv(tm, RING + 32 * j, tr[j])
        tm.r.update({'x19': RING, 'x26': TB, 'x24': TB + 480, 'x23': 960, 'x25': FOLD})
        for n in range(26):
            out['table_step'] = tm.call('tstep%d' % (n % 12))
            step_bad += getv(tm, RING + 32 * (n % 12)) != tr[n + 12]
        bm = Proc(pr); consts(bm)                    # the table-build loop writes what the step reads
        for t in (0, 29):
            i, d = divmod(t, DIG)
            cst = c[i + 1] * (1 << (KB * d)) % P
            for k, w in enumerate(limbs(cst)): bm.r['x%d' % (9 + k)] = w
            for k in range(4): bm.r['x%d' % (5 + k)] = 0
            bm.r['x20'], bm.r['x21'] = TB + 960 + 32 * t, TB + 960 * 41 + 32 * t
            n_tb = bm.call('tabbuild')
            out['table_build_entry'] = (n_tb - 1) // 40
            bad += n_tb != 40 * out['table_build_entry'] + 1
            for v in range(1, 41):
                bad += getv(bm, TB + 960 * v + 32 * t) != getv(tm, TB + 960 * v + 32 * t)
    bad += step_bad
    # emit at production routing (K = 256, no block fills)
    s1, s2 = rng.getrandbits(64), rng.getrandbits(64)
    W = walker_state(pr, {}, rt, 1 << 63, s1, s2, z, delta, 3, None)
    for _ in range(16):
        t = rng.randrange(P)
        j = W.ld(ST + 32) | (W.ld(ST + 40) << 64)
        for k in range(4): W.r['x%d' % k] = (t >> (64 * k)) & M64
        w1, w2 = khash4(t, s1), khash4(t, s2)
        d1, d2 = w1 >> 28, ((w1 & ((1 << 28) - 1)) << 8) | (w2 >> 56)
        c0 = W.ld(CURB + 8 * d1)
        out['emit'] = W.call('emit')
        slot = BLK + sb * d1 + HDR + 40 * c0
        bad += [W.ld(slot + 8 * k) for k in range(5)] != [t & M64, (t >> 64) & M64, d2 | (1 << 63), j & M64, j >> 64]
        bad += W.ld(CURB + 8 * d1) != c0 + 1
    # walk loop at production shapes on the baby walk: 1, 12 and 24 steps (no block fills)
    cg = charpoly(TW, g, P)
    trg, cur = [], TW.one
    for _ in range(12 + 24):
        trg.append(12 * cur[0][0][0] % P); cur = TW.mul(cur, g)
    wl = {}
    for nsteps in (1, 12, 24):
        Wk = walker_state(pr, {}, rt, 0, s1, s2, z, delta, 3, cg)
        for j in range(12): setv(Wk, RING + 32 * j, trg[j])
        Wk.r['x22'] = nsteps
        wl[nsteps] = Wk.call('walk')
        bad += Wk.ld(ST + 32) != nsteps
        bad += any(getv(Wk, RING + 32 * ((n_ + 12) % 12)) != trg[n_ + 12] for n_ in range(nsteps - min(nsteps, 12), nsteps))
    out['walk_step'] = wl[1] - 3
    bad += wl[12] != 3 + 12 * out['walk_step'] or wl[24] != 3 + 24 * out['walk_step'] + 1
    Wi = walker_state(pr, {}, rt, 0, s1, s2, z, delta, 3, cg)
    for j in range(12): setv(Wi, RING + 32 * j, trg[j])
    out['walk_init'] = Wi.call('walk_init')
    bad += Wi.ld(ST + 32) != 12
    # hop 1 and hop 2 at production shapes: a full walker block -> relay handler -> a full home block -> home
    hmem, rmem = {}, {}
    relay_state(rmem, rt, 3)
    home_state(hmem, 1 << 20)
    Rd = dproc(pr, rmem, 1 << 42)
    Hd = dproc(pr, hmem, 1 << 42)
    hop1, hop2 = [], []

    def r_hook(src, addr):
        base = addr - 8
        hop2.append((base, src.mem[base], src.mem.get(base + 16, 0), deliver_to(Hd, 'h_home', src.mem, base, sb)))
    Rd.on_stlr = r_hook

    def w_hook(src, addr):
        base = addr - 8
        hop1.append((base, src.mem[base], deliver_to(Rd, 'h_relay', src.mem, base, sb)))
    Wf = walker_state(pr, {}, rt, 0, s1, s2, z, delta, 3, None)
    Wf.on_stlr = w_hook
    d2s = set()
    while len(d2s) < rt['K']:                        # distinct homes, so every record's slot is known
        d2s.add(rng.getrandbits(36))
    d2s = sorted(d2s, key=lambda _: rng.random())
    hot = d2s[0]
    recs = []
    for q in range(rt['K']):                         # one full walker block for relay d1 = 5
        t = rng.randrange(P)
        recs.append([t & M64, (t >> 64) & M64, d2s[q] | (rng.getrandbits(1) << 63), q, 0])
        for k in range(5): Wf.st(BLK + sb * 5 + HDR + 40 * q + 8 * k, recs[-1][k])
    rmem[CURB + 8 * hot] = rt['K'] - 1               # the relay's block for `hot` already holds K - 1 records
    for q in range(rt['K'] - 1):
        for k in range(5): rmem[BLK + sb * hot + HDR + 40 * q + 8 * k] = q + 7 * k
    Wf.st(CURB + 8 * 5, rt['K'])
    Wf.r.update({'x6': 5, 'x16': sb, 'x17': rt['K']})
    out['callsend_w'] = Wf.call('callsend_w')
    bad += len(hop1) != 1 or Wf.ld(CURB + 8 * 5) != 0 or len(hop2) != 1 or (hop1 and hop1[0][1] != rt['K'])
    bad += (not hop2) or hop2[0][1] != rt['K'] or hop2[0][2] != 0 or rmem.get(CURB + 8 * hot) != 0
    for r_ in recs[1:]:
        d2 = r_[2] & ((1 << 36) - 1)
        slot = BLK + sb * d2 + HDR + 40 * (rmem[CURB + 8 * d2] - 1)
        bad += [rmem.get(slot + 8 * k, 0) for k in range(5)] != r_
    bad += hmem.get(HCUR, 0) + hmem.get(HCUR + 8, 0) != rt['K']
    out['home_msg_full'] = hop2[0][3] if hop2 else 0
    smem = {}
    for k, v in enumerate([1, 2, 9, 3, 4]): smem[(1 << 30) + HDR + 8 * k] = v
    for k, c_ in enumerate(flat(z)):
        for w in range(4): smem[(1 << 30) + HDR + 40 * rt['K'] + 32 * k + 8 * w] = (c_ >> (64 * w)) & M64
    for w in range(4): smem[(1 << 30) + HDR + 40 * rt['K'] + 384 + 8 * w] = (delta >> (64 * w)) & M64

    def msg(nrec, end):
        smem[1 << 30], smem[(1 << 30) + 16] = nrec, end
        return deliver_to(Rd, 'h_relay', smem, 1 << 30, sb)
    hop2.clear()
    out['relay_msg_base'] = msg(0, 0)
    out['relay_rec'] = msg(1, 0) - out['relay_msg_base']
    out['relay_full_extra'] = hop1[0][2] - out['relay_msg_base'] - rt['K'] * out['relay_rec']
    rmem[ST + REP + ENDC] = 0
    out['relay_end'] = msg(0, 1) - out['relay_msg_base']
    e2 = msg(0, 1) - out['relay_msg_base'] - out['relay_end']
    out['relay_flush_home'] = (e2 - 2) // 3
    bad += e2 != 2 + 3 * out['relay_flush_home'] or len(hop2) != 3 or any(h[2] != 1 for h in hop2)
    bad += [getv(Hd, ST + PAR + 32 * k) for k in range(12)] != flat(z) or getv(Hd, ST + DLT) != delta
    bad += hmem.get(ST + REP + READY) != 1
    # home handler: unit costs; planted overflow at cap 12, none at cap 64
    hm2 = {}
    H2 = dproc(pr, hm2, 1 << 42)
    hrec = [[rng.getrandbits(64), rng.getrandbits(64), rng.getrandbits(36) | (rng.getrandbits(1) << 63),
             rng.getrandbits(64), rng.getrandbits(62)] for _ in range(40)]
    nb = sum(1 for r_ in hrec if r_[2] >> 63 == 0)
    hs = {}
    for i, r_ in enumerate(hrec):
        for k in range(5): hs[(1 << 30) + HDR + 40 * i + 8 * k] = r_[k]
    for capv in (64, 12):
        hm2.clear(); home_state(hm2, capv)
        hs[1 << 30], hs[(1 << 30) + 16] = 40, 0
        n40 = deliver_to(H2, 'h_home', hs, 1 << 30, sb)
        bad += hm2.get(ST + 40, 0) != int(max(nb, 40 - nb) > capv)
        if capv == 64:
            got = sorted(tuple(hm2.get(base + 40 * q + 8 * k, 0) for k in range(5))
                         for base, cnt in ((BABY, hm2.get(HCUR, 0)), (GIANT, hm2.get(HCUR + 8, 0))) for q in range(cnt))
            bad += got != sorted(tuple(r_) for r_ in hrec)
            n40_ok = n40
    hs[1 << 30] = 0
    out['home_msg_base'] = deliver_to(H2, 'h_home', hs, 1 << 30, sb)
    out['home_rec'] = (n40_ok - out['home_msg_base']) // 40
    bad += n40_ok != out['home_msg_base'] + 40 * out['home_rec']
    bad += out['home_msg_full'] != out['home_msg_base'] + rt['K'] * out['home_rec']
    hs[(1 << 30) + 16] = 1
    for k, c_ in enumerate(flat(z)):
        for w in range(4): hs[(1 << 30) + HDR + 40 * rt['K'] + 32 * k + 8 * w] = (c_ >> (64 * w)) & M64
    out['home_end'] = deliver_to(H2, 'h_home', hs, 1 << 30, sb) - out['home_msg_base']
    bad += hm2.get(ST + REP + READY) != 1 or [getv(H2, ST + PAR + 32 * k) for k in range(12)] != flat(z)
    # reports: home -> relay handler (not last, then last), relay_up -> walker 0 handler, output; waits
    rp_mem = {}
    relay_state(rp_mem, rt, 2)
    R2 = dproc(pr, rp_mem, 1 << 43)
    xr = rng.randrange(R)
    hmm = {}
    home_state(hmm, 64)
    Hm = Proc(pr, hmm, None, 1 << 41); Hm.r['x27'] = ST
    reps = []
    Hm.on_stlr = lambda src, addr: reps.append(deliver_to(R2, 'h_rep', src.mem, addr - 8, HDR + 40))
    for found in (1, 0):
        for k in range(4): Hm.r['x%d' % k] = (xr >> (64 * k)) & M64 if found else 0
        Hm.r['x4'] = found
        out['report'] = Hm.call('report')
    out['rep_handler'], out['rep_handler_last'] = reps
    bad += getv(R2, ST + REP + XS) != xr or rp_mem.get(ST + REP + FND) != 1 or rp_mem.get(ST + REP + READY) != 1
    Rm = Proc(pr, rp_mem, None, 1 << 41); Rm.r['x27'] = ST
    w0 = {ST + REP + NEXP: 1}
    W0d = dproc(pr, w0, 1 << 42)
    Rm.on_stlr = lambda src, addr: deliver_to(W0d, 'h_rep', src.mem, addr - 8, HDR + 40)
    out['relay_up'] = Rm.call('relay_up')
    W0m = Proc(pr, w0, None, 1 << 41); W0m.r['x27'] = ST
    out['output'] = W0m.call('output')
    bad += sum(W0m.r['x%d' % k] << (64 * k) for k in range(4)) != xr or W0m.r['x4'] != 1
    Wt = Proc(pr, {ST: 1}); Wt.r['x1'] = ST
    out['wait'] = Wt.call('wait')

    class Stale(Proc):                               # the capped path: the word is still 0 on the first passes
        def ld_acq(self, ad):
            self.k_ = getattr(self, 'k_', 0) + 1
            return int(self.k_ >= CHARGE["wait_passes"])
    out['wait_cap'] = Stale(pr).call('wait')
    # walker flush (3 and 2 relays), setup message
    fl_msgs = []
    Wfl = walker_state(pr, {}, rt, 0, s1, s2, z, delta, 3, None)
    Wfl.on_stlr = lambda src, addr: fl_msgs.append((src.mem[addr - 8], src.mem[addr + 8]))
    n3 = Wfl.call('flush_w')
    Wfl2 = walker_state(pr, {}, rt, 0, s1, s2, z, delta, 2, None)
    Wfl2.on_stlr = lambda src, addr: None
    n2 = Wfl2.call('flush_w')
    out['flush_relay'] = n3 - n2
    out['flush_fixed'] = n2 - 2 * out['flush_relay']
    bad += len(fl_msgs) != 3 or any(f[1] != 1 for f in fl_msgs)
    w1m = {}
    W1d = dproc(pr, w1m, 1 << 42)
    W0s = walker_state(pr, {}, rt, 0, s1, s2, z, delta, 3, None)
    W0s.on_stlr = lambda src, addr: out.__setitem__('h_setup', deliver_to(W1d, 'h_setup', src.mem, addr - 8, HDR + 48))
    out['setup_send'] = W0s.call('setup_send')
    bad += getv(W1d, ST + DLT) != delta or w1m.get(ST) != s1 or w1m.get(ST + 8) != s2 or w1m.get(ST + REP + READY) != 1
    # radix passes: full-width keys (counter prefix in Python), then small digits with the zero and prefix loops
    IN, DST = 1 << 50, 1 << 52
    for small in (False, True):
        keys = [((rng.randrange(64) | rng.randrange(64) << 32) if small else rng.getrandbits(64),
                 (rng.randrange(64) | rng.randrange(64) << 32) if small else rng.getrandbits(64)) for _ in range(32)]
        rs = Proc(pr)
        for i, (k0, k1) in enumerate(keys):
            for k, v in enumerate([k0, k1, i, 0, 0]): rs.st(IN + 40 * i + 8 * k, v)
        src, dst = IN, DST
        for (wd, sh) in ((0, 0), (0, 32), (1, 0), (1, 32)):
            if small:
                rs.r.update({'x19': CNTB, 'x24': CNTB + 8 * 64})
                out['zero_call_64'] = rs.call('zero')
            else:
                rs.mem = {a_: v for a_, v in rs.mem.items() if not (CNTB <= a_ < CNTB + (1 << 35))}
            rs.r.update({'x19': CNTB, 'x24': src + 40 * len(keys), 'x1': src})
            out['radix_count_pass_32'] = rs.call('count%d%d' % (wd, sh))
            if small:
                rs.r.update({'x19': CNTB, 'x24': CNTB + 8 * 64, 'x3': 0})
                out['prefix_call_64'] = rs.call('prefix')
            else:
                run = 0
                for d_ in sorted({(rs.ld(src + 40 * i + 8 * wd) >> sh) & 0xffffffff for i in range(len(keys))}):
                    cnt = rs.ld(CNTB + 8 * d_); rs.st(CNTB + 8 * d_, run); run += cnt
            rs.r.update({'x19': CNTB, 'x20': dst, 'x22': 40, 'x24': src + 40 * len(keys), 'x1': src})
            out['radix_scatter_pass_32'] = rs.call('scatter%d%d' % (wd, sh))
            src, dst = dst, src
        bad += [(rs.ld(src + 40 * i), rs.ld(src + 40 * i + 8)) for i in range(len(keys))] != sorted(keys, key=lambda kk: (kk[1], kk[0]))
    out['radix_count'] = (out['radix_count_pass_32'] - 1) // 32
    out['radix_scatter'] = (out['radix_scatter_pass_32'] - 1) // 32
    out['zero_per_4'] = (out['zero_call_64'] - 1) // 16
    out['prefix_entry'] = (out['prefix_call_64'] - 1) // 64
    bad += out['radix_count_pass_32'] != 32 * out['radix_count'] + 1 or out['radix_scatter_pass_32'] != 32 * out['radix_scatter'] + 1
    bad += out['zero_call_64'] != 16 * out['zero_per_4'] + 1 or out['prefix_call_64'] != 64 * out['prefix_entry'] + 1
    # merge over two sorted lists drawn from a common key pool
    mg = Proc(pr)
    A2, B2 = 1 << 45, 1 << 46
    pool = [(rng.getrandbits(64), rng.getrandbits(64)) for _ in range(10)]
    ka = sorted((rng.choice(pool) for _ in range(24)), key=lambda kk: (kk[1], kk[0]))
    kb = sorted((rng.choice(pool) for _ in range(24)), key=lambda kk: (kk[1], kk[0]))
    for i, kk in enumerate(ka): mg.st(A2 + 40 * i, kk[0]); mg.st(A2 + 40 * i + 8, kk[1])
    for i, kk in enumerate(kb): mg.st(B2 + 40 * i, kk[0]); mg.st(B2 + 40 * i + 8, kk[1])
    total, steps, calls, eqx, pairs = run_merge(mg, A2, 24, B2, 24)
    seen = {(mg.ld(a_), mg.ld(a_ + 8)) for a_, _ in pairs}
    bad += seen != (set(ka) & set(kb)) or any((mg.ld(a_), mg.ld(a_ + 8)) != (mg.ld(b_), mg.ld(b_ + 8)) for a_, b_ in pairs)
    out['merge_step'], out['merge_exit'] = 19, 11
    bad += not (steps * 19 <= total <= steps * 19 + calls * 11)
    s0 = pr.labels['f12_mul']
    s1_ = min(v for v in pr.labels.values() if v > s0)
    out['f12_fp_mul_calls'] = sum(1 for op, a in pr.code[s0:s1_] if op == 'bl' and a[0] == 'fp_mul')
    out.update({"mismatches": bad, "table_step_mismatches": step_bad, "f12_scratch_slots": t12.peak,
                "program_instructions": len(pr.code)})
    return out


# ---------------------------------------------------------------- cost functions (shared by ledger and simulation)

def walk_cost(k, c):
    return 3 + c * k['walk_step'] + (c - 1) // 12


def role_costs(k, ev, wait):
    """Executed instructions per processor kind, from unit counts k and event counts ev."""
    Rn, Hn = ev['R'], ev['H']
    return {
        "walker_main": sum(walk_cost(k, c) for c in ev['walk_calls']) + ev['walk_inits'] * k['walk_init']
        + ev['F1'] * (k['callsend_w'] - 1) + sum(k['flush_fixed'] + f * k['flush_relay'] for f in ev['flushes'])
        + k['setup_send'] + 2 * wait + k['output'],
        "walker_delivery": k['h_setup'] + Rn * k['rep_handler'] + (k['rep_handler_last'] - k['rep_handler']),
        "relay_delivery_1": (ev['F1'] + 2 * Rn) * k['relay_msg_base'] + ev['N'] * k['relay_rec']
        + ev['F2'] * k['relay_full_extra'] + 2 * Rn * k['relay_end'] + Rn * (2 + ev['Hr'] * k['relay_flush_home']),
        "relay_delivery_2": Hn * k['rep_handler'] + Rn * (k['rep_handler_last'] - k['rep_handler']),
        "relay_main": Rn * (wait + k['relay_up']),
        "home_delivery": (ev['F2'] + Hn) * k['home_msg_base'] + ev['N'] * k['home_rec'] + Hn * k['home_end'],
        "home_main_transport": Hn * (wait + k['report'])}


# ---------------------------------------------------------------- multi-node simulation

def transport(seed, k):
    """2 walkers -> 4 relays -> 16 homes with the same generated listings (routing shifts select 2 + 2 bits
    of the same hash words; K = 8 records per block). Walkers run the table step on real BN254 traces."""
    TW = Tower(P, 9)
    g = TW.to_tower(G_PINNED)
    rng = random.Random(1000003 + seed)
    rt = SIMR
    pr, _ = build(rt)
    K, NR, NHp = rt['K'], rt['relays'], rt['homes_per']
    sb = stride(K)
    ms = 60 + rng.randrange(12)                       # baby j = 0..ms, giant i = 0..ms-1
    delta, s1, s2 = rng.randrange(1, R), rng.getrandbits(64), rng.getrandbits(64)
    i0, j0, k0 = rng.randrange(ms), rng.randrange(1, ms + 1), rng.randrange(12)
    xt = (delta * i0 + j0 * pow(LAM, k0, R)) % R      # planted target inside the window
    z = TW.pow(g, xt)
    cg = charpoly(TW, g, P)
    D = TW.pow(g, (-delta) % R)
    baby, cur = [], TW.one
    for j in range(ms + 1):
        baby.append(12 * cur[0][0][0] % P); cur = TW.mul(cur, g)
    giant, cur = [], z
    for i in range(ms):
        giant.append(12 * cur[0][0][0] % P); cur = TW.mul(cur, D)
    expected = {}
    for lst, trs in ((0, baby), (1, giant)):
        for idx, t in enumerate(trs):
            w1, w2 = khash4(t, s1), khash4(t, s2)
            expected.setdefault((w1 >> 62, ((w1 & 1) << 1) | (w2 >> 63)), []).append((t & ((1 << 128) - 1), lst, idx))
    runs = []
    for sched in (0, 1):
        log = {"hop1": 0, "hop1_end": 0, "hop2": 0, "hop2_end": 0, "setup": 0, "rep2": 0, "rep1": 0}
        wm = [WakeMem(), WakeMem()]
        WM = [walker_state(pr, wm[0], rt, 0, s1, s2, z, delta, NR, cg),
              walker_state(pr, wm[1], rt, 1 << 63, 0, 0, z, 0, NR, None)]
        WD = [dproc(pr, wm[w], 1 << 42) for w in range(2)]
        rmem = [WakeMem() for _ in range(NR)]
        for m_ in rmem: relay_state(m_, rt, NHp)
        RD1 = [dproc(pr, rmem[i], 1 << 42) for i in range(NR)]
        RD2 = [dproc(pr, rmem[i], 1 << 43) for i in range(NR)]
        RM = [Proc(pr, rmem[i], None, 1 << 41) for i in range(NR)]
        hmem = [[WakeMem() for _ in range(NHp)] for _ in range(NR)]
        cap = 2 * -(-(ms + 1) // (NR * NHp)) + 64
        for row in hmem:
            for m_ in row: home_state(m_, cap)
        HD = [[dproc(pr, hmem[i][j], 1 << 42) for j in range(NHp)] for i in range(NR)]
        HM = [[Proc(pr, hmem[i][j], None, 1 << 41) for j in range(NHp)] for i in range(NR)]
        for p_ in RM + [h_ for row in HM for h_ in row]:
            p_.r['x27'] = ST

        def w_hook(src, addr, w):
            base = addr - 8
            if w == 0 and base == SBLKA:
                log["setup"] += 1
                return deliver_to(WD[1], 'h_setup', src.mem, base, HDR + 48)
            d1, off = divmod(base - BLK, sb)
            assert off == 0 and 0 <= d1 < NR
            log["hop1_end" if src.mem.get(base + 16) else "hop1"] += 1
            deliver_to(RD1[d1], 'h_relay', src.mem, base, sb)

        def r_hook(src, addr, ri):
            base = addr - 8
            d2, off = divmod(base - BLK, sb)
            assert off == 0 and 0 <= d2 < NHp
            log["hop2_end" if src.mem.get(base + 16) else "hop2"] += 1
            deliver_to(HD[ri][d2], 'h_home', src.mem, base, sb)

        def h_hook(src, addr, ri):
            assert addr - 8 == SBLKA
            log["rep2"] += 1
            deliver_to(RD2[ri], 'h_rep', src.mem, SBLKA, HDR + 40)

        def m_hook(src, addr):
            assert addr - 8 == SBLKA
            log["rep1"] += 1
            deliver_to(WD[0], 'h_rep', src.mem, SBLKA, HDR + 40)
        for w in range(2):
            WM[w].on_stlr = (lambda w_: (lambda s_, a_: w_hook(s_, a_, w_)))(w)
        for i in range(NR):
            RD1[i].on_stlr = (lambda i_: (lambda s_, a_: r_hook(s_, a_, i_)))(i)
            RM[i].on_stlr = m_hook
            for j in range(NHp):
                HM[i][j].on_stlr = (lambda i_: (lambda s_, a_: h_hook(s_, a_, i_)))(i)
        # setup: walker 0 sends (delta, s1, s2); walker 1 waits, then derives its walk from its own memory
        WM[0].call('setup_send')
        WM[1].r['x1'] = ST + REP + READY
        WM[1].call('wait')
        d_recv = getv(WM[1], ST + DLT)
        D_ = TW.pow(g, (-d_recv) % R)
        WM[1].c = charpoly(TW, D_, P)
        cur = unflat([getv(WM[1], ST + PAR + 32 * q) for q in range(12)])
        for j in range(12):
            setv(WM[1], RING + 32 * j, 12 * cur[0][0][0] % P); cur = TW.mul(cur, D_)
        for j in range(12):
            setv(WM[0], RING + 32 * j, baby[j])
        left = [ms + 1 - 12, ms - 12]
        calls = []

        def walk_chunk(w, c):
            WM[w].r['x22'] = c; WM[w].call('walk'); left[w] -= c; calls.append(c)
        order = [0, 1] if sched == 0 else [1, 0]
        for w in order:
            WM[w].call('walk_init')
        if sched == 0:                                  # alternate in chunks of 12 steps, then both flushes
            while left[0] or left[1]:
                for w in (0, 1):
                    if left[w]:
                        walk_chunk(w, min(12, left[w]))
            WM[0].call('flush_w'); WM[1].call('flush_w')
        else:                                           # walker 1 walks and ends first, then walker 0
            walk_chunk(1, left[1]); WM[1].call('flush_w')
            walk_chunk(0, left[0]); WM[0].call('flush_w')
        # homes: wait for END, sort (Python stand-in for the radix listings run in mode kernels), merge, recheck
        got, pairs_seen, flags, params_ok = {}, 0, 0, True
        merge_exec = merge_bound = 0
        homes = [(i, j) for i in range(NR) for j in range(NHp)]
        for (i, j) in (homes if sched == 0 else homes[::-1]):
            hm, H_ = hmem[i][j], HM[i][j]
            H_.r['x1'] = ST + REP + READY
            H_.call('wait')
            flags += hm.get(ST + 40, 0)
            zh = unflat([getv(H_, ST + PAR + 32 * q) for q in range(12)])
            dh = getv(H_, ST + DLT)
            params_ok &= zh == z and dh == delta
            lists = []
            for lb, base in ((0, BABY), (1, GIANT)):
                rr = [tuple(hm.get(base + 40 * q + 8 * w_, 0) for w_ in range(5)) for q in range(hm.get(HCUR + 8 * lb, 0))]
                for r_ in rr:
                    got.setdefault((i, j), []).append(((r_[1] << 64) | r_[0], lb, r_[3] | (r_[4] << 64)))
                rr.sort(key=lambda r_: (r_[1], r_[0]))
                for q, r_ in enumerate(rr):
                    for w_ in range(5): hm[SORTB + (lb << 59) + 40 * q + 8 * w_] = r_[w_]
                lists.append(len(rr))
            found, xf = 0, 0
            if lists[0] and lists[1]:
                tot, steps, ncall, _, pairs = run_merge(H_, SORTB, lists[0], SORTB + (1 << 59), lists[1])
                merge_exec += tot; merge_bound += 19 * steps + 11 * ncall
                bad_steps = steps > lists[0] + lists[1]
                assert not bad_steps
                for a_, b_ in pairs:
                    pairs_seen += 1
                    jj, ii = hm[a_ + 24] | (hm[a_ + 32] << 64), hm[b_ + 24] | (hm[b_ + 32] << 64)
                    u, yv = TW.pow(g, jj), TW.mul(zh, TW.pow(TW.pow(g, (-dh) % R), ii))
                    for kk in range(12):
                        if yv == TW.frob_k(u, kk):
                            xx = (dh * ii + jj * pow(LAM, kk, R)) % R
                            if TW.pow(g, xx) == zh and not found:
                                found, xf = 1, xx
            for q in range(4): H_.r['x%d' % q] = (xf >> (64 * q)) & M64
            H_.r['x4'] = found
            H_.call('report')
        for i in (range(NR) if sched == 0 else reversed(range(NR))):
            RM[i].r['x1'] = ST + REP + READY
            RM[i].call('wait'); RM[i].call('relay_up')
        W0o = Proc(pr, wm[0], None, 1 << 41)
        W0o.r['x27'] = ST
        W0o.r['x1'] = ST + REP + READY
        W0o.call('wait'); W0o.call('output')
        xo = sum(W0o.r['x%d' % q] << (64 * q) for q in range(4))
        ok_homes = all(sorted(got.get(h_, [])) == sorted(v) for h_, v in expected.items()) and set(got) <= set(expected)
        ev = {'R': NR, 'H': NR * NHp, 'Hr': NHp, 'N': 2 * ms + 1, 'F1': log['hop1'], 'F2': log['hop2'],
              'walk_calls': calls, 'walk_inits': 2, 'flushes': [NR, NR]}
        pred = role_costs(k, ev, k['wait'])
        wake1 = all(getattr(m_, 'wake', 0) == [WAKE] for m_ in wm + rmem + [h_ for row in hmem for h_ in row])
        execd = {"walker_main": WM[0].count + WM[1].count + W0o.count,
                 "walker_delivery": WD[0].count + WD[1].count,
                 "relay_delivery_1": sum(p_.count for p_ in RD1), "relay_delivery_2": sum(p_.count for p_ in RD2),
                 "relay_main": sum(p_.count for p_ in RM), "home_delivery": sum(p_.count for row in HD for p_ in row),
                 "home_main_transport": sum(p_.count for row in HM for p_ in row) - merge_exec}
        runs.append({"schedule": ["chunks of 12 steps alternate", "walker 1 walks and ends first"][sched],
                     "messages": log, "every_record_once_at_its_home": ok_homes, "parameters_at_every_home": params_ok,
                     "overflow_flags": flags, "equal_key_pairs": pairs_seen, "output_found": W0o.r['x4'],
                     "output_equals_target": xo == xt and W0o.r['x4'] == 1,
                     "executed_by_role": execd, "cost_functions_by_role": pred,
                     "executed_equals_cost_functions": execd == pred, "each_wake_block_stored_once": wake1,
                     "merge_within_bound": merge_exec <= merge_bound,
                     "hop1_full_within_bound": log['hop1'] <= (ms + 1) // K + ms // K,
                     "hop2_full_within_bound": log['hop2'] <= (2 * ms + 1) // K,
                     "contents": sorted((h_, tuple(sorted(v))) for h_, v in got.items())})
    same = runs[0].pop("contents") == runs[1].pop("contents")
    ok = same and all(r_[f] for r_ in runs for f in ("every_record_once_at_its_home", "parameters_at_every_home",
                                                      "output_equals_target", "executed_equals_cost_functions", "each_wake_block_stored_once",
                                                      "merge_within_bound", "hop1_full_within_bound",
                                                      "hop2_full_within_bound")) and all(r_["overflow_flags"] == 0 for r_ in runs)
    return {"walkers": 2, "relays": NR, "homes": NR * NHp, "records_per_block": K, "baby_records": ms + 1,
            "giant_records": ms, "planted": {"i": i0, "j": j0, "k": k0}, "runs": runs,
            "both_schedules_same_home_contents": same, "transport_ok": ok}


# ---------------------------------------------------------------- algebra on the pinned constants

def algebra(seed):
    rng = random.Random(seed)
    TW = Tower(P, 9)
    g = TW.to_tower(G_PINNED)
    m, h = stride_params(R)
    ok_trace = True
    for _ in range(2):
        y = TW.pow(g, rng.randrange(1, R))
        s, cur = None, y
        for _ in range(12):
            s = cur if s is None else tuple(tuple(TW.f2add(a_, b_) for a_, b_ in zip(u_, v_)) for u_, v_ in zip(s, cur))
            cur = TW.frob(cur)
        fl = [c_ for i in range(2) for j in range(3) for c_ in s[i][j]]
        ok_trace &= fl[0] == 12 * y[0][0][0] % P and not any(fl[1:])
    y = TW.pow(g, rng.randrange(1, R))
    frob_is_lambda = all(TW.frob_k(y, kk) == TW.pow(y, pow(LAM, kk, R)) for kk in (1, 5, 6, 11))
    sv = [lagrange_min((1, pow(LAM, kk, R)), (0, R)) for kk in range(1, 12) if kk != 6]
    min_sup_lb = min(math.isqrt(v[0] ** 2 + v[1] ** 2) for v in sv) // 2      # sup >= Euclid / sqrt 2 > Euclid / 2
    return {"p_r_prime": is_prime(P) and is_prime(R),
            "ord_r_p_is_12": all(pow(P, kk, R) != 1 for kk in range(1, 12)) and pow(P, 12, R) == 1,
            "lambda6_is_minus1": pow(LAM, 6, R) == R - 1, "g_order_r": TW.pow(g, R) == TW.one and g != TW.one,
            "trace_is_12_Re_e0": ok_trace, "frobenius_acts_as_lambda_power": frob_is_lambda,
            "m": str(m), "h": str(h), "m_least": covers(m, R) and not covers(m - 1, R), "m_lt_r_over_2": 2 * m < R,
            "nonoverlap_sup_lower_bound_log2_floor": math.floor(math.log2(min_sup_lb) * 1000) / 1000,
            "m_log2_ceil": math.ceil(math.log2(m) * 1000) / 1000,
            "frobenius_images_disjoint": min_sup_lb > m}


PROFILE_SHA256 = "b9c9b8f4a214f4c3327f27913caa78efa7bf260e7244e7dc6ba165e2172ce708"   # tracks/bn254-gt-dlp/target-profile.json


def profile_check():
    """Load the shipped copy of the protected target profile and assert every constant used here."""
    import hashlib
    f = HERE / "target-profile.json"
    if not f.exists():
        return {"profile_found": False}
    raw = f.read_bytes()
    prof = json.loads(raw)
    gen = [int(c, 16) for c in prof["generator_g"]["coefficients"]]
    return {"profile_found": True, "sha256_matches_protected": hashlib.sha256(raw).hexdigest() == PROFILE_SHA256,
            "p_matches": int(prof["field_modulus"], 16) == P, "r_matches": int(prof["subgroup_order"], 16) == R,
            "x_matches": int(prof["parameter_x"], 16) == U, "g_matches": gen == G_PINNED,
            "fp12_basis_matches": "w^12-18*w^6+82" in prof["field_Fp12"]}


def isa_negative():
    """Every immediate and label form the generator emits must reject out-of-range values (and x18)."""
    probes = [("add", ("x0", "x1", 4096)), ("add", ("x0", "x1", -1)), ("adds_i", ("x0", "x1", 4096)),
              ("subs_i", ("x0", "x1", 4096)), ("ldp", ("x0", "x1", "x2", 512)), ("ldp", ("x0", "x1", "x2", -520)),
              ("ldp", ("x0", "x1", "x2", 12)), ("stp", ("x0", "x1", "x2", 508)), ("ldr", ("x0", "x1", 32768)),
              ("ldr", ("x0", "x1", 4)), ("str", ("x0", "x1", -8)), ("movi", ("x0", 65536)), ("movi", ("x0", -1)),
              ("lsl_i", ("x0", "x1", 64)), ("lsr", ("x0", "x1", -1)), ("eor_lsr", ("x0", "x1", "x2", 64)),
              ("add_lsl", ("x0", "x1", "x2", 64)), ("extr", ("x0", "x1", "x2", 64)), ("extr", ("x0", "x1", "x2", -1)),
              ("and_mask", ("x0", "x1", 0)), ("and_mask", ("x0", "x1", 64)), ("bc", ("gt", "l")), ("b", (7,)),
              ("cbz", ("x0", 7)), ("cbnz_l", ("x0", 7)), ("bl", (7,)), ("movz_lsl", ("x0", 65536, 32)),
              ("movz_lsl", ("x0", 1, 8)), ("sub_sp", (4097,)), ("sub_sp", (1 << 25,)), ("add_sp", (4097,)),
              ("sub_sp", (-4096,)), ("sub_sp", (-(1 << 24),)), ("add_sp", (-(1 << 40),)), ("add_sp", (-(1 << 64),)), ("csel", ("x0", "x1", "x2", "gt")), ("cset", ("x0", "lt")),
              ("stp_pre", ("x29", "x30", "sp", -520)), ("ldp_post", ("x29", "x30", "sp", 512)),
              ("lsr", ("x0", "x1", 1 << 64)), ("eor_lsr", ("x0", "x1", "x2", -1)), ("lsl_i", ("x0", "x1", -1)),
              ("mov", ("x18", "x1"))]
    good = [("add", ("x0", "x1", 4095)), ("adds_i", ("x0", "x1", 1)), ("subs_i", ("x0", "x1", 4095)),
            ("ldp", ("x0", "x1", "x2", 504)), ("ldp", ("x0", "x1", "x2", -512)), ("ldr", ("x0", "x1", 32760)),
            ("movi", ("x0", 65535)), ("lsl_i", ("x0", "x1", 63)), ("eor_lsr", ("x0", "x1", "x2", 31)),
            ("extr", ("x0", "x1", "x2", 63)), ("and_mask", ("x0", "x1", 63)), ("bc", ("ne", "l")),
            ("cbnz_l", ("x0", "l")), ("bl", ("l",)), ("movz_lsl", ("x0", 16384, 32)), ("sub_sp", (4096,)),
            ("add_sp", (4096,)), ("sub_sp", (4095 << 12,)), ("stp_pre", ("x29", "x30", "sp", -16)),
            ("ldp_post", ("x29", "x30", "sp", 16))]
    rejected = accepted = 0
    for op, a in probes:
        try:
            Program().emit(op, *a)
        except AssertionError:
            rejected += 1
    for op, a in good:
        try:
            Program().emit(op, *a); accepted += 1
        except AssertionError:
            pass
    return {"out_of_range_probes": len(probes), "rejected": rejected, "in_range_probes": len(good), "accepted": accepted,
            "all_correct": rejected == len(probes) and accepted == len(good)}


# ---------------------------------------------------------------- ledger

def lg(v):
    return Decimal(v).ln() / Decimal(2).ln()


def up(v, d):
    return (lg(v) * 10 ** d).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** d)


def ledger(k):
    getcontext().prec = 120
    m, h = stride_params(R)
    Nb, Ng = m + 1, h
    N = Nb + Ng
    K, Rn, Hr = PROD['K'], PROD['relays'], PROD['homes_per']
    Hn = Rn * Hr
    F1, F2 = Nb // K + Ng // K, N // K
    ev = {'R': Rn, 'H': Hn, 'Hr': Hr, 'N': N, 'F1': F1, 'F2': F2, 'walk_calls': [Nb - 12, Ng - 12],
          'walk_inits': 2, 'flushes': [Rn, Rn]}
    rc = role_costs(k, ev, k['wait_cap'])
    F12 = k['f12_mul']
    cap = 2 * -(-Nb // Hn)
    ent = 1 << 32
    sort_home = 2 * 4 * ((ent // 4) * k['zero_per_4'] + 1 + ent * k['prefix_entry'] + 1 + 2)
    coins = CHARGE["delta_attempts"] * (4 * CHARGE["coin_word"] + CHARGE["delta_compare"]) + CHARGE["salt_words"] * CHARGE["coin_word"]
    tables = 2 * (30 * (((1 << KB) - 1) * k['table_build_entry'] + 1 + CHARGE["table_init"]) + CHARGE["fold_table"])
    zeroing = (2 + Rn) * ((1 << 36) // 4 * k['zero_per_4'] + 1)
    setup = coins + 2 * (CHARGE["setup_products"] * F12 + CHARGE["setup_extra"]) + tables + zeroing
    rows = {"walkers": rc["walker_main"], "relay handlers": rc["relay_delivery_1"],
            "home handlers": rc["home_delivery"],
            "radix sort": N * 4 * (k['radix_count'] + k['radix_scatter']) + Hn * sort_home,
            "merge": N * k['merge_step'] + Hn * k['merge_exit'] * (CHARGE["recheck_cap"] + 2),
            "rechecks": Hn * CHARGE["recheck_cap"] * (CHARGE["recheck_products"] * F12 + CHARGE["recheck_extra"] + CHARGE["pair_expand"]),
            "reports": rc["walker_delivery"] + rc["relay_delivery_2"] + rc["relay_main"] + rc["home_main_transport"],
            "home drivers": Hn * CHARGE["home_driver"], "setup": setup,
            "node fixed": (2 + Rn + Hn) * CHARGE["node_fixed"], "output": CHARGE["output"]}
    total = sum(rows.values())
    per = k['walk_step'] + k['relay_rec'] + k['home_rec'] + 4 * (k['radix_count'] + k['radix_scatter']) + k['merge_step']
    score = up(Decimal(total) / C_EQ, 6)
    sb = stride(K)
    walker_b = 30 * (1 << KB) * 32 + 32 * 32 + Rn * sb + Rn * 8 + (1 << 20)
    relay_b = Hr * sb + Hr * 8 + (1 << 20)
    home_b = 4 * cap * 40 + (1 << 35) + (1 << 20)
    peak = 2 * walker_b + Rn * relay_b + Hn * home_b
    wk = [(WAKE, 1 << 11), (WIN, 1 << 20)]
    lay = {"walker": wk + [(RING, 384), (CONST, 64), (ST, 4096), (CURB, Rn * 8), ((1 << 41) - (1 << 16), 1 << 16),
                      ((1 << 42) - (1 << 16), 1 << 16), (SBLKA, sb), (BLK, Rn * sb), (FOLD, 1024), (TB, (1 << KB) * 960)],
           "relay": wk + [(ST, 4096), (CURB, Hr * 8), ((1 << 41) - (1 << 16), 1 << 16), ((1 << 42) - (1 << 16), 1 << 16),
                     ((1 << 43) - (1 << 16), 1 << 16), (SBLKA, sb), (BLK, Hr * sb)],
           "home": wk + [(ST, 4096), (CNTB, 1 << 35), ((1 << 41) - (1 << 16), 1 << 16), ((1 << 42) - (1 << 16), 1 << 16),
                    (SBLKA, sb), (BABY, cap * 40), (GIANT, cap * 40), (SORTB, cap * 40), (SORTB + (1 << 59), cap * 40)]}
    # success: exact coverage minus three loss terms, each bounded above and stated as 2^-e with e rounded down
    getcontext().prec = 120
    acc = Decimal(R - 1) / Decimal(2) ** 254
    q_delta = (1 - acc) ** CHARGE["delta_attempts"]
    e_delta = -lg(q_delta)
    dup_pairs = Fraction(Nb * (Nb - 1) + Ng * (Ng - 1), 2 << 250) + Fraction(11 * Ng * (Ng - 1), 2 * R)
    dup_e = Decimal(DUP_C * Nb * Nb) / Decimal(R)                   # bounds dup_pairs (checked below)
    p_dup = dup_e / Decimal(2) ** 50
    mu = Decimal(Nb) / Decimal(Hn)                                   # the larger list mean (weakest bound)
    ratio = (Decimal(cap) - Decimal(2) ** 50) / (Decimal(Nb) / Decimal(Hn))
    chern = (ratio * ratio.ln() - (ratio - 1)) * mu / Decimal(2).ln()  # -log2 of the Chernoff bound per list
    e_over = min(-lg(p_dup), chern - 1 - 72) - 1                     # overflow: repeats >= 2^50 or a list past the cap
    foreign = Decimal(cap) ** 2 * (Decimal(2) ** -128 + Decimal(2) ** -250)
    e_exh = -lg(foreign / Decimal(CHARGE["recheck_cap"]))
    fl1 = lambda v: math.floor(v * 10) / 10
    cov = Fraction(m * (12 * m + 1), R)
    cov_lb = cov - cov * cov / 2
    succ = Decimal(cov_lb.numerator) / Decimal(cov_lb.denominator) - sum(Decimal(2) ** -Decimal(str(fl1(e))) for e in (e_delta, e_over, e_exh))
    rows_lg = {a: up(Decimal(b) / C_EQ, 4) for a, b in rows.items()}
    return {"m": str(m), "h": str(h), "records": str(N), "coins": coins, "hop1_full_blocks_max": str(F1), "hop2_full_blocks_max": str(F2),
            "per_record_instructions": per, "home_cap_per_list": str(cap),
            "rows": {a: str(b) for a, b in rows.items()}, "rows_log2_up": {a: str(v) for a, v in rows_lg.items()},
            "total_instructions": str(total), "log2_W_floor_16dp": str((lg(Decimal(total) / C_EQ) * 10 ** 16).to_integral_value(rounding=ROUND_FLOOR) / Decimal(10 ** 16)),
            "log2_W_ceil_16dp": str((lg(Decimal(total) / C_EQ) * 10 ** 16).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** 16)),
            "time_log2": str(score), "score_rational": [str(int(score * 10 ** 6)), "1000000"],
            "setup_log2_up": str(up(Decimal(setup) / C_EQ, 4)),
            "walker_bytes": str(walker_b), "relay_bytes": str(relay_b), "home_bytes": str(home_b),
            "every_node_below_2^63": max(walker_b, relay_b, home_b) < 1 << 63,
            "layout_disjoint_below_2^63": all(regions_ok(v) for v in lay.values()),
            "peak_bytes_log2_up": str(up(peak, 6)),
            "coverage_lower_bound_floor": str(math.floor(Fraction(cov_lb) * 10 ** 10) / Decimal(10 ** 10)),
            "loss_exponents_floor": {"delta_sampling": fl1(float(e_delta)), "overflow": fl1(float(e_over)),
                                     "recheck_budget": fl1(float(e_exh))},
            "foreign_pairs_log2_floor": fl1(float(-lg(foreign))),
            "success_lower_bound_floor": str(math.floor(succ * 10 ** 10) / Decimal(10 ** 10)),
            "repeats_bound_ok": dup_pairs <= Fraction(DUP_C * Nb * Nb, R) < 1,
            "success_above_0.39": succ > Decimal("0.39"), "coverage_ge_0.3901": cov_lb >= TARGET}


# ---------------------------------------------------------------- toy: the whole algorithm in a BN-family field

TOY_FIELDS = {"small": -15, "large": -41}
TOY_PAIRS = {"small": 100000, "large": 30000}
TOY_SIZES = (20, 40, 60)


def ceil_dp(v, d):
    return math.ceil(v * 10 ** d) / 10 ** d


def toy_field(name):
    q, rp, a0 = bn_toy(TOY_FIELDS[name])
    m, h = stride_params(rp)
    lam = q % rp
    Bset = {(j * pow(lam, kk, rp)) % rp for j in range(m + 1) for kk in range(12)}
    mu = Fraction(h * len(Bset), rp)
    return q, rp, a0, m, h, lam, Bset, mu - mu * mu / 2


def toy_bounds():
    """Seed-independent figures of the toy: the exact Bonferroni lower bound for each field (with the actual
    size of B), and the binomial spread bounds 0.5/sqrt(n) of the run fraction and of the coverage estimate."""
    out = {}
    for f in TOY_FIELDS:
        lb = toy_field(f)[7]
        out[f] = {"bonferroni_lower_bound_floor": math.floor(lb * 10 ** 4) / 10 ** 4, "pairs": TOY_PAIRS[f],
                  "pairs_spread_at_most": ceil_dp(0.5 / math.sqrt(TOY_PAIRS[f]), 4)}
    out["run_spread_at_most"] = {n: ceil_dp(0.5 / math.sqrt(n), 3) for n in TOY_SIZES}
    return out


def toy(seed, targets, field="small", relays=4, homes_per_relay=4, cap=None):
    q, rp, a0, m, h, lam, Bset, lb = toy_field(field)
    TW = Tower(q, a0)
    rng = random.Random(10 ** 6 + seed)
    while True:
        y = tuple(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(3)) for _ in range(2))
        g = TW.pow(y, (q ** 12 - 1) // rp)
        if g != TW.one:
            break
    capd = cap or 2 * -(-(m + 1) // (relays * homes_per_relay))
    cg = charpoly(TW, g, q)
    baby = walk(cg, [12 * TW.pow(g, j)[0][0][0] % q for j in range(12)], m + 1, q)
    solved = wrong = inwin = overflow = falsep = agree_clean = clean = 0
    for _ in range(targets):
        x = rng.randrange(rp)
        z = TW.pow(g, x)
        delta = rng.randrange(1, rp)
        D = TW.pow(g, rp - delta)
        s1, s2 = rng.getrandbits(64), rng.getrandbits(64)
        init, cur = [], z
        for _ in range(12):
            init.append(12 * cur[0][0][0] % q); cur = TW.mul(cur, D)
        giant = walk(charpoly(TW, D, q), init, h, q)
        homes, flag = {}, False
        for lst, trs in ((0, baby), (1, giant)):
            for idx, t in enumerate(trs):
                Hh = homes.setdefault((khash4(t, s1) % relays, khash4(t, s2) % homes_per_relay), ([], []))
                if len(Hh[lst]) >= capd:
                    flag = True
                    continue
                Hh[lst].append((t, idx))
        found = None
        if not flag:
            for Bl, Gl in homes.values():
                Bl.sort(); Gl.sort()
                a = b = 0
                while a < len(Bl) and b < len(Gl) and found is None:
                    if Bl[a][0] != Gl[b][0]:
                        a, b = (a + 1, b) if Bl[a][0] < Gl[b][0] else (a, b + 1)
                        continue
                    kk = Bl[a][0]; a2, b2 = a, b
                    while a2 < len(Bl) and Bl[a2][0] == kk: a2 += 1
                    while b2 < len(Gl) and Gl[b2][0] == kk: b2 += 1
                    for _, j in Bl[a:a2]:
                        for _, i in Gl[b:b2]:
                            gam, gj = TW.mul(z, TW.pow(D, i)), TW.pow(g, j)
                            hit = [c_ for c_ in range(12) if gam == TW.frob_k(gj, c_)]
                            if not hit:
                                falsep += 1
                                continue
                            xx = (delta * i + j * pow(lam, hit[0], rp)) % rp
                            if found is None and TW.pow(g, xx) == z:
                                found = xx
                    a, b = a2, b2
        covered = any((x - delta * i) % rp in Bset for i in range(h))
        overflow += flag
        inwin += covered
        solved += found is not None
        wrong += found is not None and found != x
        if not flag:
            clean += 1
            agree_clean += (found == x) == covered
    # coverage estimate: fresh (x, delta) pairs, exact membership test, no group operations
    rr = random.Random(2 * 10 ** 6 + seed)
    npairs, cov = TOY_PAIRS[field], 0
    for _ in range(npairs):                          # covered iff some x - delta i (i < h) lies in B
        s, d = rr.randrange(rp), rr.randrange(1, rp)
        cov += not Bset.isdisjoint(map(rp.__rmod__, range(s, s - d * h, -d)))
    frac = cov / npairs
    spread = ceil_dp(0.5 / math.sqrt(npairs), 4)
    lbf = math.floor(lb * 10 ** 4) / 10 ** 4
    return {"field": field, "q": q, "r_prime": rp, "m": m, "h": h, "baby_set_size": len(Bset), "relays": relays,
            "homes_per_relay": homes_per_relay, "cap_per_list": capd,
            "algorithm": {"targets": targets, "covered_by_window": inwin, "solved": solved, "wrong_outputs": wrong,
                          "runs_without_overflow": clean,
                          "solved_equals_ground_truth_without_overflow": agree_clean == clean,
                          "overflow_runs": overflow, "pairs_failing_all_12_conjugates": falsep,
                          "generic_predicted": round(targets * (m + 1) * h / q, 2),
                          "run_fraction": round(inwin / targets, 4),
                          "run_fraction_spread_at_most": ceil_dp(0.5 / math.sqrt(targets), 3)},
            "coverage": {"pairs": npairs, "covered": cov, "fraction": round(frac, 5), "spread_at_most": spread,
                         "bonferroni_lower_bound_floor": lbf,
                         "fraction_minus_bound_in_spreads": math.floor((frac - lbf) / spread * 10) / 10,
                         "fraction_at_least_bound": frac >= lbf, "fraction_at_least_0.39": frac >= 0.39}}


# ---------------------------------------------------------------- stated figures (fail closed)

NUM = re.compile(r'(?<![A-Za-z0-9_.,])(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?![\d,]\d)')


def doc_lines(L, k, isa, alg):
    """The exact lines proof.md and machine.md must contain, build-filled from these values."""
    fmt = lambda n: f"{int(n):,}"
    le = L["loss_exponents_floor"]
    proof = ["**time_log2 = " + L["time_log2"] + "**",
             "Exact rational: %s/%s." % tuple(L["score_rational"]),
             "- *Memory:* at most 2^%s bytes, over nodes that each hold below 2^63 bytes." % L["peak_bytes_log2_up"],
             "- *Setup:* at most 2^%s units, included in the total." % L["setup_log2_up"],
             "Sup-norm lower bound: 2^%s; log2 m <= %s." % (alg["nonoverlap_sup_lower_bound_log2_floor"], alg["m_log2_ceil"]),
             "m = h = " + L["m"],
             "Records: (m + 1) + h = %s." % L["records"],
             "Full blocks at hop 1: at most %s." % L["hop1_full_blocks_max"],
             "Full blocks at hop 2: at most %s." % L["hop2_full_blocks_max"],
             "List cap: %s records per list." % L["home_cap_per_list"],
             "Per record: **%s** executed instructions." % fmt(L["per_record_instructions"]),
             "| L1 | delta sampling fails 128 times: (1 - (r - 1)/2^254)^128 | < 2^-%s |" % le["delta_sampling"],
             "| L2 | a list exceeds its cap (overflow) | < 2^-%s |" % le["overflow"],
             "| L3 | the home of a covering pair uses its 2^10 rechecks on foreign pairs first | < 2^-%s |" % le["recheck_budget"],
             "Expected foreign pairs at a home: below 2^-%s." % L["foreign_pairs_log2_floor"],
             "Expected repeats in all lists: at most %d(m + 1)^2/r < 1." % DUP_C,
             "node's wake block), so a wait runs at most %d passes: %d instructions, charged on every wait." % (CHARGE["wait_passes"], k["wait_cap"]),
             "- *home drivers:* 2^%d per home: the overflow test, the calls of `wait`, the 8 radix passes, `merge` and" % (CHARGE["home_driver"].bit_length() - 1),
             "Success >= 0.3901 - 2^-%s - 2^-%s - 2^-%s >= %s > 0.39." % (le["delta_sampling"], le["overflow"],
                                                                         le["recheck_budget"], L["success_lower_bound_floor"]),
             "| walker | 2 | %s |" % fmt(L["walker_bytes"]), "| relay | 2^36 | %s |" % fmt(L["relay_bytes"]),
             "| home | 2^72 | %s |" % fmt(L["home_bytes"]),
             "Peak total: 2^%s bytes (rounded up)." % L["peak_bytes_log2_up"],
             "%d probes rejected, %d in-range forms accepted." % (isa["rejected"], isa["accepted"]),
             "b9c9b8f4a214f4c3327f27913caa78efa7bf260e7244e7dc6ba165e2172ce708"]
    tb = toy_bounds()
    proof += ["| %d | %.3f |" % (n, tb["run_spread_at_most"][n]) for n in TOY_SIZES]
    proof += ["Coverage estimate, %s field: %s pairs, spread at most %.4f, Bonferroni bound %.4f." %
              (f, fmt(tb[f]["pairs"]), tb[f]["pairs_spread_at_most"], tb[f]["bonferroni_lower_bound_floor"]) for f in TOY_FIELDS]
    proof += ["| %s | %s | %s |" % (a, fmt(L["rows"][a]), L["rows_log2_up"][a]) for a in L["rows"]]
    proof += ["| **total** | **%s** | **%s** |" % (fmt(L["total_instructions"]), L["time_log2"])]
    units = [("table step", "table_step"), ("emit", "emit"), ("walk step (two calls, table step, emit, counter, branch)", "walk_step"),
             ("relay record", "relay_rec"), ("home record", "home_rec"), ("radix count, per record and pass", "radix_count"),
             ("radix scatter, per record and pass", "radix_scatter"), ("merge step", "merge_step"),
             ("walk_init (12 ring traces)", "walk_init"), ("callsend_w with send (full walker block)", "callsend_w"),
             ("relay message, fixed part", "relay_msg_base"), ("relay full-block extra (call, callsend_r, send)", "relay_full_extra"),
             ("relay END extra", "relay_end"), ("relay flush, per home", "relay_flush_home"),
             ("home message, fixed part", "home_msg_base"), ("home END extra", "home_end"),
             ("walker flush, per relay", "flush_relay"), ("walker flush, fixed", "flush_fixed"),
             ("report send", "report"), ("report handler (last report)", "rep_handler_last"),
             ("report handler (other reports)", "rep_handler"), ("relay_up", "relay_up"), ("output", "output"),
             ("wait (cap: two empty passes)", "wait_cap"), ("setup_send", "setup_send"), ("h_setup", "h_setup"),
             ("zero, per 4 counters", "zero_per_4"), ("prefix, per counter", "prefix_entry"),
             ("table build, per entry", "table_build_entry"), ("merge exit, per call", "merge_exit"),
             ("fp_mul", "fp_mul"), ("Fp12 product", "f12_mul"), ("fp_mul calls in one Fp12 product", "f12_fp_mul_calls")]
    mach = ["| %s | %s |" % (lab, fmt(k[key])) for lab, key in units]
    mach.append("| **per record** | **%s** |" % fmt(L["per_record_instructions"]))
    mach.append("Per record: %d + %d + %d + 4 x (%d + %d) + %d." % (k["walk_step"], k["relay_rec"], k["home_rec"],
                                                                  k["radix_count"], k["radix_scatter"], k["merge_step"]))
    return proof, mach


# Structural constants that may appear in proof.md and machine.md besides computed values: field and tower
# degrees, word and digit sizes, shift and register widths, routing widths, node counts as powers of two,
# block layout, encodability ranges, coverage target and the cost-model constant. Every other numeral must be
# a computed value.
STRUCTURAL = {"0", "0.5", "3", "7", "8", "13", "15", "24", "28", "32", "33", "35", "40", "41", "42", "45", "46", "48", "50",
              "51", "52", "58", "59", "60", "64", "73", "250", "255", "256", "259", "504", "512", "1000", "4095",
              "32760"}


def constants_used(L):
    """Code constants a document may quote: charges, coin total, block bytes, parameter bytes, table row
    bytes, C."""
    fmt = lambda n: f"{int(n):,}"
    vals = list(CHARGE.values()) + [DUP_C, L["coins"], stride(PROD['K']), PARB, 960, C_EQ]
    return {fmt(v) for v in vals} | {str(v) for v in vals}


def read_docs():
    names = ("submission.json", "proof.md", "machine.md")
    if not all((HERE / n).exists() for n in names):
        return None
    man = json.loads((HERE / "submission.json").read_text())
    return man, (HERE / "proof.md").read_text(), (HERE / "machine.md").read_text()


def stated_figures(L, k, isa, prof, alg, docs):
    if docs is None:
        return {"stated_figures_match": False, "note": "package files not found next to the script"}
    man, pr_, mach = docs
    st = man["claims"][0]["statement"]
    sc = man["claims"][0]["score"]
    getcontext().prec = 120
    claimed = Decimal(sc["numerator"]) / Decimal(sc["denominator"])
    fmt = lambda n: f"{int(n):,}"
    lines_p, lines_m = doc_lines(L, k, isa, alg)
    pl, ml = {x.strip() for x in pr_.splitlines()}, {x.strip() for x in mach.splitlines()}
    miss = [x for x in lines_p if x not in pl] + [x for x in lines_m if x not in ml]
    claim_need = ["2^" + L["time_log2"], "2^" + L["peak_bytes_log2_up"], "2^" + L["setup_log2_up"],
                  "success at least 0.39", "C=20,000", "below 2^63"]
    miss += [x for x in claim_need if x not in st]
    claim_ok = {L["time_log2"], L["peak_bytes_log2_up"], L["setup_log2_up"], "0.39", "20,000", "63", "0", "1", "2"}
    extra_claim = sorted(set(NUM.findall(st)) - claim_ok)
    computed = constants_used(L)
    for x in lines_p + lines_m:
        computed |= set(NUM.findall(x))
    unknown = sorted(set(NUM.findall(pr_ + "\n" + mach)) - computed - STRUCTURAL)
    rows_ok = all(up(Decimal(int(v)) / C_EQ, 4) == Decimal(L["rows_log2_up"][a]) for a, v in L["rows"].items())
    ok_static = isa["all_correct"] and all(prof.get(x) for x in ("profile_found", "sha256_matches_protected", "p_matches",
                                                                  "r_matches", "x_matches", "g_matches", "fp12_basis_matches"))
    ok = (claimed == Decimal(L["time_log2"]) and not miss and not extra_claim and not unknown and rows_ok and ok_static
          and L["success_above_0.39"] and L["repeats_bound_ok"] and L["coverage_ge_0.3901"] and L["every_node_below_2^63"] and L["layout_disjoint_below_2^63"])
    return {"claim_score_equals_ledger": claimed == Decimal(L["time_log2"]), "missing_lines": miss,
            "claim_numerals_not_computed": extra_claim, "doc_numerals_not_accounted": unknown,
            "rows_rounded_up": rows_ok, "static_checks": ok_static, "stated_figures_match": ok}


def mutation_test(k, isa, prof, alg, docs, baseline_ok):
    """The checker must reject the package when a priced count or a stated figure changes. Meaningful only
    when the unmutated package passes."""
    res = {"baseline_passes": baseline_ok}
    for name, val in (("table_step+1", k["table_step"] + 1), ("table_step=1000", 1000)):
        km = dict(k); km["table_step"] = val; km["walk_step"] = k["walk_step"] - k["table_step"] + val
        res[name] = not stated_figures(ledger(km), km, isa, prof, alg, docs)["stated_figures_match"]
    km = dict(k); km["relay_rec"] = k["relay_rec"] + 1
    res["relay_record+1"] = not stated_figures(ledger(km), km, isa, prof, alg, docs)["stated_figures_match"]
    L = ledger(k)
    if docs is not None:
        man, pr_, mach = docs
        le = L["loss_exponents_floor"]
        swaps = [("proof: L3 bound", 1, "| < 2^-%s |" % le["recheck_budget"], "| < 2^-%s |" % L["foreign_pairs_log2_floor"]),
                 ("proof: home bytes", 1, fmt_(L["home_bytes"]), fmt_(L["walker_bytes"])),
                 ("machine: relay record", 2, "| relay record | %d |" % k["relay_rec"], "| relay record | %d |" % k["home_rec"]),
                 ("claim: memory", 0, "2^" + L["peak_bytes_log2_up"], "2^" + L["time_log2"])]
        for name, which, a_, b_ in swaps:
            d = [json.loads(json.dumps(man)), pr_, mach]
            if which == 0:
                d[0]["claims"][0]["statement"] = d[0]["claims"][0]["statement"].replace(a_, b_)
            else:
                d[which] = d[which].replace(a_, b_)
            res[name] = d[which if which else 0] != docs[which if which else 0] and \
                not stated_figures(L, k, isa, prof, alg, tuple(d))["stated_figures_match"]
    return {"detected": res, "all_detected": all(res.values())}


def fmt_(n):
    return f"{int(n):,}"


def main():
    args = json.load(sys.stdin)
    mode, seed = args.get("mode", "all"), args.get("seed", 1)
    out = {}
    k = kernels(seed)
    tr = transport(seed, k)
    alg, prof, isa = algebra(seed), profile_check(), isa_negative()
    L = ledger(k)
    docs = read_docs()
    sf = stated_figures(L, k, isa, prof, alg, docs)
    mt = mutation_test(k, isa, prof, alg, docs, sf["stated_figures_match"])
    checks = {"kernels": k["mismatches"] == 0, "transport": tr["transport_ok"],
              "algebra": all(v for v in alg.values() if isinstance(v, bool)), "isa_negative": isa["all_correct"],
              "profile": all(prof.get(x) for x in ("profile_found", "sha256_matches_protected", "p_matches", "r_matches",
                                                    "x_matches", "g_matches", "fp12_basis_matches")),
              "stated_figures": sf["stated_figures_match"], "mutation_test": mt["all_detected"]}
    if mode in ("all", "kernels"):
        out["kernels"] = k
    if mode in ("all", "transport"):
        out["transport"] = tr
    if mode in ("all", "algebra"):
        out["algebra"], out["profile"], out["isa_negative"] = alg, prof, isa
    if mode in ("all", "ledger"):
        out["charges"], out["ledger"], out["mutation_test"] = CHARGE, L, mt
    out["stated_figures"] = sf
    if mode in ("all", "toy"):
        T = toy(seed, args.get("toy_targets", 40), args.get("toy_field", "small"), cap=args.get("toy_cap"))
        out["toy"] = T
        checks["toy"] = (T["algorithm"]["wrong_outputs"] == 0 and T["algorithm"]["solved_equals_ground_truth_without_overflow"]
                         and T["coverage"]["fraction_at_least_bound"])
    out["checks"] = checks
    out["all_checks_passed"] = all(checks.values())
    out["scope"] = ("every mode runs every check except the toy and exits with status 1 if one fails; kernel counts are "
                    "executed instructions of generated programs; the simulation runs the same generated "
                    "listings (8 records per block) on 22 nodes; "
                    "the toy runs the whole algorithm at toy scale; no production attack is run")
    print(json.dumps(out, indent=1))
    if not out["all_checks_passed"]:
        print("FAIL: " + ", ".join(c for c, v in checks.items() if not v), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
