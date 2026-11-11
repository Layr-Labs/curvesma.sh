"""Supporting checks for the BN254 pairing-inversion claim: 16 random-stride walks on Frobenius traces in
GT, every equal-trace pair across the walks found in one sorted list per home (recheck against all 12
conjugates, table-lookup trace step, executed message transport), then the output P = [x]P0, Q = Q0
(outstep.py). Standard library only.

Modes (input JSON on stdin):
  kernels    execute every priced program at production shapes against Python and big integers (field
             kernels mod p and mod r, table step, table build, emit, walk loop, send, setup messages, relay
             and home handlers with their END paths, flushes, reports, waits, radix passes, zero and prefix
             loops, scan) and print the counts.
  transport  a multi-node simulation: 16 walkers walk real BN254 traces with the table step, 4 relays and 16
             homes run the executed send and handler listings under two schedules; a planted cross-walk
             target is found by a home and reported over the return route to walker 0. Every executed count
             is compared with the ledger's cost functions.
  algebra    constants from the shipped protected profile, orders, trace formula, the nonzero denominators
             of every walk pair and conjugate, the exact coverage inequality for N and its failure for N - 1;
             the output step [x]P0 with counted field operations against an affine reference.
  ledger     recompute every ledger integer and the score and print the comparison of every build-filled
             line of submission.json, proof.md and machine.md (made in every mode); a mutated step cost must
             be detected. Exit status 1 on any mismatch.
  toy        the whole algorithm (16 random-stride walks, hashed routing, caps, sort, run expansion,
             12-conjugate recheck) in a BN-family toy field; solved must equal exact ground truth when no list
             overflows; the exact cross-walk coverage of the small toy field for fresh coins.
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
import outstep as OS                                 # noqa: E402
import toy as TY                                     # noqa: E402

C_EQ = 20000
LAM = P % R
PROD = dict(K=256, relays=1 << 36, homes_per=1 << 36, S1=28, A=36, B=28, C=56, F=28)
SIMR = dict(K=8, relays=4, homes_per=4, S1=62, A=63, B=62, C=63, F=62)
ST, RING, BLK, CURB, SBLKA = 1 << 33, 1 << 26, 1 << 50, 1 << 40, 1 << 45
PORTM, PORTD = ST + 2048, ST + 2064
LIST, HCUR, SCRATCH, SORTB, CNTB = 1 << 60, ST + 3072, ST + 3200, 1 << 61, 1 << 47
CNT_ENT = 1 << 43                                    # radix counters per pass (43-bit digits)
EXP_PRODUCTS = 508                                   # products in one exponentiation below 2^254 (2 x 254)
SMALL_EXP = 2 * (max(A_W).bit_length() - 1)          # products in z^a for a <= 16 (binary powering)
CHARGE = {"coin_word": 16, "draws": 2 * NW, "draw_attempts": 128, "draw_compare": 64, "salt_words": 2,
          "setup_products": 4096, "setup_extra": 1 << 24, "table_init": 1 << 10, "fold_table": 1 << 10,
          "exp_products": EXP_PRODUCTS, "recheck_cap": 1 << 10,
          "recheck_products": 3 * EXP_PRODUCTS + 2 * (SMALL_EXP + 1) + 12, "recheck_scalar": 1 << 10,
          "recheck_extra": 1 << 12, "pair_expand": 64,
          "home_driver": 1 << 12, "node_fixed": 1 << 20, "output": 1 << 16, "wait_passes": 3,
          "point_out": 1 << 22, "pairing_check": 1 << 28}


def build(rt):
    pr = Program()
    gen_fp_mul(pr, 'fp_mul', 0); gen_fp_add(pr, 'fp_add', 0); gen_fp_sub(pr, 'fp_sub', 0)
    t12 = gen_f12_mul(pr)
    for ph in range(12):
        pr.label('tstep%d' % ph); gen_tabstep(pr, ph); pr.emit('ret')
    gen_tabbuild(pr)
    gen_emit(pr, rt['S1'], rt['A'], rt['B'], rt['C'])
    gen_transport(pr, rt['K'], rt['F'])
    gen_walk(pr); gen_walk_all(pr); gen_walk_all(pr, 'walk_all_t', 24)
    for ps in range(3):
        gen_radix(pr, ps)
    gen_table_ops(pr); gen_scan(pr)
    gen_fp_mul(pr, 'fr_mul', 64)                     # the same listing with the modulus r (recheck scalars)
    return pr, t12


def setv(m_, addr, x):
    for k, w in enumerate(limbs(x)):
        m_.st(addr + 8 * k, w)


def getv(m_, addr):
    return sum(m_.ld(addr + 8 * k) << (64 * k) for k in range(4))


def consts(m_):
    setv(m_, CONST, P)
    m_.st(CONST + 32, (-pow(P, -1, 1 << 64)) & M64)
    setv(m_, CONST + 64, R)
    m_.st(CONST + 96, (-pow(R, -1, 1 << 64)) & M64)
    m_.r['x28'] = CONST


def walker_state(pr_, mem, rt, w, s1, s2, z, deltas, betas, nrelays, coeffs):
    """Walker w's state; walker 0 holds every coin, the others receive them in the setup message."""
    M_ = Proc(pr_, mem, coeffs, 1 << 41)
    for k, v in enumerate([s1, s2, 0xbf58476d1ce4e5b9, 0x94d049bb133111eb, 0, 0, BLK, stride(rt['K']), CURB,
                           w << 60, 40, rt['K'], PORTM, nrelays]):
        M_.st(ST + 8 * k, v)
    for k, c_ in enumerate(flat(z)):
        setv(M_, ST + PAR + 32 * k, c_)
    for q in range(NW):
        setv(M_, ST + DLT + 32 * q, deltas[q] if deltas else 0); setv(M_, ST + BET + 32 * q, betas[q] if betas else 0)
    M_.st(ST + REP + NEXP, nrelays); M_.st(ST + REP + SBLK, SBLKA)
    consts(M_)
    M_.r.update({'x29': ST, 'x27': ST, 'x19': RING, 'x26': TB, 'x24': TB + 480, 'x23': 960, 'x25': FOLD})
    return M_


def relay_state(m_, rt, nh):
    for k, v in enumerate([BLK, stride(rt['K']), CURB, 40, rt['K'], 0, 0, 0, PORTD]):
        m_[ST + 8 * k] = v
    m_[ST + REP + NH], m_[ST + REP + NEXP], m_[ST + REP + RBLK], m_[ST + REP + MPORT] = nh, nh, SBLKA, PORTM


def home_state(m_, cap):
    for k, v in enumerate([LIST, 0, HCUR, 40, cap, 0, SCRATCH, 0, PORTD]):
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


def run_scan(sc, base, cnt):
    """Drives the scan listing over one sorted list of cnt records; each equal-key run is expanded by the pair
    routine, which resumes the scan at the run's last record. Returns (executed count, advancing steps, calls,
    pairs)."""
    total = steps = calls = 0
    pairs = []
    last = base + 40 * (cnt - 1)
    key = lambda a_: (sc.ld(a_), sc.ld(a_ + 8))
    sc.r['x1'], sc.r['x3'] = base, last
    while cnt >= 2 and sc.r['x1'] < last:
        x1_ = sc.r['x1']
        total += sc.call('scan'); calls += 1
        steps += (sc.r['x1'] - x1_) // 40
        if sc.r['x9']:
            a0 = e_ = sc.r['x1']
            while e_ < last and key(e_ + 40) == key(a0): e_ += 40
            if e_ == a0:                             # an equal exit on unequal keys: report it as a pair
                pairs.append((a0, a0 + 40)); e_ = a0 + 40
            pairs += [(a_, b_) for a_ in range(a0, e_ + 40, 40) for b_ in range(a_ + 40, e_ + 40, 40)]
            sc.r['x1'] = e_
    return total, steps, calls, pairs


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
    for name, mod in (('fp_mul', P), ('fr_mul', R)):
        edge = [(2 * mod - 1, 2 * mod - 1), (2 * mod - 1, 1), (0, 2 * mod - 1), (mod, mod), (mod - 1, 2 * mod - 1)]
        for it in range(120):
            x, y = edge[it] if it < len(edge) else (rng.randrange(2 * mod), rng.randrange(2 * mod))
            setv(mc, A_, x); setv(mc, B_, y); mc.r['x0'], mc.r['x1'], mc.r['x2'] = O_, A_, B_
            out[name] = mc.call(name); bad += getv(mc, O_) != x * y * pow(RM, -1, mod) % mod
    YA, YB, YO = 1 << 24, (1 << 24) + 4096, (1 << 24) + 8192
    for _ in range(3):
        xa, xb = [rng.randrange(P) for _ in range(12)], [rng.randrange(P) for _ in range(12)]
        a, b = TW.to_tower(xa), TW.to_tower(xb)
        for k, v in enumerate(mem_order(a)): setv(mc, YA + 32 * k, mont(v))
        for k, v in enumerate(mem_order(b)): setv(mc, YB + 32 * k, mont(v))
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = YO, YA, YB
        out['f12_mul'] = mc.call('f12_mul')
        bad += TW.from_tower(from_mem([unmont(getv(mc, YO + 32 * k)) for k in range(12)])) != TW.schoolbook(xa, xb)
    # table step on two production walks: start y_w = z^(w+1) g^beta_w, step D_w = g^-delta_w (w = 0 and 15)
    deltas = [rng.randrange(1, R) for _ in range(NW)]
    betas = [rng.randrange(R) for _ in range(NW)]
    z = TW.pow(g, rng.randrange(1, R))
    step_bad = 0
    for w in (0, NW - 1):
        yy = TW.pow(g, (-deltas[w]) % R)
        y0 = TW.mul(TW.pow(z, A_W[w]), TW.pow(g, betas[w]))
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
            for k, w_ in enumerate(limbs(cst)): bm.r['x%d' % (9 + k)] = w_
            for k in range(4): bm.r['x%d' % (5 + k)] = 0
            bm.r['x20'], bm.r['x21'] = TB + 960 + 32 * t, TB + 960 * 41 + 32 * t
            n_tb = bm.call('tabbuild')
            out['table_build_entry'] = (n_tb - 1) // 40
            bad += n_tb != 40 * out['table_build_entry'] + 1
            for v in range(1, 41):
                bad += getv(bm, TB + 960 * v + 32 * t) != getv(tm, TB + 960 * v + 32 * t)
    bad += step_bad
    # emit at production routing (K = 256, no block fills), walker 5
    s1, s2 = rng.getrandbits(64), rng.getrandbits(64)
    W = walker_state(pr, {}, rt, 5, s1, s2, z, deltas, betas, 3, None)
    for _ in range(16):
        t = rng.randrange(P)
        j = W.ld(ST + 32) | (W.ld(ST + 40) << 64)
        for k in range(4): W.r['x%d' % k] = (t >> (64 * k)) & M64
        w1, w2 = khash4(t, s1), khash4(t, s2)
        d1, d2 = w1 >> 28, ((w1 & ((1 << 28) - 1)) << 8) | (w2 >> 56)
        c0 = W.ld(CURB + 8 * d1)
        out['emit'] = W.call('emit')
        slot = BLK + sb * d1 + HDR + 40 * c0
        bad += [W.ld(slot + 8 * k) for k in range(5)] != [t & M64, (t >> 64) & M64, d2 | (5 << 60), j & M64, j >> 64]
        bad += W.ld(CURB + 8 * d1) != c0 + 1
    # walk loop at production shapes on walk 0: 1, 12 and 24 steps (no block fills)
    D0 = TW.pow(g, (-deltas[0]) % R)
    cg = charpoly(TW, D0, P)
    trg, cur = [], TW.mul(z, TW.pow(g, betas[0]))
    for _ in range(12 + 53):
        trg.append(12 * cur[0][0][0] % P); cur = TW.mul(cur, D0)
    wl = {}
    for nsteps in (1, 12, 24):
        Wk = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 3, cg)
        for j in range(12): setv(Wk, RING + 32 * j, trg[j])
        Wk.r['x22'] = nsteps
        wl[nsteps] = Wk.call('walk')
        bad += Wk.ld(ST + 32) != nsteps
        bad += any(getv(Wk, RING + 32 * ((n_ + 12) % 12)) != trg[n_ + 12] for n_ in range(nsteps - min(nsteps, 12), nsteps))
    out['walk_step'] = wl[1] - 3
    bad += wl[12] != 3 + 12 * out['walk_step'] or wl[24] != 3 + 24 * out['walk_step'] + 1
    wa = []                                          # walk_all: q = 0 (production), q = 1, 2 (24-step chunk)
    for lab, q_ in (('walk_all', 0), ('walk_all_t', 1), ('walk_all_t', 2)):
        Wk = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 3, cg)
        for j in range(12): setv(Wk, RING + 32 * j, trg[j])
        n_ = 24 * q_ + 5
        Wk.r['x21'], Wk.r['x22'] = q_, 5
        wa.append(Wk.call(lab) - q_ * walk_cost(out, 24) - walk_cost(out, 5))
        bad += Wk.ld(ST + 32) != n_ or any(getv(Wk, RING + 32 * (i % 12)) != trg[i] for i in range(n_, n_ + 12))
    out['walk_all_fixed'], out['walk_all_call'] = wa[0], wa[1] - wa[0]
    bad += wa[2] - wa[1] != wa[1] - wa[0]
    Wi = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 3, cg)
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
    Wf = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 3, None)
    Wf.on_stlr = w_hook
    d2s = set()
    while len(d2s) < rt['K']:                        # distinct homes, so every record's slot is known
        d2s.add(rng.getrandbits(36))
    d2s = sorted(d2s, key=lambda _: rng.random())
    hot = d2s[0]
    recs = []
    for q in range(rt['K']):                         # one full walker block for relay d1 = 5
        t = rng.randrange(P)
        recs.append([t & M64, (t >> 64) & M64, d2s[q] | (rng.getrandbits(4) << 60), q, 0])
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
    bad += hmem.get(HCUR, 0) != rt['K']
    out['home_msg_full'] = hop2[0][3] if hop2 else 0
    smem = {}
    for k, v in enumerate([1, 2, 9, 3, 4]): smem[(1 << 30) + HDR + 8 * k] = v
    par = flat(z) + deltas + betas                   # the parameter area: z, every delta_w, every beta_w
    for k, c_ in enumerate(par):
        for w_ in range(4): smem[(1 << 30) + HDR + 40 * rt['K'] + 32 * k + 8 * w_] = (c_ >> (64 * w_)) & M64

    def msg(nrec, end):
        smem[1 << 30], smem[(1 << 30) + 16] = nrec, end
        return deliver_to(Rd, 'h_relay', smem, 1 << 30, sb)
    hop2.clear()
    out['relay_msg_base'] = msg(0, 0)
    out['relay_rec'] = msg(1, 0) - out['relay_msg_base']
    out['relay_full_extra'] = hop1[0][2] - out['relay_msg_base'] - rt['K'] * out['relay_rec']
    rmem[ST + REP + ENDC] = NW - 2                   # the next END is the 15th, the one after it the last
    out['relay_end'] = msg(0, 1) - out['relay_msg_base']
    bad += len(hop2) != 0
    e2 = msg(0, 1) - out['relay_msg_base'] - out['relay_end']
    out['relay_flush_home'] = (e2 - 2) // 3
    bad += e2 != 2 + 3 * out['relay_flush_home'] or len(hop2) != 3 or any(h[2] != 1 for h in hop2)
    bad += [getv(Hd, ST + PAR + 32 * k) for k in range(len(par))] != par
    bad += hmem.get(ST + REP + READY) != 1
    # home handler: unit costs; planted overflow at cap 12, none at cap 64
    hm2 = {}
    H2 = dproc(pr, hm2, 1 << 42)
    hrec = [[rng.getrandbits(64), rng.getrandbits(64), rng.getrandbits(36) | (rng.getrandbits(4) << 60),
             rng.getrandbits(64), rng.getrandbits(62)] for _ in range(40)]
    hs = {}
    for i, r_ in enumerate(hrec):
        for k in range(5): hs[(1 << 30) + HDR + 40 * i + 8 * k] = r_[k]
    for capv in (64, 12):
        hm2.clear(); home_state(hm2, capv)
        hs[1 << 30], hs[(1 << 30) + 16] = 40, 0
        n40 = deliver_to(H2, 'h_home', hs, 1 << 30, sb)
        bad += hm2.get(ST + 40, 0) != int(40 > capv)
        bad += any(LIST + 40 * capv <= a_ < LIST + 40 * 40 for a_ in hm2)    # nothing stored past the cap
        if capv == 64:
            got = [tuple(hm2.get(LIST + 40 * q + 8 * k, 0) for k in range(5)) for q in range(hm2.get(HCUR, 0))]
            bad += got != [tuple(r_) for r_ in hrec]
            n40_ok = n40
    hs[1 << 30] = 0
    out['home_msg_base'] = deliver_to(H2, 'h_home', hs, 1 << 30, sb)
    out['home_rec'] = (n40_ok - out['home_msg_base']) // 40
    bad += n40_ok != out['home_msg_base'] + 40 * out['home_rec']
    bad += out['home_msg_full'] != out['home_msg_base'] + rt['K'] * out['home_rec']
    hs[(1 << 30) + 16] = 1
    for k, c_ in enumerate(par):
        for w_ in range(4): hs[(1 << 30) + HDR + 40 * rt['K'] + 32 * k + 8 * w_] = (c_ >> (64 * w_)) & M64
    out['home_end'] = deliver_to(H2, 'h_home', hs, 1 << 30, sb) - out['home_msg_base']
    bad += hm2.get(ST + REP + READY) != 1 or [getv(H2, ST + PAR + 32 * k) for k in range(len(par))] != par
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
    # walker flush (3 and 2 relays), setup messages to walkers 1..15
    fl_msgs = []
    Wfl = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 3, None)
    Wfl.on_stlr = lambda src, addr: fl_msgs.append((src.mem[addr - 8], src.mem[addr + 8]))
    n3 = Wfl.call('flush_w')
    Wfl2 = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 2, None)
    Wfl2.on_stlr = lambda src, addr: None
    n2 = Wfl2.call('flush_w')
    out['flush_relay'] = n3 - n2
    out['flush_fixed'] = n2 - 2 * out['flush_relay']
    bad += len(fl_msgs) != 3 or any(f[1] != 1 for f in fl_msgs)
    wsm = [{} for _ in range(NW)]
    WSd = [dproc(pr, wsm[w], 1 << 42) for w in range(NW)]
    hs_counts = []

    def s_hook(src, addr):
        w_, off = divmod(addr - 8 - SBLKA, SSTR)
        assert off == 0 and 1 <= w_ < NW
        hs_counts.append(deliver_to(WSd[w_], 'h_setup', src.mem, addr - 8, HDR + 8 * SETW))
    W0s = walker_state(pr, {}, rt, 0, s1, s2, z, deltas, betas, 3, None)
    W0s.on_stlr = s_hook
    out['setup_send'] = W0s.call('setup_send')
    out['h_setup'] = max(hs_counts)
    bad += len(hs_counts) != NW - 1 or min(hs_counts) != max(hs_counts)
    for w in range(1, NW):
        bad += [getv(WSd[w], ST + DLT + 32 * q) for q in range(NW)] != deltas
        bad += [getv(WSd[w], ST + BET + 32 * q) for q in range(NW)] != betas
        bad += wsm[w].get(ST) != s1 or wsm[w].get(ST + 8) != s2 or wsm[w].get(ST + REP + READY) != 1
    # radix passes: full-width keys (counter prefix in Python), then small digits with the zero and prefix loops
    IN, DST = 1 << 50, 1 << 52
    for small in (False, True):
        keys = []
        for _ in range(32):
            kv = sum(rng.randrange(64) << lo for lo, _w in RDX) if small else rng.getrandbits(128)
            keys.append(kv)
        rs = Proc(pr)
        for i, kv in enumerate(keys):
            for k, v in enumerate([kv & M64, kv >> 64, i, 0, 0]): rs.st(IN + 40 * i + 8 * k, v)
        src, dst = IN, DST
        for ps, (lo, wd) in enumerate(RDX):
            if small:
                rs.r.update({'x19': CNTB, 'x24': CNTB + 8 * 64})
                out['zero_call_64'] = rs.call('zero')
            else:
                rs.mem = {a_: v for a_, v in rs.mem.items() if not (CNTB <= a_ < CNTB + 8 * CNT_ENT)}
            rs.r.update({'x19': CNTB, 'x24': src + 40 * len(keys), 'x1': src})
            out['radix_count_pass_32_%d' % ps] = rs.call('count%d' % ps)
            digit = lambda i: ((rs.ld(src + 40 * i) | rs.ld(src + 40 * i + 8) << 64) >> lo) & ((1 << wd) - 1)
            bad += any(rs.ld(CNTB + 8 * digit(i)) == 0 for i in range(len(keys)))
            if small:
                rs.r.update({'x19': CNTB, 'x24': CNTB + 8 * 64, 'x3': 0})
                out['prefix_call_64'] = rs.call('prefix')
            else:
                run = 0
                for d_ in sorted({digit(i) for i in range(len(keys))}):
                    cnt = rs.ld(CNTB + 8 * d_); rs.st(CNTB + 8 * d_, run); run += cnt
            rs.r.update({'x19': CNTB, 'x20': dst, 'x22': 40, 'x24': src + 40 * len(keys), 'x1': src})
            out['radix_scatter_pass_32_%d' % ps] = rs.call('scatter%d' % ps)
            src, dst = dst, src
        bad += [rs.ld(src + 40 * i) | rs.ld(src + 40 * i + 8) << 64 for i in range(len(keys))] != sorted(keys)
    for ps in range(3):
        out['radix_count%d' % ps] = (out['radix_count_pass_32_%d' % ps] - 1) // 32
        out['radix_scatter%d' % ps] = (out['radix_scatter_pass_32_%d' % ps] - 1) // 32
        bad += out['radix_count_pass_32_%d' % ps] != 32 * out['radix_count%d' % ps] + 1
        bad += out['radix_scatter_pass_32_%d' % ps] != 32 * out['radix_scatter%d' % ps] + 1
    out['zero_per_4'] = (out['zero_call_64'] - 1) // 16
    out['prefix_entry'] = (out['prefix_call_64'] - 1) // 64
    bad += out['zero_call_64'] != 16 * out['zero_per_4'] + 1 or out['prefix_call_64'] != 64 * out['prefix_entry'] + 1
    # scan over one sorted list drawn from a small key pool (equal-key runs of several lengths)
    sc = Proc(pr)
    L2 = 1 << 45
    pool = [(rng.getrandbits(64), rng.getrandbits(64)) for _ in range(10)]
    a_, b_, c_ = rng.getrandbits(64), rng.getrandbits(63), rng.getrandbits(63)
    pool += [(a_, b_), (a_, b_ + 1), (c_, a_), (c_ + 1, a_)]    # adjacent keys that share one word
    kl = sorted([rng.choice(pool) for _ in range(32)] + pool[-4:] * 2, key=lambda kk: (kk[1], kk[0]))
    for i, kk in enumerate(kl): sc.st(L2 + 40 * i, kk[0]); sc.st(L2 + 40 * i + 8, kk[1])
    total, steps, calls, pairs = run_scan(sc, L2, len(kl))
    want = {(L2 + 40 * a_, L2 + 40 * b_) for a_ in range(40) for b_ in range(a_ + 1, 40) if kl[a_] == kl[b_]}
    bad += any((sc.ld(a_), sc.ld(a_ + 8)) != (sc.ld(b_), sc.ld(b_ + 8)) for a_, b_ in pairs)
    bad += set(pairs) != want or len(pairs) != len(want)
    out['scan_step'], out['scan_exit'] = 8, 8
    bad += not (steps * 8 <= total <= steps * 8 + calls * 8)
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
    Rn, Hn, Wn = ev['R'], ev['H'], ev['W']
    return {
        "walker_main": sum(walk_cost(k, c) for c in ev['walk_calls']) + ev['walk_inits'] * k['walk_init']
        + ev['F1'] * (k['callsend_w'] - 1) + sum(k['flush_fixed'] + f * k['flush_relay'] for f in ev['flushes'])
        + k['setup_send'] + Wn * wait + k['output'] + ev.get('walk_full', 0) * (walk_cost(k, WALK_CHUNK)
        + k['walk_all_call']) + ev.get('walk_alls', 0) * k['walk_all_fixed'],
        "walker_delivery": (Wn - 1) * k['h_setup'] + Rn * k['rep_handler'] + (k['rep_handler_last'] - k['rep_handler']),
        "relay_delivery_1": (ev['F1'] + Wn * Rn) * k['relay_msg_base'] + ev['N'] * k['relay_rec']
        + ev['F2'] * k['relay_full_extra'] + Wn * Rn * k['relay_end'] + Rn * (2 + ev['Hr'] * k['relay_flush_home']),
        "relay_delivery_2": Hn * k['rep_handler'] + Rn * (k['rep_handler_last'] - k['rep_handler']),
        "relay_main": Rn * (wait + k['relay_up']),
        "home_delivery": (ev['F2'] + Hn) * k['home_msg_base'] + ev['N'] * k['home_rec'] + Hn * k['home_end'],
        "home_main_transport": Hn * (wait + k['report'])}


def recheck(TW, g, z, lam, r, par, ra, rb):
    """The home's recheck of one equal-key pair of records (walk, index): rebuild both elements from the
    parameters, find the first k < 12 with u = sigma^k(y), solve for x and keep it only if g^x = z. Returns
    (x or None, conjugate found)."""
    deltas, betas = par
    (wa, na), (wb, nb) = ra, rb
    ea, eb = (betas[wa] - deltas[wa] * na) % r, (betas[wb] - deltas[wb] * nb) % r
    u = TW.mul(TW.pow(z, A_W[wa]), TW.pow(g, ea))
    y = TW.mul(TW.pow(z, A_W[wb]), TW.pow(g, eb))
    for kk in range(12):
        if u == TW.frob_k(y, kk):
            lk = pow(lam, kk, r)
            xx = (lk * eb - ea) * pow((A_W[wa] - lk * A_W[wb]) % r, -1, r) % r
            return (xx if TW.pow(g, xx) == z else None), True
    return None, False


# ---------------------------------------------------------------- multi-node simulation

def transport(seed, k):
    """16 walkers -> 4 relays -> 16 homes with the same generated listings (routing shifts select 2 + 2 bits
    of the same hash words; K = 8 records per block). Walkers run the table step on real BN254 traces."""
    TW = Tower(P, 9)
    g = TW.to_tower(G_PINNED)
    rng = random.Random(1000003 + seed)
    rt = SIMR
    pr, _ = build(rt)
    K, NR, NHp = rt['K'], rt['relays'], rt['homes_per']
    sb = stride(K)
    lens = [20 + rng.randrange(12) for _ in range(NW)]          # records of walk w: n = 0..lens[w]-1
    deltas, betas = [rng.randrange(1, R) for _ in range(NW)], [rng.randrange(R) for _ in range(NW)]
    s1, s2 = rng.getrandbits(64), rng.getrandbits(64)
    w0, v0 = sorted(rng.sample(range(NW), 2))
    n0, n1, k0 = rng.randrange(lens[w0]), rng.randrange(lens[v0]), rng.randrange(12)
    lk = pow(LAM, k0, R)                              # planted cross-walk target: e(w0, n0) = lambda^k0 e(v0, n1)
    xt = (lk * (betas[v0] - deltas[v0] * n1) - (betas[w0] - deltas[w0] * n0)) * pow((A_W[w0] - lk * A_W[v0]) % R, -1, R) % R
    z = TW.pow(g, xt)
    Ds = [TW.pow(g, (-d) % R) for d in deltas]
    zp = [TW.one]
    for _ in range(NW):
        zp.append(TW.mul(zp[-1], z))
    traces = []
    for w in range(NW):
        cur, tr = TW.mul(zp[A_W[w]], TW.pow(g, betas[w])), []
        for _ in range(lens[w]):
            tr.append(12 * cur[0][0][0] % P); cur = TW.mul(cur, Ds[w])
        traces.append(tr)
    cps = [charpoly(TW, Ds[w], P) for w in range(NW)]
    expected = {}
    for w, trs in enumerate(traces):
        for idx, t in enumerate(trs):
            w1, w2 = khash4(t, s1), khash4(t, s2)
            expected.setdefault((w1 >> 62, ((w1 & 1) << 1) | (w2 >> 63)), []).append((t & ((1 << 128) - 1), w, idx))
    runs = []
    for sched in (0, 1):
        log = {"hop1": 0, "hop1_end": 0, "hop2": 0, "hop2_end": 0, "setup": 0, "rep2": 0, "rep1": 0}
        wm = [WakeMem() for _ in range(NW)]
        WM = [walker_state(pr, wm[0], rt, 0, s1, s2, z, deltas, betas, NR, cps[0])]
        WM += [walker_state(pr, wm[w], rt, w, 0, 0, z, None, None, NR, None) for w in range(1, NW)]
        WD = [dproc(pr, wm[w], 1 << 42) for w in range(NW)]
        rmem = [WakeMem() for _ in range(NR)]
        for m_ in rmem: relay_state(m_, rt, NHp)
        RD1 = [dproc(pr, rmem[i], 1 << 42) for i in range(NR)]
        RD2 = [dproc(pr, rmem[i], 1 << 43) for i in range(NR)]
        RM = [Proc(pr, rmem[i], None, 1 << 41) for i in range(NR)]
        hmem = [[WakeMem() for _ in range(NHp)] for _ in range(NR)]
        cap = 2 * -(-sum(lens) // (NR * NHp)) + 64
        for row in hmem:
            for m_ in row: home_state(m_, cap)
        HD = [[dproc(pr, hmem[i][j], 1 << 42) for j in range(NHp)] for i in range(NR)]
        HM = [[Proc(pr, hmem[i][j], None, 1 << 41) for j in range(NHp)] for i in range(NR)]
        for p_ in RM + [h_ for row in HM for h_ in row]:
            p_.r['x27'] = ST

        def w_hook(src, addr):
            base = addr - 8
            if SBLKA <= base < SBLKA + NW * SSTR:
                ws, off = divmod(base - SBLKA, SSTR)
                assert off == 0 and 1 <= ws < NW
                log["setup"] += 1
                return deliver_to(WD[ws], 'h_setup', src.mem, base, HDR + 8 * SETW)
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
        for w in range(NW):
            WM[w].on_stlr = w_hook
        for i in range(NR):
            RD1[i].on_stlr = (lambda i_: (lambda s_, a_: r_hook(s_, a_, i_)))(i)
            RM[i].on_stlr = m_hook
            for j in range(NHp):
                HM[i][j].on_stlr = (lambda i_: (lambda s_, a_: h_hook(s_, a_, i_)))(i)
        # setup: walker 0 sends every coin to walkers 1..15; each waits, then derives its walk from its own memory
        WM[0].call('setup_send')
        coins_ok, memo = True, {}
        gp = lambda e_: memo[e_] if e_ in memo else memo.setdefault(e_, TW.pow(g, e_))   # repeated powers of g
        for w in range(NW):
            if w:
                WM[w].r['x1'] = ST + REP + READY
                WM[w].call('wait')
            dw, bw = getv(WM[w], ST + DLT + 32 * w), getv(WM[w], ST + BET + 32 * w)
            coins_ok &= dw == deltas[w] and bw == betas[w] and WM[w].ld(ST) == s1 and WM[w].ld(ST + 8) == s2
            zw = unflat([getv(WM[w], ST + PAR + 32 * q) for q in range(12)])
            Dw = gp((-dw) % R)
            WM[w].c = charpoly(TW, Dw, P)
            cur = TW.mul(TW.pow(zw, A_W[w]), gp(bw))
            for j in range(12):
                setv(WM[w], RING + 32 * j, 12 * cur[0][0][0] % P); cur = TW.mul(cur, Dw)
        left = [n_ - 12 for n_ in lens]
        calls = []

        def walk_chunk(w, c):
            WM[w].r['x22'] = c; WM[w].call('walk'); left[w] -= c; calls.append(c)
        order = list(range(NW)) if sched == 0 else list(range(NW))[::-1]
        for w in order:
            WM[w].call('walk_init')
        if sched == 0:                                  # alternate in chunks of 12 steps, then every flush
            while any(left):
                for w in range(NW):
                    if left[w]:
                        walk_chunk(w, min(12, left[w]))
            for w in range(NW):
                WM[w].call('flush_w')
        else:                                           # walker 15 walks and ends first, then 14, ...
            for w in order:
                walk_chunk(w, left[w]); WM[w].call('flush_w')
        # homes: wait for END, sort (Python stand-in for the radix listings run in mode kernels), scan, recheck
        got, pairs_seen, flags, params_ok = {}, 0, 0, True
        scan_exec = scan_bound = 0
        homes = [(i, j) for i in range(NR) for j in range(NHp)]
        for (i, j) in (homes if sched == 0 else homes[::-1]):
            hm, H_ = hmem[i][j], HM[i][j]
            H_.r['x1'] = ST + REP + READY
            H_.call('wait')
            flags += hm.get(ST + 40, 0)
            zh = unflat([getv(H_, ST + PAR + 32 * q) for q in range(12)])
            ph = ([getv(H_, ST + DLT + 32 * q) for q in range(NW)], [getv(H_, ST + BET + 32 * q) for q in range(NW)])
            params_ok &= zh == z and ph == (deltas, betas)
            rr = [tuple(hm.get(LIST + 40 * q + 8 * w_, 0) for w_ in range(5)) for q in range(hm.get(HCUR, 0))]
            for r_ in rr:
                got.setdefault((i, j), []).append(((r_[1] << 64) | r_[0], r_[2] >> 60, r_[3] | (r_[4] << 64)))
            rr.sort(key=lambda r_: (r_[1], r_[0]))
            for q, r_ in enumerate(rr):
                for w_ in range(5): hm[SORTB + 40 * q + 8 * w_] = r_[w_]
            found, xf = 0, 0
            if len(rr) >= 2:
                tot, steps, ncall, pairs = run_scan(H_, SORTB, len(rr))
                scan_exec += tot; scan_bound += 8 * steps + 8 * ncall
                assert steps <= len(rr)
                for a_, b_ in pairs:
                    pairs_seen += 1
                    rec = lambda c_: (hm[c_ + 16] >> 60, hm[c_ + 24] | (hm[c_ + 32] << 64))
                    xx, _conj = recheck(TW, g, zh, LAM, R, ph, rec(a_), rec(b_))
                    if xx is not None and not found:
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
        pt_ok = OS.check([xo], k['fp_mul'])["output_step_ok"]
        ok_homes = all(sorted(got.get(h_, [])) == sorted(v) for h_, v in expected.items()) and set(got) <= set(expected)
        ev = {'W': NW, 'R': NR, 'H': NR * NHp, 'Hr': NHp, 'N': sum(lens), 'F1': log['hop1'], 'F2': log['hop2'],
              'walk_calls': calls, 'walk_inits': NW, 'flushes': [NR] * NW}
        pred = role_costs(k, ev, k['wait'])
        wake1 = all(getattr(m_, 'wake', 0) == [WAKE] for m_ in wm + rmem + [h_ for row in hmem for h_ in row])
        execd = {"walker_main": sum(W_.count for W_ in WM) + W0o.count,
                 "walker_delivery": sum(W_.count for W_ in WD),
                 "relay_delivery_1": sum(p_.count for p_ in RD1), "relay_delivery_2": sum(p_.count for p_ in RD2),
                 "relay_main": sum(p_.count for p_ in RM), "home_delivery": sum(p_.count for row in HD for p_ in row),
                 "home_main_transport": sum(p_.count for row in HM for p_ in row) - scan_exec}
        runs.append({"schedule": ["chunks of 12 steps alternate", "walker 15 walks and ends first"][sched],
                     "messages": log, "coins_at_every_walker": coins_ok,
                     "every_record_once_at_its_home": ok_homes, "parameters_at_every_home": params_ok,
                     "overflow_flags": flags, "equal_key_pairs": pairs_seen, "output_found": W0o.r['x4'],
                     "output_equals_target": xo == xt and W0o.r['x4'] == 1, "output_point_ok": pt_ok,
                     "executed_by_role": execd, "cost_functions_by_role": pred,
                     "executed_equals_cost_functions": execd == pred, "each_wake_block_stored_once": wake1,
                     "scan_within_bound": scan_exec <= scan_bound,
                     "hop1_full_within_bound": log['hop1'] <= sum(n_ // K for n_ in lens),
                     "hop2_full_within_bound": log['hop2'] <= sum(lens) // K,
                     "contents": sorted((h_, tuple(sorted(v))) for h_, v in got.items())})
    same = runs[0].pop("contents") == runs[1].pop("contents")
    ok = same and all(r_[f] for r_ in runs for f in ("coins_at_every_walker", "every_record_once_at_its_home",
                                                      "parameters_at_every_home", "output_equals_target",
                                                      "output_point_ok", "executed_equals_cost_functions",
                                                      "each_wake_block_stored_once", "scan_within_bound",
                                                      "hop1_full_within_bound", "hop2_full_within_bound")) \
        and all(r_["overflow_flags"] == 0 and r_["messages"]["setup"] == NW - 1 for r_ in runs)
    return {"walkers": NW, "relays": NR, "homes": NR * NHp, "records_per_block": K, "walk_records": lens,
            "planted": {"walks": [w0, v0], "indices": [n0, n1], "k": k0}, "runs": runs,
            "both_schedules_same_home_contents": same, "transport_ok": ok}


# ---------------------------------------------------------------- algebra on the pinned constants

def algebra(seed):
    rng = random.Random(seed)
    TW = Tower(P, 9)
    g = TW.to_tower(G_PINNED)
    N = walk_params(R)
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
    lk = [pow(LAM, kk, R) for kk in range(12)]
    dens = [(A_W[w] - lk[kk] * A_W[v]) % R for w in range(NW) for v in range(NW) for kk in range(12) if w != v or kk]
    # same walk pair, k != k': the coefficients -1/(a_w - lambda^k a_v) of beta_w in x_t must differ
    coef = {(w, v, kk): pow((A_W[w] - lk[kk] * A_W[v]) % R, -1, R) for w in range(NW) for v in range(w + 1, NW) for kk in range(12)}
    cpairs = [(coef[w, v, a_], coef[w, v, b_]) for w in range(NW) for v in range(w + 1, NW) for a_ in range(12) for b_ in range(a_ + 1, 12)]
    return {"p_r_prime": is_prime(P) and is_prime(R),
            "ord_r_p_is_12": all(pow(P, kk, R) != 1 for kk in range(1, 12)) and pow(P, 12, R) == 1,
            "lambda6_is_minus1": pow(LAM, 6, R) == R - 1, "g_order_r": TW.pow(g, R) == TW.one and g != TW.one,
            "trace_is_12_Re_e0": ok_trace, "frobenius_acts_as_lambda_power": frob_is_lambda,
            "walks": NW, "z_exponents": A_W, "N": str(N), "N_least": covers(N, R) and not covers(N - 1, R),
            "N_below_r": N < R, "denominators_checked": len(dens), "denominators_nonzero": all(dens),
            "beta_coefficient_pairs_checked": len(cpairs), "beta_coefficients_distinct": all(a_ != b_ for a_, b_ in cpairs)}


PROFILE_SHA256 = "649e2c45757c3fd789e168fe1b7ea58a40b924640476551586a177ebea3c1ed0"   # the protected target profile


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
            "fp12_basis_matches": "w^12-18*w^6+82" in prof["field_Fp12"],
            "p0_q0_ate_match": [int(prof["generator_P0"][c], 16) for c in "xy"] == list(OS.P0) and int(prof["pairing"]["ate_loop_count"], 16) == OS.ATE
            and tuple(tuple(int(c, 16) for c in prof["generator_Q0"][a]) for a in "xy") == OS.Q0}


# ---------------------------------------------------------------- ledger

def lg(v):
    return Decimal(v).ln() / Decimal(2).ln()


def up(v, d):
    return (lg(v) * 10 ** d).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** d)


def ledger(k):
    getcontext().prec = 120
    N = walk_params(R)
    M = NW * N
    K, Rn, Hr = PROD['K'], PROD['relays'], PROD['homes_per']
    Hn = Rn * Hr
    F1, F2 = NW * (N // K), M // K
    q = (N - 13) // WALK_CHUNK                       # full walk calls per walker; the last call has 1 to c steps
    ev = {'W': NW, 'R': Rn, 'H': Hn, 'Hr': Hr, 'N': M, 'F1': F1, 'F2': F2, 'walk_calls': [N - 12 - q * WALK_CHUNK] * NW,
          'walk_full': NW * q, 'walk_alls': NW, 'walk_inits': NW, 'flushes': [Rn] * NW}
    rc = role_costs(k, ev, k['wait_cap'])
    F12 = k['f12_mul']
    cap = 2 * -(-M // Hn)
    sort_home = 3 * ((CNT_ENT // 4) * k['zero_per_4'] + 1 + CNT_ENT * k['prefix_entry'] + 1 + 2)
    coins = CHARGE["draws"] * CHARGE["draw_attempts"] * (4 * CHARGE["coin_word"] + CHARGE["draw_compare"]) \
        + CHARGE["salt_words"] * CHARGE["coin_word"]
    tables = NW * (30 * (((1 << KB) - 1) * k['table_build_entry'] + 1 + CHARGE["table_init"]) + CHARGE["fold_table"])
    zeroing = (NW + Rn) * ((1 << 36) // 4 * k['zero_per_4'] + 1)
    setup = coins + NW * (CHARGE["setup_products"] * F12 + CHARGE["setup_extra"]) + tables + zeroing
    radix_rec = sum(k['radix_count%d' % q] + k['radix_scatter%d' % q] for q in range(3))
    pair_cost = CHARGE["recheck_products"] * F12 + CHARGE["recheck_scalar"] * k['fr_mul'] + CHARGE["recheck_extra"] + CHARGE["pair_expand"]
    rows = {"walkers": rc["walker_main"], "relay handlers": rc["relay_delivery_1"],
            "home handlers": rc["home_delivery"],
            "radix sort": M * radix_rec + Hn * sort_home,
            "scan": M * k['scan_step'] + Hn * k['scan_exit'] * (CHARGE["recheck_cap"] + 2),
            "rechecks": Hn * CHARGE["recheck_cap"] * pair_cost,
            "reports": rc["walker_delivery"] + rc["relay_delivery_2"] + rc["relay_main"] + rc["home_main_transport"],
            "home drivers": Hn * CHARGE["home_driver"], "setup": setup,
            "node fixed": (NW + Rn + Hn) * CHARGE["node_fixed"], "output": CHARGE["output"],
            "point output": CHARGE["point_out"], "pairing check": CHARGE["pairing_check"]}
    ob = OS.bounds(k['fp_mul'], F12)
    total = sum(rows.values())
    per = k['walk_step'] + k['relay_rec'] + k['home_rec'] + radix_rec + k['scan_step']
    score = up(Decimal(total) / C_EQ, 6)
    sb = stride(K)
    walker_b = 30 * (1 << KB) * 32 + 32 * 32 + Rn * sb + Rn * 8 + NW * SSTR + (1 << 20)
    relay_b = Hr * sb + Hr * 8 + (1 << 20)
    home_b = 2 * cap * 40 + 8 * CNT_ENT + (1 << 20)
    peak = NW * walker_b + Rn * relay_b + Hn * home_b
    wk = [(WAKE, 1 << 11), (WIN, 1 << 20)]
    lay = {"walker": wk + [(RING, 384), (CONST, 128), (ST, 4096), (CURB, Rn * 8), ((1 << 41) - (1 << 16), 1 << 16),
                      ((1 << 42) - (1 << 16), 1 << 16), (SBLKA, NW * SSTR), (BLK, Rn * sb), (FOLD, 1024), (TB, (1 << KB) * 960)],
           "relay": wk + [(ST, 4096), (CURB, Hr * 8), ((1 << 41) - (1 << 16), 1 << 16), ((1 << 42) - (1 << 16), 1 << 16),
                     ((1 << 43) - (1 << 16), 1 << 16), (SBLKA, sb), (BLK, Hr * sb)],
           "home": wk + [(ST, 4096), (CNTB, 8 * CNT_ENT), ((1 << 41) - (1 << 16), 1 << 16), ((1 << 42) - (1 << 16), 1 << 16),
                    (SBLKA, sb), (LIST, cap * 40), (SORTB, cap * 40)]}
    # success: exact coverage minus three loss terms, each bounded above and stated as 2^-e with e rounded down
    acc = Decimal(R - 1) / Decimal(2) ** 254                       # least acceptance probability of one attempt
    q_draw = CHARGE["draws"] * (1 - acc) ** CHARGE["draw_attempts"]
    e_draw = -lg(q_draw)
    dup_pairs = Fraction(M * (M - 1), 2) * (Fraction(1, 1 << 250) + Fraction(12, R))
    dup_e = Decimal(DUP_C * M * M) / Decimal(R)                     # bounds dup_pairs (checked below)
    p_dup = dup_e / Decimal(2) ** 50
    mu = Decimal(M) / Decimal(Hn)                                    # mean list size
    ratio = (Decimal(cap) - Decimal(2) ** 50) / mu
    chern = (ratio * ratio.ln() - (ratio - 1)) * mu / Decimal(2).ln()  # -log2 of the Chernoff bound per list
    e_over = min(-lg(p_dup), chern - 72) - 1                         # overflow: repeats >= 2^50 or a list past the cap
    foreign = Decimal(cap) * Decimal(cap - 1) / 2 * (Decimal(2) ** -128 + Decimal(2) ** -250)
    e_exh = -lg(foreign / Decimal(CHARGE["recheck_cap"]))
    fl1 = lambda v: math.floor(v * 10) / 10
    cov = Fraction(cross_triples(N), R)
    cov_lb = cov - cov * cov / 2
    succ = Decimal(cov_lb.numerator) / Decimal(cov_lb.denominator) - sum(Decimal(2) ** -Decimal(str(fl1(e))) for e in (e_draw, e_over, e_exh))
    rows_lg = {a: up(Decimal(b) / C_EQ, 4) for a, b in rows.items()}
    return {"N": str(N), "records": str(M), "coins": coins, "hop1_full_blocks_max": str(F1), "hop2_full_blocks_max": str(F2),
            "per_record_instructions": per, "radix_per_record": radix_rec, "home_cap": str(cap),
            "pair_recheck_instructions": pair_cost,
            "rows": {a: str(b) for a, b in rows.items()}, "rows_log2_up": {a: str(v) for a, v in rows_lg.items()},
            "total_instructions": str(total), "log2_W_floor_16dp": str((lg(Decimal(total) / C_EQ) * 10 ** 16).to_integral_value(rounding=ROUND_FLOOR) / Decimal(10 ** 16)),
            "log2_W_ceil_16dp": str((lg(Decimal(total) / C_EQ) * 10 ** 16).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** 16)),
            "time_log2": str(score), "score_rational": [str(Fraction(int(score * 10 ** 6), 10 ** 6).numerator),
                                                       str(Fraction(int(score * 10 ** 6), 10 ** 6).denominator)],
            "setup_log2_up": str(up(Decimal(setup) / C_EQ, 4)),
            "walker_bytes": str(walker_b), "relay_bytes": str(relay_b), "home_bytes": str(home_b),
            "every_node_below_2^63": max(walker_b, relay_b, home_b) < 1 << 63,
            "layout_disjoint_below_2^63": all(regions_ok(v) for v in lay.values()),
            "peak_bytes_log2_up": str(up(peak, 6)),
            "coverage_lower_bound_floor": str(math.floor(Fraction(cov_lb) * 10 ** 10) / Decimal(10 ** 10)),
            "loss_exponents_floor": {"coin_draws": fl1(float(e_draw)), "overflow": fl1(float(e_over)),
                                     "recheck_budget": fl1(float(e_exh))},
            "foreign_pairs_log2_floor": fl1(float(-lg(foreign))), "repeats_markov_log2_floor": math.floor(-lg(p_dup)),
            "success_lower_bound_floor": str(math.floor(succ * 10 ** 10) / Decimal(10 ** 10)),
            "repeats_bound_ok": dup_pairs <= Fraction(DUP_C * M * M, R) < 2, "output_bounds": ob,
            "output_allowances_ok": ob["point"] <= CHARGE["point_out"] and ob["pairing"] <= CHARGE["pairing_check"],
            "success_above_0.39": succ > Decimal("0.39"), "coverage_ge_0.3901": cov_lb >= TARGET}


# ---------------------------------------------------------------- stated figures (fail closed)

NUM = re.compile(r'(?<![A-Za-z0-9_.,])(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?![\d,]\d)')


def doc_lines(L, k, isa, alg):
    """The exact lines proof.md and machine.md must contain, build-filled from these values."""
    fmt = lambda n: f"{int(n):,}"
    le = L["loss_exponents_floor"]
    proof = ["**time_log2 = " + L["time_log2"] + "**",
             "Exact rational: %s/%s (in lowest terms)." % tuple(L["score_rational"]),
             "- *Memory:* at most 2^%s bytes, over nodes that each hold below 2^63 bytes." % L["peak_bytes_log2_up"],
             "- *Setup:* at most 2^%s units, included in the total." % L["setup_log2_up"],
             "Denominators checked: %s, all nonzero; coefficient pairs checked: %s, all distinct." % (
                 fmt(alg["denominators_checked"]), fmt(alg["beta_coefficient_pairs_checked"])),
             "N = " + L["N"],
             "Records: 16 N = %s." % L["records"],
             "Full blocks at hop 1: at most %s." % L["hop1_full_blocks_max"],
             "Full blocks at hop 2: at most %s." % L["hop2_full_blocks_max"],
             "List cap: %s records per home." % L["home_cap"],
             "Per record: **%s** executed instructions." % fmt(L["per_record_instructions"]),
             "| L1 | one of the 32 coin draws fails 128 times: 32 (1 - (r - 1)/2^254)^128 | < 2^-%s |" % le["coin_draws"],
             "| L2 | a list exceeds its cap (overflow) | < 2^-%s |" % le["overflow"],
             "| L3 | the home of a covering pair uses its 2^10 rechecks on foreign pairs first | < 2^-%s |" % le["recheck_budget"],
             "Expected foreign pairs at a home: below 2^-%s." % L["foreign_pairs_log2_floor"],
             "Expected repeats in all lists: at most %dM^2/r < 2." % DUP_C,
             "node's wake block), so a wait runs at most %d passes: %d instructions, charged on every wait." % (CHARGE["wait_passes"], k["wait_cap"]),
             "- *home drivers:* 2^%d per home: the overflow test, the calls of `wait`, the radix passes, `scan` and" % (CHARGE["home_driver"].bit_length() - 1),
             "Success >= 0.3901 - 2^-%s - 2^-%s - 2^-%s >= %s > 0.39." % (le["coin_draws"], le["overflow"],
                                                                         le["recheck_budget"], L["success_lower_bound_floor"]),
             "- *rechecks:* per home 2^10 pairs x (%s Fp12 products + 2^10 products mod r at the executed `fr_mul`" % fmt(CHARGE["recheck_products"]),
             "count + 2^12 + 64): u and y (each a power of g with at most %d products, a power z^(w + 1) with at" % EXP_PRODUCTS,
             "most %d and one product), 12 Frobenius maps (each at most one product), g^x (%d); the products mod r" % (SMALL_EXP, EXP_PRODUCTS),
             "So at least 2^50 repeats has probability below 2^-%d (Markov). Under H2 the distinct values in a list are" % L["repeats_markov_log2_floor"],
             "- *setup:* coins %s (32 draws x 128 attempts x (4 words x 16 + 64) + 2 salt words x 16), per walker" % fmt(L["coins"]),
             *fix_lines(k, EXP_PRODUCTS // 2, CHARGE["recheck_extra"].bit_length() - 1, CHARGE["pair_expand"]),
             "| walker | 16 | %s |" % fmt(L["walker_bytes"]), "| relay | 2^36 | %s |" % fmt(L["relay_bytes"]),
             "| home | 2^72 | %s |" % fmt(L["home_bytes"]),
             "Peak total: 2^%s bytes (rounded up)." % L["peak_bytes_log2_up"],
             "%d probes rejected, %d in-range forms accepted." % (isa["rejected"], isa["accepted"]),
             PROFILE_SHA256]
    op, om = OS.doc_lines(k["fp_mul"], k["f12_mul"], CHARGE["point_out"].bit_length() - 1,
                          CHARGE["pairing_check"].bit_length() - 1)
    tb = TY.toy_bounds()
    proof += op + ["| %d | %.3f |" % (n, tb["run_spread_at_most"][n]) for n in TY.TOY_SIZES]
    sm = tb["small"]
    proof += ["Exact coverage, small field: r' = %s, N = %d, %s cross triples, %d coin draws, Bonferroni bound %.4f." %
              (fmt(sm["r_prime"]), sm["walk_length"], fmt(cross_triples(sm["walk_length"])), TY.COV_DRAWS,
               sm["bonferroni_lower_bound_floor"])]
    proof += ["Walk length in the toy: N = %d (small), N = %d (large)." % (sm["walk_length"], tb["large"]["walk_length"])]
    proof += ["| %s | %s | %s |" % (a, fmt(L["rows"][a]), L["rows_log2_up"][a]) for a in L["rows"]]
    proof += ["| **total** | **%s** | **%s** |" % (fmt(L["total_instructions"]), L["time_log2"])]
    units = [("table step", "table_step"), ("emit", "emit"), ("walk step (two calls, table step, emit, counter, branch)", "walk_step"),
             ("relay record", "relay_rec"), ("home record", "home_rec"),
             ("radix count, pass 1", "radix_count0"), ("radix scatter, pass 1", "radix_scatter0"),
             ("radix count, pass 2", "radix_count1"), ("radix scatter, pass 2", "radix_scatter1"),
             ("radix count, pass 3", "radix_count2"), ("radix scatter, pass 3", "radix_scatter2"),
             ("scan step", "scan_step"),
             ("walk_init (12 ring traces)", "walk_init"), ("walk_all, per full walk call (MOVZ, BL, SUBS, B.NE)", "walk_all_call"),
             ("walk_all, fixed (entry, last call, exit)", "walk_all_fixed"), ("callsend_w with send (full walker block)", "callsend_w"),
             ("relay message, fixed part", "relay_msg_base"), ("relay full-block extra (call, callsend_r, send)", "relay_full_extra"),
             ("relay END extra", "relay_end"), ("relay flush, per home", "relay_flush_home"),
             ("home message, fixed part", "home_msg_base"), ("home END extra", "home_end"),
             ("walker flush, per relay", "flush_relay"), ("walker flush, fixed", "flush_fixed"),
             ("report send", "report"), ("report handler (last report)", "rep_handler_last"),
             ("report handler (other reports)", "rep_handler"), ("relay_up", "relay_up"), ("output", "output"),
             ("wait (cap: two empty passes)", "wait_cap"), ("setup_send (15 setup messages)", "setup_send"),
             ("h_setup", "h_setup"), ("zero, per 4 counters", "zero_per_4"), ("prefix, per counter", "prefix_entry"),
             ("table build, per entry", "table_build_entry"), ("scan exit, per call", "scan_exit"),
             ("fp_mul", "fp_mul"), ("fr_mul (the same listing modulo r)", "fr_mul"), ("Fp12 product", "f12_mul"),
             ("fp_mul calls in one Fp12 product", "f12_fp_mul_calls")]
    mach = ["| %s | %s |" % (lab, fmt(k[key])) for lab, key in units]
    mach.append("| **per record** | **%s** |" % fmt(L["per_record_instructions"]))
    mach.append("Per record: %d + %d + %d + (%d + %d) + (%d + %d) + (%d + %d) + %d." % (
        k["walk_step"], k["relay_rec"], k["home_rec"], k["radix_count0"], k["radix_scatter0"], k["radix_count1"],
        k["radix_scatter1"], k["radix_count2"], k["radix_scatter2"], k["scan_step"]))
    return proof, mach + om


# Structural constants that may appear in proof.md and machine.md besides computed values and code constants:
# the binomial spread factor 0.5, the exponents 24 (setup allowance 2^24) and 41 (stack top 2^41, toy field
# p(-41)), the MOVZ shift 48, the table digit width 51, the H1 exponent 250, the accumulator bounds 255, 256
# and 259, the encodability ranges 504, 4095 and 32760, and the mutated step cost 1000.
STRUCTURAL = {"0.5", "24", "41", "48", "51", "250", "255", "256", "259", "504", "1000", "4095", "32760"}


def constants_used(L):
    """Code constants a document may quote: charges, coin total, block bytes, parameter bytes, table row
    bytes, C, walks, z exponents, setup layout."""
    fmt = lambda n: f"{int(n):,}"
    vals = list(CHARGE.values()) + [DUP_C, L["coins"], stride(PROD['K']), PARB, PARW, 960, C_EQ, NW, SSTR, SETW,
                                     SMALL_EXP, CNT_ENT.bit_length() - 1, PAR, REP, HCUR - ST, SCRATCH - ST, 2048]
    vals += [a_.bit_length() - 1 for a_ in (RING, CONST, ST, CURB, SBLKA, WIN, BLK, FOLD, TB, CNTB, LIST, SORTB)]
    vals += [v for lo, wd in RDX for v in (lo, wd, lo + wd - 1)]             # radix digit bit ranges
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
    ok_static = isa["all_correct"] and all(prof.values())
    ok = (claimed == Decimal(L["time_log2"]) and not miss and not extra_claim and not unknown and rows_ok and ok_static
          and L["success_above_0.39"] and L["repeats_bound_ok"] and L["output_allowances_ok"] and L["coverage_ge_0.3901"] and L["every_node_below_2^63"] and L["layout_disjoint_below_2^63"])
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
    for key in ("relay_rec", "home_rec", "scan_step", "radix_scatter1"):
        km = dict(k); km[key] = k[key] + 1
        res[key + "+1"] = not stated_figures(ledger(km), km, isa, prof, alg, docs)["stated_figures_match"]
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
    oc = OS.check([1, R - 1, random.Random(seed).randrange(1, R)], k['fp_mul'])
    L = ledger(k)
    docs = read_docs()
    sf = stated_figures(L, k, isa, prof, alg, docs)
    mt = mutation_test(k, isa, prof, alg, docs, sf["stated_figures_match"])
    checks = {"kernels": k["mismatches"] == 0, "transport": tr["transport_ok"],
              "algebra": all(v for v in alg.values() if isinstance(v, bool)), "isa_negative": isa["all_correct"],
              "output_step": oc["output_step_ok"] and oc["identity_gives_null_and_Q0"],
              "profile": all(prof.values()),
              "stated_figures": sf["stated_figures_match"], "mutation_test": mt["all_detected"]}
    if mode in ("all", "kernels"):
        out["kernels"] = k
    if mode in ("all", "transport"):
        out["transport"] = tr
    if mode in ("all", "algebra"):
        out["algebra"], out["profile"], out["isa_negative"], out["output_step"] = alg, prof, isa, oc
    if mode in ("all", "ledger"):
        out["charges"], out["ledger"], out["mutation_test"] = CHARGE, L, mt
    out["stated_figures"] = sf
    if mode in ("all", "toy"):
        T = TY.toy(seed, args.get("toy_targets", 40), args.get("toy_field", "small"), cap=args.get("toy_cap"))
        out["toy"] = T
        checks["toy"] = (T["algorithm"]["wrong_outputs"] == 0 and T["algorithm"]["solved_equals_ground_truth_without_overflow"]
                         and T["algorithm"]["one_walk_pairs_only"]["solved_equals_ground_truth"]
                         and T["coverage"]["every_draw_at_least_bound"] and T["denominators_nonzero"]
                         and T["coverage"]["denominators_nonzero"])
    out["checks"] = checks
    out["all_checks_passed"] = all(checks.values())
    out["scope"] = ("every mode runs every check except the toy and exits with status 1 if one fails; kernel counts are "
                    "executed instructions of generated programs; the simulation runs the same generated "
                    "listings (8 records per block) on 36 nodes; "
                    "the toy runs the whole algorithm at toy scale; no production attack is run")
    print(json.dumps(out, indent=1))
    if not out["all_checks_passed"]:
        print("FAIL: " + ", ".join(c for c, v in checks.items() if not v), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
