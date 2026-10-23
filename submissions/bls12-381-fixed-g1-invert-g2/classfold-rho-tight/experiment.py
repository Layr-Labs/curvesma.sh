"""Class-folded distinguished-point rho in GT (BLS12-381, fixed-G1 pairing inversion).

Standard library only.  Reads one JSON object from stdin and prints one JSON object.
Modes (each run is independent and bounded; `all` runs the four parts in one process):
  ledger   exact integer/rational work ledger, trap-lemma bound, memory ledger
  algebra  production-field checks on real BLS12-381 GT elements (selector, Frobenius)
  toy      end-to-end algorithm in a small tower of the same shape (real Fp12 arithmetic)
  sim      exponent-domain simulation of the identical walk (walk-quality constant)
Nothing here is a production-scale run; see proof.md section 9 for the scope of each mode.
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

M64 = (1 << 64) - 1
P_FP = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORD = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
MU = Fraction(11, 10)                                    # allowed walk inefficiency (H1), unchanged
P_MAIN = Fraction(392, 1000)                             # main-event probability at S_max under H1
LN_MAIN_UP = Fraction(9951607940319401, 2 * 10 ** 16)    # ln(1000/608) rounded up (toy budget only)
T_CLAIM = float(MU) * math.log(1 / (1 - float(P_MAIN)))  # S_max^2/(2N) = mu ln(1000/608) = 0.54734
# Per-slot work other than the reference multiplication, itemised as hard bounds (proof.md section 5).
SLOT_ITEMS = (
    ('selector: 4 blst_fp_mul at 383', 4 * 383),
    ('selector: 9 modular add/sub/neg at 40', 9 * 40),
    ('selector: 11 six-limb compare+cmov at 25', 11 * 25),
    ('khash: two 6-limb words at 70 (7 mixes of 8, 6 xors, limb loads; word 1 = index and 64 DP bits, word 2 = 16)', 2 * 70),
    ('distinguished-point test: two 64-bit compares and a branch', 8),
    ('look-back: 5 binary-search steps, index byte load + six-limb compare at 26', 5 * 26),
    ('look-back: ring write of 6 key words and slot address', 8),
    ('look-back: rotation of at most 23 one-byte index entries at 4', 23 * 4),
    ('look-back: inverse-map load, two stores, rank fixup, branch', 7),
    ('look-back: window bookkeeping (slot counter, length with its not-yet-full comparison, branch)', 32),
    ('table address (w1 mod 2^14, (h+a) mod R, 576-byte stride), counters, loop control', 40),
)
C_ITEM = sum(n for _, n in SLOT_ITEMS)                   # itemised per-slot instructions
C_SLOT = 2638                                            # charged: itemised plus an explicit margin
assert C_ITEM == 2624, C_ITEM
assert C_SLOT - C_ITEM == 14, "explicit margin (0.5%), as proof.md section 5 states"

# ---------------------------------------------------------------- hashing of keys
def mix(x):
    x &= M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)

def khash(kappa, salt, nlimbs):
    """64-bit fixed mixing of the nlimbs 64-bit limbs of a class key."""
    acc = mix(salt)
    for i in range(nlimbs):
        acc = mix(acc ^ ((kappa >> (64 * i)) & M64))
    return acc

# ---------------------------------------------------------------- tower Fp2 / Fp12
def make_field(p):
    """Fp2 = Fp[u]/(u^2+1); Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u; element = 6 Fp2 coefficients of w^i."""
    xi = (1, 1)
    def f2mul(a, b):
        t0 = a[0] * b[0]; t1 = a[1] * b[1]
        return ((t0 - t1) % p, ((a[0] + a[1]) * (b[0] + b[1]) - t0 - t1) % p)
    def f2conj(a): return (a[0], (-a[1]) % p)
    def f2pow(a, e):
        r = (1, 0)
        while e:
            if e & 1: r = f2mul(r, a)
            a = f2mul(a, a); e >>= 1
        return r
    def f12mul(x, y):
        t = [[0, 0] for _ in range(11)]
        for i in range(6):
            xa, xb = x[i]
            for j in range(6):
                ya, yb = y[j]
                t0 = xa * ya; t1 = xb * yb
                c = t[i + j]
                c[0] += t0 - t1
                c[1] += (xa + xb) * (ya + yb) - t0 - t1
        out = []
        for k in range(6):
            a = t[k]
            if k + 6 < 11:
                b = t[k + 6]                       # multiply by xi = 1 + u
                out.append(((a[0] + b[0] - b[1]) % p, (a[1] + b[0] + b[1]) % p))
            else:
                out.append((a[0] % p, a[1] % p))
        return tuple(out)
    one12 = ((1, 0),) + ((0, 0),) * 5
    def f12pow(x, e):
        r = one12
        for bit in bin(e)[2:]:
            r = f12mul(r, r)
            if bit == '1': r = f12mul(r, x)
        return r
    gamma = f2pow(xi, (p - 1) // 6)
    gam = [(1, 0)]
    for _ in range(5): gam.append(f2mul(gam[-1], gamma))
    def frob(x):                                    # x -> x^p
        return tuple(f2mul(f2conj(x[i]), gam[i]) for i in range(6))
    zeta = f2mul(gamma, f2conj(gamma))
    assert zeta[1] == 0 and (zeta[0] * zeta[0] - zeta[0] + 1) % p == 0
    return dict(p=p, mul=f12mul, pow=f12pow, one=one12, frob=frob, gamma=gamma, zeta=zeta[0],
                f2mul=f2mul, f2pow=f2pow, nlimbs=(p.bit_length() + 63) // 64)

def frobk(F, x, k):
    for _ in range(k % 12): x = F['frob'](x)
    return x

def sel(F, y):
    """Class selector: 4 Fp multiplications.  Returns (kappa, kstar).  kappa = min over the 12
    first-coordinates of c1(sigma^k(y)); kstar = the k attaining it (lowest k on a tie)."""
    p = F['p']; zeta = F['zeta']; g0, g1 = F['gamma']
    a, b = y[1]
    za = zeta * a % p                      # mult 1
    d = (g0 * a + g1 * b) % p              # mults 2, 3
    zd = zeta * d % p                      # mult 4
    h2 = (za - a) % p; hd = (zd - d) % p   # zeta^2 = zeta - 1
    V = (a, d, za, zd, h2, hd, (-a) % p, (-d) % p, (-za) % p, (-zd) % p, (-h2) % p, (-hd) % p)
    # V[2j] = zeta^j a (j=0..2), V[2j+1] = zeta^j d, and V[k+6] = -V[k]
    best = V[0]; kb = 0
    for k in range(1, 12):
        if V[k] < best: best = V[k]; kb = k
    return best, kb

# ---------------------------------------------------------------- exact ledger
def isqrt_ceil(n):
    r = math.isqrt(n)
    return r if r * r == n else r + 1

def trap_bound(R, ell_min, ell_max):
    """E(l) = sum_{d<=min(l//2,R)} C(R,d) (d/R)^l 12^-d, evaluated in log2; returns (max, sum) over l."""
    lg = math.lgamma
    def lc(n, k): return (lg(n + 1) - lg(k + 1) - lg(n - k + 1)) / math.log(2)
    mx = -1e9; tot = 0.0
    for ell in range(ell_min, ell_max + 1):
        terms = []
        for d in range(1, min(ell // 2, R) + 1):
            terms.append(lc(R, d) + ell * math.log2(d / R) - d * math.log2(12))
        m = max(terms)
        s = m + math.log2(sum(2.0 ** (t - m) for t in terms))
        mx = max(mx, s); tot += 2.0 ** s
    return mx, tot

def ledger():
    getcontext().prec = 80
    N = (R_ORD - 1) // 12
    mu = MU
    # main-event probability 392/1000:  S^2 >= -2 mu N ln(1 - 392/1000) = 2 mu N ln(1000/608)
    lnv = (Decimal(P_MAIN.denominator) / Decimal(P_MAIN.denominator - P_MAIN.numerator)).ln()
    s2 = Decimal(2) * Decimal(mu.numerator) / Decimal(mu.denominator) * Decimal(N) * lnv
    S_max = int(s2.sqrt().to_integral_value(rounding="ROUND_CEILING"))
    assert Fraction(S_max) ** 2 >= 2 * mu * N * Fraction(lnv) and Fraction(S_max - 1) ** 2 < 2 * mu * N * Fraction(lnv)
    M_max = S_max + (S_max >> 13) + 1                    # probes allowed: 2^-13 of the accepted steps
    c_over = C_SLOT                                      # instructions per Fp12 multiplication slot (section 5)
    unit = Fraction(20000 + c_over, 20000)               # 11319/10000 per multiplication slot
    R = 1 << 14
    theta_bits = 80
    L_cap = 1 << 87
    D_max = 1 << 47                                      # reports cap
    ST_max = 1 << 47                                     # starts cap
    REC_max = 16                                         # recoveries cap
    add = {}
    add['pairings (g and one spare), 2*2^30 instr'] = Fraction(2 * (1 << 30), 20000)
    add['table: R*(1100 + 12) units'] = R * 1112
    add['starts: ST_max*1100 units'] = ST_max * 1100
    add['reports: D_max*2^16 instr'] = Fraction(D_max * (1 << 16), 20000)
    add['replays: REC_max*2*L_cap*(5/4)'] = Fraction(REC_max * 2 * L_cap * 5, 4)
    add['finish last trail: L_cap*unit'] = L_cap * unit
    add['G2 scalar mult + GT check (<= 2^12 units)'] = 1 << 12
    A = sum(add.values())
    W = unit * M_max + A
    # exact log2 with Decimal
    log2W = (Decimal(W.numerator) / Decimal(W.denominator)).ln() / Decimal(2).ln()
    claim = (log2W * 100).to_integral_value(rounding="ROUND_CEILING") / 100
    c100 = int(claim * 100)
    exact_le = W.numerator ** 100 <= (1 << c100) * W.denominator ** 100           # W^100 <= 2^(100*claim)
    exact_gt_prev = W.numerator ** 100 > (1 << (c100 - 1)) * W.denominator ** 100  # W^100 > 2^(100*claim - 1)
    pre = add['table: R*(1100 + 12) units'] + add['pairings (g and one spare), 2*2^30 instr']
    mem = {
        'reports_bytes_max': D_max * 96,                # key 48 B + seed 16 B + 3 pointers + colour
        'twist_table_bytes': 12 * R * 576,
        'uv_scalars_bytes': 1 << 20,
        'window_and_state_bytes': 1 << 20,
    }
    peak = sum(mem.values())
    trap_max, trap_sum = trap_bound(R, 25, 512)
    expected_dp = Decimal(S_max) * Decimal(2) ** (-theta_bits)
    return {
        'log2_r': float(Decimal(R_ORD).ln() / Decimal(2).ln()), 'log2_N_classes': float(Decimal(N).ln() / Decimal(2).ln()),
        'mu': str(mu), 'S_max': S_max, 'log2_S_max': float(Decimal(S_max).ln() / Decimal(2).ln()),
        'M_max': M_max, 'unit_per_mult': str(unit),
        'additive_items_log2': {k: float((Decimal(v.numerator if isinstance(v, Fraction) else v) / Decimal(v.denominator if isinstance(v, Fraction) else 1)).ln() / Decimal(2).ln()) for k, v in add.items()},
        'log2_additive_total': float((Decimal(A.numerator if isinstance(A, Fraction) else A) / Decimal(A.denominator if isinstance(A, Fraction) else 1)).ln() / Decimal(2).ln()),
        'log2_W_max': str(log2W)[:12], 'claim_time_log2': '%.2f' % claim,
        'claim_check_W_le_2^claim': bool(log2W <= claim),
        'exact_integer_check_W^100_le_2^%d' % c100: exact_le,
        'exact_integer_check_W^100_gt_2^%d' % (c100 - 1): exact_gt_prev,
        'slot_items_instructions': {k: n for k, n in SLOT_ITEMS},
        'slot_instructions_itemised': C_ITEM, 'slot_instructions_charged': c_over,
        'slot_margin_instructions': c_over - C_ITEM,
        'log2_preprocessing': float((Decimal(pre.numerator if isinstance(pre, Fraction) else pre) / Decimal(pre.denominator if isinstance(pre, Fraction) else 1)).ln() / Decimal(2).ln()),
        'expected_reports_log2': float(expected_dp.ln() / Decimal(2).ln()), 'report_cap_log2': 47,
        'peak_memory_log2_bytes': math.log2(peak),
        'trap_lemma_log2_max_E_l_for_l_25_to_512': trap_max, 'trap_lemma_log2_sum': math.log2(trap_sum),
        'trap_cost_exponent_log2(L_cap*sum)': 87 + math.log2(trap_sum),
        'success_main_event': '392/1000 under mu', 'success_claim': '39/100',
    }

# ---------------------------------------------------------------- algebra (production field)
def algebra(elements, seed):
    F = make_field(P_FP); p = P_FP
    rng = random.Random(seed)
    lam = p % R_ORD
    out = {'lambda_order_12': pow(lam, 12, R_ORD) == 1 and pow(lam, 6, R_ORD) == R_ORD - 1 and pow(lam, 4, R_ORD) != 1,
           'zeta_primitive_6th_root_in_Fp': (F['zeta'] ** 6) % p == 1 and (F['zeta'] ** 3) % p == p - 1,
           'gamma_hex': [format(F['gamma'][0], 'x'), format(F['gamma'][1], 'x')], 'zeta_hex': format(F['zeta'], 'x'),
           'checked_elements': 0}
    e = (p ** 12 - 1) // R_ORD
    first = True
    for _ in range(elements):
        x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
        y = F['pow'](x, e)
        assert F['pow'](y, R_ORD) == F['one'] and y != F['one']
        if first:
            assert F['frob'](y) == F['pow'](y, p) == F['pow'](y, lam)      # Frobenius is x -> x^lambda on GT
            first = False
        kap0, k0 = sel(F, y)
        for k in range(12):                                                  # candidates are the true coordinates
            assert frobk(F, y, k)[1][0] == (sel_all(F, y)[k])
        for m in range(12):
            ym = frobk(F, y, m)
            kap, km = sel(F, ym)
            assert kap == kap0 and km == (k0 - m) % 12
        t = F['pow'](y, rng.randrange(2, R_ORD))
        for m in range(12):                                                  # twisted step is class-equivariant
            ym = frobk(F, y, m); _, km = sel(F, ym)
            y1 = F['mul'](ym, frobk(F, t, (-km) % 12))
            y0 = F['mul'](y, frobk(F, t, (-k0) % 12))
            assert any(frobk(F, y0, j) == y1 for j in range(12))
        out['checked_elements'] += 1
    out['all_assertions_passed'] = True
    return out

def sel_all(F, y):
    p = F['p']; zeta = F['zeta']; g0, g1 = F['gamma']
    a, b = y[1]; za = zeta * a % p; d = (g0 * a + g1 * b) % p; zd = zeta * d % p
    h2 = (za - a) % p; hd = (zd - d) % p
    return (a, d, za, zd, h2, hd, (-a) % p, (-d) % p, (-za) % p, (-zd) % p, (-h2) % p, (-hd) % p)

# ---------------------------------------------------------------- number-theory helpers (toy)
def is_prime(n):
    if n < 2: return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0: return n == q
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def find_order(bits, rng):
    while True:
        c = rng.randrange(1 << (bits - 1), 1 << bits) // 12 * 12 + 1
        if is_prime(c): return c

def prim12(r, rng):
    while True:
        lam = pow(rng.randrange(2, r - 1), (r - 1) // 12, r)
        if pow(lam, 6, r) == r - 1 and pow(lam, 4, r) != 1: return lam

def toy_field(rp, lam, qbits=60):
    """prime q = lam (mod rp), q = 7 (mod 12), xi = 1+u neither a square nor a cube in Fq2."""
    s0 = ((7 - lam) * pow(rp, -1, 12)) % 12
    j = (1 << qbits) // (12 * rp)
    while True:
        q = lam + rp * (s0 + 12 * j)
        j += 1
        if q % 12 == 7 and is_prime(q) and pow(2, (q - 1) // 2, q) == q - 1:       # N(xi)=2 non-residue
            F = make_field_checked(q)
            if F: return F

def make_field_checked(q):
    try:
        F = make_field(q)
    except AssertionError:
        return None
    one2 = (1, 0)
    if F['f2pow']((1, 1), (q * q - 1) // 3) == one2: return None
    return F

class Walk:
    """Shared parameters of the class-folded walk.  The emulation draws the table index and the distinguished-point
    bits from independent salts, the idealised case of the production rule (where the index is the low 14 bits of the
    first word, proof.md section 5): a distinguished state is reported and never stepped from, so that reuse does not
    enter the walk statistics modelled here."""
    def __init__(self, nlimbs, R, w, db, A=8, seed=0):
        rng = random.Random(seed)
        self.R, self.w, self.A, self.db, self.nl = R, w, A, db, nlimbs
        self.s_h = rng.getrandbits(64); self.s_dp = [rng.getrandbits(64) for _ in range((db + 63) // 64)]
    def hidx(self, kappa): return khash(kappa, self.s_h, self.nl) % self.R
    def is_dp(self, kappa):
        left = self.db
        for s in self.s_dp:
            take = min(64, left); left -= take
            if khash(kappa, s, self.nl) & ((1 << take) - 1): return False
        return True

# ---------------------------------------------------------------- toy: real tower, end to end
def toy(bits, trials, seed):
    rng = random.Random(seed)
    res = []
    for t in range(trials):
        rp = find_order(bits, rng); lam = prim12(rp, rng)
        F = toy_field(rp, lam)
        q = F['p']
        assert (q % rp) == lam and pow(lam, 4, rp) != 1
        e = (q ** 12 - 1) // rp
        while True:
            x = tuple((rng.randrange(q), rng.randrange(q)) for _ in range(6))
            g = F['pow'](x, e)
            if g != F['one']: break
        a = rng.randrange(1, rp)
        z = F['pow'](g, a)
        assert F['pow'](z, rp) == F['one']
        N = (rp - 1) // 12
        res.append(toy_instance(F, g, z, a, rp, lam, N, rng.getrandbits(32), rng))
    ok = sum(r['recovered'] for r in res)
    return {'trials': trials, 'bits': bits, 'recovered_correct': ok,
            'wrong_outputs': sum(r['wrong'] for r in res),
            'within_budget': sum(1 for r in res if r['recovered'] and r['S'] <= r['S_max']),
            'mean_steps_over_sqrt_N': round(sum(r['S'] for r in res) / trials / math.sqrt(sum(r['N'] for r in res) / trials), 3),
            'replay_matched_original_path': all(r['replay_ok'] for r in res),
            'probes_total': sum(r['probes'] for r in res), 'escapes_total': sum(r['escapes'] for r in res),
            'dead_trails': sum(r['dead'] for r in res), 'budget_S_max_mean': round(sum(r['S_max'] for r in res) / trials, 1),
            'scope': 'toy tower of the production shape (q ~ 2^60, r ~ 2^%d); no production cost or success claim' % bits}

def toy_instance(F, g, z, a, rp, lam, N, wseed, rng, R=64, w=12, db=3):
    mul = F['mul']; nl = F['nlimbs']
    W = Walk(nl, R, w, db, seed=wseed)
    lamp = [pow(lam, k, rp) for k in range(12)]
    lami = [pow(lam, (-k) % 12, rp) for k in range(12)]
    U = [rng.randrange(rp) for _ in range(R)]; V = [rng.randrange(rp) for _ in range(R)]
    TW = []
    for h in range(R):
        Th = mul(F['pow'](g, U[h]), F['pow'](z, V[h]))
        fr = [Th]
        for _ in range(11): fr.append(F['frob'](fr[-1]))
        TW.append([fr[(12 - k) % 12] for k in range(12)])               # TW[h][k] = sigma^-k(T_h)
    S_max = isqrt_ceil(int(2 * MU * N * LN_MAIN_UP))
    stats = dict(probes=0, escapes=0, dead=0)

    def trail(tid, track):
        r2 = random.Random(wseed * 7919 + tid)
        al = r2.randrange(rp); be = r2.randrange(rp)
        y = mul(F['pow'](g, al), F['pow'](z, be))
        kap, k = sel(F, y)
        win = [kap]; n = 0; steps = 0; slots = 0
        while True:
            steps += 1
            if W.is_dp(kap): return ('dp', kap, y, al, be, steps)
            if slots >= 64 * (1 << db): return ('cap', kap, y, al, be, steps)
            h0 = W.hidx(kap); ok = False
            for ai in range(W.A):
                hh = (h0 + ai) % R
                slots += 1
                stats['probes'] += 1 if track is False else 0
                y2 = mul(y, TW[hh][k])
                kap2, k2 = sel(F, y2)
                if kap2 in win:
                    stats['escapes'] += 1 if track is False else 0
                    continue
                ok = True; break
            if not ok:
                stats['dead'] += 1 if track is False else 0
                return ('dead', kap, y, al, be, steps)
            if track:
                al = (al + lami[k] * U[hh]) % rp; be = (be + lami[k] * V[hh]) % rp
            y, kap, k = y2, kap2, k2
            win.append(kap)
            if len(win) > W.w: win.pop(0)
            n += 1

    store = {}; S = 0; tid = 0; recovered = False; wrong = 0; replay_ok = True
    while True:
        tid += 1
        status, kap, y, al, be, steps = trail(tid, False)
        S += steps
        if status == 'dp':
            if kap in store:
                t1 = store[kap]
                A_ = trail_track(trail, t1); B_ = trail_track(trail, tid)
                # replay reproduces the untracked path (same endpoint and key)
                replay_ok = replay_ok and A_[1] == kap and B_[1] == kap
                yA, alA, beA = A_[2], A_[3], A_[4]; yB, alB, beB = B_[2], B_[3], B_[4]
                for m in range(12):
                    if frobk(F, yA, m) == yB:
                        lm = lamp[m]; den = (lm * beA - beB) % rp
                        if den:
                            cand = ((alB - lm * alA) * pow(den, -1, rp)) % rp
                            if F['pow'](g, cand) == z:
                                recovered = True
                                if cand != a: wrong += 1
                        break
                if recovered: break
            else:
                store[kap] = tid
        if S >= 6 * S_max: break
    return dict(recovered=recovered, wrong=wrong, S=S, S_max=S_max, N=N, replay_ok=replay_ok, **stats)

def trail_track(trail, tid): return trail(tid, True)

# ---------------------------------------------------------------- sim: exponent-domain walk
def sim(bits, R, trials, seed, w=24, db=5):
    rng = random.Random(seed)
    tot = 0.0; n_first = 0; esc = 0; att = 0; dead = 0; merged = 0; coalesced = 0; within = 0
    for _ in range(trials):
        rp = find_order(bits, rng); lam = prim12(rp, rng)
        N = (rp - 1) // 12
        W = Walk(1, R, w, db, seed=rng.getrandbits(32))
        salt = rng.getrandbits(64)
        ip = [pow(lam, (-k) % 12, rp) for k in range(12)]
        a = rng.randrange(1, rp)
        U = [rng.randrange(rp) for _ in range(R)]; V = [rng.randrange(rp) for _ in range(R)]
        Mh = [(U[h] + a * V[h]) % rp for h in range(R)]
        def ssel(x):
            best = None; kb = 0; xk = x
            for k in range(12):
                hv = mix(xk ^ salt)
                if best is None or hv < best: best = hv; kb = k
                xk = xk * lam % rp
            return best, kb
        visited = {}; end_key = {}; steps = 0; first = None; tid = 0; limit = int(12 * math.sqrt(N))
        follow = None                                   # (earlier trail id) while the merged trail is followed to its end
        while (first is None or follow is not None) and steps < limit + 4 * (1 << db) * 16:
            tid += 1
            x = rng.randrange(1, rp)
            kap, k = ssel(x); win = [kap]; sl = 0
            while True:
                if first is None:
                    steps += 1
                    pv = visited.get(kap)
                    if pv is not None and pv != tid:
                        first = steps; follow = pv
                    else:
                        visited.setdefault(kap, tid)
                if W.is_dp(kap) or sl >= 16 * (1 << db) or (first is None and steps >= limit): break
                h0 = W.hidx(kap); good = False
                for ai in range(W.A):
                    hh = (h0 + ai) % R
                    att += 1; sl += 1
                    x2 = (x + ip[k] * Mh[hh]) % rp
                    kap2, k2 = ssel(x2)
                    if kap2 in win: esc += 1; continue
                    good = True; break
                if not good: dead += 1; break
                x, kap, k = x2, kap2, k2
                win.append(kap)
                if len(win) > w: win.pop(0)
            end_key[tid] = kap if W.is_dp(kap) else None
            if follow is not None:
                if tid != follow:
                    merged += 1
                    if end_key[tid] is not None and end_key[tid] == end_key.get(follow): coalesced += 1
                    follow = None
        if first is not None:
            n_first += 1; tot += first * first / (2 * N)
            if first * first / (2 * N) <= T_CLAIM: within += 1
    # control: plain birthday on N classes with the same estimator
    ctl = 0.0; cn = min(trials, 300)
    for _ in range(cn):
        rp = find_order(bits, rng); N = (rp - 1) // 12; seen = set(); s = 0
        while True:
            v = rng.randrange(N); s += 1
            if v in seen: break
            seen.add(v)
        ctl += s * s / (2 * N)
    return {'bits': bits, 'R': R, 'w': w, 'dp_bits': db, 'trials': trials, 'repeats_found': n_first,
            'mu_hat_first_inter_trail_repeat': round(tot / max(1, n_first), 3),
            'std_error_approx': round(1 / math.sqrt(max(1, n_first)), 3),
            'groups_first_repeat_within_S_max_budget': within,
            'fraction_first_repeat_within_S_max_budget': round(within / trials, 4),
            'budget_threshold_S^2/(2N)': round(T_CLAIM, 5),
            'ideal_map_fraction_1-(608/1000)^1.1': round(1 - 0.608 ** 1.1, 4), 'required_fraction': 0.3903,
            'control_plain_birthday_mu_hat': round(ctl / cn, 3), 'control_trials': cn,
            'merged_trails_followed': merged, 'merged_trails_reaching_the_same_DP': coalesced,
            'escape_rate_per_attempt': esc / max(1, att), 'predicted_two_cycle_rate_1_over_12R': 1 / (12 * R),
            'dead_trails': dead, 'ideal_birthday_value': 1.0,
            'scope': 'exponent-domain emulation of the identical walk; keys are pseudorandom hashes of the 12 conjugates'}

# ---------------------------------------------------------------- main
def main():
    cfg = json.load(sys.stdin)
    mode = cfg['mode']
    if mode == 'ledger': out = ledger()
    elif mode == 'algebra': out = algebra(cfg['elements'], cfg['seed'])
    elif mode == 'toy': out = toy(cfg['bits'], cfg['trials'], cfg['seed'])
    elif mode == 'sim': out = sim(cfg['bits'], cfg['R'], cfg['trials'], cfg['seed'], db=cfg.get('db', 5))
    elif mode == 'all':
        out = {'ledger': ledger(),
               'algebra': algebra(cfg['algebra_elements'], cfg['seed']),
               'toy': toy(cfg['toy_bits'], cfg['toy_trials'], cfg['seed']),
               'sim': sim(cfg['sim_bits'], cfg['sim_R'], cfg['sim_trials'], cfg['seed'], db=cfg['sim_db'])}
    else: raise SystemExit('unknown mode')
    print(json.dumps(out, indent=1))

if __name__ == '__main__':
    main()
