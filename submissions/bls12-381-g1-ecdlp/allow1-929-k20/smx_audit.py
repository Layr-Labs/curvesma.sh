"""SMX audit (proof.md s15): exact restart price from the shipped kernel's own Montgomery bodies (measure/side_lanes.S).
Every product and square use-site takes the canonicalizing form (Mc, Sc), so every field value of the routine outside the
bodies is canonical; the plain M/S sequences are emulated as component checks only."""
import json, random
from pathlib import Path
FP15 = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
M_OP, MC_OP, S_OP, SC_OP, C1_OP = 303, 315, 253, 265, 158
A1, A2, SB, X2C, X2 = 29, 41, 29, 26, 38
X3, X4, X8 = X2C + A1, 2 * X2C, 3 * X2C
DBL_OP = 2 * MC_OP + 5 * SC_OP + (A2 + 5 * SB + 2 * X2C + X2 + X3 + X8)      # dbl-2007-bl: 2M + 5S, all canonicalizing, 409 smalls
ADD_OP = 7 * MC_OP + 4 * SC_OP + (A1 + 8 * SB + 2 * X2C + X2 + X4)           # madd-2007-bl: 7M + 4S, all canonicalizing, 403 smalls
FERMAT_OP = 380 * SC_OP + 228 * MC_OP                                        # p-2 fixed: 380 squarings, 228 multiplications
LOOP_OP, SCAN_OP = 76 * 5, 77 * 3 + 2                                        # ladder control; first-set-bit scan
ENTRY_OP, SETUP_OP, CONV_GLUE = 656, 16, 24                                  # bounds over 2 Mc + 18 moves = 648; 12 + 3 = 15; four 6-instruction frame moves
assert ENTRY_OP >= 2 * MC_OP + 18 and SETUP_OP >= 12 + 3 and CONV_GLUE == 4 * 6
CONV_OP = SC_OP + 3 * MC_OP + 2 * C1_OP + CONV_GLUE
SMX = 76 * (DBL_OP + ADD_OP) + FERMAT_OP + LOOP_OP + SCAN_OP + CONV_OP + ENTRY_OP + SETUP_OP
assert (DBL_OP, ADD_OP, FERMAT_OP, CONV_OP, SMX) == (2364, 3668, 172520, 1550, 633787)
def audit_smx():
    import a64
    o = Path(__file__).resolve().parent
    if (o / 'files').is_dir(): o = o / 'files'
    P, RA = FP15, 1 << 384
    N0 = (-pow(P, -1, RA)) & ((1 << 64) - 1)
    tx = (o / 'measure/side_lanes.S').read_text()
    def rg(a, b):
        s = tx[tx.index(a):tx.index(b, tx.index(a))]
        return [l.split('//')[0].strip() for l in s.splitlines() if l.split('//')[0].strip() and not l.split('//')[0].strip().startswith('.')]
    mn_ = rg('// ---- inlined MULNR', '// ---- end MULNR')     # first product: a in x11..x16 (< 2p), b words paired from [x2]
    sq = rg('// ---- inlined SQR', '// ---- end SQR')          # square of x11..x16, 12 words in registers
    m1 = rg('// ---- inlined BY1', '// ---- end BY1')          # word-serial reduction of the low half (n0 in x18, p in x5..x10)
    assert (len(mn_), len(sq), len(m1)) == (295, 91, 150)
    pre = [l.split('//')[0].strip() for l in tx[:tx.index('// ---- inlined MULNR')].splitlines() if l.split('//')[0].strip()][-1]
    assert pre == 'ldp\tx17,x29,[x2]'                             # the first b-word pair load precedes the region
    rs = tx[tx.index('// ---- end BY1'):].splitlines()
    ha = [' '.join(l.split('//')[0].split()) for l in rs[1:] if l.split('//')[0].strip()][:6]
    assert [x.split()[0] for x in ha] == ['adds', 'adcs', 'adcs', 'adcs', 'adcs', 'adc']   # high-half add of the square
    ST, US = {'adds', 'adcs', 'subs', 'sbcs'}, {'adcs', 'adc', 'sbcs', 'sbc', 'csel'}
    for c in (mn_, sq, m1):
        f = False
        for ln in c:
            w = ln.split()[0]
            if w in ST: f = True
            if w in US: assert f                                    # every region sets flags before its first flag use
        for ln in c: assert 'sp' not in ln.split()[1:] and not ln.startswith(('stp', 'str'))   # no stack traffic inside the bodies
    T = ['x%d' % i for i in (11, 12, 13, 14, 15, 16)]
    U = ['x%d' % i for i in (19, 20, 21, 22, 23, 24)]
    R = ['x%d' % i for i in (1, 3, 4, 25, 26, 27)]
    PM = ['x%d' % i for i in (5, 6, 7, 8, 9, 10)]
    L1 = lambda o: ['ldp %s,%s,[sp,#%d]' % (T[2 * i], T[2 * i + 1], o + 16 * i) for i in range(3)]
    L2 = lambda o: L1(o) + ['ldp %s,%s,[sp,#%d]' % (U[2 * i], U[2 * i + 1], o + 48 + 16 * i) for i in range(3)]
    S3 = lambda o: ['stp %s,%s,[sp,#%d]' % (T[2 * i], T[2 * i + 1], o + 16 * i) for i in range(3)]
    def CH(v):
        a = v < 2
        f, g = ('adds', 'adcs') if a else ('subs', 'sbcs')
        s = U if v == 0 else (T if v == 1 else U)
        q = [f + ' %s,%s,%s' % (T[0], T[0], s[0])] + [g + ' %s,%s,%s' % (T[i], T[i], s[i]) for i in range(1, 6)]
        return q + ['adc x0,xzr,xzr'] if a else q
    S6 = ['subs %s,%s,%s' % (R[0], T[0], PM[0])] + ['sbcs %s,%s,%s' % (R[i], T[i], PM[i]) for i in range(1, 6)]
    CL = ['csel %s,%s,%s,lo' % (T[i], T[i], R[i]) for i in range(6)]
    C1 = S6 + ['sbcs xzr,x0,xzr'] + CL
    C2 = S6 + CL
    RP = ['adds %s,%s,%s' % (R[0], T[0], PM[0])] + ['adcs %s,%s,%s' % (R[i], T[i], PM[i]) for i in range(1, 5)] + ['adc %s,%s,%s' % (R[5], T[5], PM[5]), 'cmp x0,#1'] + ['csel %s,%s,%s,eq' % (T[i], T[i], R[i]) for i in range(6)]
    SM = {'a1': L2(0) + CH(0) + C1 + S3(96), 'a2': L2(0) + CH(0) + C1 + C2 + S3(96), 'sb': L2(0) + CH(2) + ['adc x0,xzr,xzr'] + RP + S3(96), 'x2': L1(0) + CH(1) + C1 + C2 + S3(96), 'xc': L1(0) + CH(1) + C1 + S3(96)}
    assert [len(SM[k]) for k in ('a1', 'a2', 'sb', 'x2', 'xc')] == [A1, A2, SB, X2, X2C]
    MG = ['mov x2,sp'] + L1(48) + ['ldp x17,x29,[x2]']           # b at [sp,#0..47] via x2, a at [sp,#48..95]
    MS = S3(96)
    SG = L1(0)
    CG = ['ldr x17,[sp,#0]', 'ldp x19,x20,[sp,#8]', 'ldp x21,x22,[sp,#24]', 'ldr x23,[sp,#40]']
    CS = ['str x17,[sp,#96]', 'stp x19,x20,[sp,#104]', 'stp x21,x22,[sp,#120]', 'str x23,[sp,#136]']
    assert len(MG + MS) == 8 and len(SG + MS) == 6 and len(CG + CS) == 8
    assert len(MG + mn_ + MS) == M_OP and len(MG + mn_ + C2 + MS) == MC_OP and len(SG + sq + m1 + ha + MS) == S_OP and len(SG + sq + m1 + ha + C2 + MS) == SC_OP and len(CG + m1 + CS) == C1_OP
    FB = 0x7000
    MD = [(P >> (64 * i)) & ((1 << 64) - 1) for i in range(6)]
    rm = lambda a, b: (a * b + (((-a * b * pow(P, -1, RA)) % RA) * P)) // RA
    r = random.Random(20261103)
    HP = (3 * P) // 2
    assert 8 * P < RA
    def em(bd, en, tr, swap=False):   # swap: product layout, b at [sp,#0] (via x2) and a at [sp,#48]
        pr, lb = a64.a64_load(en + ':\n' + '\n'.join('\t' + l for l in bd) + '\n\tret\n')
        M64 = (1 << 64) - 1
        def pt(m, a, v):
            for i in range(6): m[a + 8 * i] = (v >> (64 * i)) & M64
        for a, b, w in tr:
            m = {}
            pt(m, FB, b if swap else a); pt(m, FB + 48, a if swap else b)
            q = {f'x{i}': r.getrandbits(64) for i in range(31)}
            q.update(sp=FB, x18=N0)
            for i, v in enumerate(MD): q[f'x{5 + i}'] = v
            ou = {}
            a64.a64_run(pr, lb, en, q, m, limit=3000, out_regs=ou)
            assert sum(m[FB + 96 + 8 * i] << (64 * i) for i in range(6)) == w, (en, hex(a))
            assert ou['x18'] == N0 and ou['sp'] == FB and all(ou[f'x{5 + i}'] == v for i, v in enumerate(MD)), en   # n0, frame and modulus preserved
    ed = [0, 1, 2, P - 2, P - 1, P, P + 1, HP - 1, 2 * P - 1]
    mp = [(r.randrange(2 * P), r.randrange(2 * P)) for _ in range(14)] + [(a, b) for a in ed for b in (1, P - 1, 2 * P - 1)]
    assert all(rm(a, b) < HP for a, b in mp)
    em(MG + mn_ + MS, 'm', [(a, b, rm(a, b)) for a, b in mp], True)
    em(MG + mn_ + C2 + MS, 'mc', [(a, b, rm(a, b) % P) for a, b in mp], True)
    sp_ = [r.randrange(2 * P) for _ in range(14)] + ed
    assert all(rm(a, a) < HP for a in sp_)
    em(SG + sq + m1 + ha + MS, 's', [(a, 0, rm(a, a)) for a in sp_])
    em(SG + sq + m1 + ha + C2 + MS, 'sc', [(a, 0, rm(a, a) % P) for a in sp_])
    em(CG + m1 + CS, 'c1', [(a, 0, rm(a, 1)) for a in sp_])
    for k, hi, fn in (('a1', P, lambda a, b: (a + b) % P), ('a2', HP, lambda a, b: (a + b) % P), ('sb', P, lambda a, b: (a - b) % P), ('x2', HP, lambda a, b: (2 * a) % P), ('xc', P, lambda a, b: (2 * a) % P)):
        ee = sorted(set([0, 1, 2, P - 2, P - 1] + ([P, P + 1, HP - 1] if hi == HP else [])))
        tr = [(r.randrange(hi), r.randrange(hi)) for _ in range(10)] + [(a, b) for a in ee for b in ee]
        em(SM[k], 'q' + k, [(a, b, fn(a, b)) for a, b in tr])
    def AA(p1, q1):
        if p1 is None: return q1
        if q1 is None: return p1
        x1, y1 = p1; x2, y2 = q1
        if x1 == x2 and (y1 + y2) % P == 0: return None
        l = (3 * x1 * x1) * pow(2 * y1, -1, P) % P if x1 == x2 else ((y2 - y1) * pow(x2 - x1, -1, P)) % P
        x3 = (l * l - x1 - x2) % P
        return (x3, (l * (x1 - x3) - y1) % P)
    def AM(k, pt):
        o = None
        while k:
            if k & 1: o = AA(o, pt)
            pt = AA(pt, pt); k >>= 1
        return o
    def LD(b, Qa):
        X, Y, Z = Qa[0], Qa[1], 1
        for i in range(b.bit_length() - 2, -1, -1):
            A = X * X % P; B = Y * Y % P; C = B * B % P
            D = (2 * ((X + B) ** 2 - A - C)) % P; E = (3 * A) % P; F = E * E % P
            X, Y, Z = (F - 2 * D) % P, (E * (D - (F - 2 * D)) - 8 * C) % P, (2 * Y * Z) % P
            if (b >> i) & 1:
                W = Z * Z % P; U2 = Qa[0] * W % P; S2 = Qa[1] * W * Z % P
                H = (U2 - X) % P; HH = H * H % P; I = (4 * HH) % P; J = H * I % P
                t = (2 * (S2 - Y)) % P; V = X * I % P
                X = (t * t - J - 2 * V) % P
                X, Y, Z = X, (t * (V - X) - 2 * Y * J) % P, ((Z + H) ** 2 - W - HH) % P
        zi = pow(Z, -1, P); z2 = zi * zi % P
        return (X * z2 % P, Y * z2 * zi % P)
    qm = json.loads((o / 'params.json').read_text())
    G = (int(qm['gx'], 16), int(qm['gy'], 16))
    assert (G[1] ** 2 - G[0] ** 3 - 4) % P == 0
    for t in range(24):
        Qa = AM(r.randrange(1, 1 << 20), G)
        b = r.getrandbits(77) or 1
        assert LD(b, Qa) == AM(b, Qa), t
    for b in (1, 2, (1 << 76), (1 << 77) - 1):
        Qa = AM(r.randrange(1, 1 << 20), G)
        assert LD(b, Qa) == AM(b, Qa), b
    RR = pow(2, 384, P); RI = pow(RR, -1, P)
    mn = lambda a, b: (a * b % P) * RI % P
    e = P - 2
    assert (e.bit_length() - 1, bin(e).count('1') - 1) == (380, 228)
    for t in range(6):
        Z = r.randrange(1, P)
        ac = Z * RR % P
        for i in range(e.bit_length() - 2, -1, -1):
            ac = mn(ac, ac)
            if (e >> i) & 1: ac = mn(ac, Z * RR % P)
        assert ac == pow(Z, -1, P) * RR % P, t
    return {'SMX': SMX, 'M': M_OP, 'Mc': MC_OP, 'S': S_OP, 'Sc': SC_OP, 'C1': C1_OP, 'DBL': DBL_OP, 'ADD': ADD_OP, 'Fermat': FERMAT_OP, 'conv': CONV_OP, 'loop': LOOP_OP, 'scan': SCAN_OP, 'entry': ENTRY_OP, 'setup': SETUP_OP}
if __name__ == '__main__': print(json.dumps(audit_smx()))
