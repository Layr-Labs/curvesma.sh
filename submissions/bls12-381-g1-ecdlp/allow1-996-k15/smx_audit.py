"""SMX audit (s15)."""
import json, random
from pathlib import Path
FP15 = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
M_OP, MC_OP, S_OP, SC_OP, C1_OP = 325, 337, 257, 269, 158
A1, A2, SB, X2C, X2 = 29, 41, 29, 26, 38
X3, X4, X8 = X2C + A1, 2 * X2C, 3 * X2C
DBL_OP = M_OP + MC_OP + S_OP + 4 * SC_OP + (A2 + 5 * SB + 2 * X2C + X2 + X3 + X8)
ADD_OP = 2 * M_OP + 5 * MC_OP + 4 * SC_OP + (A1 + 8 * SB + 2 * X2C + X2 + X4)
FERMAT_OP = 380 * S_OP + 228 * M_OP
CONV_OP = S_OP + 3 * M_OP + 2 * C1_OP + 24
SMX = 76 * (DBL_OP + ADD_OP) + FERMAT_OP + 380 + 233 + CONV_OP + 656 + 16
assert (DBL_OP, ADD_OP, FERMAT_OP, CONV_OP, SMX) == (2404, 3814, 171760, 1572, 647185)
def audit_smx():
    import a64
    from pathlib import Path
    o = Path(__file__).resolve().parent
    if (o / 'files').is_dir(): o = o / 'files'
    P, RA = FP15, 1 << 384
    N0 = (-pow(P, -1, RA)) & ((1 << 64) - 1)
    tx = (o / 'side_k2.S').read_text()
    def rg(a, b):
        s = tx[tx.index(a):tx.index(b, tx.index(a))]
        return [l.split('//')[0].strip() for l in s.splitlines() if l.split('//')[0].strip()]
    lz = rg('// ---- lazy first Montgomery product', '// ---- end __mul_mont_384')
    sq = rg('// ---- inlined blst __sqr_384', '// ---- end __sqr_384')
    m1 = rg('// ---- inlined blst __mul_by_1_mont_384', '// ---- end __mul_by_1_mont_384')
    assert (len(lz), len(sq), len(m1)) == (317, 92, 150)
    rs = tx[tx.index('// ---- end __mul_by_1_mont_384'):].splitlines()
    ha = [' '.join(l.split('//')[0].split()) for l in rs[1:] if l.split('//')[0].strip()][:6]
    assert [x.split()[0] for x in ha] == ['adds', 'adcs', 'adcs', 'adcs', 'adcs', 'adc']
    ST, US = {'adds', 'adcs', 'subs', 'sbcs'}, {'adcs', 'adc', 'sbcs', 'sbc', 'csel'}
    for c in (lz, sq, m1):
        f = False
        for ln in c:
            w = ln.split()[0]
            if w in ST: f = True
            if w in US: assert f
    T = ['x%d' % i for i in (11, 12, 13, 14, 15, 16)]
    U = ['x%d' % i for i in (19, 20, 21, 22, 23, 24)]
    R = ['x%d' % i for i in (1, 3, 4, 25, 26, 27)]
    PM = ['x%d' % i for i in (5, 6, 7, 8, 9, 10)]
    L1 = lambda o: ['ldp %s,%s,[x29,#%d]' % (T[2 * i], T[2 * i + 1], o + 16 * i) for i in range(3)]
    L2 = lambda o: L1(o) + ['ldp %s,%s,[x29,#%d]' % (U[2 * i], U[2 * i + 1], o + 48 + 16 * i) for i in range(3)]
    S3 = lambda o: ['stp %s,%s,[x29,#%d]' % (T[2 * i], T[2 * i + 1], o + 16 * i) for i in range(3)]
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
    MG = ['add x2,x29,#48', 'ldp x19,x20,[x29,#0]', 'ldp x21,x22,[x29,#16]', 'ldp x23,x24,[x29,#32]', 'ldr x17,[x2]']
    MS = ['stp x11,x12,[x29,#96]', 'stp x13,x14,[x29,#112]', 'stp x15,x16,[x29,#128]']
    SG = ['ldp x11,x12,[x29,#0]', 'ldp x13,x14,[x29,#16]', 'ldp x15,x16,[x29,#32]']
    SR = ['ldp x5,x6,[x29,#144]', 'ldp x7,x8,[x29,#160]', 'ldp x9,x10,[x29,#176]']
    CG = ['ldr x17,[x29,#0]', 'ldp x19,x20,[x29,#8]', 'ldp x21,x22,[x29,#24]', 'ldr x23,[x29,#40]']
    CS = ['str x17,[x29,#96]', 'stp x19,x20,[x29,#104]', 'stp x21,x22,[x29,#120]', 'str x23,[x29,#136]']
    assert len(MG + MS) == 8 and len(SG + SR + MS) == 9 and len(CG + CS) == 8
    assert len(MG + lz + MS) == M_OP and len(MG + lz + C2 + MS) == MC_OP and len(SG + sq + SR + m1 + ha + MS) == S_OP and len(SG + sq + SR + m1 + ha + C2 + MS) == SC_OP and len(CG + m1 + CS) == C1_OP
    FB, XW = 0x7000, 0x8000
    MD = [(P >> (64 * i)) & ((1 << 64) - 1) for i in range(6)]
    rm = lambda a, b: (a * b + (((-a * b * pow(P, -1, RA)) % RA) * P)) // RA
    r = random.Random(20261103)
    HP = (3 * P) // 2
    assert 8 * P < RA
    def em(bd, en, tr):
        pr, lb = a64.a64_load(en + ':\n' + '\n'.join('\t' + l.replace('8*1', '8').replace('8*2', '16').replace('8*3', '24').replace('8*4', '32').replace('8*5', '40') for l in bd) + '\n\tret\n')
        M64 = (1 << 64) - 1
        def pt(m, a, v):
            for i in range(6): m[a + 8 * i] = (v >> (64 * i)) & M64
        for a, b, w in tr:
            m = {XW + 8 * i: 0xBEEF for i in range(5)}
            pt(m, FB, a); pt(m, FB + 48, b); pt(m, FB + 144, P)
            q = {f'x{i}': r.getrandbits(64) for i in range(31)}
            q.update(x29=FB, x30=N0, sp=0x9000, x18=XW + 32)
            for i, v in enumerate(MD): q[f'x{5 + i}'] = v
            ou = {}
            a64.a64_run(pr, lb, en, q, m, limit=3000, out_regs=ou)
            assert sum(m[FB + 96 + 8 * i] << (64 * i) for i in range(6)) == w and ou['x30'] == N0, (en, hex(a))
    ed = [0, 1, 2, P - 2, P - 1, P, P + 1, HP - 1, 2 * P - 1]
    mp = [(r.randrange(2 * P), r.randrange(2 * P)) for _ in range(14)] + [(a, b) for a in ed for b in (1, P - 1)]
    assert all(rm(a, b) < HP for a, b in mp)
    em(MG + lz + MS, 'm', [(a, b, rm(a, b)) for a, b in mp])
    em(MG + lz + C2 + MS, 'mc', [(a, b, rm(a, b) % P) for a, b in mp])
    sp = [r.randrange(2 * P) for _ in range(14)] + ed
    assert all(rm(a, a) < HP for a in sp)
    em(SG + sq + SR + m1 + ha + MS, 's', [(a, 0, rm(a, a)) for a in sp])
    em(SG + sq + SR + m1 + ha + C2 + MS, 'sc', [(a, 0, rm(a, a) % P) for a in sp])
    em(CG + m1 + CS, 'c1', [(a, 0, rm(a, 1)) for a in sp])
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
