"""Fraction-label backward kernel (section 3): derivation from the rolled base and a bit-exact model.

The batch-inversion part of pair_bwd_y0_c0 (i = B*A_(l-1), then B <- B*d with its operand reuse) is kept
verbatim; the tail after the B store is replaced by: w = M(i,i); the ID word (record step); the Y fraction
entry by three LD64B (rows 0,1; rows 2,3; T(x_Y), a zero word and y_Y); TP1 = sum_k w_k (rY_k + rX_k);
U = TP1 - T(x_Y) - T(x_X); M1 = M(y_Y,w); TP2 = sum_k M1_k r2_k with the rows of 2y_X by two LD64B;
L+ = U - TP2 and L- = U + TP2 (mod 2^192), each emitted as one record.
Y entry (256 bytes): +0 rows 0..3 of y_Y^2 R^-1, +128 T(x_Y), +152 zero, +160 y_Y (Montgomery form).
X area (512 bytes): +0 rows 0..3 of 2 y_X R^-1, +128 rows 0..3 of y_X^2 R^-1, +256 T(x_X), +280 zero.
"""
import re
FB_ARG = 'x10'


def tp(m, acc, base, off, rows, tmp, base2=None, off2=None, rows2=None):
    A3, A2, A1 = acc
    t0, t1, t2 = tmp
    s0, s1, s2, s3 = rows
    out = []
    for k, mk in enumerate(m):
        o = off + 32 * k
        out += [f' ldp {s0}, {s1}, [{base}, #{o}]', f' ldp {s2}, {s3}, [{base}, #{o + 16}]']
        if base2:
            u0, u1, u2, u3 = rows2
            o2 = off2 + 32 * k
            out += [f' ldp {u0}, {u1}, [{base2}, #{o2}]', f' ldp {u2}, {u3}, [{base2}, #{o2 + 16}]',
                    f' adds {s0}, {s0}, {u0}', f' adcs {s1}, {s1}, {u1}', f' adcs {s2}, {s2}, {u2}', f' adc {s3}, {s3}, {u3}']
        if k == 0:
            out += [f' umulh {A3}, {mk}, {s0}', f' umulh {A2}, {mk}, {s1}', f' umulh {A1}, {mk}, {s2}']
        else:
            out += [f' umulh {t0}, {mk}, {s0}', f' umulh {t1}, {mk}, {s1}', f' umulh {t2}, {mk}, {s2}',
                    f' adds {A3}, {A3}, {t0}', f' adcs {A2}, {A2}, {t1}', f' adc {A1}, {A1}, {t2}']
        out += [f' mul {t0}, {mk}, {s1}', f' mul {t1}, {mk}, {s2}', f' mul {t2}, {mk}, {s3}',
                f' adds {A3}, {A3}, {t0}', f' adcs {A2}, {A2}, {t1}', f' adc {A1}, {A1}, {t2}']
    return out


def _sections(variant):
    lines = variant.splitlines()
    def sec(a, b):
        return lines[lines.index(a + ':') + 1:lines.index(b + ':')]
    sq = sec('k_mul_store_5_y0_c0', 'fs_negative_6_y0_c0')
    sq = sq[:[i for i, l in enumerate(sq) if l.strip().startswith('b.cs fs_high')][0]]
    mu = sec('k_add_store_4_y0_c0', 'k_mul_store_5_y0_c0')
    mu = mu[:[i for i, l in enumerate(mu) if l.strip().startswith('b.cc')][0]]
    ren = lambda ls, mp: [re.sub(r'\b(x\d+)\b', lambda m: mp.get(m[1], m[1]), l) for l in ls]
    assert not any(re.search(r'\bx2[4-7]\b', l) for l in sq) and not any(re.search(r'\bx[5-8]\b', l) for l in mu)
    sqr = ren(sq, {'x9': 'x24', 'x10': 'x25', 'x3': 'x26', 'x4': 'x27'})  # i in x24..x27 -> w in x7,x8,x5,x6
    mul = ren(mu, {'x24': 'x7', 'x25': 'x8', 'x26': 'x5', 'x27': 'x6'})  # y_Y (x11..x14) * w -> x9,x10,x3,x4
    return sqr, mul


TAIL4 = lambda r, l: [f' b.cc {l}', f' adds {r[0]}, {r[0]}, #1', f' sbcs {r[1]}, {r[1]}, x0',
                      f' sbcs {r[2]}, {r[2]}, xzr', f' sbc {r[3]}, {r[3]}, x2', f'{l}:']
FLUSH = [' mov x16,x30', ' bl flush_publish_b', ' mov x30,x16']


def tp_rows(m, acc, rows, tmp, first):
    A3, A2, A1 = acc
    t0, t1, t2 = tmp
    s0, s1, s2, s3 = rows
    if first:
        out = [f' umulh {A3}, {m}, {s0}', f' umulh {A2}, {m}, {s1}', f' umulh {A1}, {m}, {s2}']
    else:
        out = [f' umulh {t0}, {m}, {s0}', f' umulh {t1}, {m}, {s1}', f' umulh {t2}, {m}, {s2}',
               f' adds {A3}, {A3}, {t0}', f' adcs {A2}, {A2}, {t1}', f' adc {A1}, {A1}, {t2}']
    return out + [f' mul {t0}, {m}, {s1}', f' mul {t1}, {m}, {s2}', f' mul {t2}, {m}, {s3}',
                  f' adds {A3}, {A3}, {t0}', f' adcs {A2}, {A2}, {t1}', f' adc {A1}, {A1}, {t2}']


def add_xrow(rows, off):
    s0, s1, s2, s3 = rows
    return [f' ldp x3, x4, [x22, #{off}]', f' adds {s0}, {s0}, x3', f' adcs {s1}, {s1}, x4',
            f' ldp x3, x4, [x22, #{off + 16}]', f' adcs {s2}, {s2}, x3', f' adc {s3}, {s3}, x4']


LO, HI = ('x10', 'x11', 'x12', 'x13'), ('x14', 'x15', 'x16', 'x17')  # one LD64B line: two rows
MUL_REN = {'x11': 'x14', 'x12': 'x15', 'x13': 'x16', 'x14': 'x17', 'x15': 'x11', 'x16': 'x12', 'x10': 'x27'}
XR2, XRX, XTX = 0, 128, 256


def tail(sqr, mul, s):
    mul = [re.sub(r'\b(x\d+)\b', lambda m: MUL_REN.get(m[1], m[1]), l) for l in mul]  # y_Y in x14..x17, M1 -> x9,x27,x3,x4
    L = list(sqr) + TAIL4(['x7', 'x8', 'x5', 'x6'], f'fw_ok{s}')
    L += [' ldr x1, [x20, #-72]!', ' add x18, x30, x1, lsr #28']
    acc, tmp = ('x24', 'x25', 'x26'), ('x27', 'x28', 'x29')
    L += [' ld64b x10, [x18]']
    L += add_xrow(LO, XRX) + tp_rows('x7', acc, LO, tmp, True)
    L += add_xrow(HI, XRX + 32) + tp_rows('x8', acc, HI, tmp, False)
    L += [' add x18, x18, #64', ' ld64b x10, [x18]']
    L += add_xrow(LO, XRX + 64) + tp_rows('x5', acc, LO, tmp, False)
    L += add_xrow(HI, XRX + 96) + tp_rows('x6', acc, HI, tmp, False)
    L += [' add x18, x18, #64', ' ld64b x10, [x18]',
          ' subs x24, x24, x10', ' sbcs x25, x25, x11', ' sbc x26, x26, x12',
          f' ldp x27, x28, [x22, #{XTX}]', f' ldr x29, [x22, #{XTX + 16}]',
          ' subs x24, x24, x27', ' sbcs x25, x25, x28', ' sbc x26, x26, x29']
    L += list(mul) + TAIL4(['x9', 'x27', 'x3', 'x4'], f'fm_ok{s}')
    acc2, tmp2 = ('x5', 'x6', 'x8'), ('x7', 'x28', 'x29')
    L += [' ld64b x10, [x22]'] + tp_rows('x9', acc2, LO, tmp2, True) + tp_rows('x27', acc2, HI, tmp2, False)
    L += [' add x18, x22, #64', ' ld64b x10, [x18]'] + tp_rows('x3', acc2, LO, tmp2, False) + tp_rows('x4', acc2, HI, tmp2, False)
    L += [' subs x9, x24, x5', ' sbcs x7, x25, x6', ' sbc x10, x26, x8',
          ' ldp x13,x14,[sp]', ' bfi x13,x7,#3,#35', ' mov x17,#32', ' ldadd x17,x15,[x13]',
          ' orr x11, x1, x9, lsr #29', ' extr x12, x9, x7, #35', ' stp x10,x12,[x15]', ' stp x14,x11,[x15,#16]',
          f' tbz x15,#12,fskip1{s}'] + FLUSH + [f'fskip1{s}:']
    L += [' adds x9, x24, x5', ' adcs x7, x25, x6', ' adc x10, x26, x8',
          ' bfi x13,x7,#3,#35', ' ldadd x17,x15,[x13]', ' orr x11, x1, #0x800000000', ' orr x11, x11, x9, lsr #29',
          ' extr x12, x9, x7, #35', ' stp x10,x12,[x15]', ' stp x14,x11,[x15,#16]',
          f' tbz x15,#12,fskip2{s}'] + FLUSH + [f'fskip2{s}:']
    return L


def derive(code):
    parts = re.split(r'(?=\.globl pair_bwd_y\d_c\d)', code)
    head, v = parts[0], parts[1]
    assert v.startswith('.globl pair_bwd_y0_c0') and len(parts) == 5
    sqr, mul = _sections(v)
    lines = v.split('\n')
    a = lines.index('k_mul_store_1_y0_c0:')
    assert lines[a + 1].strip().startswith('stp x9,x10,[x23]') and lines[a + 2].strip().startswith('stp x11,x12,[x23')
    b = lines.index('fskip2_y0_c0:')
    body = lines[:a + 3] + tail(sqr, mul, '_f') + lines[b + 1:]
    z = '\n'.join(body)
    z = z.replace('.globl pair_bwd_y0_c0\npair_bwd_y0_c0:', '.globl pair_bwd_f\npair_bwd_f:')
    pro = ' cbz x19,pb_done_y0_c0\n'
    frame = ' sub sp,sp,#48\n stp x5,x6,[sp]\n stp x7,x8,[sp,#16]\n str x9,[sp,#32]\n'
    assert z.count(pro) == 1 and z.count(frame) == 1
    z = z.replace(frame, '').replace(pro, frame + f' mov x30,{FB_ARG}\n' + pro)  # frame before the empty test
    for lab in re.findall(r'^([A-Za-z_]\w*_y0_c0):', z, re.M):
        z = re.sub(r'\b' + lab + r'\b', lab[:-6] + '_f', z)
    assert '_y0_c0' not in z and 'x20,#136' not in z and z.count('bl flush_publish_b') == 2
    for reg in ('x19',):
        t = z[z.index('pb_loop_f:'):z.index('pb_done_f:')]
        assert len(re.findall(r'\b' + reg + r'\b', t)) == 2
    tail_txt = parts[4].split('.globl flush_publish_b', 1)
    return head, z, '.globl flush_publish_b' + tail_txt[1]


P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
RR = 1 << 256; R = RR % P; RI = pow(R, -1, P); MASK = (1 << 64) - 1; ni = -pow(P, -1, RR) % RR


def mont(a, b):
    t = (a * b + (a * b * ni % RR) * P) >> 256
    return (t - P, 1) if t >= RR else (t, 0)


def ls64_regions(FB, XA):
    return ([(FB, FB + (256 << 28)), (XA, XA + 512)], [])


def frow(c):
    return [((c * pow(2, 64 * k, P) % P) << 256) // P for k in range(4)]


def t192(v):
    return (v << 192) // P


def tp_model(m, rows):
    a = 0
    for k in range(4):
        mk = (m >> (64 * k)) & MASK; r = rows[k]
        w = [(r >> (64 * j)) & MASK for j in range(4)]
        hi = [(mk * w[j]) >> 64 for j in range(3)]; lo = [(mk * w[j]) & MASK for j in (1, 2, 3)]
        a = (a + (hi[0] | hi[1] << 64 | hi[2] << 128)) % (1 << 192)
        a = (a + (lo[0] | lo[1] << 64 | lo[2] << 128)) % (1 << 192)
    return a


def tables(xX, yX, xY, yY):
    rY = frow(yY * yY * RI % P); rX = frow(yX * yX * RI % P); r2 = frow(2 * yX * RI % P)
    sig = [(rY[k] + rX[k]) % RR for k in range(4)]
    return dict(tY=t192(xY), tX=t192(xX), rY=rY, rX=rX, r2=r2, sig=sig)


def model(i_rep, xX, yX, xY, yY, yY_m=None):
    tb = tables(xX, yX, xY, yY)
    w, cw = mont(i_rep, i_rep)
    a = tp_model(w, tb['sig'])
    u = (a - tb['tY'] - tb['tX']) % (1 << 192)
    m1, cm = mont(yY * R % P if yY_m is None else yY_m, w)
    bb = tp_model(m1, tb['r2'])
    return (u - bb) % (1 << 192), (u + bb) % (1 << 192), (cw, cm)


def record(L, sign, idw, w2):
    L1, L2, L3 = L >> 128, (L >> 64) & MASK, L & MASK
    k3 = ((L2 >> 35) | (L3 << 29)) & MASK
    return (L1, k3, w2, idw | (sign << 35) | (L3 >> 29)), L2 & ((1 << 35) - 1)


def put4(mem, a, v):
    for k in range(4): mem[a + 8 * k] = (v >> (64 * k)) & MASK


def install(mem, FB, XA, xX, yX, items):
    tb = tables(xX, yX, 1, 1)
    for k in range(4):
        put4(mem, XA + 32 * k, tb['r2'][k]); put4(mem, XA + 128 + 32 * k, tb['rX'][k])
    for j in range(3): mem[XA + 256 + 8 * j] = (tb['tX'] >> (64 * j)) & MASK
    mem[XA + 280] = 0
    for xm, ym, l in items:
        e = FB + 256 * l; ty = t192(xm * RI % P); ry = frow((ym * RI % P) ** 2 * RI % P)
        for k in range(4): put4(mem, e + 32 * k, ry[k])
        for j in range(3): mem[e + 128 + 8 * j] = (ty >> (64 * j)) & MASK
        mem[e + 152] = 0
        put4(mem, e + 160, ym)


def expected(B, A, d, order, xX, yX, w2):
    n = len(order); recs = []; tails = [0] * n; errs = []
    for j in range(n - 1, -1, -1):
        xm, ym, l = order[j]
        i, c0 = mont(B, A[j]); B, c1 = mont(B, d[j])
        xY, yY = xm * RI % P, ym * RI % P
        Lp, Lm, (cw, cm) = model(i, xX, yX, xY, yY, yY_m=ym)
        tails[j] = c0 + c1 + cw + cm
        dx = (xY - xX) % P
        if dx:
            ii = pow(dx, -1, P)
            for L, s in ((Lp, -1), (Lm, 1)):
                lam = (yY + s * yX) * ii % P if s == 1 else (yY - yX) * ii % P
                x3 = (lam * lam - xX - xY) % P; e = (L - t192(x3)) % (1 << 192)
                errs.append(e - (1 << 192) if e >= 1 << 191 else e)
        for L, sg in ((Lp, 0), (Lm, 1)):
            recs.append(record(L, sg, l << 36, w2))
    return recs, tails, errs


def tail_free(n):
    return 529 * n + 21 - 14 * (n // 8)
