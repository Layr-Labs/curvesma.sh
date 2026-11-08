"""Flipped-grid leaf kernels (proof.md section 3): generator of the listings and a bit-exact model.
Streams: the k_mul and k_sqr bodies of p256_fields.s with operands in registers (register renaming only),
the carry-only tail after B.CC using the reserved constants x0 = 0xffffffff and x2 = 0xffffffff00000001.
Table lookups: a 256-bit integer v = sum_k v_k 2^(37k) (seven windows: six of 37 bits, the last of 34) selects
row v_k, slot k < 7 of a 32 TiB table (2^37 rows; row t at +256 t, slot k at +32 k, holding (c t 2^(37k) + [k = 0] add0)
mod p; slot 7 is zero and never read); the window is inserted by BFI into address bits 8..44 and slot pairs (2i, 2i+1)
are read by one LD64B of line i.
"""
import re
P = 115792089210356248762697446949407573530086143415290314195533631308867097853951
M64 = (1 << 64) - 1
RR = 1 << 256
R = RR % P
RI = pow(R, -1, P)
DELTA = RR - P
NI = -pow(P, -1, RR) % RR
def _body(src, name):
    return [l.strip() for l in src.split(name + ':\n')[1].split(' ret')[0].strip().splitlines()]
def _ren(lines, mp):
    return [re.sub('\\b(x\\d+)\\b', lambda m: mp.get(m[1], m[1]), l) for l in lines]
def tail4(r, label):
    """Carry-only REDC tail: subtract p when the raw result reached 2^256 (taken iff C = 1)."""
    return [f'b.cc {label}', f'adds {r[0]}, {r[0]}, #1', f'sbcs {r[1]}, {r[1]}, x0', f'sbcs {r[2]}, {r[2]}, xzr', f'sbc {r[3]}, {r[3]}, x2', f'{label}:']
def mul_lines(src, a, b, out, tmp, label):
    """out = REDC(a*b) < 2^256 with a, b in registers (a is consumed, b is preserved): 95 + B.CC (+4 taken tail)."""
    keep = []
    for s in _body(src, 'k_mul'):
        if s.startswith('ldp ') or s.startswith('mov x17'):
            continue
        if s.startswith('b.cc'):
            break
        keep.append(s)
    assert len(keep) == 95
    mp = {'x3': a[0], 'x4': a[1], 'x5': a[2], 'x6': a[3], 'x7': b[0], 'x8': b[1], 'x13': b[2], 'x14': b[3], 'x9': out[0], 'x10': out[1], 'x11': out[2], 'x12': out[3], 'x15': tmp[0], 'x16': tmp[1], 'x17': 'x2'}
    regs = list(a) + list(b) + list(out) + list(tmp)
    assert len(set(regs)) == 14 and 'x0' not in regs and ('x2' not in regs)
    return _ren(keep, mp) + tail4(out, label)
def sqr_lines(src, a, out, tmp, label):
    """out = REDC(a*a) < 2^256 with a in registers (consumed): 75 + B.CC (+4 taken tail). tmp = 8 registers."""
    keep = []
    for s in _body(src, 'k_sqr'):
        if s.startswith('ldp '):
            continue
        if s.startswith('adc x11, xzr, xzr'):
            break
        keep.append(s)
    assert len(keep) == 75
    orig_t = ['x1', 'x11', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17']
    mp = {'x3': a[0], 'x4': a[1], 'x5': a[2], 'x6': a[3], 'x7': out[0], 'x8': out[1], 'x9': out[2], 'x10': out[3]}
    mp.update(zip(orig_t, tmp))
    regs = list(a) + list(out) + list(tmp)
    assert len(set(regs)) == 16 and 'x0' not in regs and ('x2' not in regs)
    return _ren(keep, mp) + tail4(out, label)
WB, NW = 37, 7
def window_src(limbs, k, t):
    """Instructions leaving window k (bits 37k .. 37k+36 of the 256-bit value in limbs) in the low bits of the returned
    register (bits above the window are zero only for k = 6, BFI keeps 37 bits): (lines, register)."""
    lo, sh = divmod(WB * k, 64)
    if sh == 0:
        return [], limbs[lo]
    if sh + WB <= 64 or lo == 3:
        return [f'lsr {t}, {limbs[lo]}, #{sh}'], t
    return [f'extr {t}, {limbs[lo + 1]}, {limbs[lo]}, #{sh}'], t
def lookup_lines(b, limbs, rt, acc, c, t, addend=None, tsel=None, bsrc=None):
    """acc + c 2^256 = sum_k slot_k(row v_k) (+ addend): seven windows of the 256-bit value in `limbs`.
    b: address register; when tsel is given, b = bsrc | tsel first (table select, line bits 0).
    Otherwise b must already have line bits 6,7 = 0. rt: first LD64B register (x<rt>..x<rt+7>).
    t: temporary for the odd windows (may lie in the LD64B area). Line bits end at 3 (192).
    """
    L = []
    R8 = [f'x{rt + i}' for i in range(8)]
    assert rt % 2 == 0 and rt <= 22 and (b not in R8) and (c not in R8) and (not set(acc) & set(R8))
    assert not set(limbs) & set(R8) and t not in (b, c) and (t not in acc) and (t not in limbs)
    if tsel is not None:
        L.append(f'orr {b}, {bsrc}, #{tsel}')
    for k in range(NW):
        if k in (2, 6):
            L.append(f'orr {b}, {b}, #64')
        if k == 4:
            L.append(f'eor {b}, {b}, #192')
        pre, reg = window_src(limbs, k, t)
        L += pre + [f'bfi {b}, {reg}, #8, #{WB}']
        half = R8[4 * (k % 2):4 * (k % 2) + 4]
        if k == 0 and addend is None:
            L += [f'ldp {acc[0]}, {acc[1]}, [{b}]', f'ldp {acc[2]}, {acc[3]}, [{b}, #16]']
            continue
        L.append(f'ld64b x{rt}, [{b}]')
        src = addend if k == 0 else acc
        L += [f'adds {acc[0]}, {src[0]}, {half[0]}', f'adcs {acc[1]}, {src[1]}, {half[1]}', f'adcs {acc[2]}, {src[2]}, {half[2]}', f'adcs {acc[3]}, {src[3]}, {half[3]}']
        first_carry = k == 0 or (k == 1 and addend is None)
        L.append(f'cset {c}, hs' if first_carry else f'adc {c}, {c}, xzr')
    return L
def fold_lines(acc, c, t, label):
    """Exact fold for c<=8: 8 instructions without carry, 12 with carry."""
    t0, t1, t2, t3 = t
    return [f'lsl {t0}, {c}, #5', f'ldp {t1}, {t2}, [{t0}]', f'ldp {t3}, {t0}, [{t0}, #16]', f'adds {acc[0]}, {acc[0]}, {t1}', f'adcs {acc[1]}, {acc[1]}, {t2}', f'adcs {acc[2]}, {acc[2]}, {t3}', f'adcs {acc[3]}, {acc[3]}, {t0}'] + tail4(acc, label)
def fold_table():
    return [c * DELTA for c in range(9)]
def tp_rows(m, acc, rows, tmp, first):
    """One truncated-product row step: acc (A_lo, A_mid, A_hi) += floor(m*row/2^64) mod 2^192."""
    A3, A2, A1 = acc
    t0, t1, t2 = tmp
    s0, s1, s2, s3 = rows
    if first:
        out = [f'umulh {A3}, {m}, {s0}', f'umulh {A2}, {m}, {s1}', f'umulh {A1}, {m}, {s2}']
    else:
        out = [f'umulh {t0}, {m}, {s0}', f'umulh {t1}, {m}, {s1}', f'umulh {t2}, {m}, {s2}', f'adds {A3}, {A3}, {t0}', f'adcs {A2}, {A2}, {t1}', f'adc {A1}, {A1}, {t2}']
    return out + [f'mul {t0}, {m}, {s1}', f'mul {t1}, {m}, {s2}', f'mul {t2}, {m}, {s3}', f'adds {A3}, {A3}, {t0}', f'adcs {A2}, {A2}, {t1}', f'adc {A1}, {A1}, {t2}']
def redc(a, b):
    t = a * b + a * b * NI % RR * P >> 256
    return (t - P, 1) if t >= RR else (t, 0)
def table_entry(cst, add0, t, k):
    return (cst * t * pow(2, WB * k, P) + (add0 if k == 0 else 0)) % P if k < NW else 0
def lookup_model(cst, add0, v, addend=0):
    s = addend + sum((table_entry(cst, add0, v >> WB * k & (1 << WB) - 1, k) for k in range(NW)))
    return (s & RR - 1, s >> 256)
def fold_model(o, c):
    s = o + c * DELTA
    if s >= RR:
        s = s - RR + DELTA
    assert s < RR
    return s
YE, XE = (1024, 256)
STEP_F = (1 << 42) - 192
FLUSH_SUB = '.globl flush_publish_b\nflush_publish_b:\n and x7,x7,#0xfffffffff\n ldr x9,[sp,#16]\n ldr x10,[sp,#24]\n ldr x11,[sp,#32]\n mov x12,#4064\n sub x15,x15,x12\n str x15,[x13]\n add x9,x9,x7,lsl #3\n mov x12,#4096\n ldadd x12,x3,[x9]\n cmp x3,x10\n b.hs fpb_skip\n mov x12,#128\n stp x7,x3,[x11]\n stp x15,x12,[x11,#16]\nfpb_skip:\n ret\n'
def _uniq(lines, suf):
    labs = [l[:-1] for l in lines if l.endswith(':')]
    out = []
    for l in lines:
        for lab in labs:
            l = re.sub('\\b' + lab + '\\b', lab + suf, l)
        out.append(l)
    return out
def fwd_quad(src, s):
    """X-pairs A=2P, B=2P+1 share x_Y's seven windows: one LD64B of line k%4, row v_k of tau0 table T(k//4) of pair pair P
    gives both window-k slots [A|B]; e_A, e_B = fold(S_Y + slots); store [e_A, A_prev], A = M(A_prev, e_A), store [e_B, A],
    A = M(A, e_B) (alternating register blocks: no reload)."""
    A, B, lo, hi = (['x16', 'x17', 'x18', 'x19'], ['x4', 'x5', 'x6', 'x7'], ['x8', 'x9', 'x10', 'x11'], ['x12', 'x13', 'x14', 'x15'])
    L = ['ldp x16, x17, [x1, #32]', 'ldp x18, x19, [x1, #48]']
    for k in range(NW):
        o = ('x24', 'x25')[k % 2]
        if k % 2 == 0:
            L.append(f'ldp x24, x25, [sp, #{48 + 8 * k}]' if k + 1 < NW else f'ldr x24, [sp, #{48 + 8 * k}]')
        L += [f'add {o}, x28, {o}', f'ld64b x8, [{o}]']
        if k == 0:
            L += [f'adds {B[0]}, {A[0]}, {hi[0]}'] + [f'adcs {B[i]}, {A[i]}, {hi[i]}' for i in (1, 2, 3)] + ['cset x27, hs']
            L += [f'adds {A[0]}, {A[0]}, {lo[0]}'] + [f'adcs {A[i]}, {A[i]}, {lo[i]}' for i in (1, 2, 3)] + ['cset x26, hs']
        else:
            for acc, h, c in ((A, lo, 'x26'), (B, hi, 'x27')):
                L += [f'adds {acc[0]}, {acc[0]}, {h[0]}'] + [f'adcs {acc[i]}, {acc[i]}, {h[i]}' for i in (1, 2, 3)]
                L.append(f'adc {c}, {c}, xzr')
    L += fold_lines(A, 'x26', lo, 'fefa' + s) + fold_lines(B, 'x27', hi, 'fefb' + s)
    L += ['st64b x16, [x29]', 'add x29, x29, #64']
    L += mul_lines(src, ['x20', 'x21', 'x22', 'x23'], A, lo, ['x12', 'x13'], 'ff_ma' + s)
    L += ['st64b x4, [x29]', 'add x29, x29, #64']
    L += mul_lines(src, lo, B, ['x20', 'x21', 'x22', 'x23'], ['x12', 'x13'], 'ff_mb' + s)
    L += ['add x28, x28, x30']
    return L
WMASK = ((1 << WB) - 1) << 8
def fwd_offset_lines(o, limbs, k):
    """o = 256 v_k + 64 (k%4) + 2^45 (k//4): window k of x_Y at address bits 8..44 (AND with the 37-bit mask)."""
    lo, sh = divmod(WB * k, 64)
    sh -= 8
    if sh < 0:
        L = [f'lsl {o}, {limbs[lo]}, #{-sh}']
    elif sh + 8 + WB <= 64 or lo == 3:
        L = [f'lsr {o}, {limbs[lo]}, #{sh}']
    else:
        L = [f'extr {o}, {limbs[lo + 1]}, {limbs[lo]}, #{sh}']
    L.append(f'and {o}, {o}, #{WMASK:#x}')
    return L + [f'orr {o}, {o}, #{1 << 45:#x}'] * (k >= 4) + [f'orr {o}, {o}, #{64 * (k % 4)}'] * (k % 4 > 0)
def fwd_listing(src, U=2):
    """flip_fwd(x0 = Y count >= 1, x1 = first Y entry, x2 = batch entry 0, x3 = TB, x4 = 64 MP with MP % (2U) == 0).
    Per Y: offsets O_k = 2^45 (k//4) + 256 v_k + 64 (k%4) of the seven windows v_k of x_Y at sp+48..sp+103.
    Stores entry j = [e_j, A_(j-1)] (64 bytes) for every (Y, pair) in order and the final prefix after the last entry."""
    L = ['.globl flip_fwd', 'flip_fwd:', ' sub sp, sp, #112', ' str x30, [sp]', ' add x5, x1, x0, lsl #10', ' stp x1, x5, [sp, #8]', ' stp x3, x4, [sp, #24]', ' mov x29, x2', ' mov w0, #-1', ' mov x2, #0xffffffff00000001', ' mov x30, #0x1000000000000', ' mov x5, #32', ' ldp x20, x21, [x5]', ' ldp x22, x23, [x5, #16]', 'ff_y:', ' ldr x1, [sp, #8]', ' ldp x12, x13, [x1]', ' ldp x14, x15, [x1, #16]']
    for k in range(NW):
        L += [' ' + l for l in fwd_offset_lines(f'x{4 + k}', ['x12', 'x13', 'x14', 'x15'], k)]
    L += [' stp x4, x5, [sp, #48]', ' stp x6, x7, [sp, #64]', ' stp x8, x9, [sp, #80]', ' str x10, [sp, #96]', ' add x5, x1, #1024', ' str x5, [sp, #8]', ' ldp x28, x3, [sp, #24]', ' add x3, x29, x3', 'ff_pair:']
    for u in range(U):
        L += [l if l.endswith(':') else ' ' + l for l in _uniq(fwd_quad(src, ''), f'_u{u}')]
    L += [' cmp x29, x3', ' b.ne ff_pair', ' ldp x5, x6, [sp, #8]', ' cmp x5, x6', ' b.ne ff_y', ' stp x20, x21, [x29]', ' stp x22, x23, [x29, #16]', ' ldr x30, [sp]', ' add sp, sp, #112', ' ret']
    return '\n'.join(L) + '\n'
def label_lines(src, xo, tsel, s):
    """Records of one incidence from i (x24..x27): w = i^2, TP1, U, N1 = n_X(w) by lookup, TP2 (5 rows), L+-.
    Both signs and both members of an X-pair carry the same compact source metadata; resolution recovers the omitted choices.
    Both signs emit with ORR x11, x1, x9, lsr #28. The old plus-path MOV is not issued."""
    L = sqr_lines(src, ['x24', 'x25', 'x26', 'x27'], ['x3', 'x4', 'x5', 'x6'], ['x7', 'x8', 'x9', 'x10', 'x11', 'x12', 'x13', 'x14'], 'fw' + s)
    acc, tmp = (('x24', 'x25', 'x26'), ('x27', 'x28', 'x29'))
    LO, HI = (('x8', 'x9', 'x10', 'x11'), ('x12', 'x13', 'x14', 'x15'))
    def xrow(rows, off):
        s0, s1, s2, s3 = rows
        return [f'ldp x16, x17, [x22, #{off}]', f'adds {s0}, {s0}, x16', f'adcs {s1}, {s1}, x17', f'ldp x16, x17, [x22, #{off + 16}]', f'adcs {s2}, {s2}, x16', f'adc {s3}, {s3}, x17']
    L += ['add x18, x20, #64', 'ld64b x8, [x18]']
    L += xrow(LO, xo) + tp_rows('x3', acc, LO, tmp, True) + xrow(HI, xo + 32) + tp_rows('x4', acc, HI, tmp, False)
    L += ['add x18, x18, #64', 'ld64b x8, [x18]']
    L += xrow(LO, xo + 64) + tp_rows('x5', acc, LO, tmp, False) + xrow(HI, xo + 96) + tp_rows('x6', acc, HI, tmp, False)
    L += ['add x18, x18, #64', 'ld64b x8, [x18]', 'subs x24, x24, x8', 'sbcs x25, x25, x9', 'sbc x26, x26, x10', f'ldp x16, x17, [x22, #{xo + 128}]', f'ldp x27, x1, [x22, #{xo + 144}]', 'subs x24, x24, x16', 'sbcs x25, x25, x17', 'sbc x26, x26, x27', 'orr x1, x1, x11']
    L += lookup_lines('x16', ['x3', 'x4', 'x5', 'x6'], 8, ['x17', 'x27', 'x28', 'x29'], 'x18', 'x8', tsel=tsel, bsrc='x23')
    acc2, tmp2 = (('x4', 'x5', 'x6'), ('x3', 'x7', 'x16'))
    L += ['add x18, x20, x18, lsl #5', 'ldp x4, x5, [x18, #384]', 'ldr x6, [x18, #400]']
    L += ['add x18, x20, #256', 'ld64b x8, [x18]'] + tp_rows('x17', acc2, LO, tmp2, False) + tp_rows('x27', acc2, HI, tmp2, False)
    L += ['add x18, x18, #64', 'ld64b x8, [x18]'] + tp_rows('x28', acc2, LO, tmp2, False) + tp_rows('x29', acc2, HI, tmp2, False)
    FL = ['mov x16,x30', 'bl flush_publish_b', 'mov x30,x16']
    L += ['subs x9, x24, x4', 'sbcs x7, x25, x5', 'sbc x10, x26, x6', 'ldp x13,x14,[sp]', 'bfi x13,x7,#3,#36', 'mov x17,#32', 'ldadd x17,x15,[x13]', 'orr x11, x1, x9, lsr #28', 'extr x12, x9, x7, #36', 'stp x10,x12,[x15]', 'stp x14,x11,[x15,#16]', f'tbz x15,#12,fs1{s}'] + FL + [f'fs1{s}:']
    # Plus sign uses the same ORR as the minus sign. x1 is not written by the flush
    # body, so the previous MOV x11,x1 is redundant on every path.
    L += ['adds x9, x24, x4', 'adcs x7, x25, x5', 'adc x10, x26, x6', 'bfi x13,x7,#3,#36', 'ldadd x17,x15,[x13]', 'orr x11, x1, x9, lsr #28', 'extr x12, x9, x7, #36', 'stp x10,x12,[x15]', 'stp x14,x11,[x15,#16]', f'tbz x15,#12,fs2{s}'] + FL + [f'fs2{s}:']
    return L
def bwd_pair(src, s, ev):
    """One X-pair of the backward pass (pairs in reverse): inverse split by partial fractions, then both incidences.
    ev: even J, whose tau1..3 sit one 2^45 slot higher in the pair-pair region (tau0 is interleaved)."""
    q = (1 << 45) * ev
    L = ['ld64b x8, [x21]', 'sub x21, x21, #64', 'ldp x24, x25, [sp, #48]', 'ldp x26, x27, [sp, #64]']
    L += mul_lines(src, ['x12', 'x13', 'x14', 'x15'], ['x24', 'x25', 'x26', 'x27'], ['x4', 'x5', 'x6', 'x7'], ['x16', 'x17'], 'fbi' + s)
    L += mul_lines(src, ['x8', 'x9', 'x10', 'x11'], ['x24', 'x25', 'x26', 'x27'], ['x12', 'x13', 'x14', 'x15'], ['x16', 'x17'], 'fbb' + s)
    L += ['stp x12, x13, [sp, #48]', 'stp x14, x15, [sp, #64]', 'ldp x8, x9, [x20]', 'ldp x10, x11, [x20, #16]', 'ldp x12, x13, [x22, #416]', 'ldp x14, x15, [x22, #432]', 'subs x8, x8, x12', 'sbcs x9, x9, x13', 'sbcs x10, x10, x14', 'sbcs x11, x11, x15', f'b.cs fbd2done{s}', f'fbd2borrow{s}:', 'subs x8, x8, #1', 'adcs x9, x9, x0', 'adcs x10, x10, xzr', 'adc x11, x11, x2', f'fbd2done{s}:']
    L += mul_lines(src, ['x8', 'x9', 'x10', 'x11'], ['x4', 'x5', 'x6', 'x7'], ['x24', 'x25', 'x26', 'x27'], ['x16', 'x17'], 'fbd' + s)
    L += lookup_lines('x1', ['x4', 'x5', 'x6', 'x7'], 8, ['x16', 'x17', 'x18', 'x3'], 'x28', 'x8', addend=['x24', 'x25', 'x26', 'x27'], tsel=(1 << 45) + q, bsrc='x23')
    L += fold_lines(['x16', 'x17', 'x18', 'x3'], 'x28', ['x8', 'x9', 'x10', 'x11'], 'fbfold' + s)
    L += ['stp x16, x17, [sp, #80]', 'stp x18, x3, [sp, #96]']
    L += label_lines(src, 0, (2 << 45) + q, 'a' + s)
    L += ['ldp x24, x25, [sp, #80]', 'ldp x26, x27, [sp, #96]']
    L += label_lines(src, 256, (3 << 45) + q, 'b' + s)
    L += ['sub x22, x22, #512', 'sub x23, x23, x30']
    return L
def bwd_listing(src, U=4):
    """flip_bwd(x0 = Y count, x1 = last Y entry, x2 = last batch entry, x3 = TB + 2^47 (MP-1), x4 = XA + 512 (MP-1),
    x5 = cursor table, x6 unused, x7..x9 = flush words (as 008389f5), x10 = address of B = 1/A_last (Montgomery,
    written back at the end), x11 = 64 MP). Emits two records per incidence."""
    L = ['.globl flip_bwd', 'flip_bwd:', ' sub sp, sp, #176', ' stp x5, x6, [sp]', ' stp x7, x8, [sp, #16]', ' stp x9, x30, [sp, #32]', ' ldp x12, x13, [x10]', ' ldp x14, x15, [x10, #16]', ' stp x12, x13, [sp, #48]', ' stp x14, x15, [sp, #64]', ' sub x5, x1, x0, lsl #10', ' stp x1, x5, [sp, #112]', ' stp x3, x4, [sp, #128]', ' stp x11, x10, [sp, #144]', ' mov x21, x2', ' mov w0, #-1', ' mov x2, #0xffffffff00000001', ' mov x30, #0x800000000000', 'fb_y:', ' ldr x20, [sp, #112]', ' ldr x5, [x20, #224]', ' str x5, [sp, #8]', ' ldp x23, x22, [sp, #128]', ' ldr x19, [sp, #144]', ' sub x19, x21, x19', 'fb_pair:']
    for u in range(U):
        L += [l if l.endswith(':') else ' ' + l for l in _uniq(bwd_pair(src, '', u % 2), f'_u{u}')]
    L += [' cmp x21, x19', ' b.ne fb_pair', ' ldp x20, x5, [sp, #112]', ' sub x20, x20, #1024', ' str x20, [sp, #112]', ' cmp x20, x5', ' b.ne fb_y', ' ldr x10, [sp, #152]', ' ldp x12, x13, [sp, #48]', ' ldp x14, x15, [sp, #64]', ' stp x12, x13, [x10]', ' stp x14, x15, [x10, #16]', ' ldr x30, [sp, #40]', ' add sp, sp, #176', ' ret']
    return '\n'.join(L) + '\n' + FLUSH_SUB
def fsub_lines(a, b, out, t):
    """out = a - b (+ p on borrow): a representative below 2^256 of a - b mod p (a < 2^256, b < p); canonical if a < p."""
    m, t1, t2 = t
    return [f'subs {out[0]}, {a[0]}, {b[0]}', f'sbcs {out[1]}, {a[1]}, {b[1]}', f'sbcs {out[2]}, {a[2]}, {b[2]}', f'sbcs {out[3]}, {a[3]}, {b[3]}', f'sbc {m}, xzr, xzr', f'and {t1}, {m}, x0', f'and {t2}, {m}, x2', f'adds {out[0]}, {out[0]}, {m}', f'adcs {out[1]}, {out[1]}, {t1}', f'adcs {out[2]}, {out[2]}, xzr', f'adc {out[3]}, {out[3]}, {t2}']
def canon_lines(a, t):
    """a < 2^256 -> a mod p (subtract p by CSEL when a >= p)."""
    return [f'adds {t[0]}, {a[0]}, #1', f'sbcs {t[1]}, {a[1]}, x0', f'sbcs {t[2]}, {a[2]}, xzr', f'sbcs {t[3]}, {a[3]}, x2'] + [f'csel {a[i]}, {t[i]}, {a[i]}, hs' for i in range(4)]
MU = {N: (1 << N + 256) // P for N in (192, 448, 512)}
def barrett_lines(c, N, mb, out, t, moff=0):
    """out[4 .. 4 + N/64) = floor(2^N c / p) exactly, for canonical c (4 registers, preserved).
    q_est = floor(c mu / 2^256), mu = floor(2^(N+256)/p) = 2^N + mu' (mu' words at [mb]), is q or q - 1; the exact
    residue r = 2^N c - q_est p is formed mod 2^320 from p = 2^256 - (2^224 - 2^192 - 2^96 + 1) and q_est += [r >= p].
    out: N/64 + 4 registers (product words; no carry leaves the top one since q < 2^N); t: 9 temporaries.
    """
    m = N // 64
    assert len(out) == m + 4 and len(t) == 9
    l0, h0, l1, h1, l2, h2, l3, mu0, mu1 = t
    L = []
    for j in range(m):
        if j % 2 == 0:
            L.append(f'ldp {mu0}, {mu1}, [{mb}, #{moff + 8 * j}]' if j + 1 < m else f'ldr {mu0}, [{mb}, #{moff + 8 * j}]')
        mu = mu0 if j % 2 == 0 else mu1
        if j == 0:
            L += [f'mul {out[0]}, {c[0]}, {mu}', f'umulh {h0}, {c[0]}, {mu}', f'mul {l1}, {c[1]}, {mu}', f'umulh {h1}, {c[1]}, {mu}', f'mul {l2}, {c[2]}, {mu}', f'umulh {h2}, {c[2]}, {mu}', f'mul {l3}, {c[3]}, {mu}', f'umulh {out[4]}, {c[3]}, {mu}', f'adds {out[1]}, {h0}, {l1}', f'adcs {out[2]}, {h1}, {l2}', f'adcs {out[3]}, {h2}, {l3}', f'adc {out[4]}, {out[4]}, xzr']
        else:
            L += [f'mul {l0}, {c[0]}, {mu}', f'umulh {h0}, {c[0]}, {mu}', f'mul {l1}, {c[1]}, {mu}', f'umulh {h1}, {c[1]}, {mu}', f'mul {l2}, {c[2]}, {mu}', f'umulh {h2}, {c[2]}, {mu}', f'mul {l3}, {c[3]}, {mu}', f'adds {h0}, {h0}, {l1}', f'adcs {h1}, {h1}, {l2}', f'adcs {h2}, {h2}, {l3}', f'umulh {l1}, {c[3]}, {mu}', f'adc {l1}, {l1}, xzr', f'adds {out[j]}, {out[j]}, {l0}', f'adcs {out[j + 1]}, {out[j + 1]}, {h0}', f'adcs {out[j + 2]}, {out[j + 2]}, {h1}', f'adcs {out[j + 3]}, {out[j + 3]}, {h2}', f'adc {out[j + 4]}, {l1}, xzr']
    L += [f'adds {out[m]}, {out[m]}, {c[0]}', f'adcs {out[m + 1]}, {out[m + 1]}, {c[1]}', f'adcs {out[m + 2]}, {out[m + 2]}, {c[2]}', f'adcs {out[m + 3]}, {out[m + 3]}, {c[3]}']
    q0, q1, q2, q3, q4 = (out[4:4 + m] + ['xzr', 'xzr'])[:5]
    e1, e2, e3, t3, t4, r1, r2 = (l0, h0, l1, h1, l2, h2, l3)
    L += [f'extr {e1}, {q1}, {q0}, #32', f'extr {e2}, {q2}, {q1}, #32', f'extr {e3}, {q3}, {q2}, #32' if q3 != 'xzr' else f'lsr {e3}, {q2}, #32', f'adds {t3}, {q3}, {q0}, lsl #32', f'adc {t4}, {q4}, {e1}']
    if N == 192:
        L += [f'adds {t3}, {t3}, {c[0]}', f'adc {t4}, {t4}, {c[1]}']
    L += [f'subs {r1}, {q1}, {q0}, lsl #32', f'sbcs {r2}, {q2}, {e1}', f'sbcs {t3}, {t3}, {e2}', f'sbc {t4}, {t4}, {e3}', f'subs {t3}, {t3}, {q0}', f'sbc {t4}, {t4}, {q1}', f'sub {t4}, {t4}, {q0}', f'cmn {q0}, #1', f'sbcs xzr, {r1}, x0', f'sbcs xzr, {r2}, xzr', f'sbcs xzr, {t3}, x2', f'sbcs xzr, {t4}, xzr']
    L += [f'adcs {qq}, {qq}, xzr' for qq in out[4:4 + m]]
    return L
def _ind(lines):
    return [l if l.endswith(':') else ' ' + l for l in lines]
def ychain_f_listing(src):
    """ychain_f(x0 = lanes n >= 1, x1 = chain state of lane 0 (128 bytes per lane: xY, yY, xS, yS canonical Montgomery),
    x2 = scratch entry of lane 0): stores [d_k, A_(k-1)] (d_k = xS - xY mod p) and the final prefix after the last entry."""
    L = ['mov x19, x0', 'mov x20, x1', 'mov x21, x2', 'mov w0, #-1', 'mov x2, #0xffffffff00000001', 'mov x5, #32', 'ldp x24, x25, [x5]', 'ldp x26, x27, [x5, #16]', 'yf_loop:', 'ldp x8, x9, [x20]', 'ldp x10, x11, [x20, #16]', 'ldp x12, x13, [x20, #64]', 'ldp x14, x15, [x20, #80]']
    L += fsub_lines(['x12', 'x13', 'x14', 'x15'], ['x8', 'x9', 'x10', 'x11'], ['x12', 'x13', 'x14', 'x15'], ['x16', 'x17', 'x18'])
    L += ['stp x12, x13, [x21]', 'stp x14, x15, [x21, #16]', 'stp x24, x25, [x21, #32]', 'stp x26, x27, [x21, #48]']
    L += mul_lines(src, ['x24', 'x25', 'x26', 'x27'], ['x12', 'x13', 'x14', 'x15'], ['x3', 'x4', 'x5', 'x6'], ['x16', 'x17'], 'yfm')
    L += ['mov x24, x3', 'mov x25, x4', 'mov x26, x5', 'mov x27, x6', 'add x20, x20, #128', 'add x21, x21, #64', 'subs x19, x19, #1', 'b.ne yf_loop', 'stp x24, x25, [x21]', 'stp x26, x27, [x21, #16]', 'ret']
    return '.globl ychain_f\nychain_f:\n' + '\n'.join(_ind(L)) + '\n'
def ychain_b_listing(src):
    """ychain_b(x0 = lanes n, x1 = chain state of the last lane, x2 = its scratch entry, x3 = address of B = 1/A_last):
    Y <- Y + S for every lane (affine, canonical), B written back (Montgomery one at the end)."""
    Q = lambda a: ['x%d' % (a + i) for i in range(4)]
    T3 = ['x16', 'x17', 'x18']
    L = ['mov x19, x0', 'mov x20, x1', 'mov x21, x2', 'mov x22, x3', 'mov w0, #-1', 'mov x2, #0xffffffff00000001', 'yb_loop:', 'ldp x8, x9, [x21]', 'ldp x10, x11, [x21, #16]', 'ldp x12, x13, [x21, #32]', 'ldp x14, x15, [x21, #48]', 'ldp x24, x25, [x22]', 'ldp x26, x27, [x22, #16]']
    L += mul_lines(src, Q(12), Q(24), Q(4), ['x16', 'x17'], 'ybi')
    L += mul_lines(src, Q(8), Q(24), Q(12), ['x16', 'x17'], 'ybb')
    L += ['stp x12, x13, [x22]', 'stp x14, x15, [x22, #16]', 'ldp x8, x9, [x20, #96]', 'ldp x10, x11, [x20, #112]', 'ldp x12, x13, [x20, #32]', 'ldp x14, x15, [x20, #48]']
    L += fsub_lines(Q(8), Q(12), Q(8), T3)
    L += mul_lines(src, Q(8), Q(4), Q(24), ['x16', 'x17'], 'ybl')
    L += ['mov x8, x24', 'mov x9, x25', 'mov x10, x26', 'mov x11, x27']
    L += sqr_lines(src, Q(8), Q(3), ['x12', 'x13', 'x14', 'x15', 'x16', 'x17', 'x18', 'x1'], 'ybs')
    L += ['ldp x8, x9, [x20]', 'ldp x10, x11, [x20, #16]'] + fsub_lines(Q(3), Q(8), Q(3), T3)
    L += ['ldp x12, x13, [x20, #64]', 'ldp x14, x15, [x20, #80]'] + fsub_lines(Q(3), Q(12), Q(3), T3)
    L += canon_lines(Q(3), Q(12)) + ['stp x3, x4, [x20]', 'stp x5, x6, [x20, #16]']
    L += fsub_lines(Q(8), Q(3), Q(8), T3)
    L += mul_lines(src, Q(8), Q(24), Q(3), ['x16', 'x17'], 'ybm')
    L += ['ldp x12, x13, [x20, #32]', 'ldp x14, x15, [x20, #48]'] + fsub_lines(Q(3), Q(12), Q(3), T3)
    L += canon_lines(Q(3), Q(12)) + ['stp x3, x4, [x20, #32]', 'stp x5, x6, [x20, #48]', 'sub x20, x20, #128', 'sub x21, x21, #64', 'subs x19, x19, #1', 'b.ne yb_loop', 'ret']
    return '.globl ychain_b\nychain_b:\n' + '\n'.join(_ind(L)) + '\n'
KOFF = dict(one=0, two=32, rinv=64, m192=96, m448=128, m512=192)
def kblock():
    v = {0: 1, 32: 2, 64: RI}
    w = {}
    for off, x in v.items():
        for k in range(4):
            w[off + 8 * k] = x >> 64 * k & M64
    for N, off in ((192, 96), (448, 128), (512, 192)):
        mu = MU[N] - (1 << N)
        for k in range(N // 64):
            w[off + 8 * k] = mu >> 64 * k & M64
    return w
def yprep_listing(src):
    """yprep(x0 = lanes n, x1 = chain state of lane 0, x2 = Y entry of lane 0 (1024 bytes apart), x3 = constants block,
    x4 = IDY, x5 = W2 of lane 0): writes every Y entry field of proof.md section 3 (W2 steps by 2^36 per lane)."""
    Q = lambda a: ['x%d' % (a + i) for i in range(4)]
    T8 = ['x12', 'x13', 'x14', 'x15', 'x16', 'x17', 'x18', 'x1']
    T9 = ['x23', 'x24', 'x25', 'x26', 'x27', 'x28', 'x29', 'x30', 'x1']
    L = ['str x30, [sp, #-32]!', 'stp x4, x5, [sp, #8]', 'mov x19, x0', 'mov x20, x1', 'mov x21, x2', 'mov x22, x3', 'mov w0, #-1', 'mov x2, #0xffffffff00000001', 'yp_loop:', 'ldp x8, x9, [x20]', 'ldp x10, x11, [x20, #16]', 'stp x8, x9, [x21]', 'stp x10, x11, [x21, #16]']
    L += sqr_lines(src, Q(8), Q(3), T8, 'yps') + ['stp x3, x4, [x21, #32]', 'stp x5, x6, [x21, #48]']
    L += ['ldp x8, x9, [x20]', 'ldp x10, x11, [x20, #16]', f"ldp x12, x13, [x22, #{KOFF['one']}]", f"ldp x14, x15, [x22, #{KOFF['one'] + 16}]"]
    L += mul_lines(src, Q(8), Q(12), Q(3), ['x16', 'x17'], 'ypx') + canon_lines(Q(3), Q(8))
    L += barrett_lines(Q(3), 192, 'x22', ['x7', 'x8', 'x9', 'x10', 'x11', 'x12', 'x13'], T9, KOFF['m192'])
    L += ['stp x11, x12, [x21, #192]', 'str x13, [x21, #208]']
    L += ['ldp x8, x9, [x20, #32]', 'ldp x10, x11, [x20, #48]'] + sqr_lines(src, Q(8), Q(3), T8, 'ypy')
    L += [f"ldp x12, x13, [x22, #{KOFF['rinv']}]", f"ldp x14, x15, [x22, #{KOFF['rinv'] + 16}]"]
    L += mul_lines(src, Q(3), Q(12), Q(8), ['x16', 'x17'], 'ypr') + canon_lines(Q(8), Q(3))
    G = ['x3', 'x4', 'x5', 'x6', 'x7', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17']
    L += barrett_lines(Q(8), 448, 'x22', G, T9, KOFF['m448'])
    g = G[4:]
    for k in range(4):
        w = g[3 - k:7 - k]
        L += [f'stp {w[0]}, {w[1]}, [x21, #{64 + 32 * k}]', f'stp {w[2]}, {w[3]}, [x21, #{80 + 32 * k}]']
    L += ['ldp x8, x9, [x20, #32]', 'ldp x10, x11, [x20, #48]', f"ldp x12, x13, [x22, #{KOFF['two']}]", f"ldp x14, x15, [x22, #{KOFF['two'] + 16}]"]
    L += mul_lines(src, Q(8), Q(12), Q(3), ['x16', 'x17'], 'ypt') + canon_lines(Q(3), Q(8))
    G = ['x7', 'x8', 'x9', 'x10', 'x11', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17', 'x18']
    L += barrett_lines(Q(3), 512, 'x22', G, T9, KOFF['m512'])
    g = G[4:]
    for k in range(4):
        w = g[4 - k:8 - k]
        L += [f'stp {w[0]}, {w[1]}, [x21, #{256 + 32 * k}]', f'stp {w[2]}, {w[3]}, [x21, #{272 + 32 * k}]']
    r0, r1, r2, r3 = g[0:4]
    a = ['x3', 'x4', 'x5', 'x6', 'x7']
    L += ['add x1, x21, #384', 'stp xzr, xzr, [x1]', 'str xzr, [x1, #16]', f'stp {r1}, {r2}, [x1, #32]', f'str {r3}, [x1, #48]', f'mov {a[0]}, {r0}', f'mov {a[1]}, {r1}', f'mov {a[2]}, {r2}', f'mov {a[3]}, {r3}', f'mov {a[4]}, xzr']
    for c in range(2, 8):
        L += [f'adds {a[0]}, {a[0]}, {r0}', f'adcs {a[1]}, {a[1]}, {r1}', f'adcs {a[2]}, {a[2]}, {r2}', f'adcs {a[3]}, {a[3]}, {r3}', f'adc {a[4]}, {a[4]}, xzr', f'stp {a[1]}, {a[2]}, [x1, #{32 * c}]', f'str {a[3]}, [x1, #{16 + 32 * c}]']
    L += ['ldp x3, x4, [sp, #8]', 'stp x3, x4, [x21, #216]', 'mov x5, #0x1000000000', 'add x4, x4, x5', 'str x4, [sp, #16]', 'add x20, x20, #128', 'add x21, x21, #1024', 'subs x19, x19, #1', 'b.ne yp_loop', 'ldr x30, [sp], #32', 'ret']
    return '.globl yprep\nyprep:\n' + '\n'.join(_ind(L)) + '\n'
def tbld_listing():
    """tbld(x0 = table base, x1 = rows (multiple of 4), x2 = 16 canonical constants (32 bytes each): multipliers m_k at
    +32k, starts s_k at +256+32k): writes slot k of row t = (m_k t + s_k) mod p, canonical, for every row and slot
    (column order). Standard table: m_k = c 2^(32k), s_0 = add0, s_k = 0 (k > 0); interleaved tau0 pair: section 3."""
    L = ['mov x19, x0', 'mov x20, x1', 'mov x21, x2', 'mov w0, #-1', 'mov x2, #0xffffffff00000001', 'mov x22, #0', 'tb_slot:', 'add x23, x19, x22, lsl #5', 'add x24, x21, x22, lsl #5', 'ldp x8, x9, [x24]', 'ldp x10, x11, [x24, #16]', 'ldp x3, x4, [x24, #256]', 'ldp x5, x6, [x24, #272]', 'mov x26, x20', 'tb_loop:']
    for u in range(4):
        L += ['stp x3, x4, [x23]', 'stp x5, x6, [x23, #16]', 'add x23, x23, #256', 'adds x3, x3, x8', 'adcs x4, x4, x9', 'adcs x5, x5, x10', 'adcs x6, x6, x11', 'adc x7, xzr, xzr', 'adds x12, x3, #1', 'sbcs x13, x4, x0', 'sbcs x14, x5, xzr', 'sbcs x15, x6, x2', 'sbcs xzr, x7, xzr', 'csel x3, x12, x3, hs', 'csel x4, x13, x4, hs', 'csel x5, x14, x5, hs', 'csel x6, x15, x6, hs']
    L += ['subs x26, x26, #4', 'b.ne tb_loop', 'add x22, x22, #1', 'cmp x22, #8', 'b.ne tb_slot', 'ret']
    return '.globl tbld\ntbld:\n' + '\n'.join(_ind(L)) + '\n'
def unroll8(code, loop, done):
    """Eight uniquely labelled copies of a counted x19 loop (biased count, single-copy remainder); from 008389f5."""
    start = code.index(loop + ':\n')
    stop = code.index(done + ':\n', start)
    block = code[start + len(loop) + 2:stop]
    end = f' subs x19,x19,#1\n b.ne {loop}\n'
    assert block.endswith(end)
    body = block[:-len(end)]
    assert not re.search('\\bx19\\b', body)
    t, u, g = (loop + '_tail', loop + '_single', loop + '_eight')
    new = f' subs x19,x19,#8\n b.lo {t}\n{g}:\n'
    for k in range(8):
        b = body
        for lab in re.findall('^([A-Za-z_]\\w*):', body, re.M):
            b = re.sub('\\b' + lab + '\\b', lab + f'_u{k}', b)
        new += b
    new += f' subs x19,x19,#8\n b.hs {g}\n{t}:\n adds x19,x19,#8\n cbz x19,{done}\n{u}:\n' + body + f' subs x19,x19,#1\n b.ne {u}\n'
    return code[:start] + new + code[stop:]
def block_production(root):
    return unroll8((root / 'block.s').read_text(), 'pk_loop', 'pk_done')
DRIVER_PER_ROUND, DRIVER_FIXED = (68, 5)
def group_listing():
    """flip_group(x0 = parameter block): per round yprep, flip_fwd, k_inv, flip_bwd on the chains' current points, then
    the chain step ychain_f, k_inv, ychain_b; IDY and W2 from the round index. Its own work: 68 per round plus 5."""
    L = ['sub sp, sp, #32', 'stp x30, x0, [sp]', 'fg_round:', 'ldr x9, [sp, #8]', 'ldr x3, [x9, #160]', 'lsl x4, x3, #10', 'and x5, x4, #0xfffffff', 'lsr x6, x4, #28', 'lsl x6, x6, #45', 'ldr x7, [x9, #152]', 'orr x5, x7, x5, lsl #36', 'stp x6, x5, [x9, #56]', 'ldr x0, [x9]', 'ldr x1, [x9, #8]', 'ldp x2, x3, [x9, #40]', 'ldp x4, x5, [x9, #56]', 'bl yprep', 'ldr x9, [sp, #8]', 'ldr x0, [x9]', 'ldr x1, [x9, #40]', 'ldp x2, x3, [x9, #72]', 'ldr x4, [x9, #88]', 'bl flip_fwd', 'ldr x9, [sp, #8]', 'ldr x3, [x9]', 'ldr x4, [x9, #88]', 'mul x3, x3, x4', 'ldr x1, [x9, #72]', 'add x1, x1, x3', 'ldr x0, [x9, #96]', 'ldr x2, [x9, #32]', 'bl k_inv', 'ldr x12, [sp, #8]', 'ldr x0, [x12]', 'sub x6, x0, #1', 'ldr x1, [x12, #40]', 'add x1, x1, x6, lsl #10', 'ldr x11, [x12, #88]', 'mul x6, x0, x11', 'ldr x2, [x12, #72]', 'add x2, x2, x6', 'sub x2, x2, #64', 'ldp x4, x3, [x12, #136]', 'ldp x5, x7, [x12, #104]', 'ldp x8, x9, [x12, #120]', 'ldr x10, [x12, #96]', 'bl flip_bwd', 'ldr x9, [sp, #8]', 'ldr x0, [x9]', 'ldp x1, x2, [x9, #8]', 'bl ychain_f', 'ldr x9, [sp, #8]', 'ldr x3, [x9]', 'ldr x1, [x9, #16]', 'add x1, x1, x3, lsl #6', 'ldp x0, x2, [x9, #24]', 'bl k_inv', 'ldr x9, [sp, #8]', 'ldr x0, [x9]', 'ldp x1, x2, [x9, #8]', 'sub x4, x0, #1', 'add x1, x1, x4, lsl #7', 'add x2, x2, x4, lsl #6', 'ldr x3, [x9, #24]', 'bl ychain_b', 'ldr x9, [sp, #8]', 'ldp x3, x4, [x9, #160]', 'add x3, x3, #1', 'subs x4, x4, #1', 'stp x3, x4, [x9, #160]', 'b.ne fg_round', 'ldr x30, [sp]', 'add sp, sp, #32', 'ret']
    return '.globl flip_group\nflip_group:\n' + '\n'.join(_ind(L)) + '\n'

# Root independently derived all-one-block chain; weak REDC representatives preserved.
def inverse13(a):
 def sm(v,n,t):
  for _ in range(n):v=redc(v,v)[0]
  return redc(v,t)[0]
 t2=sm(a,1,a);t4=sm(t2,2,t2);t8=sm(t4,4,t4)
 t16=sm(t8,8,t8);t24=sm(t16,8,t8);t28=sm(t24,4,t4)
 t30=sm(t28,2,t2);t32=sm(t30,2,t2)
 v=sm(t32,32,a);v=sm(v,128,t32);v=sm(v,32,t32)
 return sm(sm(v,30,t30),2,a)
