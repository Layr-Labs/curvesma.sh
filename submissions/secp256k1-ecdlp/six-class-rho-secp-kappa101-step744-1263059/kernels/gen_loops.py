"""Emit hand-scheduled AArch64 assembly for the two counted batch loops of secp_walk.c.

usage: python3 gen_loops.py elf|macho > secp_loops.s

Field elements are four 64-bit limbs in registers. Straight-line schedules with explicit register
allocation from a pool; products are row-wise schoolbook with ADDS/ADCS chains, reduction by
2^256 = C0 (mod p), exactly the algorithm of secp_walk.c (fe_mul, fe_sqr, fe_sub, fe_canon).
"""
import sys

FLAVOR = sys.argv[1] if len(sys.argv) > 1 else 'elf'
POOL = ['x5', 'x6', 'x7', 'x8', 'x9', 'x10', 'x11', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17',
        'x19', 'x20', 'x21', 'x22', 'x23', 'x24', 'x25', 'x26', 'x27', 'x28']
CREG = 'x4'          # holds C0 = 2^32 + 977 inside the loops


class Gen:
    def __init__(self, reserved=()):
        self.out = []
        self.free = [r for r in POOL if r not in reserved]
        self.used = set()
        self.peak = 0

    def ins(self, s):
        self.out.append('\t' + s)

    def alloc(self):
        if not self.free:
            raise RuntimeError('register pool exhausted at line %d' % len(self.out))
        r = self.free.pop(0)
        self.used.add(r)
        self.peak = max(self.peak, len(POOL) - len(self.free))
        return r

    def release(self, *regs):
        for r in regs:
            if isinstance(r, list):
                self.release(*r)
            else:
                assert r not in self.free and r in POOL, r
                self.free.insert(0, r)

    # ---- memory
    def ld256(self, base, off, post=None, into=None):
        r = list(into) if into is not None else [self.alloc() for _ in range(4)]
        if post is None:
            self.ins(f'ldp {r[0]}, {r[1]}, [{base}, #{off}]')
            self.ins(f'ldp {r[2]}, {r[3]}, [{base}, #{off + 16}]')
        else:   # last access at offset 0 with post-index
            assert off == 0
            self.ins(f'ldp {r[2]}, {r[3]}, [{base}, #16]')
            self.ins(f'ldp {r[0]}, {r[1]}, [{base}], #{post}')
        return r

    def st256(self, v, base, off, post=None):
        if post is None:
            self.ins(f'stp {v[0]}, {v[1]}, [{base}, #{off}]')
            self.ins(f'stp {v[2]}, {v[3]}, [{base}, #{off + 16}]')
        else:
            assert off == 0
            self.ins(f'stp {v[2]}, {v[3]}, [{base}, #16]')
            self.ins(f'stp {v[0]}, {v[1]}, [{base}], #{post}')

    # ---- field operations
    def reduce(self, z, out=None):
        """z: 8 regs (z0..z7) holding the 512-bit product; returns 4 regs (z0..z3 updated in place, or `out`)."""
        C = CREG
        g = [self.alloc() for _ in range(4)]
        k = [self.alloc() for _ in range(4)]
        for j in range(4):
            self.ins(f'mul {g[j]}, {z[4 + j]}, {C}')
            self.ins(f'umulh {k[j]}, {z[4 + j]}, {C}')
        self.release(z[4:])
        # q1..q4 in k/g: q1 = g1 + k0, q2 = g2 + k1 + c, q3 = g3 + k2 + c, q4 = k3 + c
        self.ins(f'adds {g[1]}, {g[1]}, {k[0]}')
        self.ins(f'adcs {g[2]}, {g[2]}, {k[1]}')
        self.ins(f'adcs {g[3]}, {g[3]}, {k[2]}')
        self.ins(f'adc {k[3]}, {k[3]}, xzr')
        self.release(k[0], k[1], k[2])
        # s = z0..z3 + (g0, q1, q2, q3) + q4 * 2^256
        self.ins(f'adds {z[0]}, {z[0]}, {g[0]}')
        self.ins(f'adcs {z[1]}, {z[1]}, {g[1]}')
        self.ins(f'adcs {z[2]}, {z[2]}, {g[2]}')
        self.ins(f'adcs {z[3]}, {z[3]}, {g[3]}')
        self.ins(f'adc {k[3]}, {k[3]}, xzr')         # s4
        self.release(g)
        w0 = self.alloc()
        w1 = self.alloc()
        self.ins(f'mul {w0}, {k[3]}, {C}')
        self.ins(f'umulh {w1}, {k[3]}, {C}')
        self.release(k[3])
        self.ins(f'adds {z[0]}, {z[0]}, {w0}')
        self.ins(f'adcs {z[1]}, {z[1]}, {w1}')
        self.ins(f'adcs {z[2]}, {z[2]}, xzr')
        self.ins(f'adcs {z[3]}, {z[3]}, xzr')
        self.ins(f'csel {w0}, {C}, xzr, cs')
        d = out if out is not None else z[:4]
        self.ins(f'adds {d[0]}, {z[0]}, {w0}')
        self.ins(f'adcs {d[1]}, {z[1]}, xzr')
        self.ins(f'adcs {d[2]}, {z[2]}, xzr')
        self.ins(f'adc {d[3]}, {z[3]}, xzr')
        self.release(w0, w1)
        if out is not None:
            self.release(z[:4])
        return d

    def mul(self, A, B, keep_a=False, keep_b=False, out=None):
        """A*B mod p (weak). Frees A rows as they are consumed (unless keep_a) and B at the end (unless keep_b)."""
        z = [None] * 8
        a = list(A)
        q = self.alloc()
        # row 0 (12 instructions)
        z[0] = self.alloc(); z[1] = self.alloc(); z[2] = self.alloc(); z[3] = self.alloc(); z[4] = self.alloc()
        self.ins(f'mul {z[0]}, {a[0]}, {B[0]}')
        self.ins(f'umulh {z[1]}, {a[0]}, {B[0]}')
        self.ins(f'mul {q}, {a[0]}, {B[1]}')
        self.ins(f'umulh {z[2]}, {a[0]}, {B[1]}')
        self.ins(f'adds {z[1]}, {z[1]}, {q}')
        self.ins(f'mul {q}, {a[0]}, {B[2]}')
        self.ins(f'umulh {z[3]}, {a[0]}, {B[2]}')
        self.ins(f'adcs {z[2]}, {z[2]}, {q}')
        self.ins(f'mul {q}, {a[0]}, {B[3]}')
        self.ins(f'umulh {z[4]}, {a[0]}, {B[3]}')
        self.ins(f'adcs {z[3]}, {z[3]}, {q}')
        self.ins(f'adc {z[4]}, {z[4]}, xzr')
        if not keep_a:
            self.release(a[0])
        for i in range(1, 4):
            p = [self.alloc() for _ in range(5)]     # p0..p4 ; p4 becomes z[i+4]
            self.ins(f'mul {p[0]}, {a[i]}, {B[0]}')
            self.ins(f'umulh {p[1]}, {a[i]}, {B[0]}')
            self.ins(f'mul {q}, {a[i]}, {B[1]}')
            self.ins(f'umulh {p[2]}, {a[i]}, {B[1]}')
            self.ins(f'adds {p[1]}, {p[1]}, {q}')
            self.ins(f'mul {q}, {a[i]}, {B[2]}')
            self.ins(f'umulh {p[3]}, {a[i]}, {B[2]}')
            self.ins(f'adcs {p[2]}, {p[2]}, {q}')
            self.ins(f'mul {q}, {a[i]}, {B[3]}')
            self.ins(f'umulh {p[4]}, {a[i]}, {B[3]}')
            self.ins(f'adcs {p[3]}, {p[3]}, {q}')
            self.ins(f'adc {p[4]}, {p[4]}, xzr')
            self.ins(f'adds {z[i]}, {z[i]}, {p[0]}')
            self.ins(f'adcs {z[i + 1]}, {z[i + 1]}, {p[1]}')
            self.ins(f'adcs {z[i + 2]}, {z[i + 2]}, {p[2]}')
            self.ins(f'adcs {z[i + 3]}, {z[i + 3]}, {p[3]}')
            self.ins(f'adc {p[4]}, {p[4]}, xzr')       # z[i+4] = p4 + carry
            z[i + 4] = p[4]
            self.release(p[0], p[1], p[2], p[3])
            if not keep_a:
                self.release(a[i])
        self.release(q)
        if out is not None:
            assert list(out) == list(B) and not keep_b
            return self.reduce(z, out=out)          # B stays reserved; its registers receive the result
        if not keep_b:
            self.release(B)
        return self.reduce(z)

    def sqr(self, A, keep_a=False):
        a0, a1, a2, a3 = A
        t = {}
        # off-diagonal products
        def prod(i, j):
            lo, hi = self.alloc(), self.alloc()
            self.ins(f'mul {lo}, {A[i]}, {A[j]}')
            self.ins(f'umulh {hi}, {A[i]}, {A[j]}')
            return lo, hi
        l01, h01 = prod(0, 1); l02, h02 = prod(0, 2); l03, h03 = prod(0, 3)
        l12, h12 = prod(1, 2); l13, h13 = prod(1, 3); l23, h23 = prod(2, 3)
        # o1 = l01 ; o2 = h01 + l02 ; o3 = h02 + l03 + c ; o4 = h03 + c
        self.ins(f'adds {h01}, {h01}, {l02}')               # o2
        self.ins(f'adcs {h02}, {h02}, {l03}')               # o3
        self.ins(f'adc {h03}, {h03}, xzr')                  # o4
        self.release(l02, l03)
        # p3 = l12 ; p4 = h12 + l13 ; p5 = h13 + c
        self.ins(f'adds {h12}, {h12}, {l13}')               # p4
        self.ins(f'adc {h13}, {h13}, xzr')                  # p5
        self.release(l13)
        # o3 += p3 ; o4 += p4 + c ; o5 = p5 + l23 + c ; o6 = h23 + c
        self.ins(f'adds {h02}, {h02}, {l12}')
        self.ins(f'adcs {h03}, {h03}, {h12}')
        self.ins(f'adcs {h13}, {h13}, {l23}')               # o5
        self.ins(f'adc {h23}, {h23}, xzr')                  # o6
        self.release(l12, h12, l23)
        o1, o2, o3, o4, o5, o6 = l01, h01, h02, h03, h13, h23
        # doubling -> d1..d7
        d7 = self.alloc()
        self.ins(f'lsr {d7}, {o6}, #63')
        self.ins(f'extr {o6}, {o6}, {o5}, #63')            # d6
        self.ins(f'extr {o5}, {o5}, {o4}, #63')            # d5
        self.ins(f'extr {o4}, {o4}, {o3}, #63')            # d4
        self.ins(f'extr {o3}, {o3}, {o2}, #63')            # d3
        self.ins(f'extr {o2}, {o2}, {o1}, #63')            # d2
        self.ins(f'lsl {o1}, {o1}, #1')                    # d1
        d1, d2, d3, d4, d5, d6 = o1, o2, o3, o4, o5, o6
        # diagonals
        z = [None] * 8
        z[0] = self.alloc()
        g = [self.alloc() for _ in range(4)]
        y = [self.alloc() for _ in range(3)]
        self.ins(f'mul {z[0]}, {a0}, {a0}')
        self.ins(f'umulh {g[0]}, {a0}, {a0}')
        self.ins(f'mul {y[0]}, {a1}, {a1}')
        self.ins(f'umulh {g[1]}, {a1}, {a1}')
        self.ins(f'mul {y[1]}, {a2}, {a2}')
        self.ins(f'umulh {g[2]}, {a2}, {a2}')
        self.ins(f'mul {y[2]}, {a3}, {a3}')
        self.ins(f'umulh {g[3]}, {a3}, {a3}')
        if not keep_a:
            self.release(list(A))
        self.ins(f'adds {d1}, {d1}, {g[0]}')               # z1
        self.ins(f'adcs {d2}, {d2}, {y[0]}')               # z2
        self.ins(f'adcs {d3}, {d3}, {g[1]}')               # z3
        self.ins(f'adcs {d4}, {d4}, {y[1]}')               # z4
        self.ins(f'adcs {d5}, {d5}, {g[2]}')               # z5
        self.ins(f'adcs {d6}, {d6}, {y[2]}')               # z6
        self.ins(f'adc {d7}, {d7}, {g[3]}')                # z7
        self.release(g, y)
        z[1], z[2], z[3], z[4], z[5], z[6], z[7] = d1, d2, d3, d4, d5, d6, d7
        return self.reduce(z)

    def sub(self, A, B, keep_a=False, keep_b=False):
        """A - B mod p, B canonical; result in A's registers if A is consumed, else new registers."""
        C = CREG
        if keep_a:
            d = [self.alloc() for _ in range(4)]
        else:
            d = A
        m = self.alloc()
        self.ins(f'subs {d[0]}, {A[0]}, {B[0]}')
        self.ins(f'sbcs {d[1]}, {A[1]}, {B[1]}')
        self.ins(f'sbcs {d[2]}, {A[2]}, {B[2]}')
        self.ins(f'sbcs {d[3]}, {A[3]}, {B[3]}')
        self.ins(f'csel {m}, {C}, xzr, lo')
        self.ins(f'subs {d[0]}, {d[0]}, {m}')
        self.ins(f'sbcs {d[1]}, {d[1]}, xzr')
        self.ins(f'sbcs {d[2]}, {d[2]}, xzr')
        self.ins(f'sbc {d[3]}, {d[3]}, xzr')
        self.release(m)
        if not keep_b:
            self.release(B)
        return d

    def canon(self, V):
        C = CREG
        t = [self.alloc() for _ in range(4)]
        self.ins(f'adds {t[0]}, {V[0]}, {C}')
        self.ins(f'adcs {t[1]}, {V[1]}, xzr')
        self.ins(f'adcs {t[2]}, {V[2]}, xzr')
        self.ins(f'adcs {t[3]}, {V[3]}, xzr')
        for j in range(4):
            self.ins(f'csel {V[j]}, {t[j]}, {V[j]}, cs')
        self.release(t)
        return V


def finish(g, w_y, t_ptr, inv, last_post):
    """W <- W + T given inv = 1/dx (regs, consumed). w_y: base register pointing at W.y (x at -32).
    t_ptr: register pointing at T (x at 0, y at 32). Returns nothing; stores x3,y3 into W.
    last_post: post-index applied by the final store (None for the peeled element)."""
    ty = g.ld256(t_ptr, 32); wy = g.ld256(w_y, 0)
    dy = g.sub(ty, wy)                       # consumes ty, frees wy
    lam = g.mul(dy, inv)
    l2 = g.sqr(lam, keep_a=True)
    wx = g.ld256(w_y, -32)
    x3 = g.sub(l2, wx, keep_b=True)          # keep wx for tt
    tx = g.ld256(t_ptr, 0)
    x3 = g.sub(x3, tx)
    g.canon(x3)
    tt = g.sub(wx, x3, keep_b=True)          # consumes the kept wx
    g.st256(x3, w_y, -32)
    g.release(x3)
    uu = g.mul(lam, tt)
    wy = g.ld256(w_y, 0)
    y3 = g.sub(uu, wy)
    g.canon(y3)
    if last_post is None:
        g.st256(y3, w_y, 0)
    else:
        g.ins(f'stp {y3[2]}, {y3[3]}, [{w_y}, #16]')
        g.ins(f'stp {y3[0]}, {y3[1]}, [{w_y}], #{last_post}')
    g.release(y3)


# ------------------------------------------------------------------ functions
def sym(name):
    return ('_' + name) if FLAVOR == 'macho' else name


def header(name):
    lines = ['\t.text', '\t.p2align 2', f'\t.globl {sym(name)}']
    if FLAVOR == 'elf':
        lines.append(f'\t.type {name}, %function')
    lines.append(f'{sym(name)}:')
    return lines


SAVED = ['x19', 'x20', 'x21', 'x22', 'x23', 'x24', 'x25', 'x26', 'x27', 'x28']


def prologue(frame):
    L = [f'\tstp x19, x20, [sp, #-{frame}]!', '\tstp x21, x22, [sp, #16]', '\tstp x23, x24, [sp, #32]',
         '\tstp x25, x26, [sp, #48]', '\tstp x27, x28, [sp, #64]']
    return L


def epilogue(frame):
    return ['\tldp x27, x28, [sp, #64]', '\tldp x25, x26, [sp, #48]', '\tldp x23, x24, [sp, #32]',
            '\tldp x21, x22, [sp, #16]', f'\tldp x19, x20, [sp], #{frame}', '\tret']


def const_c0():
    return [f'\tmov {CREG}, #977', f'\tmovk {CREG}, #1, lsl #32']


def gen_pass1():
    """pass1(n, w, tp, sl): x0=n, x1=w, x2=tp, x3=sl. Slot i (64 bytes) = [dx_i | acc_{i-1}], slot n = [- | acc_{n-1}].
    acc_i = acc_{i-1} * dx_i is carried in registers F; slot 0 is not written (pass 2 uses R for walk 0)."""
    F = ['x5', 'x6', 'x7', 'x8']
    L = header('pass1') + prologue(80) + const_c0()
    g = Gen(reserved=F)
    t = g.alloc()
    g.ins(f'ldr {t}, [x2], #8')
    g.ld256(t, 0, into=F)
    wx = g.ld256('x1', 0, post=64)
    g.sub(F, wx)                       # F = dx_0 = acc_0 ; wx released by sub
    g.release(t)
    L += g.out
    L += ['\tadd x3, x3, #64', '\tsub x0, x0, #1']
    g = Gen(reserved=F)
    L.append('1:')
    t = g.alloc()
    g.ins(f'ldr {t}, [x2], #8')
    tx = g.ld256(t, 0)
    wx = g.ld256('x1', 0, post=64)
    dx = g.sub(tx, wx)
    g.release(t)
    g.st256(F, 'x3', 32)               # acc_{i-1}
    g.st256(dx, 'x3', 0, post=64)      # dx_i, and advance to slot i+1
    g.mul(dx, F, out=F)                # F <- dx * F
    L += g.out
    L += ['\tsubs x0, x0, #1', '\tb.ne 1b']
    body = len(g.out) + 2
    # acc_{n-1} for the inversion
    g2 = Gen(reserved=F)
    g2.st256(F, 'x3', 32)
    L += g2.out
    L += epilogue(80)
    return L, body, g.peak


def gen_pass2():
    """pass2(n, w, tp, sl, Rinv): replaces W_i by W_i + T_i. Slot i = [dx_i | acc_{i-1}] as written by pass1."""
    L = header('pass2') + prologue(112)
    L += ['\tldp x5, x6, [x4]', '\tldp x7, x8, [x4, #16]', '\tstp x5, x6, [sp, #80]', '\tstp x7, x8, [sp, #96]']
    L += const_c0()
    L += ['\tsub x5, x0, #1',                   # i = n-1
          '\tadd x1, x1, x5, lsl #6',           # &w[n-1]
          '\tadd x1, x1, #32',                  # -> its y limbs
          '\tadd x2, x2, x5, lsl #3',           # &tp[n-1]
          '\tadd x3, x3, x5, lsl #6',           # &slot[n-1]
          '\tmov x0, x5']                       # iterations i = n-1 .. 1
    g = Gen()
    L.append('1:')
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2], #-8')
    R = g.ld256('sp', 80)
    accv = g.ld256('x3', 32)
    inv = g.mul(R, accv, keep_a=True)
    dx = g.ld256('x3', 0, post=-64)
    R2 = g.mul(R, dx)
    g.st256(R2, 'sp', 80)
    g.release(R2)
    finish(g, 'x1', t_ptr, inv, -64)
    g.release(t_ptr)
    L += g.out
    L += ['\tsubs x0, x0, #1', '\tb.ne 1b']
    body = len(g.out) + 2
    peak = g.peak
    g = Gen()
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2]')
    inv = g.ld256('sp', 80)
    finish(g, 'x1', t_ptr, inv, None)
    g.release(t_ptr)
    L += g.out
    L += epilogue(112)
    return L, body, max(peak, g.peak)


def gen_tests():
    """Unit-test wrappers t_mul/t_sqr/t_sub/t_canon (r=x0, a=x1, b=x2) around the same emitters."""
    L = []
    for name in ('t_mul', 't_sqr', 't_sub', 't_canon'):
        L += header(name) + prologue(80) + const_c0()
        g = Gen()
        a = g.ld256('x1', 0)
        if name == 't_mul':
            b = g.ld256('x2', 0)
            r = g.mul(a, b)
        elif name == 't_sqr':
            r = g.sqr(a)
        elif name == 't_sub':
            b = g.ld256('x2', 0)
            r = g.sub(a, b)
        else:
            r = g.canon(a)
        g.st256(r, 'x0', 0)
        L += g.out
        L += epilogue(80)
    return L


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[2] == 'tests':
        print('\n'.join(gen_tests()))
        sys.exit(0)
    p1, b1, k1 = gen_pass1()
    p2, b2, k2 = gen_pass2()
    sys.stderr.write(f'pass1 loop {b1} peak {k1}; pass2 loop {b2} peak {k2}\n')
    print('\n'.join(p1 + p2))
