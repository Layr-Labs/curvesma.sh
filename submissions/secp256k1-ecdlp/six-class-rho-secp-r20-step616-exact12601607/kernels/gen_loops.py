"""Emit hand-scheduled AArch64 assembly for the two counted batch loops.

usage: python3 gen_loops.py elf|macho [tests] > secp_loops.s

Field elements are four 64-bit limbs in registers. Straight-line schedules with explicit register
allocation from a pool; p = 2^256 - C0 with C0 = 2^32 + 977.

- fe_mul: 16 MUL + 16 UMULH and 30 ADDS/ADCS/ADC in six carry chains (62), reduction 23: 85.
- fe_sqr: 10 MUL + 10 UMULH, off-diagonal sum 8, doubling 6, diagonal chain 7, reduction 23: 64.
- fused reduction with k subtrahends: 24 + 5k, so lambda^2 - x_W - x_T is 41 + 34 = 75 and
  lambda*(x_W - x3) - y_W is 62 + 29 = 91.  fe_sub: 9.
The reductions make no final fold.  All are exact except on the rare inputs of proof.md section 2: z mod p
below C0^2 + C0 (minus the subtrahends, for a fused one), a fused fold that goes negative, or a fused result
in [p, 2^256).
"""
import sys

FLAVOR = sys.argv[1] if len(sys.argv) > 1 else 'elf'
POOL = ['x5', 'x6', 'x7', 'x8', 'x9', 'x10', 'x11', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17',
        'x19', 'x20', 'x21', 'x22', 'x23', 'x24', 'x25', 'x26', 'x27', 'x28', 'x29', 'x30']
# x19-x30 are saved by the prologue and restored before RET (x29 and x30 are free inside the loops,
# which make no calls); x0-x3 hold the count and array pointers; x18 is never written
CREG = 'x4'          # holds C0 = 2^32 + 977 inside the loops

# fe_mul product: six carry chains over the columns of a*b.  Lij / Hij are the low / high words of
# a_i*b_j (column i+j / i+j+1).  ('adds'|'adcs', dst, x, y) merges two words of one column; ('adc', dst, x)
# adds only the incoming carry to a fresh high word (<= 2^64 - 2) or, last, to the top column.
MUL_PLAN = [
    # chain 1: columns 1-4, carry into H13 (column 5)
    ('adds', 's1', 'H00', 'L01'), ('adcs', 's2', 'H01', 'L02'), ('adcs', 's3', 'H02', 'L03'),
    ('adcs', 's4', 'H03', 'L13'), ('adc', 'c5a', 'H13'),
    # chain 2: columns 1-4, carry into H22
    ('adds', 'z1', 's1', 'L10'), ('adcs', 's2', 's2', 'H10'), ('adcs', 's3', 's3', 'H20'),
    ('adcs', 's4', 's4', 'H30'), ('adc', 'c5b', 'H22'),
    # chain 3: columns 2-5, carry into H23 (column 6)
    ('adds', 's2', 's2', 'L20'), ('adcs', 's3', 's3', 'L30'), ('adcs', 's4', 's4', 'L22'),
    ('adcs', 's5', 'c5a', 'c5b'), ('adc', 'c6a', 'H23'),
    # chain 4: columns 2-5, carry into H32
    ('adds', 'z2', 's2', 'L11'), ('adcs', 's3', 's3', 'H11'), ('adcs', 's4', 's4', 'H12'),
    ('adcs', 's5', 's5', 'H31'), ('adc', 'c6b', 'H32'),
    # chain 5: columns 3-6, carry into H33 (column 7)
    ('adds', 's3', 's3', 'L12'), ('adcs', 's4', 's4', 'H21'), ('adcs', 's5', 's5', 'L23'),
    ('adcs', 's6', 'c6a', 'c6b'), ('adc', 'c7', 'H33'),
    # chain 6: columns 3-6, carry into the top column
    ('adds', 'z3', 's3', 'L21'), ('adcs', 'z4', 's4', 'L31'), ('adcs', 'z5', 's5', 'L32'),
    ('adcs', 'z6', 's6', 'L33'), ('adc', 'z7', 'c7'),
]


class Gen:
    def __init__(self, reserved=()):
        self.out = []
        self.free = [r for r in POOL if r not in reserved]
        self.size = len(self.free)
        self.peak = 0

    def ins(self, s):
        self.out.append('\t' + s)

    def alloc(self):
        if not self.free:
            raise RuntimeError('register pool exhausted at line %d' % len(self.out))
        r = self.free.pop(0)
        self.peak = max(self.peak, self.size - len(self.free))
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
        """z: 8 regs holding z = L + H*2^256; returns 4 regs with t mod 2^256, where H*C0 + L = s,
        s4 = s >> 256 <= C0 and t = (s mod 2^256) + s4*C0.  That is a weak value = z (mod p) unless t
        carries out of 2^256, which needs z mod p < C0^2 + C0 < 2^65; no final fold is made."""
        C = CREG
        g0 = self.alloc()
        self.ins(f'mul {g0}, {z[4]}, {C}')
        k0 = self.alloc()
        self.ins(f'umulh {k0}, {z[4]}, {C}')
        self.release(z[4])
        g1 = self.alloc()
        self.ins(f'mul {g1}, {z[5]}, {C}')
        self.ins(f'adds {g1}, {g1}, {k0}')
        self.release(k0)
        k1 = self.alloc()
        self.ins(f'umulh {k1}, {z[5]}, {C}')
        self.release(z[5])
        g2 = self.alloc()
        self.ins(f'mul {g2}, {z[6]}, {C}')
        self.ins(f'adcs {g2}, {g2}, {k1}')
        self.release(k1)
        k2 = self.alloc()
        self.ins(f'umulh {k2}, {z[6]}, {C}')
        self.release(z[6])
        g3 = self.alloc()
        self.ins(f'mul {g3}, {z[7]}, {C}')
        self.ins(f'adcs {g3}, {g3}, {k2}')
        self.release(k2)
        k3 = self.alloc()
        self.ins(f'umulh {k3}, {z[7]}, {C}')
        self.release(z[7])
        self.ins(f'adc {k3}, {k3}, xzr')               # H*C0 = (g0, g1, g2, g3, k3)
        self.ins(f'adds {z[0]}, {z[0]}, {g0}')
        self.ins(f'adcs {z[1]}, {z[1]}, {g1}')
        self.ins(f'adcs {z[2]}, {z[2]}, {g2}')
        self.ins(f'adcs {z[3]}, {z[3]}, {g3}')
        self.ins(f'adc {k3}, {k3}, xzr')               # s4
        self.release(g0, g1, g2, g3)
        w0 = self.alloc()
        w1 = self.alloc()
        self.ins(f'mul {w0}, {k3}, {C}')
        self.ins(f'umulh {w1}, {k3}, {C}')
        self.release(k3)
        d = out if out is not None else z[:4]
        self.ins(f'adds {d[0]}, {z[0]}, {w0}')
        self.ins(f'adcs {d[1]}, {z[1]}, {w1}')
        self.ins(f'adcs {d[2]}, {z[2]}, xzr')
        self.ins(f'adc {d[3]}, {z[3]}, xzr')
        self.release(w0, w1)
        if out is not None:
            self.release(z[:4])
        return d

    def reduce_sub(self, z, subs):
        """Fused reduction and subtraction.  z: 8 regs (z = L + H*2^256); subs: (base, off, keep) for each
        canonical subtrahend in memory.  s = L + H*C0 in five limbs, each subtrahend taken off all five
        (a kept one is loaded whole and returned, else read in two halves inside its borrow chain), then the
        signed fold t = (s' mod 2^256) + s4'*C0 with s4' in [-2, C0] by MUL, SMULH and a sign word.  Returns
        (t mod 2^256, kept registers); t is congruent to z minus the subtrahends unless t leaves [0, 2^256)."""
        C = CREG
        g0 = self.alloc()
        self.ins(f'mul {g0}, {z[4]}, {C}')
        k0 = self.alloc()
        self.ins(f'umulh {k0}, {z[4]}, {C}')
        self.release(z[4])
        g1 = self.alloc()
        self.ins(f'mul {g1}, {z[5]}, {C}')
        self.ins(f'adds {g1}, {g1}, {k0}')
        self.release(k0)
        k1 = self.alloc()
        self.ins(f'umulh {k1}, {z[5]}, {C}')
        self.release(z[5])
        g2 = self.alloc()
        self.ins(f'mul {g2}, {z[6]}, {C}')
        self.ins(f'adcs {g2}, {g2}, {k1}')
        self.release(k1)
        k2 = self.alloc()
        self.ins(f'umulh {k2}, {z[6]}, {C}')
        self.release(z[6])
        g3 = self.alloc()
        self.ins(f'mul {g3}, {z[7]}, {C}')
        self.ins(f'adcs {g3}, {g3}, {k2}')
        self.release(k2)
        s4 = self.alloc()
        self.ins(f'umulh {s4}, {z[7]}, {C}')
        self.release(z[7])
        self.ins(f'adc {s4}, {s4}, xzr')
        self.ins(f'adds {z[0]}, {z[0]}, {g0}')
        self.ins(f'adcs {z[1]}, {z[1]}, {g1}')
        self.ins(f'adcs {z[2]}, {z[2]}, {g2}')
        self.ins(f'adcs {z[3]}, {z[3]}, {g3}')
        self.ins(f'adc {s4}, {s4}, xzr')                # s = L + H*C0, s4 <= C0
        self.release(g0, g1, g2, g3)
        kept = []
        for base, off, keep in subs:
            if isinstance(base, list):              # subtrahend already in registers (released)
                self.ins(f'subs {z[0]}, {z[0]}, {base[0]}')
                self.ins(f'sbcs {z[1]}, {z[1]}, {base[1]}')
                self.ins(f'sbcs {z[2]}, {z[2]}, {base[2]}')
                self.ins(f'sbcs {z[3]}, {z[3]}, {base[3]}')
                self.release(base)
            elif keep:
                b = self.ld256(base, off)
                self.ins(f'subs {z[0]}, {z[0]}, {b[0]}')
                self.ins(f'sbcs {z[1]}, {z[1]}, {b[1]}')
                self.ins(f'sbcs {z[2]}, {z[2]}, {b[2]}')
                self.ins(f'sbcs {z[3]}, {z[3]}, {b[3]}')
                kept.append(b)
            else:
                c0, c1 = self.alloc(), self.alloc()
                self.ins(f'ldp {c0}, {c1}, [{base}, #{off}]')
                self.ins(f'subs {z[0]}, {z[0]}, {c0}')
                self.ins(f'sbcs {z[1]}, {z[1]}, {c1}')
                self.ins(f'ldp {c0}, {c1}, [{base}, #{off + 16}]')
                self.ins(f'sbcs {z[2]}, {z[2]}, {c0}')
                self.ins(f'sbcs {z[3]}, {z[3]}, {c1}')
                self.release(c0, c1)
            self.ins(f'sbc {s4}, {s4}, xzr')            # s4' = s4 - borrows, in [-2, C0]
        sw = self.alloc()
        self.ins(f'asr {sw}, {s4}, #63')
        w0, w1 = self.alloc(), self.alloc()
        self.ins(f'mul {w0}, {s4}, {C}')
        self.ins(f'smulh {w1}, {s4}, {C}')
        self.release(s4)
        d = z[:4]
        self.ins(f'adds {d[0]}, {z[0]}, {w0}')
        self.ins(f'adcs {d[1]}, {z[1]}, {w1}')
        self.ins(f'adcs {d[2]}, {z[2]}, {sw}')
        self.ins(f'adc {d[3]}, {z[3]}, {sw}')
        self.release(w0, w1, sw)
        return d, kept

    def mul(self, A, B, keep_a=False, keep_b=False, out=None, subs=None):
        """A*B mod p (weak).  Products are formed just before the chain step that consumes them; an
        operand limb is released after its last product unless kept."""
        a, b = list(A), list(B)
        prods = {}
        for st in MUL_PLAN:
            for w in st[2:]:
                if w[0] in 'LH':
                    prods[w] = None
        prods['L00'] = None
        uses = {}                            # remaining products per operand limb
        for w in prods:
            i, j = int(w[1]), int(w[2])
            uses[('a', i)] = uses.get(('a', i), 0) + 1
            uses[('b', j)] = uses.get(('b', j), 0) + 1
        assert sorted(prods) == sorted(f'{h}{i}{j}' for h in 'LH' for i in range(4) for j in range(4))
        val = {}

        def product(w):
            r = self.alloc()
            i, j = int(w[1]), int(w[2])
            self.ins(f'{"mul" if w[0] == "L" else "umulh"} {r}, {a[i]}, {b[j]}')
            for key, keep, regs, k in ((('a', i), keep_a, a, i), (('b', j), keep_b, b, j)):
                uses[key] -= 1
                if uses[key] == 0 and not keep and not (out is not None and regs is b):
                    self.release(regs[k])
            return r

        val['z0'] = product('L00')
        for st in MUL_PLAN:
            op, dst = st[0], st[1]
            src = [val.pop(w) if w in val else product(w) for w in st[2:]]
            if op == 'adc':
                self.ins(f'adc {src[0]}, {src[0]}, xzr')
            else:
                self.ins(f'{op} {src[0]}, {src[0]}, {src[1]}')
                self.release(src[1])
            val[dst] = src[0]
        assert sorted(val) == [f'z{k}' for k in range(8)]
        z = [val[f'z{k}'] for k in range(8)]
        if subs:
            return self.reduce_sub(z, subs)
        if out is not None:
            assert list(out) == b and not keep_b
            return self.reduce(z, out=out)          # B stays reserved; its registers receive the result
        return self.reduce(z)

    def sqr(self, A, keep_a=False, subs=None):
        a0, a1, a2, a3 = A

        def prod(op, i, j):
            r = self.alloc()
            self.ins(f'{op} {r}, {A[i]}, {A[j]}')
            return r
        # off-diagonal O = sum_{i<j} a_i a_j 2^(64(i+j)): two chains
        l03 = prod('mul', 0, 3); l12 = prod('mul', 1, 2)
        self.ins(f'adds {l03}, {l03}, {l12}')                 # column 3 partial
        self.release(l12)
        h03 = prod('umulh', 0, 3); h12 = prod('umulh', 1, 2)
        self.ins(f'adcs {h03}, {h03}, {h12}')                 # column 4 partial
        self.release(h12)
        h13 = prod('umulh', 1, 3)
        self.ins(f'adc {h13}, {h13}, xzr')                    # carry into H13 (column 5)
        h01 = prod('umulh', 0, 1); l02 = prod('mul', 0, 2)
        self.ins(f'adds {h01}, {h01}, {l02}')                 # o2
        self.release(l02)
        h02 = prod('umulh', 0, 2)
        self.ins(f'adcs {l03}, {l03}, {h02}')                 # o3
        self.release(h02)
        l13 = prod('mul', 1, 3)
        self.ins(f'adcs {h03}, {h03}, {l13}')                 # o4
        self.release(l13)
        l23 = prod('mul', 2, 3)
        self.ins(f'adcs {h13}, {h13}, {l23}')                 # o5
        self.release(l23)
        h23 = prod('umulh', 2, 3)
        self.ins(f'adc {h23}, {h23}, xzr')                    # o6
        o1 = prod('mul', 0, 1)
        o2, o3, o4, o5, o6 = h01, l03, h03, h13, h23
        # doubling: d7..d2 by shifts; d1 = o1 << 1 is folded into the first diagonal addition
        d7 = self.alloc()
        self.ins(f'lsr {d7}, {o6}, #63')
        self.ins(f'extr {o6}, {o6}, {o5}, #63')
        self.ins(f'extr {o5}, {o5}, {o4}, #63')
        self.ins(f'extr {o4}, {o4}, {o3}, #63')
        self.ins(f'extr {o3}, {o3}, {o2}, #63')
        self.ins(f'extr {o2}, {o2}, {o1}, #63')
        # diagonals a_i^2 and the final chain
        z0 = prod('mul', 0, 0)
        g = prod('umulh', 0, 0)
        self.ins(f'adds {o1}, {g}, {o1}, lsl #1')             # z1
        self.release(g)
        diag = [('mul', 1), ('umulh', 1), ('mul', 2), ('umulh', 2), ('mul', 3)]
        zs = [o2, o3, o4, o5, o6]
        for (op, i), zr in zip(diag, zs):
            g = prod(op, i, i)
            self.ins(f'adcs {zr}, {zr}, {g}')
            self.release(g)
        g = prod('umulh', 3, 3)
        if not keep_a:
            self.release(list(A))
        self.ins(f'adc {d7}, {d7}, {g}')                      # z7
        self.release(g)
        if subs:
            return self.reduce_sub([z0, o1, o2, o3, o4, o5, o6, d7], subs)
        return self.reduce([z0, o1, o2, o3, o4, o5, o6, d7])

    def sub(self, A, B, keep_a=False, keep_b=False):
        """A - B mod p, B canonical; result in A's registers if A is consumed, else new registers."""
        C = CREG
        d = [self.alloc() for _ in range(4)] if keep_a else A
        self.ins(f'subs {d[0]}, {A[0]}, {B[0]}')
        self.ins(f'sbcs {d[1]}, {A[1]}, {B[1]}')
        self.ins(f'sbcs {d[2]}, {A[2]}, {B[2]}')
        self.ins(f'sbcs {d[3]}, {A[3]}, {B[3]}')
        if not keep_b:
            self.release(B)
        m = self.alloc()
        self.ins(f'csel {m}, {C}, xzr, lo')
        self.ins(f'subs {d[0]}, {d[0]}, {m}')
        self.ins(f'sbcs {d[1]}, {d[1]}, xzr')
        self.ins(f'sbcs {d[2]}, {d[2]}, xzr')
        self.ins(f'sbc {d[3]}, {d[3]}, xzr')
        self.release(m)
        return d


def finish(g, w_y, t_ptr, inv, last_post):
    """W <- W + T given inv = 1/dx (regs, consumed). w_y: base register pointing at W.y (x at -32).
    t_ptr: register pointing at T (x at 0, y at 32). Stores x3,y3 into W.
    last_post: post-index applied by the final store (None for the peeled element)."""
    ty = g.ld256(t_ptr, 32); wy = g.ld256(w_y, 0)
    dy = g.sub(ty, wy, keep_b=True)          # keep y_W for y3
    lam = g.mul(dy, inv)
    x3, (wx,) = g.sqr(lam, keep_a=True, subs=[(w_y, -32, True), (t_ptr, 0, False)])   # lambda^2 - x_W - x_T
    tt = g.sub(wx, x3, keep_b=True)          # x_W - x3; consumes x_W
    g.st256(x3, w_y, -32)
    g.release(x3)
    y3, _ = g.mul(lam, tt, subs=[(wy, None, False)])   # lambda*(x_W - x3) - y_W; consumes y_W
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


def prologue():
    return ['\tstp x19, x20, [sp, #-96]!', '\tstp x21, x22, [sp, #16]', '\tstp x23, x24, [sp, #32]',
            '\tstp x25, x26, [sp, #48]', '\tstp x27, x28, [sp, #64]', '\tstp x29, x30, [sp, #80]']


def epilogue():
    return ['\tldp x29, x30, [sp, #80]', '\tldp x27, x28, [sp, #64]', '\tldp x25, x26, [sp, #48]',
            '\tldp x23, x24, [sp, #32]', '\tldp x21, x22, [sp, #16]', '\tldp x19, x20, [sp], #96', '\tret']


def const_c0():
    return [f'\tmov {CREG}, #977', f'\tmovk {CREG}, #1, lsl #32']


def gen_pass1():
    """pass1(n, w, tp, sl): x0=n, x1=w, x2=tp, x3=sl. Slot i (64 bytes) = [dx_i | acc_{i-1}], slot n = [- | acc_{n-1}].
    acc_i = acc_{i-1} * dx_i is carried in registers F; slot 0 is not written (pass 2 uses R for walk 0)."""
    F = ['x5', 'x6', 'x7', 'x8']
    L = header('pass1') + prologue() + const_c0()
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
    L += epilogue()
    return L, body, g.peak + len(F)


def gen_pass2():
    """pass2(n, w, tp, sl, Rinv): replaces W_i by W_i + T_i. Slot i = [dx_i | acc_{i-1}] as written by pass1.
    R (= 1/(dx_0 ... dx_i) at iteration i) is carried in registers across iterations."""
    R = ['x25', 'x26', 'x27', 'x28']
    L = header('pass2') + prologue()
    L += [f'\tldp {R[0]}, {R[1]}, [x4]', f'\tldp {R[2]}, {R[3]}, [x4, #16]']
    L += const_c0()
    L += ['\tsub x5, x0, #1',                   # i = n-1
          '\tadd x1, x1, x5, lsl #6',           # &w[n-1]
          '\tadd x1, x1, #32',                  # -> its y limbs
          '\tadd x2, x2, x5, lsl #3',           # &tp[n-1]
          '\tadd x3, x3, x5, lsl #6',           # &slot[n-1]
          '\tmov x0, x5']                       # iterations i = n-1 .. 1
    g = Gen(reserved=R)
    L.append('1:')
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2], #-8')
    accv = g.ld256('x3', 32)
    inv = g.mul(R, accv, keep_a=True)
    dx = g.ld256('x3', 0, post=-64)
    g.mul(dx, R, out=R)                         # R <- R * dx_i
    finish(g, 'x1', t_ptr, inv, -64)
    g.release(t_ptr)
    L += g.out
    L += ['\tsubs x0, x0, #1', '\tb.ne 1b']
    body = len(g.out) + 2
    peak = g.peak + len(R)
    g = Gen(reserved=R)
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2]')
    finish(g, 'x1', t_ptr, list(R), None)      # walk 0 uses R itself as its inverse
    g.release(t_ptr)
    L += g.out
    L += epilogue()
    return L, body, max(peak, g.peak + len(R))


def gen_tests():
    """Unit-test wrappers (r=x0, a=x1, b=x2) around the same emitters; t_reduce reads 8 limbs at x1."""
    L = []
    for name in ('t_mul', 't_sqr', 't_sub', 't_reduce', 't_sqr_sub2', 't_mul_sub'):
        L += header(name) + prologue() + const_c0()
        g = Gen()
        if name == 't_reduce':
            z = g.ld256('x1', 0) + g.ld256('x1', 32)
            r = g.reduce(z)
        else:
            a = g.ld256('x1', 0)
            if name == 't_mul':
                b = g.ld256('x2', 0)
                r = g.mul(a, b)
            elif name == 't_sqr_sub2':        # r = a^2 - b - c, b at [x2], c at [x3]
                r, _ = g.sqr(a, subs=[('x2', 0, False), ('x3', 0, False)])
            elif name == 't_mul_sub':         # r = a*b - c, c at [x3]
                b = g.ld256('x2', 0)
                r, _ = g.mul(a, b, subs=[('x3', 0, False)])
            elif name == 't_sqr':
                r = g.sqr(a)
            else:
                b = g.ld256('x2', 0)
                r = g.sub(a, b)
        g.st256(r, 'x0', 0)
        L += g.out
        L += epilogue()
    return L


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[2] == 'tests':
        print('\n'.join(gen_tests()))
        sys.exit(0)
    p1, b1, k1 = gen_pass1()
    p2, b2, k2 = gen_pass2()
    sys.stderr.write(f'pass1 loop {b1} peak {k1}; pass2 loop {b2} peak {k2}\n')
    print('\n'.join(p1 + p2))
