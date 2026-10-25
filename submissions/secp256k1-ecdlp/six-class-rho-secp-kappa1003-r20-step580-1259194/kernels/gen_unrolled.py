"""Emit the two counted batch loops unrolled (pass 1 eight times, pass 2 four times), with rare exits.

usage: python3 gen_unrolled.py elf|macho > secp_loops.s     (needs 2d25756c's kernels/gen_loops.py on the path)

The emitters, walk bodies and register allocation are gen_loops.py's.  As in kaankolcu's a9297fb0, three tails
become rare exits (proof.md 2.5): a branch to a stub after RET that charges 20 milli-units (macro RCHG), runs the
tail's own instructions and branches back, so every result is bit-identical to the inline kernels.  Walk bodies
are the macros P1WALK k, p and P2WALK k, p (k numbers the copy, p holds the walk's table pointer); their stubs are
P1RARE k and P2RARE k.  A group loads two walks' table pointers with one LDP; a tail walk loads its own with LDR.
Operand lists are printed without spaces, except on branch, compare and select lines.
"""
import sys

import gen_loops as G

# 20 milli-units from the caller's x19 (the round's meter remainder), which the prologue saved at [sp]; sets no
# flags and restores x0
RCHG = ['\t.macro RCHG', '\tstr x0, [sp, #-16]!', '\tldr x0, [sp, #16]', '\tsub x0, x0, #20',
        '\tstr x0, [sp, #16]', '\tldr x0, [sp], #16', '\t.endm']


class RGen(G.Gen):
    """gen_loops.Gen whose fe_sub, plain reduction and one-subtrahend fused reduction end in rare exits."""
    def __init__(self, reserved=(), tag='1', k='\\k'):
        super().__init__(reserved)
        self.tag, self.k, self.rare, self.n = tag, k, [], 0

    def label(self, kind):
        return f'.L{kind}{self.tag}_{self.n}_{self.k}'

    def exit(self, branch, tail):
        """Emit `branch` to a stub that charges, runs `tail` and returns; returns the label after the exit."""
        r, b = self.label('r'), self.label('b')
        self.n += 1
        self.ins(f'{branch} {r}')
        self.out.append(f'{b}:')
        self.rare += [f'{r}:', '\tRCHG'] + tail + [f'\tb {b}']
        return b

    def sub(self, A, B, keep_a=False, keep_b=False):
        d = G.Gen.sub(self, A, B, keep_a, keep_b)
        tail = self.out[-3:]                     # SBCS d1, SBCS d2, SBC d3 of the second borrow chain
        assert [s.split()[0] for s in tail] == ['sbcs', 'sbcs', 'sbc']
        del self.out[-3:]
        self.exit('b.lo', tail)
        return d

    def reduce(self, z, out=None):
        d = G.Gen.reduce(self, z, out)
        if out is not None:                      # step 2 writes limbs 2 and 3 straight into the output registers
            for j in (2, 3):
                i = max(i for i, s in enumerate(self.out) if s.startswith(f'\tadcs {z[j]}, {z[j]}, '))
                self.out[i] = self.out[i].replace(f'adcs {z[j]},', f'adcs {d[j]},', 1)
        assert [s.split()[0] for s in self.out[-2:]] == ['adcs', 'adc']
        del self.out[-2:]
        self.exit('b.cs', [f'\tadcs {d[2]}, {d[2]}, xzr', f'\tadc {d[3]}, {d[3]}, xzr'])
        return d

    def reduce_sub(self, z, subs):
        d, kept = G.Gen.reduce_sub(self, z, subs)
        if len(subs) != 1:                       # x3's fold stays inline: lambda^2 makes s4 = 0 too often
            return d, kept
        fold = self.out[-7:]                     # ASR sw; MUL w0; SMULH w1; ADDS; ADCS; ADCS ..., sw; ADC ..., sw
        assert [s.split()[0] for s in fold] == ['asr', 'mul', 'smulh', 'adds', 'adcs', 'adcs', 'adc']
        del self.out[-7:]
        s4 = fold[1].split()[2].rstrip(',')
        w0, w1 = fold[1].split()[1].rstrip(','), fold[2].split()[1].rstrip(',')
        r = self.label('r')
        self.n += 1
        self.ins(f'tbnz {s4}, #63, {r}')         # s4' < 0: the signed fold, in the stub
        self.ins(f'mul {w0}, {s4}, {G.CREG}')    # s4' in [0, C0]: SMULH = UMULH and the sign word is 0
        self.ins(f'umulh {w1}, {s4}, {G.CREG}')
        self.out += fold[3:5]
        b = self.exit('b.cs', [f'\tadcs {d[2]}, {d[2]}, xzr', f'\tadc {d[3]}, {d[3]}, xzr'])
        self.rare += [f'{r}:', '\tRCHG'] + fold + [f'\tb {b}']
        return d, kept


def p1_walk(F, k='\\k'):
    g = RGen(reserved=F, tag='1', k=k)
    t = g.alloc()                                # x9 stays reserved while the table pointer (parameter p) is read
    tx = g.ld256('\\p', 0)
    wx = g.ld256('x1', 0, post=64)
    dx = g.sub(tx, wx)
    g.release(t)
    g.st256(F, 'x3', 32)                         # acc_{i-1}
    g.st256(dx, 'x3', 0, post=64)                # dx_i, and advance to slot i+1
    g.mul(dx, F, out=F)                          # F <- dx * F
    return g


def p2_walk(R, k='\\k'):
    g = RGen(reserved=R, tag='2', k=k)
    t_ptr = g.alloc()                            # x5 stays free in the body: the table pointer is parameter p
    accv = g.ld256('x3', 32)
    inv = g.mul(R, accv, keep_a=True)
    dx = g.ld256('x3', 0, post=-64)
    g.mul(dx, R, out=R)                          # R <- R * dx_i
    G.finish(g, 'x1', '\\p', inv, -64)
    g.release(t_ptr)
    return g


def macro(name, lines, params='k'):
    return [f'\t.macro {name} {params}'] + lines + ['\t.endm']


def count(lines):
    return sum(1 for s in lines if s.startswith('\t'))


def gen_pass1():
    """pass1(n, w, tp, sl): x0=n, x1=w, x2=tp, x3=sl. Slot i (64 bytes) = [dx_i | acc_{i-1}], slot n = [- | acc_{n-1}].
    acc_i = acc_{i-1} * dx_i is carried in registers F; slot 0 is not written (pass 2 uses R for walk 0)."""
    F = ['x5', 'x6', 'x7', 'x8']
    w = p1_walk(F)
    assert not any(r in s for s in w.out + w.rare for r in ('x26', 'x29'))   # x26: tp[i+1]; x29: groups
    L = RCHG + macro('P1WALK', w.out, 'k, p') + macro('P1RARE', w.rare) + G.header('pass1') + G.prologue() + G.const_c0()
    g = RGen(reserved=F, tag='1', k='w0')
    t = g.alloc()
    g.ins(f'ldr {t}, [x2], #8')
    g.ld256(t, 0, into=F)
    wx = g.ld256('x1', 0, post=64)
    g.sub(F, wx)                                 # F = dx_0 = acc_0 ; wx released by sub
    g.release(t)
    L += g.out
    group = []
    for k in range(0, 8, 2):                     # one LDP of tp[i], tp[i+1] per two walks
        group += ['\tldp x9, x26, [x2], #16', f'\tP1WALK {k}, x9', f'\tP1WALK {k + 1}, x26']
    L += ['\tadd x3, x3, #64', '\tsub x0, x0, #1', '\tlsr x29, x0, #3', '\tand x0, x0, #7', '\tcbz x0, 3f',
          '2:', '\tldr x9, [x2], #8', '\tP1WALK 8, x9', '\tsubs x0, x0, #1', '\tb.ne 2b',   # tail: (n - 1) mod 8
          '3:', '\tcbz x29, 4f',                                                         # no group of eight
          '1:'] + group + ['\tsubs x29, x29, #1', '\tb.ne 1b', '4:']
    g2 = G.Gen(reserved=F)
    g2.st256(F, 'x3', 32)                        # acc_{n-1} for the inversion
    L += g2.out
    L += G.epilogue()
    L += ['pass1_rare:'] + [f'\tP1RARE {k}' for k in range(9)] + g.rare
    return L, count(w.out), w.peak + len(F)


def gen_pass2():
    """pass2(n, w, tp, sl, Rinv): replaces W_i by W_i + T_i. Slot i = [dx_i | acc_{i-1}] as written by pass1.
    R (= 1/(dx_0 ... dx_i) at iteration i) is carried in registers across iterations."""
    R = ['x25', 'x26', 'x27', 'x28']
    w = p2_walk(R)
    assert not any(r in s for s in w.out + w.rare for r in ('x5', 'x0'))   # x5 and x0 hold the table pointers
    L = macro('P2WALK', w.out, 'k, p') + macro('P2RARE', w.rare) + G.header('pass2') + G.prologue()
    L += [f'\tldp {R[0]}, {R[1]}, [x4]', f'\tldp {R[2]}, {R[3]}, [x4, #16]']
    L += G.const_c0()
    group = []
    for k in range(0, 4, 2):                     # one LDP of tp[i-1], tp[i] per two walks; x2 = &tp[i+1]
        group += ['\tldp x0, x5, [x2, #-16]!', f'\tP2WALK {k}, x5', f'\tP2WALK {k + 1}, x0']
    L += ['\tsub x5, x0, #1',                    # i = n-1
          '\tadd x1, x1, x5, lsl #6',            # &w[n-1]
          '\tadd x1, x1, #32',                   # -> its y limbs
          '\tadd x2, x2, x0, lsl #3',            # &tp[n]
          '\tadd x3, x3, x5, lsl #6',            # &slot[n-1]
          '\tand x0, x5, #3', '\tcbz x0, 3f',
          '2:', '\tldr x5, [x2, #-8]!', '\tP2WALK 4, x5', '\tsubs x0, x0, #1', '\tb.ne 2b',  # tail: (n - 1) mod 4
          '3:', '\ttbz x3, #22, 4f',                                                     # no group of four
          '1:'] + group + ['\ttbnz x3, #22, 1b', '4:']
    g = RGen(reserved=R, tag='2', k='w0')
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2, #-8]')
    G.finish(g, 'x1', t_ptr, list(R), None)      # walk 0 uses R itself as its inverse
    g.release(t_ptr)
    L += g.out
    L += G.epilogue()
    L += ['pass2_rare:'] + [f'\tP2RARE {k}' for k in range(5)] + g.rare
    return L, count(w.out), max(w.peak + len(R), g.peak + len(R))


if __name__ == '__main__':
    p1, b1, k1 = gen_pass1()
    p2, b2, k2 = gen_pass2()
    sys.stderr.write(f'pass1 walk {b1} (group of 8: {8 * b1 + 6}) peak {k1}; '
                     f'pass2 walk {b2} (group of 4: {4 * b2 + 3}) peak {k2}\n')
    print('\n'.join(s if (s.split() or ['b'])[0][0] in 'bct' else s.replace(', ', ',') for s in p1 + p2))
