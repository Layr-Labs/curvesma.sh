"""Emit the two counted batch loops unrolled (pass 1 eight times, pass 2 four times).

usage: python3 gen_unrolled.py elf|macho > secp_loops.s     (needs 2d25756c's kernels/gen_loops.py on the path)

The field emitters, the per-walk bodies and the register allocation are those of gen_loops.py; only the loop
control changes.  Each walk body is emitted once, as an assembler macro (P1WALK, P2WALK), and invoked once
in a tail loop and U times in the main loop.  The tail loop (SUBS, B.NE) runs the (n - 1) mod U walks first.
- pass 1 runs walks 1..n-1 forward; its main loop runs (n - 1) div 8 groups of eight walks, counted in x29
  (SUBS, B.NE: 2 per group).
- pass 2 runs walks n-1..1 backward in groups of four; its main loop ends when the slot pointer x3 reaches
  slot 0, which is the slot array's only address with bit 22 clear (TBNZ x3, #22: 1 per group).  That end does
  not depend on n.
Batches have n = 1 (mod 8) walks, so the tail loops run only in rounds where walks double.
"""
import sys

import gen_loops as G


def body(emit):
    g = emit()
    return g.out, g.peak


def p1_walk(F):
    g = G.Gen(reserved=F)
    t = g.alloc()
    g.ins(f'ldr {t}, [x2], #8')
    tx = g.ld256(t, 0)
    wx = g.ld256('x1', 0, post=64)
    dx = g.sub(tx, wx)
    g.release(t)
    g.st256(F, 'x3', 32)               # acc_{i-1}
    g.st256(dx, 'x3', 0, post=64)      # dx_i, and advance to slot i+1
    g.mul(dx, F, out=F)                # F <- dx * F
    return g


def p2_walk(R):
    g = G.Gen(reserved=R)
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2], #-8')
    accv = g.ld256('x3', 32)
    inv = g.mul(R, accv, keep_a=True)
    dx = g.ld256('x3', 0, post=-64)
    g.mul(dx, R, out=R)                # R <- R * dx_i
    G.finish(g, 'x1', t_ptr, inv, -64)
    g.release(t_ptr)
    return g


def macro(name, lines):
    return [f'\t.macro {name}'] + lines + ['\t.endm']


def gen_pass1():
    """pass1(n, w, tp, sl): x0=n, x1=w, x2=tp, x3=sl. Slot i (64 bytes) = [dx_i | acc_{i-1}], slot n = [- | acc_{n-1}].
    acc_i = acc_{i-1} * dx_i is carried in registers F; slot 0 is not written (pass 2 uses R for walk 0)."""
    F = ['x5', 'x6', 'x7', 'x8']
    w = p1_walk(F)
    L = macro('P1WALK', w.out) + G.header('pass1') + G.prologue() + G.const_c0()
    g = G.Gen(reserved=F)
    t = g.alloc()
    g.ins(f'ldr {t}, [x2], #8')
    g.ld256(t, 0, into=F)
    wx = g.ld256('x1', 0, post=64)
    g.sub(F, wx)                       # F = dx_0 = acc_0 ; wx released by sub
    g.release(t)
    L += g.out
    assert not any('x29' in s for s in w.out)                       # x29 counts the groups
    L += ['\tadd x3, x3, #64', '\tsub x0, x0, #1', '\tlsr x29, x0, #3', '\tand x0, x0, #7', '\tcbz x0, 3f',
          '2:', '\tP1WALK', '\tsubs x0, x0, #1', '\tb.ne 2b',        # tail: (n - 1) mod 8 walks
          '3:', '\tcbz x29, 4f',                                     # no group of eight
          '1:'] + ['\tP1WALK'] * 8 + ['\tsubs x29, x29, #1', '\tb.ne 1b', '4:']
    g2 = G.Gen(reserved=F)
    g2.st256(F, 'x3', 32)              # acc_{n-1} for the inversion
    L += g2.out
    L += G.epilogue()
    return L, len(w.out), w.peak + len(F)


def gen_pass2():
    """pass2(n, w, tp, sl, Rinv): replaces W_i by W_i + T_i. Slot i = [dx_i | acc_{i-1}] as written by pass1.
    R (= 1/(dx_0 ... dx_i) at iteration i) is carried in registers across iterations."""
    R = ['x25', 'x26', 'x27', 'x28']
    w = p2_walk(R)
    L = macro('P2WALK', w.out) + G.header('pass2') + G.prologue()
    L += [f'\tldp {R[0]}, {R[1]}, [x4]', f'\tldp {R[2]}, {R[3]}, [x4, #16]']
    L += G.const_c0()
    L += ['\tsub x5, x0, #1',                   # i = n-1
          '\tadd x1, x1, x5, lsl #6',           # &w[n-1]
          '\tadd x1, x1, #32',                  # -> its y limbs
          '\tadd x2, x2, x5, lsl #3',           # &tp[n-1]
          '\tadd x3, x3, x5, lsl #6',           # &slot[n-1]
          '\tand x0, x5, #3', '\tcbz x0, 3f',
          '2:', '\tP2WALK', '\tsubs x0, x0, #1', '\tb.ne 2b',        # tail: (n - 1) mod 4 walks
          '3:', '\ttbz x3, #22, 4f',                                  # no group of four left
          '1:'] + ['\tP2WALK'] * 4 + ['\ttbnz x3, #22, 1b', '4:']
    g = G.Gen(reserved=R)
    t_ptr = g.alloc()
    g.ins(f'ldr {t_ptr}, [x2]')
    G.finish(g, 'x1', t_ptr, list(R), None)    # walk 0 uses R itself as its inverse
    g.release(t_ptr)
    L += g.out
    L += G.epilogue()
    return L, len(w.out), max(w.peak + len(R), g.peak + len(R))


if __name__ == '__main__':
    p1, b1, k1 = gen_pass1()
    p2, b2, k2 = gen_pass2()
    sys.stderr.write(f'pass1 walk {b1} (group of 8: {8 * b1 + 2}) peak {k1}; '
                     f'pass2 walk {b2} (group of 4: {4 * b2 + 1}) peak {k2}\n')
    print('\n'.join(p1 + p2))
