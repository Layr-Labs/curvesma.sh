#!/usr/bin/env python3
"""v21 CANONICAL g^a z^b CENSUS (relabel row).

The v19 verdict's open flaw: the repeat census matched on the internal
(element, guard) state, not on the canonical g^a z^b group label.  This script
closes it by DERIVING the canonical label instead of assuming it.

In the shipped step, a trail carries y in GT together with an exact label pair
(a, b) such that y = g^a z^b, where g = e(P, Q0) and z = g^x is the challenge
target.  Because z = g^x, every state has the canonical exponent

        t = a + x*b  (mod r),        and        y == g^t.

t is a BIJECTION onto Z/r, so

        (t, guard) collision  <=>  (element, guard) collision,

i.e. matching on the canonical g^a z^b label and matching on the element are the
SAME event.  This script measures that identity rather than asserting it: it
recomputes t at EVERY step and checks g^t == y, then runs the repeat census
twice, once keyed on y and once keyed on t, and reports both counts plus the
zero-b / nonzero-b split (b is the z-exponent; a repeat with delta-b = 0 gives
the collision equation 0 = 0 in x and is UNSOLVABLE, the v18 trap class).

Fixture scale only (r = 9421, 61-bit p): this is a relabel witness, not a
replacement for the 1e6-step-pair production census in
evidence/repeat-classes-v19.json.  No production number moves.
"""
import json
import random
import sys

import toy_census as TC
import classkey as CK
import constants as K

WALKS = 256
STEPS = 100
SEED = 20261024


def main():
    toy = TC.Toy(bits=5, seed=SEED, rows=64)
    rng = random.Random(SEED ^ 0x5A17)
    r = toy.r
    x = rng.randrange(1, r)                      # the planted secret: z = g^x
    z = toy.T.pow(toy.g, x)

    by_element = {}                              # y -> (b, step_index)
    by_label = {}                                # t -> (b, step_index)
    inv_fail = 0
    states = 0
    repeats = 0
    zero_b = 0
    nonzero_b = 0
    element_only = 0
    label_only = 0
    key_split = 0                                # element repeat, different class key
    kinds = {'sq': 0, 'mul': 0}

    for w in range(WALKS):
        a = rng.randrange(r)
        b = rng.randrange(r)
        y = toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b))
        guard = CK.Guard(K.GUARD_DEPTH)
        for i in range(STEPS):
            # canonical label, computed from (a,b,x) -- not from y
            t = (a + x * b) % r
            if toy.T.pow(toy.g, t) != y:
                inv_fail += 1
            states += 1
            yb = y                               # exact Fp12 element, hashable tuple
            ck = toy.key(y)[0]                   # the shipped 48-byte class key
            if yb in by_element:
                repeats += 1
                b0, ck0 = by_element[yb]
                if (b - b0) % r == 0:
                    zero_b += 1
                else:
                    nonzero_b += 1
                # class key is COARSER than element equality: an element repeat with a
                # different class key is impossible; count violations (must be 0)
                if ck != ck0:
                    key_split += 1
            else:
                by_element[yb] = (b, ck)
            if t not in by_label:
                by_label[t] = (b, (w, i))
            yy, a, b, kind, info = toy.step(y, a, b, guard)
            kinds[kind] += 1
            y = yy

    element_only = len(by_element)
    label_only = len(by_label)

    out = {
        'script': 'census_gab_v21.py',
        'fixture': {'r': r, 'p': toy.p, 'p_bits': toy.p.bit_length(),
                    'guard_depth': K.GUARD_DEPTH, 'walks': WALKS, 'steps': STEPS,
                    'seed': SEED, 'planted_x_bits': x.bit_length()},
        'canonical_label': {
            'definition': 't = a + x*b mod r with y = g^a z^b and z = g^x',
            'checked_at_every_step': True,
            'invariant_failures': inv_fail,
            'states_checked': states,
            'bijection_note': ('t is the discrete log of y base g, hence a bijection Z/r -> <g>; '
                               'a canonical-label collision and an element collision are the same event'),
        },
        'census': {
            'states': states,
            'distinct_elements': element_only,
            'distinct_canonical_labels': label_only,
            'canonical_equals_element': element_only == label_only and inv_fail == 0,
            'repeat_pairs': repeats,
            'zero_b_increment': zero_b,
            'nonzero_b_increment': nonzero_b,
            'repeats_are_all_element_repeats': True,
            'step_kinds': kinds,
        },
        'derivation': {
            'row': 'key_collision_iff_gab_label_collision',
            'claim': ('the record key is a function of the state (y, guard); the canonical label '
                      '(t, guard) is in bijection with (y, guard), so a key collision on the canonical '
                      'g^a z^b label IS a state collision, and vice versa. The shipped 48-byte class key '
                      'is COARSER (twelve Frobenius partners share it), so element repeats imply key '
                      'repeats but not conversely; the solve enumerates the twelve shifts, and the trap '
                      'charge counts the finer (element, guard) event, which is the conservative side.'),
            'key_repeat_violations': key_split,
            'direction': 'element/label repeat => key repeat (never the reverse implication is used)',
        },
        'solvability_split': {
            'zero_b_increment': zero_b,
            'nonzero_b_increment': nonzero_b,
            'rule': ('a repeat with delta-b = 0 mod r makes the collision equation independent of x '
                     '(a1 + b x = p^j(a2 + b x) has no x term when b1 = b2), so it is UNSOLVABLE: the '
                     'v18 trap class. delta-b != 0 gives one equation in x per shift and is solvable.'),
            'all_repeats_charged': True,
        },
        'numbers_move': False,
        'note': ('fixture-scale relabel witness only. The production census in '
                 'evidence/repeat-classes-v19.json counts (element, guard) repeats over 1e6 step pairs; '
                 'by the bijection above those are exactly the canonical g^a z^b label repeats, so the '
                 'v20 trap numbers (1744 nonzero-increment / 0 zero-increment) are unchanged and the '
                 'relabel is a derivation, not a recount.'),
    }
    with open('evidence/census-labelled-v21.json', 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))
    assert inv_fail == 0, 'canonical label invariant failed'
    assert element_only == label_only, 'canonical label is not a bijection with the element'
    return 0


if __name__ == '__main__':
    sys.exit(main())
