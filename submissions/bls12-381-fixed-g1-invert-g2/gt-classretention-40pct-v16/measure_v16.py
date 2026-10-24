#!/usr/bin/env python3
"""v16 measurements: the CORRECTED retention event, class-stream agreement vs D.

The v15 verdict (408a011b) killed the v15 predicate: retention accepted a first
joint-state meeting anywhere in 1..2D, never inside D, and its 861 "successes"
were all at lags 13..24.  The deeper error is that the joint state (element,
ring tuple) is the WRONG event: the collision solve consumes a class match.

  y1 = g^{a1} z^{b1}, y2 = g^{a2} z^{b2}, class(y1) = class(y2)
  => y2 = sigma^j(y1) for some j, and x = log_g z is recovered from
     x (b2 - p^j b1) = p^j a1 - a2, tried for j = 0..11 and VERIFIED.

A class match is a function of the ELEMENT ALONE: the class key is
min{keys12(y)}, and the ring does not enter it.  So a class match at report time
is instantaneous and cannot be destroyed by a later guard rewrite.  What the
ring CAN destroy is path-sharing after a merge, which changes the RATE of later
matches, not the usability of the match already reported.

This script therefore measures, at fixture scale, both quantities:

  A. class_stream: two trails that start on the SAME class (same element, or
     Frobenius partners) with INDEPENDENT ring histories; the class keys are
     compared k steps after the merge.  This is the durable-merge curve; it is
     reported, and it is NOT what the ledger multiplies by.

  B. class_match_rate: independent random starts; first lag at which the two
     class keys agree, and the fraction of pairs that agree at least once in a
     window.  This is the event the birthday model counts.

Depth sweep D in {12, 24, 48, 128} at R = 64 and R = 1024, with CP-95 bounds.
Output: evidence/class-retention.json
"""
import json
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import toy_census as TC
from cp95 import cp95_lower

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
DEPTHS = (12, 24, 48, 128)


def warm_ring(toy, rng, depth):
    """A ring drawn from a real walk of the same depth (independent history)."""
    y = toy.rand_gt()
    g = CK.Guard(depth)
    for _ in range(depth):
        y, _a, _b, _k, _i = toy.step(y, 1, 0, g)
    return tuple(g.slots)


def curve_one(toy, depth, pairs, steps, seed):
    """Class-stream agreement k steps after a merge at s=0."""
    rng = random.Random(seed)
    agree = [0] * (steps + 1)
    elem_agree = [0] * (steps + 1)
    for _ in range(pairs):
        y0 = toy.rand_gt()
        gA, gB = CK.Guard(depth), CK.Guard(depth)
        gA.slots = list(warm_ring(toy, rng, depth))
        gB.slots = list(warm_ring(toy, rng, depth))
        yA, yB = y0, y0
        kA0 = CK.classkey(toy.T, yA)[0]
        kB0 = CK.classkey(toy.T, yB)[0]
        agree[0] += 1 if kA0 == kB0 else 0
        elem_agree[0] += 1 if yA == yB else 0
        for s in range(1, steps + 1):
            yA, _a, _b, _c, _d = toy.step(yA, 1, 0, gA)
            yB, _e, _f, _g, _h = toy.step(yB, 1, 0, gB)
            agree[s] += 1 if CK.classkey(toy.T, yA)[0] == CK.classkey(toy.T, yB)[0] else 0
            elem_agree[s] += 1 if yA == yB else 0
    return {
        'pairs': pairs, 'steps': steps, 'guard_depth': depth, 'rows': len(toy.us),
        'class_agree_frac': [a / pairs for a in agree],
        'class_agree_cp95_lower': [float(cp95_lower(a, pairs)) for a in agree],
        'elem_agree_frac': [a / pairs for a in elem_agree],
        'class_agree_at_1': agree[1] / pairs,
        'class_agree_at_min_depth': agree[min(depth, steps)] / pairs,
        'class_agree_at_min_depth_cp95_lower': float(cp95_lower(agree[min(depth, steps)], pairs)),
        'half_life_steps': next((s for s in range(1, steps + 1)
                                 if agree[s] * 2 < pairs), None),
        'note': "both trails start at the SAME element with independent rings of "
                "the same depth; class keys (min of the twelve Frobenius keys) are "
                "a function of the element alone, so a merge is usable at s=0 and "
                "the curve measures only whether the two trails keep sharing a path.",
    }


def match_rate(toy, depth, pairs, window, seed):
    """Independent starts: fraction of pairs whose class keys agree within window."""
    rng = random.Random(seed)
    hit = 0
    first = {}
    for _ in range(pairs):
        yA, yB = toy.rand_gt(), toy.rand_gt()
        gA, gB = CK.Guard(depth), CK.Guard(depth)
        lag = None
        for s in range(1, window + 1):
            yA, _a, _b, _c, _d = toy.step(yA, 1, 0, gA)
            yB, _e, _f, _g, _h = toy.step(yB, 1, 0, gB)
            if CK.classkey(toy.T, yA)[0] == CK.classkey(toy.T, yB)[0]:
                lag = s
                break
        if lag is not None:
            hit += 1
            first[str(lag)] = first.get(str(lag), 0) + 1
    return {
        'pairs': pairs, 'window': window, 'guard_depth': depth, 'rows': len(toy.us),
        'pairs_matching': hit, 'rate': hit / pairs,
        'cp95_lower': float(cp95_lower(hit, pairs)),
        'first_match_lag_hist': dict(sorted(first.items(), key=lambda kv: int(kv[0]))),
        'note': "independent uniform starts, no planted merge: the birthday event "
                "the success model counts, at fixture scale (N = r-1 classes, "
                "S = window steps).",
    }


def main():
    os.makedirs(EV, exist_ok=True)
    out = {'depths': list(DEPTHS), 'fixture': {}, 'class_stream': {}, 'class_match': {}}
    # fixture sizes: N = r-1 classes; ring alphabet = 12 R columns
    for rows in (64, 1024):
        toy = TC.Toy(rows=rows)
        out['fixture'][str(rows)] = {'rows': rows, 'r': toy.r, 'N_classes': toy.r - 1}
        cs, cm = {}, {}
        for D in DEPTHS:
            steps = min(4 * D + 8, 256)
            pairs = 160 if rows == 64 else 120
            cs[str(D)] = curve_one(toy, D, pairs, steps, seed=1000 + D + rows)
            print("R=%4d D=%3d class_stream k=1 %.4f  k=D %.4f  half-life %s"
                  % (rows, D, cs[str(D)]['class_agree_at_1'],
                     cs[str(D)]['class_agree_at_min_depth'],
                     cs[str(D)]['half_life_steps']), flush=True)
            if D == 12:
                cm[str(D)] = match_rate(toy, D, pairs * 2, 256, seed=9000 + rows)
                print("      match_rate %.4f (%d/%d pairs)"
                      % (cm[str(D)]['rate'], cm[str(D)]['pairs_matching'],
                         cm[str(D)]['pairs']), flush=True)
        out['class_stream'][str(rows)] = cs
        out['class_match'][str(rows)] = cm
    with open(os.path.join(EV, 'class-retention.json'), 'w') as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("wrote evidence/class-retention.json", flush=True)


if __name__ == '__main__':
    main()
