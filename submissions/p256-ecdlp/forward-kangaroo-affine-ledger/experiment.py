"""Exact small-group replica of the production forward-kangaroo update.

The production walk of proof.md keeps an affine point and an explicit
non-negative scalar offset, and every update is

    point  <- point  + jump_point[i]
    offset <- offset + jump_size[i]        (mod r)

where i is a deterministic function of the canonical affine abscissa.  No
negation map, no canonicalisation and no signed offset appear anywhere, so the
identical update can be reproduced on the additive group Z_q.

This module checks, exhaustively over every target scalar, that

  * every tame/wild point collision satisfies  d_tame - d_wild = x (mod q), and
  * every distinguished-point match therefore recovers x exactly.

It is a correctness replica at toy scope.  It makes no claim about production
collision times, instruction counts or memory bounds.
"""
import json
import math
import sys

K = 8
DP_MODULUS = 4
BUDGET_FACTOR = 64


def jump_sizes(order):
    """8 jumps with mean about sqrt(order)/2 and maximum about 2*mean."""
    unit = max(1, int(round(math.isqrt(order) / 2 / 4.5)))
    return [(i + 1) * unit for i in range(K)]


def step(pos, off, order, jumps):
    j = jumps[pos % K]
    return (pos + j) % order, (off + j) % order


def is_dp(pos):
    return pos % DP_MODULUS == 0


def trail(start, off0, order, jumps, steps):
    pts = []
    offs = []
    pos, off = start, off0
    for _ in range(steps):
        pts.append(pos)
        offs.append(off)
        pos, off = step(pos, off, order, jumps)
    return pts, offs


def run(order):
    jumps = jump_sizes(order)
    steps = BUDGET_FACTOR * math.isqrt(order) + 128

    tame_pts, tame_offs = trail(order - 1, order - 1, order, jumps, steps)
    by_point = {}
    dp_table = {}
    for i, p in enumerate(tame_pts):
        by_point.setdefault(p, []).append(i)
        if is_dp(p):
            dp_table.setdefault(p, []).append(i)

    merge_pairs = 0
    scalars_with_merge = 0
    recovered = 0
    unmatched = 0
    violations = 0

    for x in range(order):
        wild_pts, wild_offs = trail(x, 0, order, jumps, steps)
        hit = False
        found = False
        for m, p in enumerate(wild_pts):
            idxs = by_point.get(p)
            if idxs is not None:
                hit = True
                for i in idxs:
                    merge_pairs += 1
                    if (tame_offs[i] - wild_offs[m]) % order != x:
                        violations += 1
            if not found and is_dp(p):
                idxs = dp_table.get(p)
                if idxs is not None:
                    found = True
                    if (tame_offs[idxs[0]] - wild_offs[m]) % order == x:
                        recovered += 1
                    else:
                        violations += 1
        if hit:
            scalars_with_merge += 1
        if not found:
            unmatched += 1

    return {
        "group": "integers modulo order under addition",
        "order": order,
        "jump_sizes": jumps,
        "steps_per_walker": steps,
        "scalars": order,
        "merge_pairs_checked": merge_pairs,
        "scalars_with_merge": scalars_with_merge,
        "scalars_recovered_at_distinguished_point": recovered,
        "scalars_without_distinguished_match": unmatched,
        "invariant_violations": violations,
        "recovery_formula": "x = (d_tame - d_wild) mod order",
        "scope": "toy group only; does not validate production cost, memory or collision-time bounds",
    }


if __name__ == "__main__":
    print(json.dumps(run(json.load(sys.stdin)["order"]), indent=2))
