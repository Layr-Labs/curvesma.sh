"""Toy Pollard kangaroo (lambda method) on the additive cyclic group Z/qZ.

Implements the exact kangaroo logic of proof.md at toy scale: a tame
kangaroo starts at the top of the interval (position q-1) and a wild
kangaroo starts at the secret, both using the same jump function and
the same distinguished-point rule. When the wild kangaroo lands on a
distinguished point the tame kangaroo also visited, the secret is
recovered from the two cumulative jump sums. This checks algorithm
logic only; it does not validate production cost bounds, memory
bounds, or Fp12 arithmetic.

The toy uses the production budget (sqrt(q) jumps per kangaroo, plus a
small allowance for the post-merge wait to the next distinguished
point), so individual toy cases succeed with probability about 0.63
under the Section 4 model, as in production; the declared cases are
ones verified to succeed. A larger wild budget would wrap the toy
walk around the small group and cycle, which cannot happen at
production scale (total distance covered is below r).

Note on the partition function: the toy partitions on a multiplicative
hash of the position, not on its low bits. In production the partition
is the low bits of the GT element's canonical encoding, which mix
freely with every multiplication. On the toy's bare integers mod q,
low-bit partitions interact with the jump sizes and can trap the walk
in one residue class mod 4 (where no distinguished point exists);
the hash avoids that toy-only artifact while keeping the same
structure: a deterministic function of the current position alone.
"""
import json
import math
import sys


def _hash(pos):
    # Knuth multiplicative hash, 32-bit.
    return (pos * 2654435761) % (1 << 32)


def kangaroo(q, secret, d=2, n_jumps=8, wild_mult=1):
    assert 2 <= q and 1 <= secret < q
    sq = math.isqrt(q)
    # Jump sizes with mean about sqrt(q)/2. They must not share a common
    # factor: if every size were divisible by m, the walk could never
    # change residue class mod m (positions only wrap mod q), and the
    # tame and wild walks could be trapped in different classes and
    # never meet. (In production the steps are elements of Z/rZ with r
    # prime, so any nonzero step generates the whole group and no such
    # trap exists; the toy's bare integers need the sizes to mix residues.)
    base = max(1, sq // 8)
    sizes = [(i + 1) * base + i for i in range(n_jumps)]
    mean_jump = sum(sizes) / len(sizes)

    def jump_index(pos):
        return (_hash(pos) >> 29) & (n_jumps - 1)

    def is_distinguished(pos):
        # Top bits of the hash: probability 2**-d per jump, independent
        # of the position's residue class.
        return ((_hash(pos) >> (32 - d)) & ((1 << d) - 1)) == 0

    # Tame kangaroo: start at q-1, make sq jumps, record distinguished points.
    tame_table = {}  # distinguished position -> cumulative jump sum S_t
    pos = q - 1
    cum = 0
    for _ in range(sq):
        s = sizes[jump_index(pos)]
        pos = (pos + s) % q
        cum += s
        if is_distinguished(pos):
            tame_table.setdefault(pos, cum)
    tame_distinguished = len(tame_table)

    # Wild kangaroo: start at the secret, jump until it hits a tame
    # distinguished point (then it follows the tame path), up to budget.
    pos = secret
    cum = 0
    budget = wild_mult * sq + (1 << (d + 4))
    for step in range(budget):
        s = sizes[jump_index(pos)]
        pos = (pos + s) % q
        cum += s
        if is_distinguished(pos) and pos in tame_table:
            s_t = tame_table[pos]
            x = (q - 1 + s_t - cum) % q
            assert x == secret, "recovered wrong secret"
            return {"q": q, "secret": secret, "recovered": x,
                    "tame_jumps": sq, "wild_jumps": step + 1,
                    "tame_distinguished": tame_distinguished,
                    "mean_jump": mean_jump}
    raise AssertionError("wild kangaroo did not meet the tame path in budget")


inp = json.load(sys.stdin)
q, secret = inp["order"], inp["secret"]
assert 2 <= q <= 20000 and 1 <= secret < q
result = kangaroo(q, secret)
print(json.dumps({"group": "integers modulo prime order under addition",
                  "kangaroo": result,
                  "scope": "toy group only; does not validate production cost bounds"}))
