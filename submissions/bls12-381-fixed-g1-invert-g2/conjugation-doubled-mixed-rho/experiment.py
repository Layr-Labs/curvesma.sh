"""Toy check of the collision algebra used in proof.md.

The walk runs on the additive group Z_r (the exponent group of a prime-order
group) with exact affine bookkeeping y = z^u * g^v, a hash-selected step from a
cyclotomic-squaring analogue (doubling), a Frobenius analogue (a fixed
multiplier) and a table addition on 1/64 of steps. Two variants are run: the
attack's walk, and a power-map-only walk with the additions removed.

For each variant the script reports how many matching (generated, stored)
element pairs were seen, how many matched through each of the two equations of
proof.md section 4, and whether the true discrete logarithm appears among the
solved values. Every genuine match must solve to the same logarithm, so a
correct run reports recovered=true and few distinct solutions.

Toy scope only: this establishes no production cost, memory or success bound.
"""
import json
import sys

cfg = json.load(sys.stdin)
r = int(cfg["order"])
chains = int(cfg.get("chains", 16))
steps = int(cfg.get("steps", 400))
if r < 5:
    raise SystemExit("order must be at least 5")
if chains < 2:
    raise SystemExit("chains must be at least 2")

A = 0x9E3779B97F4A7C15
B = 0xBF58476D1CE4E5B9
M64 = (1 << 64) - 1
TABLE = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59)
MULT = 3  # stands in for the Frobenius multiplier p


def mix(e, chain):
    x = (e * A + chain) & M64
    x ^= x >> 31
    x = (x * B) & M64
    x ^= x >> 29
    return x


def run(a, mixed):
    seen = {}
    solved = set()
    pairs = 0
    type1 = 0
    type2 = 0
    for chain in range(chains):
        u = chain + 1
        v = 1 + (chain * 7) % 97
        e = (u * a + v) % r
        for _ in range(steps):
            h = mix(e, chain)
            if mixed and h % 64 == 0:
                ci = TABLE[(h >> 6) % 16]
                e = (e + ci) % r
                v = (v + ci) % r
            elif (h >> 20) & 1:
                e = (2 * e) % r
                u = (2 * u) % r
                v = (2 * v) % r
            else:
                e = (MULT * e) % r
                u = (MULT * u) % r
                v = (MULT * v) % r
            cands = []
            for uj, vj in seen.get(e, ()):
                if (uj, vj) != (u, v):
                    cands.append((1, uj, vj))
            if e % r:
                for uj, vj in seen.get((-e) % r, ()):
                    if (uj, vj) != (u, v):
                        cands.append((2, uj, vj))
            for kind, uj, vj in cands:
                pairs += 1
                if kind == 1:
                    den = (u - uj) % r
                    if den:
                        type1 += 1
                        solved.add(((vj - v) % r) * pow(den, -1, r) % r)
                else:
                    den = (u + uj) % r
                    if den:
                        type2 += 1
                        solved.add((-(v + vj)) % r * pow(den, -1, r) % r)
            seen.setdefault(e, []).append((u, v))
    return {
        "element_pairs": pairs,
        "type1_solved": type1,
        "type2_solved": type2,
        "distinct_solutions": len(solved),
        "recovered": a in solved,
    }


a = 1234567 % r
if a == 0:
    a = 1

print(json.dumps({
    "group": "exponents modulo prime order (Z_r)",
    "order": r,
    "chains": chains,
    "steps_per_chain": steps,
    "attack_walk": run(a, True),
    "power_map_only_walk": run(a, False),
    "note": ("type1 is the same-element equation and type2 the conjugate "
             "equation of proof.md section 4. All genuine matches must solve "
             "to the one true logarithm, so distinct_solutions stays small "
             "while element_pairs grows; recovered must be true for both "
             "variants because the algebra is identical."),
    "scope": "toy group only; does not validate production cost bounds",
}))
