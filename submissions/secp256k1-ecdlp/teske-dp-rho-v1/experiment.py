"""Toy of the Teske adding walk in proof.md.

Additive group Z/n. No elliptic addition and no inversion batch. Branch
multipliers, the mix that selects them, and the 3/2 birthday budget match
the production claim. Distinguished-point filtering is off for n < 2^16:
the production modulus is larger than the toy group, so every iterate is
stored and the linear relation is checked exactly. This does not price
secp256k1 field arithmetic.
"""

import json
import math
import sys
from fractions import Fraction

WALKS = 8
GOLD = 0x9E3779B97F4A7C15
MASK64 = (1 << 64) - 1
MUL1 = 0xBF58476D1CE4E5B9

# Fixed branch multipliers. Same integers as proof.md. Not a searched table.
A = (
    9315109294132905716,
    9476356238027737019,
    17445819973253803743,
    3342237045799953233,
    2812045804143842406,
    9122178642407602516,
    13926598530136396895,
    17056811102254473403,
    11323810722078648909,
    15698505963632274993,
    9578086386675406826,
    4470105836484457736,
    17565181521087248143,
    3097285154792831098,
    12002478175040212671,
    12486765383173650733,
)
B = (
    8546391839045156921,
    14157483109116380921,
    7783733875923022258,
    5661477030581511630,
    7704659425681353120,
    15728860691533551284,
    13562443690413809913,
    17757487276652662861,
    2212146067940045771,
    14634568361653561203,
    9217311395089492335,
    9124695148773980898,
    2502434592559100122,
    7616228408119649253,
    3985732427038351031,
    2215723844435459483,
)


def ln2_upper() -> Fraction:
    acc = Fraction(0)
    for k in range(1, 61):
        acc += Fraction(1, k << k)
    return acc + Fraction(1, 61 << 60)


def branch(residue: int) -> int:
    z = (residue * GOLD) & MASK64
    z ^= z >> 30
    z = (z * MUL1) & MASK64
    z ^= z >> 27
    return z >> 60


def budget(n: int, walks: int) -> int:
    need = (2 * n * ln2_upper()).__ceil__()
    t = math.isqrt(need)
    while t * (t - 1) < need:
        t += 1
    pad = (t * 3 + 1) // 2
    return ((pad + walks - 1) // walks) * walks


def _solve(n: int, secret: int, a: int, b: int, prev) -> bool:
    a2, b2 = prev
    db = (b - b2) % n
    if math.gcd(db, n) != 1:
        return False
    got = ((a2 - a) % n) * pow(db, -1, n) % n
    return got == secret


def recover(n: int, secret: int, walks: int) -> bool:
    total = budget(n, walks)
    per = total // walks
    seen = {}
    states = []
    for w in range(walks):
        a = (w + 1) % n
        b = 1 % n
        residue = (a + b * secret) % n
        states.append([a, b, residue])
        prev = seen.get(residue)
        if prev is not None:
            if _solve(n, secret, a, b, prev):
                return True
        else:
            seen[residue] = (a, b)
    for _ in range(per):
        for st in states:
            j = branch(st[2])
            st[0] = (st[0] + A[j]) % n
            st[1] = (st[1] + B[j]) % n
            st[2] = (st[0] + st[1] * secret) % n
            prev = seen.get(st[2])
            if prev is None:
                seen[st[2]] = (st[0], st[1])
                continue
            if _solve(n, secret, st[0], st[1], prev):
                return True
    return False


def main() -> None:
    n = json.load(sys.stdin)["order"]
    if not isinstance(n, int) or isinstance(n, bool):
        raise SystemExit("order must be an integer")
    walks = WALKS if n > WALKS + 2 else max(2, n // 3)
    hits = 0
    for secret in range(n):
        if recover(n, secret, walks):
            hits += 1
    print(
        json.dumps(
            {
                "group": "integers modulo order under addition",
                "order": n,
                "walks": walks,
                "budget": budget(n, walks),
                "checked": n,
                "recovered": hits,
                "all_recovered": hits == n,
                "scope": "toy adding walk only; does not validate production instruction costs or distinguished-point tails",
            }
        )
    )


if __name__ == "__main__":
    main()
