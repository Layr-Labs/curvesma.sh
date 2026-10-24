"""Toy of the order-6 automorphism walk in proof.md.

On Z/n (n prime, n = 1 mod 6) a value plays the role of an x-coordinate and
the six multipliers {1, -1, lam, -lam, lam^2, -lam^2} play the role of
negation plus the GLV endomorphism. Canonicalisation takes the smallest
image. The budget is t0, the least t with t(t+1) >= 2*N_orb*ln2 (rational
upper bound), with no 5/4 factor, charged as total steps across 8 walks.

This revision replaces per-step basis-pair coefficient tracking with the
unit frame alpha_t = U_t * A_t, beta_t = U_t * B_t (registers A, B already
fold the start coefficients): two residue adds per step, verified on every
step against the expanded scalar. Collisions are detected only at
distinguished points (canonical value with cv == 0 mod DP_MOD); the toy
modulus 2 is the stated scale analogue of the production x mod 2^20
distinguisher. On a distinguished-point hit the relation is solved from the
two frames in O(1), without any replay. The run reports recovery at the
charged budget t0 and at 2*t0: at this size the expected DP wait (2 steps)
is a visible fraction of t0 = 15, so the second sweep shows the mechanism
converging once the delay is amortised.

The oracle applies the true scalar map value -> (alpha + secret*beta) mod n,
so all recoveries are exact checks, not statistical guesses.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction

WALKS = 8
BRANCHES = 16
DP_MOD = 2


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def inv(a: int, n: int) -> int:
    return pow(a % n, -1, n)


def ln2_upper() -> Fraction:
    """ln 2 < sum_{k=1..60} 1/(k 2^k) + 1/(61 2^60)."""
    s = Fraction(0)
    for k in range(1, 61):
        s += Fraction(1, k << k)
    return s + Fraction(1, 61 << 60)


def half_budget(n_orb: int) -> int:
    need = (2 * n_orb * ln2_upper()).__ceil__()
    t = 1
    while t * (t + 1) < need:
        t += 1
    return t


def glv_lambda(n: int) -> int:
    disc = (-3) % n
    s = pow(disc, (n + 1) // 4, n) if n % 4 == 3 else None
    roots = []
    if s is not None and s * s % n == disc:
        roots = [s, n - s]
    else:
        for x in range(2, n):
            if x * x % n == disc:
                roots = [x, n - x]
                break
    lams = []
    for rt in roots:
        for cand in ((rt - 1) * inv(2, n) % n, (-rt - 1) * inv(2, n) % n):
            if (cand * cand + cand + 1) % n == 0:
                lams.append(cand % n)
    if not lams:
        raise ValueError("no sqrt(-3) mod n")
    return min(lams)


def multipliers(lam: int, n: int) -> list[int]:
    lam2 = (lam * lam) % n
    return [1, (n - 1) % n, lam, (n - lam) % n, lam2, (n - lam2) % n]


def inv_mult(kind: int, mults: list[int], n: int) -> int:
    """Index of the inverse of mults[kind]: u * u' == 1 (mod n)."""
    for k, m in enumerate(mults):
        if (mults[kind] * m) % n == 1:
            return k
    raise ValueError(kind)


def branch_index(cv: int) -> int:
    mixed = ((cv & ((1 << 64) - 1)) * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
    return mixed >> 60


def branches(n: int) -> list[tuple[int, int]]:
    return [((2 * j + 1) % n, (3 * j + 5) % n) for j in range(BRANCHES)]


def canon(value: int, mults: list[int], n: int) -> tuple[int, int]:
    best_v, best_k = n, 0
    for kind, m in enumerate(mults):
        img = (m * value) % n
        if img < best_v:
            best_v, best_k = img, kind
    return best_v, best_k


def recover(secret: int, n: int, lam: int, mults: list[int], br: list[tuple[int, int]], budget: int) -> tuple[bool, int, int]:
    """Run WALKS walks for `budget` total steps across all walks; store only
    distinguished points; solve collisions from unit frames in O(1).
    Returns (ok, fruitless, dp_stores)."""
    dp_table: dict[int, tuple[int, int, int, int]] = {}
    inv_kinds = [inv_mult(k, mults, n) for k in range(6)]
    states: list[list[int]] = []
    fruitless = 0
    dp_stores = 0

    # State [cv, u_kind, A, B, walk_id]: canonical scalar alpha = m_U * A,
    # beta = m_U * B. At start A = w+1 and B = 1 are the RAW start
    # coefficients; the canonicalisation unit folds them through U, so no
    # separate alpha_0 storage is needed.
    def advance(st: list[int]) -> bool:
        """One walk step. True iff a DP hit solves the relation."""
        nonlocal fruitless, dp_stores
        cv, uk_old = st[0], st[1]
        a, b = br[branch_index(cv)]
        raw = (cv + a + secret * b) % n
        ncv, kind = canon(raw, mults, n)
        u = mults[kind]
        ik = inv_kinds[uk_old]
        # Frame update uses U^-1 of the CURRENT U, before composing u into U:
        #   alpha_{t+1} = u*(alpha_t + a) = u*U*(A + U^-1 * a)
        A = (st[2] + mults[ik] * a) % n
        B = (st[3] + mults[ik] * b) % n
        prod = (u * mults[uk_old]) % n
        uk_new = mults.index(prod)
        # Invariant check: canonical scalar = m_U*A + secret*m_U*B.
        alpha = (mults[uk_new] * A) % n
        beta = (mults[uk_new] * B) % n
        if (alpha + secret * beta) % n != ncv:
            raise AssertionError("unit frame drifted from canonical scalar")
        st[0], st[1], st[2], st[3] = ncv, uk_new, A, B
        if ncv % DP_MOD != 0:
            return False
        dp_stores += 1
        prev = dp_table.get(ncv)
        if prev is not None and prev[0] != st[4]:
            a2 = (mults[prev[1]] * prev[2]) % n
            b2 = (mults[prev[1]] * prev[3]) % n
            den = (beta - b2) % n
            if den != 0 and ((a2 - alpha) * inv(den, n)) % n == secret:
                return True
            fruitless += 1
        elif prev is None:
            dp_table[ncv] = (st[4], uk_new, A, B)
        return False

    for w in range(WALKS):
        raw0 = (w + 1 + secret) % n
        cv0, k0 = canon(raw0, mults, n)
        states.append([cv0, k0, w + 1, 1, w])

    steps = 0
    while steps < budget:
        for st in states:
            if steps >= budget:
                break
            steps += 1
            if advance(st):
                return True, fruitless, dp_stores
    return False, fruitless, dp_stores


def run(n: int) -> dict:
    if n < 67 or n > 2003 or n % 6 != 1 or not is_prime(n):
        raise ValueError("order must be a prime ≡ 1 (mod 6) in [67, 2003]")
    lam = glv_lambda(n)
    mults = multipliers(lam, n)
    if len(set(mults)) != 6:
        raise AssertionError("multipliers not distinct")
    orbit_bad = 0
    for v in range(1, min(n, 64)):
        if len({(m * v) % n for m in mults}) != 6:
            orbit_bad += 1
    n_orb = 1 + (n - 1) // 6
    t0 = half_budget(n_orb)
    br = branches(n)
    results = {}
    for label, budget in (("t0", t0), ("2t0", 2 * t0)):
        hits = fruitless = stores = 0
        for secret in range(n):
            ok, fr, ds = recover(secret, n, lam, mults, br, budget)
            hits += int(ok)
            fruitless += fr
            stores += ds
        results[label] = {
            "budget_steps": budget,
            "recovered": hits,
            "fruitless_relations": fruitless,
            "dp_stores": stores,
            "success": hits / n,
        }
    return {
        "group": "integers modulo order under addition, order-6 multiplier action",
        "order": n,
        "orbits": n_orb,
        "half_budget_steps": t0,
        "walks": WALKS,
        "branches": BRANCHES,
        "dp_modulus": DP_MOD,
        "checked": n,
        "sweeps": results,
        "orbit_sample_failures": orbit_bad,
        "basis": "unit frame alpha=U*A, beta=U*B; table'd unit-add updates",
        "scope": "toy group only; collision rate, frame identity and DP detection, not production instruction costs",
    }


if __name__ == "__main__":
    print(json.dumps(run(json.load(sys.stdin)["order"])))
