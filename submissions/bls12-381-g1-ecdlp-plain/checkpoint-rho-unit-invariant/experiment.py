"""Supporting checks. Not a production discrete-log run.

stdin: {"bits": int, "trials": int}

Curve checks run once. `trials` integer-group walks use a prime order n < 2^bits
with n ≡ 1 (mod 3), and the same unit matrices as the production coefficient update.
"""

from __future__ import annotations

import json
import sys

from reference import (
    apply_linear_unit,
    budget_holds,
    dlp_from_collision,
    endomorphism_matches_scalar,
    is_checkpoint,
    step_preserves_invariant,
    unit_images_match_multiplication,
    unit_multiplier,
)


def _is_prime(q: int) -> bool:
    if q < 2:
        return False
    if q % 2 == 0:
        return q == 2
    d = 3
    while d * d <= q:
        if q % d == 0:
            return False
        d += 2
    return True


def _prime_one_mod_three(bits: int) -> int:
    q = (1 << bits) - 1
    if q % 2 == 0:
        q -= 1
    while q > 2:
        if q % 3 == 1 and _is_prime(q):
            return q
        q -= 2
    raise ValueError("no prime")


def _lambda(n: int) -> int:
    for x in range(2, n):
        if (x * x + x + 1) % n == 0:
            return x
    raise ValueError("no lambda")


def _canon(s: int, n: int, lam: int) -> tuple[int, int]:
    if s % n == 0:
        return 0, 0
    best_u = 0
    best = s % n
    for u in range(6):
        img = (unit_multiplier(u, lam, n) * s) % n
        if img < best:
            best = img
            best_u = u
    return best, best_u


def solve_toy(n: int, k: int, lam: int) -> int | None:
    """Rho on Z/nZ. Both coefficient pairs use `apply_linear_unit` after the jump."""
    if k % n == 0:
        return 0
    raw_jumps = [(3 + 5 * i) % n or 1 for i in range(16)]
    jumps = []
    for j in raw_jumps:
        image, uj = _canon(j, n, lam)
        s0, s1 = apply_linear_unit(j, 0, uj, n)
        jumps.append((image, s0, s1))

    def comb(z0: int, z1: int) -> int:
        return (z0 + z1 * lam) % n

    # Start at T = [k]G, coefficients (1, 0) on T and (0, 0) on G, then canonicalise.
    point, u0 = _canon(k % n, n, lam)
    a0, a1 = apply_linear_unit(1, 0, u0, n)
    b0, b1 = apply_linear_unit(0, 0, u0, n)
    seen: dict[int, tuple[int, int, int, int]] = {}
    for _ in range(n * 4):
        if comb(a0, a1) * k % n != (point - comb(b0, b1)) % n:
            return None
        prev = seen.get(point)
        if prev is not None:
            got = dlp_from_collision(
                comb(prev[0], prev[1]),
                comb(prev[2], prev[3]),
                comb(a0, a1),
                comb(b0, b1),
                n,
            )
            if got is not None:
                return got
            image, s0, s1 = jumps[(point + 1) % len(jumps)]
            summed = (point + image) % n
            canon, u = _canon(summed, n, lam)
            a0, a1 = apply_linear_unit(a0, a1, u, n)
            b0, b1 = apply_linear_unit((b0 + s0) % n, (b1 + s1) % n, u, n)
            point = canon
            continue
        seen[point] = (a0, a1, b0, b1)
        image, s0, s1 = jumps[point % len(jumps)]
        summed = (point + image) % n
        canon, u = _canon(summed, n, lam)
        a0, a1 = apply_linear_unit(a0, a1, u, n)
        b0, b1 = apply_linear_unit((b0 + s0) % n, (b1 + s1) % n, u, n)
        point = canon
    return None


def toy_recover(bits: int, trials: int) -> dict:
    n = _prime_one_mod_three(bits)
    lam = _lambda(n)
    rows = []
    recovered = 0
    for t in range(trials):
        secret = 0 if t == 0 else 1 + ((0x5A * (t + 1)) % (n - 1))
        got = solve_toy(n, secret, lam)
        ok = got == secret
        rows.append({"secret": secret, "got": got, "ok": ok})
        if ok:
            recovered += 1
    return {
        "order": n,
        "lambda": lam,
        "maps": unit_images_match_multiplication(16, n, lam),
        "recovered": recovered,
        "trials": trials,
        "rows": rows,
    }


def main() -> None:
    data = json.load(sys.stdin)
    bits = data["bits"]
    trials = data["trials"]
    toy = toy_recover(bits, trials)
    report = {
        "unit_maps": unit_images_match_multiplication(),
        "endomorphism": endomorphism_matches_scalar(),
        "step_invariant": step_preserves_invariant(),
        "budget_holds": budget_holds(),
        "checkpoint_zero_window": is_checkpoint((1 << 117, 2)),
        "checkpoint_bit20": is_checkpoint((1 << 20, 2)),
        "toy": toy,
        "scope": (
            "curve checks are algebraic identities on the generator; "
            "toy rho is an integer group of order < 2^bits; "
            "neither measures production cost"
        ),
    }
    report["ok"] = bool(
        report["unit_maps"]
        and report["endomorphism"]
        and report["step_invariant"]
        and report["budget_holds"]
        and report["checkpoint_zero_window"]
        and not report["checkpoint_bit20"]
        and toy["maps"]
        and toy["recovered"] == trials
    )
    print(json.dumps(report))
    if not report["ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
