"""Supporting checks. Not a production discrete-log run.

stdin: {"bits": int, "trials": int}

Curve checks run once. `trials` integer-group walks use a prime order n < 2^bits
with n ≡ 1 (mod 3), and the same unit matrices as the production coefficient update.
"""

import json
import math
import sys

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
GX = 0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB
GY = 0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1
CHECKPOINT_LO = 20
CHECKPOINT_BITS = 97
LN_NUM = 5108256237659908
LN_DEN = 10**16
C = 2000
FP_MUL_INS = 383
FP_ADD_INS = 64
FP_CMP_INS = 48
SCALAR4_INS = 40
SELECT_INS = 24
LOOKUP_INS = 48
M_WALKS = 1 << 16
THETA_EXP = 97
ESCAPE_ALLOWANCE = 4096
EXCEPTIONAL_ADDS = 1 << 20
JUMPS = 1 << 20
ADDS_PER_SCALAR = 508
SCALAR_MISC_INS = 64
VERIFY_INS = 200
MEM_INS_PER_RECORD = 256
RECORD_BYTES = 256
JUMP_BYTES = 192
CODE_BYTES = 1 << 20
CHECKPOINT_SLACK = 4


def _tonelli(n, p):
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        raise ValueError("not a residue")
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    while t != 1:
        i = 1
        t2 = (t * t) % p
        while t2 != 1:
            t2 = (t2 * t2) % p
            i += 1
            if i >= m:
                raise ValueError("tonelli failed")
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p
    return r


def _cube_root(mod):
    sqrt_m3 = _tonelli((-3) % mod, mod)
    return ((-1 + sqrt_m3) * pow(2, mod - 2, mod)) % mod


def _add(p1, p2):
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if x1 == x2 and y1 == y2:
        lam = (3 * x1 * x1) * pow((2 * y1) % P, P - 2, P) % P
    else:
        lam = ((y2 - y1) % P) * pow((x2 - x1) % P, P - 2, P) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return x3, y3


def _neg(pt):
    if pt is None:
        return None
    x, y = pt
    return x, (P - y) % P


def scalar_mul(k, pt):
    k %= R
    acc = None
    base = pt
    while k:
        if k & 1:
            acc = _add(acc, base)
        base = _add(base, base)
        k >>= 1
    return acc


def _bind_roots():
    beta = _cube_root(P)
    lam = _cube_root(R)
    g = (GX, GY)
    image = ((beta * GX) % P, GY)
    lam_g = scalar_mul(lam, g)
    if lam_g == image:
        return beta, lam, lam_g
    lam_alt = (-1 - lam) % R
    lam_alt_g = scalar_mul(lam_alt, g)
    if lam_alt_g == image:
        return beta, lam_alt, lam_alt_g
    beta_alt = (-1 - beta) % P
    image_alt = ((beta_alt * GX) % P, GY)
    if lam_g == image_alt:
        return beta_alt, lam, lam_g
    if lam_alt_g == image_alt:
        return beta_alt, lam_alt, lam_alt_g
    raise RuntimeError("no consistent (beta, lambda) pairing")


BETA, LAMBDA, LAMBDA_G = _bind_roots()
assert (BETA * BETA + BETA + 1) % P == 0
assert (LAMBDA * LAMBDA + LAMBDA + 1) % R == 0
assert pow(LAMBDA, 3, R) == 1 and LAMBDA != 1


def apply_linear_unit(a, b, unit, mod=R):
    if unit == 0:
        return a % mod, b % mod
    if unit == 1:
        return (-b) % mod, (a - b) % mod
    if unit == 2:
        return (b - a) % mod, (-a) % mod
    if unit == 3:
        return (-a) % mod, (-b) % mod
    if unit == 4:
        return b % mod, (b - a) % mod
    if unit == 5:
        return (a - b) % mod, a % mod
    raise ValueError(unit)


def unit_multiplier(unit, lam=LAMBDA, mod=R):
    scales = (1, lam, (lam * lam) % mod, mod - 1, (mod - lam) % mod, (mod - (lam * lam)) % mod)
    return scales[unit]


def unit_images_match_multiplication(trials=32, mod=R, lam=LAMBDA):
    seed = 0xC0FFEE
    for _ in range(trials):
        seed = (seed * 0x9E3779B97F4A7C15 + 1) & ((1 << 128) - 1)
        a = seed % mod
        b = (seed >> 17) % mod
        scalar = (a + b * lam) % mod
        for u in range(6):
            na, nb = apply_linear_unit(a, b, u, mod)
            if (na + nb * lam) % mod != (scalar * unit_multiplier(u, lam, mod)) % mod:
                return False
    return True


def beta2_x(x):
    return (P - x - (BETA * x) % P) % P


def apply_unit_point(pt, unit):
    if pt is None or unit == 0:
        return pt
    x, y = pt
    if unit == 1:
        return (BETA * x) % P, y
    if unit == 2:
        return beta2_x(x), y
    if unit == 3:
        return x, (P - y) % P
    if unit == 4:
        return (BETA * x) % P, (P - y) % P
    if unit == 5:
        return beta2_x(x), (P - y) % P
    raise ValueError(unit)


def canonical_unit(pt):
    if pt is None:
        return 0
    x, y = pt
    xs = (x, (BETA * x) % P, beta2_x(x))
    best = min(xs)
    base = xs.index(best)
    if y > (P - y) % P:
        return base + 3
    return base


def on_curve(pt):
    if pt is None:
        return True
    x, y = pt
    return (y * y - (x * x * x + 4)) % P == 0


def jump_from_scalar(s, base=(GX, GY)):
    raw = scalar_mul(s % R, base)
    unit = canonical_unit(raw)
    s0, s1 = apply_linear_unit(s % R, 0, unit)
    return apply_unit_point(raw, unit), s0, s1


def step_state(point, a0, a1, b0, b1, jump):
    jp, s0, s1 = jump
    summed = _add(point, jp)
    unit = canonical_unit(summed)
    na0, na1 = apply_linear_unit(a0, a1, unit)
    nb0, nb1 = apply_linear_unit((b0 + s0) % R, (b1 + s1) % R, unit)
    return apply_unit_point(summed, unit), na0, na1, nb0, nb1


def is_checkpoint(pt):
    if pt is None:
        return False
    x, _y = pt
    window = (x >> CHECKPOINT_LO) & ((1 << CHECKPOINT_BITS) - 1)
    return window == 0


def dlp_from_collision(a1, b1, a2, b2, mod=R):
    da = (a1 - a2) % mod
    if da == 0:
        return None
    return ((b2 - b1) % mod) * pow(da, mod - 2, mod) % mod


def endomorphism_matches_scalar():
    g = (GX, GY)
    if not on_curve(g) or LAMBDA_G != apply_unit_point(g, 1):
        return False
    if apply_unit_point(g, 2) != apply_unit_point(LAMBDA_G, 1):
        return False
    if apply_unit_point(g, 3) != _neg(g):
        return False
    if apply_unit_point(g, 4) != _neg(LAMBDA_G):
        return False
    if apply_unit_point(g, 5) != _neg(apply_unit_point(g, 2)):
        return False
    return True


def combined_scalar(z0, z1, mod=R, lam=LAMBDA):
    return (z0 + z1 * lam) % mod


def step_preserves_invariant():
    g = (GX, GY)
    k = 9
    a0, a1, b0, b1 = 4, 1, 6, 2
    alpha = combined_scalar(a0, a1)
    beta = combined_scalar(b0, b1)
    point = scalar_mul((alpha * k + beta) % R, g)
    if point is None:
        return False
    for s in (5, 8, 13, 21):
        jump = jump_from_scalar(s, g)
        point, a0, a1, b0, b1 = step_state(point, a0, a1, b0, b1, jump)
        alpha = combined_scalar(a0, a1)
        beta = combined_scalar(b0, b1)
        expect = scalar_mul((alpha * k + beta) % R, g)
        if point != expect or not on_curve(point):
            return False
    return True


def fermat_mul_count():
    e = P - 2
    return e.bit_length() + e.bit_count() - 2


def collision_step_budget():
    n_orbits = (R - 1) // 6
    need = (2 * n_orbits * LN_NUM + LN_DEN - 1) // LN_DEN
    t = math.isqrt(need)
    if t * t < need:
        t += 1
    return t


def budget():
    r_canon = FP_MUL_INS + 2 * FP_ADD_INS + 3 * FP_CMP_INS + FP_ADD_INS + SELECT_INS
    r_coeff = 4 * SCALAR4_INS
    r_base = r_canon + r_coeff + LOOKUP_INS
    surplus_muls = max(0, fermat_mul_count() - 3)
    surplus_ins = (surplus_muls * FP_MUL_INS + M_WALKS - 1) // M_WALKS
    r_step = r_base + surplus_ins
    t = collision_step_budget()
    overhead = M_WALKS * (1 << THETA_EXP)
    h_walk = t + overhead + ESCAPE_ALLOWANCE + EXCEPTIONAL_ADDS
    h_pre = JUMPS * ADDS_PER_SCALAR
    r_pre = JUMPS * (r_canon + SCALAR_MISC_INS)
    h_ver = ADDS_PER_SCALAR
    r_ver = VERIFY_INS
    cps = CHECKPOINT_SLACK * ((h_walk + (1 << THETA_EXP) - 1) // (1 << THETA_EXP))
    r_mem = cps * MEM_INS_PER_RECORD
    h = h_walk + h_pre + h_ver
    r_ins = h_walk * r_step + r_pre + r_ver + r_mem
    return {"t": t, "R_step": r_step, "W_num": h * C + r_ins, "W_den": C,
            "H": h, "R": r_ins}


def budget_holds():
    b = budget()
    n = (R - 1) // 6
    t = b["t"]
    if t * t * LN_DEN < 2 * n * LN_NUM:
        return False
    if b["W_num"] != b["H"] * C + b["R"]:
        return False
    if b["R_step"] != 955:
        return False
    if b["W_den"] != 2000:
        return False
    return True


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
