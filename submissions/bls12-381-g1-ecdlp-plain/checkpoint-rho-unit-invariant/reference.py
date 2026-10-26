"""Six-class checkpoint rho for plain BLS12-381 G1.

Walk on orbits {±P, ±λP, ±λ²P}. One β-multiplication builds the class
representative (β²x = -x - βx). The coefficient update multiplies the
post-jump scalar by that same unit, so S = [a]T + [b]P0 is preserved.

Cost integers are defined once in `budget()` and justified in proof.md.
"""

from __future__ import annotations

import math

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
GX = 0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB
GY = 0x08B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1

# bits [0, 20) select the jump; bits [20, 117) are the checkpoint window.
JUMP_BITS = 20
CHECKPOINT_LO = 20
CHECKPOINT_HI = 117
CHECKPOINT_BITS = CHECKPOINT_HI - CHECKPOINT_LO  # 97

# ln(5/3) < 5108256237659908 / 10^16. Birthday tail e^{-t^2/(2N)} ≤ 3/5.
LN_NUM = 5108256237659908
LN_DEN = 10**16

C = 2000  # reference G1 addition price
FP_MUL_INS = 383  # measured blst_fp_mul
FP_ADD_INS = 64  # 6-limb add/sub upper bound
FP_CMP_INS = 48
SCALAR4_INS = 40  # 4-limb add/sub mod r, with conditional correction
SELECT_INS = 24
LOOKUP_INS = 48
M_WALKS = 1 << 16
THETA_EXP = 97
ESCAPE_ALLOWANCE = 4096
EXCEPTIONAL_ADDS = 1 << 20
JUMPS = 1 << 20
ADDS_PER_SCALAR = 508  # 254 doubles + 254 adds, 255-bit scalar
SCALAR_MISC_INS = 64
VERIFY_INS = 200
MEM_INS_PER_RECORD = 256
RECORD_BYTES = 256
JUMP_BYTES = 192
CODE_BYTES = 1 << 20
# Factor 4 on stored checkpoints: room for a mild canonical-x bias.
CHECKPOINT_SLACK = 4

Point = tuple[int, int] | None


def _tonelli(n: int, p: int) -> int:
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


def _cube_root(mod: int) -> int:
    """One root of X^2+X+1. The other root is -1 - this."""
    sqrt_m3 = _tonelli((-3) % mod, mod)
    return ((-1 + sqrt_m3) * pow(2, mod - 2, mod)) % mod


def _add(p1: Point, p2: Point) -> Point:
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


def _neg(pt: Point) -> Point:
    if pt is None:
        return None
    x, y = pt
    return x, (P - y) % P


def scalar_mul(k: int, pt: Point) -> Point:
    k %= R
    acc: Point = None
    base = pt
    while k:
        if k & 1:
            acc = _add(acc, base)
        base = _add(base, base)
        k >>= 1
    return acc


def _bind_roots() -> tuple[int, int, Point]:
    """Pair β with λ so that (β·x(G), y(G)) = [λ]G. Returns that point too."""
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
    raise RuntimeError("no consistent (beta, lambda) pairing on the generator")


BETA, LAMBDA, LAMBDA_G = _bind_roots()
assert (BETA * BETA + BETA + 1) % P == 0
assert (LAMBDA * LAMBDA + LAMBDA + 1) % R == 0
assert pow(LAMBDA, 3, R) == 1 and LAMBDA != 1


def apply_linear_unit(a: int, b: int, unit: int, mod: int = R) -> tuple[int, int]:
    """Multiply a + b·λ by unit index 0..5 = id, λ, λ², -1, -λ, -λ².

    λ² = -1-λ is used as an integer identity, so every entry is in {-1,0,1}
    before reduction. The same indices are used by `apply_unit_point`.
    """
    if unit == 0:
        return a % mod, b % mod
    if unit == 1:  # λ·(a+bλ) = -b + (a-b)λ
        return (-b) % mod, (a - b) % mod
    if unit == 2:  # λ²·(a+bλ) = (b-a) + (-a)λ
        return (b - a) % mod, (-a) % mod
    if unit == 3:
        return (-a) % mod, (-b) % mod
    if unit == 4:  # -λ·(a+bλ) = b + (b-a)λ
        return b % mod, (b - a) % mod
    if unit == 5:  # -λ²·(a+bλ) = (a-b) + aλ
        return (a - b) % mod, a % mod
    raise ValueError(unit)


def unit_multiplier(unit: int, lam: int = LAMBDA, mod: int = R) -> int:
    scales = (1, lam, (lam * lam) % mod, mod - 1, (mod - lam) % mod, (mod - (lam * lam)) % mod)
    return scales[unit]


def unit_images_match_multiplication(trials: int = 32, mod: int = R, lam: int = LAMBDA) -> bool:
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


def beta2_x(x: int) -> int:
    """β²x = -x - βx. One multiplication by β, then a subtraction."""
    return (P - x - (BETA * x) % P) % P


def apply_unit_point(pt: Point, unit: int) -> Point:
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


def canonical_unit(pt: Point) -> int:
    """Least x in {x, βx, β²x}; on that x, the smaller of y and -y."""
    if pt is None:
        return 0
    x, y = pt
    xs = (x, (BETA * x) % P, beta2_x(x))
    best = min(xs)
    base = xs.index(best)
    if y > (P - y) % P:
        return base + 3
    return base


def on_curve(pt: Point) -> bool:
    if pt is None:
        return True
    x, y = pt
    return (y * y - (x * x * x + 4)) % P == 0


def central_coefficient_update(
    a0: int,
    a1: int,
    b0: int,
    b1: int,
    s0: int,
    s1: int,
    unit: int,
    mod: int = R,
) -> tuple[int, int, int, int]:
    """Coefficients of φ_u(S+J).

    S = [a0+a1·λ]T + [b0+b1·λ]P0, J = [s0+s1·λ]P0.
    Add the jump in the basis (1, λ), then multiply both coefficient pairs
    by the same unit that is applied to the point. Multiplying only the
    pre-jump pairs, or adding the jump after the unit, is the scalar of a
    different point.
    """
    na0, na1 = apply_linear_unit(a0, a1, unit, mod)
    nb0, nb1 = apply_linear_unit((b0 + s0) % mod, (b1 + s1) % mod, unit, mod)
    return na0, na1, nb0, nb1


def step_state(
    point: Point,
    a0: int,
    a1: int,
    b0: int,
    b1: int,
    jump: tuple[Point, int, int],
) -> tuple[Point, int, int, int, int]:
    jp, s0, s1 = jump
    summed = _add(point, jp)
    unit = canonical_unit(summed)
    na0, na1, nb0, nb1 = central_coefficient_update(a0, a1, b0, b1, s0, s1, unit)
    return apply_unit_point(summed, unit), na0, na1, nb0, nb1


def jump_from_scalar(s: int, base: Point = (GX, GY)) -> tuple[Point, int, int]:
    """Canonical table entry for [s]base, with the scalar kept in basis (1, λ)."""
    raw = scalar_mul(s % R, base)
    unit = canonical_unit(raw)
    s0, s1 = apply_linear_unit(s % R, 0, unit)
    return apply_unit_point(raw, unit), s0, s1


def is_checkpoint(pt: Point) -> bool:
    if pt is None:
        return False
    x, _y = pt
    window = (x >> CHECKPOINT_LO) & ((1 << CHECKPOINT_BITS) - 1)
    return window == 0


def dlp_from_collision(a1: int, b1: int, a2: int, b2: int, mod: int = R) -> int | None:
    da = (a1 - a2) % mod
    if da == 0:
        return None
    return ((b2 - b1) % mod) * pow(da, mod - 2, mod) % mod


def endomorphism_matches_scalar() -> bool:
    """φ(G)=[λ]G was bound at import. Check the other units against that point."""
    g = (GX, GY)
    if not on_curve(g) or LAMBDA_G != apply_unit_point(g, 1):
        return False
    # φ([λ]G) = [λ²]G, and β²x = -x-βx must agree with a second application of β.
    if apply_unit_point(g, 2) != apply_unit_point(LAMBDA_G, 1):
        return False
    if apply_unit_point(g, 3) != _neg(g):
        return False
    if apply_unit_point(g, 4) != _neg(LAMBDA_G):
        return False
    if apply_unit_point(g, 5) != _neg(apply_unit_point(g, 2)):
        return False
    return True


def combined_scalar(z0: int, z1: int, mod: int = R, lam: int = LAMBDA) -> int:
    return (z0 + z1 * lam) % mod


def step_preserves_invariant() -> bool:
    """Four jump-then-canonicalise steps. Point must stay [α]T + [β]P0."""
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


def fermat_mul_count() -> int:
    e = P - 2
    return e.bit_length() + e.bit_count() - 2


def collision_step_budget() -> int:
    """Least t with t^2 / (2N) ≥ LN_NUM/LN_DEN, N = (r-1)/6."""
    n_orbits = (R - 1) // 6
    need = (2 * n_orbits * LN_NUM + LN_DEN - 1) // LN_DEN
    t = math.isqrt(need)
    if t * t < need:
        t += 1
    return t


def budget() -> dict[str, int]:
    r_canon = FP_MUL_INS + 2 * FP_ADD_INS + 3 * FP_CMP_INS + FP_ADD_INS + SELECT_INS
    r_coeff = 4 * SCALAR4_INS
    r_base = r_canon + r_coeff + LOOKUP_INS
    surplus_muls = max(0, fermat_mul_count() - 3)
    surplus_ins = (surplus_muls * FP_MUL_INS + M_WALKS - 1) // M_WALKS
    r_step = r_base + surplus_ins
    t = collision_step_budget()
    overhead = M_WALKS * (1 << THETA_EXP)  # m/θ total steps
    h_walk = t + overhead + ESCAPE_ALLOWANCE + EXCEPTIONAL_ADDS
    h_pre = JUMPS * ADDS_PER_SCALAR
    r_pre = JUMPS * (r_canon + SCALAR_MISC_INS)
    h_ver = ADDS_PER_SCALAR
    r_ver = VERIFY_INS
    cps = CHECKPOINT_SLACK * ((h_walk + (1 << THETA_EXP) - 1) // (1 << THETA_EXP))
    r_mem = cps * MEM_INS_PER_RECORD
    h = h_walk + h_pre + h_ver
    r_ins = h_walk * r_step + r_pre + r_ver + r_mem
    mem = cps * RECORD_BYTES + JUMPS * JUMP_BYTES + CODE_BYTES
    return {
        "N": (R - 1) // 6,
        "t": t,
        "R_canon": r_canon,
        "R_step": r_step,
        "surplus_ins": surplus_ins,
        "fermat_muls": fermat_mul_count(),
        "overhead": overhead,
        "H_walk": h_walk,
        "H_pre": h_pre,
        "R_pre": r_pre,
        "H": h,
        "R": r_ins,
        "W_num": h * C + r_ins,
        "W_den": C,
        "checkpoints": cps,
        "memory_bytes": mem,
        "pre_num": h_pre * C + r_pre,
        "pre_den": C,
    }


def budget_holds() -> bool:
    b = budget()
    n = b["N"]
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
