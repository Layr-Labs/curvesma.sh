"""Toy check of the fused cyclotomic schedule and the score inequality."""

import json
import sys


P = 43  # 43 ≡ 3 (mod 4), so u^2 + 1 is irreducible


def fp_add(a, b):
    return (a + b) % P


def fp_sub(a, b):
    return (a - b) % P


def fp_mul(a, b):
    return (a * b) % P


def f2_add(a, b):
    return (fp_add(a[0], b[0]), fp_add(a[1], b[1]))


def f2_sub(a, b):
    return (fp_sub(a[0], b[0]), fp_sub(a[1], b[1]))


def f2_mul(a, b):
    # (a0 + a1 u)(b0 + b1 u), u^2 = -1
    return (
        fp_sub(fp_mul(a[0], b[0]), fp_mul(a[1], b[1])),
        fp_add(fp_mul(a[0], b[1]), fp_mul(a[1], b[0])),
    )


def f2_sqr(a):
    return f2_mul(a, a)


def mul_by_u_plus_1(a):
    # (a0 + a1 u)(1 + u) = (a0 - a1) + (a0 + a1) u
    return (fp_sub(a[0], a[1]), fp_add(a[0], a[1]))


def sqr_fp4(a0, a1):
    """Fused norm: ret0 = a0^2 + (u+1) a1^2, ret1 = 2 a0 a1."""
    s0 = f2_sqr(a0)
    s1 = f2_sqr(a1)
    cross = f2_sub(f2_sub(f2_sqr(f2_add(a0, a1)), s0), s1)
    return (f2_add(s0, mul_by_u_plus_1(s1)), cross)


def cyc_full(a):
    """blst cyclotomic_sqr_fp12 on six Fp2 coefficients [a00,a01,a02,a10,a11,a12]."""
    a00, a01, a02, a10, a11, a12 = a
    t0 = sqr_fp4(a00, a11)
    t1 = sqr_fp4(a10, a02)
    t2 = sqr_fp4(a01, a12)

    def mix_sub(t, old):
        x = f2_sub(t, old)
        x = f2_add(x, x)
        return f2_add(x, t)

    def mix_add(t, old):
        x = f2_add(t, old)
        x = f2_add(x, x)
        return f2_add(x, t)

    r00 = mix_sub(t0[0], a00)
    r01 = mix_sub(t1[0], a01)
    r02 = mix_sub(t2[0], a02)
    twisted = mul_by_u_plus_1(t2[1])
    r10 = mix_add(twisted, a10)
    r11 = mix_add(t0[1], a11)
    r12 = mix_add(t1[1], a12)
    return (r00, r01, r02, r10, r11, r12)


def kept(elem):
    _r00, r01, r02, r10, _r11, r12 = elem
    return (r01, r02, r10, r12)


def main():
    data = json.load(sys.stdin)
    trials = data["trials"]
    seed = data["seed"]
    state = seed & 0xFFFFFFFF
    mismatches = 0
    for _ in range(trials):
        coeffs = []
        for _k in range(12):
            state = (1664525 * state + 1013904223) & 0xFFFFFFFF
            coeffs.append(state % P)
        elem = tuple((coeffs[2 * i], coeffs[2 * i + 1]) for i in range(6))
        # The four kept outputs are functions of t1 and t2 only.
        a00, a01, a02, a10, a11, a12 = elem
        t1 = sqr_fp4(a10, a02)
        t2 = sqr_fp4(a01, a12)

        def mix_sub(t, old):
            x = f2_sub(t, old)
            x = f2_add(x, x)
            return f2_add(x, t)

        def mix_add(t, old):
            x = f2_add(t, old)
            x = f2_add(x, x)
            return f2_add(x, t)

        partial = (
            mix_sub(t1[0], a01),
            mix_sub(t2[0], a02),
            mix_add(mul_by_u_plus_1(t2[1]), a10),
            mix_add(t1[1], a12),
        )
        if partial != kept(cyc_full(elem)):
            mismatches += 1
    frontier = 31447421 / 250000
    claimed = 25157 / 200
    ok = mismatches == 0 and claimed < frontier
    print(
        json.dumps(
            {
                "mismatches": mismatches,
                "claimed": claimed,
                "frontier": frontier,
                "below_frontier": claimed < frontier,
                "ok": ok,
            }
        )
    )
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
