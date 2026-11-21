"""Evidence experiment for the automorphism-aware Pollard rho submission.

Part A verifies, at full production secp256k1 parameters and with pure
standard-library arithmetic, the exact algebraic structure the attack uses:
the order-6 automorphism group (x,y)->(zeta^i x, +/-y), its multiplier lam
with [lam]P = zeta(P), and the canonicalization-with-representation
bookkeeping on random points. No discrete logarithm is solved at production
size here.

Part B runs the complete attack end to end on a small j=0 toy curve
(p=1048609, y^2=x^3+29, prime order n=1047379, p=n=1 mod 3) and reports the
measured single-shot success rate at the same budget shape used in the claim,
against the ideal birthday-tail prediction. This is evidence at toy scale
only; it does not validate production cost bounds.

Stdlib only; deterministic seeds; reads JSON input on stdin.
"""
import json
import math
import random
import sys

# ---------------------------------------------------------------------------
# generic helpers
# ---------------------------------------------------------------------------


def is_prime_trial(v):
    if v < 2:
        return False
    for d in range(2, math.isqrt(v) + 1):
        if v % d == 0:
            return False
    return True


def proot3(mod):
    """Primitive cube root of unity modulo prime mod with mod % 3 == 1."""
    for g in range(2, 64):
        v = pow(g, (mod - 1) // 3, mod)
        if v != 1:
            assert pow(v, 3, mod) == 1
            return v
    raise AssertionError("no primitive cube root of unity found")


# ---------------------------------------------------------------------------
# Part A: production secp256k1 structure checks
# ---------------------------------------------------------------------------
P_FIELD = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
R_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8


def ec_add(Pt, Qt, mod, b=7):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    x1, y1 = Pt
    x2, y2 = Qt
    if x1 == x2 and (y1 + y2) % mod == 0:
        return None
    if Pt == Qt:
        m = (3 * x1 * x1) * pow(2 * y1, -1, mod) % mod
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, mod) % mod
    x3 = (m * m - x1 - x2) % mod
    return (x3, (m * (x1 - x3) - y1) % mod)


def ec_mul(k, Pt, mod, b=7):
    acc = None
    base = Pt
    while k:
        if k & 1:
            acc = ec_add(acc, base, mod, b)
        base = ec_add(base, base, mod, b)
        k >>= 1
    return acc


def production_checks():
    p, r, G = P_FIELD, R_ORDER, (GX, GY)
    out = {}
    assert p % 3 == 1 and r % 3 == 1 and r % 2 == 1 and (r - 1) % 6 == 0
    out["p_mod_3"] = p % 3
    out["r_mod_6"] = r % 6
    zeta = proot3(p)
    zeta2 = zeta * zeta % p
    out["zeta_cubed_is_1"] = pow(zeta, 3, p) == 1 and zeta != 1
    # zeta-action preserves the curve equation
    out["zeta_preserves_curve"] = (
        (zeta * GX % p) ** 3 % p + 7) % p == GY * GY % p
    lam0 = proot3(r)
    lam1 = lam0 * lam0 % r
    zG = (zeta * GX % p, GY)
    lam, lam2 = (lam0, lam1) if ec_mul(lam0, G, p) == zG else (lam1, lam0)
    out["lam_mult_of_generator_is_zeta_image"] = ec_mul(lam, G, p) == zG
    out["lam2_mult_is_zeta2_image"] = (
        ec_mul(lam2, G, p) == (zeta2 * GX % p, GY))
    out["lam_cubed_is_1"] = pow(lam, 3, r) == 1 and lam % r != 1
    out["classes_nonzero"] = (r - 1) // 6

    def canon(pt, a, bb, Q):
        x, y = pt
        cand = [(x, 0), (zeta * x % p, 1), (zeta2 * x % p, 2)]
        xc, i = min(cand)
        if i:
            f = lam if i == 1 else lam2
            a = f * a % r
            bb = f * bb % r
        if y % 2:
            a = (-a) % r
            bb = (-bb) % r
            y = p - y
        return (xc, y), a, bb

    rng = random.Random(20261022)
    rep_ok = 0
    for _ in range(12):
        x = rng.randrange(1, r)
        Q = ec_mul(x, G, p)
        a = rng.randrange(r)
        bb = rng.randrange(r)
        R = ec_add(ec_mul(a, G, p), ec_mul(bb, Q, p), p)
        Rc, a2, b2 = canon(R, a, bb, Q)
        ok = ec_add(ec_mul(a2, G, p), ec_mul(b2, Q, p), p) == Rc
        assert ok, "canonical representation bookkeeping failed"
        rep_ok += 1
    out["canonical_rep_checks_random_points"] = rep_ok
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale end-to-end attack
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)
TOY_K = 8          # walkers
TOY_NC = 32        # r-adding walk constants


def toy_attack(trials, dp):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    assert is_prime_trial(p) and is_prime_trial(n)
    assert p % 3 == 1 and n % 3 == 1 and n % 2 == 1 and (n - 1) % 6 == 0
    # generator really has prime order n (cofactor 1)
    assert ec_mul(n, G, p, b) is None and ec_mul(1, G, p, b) == G
    assert (G[0] ** 3 + b) % p == G[1] * G[1] % p
    zeta = proot3(p)
    zeta2 = zeta * zeta % p
    lam0 = proot3(n)
    lam1 = lam0 * lam0 % n
    zG = (zeta * G[0] % p, G[1])
    lam, lam2 = (lam0, lam1) if ec_mul(lam0, G, p, b) == zG else (lam1, lam0)
    assert ec_mul(lam, G, p, b) == zG

    def canon(pt, a, bb):
        x, y = pt
        cand = [(x, 0), (zeta * x % p, 1), (zeta2 * x % p, 2)]
        xc, i = min(cand)
        if i:
            f = lam if i == 1 else lam2
            a = f * a % n
            bb = f * bb % n
        if y % 2:
            a = (-a) % n
            bb = (-bb) % n
            y = p - y
        return (xc, y), a, bb

    rng = random.Random(12345)
    consts = []
    for _ in range(TOY_NC):
        u = rng.randrange(1, n)
        Ac, c, _ = canon(ec_mul(u, G, p, b), u, 0)
        assert ec_add(ec_mul(c, G, p, b), None, p, b) == Ac
        consts.append((Ac, c))

    n_classes = (n - 1) // 6
    budget = math.ceil(1.1 * math.sqrt(n_classes))
    ideal = 1 - math.exp(-budget * budget / (2 * n_classes))

    def one(secret, seed):
        rng = random.Random(seed)
        Q = ec_mul(secret, G, p, b)

        def fresh():
            a = rng.randrange(n)
            bb = rng.randrange(n)
            R = ec_add(ec_mul(a, G, p, b), ec_mul(bb, Q, p, b), p, b)
            return list(canon(R, a, bb))

        walkers = [fresh() for _ in range(TOY_K)]
        table = {}
        steps = 0
        degenerate = 0
        while steps < budget:
            for w in walkers:
                R, a, bb = w
                j = R[0] % TOY_NC
                S = ec_add(R, consts[j][0], p, b)
                a = (a + consts[j][1]) % n
                steps += 1
                if S is None:
                    w[:] = fresh()
                    continue
                Sc, a2, b2 = canon(S, a, bb)
                xc = Sc[0]
                if dp and xc % dp != 0:
                    w[:] = [Sc, a2, b2]
                    continue
                prev = table.get(xc)
                if prev is None:
                    table[xc] = (a2, b2)
                    w[:] = [Sc, a2, b2]
                    continue
                a1, b1 = prev
                if (a1 - a2) % n == 0 and (b1 - b2) % n == 0:
                    degenerate += 1
                    w[:] = fresh()
                    continue
                if b1 == b2:
                    degenerate += 1
                    w[:] = fresh()
                    continue
                xr = (a1 - a2) * pow(b2 - b1, -1, n) % n
                if ec_mul(xr, G, p, b) == Q:
                    return True, steps, degenerate, False
                return False, steps, degenerate, True
        return False, steps, degenerate, False

    successes = 0
    verif_fail = 0
    deg_total = 0
    for ti in range(trials):
        secret = rng.randrange(1, n)
        ok, st, deg, vf = one(secret, 9000 + ti)
        successes += ok
        verif_fail += vf
        deg_total += deg
    return {
        "curve": "y^2 = x^3 + 29 over F_1048609 (j=0, p mod 3 = 1)",
        "group_order_n": n,
        "classes_n_over_6": n_classes,
        "walkers": TOY_K,
        "walk_constants": TOY_NC,
        "budget_steps_all_walkers": budget,
        "budget_shape": "ceil(1.1 * sqrt((n-1)/6))",
        "distinguished_point_mask": dp,
        "trials": trials,
        "successes": successes,
        "empirical_success": successes / trials,
        "ideal_birthday_prediction": ideal,
        "extraction_verification_failures": verif_fail,
        "degenerate_restarts": deg_total,
    }


def main():
    data = json.load(sys.stdin)
    trials = data["trials"]
    dp = data["dp"]
    out = {
        "production_checks": production_checks(),
        "toy_rho": toy_attack(trials, dp),
        "scope": (
            "Part A verifies the order-6 automorphism structure and "
            "canonicalization bookkeeping at full secp256k1 parameters. "
            "Part B demonstrates the complete attack and its single-shot "
            "success rate only on a 20-bit toy curve; it is not evidence "
            "about production cost bounds or production success constants."
        ),
    }
    print(json.dumps(out))


main()
