"""Evidence experiment for the automorphism-rho-v2 submission
(bls12-381-g1-ecdlp-plain).

Part A verifies, at full production parameters and with pure standard-library
arithmetic: the order-6 automorphism structure (zeta, lam with
[lam]P = zeta(P), the six-element multiplier group G6 = {+-1, +-lam, +-lam^2}
closed under multiplication), the class count (r-1)/6, the pairwise
class-disjoint walk-constant setup check, and - the core of the v2 design -
the lazy-representation invariant rep(canon(W)) = Lambda*(a,b) verified at
every step along short production canonical walks (canonicalization
multiplier, one mod-r addition per step, Lambda updated by a 6x6 table; the
canonical-frame representation is materialized only at collisions).

Part B runs the complete v2 attack end to end on a small j=0 toy curve with
the production-shaped machinery: canonical-state walk with lazy
representations, 2-step look-back cycle detector with fruitless/useful
discrimination and one-addition escape, distinguished points by low-bits-zero
predicate, hash-table collision extraction in the canonical frame. It
reports measured success against the ideal birthday prediction, escape
counts and rates, and extraction-verification failures. Toy evidence is
evidence at toy scope only.

Stdlib only; deterministic seeds; reads JSON input on stdin.
"""
import json
import math
import random
import sys


def is_prime_trial(v):
    if v < 2:
        return False
    for d in range(2, math.isqrt(v) + 1):
        if v % d == 0:
            return False
    return True


def proot3(mod):
    for g in range(2, 64):
        v = pow(g, (mod - 1) // 3, mod)
        if v != 1:
            assert pow(v, 3, mod) == 1
            return v
    raise AssertionError("no primitive cube root of unity found")


def ec_add(Pt, Qt, mod, b):
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


def ec_mul(k, Pt, mod, b):
    acc = None
    base = Pt
    while k:
        if k & 1:
            acc = ec_add(acc, base, mod, b)
        base = ec_add(base, base, mod, b)
        k >>= 1
    return acc


# ---------------------------------------------------------------------------
# Part A: production parameter checks
# ---------------------------------------------------------------------------
P_FIELD = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORDER = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
GEN = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
       0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
CURVE_B = 4


def production_checks():
    p, r, G, b = P_FIELD, R_ORDER, GEN, CURVE_B
    out = {}
    assert p % 3 == 1 and r % 3 == 1 and (r - 1) % 6 == 0
    out["p_mod_3"] = p % 3
    out["r_mod_6"] = r % 6
    zeta = proot3(p)
    zeta2 = zeta * zeta % p
    out["zeta_cubed_is_1"] = pow(zeta, 3, p) == 1 and zeta != 1
    out["zeta_preserves_curve"] = (
        (zeta * G[0] % p) ** 3 % p + b) % p == G[1] * G[1] % p
    lam0 = proot3(r)
    lam1 = lam0 * lam0 % r
    zG = (zeta * G[0] % p, G[1])
    lam, lam2 = (lam0, lam1) if ec_mul(lam0, G, p, b) == zG else (lam1, lam0)
    out["lam_mult_of_generator_is_zeta_image"] = ec_mul(lam, G, p, b) == zG
    out["lam2_mult_is_zeta2_image"] = (
        ec_mul(lam2, G, p, b) == (zeta2 * G[0] % p, G[1]))
    out["lam_cubed_is_1"] = pow(lam, 3, r) == 1 and lam % r != 1
    out["classes_nonzero"] = (r - 1) // 6
    g6 = [1, lam, lam2, (-1) % r, (-lam) % r, (-lam2) % r]
    out["g6_closed_under_multiplication"] = all(
        v * w % r in g6 for v in g6 for w in g6)
    gidx = {v: i for i, v in enumerate(g6)}

    def canonize(W):
        x, y = W
        cand = [(x, 0), (zeta * x % p, 1), (zeta2 * x % p, 2)]
        xc, i = min(cand)
        if y % 2:
            return (xc, p - y), i, g6[3 + i]     # multiplier -lam^i
        return (xc, y), i, g6[i]                 # multiplier +lam^i

    # lazy-representation invariant along short production walks:
    # state (C, a, b, L) with C canonical and rep(C) = L*(a, b) mod r.
    rng = random.Random(20261022)
    checked_steps = 0
    for _ in range(6):
        a0 = rng.randrange(1, r)
        b0 = rng.randrange(1, r)
        W = ec_add(ec_mul(a0, G, p, b), ec_mul(b0, G, p, b), p, b)
        C, i0, k0 = canonize(W)
        a, bb, L = a0, b0, gidx[k0]
        for _step in range(16):
            j = C[0] % 32
            # (walk constants built inline for the check: use u_j = j+1)
            Aj = ec_mul(j + 1, G, p, b)
            CAj, ij, kj = canonize(Aj)
            # rep(CAj) = kj*(j+1)
            W2 = ec_add(C, CAj, p, b)
            if W2 is None:
                break
            C2, i2, k2 = canonize(W2)
            # lazy update: a <- a + L^{-1} * c_j, b unchanged, L <- L*k2
            cinv = pow(g6[L], -1, r)
            a = (a + cinv * (kj * (j + 1) % r)) % r
            L = gidx[(g6[L] * k2) % r]
            C = C2
            assert ec_add(ec_mul(g6[L] * a % r, G, p, b),
                          ec_mul(g6[L] * bb % r, G, p, b), p, b) == C
            checked_steps += 1
    out["lazy_rep_invariant_steps_verified"] = checked_steps

    # walk-constant setup check at production scale: pairwise class-disjoint
    consts_ok = 0
    xs = []
    rng2 = random.Random(777)
    while len(xs) < 32:
        A = ec_mul(rng2.randrange(1, r), G, p, b)
        C, _, _ = canonize(A)
        if C[0] in xs or any(C[0] in (zeta * x % p, zeta2 * x % p)
                             for x in xs):
            continue
        xs.append(C[0])
        consts_ok += 1
    out["constants_pairwise_class_disjoint"] = consts_ok
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale v2 attack
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)
NC = 32


def toy_attack(trials, dbits, walkers):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    assert is_prime_trial(p) and is_prime_trial(n)
    assert p % 3 == 1 and n % 3 == 1 and (n - 1) % 6 == 0
    assert ec_mul(n, G, p, b) is None
    zeta = proot3(p)
    zeta2 = zeta * zeta % p
    lam0 = proot3(n)
    lam1 = lam0 * lam0 % n
    zG = (zeta * G[0] % p, G[1])
    lam, lam2 = (lam0, lam1) if ec_mul(lam0, G, p, b) == zG else (lam1, lam0)
    g6 = [1, lam, lam2, (-1) % n, (-lam) % n, (-lam2) % n]
    gidx = {v: i for i, v in enumerate(g6)}
    multab = [[gidx[(g6[i] * g6[j]) % n] for j in range(6)] for i in range(6)]
    inv6 = {i: gidx[pow(g6[i], -1, n)] for i in range(6)}

    def canonize(W):
        x, y = W
        cand = [(x, 0), (zeta * x % p, 1), (zeta2 * x % p, 2)]
        xc, i = min(cand)
        if y % 2:
            return (xc, p - y), i, g6[3 + i]
        return (xc, y), i, g6[i]

    rng = random.Random(4242)
    consts = []  # (canonical point, c_j) pairwise class-disjoint
    while len(consts) < NC:
        u = rng.randrange(1, n)
        A = ec_mul(u, G, p, b)
        C, i0, k0 = canonize(A)
        if any(C[0] in (c[0][0], zeta * c[0][0] % p, zeta2 * c[0][0] % p)
               for c in consts):
            continue
        cj = k0 * u % n
        assert ec_add(ec_mul(cj, G, p, b), None, p, b) == C
        consts.append((C, cj))
    linv_c = [[(g6[inv6[li]] * consts[j][1]) % n for j in range(NC)]
              for li in range(6)]

    n_classes = (n - 1) // 6
    budget = math.ceil(1.15 * math.sqrt(n_classes))
    ideal = 1 - math.exp(-budget * budget / (2 * n_classes))
    mask = (1 << dbits) - 1

    def one(secret, seed):
        rng = random.Random(seed)
        Q = ec_mul(secret, G, p, b)

        def fresh():
            a = rng.randrange(n)
            bb = rng.randrange(n)
            W = ec_add(ec_mul(a, G, p, b), ec_mul(bb, Q, p, b), p, b)
            if W is None:
                return None
            C, _, k = canonize(W)
            return [C, a, bb, gidx[k]]

        ws = [fresh() for _ in range(walkers)]
        hist = [[] for _ in range(walkers)]
        table = {}
        steps = 0
        escapes = 0
        vf = 0
        while steps < budget:
            for wi, w in enumerate(ws):
                if w is None:
                    ws[wi] = fresh()
                    if ws[wi] is None:
                        continue
                    w = ws[wi]
                C, a, bb, L = w
                h = hist[wi]
                # 2-step look-back detector on the full canonical state
                if len(h) >= 2 and h[-2][0] == C[0]:
                    pa, pb, pl = h[-2][1], h[-2][2], h[-2][3]
                    if pa == a and pb == bb and pl == L:
                        # fruitless cycle: one-addition escape
                        j2 = (C[0] + 1) % NC
                        W2 = ec_add(C, consts[j2][0], p, b)
                        steps += 1
                        escapes += 1
                        if W2 is None:
                            ws[wi] = fresh()
                            hist[wi] = []
                            continue
                        C2, _, k2 = canonize(W2)
                        a2 = (a + linv_c[L][j2]) % n
                        L2 = multab[L][gidx[k2]]
                        ws[wi] = [C2, a2, bb, L2]
                        h.append((C2[0], a2, bb, L2))
                        if len(h) > 4:
                            h.pop(0)
                        continue
                    # useful self-collision: extract in the canonical frame
                    La1, Lb1 = g6[L] * a % n, g6[L] * bb % n
                    La2, Lb2 = g6[pl] * pa % n, g6[pl] * pb % n
                    if (Lb1 - Lb2) % n != 0:
                        x_rec = (La1 - La2) * pow(Lb2 - Lb1, -1, n) % n
                        if ec_mul(x_rec, G, p, b) == Q:
                            return True, steps, escapes, vf
                    vf += 1
                    ws[wi] = fresh()
                    hist[wi] = []
                    continue
                h.append((C[0], a, bb, L))
                if len(h) > 4:
                    h.pop(0)
                j = C[0] % NC
                W = ec_add(C, consts[j][0], p, b)
                steps += 1
                if W is None:
                    ws[wi] = fresh()
                    hist[wi] = []
                    continue
                a2 = (a + linv_c[L][j]) % n
                C2, _, k2 = canonize(W)
                L2 = multab[L][gidx[k2]]
                ws[wi] = [C2, a2, bb, L2]
                # distinguished point: low dbits of canonical x are zero
                if dbits == 0 or (C2[0] & mask) == 0:
                    prev = table.get(C2[0])
                    if prev is None:
                        table[C2[0]] = (a2, bb, L2)
                        continue
                    a1, b1, l1 = prev
                    if a1 == a2 and b1 == bb and l1 == L2:
                        ws[wi] = fresh()
                        hist[wi] = []
                        continue
                    La1, Lb1 = g6[l1] * a1 % n, g6[l1] * b1 % n
                    La2, Lb2 = g6[L2] * a2 % n, g6[L2] * bb % n
                    if (Lb1 - Lb2) % n == 0:
                        ws[wi] = fresh()
                        hist[wi] = []
                        continue
                    x_rec = (La1 - La2) * pow(Lb2 - Lb1, -1, n) % n
                    if ec_mul(x_rec, G, p, b) == Q:
                        return True, steps, escapes, vf
                    vf += 1
                    ws[wi] = fresh()
                    hist[wi] = []
        return False, steps, escapes, vf

    rng2 = random.Random(31337)
    succ = 0
    esc = 0
    vf = 0
    for ti in range(trials):
        s = rng2.randrange(1, n)
        ok, st, e, v = one(s, 90000 + ti)
        succ += ok
        esc += e
        vf += v
    return {
        "curve": "y^2 = x^3 + 29 over F_1048609 (j=0, p mod 3 = 1)",
        "group_order_n": n,
        "classes_n_over_6": n_classes,
        "walkers": walkers,
        "walk_constants": NC,
        "budget_steps_all_walkers": budget,
        "budget_shape": "ceil(1.15 * sqrt((n-1)/6))",
        "distinguished_low_bits": dbits,
        "trials": trials,
        "successes": succ,
        "empirical_success": succ / trials,
        "ideal_birthday_prediction": ideal,
        "escapes": esc,
        "escape_rate_per_step": esc / (trials * budget),
        "extraction_verification_failures": vf,
    }


def main():
    data = json.load(sys.stdin)
    out = {
        "production_checks": production_checks(),
        "toy_rho_v2": toy_attack(data["trials"], data["dbits"],
                                 data["walkers"]),
        "scope": (
            "Part A verifies the order-6 automorphism structure, the "
            "six-element multiplier group, the class-disjoint constant setup "
            "check, and the lazy-representation invariant at full "
            "production parameters. Part B demonstrates the complete v2 "
            "attack (canonical-state lazy representations, 2-step look-back "
            "cycle detector with escape, low-bits distinguished points) only "
            "on a 20-bit toy curve; it is not evidence about production "
            "costs or production success constants."
        ),
    }
    print(json.dumps(out))


main()
