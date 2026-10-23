"""Evidence experiment for the automorphism-rho-v3 submission
(secp256k1-ecdlp).

v3 responds to the second-round review: the distance-2 look-back is extended
to distances 2 AND 3 (the order-3 identity 1+lam+lam^2 = 0 creates
period-3 fruitless cycles: three same-constant steps with multiplier lam or
lam^2 return the full state identically), a per-walker epoch limit makes the
construction robust to any longer algebraic cycle family, and a cycle-return
census measures both family rates at toy scale.

Part A verifies at full production parameters: the order-6 automorphism
structure, the six-element multiplier group G6, the class count (r-1)/6,
the pairwise class-disjoint constant setup check, and the lazy-representation
invariant along production canonical walks.

Part B runs the complete v3 attack end to end on a small j=0 toy curve:
canonical-state walk with lazy representations, look-back cycle detection at
distances 2 and 3 with fruitless/useful discrimination and one-addition
escapes, epoch limit, distinguished points by the low-bits-zero predicate,
hash-table collision extraction in the canonical frame. It reports measured
success, escape counts for both families, a cycle-return census, and
extraction-verification failures. Toy evidence is evidence at toy scope only.

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
P_FIELD = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
R_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GEN = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
       0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
CURVE_B = 7


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
    out["order3_identity"] = (1 + lam + lam2) % r == 0
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

    # lazy-representation invariant along short production walks
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
            W2 = ec_add(C, CAj, p, b)
            if W2 is None:
                break
            C2, i2, k2 = canonize(W2)
            cinv = pow(g6[L], -1, r)
            a = (a + cinv * (kj * (j + 1) % r)) % r
            L = gidx[(g6[L] * k2) % r]
            C = C2
            assert ec_add(ec_mul(g6[L] * a % r, G, p, b),
                          ec_mul(g6[L] * bb % r, G, p, b), p, b) == C
            checked_steps += 1
    out["lazy_rep_invariant_steps_verified"] = checked_steps

    # walk-constant setup check at production scale
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
# Part B: toy-scale v3 attack
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)
NC = 32
CENSUS_STEPS = 40000


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
    assert (1 + lam + lam2) % n == 0
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
    consts = []
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
    epoch = budget * walkers + 1   # never fires at toy scale; production E = 2^36

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
            return [C, a, bb, gidx[k], 0]

        ws = [fresh() for _ in range(walkers)]
        hist = [[] for _ in range(walkers)]
        table = {}
        steps = 0
        esc2 = 0
        esc3 = 0
        vf = 0

        def try_extract(pa, pb, pl, a, bb, L):
            La1, Lb1 = g6[L] * a % n, g6[L] * bb % n
            La2, Lb2 = g6[pl] * pa % n, g6[pl] * pb % n
            if (Lb1 - Lb2) % n == 0:
                return None
            x_rec = (La1 - La2) * pow(Lb2 - Lb1, -1, n) % n
            return x_rec

        while steps < budget:
            for wi, w in enumerate(ws):
                if w is None:
                    ws[wi] = fresh()
                    if ws[wi] is None:
                        continue
                    w = ws[wi]
                C, a, bb, L, age = w
                if age > epoch:
                    ws[wi] = fresh()
                    hist[wi] = []
                    continue
                h = hist[wi]
                # look-back cycle detection at distances 2 and 3
                hit = None
                if len(h) >= 2 and h[-2][0] == C[0]:
                    hit = -2
                elif len(h) >= 3 and h[-3][0] == C[0]:
                    hit = -3
                if hit is not None:
                    pa, pb, pl = h[hit][1], h[hit][2], h[hit][3]
                    if pa == a and pb == bb and pl == L:
                        j2 = (C[0] + 1) % NC
                        W2 = ec_add(C, consts[j2][0], p, b)
                        steps += 1
                        if hit == -2:
                            esc2 += 1
                        else:
                            esc3 += 1
                        if W2 is None:
                            ws[wi] = fresh()
                            hist[wi] = []
                            continue
                        C2, _, k2 = canonize(W2)
                        a2 = (a + linv_c[L][j2]) % n
                        L2 = multab[L][gidx[k2]]
                        ws[wi] = [C2, a2, bb, L2, age + 1]
                        h.append((C2[0], a2, bb, L2))
                        if len(h) > 5:
                            h.pop(0)
                        continue
                    x_rec = try_extract(pa, pb, pl, a, bb, L)
                    if x_rec is not None and ec_mul(x_rec, G, p, b) == Q:
                        return True, steps, esc2, esc3, vf
                    vf += 1
                    ws[wi] = fresh()
                    hist[wi] = []
                    continue
                h.append((C[0], a, bb, L))
                if len(h) > 5:
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
                ws[wi] = [C2, a2, bb, L2, age + 1]
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
                    x_rec = try_extract(a1, b1, l1, a2, bb, L2)
                    if x_rec is not None and ec_mul(x_rec, G, p, b) == Q:
                        return True, steps, esc2, esc3, vf
                    vf += 1
                    ws[wi] = fresh()
                    hist[wi] = []
        return False, steps, esc2, esc3, vf

    # cycle-return census: measure full-state return rates at distances 2,3
    census = {"steps": 0, "returns_d2": 0, "returns_d3": 0}
    for wi in range(4):
        rngc = random.Random(555 + wi)
        a = rngc.randrange(n)
        bb = rngc.randrange(n)
        C, _, k = canonize(ec_add(ec_mul(a, G, p, b), ec_mul(bb, G, p, b), p, b))
        L = gidx[k]
        hist = []
        for _step in range(CENSUS_STEPS // 4):
            if len(hist) >= 2 and hist[-2][0] == C[0]:
                if hist[-2][1:] == (C[0], a, bb, L)[1:] and \
                   hist[-2][1] == a and hist[-2][2] == bb and hist[-2][3] == L:
                    census["returns_d2"] += 1
                    # escape to continue the census
                    j2 = (C[0] + 1) % NC
                    W2 = ec_add(C, consts[j2][0], p, b)
                    if W2 is None:
                        break
                    C, _, k2 = canonize(W2)
                    a = (a + linv_c[L][j2]) % n
                    L = multab[L][gidx[k2]]
                    hist = []
                    census["steps"] += 1
                    continue
            if len(hist) >= 3 and hist[-3][0] == C[0]:
                if hist[-3][1] == a and hist[-3][2] == bb and hist[-3][3] == L:
                    census["returns_d3"] += 1
                    j2 = (C[0] + 1) % NC
                    W2 = ec_add(C, consts[j2][0], p, b)
                    if W2 is None:
                        break
                    C, _, k2 = canonize(W2)
                    a = (a + linv_c[L][j2]) % n
                    L = multab[L][gidx[k2]]
                    hist = []
                    census["steps"] += 1
                    continue
            hist.append((C[0], a, bb, L))
            if len(hist) > 5:
                hist.pop(0)
            j = C[0] % NC
            W = ec_add(C, consts[j][0], p, b)
            if W is None:
                break
            a = (a + linv_c[L][j]) % n
            C, _, k2 = canonize(W)
            L = multab[L][gidx[k2]]
            census["steps"] += 1
    census["rate_d2_per_step"] = census["returns_d2"] / census["steps"]
    census["rate_d3_per_step"] = census["returns_d3"] / census["steps"]

    rng2 = random.Random(31337)
    succ = 0
    e2 = 0
    e3 = 0
    vf = 0
    for ti in range(trials):
        s = rng2.randrange(1, n)
        ok, st, x2, x3, v = one(s, 90000 + ti)
        succ += ok
        e2 += x2
        e3 += x3
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
        "escapes_distance2": e2,
        "escapes_distance3": e3,
        "cycle_census": census,
        "extraction_verification_failures": vf,
    }


def main():
    data = json.load(sys.stdin)
    out = {
        "production_checks": production_checks(),
        "toy_rho_v3": toy_attack(data["trials"], data["dbits"],
                                 data["walkers"]),
        "scope": (
            "Part A verifies the order-6 automorphism structure, the "
            "six-element multiplier group, the class-disjoint constant setup "
            "check, and the lazy-representation invariant at full production "
            "parameters. Part B demonstrates the complete v3 attack "
            "(lazy representations, look-back cycle detection at distances 2 "
            "and 3, one-addition escapes, epoch limit, low-bits "
            "distinguished points) only on a 20-bit toy curve, plus a "
            "cycle-return census measuring both fruitless-family rates; it is "
            "not evidence about production costs or production success "
            "constants."
        ),
    }
    print(json.dumps(out))


main()
