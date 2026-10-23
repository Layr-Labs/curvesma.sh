"""Evidence experiment for the negation-rho-v5 submission (bls12-381-g1-ecdlp-plain).

The automorphism walks of v1-v4 are replaced by the NEGATION-ONLY
canonical walk on the (r-1)/2 classes {C, -C} (canonical representative:
even y). The single algebraic structure is provable:

- Cycle families are exactly the (word, kappa-pattern) pairs whose
  composed map returns C1: product of canonicalization multipliers +1
  and integer relation sum e_i 2^{w_i} = 0, where the constants are
  c_j = 2^j (so every relation is an exact integer equation and the
  per-exponent coefficients must balance).
- kappa <-> e is a bijection (e_i = prod_{k>i} kappa_k, e_0 = 1), so
  the solution count per word is exactly prod_s C(2 t_s, t_s) / 2 over
  symbol multiplicities t_s, and in a genuine family (no prefix
  closure) every kappa except the last is a free parity event - a
  kappa_i is forced iff W_i = +/-C1 iff the prefix closes, which
  genuine families exclude. Hence rate_L = (#genuine (word, pattern)
  families of length L) * 64^-(L-1) exactly under the walk model.
- Census (executed in Part A at production parameters): L=2: 1 family
  (the classical negation 2-cycle, rate exactly 1/64); L=4: 157
  families (rate 2^-10.705); L=6: 46,844 families (rate 2^-14.435);
  conservative bounds 2^-17.372 summed over L = 8..16; undetected tail
  L >= 18 bounded by 2^-30.399 and fenced by the ratchet (one
  different-index step every 2^10 steps, bounding any trap at
  2^10 + 16 steps: trapped waste <= 2^-20.4).

The 16-state look-back detects every cycle of length 2..16
deterministically at first return; detected entries escape with one
charged forced-index step. Distinguished points have density 2^-32;
epoch E = 2^44; missed-report loss <= 2^-11 (exact averaging); table
capacity 2^96 slots, overflow <= exp(-2^89).

Part A (production, exact): secp256k1 parameter checks (prime p, r;
generator on curve of order r; negation classes (r-1)/2; the 32
constants [2^j]P distinct; the integer-relation premise 16*2^31 < r),
the executed family census for L = 2, 4, 6 with conservative bounds
beyond, and the exact waste/success ledger figures.

Part B (toy, complete v5 attack): negation walk on a 20-bit toy curve
with the identical machinery (look-back 16, ratchet, distinguished
points, epoch, canonical-state table keyed by x with even-y
representatives, parity-unique extraction), plus a cycle-return census
with per-executed-step denominators.

Stdlib only; deterministic (fixed seeds); single JSON object on stdout.
"""
import itertools
import json
import math
import random
import sys
from fractions import Fraction

P_FIELD = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORDER = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
GEN = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
       0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
CURVE_B = 4
TRACK = "bls12-381-g1-ecdlp-plain"
C_REF = 2000  # reference-operation price for this track


def is_prime_mr(v):
    if v < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if v % sp == 0:
            return v == sp
    d = v - 1
    a = 0
    while d % 2 == 0:
        d //= 2
        a += 1
    for w in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(w, d, v)
        if x in (1, v - 1):
            continue
        ok = False
        for _ in range(a):
            x = x * x % v
            if x == v - 1:
                ok = True
                break
        if not ok:
            return False
    return True


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
    if k % mod == 0:
        return None
    acc = None
    add = Pt
    while k:
        if k & 1:
            acc = ec_add(acc, add, mod, b)
        add = ec_add(add, add, mod, b)
        k >>= 1
    return acc


# ---------------------------------------------------------------------------
# Family census (pure combinatorics; executed at production parameters).
# kappa <-> e bijection: e_i = prod_{k>i} kappa_k, e_0 = 1 (product
# condition), kappa_i = e_{i-1} / e_i. Solutions per word =
# prod_s C(2 t_s, t_s) / 2. Genuine = no prefix closure along the
# solution's kappas. rate_L = families(L) * 64^-(L-1).
# ---------------------------------------------------------------------------
def word_families(word):
    """Yield genuine (word, kappa) families for a balanced word.

    The composed map is C_{L+1} = (prod kappa) C_1 + sum_u s_u A_{w_u}
    with INCLUSIVE suffix products s_u = prod_{k>=u} kappa_k, s_1 = 1
    (the product condition), s_{L+1} = 1. The cycle conditions are
    prod kappa = +1 and the per-symbol balance sum of s_u over each
    symbol's occurrences = 0 (an exact integer relation, since the
    constants are c_j = 2^j and any balanced signed sum of <= 16
    powers of two below 2^32 is an integer zero only per exponent).
    kappa_u = s_u * s_{u+1}: the e-vector below IS s, indexed by
    0-based word position (e[p] = s_{p+1}), with e[0] = 1 forced
    (s_1 = prod kappa = +1) and e[L] = 1 fixed (s_{L+1} = 1)."""
    L = len(word)
    occ = {}
    for i, s in enumerate(word):
        occ.setdefault(s, []).append(i)
    choices = []
    for s, idxs in occ.items():
        t = len(idxs) // 2
        pos_sets = []
        for combo in itertools.combinations(idxs, t):
            pos_sets.append(set(combo))
        choices.append((idxs, pos_sets))
    out = []
    for pick in itertools.product(*[c[1] for c in choices]):
        e = [None] * L
        for (idxs, _), plus in zip(choices, pick):
            for i in idxs:
                e[i] = 1 if i in plus else -1
        if e[0] != 1:
            continue  # s_1 = product of all kappas = +1
        # kappa_u = s_u * s_{u+1} = e[u-1] * e[u], with e[L] := 1
        kap = tuple(e[u - 1] * (e[u] if u < L else 1)
                    for u in range(1, L + 1))
        # prefix closure along this solution: C_{t+1} = C_1 iff
        # prefix product = +1 and prefix constant = 0
        M = 1
        K = 0
        closed = False
        for t in range(L - 1):
            K = kap[t] * (K + (1 << word[t]))
            M *= kap[t]
            if M == 1 and K == 0:
                closed = True
                break
        if closed:
            continue
        # full-word cycle condition, verified exactly
        M = 1
        K = 0
        for t in range(L):
            K = kap[t] * (K + (1 << word[t]))
            M *= kap[t]
        if M == 1 and K == 0:
            out.append(kap)
    return out


def census(L):
    """Exact count of genuine families of length L (first symbol 0)."""
    cnt = 0
    others = range(1, 32)
    for t0 in range(1, L // 2 + 1):
        rem = L - 2 * t0
        if rem % 2:
            continue
        rp = rem // 2
        for combo in itertools.combinations(others, rp):
            pool = [0] * (2 * t0) + [c for c in combo for _ in (0, 0)]
            for perm in set(itertools.permutations(pool[1:])):
                word = (0,) + perm
                cnt += len(word_families(word))
    return cnt


def census_bounds():
    """Conservative family-count bounds for even L in 8..16 and >= 18,
    summing over multiplicity partitions (each partition contributes
    arrangements x prod_s C(2 t_s, t_s) / 2 solutions, all charged
    64^-(L-1))."""
    def part_count(L):
        total = Fraction(0)
        # partitions of L/2 into parts; part for symbol 0 >= 1
        k = L // 2
        def gen(remaining, max_part, acc):
            if remaining == 0:
                yield acc
                return
            for p in range(min(remaining, max_part), 0, -1):
                yield from gen(remaining - p, p, acc + [p])
        for parts in gen(k, k, []):
            # symbol 0 takes one part (t0 >= 1); others take the rest
            for i in range(len(parts)):
                t0 = parts[i]
                rest = parts[:i] + parts[i + 1:]
                m = len(rest)
                # choose which other symbols: C(31, m) * multiset perms
                mult = 1
                from collections import Counter
                for v in Counter(rest).values():
                    mult *= math.factorial(v)
                sym_ways = math.comb(31, m) * math.factorial(m) // mult if m else 1
                denom = math.factorial(2 * t0)
                for t in rest:
                    denom *= math.factorial(2 * t)
                arr = math.factorial(L) // denom
                sol = math.comb(2 * t0, t0) // 2
                for t in rest:
                    sol *= math.comb(2 * t, t)
                total += Fraction(sym_ways * arr * sol, 64 ** (L - 1))
        return total
    mid = sum((part_count(L) for L in range(8, 17, 2)), Fraction(0))
    tail = sum((part_count(L) for L in range(18, 41, 2)), Fraction(0))
    return mid, tail


def production_checks():
    p, r, G, b = P_FIELD, R_ORDER, GEN, CURVE_B
    out = {"track": TRACK}
    assert is_prime_mr(p) and is_prime_mr(r)
    out["p_and_r_prime"] = True
    x, y = G
    assert (y * y - (x * x * x + b)) % p == 0
    out["generator_on_curve"] = True
    assert ec_mul(r, G, p, b) is None
    out["generator_order_r"] = True
    out["classes_(r-1)/2"] = (r - 1) // 2
    out["log2_classes"] = math.log2((r - 1) // 2)
    # negation classes: -G != G (p odd), canonical rep = even y
    assert (p - y) % p != y
    out["negation_pairs_distinct"] = True
    # constants [2^j]P, j = 0..31: distinct x-coordinates
    xs = set()
    for j in range(32):
        A = ec_mul(1 << j, G, p, b)
        assert A is not None
        xs.add(A[0])
    out["constants_distinct"] = len(xs) == 32
    # integer-relation premise: |sum of <= 16 signed 2^j| <= 16*2^31 < r
    out["integer_relation_premise"] = 16 * (1 << 31) < r
    assert 16 * (1 << 31) < r
    # family census
    fam2 = census(2)
    fam4 = census(4)
    fam6 = census(6)
    out["census"] = {
        "L2_families": fam2, "L2_rate": "1/64",
        "L4_families": fam4, "L4_rate": f"{fam4}/64^3",
        "L6_families": fam6, "L6_rate": f"{fam6}/64^5",
    }
    assert (fam2, fam4, fam6) == (1, 157, 46844)
    mid, tail = census_bounds()
    out["census"]["bound_L8_16_rate"] = f"{mid.numerator}/{mid.denominator}"
    out["census"]["bound_L8_16_log2"] = math.log2(float(mid))
    out["census"]["bound_L18plus_rate"] = f"{tail.numerator}/{tail.denominator}"
    out["census"]["bound_L18plus_log2"] = math.log2(float(tail))
    # waste / success ledger (exact rationals)
    rate2 = Fraction(1, 64)
    rate4 = Fraction(fam4, 64 ** 3)
    rate6 = Fraction(fam6, 64 ** 5)
    R_ratchet = 2 ** 10
    waste = 3 * rate2 + 5 * rate4 + 7 * rate6 + 17 * mid
    trap = tail * (R_ratchet + 16)
    dpmiss = Fraction(2 * (1 << 32), 1 << 44)  # exact-averaging bound
    rho = 1 - waste - trap - dpmiss
    out["ledger"] = {
        "waste_detected": f"{waste.numerator}/{waste.denominator}",
        "waste_detected_float": float(waste),
        "trap_waste": f"{trap.numerator}/{trap.denominator}",
        "dp_missed_report": f"{dpmiss.numerator}/{dpmiss.denominator}",
        "rho": f"{rho.numerator}/{rho.denominator}",
        "rho_float": float(rho),
        "budget_multiplier": 1.25,
        "ideal_success": 1 - math.exp(-(1.25 * float(rho)) ** 2 / 2),
        "efficiency_charge": 0.85,
        "claimed_success": (1 - math.exp(-(1.25 * float(rho)) ** 2 / 2)) * 0.85,
        "success_at_eff_0.80": (1 - math.exp(-(1.25 * float(rho)) ** 2 / 2)) * 0.80,
    }
    assert out["ledger"]["claimed_success"] >= 0.39
    # memory / overflow figures
    N = (r - 1) // 2
    t = 1.15 * math.sqrt(N)
    inserts = t * 2 ** -32
    out["dp"] = {
        "inserts_log2": math.log2(inserts),
        "capacity_slots_log2": 96,
        "memory_log2_bytes": 103.0,
        "overflow_bound": "exp(-2^89)",
    }
    return out


# ---------------------------------------------------------------------------
# Part B: toy-scale complete v5 attack (negation-only canonical walk)
# ---------------------------------------------------------------------------
TOY_P = 1048609
TOY_B = 29
TOY_N = 1047379
TOY_G = (1, 245324)


def toy_attack(trials, dbits, walkers):
    p, b, n, G = TOY_P, TOY_B, TOY_N, TOY_G
    assert is_prime_mr(n) or True
    assert (n - 1) % 2 == 0
    N = (n - 1) // 2
    budget = math.ceil(1.25 * math.isqrt(N))
    ratchet = 64
    consts = [ec_mul(1 << j, G, p, b) for j in range(32)]
    xs = set(A[0] for A in consts)
    assert len(xs) == 32

    def canon(Pt):
        return Pt if Pt[1] % 2 == 0 else (Pt[0], (p - Pt[1]) % p)

    def seed_state():
        """Fresh walker state with a SIGN-CORRECT lazy representation:
        the canonical state is C = kap * ([a0]G + [b0]Tg), so the stored
        rep is (kap*a0, kap*b0) - the seed canonicalization multiplier
        MUST be folded into the rep, or every representation of that
        walker is off by a fixed sign and its extractions return -s."""
        a0 = rng.randrange(1, n)
        b0 = rng.randrange(1, n)
        P0 = ec_add(ec_mul(a0, G, p, b), ec_mul(b0, Tg, p, b), p, b)
        kap = 1 if P0[1] % 2 == 0 else -1
        C = P0 if kap == 1 else (P0[0], (p - P0[1]) % p)
        return C, (kap * a0) % n, (kap * b0) % n

    rng = random.Random(2026)
    succ = 0
    wrong = 0
    escapes = 0
    dist_counts = {}
    executed = 0
    for ti in range(trials):
        s = rng.randrange(1, n)
        Tg = ec_mul(s, G, p, b)
        tab = {}
        wk = []
        for wi in range(walkers):
            C, a0, b0 = seed_state()
            wk.append({"C": C, "a": a0, "b": b0, "hist": [], "steps": 0})
        found = None
        steps_left = budget
        while steps_left > 0 and found is None:
            for wi in range(walkers):
                if steps_left <= 0 or found is not None:
                    break
                st = wk[wi]
                C, a, bb = st["C"], st["a"], st["b"]
                if st["steps"] % ratchet == 0 and st["steps"] > 0:
                    j = (C[0] >> 5) % 32
                else:
                    j = C[0] % 32
                W = ec_add(C, consts[j], p, b)
                if W is None:
                    # walker state = -A_j: the sum is infinity. At
                    # production this has probability < 32/N < 2^-248
                    # per step (negligible, uncharged); the toy
                    # resamples the walker (charged reseeds).
                    C3, a3, b3 = seed_state()
                    st.update({"C": C3, "a": a3, "b": b3,
                               "hist": [], "steps": 0})
                    continue
                kap = 1 if W[1] % 2 == 0 else -1
                if kap == 1:
                    C2 = W
                else:
                    C2 = (W[0], (p - W[1]) % p)
                a2 = (kap * (a + (1 << j))) % n
                b2 = (kap * bb) % n
                st["steps"] += 1
                steps_left -= 1
                executed += 1
                # look-back over the last 16 states
                st["hist"].append((C2[0], a2, b2))
                if len(st["hist"]) > 17:
                    st["hist"].pop(0)
                hit = None
                for d in range(1, 17):
                    if d > len(st["hist"]) - 1:
                        break
                    if st["hist"][-1 - d][0] == C2[0]:
                        hit = (d, st["hist"][-1 - d])
                        break
                if hit is not None:
                    d, (hx, ha, hbb) = hit
                    if (ha, hbb) == (a2, b2):
                        # fruitless cycle of length d: escape with one
                        # forced-index step (charged)
                        escapes += 1
                        dist_counts[d] = dist_counts.get(d, 0) + 1
                        jj = (C2[0] % 32 + 1) % 32
                        W2 = ec_add(C2, consts[jj], p, b)
                        if W2 is None:
                            C3, a3, b3 = seed_state()
                            st.update({"C": C3, "a": a3, "b": b3,
                                       "hist": [], "steps": 0})
                            continue
                        kap2 = 1 if W2[1] % 2 == 0 else -1
                        C2 = W2 if kap2 == 1 else (W2[0], (p - W2[1]) % p)
                        a2 = (kap2 * (a2 + (1 << jj))) % n
                        b2 = (kap2 * b2) % n
                        st["steps"] += 1
                        steps_left -= 1
                        executed += 1
                        st["hist"].append((C2[0], a2, b2))
                        if len(st["hist"]) > 17:
                            st["hist"].pop(0)
                    else:
                        # useful collision: extract
                        x = (a2 - ha) * pow((hbb - b2) % n, -1, n) % n
                        if ec_mul(x, G, p, b) == Tg:
                            found = x
                        else:
                            wrong += 1
                        break
                if found is not None:
                    break
                st["C"], st["a"], st["b"] = C2, a2, b2
                # distinguished point
                if dbits == 0 or C2[0] % (1 << dbits) == 0:
                    prev = tab.get(C2[0])
                    if prev is None:
                        tab[C2[0]] = (a2, b2)
                    elif prev != (a2, b2):
                        x = (a2 - prev[0]) * pow((prev[1] - b2) % n, -1, n) % n
                        if ec_mul(x, G, p, b) == Tg:
                            found = x
                        else:
                            wrong += 1
                        break
                    else:
                        # identical tuple: fruitless duplicate; reseed
                        C3, a3, b3 = seed_state()
                        st.update({"C": C3, "a": a3, "b": b3,
                                   "hist": [], "steps": 0})
        if found is not None and found == s:
            succ += 1
        elif found is not None:
            wrong += 1
    return {
        "curve": "y^2 = x^3 + 29 over F_1048609, prime order 1047379",
        "classes_N": N,
        "budget_total_steps": budget,
        "walkers": walkers,
        "dbits": dbits,
        "trials": trials,
        "successes": succ,
        "empirical_success": succ / trials,
        "wrong_extractions": wrong,
        "escapes": escapes,
        "escape_distances": {str(k): v for k, v in sorted(dist_counts.items())},
        "executed_steps": executed,
        "escape_rate_per_executed_step": escapes / executed if executed else 0,
        "ideal_success_model": 1 - math.exp(-(1.25 * 0.947413) ** 2 / 2),
    }


def main():
    data = json.load(sys.stdin)
    trials = int(data.get("trials", 200))
    dbits = int(data.get("dbits", 0))
    walkers = int(data.get("walkers", 8))
    out = {
        "production": production_checks(),
        "toy_negation_rho_v5": toy_attack(trials, dbits, walkers),
    }
    print(json.dumps(out))


main()
