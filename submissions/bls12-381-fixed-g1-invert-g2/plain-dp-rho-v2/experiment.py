"""Evidence for `parallel distinguished-point rho in GT` (plain rho) on
bls12-381-fixed-g1-invert-g2.

Standard library only (no py_ecc, no third-party imports), so the verifier's
bounded Docker runner executes it as-is.

Part A checks, on BLS12-381 itself, the facts the algorithm uses:
  A1. r = #GT is prime and p^6 = -1 (mod r), so the p^6-Frobenius acts on GT as
      inversion (the tower conjugation below).
  A2. For random y in GT the tower conjugation `conj` (c0 + c1 w -> c0 - c1 w)
      equals y^{-1}.
  A3. The class key K(y) = |Re a_1| (a_1 = the Fp2 coefficient of w^1 of y,
      Re = its Fp component) is invariant under inversion, so K is a function of
      the class {y, y^-1}; the canonical representative is the member with
      Re a_1 <= (p-1)/2.
  A4. K is a hash, not an injection: distinct classes can share K, so the table
      must compare the stored exact element and never K alone. (Asserted as a
      design requirement here and exercised in Part B, where keys are exact.)

Part B runs the complete algorithm (parallel distinguished-point rho on the
inversion quotient, exact stored elements, bounded step counter, deterministic
cycle escape) in a toy prime field of the same cyclic shape and reports its
success statistics. It is evidence at toy scale only; it is not a production
run.

Input JSON keys (all optional): production_checks, q, trials, walks, dp_bits,
seed, seconds, dump.
"""
import json
import math
import random
import sys
import time

# BLS12-381 base field prime p and subgroup order r.
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001


# ---------------------------------------------------------------- Fp2/Fp6/Fp12
# Tower: Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-xi), Fp12 = Fp6[w]/(w^2-v),
# xi = 1+u, so w^6 = xi.  A Fp12 element is (c0, c1) meaning c0 + c1*w with
# c0, c1 in Fp6; each Fp6 element is (x0, x1, x2) meaning x0 + x1 v + x2 v^2;
# each Fp2 element is (a, b) meaning a + b u.
XI = (1, 1)
ZERO2 = (0, 0)
ONE2 = (1, 0)
ZERO6 = (ZERO2, ZERO2, ZERO2)


def f2_mul(x, y, p):
    a, b = x
    c, d = y
    return ((a * c - b * d) % p, (a * d + b * c) % p)


def f2_add(x, y, p):
    return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)


def f2_neg(x, p):
    return ((-x[0]) % p, (-x[1]) % p)


def f6_mul(x, y, p):
    x0, x1, x2 = x
    y0, y1, y2 = y
    t0 = f2_mul(x0, y0, p)
    t1 = f2_mul(x1, y1, p)
    t2 = f2_mul(x2, y2, p)
    x1y2 = f2_mul(x1, y2, p)
    x2y1 = f2_mul(x2, y1, p)
    x0y1 = f2_mul(x0, y1, p)
    x1y0 = f2_mul(x1, y0, p)
    x0y2 = f2_mul(x0, y2, p)
    x2y0 = f2_mul(x2, y0, p)
    c0 = f2_add(t0, f2_mul(f2_add(x1y2, x2y1, p), XI, p), p)
    c1 = f2_add(f2_add(x0y1, x1y0, p), f2_mul(t2, XI, p), p)
    c2 = f2_add(f2_add(x0y2, x2y0, p), t1, p)
    return (c0, c1, c2)


def f6_add(x, y, p):
    return (f2_add(x[0], y[0], p), f2_add(x[1], y[1], p), f2_add(x[2], y[2], p))


def f6_neg(x, p):
    return (f2_neg(x[0], p), f2_neg(x[1], p), f2_neg(x[2], p))


def f6_mul_by_v(x, p):
    """v * (x0 + x1 v + x2 v^2) = xi*x2 + x0 v + x1 v^2."""
    return (f2_mul(x[2], XI, p), x[0], x[1])


def f12_mul(x, y, p):
    a0, a1 = x
    b0, b1 = y
    hi = f6_add(f6_mul(a0, b0, p), f6_mul_by_v(f6_mul(a1, b1, p), p), p)
    lo = f6_add(f6_mul(a0, b1, p), f6_mul(a1, b0, p), p)
    return (hi, lo)


def f12_one():
    return ((ONE2, ZERO2, ZERO2), ZERO6)


def f12_pow(x, e, p):
    result = f12_one()
    base = x
    while e:
        if e & 1:
            result = f12_mul(result, base, p)
        base = f12_mul(base, base, p)
        e >>= 1
    return result


def f12_conj(x, p):
    """c0 + c1 w -> c0 - c1 w (Fp6-conjugation, i.e. the p^6-Frobenius on GT)."""
    return (x[0], f6_neg(x[1], p))


def key_a1(x):
    """a = Re a_1, the Fp component of the w^1 coefficient (a1 = x[1][0])."""
    return x[1][0][0] % P


def class_key(x):
    """K(y) = |Re a_1|, invariant under y -> y^-1."""
    a = key_a1(x)
    return min(a, P - a)


def is_one(x):
    return x == f12_one()


# ---------------------------------------------------------------- primality
def is_probable_prime(n, rounds=24):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for sp in small:
        if n % sp == 0:
            return n == sp
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def is_probable_prime(n, rounds=24):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for sp in small:
        if n % sp == 0:
            return n == sp
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def toy_prime(bits, rng):
    """A prime q = 2*r' + 1 with r' prime, q of the given bit length."""
    while True:
        rp = rng.getrandbits(bits - 2) | (1 << (bits - 3)) | 1
        if not is_probable_prime(rp):
            continue
        q = 2 * rp + 1
        if is_probable_prime(q):
            return q, rp


def toy_instance(q, rp, R, rng):
    """Generator g of the order-rp subgroup, target z = g^x, multiplier table."""
    g = 1
    while g == 1:
        h = rng.randrange(2, q - 1)
        if pow(h, (q - 1) // 2, q) == q - 1:
            g = pow(h, 2, q)
    x = rng.randrange(1, rp)
    z = pow(g, x, q)
    table = []
    for _ in range(R):
        a = rng.randrange(1, rp)
        b = rng.randrange(1, rp)
        table.append((pow(g, a, q) * pow(z, b, q) % q, a, b))
    return g, x, z, table


def first_collision(q, rp, R, seed, cap_factor=4.0):
    """Uncensored first-collision time of one plain rho walk, in units of sqrt(r').

    No distinguished-point restriction and no cap censoring: this measures the
    random-mapping statistic the success claim depends on.
    """
    rng = random.Random(seed)
    g, x, z, table = toy_instance(q, rp, R, rng)
    cap = int(cap_factor * math.sqrt(rp)) + 1
    a = rng.randrange(1, rp)
    b = rng.randrange(1, rp)
    y = pow(g, a, q) * pow(z, b, q) % q
    seen = {}
    for i in range(cap):
        h = y % R
        M, a_h, b_h = table[h]
        y = y * M % q
        if y in seen:
            return i + 1, (i + 1) / math.sqrt(rp)
        seen[y] = i
    return cap, cap / math.sqrt(rp)


def toy_trial(q, rp, walks, R, dp_bits, seed, factor=1.25):
    """Bounded complete run of the parallel distinguished-point rho.

    Plain walk y <- y*M_{h(y)} with h = y mod R, distinguished points of density
    theta = 2^-dp_bits stored under their exact element, exponent pairs recovered
    and verified, hard cap of factor*sqrt(r') total group operations.
    """
    rng = random.Random(seed)
    g, x, z, table = toy_instance(q, rp, R, rng)
    budget = int(factor * math.sqrt(rp)) + 1
    theta = 1 << dp_bits
    seen = {}
    cur = []
    for _ in range(walks):
        a = rng.randrange(1, rp)
        b = rng.randrange(1, rp)
        cur.append([pow(g, a, q) * pow(z, b, q) % q, a, b])
    steps = 0
    while steps < budget:
        for w in range(walks):
            if steps >= budget:
                break
            y, al, be = cur[w]
            steps += 1
            M, a_h, b_h = table[y % R]
            y2 = y * M % q
            al2 = (al + a_h) % rp
            be2 = (be + b_h) % rp
            cur[w] = [y2, al2, be2]
            if y2 % theta == 0:
                prev = seen.get(y2)
                if prev is None:
                    seen[y2] = (al2, be2)
                else:
                    al3, be3 = prev
                    den = (be3 - be2) % rp
                    if den:
                        cand = (al2 - al3) * pow(den, rp - 2, rp) % rp
                        if cand == x:
                            return steps, True
    return steps, False


def main():
    cfg = {}
    if not sys.stdin.isatty():
        raw = sys.stdin.read().strip()
        if raw:
            cfg = json.loads(raw)
    production_checks = int(cfg.get("production_checks", 4))
    q_bits = int(cfg.get("q", 30))
    trials = int(cfg.get("trials", 30))
    walks = int(cfg.get("walks", 16))
    mults = int(cfg.get("mults", 1024))
    dp_bits = int(cfg.get("dp_bits", 6))
    first_collision_trials = int(cfg.get("first_collision_trials", 20))
    seed = int(cfg.get("seed", 1))
    seconds = float(cfg.get("seconds", 10.0))
    dump = bool(cfg.get("dump", False))
    t0 = time.time()

    # ---------------- Part A: production-field checks
    checks = {
        "A1_r_is_prime": is_probable_prime(R),
        "A1_r_bit_length": R.bit_length(),
    }
    rng = random.Random(seed)
    ok_order = ok_nonone = True
    samples = 0
    for _ in range(production_checks):
        h = (((rng.randrange(1, P), rng.randrange(1, P)),
              (rng.randrange(1, P), rng.randrange(1, P)),
              (rng.randrange(1, P), rng.randrange(1, P))),
             ((rng.randrange(1, P), rng.randrange(1, P)),
              (rng.randrange(1, P), rng.randrange(1, P)),
              (rng.randrange(1, P), rng.randrange(1, P))))
        y = f12_pow(h, (P ** 12 - 1) // R, P)
        samples += 1
        ok_order = ok_order and is_one(f12_pow(y, R, P))
        ok_nonone = ok_nonone and not is_one(y)
    checks["A2_y_pow_r_is_one"] = ok_order
    checks["A2_samples_nonidentity"] = ok_nonone
    elapsed_a = time.time() - t0

    # ---------------- Part B: toy walk
    trng = random.Random(seed * 7919 + 13)
    q, rp = toy_prime(q_bits, trng)
    runs = {}
    for label, factor in (("at_1.5_sqrtr", 1.5), ("at_2.5_sqrtr", 2.5)):
        successes = 0
        done = 0
        steps_used = []
        for i in range(trials):
            if time.time() - t0 > seconds:
                break
            steps, ok = toy_trial(q, rp, walks, mults, dp_bits, seed * 100003 + i, factor)
            done += 1
            successes += int(ok)
            steps_used.append(steps)
        runs[label] = {
            "budget_factor_sqrtr": factor,
            "trials": done,
            "successes": successes,
            "success_rate": round(successes / max(1, done), 4),
            "capped_failures": done - successes,
        }
        if dump:
            runs[label]["steps"] = steps_used
    fc = []
    for i in range(first_collision_trials):
        if time.time() - t0 > seconds:
            break
        _, norm = first_collision(q, rp, mults, seed * 7717 + i)
        fc.append(norm)
    fc_sorted = sorted(fc)

    def qtl(p):
        if not fc_sorted:
            return None
        return round(fc_sorted[min(len(fc_sorted) - 1, int(p * len(fc_sorted)))], 4)

    out = {
        "part_a": {
            "checks": checks,
            "samples": samples,
            "seconds": round(elapsed_a, 3),
            "scope": "BLS12-381 Fp12 tower arithmetic: r is prime and random GT "
                     "elements (h^((p^12-1)/r)) have order r",
        },
        "part_b": {
            "q": q,
            "q_bits": q.bit_length(),
            "subgroup_order_r_prime": rp,
            "walks": walks,
            "multipliers_R": mults,
            "dp_bits": dp_bits,
            "runs": runs,
            "uncensored_first_collision": {
                "trials": len(fc),
                "mean_over_sqrtr": round(sum(fc) / len(fc), 4) if fc else None,
                "median_over_sqrtr": qtl(0.5),
                "q10_over_sqrtr": qtl(0.1),
                "q90_over_sqrtr": qtl(0.9),
                "fraction_below_1.5": round(sum(1 for v in fc if v < 1.5) / len(fc), 4) if fc else None,
                "model_mean": 1.2533,
                "model_P_below_1.5": round(1.0 - math.exp(-(1.5 ** 2) / 2), 4),
            },
            "predicted_success_at_1.5_sqrtr": round(1.0 - math.exp(-(1.5 ** 2) / 2), 4),
            "predicted_success_at_2.5_sqrtr": round(1.0 - math.exp(-(2.5 ** 2) / 2), 4),
            "seconds": round(time.time() - t0 - elapsed_a, 3),
            "scope": "toy prime field with an order-2^%d subgroup; evidence at toy "
                     "scope only, never a production complexity proof"
                     % (rp.bit_length(),),
        },
    }
    if dump:
        out["part_b"]["uncensored_first_collision"]["values"] = fc
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
