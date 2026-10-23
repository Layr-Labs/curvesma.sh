"""Toy run of the class walk of proof.md section 3 on a prime-order subgroup of Z_q^* with an order-12 exponent
automorphism y -> y^lam (standing in for Frobenius). Same structure as the production walk: un-canonicalised state with
conjugated jump tables, class key = least of the 12 images, jump index and distinguished-point test from a mix of the
key, period-64 checkpoint cycle detection, escape from the least-key member of the cycle by
SQUARING the state (exponents doubled). It also runs the earlier escape (ordinary table step j+1) as a control.
Measures: end-to-end recovery, merge of trails that share a class, escapes that fall back into the cycle, and the
steps-to-first-collision ratio against the random-mapping value sqrt(pi N'/2). Toy scale only: no production cost claim."""
import json
import math
import random
import sys


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def setup(order):
    r = order
    while not (r % 12 == 1 and is_prime(r)):
        r -= 1
    m = 2
    while not is_prime(m * r + 1):
        m += 1
    q = m * r + 1
    g = next(pow(h, m, q) for h in range(2, q) if pow(h, m, q) != 1)
    lam = next(pow(h, (r - 1) // 12, r) for h in range(2, r) if len({pow(pow(h, (r - 1) // 12, r), i, r) for i in range(12)}) == 12)
    return r, q, g, lam


MASK = (1 << 64) - 1


def fin(x):
    """splitmix64 finalizer."""
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9 & MASK
    x = (x ^ (x >> 27)) * 0x94D049BB133111EB & MASK
    return x ^ (x >> 31)


def mix(key):
    """Index and distinguished-test bits from the class key (one limb here, six in production): a multiply-xorshift
    chain, then two finalizers. m1 gives the jump index, m2 the distinguished test."""
    s = key * 0x9E3779B97F4A7C15 & MASK
    s ^= s >> 32
    return fin(s), fin(s + 0x9E3779B97F4A7C15 & MASK)


def make_tables(r, q, g, lam, z, rng, tab):
    H = [pow(lam, i, r) for i in range(12)]
    Hinv = [pow(h, -1, r) for h in H]
    jumps = [(rng.randrange(r), rng.randrange(r)) for _ in range(tab)]
    M = [pow(g, a, q) * pow(z, b, q) % q for a, b in jumps]
    Mc = [[pow(Mj, hi, q) for hi in Hinv] for Mj in M]
    ec = [[(a * hi % r, b * hi % r) for hi in Hinv] for a, b in jumps]
    return H, Mc, ec


def walk(r, q, y, al, be, tabs, tab, dpbits, cap, phase=0, alpha=64, escape="square", visit=None):
    """One trail from y = g^al z^be (tab a power of two). Key = least image; the jump index j and the distinguished test
    come from a mix of the key (top bits of m1; top dpbits bits of m2 zero). Each step: distinguished test first, then the fixed-period cycle rule: every alpha
    steps save the key as checkpoint C; a key equal to C (tested before C is overwritten) detects a cycle of length <= alpha,
    W is the least key since C; the trail then runs on to the class with key W and leaves it by `escape`: "square" replaces
    the table step by y -> y^2 with exponents doubled, "jplus1" uses table index j + 1. Returns (key, exps, steps, events,
    falls); a fall is an escape from a least-key member this trail has escaped from before."""
    H, Mc, ec = tabs
    lt = tab.bit_length() - 1
    t, steps, ckpt, wmin, seek, events, falls, escaped = phase, 0, None, None, None, 0, 0, set()
    while steps < cap:
        imgs = [pow(y, h, q) for h in H]
        idx = min(range(12), key=imgs.__getitem__)
        key = imgs[idx]
        if visit is not None:
            visit.append(key)
        m1, m2 = mix(key)
        if m2 >> (64 - dpbits) == 0:
            return key, (al * H[idx] % r, be * H[idx] % r), steps, events, falls
        j = m1 >> (64 - lt)
        esc = False
        if seek is not None:
            esc = key == seek
        elif ckpt is not None and key == ckpt:
            events += 1
            seek = wmin
            esc = key == seek
        elif t % alpha == 0:
            ckpt, wmin = key, key
        elif ckpt is not None:
            wmin = min(wmin, key)
        if esc:
            falls += key in escaped
            escaped.add(key)
            seek = ckpt = wmin = None
        steps += 1
        t += 1
        if esc and escape == "square":
            y, al, be = y * y % q, 2 * al % r, 2 * be % r
        else:
            if esc:
                j = (j + 1) % tab
            y = y * Mc[j][idx] % q
            al, be = (al + ec[j][idx][0]) % r, (be + ec[j][idx][1]) % r
    return None, None, steps, events, falls


def solve(r, q, g, lam, z, rng, tab=16, dpbits=3):
    tabs = make_tables(r, q, g, lam, z, rng, tab)
    classes = (r - 1) // 12
    budget = int(3.0 * classes ** 0.5) + 200
    seen, steps = {}, 0
    while steps < budget:
        al, be = rng.randrange(r), rng.randrange(r)
        y = pow(g, al, q) * pow(z, be, q) % q
        key, rec, n, _, _ = walk(r, q, y, al, be, tabs, tab, dpbits, 64 * (1 << dpbits), rng.randrange(64))
        steps += n
        if key is None:
            continue
        ca, cb = rec
        if key in seen:
            a2, b2 = seen[key]
            if b2 != cb:
                return (a2 - ca) * pow(cb - b2, -1, r) % r, steps
        else:
            seen[key] = (ca, cb)
    return None, steps


def merge_test(r, q, g, lam, rng, trials, tab=4, dpbits=6):
    """Two trails started at the same class (same group element y, different exponents and different phase of the
    period-64 counter) must finish at the same distinguished point. Cycles are frequent here (tab = 4). Run in a larger
    group (about 8e6 classes) so that random cycles longer than 64 steps, which the rule does not claim to handle, are rare."""
    z = pow(g, rng.randrange(1, r), q)
    tabs = make_tables(r, q, g, lam, z, rng, tab)
    agree = differ = capped = events = 0
    for _ in range(trials):
        y = pow(g, rng.randrange(r), q)
        out = []
        for ph in (rng.randrange(64), rng.randrange(64)):
            al, be = rng.randrange(r), rng.randrange(r)
            out.append(walk(r, q, y, al, be, tabs, tab, dpbits, 800, ph))
        events += out[0][3] + out[1][3]
        if out[0][0] is None and out[1][0] is None:
            capped += 1
        elif out[0][0] == out[1][0]:
            agree += 1
        else:
            differ += 1
    return agree, differ, capped, events


def fall_test(big, rng, walks, tab=2, steps=200):
    """Walks of `steps` steps without distinguished points in a larger group (about 8e6 classes), so that repeats
    of an escape from the same least-key member come from the rule and not from random self-collisions of the trail.
    A 2-entry table makes 2-cycles (rate 1/24 per step) frequent. Counts escapes and falls for the squaring rule and for
    the earlier table-step rule as a control."""
    r, q, g, lam = big
    out = {}
    for rule in ("square", "jplus1"):
        rnd = random.Random(rng.randrange(1 << 30))
        tabs = make_tables(r, q, g, lam, pow(g, rnd.randrange(1, r), q), rnd, tab)
        ev = fl = 0
        for _ in range(walks):
            al, be = rnd.randrange(r), rnd.randrange(r)
            _, _, _, e, f = walk(r, q, pow(g, al, q), al, be, tabs, tab, 40, steps, rnd.randrange(64), escape=rule)
            ev, fl = ev + e, fl + f
        out[rule] = (ev, fl)
    return out


def ratio_test(r, q, g, lam, rng, trials, tab=16, dpbits=3):
    """Number of distinct classes visited until the first class visited by two different trails (trails restart after a
    distinguished point or after 512 steps), divided by the random-mapping value sqrt(pi N'/2); also the same ratio in
    executed steps, which at this toy size is inflated by cycle waste (period 64 against sqrt(N') = 90). Returns
    (productive mean, 95% half-width, executed mean, reference)."""
    nprime = (r - 1) // 12 + 1
    ref = math.sqrt(math.pi * nprime / 2)
    prod, execd = [], []
    for _ in range(trials):
        z = pow(g, rng.randrange(1, r), q)
        tabs = make_tables(r, q, g, lam, z, rng, tab)
        owner, total, hit = {}, 0, False
        for trail in range(10 ** 4):
            al, be = rng.randrange(r), rng.randrange(r)
            vis = []
            n = walk(r, q, pow(g, al, q) * pow(z, be, q) % q, al, be, tabs, tab, dpbits, 512,
                     rng.randrange(64), visit=vis)[2]
            for i, key in enumerate(vis):
                if owner.setdefault(key, trail) != trail:
                    total, hit = total + i + 1, True
                    break
            if hit:
                break
            total += n + 1
        prod.append(len(owner) / ref)
        execd.append(total / ref)
    mean = sum(prod) / trials
    sd = math.sqrt(sum((v - mean) ** 2 for v in prod) / (trials - 1))
    return mean, 1.96 * sd / math.sqrt(trials), sum(execd) / trials, ref


def control_ratio(r, rng, trials, dpbits=3):
    """The same statistic for an ideal random mapping on N' points with the same distinguished-point rate and trail
    cap: it is above 1 at this size because trails are short and end at distinguished points, so it, not 1, is what the
    walk's ratio is compared with. Returns the mean and its 95% half-width."""
    n = (r - 1) // 12 + 1
    ref = math.sqrt(math.pi * n / 2)
    vals = []
    for _ in range(trials):
        f = [rng.randrange(n) for _ in range(n)]
        owner, done = {}, False
        for trail in range(10 ** 4):
            x, steps = rng.randrange(n), 0
            while True:
                if owner.setdefault(x, trail) != trail:
                    done = True
                    break
                if rng.getrandbits(dpbits) == 0 or steps >= 512:
                    break
                x, steps = f[x], steps + 1
            if done:
                break
        vals.append(len(owner) / ref)
    mean = sum(vals) / trials
    sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (trials - 1))
    return mean, 1.96 * sd / math.sqrt(trials)


def main():
    spec = json.load(sys.stdin)
    r, q, g, lam = setup(spec["order"])
    big = setup(10 ** 8)
    rng = random.Random(spec.get("seed", 1))
    ok = bad = 0
    for _ in range(spec.get("trials", 3)):
        x = rng.randrange(1, r)
        got, _ = solve(r, q, g, lam, pow(g, x, q), rng)
        if got is None:
            continue
        ok += got == x
        bad += got != x
    agree, differ, capped, events = merge_test(*big, rng, spec.get("merge_trials", 12))
    falls = fall_test(big, rng, spec.get("fall_walks", 30))
    mean, half, execd, ref = ratio_test(r, q, g, lam, rng, spec.get("ratio_trials", 100))
    cmean, chalf = control_ratio(r, rng, spec.get("ratio_trials", 100))
    print(json.dumps({"subgroup_order": r, "modulus": q, "automorphism_order": 12,
                      "recovered": ok, "wrong": bad, "merge_agree": agree, "merge_differ": differ,
                      "merge_both_capped": capped, "cycle_events_in_merge_test": events,
                      "escapes_and_falls_square": falls["square"], "escapes_and_falls_jplus1_control": falls["jplus1"],
                      "collision_ratio_mean": round(mean, 3), "collision_ratio_ci95_halfwidth": round(half, 3),
                      "collision_ratio_executed_steps_mean": round(execd, 3),
                      "random_mapping_control_ratio_mean": round(cmean, 3),
                      "random_mapping_control_ci95_halfwidth": round(chalf, 3),
                      "collision_ratio_reference_steps": round(ref, 1), "scope": "toy group only"}))
    assert bad == 0 and differ == 0


if __name__ == "__main__":
    main()
