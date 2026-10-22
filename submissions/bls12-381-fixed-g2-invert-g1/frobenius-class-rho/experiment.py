"""Toy end-to-end run of the class-walk scheme: prime-order subgroup of Z_q^* with an order-12 exponent
automorphism y -> y^lam, state kept un-canonicalised, conjugated jump tables, distinguished points.
Only a small-group sanity check; it says nothing about production costs."""
import json
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


def make_tables(r, q, g, lam, z, rng, tab):
    H = [pow(lam, i, r) for i in range(12)]
    Hinv = [pow(h, -1, r) for h in H]
    jumps = [(rng.randrange(r), rng.randrange(r)) for _ in range(tab)]
    M = [pow(g, a, q) * pow(z, b, q) % q for a, b in jumps]
    Mc = [[pow(Mj, hi, q) for hi in Hinv] for Mj in M]
    ec = [[(a * hi % r, b * hi % r) for hi in Hinv] for a, b in jumps]
    return H, Mc, ec


def walk(r, q, y, al, be, tabs, tab, dpbits, cap, phase=0, alpha=64):
    """One trail from y = g^al z^be. Each step: distinguished test first, then the fixed-period cycle rule.
    The trail counts its own steps t (starting at `phase`). Every alpha steps it saves the key as checkpoint C.
    A step whose key equals C (tested before C is overwritten) detects a cycle of length <= alpha; W is the least
    key since C. The trail then runs on to the class with key W and leaves it with jump index j + 1."""
    H, Mc, ec = tabs
    t, steps, ckpt, wmin, seek, events = phase, 0, None, None, None, 0
    while steps < cap:
        imgs = [pow(y, h, q) for h in H]
        idx = min(range(12), key=imgs.__getitem__)
        key = imgs[idx]
        if (key * 40503 >> 5) % (1 << dpbits) == 0:
            return key, (al * H[idx] % r, be * H[idx] % r), steps, events
        j = (key * 2654435761 >> 7) % tab
        if seek is not None:
            if key == seek:
                j, seek, ckpt, wmin = (j + 1) % tab, None, None, None
        elif ckpt is not None and key == ckpt:
            events += 1
            seek = wmin
            if key == seek:
                j, seek, ckpt, wmin = (j + 1) % tab, None, None, None
        elif t % alpha == 0:
            ckpt, wmin = key, key
        elif ckpt is not None:
            wmin = min(wmin, key)
        steps += 1
        t += 1
        y = y * Mc[j][idx] % q
        al, be = (al + ec[j][idx][0]) % r, (be + ec[j][idx][1]) % r
    return None, None, steps, events


def solve(r, q, g, lam, z, rng, tab=64, dpbits=3):
    tabs = make_tables(r, q, g, lam, z, rng, tab)
    classes = (r - 1) // 12
    budget = int(3.0 * classes ** 0.5) + 200
    seen, steps = {}, 0
    while steps < budget:
        al, be = rng.randrange(r), rng.randrange(r)
        y = pow(g, al, q) * pow(z, be, q) % q
        key, rec, n, _ = walk(r, q, y, al, be, tabs, tab, dpbits, 64 * (1 << dpbits), rng.randrange(64))
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
    """Two trails started at the same class (same group element y, different exponents and different phase of
    the period-64 counter) must finish at the same distinguished point. Cycles are frequent here (tab = 4)."""
    z = pow(g, rng.randrange(1, r), q)
    tabs = make_tables(r, q, g, lam, z, rng, tab)
    agree = differ = capped = events = 0
    for _ in range(trials):
        y = pow(g, rng.randrange(r), q)
        out = []
        for ph in (rng.randrange(64), rng.randrange(64)):
            al, be = rng.randrange(r), rng.randrange(r)
            out.append(walk(r, q, y, al, be, tabs, tab, dpbits, 4000, ph))
        events += out[0][3] + out[1][3]
        if out[0][0] is None and out[1][0] is None:
            capped += 1
        elif out[0][0] == out[1][0]:
            agree += 1
        else:
            differ += 1
    return agree, differ, capped, events


def main():
    spec = json.load(sys.stdin)
    r, q, g, lam = setup(spec["order"])
    rng = random.Random(spec.get("seed", 1))
    ok = bad = 0
    for _ in range(spec.get("trials", 20)):
        x = rng.randrange(1, r)
        got, _ = solve(r, q, g, lam, pow(g, x, q), rng)
        if got is None:
            continue
        ok += got == x
        bad += got != x
    agree, differ, capped, events = merge_test(r, q, g, lam, rng, spec.get("merge_trials", 200))
    print(json.dumps({"subgroup_order": r, "modulus": q, "automorphism_order": 12,
                      "recovered": ok, "wrong": bad, "merge_agree": agree, "merge_differ": differ,
                      "merge_both_capped": capped, "cycle_events_in_merge_test": events,
                      "scope": "toy group only"}))
    assert bad == 0 and differ == 0


main()
