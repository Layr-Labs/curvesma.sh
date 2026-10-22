"""Toy end-to-end run of the tame-wild class-walk Cheon attack; not evidence of full-size costs."""
import json
import random
import sys


def is_prime(n):
    return n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))


def toy_group(order_lo, d):
    r = order_lo
    while not (is_prime(r) and (r - 1) % (6 * d) == 0 and ((r - 1) // d) % 6 == 0):
        r += 1
    t = 2
    while not is_prime(t * r + 1):
        t += 1
    q = t * r + 1
    h = 2
    while pow(h, t, q) == 1:
        h += 1
    return r, q, pow(h, t, q)


def run(order_lo, d, trials):
    r, q, g = toy_group(order_lo, d)
    n_e = (r - 1) // d                  # |E|, E = d-th powers in F_r^*
    root = 2
    while pow(root, (r - 1) // 3, r) == 1:
        root += 1
    lam = pow(root, (r - 1) // 3, r)    # cube root of unity mod r, in E
    units = [(e * pow(lam, i, r)) % r for i in range(3) for e in (1, r - 1)]
    rng = random.Random(7)
    jumps = [pow(rng.randrange(2, r), d, r) for _ in range(16)]

    def canon(base, e):                 # class {+-lam^i Q} of Q = base^e
        best = min((pow(base, (u * e) % r, q), u) for u in units)
        return best[0], (best[1] * e) % r

    def trail(base, tame):
        e = pow(rng.randrange(2, r), d, r)
        pt = pow(base, e, q)
        for _ in range(200):
            cv, ec = canon(base, e)     # scalar of the canonical point, relative to base
            if cv % 4 == 0:
                return cv, ec, tame
            e = (ec * jumps[(cv // 4) % 16]) % r
            pt = pow(base, e, q)
        return None

    solved = tries = 0
    for s in range(2, 2 + trials):
        big = pow(g, pow(s, d, r), q)
        seen, v = {}, None
        while v is None:
            tries += 1
            res = trail(g, 1) if tries % 2 else trail(big, 0)
            if res is None:
                continue
            cv, ec, tame = res
            if cv in seen and seen[cv][1] != tame:
                et, ew = (ec, seen[cv][0]) if tame else (seen[cv][0], ec)
                v = (et * pow(ew, -1, r)) % r
                assert v == pow(s, d, r)
            seen[cv] = (ec, tame)
        cand = [x for x in range(1, r) if pow(x, d, r) == v and pow(g, x, q) == pow(g, s, q)]
        assert cand == [s]
        solved += 1
    return r, n_e, solved, tries


if __name__ == "__main__":
    params = json.load(sys.stdin)
    r, n_e, solved, tries = run(params["order_lo"], params["d"], params["secrets"])
    print(json.dumps({"order": r, "E_size": n_e, "classes": n_e // 6, "recovered": solved,
                      "trails": tries, "scope": "toy group; checks bookkeeping of classes, scalar ratio and d-th root step only"}))
