"""Toy Pollard-rho with automorphism equivalence classes on a small curve.

Faithful miniature of the production argument in proof.md: distinguished-point
trails over equivalence classes of the automorphism group ({+-1} or
{+-1, +-lambda, +-lambda^2} for the GLV endomorphism (x,y)->(beta*x,y) of
y^2 = x^3 + 7), canonical representatives, per-step scalar updates by the
applied automorphism, length-capped trails, and a fixed budget
t = ceil(sqrt(2 * N_classes * 0.562119)) for which the birthday bound predicts
success probability >= 0.43.

Two success metrics are reported. The birthday bound governs the number of
distinct classes VISITED; total STEPS additionally include trail-start points,
post-collision distinguished-point overshoot and trapped-trail waste. At
production scale those corrections are relatively < 2^-20 (proof.md), so the
metrics coincide; at toy scale they do not, and both are shown. Scope: this
validates the class-walk machinery (canonicalisation, scalar bookkeeping,
collision solving) and the visited-class collision statistics at toy scale.
It does not validate production instruction costs or memory bounds.
"""
import json
import math
import random
import sys


def ec_add(Pt, Qt, p):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    x1, y1 = Pt
    x2, y2 = Qt
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if Pt == Qt:
        lam = (3 * x1 * x1) * pow(2 * y1 % p, -1, p) % p
    else:
        lam = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_mul(k, Pt, p):
    R = None
    A = Pt
    while k:
        if k & 1:
            R = ec_add(R, A, p)
        A = ec_add(A, A, p)
        k >>= 1
    return R


def tonelli(n, p):
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r0 = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 1, t * t % p
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r0 = i, b * b % p, t * b * b % p, r0 * b % p
    return r0


def count_points(p):
    cnt = 1
    for x in range(p):
        rhs = (pow(x, 3, p) + 7) % p
        ls = pow(rhs, (p - 1) // 2, p)
        if ls == 1:
            cnt += 2
        elif ls == 0:
            cnt += 1
    return cnt


def factorize(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def main():
    cfg = json.load(sys.stdin)
    p = int(cfg["prime"])
    trials = int(cfg["trials"])
    seed = int(cfg["seed"])
    mode = cfg["automorphism"]
    assert mode in ("glv", "negation")
    assert p % 6 == 1

    n_group = count_points(p)
    fac = factorize(n_group)
    r = max(fac)  # largest prime factor: subgroup order
    assert fac[r] == 1 and r % 3 == 1
    cof = n_group // r

    # Generator of the order-r subgroup.
    P = None
    x = 0
    while P is None:
        rhs = (pow(x, 3, p) + 7) % p
        if pow(rhs, (p - 1) // 2, p) == 1:
            cand = ec_mul(cof, (x, tonelli(rhs, p)), p)
            if cand is not None and ec_mul(r, cand, p) is None:
                P = cand
        x += 1

    # Automorphism data: beta a primitive cube root of unity mod p,
    # lambda its eigenvalue on <P> (lambda^2 + lambda + 1 = 0 mod r).
    beta = None
    g = 2
    while beta is None:
        c = pow(g, (p - 1) // 3, p)
        if c != 1:
            beta = c
        g += 1
    beta2 = beta * beta % p
    lam = None
    for cand in range(2, r):
        if (cand * cand + cand + 1) % r == 0:
            if ec_mul(cand, P, p) == (beta * P[0] % p, P[1]):
                lam = cand
                break
    assert lam is not None
    lam2 = lam * lam % r

    if mode == "glv":
        n_classes = 1 + (r - 1) // 6
        transforms = [(1, 1, 1), (r - 1, 1, -1), (lam, beta, 1),
                      (r - lam, beta, -1), (lam2, beta2, 1), (r - lam2, beta2, -1)]
    else:
        n_classes = (r + 1) // 2
        transforms = [(1, 1, 1), (r - 1, 1, -1)]

    def canon(Pt):
        """(canonical representative of the class of Pt, scalar tau with
        rep = [tau]Pt)."""
        x0, y0 = Pt
        best = None
        best_tau = 1
        for tau, bx, sy in transforms:
            cx = bx * x0 % p
            cy = y0 if sy == 1 else (-y0) % p
            if best is None or (cx, cy) < best:
                best = (cx, cy)
                best_tau = tau
        return best, best_tau

    DP_BITS = 2
    theta_inv = 1 << DP_BITS
    L_MAX = 16 * theta_inv
    J = 32
    budget = math.ceil(math.sqrt(2 * n_classes * 0.562119))
    cap = 8 * budget

    rng_master = random.Random(seed)
    succ_steps = 0
    succ_visited = 0
    succ_cap = 0
    steps_list = []
    visited_list = []
    trapped = 0
    bad_collisions = 0

    for _ in range(trials):
        rng = random.Random(rng_master.randrange(1 << 62))
        secret = rng.randrange(r)
        Q = ec_mul(secret, P, p)
        # Per-instance adding points S_j = [u_j]P + [v_j]Q (as in production).
        U = [rng.randrange(r) for _ in range(J)]
        V = [rng.randrange(r) for _ in range(J)]
        S = [ec_add(ec_mul(U[j], P, p), ec_mul(V[j], Q, p), p) for j in range(J)]

        table = {}
        visited = set()
        total_steps = 0
        found = None
        while found is None and total_steps <= cap:
            a0 = rng.randrange(r)
            b0 = rng.randrange(r)
            R = ec_add(ec_mul(a0, P, p), ec_mul(b0, Q, p), p)
            A, B = a0, b0
            R, tau = canon(R)
            A, B = tau * A % r, tau * B % r
            steps = 0
            # A trail walks until its FIRST distinguished class (reports and
            # ends), or L_MAX steps (abandoned: trapped in a DP-free cycle).
            while True:
                visited.add(R)
                if R[0] % theta_inv == 0:  # distinguished class: report, end
                    if R in table:
                        A1, B1 = table[R]
                        if B1 != B:
                            xx = (A1 - A) * pow((B - B1) % r, -1, r) % r
                            if ec_mul(xx, P, p) == Q:
                                found = xx
                            else:
                                bad_collisions += 1
                        else:
                            bad_collisions += 1
                    else:
                        table[R] = (A, B)
                    break
                if steps >= L_MAX:
                    trapped += 1
                    break
                j = R[0] % J
                T = ec_add(R, S[j], p)
                A = (A + U[j]) % r
                B = (B + V[j]) % r
                R, tau = canon(T)
                A, B = tau * A % r, tau * B % r
                steps += 1
                total_steps += 1
                if total_steps > cap:
                    break
            # end trail
        if found is not None:
            assert found == secret
            succ_cap += 1
            steps_list.append(total_steps)
            visited_list.append(len(visited))
            if total_steps <= budget:
                succ_steps += 1
            if len(visited) <= budget:
                succ_visited += 1

    def mean(v):
        return sum(v) / len(v) if v else None

    out = {
        "scope": "toy curve only; validates class-walk machinery and "
                 "visited-class collision statistics at toy scale, not "
                 "production instruction costs or memory bounds",
        "prime": p,
        "subgroup_order": r,
        "automorphism": mode,
        "classes": n_classes,
        "trials": trials,
        "budget": budget,
        "success_rate_steps_metric": succ_steps / trials,
        "success_rate_visited_metric": succ_visited / trials,
        "successes_within_cap": succ_cap,
        "mean_steps_successes": mean(steps_list),
        "mean_visited_successes": mean(visited_list),
        "predicted_mean_visited": math.sqrt(math.pi * n_classes / 2),
        "trapped_trails": trapped,
        "bad_collisions": bad_collisions,
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
