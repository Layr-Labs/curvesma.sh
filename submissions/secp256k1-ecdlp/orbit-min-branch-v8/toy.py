"""The q25 toy of proof.md section 8 (part C of experiment.py) and the curve class both use: the toy in discrete-log
form and as a curve pipeline that computes each incidence's six leaves from one batched inverse of D."""
import math
import time


class Curve:
    """Short Weierstrass y^2 = x^3 + 7 over F_q with prime group order r."""

    def __init__(self, q, order):
        """q: field prime; order: group order. returns: None."""
        self.q, self.r = q, order

    def add(self, P, Q):
        """Affine addition. P, Q: points or None (infinity). returns: P + Q."""
        q = self.q
        if P is None:
            return Q
        if Q is None:
            return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % q == 0:
                return None
            l = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
        else:
            l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
        x = (l * l - P[0] - Q[0]) % q
        return (x, (l * (P[0] - x) - P[1]) % q)

    def neg(self, P):
        """Negation. P: point. returns: -P."""
        return None if P is None else (P[0], (-P[1]) % self.q)

    def mul(self, k, P):
        """Double-and-add. k: integer; P: point. returns: [k]P."""
        k %= self.r
        acc = None
        while k:
            if k & 1:
                acc = self.add(acc, P)
            P = self.add(P, P)
            k >>= 1
        return acc


def toy_setup(rng):
    """The q25 toy curve of package F (q = 33556867, r = 33563419, j = 0): generator, beta, lambda.
    rng: Random. returns: (Curve, G, beta, lambda)."""
    q, order = 33556867, 33563419
    E = Curve(q, order)
    while True:
        x = rng.randrange(q)
        y = pow(x ** 3 + 7, (q + 1) // 4, q)
        if (y * y - x ** 3 - 7) % q == 0:
            G = (x, y)
            break
    assert E.mul(order, G) is None
    cr = lambda m: [c for g in range(2, 50) for c in [pow(g, (m - 1) // 3, m)] if c != 1][0]
    for b in (cr(q), cr(q) ** 2 % q):
        for l in (cr(order), cr(order) ** 2 % order):
            if E.mul(l, G) == (G[0] * b % q, G[1]):
                return E, G, b, l
    raise AssertionError("no beta/lambda pair")


def toy_trial_logs(rng, rt, lam, nodes, q1, q2):
    """One toy trial in discrete-log form (proof.md 8): class of [t]G is t^6 mod r; leaf phi^t(X) + sY has logarithm
    lam^t x + s y and Q coefficient lam^t b.

    rng: Random; rt: order; lam: cube root of 1; nodes, q1, q2: shape. returns: (scalars, success, cross-node
    success, colliding pairs).
    """
    xi = rng.randrange(rt)
    sc = []
    for g in range(nodes):
        bg = rng.randrange(rt)
        xs = [(rng.randrange(rt), bg) for _ in range(q1)]
        ys = [rng.randrange(1, rt) for _ in range(q2)]
        sc.append((xs, ys))
    units = [pow(lam, t, rt) for t in range(3)]
    keys, live = [], []
    for g, (xs, ys) in enumerate(sc):
        for jx, (a, b) in enumerate(xs):
            tx = (a + b * xi) % rt
            if tx == 0 or any((c - e * w * tx) % rt == 0 for c in ys for w in units for e in (1, -1)):
                continue  # X = O or D = 0 (Y = +-phi^t(X)): the batch is a defect row (proof.md 6.4)
            live.append((g, b, tx, ys))
            keys += [pow((w * tx + e * c) % rt, 6, rt) for c in ys for w in units for e in (1, -1)]
    groups = {}
    if len(set(keys)) < len(keys):
        for g, b, tx, ys in live:
            for c in ys:
                for w in units:
                    for e in (1, -1):
                        t = (w * tx + e * c) % rt
                        groups.setdefault(pow(t, 6, rt), []).append((g, w * b % rt, t))
    succ = succ_x = False
    pairs = 0
    for v in groups.values():
        if len(v) < 2:
            continue
        for u_ in range(len(v)):
            for w_ in range(u_ + 1, len(v)):
                (ga, ba, ta), (gb, bb, tb) = v[u_], v[w_]
                pairs += 1
                u = tb * pow(ta, -1, rt) % rt
                if (bb - u * ba) % rt:
                    succ = True
                    succ_x |= ga != gb
    return (xi, sc), succ, succ_x, pairs


TOY_TAB = {}


def toy_table(E, G):
    """Fixed-base table of G: rows of [d * 256^w]G. E: curve; G: generator. returns: list of rows (cached)."""
    if G not in TOY_TAB:
        rows, P = [], G
        for _ in range(4):
            row = [None]
            for _ in range(255):
                row.append(E.add(row[-1], P))
            rows.append(row)
            P = E.add(row[-1], P)
        TOY_TAB[G] = rows
    return TOY_TAB[G]


def toy_fmul(E, tab, k):
    """[k]G by the fixed-base table. E: curve; tab: toy_table rows; k: scalar. returns: point."""
    k %= E.r
    acc, w = None, 0
    while k:
        acc = E.add(acc, tab[w][k & 255])
        k >>= 8
        w += 1
    return acc


def toy_trial_curve(E, G, beta, lam, scalars, nodes, q1, q2):
    """The same trial on the toy curve: six leaves per incidence from one inverse of D (proof.md 3.1), keys
    min(x, bx, b^2x), equal-key pairs resolved with the points and checked.

    E, G, beta, lam: toy curve; scalars: from toy_trial_logs; nodes, q1, q2: shape. returns: (success, pairs).
    """
    xi, sc = scalars
    q, rt = E.q, E.r
    tab = toy_table(E, G)
    fm = lambda k: toy_fmul(E, tab, k)
    Q = fm(xi)
    b2 = beta * beta % q
    groups = {}
    for g, (xs, ys) in enumerate(sc):
        Ypts = [fm(c) for c in ys]
        for jx, (a, b) in enumerate(xs):
            X = fm(a + b * xi)
            if X is None or any((Y[1] * Y[1] - X[1] * X[1]) % q == 0 for Y in Ypts):
                continue
            Ds = [(Y[1] * Y[1] - X[1] * X[1]) % q for Y in Ypts]
            pre = [1]
            for d in Ds:
                pre.append(pre[-1] * d % q)
            R = pow(pre[-1], -1, q)
            for iy in range(q2 - 1, -1, -1):
                inv = R * pre[iy] % q
                R = R * Ds[iy] % q
                Y = Ypts[iy]
                u, v = (Y[1] - X[1]) * inv % q, (Y[1] + X[1]) * inv % q
                for t in range(3):
                    bt = pow(beta, t, q)
                    Qt = (Y[0] * Y[0] + bt * X[0] * Y[0] + bt * bt * X[0] * X[0]) % q
                    for s, w in ((1, u), (-1, v)):
                        x = (Qt * Qt * w * w - bt * X[0] - Y[0]) % q
                        key = min(x, x * beta % q, x * b2 % q)
                        groups.setdefault(key, []).append((g, t, s, a, b, ys[iy]))
    succ, pairs = False, 0
    for v in groups.values():
        for u_ in range(len(v)):
            for w_ in range(u_ + 1, len(v)):
                pairs += 1
                pts = []
                for (g, t, s, a, b, c) in (v[u_], v[w_]):
                    lt = pow(lam, t, rt)
                    P = E.add(E.add(fm(lt * a), E.mul(lt * b % rt, Q)), fm(s * c))
                    pts.append((P, (lt * a + s * c) % rt, lt * b % rt))
                (PA, aA, bA), (PB, aB, bB) = pts
                for k in range(3):
                    for e in (1, -1):
                        img = (PA[0] * pow(beta, k, q) % q, PA[1] if e == 1 else (-PA[1]) % q)
                        if img == PB:
                            u = e * pow(lam, k, rt) % rt
                            den = (bB - u * bA) % rt
                            if den:
                                cand = (u * aA - aB) * pow(den, -1, rt) % rt
                                assert E.mul(cand, G) == Q
                                succ = True
    return succ, pairs


def mu_of(nodes, n0, rt):
    """Cross-node collision mean |Gamma| P1. nodes, n0: shape; rt: order. returns: float."""
    s = nodes * n0
    return (s * s - nodes * n0 * n0) / 2 * (6 * (rt - 1) + 1) / rt ** 2


def part_c(rng, deadline, trials, curve_trials, min_multi, min_curve):
    """Toy success near the production mu, and curve/log equivalence.

    rng: Random; deadline: experiment.py's budget function; trials, curve_trials: requested counts; min_multi,
    min_curve: counts that always run. returns: dict."""
    E, G, beta, lam = toy_setup(rng)
    rt = E.r
    res = dict(q=E.q, r=rt)
    nodes, q1, q2 = 25, 4, 4
    ok = okx = done = pairs = 0
    dl = deadline(0.70)
    while done < trials and (done < min_multi or time.time() < dl):
        _, s1, s2, pr = toy_trial_logs(rng, rt, lam, nodes, q1, q2)
        ok += s1
        okx += s2
        pairs += pr
        done += 1
    n0 = 6 * q1 * q2
    S = nodes * n0
    mu_all = S * (S - 1) / 2 * (6 * (rt - 1) + 1) / rt ** 2
    mu = mu_of(nodes, n0, rt)
    ph = ok / max(done, 1)
    p1t = (6 * (rt - 1) + 1) / rt ** 2
    gam = (S * S - nodes * n0 * n0) / 2
    cs = gam * 12 * (q1 + q2) * S * (p1t * p1t + 36 / (rt - 1) ** 2)
    res["multi_node_grid"] = dict(nodes=nodes, grid=[q1, q2], leaves_per_incidence=6, leaves=S, trials=done,
                                  success=round(ph, 5), se=round(math.sqrt(ph * (1 - ph) / max(done, 1)), 5),
                                  cross_node_success=round(okx / max(done, 1), 5),
                                  mean_colliding_pairs=round(pairs / max(done, 1), 5), mu_cross=round(mu, 6),
                                  mu_all_pairs=round(mu_all, 6), poisson=round(1 - math.exp(-mu), 5),
                                  chen_stein_b1_b2_at_this_scale=round(cs, 5))
    agree = n = curve_succ = 0
    dl = deadline(1.0)
    while n < curve_trials and (n < min_curve or time.time() < dl):
        scal, s1, _, pr = toy_trial_logs(rng, rt, lam, nodes, q1, q2)
        s2, pr2 = toy_trial_curve(E, G, beta, lam, scal, nodes, q1, q2)
        agree += (s1 == s2) and (pr == pr2)
        curve_succ += s2
        n += 1
    res["curve_pipeline_vs_logs"] = dict(trials=n, identical_outcomes=agree, curve_successes=curve_succ)
    return res
