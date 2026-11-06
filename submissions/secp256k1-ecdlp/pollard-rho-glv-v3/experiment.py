"""Toy Pollard-rho over GLV automorphism classes, v3 walk (proof.md v3).

Faithful miniature of the production algorithm in proof.md (v3):

* the walk is the production step function G of proof.md Section 3,
  implemented ONCE (make_walk) and shared by the attack and census
  modes: from the canonical representative R with partition index
  a = j(R), step to canon(R + S_a); if that class again has index a,
  bump once to canon(R' + S_{a+1}); an attempt whose entry satisfies
  S = +-R (the production exceptional branch: the sum is O or 2R)
  retries on consecutive table indices (at most 3 attempts per phase)
  and otherwise abandons the trail. Mode walk="plain" runs the v2 step
  function F (single generic addition, no bump, no retries) for
  comparison and census.
* scalar bookkeeping is the production SPLIT-COORDINATE implementation
  of proof.md Section 3.4: (A, B) are held as A = A1 + A2*lambda,
  B = B1 + B2*lambda over the integers (lambda^2 = -1 - lambda),
  updated by split increments and by the mu-action formulas, and
  REBALANCED every 36 steps by evaluating the split components
  (A1 + A2*lambda mod r) and re-splitting by Babai rounding against a
  Gauss-reduced basis of the decomposition lattice. A shadow copy
  updated by ordinary modular arithmetic serves only as a test oracle:
  the split state is asserted equal to it immediately after every step
  (before any rebalancing) and the envelope is asserted with it;
  reported scalars are the split-evaluated ones and the collision
  solve uses them.
* every canonicalisation uses the production single-multiplication
  method (beta^2*x = -(x + beta*x) mod p) and is cross-checked against
  brute-force enumeration of all six transforms.

Modes:
* "attack": distinguished-point search as in proof.md; reports success
  in the total-step metric (every step counted, repeated work
  included) and the distinct-visited-class metric, bump counts,
  affected (exceptional) steps, trapped trails and component maxima.
* "census": exact functional-graph census of the chosen step function
  on all classes: cyclic states by cycle length, 2-cycle states, and
  the v2 prediction N_classes/(6J) for the plain walk's 2-cycle states.

Scope: validates the class-walk machinery, the bump's destruction of
the systematic 2-cycle family, the split-coordinate schedule, and
visited-class collision statistics at toy scale. It does not validate
production instruction costs or memory bounds.
"""
import json
import math
import random
import sys
from fractions import Fraction


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


def gauss_reduce(b1, b2):
    """Gauss-reduce a 2D lattice basis (exact integer arithmetic)."""
    def norm2(v):
        return v[0] * v[0] + v[1] * v[1]
    for _ in range(1000):
        if norm2(b2) < norm2(b1):
            b1, b2 = b2, b1
        m = (b1[0] * b2[0] + b1[1] * b2[1])
        d = norm2(b1)
        q = (2 * m + d) // (2 * d) if m >= 0 else -((-2 * m + d) // (2 * d))
        if q == 0:
            return b1, b2
        b2 = (b2[0] - q * b1[0], b2[1] - q * b1[1])
    return b1, b2  # iteration cap; split() asserts validate the result


class Splitter:
    """GLV splitting k = k1 + k2*lambda (mod r) with small components, by
    Babai rounding against a Gauss-reduced basis of the decomposition
    lattice Lambda = {(a, b) : a + b*lambda = 0 (mod r)}."""

    def __init__(self, r, lam):
        self.r = r
        self.lam = lam
        v1, v2 = gauss_reduce((r, 0), ((-lam) % r, 1))
        assert (v1[0] + v1[1] * lam) % r == 0
        assert (v2[0] + v2[1] * lam) % r == 0
        self.v1, self.v2 = v1, v2
        det = v1[0] * v2[1] - v1[1] * v2[0]
        assert abs(det) == r
        self.det = det
        # Babai residual bound: components of the residual are at most
        # (|v1|_inf + |v2|_inf); take the safe envelope below.
        self.bound = max(abs(v1[0]), abs(v1[1]), abs(v2[0]), abs(v2[1])) * 2

    def split(self, k):
        k %= self.r
        # (c1, c2) = (k,0) - q1*v1 - q2*v2 with q the rounded exact
        # solution of [v1 v2] q = (k, 0).
        det = self.det
        q1 = Fraction(k * self.v2[1], det)
        q2 = Fraction(-k * self.v1[1], det)
        q1 = math.floor(q1 + Fraction(1, 2))
        q2 = math.floor(q2 + Fraction(1, 2))
        c1 = k - q1 * self.v1[0] - q2 * self.v2[0]
        c2 = -q1 * self.v1[1] - q2 * self.v2[1]
        assert (c1 + c2 * self.lam) % self.r == k
        assert abs(c1) <= self.bound and abs(c2) <= self.bound, (c1, c2)
        return c1, c2


def main():
    cfg = json.load(sys.stdin)
    p = int(cfg["prime"])
    trials = int(cfg.get("trials", 1))
    seed = int(cfg["seed"])
    mode = cfg.get("mode", "attack")
    walk = cfg.get("walk", "bumped")
    J = int(cfg.get("J", 32))
    theta_bits = int(cfg.get("theta_bits", 2))
    detect_cycles = bool(cfg.get("detect_cycles", False))
    assert mode in ("attack", "census")
    assert walk in ("bumped", "plain")
    assert p % 6 == 1

    n_group = count_points(p)
    fac = factorize(n_group)
    r = max(fac)
    assert fac[r] == 1 and r % 3 == 1
    cof = n_group // r

    P = None
    x = 0
    while P is None:
        rhs = (pow(x, 3, p) + 7) % p
        if pow(rhs, (p - 1) // 2, p) == 1:
            cand = ec_mul(cof, (x, tonelli(rhs, p)), p)
            if cand is not None and ec_mul(r, cand, p) is None:
                P = cand
        x += 1

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
    assert beta2 == (-1 - beta) % p
    assert (lam2 + lam + 1) % r == 0

    n_classes = 1 + (r - 1) // 6
    transforms = [(1, 1, 1), (r - 1, 1, -1), (lam, beta, 1),
                  (r - lam, beta, -1), (lam2, beta2, 1), (r - lam2, beta2, -1)]

    def canon(Pt):
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

    crosschecks = [0]

    def canon_fast(Pt):
        """Production canonicalisation: ONE multiplication by beta;
        beta^2*x = -(x + beta*x) mod p; sign from Pt's own y.
        Cross-checked against canon() on every call."""
        x0, y0 = Pt
        bx = beta * x0 % p
        b2x = (-(x0 + bx)) % p
        xs = (x0, bx, b2x)
        taus = (1, lam, lam2)
        k = min(range(3), key=lambda i: xs[i])
        if y0 > (p - 1) // 2:
            rep, tau = (xs[k], p - y0), (-taus[k]) % r
        else:
            rep, tau = (xs[k], y0), taus[k]
        ref = canon(Pt)
        crosschecks[0] += 1
        assert (rep, tau) == ref, (Pt, (rep, tau), ref)
        return rep, tau

    splitter = Splitter(r, lam)

    def mu_apply(c1, c2, tau):
        """(c1, c2) -> components of tau*(c1 + c2*lambda), exact integers,
        by the proof's formulas (lambda^2 = -1 - lambda)."""
        if tau == 1:
            return c1, c2
        if tau == r - 1:
            return -c1, -c2
        if tau == lam:            # lambda*(c1 + c2 l) = -c2 + (c1 - c2) l
            return -c2, c1 - c2
        if tau == r - lam:
            return c2, c2 - c1
        if tau == lam2:           # lambda^2*(...) = (c2 - c1) - c1 l
            return c2 - c1, -c1
        if tau == r - lam2:
            return c1 - c2, c1
        raise AssertionError(tau)

    def build_table(rng, Q):
        U = [rng.randrange(r) for _ in range(J)]
        V = [rng.randrange(r) for _ in range(J)]
        S = [ec_add(ec_mul(U[j], P, p), ec_mul(V[j], Q, p), p)
             for j in range(J)]
        assert all(s is not None for s in S)  # S_j != O, as in production
        Us = [splitter.split(u) for u in U]
        Vs = [splitter.split(v) for v in V]
        return U, V, S, Us, Vs

    def neg(Pt):
        return (Pt[0], (-Pt[1]) % p)

    def make_walk(S):
        """The production step function of proof.md Section 3 on points.
        Returns walk(R) -> (R_next, phases, affected, adds) or None to
        abandon; phases is the ordered list of (table index, alignment)
        actually applied. walk='plain' gives the v2 function F."""
        def walk_fn(R):
            a = R[0] % J
            if walk == "plain":
                T = ec_add(R, S[a], p)
                if T is None:
                    return None
                R2, tau = canon_fast(T)
                return R2, [(a, tau)], 0, 1
            affected = 0
            phases = []
            adds = 0
            # phase 1: indices a, a+1, a+2
            R2 = None
            for k in range(3):
                idx = (a + k) % J
                adds += 1
                if S[idx] == R or S[idx] == neg(R):
                    affected = 1  # production exception: S = +-R
                    continue
                T = ec_add(R, S[idx], p)
                R2, tau = canon_fast(T)
                phases.append((idx, tau))
                break
            if R2 is None:
                return None
            if R2[0] % J != a:
                return R2, phases, affected, adds
            # phase 2 (bump): indices a+1, a+2, a+3
            R3 = None
            for k2 in range(1, 4):
                idx2 = (a + k2) % J
                adds += 1
                if S[idx2] == R2 or S[idx2] == neg(R2):
                    affected = 1
                    continue
                T2 = ec_add(R2, S[idx2], p)
                R3, tau2 = canon_fast(T2)
                phases.append((idx2, tau2))
                break
            if R3 is None:
                return None
            return R3, phases, affected, adds
        return walk_fn

    def census():
        rng = random.Random(seed)
        secret = rng.randrange(r)
        Q = ec_mul(secret, P, p)
        _, _, S, _, _ = build_table(rng, Q)
        walk_fn = make_walk(S)
        # Enumerate all classes by walking the subgroup.
        reps = []
        seen = set()
        R = P
        for _ in range(r):
            if R is None:
                break
            c, _ = canon(R)
            if c not in seen:
                seen.add(c)
                reps.append(c)
            R = ec_add(R, P, p)
        assert len(reps) == n_classes - 1  # identity class excluded

        nxt = {}
        for R in reps:
            if walk == "plain":
                T = ec_add(R, S[R[0] % J], p)
                v = None if T is None else canon(T)[0]
            else:
                out = walk_fn(R)
                v = None if out is None else out[0]
            if v is not None and v in seen:
                nxt[R] = v
        from collections import deque
        indeg = {R: 0 for R in reps}
        for u, v in nxt.items():
            indeg[v] += 1
        dq = deque([u for u in reps if indeg[u] == 0])
        removed = set()
        while dq:
            u = dq.popleft()
            removed.add(u)
            v = nxt.get(u)
            if v is not None:
                indeg[v] -= 1
                if indeg[v] == 0:
                    dq.append(v)
        cyclic = [u for u in reps if u not in removed]
        seen2 = set()
        by_len = {}
        two = 0
        for u in cyclic:
            if u in seen2:
                continue
            cur, cyc = u, []
            while cur not in seen2:
                seen2.add(cur)
                cyc.append(cur)
                cur = nxt[cur]
            by_len[len(cyc)] = by_len.get(len(cyc), 0) + len(cyc)
            if len(cyc) == 2:
                two += 2
        return {
            "mode": "census", "prime": p, "subgroup_order": r,
            "classes": n_classes, "walk": walk, "J": J,
            "cyclic_states": len(cyclic),
            "two_cycle_states": two,
            "predicted_two_cycle_states_plain": (n_classes - 1) / (6 * J),
            "cycle_states_by_length": dict(sorted(by_len.items())),
            "canon_crosschecks": crosschecks[0],
        }

    if mode == "census":
        print(json.dumps(census()))
        return

    # ---------------- attack mode ----------------
    theta_inv = 1 << theta_bits
    L_MAX = 16 * theta_inv
    budget = math.ceil(math.sqrt(2 * n_classes * 0.519195))
    cap = 8 * budget
    REBAL_EVERY = 36

    rng_master = random.Random(seed)
    succ_steps = 0
    succ_visited = 0
    succ_cap = 0
    steps_list = []
    visited_list = []
    trapped = 0
    wasted_detected = 0
    trap_waste_production = 0
    bad_collisions = 0
    total_bumps = 0
    total_affected = 0
    total_additions = 0
    max_component = 0
    split_checks = 0

    for _ in range(trials):
        rng = random.Random(rng_master.randrange(1 << 62))
        secret = rng.randrange(r)
        Qg = ec_mul(secret, P, p)
        U, V, S, Us, Vs = build_table(rng, Qg)
        walk_fn = make_walk(S)

        table = {}
        visited = set()
        total_steps = 0
        found = None
        while found is None and total_steps <= cap:
            a0 = rng.randrange(r)
            b0 = rng.randrange(r)
            R = ec_add(ec_mul(a0, P, p), ec_mul(b0, Qg, p), p)
            if R is None:
                # degenerate start in the identity class: resample, as
                # in production (probability ~1/r per start)
                continue
            A1, A2 = splitter.split(a0)
            B1, B2 = splitter.split(b0)
            shadow = [a0, b0]
            R, tau = canon_fast(R)
            A1, A2 = mu_apply(A1, A2, tau)
            B1, B2 = mu_apply(B1, B2, tau)
            shadow[0] = shadow[0] * tau % r
            shadow[1] = shadow[1] * tau % r
            steps = 0
            since_rebal = 0
            envelope = 2 * splitter.bound  # start split, then one mu-map
            trail_seen = {}
            while True:
                visited.add(R)
                if detect_cycles:
                    # measuring instrument only (production relies on
                    # L_MAX): a repeated state within this trail means
                    # the trail is on a cycle; record incidence and the
                    # waste production would have burned (L_MAX - entry).
                    if R in trail_seen:
                        trapped += 1
                        wasted_detected += steps - trail_seen[R]
                        trap_waste_production += L_MAX - trail_seen[R]
                        break
                    trail_seen[R] = steps
                # split-state invariant checks (production schedule)
                split_checks += 1
                assert (A1 + A2 * lam) % r == shadow[0]
                assert (B1 + B2 * lam) % r == shadow[1]
                mc = max(abs(A1), abs(A2), abs(B1), abs(B2))
                if mc > max_component:
                    max_component = mc
                assert mc <= envelope, (mc, envelope)
                if R[0] % theta_inv == 0:  # distinguished: report, end
                    A_ev = (A1 + A2 * lam) % r
                    B_ev = (B1 + B2 * lam) % r
                    if R in table:
                        A1t, B1t = table[R]
                        if B1t != B_ev:
                            xx = (A1t - A_ev) * pow(
                                (B_ev - B1t) % r, -1, r) % r
                            if ec_mul(xx, P, p) == Qg:
                                found = xx
                            else:
                                bad_collisions += 1
                        else:
                            bad_collisions += 1
                    else:
                        table[R] = (A_ev, B_ev)
                    break
                if steps >= L_MAX:
                    trapped += 1
                    break
                out = walk_fn(R)
                if out is None:
                    trapped += 1
                    break
                (R, phases, na, adds) = out
                total_affected += na
                total_additions += adds
                if len(phases) == 2:
                    total_bumps += 1
                for (idx, tauv) in phases:
                    A1, A2 = A1 + Us[idx][0], A2 + Us[idx][1]
                    B1, B2 = B1 + Vs[idx][0], B2 + Vs[idx][1]
                    shadow[0] = (shadow[0] + U[idx]) % r
                    shadow[1] = (shadow[1] + V[idx]) % r
                    A1, A2 = mu_apply(A1, A2, tauv)
                    B1, B2 = mu_apply(B1, B2, tauv)
                    shadow[0] = shadow[0] * tauv % r
                    shadow[1] = shadow[1] * tauv % r
                steps += 1
                total_steps += 1
                since_rebal += 1
                # per full step: increment add (M + B) then mu-map (x2);
                # a bumped step does this twice: M -> 4M + 6B; 4M + 8B
                # is the safe envelope recurrence for both cases.
                envelope = 4 * envelope + 8 * splitter.bound
                # verify the post-step split state BEFORE rebalancing
                split_checks += 1
                assert (A1 + A2 * lam) % r == shadow[0]
                assert (B1 + B2 * lam) % r == shadow[1]
                mc = max(abs(A1), abs(A2), abs(B1), abs(B2))
                assert mc <= envelope, (mc, envelope)
                if since_rebal >= REBAL_EVERY:
                    # rebalance by evaluating the SPLIT components and
                    # re-splitting (production Section 3.4); the oracle
                    # only checks the evaluated values agree.
                    A1, A2 = splitter.split((A1 + A2 * lam) % r)
                    B1, B2 = splitter.split((B1 + B2 * lam) % r)
                    since_rebal = 0
                    envelope = splitter.bound
                if total_steps > cap:
                    break
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
        "scope": "toy curve only; validates the v3 class-walk (production "
                 "step function with bump and bounded retries), "
                 "split-coordinate bookkeeping with the rebalancing "
                 "schedule, and visited-class collision statistics at "
                 "toy scale, not production instruction costs or "
                 "memory bounds",
        "prime": p,
        "subgroup_order": r,
        "automorphism": "glv",
        "walk": walk,
        "classes": n_classes,
        "J": J,
        "theta_bits": theta_bits,
        "trials": trials,
        "budget": budget,
        "success_rate_steps_metric": succ_steps / trials,
        "success_rate_visited_metric": succ_visited / trials,
        "successes_within_cap": succ_cap,
        "mean_steps_successes": mean(steps_list),
        "mean_visited_successes": mean(visited_list),
        "predicted_mean_visited": math.sqrt(math.pi * n_classes / 2),
        "trapped_trails": trapped,
        "detect_cycles": detect_cycles,
        "wasted_steps_detected_cycles": wasted_detected,
        "trap_waste_if_production_cap": trap_waste_production,
        "bad_collisions": bad_collisions,
        "bumps": total_bumps,
        "affected_steps": total_affected,
        "total_additions": total_additions,
        "split_component_bound": splitter.bound,
        "max_split_component": max_component,
        "split_state_checks": split_checks,
        "canon_crosschecks": crosschecks[0],
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
