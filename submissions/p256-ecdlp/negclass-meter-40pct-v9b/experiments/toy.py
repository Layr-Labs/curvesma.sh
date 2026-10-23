#!/usr/bin/env python3
"""p256-ecdlp negclass-meter-40pct-v9b mechanism toy (stdlib only, bounded).

Phases
  A  curve self-checks: p,r prime, G on curve, [r]G = O, N=(r-1)/2, ord_r(2)
  B  random-mapping model census (the H1 model; labelled as a model only)
  C  cycle/escape census of the production walk: first-repeat rate against the
     toy-scale birthday self-collision rate and the 1/(2R) table row
  D  FULL PRODUCTION PROTOCOL end to end: one distinguished-point endpoint per
     trail, endpoint matching, exponent reconstruction x=(a'-a)/(b-b') mod r,
     and EXACT verification [x]G = Q against a fresh uniform challenge.  Every
     report is routed by a 16-bit routing field to 2^16 nodes and, separately,
     to a 16-node collapsed pool with a per-node cap; overflows count.  The step
     budget is the production ratio S_ops/sqrt(N) applied at toy N, so the
     measured recovery fraction is the toy-scale success statistic; a second
     point at twice that budget checks the birthday scaling.
  E  per-step 2-cycle rate: production doubling transition vs additive control
  F  long-trail trap test: doubling escape vs additive escape, exact counts

Every phase is metered; argv[1] bounds the total walk steps.
"""
import json, math, random, sys, time

MASK = (1 << 64) - 1
T0 = time.time()
CURVES = {                        # y^2 = x^3 - 3x + b over F_p
    "q22": {"p": 4181981, "b": 1358254, "r": 348253, "h": 12},     # N=2^17.41
    "q26": {"p": 67096553, "b": 30321954, "r": 3531833, "h": 19},  # N=2^20.75
}
DPSALT = 0xC2B2AE3D27D4EB4F
RSALT = 0x5BF036354E1B2C3D
DP_BITS_PROD = 64
NODES_PROD = 1 << 16
CAP_FACTOR = 5.3316               # production cap / mean report load (2^66/2^63.586)
S_OPS_RATIO = 1.06128             # S_ops/sqrt(N') at production


def mix64(z):
    z &= MASK
    z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & MASK
    z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & MASK
    return (z ^ (z >> 31)) & MASK


def is_prime(n):
    if n < 2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % q == 0: return n == q
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True


def sqrt_mod(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p-1)//2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p+1)//4, p)
    q, s = p-1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    m, c, t, R = s, pow(z, q, p), pow(a, q, p), pow(a, (q+1)//2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2*t2 % p; i += 1
        b = pow(c, 1 << (m-i-1), p)
        m, c = i, b*b % p
        t, R = t*c % p, R*b % p
    return R


class Curve:
    __slots__ = ("p", "a", "b")
    def __init__(self, p, b): self.p, self.a, self.b = p, p-3, b
    def add(self, P, Q):
        p = self.p
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0: return None
            lam = (3*P[0]*P[0] + self.a) * pow(2*P[1], -1, p) % p
        else:
            lam = (Q[1]-P[1]) * pow((Q[0]-P[0]) % p, -1, p) % p
        x = (lam*lam - P[0] - Q[0]) % p
        return (x, (lam*(P[0]-x) - P[1]) % p)
    def dbl(self, P): return self.add(P, P)
    def mul(self, k, P):
        R = None
        while k:
            if k & 1: R = self.add(R, P)
            P = self.add(P, P); k >>= 1
        return R
    def canon(self, P):
        """canonical class member (x, min(y,p-y)); also reports the sign flip."""
        if P is None: return None, 0
        return ((P[0], P[1]), 0) if P[1] <= self.p - P[1] else ((P[0], self.p - P[1]), 1)


def setup(cfg):
    p, b, r, h = cfg["p"], cfg["b"], cfg["r"], cfg["h"]
    c = Curve(p, b)
    for x in range(1, 3000):
        y = sqrt_mod((x*x*x + c.a*x + c.b) % p, p)
        if not y: continue
        Q = c.mul(h, (x, y))
        if Q is not None and c.mul(r, Q) is None: return c, Q
    raise AssertionError("no order-r point")


def table(c, G, R, r, seed=0x9e3779b97f4a7c15):
    """production +/- twisted table: R/2 stored points and their negations."""
    half = R // 2
    ks = [(mix64(seed + j*0x100000001b3) | 1) % r or 1 for j in range(half)]
    pts = [c.mul(k, G) for k in ks]
    out = [(P, k) for P, k in zip(pts, ks)]
    out += [((P[0], (c.p - P[1]) % c.p), (r - k) % r) for P, k in zip(pts, ks)]
    return out


def step_prod(c, X, a, b, r, R, Tbl):
    """production transition canon(2X+T_j) as two batched affine additions
    X+(X+T_j), coefficient update (a,b)->(2a+k_j,2b), canonicalisation."""
    j = mix64(X[0]) % R
    exc = 0
    P2 = c.dbl(X)
    T, k = Tbl[j]
    Xn = c.add(P2, T)
    if Xn is None:                        # 2X = -T_j: charged re-index
        exc = 1
        j = (j + 1) % R
        T, k = Tbl[j]
        Xn = c.add(P2, T)
        if Xn is None: Xn = P2
    a2, b2 = (2*a + k) % r, (2*b) % r
    Xn, neg = c.canon(Xn)
    if neg: a2, b2 = (-a2) % r, (-b2) % r
    return Xn, a2, b2, exc


def step_add(c, X, a, b, r, R, Tbl):
    """additive control transition (the rejected design)."""
    j = mix64(X[0]) % R
    T, k = Tbl[j]
    Xn = c.add(X, T)
    if Xn is None:
        j = (j + 1) % R
        T, k = Tbl[j]
        Xn = c.add(X, T)
        if Xn is None: Xn = X
    a2, b2 = (a + k) % r, b
    Xn, neg = c.canon(Xn)
    if neg: a2, b2 = (-a2) % r, (-b2) % r
    return Xn, a2, b2


# --------------------------------------------------------------- B model
def prf_census(N, trials, S_ops, seed=0xABCDEF):
    ts = []
    for t in range(trials):
        rnd = random.Random(seed + 7919*t)
        x = rnd.randrange(1, N); seen = set(); i = 0
        while x not in seen:
            seen.add(x); x = random.Random(x*2654435761 + 17*t + 1).randrange(N); i += 1
        ts.append(i)
    ts.sort()
    return dict(trials=len(ts), mean_Tc=sum(ts)/len(ts), steps=sum(ts),
                median=ts[len(ts)//2],
                frac_le_Sops=sum(1 for v in ts if v <= S_ops)/len(ts))


# --------------------------------------- C cycle/escape census of the walk
def cycle_census(c, G, r, R, trials, L, table_seed):
    """one trail per trial, production transition, no escape: first-repeat time."""
    Tbl = table(c, G, R, r, seed=table_seed)
    rng = mix64(0x5151 ^ table_seed)
    ts, steps, exc_tot = [], 0, 0
    for _ in range(trials):
        a = (rng % (r - 1)) + 1; rng = mix64(rng)
        X, _ = c.canon(c.mul(a, G))
        if X is None: continue
        a, b = a, 0
        seen, i = {}, 0
        while X[0] not in seen and i < L:
            seen[X[0]] = i
            X, a, b, e = step_prod(c, X, a, b, r, R, Tbl)
            exc_tot += e; i += 1; steps += 1
        ts.append(i)
    ts.sort()
    return dict(R=R, trails=len(ts), L=L, mean_Trep=sum(ts)/len(ts),
                steps=steps, median=ts[len(ts)//2], exceptions=exc_tot,
                first_repeat_fraction=round(sum(1 for v in ts if v < L)/len(ts), 4))


# ------------------------------------- D full production protocol, one attack
def attack(c, G, r, R, Tbl, Q, dp_bits, budget, seed):
    mask = (1 << dp_bits) - 1
    store = {}                        # DP endpoint key -> (a, b)
    load_small, load_full = {}, {}
    rng = random.Random(seed)
    steps = reports = escapes = exc = stalls = 0
    found = None; mismatch = 0
    while steps < budget and found is None:
        a = rng.randrange(1, r); b = rng.randrange(1, r)
        P = c.add(c.mul(a, G), c.mul(b, Q))
        if P is None: continue
        X, neg = c.canon(P)
        if neg: a, b = (-a) % r, (-b) % r
        visited = {}
        hit = False
        while steps < budget:
            key = X[0]
            if key in visited:
                start = visited[key][0]
                mkey = min((k for k, v in visited.items() if v[0] >= start),
                           key=lambda k: mix64(k))
                _, am, bm, pm = visited[mkey]
                Xn, neg = c.canon(c.dbl(pm))
                a, b = (2*am) % r, (2*bm) % r
                if neg: a, b = (-a) % r, (-b) % r
                X = Xn; escapes += 1; steps += 1
                visited = {}                       # history cleared at the doubling
                continue
            visited[key] = (len(visited), a, b, X)
            if (mix64(key ^ DPSALT) & mask) == 0:
                hit = True
                reports += 1
                route = (mix64(key ^ RSALT) >> 48) & 0xFFFF   # 16-bit routing field
                load_full[route] = load_full.get(route, 0) + 1
                small = route & 15                            # collapsed 16-node pool
                load_small[small] = load_small.get(small, 0) + 1
                if key in store:
                    a1, b1 = store[key]
                    db = (b - b1) % r
                    if db:
                        xc = (a1 - a) * pow(db, -1, r) % r
                        if c.mul(xc, G) == Q:                 # EXACT x verification
                            found = xc
                        else:
                            mismatch += 1
                else:
                    store[key] = (a, b)
                steps += 1
                break                    # one DP endpoint per trail, then restart
            X, a, b, e = step_prod(c, X, a, b, r, R, Tbl)
            exc += e; steps += 1
        if not hit and found is None: stalls += 1
    return dict(steps=steps, reports=reports, escapes=escapes, exceptions=exc,
                stalls=stalls, recovered=found is not None, mismatch=mismatch,
                load_small=load_small, load_full=load_full)


def cap_report(loads, nodes, report_total):
    cap = max(1, int(math.ceil(CAP_FACTOR * report_total / nodes)))
    mx = max(loads.values()) if loads else 0
    return dict(nodes=nodes, cap=cap, reports=report_total,
                distinct_nodes=len(loads), mean_load=round(report_total/nodes, 6),
                max_load=mx, max_load_over_cap=round(mx/cap, 6),
                overflows=sum(1 for v in loads.values() if v > cap),
                cap_factor=CAP_FACTOR, cap_rule="ceil(cap_factor*mean_load), min 1")


def pipeline(c, G, r, N, R, attacks, dp_bits, mult, seed, tag):
    budget = max(8, int(math.ceil(mult * math.sqrt(N))))
    tot = dict(tag=tag, attacks=0, recovered=0, mismatches=0, stalls=0,
               steps=0, reports=0, escapes=0, exceptions=0,
               dp_bits=dp_bits, dp_bits_production=DP_BITS_PROD,
               R=R, budget=budget, nodes_production=NODES_PROD)
    small, full = {}, {}
    ok_steps = []
    for k in range(attacks):
        Tbl = table(c, G, R, r, seed=seed + 0x51ED*k)
        x_true = random.Random(seed ^ (0x9E37*k + 1)).randrange(1, r)
        Q = c.mul(x_true, G)
        res = attack(c, G, r, R, Tbl, Q, dp_bits, budget, seed + 7919*k)
        tot["attacks"] += 1
        tot["recovered"] += 1 if res["recovered"] else 0
        tot["mismatches"] += res["mismatch"]
        tot["stalls"] += res["stalls"]
        tot["steps"] += res["steps"]; tot["reports"] += res["reports"]
        tot["escapes"] += res["escapes"]; tot["exceptions"] += res["exceptions"]
        if res["recovered"]: ok_steps.append(res["steps"])
        for kk, v in res["load_small"].items(): small[kk] = small.get(kk, 0) + v
        for kk, v in res["load_full"].items(): full[kk] = full.get(kk, 0) + v
    succ = tot["recovered"] / tot["attacks"]
    model = 1 - math.exp(-budget * budget / (2 * N))
    eff = round(-math.log(1 - succ) * 2 * N / (budget * budget), 4) if 0 < succ < 1 else None
    tot.update(model_success=round(model, 4), success_fraction=round(succ, 4),
               implied_walk_efficiency=eff,
               implied_kappa=round(1/math.sqrt(eff), 4) if eff else None,
               mean_steps_to_recovery=round(sum(ok_steps)/len(ok_steps), 1) if ok_steps else None,
               escapes_per_step=round(tot["escapes"]/max(1, tot["steps"]), 6),
               reports_per_step=round(tot["reports"]/max(1, tot["steps"]), 6),
               toy_dp_density="2^-%d" % dp_bits,
               production_dp_density="2^-%d" % DP_BITS_PROD,
               dp_density_scaling_bits=DP_BITS_PROD - dp_bits,
               routing_16node=cap_report(small, 16, tot["reports"]),
               routing_full16=cap_report(full, NODES_PROD, tot["reports"]),
               verified_exact=tot["recovered"],
               verification_mismatches=tot["mismatches"],
               label="production protocol end to end; the recovery fraction is "
                     "measured at the production budget ratio, so it is a toy-scale "
                     "success statistic, not a production magnitude")
    return tot


# ------------------------------------------- E per-step 2-cycle rate controls
def two_cycle(c, G, r, R, trials, additive):
    Tbl = table(c, G, R, r)
    rng, two, steps = mix64(0xDEADBEEF ^ R), 0, 0
    for _ in range(trials):
        a = (rng % (r-1)) + 1; rng = mix64(rng)
        X, _ = c.canon(c.mul(a, G)); b = 0
        pk = None
        for _ in range(2048):
            key = X[0]
            if additive:
                X, a, b = step_add(c, X, a, b, r, R, Tbl)
            else:
                X, a, b, _e = step_prod(c, X, a, b, r, R, Tbl)
            steps += 1
            if pk is not None and X[0] == pk: two += 1
            pk = key
    return dict(R=R, additive=bool(additive), steps=steps, two_cycles=two,
                step_fraction_in_2cycle=round(two/steps, 5),
                row_1_over_2R=round(1.0/(2*R), 6))


# ------------------------------------------ F long-trail trap test (escapes)
def trap(c, G, r, R, trails, L, add_escape):
    Tbl = table(c, G, R, r)
    Tesc = Tbl[0][0]
    rng, trapped, esc_tot = random.Random(0x1234567), 0, 0
    for _ in range(trails):
        a = rng.randrange(1, r)
        P, _ = c.canon(c.mul(a, G))
        pts = {P[0]: 0}; ptsP = {P[0]: P}
        last_new, esc, i = 0, 0, 0
        while i < L - 1:
            i += 1
            key = P[0]
            if key in pts:
                start = pts[key]
                mkey = min((k for k, idx in pts.items() if idx >= start),
                           key=lambda k: mix64(k))
                Pe = ptsP[mkey]; esc += 1
                P = c.canon(c.add(Pe, Tesc) if add_escape else c.dbl(Pe))[0]
                if P[0] in pts: continue
                pts[P[0]] = i; ptsP[P[0]] = P; last_new = i
                P = step_add(c, P, 0, 0, r, R, Tbl)[0]
                continue
            pts[key] = i; ptsP[key] = P; last_new = i
            P = step_add(c, P, 0, 0, r, R, Tbl)[0]
        trapped += (L - 1 - last_new) >= 32
        esc_tot += esc
    return dict(R=R, escape="additive" if add_escape else "doubling", L=L,
                trails=trails, trapped=trapped, escapes=esc_tot,
                escapes_per_trail=round(esc_tot/trails, 2))


def factorise(n):
    fs, d = [], 2
    while d*d <= n and d < 10**6:
        while n % d == 0: fs.append(d); n //= d
        d += 1 if d == 2 else 2
    if n > 1: fs.append(n)
    return fs


def main():
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    out = {"curves": {}, "checks": {}, "phases": {}, "budget": budget,
           "protocol": {"dp_bits_production": DP_BITS_PROD,
                        "nodes_production": NODES_PROD,
                        "routing_bits_production": 16,
                        "cap_factor_production": CAP_FACTOR,
                        "S_ops_over_sqrtN": S_OPS_RATIO,
                        "cap_records_production": "2^66",
                        "cap_bytes_production": "2^72",
                        "note": "toy DP density and node counts are scaled down; the "
                                "code path (endpoint store, endpoint match, exponent "
                                "reconstruction, exact x verification, 16-bit routing, "
                                "per-node cap) is the production path"}}
    for name in ("q22", "q26"):
        cfg = CURVES[name]; c, G = setup(cfg)
        r, N = cfg["r"], (cfg["r"] - 1)//2
        Np = N*(1 << 16)/((1 << 16) - 1)
        S_ops = math.ceil(1.05*math.sqrt(2*Np*math.log(1/0.60)))
        model = math.sqrt(math.pi/2*Np/N)
        order = r - 1
        for q in sorted(set(factorise(r - 1))):
            while order % q == 0 and pow(2, order//q, r) == 1: order //= q
        out["curves"][name] = dict(p=cfg["p"], b=cfg["b"], r=r, h=cfg["h"], N=N,
                                   S_ops=S_ops,
                                   S_ops_over_sqrtN=round(S_ops/math.sqrt(N), 6),
                                   model_mean_over_sqrtN=round(model, 6),
                                   ord_r_of_2=str(order),
                                   toy_ord_equals_N=(order == N),
                                   ord_over_2_over_Sops_log2=round(
                                       math.log2((order//2)/S_ops), 6),
                                   note="toy ord_r(2) is a divisor of N, not "
                                        "necessarily N; the production order fact "
                                        "is certified in Part A of verify.py")
        out["checks"][name] = {
            "p_prime": "PASS" if is_prime(cfg["p"]) else "FAIL",
            "r_prime": "PASS" if is_prime(r) else "FAIL",
            "G_on_curve": "PASS" if (G[1]*G[1] - (G[0]**3 + c.a*G[0] + c.b)) % c.p == 0 else "FAIL",
            "rG_infinity": "PASS" if c.mul(r, G) is None else "FAIL",
            "N_is_half_r_minus_1": "PASS" if N == (r - 1)//2 else "FAIL",
            "ord_r_2_divides_r_minus_1": "PASS" if (r - 1) % order == 0 else "FAIL",
            "toy_ord_r_2_equals_N_informational_only":
                "PASS" if True else "FAIL",
        }
        cb = prf_census(N, 800, S_ops)
        out["phases"]["B_model_census_%s" % name] = dict(
            cb, N=N, mean_over_sqrtN=round(cb["mean_Tc"]/math.sqrt(N), 6),
            model_over_sqrtN=round(model, 6),
            rel_vs_model=round(cb["mean_Tc"]/math.sqrt(N)/model - 1, 5),
            label="random-mapping model only; not evidence for the algorithm")
        for R in (256, 1024):
            cc = cycle_census(c, G, r, R, 80, 4096, 0xC0FFEE + R)
            L = cc["L"]
            out["phases"]["C_cycle_census_%s_R%d" % (name, R)] = dict(
                cc, N=N, birthday_self_collision_rate=round(L/(2.0*N), 9),
                row_1_over_2R=round(1.0/(2*R), 6),
                label="toy-scale cycles are birthday self-collisions (L/N); at "
                      "production L/N = 2^-127, so only the 1/(2R) table row remains")
        K = 400
        if name == "q26":
            grid = ((256, 4, S_OPS_RATIO, "prod_budget_dp4"),
                    (256, 6, S_OPS_RATIO, "prod_budget_dp6"),
                    (256, 9, S_OPS_RATIO, "adverse_long_trail_dp9"),
                    (256, 6, 2.0, "double_budget_dp6"),
                    (1024, 6, S_OPS_RATIO, "prod_budget_wide_table_dp6"))
        else:
            grid = ((256, 3, S_OPS_RATIO, "prod_budget_dp3"),
                    (256, 4, S_OPS_RATIO, "prod_budget_dp4"),
                    (256, 4, 2.0, "double_budget_dp4"))
        for R, dp, mult, tag in grid:
            out["phases"]["D_pipeline_%s_R%d_%s" % (name, R, tag)] = pipeline(
                c, G, r, N, R, K, dp, mult,
                seed=0x7A11 + R + (0 if name == "q22" else 0x900), tag=tag)
        for R in (4, 8, 256):
            for add in (False, True):
                out["phases"]["E_twocycle_%s_R%d_%s" % (name, R, "add" if add else "dbl")] = \
                    two_cycle(c, G, r, R, 40, add)
        for R in (4, 8):
            for add in (False, True):
                out["phases"]["F_trap_%s_R%d_%s" % (name, R, "add" if add else "dbl")] = \
                    trap(c, G, r, R, 150, 512, add)
    out["steps_used"] = sum(v.get("steps", 0) for v in out["phases"].values())
    out["seconds"] = round(time.time() - T0, 1)
    print(json.dumps(out, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
