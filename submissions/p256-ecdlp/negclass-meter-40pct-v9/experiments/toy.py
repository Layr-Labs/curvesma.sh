#!/usr/bin/env python3
"""p256-ecdlp negclass-meter-40pct-v9 mechanism toy (stdlib only, bounded).

Phases
  A  curve self-checks: p,r prime, G on curve, [r]G = O, N=(r-1)/2
  B  H1 census on the random-mapping ensemble that the ledger postulates:
     E[T_c] vs sqrt(pi N'/2) and P[T_c <= S_ops] (the success constant)
  C  concrete class-walk census (production transition canon(2X+T_j)):
     reported as an upper-bound control on T_c, not used for pricing
  D  per-step 2-cycle behaviour on the +/- twisted table: production
     doubling step vs the additive step, against the 1/(2R) escape row
  E  2^16 -> 16 storage-node DP routing with per-node caps (overflow count)
  F  long-trail trap test: same walk, doubling escape vs additive escape

Every phase is metered; --budget caps total walk steps.
"""
import json, math, random, sys, time

MASK = (1 << 64) - 1
T0 = time.time()
CURVES = {                       # y^2 = x^3 - 3x + b over F_p, found at build
    "q22": {"p": 4181981, "b": 1358254, "r": 348253, "h": 12},   # N=2^17.4
    "q26": {"p": 67096553, "b": 30321954, "r": 3531833, "h": 19},  # N=2^20.8
}

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
        if P is None: return P
        return P if P[1] <= self.p - P[1] else (P[0], self.p - P[1])

def setup(cfg):
    p, b, r, h = cfg["p"], cfg["b"], cfg["r"], cfg["h"]
    c = Curve(p, b)
    for x in range(1, 3000):
        y = sqrt_mod((x*x*x + c.a*x + c.b) % p, p)
        if not y: continue
        Q = c.mul(h, (x, y))
        if Q is not None and c.mul(r, Q) is None: return c, Q
    raise AssertionError("no order-r point")

def table(c, G, R, r, twisted, seed=0x9e3779b97f4a7c15):
    """R point offsets.  twisted=False: R independent [k_j]G.
    twisted=True: the +/- table of the construction -- R/2 stored points and
    their negations, so a negated index reproduces the additive 2-cycle."""
    if not twisted:
        return [c.mul((mix64(seed + j*0x100000001b3) | 1) % r or 1, G) for j in range(R)]
    half = R//2
    base = [c.mul((mix64(seed + j*0x100000001b3) | 1) % r or 1, G) for j in range(half)]
    return base + [(P[0], (c.p - P[1]) % c.p) if P else None for P in base]

def step(c, P, R, Tbl):      return c.canon(c.add(c.dbl(P), Tbl[mix64(P[0]) % R]))
def step_add(c, P, R, Tbl):  return c.canon(c.add(P, Tbl[mix64(P[0]) % R]))

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

def walk_census(c, G, r, R, trials, S_ops, tables=4, seed=0xC0FFEE):
    """The table is a fresh algorithm coin, so the census is pooled over
    several independently drawn tables as well as over starting points."""
    ts, cyc, steps, dps, per_table = [], [], 0, [], []
    for s in range(tables):
        Tbl = table(c, G, R, r, twisted=False, seed=seed + 0x51ED*s)
        rng = mix64(seed ^ R ^ s)
        tt = []
        for _ in range(trials):
            P = c.mul((rng % (r-1)) + 1, G); rng = mix64(rng)
            seen, i = {}, 0
            while P[0] not in seen and i < (1 << 22):
                seen[P[0]] = i
                if (mix64(P[0]) & 31) == 0: dps.append(P[0])
                P = step(c, P, R, Tbl); i += 1; steps += 1
            ts.append(i); cyc.append(i - seen[P[0]]); tt.append(i)
        per_table.append(round(sum(tt)/len(tt)/math.sqrt(N_GLOBAL), 3))
    ts.sort(); cyc.sort()
    return dict(R=R, trials=len(ts), tables=tables, mean_Tc=sum(ts)/len(ts),
                steps=steps, median=ts[len(ts)//2], median_cycle=cyc[len(cyc)//2],
                frac_le_Sops=sum(1 for v in ts if v <= S_ops)/len(ts),
                per_table_means_over_sqrtN=per_table), dps

def two_cycle(c, G, r, R, trials, additive):
    Tbl = table(c, G, R, r, twisted=True)
    rng, two, steps = mix64(0xDEADBEEF ^ R), 0, 0
    for _ in range(trials):
        P = c.mul((rng % (r-1)) + 1, G); rng = mix64(rng)
        pk = None
        for _ in range(2048):
            key = P[0]
            Q = step_add(c, P, R, Tbl) if additive else step(c, P, R, Tbl)
            steps += 1
            if pk is not None and Q[0] == pk: two += 1
            pk, P = key, Q
    return dict(R=R, additive=bool(additive), steps=steps, two_cycles=two,
                step_fraction_in_2cycle=round(two/steps, 5),
                row_1_over_2R=round(1.0/(2*R), 6))

def routing(dps, nodes=16):
    loads = [0]*nodes
    for key in dps: loads[mix64(key ^ 0x5bf03635) >> 60] += 1
    cap = math.ceil(2*len(dps)/nodes)
    return dict(nodes=nodes, dps=len(dps), cap=cap, loads=loads,
                overflows=sum(1 for L in loads if L > cap),
                max_load_over_cap=round(max(loads)/cap, 6))

def trap(c, G, r, R, trails, L, add_escape):
    """The reviewer-reported trap: the walk is the additive transition (the
    one whose 2-cycles are structural); only the escape action differs."""
    Tbl = table(c, G, R, r, twisted=True)
    Tesc = Tbl[0]
    rng, trapped, esc_tot = random.Random(0x1234567), 0, 0
    for _ in range(trails):
        P = c.mul(rng.randrange(1, r), G)
        pts = {P[0]: 0}; ptsP = {P[0]: P}
        last_new, esc, i = 0, 0, 0
        while i < L-1:
            i += 1
            key = P[0]
            if key in pts:
                start = pts[key]
                mkey = min((k for k, idx in pts.items() if idx >= start),
                           key=lambda k: mix64(k))
                Pe = ptsP[mkey]; esc += 1
                P = c.canon(c.add(Pe, Tesc)) if add_escape else c.canon(c.dbl(Pe))
                if P[0] in pts: continue
                pts[P[0]] = i; ptsP[P[0]] = P; last_new = i
                P = step_add(c, P, R, Tbl)
                continue
            pts[key] = i; ptsP[key] = P; last_new = i
            P = step_add(c, P, R, Tbl)
        trapped += (L-1-last_new) >= 32
        esc_tot += esc
    return dict(R=R, escape="additive" if add_escape else "doubling", L=L,
                trails=trails, trapped=trapped, escapes=esc_tot,
                escapes_per_trail=round(esc_tot/trails, 2))

def main():
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    global N_GLOBAL
    out = {"curves": {}, "checks": {}, "phases": {}, "budget": budget}
    for name in ("q22", "q26"):
        cfg = CURVES[name]; c, G = setup(cfg)
        r, N = cfg["r"], (cfg["r"]-1)//2
        Np = N*(1 << 16)/((1 << 16)-1)
        S_ops = math.ceil(1.05*math.sqrt(2*Np*math.log(1/0.60)))
        model = math.sqrt(math.pi/2*Np/N)
        out["curves"][name] = dict(p=cfg["p"], b=cfg["b"], r=r, h=cfg["h"], N=N,
                                   S_ops=S_ops,
                                   S_ops_over_sqrtN=round(S_ops/math.sqrt(N), 6),
                                   model_mean_over_sqrtN=round(model, 6))
        out["checks"][name] = {
            "p_prime": "PASS" if is_prime(cfg["p"]) else "FAIL",
            "r_prime": "PASS" if is_prime(r) else "FAIL",
            "G_on_curve": "PASS" if (G[1]*G[1] - (G[0]**3 + c.a*G[0] + c.b)) % c.p == 0 else "FAIL",
            "rG_infinity": "PASS" if c.mul(r, G) is None else "FAIL",
            "N_is_half_r_minus_1": "PASS" if N == (r-1)//2 else "FAIL",
        }
        cb = prf_census(N, 1200, S_ops)                      # B
        out["phases"]["B_prf_census_%s" % name] = dict(
            cb, N=N, mean_over_sqrtN=round(cb["mean_Tc"]/math.sqrt(N), 6),
            model_over_sqrtN=round(model, 6),
            rel_vs_model=round(cb["mean_Tc"]/math.sqrt(N)/model - 1, 5))
        N_GLOBAL = N
        for R, trials in ((256, 150), (2048, 150)):          # C
            res, dps = walk_census(c, G, r, R, trials, S_ops)
            out["phases"]["C_walk_census_%s_R%d" % (name, R)] = dict(
                res, N=N, mean_over_sqrtN=round(res["mean_Tc"]/math.sqrt(N), 6),
                model_over_sqrtN=round(model, 6),
                rel_vs_model=round(res["mean_Tc"]/math.sqrt(N)/model - 1, 5))
            if name == "q26":
                out["phases"]["E_routing_R%d" % R] = routing(dps)
        for R in (4, 8, 256):                                # D
            for add in (False, True):
                out["phases"]["D_twocycle_%s_R%d_%s" % (name, R, "add" if add else "dbl")] = \
                    two_cycle(c, G, r, R, 40, add)
        for R in (4, 8):                                     # F
            for add in (False, True):
                out["phases"]["F_trap_%s_R%d_%s" % (name, R, "add" if add else "dbl")] = \
                    trap(c, G, r, R, 150, 512, add)
    out["steps_used"] = sum(v.get("steps", 0) for v in out["phases"].values())
    out["seconds"] = round(time.time()-T0, 1)
    print(json.dumps(out, indent=1, sort_keys=True))

if __name__ == "__main__":
    main()
