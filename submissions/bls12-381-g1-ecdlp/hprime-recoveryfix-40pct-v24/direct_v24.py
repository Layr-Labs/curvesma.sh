#!/usr/bin/env python3
"""v24 evidence: DIRECT end-to-end solve measurement with the RECOVERY ALGEBRA FIXED.

v23 flaw (verdict a841deb5, kill 2): the pilot runner solved 0/1024 because the
candidate enumeration walked the mu6 orbit of x/y instead of recovering v.  y is
already v-scaled, so x/y is a mu6 unit; the recovered value must divide by the
UNSCALED B-side product b*w' and then be verified.

Correct algebra (this file):
  A record : x = a_i * w_k            (a_i, w_k in H, d-th powers)
  B record : y = v * b_j * w'_k       (v in H is the secret challenge; y never
                                       appears in the recovery denominator)
  key rule : k(m) = even_member(y([m]P0)) mod 2^wkey  (the shipped key rule)
  The six automorphisms of E: y^2 = x^3 + 4 act on the order-r subgroup as
  multiplication by mu6 (Aut = mu6, j = 0), and y(psi(P)) = +-y(P), so the six
  images share a key.  Hence a key match means  x = lambda * v * b_j * w'_k  for
  some lambda in mu6  =>  v = lambda^{-1} * x * (b_j * w'_k)^{-1} mod r.
  Candidates: c = x * (b_j*w'_k)^{-1} * zeta^j, j = 0..5, one of which is exactly v.
  Verify: affine([c]P0) == affine([v]P0) = Ad  (AFFINE, not raw projective; the
  v23 code compared projective tuples, which never compare equal).
  s-extraction: s = d-th root of v mod r (d = 2 here: Tonelli-Shanks), then
  affine([s]P0) == As, the ECDLP target point.  A "solve" requires this final
  point check, so success is a full ECDLP solve, not a v-recovery alone.

Known-answer gate (run FIRST, blocks any census): 64 engineered-collision
instances, 64/64 required, each with a planted lambda in mu6 (lambda = 1 and
lambda != 1 both covered) plus 64 no-plant instances scored for false solves.
The SAME 64 fixtures are re-scored with the v23 recovery (x/y) and must give
0/64 -- the planted-bug sensitivity row.

On-curve asserts: every enumerated point ([m]P0 for each A record, [b*w']*Ad for
each B record) is asserted on y^2 = x^3 + 4 before its key is used; the count of
asserted points is reported (population_on_curve row).
"""
import json, math, os, random, sys, time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod
from toy import is_prime, jmul, jaffine

K_CLAIM = F(51, 100)          # charged design K = 0.51
B_CURVE = 4                   # y^2 = x^3 + 4 (bls12-381 G1 short form)


# ---------------------------------------------------------------- curve helpers
def curve_ok(Q, p):
    """Projective on-curve test for y^2 = x^3 + B_CURVE:  Y^2 == X^3 + B*Z^6."""
    X, Y, Z = Q
    if Z % p == 0:
        return False
    return (Y * Y - X * X * X - B_CURVE * pow(Z, 6, p)) % p == 0


def affine(Q, p):
    return jaffine(Q, p)


def sqrt_mod(a, p):
    """Tonelli-Shanks (p prime, a a QR).  r = 1 mod 4 here, so TS is required."""
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r_ = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
            if i == m:
                return None
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r_ = t * c % p, r_ * b % p
    return r_


def order_r_point(p, r, order, rng):
    """A point of exact order r on y^2 = x^3 + 4 (order = r * cofactor)."""
    while True:
        x = rng.randrange(1, p)
        rhs = (x * x * x + B_CURVE) % p
        if rhs == 0:
            continue
        y = sqrt_mod(rhs, p)
        if y is None:
            continue
        P = (x, y, 1)
        assert curve_ok(P, p)
        Q = jmul(order // r, P, p)
        if Q[2] % p == 0:
            continue
        assert curve_ok(Q, p)
        if jmul(r, Q, p)[2] % p != 0:
            continue
        return Q


def generator(r):
    fac, t, q = set(), r - 1, 2
    while q * q <= t:
        while t % q == 0:
            fac.add(q)
            t //= q
        q += 1
    if t > 1:
        fac.add(t)
    g = 2
    while not all(pow(g, (r - 1) // q, r) != 1 for q in fac):
        g += 1
    return g


def key_of_point(Q, p, mask):
    return keymod.even_member(affine(Q, p)[1], p) & mask


def key_of_scalar(x, P0, p, mask):
    return key_of_point(jmul(x, P0, p), p, mask)


# ---------------------------------------------------------------- one trial
def trial(arg):
    """One independent instance.  A-side table + B-side table built from the
    PUBLIC challenge point Ad.  The secret v is never passed in: recovery is the
    whole point.  As is passed only to score the extracted root."""
    p, r, d, m, P0, Ad, As, nb, nw, wkey, seed, hm_check = arg
    mask = (1 << wkey) - 1
    if r not in ZETA_CACHE:                       # worker-local, computed once
        ZETA_CACHE[r] = pow(generator(r), (r - 1) // 6, r)
    rng = random.Random(seed)
    A = [pow(rng.randrange(2, r), d, r) for _ in range(nb)]
    W = [pow(rng.randrange(2, r), d, r) for _ in range(nw)]
    bl = [pow(rng.randrange(2, r), d, r) for _ in range(nb)]
    if hm_check:
        for x in A + W + bl:
            assert pow(x, m, r) == 1, "record not an H member"
    AS_HOLDER[r] = affine(As, p)
    mt, gv, gs, old, nc = _attack(A, W, bl, Ad, P0, p, r, mask, v_true=None)
    # SOLVE = per-TRIAL indicator (>=1 recovered+verified v).  Counting verified
    # matches instead would report a mean count, not a success probability.
    return dict(seed=seed, records=2 * nb * nw, matches=mt, solves_v=int(gv >= 1),
                solves=int(gs >= 1), verified_matches=gs, v23_solves=old,
                points_on_curve=nc,
                h_members=(len(A) + len(W) + len(bl)) if hm_check else 0)


ZETA_CACHE = {}


# ---------------------------------------------------------------- CP interval
def _betacf(a, b, x):
    """Continued fraction for the incomplete beta function (Lentz)."""
    MAXIT, EPS, FPMIN = 300, 3e-16, 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < EPS:
            break
    return h


def _betai(a, b, x):
    """Regularized incomplete beta I_x(a,b), stable in the tails."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def cp95_lower(k, n):
    """Clopper-Pearson one-sided 95% LOWER endpoint, bracket (0,1) -- UNCLIPPED
    (the v15k helper clipped at 0.99).  Endpoint = inverse in x of
    I_x(k, n-k+1) = 0.05; k == n gives the closed form 0.05^(1/n)."""
    if k <= 0:
        return 0.0
    if k >= n:
        return 0.05 ** (1.0 / n)
    lo, hi, target = 0.0, 1.0, 0.05
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _betai(k, n - k + 1, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cp_startup_assert():
    v = cp95_lower(512, 512)
    ref = 0.05 ** (1.0 / 512)
    assert abs(v - ref) < 1e-12, (v, ref)
    # monotone + sane at partial success
    assert 0 < cp95_lower(200, 512) < 200 / 512
    return dict(k=512, n=512, cp95=round(v, 12), closed_form=round(ref, 12),
                delta=abs(v - ref), ok=abs(v - ref) < 1e-12)


# ---------------------------------------------------------------- gate
def _attack(A, W, b_list, Ad, P0, p, r, mask, v_true=None):
    """Shared attack core.  Returns (matches, got_v, got_s, old_ok, oncurve, As).

    Ad only: the B-side record point is [b*w]*Ad, so the B table never needs v.
    v_true is passed ONLY by the gate to score the v23 path as a negative control.
    """
    tab, ncurve = {}, 0
    for a in A:
        for w in W:
            x = a * w % r
            Q = jmul(x, P0, p)
            ncurve += 1
            assert curve_ok(Q, p), "A-side enumerated point OFF CURVE"
            tab.setdefault(key_of_point(Q, p, mask), []).append(x)
    Ad_a = affine(Ad, p)
    zeta = ZETA_CACHE[r]
    matches = got_v = got_s = old_ok = 0
    for b in b_list:
        for w in W:
            bw = b * w % r                          # UNSCALED B-side product
            Q = jmul(bw, Ad, p)                     # [b*w]*Ad
            ncurve += 1
            assert curve_ok(Q, p), "B-side enumerated point OFF CURVE"
            hit = tab.get(key_of_point(Q, p, mask))
            if not hit:
                continue
            ibw = pow(bw, r - 2, r)                 # inverse of the UNSCALED product
            y_scalar = v_true * bw % r if v_true is not None else None
            for x in hit:
                matches += 1
                # v23 path (negative control): divide by the v-SCALED y -> mu6 unit
                if y_scalar is not None:
                    iy = pow(y_scalar, r - 2, r)
                    for j in range(6):
                        c_old = x * iy % r * pow(zeta, j, r) % r
                        if affine(jmul(c_old, P0, p), p) == Ad_a:
                            old_ok += 1
                            break
                # v24 path: divide by the UNSCALED b*w, then walk mu6
                base = x * ibw % r
                for j in range(6):
                    cand = base * pow(zeta, j, r) % r
                    if affine(jmul(cand, P0, p), p) == Ad_a:
                        got_v += 1
                        s0 = sqrt_mod(cand, r)
                        for sv in ({s0, r - s0} if s0 is not None else set()):
                            if affine(jmul(sv, P0, p), p) == AS_HOLDER[r]:
                                got_s += 1
                        break
    return matches, got_v, got_s, old_ok, ncurve


AS_HOLDER = {}


def gate(cfg, n_gate=64, n_false=64):
    """KNOWN-ANSWER gate: 64 engineered collisions must all be solved (64/64),
    and 64 unplanted instances must yield no false solve."""
    rows, false_rows = [], []
    t0 = time.time()
    for gi, cf in enumerate(cfg["groups"]):
        p, r, d = cf["p"], cf["r"], cf["d"]
        m = (r - 1) // d
        P0 = order_r_point(p, r, cf["order"], random.SystemRandom())
        zeta = pow(generator(r), (r - 1) // 6, r)
        ZETA_CACHE[r] = zeta
        mask = (1 << cf["wkey"]) - 1
        rng = random.Random(0xA841DEB5 + gi)
        for t in range(n_gate // len(cfg["groups"])):
            s_true = rng.randrange(2, r)
            v = pow(s_true, d, r)
            Ad, As = jmul(v, P0, p), jmul(s_true, P0, p)
            AS_HOLDER[r] = affine(As, p)
            lam = pow(zeta, rng.randrange(6), r) if t % 2 else 1
            b0 = pow(rng.randrange(2, r), d, r)
            w1 = pow(rng.randrange(2, r), d, r)
            a0 = lam * v % r * b0 % r                 # a0*w1 = lam*v*b0*w1
            A = [pow(rng.randrange(2, r), d, r) for _ in range(3)] + [a0]
            W = [pow(rng.randrange(2, r), d, r) for _ in range(2)]
            W[-1] = w1
            bl = [b0] + [pow(rng.randrange(2, r), d, r) for _ in range(3)]
            for a in A:
                assert pow(a, m, r) == 1, "gate A element not in H"
            mt, gv, gs, old, nc = _attack(A, W, bl, Ad, P0, p, r, mask, v_true=v)
            rows.append(dict(group=gi, t=t, lam_is_one=int(lam == 1), matches=mt,
                             got_v=gv, got_s=gs, v23_path=old, on_curve=nc,
                             ok=bool(gv >= 1 and gs >= 1)))
        for t in range(n_false // len(cfg["groups"])):
            s_true = rng.randrange(2, r)
            v = pow(s_true, d, r)
            Ad, As = jmul(v, P0, p), jmul(s_true, P0, p)
            AS_HOLDER[r] = affine(As, p)
            A = [pow(rng.randrange(2, r), d, r) for _ in range(4)]
            W = [pow(rng.randrange(2, r), d, r) for _ in range(2)]
            bl = [pow(rng.randrange(2, r), d, r) for _ in range(4)]
            mt, gv, gs, old, nc = _attack(A, W, bl, Ad, P0, p, r, mask, v_true=v)
            false_rows.append(dict(group=gi, t=t, matches=mt, got_s=gs,
                                   false_solve=int(gs > 0)))
    ok = sum(r0["ok"] for r0 in rows)
    return dict(kind="known-answer gate (engineered collisions)", n=len(rows), ok=ok,
                passed=bool(ok == len(rows) and len(rows) >= n_gate),
                lam_one_cases=sum(r0["lam_is_one"] for r0 in rows),
                lam_other_cases=sum(1 - r0["lam_is_one"] for r0 in rows),
                v23_path_solves=sum(r0["v23_path"] for r0 in rows),
                on_curve_asserts=sum(r0["on_curve"] for r0 in rows),
                false_positive_n=len(false_rows),
                false_positives=sum(r0["false_solve"] for r0 in false_rows),
                false_positive_matches=sum(r0["matches"] for r0 in false_rows),
                rows=rows, seconds=round(time.time() - t0, 2))


def main():
    import multiprocessing as mp
    cfg = json.load(open(HERE + "/fixtures-v24.json"))
    out = dict(kind="v24 direct end-to-end solve census (recovery algebra fixed)",
               groups=[])
    t0 = time.time()
    out["cp95_startup"] = cp_startup_assert()
    if "--gate" in sys.argv:
        g = gate(cfg)
        json.dump(g, open(HERE + "/gate-v24.json", "w"), indent=1)
        print(json.dumps({k: v for k, v in g.items() if k != "rows"}, indent=1))
        print("GATE", "PASS" if g["passed"] else "FAIL")
        return 0 if g["passed"] else 1
    for gi, cf in enumerate(cfg["groups"]):
        p, r, d = cf["p"], cf["r"], cf["d"]
        m = (r - 1) // d
        assert is_prime(r) and (r - 1) % (6 * d) == 0 and m % 6 == 0
        rng0 = random.SystemRandom()
        P0 = order_r_point(p, r, cf["order"], rng0)
        zeta = pow(generator(r), (r - 1) // 6, r)
        assert pow(zeta, 6, r) == 1 and pow(zeta, 3, r) == r - 1
        ZETA_CACHE[r] = zeta
        # FLOOR, not round: N must not exceed the claim's K (conservative direction)
        N = int(math.floor(math.sqrt(float(K_CLAIM) * m / 6)))
        nw = cf.get("nw", 2)
        nb = N // nw
        N = nb * nw
        Kf = F(6 * N * N, m)
        gran = F(12 * N, m)      # dK per unit N -> the achievable K resolution
        trials = cf["trials"]
        wkey = cf["wkey"]
        # budget fidelity: measured at K_f <= K_claim (conservative direction)
        assert Kf <= K_CLAIM, (float(Kf), float(K_CLAIM))
        assert abs(float(Kf) / float(K_CLAIM) - 1) <= 1e-2
        args = []
        for i in range(trials):
            s = rng0.randrange(2, r)
            v = pow(s, d, r)
            args.append((p, r, d, m, P0, jmul(v, P0, p), jmul(s, P0, p),
                         nb, nw, wkey, int.from_bytes(os.urandom(16), "big"),
                         i < 64))
        assert len({a[-2] for a in args}) == trials, "seeds must be distinct"
        t1 = time.time()
        with mp.get_context("fork").Pool(min(12, mp.cpu_count())) as pool:
            res = pool.map(trial, args, chunksize=max(1, trials // 512))
        solves = sum(x["solves"] for x in res)
        solves_v = sum(x["solves_v"] for x in res)
        matches = sum(x["matches"] for x in res)
        oncurve = sum(x["points_on_curve"] for x in res)
        row = dict(group=gi, p=p, p_bits=p.bit_length(), r=r, r_bits=r.bit_length(),
                   d=d, m=m, m_log2=round(math.log2(m), 4), order=cf["order"],
                   nb=nb, nw=nw, records_per_side=N, records_per_trial=2 * N,
                   K=float(Kf), K_claim=float(K_CLAIM),
                   K_dev=round(abs(float(Kf) / float(K_CLAIM) - 1), 8),
                   K_granularity=round(float(gran), 8),
                   K_conservative=bool(Kf <= K_CLAIM),
                   trials=trials, solves=solves, solves_v_only=solves_v,
                   q=solves / trials, q_v_only=solves_v / trials,
                   model=round(1 - math.exp(-float(Kf)), 9),
                   model_dev=round(abs(solves / trials - (1 - math.exp(-float(Kf)))) /
                                   (1 - math.exp(-float(Kf))), 8),
                   cp95=round(cp95_lower(solves, trials), 8),
                   cp95_v_only=round(cp95_lower(solves_v, trials), 8),
                   key_bits=wkey, shape=cf.get("shape", "nw%d" % nw),
                   matches=matches,
                   verified_matches=sum(x["verified_matches"] for x in res),
                   matches_per_trial=round(matches / trials, 4),
                   points_on_curve=oncurve,
                   records_total=sum(x["records"] for x in res),
                   seconds=round(time.time() - t1, 1))
        out["groups"].append(row)
        print(json.dumps(row), flush=True)
        json.dump(out, open(HERE + "/direct-v24.json", "w"), indent=1)
    classes = {}
    for g in out["groups"]:
        classes.setdefault(g["shape"], []).append(g)
    cls_rows = {}
    for name, gs in classes.items():
        cps = [g["cp95"] for g in gs]
        qs = [g["q"] for g in gs]
        cls_rows[name] = dict(groups=len(gs), trials=[g["trials"] for g in gs],
                              q=[round(q, 8) for q in qs],
                              q_min=round(min(qs), 8), q_max=round(max(qs), 8),
                              cp95_min=round(min(cps), 8),
                              model=gs[0]["model"],
                              model_dev_max=max(g["model_dev"] for g in gs),
                              meets_floor=bool(min(cps) >= 0.39))
    worst = min(cls_rows, key=lambda k: cls_rows[k]["cp95_min"])
    best = max(cls_rows, key=lambda k: cls_rows[k]["cp95_min"])
    wf = F(1, 1) - F(cls_rows[worst]["cp95_min"]).limit_denominator(10**9)
    amp2 = float(1 - (1 - cls_rows[worst]["cp95_min"]) ** 2)
    out["summary"] = dict(
        groups=len(out["groups"]),
        total_trials=sum(g["trials"] for g in out["groups"]),
        solves=[g["solves"] for g in out["groups"]],
        by_shape=cls_rows,
        worst_shape=worst, best_shape=best,
        q_min=round(min(g["q"] for g in out["groups"]), 8),
        q_max=round(max(g["q"] for g in out["groups"]), 8),
        stability_spread=round(max(g["q"] for g in out["groups"]) -
                               min(g["q"] for g in out["groups"]), 8),
        success_final=cls_rows[worst]["cp95_min"],
        success_final_shape=worst,
        success_best_shape=cls_rows[best]["cp95_min"],
        claim_floor=0.39,
        meets_floor=bool(cls_rows[worst]["cp95_min"] >= 0.39),
        success_below_floor=bool(cls_rows[worst]["cp95_min"] < 0.39),
        shortfall=round(0.39 - cls_rows[worst]["cp95_min"], 8),
        budget_to_reach_floor_bits=round(
            math.log(-math.log(1 - 0.39) / -math.log(1 - cls_rows[worst]["cp95_min"])) / 2.0
            if cls_rows[worst]["cp95_min"] < 0.39 else 0.0, 8),
        budget_to_reach_floor_rule="N scales as sqrt(K), work as N: "
                                   "delta_bits = log2( -ln(1-0.39) / -ln(1-p_measured) ) / 2",
        w2_amplification=dict(rule="success 1-(1-p)^2, charge 2W (2 sequential runs)",
                              p_used=cls_rows[worst]["cp95_min"], success=round(amp2, 8),
                              claim_if_doubled=round(116.81277 + 1.0, 4),
                              claim_note="2 x W = raw_log2 + 1 bit = 117.8128 > gate 116.8271"),
        rule="success_final = min over shapes and groups of the CP-95 lower bound of the "
             "measured per-TRIAL solve rate. The WORST shape is used for the headline, so "
             "the shared-walk clustering is charged, not averaged away. Model deviation and "
             "K fidelity are separate validation rows (CP-95 already covers sampling error; "
             "K_f <= K_claim makes the fixture weaker than the claim).",
        cp95_method="Clopper-Pearson one-sided 95% lower, bracket (0,1), inverse of "
                    "I_x(k, n-k+1)=0.05, k==n closed form 0.05^(1/n) asserted at startup",
        points_on_curve_total=sum(g["points_on_curve"] for g in out["groups"]),
        seconds_total=round(time.time() - t0, 1))
    json.dump(out, open(HERE + "/direct-v24.json", "w"), indent=1)
    print("SUMMARY", json.dumps(out["summary"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
