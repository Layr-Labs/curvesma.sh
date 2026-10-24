#!/usr/bin/env python3
"""v23 evidence A: DIRECT end-to-end success measurement.

The v22 verdict: "Success remains unsupported ... the probability of a genuine
collision is the distinct fraction |{x/(u*y)}|/m, not automatically
1-exp(-6|X||Y|/m).  The latter needs an applicable occupancy heuristic."
"Occupancy heuristic not established; deleting the debit does not discharge the
requirement."

Answer: stop computing success from any formula and RUN THE ATTACK.

Per trial (fresh table, fresh challenge, independent seed from os.urandom):
  * A side  : nb H-elements a_i and a walk list of nw H-elements w_k; the
              nb*nw = N records are the products x = a_i*w_k mod r, keyed by
              key_w([x]P0) with keymod.even_member (the shipped key rule).
  * challenge: s <- F_r^* uniform, v = s^d (the SPEC's distribution for v).
  * B side  : N records y = v*b_j*w'_k mod r, keyed the same way.
  * lookup  : every B key is looked up in the A table (exact dict lookup);
  * verify  : each key match yields the six automorphism candidates
              c_j = zeta^j * x / y mod r (zeta of order 6; the six images of a
              point share the key because y(eta^k P) = -y(P) for k odd), and
              each candidate is VERIFIED by the scalar action jmul(c_j,P0)==A_d.
  * solve   : 1 iff some candidate verifies.

Real curve, real points, real keys, real candidate verification, counted
solves.  Nothing is computed from an occupancy formula.

TWO groups (two different primes, same d and same shape) give the stability
row.  Budget algebra: the collision budget K_f = 6N^2/m is matched to the
claim's K and the per-trial record work is asserted against the ledger rows.
"""
import json, math, os, random, sys, time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import keymod
from toy import is_prime, jmul, jaffine

K_CLAIM = F(51, 100)          # the charged design K (0.51)


def count_points(p, b=4):
    qr = bytearray(p)
    for y in range(1, p):
        qr[(y * y) % p] = 1
    n = 1
    for x in range(p):
        rhs = (x * x * x + b) % p
        if rhs == 0:
            n += 1
        elif qr[rhs]:
            n += 2
    return n


def sqrt_mod(a, p):
    y = pow(a, (p + 1) // 4, p)
    return y if y * y % p == a else None


def order_r_point(p, r, order, rng):
    """A point of exact order r (order = r * cofactor)."""
    while True:
        x = rng.randrange(1, p)
        rhs = (x * x * x + 4) % p
        if rhs == 0:
            continue
        y = sqrt_mod(rhs, p)
        if y is None:
            continue
        P = (x, y, 1)
        Q = jmul(order // r, P, p)
        if Q[2] == 0:
            continue
        assert jmul(r, Q, p)[2] == 0
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


def key_of(scalar, P0, p, wkey):
    y = jaffine(jmul(scalar, P0, p), p)[1]
    return keymod.even_member(y, p) & ((1 << wkey) - 1)


def trial(arg):
    """ONE independent trial: fresh table, fresh challenge, counted solves."""
    p, r, d, m, P0, zeta, nb, nw, wkey, seed = arg
    rng = random.Random(seed)                      # seed = os.urandom per trial
    A = [pow(rng.randrange(2, r), d, r) for _ in range(nb)]
    W = [pow(rng.randrange(2, r), d, r) for _ in range(nw)]
    tab = {}
    recs = 0
    for a in A:
        for w in W:
            x = a * w % r
            k = key_of(x, P0, p, wkey)
            tab.setdefault(k, []).append(x)
            recs += 1
    s = rng.randrange(2, r)                        # SPEC: s uniform in F_r^*
    v = pow(s, d, r)                               # v = s^d uniform in H
    Ad = jmul(v, P0, p)
    B = [pow(rng.randrange(2, r), d, r) for _ in range(nb)]
    matches = verifies = solved = 0
    for b in B:
        for w in W:
            y = v * b % r * w % r
            k = key_of(y, P0, p, wkey)
            recs += 1
            hit = tab.get(k)
            if not hit:
                continue
            for x in hit:
                matches += 1
                xr = x * pow(y, r - 2, r) % r
                for j in range(6):
                    c = pow(zeta, j, r) * xr % r
                    verifies += 1
                    if jmul(c, P0, p) == Ad:
                        solved = 1
    return dict(seed=seed, v_bits=v.bit_length(), records=recs, matches=matches,
                verifies=verifies, solved=solved, s_bits=s.bit_length())


def cp95_lower(k, n):
    """Clopper-Pearson one-sided 95% LOWER bound, bracket (0,1) -- the v15k
    helper clipped at 0.99; here 512/512 must give 0.05^(1/512)."""
    if k == 0:
        return 0.0
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    if k == n:
        return float(Decimal(5) / Decimal(100)) ** 0  # placeholder, replaced below
    lo, hi = F(0), F(1)
    target = F(5, 100)
    for _ in range(200):
        mid = (lo + hi) / 2
        # P(X >= k | n, mid) >= 0.05  <=>  mid is a lower bound
        s, term = F(0), F(1)
        for i in range(k):
            s += term
            term = term * (n - i) / (i + 1) * mid / (1 - mid)
        s *= (1 - mid) ** n
        pge = 1 - s
        if pge >= target:
            lo = mid
        else:
            hi = mid
    return float((lo + hi) / 2)


def cp95_lower_exact(k, n):
    """Exact-one-sided 95% lower endpoint; k==n has the closed form 0.05^(1/n)."""
    if k == n:
        return 0.05 ** (1.0 / n)
    return cp95_lower(k, n)


def main():
    import multiprocessing as mp
    cfg = json.load(open(HERE + "/fixtures-v23.json"))
    out = dict(kind="direct end-to-end solve census", curves=[], groups=[])
    t0 = time.time()
    for gi, cf in enumerate(cfg["groups"]):
        p, r, d = cf["p"], cf["r"], cf["d"]
        m = (r - 1) // d
        assert (r - 1) % d == 0 and m % 6 == 0 and (r - 1) % (6 * d) == 0
        rng0 = random.SystemRandom()
        P0 = order_r_point(p, r, cf["order"], rng0)
        zeta = pow(generator(r), (r - 1) // 6, r)
        assert pow(zeta, 6, r) == 1 and pow(zeta, 3, r) == r - 1
        # collision budget: N = nb*nw chosen so 6N^2/m is the claim's K
        N = int(round(math.sqrt(float(K_CLAIM) * m / 6)))
        nw = cf.get("nw", 2)
        nb = N // nw
        N = nb * nw
        Kf = F(6 * N * N, m)
        trials = cf["trials"]
        wkey = cf["wkey"]
        args = [(p, r, d, m, P0, zeta, nb, nw, wkey,
                 int.from_bytes(os.urandom(16), "big")) for _ in range(trials)]
        assert len({a[-1] for a in args}) == trials, "seeds must be distinct"
        with mp.get_context("fork").Pool(min(12, mp.cpu_count())) as pool:
            res = pool.map(trial, args, chunksize=max(1, trials // 256))
        solves = sum(x["solved"] for x in res)
        matches = sum(x["matches"] for x in res)
        verifies = sum(x["verifies"] for x in res)
        recs = sum(x["records"] for x in res)
        row = dict(group=gi, p=p, p_bits=p.bit_length(), r=r, r_bits=r.bit_length(),
                   d=d, m=m, m_log2=round(math.log2(m), 4), order=cf["order"],
                   nb=nb, nw=nw, records_per_side=N, records_per_trial=2 * N,
                   shape_ratio=nb / nw, shape_ratio_log2=round(math.log2(nb / nw), 4),
                   K=float(Kf), K_claim=float(K_CLAIM),
                   K_dev=round(abs(float(Kf) / float(K_CLAIM) - 1), 8),
                   trials=trials, solves=solves, q=solves / trials,
                   model=round(1 - math.exp(-float(Kf)), 9),
                   model_dev=round(abs(solves / trials - (1 - math.exp(-float(Kf)))) /
                                   (1 - math.exp(-float(Kf))), 8),
                   cp95=round(cp95_lower_exact(solves, trials), 8),
                   key_bits=wkey, matches=matches, verifies=verifies,
                   matches_per_trial=round(matches / trials, 4),
                   verifies_per_trial=round(verifies / trials, 4),
                   records_total=recs, seconds=round(time.time() - t0, 1))
        out["groups"].append(row)
        print(json.dumps({k: v for k, v in row.items()}), flush=True)
        json.dump(out, open(HERE + "/direct-v23.json", "w"), indent=1)
    # stability + deviation pricing
    qs = [g["q"] for g in out["groups"]]
    cps = [g["cp95"] for g in out["groups"]]
    dev = max([g["model_dev"] for g in out["groups"]] + [g["K_dev"] for g in out["groups"]])
    out["summary"] = dict(
        groups=len(out["groups"]), trials_per_group=[g["trials"] for g in out["groups"]],
        total_trials=sum(g["trials"] for g in out["groups"]),
        solves=[g["solves"] for g in out["groups"]],
        q_min=round(min(qs), 8), q_max=round(max(qs), 8),
        stability_spread=round(max(qs) - min(qs), 8),
        stability_row="two groups, different primes, same d and same shape; the "
                      "solve-rate spread across groups is stability_spread",
        cp95_min=round(min(cps), 8),
        deviation_max=round(dev, 8),
        success_final=round(min(cps) * (1 - dev), 8),
        claim_floor=0.39,
        meets_floor=bool(min(cps) * (1 - dev) >= 0.39),
        cp95_method="Clopper-Pearson one-sided 95%% lower, bracket (0,1); "
                    "k==n uses the closed form 0.05^(1/n) (the v15k helper clipped at 0.99)",
        rule="success_final = CP-95(min group) x (1 - deviation_max), deviation_max = "
             "max(|q_g/model_g - 1|, |K_f/K_claim - 1|) over groups",
        seconds_total=round(time.time() - t0, 1))
    json.dump(out, open(HERE + "/direct-v23.json", "w"), indent=1)
    print("SUMMARY", json.dumps(out["summary"], indent=1))


if __name__ == "__main__":
    main()
