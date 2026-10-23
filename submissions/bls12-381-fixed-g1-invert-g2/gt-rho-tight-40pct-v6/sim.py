#!/usr/bin/env python3
"""Exponent-domain simulation of the same walk, with Frobenius-class quotienting.

The state is a pair (a, b) of exponents, i.e. the element g^a z^b of the order-r
subgroup.  Nothing here is a stand-in for a different walk:

  * the class of (a, b) is its twelve-element Frobenius orbit, and the class key is the
    minimum over that orbit of an independently mixed word -- the same construction as
    production, where the key is the minimum over the twelve conjugate w-coefficient
    keys.  Two states have the same key exactly when they are in the same class (up to
    the mixer), so the distinguished test, the table index and the step type are class
    functions, as in production;
  * the step is the production step with the v6 equivariant rule: the held state's
    Frobenius offset k (the argmin position, the exponent-domain twin of tower.kappa)
    selects the table twist (routing_twist + k) mod 12, so the step is a class function.
    `legacy=True` runs the v5 rule (a per-trail frozen twist j = mix(key, seg) % 12) as a
    CONTROL on the same field, table and seeds, so the cost of the defect is measured;
  * the step type is gated by f: st < f*8 selects a squaring.  f is a parameter of the
    transition rule, not just of the budget, so the f = 0 control walks a genuinely
    different map (no squaring branch at all) on the same field, the same table and the
    same trail seeds;
  * the detector is detector.Detector -- the same object toy_core.py executes, imported;
  * recovery is the production recovery: align the two endpoints by the Frobenius power
    that maps one to the other, solve x = (a2 p^j - a1)/(b1 - b2 p^j) mod r, and verify.
"""
import json, random, sys, time
from math import sqrt, exp
import constants as K
import routing
import classwalk
import detector as D


def find_field(bits, seed):
    """Prime r near 2^bits with a prime p of order 12 modulo r (ord_r(p) = 12)."""
    rng = random.Random(seed)
    while True:
        r = rng.randrange(1 << bits, 1 << (bits+1)) | 1
        if not _is_prime(r) or (r-1) % 12:
            continue
        lam = None
        for _ in range(200):
            c = rng.randrange(2, r-1)
            cand = pow(c, (r-1)//12, r)
            if cand != 1 and all(pow(cand, k, r) != 1 for k in (1, 2, 3, 4, 6)):
                lam = cand
                break
        if lam is None:
            continue
        k0 = rng.randrange(1, 200)
        for kk in range(k0, k0+4000):
            p = lam + kk*r
            if _is_prime(p):
                return p, r, (r-1)//12


def _is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d = n-1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n-1):
            continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1:
                break
        else:
            return False
    return True


def make_mix(seed):
    rng = random.Random(seed)
    M1 = rng.randrange(1 << 40, 1 << 61) | 1
    M2 = rng.randrange(1 << 40, 1 << 61) | 1
    MASK = (1 << 61) - 1
    def mix(x, y):
        h = (x*M1 + y*M2) & MASK
        h ^= h >> 29
        h = (h * M1) & MASK
        h ^= h >> 31
        return h
    return mix


class Field:
    def __init__(self, p, r, N, seed):
        self.p = p; self.r = r; self.N = N
        self.pw = [pow(p, j, r) for j in range(12)]      # Frobenius on exponents
        self.mix = make_mix(seed)

    def class_key_argmin(self, e):
        """(key, k) for the class of g^e: the minimum over the twelve images and the
        position k of that minimum -- the exponent-domain twin of tower.kappa()."""
        r = self.r; pw = self.pw; mix = self.mix
        m = None
        k = 0
        for j in range(12):
            v = mix(e*pw[j] % r, 12)      # constant 2nd arg: the orbit set is j-invariant
            if m is None or v < m:
                m = v; k = j
        return m, k

    def class_key(self, e):
        return self.class_key_argmin(e)[0]

    def table(self, R, rng):
        """R multipliers with their twelve twists: exponent pairs (u p^i, v p^i)."""
        r = self.r
        tab = []
        for _ in range(R):
            u = rng.randrange(1, r); v = rng.randrange(1, r)
            tab.append([(u*self.pw[i] % r, v*self.pw[i] % r) for i in range(12)])
        return tab


def run(field, R, groups, seed, f, dp_bits, legacy=False, tab=None, mu=K.MU_WALK):
    """One group = the production trail structure at scale: S/theta trails, each ending
    at its first report (or at the metered cutoff)."""
    r = field.r; mix = field.mix
    S = K.step_budget(mu=mu, target=K.TARGET, f=f, n=field.N)
    theta = 1 << dp_bits
    w_trail = K.TRAIL_METERS * theta                     # 16 DP lengths per trail
    n_trails = max(2, int(round(S/theta)))               # realised steps ~ S
    thresh = int(round(8*float(f)))
    out = dict(groups=groups, S=int(S), f=float(f), dp_bits=dp_bits, n_trails=n_trails,
               w_trail=w_trail, thresh=thresh, legacy=bool(legacy), recovered=0,
               recovered_events=0, hits=0, steps=0, escapes=0, aborts=0, sqr_steps=0,
               reports=0)
    for gi in range(groups):
        rng = random.Random(seed*1000003 + gi)
        secret = rng.randrange(1, r)
        table = field.table(R, rng) if tab is None else tab
        store = {}
        det = None
        sqr = 0
        got = 0
        for seg in range(n_trails):
            a = rng.randrange(1, r); b = rng.randrange(1, r)
            e = (a + b*secret) % r
            key, held0 = field.class_key_argmin(e)
            j = mix(key, seg) % 12 if legacy else 0
            det = D.Detector()
            det.start(0, key, mix(key, 1))
            t = 0
            got_this = False
            while t < w_trail:
                t += 1
                out['steps'] += 1
                key, held = field.class_key_argmin(e)
                tbl, tw, dp, st = routing.fields(key)
                if st < thresh:
                    a = 2*a % r; b = 2*b % r; e = 2*e % r; sqr += 1
                else:
                    idx = (tw + j) % 12 if legacy else classwalk.twist_index(tw, held)
                    u, v = table[tbl % R][idx]
                    a = (a + u) % r; b = (b + v) % r; e = (e + u + v*secret) % r
                key, held = field.class_key_argmin(e)
                ev = det.observe(t, key)
                if ev == 'abort':
                    out['aborts'] += 1
                    break
                if ev is not None:
                    # class-only escape: double the multiplier indexed by this class key
                    tbl, tw, dp, st = routing.fields(det.escape_key)
                    idx = (tw + j) % 12 if legacy else classwalk.twist_index(tw, held)
                    u, v = table[tbl % R][idx]
                    a = (a + u + u) % r; b = (b + v + v) % r
                    e = (e + u + u + (v + v)*secret) % r
                    key, held = field.class_key_argmin(e)
                    det.start(t, key, mix(key, t))
                if (routing.fields(key)[2] & (theta - 1)) == 0:
                    out['reports'] += 1
                    if key in store:
                        other, other_e = store[key]
                        jj = None
                        for cand in range(12):      # align the ELEMENTS, not the pairs
                            if (other_e*field.pw[cand] - e) % r == 0:
                                jj = cand
                                break
                        if jj is not None:
                            out['hits'] += 1
                            pj = field.pw[jj]
                            den = (other[1]*pj - b) % r
                            if den:
                                x = (a - other[0]*pj) * pow(den, r-2, r) % r
                                if x == secret:
                                    out['recovered_events'] += 1
                                    if not got_this:
                                        got += 1
                                        got_this = True
                    else:
                        store[key] = ((a, b), e)
                    break        # production trail ends at its first report
        if det is not None:
            out['escapes'] += det.escapes
        out['sqr_steps'] += sqr
        if got:
            out['recovered'] += 1
    return out


def pooled(rows):
    w = sum(x['groups'] for x in rows)
    p = sum(x['recovered'] for x in rows)/w
    return p, sqrt(max(p*(1-p), 1e-12)/w)


def pooled_model(fields, rows, f, mu=K.MU_WALK):
    w = sum(x['groups'] for x in rows)
    return sum((1 - exp(-x['S']**2*(1-float(f)**2)/(2*float(mu)*fields[i].N)))*x['groups']
               for i, x in enumerate(rows))/w


def equivariance_probe(field, samples=4000, seed=9091):
    """Measure the class-function property in the exponent domain.

    For each sample, take e and a Frobenius image e*p^c of it, and check that a step from
    each lands in ONE next class.  Reported for the v6 rule and for the v5 rule (a fixed
    class-derived index) on exactly the same samples and the same multiplier table.
    """
    rng = random.Random(seed)
    r = field.r; pw = field.pw
    R = 64
    t_ok = t_n = l_ok = 0
    skipped = 0
    for _ in range(samples):
        secret = rng.randrange(1, r)
        table = field.table(R, rng)
        a = rng.randrange(1, r); b = rng.randrange(1, r)
        e = (a + b*secret) % r
        c = rng.randrange(12)
        e2 = e*pw[c] % r
        key, held = field.class_key_argmin(e)
        key2, held2 = field.class_key_argmin(e2)
        if key2 != key or held2 != (held - c) % 12:
            skipped += 1            # mixer collision or kappa inconsistency: skip
            continue
        tbl, tw, _, st = routing.fields(key)
        res = {}
        for mode in ('v6', 'v5'):
            out = []
            for (ee, hh) in ((e, held), (e2, held2)):
                if st < K.SQR_THRESHOLD:
                    ne = 2*ee % r
                else:
                    idx = classwalk.twist_index(tw, hh) if mode == 'v6' else (tw + held) % 12
                    uu, vv = table[tbl % R][idx]
                    ne = (ee + uu + vv*secret) % r
                out.append(field.class_key(ne))
            res[mode] = int(out[0] == out[1])
        t_n += 1
        t_ok += res['v6']; l_ok += res['v5']
    return dict(samples=t_n, skipped=skipped, v6_agree=t_ok,
                v6_fraction=t_ok/t_n if t_n else 0.0,
                v5_agree=l_ok, v5_fraction=l_ok/t_n if t_n else 0.0)


def main():
    t0 = time.time()
    fields = []
    for bits, sd in ((20, 11), (21, 12), (22, 13), (23, 14),
                     (20, 21), (21, 22), (22, 23), (23, 24)):
        p, r, N = find_field(bits, sd)
        fields.append(Field(p, r, N, sd*7919))
    G = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    R = 64
    rows = []; ctrl = []; legs = []
    for fd in fields:
        rows.append(run(fd, R, G, seed=fd.p, f=K.FRAC_SQR, dp_bits=5))
        ctrl.append(run(fd, R, G, seed=fd.p, f=K.Fraction(0, 1), dp_bits=5))
        legs.append(run(fd, R, G, seed=fd.p, f=K.FRAC_SQR, dp_bits=5, legacy=True))
    p, se = pooled(rows); pc, sec = pooled(ctrl); pl, sel = pooled(legs)
    tot = sum(x['groups'] for x in rows); totc = sum(x['groups'] for x in ctrl)
    totl = sum(x['groups'] for x in legs)
    model = pooled_model(fields, rows, K.FRAC_SQR)
    model0 = pooled_model(fields, ctrl, K.Fraction(0, 1))
    model_l = pooled_model(fields, legs, K.FRAC_SQR)
    lens = sum(x['steps']/x['groups']/sqrt(fields[i].N)*x['groups'] for i, x in enumerate(rows))/tot
    lens0 = sum(x['steps']/x['groups']/sqrt(fields[i].N)*x['groups'] for i, x in enumerate(ctrl))/totc
    steps_over_S = sum(x['steps']/x['groups']/x['S']*x['groups'] for x in rows)/tot
    steps_over_S0 = sum(x['steps']/x['groups']/x['S']*x['groups'] for x in ctrl)/totc
    out = dict(rows=rows, control=ctrl, legacy=legs,
               pooled_p=p, pooled_se=se, pooled_groups=tot,
               control_p=pc, control_se=sec, control_groups=totc,
               legacy_p=pl, legacy_se=sel, legacy_groups=totl,
               model=model, model_control=model0, model_legacy=model_l,
               sqr_fraction=sum(x['sqr_steps'] for x in rows)/sum(x['steps'] for x in rows),
               control_sqr_fraction=sum(x['sqr_steps'] for x in ctrl)/sum(x['steps'] for x in ctrl),
               injectivity_ratio=lens/lens0 if lens0 else 0.0,
               injectivity_model=1/sqrt(1-float(K.FRAC_SQR)**2),
               ratio_success=p/pl if pl else 0.0,
               ratio_success_model=model/model_l if model_l else 0.0,
               steps_over_S=steps_over_S, steps_over_S_control=steps_over_S0,
               reports_per_group=sum(x['reports'] for x in rows)/tot,
               reports_expected=sum(x['S']/float(1 << 5)*x['groups'] for x in rows)/tot,
               legacy_over_model=pl/model_l if model_l else 0.0,
               equivariance=equivariance_probe(fields[0]),
               f=float(K.FRAC_SQR), mu=float(K.MU_WALK), target=float(K.TARGET),
               escapes=sum(x['escapes'] for x in rows), aborts=sum(x['aborts'] for x in rows),
               legacy_escapes=sum(x['escapes'] for x in legs),
               control_escapes=sum(x['escapes'] for x in ctrl),
               seconds=round(time.time()-t0, 1))
    out['pass'] = bool(p >= float(K.SUCCESS) and abs(p-model) <= 3*se
                       and abs(pc-model0) <= 3*sec
                       and abs(pl-model_l) > 3*sel)
    json.dump(out, open("sim-out.json", "w"), indent=1)
    for i, x in enumerate(rows):
        print("r~2^%2d N=%7d S=%4d trails=%3d P_f=%.4f P_legacy=%.4f sqr=%.4f"
              % (fields[i].r.bit_length(), fields[i].N, x['S'], x['n_trails'],
                 x['recovered']/x['groups'], legs[i]['recovered']/legs[i]['groups'],
                 x['sqr_steps']/x['steps']))
    print("pooled f=3/8   P = %.4f +- %.4f  (%d groups)  model %.4f" % (p, se, tot, model))
    print("pooled f=0     P = %.4f +- %.4f  (%d groups)  model %.4f" % (pc, sec, totc, model0))
    print("pooled v5 rule P = %.4f +- %.4f  (%d groups)  model %.4f"
          % (pl, sel, totl, model_l))
    print("ratio success v6/v5rule %.2f (model %.1f)  injectivity length %.4f (model %.4f)"
          % (out['ratio_success'], out['ratio_success_model'], out['injectivity_ratio'],
             out['injectivity_model']))
    print("squaring fraction f=3/8: %.4f   f=0: %.4f" % (out['sqr_fraction'], out['control_sqr_fraction']))
    print("equivariance: v6 %d/%d  v5 rule %d/%d" % (out['equivariance']['v6_agree'],
          out['equivariance']['samples'], out['equivariance']['v5_agree'],
          out['equivariance']['samples']))
    print("escapes %d aborts %d  steps/S %.4f  reports/group %.1f (expected %.1f)"
          % (out['escapes'], out['aborts'], steps_over_S, out['reports_per_group'],
             out['reports_expected']))
    print("SIM", "PASS" if out['pass'] else "FAIL", "%.1fs" % out['seconds'])
    return 0 if out['pass'] else 1


if __name__ == "__main__":
    sys.exit(main())
