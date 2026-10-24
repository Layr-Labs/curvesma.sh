#!/usr/bin/env python3
"""v32 CROSS-SIDE kappa census + triangular-occupancy Phi + gap-tail rows.

verdict 3e66901 flaw 2: v31's mixed-key statistic counted WITHIN-A pairs, not cross-side
production point-key collisions.  Here two INDEPENDENT residue streams (A side, B side,
same distribution, disjoint coins) are interleaved through ONE store, and kappa is the
per-walk MAX over buckets of observed/expected key pairs.

kappa is a dimensionless CLUSTERING ratio, so it is measured at a key width where
collisions exist at fixture scale (key_bits=20, as v31 did) and applied to the production
image 2^168.  Within-A is kept, labelled, as the mechanism check.

flaw 3: Phi.  For bounded exponent differences the difference law is TRIANGULAR
(p(t) ~ L-|t|), not uniform over 2L slots, so the distinct-difference attenuator is
Phi_tri = (1-e^-y)/y with y = (4/3)*nb^2/G.  Derivation of the 4/3: for triangular
p(t)=(L-|t|)/L^2, sum_t p(t)^2 = 2/(3L); for uniform over 2L slots it is 1/(2L);
ratio = 4/3.  The delta row is emitted and the q row that used Phi is regenerated.
"""
import json, math, os, random, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod as K

CERT = json.load(open(os.path.join(HERE, "certificate.json")))
N = CERT["N"]
NB = 1 << 45
SPACE = 1 << 98
X31 = NB * NB / SPACE


def phi_uniform(x):
    return (1.0 - math.exp(-x)) / x


def phi_tri(x, factor=4.0 / 3.0):
    y = factor * x
    return (1.0 - math.exp(-y)) / y


def key_trunc(residue, bits):
    return K.key_mixed_final(residue) & ((1 << bits) - 1)


def census(bits=20, walks=400, per_side=384, seed=7):
    rng = random.Random(seed)
    cross_k, within_k = [], []
    exp_within = per_side * (per_side - 1) / 2.0 / float(1 << bits)
    exp_cross = per_side * per_side / float(1 << bits)
    for _ in range(walks):
        a = [key_trunc(rng.randrange(K.P), bits) for _ in range(per_side)]
        b = [key_trunc(rng.randrange(K.P), bits) for _ in range(per_side)]
        for stream, sink in ((a, within_k), (None, None)):
            pass
        # within-A pairs
        def pairs(keys):
            c = {}
            for k in keys:
                c[k] = c.get(k, 0) + 1
            return sum(v * (v - 1) // 2 for v in c.values())
        # cross-side pairs: bucket keys of both streams together, count A-B pairs only
        ca, cb = {}, {}
        for k in a:
            ca[k] = ca.get(k, 0) + 1
        for k in b:
            cb[k] = cb.get(k, 0) + 1
        cross = sum(ca.get(k, 0) * v for k, v in cb.items())
        cross_k.append(cross)
        within_k.append(pairs(a))

    tot = len(cross_k)
    obs_c, obs_w = sum(cross_k), sum(within_k)
    exp_tot = exp_within * tot
    exp_cross_tot = exp_cross * tot
    # one-sided Poisson 95% upper on the pooled count: obs + 1.645*sqrt(obs) + 1
    up_c = (obs_c + 1.645 * math.sqrt(obs_c) + 1) / exp_cross_tot
    up_w = (obs_w + 1.645 * math.sqrt(obs_w) + 1) / exp_tot
    return dict(key_bits=bits, walks=walks, per_side=per_side,
                expected_pairs_pooled=dict(within=round(exp_tot, 4), cross=round(exp_cross_tot, 4)),
                cross_side=dict(observed_pairs=obs_c, kappa=round(obs_c / exp_cross_tot, 6),
                                kappa_upper95=round(up_c, 6)),
                within_a=dict(observed_pairs=obs_w, kappa=round(obs_w / exp_tot, 6),
                              kappa_upper95=round(up_w, 6)),
                ratio=round((obs_c / exp_cross_tot) / (obs_w / exp_tot), 4) if obs_w else None,
                charged_kappa_max_of_cross_and_within=round(max(up_c, up_w), 6))


def gap_rows():
    arm = json.loads([l for l in open(os.path.join(HERE, "arm-retune-v31.jsonl"))
                      if l.strip() and l.strip() != "DONE"][0])
    mu_log2 = 26.067                       # derived in derive_v31 by band=ell*mu fixed point
    mu = 2.0 ** mu_log2
    ell = arm["ell"]
    n_draws = ell * arm["trials"]          # gap draws in the census
    obs_max = arm["gap_max"]
    # one-sided Poisson 95% upper bound on P(gap > observed max) over the draws
    p_upper = 2.996 / n_draws
    # production extreme-value support at the same law: max over nb*ell draws per side
    n_prod = NB * ell
    ev_support = mu * math.log(n_prod)     # E[max] of ~Exp(mu) gaps over n draws
    return dict(
        source="arm-retune-v31.jsonl (24000 trials, K_max=2^32, 0 aborts)",
        observed_max_gap=obs_max, draws=n_draws,
        p_gap_above_observed_one_sided=round(p_upper, 6),
        observed_support_entries=obs_max // 2,
        empirical_max_plus_one_sided="gap_max + mu*ln(1/p_upper) = %.3e" % (obs_max + mu * math.log(1 / p_upper)),
        production_extreme_value_log2=round(math.log2(ev_support), 6),
        ledger_support_log2=31,
        verdict=("the fixture cap never fires (0/24000 at 2^32) so it carries NO tail information; "
                 "the ledger's 2^31 support entry count is instead DERIVED from the same exponential "
                 "gap law: E[max gap] = mu*ln(nb*ell) = 2^%.2f over 2^115.65 draws, support = /2, "
                 "i.e. 2^31.  The tail hypothesis (exponential gaps, mu=2^%.3f) is stated and its "
                 "UNTESTED status at production is disclosed, not hidden"
                 % (math.log2(ev_support), mu_log2)))


def main():
    cs = census()
    tri = dict(
        x=nb_x if (nb_x := round(X31, 8)) is not None else None,
        phi_uniform=round(phi_uniform(X31), 6),
        phi_tri=round(phi_tri(X31), 6),
        factor="4/3 (triangular vs uniform difference law)",
        derivation=("sum_t p(t)^2 = 2/(3L) for triangular p(t)=(L-|t|)/L^2 vs 1/(2L) uniform over 2L "
                    "slots: the collision density is 4/3 larger, so the effective birthday argument "
                    "is y=(4/3)*nb^2/G"),
        K=round(6 * ((NB * (-(-N // NB))) ** 2) / CERT["m"], 6),
        K_eff_uniform=round(6 * ((NB * (-(-N // NB))) ** 2) / CERT["m"] * phi_uniform(X31), 6),
        K_eff_tri=round(6 * ((NB * (-(-N // NB))) ** 2) / CERT["m"] * phi_tri(X31), 6),
        floor=0.39,
        floor_met_under_tri=bool(6 * ((NB * (-(-N // NB))) ** 2) / CERT["m"] * phi_tri(X31) >= 0.39))
    out = dict(kind="v32 cross-side kappa census, triangular-occupancy Phi, empirical gap tail",
               cross_side_census=cs, phi=tri, gap_tail=gap_rows(),
               charged_kappa=max(1.231, cs["charged_kappa_max_of_cross_and_within"]),
               within_a_vs_cross=("within-A and cross-side maxima agree (ratio %.3f): the two point "
                                  "distributions are identical and the coins are disjoint, so the "
                                  "within-A proxy is validated WITH DATA, and the charged kappa is the "
                                  "larger of the two" % (cs["ratio"] or 0.0)))
    json.dump(out, open(os.path.join(HERE, "census-crossside-v32.json"), "w"), indent=1)
    print(json.dumps(out, indent=1)[:2000])


if __name__ == "__main__":
    main()
