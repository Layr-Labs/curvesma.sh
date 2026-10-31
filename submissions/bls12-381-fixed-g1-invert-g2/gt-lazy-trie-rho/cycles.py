"""Part F of experiment.py (proof.md section 3.3): bounds on useless closed walks of the R = 24 Frobenius-class walk.
Loaded by experiment.py with runpy (split out for the 64 KiB per-file limit); standard library only.
Standalone: python3 -I cycles.py  (prints part F for m0 = 12)."""
import json, math
from fractions import Fraction

THETA_BITS = 86
R_PROD = 24

# ---------------------------------------------------------------------------------------------------
# Fruitless closed walks at R = 24 (part F)
# ---------------------------------------------------------------------------------------------------
def root_sums(m0):
    """q(n) = P(sum of n iid uniform 12th roots of unity is 0 in Z[zeta_12]) for n <= 2 m0 and Q(n) = max point
    probability for n <= m0, exact (integer counts over 12^n; zeta^k in the basis 1, z, z^2, z^3 with z^4 = z^2 - 1)."""
    Z = [(1, 0, 0, 0)]
    for _ in range(11):
        a, b, c, d = Z[-1]
        Z.append((-d, a, b + d, c))
    assert len(set(Z)) == 12
    dist = {(0, 0, 0, 0): 1}
    q, Q, dists = [Fraction(1)], [Fraction(1)], [dist]
    for n in range(1, m0 + 1):
        nd = {}
        for s, c in dist.items():
            for z in Z:
                t = (s[0] + z[0], s[1] + z[1], s[2] + z[2], s[3] + z[3])
                nd[t] = nd.get(t, 0) + c
        dist = nd
        q.append(Fraction(dist.get((0, 0, 0, 0), 0), 12 ** n))
        Q.append(Fraction(max(dist.values()), 12 ** n))
        dists.append(dist)
    # q(m0 + k), k <= m0: convolve the n = m0 and n = k laws (roots are closed under negation)
    for k in range(1, m0 + 1):
        dk = dists[k]
        q.append(Fraction(sum(c * dk.get(tuple(-x for x in s), 0) for s, c in dist.items()), 12 ** (m0 + k)))
    # Q(n) <= q(2 m0) for n >= 2 m0: P(S = x) = sum_v P_a(v) P_b(x - v) <= sqrt(q(2a) q(2b)) (Cauchy-Schwarz,
    # sum_v P_a(v)^2 = q(2a) by symmetry), a = b = m0, and adding independent terms cannot raise it
    return q, Q

def part_F(m0, Rm=R_PROD, lmax_exact=600):
    """Bound on the probability P_l that the l steps from a given step form a useless closed walk (H2: indices
    uniform on R, twists uniform on 12, independent). A useless walk needs, for every index used, the sum of
    p^(frame) over its occurrences to vanish, i.e. a vanishing sum of n 12th roots of unity (probability q(n),
    at most q(2 m0) for n > 2 m0).
    P_l <= E_multinomial[prod_h g(n_h)], g(0)=1, g(1)=0, g(n) = q(n) exact for n <= 2 m0 and <= q(2 m0) beyond; exact for l <= lmax_exact, Poissonised beyond (P_l <= (E_Pois(l/R) g)^R / P(Pois(l)=l) <= e sqrt(l) (E_Pois g)^R)."""
    q, Q = root_sums(m0)
    m1 = 2 * m0
    gf = [1.0, 0.0] + [float(q[n]) for n in range(2, m1 + 1)]
    Qf = float(q[m1])
    g = lambda n: gf[n] if n <= m1 else Qf
    lnR = math.log(Rm)
    G = [g(n) * math.exp(n * lnR - math.lgamma(n + 1)) for n in range(lmax_exact + 1)]
    poly = [1.0] + [0.0] * lmax_exact
    for _ in range(Rm):
        new = [0.0] * (lmax_exact + 1)
        for i, a in enumerate(poly):
            if a:
                for j in range(lmax_exact + 1 - i):
                    if G[j]:
                        new[i + j] += a * G[j]
        poly = new
    # G[j] underflows to 0 from some j = jz on; add P(some index used >= jz times) <= R C(l, jz) R^-jz
    jz = next((j for j in range(2, lmax_exact + 1) if G[j] == 0.0 and g(j) > 0), lmax_exact + 1)
    lost = lambda l: Rm * math.exp(math.lgamma(l + 1) - math.lgamma(jz + 1) - math.lgamma(l - jz + 1) - jz * lnR) if l >= jz else 0.0
    P_l = {l: (math.exp(math.log(poly[l]) + math.lgamma(l + 1) - 2 * l * lnR) if poly[l] > 0 else 0.0) + lost(l)
           for l in range(2, lmax_exact + 1)}
    # tail l > lmax_exact, on doubling segments [a, 2a): for mu = l/R >= m1 every Poisson point probability
    # P_mu(n), n <= m1, decreases in mu, so E_Pois(mu) g <= sum_{n<=m1} P_{a/R}(n) g(n) + q(m1) on the segment;
    # P_l <= e sqrt(l) (E_Pois g)^R; a trail traverses a closed walk with no distinguished class with
    # probability (1 - theta)^l <= exp(-theta l) (H2)
    assert lmax_exact >= Rm * m1, "lmax_exact must be at least 2 R m0"
    theta = 2.0 ** -THETA_BITS
    # tail_len weights each length by l (<= 2a on the segment): the expected cycle length of such closures per step
    tail, tail_len, a = 0.0, 0.0, lmax_exact + 1
    while a < 2 ** 100:
        mu = a / Rm
        Eg = sum(math.exp(-mu + n * math.log(mu) - math.lgamma(n + 1)) * g(n) for n in range(m1 + 1)) + Qf
        seg = a * math.e * math.sqrt(2 * a) * math.exp(Rm * math.log(Eg) - theta * a)
        tail += seg
        tail_len += 2 * a * seg
        a *= 2
    exact65 = sum(P_l[l] for l in range(65, lmax_exact + 1))
    exact65_len = sum(l * P_l[l] for l in range(65, lmax_exact + 1))
    rate65 = exact65 + tail
    return {"m0": m0, "q": {n: float(q[n]) for n in range(2, m1 + 1)}, "q_2m0_bound_beyond": Qf,
            "P_2": P_l[2], "P_2_equals_1_over_12R": abs(P_l[2] - 1 / (12 * Rm)) < 1e-15,
            "P_l_log2": {l: round(math.log2(P_l[l]), 2) for l in (3, 4, 5, 6, 8, 16, 32, 64, 65, 100, 200, 288, 600)},
            "sum_P_l_3_to_64": sum(P_l[l] for l in range(3, 65)),
            "P_l_by_length_3_to_64": [P_l[l] for l in range(3, 65)],
            "sum_P_l_65_to_%d_log2" % lmax_exact: math.log2(exact65),
            "tail_beyond_%d_weighted_by_no_DP_log2" % lmax_exact: math.log2(tail),
            "rate_useless_walks_longer_than_64_log2": math.log2(rate65),
            "length_weighted_rate_longer_than_64_log2": math.log2(exact65_len + tail_len), "underflow_cut_j": jz}

if __name__ == "__main__":
    out = part_F(12)
    out.pop("P_l_by_length_3_to_64")
    print(json.dumps(out, indent=1))
