#!/usr/bin/env python3
"""v19 challenge-faithful H-membership key census.

Fixes v18 flaw 1: the partX census compared UNSCALED points from ORDINARY prime
products.  Here both sides are H-members (H = d-th powers in Z_r^*, |H| = m,
index d), and the giant side carries the genuine challenge multiplier v = s^d
for s uniform in Z_r^*, applied as a group multiplication on the point.

Two rows are produced from one generator:
  * faithful : both lists in H, giant side scaled by v, averaged over challenges
  * control  : the v18 partX object (ordinary primes, no H filter, no v)
so the biased rate of the old census is visible next to the faithful one.
"""
import json, math, os, random, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod
from toy import jmul, jaffine, prod_primes, is_prime, A_LO, A_HI, W_LO, W_HI

P, R = keymod.P, keymod.R
D = 2 ** 20
M = (R - 1) // D
assert (R - 1) % D == 0 and M * D == R - 1


def h_elem(rng):
    """Uniform element of H = image of the d-th power map on Z_r^*.

    The map x -> x^d is a group homomorphism with kernel of size gcd(d,r-1)=d,
    so it is exactly d-to-1 onto H and uniform input gives uniform H output.
    """
    return pow(rng.randrange(2, R), D, R)


def h_list(rng, n):
    out = []
    seen = set()
    while len(out) < n:
        h = h_elem(rng)
        if h not in seen and h != 0:
            seen.add(h)
            out.append(h)
    return out


def y_of(x):
    return jaffine(jmul(x % R, (keymod.GX, keymod.GY, 1), P), P)[1]


def keys_of(xs, widths):
    """key at every width for each scalar; one point per scalar."""
    yk = [keymod.even_member(y_of(x), P) for x in xs]
    return {k: [y & ((1 << k) - 1) for y in yk] for k in widths}


def collide(ka, kb):
    ca, cb = Counter(ka), Counter(kb)
    return sum(c * cb[k] for k, c in ca.items() if k in cb)


def rows_for(xs, ys, widths):
    out = []
    pairs = len(xs) * len(ys)
    kx = keys_of(xs, widths)
    ky = keys_of(ys, widths)
    for k in widths:
        obs = collide(kx[k], ky[k])
        pred = pairs / 2.0 ** (k - 1)          # image 2^(k-1): even rep low bit is 0
        out.append(dict(width=k, image_log2=k - 1, pairs=pairs, observed=obs,
                        predicted=round(pred, 3),
                        ratio=round(obs / pred, 6) if pred >= 1 else None,
                        charged_factor=round(max(1.0, (obs + 3 * math.sqrt(obs + 1)) / pred), 6)
                        if pred >= 1 else None))
    return out


def _memb_worker(seed):
    """Module-level so multiprocessing can pickle it: H-density sample block."""
    rg = random.Random(seed)
    acc = 0
    for _ in range((1 << 20) // 16):
        x = rg.randrange(A_LO, A_HI)
        if pow(x, M, R) == 1:
            acc += 1
    return acc


def census():
    rng = random.Random(20261024)
    sizes = (64, 32)
    nw = 64
    widths = (10, 12, 14, 16, 18, 20, 22)
    n_ch = 6

    # ---------------- faithful object: both factors in H, products in H
    A = h_list(rng, 2 * max(sizes))
    W = h_list(rng, nw)
    ax, ay = A[:max(sizes)], A[max(sizes):2 * max(sizes)]
    assert not (set(ax) & set(ay))
    assert all(pow(a, M, R) == 1 for a in A) and all(pow(w, M, R) == 1 for w in W), "H filter"

    faithful = []
    for n in sizes:
        xs = [a * w % R for a in ax[:n] for w in W]
        ys = [a * w % R for a in ay[:n] for w in W]
        assert all(pow(x, M, R) == 1 for x in xs + ys), "products must stay in H"
        base = dict(records_per_side=len(xs), pairs=len(xs) * len(ys), walk=nw,
                    h_members=True, scaled=True)
        # challenge loop: v = s^d, s uniform in Z_r^*  => v uniform in H
        per_ch = []
        for _ in range(n_ch):
            s = rng.randrange(2, R)
            v = pow(s, D, R)
            assert pow(v, M, R) == 1
            gy = [v * y % R for y in ys]
            rws = rows_for(xs, gy, widths)
            per_ch.append(dict(s_bits=s.bit_length(), v_bits=v.bit_length(),
                               rows=rws))
        # average observed over challenges, per width
        avg = []
        for i, k in enumerate(widths):
            obss = [c["rows"][i]["observed"] for c in per_ch]
            pred = per_ch[0]["rows"][i]["predicted"]
            avg.append(dict(width=k, image_log2=k - 1, pairs=base["pairs"],
                            observed_mean=round(sum(obss) / len(obss), 3),
                            observed_min=min(obss), observed_max=max(obss),
                            predicted=pred,
                            ratio_mean=round((sum(obss) / len(obss)) / pred, 6) if pred >= 1 else None,
                            charged_factor=round(max(1.0, (max(obss) + 3 * math.sqrt(max(obss) + 1)) / pred), 6)
                            if pred >= 1 else None))
        base["exact_rows"] = avg
        base["challenges"] = [dict(rows=c["rows"]) for c in per_ch]
        faithful.append(base)

    # ---------------- control object: the v18 partX construction, verbatim
    ab = prod_primes(A_LO, A_HI, 2 * max(sizes))
    walk = prod_primes(W_LO, W_HI, nw)
    cx, cy = ab[:max(sizes)], ab[max(sizes):2 * max(sizes)]
    control = []
    for n in sizes:
        xs = [a * w for a in cx[:n] for w in walk]
        ys = [a * w for a in cy[:n] for w in walk]
        control.append(dict(records_per_side=len(xs), pairs=len(xs) * len(ys), walk=nw,
                            h_members=False, scaled=False,
                            exact_rows=rows_for(xs, ys, widths)))

    # ---------------- H-membership density row (band-limited, exact counts)
    json.dump(dict(faithful=faithful, control=control, d=D, m=M, widths=list(widths),
                   n_walk=nw, n_challenges=n_ch, sizes=list(sizes)),
              open(HERE + "/census-v19-partial.json", "w"), indent=1)
    nsamp = 1 << 20
    span = (A_HI - A_LO) // nsamp
    assert span > 0
    import multiprocessing as mp
    with mp.Pool(16) as pool:
        hits = sum(pool.map(_memb_worker, range(16)))
    memb = dict(band=[A_LO, A_HI], samples=nsamp, accepted=hits,
                rejected_fraction=round(1 - hits / nsamp, 8),
                expected_density_log2=-20,
                note="H has index d=2^20 in Z_r^*, so a uniform band element lies in H "
                     "with probability 1/d; ordinary primes as used by prod_primes are "
                     "NOT filtered, which is the v18 flaw")

    def kmax(rs):
        c = []
        for r0 in rs:
            v = r0.get("ratio_mean", r0.get("ratio"))
            if r0["predicted"] >= 16 and v is not None:
                c.append(v)
        return (round(max(c), 6), round(min(c), 6), len(c)) if c else (None, None, 0)

    fr_max, fr_min, fr_n = kmax([r0 for sz in faithful for r0 in sz["exact_rows"]])
    co_max, co_min, co_n = kmax([r0 for sz in control for r0 in sz["exact_rows"]])
    ctl_ratio = {}
    for sz in control:
        for r0 in sz["exact_rows"]:
            ctl_ratio[(sz["records_per_side"], r0["width"])] = r0["ratio"]
    dev = []
    for sz in faithful:
        for r0 in sz["exact_rows"]:
            c = ctl_ratio.get((sz["records_per_side"], r0["width"]))
            if c and r0["ratio_mean"] is not None:
                dev.append(dict(records_per_side=sz["records_per_side"], width=r0["width"],
                                faithful=r0["ratio_mean"], control=c,
                                delta=round(r0["ratio_mean"] - c, 6)))
    return dict(d=D, m=M, r_bits=R.bit_length(), widths=list(widths), n_walk=nw,
                n_challenges=n_ch, sizes=[sz["records_per_side"] for sz in faithful],
                faithful=faithful, control=control, membership=memb,
                faithful_ratio_max=fr_max, faithful_ratio_min=fr_min, faithful_cells=fr_n,
                control_ratio_max=co_max, control_ratio_min=co_min, control_cells=co_n,
                faithful_vs_control=dev,
                max_abs_delta=round(max(abs(x["delta"]) for x in dev), 6) if dev else None,
                transfer=dict(
                    stat="per-pair collision probability / 2^-167 (the measured ratio)",
                    argument=("conditional on the challenge, giant element v*y_j is uniform over H "
                              "because v is uniform over H and H is a group; hence for every pair "
                              "(i,j) Pr[key_i == key(v*y_j)] = 2^-167 in the measured sense, and "
                              "E_v[#collisions] = pairs/2^167 by linearity over pairs. Markov on the "
                              "challenge-random variable then applies at the list level. Residual "
                              "dependence stated: pairs sharing the same j are correlated through the "
                              "common factor v; the census averages over challenges rather than "
                              "assuming independence."),
                    s_to_d_regularity=("x -> x^d on Z_r^* has kernel size gcd(d,r-1)=d=2^20, so it is "
                                       "exactly d-to-1 onto H and s uniform => v=s^d uniform on H")),
                note="faithful = H-members and challenge-scaled giant side; control = v18 partX "
                     "object (ordinary primes, no scaling).")


if __name__ == "__main__":
    out = census()
    json.dump(out, open(HERE + "/census-v19.json", "w"), indent=1)
    print("faithful ratio max %s min %s over %d cells" % (out["faithful_ratio_max"],
                                                          out["faithful_ratio_min"],
                                                          out["faithful_cells"]))
    print("control  ratio max %s min %s over %d cells" % (out["control_ratio_max"],
                                                          out["control_ratio_min"],
                                                          out["control_cells"]))
    print("max |faithful-control| ratio delta %s" % out["max_abs_delta"])
    print("H density: %d/%d accepted, rejected fraction %s" % (out["membership"]["accepted"],
                                                               out["membership"]["samples"],
                                                               out["membership"]["rejected_fraction"]))
