#!/usr/bin/env python3
"""v22 evidence A: the POPULATION N(k), enumerated exactly.

The v21 verdict: "exact_rows ignores its m argument and constructs nw=Counter(ky)
from one challenge-scaled giant list.  This histogram sums to |Y|, not m, and is
not the defined population N(k)=|{h in H:key([h]P0)=k}|."

Here H is EXHAUSTIVELY enumerated: the d-th power subgroup of F_r^* has order
m = (r-1)/d, we list ALL m elements, compute [h]P0 for every one of them, and
build the exact histogram

    N_w(k) = |{h in H : key_w([h]P0) = k}|,     sum_k N_w(k) == m   (asserted)

for every width w.  Nothing here is a sample or an average: the whole population
is a finite list that fits in memory, and sum_k N_w(k) == m is the machine check
that the histogram IS the population.

With a REALIZED baby list B_w(k) the exact challenge-space expectation for that
fixed list is

    E_w = (|Y|/m) * sum_k B_w(k) N_w(k)

and it is compared, per challenge v in H, against the REALIZED count

    C_w(v) = sum_k B_w(k) G_w(k),   G_w = giant key histogram of {v*y_j}

so the row is a ratio E/C(v) over independent challenges, plus a weights row:
the population mass sum_{k:B(k)>0} N(k)/m against the uniform-weights value
|occupied keys| / |keyspace|.  If they disagree, the |X||Y|/D shortcut is NOT
valid on occupied keys and only the measured weights are usable.
"""
import json, math, os, random, sys
from collections import Counter
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod
from toy import jmul, jaffine, is_prime, g_curve_field

WIDTHS = tuple(range(10, 31))


def h_group(bits, dd):
    """r prime, d | r-1, m=(r-1)/d, 6|m; return (r, dd, m, [all m elements of H])."""
    r = (1 << (bits - 1)) | 1
    while True:
        if ((r - 1) % dd == 0) and ((r - 1) // dd) % 6 == 0 and is_prime(r):
            break
        r += 2
    m = (r - 1) // dd
    fac, t, q = set(), r - 1, 2
    while q * q <= t:
        while t % q == 0:
            fac.add(q)
            t //= q
        q += 1
    if t > 1:
        fac.add(t)
    g = 2
    while not all(pow(g, (r - 1) // p, r) != 1 for p in fac):
        g += 1
    step = pow(g, dd, r)                       # generator of the image subgroup H
    h, H = 1, []
    for _ in range(m):
        H.append(h)
        h = h * step % r
    assert len(set(H)) == m and h == 1
    assert all(pow(x, m, r) == 1 for x in H[::max(1, len(H) // 97)])
    assert pow(step, m // 2, r) != 1           # generator order is exactly m
    return r, dd, m, H


def _ys_block(args):
    chunk, p = args
    return [jaffine(jmul(h % p, (keymod.GX, keymod.GY, 1), p), p)[1] for h in chunk]


def ys_all(scalars, p, nproc=16):
    import multiprocessing as mp
    blocks = [(scalars[i:i + 2048], p) for i in range(0, len(scalars), 2048)]
    ctx = mp.get_context("fork")
    with ctx.Pool(min(nproc, ctx.cpu_count() - 1)) as pool:
        return [y for b in pool.map(_ys_block, blocks, chunksize=1) for y in b]


def key_list(yv, w, p):
    mask = (1 << w) - 1
    return [keymod.even_member(y, p) & mask for y in yv]


def fixture(r, m, H, p, nb, ny, challenges, widths=WIDTHS, seed=1, tag=""):
    """Exact N_w population + realized baby list + realized C_w(v) per challenge."""
    rng = random.Random(seed)
    yH = ys_all(H, p)
    yof = {h: yH[i] for i, h in enumerate(H)}
    pop, sumN = {}, {}
    for w in widths:
        c = Counter(key_list(yH, w, p))
        assert sum(c.values()) == m, (w, sum(c.values()), m)
        pop[w], sumN[w] = c, sum(c.values())
    baby, giant = rng.sample(H, nb), rng.sample(H, ny)
    rows = []
    for w in widths:
        B = Counter(key_list([yof[a] for a in baby], w, p))
        dot = sum(c * pop[w][k] for k, c in B.items())
        Dimg = 1 << (w - 1)
        Eexact = F(ny, m) * dot
        base = F(nb * ny, Dimg)
        mass = F(sum(pop[w][k] for k in B), m)          # population mass on occupied keys
        unif = F(len(B), Dimg)                          # uniform-weights value
        realized = []
        for _ in range(challenges):
            v = rng.choice(H)
            G = Counter(key_list([yof[(v * y) % r] for y in giant], w, p))
            realized.append(sum(c * G[k] for k, c in B.items()))
        pos = [c for c in realized if c]
        rows.append(dict(
            width=w, image_log2=w - 1, pairs=nb * ny, baseline=round(float(base), 6),
            density_ratio=round(m / Dimg, 3),
            sum_Nk=sumN[w], population_size=m, dot=int(dot),
            E_exact=str(Eexact), E_exact_float=round(float(Eexact), 6),
            E_over_baseline=round(float(Eexact / base), 6),
            realized=realized,
            ratio_mean=round(sum(float(Eexact / c) for c in pos) / len(pos), 6) if pos else None,
            ratio_min=round(min(float(Eexact / c) for c in pos), 6) if pos else None,
            ratio_max=round(max(float(Eexact / c) for c in pos), 6) if pos else None,
            n_zero_realized=sum(1 for c in realized if c == 0),
            mass_on_occupied=str(mass), mass_uniform=str(unif),
            weight_ratio=round(float(mass / unif), 6) if unif else None,
            distinct_B=len(B), max_mult_B=max(B.values()),
            distinct_N=len(pop[w]), max_mult_N=max(pop[w].values())))
    return dict(tag=tag, d=None, m=m, m_log2=round(math.log2(m), 3), nb=nb, ny=ny,
                challenges=challenges, curve_bits=p.bit_length(),
                h_enumerated=len(H), rows=rows)


def main():
    out = dict(identity="E_v[C] = (|Y|/m) * sum_k B(k) N(k), N(k) the EXACT enumerated population",
               enumeration="every h in H visited; sum_k N_w(k)==m asserted per width per fixture",
               fixtures=[])
    # production has m/D = 2^67.86 (population far DENSER than the keyspace), so the
    # fixtures sweep the density ratio m/2^(w-1) in {2^2..2^12}: width is set RELATIVE
    # to m, never absolutely.
    cfgs = [(14, 4, [(32, 32), (64, 64), (128, 64)]),
            (16, 4, [(64, 64), (128, 64), (128, 128)]),
            (18, 4, [(64, 64), (128, 128), (256, 128)])]
    p = g_curve_field(26)
    for bits, dd, shapes in cfgs:
        r, dd, m, H = h_group(bits, dd)
        mb = round(math.log2(m))
        widths = tuple(mb - k for k in (14, 12, 10, 8, 6, 4, 2) if mb - k >= 5)
        print("enum H: r bits=%d m=%d (2^%.2f) |H|=%d widths=%s"
              % (r.bit_length(), m, math.log2(m), len(H), widths), flush=True)
        for nb, ny in shapes:
            fx = fixture(r, m, H, p, nb, ny, 3, widths=widths,
                         seed=bits * 1000 + nb * 7 + ny, tag="r%d_nb%d_ny%d" % (bits, nb, ny))
            fx.update(field_r_bits=r.bit_length(), d=dd)
            out["fixtures"].append(fx)
            print("  fixture %s rows=%d" % (fx["tag"], len(fx["rows"])), flush=True)
            json.dump(out, open(HERE + "/population-v22.json", "w"), indent=1)
    cells = [(f, rw) for f in out["fixtures"] for rw in f["rows"] if rw["ratio_mean"]]
    ratios = [rw["ratio_mean"] for _, rw in cells]
    wr = [rw["weight_ratio"] for _, rw in cells if rw["weight_ratio"]]
    buckets = {}
    for _, rw in cells:
        b = int(round(math.log2(rw["density_ratio"])))
        buckets.setdefault(b, []).append(rw)
    trend = [dict(density_ratio_log2=b,
                  ratio_mean=round(sum(x["ratio_mean"] for x in v) / len(v), 4),
                  ratio_min=round(min(x["ratio_mean"] for x in v), 4),
                  ratio_max=round(max(x["ratio_mean"] for x in v), 4),
                  weight_ratio_mean=round(sum(x["weight_ratio"] for x in v) / len(v), 6),
                  cells=len(v)) for b, v in sorted(buckets.items())]
    out["summary"] = dict(
        fixtures=len(out["fixtures"]), cells_scored=len(cells), trend=trend,
        ratio_min=round(min(ratios), 4), ratio_max=round(max(ratios), 4),
        ratio_mean=round(sum(ratios) / len(ratios), 4),
        weight_ratio_min=round(min(wr), 4), weight_ratio_max=round(max(wr), 4),
        weight_ratio_cells=len(wr),
        sum_Nk_asserted_all=True,
        notes=("N(k) is enumerated, never sampled: every h in H contributes exactly one key at "
               "every width and sum_k N(k)==m is asserted per width.  ratio = E/C(v), E the "
               "population expectation for the FIXED realized baby list, C(v) the realized "
               "cross-side count of an independent challenge v.  weight_ratio = [sum_{k:B(k)>0} "
               "N(k)/m] / [|occupied|/2^(w-1)]: 1.0 would make the uniform-weights shortcut exact."))
    json.dump(out, open(HERE + "/population-v22.json", "w"), indent=1)
    print("SUMMARY", json.dumps(out["summary"], indent=1))


if __name__ == "__main__":
    main()
