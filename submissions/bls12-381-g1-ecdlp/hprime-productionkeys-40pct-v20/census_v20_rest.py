#!/usr/bin/env python3
"""v20 part 2: seed-coin gap census, two-geometry coverage census, routing-load
census, and the band-limited H-PRIME control for the key census.

Run after `python3 census_v20.py key`.
"""
import json, math, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import census_v20 as C
import keymod
from toy import is_prime, grid_census, prime_field, g_field, hprimes_range

P, R = keymod.P, keymod.R
D = 2 ** 20
M = (R - 1) // D


# ------------------------------------------------------------------ gaps
def hprime_gaps(p, dd, seed, n_gap, start_lo, start_hi):
    """Consecutive H-primes after a seed-derived start: gaps on the ACTUAL list.

    The start point is a coin of the construction (hash of the seed), so the
    per-seed gap maximum is a seed-random variable - which is the object the
    verdict asks for, instead of a fixed-list allowance.
    """
    m = (p - 1) // dd
    rng = random.Random("gap|%d|%d|%d" % (p, dd, seed))
    start = rng.randrange(start_lo, start_hi) | 1
    q, out = start, []
    while len(out) < n_gap + 1:
        if is_prime(q) and pow(q, m, p) == 1:
            out.append(q)
        q += 2
    gaps = [out[i + 1] - out[i] for i in range(n_gap)]
    mean = sum(gaps) / len(gaps)
    mx = max(gaps)
    # KS against the geometric law with the measured mean
    srt = sorted(gaps)
    ks = 0.0
    n = len(srt)
    for i, g in enumerate(srt):
        cdf = 1.0 - math.exp(-g / mean)
        ks = max(ks, abs(cdf - (i + 1) / n), abs(cdf - i / n))
    return dict(seed=seed, start_bits=start.bit_length(), n=n, mean_gap=round(mean, 3),
                mean_gap_bits=round(math.log2(mean), 4), max_gap=mx,
                max_gap_bits=round(math.log2(mx), 4), max_over_mean=round(mx / mean, 4),
                ks=round(ks, 5), ks_cut=round(1.36 / math.sqrt(n), 5), ks_pass=ks < 1.36 / math.sqrt(n),
                first=out[0], last=out[-1])


def _gap_seed(task):
    bits, dd, seed, n_gap = task
    p = _gap_field(bits, dd)
    return dict(scale_bits=bits, d=dd, p=p,
                **hprime_gaps(p, dd, seed, n_gap, 1 << (bits - 2), 1 << (bits - 1)))


def gap_part(seeds=64):
    import multiprocessing as mp
    N_GAP = 401
    tasks = [(bits, 8, s, N_GAP) for bits in (36, 38) for s in range(seeds)]
    ctx = mp.get_context("fork")
    with ctx.Pool(min(16, ctx.cpu_count() - 1)) as pool:
        rows = pool.map(_gap_seed, tasks, chunksize=2)
    by = {}
    for r in rows:
        by.setdefault((r["scale_bits"], r["d"]), []).append(r)
    agg = []
    THRESH = 2 ** 7.7               # production requirement: max gap <= 2^7.7 * mean
    for (bits, dd), rs in sorted(by.items()):
        mx = max(r["max_over_mean"] for r in rs)
        agg.append(dict(scale_bits=bits, d=dd, seeds=len(rs),
                        mean_gap_bits=[min(r["mean_gap_bits"] for r in rs), max(r["mean_gap_bits"] for r in rs)],
                        max_over_mean_max=mx, max_over_mean_median=sorted(r["max_over_mean"] for r in rs)[len(rs) // 2],
                        ks_pass_fraction=round(sum(r["ks_pass"] for r in rs) / len(rs), 4),
                        production_ratio_bound=round(math.log2(THRESH), 4),
                        seeds_exceeding_production_ratio=sum(1 for r in rs if r["max_over_mean"] > THRESH),
                        seeds_exceeding_margin=sum(1 for r in rs if r["max_over_mean"] > THRESH / 16)))
    fails = sum(1 for r in rows if r["max_over_mean"] > THRESH)
    return dict(rows=rows, aggregate=agg, thresh_ratio=THRESH,
                thresh_ratio_log2=round(math.log2(THRESH), 4), seeds=len(rows),
                scale_pairs=len(agg), failures=fails,
                empirical_fail_rate=round(fails / len(rows), 6),
                rule_of_three=round(3 / len(rows), 6),
                rule_of_three_log2=round(math.log2(3 / len(rows)), 4),
                note=("Each row is ONE independent seed = one construction coin: the start of "
                      "the consecutive H-prime band is seed-derived, so the gap maximum is a "
                      "seed-random variable. Failures count seeds whose max/mean ratio exceeds "
                      "the production requirement 2^7.7; the empirical rate and the rule-of-three "
                      "upper bound 3/128 are the debits used in the ledger."))


def _gap_field(bits, dd):
    x = (1 << (bits - 1)) | 1
    while True:
        if (x - 1) % (6 * dd) == 0 and is_prime(x):
            return x
        x += 2


# ------------------------------------------------------- band-limited control
def bandlimited_control(dd=16, nb=64, nw=64, widths=(10, 12, 14, 16, 18, 20, 22, 24)):
    """Production curve, SMALL index d so H-primes are affordable: the attack's
    actual objects (band-limited primes inside H) against uniform H elements,
    same width grid.  Answers 'you sampled uniform H factors, not the attack's
    band-limited H-prime factors' by measuring both."""
    m = (R - 1) // dd
    assert (R - 1) % dd == 0
    rng = random.Random(20261025)
    lo, hi = 1 << 77, 1 << 80
    primes, q = [], lo | 1
    need = 2 * nb + nw
    while len(primes) < need and q < hi:
        if is_prime(q) and pow(q, m, R) == 1:
            primes.append(q)
        q += 2
    assert len(primes) == need, (len(primes), need)
    bl_ax, bl_ay, bl_w = primes[:nb], primes[nb:2 * nb], primes[2 * nb:]
    bl_x = [a * w % R for a in bl_ax for w in bl_w]
    bl_y = [a * w % R for a in bl_ay for w in bl_w]
    assert all(pow(x, m, R) == 1 for x in bl_x)

    def hel(n):
        o = []
        for _ in range(n):
            o.append(pow(rng.randrange(2, R), dd, R))
        return o
    u_ax, u_ay, u_w = hel(nb), hel(nb), hel(nw)
    un_x = [a * w % R for a in u_ax for w in u_w]
    un_y = [a * w % R for a in u_ay for w in u_w]

    pairs = len(bl_x) * len(bl_y)
    yb_x, yb_y = C.ys_all(bl_x), C.ys_all(bl_y)
    yu_x, yu_y = C.ys_all(un_x), C.ys_all(un_y)
    bl_rows = C.rows_for(yb_x, yb_y, widths, pairs)
    un_rows = C.rows_for(yu_x, yu_y, widths, pairs)
    return dict(d=dd, m=m, band=[lo, hi], nb=nb, nw=nw, records_per_side=len(bl_x),
                pairs=pairs, widths=list(widths),
                bandlimited=bl_rows, uniform_h=un_rows,
                max_ratio_bandlimited=max((r["ratio"] for r in bl_rows if r["ratio"]), default=None),
                max_ratio_uniform=max((r["ratio"] for r in un_rows if r["ratio"]), default=None),
                note=("d = 2^4 here so H-primes are affordable on the PRODUCTION curve: this row "
                      "is a distribution control for the key map, not a production-index measurement."))


# ------------------------------------------------------------ routing loads
ROUTE_WIDTHS = (9, 12, 15)


def routing_part(records=4096, walk=64):
    """Fixture routing table: measured per-class loads against the uniform
    formula, and prefix-equal pair counts against pairs/2^rb."""
    rng = random.Random(99)
    nb = records // walk
    A = C.h_list(rng, 2 * nb)
    W = C.h_list(rng, walk)
    xs = [a * w % R for a in A[:nb] for w in W]
    ys = [a * w % R for a in A[nb:] for w in W]
    yx, yy = C.ys_all(xs), C.ys_all(ys)
    kx = [keymod.key_scalar(y, P) for y in yx]
    ky = [keymod.key_scalar(y, P) for y in yy]
    rows = []
    for rb in ROUTE_WIDTHS:
        cls_x = [k >> (keymod.KEY_BITS - rb) for k in kx]
        cls_y = [k >> (keymod.KEY_BITS - rb) for k in ky]
        cx = {}
        for c in cls_x:
            cx[c] = cx.get(c, 0) + 1
        cy = {}
        for c in cls_y:
            cy[c] = cy.get(c, 0) + 1
        ncls = 1 << rb
        mean_x = len(cls_x) / ncls
        pref = sum(c * cy.get(k, 0) for k, c in cx.items())
        pred_pairs = len(cls_x) * len(cls_y) / ncls
        rows.append(dict(route_bits=rb, classes=ncls, records_side=len(cls_x),
                         mean_load_formula=round(mean_x, 6), max_load_measured=max(cx.values()),
                         max_over_formula=round(max(cx.values()) / mean_x, 6) if mean_x else None,
                         prefix_pairs=pref, prefix_predicted=round(pred_pairs, 6),
                         prefix_ratio=round(pref / pred_pairs, 6) if pred_pairs >= 1 else None,
                         ucb_prefix_ratio=round(max(1.0, (pref + 3 * math.sqrt(pref + 1)) / pred_pairs), 6)
                         if pred_pairs >= 1 else None))
    return dict(rows=rows, route_bits_production=keymod.ROUTE_BITS,
                sort_bits_production=keymod.SORT_BITS,
                max_load_over_formula=max(r["max_over_formula"] for r in rows),
                max_prefix_ratio=max((r["ucb_prefix_ratio"] for r in rows if r["ucb_prefix_ratio"]), default=None),
                note=("Route widths are scaled down so the fixture has measurable occupancy "
                      "(O(1) records per class); the ratio max/mean and the prefix-pair ratio "
                      "are the measured deviation factors charged at production."))


# ------------------------------------------------------- coverage, 2 scales
COV_D = (8, 64)


def coverage_part():
    scales = [dict(bits=23, nb=8), dict(bits=25, nb=16)]
    cells = []
    for sc in scales:
        r = prime_field(sc["bits"])
        for dd in COV_D:
            c = grid_census(r, dd, 0.51, nb=sc["nb"])
            mm = c["m"]
            rowB = (c["predicted_distinct_eff"] - c["distinct_covered"]) / mm
            c["scale"] = dict(bits=sc["bits"], nb=sc["nb"])
            c["rowB"] = rowB
            c["rowA"] = (c["predicted_distinct"] - c["distinct_covered"] - c["diagonal_repetitions"]) / mm
            cells.append(c)
    # stability of |rowB| across the two scales at matched d
    stab = []
    for dd in COV_D:
        a = [c for c in cells if c["d"] == dd and c["scale"]["bits"] == 23][0]
        b = [c for c in cells if c["d"] == dd and c["scale"]["bits"] == 25][0]
        stab.append(dict(d=dd, small_abs=abs(a["rowB"]), large_abs=abs(b["rowB"]),
                         ratio=round(abs(b["rowB"]) / abs(a["rowB"]), 6) if a["rowB"] else None,
                         bounded_by=round(max(abs(a["rowB"]), abs(b["rowB"])), 8),
                         stable=abs(b["rowB"]) <= 2 * abs(a["rowB"]) + 1e-9))
    return dict(cells=cells, stability=stab, K_target=0.51, d_values=list(COV_D),
                max_abs_rowB=max(abs(c["rowB"]) for c in cells),
                max_abs_rowB_scale2=max(abs(c["rowB"]) for c in cells if c["scale"]["bits"] == 25),
                max_abs_rowA=max(abs(c["rowA"]) for c in cells),
                cells_count=len(cells), bases_per_side=[c["nb"] for c in cells],
                stable_count=sum(1 for s in stab if s["stable"]),
                note=("Two geometric scales, same SHAPE: list sizes doubled and the field "
                      "doubled in bits, so K_cand is matched and pool:walk geometry is "
                      "preserved. |rowB| is reported at both scales; no cross-scale "
                      "extrapolation is claimed beyond the stability table."))


if __name__ == "__main__":
    out = {}
    print("gaps...", flush=True)
    out["gap"] = gap_part(64)
    print("gap aggregate", [(a["scale_bits"], a["max_over_mean_max"], a["ks_pass_fraction"],
                             a["seeds_exceeding_production_ratio"]) for a in out["gap"]["aggregate"]], flush=True)
    print("coverage...", flush=True)
    out["coverage"] = coverage_part()
    print("coverage", [(c["d"], c["scale"]["bits"], round(c["rowB"], 8)) for c in out["coverage"]["cells"]], flush=True)
    print("routing...", flush=True)
    out["routing"] = routing_part()
    print("band-limited control...", flush=True)
    out["bandlimited"] = bandlimited_control()
    json.dump(out, open(HERE + "/census-rest-v20.json", "w"), indent=1)
    print("wrote census-rest-v20.json")
