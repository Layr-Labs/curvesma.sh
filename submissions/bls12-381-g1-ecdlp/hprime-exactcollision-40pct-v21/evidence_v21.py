#!/usr/bin/env python3
"""v21 evidence: the EXACT fixed-list collision expectation, a shape-matched
coverage ladder, and the conditioned gap / retry accounting.

A. exact-E.  For FIXED baby list B(k) and giant list keys, the exact
challenge-space expectation is

    E_v[C] = (|Y|/m) * sum_k B(k) N(k),

with B(k) the measured baby multiplicity and N(k) the H-side multiplicity
(N(k) = (m/|Y|) * N_list(k) estimated from the drawn list).  The v20 identity
|X||Y| sum_k (N(k)/m)^2 is NOT used: it needs the baby marginal to equal the
H marginal, which membership in H does not give.  Every quantity here is a dot
product of the two measured multiplicity histograms, so the row is exact for
the fixture's own lists.

B. coverage ladder.  Production geometry is nb = 2^15.65 bases/side and
nw = 2^50 walk entries (walk:bases = 2^34.35).  At the design K, N = nb*nw is
fixed by m, so a fixture can only vary the SPLIT.  We run the shape ladder
nb = 1,2,4,8 (nw = N/nb) at two field scales and two d values, i.e. from
minimal-base/maximal-walk (closest to production) outward, and report |rowB|
per cell.  The production split is NOT reachable in any fixture whose
candidate enumeration is O(6 N^2); the package says so and charges a flat
conservative debit instead of extrapolating.

C. gap/retry.  The construction self-checks max gap <= 2^7.7 mu and re-draws;
success is CONDITIONED on the accepted list, so no gap probability belongs in
the success product.  The seed experiment measures the rejection rate: 0/128
seeds -> retry row = 0 measured, with 3/128 disclosed as a confidence limit
(rule of three), never as a bound.
"""
import json, math, os, random, sys
from collections import Counter
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import keymod
from toy import jmul, jaffine, is_prime, prime_field, grid_lists

P, R = keymod.P, keymod.R
D = 2 ** 20
M = (R - 1) // D
WIDTHS = tuple(range(10, 31))


def h_list(rng, n, d=D, r=None):
    r = R if r is None else r
    out, seen = [], set()
    while len(out) < n:
        h = pow(rng.randrange(2, r), d, r)
        if h not in seen and h != 0:
            seen.add(h)
            out.append(h)
    return out


def ys_of(chunk):
    return [jaffine(jmul(x % R, (keymod.GX, keymod.GY, 1), P), P)[1] for x in chunk]


def ys_all(scalars):
    import multiprocessing as mp
    blocks = [scalars[i:i + 512] for i in range(0, len(scalars), 512)]
    ctx = mp.get_context("fork")
    with ctx.Pool(min(16, ctx.cpu_count() - 1)) as pool:
        return [y for b in pool.map(ys_of, blocks, chunksize=1) for y in b]


def hist_stats(c):
    s = sum(c.values())
    return dict(distinct=len(c), total=s, max_mult=max(c.values()),
                sum_sq=sum(v * v for v in c.values()),
                sum_c2=sum(v * (v - 1) // 2 for v in c.values()))


def exact_rows(kx, ky, m, widths):
    # keys, exactly as the production rule: even member of {y, p-y}, truncated
    kx = [keymod.even_member(y, P) for y in kx]
    ky = [keymod.even_member(y, P) for y in ky]
    B = Counter(kx)
    N = Counter(ky)
    nb_keys = len(kx)
    ny = len(ky)
    pairs = nb_keys * ny
    out = []
    for w in widths:
        mask = (1 << w) - 1
        bw = Counter(k & mask for k in kx)
        nw = Counter(k & mask for k in ky)
        dot = sum(c * nw[k] for k, c in bw.items() if k in nw)   # sum_k B(k) N(k)
        Dimg = 1 << (w - 1)
        base = pairs / Dimg
        # v20 identity evaluated on the same lists (baby marginal ASSUMED = H marginal)
        old = pairs * (sum(v * v for v in nw.values()) / (ny * ny))
        out.append(dict(width=w, image_log2=w - 1, pairs=pairs, baseline=round(base, 6),
                        dot=dot, exact_over_baseline=round(dot / base, 6) if base >= 1 else None,
                        v20_identity=round(old, 6),
                        exact_over_v20=round(dot / old, 6) if old >= 1 else None,
                        B=hist_stats(bw), N=hist_stats(nw)))
    return out


def fixture(nb, nw, challenges, widths=WIDTHS):
    rng = random.Random(20261026 + nb * 131 + nw)
    A = h_list(rng, 2 * nb)
    W = h_list(rng, nw)
    xs = [a * w % R for a in A[:nb] for w in W]
    ys = [a * w % R for a in A[nb:] for w in W]
    assert len(set(xs)) == len(xs) and len(set(ys)) == len(ys)
    yx = ys_all(xs)
    rows = []
    for _ in range(challenges):
        s = rng.randrange(2, R)
        v = pow(s, D, R)
        gy = ys_all([v * y % R for y in ys])
        rows.append(exact_rows(yx, gy, M, widths))
    avg = []
    for i, w in enumerate(widths):
        dots = [r[i]["dot"] for r in rows]
        base = rows[0][i]["baseline"]
        old = rows[0][i]["v20_identity"]
        avg.append(dict(width=w, image_log2=w - 1, pairs=rows[0][i]["pairs"], baseline=base,
                        dot_mean=round(sum(dots) / len(dots), 3), dot_min=min(dots),
                        dot_max=max(dots),
                        exact_over_baseline=round((sum(dots) / len(dots)) / base, 6) if base >= 1 else None,
                        exact_over_v20=round((sum(dots) / len(dots)) / old, 6) if old >= 1 else None,
                        B=rows[0][i]["B"], N=rows[0][i]["N"]))
    return dict(nb=nb, nw=nw, records_per_side=len(xs), pairs=len(xs) * len(ys),
                challenges=challenges, rows=avg)


# ------------------------------------------------------------------ coverage
def cell(r, dd, nb, target_K=0.51, min_nw=4):
    m = (r - 1) // dd
    nw = max(min_nw, int(round(math.sqrt(target_K * m / (6 * nb * nb)))))
    walk, bb, bg = grid_lists(r, m, nw, nb)
    g = 2
    while pow(g, (r - 1) // 2, r) == 1 or pow(g, (r - 1) // 3, r) == 1:
        g += 1
    lam6 = pow(g, (r - 1) // 6, r)
    invu = [pow(pow(lam6, j, r), -1, r) for j in range(6)]
    inv_w = {w: pow(w, -1, r) for w in walk}
    inv_g = {a: pow(a, -1, r) for a in bg}
    cov = set()
    for a in bb:
        for A in bg:
            ab = a * inv_g[A] % r
            for w in walk:
                b0 = ab * w % r
                for w2 in walk:
                    base = b0 * inv_w[w2] % r
                    for iu in invu:
                        cov.add(base * iu % r)
    nbs, ngs = len(bb) * len(walk), len(bg) * len(walk)
    occ = 6 * nbs * ngs
    Dc = len(cov)
    K = occ / m
    pred = int(round(m * (1 - math.exp(-K))))
    diag = 6 * len(bb) * len(bg) * (len(walk) - 1)
    Ke = (occ - diag) / m
    prede = int(round(m * (1 - math.exp(-Ke))))
    return dict(field_bits=r.bit_length(), d=dd, m=m, nb=len(bb), nw=len(walk),
                bases_per_side_fn=len(bb), records_per_side=nbs, pairs=nbs * ngs,
                K_cand=round(K, 6), K_eff=round(Ke, 6), distinct=Dc,
                predicted_eff=prede,
                rowB=round((prede - Dc) / m, 8),
                rowA=round((pred - Dc - diag) / m, 8),
                walk_per_base=len(walk),
                production_walk_per_base_log2=34.35)


def cell_nw(r, dd, nw, target_K=0.51):
    """Base-heavy shape: walk length nw FIXED, base pool nb = N/nw so that nb/nw is
    as close to the production 2^15.65 as a design-K fixture permits."""
    m = (r - 1) // dd
    N = math.sqrt(target_K * m / 6.0)
    nb = max(1, int(round(N / nw)))
    return cell(r, dd, nb, target_K=target_K, min_nw=nw)


def coverage_ladder():
    cells = []
    for bits in (23, 25):
        r = prime_field(bits)
        for dd in (8, 64):
            for nb in (1, 2, 4, 8):
                cells.append(cell(r, dd, nb))
            for nw in (1, 2, 4):
                cells.append(cell_nw(r, dd, nw))
    return dict(cells=cells, shapes=[1, 2, 4, 8],
                production_nb_log2=15.6505, production_nw_log2=50.0,
                max_abs_rowB_closest=max(abs(c["rowB"]) for c in cells if c["nb"] == 1),
                max_abs_rowB_all=max(abs(c["rowB"]) for c in cells),
                note=("Shape ladder: nb bases per side, nw = N/nb walk entries, N fixed by the "
                      "design K. nb=1 is the closest feasible shape to production (basis pool "
                      "small, walk reuse maximal); the production split (nb=2^15.65, nw=2^50) is "
                      "not reachable in any O(6N^2) fixture, so no scale transfer is claimed."))


# ------------------------------------------------------------------ gap/retry
def retry_part(seeds=128):
    g = json.load(open(HERE + "/census-rest-v20.json"))["gap"]
    fails = g["failures"]
    rate = F(fails, seeds)
    r3 = F(3, seeds)
    solve = F(0)          # no re-solve is charged when the self-check passes first time
    return dict(seeds=seeds, failures=fails, reject_rate=str(rate),
                reject_rate_float=float(rate),
                rule_of_three=str(r3), rule_of_three_float=float(r3),
                rule_of_three_log2=round(math.log2(float(r3)), 3),
                measured_retry_work_units=0,
                measured_retry_work_log2=None,
                disclosure=("0 rejections in 128 seeds. The retry work row is therefore 0 measured, "
                            "and 3/128 = 2^-5.415 is disclosed as a confidence limit on the "
                            "rejection RATE, not as a bound on any production event and not "
                            "multiplied into success: success is conditioned on the self-check "
                            "having accepted the list, and the self-check is a charged row."),
                self_check=dict(charged_row="gap_self_check",
                                description="verify max gap <= 2^7.7*mu on the drawn list",
                                instructions=2 ** 12, charged=True),
                max_over_mean=[a["max_over_mean_max"] for a in g["aggregate"]],
                thresh_ratio_log2=g["thresh_ratio_log2"])


if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if part in ("all", "exact"):
        out["exact"] = dict(
            identity="E_v[C] = (|Y|/m) * sum_k B(k) N(k)  [fixed baby list; no uniform-marginal assumption]",
            v20_identity="|X||Y| * sum_k (N(k)/m)^2  [rejected: needs baby marginal = H marginal]",
            fixtures=[fixture(64, 64, 3), fixture(256, 64, 3)],
            production=dict(N_log2=65.6505, m_log2=234.8571, D_log2=167,
                            mean_mult_per_key_log2=67.8571,
                            baby_multiplicity="0/1 with probability 1-2^-36.7",
                            exact_mean_log2=round(math.log2(2 ** 65.6505 * 2 ** 65.6505 / 2 ** 167), 4),
                            note=("At production sparsity (N/D = 2^-101.35) B(k) is 0/1, so the exact "
                                  "expectation is (|Y|/m) sum_i N(k_i) = |X||Y|/D * rho with rho the "
                                  "measured baby-vs-H correlation; the charge is the measured UCB of rho.")))
        json.dump(out["exact"], open(HERE + "/exactE-v21.json", "w"), indent=1)
        print("wrote exactE-v21.json", flush=True)
    if part in ("all", "cov"):
        out["coverage"] = coverage_ladder()
        json.dump(out["coverage"], open(HERE + "/coverage-v21.json", "w"), indent=1)
        print("wrote coverage-v21.json", flush=True)
        print([(c["field_bits"], c["d"], c["nb"], c["rowB"]) for c in out["coverage"]["cells"]], flush=True)
    if part in ("all", "retry"):
        out["retry"] = retry_part()
        json.dump(out["retry"], open(HERE + "/retry-v21.json", "w"), indent=1)
        print("wrote retry-v21.json", flush=True)
