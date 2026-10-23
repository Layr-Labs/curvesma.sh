"""Toy evidence runner for sixclass-meter-40pct-v13b.

Modes
  --offline            full captured suite -> evidence/*.json + sweep-summary.json
  --check              rerun a small fixed suite and print the pass lines
  (stdin JSON)         one hosted-style run, printed as JSON

Every number below comes from running the production mechanism: the walks use
the two production detectors (a two-step key compare and a 64-step checkpoint
with the traversal-minimum rule), escapes double the canonical representative of
the minimum-fingerprint class member, the history is cleared at every escape, and
a "recovered" trial is one where a distinguished-point key match produced an
exponent relation whose candidate x was verified as [x]P0 == T.

The census keeps one global class->owner map.  It is instrumentation for the
all-arrival statistic T_c: the walk never reads it, and no escape is decided by
it.  Triples are floats here; the production bound is recomputed in ledger.py and
certificate.json with integers and Fractions only.
"""
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / 'scripts'))
import toy_core as T

BITS = {'q25': 25, 'q29': 29}
MASK64 = (1 << 64) - 1


def make(bits_name, seed=1501):
    cur = T.find_curve(BITS[bits_name], seed)
    rnd = random.Random(seed + 1)
    base = T.base_point(cur, rnd)
    assert cur.is_on(base) and T.check_glv(cur, rnd, 20)
    return cur, base, (cur.r - 1) // 6


def s_ops(N, R, model=0.42):
    return math.sqrt(2 * N * R / (R - 1) * math.log(1 / (1 - model)))


# ------------------------------------------------------------------- census
def census_trial(cur, tbl, R, P0, rnd, ntrails, Tt, dp_bits=0, node_cap=None,
                 node_bits=0, stop_on_collision=True, max_steps=10 ** 6):
    """Production walks; T_c is the first arrival at another trail's class."""
    trails, owner = [], {}
    tot = 0
    for i in range(ntrails):
        a = rnd.randrange(1, cur.r)
        b = rnd.randrange(1, cur.r)
        trails.append(T.start_trail(cur, P0, Tt, a, b))
    def arrive(i, key):
        o = owner.setdefault(key, i)
        return o != i
    for i, tr in enumerate(trails):
        arrive(i, tr.K)
    t_c = None
    reports = []
    aborts = 0
    while tot < max_steps:
        for i, tr in enumerate(trails):
            if tr.step(cur, tbl, R) == 2:
                a = rnd.randrange(1, cur.r)
                b = rnd.randrange(1, cur.r)
                trails[i] = T.start_trail(cur, P0, Tt, a, b)
                aborts += 1
                tot += 1
                continue
            tot += 1
            if t_c is None and arrive(i, tr.K):
                t_c = tot
                if stop_on_collision:
                    return t_c, tot, tr, trails, reports, aborts
            if dp_bits and (tr.K & ((1 << dp_bits) - 1)) == 0:
                reports.append((tr.K >> (dp_bits + 4)) & ((1 << node_bits) - 1)
                               if node_bits else 0)
    return t_c, tot, None, trails, reports, aborts


def run_census(bits_name, R, trials, seed, ntrails=16):
    cur, P0, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    tbl = T.build_table(cur, R, P0, rnd)
    S = s_ops(N, R)
    vals, e2, ec, steps_tot, aborts = [], 0, 0, 0, 0
    for t in range(trials):
        r = random.Random(seed * 100003 + t)
        Tt = cur.mul(r.randrange(1, cur.r), P0)
        tc, tot, _, trails, _, ab = census_trial(cur, tbl, R, P0, r, ntrails, Tt)
        aborts += ab
        vals.append(tc if tc is not None else tot)
        e2 += sum(x.esc2 for x in trails)
        ec += sum(x.escchk for x in trails)
        steps_tot += vals[-1]
    hit = sum(1 for v in vals if v <= S)
    return dict(kind='census', curve=bits_name, R=R, classes=N, order=cur.r,
                trials=trials,
                mean_tc=round(sum(vals) / trials, 3),
                mean_over_sqrtn=round(sum(vals) / trials / math.sqrt(N), 6),
                model_mean_over_sqrtn=round(1.2533 * math.sqrt(R / (R - 1)), 6),
                s_ops=round(S, 3), p_within_s_ops=round(hit / trials, 5),
                model_p=0.420, escapes_2cycle=e2, escapes_checkpoint=ec,
                two_cycle_per_step=round(e2 / steps_tot, 8),
                exceptional_additions=aborts,
                model_two_cycle=round(1 / (6 * R), 8),
                frontier_two_cycle=round(1 / (6 * R) + 2 ** -20, 8))


# ----------------------------------------------------------------- pipeline
def solve_pair(cur, P0, Tt, r1, r2):
    """x = (u beta2 - beta1)/(alpha1 - u alpha2); u=1 is the stated formula."""
    a1, b1 = r1[0], r1[1]
    a2, b2 = r2[0], r2[1]
    for u in cur.U:
        g = (a1 - u * a2) % cur.r
        if g == 0:
            continue
        x = (u * b2 - b1) * pow(g, cur.r - 2, cur.r) % cur.r
        if cur.mul(x, P0) == Tt:
            return x, u
    return None, None


def pipeline_trial(cur, tbl, R, P0, rnd, ntrails=16, dp_bits=4, cap_mult=8.0,
                   node_bits=0, node_cap=None):
    N = (cur.r - 1) // 6
    xt = rnd.randrange(1, cur.r)
    Tt = cur.mul(xt, P0)
    tab, trails, owner, tid = {}, [], {}, []
    tot = 0
    e2 = ec = 0
    maxlen = 40 * (1 << dp_bits)
    t_c = None
    matches = deg = 0
    aborts = 0
    units = {}
    serial = 0
    for i in range(ntrails):
        a = rnd.randrange(1, cur.r)
        b = rnd.randrange(1, cur.r)
        trails.append(T.start_trail(cur, P0, Tt, a, b))
        tid.append(serial)
        serial += 1
    while tot < cap_mult * s_ops(N, R):
        for i, tr in enumerate(trails):
            if tr.step(cur, tbl, R) == 2:
                a = rnd.randrange(1, cur.r)
                b = rnd.randrange(1, cur.r)
                e2 += tr.esc2
                ec += tr.escchk
                trails[i] = T.start_trail(cur, P0, Tt, a, b)
                tid[i] = serial
                serial += 1
                aborts += 1
                tot += 1
                continue
            tot += 1
            o = owner.setdefault(tr.K, tid[i])
            if t_c is None and o != tid[i]:
                t_c = tot
            if (tr.K & ((1 << dp_bits) - 1)) == 0:
                rec = (tr.alpha, tr.beta, tot)
                lst = tab.setdefault(tr.K, [])
                for old in lst:
                    matches += 1
                    x, u = solve_pair(cur, P0, Tt, rec, old)
                    if x is None:
                        deg += 1
                        continue
                    units[u] = units.get(u, 0) + 1
                    return dict(ok=True, xok=(x == xt), steps=tot, t_c=t_c,
                                lag=(tot - t_c) if t_c else None,
                                esc2=e2 + sum(t.esc2 for t in trails),
                                escch=ec + sum(t.escchk for t in trails),
                                matches=matches, deg=deg, units=units,
                                aborts=aborts,
                                moved=sum(1 for t in trails if t.esc2 + t.escchk))
                lst.append(rec)
            if tr.steps >= maxlen:
                a = rnd.randrange(1, cur.r)
                b = rnd.randrange(1, cur.r)
                e2 += tr.esc2
                ec += tr.escchk
                trails[i] = T.start_trail(cur, P0, Tt, a, b)
                tid[i] = serial
                serial += 1
    return dict(ok=False, steps=tot, t_c=t_c, lag=None,
                esc2=e2 + sum(t.esc2 for t in trails),
                escch=ec + sum(t.escchk for t in trails),
                matches=matches, deg=deg, units=units, aborts=aborts,
                moved=sum(1 for t in trails if t.esc2 + t.escchk))


def run_pipeline(bits_name, R, trials, seed, dp_bits=4, cap_mult=8.0):
    cur, P0, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    tbl = T.build_table(cur, R, P0, rnd)
    S = s_ops(N, R)
    ok = 0
    coll = 0
    steps, lag, esc2, escch, matches, deg, moved = [], [], 0, 0, 0, 0, 0
    units = {}
    for t in range(trials):
        r = pipeline_trial(cur, tbl, R, P0, random.Random(seed * 7919 + t),
                           dp_bits=dp_bits, cap_mult=cap_mult)
        ok += 1 if r['ok'] else 0
        coll += 1 if r['t_c'] is not None and r['t_c'] <= S else 0
        steps.append(r['steps'])
        esc2 += r['esc2']
        escch += r['escch']
        matches += r['matches']
        deg += r['deg']
        moved += r['moved']
        if r['ok'] and r['lag'] is not None:
            lag.append(r['lag'])
        for u, n in r['units'].items():
            units[str(u)] = units.get(str(u), 0) + n
    hit = sum(1 for s in steps if s <= S)
    lag.sort()
    return dict(kind='pipeline', curve=bits_name, R=R, classes=N, trials=trials,
                seed=seed, dp_bits=dp_bits, s_ops=round(S, 3),
                collisions_within_s_ops=coll,
                recoveries_verified=ok, p_recovered_within_s_ops=round(hit / trials, 5),
                p_recovered_overall=round(ok / trials, 5),
                mean_steps_over_sqrtn=round(sum(steps) / trials / math.sqrt(N), 6),
                mean_lag_total_steps=round(sum(lag) / len(lag), 3) if lag else None,
                mean_lag_per_trail=round(sum(lag) / len(lag) / 16, 3) if lag else None,
                max_lag_total_steps=max(lag) if lag else None,
                lag_units=2 ** dp_bits,
                escapes_2cycle=esc2, escapes_checkpoint=escch,
                trails_that_escaped=moved, matches=matches, degenerate=deg,
                unit_histogram=units, unit_one_fraction=round(units.get('1', 0) / max(matches - deg, 1), 4))


# ------------------------------------------------------- 2^16 node routing
MIX = 0x9E3779B97F4A7C15
MIX2 = 0xBF58476D1CE4E5B9


def routing_word(m):
    """The toy's routing word: a fixed 64-bit mix of the class key.

    Production routes on bits 128..143 of the 192-bit fingerprint; a 25-bit toy
    key has no such bits, so the same role is played by a deterministic mix of
    the class key.  It is a function of the class alone, which is what merging
    and the uniformity clause need.
    """
    w = (m * MIX) & MASK64
    w ^= w >> 29
    w = (w * MIX2) & MASK64
    w ^= w >> 32
    return w


def run_routing16(bits_name, samples, seed, node_bits=16):
    """Uniformity and capacity of the routing word over 2^16 nodes."""
    cur, P0, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    p, beta = cur.p, cur.beta
    load = [0] * (1 << node_bits)
    for _ in range(samples):
        x = rnd.randrange(p)
        x2 = beta * x % p
        x3 = (-x - x2) % p
        m = x if (x <= x2 and x <= x3) else (x2 if (x2 <= x3) else x3)
        load[(routing_word(m) >> (64 - node_bits)) & ((1 << node_bits) - 1)] += 1
    mean = samples / (1 << node_bits)
    mx = max(load)
    # multiplicative Chernoff, union over 2^16 nodes, bound 2^-10
    delta = 1.0
    while math.exp(-delta * delta * mean / 3) > 2 ** -26:
        delta *= 1.05
    cap = math.ceil(mean * (1 + delta))
    over = sum(1 for v in load if v > cap)
    used = sum(1 for v in load if v)
    chi = sum((v - mean) ** 2 / mean for v in load)
    return dict(kind='routing16', curve=bits_name, nodes=1 << node_bits,
                samples=samples, mean_per_node=round(mean, 4), max_load=mx,
                max_over_mean=round(mx / mean, 4), cap=cap,
                cap_over_mean=round(cap / mean, 4), overflows=over,
                nodes_used=used, chi_square=round(chi, 2),
                chi_square_dof=(1 << node_bits) - 1,
                chernoff_delta=round(delta, 4))


# ------------------------------------------------------------- long trails
def run_traps(bits_name, trails, length, seed, R=4):
    """Long trails with many doubling escapes and no distinguished points.

    A trail is counted when it arrives at a class it visited in an *earlier*
    escape epoch: that is the non-fruitless return a composite cycle would
    produce.  The detection is the production detection, so the escape rate is
    the production rate, not the seen-set rate.
    """
    cur, P0, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    Tt = cur.mul(random.Random(seed + 2).randrange(1, cur.r), P0)
    tbl = T.build_table(cur, R, P0, rnd)
    nf = esc = 0
    for t in range(trails):
        r = random.Random(seed * 104729 + t)
        tr = T.start_trail(cur, P0, Tt, r.randrange(1, cur.r), r.randrange(1, cur.r))
        epoch, ever, cross = {tr.K}, {tr.K}, 0
        while tr.steps < length:
            before = tr.esc2 + tr.escchk
            tr.step(cur, tbl, R)
            if (tr.esc2 + tr.escchk) > before:
                epoch = {tr.K}              # a new epoch starts at the escape
                ever.add(tr.K)
            else:
                if tr.K in ever and tr.K not in epoch:
                    cross += 1
                epoch.add(tr.K)
                ever.add(tr.K)
        nf += 1 if cross else 0
        esc += tr.esc2 + tr.escchk
    return dict(kind='traps', curve=bits_name, R=R, trails=trails,
                length=length, escapes=esc,
                escapes_per_trail=round(esc / trails, 1),
                trails_with_cross_epoch_returns=nf,
                predicted_self_collision=round(
                    1 - math.exp(-length * length / (2 * N)), 5))


# ---------------------------------------------------------------- merge test
def run_merge(bits_name, trials, steps, seed, R=8):
    """Six unit partners of one point must walk one class sequence."""
    cur, P0, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    tbl = T.build_table(cur, R, P0, rnd)
    bad = 0
    seq_checks = 0
    esc = 0
    for t in range(trials):
        Q = cur.mul(rnd.randrange(1, cur.r), P0)
        seqs = []
        for u in cur.U:
            tr = T.Trail(cur, cur.mul(u, Q), 0, 1)
            s = []
            for _ in range(steps):
                tr.step(cur, tbl, R)
                s.append(tr.K)
            seqs.append(s)
            esc += tr.esc2 + tr.escchk
        if any(s != seqs[0] for s in seqs[1:]):
            bad += 1
        seq_checks += steps * 5
    return dict(kind='merge', curve=bits_name, R=R, trials=trials,
                steps_per_trail=steps, class_sequence_checks=seq_checks,
                merge_failures=bad, escapes=esc)


# --------------------------------------------------------- checkpoint test
def run_checkpoint(bits_name, trails, length, seed, R=4):
    """Small R: the checkpoint window fires often; report both detectors."""
    cur, P0, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    tbl = T.build_table(cur, R, P0, rnd)
    Tt = cur.mul(rnd.randrange(1, cur.r), P0)
    e2 = ec = steps = 0
    per_trail = []
    for t in range(trails):
        r = random.Random(seed * 65537 + t)
        tr = T.start_trail(cur, P0, Tt, r.randrange(1, cur.r), r.randrange(1, cur.r))
        while tr.steps < length:
            tr.step(cur, tbl, R)
        e2 += tr.esc2
        ec += tr.escchk
        steps += tr.steps
        per_trail.append(tr.esc2 + tr.escchk)
    return dict(kind='checkpoint', curve=bits_name, R=R, trails=trails,
                length=length, steps=steps, escapes_2cycle=e2,
                escapes_checkpoint=ec,
                checkpoint_share=round(ec / max(e2 + ec, 1), 4),
                escapes_per_trail=round((e2 + ec) / trails, 2),
                min_escapes_per_trail=min(per_trail))


# ---------------------------------------------------------------- index rule
def run_index_rule(n_random, n_adv, n_orbit, seed, gap=40, guard=8):
    """The production 8-product rule on the production prime and beta."""
    p = int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f', 16)
    beta = int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee', 16)
    rnd = random.Random(seed)
    err_max = [0, 0, 0]
    mism_no_fb = mism_fb = fb_random = 0
    cheap_wrong_random = 0
    worst = 0
    wraps = 0
    for _ in range(n_random):
        x = rnd.randrange(1, p)
        a, fb, comp, u, tie, wrap = T.truncated_index(p, beta, x, gap, guard)
        ex = T._exact_argmin(u)
        for i in range(3):
            d = abs(comp[i] - u[i])
            if d > (1 << 63):
                wraps += 1
            e = min(d, (1 << 64) - d)          # cyclic error in ulps
            if e > err_max[i]:
                err_max[i] = e
            if e > worst:
                worst = e
        if a != ex:
            if fb:
                mism_fb += 1
            else:
                mism_no_fb += 1
        if T._exact_argmin(comp) != ex:
            cheap_wrong_random += 1
        fb_random += 1 if fb else 0
    # adversarial near-ties: x*(beta-1) = delta, x*(beta^2-1) = delta, beta*x*(beta-1) = delta
    inv1 = pow((beta - 1) % p, p - 2, p)
    inv2 = pow((beta * beta - 1) % p, p - 2, p)
    inv3 = pow(beta * (beta - 1) % p, p - 2, p)
    adv_n = adv_fb = adv_mism_no_fb = adv_mism_fb = 0
    adv_cheap_wrong = 0
    adv_tie_checked = 0
    hist = {}
    for _ in range(n_adv):
        d = rnd.randrange(0, 401)
        dx = d << 192                      # d ulps of the 64-bit fraction
        which = rnd.randrange(3)
        mul = (inv1, inv2, inv3)[which]
        x = dx * mul % p
        if x == 0:
            continue
        a, fb, comp, u, tie, wrap = T.truncated_index(p, beta, x, gap, guard)
        ex = T._exact_argmin(u)
        su = sorted(u)
        if abs((su[1] - su[0]) - d) > 1:
            continue
        adv_n += 1
        if d <= gap:
            adv_tie_checked += 1
        hist[d // 20] = hist.get(d // 20, 0) + 1
        if fb:
            adv_fb += 1
        if T._exact_argmin(comp) != ex:
            adv_cheap_wrong += 1
        if a != ex:
            if fb:
                adv_mism_fb += 1
            else:
                adv_mism_no_fb += 1
    # orbit invariance: all three orbit entries must select one canonical x
    orb_mism = orb_fb = 0
    for _ in range(n_orbit):
        x = rnd.randrange(1, p)
        vals = [x, beta * x % p, (-x - beta * x) % p]
        chosen = []
        for xi in vals:
            a2, fb2, comp2, u2, _t, _w = T.truncated_index(p, beta, xi, gap, guard)
            orb_fb += 1 if fb2 else 0
            orbit_i = [xi, beta * xi % p, (-xi - beta * xi) % p]
            idx = T._exact_argmin(u2) if fb2 else a2
            chosen.append(orbit_i[idx])
        if len(set(chosen)) != 1:
            orb_mism += 1
    return dict(kind='index_rule', prime_bits=p.bit_length(), gap=gap, guard=guard,
                random_inputs=n_random, random_mismatches_no_fallback=mism_no_fb,
                random_mismatches_with_fallback=mism_fb,
                random_fallbacks=fb_random,
                fallback_rate_measured=round(fb_random / n_random, 8),
                max_abs_error_ulps=err_max, worst_error_ulps=worst,
                gross_wraps=wraps,
                adversarial_inputs=adv_n, adversarial_within_gap=adv_tie_checked,
                adversarial_fallbacks=adv_fb,
                adversarial_mismatches_no_fallback=adv_mism_no_fb,
                adversarial_mismatches_with_fallback=adv_mism_fb,
                adversarial_cheap_rule_would_be_wrong=adv_cheap_wrong,
                random_cheap_rule_would_be_wrong=cheap_wrong_random,
                orbit_checks=n_orbit * 3, orbit_fallbacks=orb_fb,
                orbit_mismatches=orb_mism)


def run_order2(max_m):
    """2^m is never in the unit group U modulo r, for every m <= max_m."""
    r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)
    lam = int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72', 16)
    U = {1, lam * lam % r, lam, r - 1, (r - lam) % r, (r - lam * lam) % r}
    x = 1
    hits = [m for m in range(1, max_m + 1) if (x := x * 2 % r) in U]
    return dict(kind='order2', max_m=max_m, unit_group_size=len(U),
                exponents_in_U=hits, v2_of_r_minus_1=6,
                r_minus_1_small_factors=[2, 2, 2, 2, 2, 2, 3, 149, 631])


def offline():
    ev = HERE / 'evidence'
    ev.mkdir(exist_ok=True)
    cfg = [('q25', 256, 2400), ('q25', 2048, 2400), ('q29', 256, 1000),
           ('q29', 2048, 1000)]
    points = []
    for i, (c, R, tr) in enumerate(cfg):
        row = run_census(c, R, tr, 1501 + i)
        points.append(row)
        (ev / f'census-{i}.json').write_text(json.dumps(row) + '\n')
        print(json.dumps(row), flush=True)
    pipe = run_pipeline('q25', 2048, 4000, 1601, dp_bits=3)
    (ev / 'pipeline.json').write_text(json.dumps(pipe) + '\n')
    print(json.dumps(pipe), flush=True)
    pipe1 = run_pipeline('q25', 2048, 4000, 1601, dp_bits=1)
    (ev / 'pipeline-dp1.json').write_text(json.dumps(pipe1) + '\n')
    print(json.dumps(pipe1), flush=True)
    rt = run_routing16('q25', 1 << 22, 1651)
    (ev / 'routing16.json').write_text(json.dumps(rt) + '\n')
    print(json.dumps(rt), flush=True)
    lt = run_traps('q29', 150, 2048, 1701, R=4)
    lt8 = run_traps('q29', 150, 2048, 1702, R=8)
    (ev / 'traps-4.json').write_text(json.dumps(lt) + '\n')
    (ev / 'traps-8.json').write_text(json.dumps(lt8) + '\n')
    print(json.dumps(lt), flush=True)
    print(json.dumps(lt8), flush=True)
    mg = run_merge('q25', 150, 200, 1751, R=8)
    (ev / 'merge.json').write_text(json.dumps(mg) + '\n')
    print(json.dumps(mg), flush=True)
    ck = run_checkpoint('q25', 200, 600, 1761, R=4)
    (ev / 'checkpoint.json').write_text(json.dumps(ck) + '\n')
    print(json.dumps(ck), flush=True)
    idx = run_index_rule(1500000, 600000, 200000, 1801)
    (ev / 'indexrule.json').write_text(json.dumps(idx) + '\n')
    print(json.dumps(idx), flush=True)
    o2 = run_order2(1 << 22)
    (ev / 'order2.json').write_text(json.dumps(o2) + '\n')
    print(json.dumps(o2), flush=True)
    summary = dict(parameters=dict(R=2048, dp_bits=69, nodes=2 ** 16,
                                   production_R=2 ** 20, checkpoint_period=64,
                                   dp_bits_toy_pipeline=3),
                   points=points, pipeline=pipe, pipeline_dp1=pipe1, routing16=rt,
                   traps=[lt, lt8], merge=mg, checkpoint=ck,
                   index_rule=idx, order2=o2)
    (HERE / 'sweep-summary.json').write_text(json.dumps(summary, indent=1))


def check():
    out = []
    out.append(('census', run_census('q25', 2048, 1200, 2101)))
    out.append(('pipeline', run_pipeline('q25', 2048, 800, 2401, dp_bits=3)))
    out.append(('index', run_index_rule(60000, 30000, 8000, 2201)))
    out.append(('traps', run_traps('q29', 40, 1024, 2301, R=4)))
    out.append(('merge', run_merge('q25', 40, 150, 2351, R=8)))
    out.append(('checkpoint', run_checkpoint('q25', 60, 400, 2361, R=4)))
    out.append(('routing16', run_routing16('q25', 1 << 18, 2371)))
    print(json.dumps(dict(out)))


if __name__ == '__main__':
    if '--offline' in sys.argv:
        offline()
    elif '--check' in sys.argv:
        check()
    else:
        cfg = json.loads(sys.stdin.read())
        c = cfg.get('curve', 'q25')
        R = cfg.get('R', 2048)
        tr = cfg.get('trials', 1200)
        sd = cfg.get('seed', 1)
        rows = [run_census(c, R, max(tr // 4, 1), sd * 13 + j) for j in range(4)]
        hit = sum(round(x['p_within_s_ops'] * x['trials']) for x in rows)
        n = sum(x['trials'] for x in rows)
        out = dict(census=rows[0],
                   census_pooled=dict(seeds=[sd * 13 + j for j in range(4)],
                                      trials=n,
                                      p_within_s_ops=round(hit / n, 5),
                                      mean_over_sqrtn=round(
                                          sum(x['mean_tc'] * x['trials'] for x in rows)
                                          / n / math.sqrt(rows[0]['classes']), 6)))
        if cfg.get('pipeline'):
            out['pipeline'] = run_pipeline(c, cfg.get('R', 2048),
                                           max(cfg.get('trials', 1200) // 4, 120),
                                           cfg.get('seed', 1) + 500,
                                           dp_bits=cfg.get('dp_bits', 3))
        print(json.dumps(out))
