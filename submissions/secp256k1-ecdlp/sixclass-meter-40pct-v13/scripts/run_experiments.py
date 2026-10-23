"""Toy evidence runner for sixclass-meter-40pct-v13.

Modes
  --offline            full captured suite -> evidence/*.json + sweep-summary.json
  --check              rerun a small fixed suite and print the pass lines
  (stdin JSON)         one hosted-style run, printed as JSON

The walk, the class key, the twisted table, the doubling escape and the node
routing are the production constructions at toy scale. Triples are floats: the
production bound is proved in ledger.py with integers/Fractions only.
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


def make(bits_name, seed=1501):
    cur = T.find_curve(BITS[bits_name], seed)
    rnd0 = random.Random(seed + 1)
    base = T.base_point(cur, rnd0)
    best = None
    for lam in (cur.lam, cur.lam * cur.lam % cur.r):
        c2 = T.Curve(cur.p, cur.r, cur.b, cur.beta, lam)
        if all(c2.phi(P) == c2.mul(lam, P) for P in
               [c2.mul(rnd0.randrange(1, cur.r), base) for _ in range(8)]):
            best = c2
            break
    assert best is not None, 'GLV orientation'
    cur = best
    N = (cur.r - 1) // 6
    return cur, base, N


def census(cur, tbl, R, base, rnd, ntrails=16, max_steps=400000,
           dp_bits=0, node_bits=0):
    """Total arrivals until two distinct trails share a class (all arrivals)."""
    logR = R.bit_length() - 1
    st, owner, reports = [], {}, []
    for i in range(ntrails):
        P = cur.mul(rnd.randrange(1, cur.r), base)
        k, Pc = cur.exact_key(P)
        st.append(dict(P=Pc, key=k, ring=[k], seen={k}))
        owner.setdefault(k, i)
    steps = dbl = 0
    while steps < max_steps:
        for i, s in enumerate(st):
            xc = s['key'][0]
            s['P'] = cur.add(s['P'], tbl[xc & (R - 1)][(xc >> (logR + 8)) % 6])
            steps += 1
            nk, s['P'] = cur.exact_key(s['P'])
            if dp_bits and (nk[0] & ((1 << dp_bits) - 1)) == 0:
                reports.append((nk[0] >> (dp_bits + 8)) & ((1 << node_bits) - 1)
                               if node_bits else 0)
            if nk in owner and owner[nk] != i:
                return steps, dbl, reports
            if nk in s['seen']:
                s['P'] = cur.dbl(s['P'])
                steps += 2
                dbl += 1
                nk, s['P'] = cur.exact_key(s['P'])
                s['seen'] = {nk}
                s['ring'] = [nk]
                s['key'] = nk
                owner.setdefault(nk, i)
                continue
            s['seen'].add(nk)
            s['ring'].append(nk)
            if len(s['ring']) > 8:
                s['ring'].pop(0)
            s['key'] = nk
            owner.setdefault(nk, i)
    return max_steps, dbl, reports


def s_ops(N, R):
    npr = N * R / (R - 1)
    return math.sqrt(2 * npr * math.log(1 / 0.58))


def run_census(bits_name, R, trials, seed, ntrails=16):
    cur, base, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    tbl = T.build_table(cur, R, base, rnd)
    S = s_ops(N, R)
    vals, dbl, steps_tot = [], 0, 0
    for t in range(trials):
        v, d, _ = census(cur, tbl, R, base, random.Random(seed * 100003 + t), ntrails)
        vals.append(v)
        dbl += d
        steps_tot += v
    hit = sum(1 for v in vals if v <= S)
    return dict(kind='census', curve=bits_name, R=R, classes=N, order=cur.r,
                trials=trials,
                mean_tc=round(sum(vals) / trials, 3),
                mean_over_sqrtn=round(sum(vals) / trials / math.sqrt(N), 6),
                model_mean_over_sqrtn=round(1.2533 * math.sqrt(R / (R - 1)), 6),
                s_ops=round(S, 3), p_within_s_ops=round(hit / trials, 5),
                model_p=0.420, two_cycle_per_step=round(dbl / steps_tot, 8),
                model_two_cycle=round(1 / (6 * R), 8))


def run_pipeline(bits_name, R, trials, seed, dp_bits=6, node_bits=4):
    """16 storage nodes, per-node cap at twice the mean load, 0 overflows."""
    cur, base, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    tbl = T.build_table(cur, R, base, rnd)
    rep, succ, esc, steps_tot = [], 0, 0, 0
    S = s_ops(N, R)
    for t in range(trials):
        v, d, rr = census(cur, tbl, R, base,
                          random.Random(seed * 7919 + t), 16, 400000,
                          dp_bits=dp_bits, node_bits=node_bits)
        rep.extend(rr)
        esc += d
        steps_tot += v
        succ += 1 if v <= S else 0
    stage = dict(curve=bits_name, R=R, classes=N, trials=trials, seed=seed,
                 dp_bits=dp_bits, node_bits=node_bits, reports=len(rep),
                 recovered=succ, escapes=esc, s_ops=round(S, 3),
                 mean_tc_over_sqrtn=round(steps_tot / trials / math.sqrt(N), 6),
                 p_within_s_ops=round(succ / trials, 5))
    mean = len(rep) / (1 << node_bits)
    cap = max(2, math.ceil(2 * mean))
    load = {}
    over = 0
    for nd in rep:
        load[nd] = load.get(nd, 0) + 1
        if load[nd] > cap:
            over += 1
    stage.update(mean_per_node=round(mean, 2), node_cap=cap, overflows=over,
                 max_load_over_cap=round(max(load.values()) / cap, 4),
                 nodes_used=len(load))
    return stage


def _one_long_trail(cur, tbl, R, base, rnd, length, mode, add_pt):
    """Escapes only; counts returns to a class from an EARLIER escape epoch."""
    logR = R.bit_length() - 1
    P = cur.mul(rnd.randrange(1, cur.r), base)
    key, P = cur.exact_key(P)
    epoch, ever = {key}, {key}
    steps = esc = cross = 0
    while steps < length:
        xc = key[0]
        P = cur.add(P, tbl[xc & (R - 1)][(xc >> (logR + 8)) % 6])
        steps += 1
        nk, P = cur.exact_key(P)
        if nk in epoch:
            if mode == 'doubling':
                P = cur.dbl(P)
                steps += 2
            else:
                P = cur.add(P, add_pt)
                steps += 1
            esc += 1
            nk, P = cur.exact_key(P)
            epoch = {nk}
            ever.add(nk)
            key = nk
            continue
        if nk in ever:
            cross += 1
        epoch.add(nk)
        ever.add(nk)
        key = nk
    return esc, cross


def run_longtrails(bits_name, trails, length, seed, R=4):
    """No DPs, many escapes per trail: doubling vs the rejected additive control."""
    cur, base, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    add_pt = cur.mul(random.Random(seed + 1).randrange(1, cur.r), base)
    out = {}
    for mode in ('doubling', 'additive'):
        tbl = T.build_table(cur, R, base, rnd)
        nf = esc = 0
        for t in range(trails):
            e, c = _one_long_trail(cur, tbl, R, base,
                                   random.Random(seed * 104729 + t), length,
                                   mode, add_pt)
            nf += 1 if c else 0
            esc += e
        out[mode] = dict(trails_with_nonfruitless_returns=nf, escapes=esc,
                         escapes_per_trail=round(esc / trails, 1))
    return dict(kind='long_trails', curve=bits_name, R=R, trails=trails,
                length=length, predicted_self_collision=round(
                    1 - math.exp(-length * length / (2 * N)), 5), **out)


def run_index_rule(bits_name, n_random, n_orbit, seed, k=8, gap=1):
    """0 mismatches for the truncated rule with fallback, plus orbit invariance."""
    cur, base, N = make(bits_name, 1501)
    rnd = random.Random(seed)
    mism = fall = 0
    for _ in range(n_random):
        P = cur.mul(rnd.randrange(1, cur.r), base)
        e, _ = cur.exact_key(P)
        t, _, fb = cur.truncated_key(P, k, gap)
        fall += 1 if fb else 0
        mism += 1 if e != t else 0
    orb_mism = orb_fb = 0
    U = [1, cur.lam, cur.lam * cur.lam % cur.r, cur.r - 1,
         (cur.r - cur.lam) % cur.r, (cur.r - cur.lam * cur.lam) % cur.r]
    for _ in range(n_orbit):
        P = cur.mul(rnd.randrange(1, cur.r), base)
        base_key, _, _ = cur.truncated_key(P, k, gap)
        for u in U:
            Q = cur.mul(u, P)
            e, _ = cur.exact_key(Q)
            t, _, fb = cur.truncated_key(Q, k, gap)
            orb_fb += 1 if fb else 0
            if e != t or t != base_key:
                orb_mism += 1
    return dict(kind='index_rule', curve=bits_name, k_bits=k, gap=gap,
                random_inputs=n_random, random_mismatches=mism,
                random_fallbacks=fall, orbit_checks=n_orbit * 6,
                orbit_fallbacks=orb_fb, orbit_mismatches=orb_mism)


def run_order2(max_m):
    """2^m is never in the unit group U modulo r, for every m <= max_m."""
    r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)
    lam = int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72', 16)
    U = {1, lam, lam * lam % r, r - 1, (r - lam) % r, (r - lam * lam) % r}
    x = 1
    hits = [m for m in range(1, max_m + 1) if (x := x * 2 % r) in U]
    return dict(kind='order2', max_m=max_m, unit_group_size=len(U),
                exponents_in_U=hits, v2_of_r_minus_1=6,
                r_minus_1_small_factors=[2, 2, 2, 2, 2, 2, 3, 149, 631])


def offline():
    ev = HERE / 'evidence'
    ev.mkdir(exist_ok=True)
    cfg = [('q25', 256, 6000), ('q25', 2048, 6000), ('q29', 256, 1200),
           ('q29', 2048, 1200)]
    points = []
    for i, (c, R, tr) in enumerate(cfg):
        row = run_census(c, R, tr, 1501 + i)
        points.append(row)
        (ev / f'census-{i}.json').write_text(json.dumps(row) + '\n')
        print(json.dumps(row), flush=True)
    pipe = run_pipeline('q25', 2048, 3000, 1601)
    (ev / 'pipeline.json').write_text(json.dumps(pipe) + '\n')
    print(json.dumps(pipe), flush=True)
    lt = run_longtrails('q29', 150, 2048, 1701, R=4)
    lt8 = run_longtrails('q29', 150, 2048, 1702, R=8)
    (ev / 'longtrails-4.json').write_text(json.dumps(lt) + '\n')
    (ev / 'longtrails-8.json').write_text(json.dumps(lt8) + '\n')
    print(json.dumps(lt), flush=True)
    print(json.dumps(lt8), flush=True)
    idx = run_index_rule('q25', 2000000, 200000, 1801)
    (ev / 'indexrule.json').write_text(json.dumps(idx) + '\n')
    print(json.dumps(idx), flush=True)
    o2 = run_order2(1 << 22)
    (ev / 'order2.json').write_text(json.dumps(o2) + '\n')
    print(json.dumps(o2), flush=True)
    summary = dict(parameters=dict(R=2048, dp_bits=69, nodes=2 ** 16,
                                   production_R=2 ** 20),
                   points=points, pipeline=pipe, longtrails=[lt, lt8],
                   index_rule=idx, order2=o2)
    (HERE / 'sweep-summary.json').write_text(json.dumps(summary, indent=1))


def check():
    out = []
    r = run_census('q25', 2048, 1200, 2101)
    out.append(('census_q25_2048', r))
    r2 = run_index_rule('q25', 30000, 5000, 2201)
    out.append(('index_rule', r2))
    lt = run_longtrails('q29', 150, 2048, 2301, R=4)
    out.append(('long_trails', lt))
    p = run_pipeline('q25', 2048, 800, 2401)
    out.append(('pipeline16', p))
    print(json.dumps(dict(out)))


if __name__ == '__main__':
    if '--offline' in sys.argv:
        offline()
    elif '--check' in sys.argv:
        check()
    else:
        cfg = json.loads(sys.stdin.read())
        out = dict(census=run_census(cfg.get('curve', 'q25'), cfg.get('R', 2048),
                                     cfg.get('trials', 1200), cfg.get('seed', 1),
                                     cfg.get('parity', 16)))
        if cfg.get('pipeline'):
            out['pipeline'] = run_pipeline(cfg.get('curve', 'q25'),
                                           cfg.get('R', 2048),
                                           max(cfg.get('trials', 1200) // 4, 100),
                                           cfg.get('seed', 1) + 500,
                                           dp_bits=cfg.get('dp_bits', 6))
        print(json.dumps(out))
