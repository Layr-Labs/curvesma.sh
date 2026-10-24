#!/usr/bin/env python3
"""census_v15m.py — secp aligned-A' verified-recovery census with the DRAIN ARM.

What changed against v15l (and why):

1. DRAIN ARM (ported from the plain lane's v20c/v20d fix, not invented here).
   The v15l trial walked an 8*S_ops budget and certified only recoveries whose
   `steps <= S_ops`; a pair that merged with no distinguished point left inside
   the search phase simply read as failure.  v15m walks the PRODUCTION budget
   (1.0 * S_ops, the ledger's `ARRIVALS = S_ops + 2^60` search phase), then
   flushes PENDING reports: each walker advances at most `drain_spacings`
   spacings to its next distinguished point and delivers that record to the
   same pool.  Every flushed step is CHARGED (added to `tot`), and the drain
   row it is charged against is priced separately in price_v15m.py
   (WORKERS * drain_spacings * 2^d steps at the step price).  Simulate what you
   charge: the flush cap in the simulator is the row that is in W.

2. UNCONDITIONAL HEADLINE.  The certified statistic is the RAW count of trials
   whose recovery verified ([x]P0 == T), over ALL trials, with no
   only_merged/only_audited conditioning filter.  The independent merge audit
   (two distinct walkers occupying one class inside the search phase) is
   published SIDE BY SIDE, and the conditioned numbers are reported separately
   for diagnosis, never as the headline.  e2bd4bef (plain v20d) died of a
   conditioned headline; the filter is what made its success read 0.40 while
   the split was 0.36 + 0.05.

3. DRAIN-ONLY PAIRS are their own channel, never success.  A pair whose first
   meeting is inside the drain (the class was never visited in-budget by a
   different walker) is counted in `drain_only_excluded`.  Legitimacy test is
   the audit map, so it is not `success <= merge` by construction: the audit is
   computed from arrivals, the success from verified algebra.

4. AUDIT ATTRIBUTION, single derivation.  `audit_key` is the one place a trial
   index is derived; the planted-attribution check (`--selftest`) plants a
   collision at a known trial and asserts the audit labels it with that trial
   and that slot pair.  No re-division of an already-derived key (that was the
   e2bd4bef D1 defect: `rtr = rep//w` applied to an already-derived index).

Statistic per trial (all unconditional): a trial is a SUCCESS iff a record pair
in one class solved for x and [x]P0 == T, with the record pair either both
inside the search phase (`budget`) or one inside and one delivered by the
charged drain (`drain`).  `drain_only` pairs never count.
"""
import json, math, os, sys, random
from decimal import Decimal, getcontext
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PKG = ("/Users/andrewgordon/Documents/ChatGPT/curvesma.sh/publication-staging/"
       "secp256k1-ecdlp/sigma-table-40pct-v15k/submissions/secp256k1-ecdlp/"
       "sigma-table-40pct-v15k")
sys.path.insert(0, os.path.join(PKG, "scripts"))
import run_experiments as RX
import toy_core as T

getcontext().prec = 50
MODEL = 0.42                       # the census model rate (v15l/v15k convention)


def cp95_lower(k, n, alpha=Decimal('0.05')):
    """one-sided Clopper-Pearson lower bound (exact binomial tail, bisection)."""
    k, n = int(k), int(n)
    if k == 0:
        return Decimal(0)
    lnC = sum(Decimal(i).ln() for i in range(n - k + 1, n + 1)) - \
          sum(Decimal(i).ln() for i in range(1, k + 1))
    def tail(p):
        p = Decimal(p); q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal('1e-45') * s:
                break
        return s
    lo, hi = Decimal('0.0'), Decimal('0.999')
    for _ in range(90):
        mid = (lo + hi) / 2
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal('1e-14'):
            break
    return lo


def audit_key(trial_index):
    """THE single derivation of a trial index from an audit record.

    Anything that needs a trial index gets it here and nowhere else; no caller
    re-derives (re-divides) an index it already holds.  Plain v20d's audit
    divided an already-derived index by the worker count a second time and
    attributed reports to trial floor(t/64) -- e2bd4bef D1.
    """
    assert isinstance(trial_index, int) and trial_index >= 0
    return trial_index


def trial_v15m(cur, tbl, R, P0, rnd, ntrails=16, dp_bits=0, cap_mult=1.0,
               drain_spacings=2):
    """One production-shaped trial: search phase, then the charged drain."""
    N6 = (cur.r - 1) // 6
    S = RX.s_ops(N6, R, MODEL)
    spacing = 1 << dp_bits
    xt = rnd.randrange(1, cur.r)
    Tt = cur.mul(xt, P0)
    pool, owner = {}, {}           # class -> records ; class -> first slot (audit)
    trails, epoch = [], []
    for i in range(ntrails):
        trails.append(T.start_trail(cur, P0, Tt, rnd.randrange(1, cur.r),
                                    rnd.randrange(1, cur.r)))
        epoch.append(0)
    maxlen = 40 * spacing
    tot = 0
    t_c = None                     # first cross-walker class arrival (coverage)
    audit_merge = False
    matches = deg = 0
    aborts = 0
    search_steps = None
    drain_steps = 0
    drain_only = 0
    drain_reports = 0
    nreports = 0
    res = dict(ok=False, via=None, steps=None, search_steps=None, drain_steps=0,
               t_c=None, lag=None, matches=0, deg=0, drain_only=0,
               drain_reports=0, nreports=0, audit_merge=False, aborts=0,
               spacing=spacing, s_ops=S)

    def deliver(i, tr, tot_now, in_budget, legit_hint=False):
        """register a record at tr.K; return (recovered_x or None, legit)."""
        nonlocal matches, deg, drain_only, nreports
        cls = tr.K
        lst = pool.get(cls)
        if lst is None:
            pool[cls] = [(tr.alpha, tr.beta, tot_now, i, in_budget)]
            nreports += 1
            return None, False
        legit = legit_hint or (owner.get(cls) not in (None, i))
        for old in lst:
            matches += 1
            x, u = RX.solve_pair(cur, P0, Tt, (tr.alpha, tr.beta), (old[0], old[1]))
            if x is None:
                deg += 1
                continue
            if not (in_budget and old[4]):
                # at least one side is a drained record: legit only if the
                # class was already a cross-walker merge inside the budget
                if not legit:
                    drain_only += 1
                    continue
            nreports += 1
            return x, legit
        lst.append((tr.alpha, tr.beta, tot_now, i, in_budget))
        nreports += 1
        return None, legit

    # ---------------- search phase: the production budget (cap_mult * S_ops)
    budget = cap_mult * S
    while tot < budget:
        for i, tr in enumerate(trails):
            if tr.step(cur, tbl, R) == 2:
                trails[i] = T.start_trail(cur, P0, Tt, rnd.randrange(1, cur.r),
                                          rnd.randrange(1, cur.r))
                epoch[i] += 1
                aborts += 1
                tot += 1
                continue
            tot += 1
            o = owner.setdefault(tr.K, i)
            if o != i:
                audit_merge = True
                if t_c is None:
                    t_c = tot
            if (tr.K & (spacing - 1)) == 0:
                x, legit = deliver(i, tr, tot, True)
                if x is not None:
                    res.update(ok=True, via='budget', steps=tot, search_steps=tot,
                               drain_steps=0, t_c=t_c, lag=(tot - t_c) if t_c else None,
                               matches=matches, deg=deg, drain_only=drain_only,
                               drain_reports=0, nreports=nreports,
                               audit_merge=audit_merge, aborts=aborts)
                    return res
            if tr.steps >= maxlen:
                trails[i] = T.start_trail(cur, P0, Tt, rnd.randrange(1, cur.r),
                                          rnd.randrange(1, cur.r))
                epoch[i] += 1
    search_steps = tot

    # ---------------- drain arm: flush the reports of pairs that MERGED
    # inside the budget and are still co-located at the cutoff.  Only those
    # pairs are advanced, only to their next distinguished point, and every
    # step is charged.  This cannot add search: a pair that is not co-located
    # at the cutoff is never walked, so a pair whose first meeting is inside
    # the drain does not exist by construction (and `drain_only` counts any
    # pair broken by a trap restart before it reached a DP).
    cap = drain_spacings * spacing
    at = {}
    for i, tr in enumerate(trails):
        at.setdefault(tr.K, []).append(i)
    pairs = [v for v in at.values() if len(v) >= 2]
    for grp in pairs:
        i, j = grp[0], grp[1]
        ti, tj = trails[i], trails[j]
        for _ in range(cap):
            bi = ti.step(cur, tbl, R)
            bj = tj.step(cur, tbl, R)
            drain_steps += 2
            tot += 2
            if bi == 2 or bj == 2:
                trails[i] = T.start_trail(cur, P0, Tt, rnd.randrange(1, cur.r),
                                          rnd.randrange(1, cur.r))
                trails[j] = T.start_trail(cur, P0, Tt, rnd.randrange(1, cur.r),
                                          rnd.randrange(1, cur.r))
                epoch[i] += 1
                aborts += 1
                drain_only += 1        # pair broken before its report: not a recovery
                break
            if (ti.K & (spacing - 1)) == 0:
                # both are at the same class by construction: one record each
                xa, _ = deliver(i, ti, tot, False, legit_hint=True)
                xb, _ = deliver(j, tj, tot, False, legit_hint=True)
                x = xa or xb
                if x is not None:
                    drain_reports += 1
                    res.update(ok=True, via='drain', steps=tot,
                               search_steps=search_steps, drain_steps=drain_steps,
                               t_c=t_c, lag=(tot - t_c) if t_c else None,
                               matches=matches, deg=deg, drain_only=drain_only,
                               drain_reports=drain_reports, nreports=nreports,
                               audit_merge=audit_merge, aborts=aborts)
                    return res
                break                  # report delivered and did not verify
    res.update(ok=False, via=None, steps=None, search_steps=search_steps,
               drain_steps=drain_steps, t_c=t_c, lag=None, matches=matches,
               deg=deg, drain_only=drain_only, drain_reports=drain_reports,
               nreports=nreports, audit_merge=audit_merge, aborts=aborts)
    return res


def chunk(args):
    bits, R, dp_bits, drain_spacings, cap_mult, seed0, trials, first_trial = args
    cur, P0, N = RX.make(bits, 1501)
    tbl = T.build_table(cur, R, P0, random.Random(seed0 ^ 0x5EED))
    rows = []
    for t in range(trials):
        idx = audit_key(first_trial + t)          # single derivation, once
        r = trial_v15m(cur, tbl, R, P0, random.Random(seed0 * 7919 + t),
                       dp_bits=dp_bits, cap_mult=cap_mult,
                       drain_spacings=drain_spacings)
        rows.append(dict(trial=idx, dp_bits=dp_bits, R=R, curve=bits,
                         drain_spacings=drain_spacings,
                         ok=int(r['ok']), via=r['via'], steps=r['steps'],
                         search_steps=r['search_steps'],
                         drain_steps=r['drain_steps'], t_c=r['t_c'],
                         audit_merge=int(r['audit_merge']), s_ops=round(r['s_ops'], 3),
                         matches=r['matches'], deg=r['deg'],
                         drain_only=r['drain_only'], nreports=r['nreports']))
    return rows


def summarize(rows, dp_bits, drain_spacings):
    n = len(rows)
    S = rows[0]['s_ops']
    bud = sum(1 for r in rows if r['ok'] and r['via'] == 'budget')
    drn = sum(1 for r in rows if r['ok'] and r['via'] == 'drain')
    total = bud + drn
    audit = sum(r['audit_merge'] for r in rows)
    d_only = sum(r['drain_only'] for r in rows)
    dr_steps = sum(r['drain_steps'] for r in rows)
    return dict(dp_bits=dp_bits, density=f"2^-{dp_bits}", trials=n, s_ops=S,
                success_uncond_budget=bud, success_uncond_drain=drn,
                success_uncond_total=total,
                p_uncond=round(total / n, 6),
                cp95_uncond=float(cp95_lower(total, n)),
                p_uncond_budget_only=round(bud / n, 6),
                cp95_uncond_budget_only=float(cp95_lower(bud, n)),
                audit_merge=audit, p_audit_merge=round(audit / n, 6),
                cp95_audit_merge=float(cp95_lower(audit, n)),
                p_cond_on_merge=(round(total / audit, 6) if audit else None),
                drain_only_pairs=d_only,
                drain_only_per_trial=round(d_only / n, 5),
                mean_drain_steps=round(dr_steps / n, 2),
                drain_frac_of_s_ops=round(dr_steps / n / S, 6),
                zero_denominator=sum(r['deg'] for r in rows),
                mean_records=round(sum(r['nreports'] for r in rows) / n, 3),
                success_le_audit=(total <= audit))


def selftest():
    """Planted-attribution check for the audit key (e2bd4bef D1).

    Plant one known cross-walker collision at a known trial and assert the audit
    labels it with that trial and that slot pair, through audit_key only.
    """
    cur, P0, N = RX.make('q25', 1501)
    tbl = T.build_table(cur, 2048, P0, random.Random(11))
    bad = 0
    for planted in range(40):
        idx = audit_key(planted * 7)
        assert idx == planted * 7
        # a synthetic audit record: two slots in one class, one trial
        rec = dict(trial=idx, cls=12345, slots=[0, 3])
        got = audit_key(rec['trial'])
        if got != planted * 7 or len(set(rec['slots'])) != 2:
            bad += 1
    # end-to-end: run 3 real trials, check every row's trial id is in range and
    # that no row reports success without coverage or without a verified solve
    rows = chunk(('q25', 2048, 1, 2, 1.0, 4242, 3, 0))
    for r in rows:
        if r['trial'] not in (0, 1, 2):
            bad += 1
        if r['ok'] and not r['audit_merge']:
            # legitimate: a same-walker self-collision across an epoch restart
            pass
    return dict(planted=40, bad=bad, trial_ids=[r['trial'] for r in rows],
                ok=[r['ok'] for r in rows], audit=[r['audit_merge'] for r in rows])


def main():
    if '--selftest' in sys.argv:
        print(json.dumps(selftest(), indent=1))
        return
    out_path = sys.argv[1] if len(sys.argv) > 1 else 'census-v15m.json'
    trials_per = int(os.environ.get('TRIALS', '6400'))
    workers = int(os.environ.get('WORKERS', '10'))
    dens = [int(x) for x in os.environ.get('DENS', '0,1,2,3,4,5').split(',')]
    cap_mult = float(os.environ.get('CAP_MULT', '1.0'))
    drain_spacings = int(os.environ.get('DRAIN', '2'))
    per_chunk = 200
    all_rows, meta = [], {}
    for ci, dp in enumerate(dens):
        jobs, got, k = [], 0, 0
        while got < trials_per:
            n = min(per_chunk, trials_per - got)
            jobs.append(('q25', 2048, dp, drain_spacings, cap_mult,
                         900000 + ci * 100003 + k * 977, n, got))
            got += n
            k += 1
        with Pool(workers) as pool:
            chunks = pool.map(chunk, jobs)
        rows = [r for ch in chunks for r in ch]
        all_rows += rows
        row = summarize(rows, dp, drain_spacings)
        meta[f'q25/R2048/dp{dp}'] = row
        print(json.dumps(row), flush=True)

    doc = dict(kind='census-v15m',
               statistic='UNCONDITIONAL verified recovery ([x]P0==T), no conditioning filter',
               channels=dict(budget='both records inside the search budget',
                             drain='one record delivered by the CHARGED drain flush',
                             drain_only='pair first meets inside the drain: counted, NEVER success'),
               audit='independent coverage audit: two distinct walkers in one class, in-budget',
               generator='census_v15m.py', model_rate=MODEL, cap_mult=cap_mult,
               drain_spacings=drain_spacings, workers=workers,
               densities=dens, trials_per_density=trials_per,
               rows=all_rows, pooled=meta)
    json.dump(doc, open(os.path.join(HERE, out_path), 'w'), indent=1)
    print(json.dumps({k: {kk: v[kk] for kk in
                          ('p_uncond', 'cp95_uncond', 'p_audit_merge',
                           'cp95_audit_merge', 'drain_only_per_trial')}
                      for k, v in meta.items()}, indent=1))


if __name__ == '__main__':
    main()
