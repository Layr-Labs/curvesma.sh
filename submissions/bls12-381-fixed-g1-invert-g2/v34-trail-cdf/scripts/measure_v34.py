#!/usr/bin/env python3
"""v34 measure: the tail on the RIGHT OBJECT, the detection denominator on ALL
class merges, and mu derived from the trace.

WHAT v33 GOT WRONG (verdict-v33, objections 1-3):

  (1) TAIL.  v33 charged q^R with R = min(cap, S - t_launch), the AGGREGATE
      remaining budget, and q^R done PER STEP.  Both are modelling assumptions:
      routing.dp_hit reads a field of the CLASS KEY, so a step's DP outcome is a
      function of the class the step lands on, and revisiting a class repeats its
      DP result.  v34 measures the trail's lifetime directly: T = executed steps
      from a launch to that trail's FIRST DP report, right-censored at the
      per-trail cap or at the budget halt, on the per-trail population of the
      shipped protocol.  The charged row is P(T > R_trail), R_trail the PER-TRAIL
      lifetime cap from the schedule -- `self.cap = cap_mult * self.theta`
      (protocol.Cfg, cap_mult=128), enforced by `limit = cfg.cap` in the step
      loop.  No q, no per-step i.i.d. claim, no aggregate R.

  (2) DETECTION.  v33 counted a merge only after a DP hit AND a store lookup, so
      its denominator was P(recover | shared report).  v34 defines the merge
      event WITHOUT the store: a CLASS-MERGE is a step of trail j landing on a
      class key already visited by a DIFFERENT trail (first such arrival for j).
      Trajectory merges that never reach a DP are IN the denominator.  The
      outcome is resolved on the class map built from the same run: DETECTED iff
      the class's forward orbit reaches a DP class inside the follower's
      remaining lifetime (cap AND budget); the twelve-shift solve is re-run and
      verified against the planted z on a subsample.  The class map is
      DETERMINISTIC, so events at one class share an outcome and pooling them as
      binomials is unjustified: the per-run CP is reported and the MINIMUM
      per-run CP-95 lower bound is charged.

  (3) mu.  v33 imported mu = 11/10.  v34 derives it from the trace:
      mu = completed / (launches * theta), the executed work per launch in units
      of the report interval.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy as TY                  # noqa: E402
import protocol as PR             # noqa: E402
import classkey as CK             # noqa: E402
import routing                    # noqa: E402
from cp95 import cp95_lower       # noqa: E402

THETA_BITS = 5
THETA = 1 << THETA_BITS
CAP_MULT = 128
CAP = CAP_MULT * THETA            # the per-trail lifetime cap, 4096 steps
N_RUNS = 24
BUDGET = 20000
EV = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'evidence'))
SAMPLE_TRAILS = 40
VERIFY_MAX = 256


# ---------------------------------------------------------------------------
# 1. the class map and the DP distance (the object the tail is measured on)
# ---------------------------------------------------------------------------
def class_map(toy, states, theta_bits=THETA_BITS, max_nodes=4000):
    """Successor class + DP flag for every class reachable from the run's states.

    step_v26 is class-equivariant, so the successor CLASS is a function of the
    class alone; a representative (y,a,b) taken from the run is exact.
    """
    succ, states = {}, dict(states)
    todo = list(states.items())
    while todo:
        m, (y, a, b) = todo.pop()
        if m in succ:
            continue
        y2, a2, b2, _k, _i = CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U, toy.r,
                                         toy.tbits)
        m2 = toy.key(y2)[0]
        succ[m] = m2
        if m2 not in states and len(states) < max_nodes:
            states[m2] = (y2, a2, b2)
            todo.append((m2, (y2, a2, b2)))
    dp = {m: routing.dp_hit(m, theta_bits) for m in succ}
    # DP distance: 0 on a DP class, d+1 along the deterministic map, None if the
    # orbit never reaches a DP class (a DP-free cycle, or a path into one).
    dd = {m: (0 if dp[m] else None) for m in succ}
    for _ in range(len(succ) + 2):
        moved = 0
        for m, n in succ.items():
            if dd[m] is None and dd.get(n) is not None:
                dd[m] = dd[n] + 1
                moved += 1
        if not moved:
            break
    return dict(succ=succ, dp=dp, dd=dd, states=states)


def dstar(events, cmap):
    """Resolve every class-merge event against the class map."""
    out = []
    for e in events:
        d = cmap['dd'].get(e['mkey'])
        det = (d is not None) and (d <= e['left'])
        row = dict(e)
        row['dp_dist'] = d
        row['detected'] = bool(det)
        out.append(row)
    return out


def reverify(toy, x, z, events, limit=VERIFY_MAX, sample=150):
    """End-to-end check on a subsample: walk the follower's own (y,a,b) forward to
    the terminal DP class, run the shipped twelve-shift solve against the labels
    the store holds for that class, and verify the candidate against z."""
    gpow = [pow(toy.p, j, toy.r) for j in range(12)]
    ok = bad = skipped = 0
    rows = []
    pool = [e for e in events if e['dp_dist'] is not None
            and e['dp_dist'] <= limit][:sample]
    for e in pool:
        d = e['dp_dist']
        y, a, b = e['state']
        for _ in range(d):
            y, a, b, _k, _i = CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U, toy.r,
                                          toy.tbits)
        rec = e['store'].get(toy.key(y)[0])
        if rec is None:
            skipped += 1
            continue
        hit = None
        for j, xx in CK.solve_candidates(gpow, None, a, b, rec[0], rec[1], toy.r):
            if toy.T.pow(toy.g, xx) == z:
                hit = j
                break
        ok += 1 if hit is not None else 0
        bad += 0 if hit is not None else 1
        rows.append(dict(trail=e['trail'], dp_dist=d, j=hit))
    return dict(checked=ok + bad, recovered=ok, failed=bad,
                skipped=len(events) - len(pool), sample_cap=sample,
                dp_dist_cap=limit, rows=rows[:20])


# ---------------------------------------------------------------------------
# 2. the per-trail lifetime-to-DP CDF
# ---------------------------------------------------------------------------
def exact_T(rows, cmap):
    """EXACT per-trail lifetime-to-DP.  The class map is deterministic, so a trail
    sitting on class c at s steps has T = s + dp_dist(c) exactly: no censoring, no
    estimator.  The only randomness is the launch draw (y = g^a z^b), i.e. which
    basin the trail starts in.  dp_dist(c) = None means the orbit never reaches a
    DP class: T = infinity (the trail never reports, however long it runs)."""
    out = []
    for r in rows:
        d = cmap['dd'].get(r.get('mkey')) if r.get('mkey') is not None else None
        if r['dp_step'] is not None:
            t = r['dp_step']
        elif r.get('mkey') is None:
            t = None
        else:
            t = None if d is None else r['steps'] + d
        out.append(dict(r, T=t, never=(t is None)))
    return out


def cp95_upper(k, n):
    """One-sided 95% Clopper-Pearson UPPER bound, by the binomial symmetry
    cp_upper(k,n) = 1 - cp_lower(n-k, n).  Charged on a PENALTY (delta)."""
    return float(1 - cp95_lower(n - k, n))


def cdf_rows(rows):
    """Empirical survival of the exact T at every quota of the schedule."""
    out = []
    for mult in (1, 2, 4, 8, 16, 32, 64, 128):
        q = mult * THETA
        over = sum(1 for r in rows if r['T'] is None or r['T'] > q)
        out.append(dict(quota=q, theta_mult=mult, n=len(rows), over_q=over,
                        surv_qgt=over / len(rows), cp95_upper=cp95_upper(over, len(rows)),
                        max_finite_T=max((r['T'] for r in rows if r['T'] is not None),
                                         default=None)))
    return out


# ---------------------------------------------------------------------------
# 3. the runs
# ---------------------------------------------------------------------------
def main():
    os.makedirs(EV, exist_ok=True)
    trails, per_run, events_by_seed = [], [], {}
    for i in range(N_RUNS):
        seed = 1000 + 7 * i
        toy = TY.Toy(bits=5, seed=seed)
        x = toy.rng.randrange(1, toy.r)
        z = toy.T.pow(toy.g, x)
        states, tl, cl = {}, [], []
        orig = CK.step_v26

        def patched(Tl, y, a, b, TBL, U, r, tbits=None, _st=states, _o=orig):
            m = toy.key(y)[0]
            if m not in _st:
                _st[m] = (y, a, b)
            return _o(Tl, y, a, b, TBL, U, r, tbits)

        CK.step_v26 = patched
        try:
            res = PR.run_protocol(toy, z, PR.Cfg(THETA_BITS, BUDGET),
                                  trail_log=tl, collision_log=cl,
                                  stop_on_success=False)
        finally:
            CK.step_v26 = orig
        cmap = class_map(toy, states)
        ev = dstar(cl, cmap)
        for e in ev:
            e['store'] = res['store'].d
        rv = reverify(toy, x, z, ev)
        clean = [{k: v for k, v in e.items() if k not in ('state', 'store')}
                 for e in ev]
        tl = exact_T(tl, cmap)
        events_by_seed[seed] = clean
        trails.extend(dict(r, seed=seed) for r in tl)
        per_run.append(dict(seed=seed, launches=res['launches'],
                            completed=res['completed'], consumed=res['consumed'],
                            merged=res['merged'], detected=res['detected'],
                            reports=res['reports'], halt=res['halt'],
                            executed_equals_S=res['executed_equals_S'],
                            consumed_equals_executed=res['consumed_equals_executed'],
                            launches_le_bound=res['launches_le_bound'],
                            mu=res['completed'] / (res['launches'] * THETA),
                            n_trails=len(tl),
                            n_dp=sum(1 for r in tl if r['dp_step'] is not None),
                            n_cap=sum(1 for r in tl if r['exit'] == 'cap'),
                            n_halt=sum(1 for r in tl if r['exit'] == 'halt'),
                            distinct_classes=len(cmap['succ']), verify=rv))
        print('run %2d seed %d trails %3d dp %3d events %3d det %3d reports %2d classes %3d'
              % (i, seed, len(tl), per_run[-1]['n_dp'], len(clean),
                 sum(1 for e in clean if e['detected']), res['reports'],
                 len(cmap['succ'])), flush=True)

    rows = cdf_rows(trails)
    never = sum(1 for r in trails if r['T'] is None)
    surv_cap = never / len(trails)
    delta_cp_upper = cp95_upper(never, len(trails))
    per_run_merge = []
    for seed, ev in events_by_seed.items():
        n = len(ev)
        k = sum(1 for e in ev if e['detected'])
        per_run_merge.append(dict(seed=seed, events=n, detected=k,
                                  p=(k / n) if n else None,
                                  cp95_lower=float(cp95_lower(k, n)) if n else None))
    tot = sum(r['events'] for r in per_run_merge)
    tot_det = sum(r['detected'] for r in per_run_merge)
    shared_matches = sum(r['merged'] for r in per_run)
    valids = [r['cp95_lower'] for r in per_run_merge if r['cp95_lower'] is not None]
    min_cp = min(valids) if valids else None
    pooled_cp = float(cp95_lower(tot_det, tot)) if tot else None
    p_shared = shared_matches / tot if tot else None
    cls, det_cls = {}, {}
    for ev in events_by_seed.values():
        for e in ev:
            cls[e['mkey']] = cls.get(e['mkey'], 0) + 1
            if e['detected']:
                det_cls[e['mkey']] = det_cls.get(e['mkey'], 0) + 1
    reverify_all = dict(
        checked=sum(r['verify']['checked'] for r in per_run),
        recovered=sum(r['verify']['recovered'] for r in per_run),
        failed=sum(r['verify']['failed'] for r in per_run),
        skipped=sum(r['verify']['skipped'] for r in per_run))

    tail = dict(
        theta=THETA, cap=CAP, cap_mult=CAP_MULT, budget=BUDGET,
        cap_line='protocol.Cfg.__init__: self.cap = cap_mult * self.theta (cap_mult=128)',
        cap_enforced_line='protocol.run_protocol: limit = cfg.cap; '
                          'if tr.steps < limit: still.append(tr) else exit=cap',
        n_trails=len(trails), n_runs=len(per_run),
        n_dp_reports=sum(1 for r in trails if r['dp_step'] is not None),
        n_censored_cap=sum(1 for r in trails if r['exit'] == 'cap'),
        n_censored_halt=sum(1 for r in trails if r['exit'] == 'halt'),
        cdf=rows, p_T_gt_cap=surv_cap, delta_cp95_upper=delta_cp_upper,
        n_never=never,
        charged=delta_cp_upper,
        charged_kind='CP-95 UPPER bound on the empirical P(T > R_trail), '
                     'R_trail = 128*theta = the per-trail cap; T is EXACT (the '
                     'class map is deterministic, so dp_dist resolves every trail '
                     'without censoring); upper bound charged because delta is a '
                     'PENALTY in the sizing law',
        old_v33_row=dict(form='mean of q^R with R = min(cap, S - t_launch) AGGREGATE',
                         value=3.790e-4,
                         defect='wrong random object: R is not this trail\'s '
                                'lifetime and the steps are not independent'),
        mu_per_run=[r['mu'] for r in per_run],
        mu_derived=sum(r['mu'] for r in per_run) / len(per_run),
        mu_line='mu = completed / (launches * theta) per run, averaged',
        mean_distinct_classes=sum(r['distinct_classes'] for r in per_run) / len(per_run),
        mean_distinct_per_trail=sum(r['distinct'] for r in trails) / len(trails),
        sample_trails=[{k: r[k] for k in ('seed', 'tid', 'steps', 'dp_step', 'exit',
                                          'distinct')} for r in trails[:SAMPLE_TRAILS]],
        per_run=per_run,
    )
    merge = dict(
        event_definition='a step of trail j landing on a class key already visited by '
                         'a DIFFERENT trail (first such arrival for j); NO store and '
                         'NO DP requirement, so misses are in the denominator',
        events_total=tot, detected_total=tot_det,
        p_detect_all_merges=(tot_det / tot) if tot else None,
        pooled_cp95_lower=pooled_cp, min_per_run_cp95_lower=min_cp,
        charged=min_cp,
        charged_kind='MIN over runs of the per-run CP-95 lower bound (pooling '
                     'refused: merge destination classes are shared, so events are '
                     'not exchangeable binomials)',
        per_run=per_run_merge,
        distinct_merge_classes=len(cls),
        top_classes=sorted(cls.items(), key=lambda kv: -kv[1])[:10],
        distinct_detected_classes=len(det_cls),
        old_v33_object=dict(merged=shared_matches,
                            note='merges that reached a DP AND a store record: '
                                 'P(recover | shared report), the selection-biased '
                                 'object v33 charged'),
        selection_bias_delta=((tot_det / tot) - p_shared) if tot else None,
        p_shared_report=p_shared,
        reverify=reverify_all,
    )
    json.dump(tail, open(os.path.join(EV, 'tail34.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    json.dump(merge, open(os.path.join(EV, 'merge34.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    print('\nTRAILS %d  dp %d  cap-censored %d  halt-censored %d'
          % (tail['n_trails'], tail['n_dp_reports'], tail['n_censored_cap'],
             tail['n_censored_halt']))
    for r in rows:
        print('  quota %5d  P(T>q) %.6f  CP95-up %.6f  over %d/%d  max finite T %s'
              % (r['quota'], r['surv_qgt'], r['cp95_upper'], r['over_q'], r['n'],
                 r['max_finite_T']))
    print('CHARGED P(T>%d) point %.6f  CP95-upper %.6f  never %d/%d  (v33 q^R 3.790e-04)'
          % (CAP, surv_cap, delta_cp_upper, never, len(trails)))
    print('MU derived %.6f  (v33 imported 1.100), per-run %s'
          % (tail['mu_derived'], ['%.3f' % v for v in tail['mu_per_run'][:6]]))
    print('MERGES all %d det %d p %.6f pooled CP-lo %.6f MIN per-run CP-lo %s'
          % (tot, tot_det, merge['p_detect_all_merges'], pooled_cp, min_cp))
    print('  old object p_shared %.6f   selection-bias delta %.6f'
          % (p_shared, merge['selection_bias_delta']))
    print('  distinct merge classes %d   reverify %s' % (len(cls), reverify_all))


if __name__ == '__main__':
    main()
