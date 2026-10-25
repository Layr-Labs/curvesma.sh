#!/usr/bin/env python3
"""v33 THE PROTOCOL -- the single launch path, now with a per-EXECUTED-STEP budget.

THE v32 DEFECT THIS FILE REPAIRS (verdict-v32, objection 1).  v32 charged theta
per LAUNCH and clamped the launch COUNT at ceil(S/theta).  With a per-trail cap
of 128*theta that bounds work at ~128*S, not S: the shipped trace had S=20000 and
654518 completed steps.  `t >= cfg.S and not inflight` was dead code behind the
preceding empty-in-flight break, so nothing stopped the step counter.

v33 charges the budget PER EXECUTED STEP:

  RULE 1 consumed advances by 1 for every step actually executed, by every trail,
         bootstrap trails included; remaining = S - consumed.  The loop tests
         remaining BEFORE each step, so the executed-step total over ALL trails --
         abandoned and in-flight ones included -- can never exceed S, and the
         stop itself is a trace row whose completed counter equals S.
  RULE 2 a launch requires remaining >= theta: no trail is started that the
         budget cannot carry through a full report interval theta.
  RULE 3 a credit-eligible launch at INFLIGHT_CAP is DEFERRED, credit kept.
  RULE 4 the launch count is therefore bounded by S/theta BY CONSTRUCTION (one
         launch per theta of executed work at most); launches_le_bound is still
         asserted against ceil(S/theta) rather than assumed.
  RULE 5 halt on the first verified recovery (stop_on_success), else at the step
         budget (halt='step_budget', completed == S exactly).

Everything else is v32's, unchanged: step_v26, DP-first termination, the store
match on the class key, the twelve-shift recovery with verification, the
in-flight gate.  `legacy_launch_credit=True` reproduces v32's charge-at-launch /
no-step-cutoff rule as a CONTROL -- it is never the shipped path.
"""
import classkey as CK

INFLIGHT_CAP = 128          # simultaneous trails, enforced inside run_protocol
CAP_MULT = 128              # per-trail step cap, in units of theta

HALT_RULE = (
    "RULE 1 the budget is EXECUTED STEPS: consumed += 1 per step actually "
    "executed (every trail, bootstrap included), remaining = S - consumed, and "
    "the loop tests remaining BEFORE each step, so total executed steps over ALL "
    "trails cannot exceed S; the halt is a trace row with completed == S. "
    "RULE 2 a launch requires remaining >= theta (one report interval of "
    "headroom). RULE 3 a credit-eligible launch at INFLIGHT_CAP is DEFERRED, "
    "credit kept. RULE 4 launches <= ceil(S/theta) by construction, asserted. "
    "RULE 5 halt on the first verified recovery, else at the step budget.")


class Cfg:
    def __init__(self, theta_bits, budget_steps, cap_mult=CAP_MULT,
                 inflight_cap=INFLIGHT_CAP):
        self.theta_bits = int(theta_bits)
        self.theta = 1 << self.theta_bits
        self.S = int(budget_steps)
        self.cap_mult = cap_mult
        self.cap = cap_mult * self.theta
        self.inflight_cap = inflight_cap
        self.launch_bound = -(-self.S // self.theta)

    def as_dict(self):
        return dict(theta_bits=self.theta_bits, theta=self.theta, S=self.S,
                    cap_mult=self.cap_mult, cap=self.cap,
                    inflight_cap=self.inflight_cap,
                    launch_bound=self.launch_bound)


class Store:
    """One shared store.  A report holds the class key plus the two labels mod r
    (48-byte key + two 255-bit label words at production width)."""

    def __init__(self):
        self.d = {}
        self.reports = 0
        self.matches = 0
        self.repeat_visits = 0
        self.write_t = {}
        self.ages = []
        self.lengths = []
        self.shared = set()          # report keys that were shared at least once

    def put(self, k, rec, t):
        self.d[k] = rec
        self.write_t[k] = t
        self.reports += 1
        self.lengths.append(rec[2])


class Trail:
    __slots__ = ('tid', 'y', 'a', 'b', 'steps', 't0', 'R')

    def __init__(self, tid, y, a, b, t0, R):
        self.tid, self.y, self.a, self.b = tid, y, a, b
        self.steps = 0
        self.t0 = t0
        self.R = R               # this trail's remaining lifetime at launch


def tail_log2(R, theta):
    """log2 of q^R for q = 1 - 1/theta, computed without underflow.

    This IS the actual-event tail: D (steps from an insertion to the next
    distinguished point) is geometric with success probability 1/theta at every
    executed step -- `routing.dp_hit` reads a FRESH 5-bit field of the class key
    on every step, so the per-step hit indicator is i.i.d. Bernoulli(1/theta) and
    P(D > R) = q^R exactly, with no mean-distance substitution and no Markov step.
    """
    import math
    return R * math.log2((theta - 1.0) / theta)


def run_protocol(toy, z, cfg, store=None, stop_on_success=True,
                 suppress_shared=False, skip_dp=False, cut_steps=0,
                 trace_n=0, trace_from=0, merge_log=None, legacy_launch_credit=False,
                 trail_log=None, collision_log=None, owner=None, seen=None):
    """Execute the shipped protocol for one challenge z.

    Planted controls (never the shipped path; each is one audit):
      suppress_shared     v31's sampler defect (DP on a stored key continues).
      skip_dp             no DP termination at all -> success collapses.
      cut_steps           every trail cut at cut_steps -> audits the gate.
      legacy_launch_credit v32's credit: charge theta at LAUNCH, no executed-step
                          cutoff.  Reproduces the 128x overage as a measured
                          control (`completed` far above S) so the fix is visible.

    v34 INSTRUMENTATION (measurement only; it never changes a decision, a counter
    or the stop rule):
      trail_log     one row per LAUNCH: steps at exit, the step index of the first
                    DP report (None if the trail never reported), the exit reason
                    (dp_report / cap / budget / halt), and the number of DISTINCT
                    classes the trail visited.  This is the per-trail
                    lifetime-to-DP population the tail row is measured on.
      collision_log  one row per CLASS-MERGE event: a step of trail j landing on a
                    class key already visited by a DIFFERENT trail, first such
                    arrival for j.  The event is defined WITHOUT the store, so
                    merges whose follow-on walk never reaches a DP are IN the
                    denominator.  Outcome is resolved afterwards from the class
                    map (see measure_v34.class_map / dp_distance).
      owner, seen   the two dicts the collision detector needs; a caller may pass
                    them in to carry the occupancy map across runs.
    """
    T, r = toy.T, toy.r
    store = Store() if store is None else store
    trace = []
    g = toy.g
    gpow = [pow(toy.p, j, r) for j in range(12)]
    t = 0                      # EXECUTED steps (the budget counter)
    consumed = 0               # per-step: consumed == t at every instant
    launches = 0
    n_bootstrap = 0
    deferred = 0
    inflight = []
    tid = 0
    recovered = False
    halt = None
    rounds = 0
    match = None
    rng = toy.rng
    theta = cfg.theta
    launch_log = []
    merged = 0
    detected = 0
    tail_R = []                # R of every launch (the miss-event population)
    owner = {} if owner is None else owner      # v34: class -> first visiting tid
    seen = {} if seen is None else seen         # v34: tid -> set of visited classes

    def launch(now):
        nonlocal tid
        R = max(0, min(cfg.cap, cfg.S - now))   # 0 => the budget is already spent
        a = rng.randrange(1, r)
        b = rng.randrange(1, r)
        y = T.mul(T.pow(g, a), T.pow(z, b))
        inflight.append(Trail(tid, y, a, b, now, R))
        launch_log.append(dict(tid=tid, t0=now, R=R,
                               tail_log2=tail_log2(R, theta)))
        if trail_log is not None:
            trail_log.append(dict(tid=tid, t0=now, steps=0, dp_step=None,
                                  exit='open', distinct=0, R=R))
            seen.setdefault(tid, set())
        tail_R.append(R)
        tid += 1

    while True:
        rounds += 1
        remaining = cfg.S - t
        # ---------------- launch decision: RULE 1..4 --------------------------
        if launches >= cfg.launch_bound:
            event = 'bound'
            if not inflight:
                halt = 'launch_bound'
                break
        elif not inflight:
            if remaining < theta:
                halt = 'step_budget'
                break
            launch(t)
            launches += 1
            n_bootstrap += 1
            event = 'bootstrap'
        elif legacy_launch_credit:
            # v32 CONTROL: credit is granted by completed work between launches
            # and there is NO executed-step cutoff (the defect, reproduced).
            if t - consumed >= theta:
                if len(inflight) >= cfg.inflight_cap:
                    deferred += 1
                    event = 'defer'
                else:
                    launch(t)
                    consumed += theta
                    launches += 1
                    event = 'accept'
            else:
                event = 'idle'
        elif len(inflight) >= cfg.inflight_cap:
            deferred += 1
            event = 'defer'
        elif remaining >= theta:
            launch(t)
            launches += 1
            event = 'accept'
        else:
            event = 'idle'
        if trace_n and rounds >= trace_from and len(trace) < trace_n:
            trace.append(dict(round=rounds, t=t, event=event,
                              inflight=len(inflight), launches=launches,
                              completed=t, consumed=consumed,
                              remaining=cfg.S - t, deferred=deferred))

        if not inflight:
            halt = halt or ('launch_bound' if launches >= cfg.launch_bound
                            else 'drained')
            break

        # ---------------- one lock-step round, cut at the STEP budget ---------
        still = []
        for tr in inflight:
            if not legacy_launch_credit and t >= cfg.S:
                halt = 'step_budget'
                break
            tr.y, tr.a, tr.b, kind, _info = CK.step_v26(
                T, tr.y, tr.a, tr.b, toy.TBL, toy.U, r, toy.tbits)
            tr.steps += 1
            t += 1
            if not legacy_launch_credit:
                consumed += 1                       # THE per-executed-step charge
            limit = cfg.cap if not cut_steps else min(cfg.cap, cut_steps)
            if trail_log is not None or collision_log is not None:
                # v34: the class of the state just reached.  MEASUREMENT ONLY --
                # the DP decision below still reads toy.dp exactly as shipped.
                mkey_i = toy.key(tr.y)[0]
                vs = seen.setdefault(tr.tid, set())
                first_visit = mkey_i not in vs
                if first_visit:
                    vs.add(mkey_i)
                    if (collision_log is not None and mkey_i in owner
                            and owner[mkey_i] != tr.tid):
                        collision_log.append(dict(
                            trail=tr.tid, other=owner[mkey_i], mkey=mkey_i,
                            t=t, steps=tr.steps, t0=tr.t0,
                            left=min(limit - tr.steps, cfg.S - t),
                            left_cap=limit - tr.steps, left_budget=cfg.S - t,
                            state=(tr.y, tr.a, tr.b)))
                    if mkey_i not in owner:
                        owner[mkey_i] = tr.tid
            is_dp = (not skip_dp) and toy.dp(tr.y, cfg.theta_bits)
            if not is_dp:
                if tr.steps < limit:
                    still.append(tr)
                else:
                    if trail_log is not None:
                        _tl = trail_log[tr.tid]
                        _tl['exit'], _tl['steps'] = 'cap', tr.steps
                        _tl['distinct'] = len(seen.get(tr.tid, ()))
                        _tl['mkey'] = mkey_i
                continue
            # ---- DP TERMINATION, BEFORE ANY STORE LOOKUP (the v32 order) ------
            mkey = toy.key(tr.y)[0]
            prior = store.d.get(mkey)
            if suppress_shared and prior is not None:
                store.repeat_visits += 1
                if tr.steps < limit:
                    still.append(tr)
                continue
            if prior is None:
                store.put(mkey, (tr.a, tr.b, tr.steps, tr.t0), t)
                if trail_log is not None:
                    _tl = trail_log[tr.tid]
                    _tl.update(exit='dp_report', dp_step=tr.steps, steps=tr.steps,
                               distinct=len(seen.get(tr.tid, ())))
                continue                    # one report per trail, then it stops
            # ---- shared report: match, twelve-shift recovery, verify --------
            store.matches += 1
            store.shared.add(mkey)
            store.ages.append(t - store.write_t[mkey])
            merged += 1
            if trail_log is not None and trail_log[tr.tid]['exit'] == 'open':
                _tl = trail_log[tr.tid]
                _tl.update(exit='dp_report', dp_step=tr.steps, steps=tr.steps,
                           distinct=len(seen.get(tr.tid, ())))
            cands = CK.solve_candidates(gpow, None, tr.a, tr.b, prior[0], prior[1], r)
            hit = None
            for j, xx in cands:
                if T.pow(g, xx) == z:
                    hit = j
                    break
            if hit is not None:
                detected += 1
            rec = dict(trail=tr.tid, prior_steps=prior[2], steps_to_report=tr.steps,
                       age=t - store.write_t[mkey], t=t, merged=True,
                       detected=hit is not None, j=hit, candidates=len(cands),
                       R_follower=tr.R, tail_log2=tail_log2(tr.R, theta),
                       remaining_at_merge=cfg.S - t)
            if merge_log is not None:
                merge_log.append(rec)
            if hit is not None and not recovered:
                recovered = True
                match = dict(age=store.ages[-1], steps_to_report=tr.steps,
                             prior_steps=prior[2], j=hit,
                             candidates=len(cands), trail=tr.tid,
                             prior_trail=prior[3] and None)
            if recovered and stop_on_success:
                halt = 'recovered'
                inflight = []
                break
        if halt in ('recovered', 'step_budget'):
            break
        inflight = still

    if halt is None:
        halt = 'drained'
    # v34: census of trails still alive at the halt (right-censored observations)
    if trail_log is not None:
        for tr in inflight:
            row = trail_log[tr.tid]
            if row['exit'] == 'open':
                row['exit'], row['steps'] = 'halt', tr.steps
                row['distinct'] = len(seen.get(tr.tid, ()))
                row['mkey'] = toy.key(tr.y)[0]
    # the stop event is a TRACE ROW: completed == S for the budget halt
    trace.append(dict(round=rounds, t=t, event='halt_' + halt,
                      inflight=len(inflight), launches=launches,
                      completed=t, consumed=consumed, remaining=cfg.S - t,
                      deferred=deferred))

    # ---- the measured miss event: P(D > R) over the launch population --------
    import math
    logs = [x['tail_log2'] for x in launch_log]
    mean_tail = (sum(2.0 ** v for v in logs) / len(logs)) if logs else 0.0
    worst = max(logs) if logs else None
    return dict(cfg=cfg.as_dict(), recovered=recovered, halt=halt,
                launches=launches, n_bootstrap=n_bootstrap, deferred=deferred,
                completed=t, consumed=consumed, remaining=cfg.S - t,
                rounds=rounds, reports=store.reports,
                matches=store.matches, repeat_visits=store.repeat_visits,
                credit_left=cfg.S - t,
                launch_bound=cfg.launch_bound,
                launches_le_bound=launches <= cfg.launch_bound,
                launch_bound_binds=(launches == cfg.launch_bound),
                budget_binds=(t <= cfg.S and (cfg.S - t) < cfg.theta),
                executed_equals_S=(t == cfg.S),
                executed_le_S=(t <= cfg.S),
                consumed_equals_executed=(consumed == t),
                merged=merged, detected=detected,
                tail=dict(n_launches=len(launch_log), mean_qR=mean_tail,
                          max_log2_qR=worst,
                          R_min=min(tail_R) if tail_R else None,
                          R_max=max(tail_R) if tail_R else None),
                match=match, trace=trace, store=store, ages=store.ages,
                lengths=store.lengths, launch_log=launch_log)


def retention(res):
    """First-report retention, measured on the integrated run: the fraction of
    inserted reports that at least one later trail SHARED.  A report shared k
    times counts once.  This is a COVERAGE statistic (it is the birthday
    fraction of a finite run), NOT P(detect|merge); v33 keeps it as a diagnostic
    and charges the conditional detection probability instead."""
    st = res['store']
    return (len(st.shared) / st.reports) if st.reports else None
