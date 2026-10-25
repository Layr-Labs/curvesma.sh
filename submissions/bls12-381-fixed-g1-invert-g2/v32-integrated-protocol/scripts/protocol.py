#!/usr/bin/env python3
"""v32 THE PROTOCOL -- the single launch path.  Every walk-dependent number in the
package comes out of `run_protocol`; there is no replacement loop, no sampler fork
and no harness-side scheduler.

PER TRAIL, IN ORDER (verdict v31 defect 1):
  1. step_v26 walk (`classkey.step_v26`, the family-free squared step)
  2. distinguished-point termination, tested BEFORE any store lookup
  3. report insertion / report matching on the class key
  4. the twelve-shift recovery (`classkey.solve_candidates`) and its verification
  5. credit accounting and the in-flight gate, bootstrap launches included

CREDIT (verdict v31 defect 1, second half).  v31's RULE 2 launched a bootstrap
trail for free, so ceil(S/theta) was not an enforced launch bound: trails ending
before theta could re-bootstrap for ever.  v32 charges EVERY launch -- bootstrap
included -- from the same S budget, and clamps the launch COUNT at
launch_bound = ceil(S/theta).  Both ride the same counter.
"""
import classkey as CK

INFLIGHT_CAP = 128          # simultaneous trails, enforced inside run_protocol
CAP_MULT = 128              # per-trail step cap, in units of theta

HALT_RULE = (
    "RULE 1 credit = completed - consumed; a launch needs credit >= theta. "
    "RULE 2 an empty in-flight set triggers a bootstrap launch and that launch "
    "CONSUMES theta from the same S budget (v32; v31 charged nothing). "
    "RULE 3 a credit-eligible launch that would exceed INFLIGHT_CAP is DEFERRED; "
    "the credit is kept and re-offered next round. "
    "RULE 4 total launches are clamped at launch_bound = ceil(S/theta), bootstrap "
    "included. "
    "RULE 5 halt on the first verified recovery, or when the launch bound is "
    "reached and the in-flight set has drained.")


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
    __slots__ = ('tid', 'y', 'a', 'b', 'steps', 't0')

    def __init__(self, tid, y, a, b, t0):
        self.tid, self.y, self.a, self.b = tid, y, a, b
        self.steps = 0
        self.t0 = t0


def run_protocol(toy, z, cfg, store=None, stop_on_success=True,
                 suppress_shared=False, skip_dp=False, cut_steps=0,
                 trace_n=0, trace_from=0):
    """Execute the shipped protocol for one challenge z.

    Planted controls (never the shipped path; each is one audit):
      suppress_shared  v31's sampler defect: a DP whose key is already stored does
                       `continue` instead of terminating and matching, so the
                       shared-report event is suppressed and repeat visits counted.
      skip_dp          no DP termination at all  -> audit (a), success collapses.
      cut_steps        every trail is cut at cut_steps  -> audit (c), launches
                       outrun the budget and the launch bound must bind.
    """
    T, r = toy.T, toy.r
    store = Store() if store is None else store
    trace = []
    g = toy.g
    gpow = [pow(toy.p, j, r) for j in range(12)]
    t = 0
    consumed = 0
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

    def launch(now):
        nonlocal tid
        a = rng.randrange(1, r)
        b = rng.randrange(1, r)
        y = T.mul(T.pow(g, a), T.pow(z, b))
        inflight.append(Trail(tid, y, a, b, now))
        tid += 1

    while True:
        rounds += 1
        # ---------------- launch decision: RULE 1..4 --------------------------
        if launches >= cfg.launch_bound:
            event = 'bound'
            if not inflight:
                halt = 'launch_bound'
                break
        elif not inflight:
            launch(t)                                  # RULE 2, charged below
            consumed += cfg.theta
            launches += 1
            n_bootstrap += 1
            event = 'bootstrap'
        elif t - consumed >= cfg.theta:
            if len(inflight) >= cfg.inflight_cap:
                deferred += 1
                event = 'defer'
            else:
                launch(t)
                consumed += cfg.theta
                launches += 1
                event = 'accept'
        else:
            event = 'idle'
        if trace_n and rounds >= trace_from and len(trace) < trace_n:
            trace.append(dict(round=rounds, t=t, event=event,
                              inflight=len(inflight), launches=launches,
                              completed=t, consumed=consumed, credit=t - consumed,
                              deferred=deferred))

        if not inflight:
            halt = halt or ('launch_bound' if launches >= cfg.launch_bound
                            else 'drained')
            break
        if t >= cfg.S and not inflight:
            halt = 'budget'
            break

        # ---------------- one lock-step round --------------------------------
        still = []
        for tr in inflight:
            tr.y, tr.a, tr.b, kind, _info = CK.step_v26(
                T, tr.y, tr.a, tr.b, toy.TBL, toy.U, r, toy.tbits)
            tr.steps += 1
            t += 1
            is_dp = (not skip_dp) and toy.dp(tr.y, cfg.theta_bits)
            limit = cfg.cap if not cut_steps else min(cfg.cap, cut_steps)
            if not is_dp:
                if tr.steps < limit:
                    still.append(tr)
                continue
            # ---- DP TERMINATION, BEFORE ANY STORE LOOKUP (the v32 order) -----
            mkey = toy.key(tr.y)[0]
            prior = store.d.get(mkey)
            if suppress_shared and prior is not None:
                store.repeat_visits += 1
                if tr.steps < limit:
                    still.append(tr)
                continue
            if prior is None:
                store.put(mkey, (tr.a, tr.b, tr.steps, tr.t0), t)
                continue                    # one report per trail, then it stops
            # ---- shared report: match, twelve-shift recovery, verify --------
            store.matches += 1
            store.shared.add(mkey)
            store.ages.append(t - store.write_t[mkey])
            cands = CK.solve_candidates(gpow, None, tr.a, tr.b, prior[0], prior[1], r)
            for j, xx in cands:
                if T.pow(g, xx) == z:
                    if not recovered:
                        recovered = True
                        match = dict(age=store.ages[-1], steps_to_report=tr.steps,
                                     prior_steps=prior[2], j=j,
                                     candidates=len(cands), trail=tr.tid,
                                     prior_trail=prior[3] and None)
                    break
            if recovered and stop_on_success:
                halt = 'recovered'
                inflight = []
                break
        inflight = still
        if halt == 'recovered':
            break

    return dict(cfg=cfg.as_dict(), recovered=recovered, halt=halt,
                launches=launches, n_bootstrap=n_bootstrap, deferred=deferred,
                completed=t,
                consumed=consumed, rounds=rounds, reports=store.reports,
                matches=store.matches, repeat_visits=store.repeat_visits,
                credit_left=t - consumed,
                launch_bound=cfg.launch_bound,
                launches_le_bound=launches <= cfg.launch_bound,
                launch_bound_binds=(launches == cfg.launch_bound),
                match=match, trace=trace, store=store, ages=store.ages,
                lengths=store.lengths)


def retention(res):
    """First-report retention, measured on the integrated run: the fraction of
    inserted reports that at least one later trail SHARED.  A report shared k
    times counts once, so the statistic is a probability, not a visit count --
    which is exactly what v31's dedup-first sampler could not measure."""
    st = res['store']
    return (len(st.shared) / st.reports) if st.reports else None
