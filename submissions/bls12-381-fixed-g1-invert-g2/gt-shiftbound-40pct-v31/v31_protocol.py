#!/usr/bin/env python3
"""v31 SHIPPED LAUNCH PROTOCOL.  This module is the algorithm, not a harness.

v30's verdict, flaw 2: "The success measurements do not execute the declared
launch protocol.  v30_measure._batch immediately starts replacement trails,
whereas scheduler_trace tests the cadence and atomic cap using trails that always
last the maximum length...  Under the literal one-launch-per-theta-completed-steps
rule, an initial trail reporting before theta can leave no live trail and
insufficient completed steps to authorize another launch; no handling of this
case is specified."

v31 ships ONE loop, `run_protocol`, with three named rules and no other launch
path.  Every other module in this package (ledger31, the regen, the audit) imports
INFLIGHT_CAP and HALT_RULE from HERE, so the bound is a property of the shipped
object:

  RULE 1 (cadence)    a launch is authorized when credit = completed aggregate
                      steps not yet consumed >= LAUNCH_PERIOD.  A launch consumes
                      exactly LAUNCH_PERIOD credit.
  RULE 2 (bootstrap)  if NO trail is live, a launch is authorized immediately and
                      consumes no credit.  This is the case the reviewer named: a
                      trail that reports before theta completes cannot deadlock
                      the schedule.
  RULE 3 (gate)       if a launch is authorized but len(inflight) == INFLIGHT_CAP,
                      the launch is DEFERRED: the credit is retained (not spent,
                      not lost) and the request is retried at the next completion
                      event.  len(inflight) can therefore never exceed
                      INFLIGHT_CAP.  This line is the bound.

HALT RULE (ledger31.stop_loss31 quotes this string): stop authorizing new launches
when completed >= S - reserve, reserve = INFLIGHT_CAP * CAP_MULT * theta.

The trace below is the executable evidence for the gate: it records the tick at
which in-flight first reaches 128 and the launch request that is refused at 128.
"""
import json
import math
import os

THETA = 8                      # scaled clock: theta = 8 steps
CAP_MULT = 128                 # trail cap = 128 * theta own steps
INFLIGHT_CAP = 128             # RULE 3, enforced here in the shipped loop
LAUNCH_PERIOD = THETA
HALT_RULE = ("stop authorizing launches when completed aggregate steps >= S - "
             "reserve, reserve = INFLIGHT_CAP * CAP_MULT * theta")
RESERVE_STEPS = INFLIGHT_CAP * CAP_MULT * THETA


def run_protocol(rounds=4000, inflight_cap=None, trail_len=None,
                 launch_period=None, gate=True, stop_at=None,
                 trace_from=None, trace_max=14):
    """The shipped loop.  Returns counters + the literal trace lines.

    trail_len is the OWN-step cap of a trail (worst case: every trail runs to its
    cap, which is the adversarial case the reserve must cover).  gate=False is a
    control only: it removes RULE 3 and is never a shipped configuration.
    Every default is read from the module globals at CALL time, so patching
    v31_protocol.INFLIGHT_CAP changes the shipped path (control (b)).
    """
    inflight_cap = INFLIGHT_CAP if inflight_cap is None else inflight_cap
    launch_period = LAUNCH_PERIOD if launch_period is None else launch_period
    trail_len = CAP_MULT * THETA if trail_len is None else trail_len
    inflight = []
    completed = 0
    consumed = 0
    accepted = 0
    deferred = 0
    denied_at_cap = 0
    bootstrap = 0
    best = 0
    first_cap_tick = None
    trace = []
    denied_lines = []
    for tick in range(rounds):
        credit = completed - consumed
        want = False
        why = None
        if not inflight:
            want, why = True, 'bootstrap'          # RULE 2
        elif credit >= launch_period:
            want, why = True, 'credit'             # RULE 1
        if want and gate and len(inflight) >= inflight_cap:
            deferred += 1                          # RULE 3
            denied_at_cap += 1
            if len(denied_lines) < trace_max:
                denied_lines.append(
                    't=%d inflight=%d completed=%d credit=%d launch#%d DENIED '
                    'gate: len(inflight)==%d (deferred=%d, credit retained)'
                    % (tick, len(inflight), completed, credit, accepted + 1,
                       inflight_cap, deferred))
            if trace_from is not None and tick >= trace_from and len(trace) < trace_max:
                trace.append(denied_lines[-1])
        elif want:
            if why == 'credit':
                consumed += launch_period
            else:
                bootstrap += 1
            inflight.append(trail_len)
            accepted += 1
            if len(inflight) == inflight_cap and first_cap_tick is None:
                first_cap_tick = tick
                if len(trace) < trace_max:
                    trace.append('t=%d inflight=%d->%d launch#%d ACCEPT '
                                 'completed=%d credit=%d (first time in-flight '
                                 'reaches %d)'
                                 % (tick, len(inflight) - 1, len(inflight),
                                    accepted, completed, credit, inflight_cap))
        # every live trail executes one own step; completions bank credit
        completed += len(inflight)
        still = []
        for x in inflight:
            if x - 1 > 0:
                still.append(x - 1)
        inflight = still
        best = max(best, len(inflight))
        if stop_at is not None and completed >= stop_at:
            break
    return dict(max_inflight=best, accepted=accepted, deferred=deferred,
                denied_at_cap=denied_at_cap, bootstrap=bootstrap,
                completed=completed, rounds=rounds, inflight_cap=inflight_cap,
                trail_len=trail_len, launch_period=launch_period, gate=gate,
                first_cap_tick=first_cap_tick, trace=trace,
                denied_lines=denied_lines)


def controls():
    """(b) enforcement lives in the object; (c) the 128 boundary is exercised."""
    ship = run_protocol()
    loosened = run_protocol(inflight_cap=INFLIGHT_CAP + 1)
    bare = run_protocol(gate=False)
    small = run_protocol(inflight_cap=INFLIGHT_CAP - 1)
    return dict(ship=ship, loosened=loosened, bare=bare, small=small)


def main():
    os.makedirs('evidence', exist_ok=True)
    c = controls()
    out = dict(rule='ONE shipped loop: RULE 1 cadence, RULE 2 bootstrap, RULE 3 '
                    'atomic in-flight gate; no other launch path exists',
               shipped=dict(inflight_cap=INFLIGHT_CAP, cap_mult=CAP_MULT,
                            theta=THETA, trail_cap_steps=CAP_MULT * THETA,
                            launch_period=LAUNCH_PERIOD,
                            reserve_steps=RESERVE_STEPS,
                            reserve_log2=math.log2(RESERVE_STEPS),
                            halt_rule=HALT_RULE),
               traces={k: v for k, v in c.items()},
               gate_lines=c['ship']['denied_lines'][:6],
               hit_lines=c['ship']['trace'][:6])
    json.dump(out, open('evidence/v31-protocol.json', 'w'), indent=1,
              sort_keys=True)
    s = c['ship']
    print('shipped: max_inflight %d accepted %d deferred %d denied_at_cap %d '
          'first_cap_tick %s' % (s['max_inflight'], s['accepted'], s['deferred'],
                                 s['denied_at_cap'], s['first_cap_tick']))
    for line in s['trace'][:4]:
        print('  ' + line)
    for line in s['denied_lines'][:3]:
        print('  ' + line)
    print('patched cap 129 -> max %d ; gate removed -> max %d ; cap 127 -> max %d'
          % (c['loosened']['max_inflight'], c['bare']['max_inflight'],
             c['small']['max_inflight']))


if __name__ == '__main__':
    main()
