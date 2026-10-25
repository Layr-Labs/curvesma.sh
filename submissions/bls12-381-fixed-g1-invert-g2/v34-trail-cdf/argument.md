# Argument: the tail on the trail, the denominator on every merge, mu and the constants derived (v34)

v33 (ticket 90701ee2) was **fully judged and rejected**. The judge accepted two of
its repairs by name -- the per-executed-step budget and the single C_REF division
-- and named three surviving objections plus one accounting gap. v34 repairs
exactly those four. `run_protocol` is still the only execution path: DP-first
termination, the twelve-shift recovery with verification, the in-flight gate, the
instruction-counted pairing row and the one-division ledger are unchanged.

Claim **127.4831** = 1274831/10000, the least 4-decimal log2 c with
2^c >= W/C_REF, where W is the exact INSTRUCTION total of
`evidence/ledger34.json` (W/C_REF = 2^127.483093 target ops after the one
division). Emitted by `scripts/regen_v34.py`; no claim literal is typed.

## 1. The tail is measured on the TRAIL, not on a schedule model (objection 1)

v33 charged `q^R` with `R = min(cap, S - t_launch)` -- the AGGREGATE remaining
budget -- once per step, with `q = 1 - 1/theta`. All three premises are false:
`routing.dp_hit` reads a field of the CLASS KEY, so a step's DP outcome is a
function of the class it lands on and a revisit repeats the earlier answer; R was
not this trail's lifetime; and launches happen per round, not per theta of
executed work.

v34 measures the object directly. `run_protocol` gains a `trail_log` (one row per
launch: steps at exit, step index of the first DP report, exit reason, distinct
classes visited) and a `collision_log`; neither touches a decision, a counter or
the stop rule. Then, because `step_v26` is class-equivariant, the successor CLASS
is a deterministic function of the class, so every trail's lifetime is resolved
EXACTLY: a trail sitting on class c after s steps has

    T = s + dp_dist(c),   dp_dist(c) = 0 on a DP class, 1 + dp_dist(succ(c)) otherwise,
                          dp_dist = infinity on a DP-free cycle.

No censoring, no estimator, no i.i.d.-per-step claim. The randomness is the launch
draw (`y = g^a z^b`), i.e. which basin the trail starts in. On 24 runs
(S = 20000, theta = 32, cap = 128*theta):

    trails            9844      (7981 report, 1790 never)
    largest finite T  61        (the survival is FLAT above q = 32)
    P(T > 128*theta)  1790/9844 = 0.181837  pooled
                      CP-95 upper bound 0.188351
                      worst run 0.806646 (per-run never-rates span 0.002..0.753)
    mu = completed/(launches*theta) = 1.906489 (per-run 0.896..3.676)

The per-trail cap is the schedule line `self.cap = cap_mult * self.theta`
(`protocol.Cfg`, `cap_mult = 128`, 4096 steps), enforced by `limit = cfg.cap` in
the step loop. No aggregate R appears anywhere in v34.

**What is charged, and why delta = 0.** The lifetime loss is real, but it is not a
separate multiplicative event: a merge whose orbit never reaches a DP is exactly
what makes P(detect|merge) < 1, and v34 measures that denominator directly
(section 2). Charging it twice would count one loss twice. So delta is deleted --
the verdict's own first remedy ("delete the delta factor entirely and make the
protocol's stop rule the only lifetime statement") -- and the loss is charged
once, on the merge. Both factorizations are printed:

    frame A  delta = 0,        p_det = merge-all rate at its cluster lower bound = 0.495525  <-- CHARGED
    frame B  delta = 0.806646, p_det = detection given a report = 1.000000     -> effective 0.193354  INFEASIBLE

Frame B is not merely conservative, it is impossible: (1 - delta) = 0.1934 < 0.392,
so no S reaches the target at all. That is a statement about the degenerate
785-class fixture (the worst seed puts three trails in four on DP-free cycles),
not about production, and it is disclosed rather than hidden inside a mean.
Frame A is the same measurement with the loss carried where it occurs, and it is
the frame that is both feasible and faithful. The two factorizations agree closely
where both are defined (0.7887 vs 0.8066 unpenalised), which is the cross-check
that no loss is invented or dropped.

## 2. Merges are counted WITHOUT the store (objection 2)

v33 counted a merge only after a DP hit AND a successful store lookup, so
1054/1054 measured P(recover | shared report). v34 defines the event with no store
and no DP requirement:

> a CLASS-MERGE is a step of trail j landing on a class key already visited by a
> DIFFERENT trail (the first such arrival for j).

Undetected trajectory merges are therefore inside the denominator. Each event is
resolved on the same run's class map: DETECTED iff `dp_dist(class) <= left`, where
`left = min(cap - steps, S - t)` is the follower's remaining lifetime; a
150-event-per-run subsample is re-verified end to end (forward walk from the
follower's own (y,a,b) to the terminal DP class, twelve-shift solve against the
stored labels, `g^xx == z`): **3600 checked, 3600 recovered, 0 failed**.

    merge events          157076
    detected              123893
    P(detect|merge)       0.788746      pooled CP-95 lower 0.787044
    per-run rates         24 values in [0.0487, 1.0000]
    min per-run CP-lo     0.043054      (pessimistic diagnostic, NOT charged)
    CLUSTER lower bound   0.495525      <-- CHARGED
    old shared-report     0.048849      selection bias 0.739897 (charged as a row)

**Why the cluster bound, not a pooled binomial.** The class map is deterministic
and destination classes are shared, so the events are not exchangeable Bernoulli
draws; a pooled CP would be decoration, which is what the judge called v33's. The
run is the cluster and the estimand is the run-LEVEL MEAN (production is one draw
from the same seed law, not the best or the worst run). The charged value is the
normal-approximation lower bound on the mean of the 24 per-run rates
(0.638049 - 1.96*se = 0.495525). The minimum per-run CP-lo is printed beside it so
the pessimistic end of the same measurement is visible.

## 3. mu is derived, not imported (objection 3)

v33 charged mu = 11/10 as an imported premise. v34 measures it on the same traces:

    mu = completed / (launches * theta) = 1.906489

(mean over the 24 runs; `launches` and `completed` are the shipped counters). An
audit perturbs the launch count by +10 percent and shows the row move by exactly
1/1.1, so the number is a function of the trace, not a constant. It replaces 1.1 in
the sizing law, and it is larger, i.e. a debit. The traces also expose an honest
wrinkle worth stating: 22 runs halt at the step budget with `completed == 20000`
and 2 halt EARLIER at the launch bound (completed 17912 and 17932) once `launches`
reaches `ceil(S/theta) = 625` and the population drains. So `completed == S` is not
a universal of the protocol; `completed <= S` is, and equality holds on every run
that the budget itself stops. Both are recorded per run in `evidence/tail34.json`.

## 4. The two constants are derived (accounting objection)

**FINAL_EXP.** v33 asserted 300 Fp12 multiplications with no chain. v34 counts the
chain from the shipped curve parameter:

    easy part   f^((p^6-1)(p^2+1)): 1 Fp12 inversion + 2 Fp12 muls
                (the p^2 Frobenius is conjugation in the tower)
    inversion   54 + 36 + 15 Fp muls down the tower + 1 Fp inversion
                (Fermat, counted on p-2: 380 squarings + hamming(p-2)-1 muls)
    hard part   3 exponentiations by |x| = 0xd201000000010000
                (64 bits, Hamming weight 6) -> 63 squarings + 5 muls each = 204
                plus the 12 Fp12 muls of the standard hard-part chain
    total       216 + 15.204 = 231.204 -> 232 Fp12 multiplications

The count is emitted by `ledger34.final_exp_chain()` and printed with its bit
length and Hamming weight, so it can be re-derived from the parameter alone.

**G2 / 54.** The divisor is an identity, not an assertion:

    Fp12 = Fp6[w]/(w^2 - v)   1 Fp12 mul = 3 Fp6 muls (Karatsuba)
    Fp6  = Fp2[v]/(v^3 - xi)  1 Fp6  mul = 6 Fp2 muls
    Fp2  = Fp[u]/(u^2 - beta) 1 Fp2  mul = 3 Fp  muls
    => 3 * 6 * 3 = 54 Fp muls per Fp12 mul, printed as `fp12_tower`

and the additions and loop control of the same G2 chain are charged in
INSTRUCTIONS at the organiser prices instead of disappearing inside the division:
7112 Fp additions (254 DBL x 16 + 127 ADD x 24) x 20 = 142240, plus 255 loop
iterations x 32 = 8160, plus 64 x 4 x 20 = 5120 inside the Miller loop. The
pairing row is 75271076 instructions = 2^26.165592, printed with the chain.

## 5. Rows (instructions, and target ops after the one division)

    row                 log2 instructions    what it is
    w_effective         2^141.770805         S * step price, S = 139528627497097697523360967773670271203
    abandon_row         2^50.498678          cycle-abandon charge
    aux_sum             2^71.690182          stores, starts, solve, verify
    pairing_row         2^26.165592          3755.8 Fp12 mul = Miller 1533 + fexp 232 + G2 225.8 + verify 1765
    overflow_row        0                    reports capped at ceil(S/theta) literally
    stop_loss_row       2^109.056976         in-flight population + abandoned work

    W = 4756251854796248714126699096368027969999584 instructions
    W_ref = W / C_REF = 2^127.483093          <-- THE ONE DIVISION LINE
    claim = least 4-decimal c with 2^c >= W_ref = 127.4831

S solves `(1-delta)(1-exp(-S^2(1-f^2)p_det/(2 mu N))) >= 0.392` with delta = 0,
p_det = 0.495525, mu = 1.906489, f = 3/8, N = (r-1)/12. Reports =
2*ceil(S/theta) = 230830751123466.

## 6. Scope, sensitivity, and what is NOT claimed

* The fixture is 785 classes at theta = 2^5, so theta ~ sqrt(N) there; production is
  theta = 2^80 against sqrt(N) = 2^125.6. The fixture measures the MECHANISM and the
  per-trail lifetime distribution, not the production birthday regime. No production
  run and no production mixing theorem is claimed.
* Both measured penalties are charged at their FIXTURE values, which over-state the
  production values in the same direction, so the claim is conservative there. The
  sensitivity row, with v33's imported premises (delta = 0, mu = 1.1, p_det = the
  pooled conditioned 0.997162), gives **126.582**; the honest charge of this
  measurement gives **127.4831**. The 0.90-bit difference is exactly the price of
  the three repairs: mu 1.1 -> 1.906489 and p_det 0.997162 -> 0.495525 (the
  merge-all cluster lower bound).
* The claim is behind this account's other g1-invert work and the live board
  minimum; it is submitted because the repaired measurement is the point.
* `evidence/` is regenerated by `scripts/measure_v34.py` (measured rows),
  `scripts/finalize34.py` (the cluster charges), `scripts/regen_v34.py --write`
  (ledger + certificate) and `scripts/audit34.py` (the four controls). `--check`
  re-emits and agrees; `--control-ledger` runs the one-division control.
