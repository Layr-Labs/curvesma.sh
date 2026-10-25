# g1-invert v33 verdict -- ticket 90701ee2 (REJECTED, unscored/invalid)

Ladder was clean for the first time in this lane: commit 98302d9 -> Challenge tests
OK (run 40151590345) -> Submissions add-only OK (40151589938) -> workflow-prepare
intake OK -> validating -> initial experiments OK -> initial judge 2 m 8 s OK ->
deterministic finalization OK -> verdict. No intake kill.

The judge ACCEPTED two of the four repairs by name:

* "the revised cutoff bounds executed steps by S"                     (objection 1)
* "the displayed ledger's single division gives log2(W/20000)
   approximately 126.581950, whose four-decimal ceiling is 126.582"   (objection 4)

Three objections survive, plus an accounting objection.

## Full rejectionReason, verbatim

The complete claim is not supported under exploratory-v1. The reduction is sound when z0=e(P,Q0), the revised cutoff bounds executed steps by S, and the displayed ledger's single division gives log2(W/20000) approximately 126.581950, whose four-decimal ceiling is 126.582. These repairs do not establish the production success or complete resource claims.

1. The claimed exact geometric tail has false premises. routing.dp_hit is a deterministic function of the class key, not fresh independent randomness at each step; revisiting a key repeats its DP result. Even granting an idealized hashing heuristic, R=min(cap,S-t_launch) is remaining aggregate work, not the number of steps allocated to that trail: concurrent trails share the remaining budget. Moreover, run_protocol launches on each eligible round, not once per theta executed steps. Its explicit launch-count cap does not establish the launch spacing used in the geometric-sum argument. Averaging these launch-time q^R values therefore does not bound the claimed conditional miss probability after a merge. The fixture's computed tail values and their monotonicity under rescaling R do not repair these event and scheduling errors.

2. The detection experiment conditions on an already observed shared report. In run_protocol, merged is incremented only after a DP hit and a successful lookup of a previously stored key. Thus 1054/1054 measures successful scalar recovery conditional on a shared report, not the probability that an arbitrary trajectory merge reaches a detectable report before termination. Undetected trajectory merges never enter its denominator. Pooling these matches from seven runs also has no supplied binomial-independence justification. Computing a Clopper-Pearson number cannot correct the selection problem or establish its applicability at production parameters.

3. The production birthday formula, including its quantitative mu=1.1 allowance, lacks relevant support for this specific walk. The strongest supplied response is the acknowledged limitation that the degenerate 785-class fixture tests mechanisms rather than production collision behavior. That is an appropriate limitation, but it leaves the required production success probability of at least 0.39 unsupported; charging the fixture's CP lower bound does not bridge that gap. This objection does not demand an unconditional mixing theorem or a full-scale attack run.

Separately, regenerating the ledger establishes arithmetic consistency, not complete instruction accounting. For example, FINAL_EXP=300 is asserted without a supplied exponentiation chain or operation derivation, and converting G2 field-multiplication counts through a factor of 54 does not justify the complete scalar instruction cost of its arithmetic, additions and control. These accounting gaps remain unresolved, while the success and tail defects already independently prevent acceptance of the original claim as written.

runtime 128.4 s, 1 attempt, gpt-6-astra, 6095 completion tokens; evidence_refs
include submission.json, argument.md, scripts/protocol.py, scripts/routing.py,
scripts/measure_v33.py, evidence/detect33.json, evidence/tail33.json,
scripts/ledger33.py, evidence/ledger33.json.

## Surviving flaws (my read)

1. TAIL IS STILL THE WRONG RANDOM OBJECT (hard, and it is a fair hit). Three
   distinct defects: (i) dp_hit is a DETERMINISTIC function of the class key, so
   per-step independent Bernoulli(1/theta) is a modelling assumption, not a fact
   of the code; (ii) R = min(cap, S - t_launch) is the AGGREGATE remaining budget
   while other trails are still spending it, so R is not this trail's allocated
   lifetime; (iii) launches happen on every eligible round, not once per theta of
   executed work, so the geometric-sum spacing (and hence 0.582*theta/S) is not
   derived from the shipped schedule.
   FIX: stop deriving delta from a schedule model. Either (a) delete the delta
   factor entirely and make the protocol's stop rule the only lifetime statement
   (a trail halts at DP, at the per-trail cap, or at the global budget -- and the
   budget stop is not a per-merge event at all), or (b) measure the tail from the
   ACTUAL per-trail allocated lifetime: give each trail an explicit step quota
   q_i <= cap written at launch and charged from the budget at launch, then
   P(D > q_i) is a statement about that trail alone and the fixture measures it
   directly (Kaplan-Meier on trail lengths, per-quota rows), no i.i.d. claim
   needed beyond the DP test's hash behaviour on FRESH classes, and state the
   revisit caveat as the only remaining assumption with its measured rate.
   Also: if the launch schedule is what carries the bound, instrument it -- log
   every launch time and prove the spacing distribution, don't assert theta.

2. P(detect|merge) IS SELECTION-BIASED (hard, fair). Merges are only counted when
   a DP hit AND a stored key match both happened; a trajectory merge that never
   reaches a DP is invisible, so 1054/1054 is P(recover | shared report), not
   P(merge -> detectable before termination). Pooling 7 runs without an
   independence argument is a second defect, and a CP bound computed on a
   non-random sample is decoration.
   FIX: define the merge event WITHOUT the store: plant collisions (two trails
   started from labels constrained to meet at a known class, or a same-class
   pair drawn from the state space), then measure the fraction whose FOLLOW-ON
   walk reaches a DP (anywhere, budget permitting) and is recovered -- that
   denominator includes the misses. Report it per run with a run-level
   (cluster) interval, not a pooled binomial: n_runs >= 30 with per-run rates, or
   a hierarchical/Clopper-Pearson on the run-level mean. Charge the lower bound.

3. PRODUCTION SUCCESS STILL UNSUPPORTED (hard, structural). The judge says
   plainly it needs neither a mixing theorem nor a full-scale run -- but the
   degenerate fixture cannot carry 0.39, and mu=1.1 has no independent support.
   FIX: build a NON-degenerate fixture in the same protocol code (larger N at
   production-like theta/sqrt(N) ratio: theta << sqrt(N), which the v33 note
   itself concedes the 785-class fixture violates), measure the birthday
   constant and the walk-quality factor q there, and charge the measured q in
   place of mu=1.1. If no such fixture is possible in the budget, the honest move
   is to claim NO success probability at all: submit the cost bound as a
   conditional "if success >= 0.39 then W" and drop the mu term from the claim
   arithmetic (a claim that does not need the success premise cannot be rejected
   for lacking it) -- or lower the target to what the fixture supports and say so.

4. INSTRUCTION ACCOUNTING INCOMPLETE (medium, cheap): FINAL_EXP=300 has no
   supplied exponentiation chain, and the G2 row divides Fp-mul counts by 54
   without charging additions and control. FIX: give the final-exponentiation
   chain explicitly (easy part: 1 inversion + ~10 Fp12 mul; hard part: 3 x 64-bit
   exponentiations with the DBL/ADD counts written out) and price the G2 scalar
   multiplication as Fp2/Fp arithmetic with explicit add and carry rows instead
   of a single /54 division.

## What held

* Step budget: accepted. completed == S traced, consumed == executed, launches
  bounded by construction; the v32 rule measured at 32.70x on the same seed.
* Unit chain: accepted. One division, W = 2^140.869663 instructions,
  W/C_REF = 2^126.581950, claim 126.582 = 63291/500; control (d) shifts by
  log2(20000).
* Intake: title 163 ASCII chars, score {"numerator":"63291","denominator":"500"}
  as strings, 28 packed files / 158 KB, ladder exit 0 on every rung.
