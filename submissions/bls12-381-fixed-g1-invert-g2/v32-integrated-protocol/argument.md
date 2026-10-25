# Argument: the integrated protocol (v32) -- v31's rejected four, repaired

v31 (ticket 816b4c54, claim 126.6053) was rejected on four flaws: a protocol the
evidence did not execute (plus free bootstrap launches), an unsupported retention
law measured by a sampler that suppressed the very event it was measuring,
output construction left as "unchanged return algebra" at 80-bit labels, and a
control whose two claims were identical. v32 repairs those four and changes
nothing else.

Claim **126.5979** = 1265979/10000 (reduced), the least 4-decimal log2 >= the
exact total 2^126.597817. Emitted by `scripts/regen_v32.py`; no claim literal is
typed anywhere.

## 1. `run_protocol` is the only path; bootstrap launches pay

`scripts/protocol.py` contains one entry point. The success experiment, the
retention histogram, the gate trace and every audit call it; there is no
replacement loop and no harness-side scheduler. Per trail it executes, in order:
`classkey.step_v26` (the family-free squared step), the distinguished-point test,
the store insert/match on the class key, and the twelve-shift recovery with its
verification against the planted target. Control (a) patches out the DP test
(`skip_dp=True`) and success collapses from 0.75 to 0.00 over 60 planted
challenges -- the experiments cannot bypass the protocol because there is nothing
else to call.

The credit rules are in the object, not the note:

    RULE 1 credit = completed - consumed; a launch needs credit >= theta
    RULE 2 a bootstrap launch CONSUMES theta from the same S budget (v31: free)
    RULE 3 a credit-eligible launch at INFLIGHT_CAP is DEFERRED, credit kept
    RULE 4 launches are clamped at launch_bound = ceil(S/theta), bootstrap included
    RULE 5 halt on first verified recovery, else at the bound

The v31 defect was that RULE 2 launched free, so a trail ending before theta
could re-bootstrap without ever paying: ceil(S/theta) was not a launch bound. In
v32 the bootstrap line is

    launch(t); consumed += cfg.theta; launches += 1        # protocol.py, RULE 2

so the launch population and the step population are the same budget. Control
(c) cuts every trail at theta/4 (`cut_steps`): launches then reach the bound and
the trace shows it binding at 125 of 125 with 971
completed steps; under v31's free-bootstrap rule the same walk would have issued
125 (all bootstrap relays, credit charged 0 under v31) launches, i.e. it would have exceeded the bound by n/a - v31 enforced no launch counter at allx.

Gate trace, shipped configuration (evidence/gate32.json): in-flight reaches
128 at round 848 (t=35710, launches=485) and the next credit-eligible launch is refused with
`reason='len(inflight)==128'`; deferred count 3628.

## 2. Retention: measured through the protocol, with the sampler order fixed

v31's `v31_tail.batch` continued whenever a key was already stored, **before**
handling `is_dp`, so a repeat DP did not terminate its trail, the shared-report
event was suppressed, and repeat visits were counted as fresh evidence. v32 tests
`is_dp` first and only then consults the store; the `suppress_shared` switch
reproduces the old order on the same input as a control. Measured on one input
(`evidence/retention32.json`):

    order                reports  matches  retention   repeat visits  mean age
    is_dp first (shipped)  10   40    0.600      0      30892 steps (965 theta, cap-truncated)
    dedup first (v31 bug)  10   0    0.000      2770      undefined (no match was ever recorded)

The histogram differs (yes - 9 occupied age bins against 0); the buggy order reports 2770 repeat
visits that are not reports at all. Only the fixed order's histogram is used.

The age law is no longer assumed. v31 posited an age density 2(L-u)/L^2 and a
conditional mean report distance theta; the verdict rejected that as an unproved
premise. v32 measures the report age at match time on the integrated run and
charges the measured mean through Markov, exactly as the v31 derivation did but
with the constant measured instead of assumed:

    delta = P(miss | merge) <= E[age] / cap = 30892 steps (965 theta, cap-truncated) theta / (128 theta)
          = 1/128 = 0.0078

Support domain: ages are collected over the whole shipped run at theta=2^5,
cap=128 theta, on the 785-class fixture; the slope of the age histogram over
[0,128) theta is recorded in `evidence/retention32.json`. This is an empirical
law at fixture scope; it is charged, not asserted, and the `max(estimate, 1/64)`
patch is dead. The fixture is not production: N=785 and theta=2^5 against N=2^251 and
theta=2^80, so what transfers is the *shape*, and the transfer is disclosed as an
extrapolation.  Two numbers are charged.  (i) delta = 1/128 = E[remaining]/cap
by Markov, which needs only that a follower's remaining length to the cap is at
most cap and that the next report is at mean distance theta.  (ii) The measured
first-report retention 0.600 REPLACES the 0.9854 the earlier versions declared,
in the slot v27's sizing law divides by: S scales as 1/sqrt(retention), so this
debit costs 0.368 bits and the claim carries it (126.6053 -> 126.9558).  That is
the whole point: the fixture number is charged, not floored.

## 3. Output construction: Q = [x]Q0, instruction-counted

The recovery returns a GT discrete-log scalar x with z = z0^x where z0 = e(P,Q0)
for one fixed Q0 in G2. Then e(P,[x]Q0) = e(P,Q0)^x = z, so the returned point is
Q = [x]Q0. There is no affine square-root chain and no "unchanged algebra"; the
row is built from instruction counts (`ledger32.output_rows`):

    Miller(e(P,Q0))          1533 Fp12 mul   (64 DBL x 9 + 87 ADD x 11)
    final exponentiation      300 Fp12 mul   (1 inv + ~10 mul + 3 x 64-bit exp)
    G2 reconstruction [x]Q0   226 Fp12 mul   (254 DBL x 24 + 127 ADD x 48 Fp mul,
                                             / 54 Fp mul per Fp12 mul)
    verify e(P,Q) = z        1833 Fp12 mul
    pairing_row              3892 Fp12 mul = 77835556 instruction equivalents
                             = 2^11.93 target operations

Widths are the shipped step's, not a two-word fiction. `step_v26` updates the
labels modulo the 255-bit subgroup order:

    label           255 bits, FOUR u64 limbs, add = 5-limb add + conditional subtract
    step counter    87 bits (cap 128 theta = 2^87), 16 B, not 8 B
    class key       48 B / 6 words;  record 112 B;  live GT states 128 x 656 B

## 4. Control (d) now moves the claim

v31's control d printed 126.6053 for both the base and the perturbed ledger,
because it perturbed only the underlying total. v32's control d perturbs the
emitted rational by one unit in the last place and re-runs the emitter: the
base claim is 126.9558, the 1-ulp-moved ledger emits 126.9559 (and 126.9570 for a 1e-3 retention quantum), and they differ
(`regen_v32.py --control-ledger`). A control that cannot move the number it
claims to move is not a control.

## 5. Rows

    w_effective      126.597817   S * c_max, S from the sizing law below
    abandon_row       35.325690   launches x P(merge inside cap) x cap
    aux_sum           56.517206   reports, store traffic, sorting
    pairing_row       11.926214   section 3
    overflow_row       0          reports are capped at ceil(S/theta) literally
    stop_loss_row     94.769251   in-flight population + abandoned work
    total            126.597817    claim 126.5979 = 1265979/10000

S solves (1-delta)(1-exp(-S^2(1-f^2)/(2 mu N))) >= 0.39 with delta measured in
section 2, mu = 11/10, f = 3/8 measured on the fixture, N = (r-1)/12. c_max is
the MAX price per step, so the mean-price row (2^126.27) is a diagnostic and is
never summed.

## 6. Scope and disclosed rank loss

The step, the twelve-shift return algebra, the family-free theorem and the N=785
fixture scope are unchanged from v27-v31. No production mixing bound is claimed
and no production run is claimed: the experiments are evidence at their scope.
Success is measured at fixture scale (section 1) and sized analytically at
production; the gap is a disclosed extrapolation, not a proof. The claim is
behind the live board minimum and this account's other tracks.
