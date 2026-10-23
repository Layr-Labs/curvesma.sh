# Argument: confirmed detector, replay records, node-wise debit (v9d)

Claim: total classical work <= 12939/100 = 2^129.39 target-group-operation
equivalents at C=1000, success >= 0.39, advice 0, peak memory < 2^73 B.
Live head at dispatch: 129.46 (`5fa78d62`). Supersedes v9c (`5e465e2`, rejected).

## 1. Target, normalization, algorithm

NIST P-256, p = 2^256-2^224+2^192+2^96-1, a = -3, cofactor 1, order r prime. C=1000
is the borrowed secp256k1 unit; P-256 arithmetic is charged by its own schedule
(section 2). H=0, so score = log2(T/1000) with T in instructions.

Construction inherited unchanged: negation classes {P,-P}, N = (r-1)/2 classes,
key = x, canonical member (x,|y|), transition canon(2X_i+T_j) with
j = H(X_i) mod 2^16 and T_j from a +/- twisted table, an escape that doubles the
minimum-fingerprint class of the detected cycle, recovery x = (a'-a)/(b-b') mod r
at a class collision, verified by [x]P = Q.

## 2. Complete step price (3087 instructions, re-itemised)

  two affine additions, core     4M+2S+12A                     1484
  batch-inversion products       2 slots x 3M                  1320
  batched Fermat inversion share 2 x ceil(51856/2^16)             2
  coefficient update             (2a+k, 2b) mod r                63
  canonicalisation of y          p-y, cmp, csel, store           14
  degeneracy tests               x2==x1 for both adds            12
  bounded detection              2 fp compares + ring + MIN-fp
                                 select + prime-67 boundary +
                                 running-min compare EVERY step  16
  coefficient negation           2 x (negate 16 + select/store 8) 48
  control (itemised, 6 rows)     hash/DP/route/cap/load          128
  total                                                         3087

The detector row is 16 instructions with the running-minimum comparison charged on
EVERY iteration (v9c amortised it by 1/64; the reviewer was right that attack()
compares every step). Six named lines, no amortisation, no division by 64. The
window is prime (67), so no per-window division is charged at all.

## 3. Ledger (P4: one constant, imported)

KAPPA = 11/10 in `experiments/protocol.py`, imported by toy and ledger:
S_ops = ceil(KAPPA sqrt(2 N' ln(1/0.60))) = 2^127.652966, budget ratio 1.111844.

  phase 1      = 3.087 (S_ops + 2^60)           = 2^129.279171
  reports      = 2^63.652966 x 128 instructions = 2^60.687
  confirm row  = 2^63.652966 matches x 24       = 2^58.272
  replay row   = 2^10 x 2^63 x 3087             = 2^74.626
  setup        = 2^16 table entries + node init = 2^25.045
  total                                        = 2^129.279171
  stress (every row x 6/5)                     = 2^129.542206
  claim                                        = 129.39

New rows, both itemised and both inside the total:

* **Confirm row (O1a).** Every fingerprint match is confirmed by FULL POINT
  EQUALITY (x equality and y equality, 12 instructions each = 24). At DP density
  2^-64 over 2^127.652966 steps there are 2^63.652966 fingerprint matches, so the
  confirm row is 2^68.2 instructions = 2^58.272 units, 2^-71 of phase 1. A false
  match therefore costs a confirm and is ignored; it can never trigger an escape.
* **Replay row (O2).** Endpoint records hold no coefficients (section 7). On a
  match the owner trail is replayed from its start: expected replay length is
  trail_mean/2 = 2^63 steps at the step price 3087. Match events follow the
  birthday law, reports^2/(2N) = 0.41 expected; the row is charged at the
  conservative bound 2^10 events, i.e. 2^84.6 instructions = 2^74.626 units,
  2^-54.7 of phase 1.

Escapes are arrivals: a doubling of at most 1620 instructions consumes a step
slot priced at 3087, so every escape is charged inside S_ops.

## 4. Reports, capacity, routing

DP width 64, so expected reports = 2^63.652966 <= 2^64 = cap/4. Capacity is
2^16 nodes x 2^50 records x 26 B = 2^70.70 B < 2^73 B (section 7).

Field spec unchanged from v9c: key K; the DP predicate at width b tests
mix64((K mod 2^64) ^ DPSALT) mod 2^b == 0 and reads ONLY bits 0..63; the router
computes mix64(((K >> 128) mod 2^64) ^ RSALT) >> 48 and reads ONLY bits
128..191. Production keys are the 256-bit x-coordinate. The toy now uses 128-bit
class keys, so the fingerprint (64 bits) is strictly shorter than the key, as in
production, and a fingerprint match no longer implies point equality.

Evidence at production width (phase G, unchanged): 150000 DP-conditional 256-bit
keys reach 58884 of 65536 nodes, max load 11 against expected 2.2888, union bound
2^-30.74 at tolerance 24; 1.5M keys attain all 65536 route values; the v9b spec
reaches 1 node with load 150000 on the same keys.

## 5. Detector (O1): confirm, minimum fingerprint, prime window

Three defects were named; each is repaired in code and tested.

(a) **Full point equality.** `points_equal(X, Y)` compares x and y. Both escape
branches call it; a fingerprint match that fails the confirm increments
`false_matches` and the walk continues. Charged in the ledger confirm row.

(b) **Minimum-fingerprint select.** The 2-cycle branch takes
`m = r1 if (fp1, key1) <= (fp2, key2) else r2`: the minimum fingerprint of the
cycle, tie-broken by the full key, exactly as specified - not the x-coordinate.

(c) **Prime window, randomised phase.** The checkpoint window is 67 steps, not
64, and its phase is drawn per trail from the start-state hash
(`phase = mix64(key_start ^ const) mod 67`). Because gcd(67, L) = 1 for every
L <= 66, the phase-shifted boundaries sweep all residues mod L, so every short
cycle's minimum is captured by some window and its next recurrence matches the
stored running minimum. Fixed 64-step alignment is the trap the reviewer named:
for L | 64 the boundary offset never advances.

(d) **Toy keys exceed 64 bits.** Toy class keys are 128 bits: the low field is
mix64(x) (read only by the DP predicate), the high field is mix64(y_canon), which
carries the canonical sign/parity (read only by the router). Toy state space is
2^26, so the toy false-match expectation is 2^-64 x 2^26 = 2^-38 (none occur),
while production is 2^-64 x 2^127.65 = 2^63.65 confirms; every production confirm
is charged. The mechanism itself is exercised deterministically by the
false-match probe (phase H): equal 64-bit fingerprint on different points is
rejected and counted, a genuine repeat is confirmed and escapes.

**Exhaustive adversarial test (phase H).** All L in 2..70 x all 67 phases:
69 x 67 = 4623 cases, and the shipped rule escapes in every one. The control
implements the fixed-64-window boundary-only rule and traps all 64 phases of each
L in {2,4,8,16,32,64}; the prime-67 rule traps zero of them.

## 6. Fixed-cap refusal and the node-wise debit (O3)

Each node cap is fixed BEFORE the run from the DECLARED expected report count:
cap = max(1, ceil(5.087772 x declared reports / nodes)). A report finding its
node full is REFUSED, counted, flags the trial cap-loss and is not stored,
matched or overwritten. Every pipeline cell runs that path with 0 cap violations.

**Debit arithmetic, corrected.** The event that loses a trial is ONE node
reaching its cap, not the total report count reaching 2^66. Per-node mean is
2^63.652966/2^16 = 2^47.652966 and the cap is 2^50 = 5.087772x the mean
(delta = 4.087772). Chernoff:
  per node      log2 P <= -1.337456813887183e15 = -2^50.248414
  union x 2^16  log2 P <= -1.337456813887183e15 + 16 = -(2^50.248414 - 16)
This is the number used in the cap-loss and refusal rows: P_cap = P_refuse
<= 2^-1.337e15. The global 2^66-report bound (log2 <= -8.765e19) is retained only
as a diagnostic; it does not bound the single-node event.

## 7. Endpoint record and replay (O2)

The 64-byte record of v9c could not hold a class key plus two 256-bit
coefficients (88 B before indexing). Records are now **26 B**:
route16 | dp-fp64 | owner-trail-id64 | step64 = 208 bits.

A 40-bit trail id cannot name the 2^63.653 trails a 2^127.653-step run produces
(40 + 64 = 104 address bits < 127.653), so the id is 64 bits; that is the only
change from the requested 24 B layout, and it is charged. Trail starts are
derived deterministically, `a0 = PRF1(trail_id)`, `b0 = PRF2(trail_id)`, so no
start table is stored. Trail length is capped at 2^64-1 steps (restart as an
arrival), so the step field cannot overflow.

On a class-key match the owner trail is replayed from its start to the recorded
step, reconstructing (X,a,b) exactly (the transition is deterministic given X),
then the collision is solved and verified by [x]P = Q. The toy executes exactly
this path: per-attack start table, replay on every match, `replays` and
`replay_steps` counted, 0 verification mismatches. Replay cost is charged in the
ledger replay row.

Memory: 2^66 records x 26 B = 2^70.700 B; detector state 2^16 trails x 4 entries
x 136 B = 2^25.09 B; peak 2^70.70 B < 2^73 B.

## 8. Success arithmetic with named debits

success >= 0.39 - P_halt - P_cap - P_refuse, with
  P_cap    <= 2^-1.337456813887183e15       (node-wise Chernoff, section 6)
  P_refuse <= P_cap                         (same event: a full node)
  P_halt   <= extras/headroom = 2^-64.95    (Markov; extras = report + confirm +
            replay + setup rows, headroom = 2^129.39 - 2^129.279171)
No successful trial is truncated: trial cost is a hard function of the step
counter. Modelled success at S_ops is 1-0.60^1.21 = 0.4610; the toy measures walk
efficiency 0.86-0.92 in the refusal-on production cells, giving 0.42-0.44; claim
0.39.

## 9. Toy (production protocol end to end)

q22/q26, 400 attacks per cell (3200 in the volume cell), fresh challenge and table
per attack; 128-bit class keys, prime-67 randomised-phase detector with full-point
confirms, minimum-fingerprint select, replay-based recovery, fixed-cap refusal.
Refusal ON / same seeds OFF reported for every cell; adverse cells (dp9, half-mean
cap, 16-node saturation) kept as declared stress. Recovery at the production budget
ratio is reported for each cell; the q26 dp6 refusal-on cell is the headline
statistic. No verification mismatch occurs in any cell.

## 10. Corrections and limits

v9d repairs, against `5e465e2`: (O1) full-point confirm charged, minimum-fingerprint
select, prime-67 randomised phase, 128-bit toy keys, exhaustive short-cycle test
with the fixed-window trap control; (O2) every-step running-min row, 26 B endpoint
record, replay reconstruction charged; (O3) node-wise union bound in the cap-loss
row. Cost: none to the claim (new rows are 2^-54.7 and 2^-71 of phase 1); claim
stays 129.39 and step stays 3087. ord_r(2) = N is even, so 2^(N/2) = -1 and the
shortest sign-ambiguous gap is N/2 = 2^126.414 x S_ops. H1 carries the success
bound and is supported by the toy census, not proved. The replay bound (2^10
events) is 2^11 above the birthday expectation. Not a production run, not a
measured kernel, not a practical attack.
