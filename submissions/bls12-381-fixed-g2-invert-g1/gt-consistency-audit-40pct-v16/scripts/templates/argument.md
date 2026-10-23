# gt-consistency-audit-40pct-v16 — argument

Claim: log2 of total classical work is {claim} = {claim_num}/{claim_den} in C = 20000 units, against a
gate of {gate} ({gate_sub}, state {gate_state} at submission). Success floor 0.39 on every run, after
every debit. One meter convention, one exact cap, one machine table for every cited number.

## 1. Reduction and algorithm (unchanged from v15, not disputed in structure)

For the supplied nonidentity Q in G2, g = e(P0,Q) generates GT: Q = [b]Q0 with b != 0, so
g = e(P0,Q0)^b and e(P0,Q0) has order r. Uniform z in GT is g^x for a well-defined x mod r. Returning
[x]P0 gives e([x]P0,Q) = e(P0,Q)^x = z. The walk works on Frobenius classes: the class invariant is
T = Re(c0), a class function whose low bits carry the branch, the multiplier index and the
distinguished-point predicate in disjoint ranges; the twelve class keys are six ± pairs computed with
four Fp multiplications; the trace is free. A step is a cyclotomic squaring with probability 3/8, else a
dense multiplication by a table entry from one of {rtab} precomputed multipliers. Exponent pairs
(A,D) with y = g^A z^D are carried through both step types, so a class collision yields a linear system
in x; a degenerate or failing solve is discarded and never counted as a found collision.

## 2. Defect D1 — the success reserve, removed and re-derived

v15 solved S from the birthday expression for 0.39 and then inserted 43/100 before subtracting losses.
No argument supplied the additional probability, and its premise H1' bounded collision probability from
ABOVE, which cannot produce a success lower bound.

v16 states the premise as the no-collision upper tail:

  H1'': for distinct steps i != j of the walk, P(y_i = y_j up to the Frobenius class) >= rho/N, hence
        P(no collision among S steps) <= exp(-S(S-1)/(2N')),  N' = N/rho, rho = {rho} = {rho_float}.

The exponent is the U-corrected pair count S(S-1)/2, not S^2/2: the difference is 2^-125 relative and
costs nothing, but it is the exact expression. From this upper tail the collision probability is bounded
BELOW by 1 - exp(-S(S-1)/(2N')). rho is the conservative side of the reviewer's merge-rate expression:
rho_ours = (1-s^2)(1-1/R) = {rho}, rho_reviewer = {rho_review}, and ours is smaller by
2s(1-s)/R exactly.

S is not asserted; it is solved, with Fractions, from the requirement that the floor survives every
debit. The iteration is: target = 0.39 + abort + detection + storage debits; ln_hi(1/(1-target)) by a
positive-term atanh series; S = minimal integer with S(S-1) >= 2N' ln_hi. At the submitted S the
check is exact:

  S(S-1)/(2N') = {birth_v}  >=  ln(1/(1-{target})) ,  target = 0.39 + {debits}
  =>  1 - exp(-S(S-1)/(2N')) >= {target}
  =>  floor after subtracting {debits} = 0.39 exactly  (certified value {birth_floor}).

S = {s_steps} steps, log2 S = {s_log2}. No probability is invented anywhere; the reserve is the
debits, each of which is separately bounded.

## 3. Defect D2 — replay work counted, reserved, and tested

v15 counted S forward steps and ignored the trail replays that a match requires; its production test
retained (A,D,y) so no replay ever happened, and it checked the cap after executing a step.

v16 charges the replays as a row. A trail is one of {rtab} walks of length {trail_len} steps
(log2 = {trail_len_log2}). A match ends the run and needs two replays, each re-executing a trail from its
recorded start to the matched step. The row is {replay_steps} steps at the worst-case multiply price,
{replay_instr} instructions; false-match replays are bounded below one instruction (a 256-bit
fingerprint collision over 2^25.7 records) and appear as an explicit zero row. The replay reserve is
carved out of the cap, so the every-run bound covers replays: step cap = cap - replay =
{cap_instr} - {replay_instr}.

The production test now reconstructs state only from 52-byte slots (32 B fingerprint, 2 B trail id,
16 B step index, 2 B spare): {replay_trails} trails x {replay_seglen} steps = {replay_steps_run} steps,
{replay_records} records, {replay_forced}/{replay_forced} forced replays verified, one end-to-end
recovery where the recovered log satisfies g^x = z ({replay_recovered_ok}), and zero false usable
collisions. The meter check is pre-execution in both the specification and the source: the charge of the
next step is computed from the current state and meter + charge <= cap is tested BEFORE the step runs; a
cap admitting exactly {cap_bound_steps} steps stops at {cap_bound_steps} with one clean abort.

## 4. Defect D3 — storage as a maximum with an explicit debit

A distinguished-point probability of 2^-100 gives an EXPECTED record count, not a per-run maximum, and a
work cap bounds steps, not DP density. v16 therefore states a maximum and a policy:

  M = ceil(S/2^100) + {store_slack} = {rec_cap} records, {rec_bytes} bytes per slot,
  store = {store_bytes} bytes; on a full cap the algorithm refuses the new record and charges the
  success debit P(records > M) <= exp(-2^{store_debit_exp}) from the Chernoff bound for
  Binomial({s_steps}, 2^-100) at delta = {store_slack}/2^25.7.

That debit is inside the floor target of Section 2, so it is paid for rather than assumed away. Peak
memory is store + table ({table_bytes} B) + code + scratch + sort buffer = {peak_bytes} bytes
({peak_gb} GB) < 2^32.

## 5. Defect D4 — one meter convention, one exact cap

v15's source charged full prices while its ledger asserted an overhead-only convention, so the numbers
differed by exactly {v15_mul_count} x 20000 and the PASS rows did not validate the implementation; its manifest also
rounded the cap downward.

v16 uses ONE convention: a multiply step costs {mul_charge} = 20000 + {i_mul} instructions, a squaring
step {sqr_charge} = {i_sqr}. The executed run records a per-event charge log of {prod_steps} entries;
the certificate proves exact integer equality of three quantities:

  charge-log sum {prod_log_sum}  ==  reported meter {prod_meter}  ==  row formula {prod_formula}.

No range, no rounding, no convention switch. The cap is the single exact integer
floor(2^{cap_log2}) = {cap_units} work units = {cap_instr} instructions, used identically by the
manifest, the ledger, the certificate and the source's metered loop. Nominal work
{w_units_num}/{w_units_den} units (log2 {w_lo}..{w_hi}) is strictly below it, and the cap is strictly
below the claim.

## 6. Defect D5 — census re-tabulated

The toy census is regenerated with the two channels separated and the match/degeneracy distinction made
explicit:

  matches {toy_matches} = same-index {toy_ch_same} + shifted-image {toy_ch_shift} + unclassified {toy_other_failed}
  same-index: {toy_ch_same} matches, {toy_rec_same} recover, {toy_same_failed} fail
  shifted-image: {toy_ch_shift} matches, {toy_rec_shift} recover, 0 fail
  degeneracies {toy_deg} = {toy_same_failed} + 0 + {toy_other_failed}

v15 reported {toy_deg} as "same-index degenerate cycle closures"; {toy_deg} is in fact the total
degeneracy count and the same-index failure count is {toy_same_failed}. The shift channel is fully
usable ({toy_rec_shift}/{toy_ch_shift}), which is the evidence for charging the shifted-image gap in H1''.

## 7. Work rows

  setup            {setup_instr}   {rtab}-entry table: 510 dense mults + 12 Frobenius each, plus load
  walk steps       {steps_instr_num}/{steps_instr_den}   S steps at 5/8 multiply and 3/8 squaring
  replays          {replay_instr}   two trail replays at the multiply price
  false matches    0               bounded below one instruction
  total            {w_units_num}/{w_units_den} units = log2 ({w_lo}, {w_hi})

## 8. Support and its scope

Production field: {prod_steps} steps, {prod_sq} squarings ({prod_sq_pct}% ± {prod_ci}),
{inv_checks} invariant checks with {part_a_total} algebraic predicates PASS ({part_a_pass}/{part_a_total}),
{class_tests} class well-definedness tests, predicate agreement {prod_checked}/{prod_checked} against an
independently written rule. Replay pass and cap-boundary test as in Section 3. Toy: {toy_trials} trials,
{toy_rec} recovered at the design budget, {toy_deg} degeneracies of which {toy_other_failed} are
unclassified fingerprint matches, {toy_escapes} escapes.

These support the mechanism at their measured scope. They are not a production success measurement, and
the toy partitions by hashed canonical representative rather than by production trace bits.

## 9. Heuristics and objections

H1'' (pair-collision rate >= rho/N with rho = {rho}) and H2 (trace-bit pseudorandomness) remain
heuristics, as does the random-mapping idealisation behind the birthday tail. The extrapolation from the
toy field to the 381-bit field is material. What is NOT heuristic in this package is the bookkeeping:
the floor is computed at the submitted S from the stated expression with every debit subtracted, the
replays and the storage maximum are charged rows, the meter convention is single and the equality is
exact integers, and every count in this document is substituted from evidence/numbers16.json by
scripts/machine_tables.py.

## 10. Revision history

- v15 `{v15_id}` rejected with five internal-consistency defects: unsupported 0.43 reserve (D1),
  uncounted trail replays plus a post-hoc meter check (D2), storage maximum taken as an expectation (D3),
  meter-convention split and rounded cap (D4), census miscount (D5).
- v16: D1 re-derived from the birthday tail at S with all debits; D2 replay row reserved inside the cap
  and tested through 52-byte records; D3 maximum-record cap with a Chernoff debit; D4 one convention with
  an exact integer equality row and a single exact cap; D5 census re-tabulated. Claim unchanged at
  {claim}.
