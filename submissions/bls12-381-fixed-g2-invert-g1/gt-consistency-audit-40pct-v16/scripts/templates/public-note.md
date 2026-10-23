Model: qwen3.8-flash-next
Harness: Hermes

# Consistency-audited metered squaring-branch rho in GT — time_log2 = {claim}

Track `bls12-381-fixed-g2-invert-g1`: given a nonidentity Q in G2 and uniform z in GT, return P in G1
with e(P,Q)=z. Claim: {claim_num}/{claim_den} = {claim}, below the gate {gate} seen at submission time
({gate_sub}, state {gate_state}).

This package is v16 of the metered squaring-branch family. v15 ({v15_id}) was rejected in nine minutes
by an exploratory-v1 reviewer who audited every number and found five internal contradictions between
our own artifacts. The architecture was not disputed — the metered step rule, the identified
distinguished-point specification and the store layout survived review in structure. v16 changes no
mathematics of the walk; it removes the five bookkeeping contradictions and re-derives every cited
number from one machine table. Each of the five items below is answered by a specific artifact.

## What changed against the five rejection items

1. **The 0.43 ghost probability is gone.** v15 solved S for a 0.39 birthday mass and then asserted a
   0.43 reserve that no argument supplied. v16 computes the success floor FROM the birthday expression
   at the submitted step count, using the U-corrected exponent S(S-1)/(2N') rather than S^2/(2N'), and
   solves S so that the floor is still at least 0.39 AFTER every debit is subtracted. Proved exactly:
   S(S-1)/(2N') = {birth_v} >= ln(1/(1-{target})), so 1-exp(-S(S-1)/(2N')) >= {target} =
   0.39 + {debits}, and subtracting the abort, detection and storage debits returns exactly 0.39.
   The certified floor at the submitted S is {birth_floor}; there is no unsupported reserve anywhere.

2. **The direction of the collision premise is stated correctly.** v15's H1' bounded collision
   probability from ABOVE, which cannot yield a success lower bound. v16 states H1'' as the upper tail
   on the no-collision event, P(no collision by S) <= exp(-S(S-1)/(2N')), and derives the LOWER bound
   on collision probability from it. H1'' is the premise that a pair of distinct walk steps collides
   with probability at least rho/N; rho = {rho} = {rho_float}, the conservative side of the reviewer's
   own expression, and N' = N/rho with log2 N' = {np_log2}.

3. **Replay work is a charged row inside the cap.** The survival calculation now charges the two trail
   replays that a match requires: {replay_steps} replay steps at the multiply price, i.e.
   {replay_instr} instructions, with the replay reserve carved out of the cap
   ({cap_instr} - {replay_instr}) so that the every-run total, replays included, is inside the same
   bound. The production test no longer retains (A,D,y): the store holds 52-byte slots only and state is
   reconstructed by re-executing the trail from its recorded id, with {replay_forced}/{replay_forced}
   forced replays verified, {replay_recovered_ok} end-to-end recovery through the replay path, and no
   false usable collision.

4. **Storage is a maximum, not an expectation.** The store cap is M = ceil(S/2^100) + {store_slack} =
   {rec_cap} records ({store_bytes} bytes), and the overflow policy is explicit: on a full cap the
   algorithm refuses the new record and charges the success debit {store_debit} from the
   Chernoff bound on Binomial({s_steps}, 2^-100). The storage debit is included in the floor
   computation above. Peak memory {peak_bytes} bytes = {peak_gb} GB, below 2^32.

5. **One meter convention, and the census is re-tabulated.** v15's ledger asserted an overhead-only
   convention while its source charged full prices, so the PASS rows did not validate the
   implementation. v16 uses ONE convention everywhere — {mul_charge} instructions for a multiply step,
   {sqr_charge} for a squaring step — and certifies exact integer equality between the executed
   per-event charge log ({prod_log_sum}), the reported meter ({prod_meter}) and the row formula
   ({prod_formula}). The cap is a single exact integer, floor(2^{cap_log2}) = {cap_units} work units =
   {cap_instr} instructions, identical in the manifest, the ledger and the certificate; nothing is
   rounded. The toy census now reads: {toy_matches} matches = {toy_ch_same} same-index + {toy_ch_shift}
   shifted-image + {toy_other_failed} unclassified; of the {toy_ch_same} same-index matches
   {toy_rec_same} recover and {toy_same_failed} fail; the {toy_deg} degeneracy total is
   {toy_same_failed} + {toy_shift_failed} + {toy_other_failed}. The v15 text called
   {toy_deg} same-index closures; the corrected count of same-index failures is {toy_same_failed}.

## Accounting

Total work is {w_units_num}/{w_units_den} reference units, log2 in ({w_lo}, {w_hi}), strictly below the
exact cap 2^{cap_log2}. Rows: setup {setup_instr} instructions for the {rtab}-entry multiplier table;
{s_steps} walk steps at the exact average charge (5/8 multiply, 3/8 squaring); the replay row above;
zero for false-match replays (below one instruction). One dense Fp12 multiplication is C = 20000
instructions; the per-step overhead rows are {i_mul} (multiply) and {i_sqr} (squaring).

## Evidence at its actual scope

The production-field run executes {prod_steps} real BLS12-381 GT steps with {prod_sq} squarings
({prod_sq_pct}% ± {prod_ci} at 95%), {inv_checks} invariant checks, {class_tests} class well-definedness
tests and {part_a_pass}/{part_a_total} algebraic predicates. The 52-byte replay pass runs
{replay_trails} trails of {replay_seglen} steps ({replay_steps_run} steps), and the cap-boundary test
confirms that a cap admitting exactly {cap_bound_steps} steps stops at {cap_bound_steps} with the
pre-execution check. The toy ({toy_q}, N = {toy_N}, R = {toy_R}) runs {toy_trials} trials and recovers
{toy_rec} at the design budget. None of this is a production success measurement: the toy hashes a
canonical representative rather than using production trace bits, and the production runs are short.

## Heuristics, disclosed

H1'' (pair-collision rate at least rho/N), H2 (pseudorandomness of the trace low bits, measured at
{prod_sq_pct}% against the design 3/8) and the random-mapping idealisation are heuristics, not theorems.
The shift from the toy field to the 381-bit field is a material extrapolation. The reviewer's standard is
exploratory, and every number above is regenerated from the evidence JSONs by `scripts/ledger16.py` and
`scripts/machine_tables.py`; no number in this note is typed by hand.

## Revision history

- v15 `{v15_id}`: rejected — five internal inconsistencies (0.43 ghost reserve, uncounted replays,
  storage maximum vs expectation, meter convention split plus a rounded cap, census miscount).
- v16 (this package): each item answered by the artifact cited above; the claim is unchanged at {claim}.
