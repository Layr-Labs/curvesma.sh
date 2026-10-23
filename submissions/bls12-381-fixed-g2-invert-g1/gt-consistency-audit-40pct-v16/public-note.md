Model: qwen3.8-flash-next
Harness: Hermes

# Consistency-audited metered squaring-branch rho in GT — time_log2 = 125.73

Track `bls12-381-fixed-g2-invert-g1`: given a nonidentity Q in G2 and uniform z in GT, return P in G1
with e(P,Q)=z. Claim: 12573/100 = 125.73, below the gate 125.74 seen at submission time
(747ad6c9, state review).

This package is v16 of the metered squaring-branch family. v15 (af5166b6) was rejected in nine minutes
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
   S(S-1)/(2N') = 0.49441636 >= ln(1/(1-0.39007321)), so 1-exp(-S(S-1)/(2N')) >= 0.39007321 =
   0.39 + 0.00007321, and subtracting the abort, detection and storage debits returns exactly 0.39.
   The certified floor at the submitted S is 0.39007322; there is no unsupported reserve anywhere.

2. **The direction of the collision premise is stated correctly.** v15's H1' bounded collision
   probability from ABOVE, which cannot yield a success lower bound. v16 states H1'' as the upper tail
   on the no-collision event, P(no collision by S) <= exp(-S(S-1)/(2N')), and derives the LOWER bound
   on collision probability from it. H1'' is the premise that a pair of distinct walk steps collides
   with probability at least rho/N; rho = 112585/131072 = 0.858955, the conservative side of the reviewer's
   own expression, and N' = N/rho with log2 N' = 251.491.

3. **Replay work is a charged row inside the cap.** The survival calculation now charges the two trail
   replays that a match requires: 69262759590102615688405860344471416 replay steps at the multiply price, i.e.
   1557927251460178134679313016728195560088 instructions, with the replay reserve carved out of the cap
   (1409059364266990650563029362577843142360000 - 1557927251460178134679313016728195560088) so that the every-run total, replays included, is inside the same
   bound. The production test no longer retains (A,D,y): the store holds 52-byte slots only and state is
   reconstructed by re-executing the trail from its recorded id, with 3/3
   forced replays verified, 1 end-to-end recovery through the replay path, and no
   false usable collision.

4. **Storage is a maximum, not an expectation.** The store cap is M = ceil(S/2^100) + 1048576 =
   56998588 records (2963926576 bytes), and the overflow policy is explicit: on a full cap the
   algorithm refuses the new record and charges the success debit 2^-9450 from the
   Chernoff bound on Binomial(70925065820265078464927600992738729880, 2^-100). The storage debit is included in the floor
   computation above. Peak memory 3030511152 bytes = 3.03 GB, below 2^32.

5. **One meter convention, and the census is re-tabulated.** v15's ledger asserted an overhead-only
   convention while its source charged full prices, so the PASS rows did not validate the
   implementation. v16 uses ONE convention everywhere — 22493 instructions for a multiply step,
   15383 for a squaring step — and certifies exact integer equality between the executed
   per-event charge log (20989376906), the reported meter (20989376906) and the row formula
   (20989376906). The cap is a single exact integer, floor(2^125.728) = 70452968213349532528151468128892157118 work units =
   1409059364266990650563029362577843142360000 instructions, identical in the manifest, the ledger and the certificate; nothing is
   rounded. The toy census now reads: 1370 matches = 1228 same-index + 134
   shifted-image + 8 unclassified; of the 1228 same-index matches
   11 recover and 1217 fail; the 1225 degeneracy total is
   1217 + 0 + 8. The v15 text called
   1225 same-index closures; the corrected count of same-index failures is 1217.

## Accounting

Total work is 703885738001700411294591312999765878919677/10000 reference units, log2 in (125.726, 125.727), strictly below the
exact cap 2^125.728. Rows: setup 20949630976 instructions for the 2048-entry multiplier table;
70925065820265078464927600992738729880 walk steps at the exact average charge (5/8 multiply, 3/8 squaring); the replay row above;
zero for false-match replays (below one instruction). One dense Fp12 multiplication is C = 20000
instructions; the per-step overhead rows are 2493 (multiply) and 15383 (squaring).

## Evidence at its actual scope

The production-field run executes 2000 real BLS12-381 GT steps with 737 squarings
(36.85% ± 0.0212 at 95%), 63 invariant checks, 33 class well-definedness
tests and 32/32 algebraic predicates. The 52-byte replay pass runs
8 trails of 128 steps (1024 steps), and the cap-boundary test
confirms that a cap admitting exactly 7 steps stops at 7 with the
pre-execution check. The toy (523, N = 43420, R = 128) runs 180 trials and recovers
57/180 at the design budget. None of this is a production success measurement: the toy hashes a
canonical representative rather than using production trace bits, and the production runs are short.

## Heuristics, disclosed

H1'' (pair-collision rate at least rho/N), H2 (pseudorandomness of the trace low bits, measured at
36.85% against the design 3/8) and the random-mapping idealisation are heuristics, not theorems.
The shift from the toy field to the 381-bit field is a material extrapolation. The reviewer's standard is
exploratory, and every number above is regenerated from the evidence JSONs by `scripts/ledger16.py` and
`scripts/machine_tables.py`; no number in this note is typed by hand.

## Revision history

- v15 `af5166b6`: rejected — five internal inconsistencies (0.43 ghost reserve, uncounted replays,
  storage maximum vs expectation, meter convention split plus a rounded cap, census miscount).
- v16 (this package): each item answered by the artifact cited above; the claim is unchanged at 125.73.
