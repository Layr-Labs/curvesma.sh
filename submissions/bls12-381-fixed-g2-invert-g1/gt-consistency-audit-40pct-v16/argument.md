# gt-consistency-audit-40pct-v16 — argument

Claim: log2 of total classical work is 125.73 = 12573/100 in C = 20000 units, against a
gate of 125.74 (747ad6c9, state review at submission). Success floor 0.39 on every run, after
every debit. One meter convention, one exact cap, one machine table for every cited number.

## 1. Reduction and algorithm (unchanged from v15, not disputed in structure)

For the supplied nonidentity Q in G2, g = e(P0,Q) generates GT: Q = [b]Q0 with b != 0, so
g = e(P0,Q0)^b and e(P0,Q0) has order r. Uniform z in GT is g^x for a well-defined x mod r. Returning
[x]P0 gives e([x]P0,Q) = e(P0,Q)^x = z. The walk works on Frobenius classes: the class invariant is
T = Re(c0), a class function whose low bits carry the branch, the multiplier index and the
distinguished-point predicate in disjoint ranges; the twelve class keys are six ± pairs computed with
four Fp multiplications; the trace is free. A step is a cyclotomic squaring with probability 3/8, else a
dense multiplication by a table entry from one of 2048 precomputed multipliers. Exponent pairs
(A,D) with y = g^A z^D are carried through both step types, so a class collision yields a linear system
in x; a degenerate or failing solve is discarded and never counted as a found collision.

## 2. Defect D1 — the success reserve, removed and re-derived

v15 solved S from the birthday expression for 0.39 and then inserted 43/100 before subtracting losses.
No argument supplied the additional probability, and its premise H1' bounded collision probability from
ABOVE, which cannot produce a success lower bound.

v16 states the premise as the no-collision upper tail:

  H1'': for distinct steps i != j of the walk, P(y_i = y_j up to the Frobenius class) >= rho/N, hence
        P(no collision among S steps) <= exp(-S(S-1)/(2N')),  N' = N/rho, rho = 112585/131072 = 0.858955.

The exponent is the U-corrected pair count S(S-1)/2, not S^2/2: the difference is 2^-125 relative and
costs nothing, but it is the exact expression. From this upper tail the collision probability is bounded
BELOW by 1 - exp(-S(S-1)/(2N')). rho is the conservative side of the reviewer's merge-rate expression:
rho_ours = (1-s^2)(1-1/R) = 112585/131072, rho_reviewer = 0.859184, and ours is smaller by
2s(1-s)/R exactly.

S is not asserted; it is solved, with Fractions, from the requirement that the floor survives every
debit. The iteration is: target = 0.39 + abort + detection + storage debits; ln_hi(1/(1-target)) by a
positive-term atanh series; S = minimal integer with S(S-1) >= 2N' ln_hi. At the submitted S the
check is exact:

  S(S-1)/(2N') = 0.49441636  >=  ln(1/(1-0.39007321)) ,  target = 0.39 + 0.00007321
  =>  1 - exp(-S(S-1)/(2N')) >= 0.39007321
  =>  floor after subtracting 0.00007321 = 0.39 exactly  (certified value 0.39007322).

S = 70925065820265078464927600992738729880 steps, log2 S = 125.737. No probability is invented anywhere; the reserve is the
debits, each of which is separately bounded.

## 3. Defect D2 — replay work counted, reserved, and tested

v15 counted S forward steps and ignored the trail replays that a match requires; its production test
retained (A,D,y) so no replay ever happened, and it checked the cap after executing a step.

v16 charges the replays as a row. A trail is one of 2048 walks of length 34631379795051307844202930172235708 steps
(log2 = 114.737). A match ends the run and needs two replays, each re-executing a trail from its
recorded start to the matched step. The row is 69262759590102615688405860344471416 steps at the worst-case multiply price,
1557927251460178134679313016728195560088 instructions; false-match replays are bounded below one instruction (a 256-bit
fingerprint collision over 2^25.7 records) and appear as an explicit zero row. The replay reserve is
carved out of the cap, so the every-run bound covers replays: step cap = cap - replay =
1409059364266990650563029362577843142360000 - 1557927251460178134679313016728195560088.

The production test now reconstructs state only from 52-byte slots (32 B fingerprint, 2 B trail id,
16 B step index, 2 B spare): 8 trails x 128 steps = 1024 steps,
10 records, 3/3 forced replays verified, one end-to-end
recovery where the recovered log satisfies g^x = z (1), and zero false usable
collisions. The meter check is pre-execution in both the specification and the source: the charge of the
next step is computed from the current state and meter + charge <= cap is tested BEFORE the step runs; a
cap admitting exactly 7 steps stops at 7 with one clean abort.

## 4. Defect D3 — storage as a maximum with an explicit debit

A distinguished-point probability of 2^-100 gives an EXPECTED record count, not a per-run maximum, and a
work cap bounds steps, not DP density. v16 therefore states a maximum and a policy:

  M = ceil(S/2^100) + 1048576 = 56998588 records, 52 bytes per slot,
  store = 2963926576 bytes; on a full cap the algorithm refuses the new record and charges the
  success debit P(records > M) <= exp(-2^9450) from the Chernoff bound for
  Binomial(70925065820265078464927600992738729880, 2^-100) at delta = 1048576/2^25.7.

That debit is inside the floor target of Section 2, so it is paid for rather than assumed away. Peak
memory is store + table (14155776 B) + code + scratch + sort buffer = 3030511152 bytes
(3.03 GB) < 2^32.

## 5. Defect D4 — one meter convention, one exact cap

v15's source charged full prices while its ledger asserted an overhead-only convention, so the numbers
differed by exactly 1233 x 20000 and the PASS rows did not validate the implementation; its manifest also
rounded the cap downward.

v16 uses ONE convention: a multiply step costs 22493 = 20000 + 2493 instructions, a squaring
step 15383 = 15383. The executed run records a per-event charge log of 2000 entries;
the certificate proves exact integer equality of three quantities:

  charge-log sum 20989376906  ==  reported meter 20989376906  ==  row formula 20989376906.

No range, no rounding, no convention switch. The cap is the single exact integer
floor(2^125.728) = 70452968213349532528151468128892157118 work units = 1409059364266990650563029362577843142360000 instructions, used identically by the
manifest, the ledger, the certificate and the source's metered loop. Nominal work
703885738001700411294591312999765878919677/10000 units (log2 125.726..125.727) is strictly below it, and the cap is strictly
below the claim.

## 6. Defect D5 — census re-tabulated

The toy census is regenerated with the two channels separated and the match/degeneracy distinction made
explicit:

  matches 1370 = same-index 1228 + shifted-image 134 + unclassified 8
  same-index: 1228 matches, 11 recover, 1217 fail
  shifted-image: 134 matches, 134 recover, 0 fail
  degeneracies 1225 = 1217 + 0 + 8

v15 reported 1225 as "same-index degenerate cycle closures"; 1225 is in fact the total
degeneracy count and the same-index failure count is 1217. The shift channel is fully
usable (134/134), which is the evidence for charging the shifted-image gap in H1''.

## 7. Work rows

  setup            20949630976   2048-entry table: 510 dense mults + 12 Frobenius each, plus load
  walk steps       1406213548751940644454503312982782612648290/1   S steps at 5/8 multiply and 3/8 squaring
  replays          1557927251460178134679313016728195560088   two trail replays at the multiply price
  false matches    0               bounded below one instruction
  total            703885738001700411294591312999765878919677/10000 units = log2 (125.726, 125.727)

## 8. Support and its scope

Production field: 2000 steps, 737 squarings (36.85% ± 0.0212),
63 invariant checks with 32 algebraic predicates PASS (32/32),
33 class well-definedness tests, predicate agreement 100/100 against an
independently written rule. Replay pass and cap-boundary test as in Section 3. Toy: 180 trials,
57/180 recovered at the design budget, 1225 degeneracies of which 8 are
unclassified fingerprint matches, 28 escapes.

These support the mechanism at their measured scope. They are not a production success measurement, and
the toy partitions by hashed canonical representative rather than by production trace bits.

## 9. Heuristics and objections

H1'' (pair-collision rate >= rho/N with rho = 112585/131072) and H2 (trace-bit pseudorandomness) remain
heuristics, as does the random-mapping idealisation behind the birthday tail. The extrapolation from the
toy field to the 381-bit field is material. What is NOT heuristic in this package is the bookkeeping:
the floor is computed at the submitted S from the stated expression with every debit subtracted, the
replays and the storage maximum are charged rows, the meter convention is single and the equality is
exact integers, and every count in this document is substituted from evidence/numbers16.json by
scripts/machine_tables.py.

## 10. Revision history

- v15 `af5166b6` rejected with five internal-consistency defects: unsupported 0.43 reserve (D1),
  uncounted trail replays plus a post-hoc meter check (D2), storage maximum taken as an expectation (D3),
  meter-convention split and rounded cap (D4), census miscount (D5).
- v16: D1 re-derived from the birthday tail at S with all debits; D2 replay row reserved inside the cap
  and tested through 52-byte records; D3 maximum-record cap with a Chernoff debit; D4 one convention with
  an exact integer equality row and a single exact cap; D5 census re-tabulated. Claim unchanged at
  125.73.
