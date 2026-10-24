# p256-ecdlp v11 — verdict record (orbit-row plain-lineage PORT)

## Tickets
| ticket | claim | status | note |
|---|---|---|---|
| 27d28f6 (27d28f66-…) | 128.537 | rejected (unscored/invalid) | first fire of the port |
| 1187a39 (1187a39d-…) | 128.537 | validating at 16:0xZ | rev-2: contradicted assertions fixed, overclaims withdrawn |

## Verdict 1 (27d28f66, verbatim summary of the four objection classes)
1. **Production success unsupported.** The 0.682592 CP lower bound is a small
   additive-group simulator row, not uniform P-256 challenges; the coordinate-based
   production walk's mixing is not shown to match the scalar-class model; two small
   orders do not establish the transfer; the package expressly disclaims the
   DP-density axis (2^-5 → 2^-60) and supplies no alternative argument. The
   restart debit 3/6400 pooled matched lag-1/lag-3 arms sharing starts, and with zero
   lag-3 firings those are not 6,400 independent restart trials.
2. **Computation bound not justified in the protected units.** The leading 1 in
   STEP_PRICE treats a P-256 addition as the borrowed reference operation; C=1000 is
   a normalization, not an implementation bound; the 180-instruction Solinas premium
   and M_insn=97 are asserted, not derived; row-sum checks cannot establish
   log2(W) ≤ 128.537.
3. **Report/retention accounting.** The ledger budgets ~2^68.20 reports against
   2^67 declared record slots; no retention analysis.
4. **Assertions contradicted by the supplied evidence.** (a) the claim called P-256
   "j = 1728"; with a ≠ 0 and b ≠ 0 the correct statement is j ∉ {0, 1728};
   (b) the no-escape control was declared 1/3200 while census-design.json reports
   1933/3200; (c) the claim declared R=64 while the ledger says R=16384;
   (d) absorption was called exhaustive when the script samples 2,048 starting
   classes × 8 secrets with 900-step truncation.

## What rev-2 (1187a39) changed — zero scalar change (claim 128.537 unchanged)
- j statement corrected to j ∉ {0, 1728} and the invariant is now COMPUTED and
  asserted in `ledger_v11.py` (`j_invariant` field in `evidence/ledger-out.json`).
- control-arm number read from `evidence/census-design.json` (1933/3200), not typed.
- fixture R read from the ledger's `fixture_R` (16384) everywhere; the "R=64"
  literal is gone.
- absorption re-described as sampled (16,384 walks/size, 900-step cap), not exhaustive.
- restart denominator = the lag-1 armed block alone (3/3200), matching the matched-pair
  design and the zero lag-3 firing count.
- record capacity raised to 2^69 slots (nodes × cap 2^53) against 2^68.20 reports, and
  the row is emitted (`report_slots_log2`, `records_log2`).
- step-price sensitivity is now a machine-emitted row: charging a P-256 addition as
  TWO reference units gives claim 129.379 (+0.842 bits) and the rank posture is
  unchanged (no rank either way).
- success re-framed as a DECLARED HEURISTIC (H1: the collision law is a function of
  x = visits/sqrt(N) and of the report density in the measured shape), with the
  DP-density transfer named as the open premise instead of being called a stress
  test that conservatively bounds production.

## Surviving structural objection (not fixable in this lane)
Objection class 1 is the same premise that killed secp v17 (verdicts 40b720e2,
b86f973b, cd18d785) and it is not repaired here: the port has no production-density
argument, and the uncounted DP-density cost is a ~2^55 gap in the worst reading. The
P-256 track's frontier (127.237678) is also 1.30 bits below this claim, so the ticket
is a baseline, not a rank move.

## Numbers (all from evidence/ledger-out.json)
- STEP_INSN 263 (P-256 recount: 83 overhead + 180 Solinas premium), STEP_PRICE 1.263.
- log2 S0 127.49244191 (N = (r-1)/2 = 2^255.000000, orbit size 2).
- executed shape x = 1.600202 → log2 S_OPS 128.19973435; log2 W 128.53670650; claim 128.537.
- success 0.681319 = CP95 0.682592 − debits 3.360e-04 − restart 9.375e-04 (floor met).
- order check: CP-95 0.609113 at 255.7× fewer classes (floor met at both orders).
- memory log2 bytes 74.585; preprocessing (table build) log2 43.60.
