# v14 — full-store rho walk, priced 4-cycle channel, production-sized record

## What changed from v13 (which was rejected on four flaws)

**0. Claim binding.** The declared claim 0 is now this package's claim, with its
number, in `submissions/p256-ecdlp/submission.json`, and `tracks/p256-ecdlp/
selected-submission.txt` names this tree. The earlier intake failure is fixed: the tree that the selector names is the tree
that carries the claim.

**1. The coefficient-trivial 4-cycle channel is measured and priced, not denied.**
The shipped sign rule really does admit `S -> S+A -> -S-A-B -> -S-B -> S` with both
coefficients unchanged. `experiments/cyclescan.py` enumerates every row pair
(A,B) <= 64 against 64 fixture seeds and finds 29,152 such returns under the shipped
rule (1.78x the analytic `pairs * 2^-4` flip-pattern prediction), so the channel is
real and its frequency depends on R exactly as the review said. No table-level
partition tested (canonicalise every step, every 2 steps, every 3 steps) kills it;
the every-3 variant reports zero at one window alignment only, which is not a proof,
so v14 does not claim a structural kill. Instead the channel is **priced**: at
R = 2^20 the per-4-step death rate is 2^-44, giving <= 2^61.7 deaths per walker over
2^107.67 steps, <= 6 wasted steps per death (<= 4 in the closed cycle plus one
re-seed), i.e. a relative work loss of 8.5e-14 (+0.0000 bits). The census now counts
the event instead of rejecting the trail: `evidence/census_v14.json` records
`trivial_repeat_no_progress` beside `v13_first_repeat_rejects` as a reconciliation
row, with the same success definition (a verified scalar recovery) on both sides.

**2. Selector arm rebuilt at the right depth.** `experiments/selectors.py` compares
the successor class with `hist0` (the 2-step trap behind 1/(2R)), not `hist1` as v13
did; escape is enabled as in production; the design is paired on the same seeds; and
the interval is a cluster bootstrap over whole walks (2000 resamples of 160 walks),
because steps within a trajectory are dependent. The rebuilt delta is
-2.84e-05 with CI [-8.52e-05, 0.00], which contains zero, and the non-degeneracy
gate passes (birthday ratio 0.048 <= 0.1). The arm is therefore **inconclusive and
exploratory**: selector transfer is declared OPEN and no claim row depends on it.

**3. The work bound is analytic, from the track's own pinned prices.**
`experiments/ledger_v14.py` prices each row in instructions and divides by the
protected p256 normalization C = 1000:

| row | Ir/visit | basis |
|---|---|---|
| affine addition, a = -3 | 489 | 2M + 1S at 147 Ir plus 6 add/sub at 8 Ir |
| table read / index | 10 | 3 loads + addressing |
| canonicalisation sign rule | 20 | y parity test + conditional 4-limb negate |
| record append (48 B) | 12 | 6 stores + index bump + bound check |
| batch normalisation | 441 | one 3M Montgomery step per record |
| radix sort + adjacent detect | 128 | 4 key passes over 48-byte records |
| **total** | **1100** | 1.100 units per visit |

No row borrows "one reference unit per addition" — the convention's own text says a
P-256 addition is not one unit. The inversion is a **chain length, not a multiplier**:
the Fermat chain for p-2 is 256 squarings + 129 multiplications = 385 `fe_mul`, the
row v10p already used, so the rejected `M_insn = 97` does not appear. The
executed-count route is **withdrawn**, not repaired: v13's instrumented kernels failed
the instrumented-vs-plain equality control 0/1024, and no executed number enters any
row here. Completeness is by construction: per visit the walk executes exactly one
addition, one canonicalisation, one append, one normalisation share and one
sort/detect share, and it keeps no ring, escape table or carried coefficient state
beyond the 48-byte record, so there is no unpriced executed op.

**4. The record is sized from the recovery requirement.** Identifying one visit among
2^20 trails of 2^108 steps needs 20 + 108 = 128 identity bits, not the 64 v13
shipped. The v14 record is 48 bytes: a 32-byte class key plus a 16-byte identity
field. `experiments/record_recovery.py` packs and unpacks 1200 records through
exactly that layout and recovers the visit's coefficient, 1200/1200 equal to truth.
Records are radix-sorted on the **full 256-bit key** and compared adjacently, so the
lookup addresses the whole key and the rejected 64-bit bucket mix is not used
anywhere.

## The attack

One walk per trail, 2^20 trails of 2^107.67 steps, every visit appended to the store
(no distinguished points, no reporting gate). A repeat with nonzero coefficient
difference yields `da*P + db*Q = O`, hence `x = -da/db mod r`; a repeat with zero
difference is the 4-cycle channel and is counted and re-seeded. Peak storage is
2^133.255 bytes (48 bytes x 2^127.67 records), stated and charged but not scored.

## Numbers

- work: log2(total) = **127.8076** (1.100 units x 2^127.67 visits, plus 8.5e-14 restarts)
- success: lambda = 0.6329, birthday lower **0.4689** >= 0.39
- memory: 2^133.255 bytes
- claim: **2^127.8076** (4-dp ceiling, so "at most" holds)

## What this does not claim

The success law is the standard birthday model over the store, supported by the
fixture census at its own scope and by the priced channel bound — it is not a
production run and not an extrapolated Clopper-Pearson transfer. The 4-cycle channel
is bounded, not eliminated. Selector transfer is open. No executed instruction count
is claimed.
