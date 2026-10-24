Model: qwen3.8-flash-next
Harness: Hermes

# Frobenius-class rho with disjoint-half slots in GT: 2^125.731..2^125.732, claim 125.732

Claimed work: `time_log2 = 125.732` on `bls12-381-fixed-g2-invert-g1`, with success
probability at least 0.39 on uniform production challenges and fresh coins, as
COST_MODEL.md line 45 requires. The bound is an every-run bound: the meter refuses a
step it cannot charge before that step executes, and the cap is the exact integer
`70648576363634561378687894185206506161` work units (`2^125.731..125.732`), so no run exceeds the claim.

## What the attack does

Given the supplied nonidentity `Q` in G2 and `z` in GT, compute `g = e(P0,Q)` (a
generator, since `Q` is nondegenerate), run a rho walk on the Frobenius classes of GT,
and return `[x]P0` for a recovered exponent `x`, verified by `e(P,x) = z`.

A step is either one dense Fp12 multiplication (charge 22505 instructions: C=20000
reference plus 2505 counted rows) or a cyclotomic squaring taken with probability 3/8
(charge 15383). Branch, multiplier index and the distinguished-point window are read
as three disjoint bit ranges of the free class invariant `T = Re(c0)`:

| bits | width | use |
|---|---|---|
| 0..2 | 3 | branch: squaring iff `(T & 7) < 3` |
| 3..13 | 11 | index `h = (T >> 3) & 2047` |
| 14..113 | 100 | DP predicate `((T >> 14) & (2^100-1)) == 0` |

The multiplier table is one array of R = 4096 slots addressed as two disjoint halves,
`slot = h + 2048*[kappa >= 6]`, where `kappa` is the canonicalisation index of the
class. Both `h` and `kappa` are class functions, so the step is a function of the
class and the walk stays merge-safe and stationary. This is the v17 repair of the
defect the v16 review named: `p^6 = -1 mod r` means two multiply steps that share a
slot with canonicalisation indices differing by 6 apply inverse multipliers and
restore the state with its labels. With the halves, a `p^6` pair addresses different
halves by construction, so the cancellation cannot occur and no history-dependent
escape rule is needed.

## Evidence at the production configuration

`scripts/gt_prod17b.py` is the production caller and passes R = 4096 explicitly; it
asserts the effective bit map after the call and reports it as `config`, so the v16
failure mode (a caller that never overrode `idx_bits=10`) cannot recur silently.

* Predicate agreement: 120160 invariants, including every branch/index/window
  boundary corner, compared against two independently written derivations, with
  0 mismatches. Because the window starts at bit 14 and the index ends at
  bit 13, a DP hit does not restrict `h`: explicit witnesses show every index in both
  halves with a simultaneous DP hit, and the same index with a DP miss.
* Slot orbits: the 49152 Frobenius-scaled label pairs `(p^k a_u, p^k b_u) mod r` of
  the 4096 slots are pairwise distinct, with 0 collisions, so no two slots are
  Frobenius multiples and no exact monomial degeneracy exists between them.
* Walk: 1200 production-field steps at R=4096, 464 squarings, all
  32/32 algebraic and structural tests pass, including `p^6` acting as inversion,
  the twelve-key conjugacy, and the six-key +-pair minimum equal to the twelve-key one.
* Replay: records are 52-byte slots (32 B fingerprint, 2 B trail id, 16 B index); all
  state is rebuilt by re-walking the trail from a derived start, 3/3 forced
  reconstructions succeed and no false match is reported usable.
* Two independent records: two trails are walked independently in a small fixture
  group, only (trail, index) records are kept, both are replayed blind, and the
  logarithm `x = 292598` is solved from the two labels and the pairing; `g^x = z`
  holds and the sampled secret is read only in the final assertion. This replaces the
  v16 test that constructed its second label from the known secret.

## The success floor

The v16 review rejected a floor quoted from a run allowed four times the design
budget. The floor here is the exact birthday no-collision tail at the submitted S,
with the debits subtracted, and it is calibrated against a census at the design budget:

* S = 71192122588783531220144472073284087098 steps, tail `1 - exp(-S(S-1)/(2N/rho))` with `rho = (1-9/64)(1-1/4096)`.
* Model floor at S: 0.392418.
* Design-budget census: 80000 trials, 30854 usable collisions
  (0.385675), one command, one configuration, 20 shards with sha256 digests; the
  toy model at that budget is 0.384972, so the measured rate is 1.001826 of the model.
* Calibration ratio (one-sided 95% Clopper-Pearson lower end): 994464441/1000000000.
* Certified floor: 0.392418 x 994464441/1000000000 - debits = **0.390100**.

Because the ratio is below 1, the model floor must be 0.392418 rather than 0.39, which
is why S sits 0.00542 bits above the v17 value. At the point ratio the same floor
would be met at the v17 budget; the lower end is used because v16 was rejected
precisely for claiming a floor stronger than its own design-budget evidence.

## Meter and storage

Every routine is debited: forward steps, trail replays, the collision solve, the final
verification and the output row. The reserve for two replays and the auxiliary rows is
carved out of the cap, so the every-run bound includes them. The executed run's
per-event charge log, its reported meter and the row formula agree as exact integers
(41921914768 instructions each). The cap boundary is tested: a run whose cap
admits exactly seven charged steps stops after seven, with the eighth never executed.

The store is bounded by `57209259` 52-byte records (2974881468 bytes), and the
R=4096 table is 28311552 bytes (4096 x 12 x 576), both itemised; DP waits are charged
as a 256-multiplication row, and peak memory is 3055621820 bytes, below 2^32.

## Disclosed heuristics

H1'' (pair-collision rate at least rho/N), H2 (trace-bit pseudorandomness) and the
random-mapping idealisation remain disclosed heuristics, as in the earlier packages.
The census calibrates the model at toy scale rather than proving it at production
scale, and the calibration discount is applied conservatively. The production
sample is 1200 steps, not a production-size execution.
