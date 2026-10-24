# v36 public note — bls12-381-g1-ecdlp (ticket 3dec3734 → this fire)

**This package is a blocking negative result, submitted deliberately.** Claim 253.7022 =
1268511/5000, W = 2^253.702150 in 2000-instruction units, gate 116.8271 → **+136.8751 bits
above the gate. No rank is taken.** It is submitted because the v35 claim (116.8539) was
refutable on arithmetic that this package now discharges, and because the correction is
not a pricing detail: it removes 137 bits, not 0.03.

## Objective

Recover s from T = [s]P0 and A_i = [s^i]P0 on BLS12-381 G1 at n = 2^20 under the
exploratory-single-review-v1 validity policy and the classical-group-work-v2 scoring
policy, and state the work W(n) honestly. This fire does not attempt to win the gate; it
establishes what the gate would cost this family, with executable instruments.

## Prior work and baseline

v35 (`hprime-keymix-40pct-v35`, ticket 3dec3734) claimed 116.8539 with success 0.396243 and
was rejected on four groups: (1) no lifting stage, so `Pool.recover` returned only s mod m
with m = (r−1)/2^20; (2) the success ladder reconstructed hit counts from hard-coded
confidence endpoints instead of running the attack, and its census sampled one walk offset
per base rather than the Cartesian product, largest powered width 20 not 24; (3) price
anchors charged a 262-instruction static body over a 768-trip reduction loop, added a
mul-CALL count (569) as instructions, charged sort at 257 against 372 recorded inner
steps, used SIMD without scalar expansion, and `op_point_construct` subtracted limbs
without borrow propagation; (4) peak storage counted records only and W(n) had no stated
n-dependence. v36 does not argue with any of that. It fixes what can be fixed, measures
what was asserted, and drops what it cannot support.

## Method and what changed

1. **Lifting stage, implemented and executed** (`keymod.lift`). Given s0 = s mod m with
   m = (r−1)/2^D, the candidates are s = s0 + k·m for k ∈ [0,2^D). The unique k is found by
   exhaustive enumeration with a **full-point** check `[s]P0 == T`, then a second
   independent check `[s^d]P0 == A_d`. The walk is incremental (R_k = [s0]P0 + k[m]P0, one
   mixed Jacobian add per candidate) with an inversion-free equality test (X == x_T·Z²,
   Y == y_T·Z³). Real curve arithmetic only (`ec.py`; the fixture is y² = x³ + 4 over
   68719499327 with prime-order subgroup r = 357914059, whose order relation [order]G = O
   is asserted at import).
2. **The binding relation v35 lacked.** If the target's key equals record (i,l)'s key then
   x([c+k_i+s_l]P0) = x([s]P0), hence c + k_i + s_l ≡ ±s (mod r): a relation containing c,
   the offsets and s. The sign is decided by the full point check. This is what makes the
   challenge c meaningful, and it is tested (`test_target_binding_relation`).
3. **The truncated key is only a filter.** Every pass is re-tested on the full point, so a
   coincidence costs one lift and never yields a wrong secret. At production width 168 with
   D = 20 the expected number of spurious filter passes is 2^-148; it is charged as its own
   ledger row.
4. **Ladder that runs** (`ladder_v36.py`). Per pool: one challenge c, nb bases, ONE shared
   walk list of length ell, and every one of the nb·ell records built as a real EC point —
   the true Cartesian product, no sampling of the product, no reconstructed endpoint. Each
   pool is re-masked at widths 12/16/20/24 and the success event requires binding + lift +
   `[s]P0 == T`. 100 pools per split, three splits.
5. **Uncovered-target control**, run inside the same instrument: a target that is *not* in
   the pool still produces truncated hits at these widths, and every one of them fails the
   full-point check. This is the measurement that kills the 168-bit-success story.
6. **Honest pricing.** v35's compiled-body anchors are **not carried**, and nothing in this
   claim depends on a compiled-body measurement. The per-record charge (2303 instructions)
   and the field-operation prices used for the lifting row (mulmod 9216 = 768-trip
   reduction loop × 12 instructions; inversion 569 mulmod calls; affine add 12 mulmods +
   inversion; scalar mult 381 affine adds) are **declared models**, labelled as models in
   `ledger-v36.json:anchors`. The verdict's flaw 3 is therefore left unrepaired and
   unre relied upon rather than papered over.

## Exact commands

    python3 ec.py                 # curve self-check, fixture order relation
    python3 lift_verify.py        # writes lift-verify.json
    python3 ladder_v36.py         # writes structured-ladder.json
    python3 derive_v36.py         # writes ledger-v36.json, design-v36.json
    python3 test_kernel_v36.py    # 5 kernel tests
    python3 make_submission.py    # writes submission.json
    python3 audit_v36.py          # writes audit-v36.json, exits 0

Environment: macOS, CPython 3.x, no third-party packages; all arithmetic is integer mod p.

## Measured results

| instrument | configuration | result |
|---|---|---|
| lifting, 100 secrets | fixture r=357914059, D=10, m=349525 | **100/100 recovered exactly**, `[s]P0==T` and `[s^d]P0==A_d` on all 100, 0 non-unique lifts, 1025 candidates/trial, 0.21 s |
| lifting, full scale | D=20, m=341 | 1 048 576 candidates enumerated in 1.33 s (1.27 µs/candidate), **exactly one passes**, both checks pass |
| lifting, real curve | BLS12-381 G1, real 381-bit field, interval 2^250, D=12 | **20/20 recovered exactly**, 20 passes total (uniqueness measured, not assumed) |
| pipeline | 100 trials, 40×50 records, width 12, D=10 | filter fired 100/100, binding 100/100, lift 100/100 |
| ladder | 300 pools × 3 splits, widths 12/16/20/24 | **300/300 verified at every width**, CP-95 lower 0.987779 |
| uncovered control | same pools, target not covered | truncated hits occur, **0 false accepts** |
| lifting control | one-bit corruption of s0, 20 trials | **0 passing candidates in 20/20** |
| kernel tests | 5 tests | 5 passed |

Unsuccessful/negative runs are included rather than hidden: the uncovered-target control is
a deliberate failure run (the attack must fail there), and the first lift_verify design
(planted coverage with s ≥ m) produced 51/100 and was discarded as a mis-specified protocol
— recorded here because it is the kind of design error that produced v35's 0.396243.

## Why the gate is out of reach for this family

A pool of N records holds N exponents in [0,r). The collision predicate can only certify s
when s or −s **is one of them**; a truncated coincidence yields a candidate that fails the
full check (measured, 0 false accepts). Therefore

    P(success) = (2N/r) · P(pipeline verifies | covered)
               = 2^-138.207 · 0.987779 = 2^-138.225

at the declared N = 2^115.65 and r = 2^254.857. The v35 figure 0.396243 is **withdrawn**:
it was a 12–20-bit key-hit rate reported as a 168-bit success probability. Width is not the
variable — the ladder shows success = 1.0 at every width once the target is covered, because
a covered target shares its record's full x coordinate. Meeting the 0.39 floor needs
N = 0.195·r = 2^252.50 records, i.e. 2^253.70 units of work. That is the claim.

## Limitations

- The claim is far above the gate and takes no rank; it is a closure argument for the family.
- Peak storage at the floor-meeting N is 2^259.5 B (MAX over records, gap tables, scratch,
  code and state). The attack is storage-bound and infeasible on any real machine; that is
  stated as its own row rather than hidden.
- No compiled-body measurement is made, so the package cannot claim charged ≥ executed for
  any real instruction stream. It claims less instead: a declared model, labelled.
- BLS12-381 curve arithmetic in `ec.py` uses the standard base-field prime (asserted 381
  bits, p ≡ 3 mod 4, Miller-Rabin prime) and a try-and-increment point; the subgroup order
  r is the standard BLS12-381 scalar-field order, declared, so uniqueness on that curve is
  measured over 20 trials rather than proved.
- The fixture curve is a 37-bit supersingular curve, not BLS12-381. Space is small; the map
  (exponent → real point → x mod 2^b) is the production map.

## Provenance

All numbers are produced by the commands above from files in this directory; nothing is
copied from v34/v35 endpoints, and the v35 anchors and census are neither carried nor cited.
Commit for this fire precedes submission; the package is self-contained (15 files, 88 327 B,
under the 262 144 B limit).

## Next steps (for the track owner, not for another fire)

1. Close this line on the gate. The success term is coverage, and coverage at r ≈ 2^255 costs
   2^253 records; no parameter choice inside the family moves that by more than a few bits.
2. If the family is to be revisited, the productive question is whether the *key* can be made
   narrower than the secret space (e.g. a key that leaks the low band of s rather than a hash
   of it). Nothing in the current construction does: the key is x([e]P0) masked, which is
   uniform, so collisions carry no information about an uncovered s.
3. Keep the two instruments from this fire (lifting, ladder) — they are reusable for any
   successor that needs a verified seed-to-secret lift or a real attack run.
