# secp256k1-ecdlp: parallel ρ on classes of six points, truncated at success 0.39

Claim: W ≤ 2^127.311835 target-group-operation equivalents (C = 1,000), success
probability ≥ 0.39 over uniform x and fresh algorithm coins. Scored value **127311835/1000000**.
`production.py` recomputes every number below from the protected constants.
`experiment.py` runs this walk end to end on toy j=0 curves.

## 1. Facts (verified in `production.py`)

- n is prime and n ≡ 1 (mod 3); p ≡ 1 (mod 3).
- β = 0x851695d4…6afa40 and λ = 0xac9c52b3…b51283ce are primitive cube roots of unity
  in F_p and F_n. They satisfy φ(P) = (βx, y) = [λ]P and φ²(P) = [λ²]P on the
  protected generator P.
- The automorphism group ⟨−1, φ⟩ has order 6. The class of R ≠ O is
  {±R, ±φR, ±φ²R}, whose x-coordinates are exactly {x, βx, −x−βx}, each appearing
  with both signs of y.
- **Canonical representative:** the x-coordinate is the minimum of the three
  canonical integers above, and y is chosen even. Computing it costs one F_p
  multiplication plus additions, comparisons and one conditional negation.

## 2. Algorithm (van Oorschot–Wiener parallel collision search on classes)

1. Choose r = 2048 random pairs (c_i, d_i) and set R_i = [c_i]P + [d_i]Q.
2. Run 2^32 parallel walks in lockstep. Each round's 2^32 additions are batched
   in groups of 2^12 using Montgomery's simultaneous inversion; these are exactly
   the reference operation. A walk starts at canon([a]P + [b]Q) for fresh random
   (a, b). One step is
   X ← canon(X + R_{h(X)}), where h(X) uses 11 bits of the canonical x.
   - **2-cycle look-ahead:** if h(canon(X + R_i)) = i, use i+1 instead (and so on).
     This is the standard fix for fruitless 2-cycles of the negation map.
   - **Longer fruitless cycles** (probability O(1/r²) per step) are caught by a
     periodic check every 2^10 steps. The walk keeps the minimum x over a window
     and restarts if a cycle is found.
3. **Distinguished points:** X is distinguished if the low 40 bits of its canonical x
   are zero. Store the point with its trail's start seed and step count. Then
   restart the trail from a fresh random start. Trails longer than 20·2^40 steps
   are abandoned; this happens with probability e^−20 per trail and is charged.
4. Stop after the budget T of additions in total, or on the first DP collision
   between two distinct trails. Recompute both trails from their seeds while
   tracking coefficients, which gives a[P] + b[Q] = a'[P] + b'[Q] for the same
   canonical point. Output x = (a − a')/(b' − b) mod n, and check [x]P = Q.
   If b = b' the relation is useless; this has probability ~6/n.

**Budget:** T = ⌈1.04 · √(2 ln(1/0.61) · n/6)⌉ = 143650071441869837743865720695360736901
≈ 2^126.756.

## 3. Success probability (heuristic: r-adding walk on classes ≈ random mapping)

There are n/6 classes. With T points visited across all trails, a uniform random
mapping gives P(collision) ≈ 1 − exp(−T²/(2·n/6)). The bare √(2 ln(1/0.61)·n/6)
gives exactly 0.39, and the 1.04 factor raises this to **0.414**. The slack covers:
- look-ahead additions (≈ 1/r = 0.05%);
- abandoned or cycling trails (≤ 2^10 steps lost per fruitless cycle, at O(1/r²)
  per step);
- the known small inefficiency of r-adding walks with r = 2048 (a factor of order
  1 + 1/(2r)).

Collisions occurring before the budget are detected by finishing each live trail
to its next DP after the budget. That tail is charged separately as
2^32 · 2^40 steps, so it is not taken from the slack.

**Toy evidence** (`experiment.py`). Each run counts every look-ahead and restart
addition against the same budget formula. All points are stored, so detection is
immediate; the DP tail is charged analytically above.

| p bits | r | walks | trials | success within budget | predicted (random mapping) |
|---|---|---|---|---|---|
| 24 | 32 | 8 | 400 | 0.4225 | 0.414 |
| 28 | 64 | 16 | 200 | 0.420 | 0.414 |
| 24 | 2048 | 64 | 131 | 0.496 | 0.414 |

Look-ahead overhead was 3.1%, 1.6% and 0.05% of additions respectively (≈ 1/r).
Every recovered logarithm was correct. The effects of a small r are visible and
still stay within the slack.

## 4. Cost ledger (C = 1,000; secp256k1 field multiplication = 179 instructions, organizer measurement)

**Per step: E = 1 reference addition + 470 instructions = 1.47 units.** The non-reference
work is **measured**, not estimated (see `kernels/`). A full class-ρ step and the bare
batched-addition pipeline were both built from the pinned libsecp256k1 field code and
counted with Callgrind. The organizer's own image, compiler flags and pinned checkouts
were used, and the same rig first reproduced the organizer's reference totals exactly.

| Item | Instructions |
|---|---:|
| Batched incomplete affine addition (the reference operation) | 1 unit |
| Measured full step (1,512.7) − measured add-only pipeline (1,061.1): canonicalisation (β·x, −x−βx, three normalisations, min of three, even y), index, look-ahead and DP tests, state store | 452 |
| One inversion (measured ≈ 8.7k) per batch of 2^12 additions | 3 |
| Allowance for address arithmetic, DP store branch and rare look-ahead redo | 15 |
| **Non-reference total** | **470** |

| Component | Units |
|---|---:|
| Steps T · 1.47 | ≈ 2^127.3118 |
| DP tail 2^72 · 1.47 | ≈ 2^72.6 |
| Trail starts (DPs + walks) × 2,000 (two scalar multiplications) | ≈ 2^97.7 |
| DP storage, ≈ 2^86.7 records × 600 instructions | ≈ 2^86.3 |
| Fruitless-cycle checks, 200 instructions per 2^10 steps | ≈ 2^114.4 |
| Step table (2048 × 2 scalar multiplications) | ≈ 2^22 |
| Final trail recomputation (2 × 20 × 2^40 steps) + solve/verify | ≈ 2^45.8 |
| Decode, constants, output | 2^20 |
| **W total** | **5406557745548889748140637652810448373363609/25600** |

log2 W = 127.3118341448…, rounded up to **127.311835**.

**Memory:** ≈ 2^92.3 bytes (DP table at 48 B/record, plus walk states). **Preprocessing:**
the step table, ≈ 2^22, included. **Advice:** none.

## 5. Heuristics and limitations

- **H1 (standard Pollard-ρ heuristic).** The r-adding walk on automorphism classes,
  with 2-cycle look-ahead, behaves like a random mapping up to the 4% slack. This is
  supported by the toy runs above and by the literature: Teske on r-adding walks;
  Duursma–Gaudry–Morain on automorphism classes; Bernstein–Lange–Schwabe,
  "On the correct use of the negation map", on fruitless cycles. It is not proven.
- **H2.** Fruitless cycles longer than 2 are rare enough that periodic detection
  every 2^10 steps wastes far less than the slack.
- The per-step instruction count is a Callgrind measurement of a faithful C kernel (`kernels/secp_step.c`, run on random field inputs). A production implementation may differ slightly; the 15-instruction allowance and the separately charged cycle checks cover small differences.
- This is an analytic bound under exploratory review, not a practical attack.

References: van Oorschot–Wiener, J. Cryptology 1999; Gallant–Lambert–Vanstone,
CRYPTO 2001; Duursma–Gaudry–Morain, ASIACRYPT 1999; Bernstein–Lange–Schwabe, PKC 2011.
The construction above is self-contained.
