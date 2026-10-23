# secp256k1-ecdlp: parallel ρ on classes of six points, truncated at success 0.39

Claim: W ≤ 2^127.313705 target-group-operation equivalents (C = 1,000), success
probability ≥ 0.39 over uniform x and fresh algorithm coins. Scored value **127313705/1000000**.
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
   - **2-cycle look-ahead (exact rule):** at canonical X, let j0 = h(X). Use the
     first j = j0, j0+1, … (mod r) with h(canon(X + R_j)) ≠ j. This is still a
     deterministic function of the class.
   - **Longer fruitless cycles** (cycle length ≥ 4) are caught by a periodic check.
     Every 2^10 steps the walk records its canonical x. If that x recurs within the
     next 2^11 steps, the walk is in a cycle: the trail is abandoned and restarted
     from a fresh random start. Every addition spent counts against the budget T,
     and every restart is charged 2,000 units.

**Lemma (negation 2-cycles are removed).** Let S = canon(X + R_j) = μ(X + R_j)
with μ ∈ μ6, and suppose the step from S uses index j'. A 2-cycle means
S + R_j' = νX for some ν ∈ μ6, i.e. (ν − μ)X = μR_j + R_j'.
- If ν ≠ μ, then ν − μ is a unit mod the group order, so X is one of at most
  36·r² fixed classes. That is probability ≈ 2^−200 per step.
- If ν = μ, the equation needs R_j' = −μR_j. For random table points this forces
  j' = j and μ = −1 (any other coincidence has probability ≤ 6r²/#group), giving
  S + R_j = −X.
This is the family the reviewer identified in other submissions, occurring at a
rate of about 1/(2r) per step. It needs the walk at S to use index j. If X used
j = h(X) (k = 0), the look-ahead had already rejected j, because h(S) = j would be
required. If X itself skipped to j, then reusing j at S also needs h(canon(S + R_j))
= h(X) ≠ j to pass. Either way such a 2-cycle needs two independent index
coincidences, so its rate is O(1/r²). The λ-type 3-cycle
X → λ(X+R_j) → … → X also needs the same index on consecutive steps and is excluded
the same way. The remaining fruitless cycles (≥ 4 steps, plus the O(1/r²)
leftovers) occur at a rate we bound by **2^−20 per step**. The textbook negation-map
4-cycle rate is (r−1)/(4r³) ≈ 2^−24 at r = 2048. In the toy runs, restarts per step
were 3·10^−5 at r = 32 and 7·10^−6 at r = 64, below 1/(4r²) in both cases.

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
gives exactly 0.39, and the 1.04 factor raises this to **0.414**. The success bound
holds as long as the walk's collision efficiency is at least c with
1 − exp(−c²·1.04²·0.4943) ≥ 0.39, i.e. **c ≥ 0.962**. The losses are explicit and add
up to under 0.4%:
- **r-adding walk:** Brent–Pollard's analysis gives a collision-rate factor of
  √(1 − 1/r) = 0.99976 for r = 2048. Bernstein–Lange use the same analysis for
  r-adding walks, and Bernstein–Lange–Schwabe apply it to negation-map class walks.
- **look-ahead additions:** about 1/r = 0.05% of T, plus O(1/r²) for multiple skips.
- **fruitless-cycle waste:** at most 2^−20 per step × at most 2^11 wasted steps =
  0.05%.
- **abandoned long trails:** probability e^−20 each.

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
| Fruitless-cycle restarts: T · 2^−20 × 2,000 units | ≈ 2^117.7 |
| DP hash table initialisation: 2^88 slots × 6 words | ≈ 2^80.6 |
| Decode, constants, output | 2^20 |
| **W total** | **1732343010941939147426954689538026451428897517/8192000** |

log2 W = 127.3137046030…, rounded up to **127.313705**.

**Memory:** ≈ 2^93.6 bytes (DP hash table of 2^88 slots × 48 B, plus walk states). **Preprocessing:**
the step table (≈ 2^22) plus clearing the DP table, ≈ 2^80.6, all included. **Advice:** none.

## 5. Heuristics and limitations

- **H1 (standard Pollard-ρ heuristic).** The r-adding walk on automorphism classes,
  with 2-cycle look-ahead, behaves like a random mapping up to the 4% slack. This is
  supported by the toy runs above and by the literature: Teske on r-adding walks;
  Duursma–Gaudry–Morain on automorphism classes; Bernstein–Lange–Schwabe,
  "On the correct use of the negation map", on fruitless cycles. It is not proven.
- **H2.** The rate of fruitless cycles left after the look-ahead is at most 2^−20 per step. The textbook rate is ≈ 2^−24, and the toy measurements sit below 1/(4r²).
- The per-step instruction count is a Callgrind measurement of a faithful C kernel (`kernels/secp_step.c`, run on random field inputs). A production implementation may differ slightly; the 15-instruction allowance and the separately charged cycle checks cover small differences.
- This is an analytic bound under exploratory review, not a practical attack.

References: van Oorschot–Wiener, J. Cryptology 1999; Gallant–Lambert–Vanstone,
CRYPTO 2001; Duursma–Gaudry–Morain, ASIACRYPT 1999; Bernstein–Lange–Schwabe, PKC 2011.
The construction above is self-contained.
