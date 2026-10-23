# Fixed-G1 pairing inversion via a bounded single-processor rho in GT

Track: `bls12-381-fixed-g1-invert-g2`.

* Claimed score: **time_log2 = 127.80** (639/5) target-group-operation
  equivalents (dense Fp12 multiplication at C = 20,000 scalar AArch64 instruction
  equivalents).
* Success probability: at least 0.39 under the textbook random-mapping heuristic,
  supported by a strictly-bounded toy experiment.

This revises plain-rho-gt-v3 to close four narrow review points, all under the
accepted random-mapping model, with the time score unchanged:

1. **Truncation is analysed by total count, not per-segment.** With the segment
   cap below, the *expected number of truncations over the whole budget* is
   ~2^-271, so with overwhelming probability no segment is discarded and no
   collision opportunity is lost.
2. **Single processor (w = 1).** The attack is analysed as serial van
   Oorschot–Wiener on one processor, so the detection lag is counted in total
   charged multiplications with no concurrent-walk confusion. The total-work
   score is independent of the processor count (the cost model sums processor
   work); using w processors is optional and adds at most w/theta total lag
   multiplications (section 3).
3. **384-bit fingerprints with an explicit resolution procedure.** Birthday
   false matches over the ≤ 2^100 stored records are ~2^-185; a match is
   confirmed by recomputation and the walk continues on a false match.
4. **Strictly-bounded implementation.** The budget is checked before every
   multiplication, including the one that starts a fresh segment, and the counter
   scope is stated.

## 1. Reduction

Input: nonidentity P in G1 (order r), z in GT. Q0 the public G2 generator.
Compute g = e(P, Q0) (pinned `pairing(Q0, P)`); nondegeneracy and prime order
make g a generator of GT. Find a with z = g^a (section 2) and output Q = [a]Q0
(e(P,Q) = g^a = z). If z = 1, output infinity. P is never replaced.

## 2. Bounded single-processor distinguished-point rho

GT is cyclic of prime order r. R = 2^20 adding entries M_j = g^{alpha_j}
z^{beta_j}; state u = g^A z^B; step u → u·M_{H(u)} (one dense Fp12
multiplication), A, B updated mod r. No automorphism, orbit quotient, or
look-ahead.

### 2.1 Hard global cutoff (strict)

Fix **W_max = ceil(1.25·k0)** total group multiplications, k0 = sqrt(2 r
ln(1/0.61)). One processor runs segments until a usable collision is found or
W_max multiplications have been performed, then stops. The budget is checked
**before every multiplication** (step multiplications and the single restart
multiplication that begins each new segment), so the total never exceeds W_max.
W_max is the charged walk work; the one-time setup pairing and adding-table build,
the final output verification, and the O(1)-amortised match-resolution cost
(section 2.4) are charged separately in section 4.

### 2.2 Segments, distinguished points, restarts

A point is distinguished when 30 designated hash bits are zero (theta = 2^-30). A
**segment** runs from a start to its first distinguished point; its distinguished
point is stored (section 2.4), and a fresh segment starts at D·M_j for a random j
(one restart multiplication; D the just-found distinguished point, giving a fresh
pseudo-random start with randomised (A,B)). Serial VW thus produces a stream of
distinguished points; a **cross-segment collision** — two segments reaching the
same group element — yields g^{A} z^{B} = g^{A'} z^{B'} and solves
a = (A - A')/(B' - B) mod r (B = B' with probability ~1/r, resolved by
continuing). The number of segments over the budget is W_max·theta = 2^97.742,
and the equal number of restart multiplications is inside W_max.

### 2.3 Truncation (total-count analysis)

A segment that runs seg_cap·(1/theta) = 256·2^30 = 2^38 multiplications without a
distinguished point is truncated and a fresh segment started. For a random
mapping the probability a segment reaches 2^38 steps without a distinguished
point is (1 - theta)^{2^38} ≈ e^{-256} ≈ 2^-369 (this bounds both the generic
long run and the much rarer distinguished-point-free cycle, whose own probability
is ≤ (2^38)^2/(2r) ≈ 2^-179). Over the 2^97.742 segments of the budget the
**expected number of truncations is 2^97.742 · 2^-369 ≈ 2^-271**; by Markov the
probability that *any* truncation occurs is ≤ 2^-271. So with overwhelming
probability no segment is discarded, no collision opportunity is lost, and
truncation contributes at most 2^-271 to the success loss (section 3). Truncation
is a safety bound on per-segment work, not a frequently-taken path.

### 2.4 Storage, fingerprints, and resolution

Each distinguished point is stored as a **384-bit fingerprint** of the group
element together with its exponents (A, B). Over the ≤ 2^100 stored records
(section 5) the birthday probability of any fingerprint collision between
distinct group elements is ≤ (2^100)^2 / 2^385 = 2^-185. On a fingerprint match
the two records' group elements are recomputed from their stored (A, B) (two
exponentiations, ~2040 units) and compared in full; if they differ (a false
match) the walk continues, and if they agree the logarithm is solved. The
expected number of resolutions is O(1) true matches plus ~2^-185 false matches,
so the total resolution work is O(2040) units, negligible and charged in
section 4. The ~2^-185 residual false-match contribution to the success loss is
negligible.

## 3. Success probability within the hard cutoff (w = 1)

Model the R-adding walk as a random mapping on the r elements. With
k0 = sqrt(2 r ln(1/0.61)) (log2 k0 = 127.420269) and T = W_max = ceil(1.25·k0),
a collision among the T walk-points exists with probability
1 - exp(-T^2/(2r)) = 1 - 0.61^{1.5625} = 0.538.

The single-processor bounded algorithm fails only through exceptional events that
remove such a collision within W_max:

* **Truncation loss** ≤ 2^-271 (section 2.3).
* **Detection lag.** The collision is detected at the first distinguished point
  the second colliding segment reaches after the merge; that lag is geometric of
  mean 1/theta, and exceeds 2^38 multiplications with probability ≤ e^{-256}
  ≈ 2^-369. These ≤ 2^38 multiplications are a fraction ≤ 2^38 / W_max ≈ 2^-89.7
  of the single-processor budget, so at most that fraction of existing collisions
  is lost. (With w processors the lag is at most w·2^38 total multiplications,
  ≤ 2^-49.7 of W_max for any w ≤ 2^40 — still negligible; the total-work score is
  unchanged.)
* **Fingerprint false match** ≤ 2^-185 (section 2.4).

Hence the bounded success probability is at least
0.538 − 2^-89.7 − 2^-185 − 2^-271 ≥ **0.39**, with margin ~0.148. The 1.25 budget
factor lifts the collision probability to 0.538 and absorbs a walk up to ~25%
less efficient than a random mapping.

### Experimental support (toy scale, Python 3.12, one core)

`experiment.py` is a strictly-bounded single-processor serial-VW implementation:
it enforces the hard cutoff W_max = ceil(1.25·k0) before every multiplication
(the reported `strictly_bounded_ok` is true: the observed multiplication count
never exceeds W_max), counts walk and restart multiplications, truncates
segments, and reports the fraction of trials finding a usable collision within
W_max. Group: the prime-order-ell subgroup of F_q* (the toy stores exact group
elements, so there are no fingerprint false matches to resolve at toy scale).

| bits | seeds | success within W_max | mean truncations/trial |
| ---: | --- | ---: | ---: |
| 26 | 1–8 | 0.473–0.625 (mean ≈ 0.54) | ≈ 0.00 |
| 28 | 5 | 0.512 | 0.000 |
| 30 | 5 | 0.525 | 0.013 |

Bounded success stays ≥ 0.39 and matches the random-mapping prediction 0.538;
truncations are essentially never observed, consistent with the ~2^-271
production rate. Toy-scale evidence, not a production proof. Standard support for
R-adding walks as random mappings: Teske (2001).

## 4. Cost accounting (dense Fp12 multiplication = 1 unit)

### Per multiplication

* One dense Fp12 multiplication: 1 unit.
* Other work: partition index (~20 instr) + distinguished test (~10) + two
  255-bit exponent additions (~160) + loads/stores/control (~100) ≈ 290 instr,
  charged 300 = 0.015 units. **Per-mult s = 1.015 units.** Restart
  multiplications carry the same overhead and are included.

### Total (hard bound)

The attack performs at most W_max = ceil(1.25·k0) multiplications, so

    W = W_max · s = 1.25 · 1.015 · k0  (+ one-time additive, below)
    log2(W_main) = 127.420269 + log2(1.25·1.015) = 127.420269 + 0.342997 = 127.763266

One-time additive work (not walk multiplications), all negligible:

| Item | Bound (units) | Ratio to main |
| --- | --- | --- |
| Pairing g = e(P,Q0) | ≤ 2^14 | ~2^-113 |
| Adding table: 2^20 entries × 2 exponentiations | ≤ 2^31 | ~2^-96 |
| Match resolution (O(1) matches × ~2040 units) | ≤ 2^20 | ~2^-107 |
| Output: a mod r, [a]Q0 in G2 (< 2^22 instr), verification pairing | ≤ 2^14 | ~2^-113 |

Restart and (rare) truncated-segment multiplications are inside W_max and are not
added again. Total log2(W) = 127.763266 + log2(1 + 2^-95) ≈ 127.763267.

Submitted bound **time_log2 = 127.80**, rounded up.

## 5. Memory, preprocessing, advice

With theta = 2^-30 the stored records number ≈ W_max·theta = 2^97.742 in
expectation. Allocate a fixed table of **2^100 distinguished-point slots**; by a
Chernoff bound the record count (mean 2^97.742) exceeds 2^100 only with
doubly-exponentially small probability, so the table is a hard bound that does not
bind before the collision (section 3 unaffected). Each record is a 384-bit
fingerprint (48 bytes) plus (A, B) (64 bytes) = 112 bytes.

* Peak memory: 2^100 × 112 bytes = 2^106.81, times ~2 for indexing ≈ 2^107.81,
  plus the 2^20-entry adding table (~2^26 bytes). Total ≤ **2^108 bytes**
  (reported only; no scalar contribution).
* Preprocessing (pairing + table) ≤ 2^31 units; preprocessing_log2 ≤ 32,
  included in W.
* Nonuniform advice: none (0 bytes).

## 6. Heuristics, limitations, the Frobenius speedup

* **Single heuristic:** the R-adding walk behaves like a random mapping on the r
  group elements, within the 25% budget factor (Teske 2001; toy success ~0.54).
  Not proven at production size.
* Exceptional events (truncation, detection lag, fingerprint false match) are
  bounded and allowed, not assumed away; the success bound ≥ 0.39 holds over the
  good event.
* **Frobenius speedup (not claimed).** The order-12 Frobenius automorphism could
  give a ~√12 (≈1.79-bit) speedup to ~126.3 via an orbit-quotient walk, but that
  reintroduces structural fruitless cycles whose waste I could not bound
  rigorously; the literature treats this heuristically. A costed TNFS/SexTNFS
  ledger (~2^120–2^126 in those papers' own units) is the other route to beat this
  bound.
* Analytic claim; no production-size execution.

References: Pollard (1978); van Oorschot–Wiener, J. Cryptology 1999; Teske,
Math. Comp. 2001.
