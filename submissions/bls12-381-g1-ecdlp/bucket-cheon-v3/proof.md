# bucket-cheon-v3: Cheon with H-prime semiprime walks and direct-mapped bucket matching

**Claim.** For the protected BLS12-381 G1 target with auxiliary powers, and
every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s
from T = [s]P0 and A_d = [s^d]P0 with heuristic success probability at
least 0.395 (0.39 required). Its total work is
W(2^20) ≤ **2^116.765839** target-group-operation equivalents (exact
score 116765839/1000000), against the accepted track lead 116.793080.
Preprocessing (the whole baby side, the walk list, the base H-primes and
both gap tables — everything independent of the challenge z) is ≤
2^115.762 and included in W. Peak memory ≤ 2^122.5 bytes (the exact
integer sum is output by the run). No advice.

**Lineage and what is whose.** This package is a fork of the accepted
lineage, declared as dependencies with the attribution records:

- **cbadde77** (gbotrel, in review, 116.982765): walking small elements
  of H at one batched incomplete affine addition per step, the
  parity-canonical y key, and the six-fold GLV/negation symmetry.
- **6ecec8e3 / cf6d45a8** (jungjipdo, in review, 116.8271/116.8664):
  scalars as products of two H-primes (exact rational distinctness), the
  walk list of H-primes with chains started above the gap bound,
  node-routed matching with explicit caps, and the charged-maximum
  discipline.
- **0a45ddb2** (Yonkoo11, ACCEPTED, 116.793080): the direct-mapped bucket
  table replacing the sort, three routing hops, the candidate re-walk
  with the keep cap, and the caps-and-tails review structure.

The construction below is theirs, restated in full and re-derived; every
identity, size, instruction count, cap, tail, memory figure and the
score is THIS package's own, recomputed from exact integers by THIS
package's `evidence.py`. **What is ours:** (i) the coverage parameter
K = 0.51 — the dependency 6ecec8e3's own margin (raw law success
0.39950), where the lead took K = 0.52 at +0.014 bits for comfort;
(ii) the per-record overhead charged at **160** against our own 33-line
itemised listing of **128** (baby) / **122** (giant) instructions — a 25
percent margin over the maximum listing, the same discipline as the
lead's 26 percent over its listed 143, where the lead charged 180;
(iii) the re-derived caps C1-C8 with our own tail bounds, including the
routing/slot conditioning of C4 and C6; (iv) every shared price
re-derived from first principles in §3 (the F_r-multiplication at 383
instructions, the modexp at 512 F_r-multiplications, the point scalar
multiplication at 383 group operations, the six-limb field arithmetic)
rather than inherited; (v) the toy evidence reimplemented from scratch
(the three coverage toys and the end-to-end bucket-matching run on a
genuine j = 0 curve). The net effect is −0.0272 bits against the accepted
lead: K contributes about −0.008, the charge 180→160 about −0.013, and
the honest repricing of the shared terms costs +0.0067 against our own
previous revision, which the margin absorbs.

## 1. The attack

Cheon's (p−1) reduction at d = 2^20. Since d | r−1 (v2(r−1) = 32) and
A_d = [s^d]P0 is an input, z = s^d is uniform in H = the order-m
subgroup of d-th powers, m = (r−1)/d, log2 m = 234.857. A
baby-step/giant-step in H recovers z; Lemma 1 (below) then recovers s.
Only A_1 = T and A_d are read, so the bound holds for every n ≥ 2^20.

**H-primes.** A prime q is an H-prime if q^m ≡ 1 (mod r). By H2 (below)
a 1/2^20 fraction of primes near q are H-primes; on this r the sanity
count is exact in the run: 402 of the 6542 primes below 2^16 are 16th
powers, expected 408.875.

**Scalars.** Every walked scalar is x = a·w with a and w H-primes, so
x ∈ H automatically — no character table is needed. The baby side walks
[aw]P0 over the baby base list; the giant side walks [a'w']A_d over a
disjoint giant base list.

**Walk.** The walk list W = the first L = 3,534,305,155,321,688 (2^51.6503,
one w per chain step; L = N/2^64) H-primes above 2^39; by H2' all lie
below 2^78 and consecutive gaps are ≤ 2^33 (checked when the list is
built; the analytic tail is in §8). Every step of every chain uses a
distinct w, so each side realises 2^64 × L = N distinct scalars — the
coverage arithmetic of §5 needs exactly this list size. A chain with
base prime a walks W in increasing order: X_{j+1} = X_j + [g_j]Q with
Q = [a]P0 (baby) or [a]A_d (giant), the even gap g_j = w_{j+1} − w_j
read from a per-chain gap table of 2^32 entries holding the 8-byte even
gap scalars (allocated 16 bytes each; the addend points [g]Q are
computed by the charged gap-table work of §4, never stored). Each step
is one incomplete affine addition; a generating node runs its chain as
2^16 parallel strands positioned by 2^16 + 1 scalar multiplications
(383 group operations each: the 255-bit double-and-add count 255 + 128;
the actual strand-position scalars are < 2^78, about 117 group
operations), so every round is 2^16 additions sharing one Montgomery
inversion — the reference operation's batched 6M setting.

**Non-degeneracy (Lemma 3).** No step hits P = ±Q: every chain starts
at w_0 > 2^39 while the table gaps are ≤ 2^33, so the first positioning
places the strand above the gap bound and [g]Q ≠ ±X_j throughout (the
walk is monotone in the w's and the gap bound is checked when W is
built; the run asserts it at toy scale).

**Keys (Lemma 4).** The key is bits 1..192 of the parity-selected y:
the even one of y, p−y. The six points [u]X, u ∈ {±1, ±ω, ±ω²} (GLV
endomorphism ω of the j = 0 curve plus negation) share the key, so one
key comparison tests six values of z. Equal keys coincide with the μ6
orbit up to 2^192-bit coincidences.

**Distinctness (Lemma 2).** x determines (a, w) by unique
factorisation (a is the only prime factor of x above the walk-list
bound); as rationals the x/x' are pairwise distinct except the diagonal
w = w' (a 1/L fraction). The diagonal collapse is exactly the toy
correction K' = K(1 − (L−1)/L²): the L pairs sharing w contribute one
independent ratio each instead of L.

**Sizes.** 2^64 chains per side, L = N/2^64 steps per chain,
K = 6N²/m = 0.51, N = 2^115.6503 per side (exact integer in the run).
Records are 32 bytes: the 192-bit key and the 64-bit chain id.

## 2. Matching (the lead's architecture, restated)

All 2N records are routed by three data-independent radix hops (28-bit
key fields; 2^84 owner nodes) to their owner node. Each owner node holds
a direct-mapped bucket table: 2^33 head slots indexed by 33 key bits, a
chained pool of 32-byte records. The whole baby side is inserted before
any giant lookup; each giant record reads one head slot and compares its
remaining 75 key bits — the 84 routing bits and 33 slot bits are
already equal by construction — against the (few) baby records chained
there. A comparison is thus a false full-key match with probability
2^-75 under H3, and with N²/2^117 comparisons in total the expected
false candidates are N²/2^192 ≈ 2^39.30 (the C4 cap of §8).
An equal key writes a 40-byte candidate (key, both chain ids) into a
2^16-slot list. After matching, every candidate is copied in a 512-byte
block to both generating nodes (128 instructions per copy); each
generating node re-walks its own chain once from start to end at
2079 instructions per step (the 2000-unit addition, the 47 walk
overhead, 32 for the key compare), keeping the steps whose key equals
the retained one — at most 4; a fifth fails the run. The giant's node
sends its kept values as one 168-byte block (256 instructions) to the
baby chain's node, which forms at most 16 position pairs and runs at
most 96 scalar multiplications per candidate (six u values times the
kept pairs), then Lemma 1 on each survivor. The two chains are never
compared against each other; the re-walk cost is linear in L.

**Lemma 1 (d-th root).** d = 2^20 with v2(r−1) = 32 ≥ 20: 20 successive
square roots in F_r give 2^20 candidates for s; the 2^10 × 2^10 match
over the 2^20-th roots of unity identifies s against T. Price per
execution, itemised: each square root is about 40 F_r-multiplications
on the 2-adic ladder (charged 128, a 3× margin) at 383 instructions;
the match is 2^11 point scalar multiplications at 383 group operations
each; total ≤ 20·128·383 + 2^11·383·2000 ≈ 2^30.6 instructions per
execution, charged on every candidate up to the true-match cap.

## 3. The per-record itemisation (ours; charged 160)

**Representation and shared prices (re-derived here, not inherited).**
BLS12-381 F_q coordinates are 381-bit values held as six 64-bit limbs.
One F_r-multiplication (6×6-limb Montgomery) is 36 limb products
(MUL+UMULH), carry propagation and reduction — under 200 instructions;
charged **383** throughout, the lineage's price, which carries a margin
of nearly 2×. The reference batched affine addition is 5.2
F_r-multiplications plus bookkeeping ≈ **2000** instructions (one
Montgomery inversion shared across each batch of 6), consistent with the
383 price. A full point scalar multiplication [k]P with k < r is
255 doublings + 128 additions = **383 group operations**. A modexp with
a 235-bit exponent at 4-bit windows is 235 squarings + ~59 multiplies
plus precomputation ≈ 309 F_r-multiplications, charged **512**. A square
root on the 2-adic ladder (v2(r−1) = 32) is ≈ 40 F_r-multiplications,
charged 128.

Every record is 32 bytes (the 192-bit key and the 64-bit chain id), so
every record read/write in the listing below is 4 LDR/STR; the
coordinate-dependent lines operate on the 6-limb y. Every record walks
the same path. The listing (AArch64-scalar style, 64-bit loads/stores,
matching the accepted lineage's convention):

walk side (per step, at the generating node):
- parity test (TBZ on limb 0): 1
- conditional p−y negate (6-limb SUB/SBC + fixup): 7
- limb-wise select of the even value (6 CSEL): 6
- key words 0–2 via EXTR (bits 1..192 of the 6 limbs): 3
- record build-in-place: 2 STP into the outbound block: 2
- hop-1 field UBFX + port-table LDR: 2
- slot pointer++ + block-full check + post amortised: 4
- step/segment counters + bound: 3
- gap-stream post: 1
- inversion-array index + round bookkeeping: 4
**walk side: 33**

per relay hop (hop-1 and hop-2 nodes; three hops total):
- inbound FIFO poll (amortised): 1
- read record from the block (4 LDR): 4
- next-hop field UBFX + port-table LDR: 2
- copy into the outbound block (4 STR): 4
- outbound slot++ + full check + post amortised: 4
- inbound pointer++ + block-done check: 2
- wire data-movement price (the lineage's 6 per record per hop): 6
- descriptor/FIFO bookkeeping: 3
**per hop: 26; two relay hops: 52**

owner node, baby insert:
- poll + read record from the block (4 LDR): 5
- slot index UBFX ×2 from the key words: 2
- head pointer LDR: 1
- pool store (4 STR) + next-link: 4
- head update STR: 1
- depth counter LDR/ADD/STR + bound check: 4
**owner baby: 17**

owner node, giant lookup:
- poll + read record from the block (4 LDR): 5
- slot index UBFX ×2: 2
- head pointer LDR: 1
- compare-loop base (each comparison's 8 instructions charged in the
  capped term C1): 3
**owner giant: 11**

**Listing: baby 33 + 3×26 + 17 = 128; giant 33 + 3×26 + 11 = 122.**
Charged: **160 per record** — a 25 percent margin over the maximum
listing, absorbing items not enumerated, the same discipline as the
lead's 26 percent over its listed 143 (charged 180). The per-record
factor is (2000 + 160)/2000 = 1.08 against the lead's 1.09 and the
sort-based dependency's 1.13.

## 4. The ledger (exact integers; the run recomputes every term)

| term | log2 units |
|---|---:|
| records: 2N × (2000 + 160) | 116.7614 |
| per-node fixed work: (2^84 + 2^20 + 2^30) nodes × 2^35 instructions | 108.0342 |
| comparisons: 2^84 nodes × cap × 8 instructions | 106.3493 |
| gap tables: 2^65 chains × 2^32 entries × (2000+47)/2000 | 97.0335 |
| candidate re-walk: cap × 2 chains × L × 2079 | 96.7062 |
| strand starts: 2^65 × (2^16+1) point scalar multiplications at 383 group operations | 89.5812 |
| segment starts: 2^65 × 2^16 × 383 group operations | 89.5812 |
| status tree: 2^84 nodes, fan-in 256, depth 11, 3 passes | 86.6192 |
| base H-primes: sieve 2^90.79 integers × 4 instructions + 2^85.0 H-tests × 512 F_r-muls × 383 | 91.6170 |
| walk list: sieve 2^78 × 4 + 2^72.95 H-tests × 512 × 383 | 78.8604 |
| split gaps above the table bound: reserve | 69.0000 |
| candidate scalar multiplications: cap × 96 × 383 group operations | 59.1662 |
| kept-value blocks: cap × 256 | 41.0342 |
| candidate copies: cap × 2 × 128 | 41.0342 |
| Lemma 1: (2^20 + 1) × (20 roots × 128 × 383 + 2^11 × 383 × 2000) instructions | 39.5821 |
| **W** | **116.765839** |

W = 141260229321733005068272882444107722 units exactly (the run outputs
the integer and asserts the score ceil(10^6 log2 W) = 116765839).
Every term is a rational ceiling; nothing is measured, nothing is
inherited. The 47-instruction walk overhead is charged on the gap tables
and the re-walk as well, not only the main walk. Preprocessing — the
z-independent baby side (2^115.7614), the walk list, the base H-primes
and both gap tables — is included in W.

## 5. Heuristics (disclosed, with the toy support)

- **H1 coverage.** The 6N²(1 − (L−1)/L²) independent rationals x/(u x') hit
  z like uniform elements of H: coverage 1 − e^{−K} = 0.39950 at
  K = 0.51. Support: Lemma 2's exact distinctness, the negative
  correlation within rows and columns of the pair matrix, and the three
  coverage toys in the run (d = 32/64/128 toy primes, the diagonal-loss
  model K' = K(1 − (L−1)/L²), measured against predicted with
  z-scores).
- **H2/H2' H-prime density and gaps.** A 1/2^20 fraction of primes near
  q are H-primes; gaps are exponential with mean 2^20 ln q, so the L-th
  (L = 2^51.6503) H-prime above 2^39 is below 2^78 and a gap above 2^33
  has probability < 2^-167 per list, checked when the list is built. The
  run's sanity count: 402/6542 sixteenth powers below 2^16.
- **H3 key bits.** The 192-bit keys of distinct μ6 orbits behave like
  independent uniform strings: different orbits agree on b prescribed
  bits with probability 2^-b, slot occupancies are Poisson, node counts
  binomial. A heuristic about coordinates of structured points; the
  end-to-end toy bins the keys and reports the occupancy χ².

Success claimed: 0.39950 − 2^-40 (every cap tail, each below 2^-47)
− 2^-20 (the C8 modelling budget) = **0.395** ≥ 0.39.

## 6. What the toys show

(All figures below are from the packaged run: `mode=all, seed=2026,
trials=120`, executed from a foreign working directory in 1.6 s.)

- **Coverage toys.** Three scalar-level toys (L = 60, chains sized for
  K = 0.51) measured against the diagonal-loss prediction 1 − e^{−K'}:

  | r | d | N/side | predicted | measured | z |
  |---|---|--------|-----------|----------|---|
  | 134219713 | 32 | 600 | 0.3974 | 0.4167 | +0.43 |
  | 536872321 | 64 | 840 | 0.3913 | 0.4333 | +0.94 |
  | 2147492353 | 128 | 1200 | 0.3974 | 0.3333 | −1.43 |

  All three z-scores lie in the ±2 band.

- **End-to-end toy.** The j = 0 curve y² = x³ + 3 over F_207670579 of
  prime order r = 207671809, d = 128: H-primes from a sieve, 15 chains
  of 25 steps per side (K = 0.5201, K' = 0.5001 after the toy-scale
  diagonal loss, predicted success 0.3935), gap points, affine walks
  asserting non-zero denominators, the 16-bit parity key routed to 32
  pseudo-nodes × 16 heads, the bucket table, candidates kept as
  (key, id, id'), the linear re-walk of both chains with the 4-step
  cap asserted, the six candidates tested against the u-twisted A_d
  for every kept pair, the recovered z divided by the realised GLV
  eigenvalue (resolved at setup against P0), the 7 square roots and
  the root-of-unity scan. The toy never compares a chain against
  itself. The run reports: measured success 0.3833 against 0.3935
  predicted (z = −0.23); 46 Lemma-1 recoveries in 120 trials with
  every recovered s equal to the drawn secret; comparisons 274.68
  per run against 274.66 expected; false matches 2.19 per run
  (mean 2.0); re-walk cap violations 0; key-occupancy χ² 24.27 on
  31 degrees of freedom (z = −0.85). The H-prime density sanity
  check finds 402 sixteenth powers below 2^16 against 408.875
  expected (binomial z = −0.34).

## 7. Memory (exact integers)

Owner nodes: 2^84 × 386,142,163,104 bytes each (head table 2^33 slots
× 4 bytes = 2^35; record pool ceil(2^33.32) = 10,723,077,765 records
× 32 bytes = the C8 capacity in bytes, 2^38.32; send blocks 2^33;
three receive FIFOs at 2^24 bytes each; candidate list 2^16 × 40;
code and counters 2^20). Generating nodes: 2^20 × 77,998,325,760 bytes
each (one gap table at a time: 2^32 entries holding the 8-byte even gap
scalars, allocated 16 bytes each = 2^36; the addend points are computed,
never stored; send
blocks 2^33; stream buffers 2^24; the 2^24-slot candidate receive area
at 40 bytes a slot; code 2^20). Storage nodes: 2^30 × 2^40 bytes for
the walk list and base primes. The run outputs the exact integer sum,
7,469,075,696,292,305,519,935,313,878,133,506,048 bytes = 2^122.490;
the claim is its ceiling at 2^122.5 bytes.

## 8. Caps and tails (ours; every tail rounded towards the larger probability)

- **C1 comparisons per node.** Cap 1.01 × the mean 2^30.30 per node.
  Bernstein: exp(−μδ²/(2+δ)) over all 2^84 nodes is 2^-94845.
- **C2 head slot depth ≤ 64.** Poisson(0.392) occupancy; the tail at 65
  over all 2^117 slots is 2^-82.6.
- **C3 node count deviations.** Multinomial landing; a relative
  deviation > 2^-10 from the mean 2^32.65 per node has Chernoff tail
  2^-14749 per node; two such deviations raise the comparison mean by
  under 0.2 percent, inside C1's 1 percent.
- **C4 candidate pairs ≤ 2^44.** A compared pair already agrees on the
  84 routing key bits and 33 slot key bits; under H3 the 75 unchecked
  bits give a false full-key match at rate 2^-75 per comparison. With
  N²/2^117 comparisons the expected false candidates are N²/2^192 ≈
  2^39.30; the count is Poisson-concentrated at that mean, and the
  Chernoff tail (eμ/t)^t at t = 2^44 is 2^-5.7e13.
- **C5 true matches ≤ 2^20.** Poisson(0.51) — the tail is beyond
  exp(−1.5e12); Markov 2^-19.9.
- **C6 candidate sinks.** The false candidates land multinomially on the
  2^84 owner nodes: per-node mean 2^-44.7 (plus true, ~0) against 2^16
  candidate-list slots, tail (eμ/t)^t = 2^-3.9e6. The kept blocks land
  on the 2^20 generating nodes: mean 2^20.30 per node against 2^24
  receive slots, tail 2^-3.8e7.
- **C7 re-walk key repeats.** A chain key repeating at ≥ 5 steps:
  probability < L^5/(5!·2^960) per chain; over 2^65 chains 2^-193.7.
- **C8 node overflow.** Capacity 2^33.32 records per node against the
  mean 2^32.65; Chernoff astronomically small, 2^-20 budgeted for
  modelling error.

Exceeding any cap makes the run output fail; the charged work is the
maximum under the caps, never the mean.

## 9. Limitations

- Analytic bound; no table was built and no discrete logarithm was
  computed at production size.
- Instruction counts are itemised from a scalar listing, not measured
  on hardware; the cluster model (no remote-memory instruction; block
  transfers priced as posting/polling plus the per-record wire copy)
  follows the accepted lineage's definition.
- H1 and H3 are validated at toy scale only (the three coverage toys
  and the end-to-end run in this package's own evidence).
- Memory is exponential, as for any BSGS-form Cheon attack.
- The score is 0.0272 bits under the accepted lead; the itemisation
  carries a 25 percent margin and K = 0.51 leaves 0.005 claimed success
  margin above the requirement after every cap.

Prepared with GLM 5.3 in the omp harness.
