# bucket-cheon-v11: Cheon with H-prime semiprime walks and direct-mapped bucket matching

**Claim.** For the protected BLS12-381 G1 target with auxiliary powers, and
every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s
from T = [s]P0 and A_d = [s^d]P0 with heuristic success probability at
least 0.3907 (0.39 required). Its total work is
W(2^20) ≤ **2^116.740733** target-group-operation equivalents (exact
score 116740733/1000000), against the accepted track lead 116.793080 and
the pending fork entry 116.761224. Preprocessing (the whole baby side,
the walk list, the base H-primes and both gap tables — everything
independent of the challenge z — the COMPLETE cold-attack
composition: the baby records, the walk list, the base H-primes, both
gap tables, the baby side's positioning share, and ALL per-node setup
including the owner head tables and depth counters) is ≤ 2^115.7432
(the run outputs the exact integer sum) and included in W. Peak memory ≤ 2^122.81 bytes (the
exact integer sum is output by the run). No advice.

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
K = F(4980,10000) = 0.498 exactly (raw law success 0.39226), where the
lead took K = 0.52 at +0.014 bits for comfort; (ii) the per-record
overhead charged at **148** against our own 33-line itemised listing of
**129** (baby) / **122** (giant) instructions — a 14.7 percent margin
over the maximum listing, the same discipline as the lead's 26 percent
over its listed 143, where the lead charged 180; (iii) the re-derived
caps C1-C9 with our own tail bounds, including the routing/slot
conditioning of C4 and C6, C1 in the bounded-increment Bernstein form
with the C3 deviation carried inside it, and C5's two failure channels;
(iv) every shared price re-derived from first principles in §3 as a
worst-case bound with its algorithm (the F_r-multiplication at 383
instructions, the modexp at 512 F_r-multiplications, the inversion at
2^16 instructions by binary extended Euclid, the point scalar
multiplication at 662 reference additions by projective 4-bit windows
INCLUDING the window table's affine normalisation at the full standard
batch-inversion cost, the square root at 2048 F_r-multiplications by
Tonelli–Shanks, the six-limb field arithmetic) rather than inherited;
(v) the toy evidence reimplemented from scratch (the three coverage
toys, the L- and base-stratified exact-coverage census of §6, and the end-to-end
bucket-matching run on a genuine j = 0 curve); (vi) the explicit walk
schedule of §1 — the retained value-indexed gap table [g]Q built per
chain by the dependent affine recurrence with one own field inversion
per entry, its build charged as the gap-tables term, its storage as the
generating-node memory, its per-step reads in the walk overhead, and
the rare > 2^33 gaps split against the reserve; (vii) N as an exact
multiple of 2^64, so 2^64 · L = N is an integer identity. The net
effect is −0.0523 bits against the accepted lead and −0.0205 against
the pending fork entry: K 0.52→0.498 contributes about −0.045 across
the record, comparison and re-walk terms, the charge 180→148 about
−0.008, and the honest scalar-mul repricing 640→662 adds back +0.031
on the strand, segment, candidate and Lemma-1 terms (which are
subdominant, so the score moves by under 0.001 net of them).

**This revision (v11) answers the v10 review's three objections:**
(1) the preprocessing claim is now the COMPLETE z-independent
composition — the baby side, the walk list, the base H-primes, both
gap tables, half the positioning (the baby side's chains) and all
per-node setup, the owner head-table and depth-counter initialisation
included (it lives inside the ledger's per-node fixed work: 2^35
instructions per node cover the 2^33 head-slot stores, the 2^33
depth-counter words, the FIFOs, the candidate lists and code —
nothing the insert algorithm reads is left uninitialised and
uncharged); (2) the census fit is reported with BOTH parameters —
|dev| ≈ 0.235/L + 0.000514 — so the model's L→∞ floor (5.14e-4) is
stated, not silently discarded, and the 0.0015 allowance is justified
as 2.9× that floor; (3) the base-prime-count extrapolation (two toy
bases against production's 2^64) is controlled two ways: provably,
cross-base-pair ratio coincidences are mod-r coincidences of unequal
rationals whose expected production count is 2^128·36/m ≈ 2^-69, so
the base-pair blocks are independent and the union coverage is the
law; and empirically, census strata at nb = 2, 4 and 8 base primes
show the residual does not grow with nb.

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
x ∈ H automatically — no character table is needed. The base H-primes
a are searched from 2^80 upward over the fixed-point interval
[2^80, F] whose length X = F − 2^80 solves X = 2^66·2^20·ln(F)
(X ≈ 2^92.0, F ≈ 2^92.0): since 1/ln(q) ≥ 1/ln(F) on the interval,
the expected supply is at least X/(2^20·ln(F)) = 2^66 — twice the
2^65 requirement — with the shortfall a Chernoff tail under H2''
(C9, §8). The baby side walks
[aw]P0 over the baby base list; the giant side walks [a'w']A_d over a
disjoint giant base list.

**Walk.** The walk list W = the first L = 3,492,477,585,540,061
(2^51.6332, one w per chain step; L = N/2^64 exactly) H-primes above
2^39; by H2' all lie below 2^78 and consecutive gaps are ≤ 2^33 w.h.p.
(checked when the list is built; the analytic tail is in §8). Every
step of every chain uses a distinct w, so each side realises 2^64 × L
= N distinct scalars — the coverage arithmetic of §5 needs exactly
this list size. A chain with base prime a walks W in increasing order:
X_{j+1} = X_j + [g_j]Q with Q = [a]P0 (baby) or [a]A_d (giant), the
even gap g_j = w_{j+1} − w_j.

**The gap table (explicit schedule).** Per chain, the generating node
builds and RETAINS the value-indexed point table [g]Q for every even
g < 2^33 — 2^32 affine points at 96 bytes (two 381-bit coordinates),
one chain's table at a time. The build is a dependent affine
recurrence on entries T_i = [2(i+1)]Q: T_0 = [2]Q is one doubling of
Q, T_1 = T_0 + [2]Q is one doubling of T_0 (the formula's two
exceptional cases, each priced as one reference addition), and for
i ≥ 1, T_{i+1} = T_i + [2]Q is a general incomplete addition —
T_i = ±[2]Q would need 2(i+1) ≡ ±2 (mod r), i.e. i = 0 or i+1 ≈ 2^253,
both outside the general-addition range. Because each entry depends
on the previous one, the inversions CANNOT be batched across the
build: each entry pays its own affine addition (2000 instructions)
plus one field inversion (BINARY EXTENDED Euclid, charged 2^16 =
65,536 instructions with the full four-value derivation of §3), plus
the 47-instruction walk overhead — (2000 + 65536 + 47)/2000 ≈ 33.79
reference additions per entry, the gap-tables term of §4.
Every walk step then READS its addend [g_j]Q from the table by gap
value (the read is charged in the 47-instruction walk overhead); a
gap above 2^33, whose H2' probability is < 2^-167 per gap, splits
into two table entries and a second addition (the split-gaps reserve
of §4). The chain runs as 2^16 lockstep segments positioned by
2^16 + 1 scalar multiplications (662 reference additions each by
the projective 4-bit-window worst case of §3), all segments reading
the same retained table, one Montgomery inversion per round of 2^16
additions — the reference operation's batched 6M setting (the
lockstep rounds are independent, so their inversions DO batch). A
candidate re-walk rebuilds its two chains' tables first, at the same
per-entry price (charged in the re-walk term of §4). The addend
points are therefore available at every one of the L steps with no
uncharged recomputation: the schedule's total work is the table
build (2^32 dependent additions with own inversions per chain) plus
the walk itself (L batched additions per chain), both charged in §4,
and its storage is the 2^32 × 96-byte retained table, charged in §7.

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

**Sizes.** N is the smallest exact multiple of 2^64 with
6N²/m ≥ K = 4980/10000, so every chain has exactly L = N/2^64 steps
and 2^64 × L = N holds as an integer identity; the realised 6N²/m is
0.498 × (1 + 2^-54), and every coverage figure below is stated at the
lower bound K = 0.498. N = 64424840203624563870472181515288576 =
2^115.6332 per side (exact integer in the run).
The wire record is 32 bytes: the 192-bit key and the 64-bit chain id;
the owner node's pool copy adds a 34-bit pool link (40 bytes total).

## 2. Matching (the lead's architecture, restated)

All 2N records are routed by three data-independent radix hops (28-bit
key fields; 2^84 owner nodes) to their owner node. Each owner node holds
a direct-mapped bucket table: 2^33 head slots indexed by 33 key bits,
each head a 5-byte pool index (the pool holds up to ceil(2^33.32)
records, beyond a 32-bit index), a chained pool of 40-byte linked
records (the 32-byte wire payload plus the 34-bit link). The whole baby
side is inserted before
any giant lookup; each giant record reads one head slot and compares its
remaining 75 key bits — the 84 routing bits and 33 slot bits are
already equal by construction — against the (few) baby records chained
there. A comparison is thus a false full-key match with probability
2^-75 under H3, and with N²/2^117 comparisons in total the expected
false candidates are N²/2^192 ≈ 2^39.27 (the C4 cap of §8).
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
execution, itemised: each square root is exponentiation by the odd
part of (r−1)/2 (a 227-bit exponent at 4-bit windows: 227 squarings +
57 multiplies + precomputation) plus the 2-adic correction ladder
(~32 further squarings) ≈ 320 F_r-multiplications, charged **512** like
a modexp; the match is 2^11 point scalar multiplications at 662
reference additions each; total ≤ 20·512·383 + 2^11·662·2000 ≈ 2^31.3
instructions per execution, charged on every candidate up to the
true-match cap.

## 3. The per-record itemisation (ours; charged 148)

**The organiser's normalisation clause (why instruction-level pricing is
the sanctioned unit convention).** The track's `COST_MODEL.md` fixes the
score semantics: W is a count of *target-group-operation equivalents*,
one reference operation = one unit, and the priced quantities are
defined by rounded constants — in the organiser's own words, the
prices are "**rounded organizer normalization constants** grounded in
the pinned measurements; **they are not CPU timings or exact measured
instruction counts**". The same file sets `candidate_prices_allowed:
false`: a candidate may not re-price the reference operation at its
own implementation's instruction count. Two consequences, which this
package follows: (a) the reference addition's 2000-instruction budget
is a *normalisation constant*, so non-group work must be converted to
units through it — the per-record overhead below is charged as
(2000 + 148)/2000 = 1.074 reference operations per record, i.e. actual
itemised instructions divided by the pinned reference budget, never
re-scaled to a cheaper "measured" mul; (b) where a sub-routine cannot
batch inversions (the gap-table build, one scalar multiplication),
its field operations are priced in F_r-multiplication units at the
pinned 383 instructions per F_r-multiplication and the pinned
65,536-instruction inversion, again through the 2000-instruction
reference — the unit convention of the normative files, applied
uniformly. What the clause does NOT do is fix the sub-routine costs
themselves: those are derived here, itemised, as worst-case bounds.

**Representation and shared prices (re-derived here, not inherited).**
BLS12-381 F_q coordinates are 381-bit values held as six 64-bit limbs.
One F_r-multiplication (6×6-limb Montgomery) is 36 limb products
(MUL+UMULH), carry propagation and reduction — under 200 instructions;
charged **383** throughout, the lineage's price, which carries a margin
of nearly 2×. One field **inversion** by BINARY EXTENDED Euclid on 381-bit
values. The algorithm maintains FOUR six-limb values: the two GCD
operands (a, b) and the two Bézout coefficients (u1, u2). Per
iteration: two conditional 6-limb compares with subtracts (~24 limb
operations), two 6-limb right shifts (~12), the coefficient halving
with its parity fixup — a conditional 6-limb add of r before the
shift (~14) — and control (~6): ≈ 56 limb operations. At most
2 × 381 = 762 iterations: ≤ 42,672 instructions; charged
**2^16 = 65,536**. (The earlier revisions priced only the operand
updates, omitting the coefficient maintenance.) It is priced
EXPLICITLY wherever it cannot be batched (the gap-table build) and
amortised at 1/6 per addition (≈ 10,923) where it can (the lockstep
walk rounds). A point **doubling** in affine form is the same shape as
the reference addition — 2 F_r-multiplications plus an inversion share
— and is priced as one reference addition. The reference **batched
affine addition** closes at the honest F_r prices: at six limbs a
schoolbook Montgomery multiplication is under 200 instructions (the
383 charged elsewhere is a deliberate 2× overcharge), so one addition
is 3M + 1S ≈ 800, field adds/subs and control ≈ 170, and the
Montgomery-batched inversion share at the lockstep batch of 2^16
(3 muls per value plus one inversion per batch) ≈ 603: ≈ 1573
instructions, covered by the lineage's **2000**. This 2000 normalises
additions and doublings; anything that cannot batch (a lone
inversion, a scalar multiplication) is priced in its own units below,
never silently as 2000. A point **scalar multiplication** [k]P with
k < r, fully itemised at the 383 price: 4-bit fixed windows in
Jacobian coordinates — 252 doublings × (4M + 4S), 63 window
selections × mixed addition (8M + 3S), a 15-entry precomputation ×
(8M + 3S): 2874 F_r-multiplications = 1,100,742 instructions;
non-multiplication field work (the coordinate adds, subs and
compares of the same formulas) ≈ 150 per operation × 330 = 49,500;
window indexing 504; the final affine conversion 3M = 1149 plus one
inversion 65,536; **the 15-entry window table normalised to affine
entries by ONE batched inversion, at the standard algorithm's full
cost: the prefix/backward Montgomery batch inverts n values with
3(n−1) = 42 F_r-multiplications plus the one inversion, and each
entry's affine rescale x = X/z², y = Y/z³ then costs 4 further
multiplications (z², z³ and the two coordinate divisions) = 60 more:
102 F_r-multiplications + 65,536 = 39,066 + 65,536 = 104,662
instructions** (mixed additions require affine addends, so this
normalisation is on the critical path of every scalar
multiplication); control ≈ 1000. Total ≈ **1,322,593 instructions =
661.3 reference additions; charged 662**. (v8 charged 640 without
the table normalisation; v9 charged 648 with an undercounted batch —
29 multiplications is the prefix pass alone; the standard backward
pass and the per-entry rescale bring the honest figure to 102, and
the charge to 662.)
A **modexp** with a 235-bit exponent at 4-bit windows is
235 squarings + ~59 multiplies + precomputation ≈ 309
F_r-multiplications, charged **512**. A **square root** in F_r by Tonelli–Shanks at v2(r−1) = 32 is the
odd-part exponentiation (227-bit exponent at 4-bit windows: 227
squarings + 57 multiplies + precomputation ≈ 310 F_r-multiplications)
plus at most 32 correction rounds, each a worst-case 32-step order
search and update (≤ 34 F_r-multiplications): at most ≈ 1400
F_r-multiplications; charged **2048**.

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
- pool store (5 STR: 40-byte linked record) + next-link: 5
- head update STR: 1
- depth counter LDR/ADD/STR + bound check: 4
**owner baby: 18**

owner node, giant lookup:
- poll + read record from the block (4 LDR): 5
- slot index UBFX ×2: 2
- head pointer LDR: 1
- compare-loop base (each comparison's 8 instructions charged in the
  capped term C1): 3
**owner giant: 11**

**Listing: baby 33 + 3×26 + 18 = 129; giant 33 + 3×26 + 11 = 122.**
Charged: **148 per record** — a 14.7 percent margin over the maximum
listing, absorbing items not enumerated, the same discipline as the
lead's 26 percent over its listed 143 (charged 180). The per-record
factor is (2000 + 148)/2000 = 1.074 against the lead's 1.09 and the
sort-based dependency's 1.13.

## 4. The ledger (exact integers; the run recomputes every term)

| term | log2 units |
|---|---:|
| records: 2N × (2000 + 148) | 116.7362 |
| per-node fixed work: (2^84 + 2^20 + 2^30) nodes × 2^35 instructions | 108.0342 |
| comparisons: 2^84 nodes × cap 1,304,364,425 × 8 instructions | 106.3149 |
| gap tables: 2^65 chains × (2^32+1) entries × (2000+65536+47)/2000 (dependent build, own inversion per entry) | 102.0786 |
| candidate re-walk: cap × 2 chains × (L × 2079 + the table rebuild) | 96.6891 |
| strand starts: 2^65 × (2^16+1) point scalar multiplications at 662 reference additions | 90.3707 |
| segment starts: 2^65 × 2^16 × 662 | 90.3707 |
| status tree: 2^84 nodes, fan-in 256, depth 11, 3 passes | 86.6192 |
| base H-primes: sieve the fixed-point interval (X = 2^66·2^20·ln(2^80+X) ≈ 2^92.0, supply ≥ 2^66) + ~2^86 H-tests × 512 × 383 | 92.6173 |
| walk list: sieve 2^78 × 4 + 2^72.95 H-tests × 512 × 383 | 78.8604 |
| split gaps above the table bound: reserve | 69.0000 |
| candidate scalar multiplications: cap × 96 × 662 | 59.9556 |
| kept-value blocks: cap × 256 | 41.0342 |
| candidate copies: cap × 2 × 128 | 41.0342 |
| Lemma 1: (2^20 + 1) × (20 roots × 2048 × 383 + 2^11 × 662 × 2000) instructions | 40.3790 |
| **W** | **116.740733** |

W = 138823277198598325615651887091741570 units exactly (the run outputs
the integer and asserts the score ceil(10^6 log2 W) = 116740733).
Every term is a rational ceiling; nothing is measured, nothing is
inherited. The 47-instruction walk overhead is charged on the gap tables
and the re-walk as well, not only the main walk. Preprocessing — the
COMPLETE z-independent composition: the baby records (2^115.7362), the
gap tables, the walk list and the base H-primes (together
2^102.0806), the baby side's share of the strand and segment
positioning (2^90.3707), and the per-node setup (2^108.0342: per owner
node 2^35 instructions cover the 2^33 head-slot stores, the 2^33
depth-counter words, the FIFO setup, the candidate lists and code —
the initialisation the insert algorithm reads) — sums to 2^115.74319,
claimed at its computed ceiling 2^115.7432, and is included in W.
(v9 stated a truncated partial sum; v10 stated a partial composition
whose ceiling left only 2^99.7 of rounding slack the head-table
initialisation would exceed; this revision states the complete
composition, so the claim is a genuine cold-attack bound.)

## 5. Heuristics (disclosed, with the toy support)

- **H1 coverage.** The 6N²(1 − (L−1)/L²) independent rationals x/(u x') hit
  z like uniform elements of H: coverage 1 − e^{−K} = 0.39226 at
  K = 0.498. Support: Lemma 2's exact distinctness, the negative
  correlation within rows and columns of the pair matrix, the three
  coverage toys in the run (d = 32/64/128 toy primes, the diagonal-loss
  model K' = K(1 − (L−1)/L²), measured against predicted with
  z-scores), and — the load-bearing quantitative support — the
  L- and base-stratified exact-coverage census of §6: on twenty-six
  genuine toy primes in five strata (twenty at two base primes with
  L = 74 to 295, six more at four and eight base primes), the exact
  covered set is enumerated, its density matches the same
  diagonal-corrected law with a residual whose mean absolute value
  shrinks monotonically in L, and the two-parameter fit
  |dev| ≈ 0.235/L + 0.000514 carries an explicit L→∞ floor of
  5.14e-4 with rms per-prime scatter 0.00136.
- **H2/H2' H-prime density and gaps.** A 1/2^20 fraction of primes near
  q are H-primes; gaps are exponential with mean 2^20 ln q, so the L-th
  (L = 2^51.6332) H-prime above 2^39 is below 2^78 and a gap above 2^33
  has probability < 2^-167 per list, checked when the list is built. The
  run's sanity count: 402/6542 sixteenth powers below 2^16.
- **H2'' membership independence.** The H-primality events of the
  sieved primes are independent, so the base-prime supply concentrates
  at its mean and the C9 Chernoff tail applies. A heuristic about
  prime membership in a number-theoretic subgroup; the sanity count
  above is consistent with it but does not prove it.
- **H3 key bits.** The 192-bit keys of distinct μ6 orbits behave like
  independent uniform strings: different orbits agree on b prescribed
  bits with probability 2^-b, slot occupancies are Poisson, node counts
  binomial. A heuristic about coordinates of structured points; the
  end-to-end toy bins the keys and reports the occupancy χ².

Success claimed: 0.39226 − 2^-40 (every cap tail, each below 2^-47)
− 2^-20 (the C8 modelling budget) − 0.0015 (the census-informed
allowance) = **0.39075, claimed at 0.3907** ≥ 0.39. The allowance is
justified from the census's two-parameter fit, both parameters
reported: |dev| ≈ 0.235/L + 0.000514, so the model's L→∞ floor is
5.14e-4 and the 0.0015 allowance is 2.9× that floor. The per-prime
scatter around the fit (rms 0.00136, worst single prime 0.0067 at
L = 104) is toy-scale mod-r coincidence luck whose expected count
does not grow with the modulus; the fit's floor is the systematic
part, and the allowance triples it.

## 6. What the toys and the census show

(All figures below are from the packaged run: `mode=all, seed=2026,
trials=120`, executed from a foreign working directory in 4.5 s.)

- **Coverage toys.** Three scalar-level toys (L = 60, chains sized for
  K = 0.498) measured against the diagonal-loss prediction 1 − e^{−K'}:

  | r | d | N/side | predicted | measured | z |
  |---|---|--------|-----------|----------|---|
  | 134219713 | 32 | 600 | 0.3974 | 0.4167 | +0.43 |
  | 536872321 | 64 | 840 | 0.3913 | 0.4333 | +0.94 |
  | 2147492353 | 128 | 1200 | 0.3974 | 0.3333 | −1.43 |

  All three z-scores lie in the ±2 band.

- **The L- and base-stratified exact-coverage census (the load-bearing
  H1 support).** For twenty-six genuine primes r_t with
  6·d_t | r_t − 1 (d_t = 64 or 128), the census builds the SAME
  semiprime structure as production — disjoint base lists of H-primes
  times one shared walk list of L entries, with L/N sized at the
  claim's K — and then enumerates the covered set EXACTLY:
  C = {b·u/g mod r_t : b ∈ B, g ∈ G, u ∈ μ6}, all of which lies in H.
  The exact density |C|/m_t is compared with the diagonal-corrected law
  1 − e^{−K(1−(L−1)/L²)} at the realised K' (0.489 to 0.495). Twenty
  primes run at two base primes per side in four growing strata
  (r_t ≈ 2^25 to 2^28, L = 74 to 295); six more hold the modulus scale
  while growing the base-prime count (three at nb = 4, three at
  nb = 8).

  **The wrap regime, stated honestly.** At production the scalars do
  NOT wrap: a < 2^93 and w < 2^78, so a·w < 2^171 < r ≈ 2^255, and
  Lemma 2's rational distinctness therefore holds at the scalar level
  exactly — distinct semiprime scalars are distinct residues. The
  census toys are smaller and some rows' semiprime products exceed
  r_t; the run reports `products_wrap` per row, and the wrapping and
  non-wrapping subsets agree (mean |dev| 0.0023 against 0.0024),
  so the regime distinction does not drive the residual. No dominance
  between the regimes is claimed — both are measured, both are
  reported.

  **The base-prime count, from two to 2^64.** The census cannot run
  2^64 base primes; two arguments control the extrapolation.
  Provable: a cross-base-pair ratio coincidence b·u/g = b'·u'/g'
  with disjoint base lists is a mod-r coincidence of unequal
  rationals (unique factorisation forbids integer equality), so its
  probability per cross base-pair is ~36/m; at production the base
  pairs number 2^64 × 2^64 = 2^128 and the expected cross coincidences
  are 2^128 · 36/m ≈ 2^-69 — the base-pair blocks are independent,
  and the union coverage over 2^128 blocks is the law. Empirical: the
  nb-strata residuals do not grow with nb —

  | nb (base primes/side) | n | mean L | mean abs dev | worst |
  |---|---|-----|------|-------|
  | 2 | 20 | 151 | 0.00237 | 0.00671 |
  | 4 | 3 | 84 | 0.00137 | 0.00198 |
  | 8 | 3 | 59 | 0.00077 | 0.00159 |

  — if anything they shrink (more base pairs average away per-prime
  coincidence luck).

  **The residual: both fit parameters reported.** Per L-stratum at
  nb = 2 (L range, mean |dev|, worst |dev|):

  | stratum | n | L | mean abs dev | worst |
  |---|---|-----|------|-------|
  | r_t ≈ 2^25 | 6 | 74–104 | 0.00361 | 0.00671 |
  | r_t ≈ 2^26 | 5 | 104–148 | 0.00244 | 0.00349 |
  | r_t ≈ 2^27 | 5 | 148–209 | 0.00170 | 0.00388 |
  | r_t ≈ 2^28 | 4 | 209–295 | 0.00127 | 0.00211 |

  The mean absolute deviation shrinks monotonically as L grows,
  across a fourfold L range, with signed means at +0.0013 or below in
  every stratum (no systematic bias; overall 16 positive, 10
  negative). The least-squares fit of |dev| against 1/L over the
  twenty nb = 2 primes is **|dev| ≈ 0.235/L + 0.000514 — both
  parameters reported**: the slope dominates at toy L (it predicts
  0.0037 at L = 74 and 0.0013 at L = 295, matching the strata), and
  the intercept is the model's L→∞ floor, 5.14e-4, with rms per-prime
  scatter 0.00136 around the fit. The claim's 0.0015 allowance is
  **2.9× the fitted floor**. (v9's census paragraph wrongly described
  the wrapping toys as the production regime with a 0.0008 allowance;
  v10 showed the shrink but reported only the fit's slope, discarding
  the intercept; this revision reports both parameters and sizes the
  allowance from the floor.)

  The census is scalar-level, exact, and reproducible from the packaged
  run; it measures the coverage law itself at the claim's K rather
  than only the end-to-end pipeline.

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

  The toy's gap handling matches the claimed production schedule in
  the load-bearing respect: it builds a per-chain table of the addend
  points [g]Q and RETAINS it for the whole chain, reading each step's
  addend from it, exactly as §1's schedule retains the value-indexed
  table (the v3 text wrongly said the points were never stored, which
  the toy contradicted; this revision states what both actually do). At the toy's
  24-gap scale the table is built by direct scalar multiplication where
  production builds the value-indexed 2^32-entry table by the dependent
  recurrence of §1 — a toy-scale artifact; the production build cost is the
  gap-tables term of §4 and the toy's build is not charged against any
  production figure.

## 7. Memory (exact integers)

Owner nodes: 2^84 × 480,516,719,816 bytes each (head table 2^33 slots
× 5-byte pool indices = 2^35.7; record pool ceil(2^33.32) =
10,723,077,765 linked records × 40 bytes — the 32-byte wire payload
plus a 34-bit pool link, the C8 capacity — = 2^38.65; send blocks 2^33;
three receive FIFOs at 2^24 bytes each; candidate list 2^16 × 40;
code and counters 2^20). Generating nodes: 2^20 × 421,595,709,440
bytes each (the retained value-indexed gap table [g]Q, one chain's at
a time: 2^32 affine points at 96 bytes = 2^38.58, per the explicit
schedule of §1; send blocks 2^33; stream buffers 2^24; the 2^24-slot
candidate receive area at 40 bytes a slot; code 2^20). Storage nodes:
2^30 × 2^40 bytes for the walk list and base primes. The run outputs
the exact integer sum,
9,294,545,109,473,455,968,826,126,053,901,598,720 bytes = 2^122.806;
the claim is its ceiling at 2^122.81 bytes.

## 8. Caps and tails (ours; every tail rounded towards the larger probability)

- **C1 comparisons per node.** Cap 1,304,364,425 = 1.01 × the mean
  1.2914 × 2^30 per node. The comparison count of one node is a sum of
  bounded increments — each giant lookup contributes the occupancy of
  the one head slot it reads, an increment in [0, 64] (the C2 cap; a
  deeper slot fails the run) — so Bernstein's inequality for bounded
  increments applies with M = 64 and variance bound 64μ:
  P(X − μ ≥ t) ≤ exp(−t²/(2(64μ + 64t/3))). The tail is evaluated at
  the WORST mean inside the C3 envelope (see C3): μ_worst = μ(1 + 2^-9)
  against the cap 1.01μ, i.e. t = 0.008μ: the per-node tail is 2^-938.2
  and the union over all 2^84 nodes is **2^-854.2**. (v8 priced C1 as a
  unit-Bernoulli Bernstein exp(−μδ²/(2+δ)), which ignored the bounded
  increments; this revision carries the correct form, and the C3
  shifted mean inside it rather than as a separate remark.)
- **C2 head slot depth ≤ 64.** Poisson(0.388) occupancy; the tail at 65
  over all 2^117 slots is 2^-82.6.
- **C3 node count deviations.** Multinomial landing; a relative
  deviation > 2^-10 from the mean 2^31.63 baby records (2^32.63 total)
  per node has Chernoff tail 2^-3055 per node, 2^-2971 over all
  nodes. The deviation's COST is not charged here: C1's tail is
  computed at the worst mean within this envelope, so the cap bounds
  the shifted mean itself and C3 prices only the envelope's breach.
- **C4 candidate pairs ≤ 2^44.** A compared pair already agrees on the
  84 routing key bits and 33 slot key bits; under H3 the 75 unchecked
  bits give a false full-key match at rate 2^-75 per comparison. With
  N²/2^117 comparisons the expected false candidates are N²/2^192 ≈
  2^39.27; the count is Poisson-concentrated at that mean, and the
  Chernoff tail (eμ/t)^t at t = 2^44 is 2^-5.7e13.
- **C5 true matches ≤ 2^20.** Two failure channels. (i) Off-diagonal:
  the count is Poisson(0.498)-concentrated; the tail at 2^20 is
  (eμ/t)^t = exp(−1.42e7) = 2^-2.0e7. (ii) Diagonal challenges: for a
  baby base a and giant base a' with z = a/a' (over the six twists),
  every one of the L shared-w positions is a true match, so the run
  aborts at the caps (the same single event that would break C4 and
  C7); the union over all 2^64 × 2^64 base pairs and six twists has
  probability 6·2^128/m = 2^-104.27. Both channels are far below
  2^-40 and are charged inside the common 2^-40 cap deduction.
- **C6 candidate sinks.** The false candidates land multinomially on the
  2^84 owner nodes: per-node mean 2^-44.73 (plus true, ~0) against 2^16
  candidate-list slots, tail (eμ/t)^t = 2^-3.9e6. The kept blocks land
  on the 2^20 generating nodes: mean 2^20.30 per node against 2^24
  receive slots, tail 2^-3.8e7.
- **C7 re-walk key repeats.** The retained key already occurs (the
  candidate matched on it), so ≥ 5 kept steps need 4 further positions
  to hit it: under H3 that is 2^-768 per 4-tuple and the L^5/5!
  ordered tuples give ≤ L^5/(5!·2^768) per chain; over 2^65 chains
  2^-451.7.
- **C8 node overflow.** Capacity 2^33.32 baby records per node against the
  mean 2^31.63; Chernoff astronomically small, 2^-20 budgeted for
  modelling error.
- **C9 base-prime supply.** The search interval [2^80, F] has length
  X = F − 2^80 solving the fixed point X = 2^66·2^20·ln(F)
  (X ≈ 2^92.0); the integrand bound 1/ln(q) ≥ 1/ln(F) gives an
  expected supply of at least X/(2^20·ln F) = 2^66 against the 2^65
  requirement. The shortfall is a Chernoff tail at δ = 1/2:
  2^-1.3e19, under H2'' (the prime-membership indicators are
  independent across the sieved primes).

Exceeding any cap makes the run output fail; the charged work is the
maximum under the caps, never the mean.

## 9. Limitations

- Analytic bound; no table was built and no discrete logarithm was
  computed at production size.
- Instruction counts are itemised from a scalar listing, not measured
  on hardware; the cluster model (no remote-memory instruction; block
  transfers priced as posting/polling plus the per-record wire copy)
  follows the accepted lineage's definition. The unit convention —
  non-group work divided by the organiser's pinned 2000-instruction
  reference budget — is the normative files' own (§3).
- H1 and H3 are validated at toy scale only (the three coverage toys,
  the L- and base-stratified exact-coverage census and the end-to-end
  run in this package's own evidence). The census's production
  extrapolation is a measured two-parameter fit (0.235/L + 0.000514)
  over a fourfold L range and a 2-to-8 base-prime range, not a proof;
  the 0.0015 allowance is 2.9× the fit's L→∞ floor.
- Memory is exponential, as for any BSGS-form Cheon attack.
- The score is 0.0523 bits under the accepted lead and 0.0205 under the
  pending fork entry; the itemisation carries a 14.7 percent margin and
  K = 0.498 leaves 0.0007 claimed success margin above the
  requirement after every cap and the census-informed allowance.

Prepared with GLM 5.3 in the omp harness.
