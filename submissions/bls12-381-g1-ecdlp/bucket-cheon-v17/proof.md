# bucket-cheon-v17: Cheon with H-prime semiprime walks and direct-mapped bucket matching

**Claim.** For the protected BLS12-381 G1 target with auxiliary powers, and
every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s
from T = [s]P0 and A_d = [s^d]P0 with heuristic success probability at
least 0.3913 (0.39 required). Its total work is
W(2^20) ≤ **2^116.749266** target-group-operation equivalents (exact
score 116749266/1000000), against the accepted track lead 116.793080 and
the pending fork entry 116.761224. Preprocessing (the whole baby side,
the walk list, the base H-primes and both gap tables — everything
independent of the challenge z — the COMPLETE cold-attack
composition: the baby records, the walk list, the base H-primes, both
gap tables, the baby side's positioning share, and ALL per-node setup
including the packed owner head words) is ≤ 2^115.7517
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
K = F(503,1000) = 0.503 exactly (raw law success 0.39529), where the
lead took K = 0.52 at +0.014 bits for comfort; (ii) the per-record
overhead charged at **150** against our own 33-line itemised listing of
**131** (baby) / **123** (giant) instructions — a 14.5 percent margin
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
effect is −0.0438 bits against the accepted lead and −0.0120 against
the pending fork entry: K 0.52→0.503 contributes about −0.037 across
the record, comparison and re-walk terms, the charge 180→150 about
−0.010, and the honest scalar-mul repricing 640→662 adds back +0.031
on the strand, segment, candidate and Lemma-1 terms (which are
subdominant, so the score moves by under 0.002 net of them).

**This revision (v13) answers the v12 review's three objections:**
(1) C2's tail is now the true failure probability: the run fails
when a slot's 64th record ARRIVES, i.e. X ≥ 64, so the tail is the
SUM of the Poisson terms from 64 upward — v12 computed the single
term exp(−λ)λ^64/64! (the probability of exactly 64), understating
the all-slots quantity 4.15052e-81 as 4.12577e-81; the run now sums
the series and reports 5.70e-81 at K = 503/1000, displayed rounded UP
as 10^-80.2;
(2) C1's worst mean now uses the C3 envelope correctly on the
comparison DRIVER: a relative deviation ε = 2^-10 on each side's
record count permits the baby × giant product to grow by
(1+ε)² = 1 + 2^-9 + 2^-20 — v12 used 1 + 2ε, omitting the ε² term;
μ_worst = μ(1+2^-9+2^-20) against the 1.01μ cap, per-node tail
2^-962.4, all-nodes 2^-878.4; (3) the success allowance is no longer
derived from the census fit — it is the WORST single-prime absolute
deviation observed among the census's thirteen L ≥ 148 primes (the
scale closest to production; measured 0.0038815 at K = 503/1000),
rounded up to 0.0039 and subtracted in full: the claim expenses the
production coverage as if production were one more census prime at
that scale drawing the worst luck of all thirteen. No extrapolation
factor is applied; the two-parameter fit (both parameters reported)
remains in the census block as trend evidence only. K moves from
0.498 to 0.503 so that the raw law 0.39529 clears the 0.39 bar with
this worst-case allowance subtracted in full.

**This revision (v17) restates C6 as a cluster model.** The v16 review caught that candidates are Cartesian products of key-sharing clusters, not independent balls — the multinomial C6 tail did not follow from H3 — and gave the example (63 babies × 1041 giants sharing a key: 65,583 candidates at one owner, the specific-set probability 2^-211,860 at a fixed owner). v17 computes the honest union over cluster shapes: the dominant channel is exactly that corner shape at 2^-93,338, with the true-class channel at 2^-1,019,000 and the withdrawn ball channel at 2^-3,884,748; every channel is far below the 2^-40 budget. The earlier consistency pass follows.
that v13's DISPLAYED bounds went stale when K moved from 0.498 to
0.503: the C2 display (4.15e-81, shown as 10^-80.3) and the C1 union
(2^-868.9) were the K = 0.498 values, smaller than the run's actual
outputs (5.70e-81 and 2^-878.4) — contradicting the claim that every
tail is rounded towards the larger probability — and the census block
in the run output still carried the retired 0.0015/2.9×-floor
allowance framing alongside the new one. v14 regenerates every
displayed bound from the run at K = 503/1000 (C2: per slot 3.43e-116,
all slots 5.70e-81, displayed 10^-80.2; C1: per-node 2^-962.4,
all-nodes 2^-878.4; C6 mean 2^-44.72) and the census block's allowance
justification now reports only the worst-observed anchor (measured
0.0038815 at L ≥ 148, rounded up to 0.0039). v15 completed the pass inside the claim, and v16 finished it by regenerating the dependency-declaration evidence text, which still quoted v10-era cap tails:
the two declared-support sentences that still quoted retired figures
(the statement's 0.0015/2.9×-floor allowance clause and the census
sign count) and the lineage narrative's v13-era numbers are
regenerated from the run. No mathematical claim changed: the same
ledger, the same K, the same allowance, the same score 116.749266.

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

**Walk.** The walk list W = the first L = 3,509,966,315,759,140
(2^51.6441, one w per chain step; L = N/2^64 exactly) H-primes above
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
6N²/m ≥ K = 503/1000, so every chain has exactly L = N/2^64 steps
and 2^64 × L = N holds as an integer identity; the realised 6N²/m is
0.503 × (1 + 2^-53), and every coverage figure below is stated at the
lower bound K = 0.503. N = 64747450334150064562029684053770240 =
2^115.6404 per side (exact integer in the run).
The wire record is 32 bytes: the 192-bit key and the 64-bit chain id;
the owner node's pool copy adds a 34-bit pool link (40 bytes total).

## 2. Matching (the lead's architecture, restated)

All 2N records are routed by three data-independent radix hops (28-bit
key fields; 2^84 owner nodes) to their owner node. Each owner node holds
a direct-mapped bucket table: 2^33 head slots indexed by 33 key bits,
each head a 5-byte PACKED word — a 34-bit pool index and a 6-bit depth
counter in 40 bits (the pool holds up to ceil(2^33.32) records, beyond
a 32-bit index; the C2 depth counter lives inside the head word, so
no separate counter array exists, and a slot holding 63 records
refuses its 64th — the run fails), a chained pool of 40-byte linked
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

## 3. The per-record itemisation (ours; charged 150)

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
(2000 + 150)/2000 = 1.075 reference operations per record, i.e. actual
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
- head word LDR (34-bit index + 6-bit depth packed in 40 bits): 1
- depth UBFX + increment + bound check (a 64th record fails the run): 3
- old-index UBFX into the pool link: 1
- new head word build (2 BFI: depth + index): 2
- pool store (5 STR: 40-byte linked record) + next-link: 5
- head word update STR: 1
**owner baby: 20**

owner node, giant lookup:
- poll + read record from the block (4 LDR): 5
- slot index UBFX ×2: 2
- head word LDR: 1
- head-word index UBFX: 1
- compare-loop base (each comparison's 8 instructions charged in the
  capped term C1): 3
**owner giant: 12**

**Listing: baby 33 + 3×26 + 20 = 131; giant 33 + 3×26 + 12 = 123.**
Charged: **150 per record** — a 14.5 percent margin over the maximum
listing, absorbing items not enumerated, the same discipline as the
lead's 26 percent over its listed 143 (charged 180). The per-record
factor is (2000 + 150)/2000 = 1.075 against the lead's 1.09 and the
sort-based dependency's 1.13. (The depth counter is packed in the
head word, so the listing pays its maintenance — extract, increment,
bound, rebuild — but no separate array is read, written, initialised
or stored.)

## 4. The ledger (exact integers; the run recomputes every term)

| term | log2 units |
|---|---:|
| records: 2N × (2000 + 150) | 116.7447 |
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
| **W** | **116.749266** |

W = 139646752555436848684353533493429875 units exactly (the run outputs
the integer and asserts the score ceil(10^6 log2 W) = 116749266).
Every term is a rational ceiling; nothing is measured, nothing is
inherited. The 47-instruction walk overhead is charged on the gap tables
and the re-walk as well, not only the main walk. Preprocessing — the
COMPLETE z-independent composition: the baby records (2^115.7447), the
gap tables, the walk list and the base H-primes (together
2^102.0806), the baby side's share of the strand and segment
positioning (2^90.3707), and the per-node setup (2^108.0342: per owner
node 2^35 instructions cover the 2^33 zeroed head-word stores, the
FIFO setup, the candidate lists and code —
the initialisation the insert algorithm reads — the 2^33 zeroed head
words carry the depth inline) — sums to 2^115.75165, claimed at its
computed ceiling 2^115.7517, and is included in W.
(v9 stated a truncated partial sum; v10 stated a partial composition
whose ceiling left only 2^99.7 of rounding slack the head-table
initialisation would exceed; this revision states the complete
composition, so the claim is a genuine cold-attack bound.)

## 5. Heuristics (disclosed, with the toy support)

- **H1 coverage.** The 6N²(1 − (L−1)/L²) independent rationals x/(u x') hit
  z like uniform elements of H: coverage 1 − e^{−K} = 0.39529 at
  K = 0.503. Support: Lemma 2's exact distinctness, the negative
  correlation within rows and columns of the pair matrix, the three
  coverage toys in the run (d = 32/64/128 toy primes, the diagonal-loss
  model K' = K(1 − (L−1)/L²), measured against predicted with
  z-scores), and — the load-bearing quantitative support — the
  L- and base-stratified exact-coverage census of §6: on twenty-eight
  genuine toy primes in six strata (twenty-two at two base primes
  with L = 74 to 296, six more at four and eight base primes), the exact
  covered set is enumerated, its density matches the same
  diagonal-corrected law with a residual whose mean absolute value
  shrinks monotonically in L, and the two-parameter fit
  |dev| ≈ 0.232/L + 0.000437 carries an explicit L→∞ floor of
  4.37e-4 with rms per-prime scatter 0.00128.
- **H2/H2' H-prime density and gaps.** A 1/2^20 fraction of primes near
  q are H-primes; gaps are exponential with mean 2^20 ln q, so the L-th
  (L = 2^51.6441) H-prime above 2^39 is below 2^78 and a gap above 2^33
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

Success claimed: 0.39529 − 2^-40 (every cap tail, each below 2^-47)
− 2^-20 (the C8 modelling budget) − 0.0039 (the census-worst
allowance) = **0.39139, claimed at 0.3913** ≥ 0.39. The allowance is
the WORST single-prime absolute deviation observed among the census's
thirteen L ≥ 148 primes — the scale closest to production; measured
0.0038815, rounded up to 0.0039 — subtracted in full: the claim
expenses the production coverage as if production were one more
census prime at that scale drawing the worst luck of all thirteen. No
extrapolation factor is applied. The census's monotone shrink across
its strata (and the two-parameter fit, both parameters reported, as
trend evidence only) indicates this overstates the production
residual; the claim does not rely on that indication.

## 6. What the toys and the census show

(All figures below are from the packaged run: `mode=all, seed=2026,
trials=120`, executed from a foreign working directory in 4.5 s.)

- **Coverage toys.** Three scalar-level toys (L = 60, chains sized for
  K = 0.503) measured against the diagonal-loss prediction 1 − e^{−K'}:

  | r | d | N/side | predicted | measured | z |
  |---|---|--------|-----------|----------|---|
  | 134219713 | 32 | 600 | 0.3974 | 0.4167 | +0.43 |
  | 536872321 | 64 | 840 | 0.3913 | 0.4333 | +0.94 |
  | 2147492353 | 128 | 1200 | 0.3974 | 0.3333 | −1.43 |

  All three z-scores lie in the ±2 band.

- **The L- and base-stratified exact-coverage census (the load-bearing
  H1 support).** For twenty-eight genuine primes r_t with
  6·d_t | r_t − 1 (d_t = 64 or 128), the census builds the SAME
  semiprime structure as production — disjoint base lists of H-primes
  times one shared walk list of L entries, with L/N sized at the
  claim's K — and then enumerates the covered set EXACTLY:
  C = {b·u/g mod r_t : b ∈ B, g ∈ G, u ∈ μ6}, all of which lies in H.
  The exact density |C|/m_t is compared with the diagonal-corrected law
  1 − e^{−K(1−(L−1)/L²)} at the realised K' (0.486 to 0.502).
  Twenty-two primes run at two base primes per side in four growing
  strata (r_t ≈ 2^25 to 2^28, L = 74 to 296); six more hold the
  modulus scale while growing the base-prime count (three at nb = 4,
  three at nb = 8).

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

  **The base-prime count, from two to 2^64 — stated honestly.** The
  census cannot run 2^64 base primes, and the extrapolation is
  supported EMPIRICALLY, not proven. What unique factorisation gives
  is only the absence of STRUCTURED collisions: distinct semiprime
  rationals never collide as rationals, so beyond the diagonal
  correction no multiplicative identity forces two ratios together
  (dependencies like (a/b)(c/d) = (a/d)(c/b) among base ratios
  create correlated VALUES, not forced collisions). The ACCIDENTAL
  mod-r collisions — at the binom(B,2)·(6L²)²/m scale across block
  pairs — are the uniform-draw assumption (H1) itself: the coverage
  law 1 − e^{−K} prices exactly these, counting unique hits among the
  6N² draws at the n²/2m birthday rate, which is the substance of
  the law, not a correction to it. The empirical control: the
  nb-strata residuals do not grow with nb —

  | nb (base primes/side) | n | mean L | mean abs dev | worst |
  |---|---|-----|------|-------|
  | 2 | 22 | 160 | 0.00218 | 0.00592 |
  | 4 | 3 | 84 | 0.00143 | 0.00217 |
  | 8 | 3 | 59 | 0.00079 | 0.00159 |

  — if anything they shrink (more base pairs average away per-prime
  coincidence luck).

  **The residual: both fit parameters reported.** Per L-stratum at
  nb = 2 (L range, mean |dev|, worst |dev|):

  | stratum | n | L | mean abs dev | worst |
  |---|---|-----|------|-------|
  | r_t ≈ 2^25 | 6 | 74–105 | 0.00331 | 0.00592 |
  | r_t ≈ 2^26 | 5 | 105–148 | 0.00250 | 0.00357 |
  | r_t ≈ 2^27 | 5 | 148–210 | 0.00171 | 0.00388 |
  | r_t ≈ 2^28 | 6 | 210–296 | 0.00118 | 0.00209 |

  The mean absolute deviation shrinks monotonically as L grows,
  across a fourfold L range, with signed means at +0.0013 or below in
  every stratum (no systematic bias; overall 18 positive, 10
  negative). The least-squares fit of |dev| against 1/L over the
  twenty-two nb = 2 primes is **|dev| ≈ 0.232/L + 0.000437 — both
  parameters reported**: the slope dominates at toy L (it predicts
  0.0036 at L = 74 and 0.0012 at L = 296, matching the strata), and
  the intercept is the model's L→∞ floor, 4.37e-4, with rms per-prime
  scatter 0.00128 around the fit. The fit is TREND EVIDENCE ONLY.

  **The allowance: the worst observed, subtracted in full.** The
  claim's success allowance is NOT derived from the fit. It is the
  worst single-prime absolute deviation observed among the thirteen
  census primes at L ≥ 148 — the scale closest to production —
  measured 0.0038815 at K = 503/1000 and rounded up to **0.0039**,
  subtracted in full: the claim expenses the production coverage as
  if production were one more census prime at that scale drawing the
  worst luck of all thirteen. No extrapolation factor is applied and
  no model of the production residual is asserted; the monotone
  shrink visible in the strata indicates the allowance overstates
  the production residual, and the claim does not rely on that
  indication. (The allowance history: v9 0.0008 with a wrong regime
  statement; v10/v11 0.0015 tied to a fit whose intercept v10
  discarded; this revision anchors the allowance to the worst
  observed deviation at the closest scale and lifts K from 0.498 to
  0.503 so the claim clears the 0.39 bar with that allowance
  subtracted in full.)

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
× 5-byte PACKED head words — a 34-bit pool index and a 6-bit depth
counter in 40 bits, so the C2 depth counter lives inside the
already-counted head word and no separate counter array exists —
= 2^35.7; record pool ceil(2^33.32) =
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

- **C1 comparisons per node.** Cap 1,317,460,454 = 1.01 × the mean
  1.3044 × 2^30 per node. The comparison count of one node is a sum of
  bounded increments — each giant lookup contributes the occupancy of
  the one head slot it reads, an increment in [0, 63] (the C2 cap; a
  64th record in the slot fails the run) — so Bernstein's inequality
  for bounded increments applies with M = 63 and variance bound 63μ:
  P(X − μ ≥ t) ≤ exp(−t²/(2(63μ + 63t/3))). The tail is evaluated at
  the WORST mean inside the C3 envelope (see C3): a relative deviation
  ε = 2^-10 on each side's count permits the baby × giant product to
  grow by (1+ε)² = 1 + 2^-9 + 2^-20, so μ_worst = μ(1+2^-9+2^-20)
  against the cap 1.01μ, i.e. t ≈ 0.008μ (v12 used 1 + 2ε, omitting
  the ε² term): the per-node tail is 2^-962.4 and the union over all
  2^84 nodes is **2^-878.4**. (v8 priced C1 as a unit-Bernoulli
  Bernstein exp(−μδ²/(2+δ)), which ignored the bounded increments;
  v9-v12 carried the bounded form with an incomplete envelope; this
  revision carries the correct form and the full envelope.)
- **C2 head slot depth ≤ 63.** Poisson(0.390) occupancy; the run fails
  when a slot's 64th record ARRIVES, so the tail is P(X ≥ 64) — the
  SUM of the Poisson terms from 64 up (v12 computed the single term
  at 64, understating it): per slot 3.43e-116, over all 2^117 slots
  5.70e-81, displayed rounded UP as **10^-80.2** (2^-266.3). The depth
  counter is the 6-bit field of the packed head word; a separate
  counter array would add 2^33 bytes per owner node and is NOT
  allocated.
- **C3 COMBINED node count deviations.** The 2N records land
  multinomially on the 2^84 nodes; the COMBINED count (baby + giant)
  has mean 2N/2^84 = 2^32.71 per node, and a relative deviation
  > 2^-10 has Chernoff tail 2^-3070 per node, 2^-2986 over all nodes
  (the run derives 2^-3070.34 and 2^-2986.34 at the combined mean;
  both displayed rounded towards the larger probability). The
  deviation's COST is not charged here and no per-side event is
  claimed: C1 consumes the combined envelope through AM-GM
  (B·G ≤ ((B+G)/2)² ≤ (1+ε)²E[B]E[G], the sides having equal means),
  so one combined event bounds the comparison driver, the cap bounds
  the shifted mean itself, and C3 prices only that envelope's breach.
- **C4 candidate pairs ≤ 2^44.** A compared pair already agrees on the
  84 routing key bits and 33 slot key bits; under H3 the 75 unchecked
  bits give a false full-key match at rate 2^-75 per comparison. With
  N²/2^117 comparisons the expected false candidates are N²/2^192 ≈
  2^39.28; the count is Poisson-concentrated at that mean, and the
  Chernoff tail (eμ/t)^t at t = 2^44 is 2^-5.7e13.
- **C5 true matches ≤ 2^20.** Two failure channels. (i) Off-diagonal:
  the count is Poisson(0.503)-concentrated; the tail at 2^20 is
  (eμ/t)^t = exp(−1.42e7) = 2^-2.0e7. (ii) Diagonal challenges: for a
  baby base a and giant base a' with z = a/a' (over the six twists),
  every one of the L shared-w positions is a true match, so the run
  aborts at the caps (the same single event that would break C4 and
  C7); the union over all 2^64 × 2^64 base pairs and six twists has
  probability 6·2^128/m = 2^-104.27. Both channels are far below
  2^-40 and are charged inside the common 2^-40 cap deduction. The
  owner's 2^16-slot candidate list binds before the 2^20 budget on
  the true channel: the true class's product passing 2^16 pairs has
  tail (eμ/t)^t = 2^-1,019,000 (C6's true-class channel); the 2^20
  budget remains the charge cap for the re-walk and
  scalar-multiplication terms.
- **C6 candidate sinks, CLUSTER MODEL.** Candidates are NOT
  independent balls: a full-key class of b baby records and g giant
  records at one owner contributes the Cartesian product b·g
  candidates there. The owner list (2^16 slots) overflows through
  key-sharing clusters, and the v12-v16 multinomial ball model
  (2^-3.88e6 per owner) is WITHDRAWN. Under H3 a specific (b+g)-set
  shares one full 192-bit key with probability 2^-(84 + 192(b+g−1))
  at a fixed owner (the form of the v16 review's example, 63×1041
  records giving 2^-211,860); summing over all 2^84 owners cancels
  the 84, so the expected number of such classes anywhere is
  C(N,b)·C(N,g)·2^-192(b+g−1) — already a union over all record sets
  and owners. With b capped at 63 by C2, the dominant overflowing
  shape is the largest b with the smallest overflowing g: (63, 1041),
  b·g = 65,583 ≥ 2^16, at **2^-93,338** (the run computes the
  exponent term-wise; lgamma differences at N = 2^115.6 cancel
  catastrophically in float). The true-match class (the challenge's
  own key, Poisson(0.503) pairs) overflows the list only past 2^16
  pairs: (eμ/t)^t = 2^-1,019,000. The withdrawn ball-concentration
  channel was 2^-3,884,748 — the cluster model is the honest, weaker,
  still-binding bound, displayed **2^-9.2e4**, and every channel
  sits far below the 2^-40 budget. The kept-block sink at the 2^20
  generating nodes (mean 2^20.30 per node against 2^24 receive
  slots) keeps its multinomial bound 2^-3.8e7: the kept blocks are
  spread one-per-candidate-chain across nodes (a class of b babies
  and g giants sends at most max(b, g) ≤ 2^16 blocks to any one
  node), and the extreme cluster channel there is 2^-1.66e9, below
  the concentration bound.
- **C7 re-walk key repeats.** The retained key already occurs (the
  candidate matched on it), so ≥ 5 kept steps need 4 further positions
  to hit it: under H3 that is 2^-768 per 4-tuple and the L^5/5!
  ordered tuples give ≤ L^5/(5!·2^768) per chain; over 2^65 chains
  2^-451.6 (the run derives 2^-451.69, displayed rounded towards the
  larger probability).
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
  extrapolation is empirical, not a proof: the two-parameter fit
  (0.232/L + 0.000437, both parameters reported) is trend evidence
  over a fourfold L range and a 2-to-8 base-prime range, and the
  0.0039 success allowance is the worst single-prime deviation
  observed at the closest scale (L ≥ 148), subtracted in full rather
  than extrapolated.
- Memory is exponential, as for any BSGS-form Cheon attack.
- The score is 0.0438 bits under the accepted lead and 0.0120 under the
  pending fork entry; the itemisation carries a 14.5 percent margin and
  K = 0.503 leaves 0.0013 claimed success margin above the
  requirement after every cap and the census-worst allowance.

Prepared with GLM 5.3 in the omp harness.
