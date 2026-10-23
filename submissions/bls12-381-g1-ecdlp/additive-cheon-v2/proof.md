# Additive Cheon: one addition per baby/giant step, in a fixed-depth key store

**Claim.** For every n with 2^20 ≤ n < r, the deterministic algorithm below recovers s with heuristic
success probability 0.4055 (≥ 0.39 required). Its total work is W(n) = W(2^20) ≤ **2^116.889745**
target-group-operation equivalents.
- Preprocessing, meaning the challenge-independent work, is ≤ 2^115.898 and is included in W. It
  consists of the baby walk with its key store, the per-prime tables, the characters and the sieves.
- Peak memory is ≤ 2^122.88 bytes. There is no advice.
- The algorithm reads only A_1 = T and A_d, with d = 2^20.

The idea is Cheon's (p−1)-attack, run as a baby-step/giant-step (BSGS) inside the order-m subgroup H of
d-th powers, with one change. The baby and giant sets are **small integers that lie in H**, not geometric
progressions ζ^i. Small integers listed in increasing order can be walked with **one batched affine
addition per element**. A geometric step needs a scalar multiplication, which is 2 additions with the
Kozaki–Kutsuma–Matsuo (KKM) c = 3 tables.

This document is a revision of the published package of a submission reviewed valid at 2^116.982765
(2026-10-23T14:42Z). The algorithm, coverage argument, toy evidence and reference-addition convention are
unchanged. What changes is §4: the key store is one flat slot array in a single extended-memory heap, its
two page-table levels carry two-word entries so every stored pointer holds the bits it must hold, and
every access charges a **bounded** number of probes and the full multiword address arithmetic, so no part
of W is an average and no address is narrower than the space it indexes. Each change is listed in §4.6.

## 1. Facts about the fixed target

These are checked by `evidence.py` and by direct computation.
- r ≈ 2^254.857 is prime and v₂(r−1) = 32, so d = 2^20 divides r−1.
- H = {x^d : x ∈ F_r^*} = ker(x ↦ x^m) has order m = (r−1)/d ≈ 2^234.857, and 6 | m.
- **The hidden value is uniform on H.** t := s^d is uniform on H, because s is uniform on F_r^* and
  x ↦ x^d is d-to-1 onto H. A_d = [t]P0 is public.
- **The symmetry group sits inside H.** Let U = {±1, ±λ, ±λ²}, where λ = x² − 1 is the GLV eigenvalue
  (λ² + λ + 1 = r) with [λ](x, y) = (βx, y). U ⊂ H because its order 6 divides m.
- **One key identifies a whole orbit.** All six points [u]Q, for u ∈ U, have y-coordinate ±y(Q),
  because λ fixes y. Conversely, the points with y-coordinate ±y₀ are exactly the six (βʲx₀, ±y₀). So
  canonical(Q) := min(y, p − y) identifies the U-orbit of Q.
- **Membership in H is a character test.** χ(x) := log_ζ(x^m) ∈ Z/d is a character with kernel H, for
  a fixed generator ζ of μ_d. So a·b ∈ H ⇔ χ(a) + χ(b) ≡ 0 (mod d).
- **2 ∉ H.** Indeed 2^m ≠ 1 (mod r); 2^m has order 2^16.
- **The store's slot index is uniform.** The canonical value c = min(y, p−y) of a point with y uniform
  on [0, p) takes each value in [1, (p−1)/2] with probability 2/p, and the value 0 with probability
  1/p. Each half-open interval [j·2^263, (j+1)·2^263) below (p−1)/2 therefore contains exactly 2^263
  admissible c, so j = ⌊c/2^263⌋ is uniform on j ∈ [1, S−1] with S = ⌊((p−1)/2)/2^263⌋ + 1 =
  2^116.700672, and j = 0 carries the single canonical value 0, whose weight 1/p is negligible against
  the 2^263·2/p of every other slot. This is the slot index of §4.2; the store's load is therefore
  0.50584 and not an assumed figure.

## 2. Algorithm

Parameters: K = 0.52 as the coverage exponent, a discard allowance δ = 0.01 for the store (§3, H3),
X1 = 2^76, X2 = 2^70 < X1/2, GMAX = 2^30, N_B = ⌈√(K·m·c_G/(6(1−δ)c_B))⌉ ≈ 2^115.717 and
N_G = ⌈K·m/(6(1−δ)N_B)⌉ ≈ 2^115.626, with c_B, c_G the per-element costs of §4.

1. **Precompute** the characters, sorted into class lists.
   - Compute χ(q) for the primes q ≤ X2 and for the primes used in (X1/2, X1]. Each costs one 235-bit
     exponentiation in F_r and a Pohlig–Hellman logarithm in μ_{2^20}.
   - Sieve χ(b) for all b ≤ X2 by multiplicativity.
   - Stream each b into its class list L_c = {b ≤ X2 : χ(b) = c}, sorted in increasing order.
2. **Baby walk.** For fresh primes a ∈ (X1/2, X1]:
   - set Q_a = [a]P0 and build T[k] = [k]Q_a for k ≤ GMAX (T[2] by doubling, then T[k] = T[k−1] + Q_a);
   - walk b over L_{−χ(a)} in increasing order, with [ab]P0 = [ab_prev]P0 + T[b − b_prev];
   - for a gap greater than GMAX, first add T[GMAX] repeatedly;
   - insert the key canonical([ab]P0) with locator (a, b) into the key store of §4.2.

   Stop after N_B elements. Every element a·b lies in H. They are distinct integers below r, because a
   is the unique prime factor of a·b above X2. Two elements can share a key only if
   a′b′ ≡ u·ab (mod r) with u ∈ U \ {±1}: for u = ±1 both products lie in (0, 2^146) ⊂ (0, r/2), so
   equality of keys forces the same integer and the same pair, while for u = ±λ, ±λ² the congruence
   λ·ab mod r lands on a given a′b′ < 2^146 with probability 2^-108.9, i.e. N_B · 2^-108.9 ≈ 2^6.8 times
   over the whole baby set. Such a pair is two distinct points of one U-orbit, which is exactly the
   relation the collision test accepts; the second record lands in the next free slot of the same
   window, and both locators are valid, so the effect is a 2^-108.9 loss of coverage, not an error.
   Inserts therefore never test for duplicate keys.

   A key whose window and whose block's overflow area are both full is **discarded**: it is not stored,
   and it adds no coverage. At most a δ = 0.01 fraction of the N_B inserts is lost this way (§3, H3).
3. **Giant walk.** Do the same with disjoint primes a′ and base A_d, which gives the points [a′b′·t]P0.
   Look up the canonical key of each point. Stop after N_G elements, or earlier on success. Lookups
   insert nothing, so nothing is discarded on this side.
4. **Collision.** A key match means [g·t]P0 = [u·ab]P0 for some u ∈ U, where g = a′b′. Test the six
   candidates t = u·ab/g by checking [t]P0 = A_d.
5. **Recover s.** Take one d-th root s0 of t: twenty square roots, keeping at each step the root that is
   still a 2^j-th power. Then s = s0·ω^i for some i < 2^20, where ω generates μ_d. Test [s0·ω^i]P0 = T
   for each i and output s.

**Batching and representation.** Each prime's walk, and its table build, is sequential. The 2^16
primes in flight are therefore advanced in lockstep. Every round performs 2^16 independent incomplete
affine additions, one per prime, sharing a single Montgomery batch inversion; this is the reference
operation's batched setting. The inversion costs ≈ 100 units per round, about 0.0015 units (3
instructions) per element, which fits inside the walk-bookkeeping budget. All stored points, meaning
tables, walk states and Q_a, are affine: 96 bytes each.
- Q_a = [a]P0 (or [a]A_d) is computed in Jacobian coordinates and normalised by one inversion, charged
  in the per-prime scalar budget.
- T[2] is one affine doubling, also charged there.

**Incomplete additions.** A walk step adds [b_prev]Q and [gap]Q, with Q ≠ O and b = b_prev + gap < r. The
x-coordinates can match only if b_prev ≡ ±gap (mod r).
- b_prev ≡ −gap would force b ≡ 0, which is impossible.
- b_prev = gap means b = 2·b_prev in the same class list, which would need χ(2) = 0. Since 2 ∉ H, this
  never happens.
- *Overflow chunks.* [b_prev + i·GMAX]Q + [GMAX]Q collides only when b_prev = GMAX and i = 0. That is at
  most once per prime; it is done by a doubling, charged to the per-prime scalar budget. The remainder
  addition is safe because b_prev + j·GMAX > GMAX ≥ rem for j ≥ 1.
- *Table build.* T[k] = [k−1]Q + Q is safe for k ≥ 3.
- Giant points are multiples of Q = [a′t]P0 ≠ O, so the same argument covers them.

## 3. Success probability

Steps 4–5 are deterministic. Let B_kept be the baby keys the store placed, of size N_B(1−δ) at least.
The algorithm succeeds exactly when t ∈ S := U·B_kept·G^{-1}. Because t is uniform on H,
Pr[success] = |S|/m.

- **H1 (coverage).** Apart from the structural loss below, the 6·N_B·N_G products u·b·g^{-1} coincide no
  more often than random products would. Then |S|/m ≈ 1 − e^{−(1−δ)K} with (1−δ)K = 0.52, i.e. 0.4055.
  Two kinds of coincidence are ruled out exactly:
  - A pair with the same baby prime a and the same giant prime a′ that coincides with another such pair
    would need an integer identity below r, which is impossible even with u = −1.
  - Products built from different primes cannot be equal as rationals, because a is determined.

  *Structural loss.* Rationals repeat when b/b′ = c/c′, which happens through a shared factor
  g = gcd(b, b′) and a cofactor k with χ(k) = χ(g). The number of repeats per pair is at most
  2·E[gcd]/d ≈ 2·(6/π²)(ln X2 + 0.15)/d ≈ 2^−14.1, where only this one quantity — a property of the
  integer factorisation of the b's, which does not enter W — is an average. That changes the success
  probability by less than 10⁻⁴. The margin is (1−δ)K = 0.52 against the 0.494 needed for 0.39.
- **H2 (equidistribution).** The values χ(a) of the primes used are spread evenly over Z/d. So the
  average length of the class lists actually walked is ≈ X2/d = 2^50; the lists' total length is
  exactly X2. This affects only the number of primes needed: at the average class length X2/d = 2^50,
  N_B + N_G elements need about 2^66.672 primes, and the ledger charges N_A = 2(N_B + N_G)d/X2 = 2^67.672,
  twice that, as the H2 allowance. There are 2^69.3 primes available in (X1/2, X1], so the charged count
  still has 3× slack. H2 does not affect success.
- **H3 (store truncation).** The store places a key in its 8-slot window if any of those slots is free,
  else in one of its block's 4 overflow records, else discards it (§4.2). Every access is charged at
  the worst case, 8 window probes plus 4 overflow comparisons, so no average probe count enters W. The
  only quantity that enters is the fraction of keys discarded, charged at δ = 0.01 and paid for by
  raising the coverage exponent from 0.52 to 52/99, so that the exponent actually realised is
  0.52·(1 ± 0) exactly:
  - 6·N_B·N_G/m ≥ 52/99 = 0.525253... by the choice of N_B and N_G, so
    (1−δ)·6·N_B·N_G/m ≥ 0.52 and Pr[success] ≥ 1 − e^{−0.52} = 0.40548.
  - Evidence at production load: on 2^19-slot tables at the production load 0.50584, the fraction of
    inserts whose whole 8-slot window is occupied is 0.005149, and with the 4 overflow records the
    discarded fraction is 0.000758 = 2^-10.37 (mean of 3 seeds). The charge 0.01 is 13.2× that
    measurement and 1.94× the window-full fraction, so it also covers the case where the overflow area is
    full everywhere and the overflow absorbs nothing. The rule is local, so the fraction depends on the
    load, D and M, not on the slot count; the measurement is made at D = 8, M = 4, load 0.50584, which
    are the production values.
  - Discarded keys only remove coverage. A lookup finds every key that was stored: it probes the same
    eight window slots in the same order and then the same four overflow records, and no slot is ever
    freed.

Evidence: `evidence.py` with the declared inputs. These results hold at toy scale only.

- **Coverage.** Exact |S|/m for a toy prime r ≈ 2^28 with d = 2^8, K = 0.52 and prediction 0.4056:
  - random sets: 0.4061;
  - smallest integers in H: 0.4069;
  - product construction at X2/d = 4, 16, 64 and 256: 0.4038, 0.4020, 0.4028 and **0.4056**.

  The measured shortfalls from the prediction are 0.0018, 0.0036, 0.0028 and 0.0000. They are small
  (≤ 0.0036, well inside the 0.0155 margin) and shrink as X2/d grows. The gcd formula is only a
  conservative upper bound on the loss: it is not a model of these deficits, and it grows with ln X2.
  The class-length spread also narrows as X2/d grows (0.76–1.22 × mean at X2/d = 256). Production has
  X2/d = 2^50.
- **Store.** The fixed-depth window rule measured at the production load, as listed under H3: window-full
  fraction 0.005149, discarded fraction 0.000758 with M = 4 (0.002205 with M = 2, 0.005149 with no
  overflow).
- **Toy attack.** The whole algorithm on a 24-bit j = 0 curve, with d = 2^6 and GMAX = 16, so that the
  overflow path runs (1.3·10⁶ chunk additions; the doubling case b_prev = GMAX did not occur).
  - 1,500 trials gave success 0.3953 against 0.4095 predicted, about 1.1σ low, at X2/d = 16 where a
    small toy-scale loss is plausible. The `conservative_gcd_bound_toy_scale` field in the output (0.366)
    applies at toy scale only; at production scale the same bound is 2^−14.1.
  - Every recovered scalar was correct.
  - Earlier runs of 500 trials on seeds 1–3 gave 0.404, 0.420 and 0.450.
  - Its key table is an exact dictionary, so it tests the walk and the coverage argument; the store's
    window rule is measured separately, above.

## 4. Cost (unit: one batched incomplete affine G1 addition = C = 2,000 instructions)

Cost is H + R/C, where H counts reference additions and R counts every other scalar AArch64 instruction
once: addressing, loads and stores, comparisons, control. Nothing is charged twice and nothing is
omitted: the reference unit is the 6M batched incomplete affine addition of the protected convention, and
its internal arithmetic is never re-charged.

### 4.1 The walk

The two walks perform N_B + N_G ≈ 2^116.672 reference additions, one per element, and each element also
carries the per-element instruction budget of §4.3. So the walk work is N_B·c_B + N_G·c_G, with
c_B = 1 + (26+29+106+80)/2000 = 1.1205 and c_G = 1 + (26+29+253+80)/2000 = 1.194.

### 4.2 The key store: a fixed-depth window with a bounded overflow

The baby side stores N_B keys, of which at least N_B(1−δ) are placed. The store is a single flat array
of blocks, each holding 1024 primary slots followed by 4 overflow records: 1028 records of 72 bytes,
74,016 bytes per block, an exact multiple of 72, so no record straddles a block. There are
⌈S/1024⌉ = 2^106.700672 blocks, 2^122.876 bytes in one extended-memory heap. A slot holds the 9-word
record: the key, which is the six limbs of the canonical y in words k0…k5, and the locator (a, b) packed
into three words (a ≤ 2^76 and b ≤ 2^70 take 146 bits). Bit 63 of k5 is the occupancy flag: it is free
because c < p/2 < 2^380.71, so k5 < 2^60.71. A slot is free iff that bit is 0. All six points of a
U-orbit share the key, and the key comparison is exact.

**Slot index.** j = ⌊c/2^263⌋ = the top 118 bits of c. By §1, j is uniform on [1, S−1] with a negligible
weight at 0, and the load after N_B inserts is N_B/S = 0.50584. No mixing or modulo reduction is needed:
S is by construction the number of distinct prefix values, so there is no bucket skew. The home slot j
lies in block P = j ≫ 10 at in-block offset o = j mod 1024.

**Fixed-depth probing, bounded overflow, and the truncation rule.** An insert probes exactly the eight
slots (o + t) mod 1024 of its own block, t = 0, 1, …, 7, in that order — the window wraps to the block's
start rather than crossing into the next block, so every probe of an access stays under one page-table
translation. The first slot whose flag is 0 takes the record. Duplicate keys are not tested for, since a
repeat can only be a U-orbit partner and both locators are valid (§2). If all eight slots are occupied
the record goes to the first free one of the block's four overflow records, at in-block offsets
1024…1027. **If those four are occupied too, the record is discarded** and its element adds no coverage;
§3 H3 charges at most δ = 0.01 of the inserts that way and shows success still clears 0.39 with margin.
A lookup probes the same eight slots in the same order, comparing the full 381-bit key at each occupied
one, and then, if no slot matched, the block's four overflow records; nothing is ever freed, so a key
that is in the store is found.

*Worst case per access, and why no average appears in W.* The probe sequence is data dependent, so the
charge is the worst case over all inputs: **eight window probes plus four overflow comparisons, twelve
slot visits, every time**. An insert tests only the flag word of each visit (a store needs no key
comparison) and writes one record; a lookup compares the full key at each of the twelve visits. Because
the bound is fixed by D and M and not by the occupancy pattern, W contains no probe count that depends
on the key set; the occupancy enters only through δ, which is charged and paid for by raising K.

**Address translation.** A slot address is 123 bits wide, because the heap is 2^122.876 bytes, so neither
a block's byte address nor a level-2 table's byte address fits in one 64-bit word. Both page-table levels
therefore carry **two-word (16-byte) entries**: 128 bits per entry, of which the pointer field uses 123
bits at level 2 (5 spare) and 111 bits at level 1 (17 spare). The block number P = j ≫ 10 is 107 bits,
split 54 + 53 across the two levels, so each index is at most 54 bits and fits in one register. The
convention counts instructions, and a load whose address is held in a register pair costs one instruction
per word it fetches; what the multiword arithmetic costs — the second word of each entry, the carry into
the high word, the shift/or that builds each index — is itemised here and in §4.3 and charged, so no
address is charged as narrower than it is:

- the block number P = j ≫ 10, 107 bits in two words (4 instructions; the slot index j itself is 118 bits
  in two words, built by shift/shift/or per word, 6);
- the level-1 index P ≫ 53 = (P_hi ≪ 11) | (P_lo ≫ 53) (3); its entry address T1base + 16·index (shift,
  add, 2) and the entry loaded as two words (2). T1 has 2^54 entries of 16 bytes, 2^58 bytes, and its own
  base is a 58-bit address, so T1 is named by a single 64-bit pointer. The entry holds the base of one
  level-2 table: 111 bits, 17 spare;
- the level-2 index P mod 2^53 (1); its entry address is the two-word table base + 16·index (shift, add,
  adc, 3) and the entry is loaded as two words (2). Level 2 needs ⌈2^106.700672/2^53⌉ = 2^53.700672
  tables of 2^53 entries of 16 bytes, 2^57 bytes each and 2^110.700672 bytes in total. The entry holds
  the block's base address 74,016·P < 2^122.876: 123 bits, 5 spare;
- the slot's offset in its block, o = j mod 1024 (1) and 72·o (3), and the record address, the two-word
  block base plus that offset (add, adc, 2).

The two-word record address is then the operand pair of the record's own loads and stores, and each of
the twelve slot visits advances it by 72 with add and adc (2 instructions) plus the in-block wrap (a
compare and, when it wraps, the two-word reset of the address to the block base, charged at 3 per visit).

The page tables are ordinary memory: 2^58 + 2^110.700672 ≈ 2^110.701 bytes, charged at 2^110.71 in the
peak memory of §4.5 and zeroed once at start-up (§4.4: 2^55 + 2^108 words). A dependent load chain of
depth two is what the machine's own address translation does; the difference is only that this table is
sized for 2^122.876 bytes instead of 2^48. There is no store-and-forward: the record is not copied
between nodes, because the heap is not distributed. Each level brings in two words, and the record brings
in its own nine; those words plus the multiword index arithmetic itemised above are the complete transfer
and addressing cost of one access. The overflow records of a block sit at a fixed offset from the same
base, so they need no second translation: only the offset construction (add, adc, 2 instructions) that
the itemisation charges.

**Why the store may be one heap.** The protected convention counts scalar instructions, not capacity,
bandwidth or elapsed time, and it gives memory no contribution and no tie-break; a data structure must
account for its implementation, which is what this section does. Nothing in the track rules requires a
2^122.876-byte store to be spread over processors, and spreading it would add store-and-forward work that
the algorithm does not need: the published package's 12 hops at 18 instructions were the price of a
sharded store it had chosen, not a property of Cheon's attack. What the implementation must do is
address, transfer and compare, and all three are charged here: the two-word index arithmetic, the
two-level translation with its four dependent entry words, the twelve slot visits, the nine record words
and the full key comparison, with the page tables counted in peak memory and zeroed in the ledger. The
store's size changes the depth of the translation and the width of every address, both of which are
charged, and nothing else.

### 4.3 Per-element instruction budgets

In instructions on top of the reference addition. Each budget is the itemised count, plus a small
slack.

| Block | Instructions | Itemised |
|---|---:|---|
| canonical y, both walks | 26 | p − y (7), 6-limb compare (7), 6-limb select (6), move the key words into registers (6) |
| addressing and page-table walk, both paths | 29 | slot index j, 118 bits in two words (6); block number P = j ≫ 10 in two words (4); level-1 index built as (P_hi ≪ 11) \| (P_lo ≫ 53) (3); level-1 entry address T1base + 16·index (2); level-1 entry, two words (2); level-2 index P mod 2^53 (1); level-2 entry address over a two-word table base (3); level-2 entry, two words (2); in-block slot offset o = j mod 1024 (1) and byte offset 72·o (3); record address, two-word block base plus 72·o (2) |
| insert slot visits, baby | 104, budgeted 106 | 8 window probes at 8 each = 64 (two-word advance by 72 and in-block wrap, 5; flag word load, 1; flag test, 1; branch, 1), overflow base over the two-word block base (2), 4 overflow flag probes at 5 each = 20 (two-word advance, 2; flag load, 1; flag test, 1; branch, 1), store the 9-word record (9), class-list advance and gap handling (9) |
| lookup slot visits, giant | 251, budgeted 253 | 8 window probes at 21 each = 168 (two-word advance by 72 and in-block wrap, 5; flag test, 1; 6 key words loaded and compared with branches, 15), overflow base (2), 4 overflow records at 18 each = 72 (two-word advance, 2; flag test, 1; 6 key words compared, 15), class-list advance and gap handling (9) |
| walk bookkeeping, both walks | 80 | read b and the gap, overflow test (6); address T[gap] and load the 96-byte affine point (15); load and store the 96-byte walk state (24); the batch-inversion prefix and suffix buffers, memory traffic distinct from the 3M of arithmetic inside the reference unit (24); round control, counters, loop (7); the rare chunk path above GMAX (4, amortised since Σ gaps ≤ X2 per prime) |
| **per element** | **241 baby / 388 giant** | 26 + 29 + 106 + 80 and 26 + 29 + 253 + 80 |

Twelve slot visits appear in both access paths: eight window probes and four overflow records, the
worst case of the rule in §4.2. An insert visits each slot only for its flag word, an insert needs no key
compare, and a lookup compares the full key at every visit; neither path contains a term that depends on
how many slots happen to be occupied. Each visit carries the two-word advance of the record address and
its in-block wrap, charged at the worst case on every visit. These give c_B = 1 + 241/2000 = 1.1205 and
c_G = 1 + 388/2000 = 1.194.

### 4.4 The ledger

Recomputed exactly by `evidence.py` (which prints every term and the exact rational score):

| Term | log₂ (units) |
|---|---:|
| baby walk N_B·c_B | 115.882 |
| giant walk N_G·c_G | 115.882 |
| gap-overflow additions (deterministic: Σ gaps ≤ X2 per prime) | ≤ 107.672 |
| per-prime tables T[k] (2^30 additions × 2^67.672 primes) | 97.672 |
| per-prime scalar multiplications (3 × 1,000 units per prime) | 79.223 |
| χ for primes ≤ 2^76 (100 units = 2·10⁵ instructions each) | 77.258 |
| prime sieve to X1 (64 instructions per integer); χ-sieve and class lists to X2 (200 per b) | 71.103 |
| store and page-table initialisation (1028 records × 2^106.700672 blocks × 9 + 2^108 + 2^55 page-table words) | 108.911 |
| collision checks, d-th root, μ_d search, input checks | ≤ 31.551 |
| **total W(2^20)** | **116.889745** |

- A scalar multiplication is charged at 1,000 units for a ≤ 76-bit scalar and 3,000 units for a 255-bit
  scalar. These are generous Jacobian ladders at ≤ 5 units per step, plus one inversion at 100 units.
- log₂ W = 116.889744940. The score is rounded up: 116.889745 = 23377949/200000. `evidence.py` proves
  the rounding with exact integers: writing W = 2^116·f with 1 ≤ f < 2 and g = ⌊f^1000·2^40⌋, the
  integers satisfy 2^929744 < g^1000 and (g+1)^1000 ≤ 2^929745, so 2^116889744 < W^10^6 < 2^116889745
  and the sixth decimal is a ceiling, not an estimate.
- Every record's nine words are written once to zero the flags before the walk — 1028 per block,
  2^106.700672 blocks — plus 2^108 page-table words for level 2 and 2^55 for level 1 (two words per
  16-byte entry), and all of it is charged rather than assumed free.
- The store's discard allowance is inside K: N_B and N_G are chosen against K/(1−δ) = 52/99, so the
  ledger above already pays for it.

### 4.5 Memory

| Component | Bytes |
|---|---:|
| block array, 2^106.700672 blocks × 74,016 | 2^122.876 |
| page tables, 2^58 (level 1) + 2^110.700672 (level 2), two-word entries | 2^110.701 |
| class lists, X2 × 9 | 2^73.2 |
| walk state and affine tables for 2^16 primes in flight, 2^16·2^30·96 | 2^52.6 |
| **peak, sum** | **≤ 2^122.88** |

Memory is not scored and carries no tie-break. Any BSGS form of Cheon stores one side, about 2^115.72
entries; the four overflow records per block add 2^−8.0 of the array and the two-word page tables
2^−12.2. The sum is 2^122.8765.

### 4.6 What changed from the published package

1. **Store organisation and address width.** The published package sharded the table over 2^60 nodes
   and charged a radix-32 store-and-forward of 18 instructions per hop for 12 hops. That organisation is
   not a property of the algorithm: the keys can equally be held in one address space. This package uses
   one flat block array addressed by a two-level radix page table whose entries are **two words (16
   bytes) at both levels**, so the level-2 entry really holds a 123-bit block base and the level-1 entry a
   111-bit level-2 table base; the two-word index build, the second word of each entry and the record's
   two-word address and two-word advance are itemised in §4.3. That is 29 instructions of addressing and
   the record's own 9 words per access, against 12 hops of 9+9 words plus the final read, and the page
   tables are counted in memory and zeroed in the ledger.
2. **Bounded probing.** The published package charged the *average* probe count of linear probing (and
   this revision's predecessor the same load-factor average, 2.42 probes per access), which does not bound
   a data-dependent probe sequence. Both are replaced by the fixed-depth window of §4.2: eight probes, four
   overflow records, and an explicit discard rule whose fraction δ = 0.01 is charged and paid for by
   raising the coverage exponent from 0.52 to 52/99. No probe count in W is an average any more.
3. **Budgets are itemised counts.** The published per-element envelope was 40 + 300 + 150 = 490
   instructions; here it is 26 + 29 + 106 + 80 = 241 for a baby element and 26 + 29 + 253 + 80 = 388 for a
   giant one, with the shared addressing row, the twelve slot visits of the bounded rule and the two-word
   page-table entries itemised inside those rows.
4. **Unbalanced sets.** With c_G > c_B the optimal split is N_B > N_G by the ratio c_G/c_B, which the
   ledger solves exactly; the published package used N_B = N_G.
5. **Exact rounding.** The score is computed from the exact rational W and an exact integer proof that
   W^10^6 lies strictly between 2^116889744 and 2^116889745, so the upward rounding at the sixth decimal
   cannot slip.
6. **Everything else is unchanged**: K = 0.52, X1, X2, GMAX, the class-list construction, the character
   test, the walk, the collision test, the root recovery, the toy experiments and their results.

## 5. W(n) and applicability

The algorithm uses only A_1 and A_{2^20}, so W(n) = 2^116.889745 for every 2^20 ≤ n < r. Nothing is
claimed for n < 2^20. More generally, for any d | r−1 with 6 | (r−1)/d and d ≤ n:
- W ≈ 2·√(c_B·c_G)·√(K(r−1)/(6(1−δ)d)), plus the lower-order terms above;
- X2 = d·2^50.

d = 2^20 is the largest divisor of r−1 that does not exceed the supplied n = 2^20 (2^19·3 = 1,572,864
exceeds it), so this is the best d the powers sequence allows.

## 6. Context (unscored estimates)

| Method | log₂ W |
|---|---:|
| Pollard rho on G1, order-6 classes: 0.5·log₂(πr/12) additions, plus overhead | ≥ 126.5 |
| Cheon BSGS, double-and-add scalar multiplication per step (≈ 383 additions) | ≈ 125.2 |
| Cheon BSGS with KKM c = 3 tables (2 additions plus one F_r multiplication per step); same K and per-element budgets as above, computed by `evidence.py` | 117.862 |
| **this work** | **116.889745** |

- The gain over KKM-accelerated Cheon is a constant factor, about 0.97 bit (117.862 against 116.889745).
- **Relation to the generic bound.** Boneh–Boyen's generic lower bound Ω(√(ε·r/d)) has unstated
  constants. Here the walks use N_B + N_G ≈ 2^116.672 additions, and the generic coverage is
  N_B·N_G/m = K/(6(1−δ)) ≈ 0.0875. The factor 6 that lifts success to 1 − e^{−K} comes from the cheap
  action of U, through negation and the GLV endomorphism. The model does not capture that action, so
  there is no conflict with the bound.
- **Context, not a premise.** KKM (Pairing 2007, paraphrased) set aside a table of small-integer
  multiples as hard to use, because such integers are scattered among the exponents. The product
  construction a·b, with χ(a) + χ(b) ≡ 0, is one way to use them. We found no prior description of it,
  but nothing in the claim depends on novelty.

**Disclosed heuristics:**
- H1, with its structural loss bounded;
- H2;
- H3, the store's truncation rule with its discard fraction charged at δ = 0.01 and its per-access cost
  charged at the worst case of §4.2.

No other assumptions. No part of W is an average, a rate or a load-factor estimate: every probe, every
slot visit and every comparison is counted at the worst case the rule allows.
