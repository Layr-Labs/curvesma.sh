# Stored x-only leaf birthday search on secp256k1: six leaves per inverse, a Frobenius class key, metered rare-carry branches (125.252273153078)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; on every run its total work is at most
`time_log2 = 125.252273153078 = 125252273153078/10^12` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit). Peak memory is at most 2^135.9988 bytes; preprocessing, at most 2^110.1472 units, is included;
there is no advice. `certificate.py` checks these declared figures (`DECLARED`): the work, memory and preprocessing
figures are the least multiples of 10^-12, 10^-4 and 10^-4 at or above their exact values, the work per leaf 366.7729
the least multiple of 10^-4 at or above R/S, and 0.390000001 is below the success bound.

**Reviewer's guide.** The bound rests on seven things.
- **A short class key (2.2, 3.3).** With secp256k1's Frobenius `π = A + Bω` (`A² - AB + B² = p`), `B*β ≡ A - B` and
  `B*β² ≡ -A (mod p)`, so the class key `K = min(Bx mod p, (A-B)x mod p, -Ax mod p)` costs one product of x by a two-word constant and
  one fold, instead of a four-word product by β^j on two leaves in three. A leaf whose fold carries out of 2^256 is
  dropped (below 2^-91.4 per leaf); otherwise its key is exact outside a set of about 2^-221 per leaf. KEY performs no runtime
  window tests: a rare wrong index gives a leaf the key of another member of its own class, a charged missed collision, never
  a false match.
- **Six leaves per inverse (3.1).** `D = y_Y^2 - y_X^2 = x_Y^3 - x_X^3` is the product of the three denominators
  `x_Y - β^t x_X`, so one batched inverse of D gives the six leaves `x(φ^t(X) ± Y)`, at 1403 instructions per pair; the
  run checks every leaf of its batches of 1, 2, 5 and 9 incidences against exact affine arithmetic.
- **Rare carries by branch, metered (3.4).** SHX, FDX, F5, FS and FADDW branch over the rest of their carry or borrow
  chain unless its carry or borrow occurred. The skipped instructions would only add or subtract 0, so every result is
  bit-identical to the inline chain's; only the flags differ, and no path reads them before setting them (checked). A
  fixup costs 3 or 4 instructions and adds C0 to a stack word. It needs a row or result word in a set of about 2^-32 of
  all words (2^-30 for folds), so a node's mean fixup count is below 2^28.35. A node starts a batch only while its count
  is at most `2^45 - 46q2`: the ledger charges `4*2^45` per node, and Markov's inequality charges the cap's failure
  (below 2^-16.6 per leaf).
- **No uncharged transfer or bookkeeping (2.1).** A node's per-record steps are instructions of its own processor,
  charged at their executed counts: the emit inside `keys` (10 per record) at hop 1, `h2s` (10.25) at hop 2, and `bkt`,
  `pass4` and `scan_ranges` (17.546875), the unrolled loops' control included (section 4). Each message also pays 8 at each
  hop for its link: its delivery as `dlv` in `route.s` (6), which the run executes for every message of its 4 × 4
  machine, and its port's in-flight count, which DSB waits on (2). So the hops are charged 18 and 18.25. Links keep no
  flags, rings or counters. Every address, cursor, slot count and cap test is a counted instruction of the sender or of
  the node that later reads the data.
- **Success (section 6):** Chen–Stein on cross-node leaf pairs with b3 = 0 (6.2; its one non-routine step is the
  non-proportionality claim, for six leaf types), usefulness (6.3), the defect table (6.4), whose batch row comes from
  the arithmetic-site table of 3.4 and whose cap rows bound every pair of leaves, dependent ones included, and the final
  inequality (6.5). `certificate.py` evaluates the terms with exact Fractions and prints the bound 0.3900000039 (the
  run's `success_lower_bound`).
- **Work (section 7):** `certificate.py` recomputes the ledger exactly; the run checks each kernel call's count against
  its constants (`counts_match_certificate`).
- **Machine (2.1; scope in section 9):** per classical-group-work-v2, loads, stores, addressing, control and memory
  transfers count (the nodes', and each delivery's as `dlv`), instructions rather than cycles (a stalled store, load or
  WFE counts once), all processors' work with coordination (barriers, count exchanges), and memory has no scalar
  contribution (section 5).
- **Not relied on:** any heuristic, any prior acceptance or score, or extrapolation from the toy or the 4 × 4 run.

## 1. Overview

Each node computes the **leaves** `x(φ^t(X) ± Y)`, t = 0, 1, 2, for independent points X and Y, six per pair from one
batched inverse: Galbraith, Wang and Zhang (ePrint 2015/605) share one inverse between `X ± Y`, and φ(x, y) = (βx, y)
extends this to the three twists of X. Keyed by min(Bx mod p,(A-B)x mod p,-Ax mod p), except for the explicitly charged classifier and arithmetic exceptions, each leaf
goes in two bulk-synchronous hops to a home node, where equal keys give the logarithm. This is an analytic bound, not a
practical break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each with one scalar AArch64 processor and a private memory
below 2^64 bytes. The ledger counts the instructions the processors execute and, per message, its delivery and its
in-flight count (below); storage and wiring carry no price. The processor implements the base scalar AArch64 instruction
set with FEAT_RNG (`MRS Xt, RNDR`, section 7) and FEAT_LSE (`LDADD`, `STADD`, sections 3.4 and 4); every counted step is
one such instruction, whatever memory accesses it performs.
- **Ports and links.** A node has a column port and a row port, 48-byte blocks of device registers in its address space.
  To send, it stores a destination index (a row of its column or a column of its row, its own included) at +0, a
  destination address at +8 and four data words at +16..+47; the STP to +32 completes the message. The link delivers
  it: it writes the four words at that address of the destination's memory and reports the write back to its port
  (the in-flight count's decrement below). It keeps no flags, rings or counters; its only state is the message it
  carries. A store that finds the link busy stalls until the link takes it; a stalled instruction counts once.
- **Delivery is charged.** A delivery is charged as the work of `dlv` (`route.s`), 6 instruction-equivalents per
  message at each hop: three LDPs read the message's 48 bytes (index, address, data) from the link's block, two STPs
  write the four data words at the address, and RET. Performed by a processor, a delivery is this code run as a loop,
  with B in place of RET: its first LDP waits until the block holds a complete message, as a sender's store waits for
  a busy link, and a stalled instruction counts once. So 6 per message is all of a delivery's work, with no call or
  poll; the destination index selects the link's output, which is wiring. Part B of the experiment performs every
  delivery by executing `dlv` on the message, with the destination's memory as its address space, and counts 6 each.
- **Completion.** Each port counts its messages in flight: the count rises when the link takes a message and falls when
  the destination's write completes. These two updates are charged as 2 instruction-equivalents per message (sections 4
  and 7). DSB completes once the node's own stores are written and its ports' counts are zero. A node executes DSB
  before joining a barrier, so after a barrier every message sent before it is in memory.
- **Barriers.** A combining tree of epoch writes (the children of `(i, j)` are `(2i+1, j)`, `(2i+2, j)` and, for i = 0,
  `(0, 2j+1)` and `(0, 2j+2)`): a node waits for its children's arrival words, executes DSB and writes its parent's,
  then waits for its release word and writes its children's. A wait is LDAXR, CMP, B.EQ and WFE: a store into the word
  clears the node's exclusive monitor and wakes the WFE, each word has its own reservation granule, and nothing else
  wakes it.

### 2.2 Leaves and records

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145289174` and
`q2 = 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the
x-coordinates of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and
`S = 2^72*n0 = 2^126.699314` in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and
`n = λ^t*b_g`. A point Z ≠ O has the class `{±Z, ±φ(Z), ±φ²(Z)}`, its x-coordinates `{x, βx, β²x}` (mod p) have the
class key `K=min(Bx mod p,(A-B)x mod p,-Ax mod p)` (limbs K0..K3), where B is the Frobenius coordinate of3.3. Multiplication by beta permutes these three values, so the minimum is shared by the six-point class. Equality of complete exact keys implies equality of classes because all selected constants are B times powers of beta. A record is `[K2, K3, w2, w3]` with `w2 = i | k<<36`
and `w3 = d2 | c<<36 | l<<39`, where `c = 1 + 2t + [s = -1]` in 1..6 is the leaf type; its digits are `d1 = K0 mod 2^36`
(hop-1 row), `d2 = K1 mod 2^36` (hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3 mod 2^48`. At its home a record's source
is `(w2 mod 2^36, j)`, where j is the index of the slot it arrived in (2.3).

### 2.3 Phases

A pair slot holds `cap_p = 2^22 - 2^12` records of 32 bytes. A node in row i writes at hop 1 only into slot i of R1, at
`R1 + 32*cap_p*i`, at each node of its column, and a node in column j writes at hop 2 only into slot j of R2 at each
node of its row. Each port's cursor table holds, per destination, the sender's cursor in its slot there; the slot ends
at the same address at every destination, so the sender keeps that end in a register; `tinit` sets each of the 2^36
entries to the slot's base.

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators with their table
   entries, header constants, the two cursor tables, the d4 counters (`zero4`) and the d3 cursor table (`init3`). Then Q
   is broadcast in 2 hops.
1. **Hop 1:** the batch kernels (section 3) run for each X_k, and the emit advances the cursor of d1 by LDADD and writes
   each leaf through the column port to `(d1, j)` at the old cursor; a leaf whose old cursor is at or past the slot end
   is dropped. Then `fin` writes the final cursor of each row d, clamped at the slot end (drops advance cursors too),
   into entry i of d's table E1, and a barrier follows.
2. **Hop 2:** `h2s` reads slot i of R1 up to `E1[i]`, for every i, and writes each record through the row port to
   `(own row, d2)` by the same cursor rule; `fin` writes each home's final cursor into its `E2[j]`; a barrier follows.
3. **Grouping:** `bkt` reads each slot of R2 up to its E2 entry and appends the record's address and d4 to its d3 bucket
   of `cap_3` entries, which ends at a multiple of 2^38, or drops the record if the bucket is full; a stored record
   counts its d4. `prefix4` forms the d4 cursors, `pass4` writes each entry's address into the pointer array at its d4
   cursor, bucket by bucket in d3 order, which leaves each d4 counter at its range's end, and `scan_ranges` lists,
   within each d4 range, adjacent records with equal K2 as candidates, up to 2^28; a full list ends it.
4. **Resolution,** in four sub-phases: a home `(h_r, h_c)` sends two requests per candidate to the sources `(g_r, g_c)`
   via the relay `(g_r, h_c)`, and replies of at most 8 records (one leaf's scalars and points) return via `(h_r, g_c)`.
   Each sub-phase is a count exchange: senders write each receiver their count and the epoch; after a barrier each
   receiver returns each sender's base, in sender order, or a drop mark if its records would pass the cap (2^29 requests
   at a source, 2^30 records at a relay); after a barrier senders write the records at their bases; a barrier ends it.
   The home recomputes `L = [m]P + [n]Q` and `L'` from the records' leaf types, finds the unit u with `L' = uL` if any,
   and, if `n' - u*n ≠ 0`, outputs `x = (u*m - m')/(n' - u*n) mod r` after checking `[x]P = Q`.

**Barriers:** 21 (program delivery 2, Q broadcast 2, the ends of phases 0 to 3 4, phase 4 12, output 1), within the
allowance of 32. **Caps:** arithmetic fixups (3.4), j = 2 key paths (3.2), pair slots, d3 buckets, candidates, requests,
relay records and coin retries are capped or bounded in total, so the ledger holds on every run.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's. With `D_l = y(Y_l)^2 - y(X)^2`, the forward pass stores `D_l` and `acc_(l-1)` in entry
l and forms `acc_l = acc_(l-1)*D_l` (`acc_(-1) = 1`); fe_inv gives `R = 1/acc_(q2-1)`; the backward pass, for
`l = q2-1..0`, sets `inv = R*acc_(l-1)`, `R <- R*D_l`, `u = inv*(y_Y - y_X)` and `v = inv*(y_Y + y_X)`, and for t = 0,
1, 2 writes

    x(φ^t(X) + Y) = Q_t²u² - β^t x_X - x_Y,    x(φ^t(X) - Y) = Q_t²v² - β^t x_X - x_Y,
    Q_t = x_Y² + β^t x_X x_Y + β^(2t) x_X².

**Why.** As R = 1/acc_l before step l, inv = 1/D_l. The curve gives `y² = x³ + 7`, so `D = x_Y³ - x_X³`, and as β³ = 1,
`(x_Y - β^t x_X)*Q_t = x_Y³ - x_X³ = D`. The slope of `φ^t(X) ± Y` with `φ^t(X) = (β^t x_X, y_X)` is
`(±y_Y - y_X)/(x_Y - β^t x_X) = (±y_Y - y_X)*Q_t/D`, which is `Q_t*u` or `-Q_t*v`; its square, less the two
x-coordinates, is the leaf. D ≠ 0 unless `x_Y³ = x_X³`, that is `Y = ±φ^t(X)` (3.4).

**Arithmetic per pair**, after W2S (`inv`, the new R) and W2L (u and v, on the shared rows of inv):
- u² and v² (SQ5 and F5: a square and its fold);
- `Q_t = (x_Y + β^t x_X/2)² + (3/4)β^(2t) x_X²` for t = 0, 1 (FADDW, SQ5, AC4 adds the constant, F5), and
  `Q_2 = 3x_Y² - Q_0 - Q_1` (two FS), as `1 + β + β² = 0`;
- for each t: Q_t² (SQ5, F5), the two products Q_t²*u² and Q_t²*v² on the shared rows of Q_t² (the pair schedule of
  W2L), and from each product's six-word sum `Z` the leaf `Z - s'' + 2^256`, folded (SB6, FDX), where
  `s'' = β^t x_X + (x_Y + C0)` (FADDW), so that the fold's 2^256 ≡ C0 cancels.

The header holds `x_X, y_X, βx_X, β²x_X, x_X/2, βx_X/2, 3x_X²/4, 3β²x_X²/4`; a Y entry holds `y², y, x, 3x², x + C0`;
the six leaves overwrite the 256-byte entry in the order `c - 1 = 2t + [s = -1]`. The arithmetic is F's row schedule;
the rare carries and borrows of SHX, FDX, F5, FS and FADDW are taken by branches over fixups that give the inline
chains' results (3.4), and its only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks entries downward and processes their six leaves. KEY_APPROX uses resident B and D=A-B to approximate their products' top words; the exact branch rule is proved in KEY_LEMMA.md. Branch-specialized KPC multiplies by B or D without per-leaf constant loads. The j=1 continuation is out of line. The j=2 continuation uses KNEG; its negation constants -C0 and -1 are initialized once in x11 and x15. A fold carry drops that leaf before emission and before incrementing its j=2 counter.

The unchanged emit costs10 instructions: compute d1 and cursor address, LDADD by32, CMP/B.HS against slot end (5); then form w3 and write three STPs (5). Pointer/tag layout remains bit-for-bit the original36/3/25-bit metadata: entries256bytes, base2^34 aligned, pointer shift31. No phased routing or35-bit layout is part of this candidate.

Exact `keys` instruction count for q2 incidences is
`410*q2+11+12*j2-j0-5*slotdrops-10*F0-11*F1-4*F2`.
Here j0 counts every classifier-index0 leaf, including fold drops; j2 counts only non-fold-dropped index2 leaves; Fi counts fold drops on branch i. Normal leaf costs are67,68,80 for branches0,1,2; folded-away leaves cost57,57,64. Loop SUB/TBZ adds2 per incidence and fixed setup/return adds11. The ledger conservatively ignores every saving term and charges410*q2+11 plus12 per actually counted j2 leaf. The protected emulator tests the exact formula, including planted fold drops and full slots.

The node admits a batch only when its j2 count is at most cap_j-6q2, where cap_j=ceil(n0/3)+2^42. Use the unconditional ideal-point domination J and variance proof in KEY_LEMMA.md, with the retained larger d_wj. Computational defects are union-charged separately; do not condition the variance argument on cleanliness. Fold drops only reduce the count.

All other leaf arithmetic is unchanged: fwd99*q2+5, bwd1403*q2+5 plus their metered3/4-instruction fixups; fe_inv23646. Kernels may overwrite x0..x29 and flags, preserve sp, and preserve x30 as specified by the original routines. The driver keeps no registers across calls. Entry-pointer termination remains TBZ bit33 with q2=2^25; no enlarged batch is used here.

### 3.3 The class key (macros KEY, KCOM and KNEG)

**Constants.** With the trace `τ = p + 1 - r`, `4p = τ² + 3B²` for the positive `B = 0xe4437ed6010e88286f547fa90abfe4c3`
(below 2^128), and `A = (τ + B)/2 = 2^128 + a'` (a' below 2^125) gives `A² - AB + B² = p`: `π = A + Bω` is secp256k1's
Frobenius in Z[ω]. Hence `((A - B)/B)² + (A - B)/B + 1 = p/B² ≡ 0 (mod p)`, so `(A - B)/B` is a cube root of 1 mod p
other than 1; it is the package's β, and `B*β ≡ A - B`, `B*β² = (A - B)²/B ≡ -A (mod p)`. So `c_j = B*β^j mod p` is
B (128 bits), A - B (126 bits) or `p - A` for j = 0, 1, 2 (the run checks `frobenius_norm_p_and_trace` and
`key_constants_are_B_beta_j`).

**The key and classifier.** The exact desired index minimizes c_j*x modulo p, not beta^j*x. The desired class key is their minimum; the classifier selects it outside the charged window set. KEY_APPROX uses five instructions per128-bit constant; comparisons and branches select the product continuation. KEY_LEMMA.md gives the complete global weak-input approximation, classifier correctness outside explicit windows, fixed bad-set count at most352*2^192, j2-meter domination and unconditional variance proof. This count is smaller than the retained d_wj allowance; no probability improvement is claimed from the spare margin.

Wrong indices still return a selected B*beta^j multiple of the same x, so complete-key equality cannot join distinct six-point classes. Witness misses are charged by d_wj. For capacity and intruder analysis, the unconditional ideal-key coupling in KEY_LEMMA.md is essential: the global LT_BAD/KEY_BAD charge excludes noncanonical leaves/outputs, and wrong-index competitors are included in the fixed ideal map. No unproved representative-independence assumption is used.

**KPC/KCOM (40) and KNEG (50).** KCOM multiplies x by c_j's two words in two rows of 12 and adds them (5): `z = c_j*x`,
below 2^384. Its fold (10) forms `t = (z mod 2^256) + floor(z/2^256)*C0`, below `2^256 + 2^161`, and the carry c of t
out of 2^256, which its B.HS tests. KNEG forms `z = a'*x + 2^128*x = A*x`, below 2^385, from the two rows, their sum and
four adds of x at limb 2 (29 + 5), folds words 4..6 the same way (11; the seventh word is at most 1), tests c by a B.HS,
and returns `p - (t mod 2^256)` using the resident low word -C0 and all-ones word: SUBS, two SBCS and SBC cost4. The two constant-forming instructions are removed from every j2 leaf; the extra setup instruction is included in KEY_FIX=11. **Drops.** If c = 1 the leaf is dropped: KCOM's B.HS goes to the end of the emit and KNEG's past the j = 2 count
and the emit, so no record is written, no cursor moves and the count does not rise. Then `t mod 2^256 < 2^162`, and the
path value, `c_j*x mod p` (KCOM) or `A*x mod p` (KNEG), is `(t mod 2^256) + C0 < 2^162 + C0`: a path drops only leaves
whose key K is below `2^162 + C0` (KCOM) or above `p - 2^162 - C0` (KNEG), at most `3(2^162 + C0)` leaf values over the
three paths, as `x -> c_j*x mod p` is a bijection; 6.4 charges them as the per-leaf row `d_kd`, below 2^-91.4.
**Exactness.** If c = 0, t is congruent to the exact value and below 2^256, so the result is exact unless `t >= p`,
which needs the exact value of `c_j*x` (KCOM) or of `A*x` (KNEG) mod p to be below C0: then KCOM returns K + p with K
below C0, and KNEG returns K + C0 with K above `p - C0` (t is never 0, as x ≢ 0: no point has x = 0). So a path returns
`c_j*x mod p` exactly, or drops the leaf, unless that value lies in **BAD**, `{K < C0}` for KCOM and `{K > p - C0}` for
KNEG, C0 values per path; at most 3·C0 leaf values meet BAD on any path (the key row of 3.4). The run builds path values
below C0 (K above `p - C0` for KNEG), in `[C0, 2^128)`, where the fold carries but for a fraction of about 2^-30, in
`[2^128, 2^162)`, where it may, just past that range and at random, and checks that each path drops the leaf exactly
when its fold carries, deviates exactly when its bit-level predicate holds (c = 0 with `t >= p`), and that every drop
and deviation lies in its set (`site_checks.key_paths`).

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts `B_(k+1)
= SHX(B_k)`; RF and RA sum `Z = Σ a_k*B_k < 2^322` exactly (each high word is at most 2^64 - 2), and a fold (FDX)
computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ5 forms the square of a weak a as `V = L +
H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), AC4 may add a constant below p (top limb at most
C0 + 2), and F5 folds V the same way. SB6 forms `Z' = Z - s'' + 2^256` exactly (top word at most 4). Only these steps drop
a carry or a borrow: SHX's last ADCS, the final ADC of FDX and of F5, and FS when its subtrahend is a fold's output at
or above p, besides the key paths of 3.3. Without a dropped carry each returns a value below 2^256 congruent to the
exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; the fused fold of
  `Z' < 2^322 + 2^256` likewise needs the leaf in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`, or Q_t in
  `[C0, C0² + 3C0)` with the constant added.
- **FS** `a - b` with `a < p` and `b < 2^256`: wrong iff `b - a > p`, which needs `b >= p`; the subtrahends Q_0 and Q_1
  are fold outputs, at or above p only when the true value is below C0.

A fold returns a value `>= p` only when the true value is below C0. That matters only where a canonical value is needed:
the leaf entering KEY, the emitted key K, and FS's subtrahend (above). Every other operation is exact on its inputs: RF
and RA; SB6; the key paths' rows and sums (exact integers below 2^385; only their folded output can deviate, 3.3); FS
and FADDW on canonical inputs; FADDW of a canonical value and `x_Y + C0` (the sum is
below `p + 2^256`, so one conditional subtraction of p leaves it below 2^256); KEY's integer arithmetic (3.3 covers its
index); the emit; and F's `fe_mul_f` and `fe_sqr_f`, which fold the 512-bit product, fold its fifth limb again and add
C0 exactly when that last sum carries past 2^256, so fe_inv (255 squarings, 15 multiplications) is exact for every input
below 2^256. Generator construction, the table entries, the header values, the resolver and `[x]P = Q` run on
`fe_mul_f`, `fe_sqr_f` and fe_inv, each coordinate brought below p (an ADDS chain of C0 detects `t >= p`, CSELs keep
`t - p`) before it is subtracted, compared or stored, with exact arithmetic mod r for the solve. So they have no sites;
the tables hold exact canonical values, and `x_Y + C0 < 2^256` as `x_Y < p`.

**Distribution.** Given X ≠ O and every other generator, `Y_l` is uniform on the r - 1 points of G \ {O}. A value of
x_Y comes from at most 2 points, and of y_Y from at most 3. So `D_l = x_Y³ - x_X³`, `acc_l = acc_(l-1)*D_l`,
`R_l = 1/acc_l` and `inv_l = 1/D_l` take each value with probability at most `6/(r-1)` (the rest are bijections while
no `D_m = 0` for m ≤ l, the degenerate case below); `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)` (as
`D = (y_Y - y_X)(y_Y + y_X)`) at most `3/(r-1)`; u² and v² at most `6/(r-1)`; Q_t, quadratic in x_Y, at most `4/(r-1)`;
Q_t² at most `8/(r-1)`; and each leaf at most `2/(r-1)`, as `φ^t(X) ± Y_l` is uniform on r - 1 points. A row
`B_k ≡ v*2^(64k)` determines its operand v, so at most 2^97 values of v make a given row shift carry.

| site, at incidence l | runs per batch | fresh generator: operand | bad values | points per value |
| --- | --- | --- | --- | ---: |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | `Y_(l-1)`: `acc_(l-1)` | 2^97 | 6 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | `Y_l`: `acc_l` | 2^66*C0 | 6 |
| bwd: row shifts of `R_l` | 3q2 | `Y_l`: `1/acc_l` | 2^97 | 6 |
| bwd: fold of `inv_l = R_l*acc_(l-1)`, l ≥ 1 | q2 - 1 | `Y_l`: `1/D_l` | 2^66*C0 | 6 |
| bwd: fold of `R_(l-1) = R_l*D_l`, l ≥ 1 | q2 - 1 | `Y_(l-1)`: `1/acc_(l-1)` | 2^66*C0 | 6 |
| bwd: row shifts of `inv_l` | 3q2 | `Y_l`: `1/D_l` | 2^97 | 6 |
| bwd: folds of u and v | 2q2 | `Y_l`: u, v | 2^66*C0 | 3 |
| bwd: square folds of u² and v² | 2q2 | `Y_l`: u², v² | C0² + 2C0 | 6 |
| bwd: square folds of Q_0 and Q_1 | 2q2 | `Y_l`: Q_t | C0² + 3C0 | 4 |
| bwd: FS forming Q_2 | 2q2 | `Y_l`: Q_0, Q_1 | C0 | 4 |
| bwd: square folds of Q_t², t = 0, 1, 2 | 3q2 | `Y_l`: Q_t² | C0² + 2C0 | 8 |
| bwd: row shifts of Q_t² | 9q2 | `Y_l`: Q_t² | 2^97 | 8 |
| bwd: fused folds of the six leaves | 6q2 | `Y_l`: the leaves | (2^66 + 2)*C0 | 2 |
| keys: key paths of the six leaves (KCOM or KNEG; a drop is not a deviation) | 6q2 | `Y_l`: the leaves | 3*C0 | 2 |

A run fails with probability at most its bad count times its points per value over r - 1. At l = 0 three row shifts
and three folds cannot fail, and are not in the table:
- fwd's rows of `acc_(-1) = 1` are 1, 2^64 and 2^128, and its fold of `D_0*1` has `Z < 2^256`;
- bwd's fold of `inv_0 = R_0*1` has `Z = R_0 < 2^256`;
- bwd's last update `R_0*D_0 ≡ 1` folds to 1 or p + 1, as nothing else congruent to 1 lies below `2^256 + 2^66*C0`. It
  is left in x25..x28 after the last pass and never read.

The degenerate cases are X = O (probability 1/r) and `D_l = 0`, that is `Y_l = ±φ^t(X)` for some t (at most `6/(r-1)`
per incidence), which zeroes acc and spoils the batch; when `D_l ≠ 0`, `Q_t = Π_(t'≠t)(x_Y - β^(t') x_X) ≠ 0`. The
first site of a batch whose predicate holds has exact operands, so a batch has a defect with probability at most
`d_batch`: the table's sum plus `6q2/(r-1) + 1/r` (`SITES` in `certificate.py`; 6.4).

**Rare carries by branch.** In SHX the carry out of its first ADDS, in FDX and F5 the carry into limb 2, in FS the
borrow of its C0 subtraction and in FADDW the carry of its C0 addition are rare. Each macro branches over its **fixup**
(B.LO; B.HS for FS) unless that carry or borrow occurred. The fixup is the rest of the inline chain, two instructions
(ADCS and ADCS or ADC) or three for FS and FADDW, and `STADD x4, [sp]`. Where the branch is taken the skipped
instructions would add or subtract 0, so on every input every register equals the inline chain's, and so does every
deviation of the table above; only the flags differ, and every path from the join sets them before reading them or
returns (`fixups.flags_dead_at_every_join`; kernels may leave the flags changed). A fixup adds 3 instructions (SHX, FDX,
F5) or 4 (FS, FADDW) and C0 to the word at [sp], which opens a 16-byte frame pushed as zero by `fwd` (STR) and by
`bwd` beside its saved x30 (STP); each returns the word in x0, at most `46q2*C0 < 2^63` per call.
The driver adds x0/C0 to the node's fixup count N and starts a batch only while `N <= 2^45 - 46q2` (46 sites per
incidence: 18 SHX, 11 FDX, 7 F5, 4 FS and 6 FADDW). So N is at most 2^45 on every run, and section 7 charges `4*2^45`
per node.

**How often.** A fixup needs a word in a small set:
- SHX's carry needs the row B it shifts to have its low word above `2^64 - C0`, as the word it adds is below C0: at most
  `C0*2^192` rows.
- FDX's carry into limb 2 leaves its output T with `T mod 2^128 < F*C0 <= 2^66*C0` and `T >= 2^128`, where `F =
  floor(Z/2^256) <= 2^66` (for the fused fold too, as `Z' < 2^322 + 2^256`): at most `2^194*C0` words.
- F5's carry into limb 2 leaves `T mod 2^128 < 2^65`, as its fifth word is at most C0 + 2, so that word times C0 has a
  high word of at most 1: at most 2^193 words.
- FS's borrow leaves an output whose low word is at least `2^64 - C0`, and FADDW's carry one whose low word is below C0:
  at most `C0*2^192` words each.

While no site of the batch up to and including this one has deviated, the row or output is a word below 2^256 congruent
to the site's operand value (times `2^(64k)` for the k-th row), which therefore determines that value. So, given X and
every other generator, a site takes its fixup with probability at most its set's size times the points of `Y_l` per
operand value over r - 1, with the fresh generator and multiplicities of the table above; `y_Y ± y_X` have 3 points per
value, and `x_Y + β^t x_X/2` and `s''` have 2. At l = 0 the six sites the table leaves out cannot take a fixup either:
rows 1, 2^64 and 2^128 have a low word of at most 1, folds of values below 2^256 add nothing, and `R_0*D_0` folds to 1
or p + 1, neither of which is at least 2^128 with its low 128 bits below `2^66*C0`. Summed over a batch (`FIXUPS` and
`E_FIX` in `certificate.py`), the expected number of fixups a batch takes before its first deviation, and so of a batch
without a defect, is at most `E_FIX`, and over a node's batches at most `q1*E_FIX < 2^28.35`.

| fixup sites, at incidence l | runs per batch | words | points per value |
| --- | --- | --- | ---: |
| fwd: FS forming `D_l` | q2 | `C0*2^192` | 6 |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | `C0*2^192` | 6 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | `2^194*C0` | 6 |
| bwd: row shifts of `R_l` | 3q2 | `C0*2^192` | 6 |
| bwd: folds of `inv_l` and `R_(l-1)`, l ≥ 1 | 2(q2 - 1) | `2^194*C0` | 6 |
| bwd: FS of `y_Y - y_X`, and FADDW of `y_Y + y_X` | q2 each | `C0*2^192` | 3 |
| bwd: row shifts of `inv_l` | 3q2 | `C0*2^192` | 6 |
| bwd: folds of u and v | 2q2 | `2^194*C0` | 3 |
| bwd: square folds of u² and v² | 2q2 | 2^193 | 6 |
| bwd: FADDW of `x_Y + β^t x_X/2`, t = 0, 1 | 2q2 | `C0*2^192` | 2 |
| bwd: square folds of Q_0 and Q_1 | 2q2 | 2^193 | 4 |
| bwd: two FS forming Q_2 | 2q2 | `C0*2^192` | 4 |
| bwd: square folds of Q_t² | 3q2 | 2^193 | 8 |
| bwd: row shifts of Q_t² | 9q2 | `C0*2^192` | 8 |
| bwd: FADDW of `s''` | 3q2 | `C0*2^192` | 2 |
| bwd: fused folds of the six leaves | 6q2 | `2^194*C0` | 2 |

**The meter row.** A node stops early only if N exceeds `2^45 - 46q2` at the start of some batch. On the event that none
of its batches has a defect, N is that clean count, so by Markov's inequality a node stops with probability at most
`q1*E_FIX/(2^45 - 46q2) + q1*d_batch < 2^-16.6` (`d_f`, a node-local row of 6.4).

## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (10; 5 if the slot is full), and `dlv` (6) | 16 | 18 |
| hop 2: `h2s` (10; 9 if the slot is full; loop 0.25), and `dlv` (6) | 16.25 | 18.25 |
| bucket pass `bkt` (10; 8 if the bucket is full; loop 0.25) | 10.25 | 10.25 |
| local pass `pass4` (4; loop 0.25) | 4.25 | 4.25 |
| `scan_ranges`: 3.046875 per record (paired CMP/CCMP tests), plus 8 per listed candidate and at most 26 per d4 range | 3.046875 | 3.046875 |
| **total** | **49.796875** | **53.796875** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlv`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. `h2s`, `bkt` and `pass4` take a slot's or bucket's records in blocks
of 1, 2 and 4, as the low bits of its count say, then 8 at a time; the 8-record loop's CMP and B.LO, 2 per pass, are
charged as 0.25 per record (a slot of n records runs floor(n/8) passes), and the tests and set-up, 11, 12 and 11 per
slot or bucket, are per-node rows. **Cursors.** The emit and `h2s` advance a destination's cursor with one LDADD, which
returns the old cursor: the record goes there unless that cursor is at or past the sender's slot end, which is the same
at every destination and held in a register. On a drop `h2s` stores the end back, so a hop-2 cursor stays at most its
slot end; a hop-1 cursor can pass its end, but only by the node's own leaves, 32 bytes each, so it stays below 2^63.2
(checked), and `fin` clamps each final cursor at the slot end (CMP, CSEL). `h2s` and `keys` add 64 to their port
argument and address the register block at -64, -48 and -32. **Buckets.** Bucket t's `cap_3` entries end at `R1 +
(t+1)*2^38`, so the cursor LDADD returns to `bkt` has its low 38 bits zero exactly when the bucket is full: TST and B.NE
send a bucket with room to the store, and on a full bucket a STADD undoes the advance, so after each record a cursor is
at most its bucket's end. `bkt` steps its record pointer with a pre-indexed LDP and counts d4 with one STADD; `pass4`
moves its d4 cursors with register-offset LDR and STR. **Scan** (`scan_ranges`, from submissions 80e9e8af and 1efcef88).
`prefix4` sets each d4 counter to its range's start and `pass4` advances it per pointer it writes, so after pass4 the
counter is the range's exclusive end. `scan_ranges` visits the 2^48 ranges in order and never compares across one. In a
range of n >= 1 records it loads the first record, handles one odd comparison when needed, then compares adjacent pairs using CMP A,B; CCMP B,C,#4,NE; B.EQ. This branches exactly when A=B or B=C. The slow path rechecks both equalities, appends in original order and checks capacity separately before either store. Therefore the candidate list and abort position are unchanged.

Three key registers rotate and four pointer-register pairs rotate. Initial copies put the range's first key and address in every possible suffix-start register. The pair remainder r=((n-1)//2) mod32 selects entry at pair32-r (r=0 selects the end test), using five TBZs and one additional branch when r is odd. Each pair uses the previous key and record pointer, loads two new pointers and their keys, and retains the second key/address for the next pair. All potential entry registers initially hold the same first key/address. Since32 mod3=2 and32 mod4=0, the final key is x7 and final pointer x16, so MOV x5,x7 normalizes the next full block. The hot block costs32*(LDP+2LDR+CMP+CCMP+B.EQ)=192 instructions plus MOV,CMP,B.LO:195 per64 records.

Write n-1=64q+2r+epsilon with0<=r<=31 and epsilon in{0,1}. A nonmatching nonempty range costs25+6epsilon+(r mod2)+6r+195q including range control. Subtracting195n/64 gives1405/64+189epsilon/64-3r/32+(r mod2). Its maximum is1652/64=25.8125 at epsilon=1,r=1, bounded by26. Empty ranges cost5. The slow path adds7 instructions for only the first equality,8 for only the second,11 for both; an odd-prefix hit adds4. Charge8 per listed candidate plus one additional8 for a full-list append attempt. Every append has its own capacity check and preserves original adjacent-pair ordering, even when both equalities hold. Global guard and RET cost2.

Thus a node uses at most195n/64+26*2^48+8*(2^28+1)+2 scanner instructions. This rational inequality holds for each node separately, so summing yields195S/64 without an uncharged ceiling. Production S is divisible by64. The remainder suffix consumes exactly2r records after the optional odd comparison; every full block consumes64, and the exclusive end test prevents out-of-range loads. Empty ranges read no record, and no comparison crosses a range boundary. Candidate identities, order, capacity-abort behavior, and the original partial-match success/cost analysis are unchanged.

A hop's receiving
node executes nothing on arrival; it reads the records in the next phase, with the LDPs of `h2s` (hop 1) or `bkt` (hop
2). Per-node costs are rows of section 7, and every address a kernel loads or stores lies in the regions of section 5.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 key bits. The records with one key form one run unless an **intruder**, a record
with another key and the same 144 bits, sits inside (6.4); a run lies in one d4 range with one K2, so its adjacent pairs
are all listed unless the candidate list is full. `scan_ranges` lists every adjacent pair of a range with equal K2,
which agrees on d1, d2, K2 and d4, 184 key bits: besides a run's pairs it can list partial matches of different classes,
which the resolver rejects (no unit relates their points) and 6.4 charges to the candidate, request and relay caps
through KB(184). Writing the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`, an
equivalence relation, so a run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums the regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.0080 used) | code, stack, the five barrier words (one reservation granule each), header, X table (eight header values and y² per X, 288 bytes), Y table (160 bytes per Y), scalars, batch entries (2^33 bytes at a multiple of 2^34), the two ports with cursor tables of 2^39 bytes each, E1 and E2 (2^41 each), the d3 cursor table, the d4 counters (2^51) and the phase-4 buffers and tables (`FIXED` in `certificate.py`) | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | hop-1 pair slots | phases 1-2 |
| the same range, from phase 3 | d3 buckets in 2^24 strides of 2^38 bytes (`16*cap_3 = 2^37.7805` bytes at the end of each) and the pointer array `8*2^24*cap_3`, together 2^62.5155 | phase 3 |
| R2 after R1, ending at 2^63.998943 | hop-2 pair slots (phase 4 reads candidates' records at their addresses) | phases 2-4 |

So every node stays below 2^64 bytes, and the peak, in phase 2, is `2^72*(fixed + |R1| + |R2|)`, at most 2^135.9988
bytes (log2 ≈ 135.9987677). A remote item is reached only by a message.

## 6. Success (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494345441` (the run's `mu`).
Same-node pairs are left out; their collisions only add to success.

### 6.2 Chen–Stein on cross-node pairs

Let `W = Σ_(α∈Γ) I_α`, and let `B_α` be the pairs that share a generator point with α. I_α depends only on α's four
generators, distinct because they come from two nodes, so I_α is independent of `{I_β : β ∉ B_α}` (b3 = 0). A generator
lies in 6q2 (an X) or 6q1 (a Y) leaves, so `|B_α| <= 12(q1 + q2)*S`.

For a unit u, the event `L' = u*L` with `L = φ^t(X) + s*Y` and `L' = φ^(t')(X') + s'*Y'` is the linear relation
`λ^(t')x' + s'y' - uλ^t x - us*y = 0` in the generators' logarithms; I_α is at most the sum of these 6 events. Take
α ≠ β in Γ sharing a generator, and units u and v. **The two relations are not proportional:** proportional relations
have the same support, so each node contributes exactly one X and one Y to both, and β's leaf in each node uses both:
β is `{φ^(t2)(X) + s2*Y, φ^(t2')(X') + s2'*Y'}` and its event is `L2' = v*L2`. A factor c with `λ^(t2') = c*λ^(t')` and
`s2' = c*s'` is both ±1 and a power of λ, so c = 1 (-1 is not: every power of λ cubes to 1, and (-1)³ = -1),
`t2' = t'` and `s2' = s'`; then `vλ^(t2) = uλ^t` and `v*s2 = u*s` force `λ^(t2 - t) = s2/s`, so `t2 = t`, `s2 = s` and
β = α. So some 2 × 2 minor is invertible: given all
other generators, the relations fix two generators, each uniform on at least r - 1 points, and
`E[I_α I_β] <= 36/(r-1)²`. The Chen–Stein method (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989), Theorem 1)
gives `|P(W = 0) - e^-μ| <= b1 + b2 <= |Γ|*12(q1 + q2)S*(P1² + 36/(r-1)²)`, which the certificate checks is below
2^-95.

### 6.3 Usefulness

A cross-node colliding pair is **detected** if both records carry their class's key (no batch defect, no wrong class
index) and are emitted (no key-fold drop), both survive the caps, no intruder enters their run, and the candidate,
request and relay caps are not reached first. Its run is then resolved, a useful adjacent pair (section 4) gives x, and
`[x]P = Q` is checked before output. The pair is useful iff `n' ≠ u*n`, where `n = λ^t*b_g` and `n' = λ^(t')*b_(g')` for
its two nodes; given all points, this fails with probability exactly 1/r. So `P(fail) <= P(W = 0) + E #{α : I_α = 1, α
useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (a batch defect of 3.4, a wrong class index, a key-fold drop, either meter) is a function
of node g's coins, and L' is uniform given them, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs at most
`μ(1 + 1/r)` times its per-leaf rate.

**Caps and intruders.** Condition on α's four generators. A leaf with a generator outside them, a **free** generator,
is then uniform on at least r - 1 points, so at most `6*2^(256-f)/(r-1)` (`KB(f)`) of its mass falls on given values of
f key bits. The ideal point-only key map in KEY_LEMMA.md has at most six point preimages per value, including wrong-index points. All ideal variables below are defined under the original generator distribution, never conditioned on global cleanliness; the explicitly stated witness-generator conditioning is retained; the global exceptional-set charge couples retained actual records to these ideal records. The
12 leaves of α's two incidences are fixed; the six leaves of an incidence through one of α's generators and a free one
are functions of that free generator.
- **Pair slots and d3 buckets** count single leaves. A pair slot has mean at most `m_p = n0*KB(36) + 1` from the free
  leaves (at hop 1 it holds one node's leaves with a given d1, at hop 2 one column's leaves with a given (d1, d2), the
  same bound), and a d3 bucket `m_3 = S*KB(96) + 1`; each cap leaves room for the at most 12 fixed leaves. Two leaves'
  indicators have covariance 0 without a common free generator or with a common free Y; at most `1/(4(r-1)²)` with a
  common free X (at most `KB(f)/(r-1)` when one of them also uses α's Y); at most `1/(r²(r-1))` times the two cells'
  sizes for two leaves of one incidence with both generators free, whose pair is uniform on r(r - 1) pairs (an
  invertible linear image of (X, Y), with Y ≠ O); and at most one leaf's mean for two leaves of one incidence through
  α's generators, at most 60(q1 + q2) ordered pairs. So a count's variance is at most its mean plus
  `1 + 60(q1 + q2)*KB(f)`, and Chebyshev bounds the drop of a given record by `(m + 1 + 60(q1 + q2)KB(f))/(c - 12 - m)²`
  for a cap c and mean bound m (`d_p`, `d_3`).
- **Intruders and the candidate, request and relay loads** count pairs of leaves, jointly with I_α. For a pair
  π = {A, B} ≠ α with event E_π (keys equal on f bits; at the candidate list and the relay, also A's home or home column
  is α's), `E[I_α*1(E_π)] <= ρ*(6/r)*P*(π)`, where `P*(π)` is E_π's probability for independent uniform leaves and
  `ρ = (r/(r-1))^4` covers the exclusion Y ≠ O. (i) If one of them, say A, has a free generator that the other lacks,
  condition on everything else: A is uniform on at least r - 1 points while I_α, B and α's home are fixed, so the f-bit
  match costs at most `(r/(r-1))KB(f)` and the home condition falls on B. If B has a free generator, B is uniform given
  α's generators, and its home is α's with probability at most `(r/(r-1))KB(72)` (KB(36) for a home column); if B is a
  partner of L (or L'), (L, B) (or (L', B)) is uniform up to ρ over that node's coins, given which `P(I_α) <= 6/r`; if
  B is L or L', π is a pair through L (below). (ii) Otherwise A and B lie in
  one incidence or in α's incidences. If their incidence is free of α's generators, condition on α's: (A, B) is an
  invertible linear image of the incidence's two generators, so uniform up to ρ. If both lie in node g of α and their
  incidence uses α's X or Y there (α's own incidence included), E_π (with α's home, which on I_α is L's) is a function
  of g's coins, given which L' is uniform, so `P(I_α | g's coins) <= 6/r`; over g's coins (A, B) is uniform up to ρ, and
  L is uniform given (A, B) through a generator they lack, unless A is L (the pairs through L, and likewise L',
  are counted without the home condition, as both are at α's home) or A and B are both partners of L (20 pairs, 10 per
  incidence of α, also counted without it). If A is in L's incidence and B in L''s, `(L, A)` and `(L', B)` are uniform
  up to ρ, so each of the six events `L' = uL` meets E_π with probability at most `ρ*P*(π)/r`.

  Hence the intruder rate is at most `ρ(S + 12)KB(144)` per colliding pair; the candidate list at α's home has mean
  `other_h <= ρ(S²*KB(72)*KB(184) + (2S + 20)*KB(184))` (pairs at the home equal on d1, d2, K2 and d4, the 184 key
  bits `scan_ranges` matches; L and L' are there); a source receives at most `other_g <= 2ρ*n0*S*KB(184)` requests; a
  relay holds at most `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(184) + (2S + 20)*KB(184))` records (a request is one record and a reply at
  most 8).
  Markov with the caps less α's own entries (4 in the candidate list, 8 requests, 16 records) gives `d_cl`, `d_rq` and
  `d_rl`. Leaves of a defective batch elsewhere are not uniform; the whole-run row "a computational defect in any batch"
  excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-126 |
| wrong class index (3.3): the new fixed bad-set union, bounded conservatively by the retained old allowance | `d_wj` | < 2^-54.29 |
| j = 2 meter, Chebyshev with slack `2^42 - 6q2 - 1 - n0*d_wj` | `d_j` | < 2^-31 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^45 - 46q2` (3.4) | `d_f` | < 2^-16.6 |
| key fold carries: the leaf is dropped (3.3) | `d_kd` | < 2^-91.4 |
| pair slot at hop 1, and again at hop 2 (Chebyshev, `cap_p`) | `d_p`, counted twice | < 2^-19.9 each |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `ρ(S + 12)KB(144)` | `lam_int` | < 2^-14.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-25.4; < 2^-28; < 2^-23.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.9 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.9 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.390000001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 145289174, where the left side is at least 0.3900000039, and
the target is 10^-9 above the rules' 0.39.


## Remaining proof sections

Sections7-9, the exact ledger, experimental scope and limitations continue in LEDGER_EVIDENCE.md. The replacement class-key and ideal-coupling proofs are in KEY_LEMMA.md. These files together form this submission's complete argument.
