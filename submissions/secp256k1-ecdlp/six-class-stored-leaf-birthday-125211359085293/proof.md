# Stored x-only leaf birthday search on secp256k1: six leaves per inverse, a fraction class key, metered rare-carry branches (125.211359085293)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; on every run its total work is at most
`time_log2 = 125.211359085293 = 125211359085293/10^12` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit). Peak memory is at most 2^135.9988 bytes; preprocessing, at most 2^110.2061 units, is included;
there is no advice. `certificate.py` checks these declared figures (`DECLARED`): the work, memory and preprocessing
figures are the least multiples of 10^-12, 10^-4 and 10^-4 at or above their exact values, the work per leaf 355.7798
the least multiple of 10^-4 at or above R/S, and 0.390000001 is below the success bound.

**Reviewer's guide.** The bound rests on seven things.
- **A fraction class key (2.2, 3.3).** The key is the top 176 bits of the class's canonical x-coordinate
  `v* = min(x, βx, β²x)`. A leaf whose own x is canonical (KEY's index j = 0, a third of leaves) reads it from x with no
  arithmetic; j = 1 and j = 2 extend KEY's top column of `xβ/p` by two more columns of a truncated product, and j = 2
  subtracts from x. Each path's 192-bit value is within 7 below (j = 1) or 10 above (j = 2) the exact fraction
  `floor(2^192 v*/p)`, so a test of its low 16 bits drops every leaf near a boundary, and every other member of exact
  index emits exactly `v* >> 80`. The drops are a charged row (below 2^-9.6 per leaf, a bound about 10 times the
  expected rate). KEY has no windows: a rare wrong index gives a leaf another member's fraction, a charged missed
  collision; a class sharing the 176 bits gives a false candidate or an intruder, both charged through the mass bound
  KB.
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
  charged at their executed counts: the emit inside `keys` (11 per record) at hop 1, `h2s` (10.25) at hop 2, and `bkt`,
  `pass4` and `scan_ranges` (17.546875), the unrolled loops' control included (section 4). Each message also pays 8 at
  each hop for its link: its delivery as `dlv` in `route.s` (6), which the run executes for every message of its 4 × 4
  machine, and its port's in-flight count, which DSB waits on (2). So the hops are charged 19 and 18.25. Links keep no
  flags, rings or counters. Every address, cursor, slot count and cap test is a counted instruction of the sender or of
  the node that later reads the data.
- **Success (section 6):** Chen–Stein on cross-node leaf pairs with b3 = 0 (6.2; its one non-routine step is the
  non-proportionality claim, for six leaf types), usefulness (6.3), the defect table (6.4), whose batch row comes from
  the arithmetic-site table of 3.4 and whose cap rows bound every pair of leaves, dependent ones included, and the final
  inequality (6.5). `certificate.py` evaluates the terms with exact Fractions and prints the bound 0.3900000047 (the
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
extends this to the three twists of X. Keyed by the top 176 bits of its class's canonical x, each leaf
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

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145590424` and `q2
= 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the x-coordinates
of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and `S = 2^72*n0 = 2^126.702303`
in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and `n = λ^t*b_g`. A point Z ≠ O
has the class `{±Z, ±φ(Z), ±φ²(Z)}`, its x-coordinates `{x, βx, β²x}` (mod p) have the canonical member `v* = min(x, βx
mod p, β²x mod p)`, and the **key** is `v* >> 80`, its top 176 bits (3.3), shared by the class. A leaf's key path
computes a 192-bit value with limbs `(lo, mid, hi)` whose top 176 bits are the key when its index is exact, unless it
drops the leaf (3.3). A record is `[K2, K3, w2, w3]` with `K2 = (hi mod 2^28) | (mid mod 2^36)<<28`, `K3 = lo`, `w2 = i
| k<<36` and `w3 = d2 | c<<36 | l<<39`, where `c = 1 + 2t + [s = -1]` in 1..6 is the leaf type; its digits are `d1 = hi
>> 28` (hop-1 row), `d2 = mid >> 28` (hop-2 column), `d3 = K2 mod 2^24` (hi's low 24 bits) and `d4 = K3 >> 16` (lo's top
48). They are 144 distinct key bits, and with K2 a record carries all 176. At its home a record's source is `(w2 mod
2^36, j)`, where j is the index of the slot it arrived in (2.3).

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
   at a source, 2^34 records at a relay); after a barrier senders write the records at their bases; a barrier ends it.
   The home recomputes `L = [m]P + [n]Q` and `L'` from the records' leaf types, finds the unit u with `L' = uL` if any,
   and, if `n' - u*n ≠ 0`, outputs `x = (u*m - m')/(n' - u*n) mod r` after checking `[x]P = Q`.

**Barriers:** 21 (program delivery 2, Q broadcast 2, the ends of phases 0 to 3 4, phase 4 12, output 1), within the
allowance of 32. **Caps:** arithmetic fixups (3.4), key paths (3.2), pair slots, d3 buckets, candidates, requests,
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

`keys` walks the entries downwards and takes each entry's six leaves. KEY (3.3) gives the class index j: its B.HS takes
j = 2 and a CBNZ takes j = 1 to their blocks, and j = 0 goes on to the emit with `(lo, mid, hi) = x >> 64`, x's limbs
1..3 already in x21..x23. The j = 1 block computes A1 (3.3) into x21..x23, adds 43 to the key-path count in x18, drops
the leaf if its window holds and branches back to the emit; the j = 2 block forms `A2 = -(x >> 64) - A1 mod 2^192`
there, adds 47 and tests its window. The **emit** forms `d1 = hi >> 28` and its cursor's address in the table at port +
64, advances the cursor by 32 with LDADD, which returns the old one, and drops the leaf if that is at or past the slot
end in x12 (5 instructions); otherwise it forms w3 from `d2 = mid >> 28`, inserts mid's low 36 bits above hi's low 28
(BFI) to form K2, and writes the message with three STPs (6). The leaf pointer, advanced 32 bytes per leaf by a
post-indexed LDP, gives w3's high bits: entries are 256 bytes at a 2^34-aligned base, so `pointer << 31` is `c<<36 |
l<<39`. **Meter.** The exact index is uniform on {0, 1, 2} for a leaf other than O, and KEY's index equals it outside a
set of probability at most `d_wj` per leaf (3.3), so a leaf adds 43, or 47, with probability at most `1/3 + d_wj` each.
A node starts a batch only while its count is at most `cap_w - 47*6q2`, `cap_w = ceil(30*n0) + 47*2^42`. Two leaves of
one incidence are an invertible linear image of (X, Y), and leaves sharing one generator are independent given it up to
one excluded point, so a node's count has mean at most `90*n0*(1/3 + d_wj)` and variance at most `(43² + 47²)(1/3 +
d_wj)*n0 + 47²` (the dependent pairs' covariances are of order 47²/r each) (6.4's meter row, while no batch of the node
is defective; a dropped leaf's path ran and is counted). Executed counts per incidence (fixed part in parentheses):
`fwd` 99 (5) and `bwd` 1403 (5), each plus 3 or 4 per fixup (3.4), and `keys` 188 for six leaves on the j = 0 path,
including their six emits of 11 (7), plus the count it adds to x18; fe_inv (`secp_fe.s`) is 23646. Per leaf `keys`
executes two LDPs, KEY (16), its B.HS and CBNZ and the emit: 31, plus the loop's SUB and TBZ per six leaves; a leaf
takes 43 more on the j = 1 path, 47 more on the j = 2 path (which skips the CBNZ), and a window drop skips the emit and
the branch back (12 fewer). Kernels may overwrite x0-x29 and the flags and preserve sp; `bwd` and fe_inv save and
restore x30, which `fwd`, `keys` and the transport kernels do not write; `fwd` and `bwd` return their fixup word in x0
(3.4); the driver keeps no registers across calls and loads the key-path count into x18 for `keys`. `bwd` keeps R in
x25..x28 from pass to pass; `bwd` and `keys` end with a TBZ on bit 33 of the entry pointer: the 2^25 entries of 256
bytes fill the lower 2^33 bytes of their 2^34-aligned block, so bit 33 is clear at every entry and set one entry below
the block.

### 3.3 The class key (KEY and the paths J1 and J2)

**Fractions.** For a leaf value x in [1, p) let `v_j = β^j x mod p` and `T(v) = floor(2^192 v/p)`. The canonical member
`v* = min v_j` has index j = argmin, and the class's key is `v* >> 80`. As `2^192 v/p = v/2^64 + C0*v/(2^64 p)` with the
second term in [0, 2^-31), `T(v) = (v >> 64) + e0(v)`, `e0(v)` in {0, 1}; write `z = v* >> 64`.

**KEY (16 instructions).** This is F's canonical-index row without its windows. `G1..G5` are limbs 1..5 of
`G = floor(β*2^384/p)`, F0 is the top limb of x, `F1 = Σ_i [lo(x_i*G_(5-i)) + hi(x_i*G_(4-i))] mod 2^64`, `σ = F0 + F1`,
`d = F0 - F1` and `e = min(F0, F1) + σ`. Then j = 2 if e carries, else 1 if F0 ≥ F1, else 0; KEY sets x29 to 1 when
F0 ≥ F1, and the B.HS after it takes j = 2 to its block. The **window set** is where (W1) `F1 >= 2^64 - 16`, (W2)
`σ >= 2^64 - 17`, (W3) `e mod 2^64 >= 2^64 - 33` or (W4) `d < 16`; it only describes where j may be wrong.

**Exactness of the index.** Let `V_i = 2^64*v_i/p` for `v0 = x`, `v1 = βx mod p` and `v2 = β²x mod p`, which are
distinct and sum to p or 2p. Then `F0 = V0 - a` with `0 <= a < 2`, and `a < 1` unless `x mod 2^192 > 2^192 - C0`. Bits
320..383 of `x*G` equal `F1 + c` and `floor(V1 - ε)` mod 2^64, with `0 <= c <= 8` (the terms left out of F1 sum below
8.01*2^320 mod 2^384) and `0 <= ε < 2^-64` (as `0 <= β*2^384/p - G < 1`). So outside (W1) `F1 = V1 - b` with `0 <= b <
10`, and outside (W1) and (W2) `2^64 - σ = V2 + a + b`. Hence j = 2 gives `V0, V1 >= V2`; j = 1 gives `V0 - V1 > 16 -
10` and `V2 - V1 > 34 - 2 - 20` by (W4) and (W3); j = 0 gives `V2 - V0 > 34 - 4 - 10` and `V1 - V0 >= 1 - a > 0` when `a
< 1`. Outside the window set j is therefore `argmin v_j`, except when j = 0 and `x mod 2^192 > 2^192 - C0`. **A wrong
index.** KEY itself drops no leaf (a key path can, below). Where j is wrong the leaf takes another member's path: its
key may differ from the class's, a missed collision charged as the row `d_wj` of 6.4, and may equal another class's key,
a false candidate that the mass bound below covers. **Rarity.** (W1) puts v1 within 2^196 of 0 or p; (W2) and (W4) put
`(1±β)x mod p` within 2^197 of 0; (W3) puts `(2+β)x` or `(1+2β)x mod p` within 2^198 of 0. As these maps are bijections,
at most `13*2^197` values of x lie in the window set and at most `2^64*C0` more give the wrong index; a uniform point's
x takes each value with probability at most `2/(r-1)`, so `d_wj <= (13*2^197 + 2^64*C0)*2/(r-1) < 2^-54.29` per leaf.

**j = 0.** The leaf is v*, and the emit reads `x >> 64 = z`: its top 176 bits are exactly `v* >> 80`.

**j = 1 (J1, 43).** `G' = floor(β*2^448/p)` has limbs g0..g6 with `g1..g6 = G0..G5` and `g0 = 0x03fde1631c4b80af`;
x7..x11 hold G1..G5 and x14, x15 hold g0, g1. KROW0 and three KROWs add, for i = 0..3, the low half of `x_i*g_(4-i)` and
the high half of `x_i*g_(3-i)` (column 4) and the low half of `x_i*g_(5-i)` and the high half of `x_i*g_(4-i)` (column
5), with their carries into a third word, to which an ADD brings KEY's F1 (column 6's own terms mod 2^64): `A1 = (c4 +
2^64*c5 + 2^128*F1) mod 2^192` for the exact column sums c4 and c5. **Bound.** As `x*G'/2^448 = xβ/p - xδ/2^448` with δ
in [0, 1), the fraction `{x*G'/2^448}` lies within 2^-192 below `v1/p` (mod 1), so `floor(2^192 {x*G'/2^448})` is
`T(v1)` or `T(v1) - 1` mod 2^192 (the latter wraps when `v1 < p*2^-192`); it equals columns 4 to 6 of `x*G'` with the
carry from columns 0 to 3. A1 leaves that carry out: columns 0 to 3 have 1, 3, 5 and 7 terms below 2^64, so their
carries out are at most 0, 2, 4 and 6. Hence `A1 ≡ T(v1) - e1 (mod 2^192)` with e1 in [0, 7], for every x. For j = 1,
`v1 = v*`, so `A1 ≡ z + d (mod 2^192)` with d in [-7, 1]. **Window.** J1 drops the leaf if `(A1 + 7) mod 2^16 < 8`, that
is `A1 mod 2^16` in `{0} ∪ [2^16 - 7, 2^16)`; otherwise `A1 - d` lies in A1's block of 2^16, so `z = A1 - d` and the
emitted `A1 >> 16` is `v* >> 80`.

**j = 2 (J2, 47).** As `v0 + v1 + v2` is p or 2p, `2^192(v0 + v1 + v2)/p ≡ 0 (mod 2^192)`. Writing `2^192 v/p = T(v) +
f(v)`, `T(v2) = -T(x) - T(v1) - k mod 2^192` with `k = ceil(f(x) + f(v1))` in {0, 1, 2}, as `f(v2)` is in [0, 1) and the
three multiples sum to a multiple of 2^192. J2 forms `A2 = -(x >> 64) - A1 mod 2^192` (ADDS, ADCS, ADC, then SUBS, SBCS
and SBC from XZR); as `x >> 64 = T(x) - e0(x)` and `A1 ≡ T(v1) - e1`, `A2 - T(v2) ≡ e0(x) + e1 + k (mod 2^192)`, in [0,
10]. For j = 2, `v2 = v*`, so `A2 ≡ z + d` with d in [0, 11]. J2 drops the leaf if `A2 mod 2^16 < 16`; otherwise `A2 -
d` lies in A2's block, so `A2 >> 16 = v* >> 80`.

So a leaf whose index is exact emits its class's key or is dropped; the run checks this on random and window-edge
classes and reports the largest errors it sees, which stay inside these bounds. **Drops.** A j = 1 drop needs `T(v1) mod
2^16` in the 15 residues `[2^16 - 7, 2^16) ∪ [0, 7]`, and a j = 2 drop `T(v2) mod 2^16` in the 26 residues `[2^16 - 10,
2^16) ∪ [0, 15]`. As `x -> β^j x mod p` is a bijection and each value of T is taken by at most 2^64 values (p < 2^256),
at most `41*2^240` values of x qualify, each of probability at most `2/(r-1)`: the row `d_win < 2^-9.64` of 6.4, which
counts every x, not only those of index j (about 2^-13 is expected: 8 or 16 low residues of 2^16 on two thirds of
leaves). **Mass bound.** Whatever index KEY returns, the emitted key takes given values on f of its 176 bits with
probability at most `KB(f) = 6*2^(256-f)*(1 + 2^-13)/(r-1)`: on the j = 0 path the bits are x's own, at most `2^(256-f)`
values; on the j = 1 and j = 2 paths, through the bijection `x -> v_j`, T(v_j) lies in at most `2^(176-f)` blocks of
2^16 widened by the path's error (7 or 10), at most `2^(256-f)(1 + 10*2^-16)` values; summed, at most `3(1 +
2^-13)2^(256-f)` values of x, each 2 points.

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts
`B_(k+1) = SHX(B_k)`; RF and RA sum `Z = Σ a_k*B_k < 2^322` exactly (each high word is at most 2^64 - 2), and a fold
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ5 forms the square of a weak a as `V =
L + H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), AC4 may add a constant below p (top limb at
most C0 + 2), and F5 folds V the same way. SB6 forms `Z' = Z - s'' + 2^256` exactly (top word at most 4). Only these
steps drop a carry or a borrow: SHX's last ADCS, the final ADC of FDX and of F5, and FS when its subtrahend is a
fold's output at or above p. Without a dropped carry each returns a value below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; the fused fold of
  `Z' < 2^322 + 2^256` likewise needs the leaf in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`, or Q_t in
  `[C0, C0² + 3C0)` with the constant added.
- **FS** `a - b` with `a < p` and `b < 2^256`: wrong iff `b - a > p`, which needs `b >= p`; the subtrahends Q_0 and Q_1
  are fold outputs, at or above p only when the true value is below C0.

A fold returns a value `>= p` only when the true value is below C0. That matters only where a canonical value is needed:
the leaf entering KEY and the key paths, and FS's subtrahend (above). Every other operation is exact on its inputs: RF
and RA; SB6; the key paths' integer rows and sums (3.3); FS
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
| hop 1: the emit in `keys` (11; 5 if the slot is full), and `dlv` (6) | 17 | 19 |
| hop 2: `h2s` (10; 9 if the slot is full; loop 0.25), and `dlv` (6) | 16.25 | 18.25 |
| bucket pass `bkt` (10; 8 if the bucket is full; loop 0.25) | 10.25 | 10.25 |
| local pass `pass4` (4; loop 0.25) | 4.25 | 4.25 |
| `scan_ranges`: 195/64 per record (K2 only, pointer pairs, CMP and CCMP per pair, 64 records per end test), plus 8 per listed candidate and at most 26 per d4 range | 3.046875 | 3.046875 |
| **total** | **50.796875** | **54.796875** |

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
(t+1)*2^38`, so the cursor LDADD returns to `bkt` has its low 38 bits zero exactly when the bucket is full: TST and
B.NE send a bucket with room to the store, and on a full bucket a STADD undoes the advance, so after each record a
cursor is at most its bucket's end. `bkt` steps its record pointer with a pre-indexed LDP, takes `d4 = K3 >> 16` (LSR)
and counts it with one STADD; `pass4` moves its d4 cursors with register-offset LDR and STR. **Scan** (`scan_ranges`:
the range layout of submissions 80e9e8af and 1efcef88, the paired form of 892a3df8 and b00400fa). `prefix4` sets each
d4 counter to its range's start and `pass4` advances it per pointer it writes, so after pass4 the counter is the
range's exclusive end. `scan_ranges` visits the 2^48 ranges in order and never compares across one. In a range of n >=
1 records, `n - 1 = 64q + 2r + o` with r < 32 and o in {0, 1}, it loads the first record's K2, makes one single
comparison (two LDRs, CMP, B.EQ) when o = 1, enters a block of 32 pairs through a five-level TBZ tree (and a B when r
is odd) for the r remaining pairs, and then runs the block q times. A pair loads two pointers with one LDP and their
K2s, compares the older adjacent pair with CMP and the newer with `CCMP ...,#4,NE` (Z when either is equal) and
branches once: 6 per two records, and 195 per 64 with the block's end test (MOV, CMP, B.LO). Three key registers and
four pointer pairs rotate through the block; the entry copies the range's first record into every register a remainder
can start from, so each comparison is between a record and the one before it in the range, and as 32 mod 3 = 2 and 32
mod 4 = 0 one MOV realigns the keys for the next block, whose CMP with the range's end (x18) stops it there, reading no
pointer past the range. With range control a range without candidates costs `25 + 6o + (r mod 2) + 6r + 195q`, which
exceeds `(195/64)n` by `1405/64 + 189o/64 - 3r/32 + (r mod 2)`, at most 25.8125 (o = 1, r = 1): 26 per range is
charged; an empty range costs 5. A hit re-tests both comparisons on its out-of-line path (CMP, B.NE each) and appends
each equal pair in order after its list-end test (CMP, B.HS, STP), then branches back: 7, 8 or 11 for the first, the
second or both, 4 for the single comparison, so at most 8 per listed candidate; the attempt that finds the list full
ends the scan. So a node's scan is at most `(195/64)n + 26*2^48 + 8(2^28 + 1) + 2`, the row of section 7 (S is a
multiple of 64). A hop's receiving node executes nothing on arrival; it reads the records in the next phase, with the
LDPs of `h2s` (hop 1) or `bkt` (hop 2). Per-node costs are rows of section 7, and every address a kernel loads or
stores lies in the regions of section 5.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 key bits. The records of one class form one run unless an **intruder**, a record
with the same 144 bits that does not carry the run's key, sits inside (6.4); a run lies in one d4 range with one K2, so
its adjacent pairs are all listed unless the candidate list is full. `scan_ranges` lists every adjacent pair of a range
with equal K2, which agrees on d1, d2, K2 and d4, all 176 key bits: besides a run's pairs it can list pairs of different
classes, which the resolver rejects (no unit relates their points) and 6.4 charges to the candidate, request and relay
caps through KB(176). Writing the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`,
an equivalence relation, so a run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums the regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.0087 used) | code, stack, the five barrier words (one reservation granule each), header, X table (eight header values and y² per X, 288 bytes), Y table (160 bytes per Y), scalars, batch entries (2^33 bytes at a multiple of 2^34), the two ports with cursor tables of 2^39 bytes each, E1 and E2 (2^41 each), the d3 cursor table, the d4 counters (2^51) and the phase-4 buffers and tables (`FIXED` in `certificate.py`) | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | hop-1 pair slots | phases 1-2 |
| the same range, from phase 3 | d3 buckets in 2^24 strides of 2^38 bytes (`16*cap_3 = 2^37.7828` bytes at the end of each) and the pointer array `8*2^24*cap_3`, together 2^62.5161 | phase 3 |
| R2 after R1, ending at 2^63.998943 | hop-2 pair slots (phase 4 reads candidates' records at their addresses) | phases 2-4 |

So every node stays below 2^64 bytes, and the peak, in phase 2, is `2^72*(fixed + |R1| + |R2|)`, at most 2^135.9988
bytes (log2 ≈ 135.9987678). A remote item is reached only by a message.

## 6. Success (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.496397569` (the run's `mu`).
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
index) and are emitted (no window drop), both survive the caps, no intruder enters their run, and the candidate,
request and relay caps are not reached first. Its run is then resolved, a useful adjacent pair (section 4) gives x, and
`[x]P = Q` is checked before output. The pair is useful iff `n' ≠ u*n`, where `n = λ^t*b_g` and `n' = λ^(t')*b_(g')` for
its two nodes; given all points, this fails with probability exactly 1/r. So `P(fail) <= P(W = 0) + E #{α : I_α = 1, α
useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (a batch defect of 3.4, a wrong class index, a window drop, either meter) is a function
of node g's coins, and L' is uniform given them, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs at most
`μ(1 + 1/r)` times its per-leaf rate.

**Caps and intruders.** Condition on α's four generators. A leaf with a generator outside them, a **free** generator,
is then uniform on at least r - 1 points, so at most `KB(f)` of its mass falls on given values of f key bits,
whatever index KEY returns (3.3's mass bound). The
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
  `other_h <= ρ(S²*KB(72)*KB(176) + (2S + 20)*KB(176))` (pairs at the home equal on d1, d2, K2 and d4, the 176 key bits
  `scan_ranges` matches; L and L' are there); a source receives at most `other_g <= 2ρ*n0*S*KB(176)` requests; a relay
  holds at most `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(176) + (2S + 20)*KB(176))` records (a request is one record and a
  reply at most 8). Markov with the caps less α's own entries (4 in the candidate list, 8 requests, 16 records) gives
  `d_cl`, `d_rq` and `d_rl`. Leaves of a defective batch elsewhere are not uniform; the whole-run row "a computational
  defect in any batch" excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-126 |
| wrong class index (3.3): KEY's window set, or x mod 2^192 > 2^192 - C0 | `d_wj` | < 2^-54.29 |
| key-path meter, Chebyshev with slack `47(2^42 - 6q2 - 1) - 90*n0*d_wj` | `d_w` | < 2^-30 |
| fixup meter, Markov on the mean `q1*E_FIX` with the cap `2^45 - 46q2` (3.4) | `d_f` | < 2^-16.6 |
| a key window drops the leaf (3.3) | `d_win` | < 2^-9.6 |
| pair slot at hop 1, and again at hop 2 (Chebyshev, `cap_p`) | `d_p`, counted twice | < 2^-19.9 each |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `ρ(S + 12)KB(144)` | `lam_int` | < 2^-14.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-17.4; < 2^-19.9; < 2^-19.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.9 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.9 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.390000001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 145590424, where the left side is at least 0.3900000047, and
the target is 10^-9 above the rules' 0.39.

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 2^20 per point:** a complete Jacobian
  double-and-add, exceptional-case tests, two inversions with the affine conversion, all on the exact arithmetic of 3.4:
  `fe_mul_f` and `fe_sqr_f` at 115 and 82 instructions as executed, plus 4 and 3 for the call as in fe_inv, and an
  addition or subtraction with its operands brought below p at 27; 665,910 for the 256 doublings, 256 mixed additions
  and two conversions. The scalar's bit tests, the pointer step and the return add at most 6 per iteration, 2,048 in
  all. A Y's table entry (y², 3x², x + C0 and ten STPs, 239) or an X's mixed addition `R_g + [a_k]P` with its eight
  header values, y² and eighteen STPs (2,379) brings the total to at most 670,337, within the allowance (1.5642x).
- **Node set-up, 2^20:** a program of at most 2^16 bytes (4,313 kernel instructions, generator construction, resolver,
  phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 1,538,117:** two requests and two replies of at most 8 records over 2 hops at 40, two exact additions
  with inversion, the six-unit test, the solve mod r and `[x]P = Q`, on the same exact arithmetic; 769,058, charged at
  twice that plus 1 (as in 1fbd26e8).
- **Phase-4 count exchange, 2^41 per node and sub-phase:** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 88 per message (its count,
  the count message, the offset reply and the write, at most 22 each with in-flight count and delivery) for at most
  8*2^29 messages (2^29 requests of up to 8 records): below 2^40.5.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^10 per node:** 11 kernel calls between barriers, each at most 33 (eight argument registers of at
  most four MOVZ or MOVK each, and BL), the E-entry addresses and `scan_ranges`' arguments: below 2^9.
- **d4 prefix sums, 10 per four counters:** two LDPs, two STPs, four ADDs, SUBS and B.NE; the RET is the row's +1.
- **Header 32** (G1..G5, g0 and g1: seven words of MOVZ and three MOVK, three STPs and an STR),
  **batch set-up 128** (X's eight header values copied into H by 16 LDPs and 16 STPs, y_X² into registers, acc = 1,
  pointers, the index word, the slot end for `keys`, both meters' updates and tests and the loop: at most 119, itemized
  in `certificate.py`), **inversion call 16, request served 2^10, output 2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on
  the failure flag, plus a range test. RNDR is this machine's coin source: it returns 64 fresh uniform bits, the
  algorithm coins of the rules' probability space, and never reports failure, so the charged B.EQ is never taken and
  coins add no failure row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, the barriers and the output. A record row uses
S, as each leaf gives at most one record.

| row (P: preprocessing) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1 + q2 + 1)*64 (P) | 105.4166 |
| Y generators with their table entries: q2*2^20 (P) | 117.0000 |
| R_g and X generators with their header values: (q1 + 1)*2^20 | 119.1173 |
| header: 32 (P) | 77.0000 |
| cursor tables: 2*(3*2^36 + 1) (P) | 110.5850 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| d4 counters: 2^48 + 1 (P) | 120.0000 |
| d3 cursor table: 4*2^24 + 1 (P) | 98.0000 |
| forward: q1*(99*q2 + 5) | 130.7467 |
| inversion: q1*(23646 + 16) | 113.6476 |
| backward: q1*(1403*q2 + 5) | 134.5716 |
| arithmetic fixups: 2^45*4 | 119.0000 |
| key pass without emits: q1*(122*q2 + 7) | 131.0481 |
| j = 1 and j = 2 key paths: cap_w | 131.6095 |
| batch set-up: q1*128 | 106.1173 |
| hop 1, the emits, in-flight counts and deliveries: q1*114*q2 | 130.9502 |
| hop 2, `h2s`, in-flight counts and deliveries: S*18.25, plus 11*2^36 + 3 per node | 130.8921 |
| final cursors, in-flight counts and deliveries, both hops: 2*(17*2^36 + 3) | 113.0875 |
| bucket pass `bkt`: S*10.25, plus 12*2^36 + 4 per node | 130.0599 |
| d4 prefix sums: 2^48/4*10 + 1 | 121.3219 |
| local pass `pass4`: S*4.25, plus 11*2^24 + 1 per node | 128.7898 |
| scan_ranges: S*195/64, plus 26*2^48 + 2 + 8*(2^28 + 1) per node | 128.4233 |
| candidates: 2^28*1538117 | 120.5527 |
| requests served: 2^29*2^10 | 111.0000 |
| phase-4 count exchanges: 4*2^41 | 115.0000 |
| barriers: 32*2^82 | 87.0000 |
| phase control: 2^10 | 82.0000 |
| output | 20.0000 |
| **total R** | **135.177143** |

`R = 49246444704129254176261145737992665563136` instructions, at most 355.7798 per leaf (the run's `per_leaf`), and `W =
R/1000` with log2 W ≈ 125.2113590853. **The claim, 125.211359085293, is the least multiple of 10^-12 with 2^claim >=
W**; the certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds on ln 2. Rows marked P are
preprocessing: their sum over 1000 has log2 ≈ 110.20607 and is declared as 2^110.2061, the least multiple of 10^-4 at or
above it (checked).

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "trials": 6000, "curve_trials": 100, "seconds": 15}` (`seconds` 5 to
15); it finds the package files beside itself, under files/ or under /snapshot/files. Budget-limited checks run
guaranteed minimums (2,000 random and 600 near-tie class-rule inputs, 3,000 class-key classes, 400 toy trials, 10 curve
checks), then their share of the remaining budget. **Part A** runs the leaf kernels and fe_inv on production batches of
1, 2, 5 and 9 incidences: every leaf against exact affine addition of `φ^t(X)` and `±Y`, each port message and cursor
against an exact model of the emit, every emitted key of a leaf outside KEY's window set against its class's `v* >> 80`,
and KEY's index against the exact class rule, which it matches outside the window set (`class_rule`; a batch of
window-set leaves checks their messages too). On at least 3,000 classes, a third at window edges (some with `v* mod 2^64
> 2^64 - C0 - 5`), every member whose index is exact emits `v* >> 80` or is dropped by its window, with the path errors
inside 0..7 and 0..10 (`class_keys`). The `fwd` and `bwd` calls of these batches, of the planted-point batch and of the
forced copy below count their formula plus 3 or 4 per fixup they execute, and return C0 times their fixups
(`fwd_bwd_formulas_hold`, `fixup_batch`, `forced_fixups`). A batch of planted points takes chosen fixups in both passes
(fwd's FS and the next incidence's row shift, bwd's FS and three FADDWs) with exact leaves and counts (`fixup_batch`). A
copy of `leaf.s` whose fixup branches are never taken runs all 46 fixups in place on 2 incidences, as inline chains: its
leaves are exact and it adds 148 instructions and 46 fixups per incidence (`forced_fixups`). A batch of planted leaves
whose key windows hold drops exactly those six leaves, three on each of the paths j = 1 and 2, with the emit model's
messages, cursors, key-path count and count (`window_drop_batch`, whose `formula`, the terms of `keys_formula`, equals
its count). For each site type of 3.4 (`site_checks` in `sites.py`: row shift, fold, fused fold, the five square folds,
FS on a fold output), inputs built to meet its predicate deviate exactly when it holds, other boundary and random inputs
are exact, and all match its bit-level model; each fixup runs exactly when its carry or borrow occurs, adds C0 to the
word and 3 or 4 to the count, and its row or output lies in its set of 3.4 (`fixups`); each key path, on random leaves
and on members of random and window-edge classes, matches its bit-exact model, adds 43 or 47 to the count and drops a
leaf exactly when its window holds, and a member whose index is exact emits `v* >> 80` unless dropped (`key_paths`); FS
and FADDW on canonical inputs (with fixups), FADDW of a canonical and a weak value, `fe_mul_f` and `fe_sqr_f` are exact.
A batch with `Y_l = φ(X)` and `Y_m = -X` gets acc = 0 and loses all its other leaves; the last R update is ≡ 1 and
unread; the SHX, FDX, F5, FS and FADDW invocations per pass match the tables (`sites_match_certificate`), and every path
from a fixup's join sets the flags before reading them or returns (`fixups.flags_dead_at_every_join`). **Part B** runs a
4 × 4 machine through `leaf.s` and `route.s` (`tinit`, `init3`, `zero4`; `keys` on 4 planted incidences of six leaves
per node, d1, d2 < 4 and d4 < 16, some classes shared; `fin`; `h2s`; `fin`; `bkt`, `prefix4`, `pass4` and `scan_ranges`
at each home), with pair slots of 7 and buckets of 8, delivering each message by executing `dlv` on it (`counts.dlv` is
6 × `deliveries`). Every call's count equals its formula from the certificate's constants (for `scan_ranges`, its exact
path model, within the certificate's per-record and per-range bound); messages, slots, bucket entries, pointer order and
candidates equal an exact model; full slots, a full bucket (its LDADD undone, its cursor left at the bucket's end) and a
full list are forced; `unrolled_paths` runs `h2s`, `bkt` and `pass4` on slots and buckets of 0 to 19 records (every
block combination, up to two 8-record passes), some filling, and `scan_ranges` alone on ranges of 0 to 129 records
(every remainder of pairs with either parity, up to two 64-record blocks) with K2 from three values or random, on ranges
all equal, with a full list at a pair's first and at its second append, and with no range: candidates and counts equal
the path model, and every instruction of `scan_ranges` runs. **Part C** (`toy.py`) runs the q25 toy in log form (the
class of `[t]G` is `t^6 mod r`, φ is `[λ]`): 25 nodes of 4 × 4 incidences with six leaves each, and a curve pipeline
computing each incidence's six leaves from one batched inverse of D by 3.1's formulas, which must match the log form;
the toy keys a class by its canonical x; the production key's 176 bits group the same leaves outside the charged
wrong-index and window-drop sets, and besides pair classes that share those bits, which 6.4 charges as candidates and
intruders. **Part D** runs `certificate.py`, and `counts_match_certificate` compares its constants with parts A and B.

**Census.** `census-toy.json` is b3716a20's, reproduced exactly by this package's code: an off-runner run (`seconds`
1500, above the hosted maximum of 15; each toy count reaches its request). Over 100000 toy trials, success is 0.3953
with all pairs and 0.39027 with cross-node pairs only, against Poisson 0.38997 at the cross-node μ = 0.494252 and a toy
Chen–Stein bound of 0.0407, and the curve pipeline matches the log form in 1000 of 1000 trials. The toy and the 4 × 4
run are evidence at their own scope only; a declared run's toy figure has its own standard error.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1, whose links keep no state beyond the message in flight. The sender's
  instructions (including those that form the address), each delivery (as `dlv`) and the ports' in-flight counts are
  charged. Each hop is one link from the sender's port to its destination; a network that forwarded a message through
  further nodes would add their work at each extra hop, and that machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's meter tests, and phase-4 exchanges), the in-flight counts, and the production attack, machine and
  links.
- **Provenance**, not a premise: package F (our 8f733826) for the row products (RF, RA, and the pair schedule of W2S and
  W2L) and FS's CSEL form, with kaankolcu's c8a6819e, its base, for the kernel macro set and `kernels/a64.py` (whose
  extensions its docstring lists), aburan28's 02feb284 for fe_inv and the F0 and F1 rows with their error analysis,
  kaankolcu's 4595bad7 for SB6's subtraction of a stored offset, and aburan28's fce6e7ef for charging link-side counter
  updates as instruction-equivalents (the in-flight count); submissions 3828873f, 6ef3e7b3, d637fb0a, 9798f5f4 and
  d9e3dbf4 for STADD in `h2s` and `bkt`, paired scan comparisons, `prefix4`, the generator and candidate allowances and
  a claim on a decimal grid; submissions 472226d8, 4fdb8a9c and e1ecd7cc for a success target nearer 0.39 and a
  12-decimal claim; submissions 80e9e8af and 1efcef88 for `scan_ranges`' range layout (K2-only comparisons within d4
  ranges, a TBZ entry tree, the 2^28 and 2^29 caps; its first form was developed on 3828873f), and kaankolcu's 892a3df8
  and b00400fa for its paired CMP/CCMP block of 64 records with the bound of 26 per range, taken unchanged; aburan28's
  1fbd26e8 for charging a candidate at twice its itemized 769,058 plus 1; submission f9d5a673 (ours) for the LDADD
  cursors and the mask-tested buckets, and b3716a20 (ours) for the six-leaf `fwd` and `bwd`, the site table and the
  success argument, all unchanged here but for the key, the scan and the fixups; kaankolcu's 807e4bf5 for rare exits
  whose stubs run the original carry instructions and are charged, with a9297fb0's per-walk meter, which our 8f733826
  used in SH and FD and our 34a30cf5 in fe_sub; our 502d7155, eb0aea59 and 5b6cd013 for the windowless KEY, a class key
  with the range scan, and the metered fixups, kept here; Galbraith, Wang and Zhang for X + Y and X - Y from one
  inverse; Arratia, Goldstein and Gordon for Chen–Stein. The six-leaf sharing of 3.1 and the fixups of 3.4 are ours
  (b3716a20, 5b6cd013); new here are the fraction key of 3.3 (the j = 0 key read from x, the extension of KEY's column
  for j = 1 and 2 with its error bounds and windows, the drop and mass bounds and the weighted key-path meter) and the
  record layout of 2.2 that carries all 176 key bits.
