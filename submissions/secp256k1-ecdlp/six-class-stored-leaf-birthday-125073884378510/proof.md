# Stored x-only leaf birthday search on secp256k1: six leaves per inverse, a two-fraction class label, metered rare-carry branches, k-record block messages, tbz slot test, fwd tbz (125.073884378510)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; on every run its total work is at most
`time_log2 = 125.073884378510 = 125073884378510/10^12` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit). Peak memory is at most 2^135.9989 bytes; preprocessing, at most 2^110.1508 units, is included;
there is no advice. `certificate.py` checks these declared figures (`DECLARED`): the work, memory and preprocessing
figures are the least multiples of 10^-12, 10^-4 and 10^-4 at or above their exact values, the work per leaf 324.1129
the least multiple of 10^-4 at or above R/S, and 0.390000001 is below the success bound.

**Reviewer's guide.** The bound rests on seven things.
- **A two-fraction class label (2.2, 3.3).** A class's label is the top 100 bits of the fractions `v*/p` and `(βv* mod
  p)/p`, `v* = min(x, βx, β²x)`. Each leaf has its x's three fractions cheaply: x's top bits; KEY's top column of `xβ/p`
  and one more, summed as words shifted right by 3 (15 instructions, no carry chain), within 68 below at 128 bits; and
  the third from the first two, within [−1, 70]. KEY's index j picks the pair (j, j + 1): 41 instructions for j = 0, 46
  for j = 1, 2. A member whose computed fraction lies near a block boundary, or whose x lies just below a multiple of
  2^156, may carry another label: a charged missed-collision row below 2^-19.8 per leaf; no leaf is dropped. The mass
  bound KB is a lattice count: the pairs `(v, βv mod p)` lie on a lattice with the short basis
  `(B, A - B), (A - B, -A)` of secp256k1's Frobenius `π = A + Bω`, so a box of label values holds at most `1 + 2^-19`
  times its share.
- **Six leaves per inverse (3.1).** `D = y_Y^2 - y_X^2 = x_Y^3 - x_X^3` is the product of the three denominators
  `x_Y - β^t x_X`, so one batched inverse of D gives the six leaves `x(φ^t(X) ± Y)`, at 1397 instructions per pair (LT
  adds `T = (p - β^t x_X) + (p - x_Y)` from stored negatives); the run checks every leaf of its batches of 1, 2, 5 and 9
  incidences against exact affine arithmetic.
- **Rare carries by branch, metered (3.4).** SHX, FDX, F5, FS and FADDW branch over the rest of their carry or borrow
  chain unless its carry or borrow occurred. The skipped instructions would only add or subtract 0, so every result is
  bit-identical to the inline chain's; only the flags differ, and no path reads them before setting them (checked). A
  fixup costs 3 or 4 instructions and adds C0 to a stack word. It needs a row or result word in a set of about 2^-32 of
  all words (2^-30 for folds), so a node's mean fixup count is below 2^28.35. A node starts a batch only while its count
  is at most `2^45 - 43q2`: the ledger charges `4*2^45` per node, and Markov's inequality charges the cap's failure
  (below 2^-16.6 per leaf).
- **No uncharged transfer or bookkeeping (2.1, 4).** Per-record steps are instructions of the node's own processor at
  their executed counts: the emit in `keys` (7 per record) at hop 1, `h2s` (7.25) at hop 2, and `bkt`, `pass4` and
  `scan_ranges` (14.546875), loop control included; one AND, BIC or BFI forms each record's table address and
  destination word. Records go into their destination's **record block** in the sender's port; after k = 128, a flush
  (7) stores its header, making it one message, or drops it. Every message pays at each hop its delivery as executed
  `route.s` code, `dlvk` (514: per record the two LDPs and two STPs of the one-record `dlv`, plus a header LDP and RET)
  or `dlv` (6), and its port's in-flight count, which DSB waits on (2): 7 + 523/128 and 7.25 + 523/128 per record.
  Blocks are kept or dropped whole, as if each record were sent alone (the keep lemma, section 4). Links keep no flags,
  rings or counters; every address, cursor, slot count and cap test is a counted instruction of the sender or of the
  node that later reads the data.
- **Success (section 6, in `success.md`):** Chen–Stein on cross-node leaf pairs with b3 = 0 (6.2; its one non-routine
  step is the
  non-proportionality claim, for six leaf types), usefulness (6.3), the defect table (6.4), whose batch row comes from
  the arithmetic-site table of 3.4 and whose cap rows bound every pair of leaves, dependent ones included, and the final
  inequality (6.5). `certificate.py` evaluates the terms with exact Fractions and prints the bound 0.3900000051 (the
  run's `success_lower_bound`).
- **Work (section 7):** `certificate.py` recomputes the ledger exactly; the run checks each kernel call's count against
  its constants (`counts_match_certificate`).
- **Machine (2.1; scope in section 9):** per classical-group-work-v2, loads, stores, addressing, control and transfers
  count (each delivery as `dlv` or `dlvk`), instructions not cycles (a stalled instruction counts once), all processors'
  work with coordination; memory has no scalar contribution (section 5).
- **Not relied on:** any heuristic, any prior acceptance or score, or extrapolation from the toy or the 4 × 4 run.

## 1. Overview

Each node computes the **leaves** `x(φ^t(X) ± Y)`, t = 0, 1, 2, for independent points X and Y, six per pair from one
batched inverse: Galbraith, Wang and Zhang (ePrint 2015/605) share one inverse between `X ± Y`, and φ(x, y) = (βx, y)
extends this to the three twists of X. Labelled by two fractions of its class's canonical x, each leaf goes in two
bulk-synchronous hops to a home node, where equal keys give the logarithm. This is an analytic bound, not a practical
break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each with one scalar AArch64 processor and a private memory
below 2^64 bytes. The ledger counts the instructions the processors execute and, per message, its delivery and its
in-flight count (below); storage and wiring carry no price. The processor implements the base scalar AArch64 instruction
set with FEAT_RNG (`MRS Xt, RNDR`, section 7) and FEAT_LSE (`LDADD`, `STADD`, sections 3.4 and 4); every counted step is
one such instruction, whatever memory accesses it performs.
- **Ports and links.** A node has a column port and a row port, device registers in its address space. A port's
  48-byte **message block** sends one record: the node stores at +0 a destination word, whose bits 3..38 are the
  destination index (a row of its column or a column of its row, its own included; the link ignores the other bits), a
  destination address at +8 and four data words at +16..+47, and the STP to +32 completes the message. Its **record
  blocks**, one per destination index d at `PB + 2^13 d` (PB = 2^49 for the column port, 2^50 for the row port), each
  hold a header (destination word, address) at +0 and k = 128 record slots, slot n at +32n; the STP of the header
  completes a k-record message. The link delivers a message: it writes its data (four words, or the k records in
  slot order) at that address of the destination's memory and reports the write to its port (the in-flight count's
  decrement below). It keeps no flags, rings or counters; its only state is the message it carries. A
  store that completes a message stalls until the link has taken it, which for a block means read it, and a stalled
  instruction counts once: a port holds at most one message, and no store refills a block the link has not read.
- **Delivery is charged.** A delivery is charged as the work of `route.s` code at each hop: `dlv`, 6
  instruction-equivalents per one-record message (three LDPs read its 48 bytes from the link's block, two STPs write the
  four data words at the address, RET), or `dlvk`, `4k + 2 = 514` per k-record message (an LDP of the header, per record
  two LDPs from its slot and two post-indexed STPs at the address, RET). Per record `dlvk` does `dlv`'s data work.
  Performed by a processor for each kind of message, a delivery is this code run as a loop, with B in place of RET
  (`dlvk` returns x0 to its start): its first LDP waits until x0, the link's window on the message, holds a complete
  one, as a sender's store waits for the link, and a stalled instruction counts once. So 6 or 514 per message is all of
  a delivery's work, with no call or poll; the destination index selects the link's output, which is wiring.
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

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145289358` and `q2
= 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the x-coordinates
of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and `S = 2^72*n0 = 2^126.699316`
in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and `n = λ^t*b_g`. A point Z ≠ O
has the class `{±Z, ±φ(Z), ±φ²(Z)}`, its x-coordinates `{x, βx, β²x}` (mod p) have the canonical member `v* = min(x, βx
mod p, β²x mod p)`, and its **label** is `(comp(v*), comp(βv* mod p))`, `comp(v) = floor(2^100 v/p)`: two 100-bit
components shared by the class (3.3). A record is `[K2, K3, w2, w3]` with K2 and K3 the top 64 bits of the first and the
second component, `w2 = i | k<<36` and `w3 = d2 | c<<36 | l<<39`, where `c = 1 + 2t + [s = -1]` in 1..6 is the leaf
type; its digits are d1 and d2, the components' low 36 bits (hop-1 row, hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3
mod 2^48`. At its home a record's source is `(w2 mod 2^36, j)`, where j is the index of the slot it arrived in (2.3).

### 2.3 Phases

A pair slot holds `cap_p = 2^22 - 2^12` records of 32 bytes, a multiple of k. A node in row i writes at hop 1 only into
slot i of R1, at `R1 + 32*cap_p*i`, at each node of its column, and a node in column j writes at hop 2 only into slot j
of R2 at each node of its row. Each port has two tables of 2^36 entries, one per destination d: the **LC** entry is the
next slot of d's record block (`init3` sets it to slot 1) and the **RC** entry the sender's cursor in its pair slot at d
(`tinit` sets it to the slot's base); the slot ends at the same address at every destination, so the sender keeps that
end in a register. The column port's LC table is at address 0 and the row port's at 2^39, so d's entry, 8d or 2^39 + 8d,
is also a destination word for d; the RC tables lie 2^43 above them (ORR of 2^43).

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators with their table
   entries, header constants, the LC and RC tables, the d4 counters (`zero4`) and the d3 cursor table (`init3`). Then Q
   is broadcast in 2 hops.
1. **Hop 1:** the batch kernels (section 3) run for each X_k; the emit writes each leaf's record into d1's block in the
   column port, and after slot k the block goes to `(d1, j)` at its RC cursor, or is dropped if that cursor is at or
   past the slot end (section 4). Then `fin` sends each row d's final block by the same test and writes its final
   cursor, clamped at the slot end (drops advance cursors too), into entry i of d's table E1; a barrier follows.
2. **Hop 2:** `h2s` reads slot i of R1 up to `E1[i]`, for every i, and writes each record through the row port's blocks
   to `(own row, d2)` by the same rule; `fin` sends the final blocks and writes each home's final cursor into its
   `E2[j]`; a barrier follows.
3. **Grouping:** `bkt` reads each slot of R2 up to its E2 entry and appends the record's address and its d4 counter's
address to its d3 bucket
   of `cap_3` entries, which ends at a multiple of 2^38, or drops the record if the bucket is full; a stored record
   counts its d4. `prefix4` forms the d4 cursors, `pass4` writes each entry's address into the pointer array at its d4
   cursor, bucket by bucket in d3 order, which leaves each d4 counter at its range's end, and `scan_ranges` lists,
   within each d4 range, adjacent records with equal K2 as candidates, up to 2^28; a full list ends it.
4. **Resolution,** in four sub-phases: a home `(h_r, h_c)` sends two requests per candidate to the sources `(g_r, g_c)`
   via the relay `(g_r, h_c)`, and replies of at most 8 records (one leaf's scalars and points) return via `(h_r, g_c)`.
   Each sub-phase is a count exchange: senders write each receiver their count and the epoch; after a barrier each
   receiver returns each sender's base, in sender order, or a drop mark if its records would pass the cap (2^29 requests
   at a source, 2^34 records at a relay); after a barrier senders write the records at their bases; a barrier ends it. A
   count goes to the sender's entry of the receiver's count table and a base to the receiver's entry of the sender's
   base table (fixed addresses); a drop mark makes the sender skip that receiver. A destination word is 8 times the
   destination's index (LSL, or ADD of 8 when stepping); a request carries its home's indices and its candidate's list
   position, which its reply returns. The home recomputes `L = [m]P + [n]Q` and `L'` from the records' leaf types, finds
   the unit u with `L' = uL` if any, and, if `n' - u*n ≠ 0`, outputs `x = (u*m - m')/(n' - u*n) mod r` after checking
   `[x]P = Q`.

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
  W2L), and from each product's six-word sum `Z` the leaf `Z + T`, folded (AD6, FDX), where `T = (p - β^t x_X) + (p -
  x_Y)` is the sum of two stored negatives (ADDS, three ADCS and CSET of its carry), `T < 2p`.

The header holds `p - x_X, y_X, p - βx_X, p - β²x_X, x_X/2, βx_X/2, 3x_X²/4, 3β²x_X²/4`; a Y entry holds `y², y, x, 3x²,
p - x`; the six leaves overwrite the 256-byte entry in the order `c - 1 = 2t + [s = -1]`. The arithmetic is F's row
schedule; the rare carries and borrows of SHX, FDX, F5, FS and FADDW are taken by branches over fixups that give the
inline chains' results (3.4), and its only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks the entries downwards and takes each entry's six leaves. For each leaf it forms KEY's column F1 and the
next column's shifted sum s (3.3), then decides KEY's index j: with `σ = F0 + F1`, if F0 ≥ F1 (CMP, B.HS) it branches to
the block JH, where j = 2 if `F1 + σ` carries (ADDS, B.HS to JT) and j = 1 otherwise; else j = 2 if `F0 + σ` carries
(ADDS, B.HS to JT) and j = 0 otherwise, which goes on inline. The j = 0 path forms `c6 = F1 + (s >> 61)` and emits the
pair (x's, A's); JH forms c6 and `N' = ~((x >> 128) + A)` (ADDS, ADC, MVN) and emits (A's, N''s); JT forms N' (ADDS,
ADC, ADD of `s >> 61`, MVN) and emits (N''s, x's). The **emit** forms 8*d1 by one AND (j = 0, 1) or BIC (j = 2) of the
register x1 = (2^36 - 1) << 3 with the word holding d1 (3.3: x2, s or the sum's low word) shifted right by 25, 22 or 25
bits: with the column LC table at address 0 this is both d1's LC entry and the destination word. LDADD advances the
entry by 32 and returns the record's slot; the emit forms w3 from d2, writes the record with two STPs and tests for slot
k (one TBZ of bit 12, section 4): 7 instructions. After slot k, FLUSH (7, or 6 for a dropped block) ends the block (section 4). The leaf
pointer, advanced 32
bytes per leaf by a post-indexed LDP, gives w3's high bits: entries are 256 bytes at a 2^34-aligned base, so `pointer <<
31` is `c<<36 | l<<39`. JH and JT add 5 to the key-path count in x18 and branch back. **Meter.** The exact index is
uniform on {0, 1, 2} for a leaf other than O, and KEY's index equals it outside a set of probability at most `d_wj` per
leaf (3.3), so a leaf takes JH or JT with probability at most `2/3 + 2d_wj`. A node starts a batch only while its count
is at most `5(cap_j - 6q2)`, `cap_j = ceil(2n0/3) + 2^42`. Two leaves of one incidence are an invertible linear image of
(X, Y), and leaves sharing one generator are independent given it up to one excluded point, so a node's number of such
leaves has mean at most `n0(2/3 + 2d_wj)` and variance at most `(2/9 + 2d_wj)n0 + 1` (the dependent pairs' covariances
are O(1/r) each; 6.4's meter row, while no batch of the node is defective). Executed counts per incidence (fixed part in
parentheses): `fwd` 98 (5) and `bwd` 1397 (5), each plus 3 or 4 per fixup (3.4), and `keys` 248 (9) for six leaves on
the j = 0 path, including their six emits of 7, plus the count it adds to x18 and its flushes; fe_inv (`secp_fe.s`)
is 23646. Per leaf `keys` executes two LDPs, F1 (11), the column (15), the index's ADD, CMP, B.HS, ADDS and B.HS, the
ADD of c6 and the emit: 41, plus the loop's SUB and TBZ per six leaves; a leaf on JH or JT executes 5 more. Kernels
may overwrite x0-x29 and the flags and preserve sp; `bwd` and fe_inv save and restore x30, which `fwd`, `keys` and
the transport kernels do not write; `fwd` and `bwd` return their fixup word in x0 (3.4); the driver keeps no
registers across calls and loads the key-path count into x18 for `keys`. `bwd` keeps R in x25..x28 from pass to pass; `fwd` ends with a TBZ on bit 33 of the entry pointer, ascending from
the 2^34-aligned base to the block end, while `bwd` and `keys` end with a TBZ on bit 33 descending: the 2^25 entries
of 256 bytes fill the lower 2^33 bytes of the block, so bit 33 is clear at every entry, set past the last one, and
set one entry below the block.

### 3.3 The class label (the index, the column, JH and JT)

**Fractions and the label.** For v in [0, p) let `T(v) = floor(2^128 v/p)` and `comp(v) = floor(2^100 v/p) = T(v) >>
28`. For a leaf value x in [1, p) let `v_j = β^j x mod p`; the canonical member `v* = min v_j` has index j = argmin, and
the class's **label** is `(comp(v*), comp(βv* mod p))`: the member `β^m x` has index `j - m` and the same pair. As
`2^128 v/p = v/2^128 + C0*v/(2^128 p)` with the second term in [0, 2^-95), `T(v) = (v >> 128) + e0(v)`, e0 in {0, 1};
and `comp(v) = v >> 156` unless `v mod 2^156 > 2^156 - 2^34`.

**The index.** This is F's canonical-index row without its windows. `G = floor(β*2^384/p)` has limbs G0..G5 (x16,
x7..x11), F0 is the top limb of x, `F1 = Σ_i [lo(x_i*G_(5-i)) + hi(x_i*G_(4-i))] mod 2^64` (column 5 of `x*G` without
its carry-in), `σ = F0 + F1`, `d = F0 - F1` and `e = min(F0, F1) + σ`. Then j = 2 if e carries, else 1 if F0 ≥ F1, else
0; 3.2's branches decide exactly this without forming the minimum. The **window set** is where (W1) `F1 >= 2^64 - 16`,
(W2) `σ >= 2^64 - 17`, (W3) `e mod 2^64 >= 2^64 - 33` or (W4) `d < 16`; it only describes where j may be wrong.

**Exactness of the index.** Let `V_i = 2^64*v_i/p` for `v0 = x`, `v1 = βx mod p` and `v2 = β²x mod p`, which are
distinct and sum to p or 2p. Then `F0 = V0 - a` with `0 <= a < 2`, and `a < 1` unless `x mod 2^192 > 2^192 - C0`. Bits
320..383 of `x*G` equal `F1 + c` and `floor(V1 - ε)` mod 2^64, with `0 <= c <= 8` (the terms left out of F1 sum below
8.01*2^320 mod 2^384) and `0 <= ε < 2^-64` (as `0 <= β*2^384/p - G < 1`). So outside (W1) `F1 = V1 - b` with `0 <= b <
10`, and outside (W1) and (W2) `2^64 - σ = V2 + a + b`. Hence j = 2 gives `V0, V1 >= V2`; j = 1 gives `V0 - V1 > 16 -
10` and `V2 - V1 > 34 - 2 - 20` by (W4) and (W3); j = 0 gives `V2 - V0 > 34 - 4 - 10` and `V1 - V0 >= 1 - a > 0` when `a
< 1`. Outside the window set j is therefore `argmin v_j`, except when j = 0 and `x mod 2^192 > 2^192 - C0`. **A wrong
index.** Where j is wrong the leaf emits another member's pair: a missed collision charged as the row `d_wj` of 6.4, or,
if it equals another class's label, a false candidate, which the mass bound below covers. **Rarity.** (W1) puts v1
within 2^196 of 0 or p; (W2) and (W4) put `(1±β)x mod p` within 2^197 of 0; (W3) puts `(2+β)x` or `(1+2β)x mod p` within
2^198 of 0. As these maps are bijections, at most `13*2^197` values of x lie in the window set and at most `2^64*C0`
more give the wrong index; a uniform point's x takes each value with probability at most `2/(r-1)`, so `d_wj <=
(13*2^197 + 2^64*C0)*2/(r-1) < 2^-54.29` per leaf.

**The column (15).** Column 4 of `x*G` has the eight words `lo(x_i*G_(4-i))` and `hi(x_i*G_(3-i))`. The kernel adds
`hi(x_0*floor(G3/8))` (x15 = G3 >> 3) and the other seven words shifted right by 3 (ADD with LSR #3) into s: each term
is below 2^61, so `s < 2^64`. Then `A = 2^64*F1 + 8s mod 2^128`, with top word `c6 = F1 + (s >> 61)`. **Bound.** As
`x*G/2^384 = xβ/p - xδ/2^384` with δ in [0, 1), `floor(2^128 {x*G/2^384})` is `T(v1)` or `T(v1) - 1` mod 2^128, and it
is columns 4 and 5 of `x*G` with the carry from columns 0 to 3. A leaves out that carry (columns 0 to 3 have 1, 3, 5 and
7 words, so their carries out are at most 0, 2, 4 and 6) and the bits the shifts drop: `Σ (t mod 8)` over seven words,
at most 49, and `hi(x_0 G3) - 8 hi(x_0 floor(G3/8))`, at most 12 as `G3 mod 8 = 5`. Hence `T(v1) - A ≡ e1 (mod 2^128)`,
e1 in [0, 68], for every x. **The third fraction.** As `v0 + v1 + v2` is p or 2p, `T(v2) ≡ -T(x) - T(v1) - k (mod
2^128)` with `k = ceil(f0 + f1)` in {0, 1, 2}, f0 and f1 the fractional parts of `2^128 x/p` and `2^128 v1/p` (the
fractional part of `2^128 v2/p` is in [0, 1) and the three multiples sum to a multiple of 2^128). `N' = ~((x >> 128) +
A) = -(x >> 128) - A - 1 mod 2^128`, and `x >> 128 = T(x) - e0(x)`, so `N' - T(v2) ≡ e0(x) + e1 + k - 1`, in [−1, 70].

**The paths.** A leaf of index j emits `(x >> 156, A >> 28)` for j = 0, `(A >> 28, N' >> 28)` for j = 1 and
`(N' >> 28, x >> 156)` for j = 2: its fields are each component's top 64 bits (K2, K3) and low 36 (d1, d2), so for j = 0
K2 = x3, d1 = x2 >> 28, K3 = c6 and d2 = bits 25..60 of s; for j = 1 K2 = c6, d1 = bits 25..60 of s, K3 = N''s top word
and d2 its bits 28..63 (EOR of the sum's bits 28..63 with 2^36 - 1); for j = 2 K2 = N''s top word, d1 its bits 28..63,
K3 = x3, d2 = x2 >> 28. With the index exact, `v_j = v*`, so the pair is the class's label unless a computed value's
error carries it across a block of 2^28 or x's own bits differ from comp: the **label-error set** D is where `T(v1) mod
2^28 < 70`, `T(v2) mod 2^28 >= 2^28 - 71` or `= 0`, or `x mod 2^156 > 2^156 - 2^34`. Outside D a leaf whose index is
exact emits its class's label. As each value of T is taken by at most 2^128 values and `x -> v_j` is a bijection, D has
at most `142*2^228 + 2^134` values of x, each of probability at most `2/(r-1)`: the row `d_lab < 2^-19.85` of 6.4, a
missed collision like `d_wj`; no leaf is dropped.

**Mass bound.** secp256k1's Frobenius is `π = A + Bω` with `B = 0xe4437ed6010e88286f547fa90abfe4c3` and
`A = 0x114ca50f7a8e2f3f657c1108d9d44cfd9`: `A² - AB + B² = p`, `B*β ≡ A - B` and `(A - B)*β ≡ -A (mod p)` (checked by
the run and the certificate). So `Λ = {(v, u) ∈ Z²: u ≡ βv (mod p)}` has the basis `(B, A - B), (A - B, -A)` of
determinant -p, and as `A > B > 0` its half-open parallelogram P spans [0, A) in v and (-A, A - B) in u. For each point
λ of Λ in a box `[v0, v0 + a) × [u0, u0 + b)` of integer sides a and b, `λ + P` lies in `[v0, v0 + a - 1 + A) × (u0 - A,
u0 + b - 1 + A - B)`; these translates are disjoint, as P tiles the plane under Λ, and each has area p, so the box holds
at most `(a + A)(b + 2A - B)/p` points of Λ. Fixing f of a label's 200 bits leaves `2^(200-f)` pairs of component
values; a component value is an interval of at most 2^156 values of v (of u = βv mod p for the second), widened at each
end by 72 units of `p/2^128` for a computed component (one that wraps mod 2^128 may be two end pieces of [0, p), which
count as one interval as `pZ² ⊂ Λ`). So on each path at most `2^(200-f)(2^156 + w + A)(2^156 + w + 2A - B)/p < (1 +
2^-19)2^-f p` values of `v_j`, and of x, give those bits, with w the widening (the certificate checks the factor),
whatever index KEY returns. With three paths and two points per x, the emitted pair takes given values on f of its bits
with probability at most `KB(f) = 6*2^(256-f)(1 + 2^-19)/(r-1)`.

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts
`B_(k+1) = SHX(B_k)`; RF and RA sum `Z = Σ a_k*B_k < 2^322` exactly (each high word is at most 2^64 - 2), and a fold
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ5 forms the square of a weak a as `V =
L + H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), AC4 may add a constant below p (top limb at
most C0 + 2), and F5 folds V the same way. AD6 forms `Z' = Z + T` exactly (top word at most 4). Only these
steps drop a carry or a borrow: SHX's last ADCS, the final ADC of FDX and of F5, and FS when its subtrahend is a
fold's output at or above p. Without a dropped carry each returns a value below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; the fused fold of
  `Z' < 2^322 + 2^257` likewise needs the leaf in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`, or Q_t in
  `[C0, C0² + 3C0)` with the constant added.
- **FS** `a - b` with `a < p` and `b < 2^256`: wrong iff `b - a > p`, which needs `b >= p`; the subtrahends Q_0 and Q_1
  are fold outputs, at or above p only when the true value is below C0.

A fold returns a value `>= p` only when the true value is below C0. That matters only where a canonical value is needed:
the leaf entering KEY and the key paths, and FS's subtrahend (above). Every other operation is exact on its inputs: RF
and RA; T and AD6; the label's column sums and path additions (3.3); FS and FADDW on canonical inputs; KEY's integer
arithmetic (3.3 covers its index); the emit; and F's `fe_mul_f` and `fe_sqr_f`, which fold the 512-bit product, fold its
fifth limb again and add C0 exactly when that last sum carries past 2^256, so fe_inv (255 squarings, 15 multiplications)
is exact for every input below 2^256. Generator construction, the table entries, the header values, the resolver and
`[x]P = Q` run on `fe_mul_f`, `fe_sqr_f` and fe_inv, each coordinate brought below p (an ADDS chain of C0 detects `t >=
p`, CSELs keep `t - p`) before it is subtracted, compared or stored, with exact arithmetic mod r for the solve. So they
have no sites; the tables hold exact canonical values (`p - x` too: no point has x = 0, as 7 is not a square mod p).

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
`bwd` beside its saved x30 (STP); each returns the word in x0, at most `43q2*C0 < 2^63` per call.
The driver adds x0/C0 to the node's fixup count N and starts a batch only while `N <= 2^45 - 43q2` (43 sites per
incidence: 18 SHX, 11 FDX, 7 F5, 4 FS and 3 FADDW). So N is at most 2^45 on every run, and section 7 charges `4*2^45`
per node.

**How often.** A fixup needs a word in a small set:
- SHX's carry needs the row B it shifts to have its low word above `2^64 - C0`, as the word it adds is below C0: at most
  `C0*2^192` rows.
- FDX's carry into limb 2 leaves its output with low 128 bits below `F*C0 <= (2^66 + 1)*C0` and the rest nonzero, where
  `F = floor(Z/2^256) <= 2^66 + 1` (`Z < 2^322`, or `Z' < 2^322 + 2^257` for the fused fold): at most
    `F_w = (2^66 + 1)*2^128*C0` words.
- F5's carry into limb 2 leaves its output's low 128 bits below 2^65, as its fifth word is at most C0 + 2, so that word
  times C0 has a high word of at most 1: at most 2^193 words.
- FS's borrow leaves an output whose low word is at least `2^64 - C0`, and FADDW's carry one whose low word is below C0:
  at most `C0*2^192` words each.

While no site of the batch up to and including this one has deviated, the row or output is a word below 2^256 congruent
to the site's operand value (times `2^(64k)` for the k-th row), which therefore determines that value. So, given X and
every other generator, a site takes its fixup with probability at most its set's size times the points of `Y_l` per
operand value over r - 1, with the fresh generator and multiplicities of the table above; `y_Y ± y_X` have 3 points per
value, and `x_Y + β^t x_X/2` has 2. At l = 0 the six sites the table leaves out cannot take a fixup either: rows 1, 2^64
and 2^128 have a low word of at most 1, folds of values below 2^256 add nothing, and `R_0*D_0` folds to 1 or p + 1,
neither of which is at least 2^128 with its low 128 bits below `2^66*C0`. Summed over a batch (`FIXUPS` and `E_FIX` in
`certificate.py`), the expected number of fixups a batch takes before its first deviation, and so of a batch without a
defect, is at most `E_FIX`, and over a node's batches at most `q1*E_FIX < 2^28.35`.

| fixup sites, at incidence l | runs per batch | words | points per value |
| --- | --- | --- | ---: |
| fwd: FS forming `D_l` | q2 | `C0*2^192` | 6 |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | `C0*2^192` | 6 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | `F_w` | 6 |
| bwd: row shifts of `R_l` | 3q2 | `C0*2^192` | 6 |
| bwd: folds of `inv_l` and `R_(l-1)`, l ≥ 1 | 2(q2 - 1) | `F_w` | 6 |
| bwd: FS of `y_Y - y_X`, and FADDW of `y_Y + y_X` | q2 each | `C0*2^192` | 3 |
| bwd: row shifts of `inv_l` | 3q2 | `C0*2^192` | 6 |
| bwd: folds of u and v | 2q2 | `F_w` | 3 |
| bwd: square folds of u² and v² | 2q2 | 2^193 | 6 |
| bwd: FADDW of `x_Y + β^t x_X/2`, t = 0, 1 | 2q2 | `C0*2^192` | 2 |
| bwd: square folds of Q_0 and Q_1 | 2q2 | 2^193 | 4 |
| bwd: two FS forming Q_2 | 2q2 | `C0*2^192` | 4 |
| bwd: square folds of Q_t² | 3q2 | 2^193 | 8 |
| bwd: row shifts of Q_t² | 9q2 | `C0*2^192` | 8 |
| bwd: fused folds of the six leaves | 6q2 | `F_w` | 2 |

**The meter row.** A node stops early only if N exceeds `2^45 - 43q2` at the start of some batch. On the event that none
of its batches has a defect, N is that clean count, so by Markov's inequality a node stops with probability at most
`q1*E_FIX/(2^45 - 43q2) + q1*d_batch < 2^-16.6` (`d_f`, a node-local row of 6.4).

## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (7), and per block of k = 128 records its flush (7) and `dlvk` (514) | 11.0703125 | 11.0859375 |
| hop 2: `h2s` (7; loop 0.25), and per block its flush (7; 8 if dropped) and `dlvk` (514) | 11.3203125 | 11.3359375 |
| bucket pass `bkt` (8; 7 if the bucket is full; loop 0.25) | 8.25 | 8.25 |
| local pass `pass4` (3; loop 0.25) | 3.25 | 3.25 |
| `scan_ranges`: 195/64 per record (K2 only, pointer pairs, CMP and CCMP per pair, 64 records per end test), plus 8 per listed candidate and at most 26 per d4 range | 3.046875 | 3.046875 |
| **total** | **36.9375** | **36.96875** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlvk`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. A destination's full blocks number at most its records over k, and a
dropped block costs 6 or 8 and sends nothing, so 523 per k records bounds each hop's blocks; final blocks are `fin`'s
rows (section 7). `h2s`, `bkt` and `pass4` take a slot's or bucket's records in groups
of 1, 2 and 4, as the low bits of its count say, then 8 at a time; the 8-record loop's CMP and B.LO, 2 per pass, are
charged as 0.25 per record (a slot of n records runs floor(n/8) passes), and the tests and set-up, 11, 12 and 11 per
slot or bucket, are per-node rows. **Addresses.** One instruction forms each per-record table address: the emit's AND or
BIC (3.2); `h2s`'s BFI of d2 (w3's low 36 bits) into bits 3..38 of x10, which holds 2^39 elsewhere (d2's row LC entry
and destination word); the flushes' ORR of 2^43 (the RC entry); `bkt`'s BFIs of K2's low 24 and K3's low 48 bits into
bits 3..26 and 3..50 of the d3 table's base (a multiple of 2^27) and the d4 counters' (2^51). **Blocks and cursors.**
The emit and `h2s` take a record's slot by one LDADD of 32 on its LC entry, which holds the absolute slot address:
block base (low 13 bits zero) plus offset 32n below 2^13, so bit 12 is set exactly after slot k (32k = 2^12). One TBZ
of bit 12 of LDADD's old value finds it (checked); there the **flush** returns
the LC entry to slot 1 (SUB, STR), forms the RC entry (ORR) and advances the RC cursor by 32k with LDADD, which returns
the old cursor c. If c is below the sender's slot end, which is the same at every destination and held in a register,
the STP of the header (destination word, c) sends the block; otherwise it is dropped (CMP and B.HS, or B.LO in `h2s`),
and `h2s` stores the end back, so a hop-2 cursor stays at most its slot end; a hop-1 cursor can pass its end, but only
by the node's own leaves, 32 bytes each, so it stays below 2^63.2 (checked). `fin` takes the L = (LC - slot 1)/32
records left in each destination's block, sends it by the header (8d, c) if L > 0 and c is below the end, and writes
the final cursor `min(c + 32L, end)` (CMP, CSEL). **Keep lemma.** Cursors start at the slot base and move by 32k per
block, and `cap_p` is a multiple of k, so c and the end are congruent mod 32k, and c below the end gives `c + 32k <=
end`: a block is kept or dropped whole, and as its records take the cursors c, c + 32, ..., the records kept, their
slot order and the final cursor equal those of sending each record alone at its own cursor (the one-record rule). A
final block's k - L slots after its L records hold older records or none and land in `[c + 32L, c + 32k)`, inside the
sender's slot and at or past the final cursor, where every later read stops. **Buckets.** Bucket t's `cap_3` entries
end at `R1 + (t+1)*2^38`, so the cursor LDADD returns to `bkt` has its low 38 bits zero exactly when the bucket is
full: TST and B.NE send a bucket with room to the store, and on a full bucket a STADD undoes the advance, so after each
record a cursor is at most its bucket's end. `bkt` steps its record pointer with a pre-indexed LDP, stores the record's
address with its d4 counter's and counts d4 with one STADD; `pass4` takes an entry's pointer-array cursor from that
counter with one LDADD of 8. **Scan** (`scan_ranges`: the range layout of submissions 80e9e8af and 1efcef88, the paired
form of 892a3df8 and b00400fa). `prefix4` sets each d4 counter to its range's start and `pass4` advances it per pointer
it writes, so after pass4 the counter is the range's exclusive end. `scan_ranges` visits the 2^48 ranges in order and
never compares across one. In a range of n >= 1 records, `n - 1 = 64q + 2r + o` with r < 32 and o in {0, 1}, it loads
the first record's K2, makes one single comparison (two LDRs, CMP, B.EQ) when o = 1, enters a block of 32 pairs through
a five-level TBZ tree (and a B when r is odd) for the r remaining pairs, and then runs the block q times. A pair loads
two pointers with one LDP and their K2s, compares the older adjacent pair with CMP and the newer with `CCMP ...,#4,NE`
(Z when either is equal) and branches once: 6 per two records, and 195 per 64 with the block's end test (MOV, CMP,
B.LO). Three key registers and four pointer pairs rotate through the block; the entry copies the range's first record
into every register a remainder can start from, so each comparison is between a record and the one before it in the
range, and as 32 mod 3 = 2 and 32 mod 4 = 0 one MOV realigns the keys for the next block, whose CMP with the range's
end (x18) stops it there, reading no pointer past the range. With range control a range without candidates costs `25 +
6o + (r mod 2) + 6r + 195q`, which exceeds `(195/64)n` by `1405/64 + 189o/64 - 3r/32 + (r mod 2)`, at most 25.8125 (o =
1, r = 1): 26 per range is charged; an empty range costs 5. A hit re-tests both comparisons on its out-of-line path
(CMP, B.NE each) and appends each equal pair in order after its list-end test (CMP, B.HS, STP), then branches back: 7,
8 or 11 for the first, the second or both, 4 for the single comparison, so at most 8 per listed candidate; the attempt
that finds the list full ends the scan. So a node's scan is at most `(195/64)n + 26*2^48 + 8(2^28 + 1) + 2`, the row of
section 7 (S is a multiple of 64). A hop's receiving node executes nothing on arrival; it reads the records in the next
phase, with the LDPs of `h2s` (hop 1) or `bkt` (hop 2). Per-node costs are rows of section 7, and every address a
kernel loads or stores lies in the regions of section 5.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 label bits. The records of one class form one run unless an **intruder**, a record
with the same 144 bits that does not carry the run's label, sits inside (6.4); a run lies in one d4 range with one K2,
so its adjacent pairs are all listed unless the candidate list is full. `scan_ranges` lists every adjacent pair of a
range with equal K2, which agrees on d1, d2, K2 and d4, 184 label bits: besides a run's pairs it can list pairs of
different classes, which the resolver rejects (no unit relates their points) and 6.4 charges to the candidate, request
and relay caps through KB(184). Writing the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i =
n_j/u_j`, an equivalence relation, so a run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums the regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.5912 used) | code, stack, the five barrier words (one reservation granule each), header, X table (eight header values and y² per X, 288 bytes), Y table (160 bytes per Y), scalars, batch entries (2^33 bytes at a multiple of 2^34), the ports' message blocks (at 2^49 - 64 and 2^49 - 128), their LC tables at 0 and 2^39 and RC tables at 2^43 and 2^43 + 2^39 (2^39 bytes each), the d3 cursor table at 2^40 (2^27 bytes), the record blocks at 2^49 and 2^50 (2^49 bytes each), the d4 counters at 2^51 (2^51 bytes), E1 and E2 (2^41 each) and the phase-4 buffers and tables (`FIXED` in `certificate.py`); all but the aligned tables, 2^43.51 bytes, lie between 2^44 and 2^49 | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | hop-1 pair slots | phases 1-2 |
| the same range, from phase 3 | d3 buckets in 2^24 strides of 2^38 bytes (`16*cap_3 = 2^37.7806` bytes at the end of each) and the pointer array `8*2^24*cap_3`, together 2^62.5155 | phase 3 |
| R2 after R1, ending at 2^63.998943 | hop-2 pair slots (phase 4 reads candidates' records at their addresses) | phases 2-4 |

So every node stays below 2^64 bytes, and the peak, in phase 2, is `2^72*(fixed + |R1| + |R2|)`, at most 2^135.9989
bytes (log2 ≈ 135.9988560). A remote item is reached only by a message.

## 6. Success (no heuristic)

In `success.md`: the model (6.1), Chen–Stein on cross-node pairs (6.2), usefulness (6.3), the defect rows (6.4) and the
bound (6.5).

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 676,627 per point:** a complete Jacobian
  double-and-add, exceptional-case tests, two inversions with the affine conversion, all on the exact arithmetic of 3.4:
  `fe_mul_f` and `fe_sqr_f` at 115 and 82 instructions as executed, plus 4 and 3 for the call as in fe_inv, and an
  addition or subtraction with its operands brought below p at 27; 665,910 for the 256 doublings, 256 mixed additions
  and two conversions; the exceptional-case tests (Z and H brought below p and tested for zero per iteration, 13 each,
  and 64 for the paths they select) add 6,720. The scalar's bit tests and the pointer step add at most 6 per iteration,
  and with set-up and return 1,600 in all. A Y's table entry (y², 3x², `p - x` by two LDPs of p's words, SUBS and three
  SBCS, and ten STPs, 240) or an X's mixed addition `R_g + [a_k]P` with its eight header values (three of them
  negatives, 6 each), y² and eighteen STPs (2,397) brings the total to at most 676,627, the allowance (charged at its
  itemization as in aburan28's 7e8d13c7 and ecc3ec1b).
- **Node set-up, 2^20:** a program of at most 2^16 bytes (4,676 kernel instructions, generator construction, resolver,
  phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD of 8 to its destination word, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 1,538,117:** two requests and two replies of at most 8 records over 2 hops at 40, two exact additions
  with inversion, the six-unit test, the solve mod r and `[x]P = Q`, on the same exact arithmetic; 769,058, charged at
  twice that plus 1 (as in 1fbd26e8).
- **Phase-4 count exchange, 2^41 per node and sub-phase:** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 88 per message (its count,
  the count message, the offset reply and the write, at most 22 each with in-flight count and delivery) for at most
  8*2^29 messages (2^29 requests of up to 8 records): below 2^40.5.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^10 per node:** 13 kernel calls between barriers, each at most 33 (eight argument registers of at
  most four MOVZ or MOVK each, and BL), the E-entry addresses and `scan_ranges`' arguments: below 2^9.
- **d4 prefix sums, 10 per four counters:** two LDPs, two STPs, four ADDs, SUBS and B.NE; the RET is the row's +1.
- **Header 27** (G0..G5: six words of MOVZ and three MOVK, and three STPs),
  **batch set-up 128** (X's eight header values copied into H by 16 LDPs and 16 STPs, y_X² into registers, acc = 1,
  pointers, the index word, the slot end for `keys`, both meters' updates and tests and the loop: at most 115, itemized
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
| coins: (q1 + q2 + 1)*64 (P) | 105.4141 |
| Y generators with their table entries: q2*676627 (P) | 116.3680 |
| R_g and X generators with their header values: (q1 + 1)*676627 | 118.4824 |
| header: 27 (P) | 76.7549 |
| LC and RC tables: 2*(4*2^36 + 1) + 2*(3*2^36 + 1) (P) | 111.8074 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| d4 counters: 2^48 + 1 (P) | 120.0000 |
| d3 cursor table: 4*2^24 + 1 (P) | 98.0000 |
| forward: q1*(98*q2 + 5) | 130.7291 |
| inversion: q1*(23646 + 16) | 113.6446 |
| backward: q1*(1397*q2 + 5) | 134.5625 |
| arithmetic fixups: 2^45*4 | 119.0000 |
| key pass without emits: q1*(206*q2 + 9) | 131.8009 |
| j = 1 and j = 2 key paths: cap_j*5 | 128.4366 |
| batch set-up: q1*128 | 106.1144 |
| hop 1, the emits and full blocks (flush, in-flight count, `dlvk`): S*(7 + 523/128) | 130.1700 |
| hop 2, `h2s` and full blocks: S*(7.25 + 523/128), plus 11*2^36 + 4 per node | 130.2022 |
| final blocks and cursors, in-flight counts and deliveries, both hops: 2*(541*2^36 + 2) | 118.0795 |
| bucket pass `bkt`: S*8.25, plus 12*2^36 + 4 per node | 129.7437 |
| d4 prefix sums: 2^48/4*10 + 1 | 121.3219 |
| local pass `pass4`: S*3.25, plus 11*2^24 + 2 per node | 128.3998 |
| scan_ranges: S*195/64, plus 26*2^48 + 2 + 8*(2^28 + 1) per node | 128.4205 |
| candidates: 2^28*1538117 | 120.5527 |
| requests served: 2^29*2^10 | 111.0000 |
| phase-4 count exchanges: 4*2^41 | 115.0000 |
| barriers: 32*2^82 | 87.0000 |
| phase control: 2^10 | 82.0000 |
| output | 20.0000 |
| **total R** | **135.039669** |

`R = 44770389438292344892248225799683571712000` instructions, at most 324.1129 per leaf (the run's `per_leaf`), and `W =
R/1000` with log2 W ≈ 125.0738843785. **The claim, 125.073884378510, is the least multiple of 10^-12 with 2^claim >=
W**; the certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds on ln 2. Rows marked P are
preprocessing: their sum over 1000 has log2 ≈ 110.15074 and is declared as 2^110.1508, the least multiple of 10^-4 at or
above it (checked).

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "trials": 6000, "curve_trials": 100, "seconds": 15}` (`seconds` 5 to
15). Budget-limited checks run guaranteed minimums (2,000 random and 600 near-tie class-rule inputs, 3,000 class-label
classes, 400 toy trials, 10 curve checks), then their share of the remaining budget. **Part A** runs the leaf kernels
and fe_inv on production batches of 1, 2, 5 and 9 incidences (`fwd` top-justified, its entries copied to the block
base for `bwd`/`keys`): every leaf against exact affine addition of `φ^t(X)` and `±Y`, each record's slot and LC entry against an exact model of the emit, every emitted label of a leaf outside KEY's
window set and D against its class's label, and KEY's index against the exact class rule, which it matches outside the
window set and the C0 case of 3.3 (`class_rule`; a batch of window-set leaves checks their records too). On at least
3,000 classes, a third built near a block boundary of T or of x's own bits, every member whose index is exact carries
the class label unless it lies in D, and the computed fractions' errors stay in [0, 68] and [−1, 70] (`class_labels`).
The `fwd` and `bwd` calls of these batches, of the planted-point batch and of the forced copy below count their formula
plus 3 or 4 per fixup they execute and return C0 times their fixups (`fwd_bwd_formulas_hold`). Planted points take
chosen fixups in both passes (fwd's FS and the next incidence's row shift, bwd's FS and two FADDWs) with exact leaves
and counts (`fixup_batch`). A copy of `leaf.s` whose fixup branches are never taken runs all 43 fixups in place on 2
incidences, as inline chains: exact leaves, 136 more instructions and 43 fixups per incidence (`forced_fixups`). 18
leaves near D's boundaries, at least 9 in D, give the emit model's records, key-path count and count
(`boundary_batch`). For each site type of 3.4 (`site_checks` in `sites.py`: row shift, fold, fused fold, the five
square folds, FS on a fold output), inputs built to meet its predicate deviate exactly when it holds, other boundary and
random inputs are exact, and all match its bit-level model; each fixup runs exactly when its carry or borrow occurs,
adds C0 to the word and 3 or 4 to the count, and its row or output lies in its set of 3.4 (`fixups`); FS and FADDW on
canonical inputs, `fe_mul_f` and `fe_sqr_f` are exact. A batch with `Y_l = φ(X)` and `Y_m = -X` gets acc = 0 and loses
all its other leaves; the last R update is ≡ 1 and unread; the SHX, FDX, F5, FS and FADDW invocations per pass match the
tables (`sites_match_certificate`), and every path from a fixup's join sets the flags before reading them or returns
(`fixups.flags_dead_at_every_join`). **Part B** (`machine.py`) runs a 4 × 4 machine through `leaf.s` and `route.s`
(`init3`, `tinit`, `zero4`; `keys` on 4 planted incidences of six leaves per node, labels with d1, d2, d3 < 4 and d4 <
16 at production field widths, each class's v* the point of 3.3's lattice nearest a chosen pair of component blocks,
some classes shared; `fin`; `h2s`; `fin`; `bkt`, `prefix4`, `pass4` and `scan_ranges` at each home), with pair slots of
k records (so its blocks are final blocks) and buckets of 8, delivering each message by executing `dlvk` or `dlv`,
which must return x0 to the link's window.
`block_paths` forces the other paths: `keys` with 390, 140, 128, 2 and 0 records to five destinations, the first three
on the j = 0, 1 and 2 paths, cursors starting 1, 2, 1, 1 and 3 blocks below the end (blocks sent and dropped, a cursor
past the end, final blocks sent with new, older and unwritten slots, held, and empty), then `fin`; `h2s` on slots of 0
to 19 records and two of 160, with 300, 140, 50 and 20 records to four full-width columns (blocks sent, and dropped with
the end stored back). Every call's count equals its formula from the certificate's constants (`scan_ranges`: its exact
path model, within the certificate's bound); messages, LC and RC entries, block, pair and bucket slots, pointer order
and candidates equal an exact model; delivered slots equal the one-record rule's, with every delivery store inside the
sender's slot; `unrolled_paths` runs `bkt` and `pass4` (full-width d3, d4) on slots and buckets of 0 to 19 records
(every group combination, up to two 8-record passes), some filling (a full bucket's LDADD undone, its cursor left at its
end), and `scan_ranges` alone on ranges of 0 to 129 records (every remainder of pairs with either parity, up to two
64-record blocks) with K2 from three values or random, on all-equal ranges, with a full list at a pair's first and
second append, and with no range: candidates and counts equal the path model, and every instruction of `scan_ranges`
runs. **Part C** (`toy.py`) runs the q25 toy in log form (the class of `[t]G` is `t^6 mod r`, φ is `[λ]`): 25 nodes of 4
× 4 incidences with six leaves each, and a curve pipeline computing each incidence's six leaves from one batched inverse
of D by 3.1's formulas, which must match the log form; the toy keys a class by its canonical x, and the production
label's 184 matched bits group the same leaves outside the charged wrong-index and label-error sets, besides pair
classes sharing those bits, which 6.4 charges as candidates and intruders. **Part D** runs `certificate.py`, and
`counts_match_certificate` compares its constants with parts A and B.

The toy and the 4 × 4 run are evidence at their own scope only; a run's toy figure has its own standard error.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1, whose links keep no state beyond the message in flight. The sender's
  instructions (including those that form the address), each delivery (as `dlv` or `dlvk`) and the ports' in-flight
  counts are charged. Each hop is one link from the sender's port to its destination; a network that forwarded a message
  through further nodes would add their work at each extra hop, and that machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's meter tests, and phase-4 exchanges), the in-flight counts, and the production attack, machine and
  links.
- **Provenance**, not a premise: package F (our 8f733826) for the row products (RF, RA, the pair schedule of W2S and
  W2L) and FS's CSEL form, with kaankolcu's c8a6819e, its base, for the macro set and `kernels/a64.py` (extensions in
  its docstring); aburan28's 02feb284 for fe_inv and the F0 and F1 rows with their error analysis (3.3's column extends
  F1); kaankolcu's 4595bad7 for a stored negation added in place of a subtraction (T's negatives); aburan28's fce6e7ef
  for charging link-side counter updates (the in-flight count); 3828873f, 6ef3e7b3, d637fb0a, 9798f5f4 and d9e3dbf4 for
  STADD in `h2s` and `bkt`, paired scan comparisons, `prefix4`, the generator and candidate allowances and a
  decimal-grid claim; 472226d8, 4fdb8a9c and e1ecd7cc for a success target nearer 0.39 and a 12-decimal claim; 80e9e8af
  and 1efcef88 for `scan_ranges`' range layout (K2-only comparisons in d4 ranges, the TBZ entry tree, the 2^28 and 2^29
  caps; first on 3828873f) and kaankolcu's 892a3df8 and b00400fa for its paired CMP/CCMP block of 64 records with 26 per
  range, unchanged; aburan28's 1fbd26e8 for a candidate at twice its itemized 769,058 plus 1, and 7e8d13c7 and ecc3ec1b
  for the generator charged at its itemization; kaankolcu's 807e4bf5 for charged rare exits running the original carry
  instructions, with a9297fb0's per-walk meter; ours: 10c5d938, the base (T, the address forms of section 4), and its
  ancestors bbe1cd85 (the two-fraction label of 3.3 and its lattice count), a3ca3235, 5b6cd013, eb0aea59, 502d7155,
  f9d5a673 and b3716a20 (six leaves per inverse, the success argument); Galbraith, Wang and Zhang for X ± Y from one
  inverse; Arratia, Goldstein and Gordon for Chen–Stein. Base: our e3a9665c (TBZ slot test on e59be87's k-record
  messages).
- **This package.** New here is the fwd loop-end TBZ: one TBZ of bit 33 replaces SUBS and B.NE at the end of `fwd`
  (3.2), saving one instruction per incidence; the run and `certificate.py` re-verify the changed code.
