# Stored x-only leaf birthday search on secp256k1: six leaves per inverse labelled by truncated fraction products, grouped hash-chain abort, unequal rung ladder (124.677272491447240)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; its expected total work is at most
`time_log2 = 124.677272491447240 = 124677272491447240/10^15` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit), and on every run at most R_worst = 2^188.698089 (the certificate's least 6-decimal upper bound; disclosed: all 114 groups run, every hash walk
visits its whole chain of at most HOME_MAX links). Peak memory is at most 2^135.9989 bytes; preprocessing, at most
2^118.0361 units, is included; there is no advice.
`certificate.py` checks these declared figures (`DECLARED`): the work, memory and preprocessing
figures are the least multiples of 10^-15, 10^-4 and 10^-4 at or above their exact values, the work per leaf 246.1820
the least multiple of 10^-4 at or above E/S, and 0.390000001 is below the success bound. The bound rests on:
- **Early abort with E-priced group costs (2.3, 6.6, 7).** Group g's leaf and route work is exact given g and runs
  iff T > g - 1; its main-pair hash-walk work has mean bounded uniformly given the history H_g = (X_<g, Y, b, Q), over
  which the stopping indicator is fixed, so E[work] = E_pre + sum_g PHI[g] W_g + same-X/K2/list/defect rows +
  worst-case rest by history-conditional domination with a validity split (no independence is claimed); PHI[g] is
  the Chen-Stein prefix upper tail plus prefix defect rows (certificate `PHI`, nonincreasing, PHI[0] = 1); 3.2's
  key paths are E-priced the same way (6.6). Caps, preprocessing, resolution, barriers and control are worst-case. Stopping needs a resolver-confirmed x ([x]P = Q
  checked), which latches; E1-E7 give the same success event as the ungrouped run.
- **Order-free grouped detection by hash chains (4).** `zero4096` zeroes the 2^57 head words once; `hchg` reads each
  slot from its read cursor to E2 and advances the cursor, links the new record to the old head through the packed
  next-array, walks the K2 chain (R2 addresses) to the first full match or the chain end (uncapped), verifies d4,
  lists unless full, and makes the record the head. Same-key records form arrival-order chains, so every colliding
  pair with no intruder between its records is listed on every run.
- **Six leaves per inverse, never formed (3.1).** One batched inverse of `D = y_Y^2 - y_X^2` gives `u = 1/(y_Y +
  y_X)` and `v = 1/(y_Y - y_X)`; with E, F, G of 3.1 and w = u² or v², the leaves `x(φ^t(X) ± Y) = (Ew - b) + β^t(Fw -
  a) + β^(2t)Gw` are a Fourier transform over the cube roots of unity, and `bwd` stores two 128-bit fractions of each
  component: 1285.5 instructions per pair.
- **Fractions by truncated products (3.3).** On the rows `c*2^(64k) mod p` of c = E, F or G and a multiplier M, columns
  2 and 3 of `Σ M_k*row_k` (16 single-word products) give `T(cM mod p)` less at most 64; M = w and βw give a
  component's fractions. Fractions mod 1 add and β rotates them, so a leaf's are sums of its components' (with
  correction words for the linear terms), within [−132, 194], [−132, 194] and [−389, 264] at 128 bits; KEY's index j on
  their top words picks the pair (j, j + 1). A leaf in the window set W or the label-error set D' may carry another
  label: missed-collision rows below 2^-58.47 and 2^-16.65 per leaf; no leaf is dropped. A lattice count from
  secp256k1's Frobenius `π = A + Bω` bounds a box of label values at `1 + 2^-17` times its share.
- **Rare carries by branch, metered (3.4).** SHX, SHXT, FDX, F5, FS and FADDW branch over the rest of a carry or borrow
  chain unless it occurred: results are bit-identical to the inline chain's, and no path reads the flags before setting
  them (checked). A fixup costs 3 or 4 and adds C0 to a stack word; a node's mean fixup count is below 2^28.68 and its
  cap `2^45 - 33q2`, so the ledger charges `4*2^45` per node and Markov's inequality the cap's failure (below 2^-16.32
  per leaf).
- **No uncharged transfer or bookkeeping (2.1, 4).** Per-record work is executed code: the emit (7) at hop 1, `h2s`
  (7.25) at hop 2, `hchg` 7.25 fixed per record plus expected chain walks (6.5 per expected visit) bounded by
  history-conditional KB(129) pair counting (same-X pairs carved out globally, defective inserts allowed globally),
  the j = 1, 2 key paths of 3.2 bounded the same way (6.6), K2 verifications and listings
  charged globally, and one head-table zeroing plus one read-cursor init charged in preprocessing (next words need no
  zeroing: each is written at its record's insert, before any read). A block of k = 256 records is one message,
  charged as executed `dlvk` (1026) plus its in-flight count (2), kept or dropped whole (the keep lemma, 4, also at
  cap_pg). Links keep no flags, rings or counters.
- **Success, work, machine.** Chen–Stein on cross-node leaf pairs, usefulness, the defect rows and the final inequality
  (6, in `success.md`), evaluated with exact Fractions (grouped bound 0.3900000015, the run's `success_lower_bound`);
  the E-ledger recomputed exactly, each kernel call's count checked against its constants (7,
  `counts_match_certificate`); loads, stores, addressing, control and transfers count, instructions not cycles (2.1,
  5; scope in 9).
- **Not relied on:** any heuristic, any prior acceptance or score, or extrapolation from the 4 × 4 run.

## 1. Overview

Each node computes the **leaves** `x(φ^t(X) ± Y)`, t = 0, 1, 2, for independent points X and Y, six per pair from one
batched inverse: Galbraith, Wang and Zhang (ePrint 2015/605) share one inverse between `X ± Y`, and φ(x, y) = (βx, y)
extends this to the three twists of X. The six leaves are held as their Fourier components over the cube roots of
unity and labelled without being formed, by truncated fraction products of the components: by two fractions of its
class's canonical x, each leaf goes in two
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
  48-byte message block sends one record: a destination word at +0 (its bits 3..38 select a row of the node's column
  or a column of its row, its own included; the link ignores the other bits), a destination address at +8 and four
  data words at +16..+47, completed by the STP to +32. Its record blocks, one per destination index d at
  `PB + 2^14 d` (PB = 2^50 column, 2^51 row), hold a header (destination word, address) and k = 256 slots of 32
  bytes, slot n at +32n; the header's STP completes a k-record message. The link writes a message's data at that
  address of the destination and reports the write to its port; it keeps no flags, rings or counters. A store that
  completes a message stalls until the link has taken it, and a stalled instruction counts once.
- **Delivery is charged.** A delivery is `route.s` code run by a processor as a loop: `dlv`, 6 per one-record
  message (three LDPs read its 48 bytes, two STPs write the four data words, RET), or `dlvk`, `4k + 2 = 1026` per
  k-record message (header LDP, then per record two LDPs from its slot and two post-indexed STPs, RET); `dlvk` does
  `dlv`'s data work per record. Its first LDP waits until the link's window holds a complete message, and a stalled
  instruction counts once, so 6 or 1026 is all of a delivery's work, with no call or poll; the destination index
  selects the link's output, which is wiring.
- **Completion.** A port's in-flight count rises when the link takes a message and falls when the destination's
  write completes; both updates are charged, 2 per message (sections 4 and 7). DSB completes once the node's stores
  are written and its ports' counts are zero, and a node executes DSB before joining a barrier, so after a barrier
  every message sent before it is in memory.
- **Barriers.** A combining tree of epoch writes (the children of `(i, j)` are `(2i+1, j)`, `(2i+2, j)` and, for
  i = 0, `(0, 2j+1)` and `(0, 2j+2)`): a node waits for its children's arrival words, executes DSB, writes its
  parent's, then waits for its release word and writes its children's. A wait is LDAXR, CMP, B.EQ and WFE: a store
  into the word clears the node's exclusive monitor and wakes the WFE, each word has its own reservation granule,
  and nothing else wakes it.

### 2.2 Leaves and records

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145307658` and `q2
= 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the x-coordinates
of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and `S = 2^72*n0 = 2^126.699497`
in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and `n = λ^t*b_g`. A point Z ≠ O
has the class `{±Z, ±φ(Z), ±φ²(Z)}`, its x-coordinates `{x, βx, β²x}` (mod p) have the canonical member `v* = min(x, βx
mod p, β²x mod p)`, and its **label** is `(comp(v*), comp(βv* mod p))`, `comp(v) = floor(2^100 v/p)`: two 100-bit
parts shared by the class (3.3). A record is `[K2, K3, w2, w3]` with K2 and K3 the top 64 bits of the first and the
second part, `w2 = i | k<<36` and `w3 = d2 | c<<36 | l'<<39`, where c = 2 and l' = l for s = +1, and c = 5 and
`l' = l - 1 mod 2^25` for s = −1 (3.2): w3 gives the sign, not t. Its digits are d1 and d2, the parts' low 36 bits
(hop-1 row, hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3 mod 2^46`. At its home a record's source is `(w2 mod 2^36,
j)`, where j is the index of the slot it arrived in (2.3).

### 2.3 Phases

A pair slot holds `cap_p = 2^22 - 2^12` records of 32 bytes, a multiple of k. A node in row i writes at hop 1 only into
slot i of R1, at `R1 + 32*cap_p*i`, at each node of its column, and a node in column j writes at hop 2 only into slot j
of R2 at each node of its row. Each port has two tables of 2^36 entries, one per destination d: the **LC** entry is the
next slot of d's record block (`init3` sets it to slot 1) and the **RC** entry the sender's cursor in its pair slot at d
(`tinit` sets it to the slot's base); the slot ends at the same address at every destination, so the sender keeps that
end in a register. The column port's LC table is at address 0 and the row port's at 2^39, so d's entry, 8d or 2^39 + 8d,
is also a destination word for d; the RC tables lie 2^43 above them (ORR of 2^43).

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators with their table
  entries, the rows of β, the row port's LC/RC tables, the hash-head zeroing (`zero4096` over 2^57 words) and
  the hash-read cursors (`init3` of each slot's base into its 8-byte entry at 2^48). Then Q is broadcast in 2 hops
  and R_g is generated.
1-4. **Groups** g = 1..114, rung g starting at round(q1*(g/114)^0.65) (unequal sizes; 6.6, `changes.md`), each over a
  contiguous range of X's (all nodes in lockstep; the group loop tests the
  latched stop flag first and halts the run if set). Per group: generate its X's; run the batch kernels (section 3)
  for them; hop 1 into fresh per-group slots of cap_pg = cap_p/8 records (column LC/RC reset by `init3`/`tinit`);
  hop 2 into the cumulative R2 slots of cap_p (row cursors persist); hash-chain detection over the records appended
  since the previous group call (`hchg` reads each slot from its read cursor to E2 and stores E2 back, section 4);
  the four resolution sub-phases for the group's new candidates; then the stop test (below).
  Barriers separate every phase of every group.

1. **Hop 1 (per group):** the batch kernels run for each X_k; the emit writes each leaf's record into d1's block
   in the column port, and after slot k the block goes to `(d1, j)` at its RC cursor, or is dropped if that cursor is
   at or past the slot end (section 4). `fin` sends each row d's final block by the same test and writes its final
   cursor, clamped at the slot end (drops advance cursors too), into entry i of E1.
2. **Hop 2 (per group):** `h2s` reads slot i of R1 up to `E1[i]` for every i and writes each record through the row
   port's blocks to `(own row, d2)` by the same rule into the cumulative R2; `fin` sends the final blocks and writes
   each home's final cursor into its `E2[j]`.
3. **Grouping (per group):** `hchg` reads each slot of R2 from its read cursor to E2, stores E2 back, and hash-joins
   the new records: bucket K2 mod 2^57, R2-address chains via the packed next-array (one 8-byte word per slot
   position, at (A>>2)+C), uncapped walks, d4 verified on a K2 hit, the pair listed unless the cumulative list (up
   to 2^28) is full, then prepended. Per-group [cursor, E2) ranges partition each slot (every record read and
   inserted once); the next-array covers every receive pattern up to the slot caps by construction. No buckets,
   passes, snapshots or scans (section 4, 6.6).
4. **Resolution (per group),** four sub-phases: a home `(h_r, h_c)` sends two requests per new candidate to the
   sources `(g_r, g_c)` via the relay `(g_r, h_c)`; replies of at most 8 records (one leaf's scalars and points)
   return via `(h_r, g_c)`. Each sub-phase is the base count exchange (counts, bases in sender order, drop marks at
   2^29 requests at a source and 2^34 records at a relay, three barriers). The source serves a request from its
   tables (X_k's affine point from the retained all-groups X affine store, Y_l from the Y table, the scalars) inside
   2^10, and old groups serve identically. As w3 gives s but not t, the home forms each record's three possible
   leaves `L = φ^t(X_k) + s*Y_l = [m]P + [n]Q`, t = 0, 1, 2, tests the nine pairs in a fixed order for a unit u with
   `L' = uL`, and from the first with `n' - u*n ≠ 0` outputs `x = (u*m - m')/(n' - u*n) mod r` after checking
   `[x]P = Q`.
5. **Stop test (per group):** if any home resolved an x this group, the stop flag (latched at a fixed address, set
   by the resolving home) broadcasts over 2 hops and every node halts after the barrier; the latched x is the
   output. The test is 5 instructions per group per node (LDR, CMP, B, group-loop ADD/STR), charged in W_g.

**Barriers:** at most 1943 (pre 4, per group at most 17, output 1), within the allowance of 2048. **Caps:** arithmetic fixups (3.4),
pair slots (cap_pg per group at hop 1, cap_p cumulative at hop 2), hash walks (uncapped, at most
HOME_MAX visits per record since a chain holds one home's records of one fold), candidates, requests, relay records
and coin retries are capped or bounded in total, so the worst-case ledger holds on every run and the E-ledger (7)
averages only group-granular executed work.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's (q2 even: both loops take two incidences per pass). With `D_l = y(Y_l)^2 - y(X)^2`, the
forward pass stores `D_l` and `acc_(l-1)` in entry l and forms `acc_l = acc_(l-1)*D_l` (`acc_(-1) = 1`); fe_inv gives
`R = 1/acc_(q2-1)`; the backward pass, for `l = q2-1..0`, sets `inv = R*acc_(l-1)`, `R <- R*D_l`, `u = inv*(y_Y - y_X)`
and `v = inv*(y_Y + y_X)`, and stores twelve **fractions**, two of each **component** `α' = Ew`, `β'' = Fw` and `γ = Gw`
(3.3), with `a = x_X`, `b = x_Y`, `E = b^4 + 2a^3 b`, `F = a^4 + 2ab^3`, `G = 3a^2 b^2` and w = u² (s = +1) or v²
(s = −1). For t = 0, 1, 2 the leaves, like the components never formed, are

    x(φ^t(X) + sY) = (α' - b) + β^t (β'' - a) + β^(2t) γ   (mod p).

**Why.** As R = 1/acc_l before step l, inv = 1/D_l, so `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)`. The curve gives
`y² = x³ + 7`, so `D = b³ - a³`, and as β³ = 1, `(b - β^t a)*Q_t = D` with `Q_t = b² + β^t ab + β^(2t) a²`. The slope of
`φ^t(X) ± Y` with `φ^t(X) = (β^t a, y_X)` is `(±y_Y - y_X)/(b - β^t a) = (±y_Y - y_X)*Q_t/D`, which is `Q_t*u` or
`-Q_t*v`; its square, less the two x-coordinates, is the leaf `Q_t²w - β^t a - b`, and expanding with `β^(3t) = 1` gives
`Q_t² = E + β^t F + β^(2t) G`, hence the display. D ≠ 0 unless `b³ = a³`, that is `Y = ±φ^t(X)` (3.4).

**Arithmetic per pair**, after W2S (`inv`, the new R) and W2L (u and v, on the shared rows of inv): u² and v² (SQ13 and
F5: a square and its fold); βu² and βv² (PAIRB: products on the tabled rows `β*2^(64k) mod p`, folded; limb pairs 0, 1, 2 cached in (x29,x30), (x27,x28), (x7,x8), only pair 3 staged); E, F and G,
each a product of a Y-entry word (b, b³ or 3b²) on the four tabled rows `c*2^(64k) mod p`, k = 0..3, of an X constant c
(2a³, 2a or a²), plus b⁴ or a⁴ for E and F (AD4), folded (FDX); then for each of E, F and G its rows (SHX, SHX, then
SHXT, which forms only the words 1..3 of its row that the products read) and four **truncated products** (TR0, TRK;
3.3) with the limbs of u², v², βu² and βv² (β pairs 0, 1, 2 read from the cache, only pair 3 reloaded), finalized, with the correction words for α' and β'', and stored. The three-pair register allocation is 66043fba's (section 9); the run checks every stored fraction bit for bit against the model of `bwd` (8). The header
H (bwd's caller's sp) holds the rows of 2a³, 2a and a², y_X, a⁴, `-T(a)` and `-T(βa)`; the rows of β lie at H − 144,
below bwd's 16-byte frame, written once per node; a Y entry (192 bytes) holds `y², y, b, b³, 3b², b⁴`; batch entry l
holds D_l (+0) and acc_(l-1) (+32) from fwd, stages R' at +224 and the pair `(βu²_3, βv²_3)` at +112, receives the
fractions (V0, V1) of α'−, β''−, γ−, γ+, β''+ and α'+ at +0, +32, +64, +96, +128 and +160 (γ's V1 as (s, hi), its low
word 8s), and holds the correction words `-T(b)` and `-T(βb)` at +192, written with the Y table and never overwritten.
The rare carries of SHX, SHXT, FDX, F5, FS and FADDW are taken by branches over fixups that give the inline chains'
results (3.4), and the only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks the entries downwards, one sign at a time. It loads the sign's six fractions (6 LDPs) and forms its three
leaves' W0 and W1 by two-word additions and subtractions (28 in all): with a, b and g the fractions of α, β' and γ
(3.3), `Sb = b0 + b1` and `Sg = g0 + g1`, they are `a0 + b0 + g0` and `a1 + b1 + g1` for t = 0, `a0 + b1 - Sg` and `a1 +
g0 - Sb` for t = 1, and `a0 + g1 - Sb` and `a1 + b0 - Sg` for t = 2, γ's V1 entering as a shifted operand (LSL #3). For
each leaf it then decides KEY's index j from the top words F0 of W0 and F1 of W1: with `σ = F0 + F1`, if F0 ≥ F1 (CMP,
B.HS) it branches to the block JH, where j = 2 if `F1 + σ` carries (CMN, B.HS to JT) and j = 1 otherwise; else j = 2 if
`F0 + σ` carries (CMN, B.HS to JT) and j = 0 otherwise, which goes on inline. The j = 0 path emits the pair (W0, W1); JH
forms `N' = ~(W0 + W1)` (ADDS, ADC, MVN) and emits (W1, N'); JT forms N' and emits (N', W0). The **emit** forms 8*d1 by
one AND (j = 0, 1) or BIC (j = 2) of the register x1 = (2^36 - 1) << 3 with the low word of the pair's first fraction
shifted right by 25 bits: with the column LC table at address 0 this is both d1's LC entry and the destination word.
LDADD advances the entry by 32 and returns the record's slot; the emit forms w3 from d2 (LSR #28, or for N' EOR of the
shifted sum with 2^36 - 1) and the entry pointer, writes the record with two STPs and tests for slot k (one TBZ of bit 13, section 4): 7
instructions. After slot k, FLUSH (8, or 7 for a dropped block) ends the block (section 4). The entry pointer x3, moved
by one post-indexed LDP per sign (−96, then −160), is entry + 64 at a sign +'s emits and entry − 96 at a sign −'s;
entries are 256 bytes at a 2^34-aligned base, so `x3 << 31` mod 2^64, which the emit adds to d2, is `2<<36 | l<<39` or
`5<<36 | (l - 1 mod 2^25)<<39`. JH and JT branch back. **Key paths.** The exact index is uniform on {0, 1, 2} for a leaf other
than O, and a leaf takes JH or JT only if its exact index is 1 or 2 or it lies in W (3.3): an event of its value, of
probability at most `2/3 + 2d_wj`. Their work (4 more per such leaf) is E-priced in 6.6 by the same history-conditional
domination as the walks; there is no key-path meter, and the worst case charges every leaf (7). Executed counts per incidence, q2 even (fixed part in parentheses): `fwd` 98 (5) and `bwd` 1285.5 (5), each
plus 3 or 4 per fixup (3.4), and `keys` 140.5 (5) for six leaves on the j = 0 path, including their six emits of 7, plus 4
per leaf on JH or JT and its flushes; fe_inv (`secp_fe.s`) is 22549. Per sign `keys` executes 6 LDPs, the sums (28)
and per leaf the index's ADD, CMP, B.HS, CMN and B.HS and the emit: 70, plus a TBZ per two incidences; a leaf on JH or
JT executes 4 more (ADDS, ADC, MVN of N', B back). Kernels may
overwrite x0-x29 and the flags and preserve sp; `bwd` and fe_inv save and restore x30, which `keys`, `fwd` and the
transport kernels do not write; `fwd` and `bwd` return their fixup word in x0 (3.4); the driver keeps no registers
across calls and sets sp to H for `bwd`. `bwd` keeps R in x25..x28 from
pass to pass; `fwd` takes two incidences per pass (SUBS of 2, B.NE), and `bwd` and `keys` end each pass with a TBZ on
bit 33 of the entry pointer: the 2^25 entries of 256 bytes fill the lower 2^33 bytes of their 2^34-aligned block, so
bit 33 is clear at every entry and set one entry below the block (both loop ends test bit 33 and the 36 emit TBZs test
bit 13: `loop_ends_on_bit_33`). In `keys` each leaf's j = 1 and j = 2 blocks (LJX) sit after the loop's RET and branch
back to the leaf's join (JPATH 4).

### 3.3 The class label (truncated fraction products, the index, the paths)

**Fractions and the label.** For v in [0, p) let `T(v) = floor(2^128 v/p)` and `comp(v) = floor(2^100 v/p) = T(v) >>
28`. For a leaf value x in [1, p) let `v_m = β^m x mod p`; the canonical member `v* = min v_m` has index j = argmin, and
the class's **label** is `(comp(v*), comp(βv* mod p))`: the member `β^m x` has index `j - m` and the same pair. As
`2^128 v/p = v/2^128 + C0*v/(2^128 p)` with the second term in [0, 2^-95), `T(v) = (v >> 128) + e0(v)`, e0 in {0, 1}.
For a field element c let `T_m(c) = T(β^m c mod p)`.

**Truncated products.** Let c be a stored E, F or G, its rows `R_k ≡ c*2^(64k)` (mod p) below 2^256 (`R_0 = c`,
`R_(k+1) = SHX(R_k)`), and M below 2^256 a multiplier with limbs `M_k`. Then `Z = Σ M_k R_k ≡ cM` (mod p) and
`Z < 2^322`, and for `z = Z mod p`, `T(z) ≡ floor(Z/2^128) + e` (mod 2^128) with e in {0, 1}, as
`2^128 Z/p - Z/2^128 = Z*C0/(2^128 p) < 2^-28`. TR0 and TRK keep the words of columns 2 and 3 of each `M_k R_k`:
`f = Σ lo(M_k r_k3) + hi(M_k r_k2) mod 2^64` and `s = Σ (lo(M_k r_k2) >> 3) + (hi(M_k r_k1) >> 3)`, below 2^64 as its
eight terms are each below 2^61, then `A = f*2^64 + 8s mod 2^128` (ADD of `s >> 61` to f). `floor(Z/2^128) - A ≡ c_low +
d` (mod 2^128): `c_low <= 7` is the carry out of columns 0 and 1 (each row's part of them is below 2^129) and `d <= 56`
the eight 3-bit parts the shifts drop. So **`T(z) - A ∈ [0, 64]`** (mod 2^128), for every such `R_k` and M. With M = w
it estimates `T_0(cw)`, the component's V0; with `M = βw` (PAIRB) it estimates `T(βcw mod p) = T_1(cw)`, its V1.

**Corrections and sums.** A sign's leaves are `L_t = α + β^t β' + β^(2t) γ` with `α = α' - b` and `β' = β'' - a` (3.1).
`bwd` adds `-T_m(b) mod 2^128` to α''s `V_m` and `-T_m(a)` to β'''s, m = 0, 1; as `T_m(c - y) ≡ T_m(c) - T_m(y) -
{0, 1}`, they estimate α's and β''s fractions. `keys` takes `U_m = V_m` (m = 0, 1) and `U_2 = -U_0 - U_1`; the values
`β^m c mod p` sum to 0, p or 2p, so `T_0 + T_1 + T_2 ≡ -k` with k in {0, 1, 2}, and U_2 errs by
`-(T_0 - U_0) - (T_1 - U_1) - k`. As `β^m L_t` is the sum of `β^m α`, `β^(m+t) β'` and `β^(m+2t) γ` in F_p,
`T(v_m) ≡ T_m(α) + T_(m+t)(β') + T_(m+2t)(γ) + κ`, κ in {0, 1, 2}; so `W_m = U_m(α) + U_(m+t)(β') + U_(m+2t)(γ) mod
2^128` errs by its three terms' errors plus κ, and `N' = ~(W_0 + W_1)` by `-(T(v_0) - W_0) - (T(v_1) - W_1) - k + 1`:

| estimate | T − estimate, mod 2^128 | from |
| --- | --- | --- |
| γ's `U_0`, `U_1`; γ's `U_2` | [0, 64]; [−130, 0] | e + c_low + d; `-(T_0 - U_0) - (T_1 - U_1) - k` |
| α's and β''s `U_0`, `U_1`; β''s `U_2` | [−1, 64]; [−130, 2] | less {0, 1}; as γ's |
| `W_0` for t = 0, 1, 2 | [−2, 194], [−132, 130], [−131, 132] | `U_0(α) + U_t(β') + U_2t(γ)`, plus κ |
| `W_1` for t = 0, 1, 2 | [−2, 194], [−131, 132], [−132, 130] | `U_1(α) + U_(1+t)(β') + U_(1+2t)(γ)`, plus κ |
| `N'` for t = 0, 1, 2 | [−389, 5], [−263, 264], [−263, 264] | the two rows above, less k, plus 1 |

Over t: `T(v_0) - W_0 ∈ [-132, 194]`, `T(v_1) - W_1 ∈ [-132, 194]` and `T(v_2) - N' ∈ [-389, 264]` (indices of U mod
3). The certificate derives this table from the three terms and ties D''s 1305 residues and the 654-unit widening to it.

**The index.** KEY decides j from the top words F0 of W_0 and F1 of W_1 (3.2): with `σ = F0 + F1 mod 2^64`, j = 2 if
`min(F0, F1) + σ` carries, else 1 if F0 ≥ F1, else 0. Let `X_m = 2^64 v_m/p`; `X_0 + X_1 + X_2` is 2^64 or 2^65. The
**window set** W is where (W1) X_0 or X_1 is below 1 or above `2^64 - 1`, (W2) X_2 is below 1 or above `2^64 - 4`, (W3)
`|X_0 - X_1| < 2`, or (W4) `X_0 - X_2` or `X_1 - X_2` lies in (−1, 4). **Exactness.** Outside (W1),
`2^64 - 1 <= T(v_m) <= 2^128 - 2^64` for m = 0, 1, so `W_m = T(v_m) - e_m` without wrap, and as
`2^128 v_m/p - T(v_m) ∈ [0, 1)`, both `a' = X_0 - F0` and `b' = X_1 - F1` lie in `[-132*2^-64, 1 + 195*2^-64)`. Outside
(W2), `σ = 2^64 - X_2 - a' - b'`. Hence F0 ≥ F1 iff `X_0 - X_1 >= a' - b'`, which lies in (−2, 2); `F1 + σ` carries iff
`X_1 - X_2 >= a' + 2b'`, and `F0 + σ` iff `X_0 - X_2 >= 2a' + b'`, both in (−1, 4). Outside (W3) and (W4) each test
has the sign of the exact difference: F0 ≥ F1 means v_0 > v_1, and then the carry of `F1 + σ` decides `v_2 < v_1`;
otherwise v_0 < v_1, and the carry of `F0 + σ` decides `v_2 < v_0`. So outside W, j = argmin v_m. **A wrong index.**
Where j is wrong the leaf emits another member's pair: a missed collision charged as the row `d_wj` of 6.4, or, if it
equals another class's label, a false candidate, which the mass bound below covers. **Rarity.** `X_m < 1` means
`v_m < p/2^64 < 2^192`, and `|X_i - X_j| < c` means `(β^i - β^j)x mod p` is below `c*2^192` or above `p - c*2^192`. As
these maps are bijections, (W1) holds for at most `4*2^192` leaf values, (W2) for `5*2^192`, (W3) for `4*2^192` and (W4)
for `10*2^192`; a uniform point's x takes each value with probability at most `2/(r-1)`, so
`d_wj <= 46*2^192/(r-1) < 2^-58.47` per leaf.

**The paths.** A leaf of index j emits `(W_0, W_1)` for j = 0, `(W_1, N')` for j = 1 and `(N', W_0)` for j = 2: of each
fraction, bits 28..63 are its digit (d1, d2) and bits 64..127 its top word (K2, K3). With the index exact, `v_j = v*`
and `v_(j+1) = βv*` (indices mod 3), so the pair is the class's label unless a computed fraction's error carries it
across a block of 2^28: the **label-error set** D' is where `T(v_0) mod 2^28` or `T(v_1) mod 2^28` lies in
`[0, 194) ∪ [2^28 - 132, 2^28)`, or `T(v_2) mod 2^28` in `[0, 264) ∪ [2^28 - 389, 2^28)` (a wrap mod 2^128 is at these
residues too). Outside W and D' a leaf of a batch without a deviation (3.4) emits its class's label, whatever the
representatives. As each value of T is taken by at most 2^128 values of v and `x -> v_m` is a bijection, D' has at
most `1305*2^228` values of x, each of probability at most `2/(r-1)`: the row `d_lab < 2^-16.65` of 6.4, a missed
collision like `d_wj`; no leaf is dropped.

**Mass bound.** secp256k1's Frobenius is `π = A + Bω` with `B = 0xe4437ed6010e88286f547fa90abfe4c3` and
`A = 0x114ca50f7a8e2f3f657c1108d9d44cfd9`: `A² - AB + B² = p`, `B*β ≡ A - B` and `(A - B)*β ≡ -A (mod p)` (checked by
the run and the certificate). So `Λ = {(v, u) ∈ Z²: u ≡ βv (mod p)}` has the basis `(B, A - B), (A - B, -A)` of
determinant -p, and as `A > B > 0` its half-open parallelogram P spans [0, A) in v and (-A, A - B) in u. For each point
λ of Λ in a box `[v0, v0 + a) × [u0, u0 + b)` of integer sides a and b, `λ + P` lies in `[v0, v0 + a - 1 + A) × (u0 - A,
u0 + b - 1 + A - B)`; these translates are disjoint, as P tiles the plane under Λ, and each has area p, so the box holds
at most `(a + A)(b + 2A - B)/p` points of Λ. Fixing f of a label's 200 bits leaves `2^(200-f)` pairs of component
values; a component value is an interval of at most 2^156 values of v (of u = βv mod p for the second), widened at each
end by its fraction's error range, at most 654 units of `p/2^128` in all (N''s [−389, 264]), for a computed fraction
(one that wraps mod 2^128 may be two end pieces of [0, p), which count as one interval as `pZ² ⊂ Λ`). So on each path at
most `2^(200-f)(2^156 + w + A)(2^156 + w + 2A - B)/p < (1 + 2^-17)2^-f p` values of `v_j`, and of x, give those bits,
with w the widening (the certificate checks the factor), whatever index KEY returns. With three paths and two points per
x, the emitted pair takes given values on f of its bits with probability at most `KB(f) = 6*2^(256-f)(1 + 2^-17)/(r-1)`.

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts
`B_(k+1) = SHX(B_k)`; the rows of an X constant and of β are tabled exactly. RF and RA sum `Z = Σ a_k*B_k < 2^322`
exactly (each high word is at most 2^64 - 2), AD4 adds b⁴ or a⁴ (below p) for E and F, `Z' < 2^322 + 2^256`, and a fold
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ13 forms the square of a weak a as
`V = L + H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), and F5 folds V the same way. Only these
steps drop a carry: SHX's last ADCS (SHXT is SHX without its low word's MUL) and the final ADC of FDX and of F5. Without
a dropped carry each returns a value below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; of
  `Z' < 2^322 + 2^256` (E, F, G), in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`.

A fold returns a value `>= p` only when the true value is below C0: a multiplier stored so is still its value mod p,
and 3.3's bound holds for every multiplier below 2^256; an E, F or G stored so has its low 192 bits above `2^192 -
2^33`, in the bad set of its first row shift, which the table below charges. Every other operation is exact on its
inputs: RF and RA; AD4; FS and FADDW, whose operands (`y_Y²` and `y_X²` in `fwd`, `y_Y` and `y_X` in `bwd`) are
canonical table values; the truncated products (f wraps mod 2^64 by design, s stays below 2^64); the corrections, sums
and path additions (3.3); KEY's integer arithmetic (3.3 covers its index); the emit; and F's `fe_mul_f` and `fe_sqr_f`,
which fold the 512-bit product, fold its fifth limb again and add C0 exactly when that last sum carries past 2^256, so
fe_inv (255 squarings, 15 multiplications) is exact for every input below 2^256. Generator construction, the table
entries, the header rows and values, the correction words, the resolver and `[x]P = Q` run on `fe_mul_f`, `fe_sqr_f` and
fe_inv, each coordinate brought below p (an ADDS chain of C0 detects `t >= p`, CSELs keep `t - p`) before it is
subtracted, compared or stored, with exact arithmetic mod r for the solve. So they have no sites; the tables hold exact
canonical values (no point has x = 0, as 7 is not a square mod p, so a, b ≠ 0).

**Distribution.** Given X ≠ O and every other generator, `Y_l` is uniform on the r - 1 points of G \ {O}. A value of
x_Y comes from at most 2 points, and of y_Y from at most 3. So `D_l = x_Y³ - x_X³`, `acc_l = acc_(l-1)*D_l`,
`R_l = 1/acc_l` and `inv_l = 1/D_l` take each value with probability at most `6/(r-1)` (the rest are bijections while
no `D_m = 0` for m ≤ l, the degenerate case below); `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)` (as
`D = (y_Y - y_X)(y_Y + y_X)`) at most `3/(r-1)`; u² and v², and βu² and βv² (bijections of them), at most `6/(r-1)`;
and E, F and G, of degree 4, 3 and 2 in x_Y, at most `8/(r-1)`. A row `B_k ≡ v*2^(64k)` determines its operand v, so
at most 2^97 values of v make a given row shift carry. Each site's fresh generator is `Y_l`, except `Y_(l-1)` for fwd's
row shifts of `acc_(l-1)` and bwd's fold of `R_(l-1)`; the fixup words are those of **How often** below:

| site, at incidence l | runs per batch | bad values | fixup words | points per value |
| --- | --- | --- | --- | ---: |
| fwd: FS forming `D_l` | q2 | none | `C0*2^192` | 6 |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | 2^97 | `C0*2^192` | 6 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | 2^66*C0 | `F_w` | 6 |
| bwd: row shifts of `R_l` | 3q2 | 2^97 | `C0*2^192` | 6 |
| bwd: fold of `inv_l = R_l*acc_(l-1)`, l ≥ 1 | q2 - 1 | 2^66*C0 | `F_w` | 6 |
| bwd: fold of `R_(l-1) = R_l*D_l`, l ≥ 1 | q2 - 1 | 2^66*C0 | `F_w` | 6 |
| bwd: FS of `y_Y - y_X`, and FADDW of `y_Y + y_X` | q2 each | none | `C0*2^192` | 3 |
| bwd: row shifts of `inv_l` | 3q2 | 2^97 | `C0*2^192` | 6 |
| bwd: folds of u and v | 2q2 | 2^66*C0 | `F_w` | 3 |
| bwd: square folds of u² and v² | 2q2 | C0² + 2C0 | 2^193 | 6 |
| bwd: folds of βu² and βv² | 2q2 | 2^66*C0 | `F_w` | 6 |
| bwd: folds of E, F, G | 3q2 | (2^66 + 1)*C0 | `F_w` | 8 |
| bwd: row shifts of E, F, G (the third by SHXT) | 9q2 | 2^97 | `C0*2^192` | 8 |

A run fails with probability at most its bad count times its points per value over r - 1. At l = 0 three row shifts
and three folds can neither fail nor take a fixup, and are not counted:
- fwd's rows of `acc_(-1) = 1` are 1, 2^64 and 2^128 (low words at most 1), and its fold of `D_0*1` has `Z < 2^256`;
- bwd's fold of `inv_0 = R_0*1` has `Z = R_0 < 2^256`;
- bwd's last update `R_0*D_0 ≡ 1` folds to 1 or p + 1, as nothing else congruent to 1 lies below `2^256 + 2^66*C0`;
  neither is at least 2^128 with its low 128 bits below `2^66*C0`. It is left in x25..x28 after the last pass, unread.

The degenerate cases are X = O (probability 1/r) and `D_l = 0`, that is `Y_l = ±φ^t(X)` for some t (at most `6/(r-1)`
per incidence), which zeroes acc and spoils the batch. The first site of a batch whose predicate holds has exact
operands, so a batch has a defect with probability at most `d_batch`: the table's sum plus `6q2/(r-1) + 1/r` (`SITES`
in `certificate.py`; 6.4).

**Rare carries by branch.** In SHX and SHXT the carry out of the first ADDS, in FDX and F5 the carry into limb 2, in FS
the borrow of its C0 subtraction and in FADDW the carry of its C0 addition are rare. Each macro branches over its
**fixup** (B.LO; B.HS for FS) unless that carry or borrow occurred. The fixup is the rest of the inline chain, two
instructions (ADCS and ADCS or ADC) or three for FS and FADDW, and `STADD x4, [sp]`. Where the branch is taken the
skipped instructions would add or subtract 0, so on every input every register equals the inline chain's, and so does
every deviation of the table above; only the flags differ, and every path from the join sets them before reading them
or returns (`fixups.flags_dead_at_every_join`; kernels may leave the flags changed). A fixup adds 3 instructions (SHX,
SHXT, FDX, F5) or 4 (FS, FADDW) and C0 to the word at [sp], which opens a 16-byte frame pushed as zero by `fwd` (STR)
and by `bwd` beside its saved x30 (STP); each returns the word in x0, at most `33q2*C0 < 2^63` per call. The driver adds
x0/C0 to the node's fixup count N and starts a batch only while `N <= 2^45 - 33q2` (33 sites per incidence: 18 SHX, 3
of them SHXT, 10 FDX, 2 F5, 2 FS and 1 FADDW). So N is at most 2^45 on every run, and section 7 charges `4*2^45` per
node.

**How often.** A fixup needs a word in a small set:
- A row shift's carry needs the row it shifts to have its low word above `2^64 - C0`: at most `C0*2^192` rows.
- A fold's carry into limb 2 leaves its output with low 128 bits below `F*C0`, `F = floor(Z/2^256) <= 2^66 + 1` (for
  E, F, G with their constant), so at most `F_w = (2^66 + 1)*2^128*C0` words; a square fold's leaves its low 128
  bits below 2^65, so at most 2^193 words.
- FS's borrow and FADDW's carry leave an output whose low word is in a set of at most `C0*2^192` words.

While no site of the batch up to and including this one has deviated, the row or output is a word below 2^256
congruent to the site's operand value (times `2^(64k)` for the k-th row), which determines that value. So, given X
and every other generator, a site takes its fixup with probability at most its set's size times the points of `Y_l`
per operand value over r - 1, with the table's fresh generator and multiplicities. Summed over a batch (`FIXUPS` and
`E_FIX` in `certificate.py`), the expected number of fixups a batch takes before its first deviation, and so of a
batch without a defect, is at most `E_FIX`, and over a node's batches at most `q1*E_FIX < 2^28.68`.

**The meter row.** A node stops early only if N exceeds `2^45 - 33q2` at the start of some batch. On the event that
none of its batches has a defect, N is that clean count, so by Markov's inequality a node stops with probability at
most `q1*E_FIX/(2^45 - 33q2) + q1*d_batch < 2^-16.32` (`d_f`, a node-local row of 6.4).


## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (7), and per block of k = 256 records its flush (8) and `dlvk` (1026) | 11.0390625 | 11.046875 |
| hop 2: `h2s` (7; loop 0.25), and per block its flush (8; 9 if dropped) and `dlvk` (1026) | 11.2890625 | 11.296875 |
| hash join `hchg` fixed (7; loop 0.25) | 7.25 | 7.25 |
| hash walks (expected; 6 per visit, one B after an even walk: 6.5 per expected visit) | ~3.2 | E-rows |
| K2 verifications and listings (expected; global pair counting) | ~0.3 | E-rows |
| **fixed total** | **29.578125** | **29.59375** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlvk`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. A destination's full blocks number at most its records over k, and a
dropped block costs 7 or 9 and sends nothing, so 1036 per k records bounds each hop's blocks; final blocks are `fin`'s
rows (section 7). `h2s` and `hchg` take a slot's records in groups
of 1, 2 and 4, as the low bits of its count say, then 8 at a time; the 8-record loop's CMP and B.LO, 2 per pass, are
charged as 0.25 per record (a slot of n records runs floor(n/8) passes), and the tests and set-up, 11 and 13 per
slot (`hchg` loads the slot's read cursor and stores E2 back), are per-node rows. **Addresses.** One instruction
**Blocks and cursors.** The emit and `h2s` take a record's slot by one LDADD of 32 on its LC entry, which
holds the absolute slot address (block base with low 13 bits zero, offset 32n below 2^14), so one TBZ of bit 13 of
the old value finds slot k (32k = 2^13, checked). The flush resets the LC entry to slot 1, forms the RC entry (ORR)
and advances the RC cursor by 32k with LDADD, returning the old cursor c; the block is sent by the STP of the header
(destination word, c) when c is below the sender's slot end, which is the same at every destination and held in a
register, and dropped otherwise, and `h2s` stores the end back so a hop-2 cursor stays at most its slot end (a hop-1
cursor can pass its end only by the node's own leaves, 32 bytes each, staying below 2^63.2, checked). `fin` sends the
L = (LC - slot 1)/32 records left in a block by the header (8d, c) if L > 0 and c is below the end, and writes
`min(c + 32L, end)`. **Keep lemma.** Cursors start at the slot base and move by 32k per block with `cap_p` a multiple
of k, so c and the end are congruent mod 32k and c below the end gives `c + 32k <= end`: a block is kept or dropped
whole, and the records kept, their slot order and the final cursor equal those of sending each record alone at its own
cursor (the one-record rule). A final block's unused slots land in `[c + 32L, c + 32k)`, at or past the final cursor,
where every later read stops. At hop 2 the LC and RC entries persist across groups, so the next group keeps filling
the block from slot L + 1 and rewrites `[c, c + 32L)` byte-identically before its own records at
`[c + 32L, c + 32k)`; the lemma holds at `cap_pg` too (a multiple of k, checked). **Hash table and join.** Each node
owns 2^57 head words at 2^62, zeroed once by `zero4096` (2050 per 4096 words), and a packed next-array at 2^61 with
one 8-byte next word per R2 slot position (section 5). Bucket `K2 mod 2^57`: same K2 always lands in the same bucket,
so same-key records meet with no ordering step. `hchg` reads each slot from its read cursor to its E2 cursor and
stores E2 back (the per-slot 13), so the per-group ranges partition each slot and every record is read exactly once
however the slots' previous ends differ; cross-group state otherwise travels only through the cumulative table. Per
record (7): LDP; BFI of the bucket into x10 (the head's address); LDR of the head; ADD of C to A >> 2 (the new
record's next word); STR of the head into it before the walk, which reads only priors' next words; CBZ to the final
STR on an empty bucket; that STR of A into the head. The walk is unrolled by two, the prior alternating between x11
and x21: a visit loads the prior's K2 and, by ADD of C to its address >> 2, its next word, then CMP, B.EQ and CBZ or
CBNZ (6), and an even chain end takes one B to the final STR (6.5 per visit bounds a walk). A K2 hit verifies d4 (5)
and rejoins the walk on a mismatch; a full hit lists (4, or 3 from the second half) unless the list is full (2); the
matched visit skips its CBZ or CBNZ (the run counts 5, the ledger 6). Uncapped, so every new record links to its most
recent same-key predecessor on every run. The next-array fits every receive pattern up to the slot caps by
construction, next words need no zeroing (written before any read), records stay intact in R2 for phase 4, and a
receiving node executes nothing on arrival: it reads the records with the LDPs of `h2s` or `hchg`. (`bkt`, `pass4`
and `scan_ranges` stay for part B's base-path evidence; the grouped run never executes them.)

**Chains replace the sorted run.** A home's records share d1 and d2; records with equal (K2, d4) — the 182 label bits
— land in one bucket and form one arrival-order chain, each member linked to its most recent same-key predecessor.
The records of one class form one connected chain unless an **intruder**, a record with the same 142 bits that does not
carry the class's label, arrives between two members (6.4, 6.6); every listed pair agrees on d1, d2, K2 and d4, so
besides a class's edges it can list pairs of different classes, which the resolver rejects (no unit relates their
points) and 6.4 charges to the candidate, request and relay caps through KB(182). Writing the members `u_i*Z` with Q
coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`, an equivalence relation, so a chain with a useful pair
has a useful chain edge.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums the regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.3297 used) | code, stack, the five barrier words (one reservation granule each), header H with the rows of β below it, X table of the running group (512 bytes per X), X affine store (64 bytes per X, all groups), Y table (192 bytes per Y), scalars, batch entries (2^33 bytes at a multiple of 2^34), the ports' message blocks (at 2^48 - 64 and 2^48 - 128), their LC tables at 0 and 2^39 and RC tables at 2^43 and 2^43 + 2^39 (2^39 bytes each), the hash-read cursors at 2^48 (2^39 bytes), the record blocks at 2^50 and 2^51 (2^50 bytes each), E1 and E2 (2^41 each) and the phase-4 buffers and tables (`FIXED` in `certificate.py`, whose d3/d4 tables stay reserved but unused by the grouped run); all but the aligned tables lie between 2^44 and 2^48 | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | per-group hop-1 slots of cap_pg at 2^52 (`2^59.9986`), the hash next-array at 2^61 (`2^60.9986`, R2 address A's word at (A>>2) + C), the hash heads at 2^62 (`2^60`: a head's address is 2^62 with the bucket in bits 3..59); disjoint, inside R1 (checked), together 2^61.9989 | per group |
| R2 after R1, ending at 2^63.998943 | cumulative hop-2 pair slots (phase 4 reads candidates' records at their addresses) | groups 1-32 |

So every node stays below 2^64 bytes, and the peak is `2^72*(fixed + |R1| + |R2|)`, at most 2^135.9989
bytes (log2 ≈ 135.9988120). A remote item is reached only by a message.

## 6. Success (no heuristic)

In `success.md`: the model (6.1), Chen–Stein on cross-node pairs (6.2), usefulness (6.3), the defect rows (6.4), the
bound (6.5) and the grouped hash-chain abort (6.6: stopping-time bound, chain detection E1-E7, uncapped walks with
history-conditional walk means and a defective-run allowance).

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose itemization
`certificate.py` states in its comments and checks: generator construction 645,410 per point (645,414 for a grouped
X, whose four affine-store STPs the X row charges; a Y's table entry and no mixed addition make it 642,930,
`GEN_MULY`), node set-up 2^20, a message outside the hops 40 per hop, the candidate 1,557,745, the phase-4 count
exchange 2^41 per node and sub-phase, the barrier 131 per node (2^10 per node charged), phase control 2^16 per node,
`zero4096` 2050 per 4096 words, `init3` 4 per slot, the rows of β 72, batch set-up 160 (145 itemized), the inversion
call 16, a request served 2^10 and the output 2^20. Coins cost 64 per scalar: two 4-word `MRS Xt, RNDR` draws
(FEAT_RNG) with a B.EQ per draw and a range test; RNDR returns 64 fresh uniform bits, the algorithm coins of the
rules' probability space, and never reports failure, so the charged branch is never taken and coins add no failure
row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, the barriers and the output. A record row uses
its group's records; E-averaged rows are sums over groups weighted by PHI. With q1g the grouped q1 (145307571):

| row (P: preprocessing; E: E-averaged; W: worst-case) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1g + q2 + 1)*64 (P) | 105.4143 |
| Y generators with their table entries: q2*642930 (P) | 116.2943 |
| R_g generator | 91.2999 |
| rows of β: 72 (P) | 78.1699 |
| row LC and RC tables: (4*2^36 + 1) + (3*2^36 + 1) (P) | 110.8074 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| hash-head zeroing: zero4096 2^57/4096*2050 + 1 (P) | 128.0014 |
| hash-read-cursor init: init3 4*2^36 + 1 (P) | 110.0000 |
| leaf work, E-averaged: sum_g PHI[g] (X gens, fwd, inv, bwd, keys w/o emits, batch set-up) (E) | 134.4279 |
| arithmetic fixups: 2^45*4 (W) | 119.0000 |
| j = 1 and j = 2 key paths, E-averaged: sum_g PHI[g] 4 s_g PJ + DEFJ (E; every leaf, 128.6995, is the worst case) | 127.8946 |
| route work, E-averaged: sum_g PHI[g] (emits, blocks, h2s, fin, column tables) (E) | 130.9959 |
| group work, E-averaged: sum_g PHI[g] (hchg fixed, main-pair walks, stop test) (E) | 129.9257 |
| same-X walks: global pair counting (E) | 29.5699 |
| K2 verifications: global pair counting (E) | 124.8909 |
| listings: global pair counting (E) | 75.9840 |
| defective-run rows: DEFV + DEFK2 + DEFL + DEFJ (E) | 63.3670 |
| candidates: 2^28*1491011 (W) | 120.5079 |
| requests served: 2^29*2^10 (W) | 111.0000 |
| phase-4 count exchanges: 114*4*2^41 (W) | 121.8329 |
| barriers: 2048*2^82 (W) | 93.0000 |
| phase control: 2^16 (W) | 88.0000 |
| output + stop broadcast (W) | 85.1548 |
| **expected total E** | **134.643057** |

`E` is the exact `E_total` of `certificate.py` (expected instructions), at most 246.1820 per leaf (the run's
`per_leaf`), and `W_E = E/1000` with log2 W_E ≈ 124.6772724914. **The claim, 124.677272491447240, is the least multiple
of 10^-15 with 2^claim >= W_E** (the base claim stays on 10^-12); the certificate encloses 2^c by an exact series for
e^(f*ln 2) with interval bounds on ln 2, both ways (2^(claim - 10^-15) < W_E), from an 80-digit decimal start.
Rows marked P are preprocessing: their sum over 1000 has log2 ≈ 118.0360707 and is declared as 2^118.0361,
the least multiple of 10^-4 at or above it (checked). The worst case is at most 2^188.698089, the least multiple of
10^-6 at or above log2 R_worst, both ways certified (2^188.698088 < R_WORST <= 2^188.698089): all groups run, with
uncapped walks to HOME_MAX visits and verification and listing on every record. It is disclosed, not scored.

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "seconds": 12}` (`seconds` 5 to 15). Budget-limited checks run
The declared inputs are seed 1 and 12 s (`seconds` 5 to 15). Budget-limited checks run their guaranteed
minimums and then their share of the remaining budget. **Part A** runs `leaf.s` and fe_inv on production batches of
2, 4, 6 and 10 incidences and checks every stored fraction against a bit-exact model of `bwd` and its bound, the
truncated-product terms, 3.1's identity against exact affine addition, the records and slots against a bit-exact model
of `keys` and the emit, the class rule and class labels outside W and D', the fwd/bwd formula counts, planted fixups,
a copy of `leaf.s` with every fixup taken inline, and each site type of 3.4 against its bit-level model
(`fraction_mismatches`, `fractions_outside_bounds`, `tfp_terms`, `identity_mismatches`, `wrong_index`,
`keys_not_canonical`, `class_rule`, `class_labels`, `window_batch`, `boundary_batch`, `fwd_bwd_formulas_hold`,
`fixup_batch`, `forced_fixups`, `site_checks`, `sites_match_certificate`, `fixups`,
`fixups.flags_dead_at_every_join`). **Part B** (`machine.py`) executes `leaf.s` and `route.s` on a 4 x 4 machine,
delivering every message by executing `dlvk` or `dlv`, and checks the block, cursor, bucket, pass and scan paths and
the grouped hash-chain run (stop, full and full-list cases) against the same exact models (`block_paths`,
`unrolled_paths`, `two_group_finals`, `grouped_paths`, `old_x_served`). **Part D** runs `certificate.py`;
`D_certificate.checks` lists every check with its value and `counts_match_certificate` compares its constants with
parts A and B.

The 4 × 4 run is evidence at its own scope only.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1 (links keep no state beyond the message in flight; the sender's
  instructions, each delivery as `dlv` or `dlvk` and the in-flight counts are charged); a forwarding network would
  add work per extra hop and is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's fixup-meter test, stop test and phase-4 exchanges), the in-flight counts, and the production attack,
  machine and links. The toy stop fires on a listed candidate; production stops only on a resolver-confirmed x, a
  stronger rule, so the E-bound (which charges confirmed stopping) is conservative over the toy.
- **Provenance**, not a premise: the base package 26e42f93 by kaankolcu and its credits stand (d91b6dc5,
  bradcush's 25e43d87 and 4598c05a, e59be87d, bbe1cd85, 10c5d938, 02feb284, 80e9e8af, 1efcef88, 892a3df8,
  b00400fa, 1fbd26e8, 7e8d13c7, 807e4bf5, a9297fb0, 66043fba, b63af001, 7c8264bf, 03438c81, 266e82b9, 46385d36,
  aburan28's 79854ee2, pepedesigner's b505cb4e and p256 1772672b); Galbraith, Wang and Zhang for X ± Y from one
  inverse; Arratia, Goldstein and Gordon for Chen–Stein. New vs 26e42f93: the unequal abort ladder (2.3, 6.6,
  `changes.md`). The run and `certificate.py` re-verify the code.
