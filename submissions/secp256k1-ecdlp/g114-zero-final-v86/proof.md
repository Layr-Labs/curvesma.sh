# Stored x-only leaf birthday search on secp256k1: six leaves per inverse, grouped hash-chain abort (124.797251154654)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; its expected total work is at most
`time_log2 = 124.797251154654 = 124797251154654/10^12` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit), and on every run at most R_worst = 2^188.604980 (disclosed: all 114 groups run, every hash walk
visits its whole chain of at most HOME_MAX links). Peak memory is at most 2^135.9989 bytes; preprocessing, at most
2^118.0361 units, is included; there is no advice.
`certificate.py` checks these declared figures (`DECLARED`): the work, memory and preprocessing
figures are the least multiples of 10^-12, 10^-4 and 10^-4 at or above their exact values, the work per leaf 267.5303
the least multiple of 10^-4 at or above E/S, and 0.390000001 is below the success bound.

**Reviewer's guide.** **New against d91b6dc5** (124.988679192098, same line): 114 X-groups with early abort at the
first resolving group (2.3, 6.6, 7); per-group fresh hop-1 slots (cap_p/8), cumulative hop-2 slots with per-slot
hash-read cursors, and cumulative hash chains over 2^57 zeroed head words keyed by K2's low 57 bits with uncapped
walks (4); new kernels hchg, zero4096 (4; leaf.s untouched); E-work ledger with prefix non-stopping bounds PHI (7,
`E_TOTAL`); grouped q1 = 145307752 re-searched least. **Repair note:** hash chains (not v1's unmerged runs) link
each record to its most recent same-key predecessor with no record order; v4 added read cursors, domination and
uncapped walks; v5 links chains through a packed next-array (one 8-byte next word per R2 slot position, every
receive pattern covered by construction) and allows defective-run hash work globally (DEFV/DEFK2/DEFL) with
validity-split domination. 3.1, 3.3,
the resolver and the success rows are d91b6dc5's, re-evaluated at prefixes where stated. The run prints every
certificate check with its value (`D_certificate.checks`). The bound rests on:
- **Early abort with E-priced group costs (2.3, 6.6, 7).** Group g's leaf and route work is exact given g and runs
  iff T > g - 1; its main-pair hash-walk work has mean bounded uniformly given the history H_g = (Q, all node b values, X_<g, Y), over
  which the stopping indicator is fixed, so E[work] = E_pre + sum_g PHI[g] W_g + same-X/K2/list/defect rows +
  worst-case rest by history-conditional domination with a validity split (no independence is claimed); PHI[g] is
  the Chen-Stein prefix upper tail plus prefix defect rows (certificate `PHI`, nonincreasing, PHI[0] = 1). Caps,
  preprocessing, resolution, barriers and control are worst-case. Stopping needs a resolver-confirmed x ([x]P = Q
  checked), which latches; E1-E7 give the same success event as the ungrouped run.
- **Order-free grouped detection by hash chains (4).** `zero4096` zeroes the 2^57 head words once; `hchg` reads each
  slot from its read cursor to E2 and advances the cursor, then walks the K2 chain (R2 addresses, packed next-array)
  to the first full match or the chain end (uncapped), verifies d4, lists unless full, and prepends. Same-key records
  form arrival-order chains, so every colliding pair with no intruder between its records is listed on every run.
- **Six leaves per inverse, never formed (3.1).** One batched inverse of `D = y_Y^2 - y_X^2 = x_Y^3 - x_X^3` gives
  `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)`. With `a = x_X`, `b = x_Y`, `E = b^4 + 2a^3 b`, `F = a^4 + 2ab^3`, `G =
  3a^2 b^2` and w = u² or v², the leaves `x(φ^t(X) ± Y) = (Ew - b) + β^t(Fw - a) + β^(2t)Gw` are a Fourier transform
  over the cube roots of unity, and `bwd` stores only the six **components** Ew, Fw and Gw: 1206.5 instructions per
  pair.
- **A two-fraction class label from the components (2.2, 3.3).** The label is the top 100 bits of `v*/p` and
  `(βv* mod p)/p`, `v* = min(x, βx, β²x)`. Fractions mod 1 add and β rotates a component's three fractions, so a leaf's
  fractions are sums of its components' (with correction words for the linear terms), within [−73, 73], [−73, 206] and
  [−212, 146] at 128 bits; KEY's index j on their top words picks the pair (j, j + 1). A leaf in the window set W or the
  label-error set D' may carry another label: missed-collision rows below 2^-58.47 and 2^-17.38 per leaf; no leaf is
  dropped. A lattice count from secp256k1's Frobenius `π = A + Bω` bounds a box of label values at `1 + 2^-18` times its
  share.
- **Rare carries by branch, metered (3.4).** SHX, FDX, F5, FS and FADDW branch over the rest of a carry or borrow chain
  unless it occurred: results are bit-identical to the inline chain's, and no path reads the flags before setting them
  (checked). A fixup costs 3 or 4 and adds C0 to a stack word; a node's mean fixup count is below 2^29.15 and its cap
  `2^45 - 37q2`, so the ledger charges `4*2^45` per node and Markov's inequality the cap's failure (below 2^-15.85 per
  leaf).
- **No uncharged transfer or bookkeeping (2.1, 4).** Per-record work is executed code: the emit (7) at hop 1, `h2s`
  (7.25) at hop 2, grouping 11.25 fixed per record (`hchg`: LDP, AND-mask, head address/load/save, CBZ, head-address
  recompute, next-address pair, next store, head store, loop eighth) plus expected chain walks bounded by
  history-conditional KB(129) pair counting (same-X pairs carved out globally, defective inserts allowed globally),
  K2 verifications and listings charged globally, and one head-table zeroing plus one read-cursor init charged in
  preprocessing (next words need no zeroing: each is written at its record's insert, before any read). Records fill
  per-destination record blocks of
  k = 256 in the sender's port; a full block is one message, charged its delivery as executed `route.s` `dlvk` (1026)
  and its port's in-flight count (2), and kept or dropped whole, as if each record were sent alone (the keep lemma, 4,
  also at cap_pg). Links keep no flags, rings or counters.
- **Success, work, machine.** Chen–Stein on cross-node leaf pairs, usefulness, the defect rows and the final inequality
  (6, in `success.md`), evaluated with exact Fractions (grouped bound 0.3900000029, the run's `success_lower_bound`);
  the E-ledger recomputed exactly, each kernel call's count checked against its constants (7,
  `counts_match_certificate`); loads, stores, addressing, control and transfers count, instructions not cycles, and
  memory has no scalar contribution (2.1, 5; scope in 9).
- **Not relied on:** any heuristic, any prior acceptance or score, or extrapolation from the 4 × 4 run.

## 1. Overview

Each node computes the **leaves** `x(φ^t(X) ± Y)`, t = 0, 1, 2, for independent points X and Y, six per pair from one
batched inverse: Galbraith, Wang and Zhang (ePrint 2015/605) share one inverse between `X ± Y`, and φ(x, y) = (βx, y)
extends this to the three twists of X. The six leaves are held as their Fourier components over the cube roots of
unity and labelled without being formed: by two fractions of its class's canonical x, each leaf goes in two
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
  blocks**, one per destination index d at `PB + 2^14 d` (PB = 2^50 for the column port, 2^51 for the row port), each
  hold a header (destination word, address) at +0 and k = 256 record slots, slot n at +32n; the STP of the header
  completes a k-record message. The link delivers a message: it writes its data (four words, or the k records in
  slot order) at that address of the destination's memory and reports the write to its port (the in-flight count's
  decrement below). It keeps no flags, rings or counters; its only state is the message it carries. A
  store that completes a message stalls until the link has taken it, which for a block means read it, and a stalled
  instruction counts once: a port holds at most one message, and no store refills a block the link has not read.
- **Delivery is charged.** A delivery is charged as the work of `route.s` code at each hop: `dlv`, 6
  instruction-equivalents per one-record message (three LDPs read its 48 bytes from the link's block, two STPs write the
  four data words at the address, RET), or `dlvk`, `4k + 2 = 1026` per k-record message (an LDP of the header, per record
  two LDPs from its slot and two post-indexed STPs at the address, RET). Per record `dlvk` does `dlv`'s data work.
  Performed by a processor for each kind of message, a delivery is this code run as a loop, with B in place of RET
  (`dlvk` returns x0 to its start): its first LDP waits until x0, the link's window on the message, holds a complete
  one, as a sender's store waits for the link, and a stalled instruction counts once. So 6 or 1026 per message is all of
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

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145307752` and `q2
= 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the x-coordinates
of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and `S = 2^72*n0 = 2^126.699499`
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
slot i of R1, at `R1 + 32*cap_pg*i`, at each node of its column, and a node in column j writes at hop 2 only into slot j
of R2, whose physical stride is `32*(cap_p+k)` with logical cap `cap_p` (guard-repair.md). Each port has two tables of 2^36 entries, one per destination d: the **LC** entry is the
next slot of d's record block (`init3` sets it to slot 1) and the **RC** entry the sender's cursor in its pair slot at d
(`tinit` sets it to the slot's base); the slot ends at the same address at every destination, so the sender keeps that
end in a register. The column port's LC table is at address 0 and the row port's at 2^39, so d's entry, 8d or 2^39 + 8d,
is also a destination word for d; the RC tables lie 2^43 above them (ORR of 2^43).

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators with their table
  entries, keys' block constants, the row port's LC/RC tables, the hash-head zeroing (`zero4096` over 2^57 words) and
  the hash-read cursors (`init3` of each slot's base into its 8-byte entry at 2^48). Then Q is broadcast in 2 hops
  and R_g is generated.
1-4. **Groups** g = 1..114, each over a contiguous range of X's (all nodes in lockstep; the group loop tests the
  latched stop flag first and halts the run if set). Per group: generate its X's; run the batch kernels (section 3)
  for them; hop 1 into fresh per-group slots of cap_pg = cap_p/8 records (column LC/RC reset by `init3`/`tinit`);
  hop 2 into the cumulative R2 slots of cap_p (row cursors persist); hash-chain detection over the records appended
  since the previous group call (`hchg` reads each slot from its read cursor to E2 and stores E2 back, section 4);
  the four resolution sub-phases for the group's new candidates; then the stop test (below).
  Barriers separate every phase of every group.

1. **Hop 1 (per group):** the batch kernels run for each X_k of the group; the emit writes each leaf's record into
   d1's block in the column port, and after slot k the block goes to `(d1, j)` at its RC cursor, or is dropped if
   that cursor is at or past the slot end (section 4). Then `fin` sends each row d's final block by the same test
   and writes its final cursor, clamped at the slot end (drops advance cursors too), into entry i of d's table E1.
2. **Hop 2 (per group):** `h2s` reads slot i of R1 up to `E1[i]`, for every i, and writes each record through the row
   port's blocks to `(own row, d2)` by the same rule into the cumulative R2; `fin` sends the final blocks and writes
   each home's final cursor into its `E2[j]`.
3. **Grouping (per group):** `hchg` reads each slot of R2 from its read cursor to E2, stores E2 back, and
   hash-joins the new records into the cumulative table: bucket K2 mod 2^57, R2-address chains via the packed
   next-array (one 8-byte word per slot position, at (A>>2)+C), uncapped walks, d4 verified on a K2 hit, the pair
   listed unless the cumulative list (up to 2^28) is full, then prepended. Per-group [cursor, E2) ranges partition
   each slot (every record read and inserted once); the next-array covers every receive pattern up to the slot caps
   by construction, with no bump and no capacity check. No buckets, passes, snapshots or scans (section 4, 6.6).
4. **Resolution (per group),** in four sub-phases: a home `(h_r, h_c)` sends two requests per new candidate to the
  sources `(g_r, g_c)` via the relay `(g_r, h_c)`, and replies of at most 8 records (one leaf's scalars and
  points) return via `(h_r, g_c)`.
  Each sub-phase is the base count exchange (counts, bases in sender order, drop marks at 2^29 requests at a source
  and 2^34 records at a relay, three barriers). The source serves a request from its tables — X_k's affine point
  from the retained all-groups X affine store (section 5), Y_l from the Y table, and the scalars — copied into the
  reply (under 100 instructions of address forms, loads and stores, inside 2^10); old groups serve identically.
  As w3 gives s but not t, the home forms each record's three possible
  leaves `L = φ^t(X_k) + s*Y_l = [m]P + [n]Q`, t = 0, 1, 2, tests the nine pairs in a fixed order for a unit u with
  `L' = uL`, and with the first it finds that has `n' - u*n ≠ 0` outputs `x = (u*m - m')/(n' - u*n) mod r` after
  checking `[x]P = Q`.
5. **Stop test (per group):** if any home resolved an x this group, the stop flag (latched at a fixed address, set by
  the resolving home) broadcasts over 2 hops and every node halts after the barrier; the latched x is the output. The
  test is 5 instructions per group per node (LDR, CMP, B, group-loop ADD/STR), charged in W_g.

**Barriers:** at most 1943 (pre 4, per group at most 17, output 1), within the allowance of 2048. **Caps:** arithmetic fixups (3.4),
key paths (3.2), pair slots (cap_pg per group at hop 1, cap_p cumulative at hop 2), hash walks (uncapped, at most
HOME_MAX visits per record since a chain holds one home's records of one fold), candidates, requests, relay records
and coin retries are capped or bounded in total, so the worst-case ledger holds on every run and the E-ledger (7)
averages only group-granular executed work.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's (q2 even: both loops take two incidences per pass). With `D_l = y(Y_l)^2 - y(X)^2`, the
forward pass stores `D_l` and `acc_(l-1)` in entry l and forms `acc_l = acc_(l-1)*D_l` (`acc_(-1) = 1`); fe_inv gives
`R = 1/acc_(q2-1)`; the backward pass, for `l = q2-1..0`, sets `inv = R*acc_(l-1)`, `R <- R*D_l`, `u = inv*(y_Y - y_X)`
and `v = inv*(y_Y + y_X)`, and writes six **components**: with `a = x_X`, `b = x_Y`, `E = b^4 + 2a^3 b`,
`F = a^4 + 2ab^3`, `G = 3a^2 b^2` and w = u² (s = +1) or v² (s = −1), `α' = Ew`, `β'' = Fw` and `γ = Gw`. For t = 0,
1, 2 the leaves, never formed, are

    x(φ^t(X) + sY) = (α' - b) + β^t (β'' - a) + β^(2t) γ   (mod p).

**Why.** As R = 1/acc_l before step l, inv = 1/D_l, so `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)`. The curve gives
`y² = x³ + 7`, so `D = b³ - a³`, and as β³ = 1, `(b - β^t a)*Q_t = D` with `Q_t = b² + β^t ab + β^(2t) a²`. The slope of
`φ^t(X) ± Y` with `φ^t(X) = (β^t a, y_X)` is `(±y_Y - y_X)/(b - β^t a) = (±y_Y - y_X)*Q_t/D`, which is `Q_t*u` or
`-Q_t*v`; its square, less the two x-coordinates, is the leaf `Q_t²w - β^t a - b`, and expanding with `β^(3t) = 1` gives
`Q_t² = E + β^t F + β^(2t) G`, hence the display. D ≠ 0 unless `b³ = a³`, that is `Y = ±φ^t(X)` (3.4). A leaf's label
fractions are sums of its components' (3.3), so no leaf is formed.

**Arithmetic per pair**, after W2S (`inv`, the new R) and W2L (u and v, on the shared rows of inv): u² and v² (SQ13 and
F5: a square and its fold); E, F and G, each a product of a Y-entry word (b, b³ or 3b²) on the four tabled rows
`c*2^(64k) mod p`, k = 0..3, of an X constant c (2a³, 2a or a²), plus b⁴ or a⁴ for E and F (AD4), folded (FDX); then,
for each of E, F and G, its products by u² and by v² on its shared rows (row shifts by SHX, the pair schedule of W2L),
each folded and stored. The header H (bwd's caller's sp) holds the rows of 2a³, 2a and a², y_X and a⁴; a Y entry (192
bytes) `y², y, b, b³, 3b², b⁴`; batch entry l receives α'+, β''+, γ+, α'−, β''− and γ− at +160, +128, +96, +64, +32 and
+0 and stages R' at +224, and holds the correction words `-T(b)` and `-T(βb)` at +192 (3.3), written with the Y table
and never overwritten. The arithmetic is F's row schedule; the rare carries of SHX, FDX, F5, FS and FADDW are taken by
branches over fixups that give the inline chains' results (3.4), and its only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks the entries downwards, one sign at a time. For each of the sign's three components it forms U0 (its top
two words) and F1C's U1 (KEY's column F1 and the next column's shifted sum s, 3.3; 27 instructions); it adds the
correction words `-T(b)` and `-T(βb)` (entry +192) to α''s U0 and U1, and `-T(a)` and `-T(βa)` (K + 64) to β'''s, and
forms the three leaves' W0 and W1 by two-word additions and subtractions (28 in all, as `U2 = -U0 - U1` enters as a
subtraction of `U0 + U1`). For each leaf it then decides KEY's index j from the top words F0 of W0 and F1 of W1: with
`σ = F0 + F1`, if F0 ≥ F1 (CMP, B.HS) it branches to the block JH, where j = 2 if `F1 + σ` carries (CMN, B.HS to JT)
and j = 1 otherwise; else j = 2 if `F0 + σ` carries (CMN, B.HS to JT) and j = 0 otherwise, which goes on inline. The
j = 0 path emits the pair (W0, W1); JH forms `N' = ~(W0 + W1)` (ADDS, ADC, MVN) and emits (W1, N'); JT forms N' and
emits (N', W0). The **emit** forms 8*d1 by one AND (j = 0, 1) or BIC (j = 2) of the register x1 = (2^36 - 1) << 3 with
the low word of the pair's first fraction shifted right by 25 bits: with the column LC table at address 0 this is both
d1's LC entry and the destination word. LDADD advances the entry by 32 and returns the record's slot; the emit forms w3
from d2 (LSR #28, or for N' EOR of the shifted sum with 2^36 - 1) and the entry pointer, writes the record with two STPs
and tests for slot k (one TBZ of bit 13, section 4): 7 instructions. After slot k, FLUSH (8, or 7 for a dropped block)
ends the block (section 4). The entry pointer x3, stepped 32 bytes per component by post-indexed LDPs, is entry + 64 at
a sign +'s emits and entry − 96 at a sign −'s; entries are 256 bytes at a 2^34-aligned base, so `x3 << 31` mod 2^64,
which the emit adds to d2, is `2<<36 | l<<39` or `5<<36 | (l - 1 mod 2^25)<<39`. JH and JT add 5 to the key-path count
in x18 and branch back. **Meter.** The exact index is uniform on {0, 1, 2} for a leaf other than O, and a leaf takes JH
or JT only if its exact index is 1 or 2 or it lies in W (3.3): an event of its value, of probability at most `2/3 +
2d_wj`, whose leaves bound the meter's count. A node starts a batch only while its count is at most `5(cap_j - 6q2)`,
`cap_j = ceil(2n0/3) + 2^42` with the base n0 (fixed; grouped means are prefix values, looser and sound). Two leaves of one incidence are an invertible linear image of (X, Y), and leaves sharing
one generator are independent given it up to one excluded point, so a node's number of such leaves has mean at most
`n0(2/3 + 2d_wj)` and variance at most `(2/9 + 2d_wj)n0 + 1` (the dependent pairs' covariances are O(1/r) each; 6.4's
meter row, while no batch of the node is defective). Executed counts per incidence, q2 even (fixed part in
parentheses): `fwd` 98 (5) and `bwd` 1206.5 (5), each plus 3 or 4 per fixup (3.4), and `keys` 328.5 (11) for six
leaves on the j = 0 path, including their six emits of 7, plus the count it adds to x18 and its flushes; fe_inv
(`secp_fe.s`) is 22549. Per sign `keys` executes 11 LDPs (the components, the correction words, and G3 >> 3 and G0
reloaded), three F1C (81), the corrections (8), the sums (28) and per leaf the index's ADD, CMP, B.HS, CMN and B.HS
and the emit: 164, plus a TBZ per two incidences; a leaf on JH or JT executes 5 more. Kernels may overwrite x0-x29 and
the flags and preserve sp; `bwd`, `keys` and fe_inv save and restore x30 (`keys` at K + 96), which `fwd` and the
transport kernels do not write; `fwd` and `bwd` return their fixup word in x0 (3.4); the driver keeps no registers
across calls, loads the key-path count into x18 for `keys` and sets sp to H for `bwd` and to K for `keys`, the block
of G1..G5, G0 (+0..+47), G3 >> 3 and G0 (+48) and X's correction words (+64). `bwd` keeps R in x25..x28 from pass to
pass; `bwd` and `keys` end with a TBZ on bit 33 of the entry pointer: the 2^25 entries of 256 bytes fill the lower
2^33 bytes of their 2^34-aligned block, so bit 33 is clear at every entry and set one entry below the block (both loop
ends test bit 33: `loop_ends_on_bit_33`).

### 3.3 The class label (fractions from components, the index, the paths)

**Fractions and the label.** For v in [0, p) let `T(v) = floor(2^128 v/p)` and `comp(v) = floor(2^100 v/p) = T(v) >>
28`. For a leaf value x in [1, p) let `v_m = β^m x mod p`; the canonical member `v* = min v_m` has index j = argmin, and
the class's **label** is `(comp(v*), comp(βv* mod p))`: the member `β^m x` has index `j - m` and the same pair. As
`2^128 v/p = v/2^128 + C0*v/(2^128 p)` with the second term in [0, 2^-95), `T(v) = (v >> 128) + e0(v)`, e0 in {0, 1}.
For a stored component c below 2^256 let `c̄ = c mod p` and `T_m(c) = T(β^m c̄ mod p)`.

**A component's fractions.** `keys` forms `U_0(c) = c >> 128`, `U_1(c) = 2^64*F1 + 8s mod 2^128` (F1C below) and `U_2 =
-U_0 - U_1 mod 2^128`. The table below gives each estimate's error T − estimate mod 2^128: for c ≥ p, c̄ < C0, so `T_0 =
0` and `U_0 = 2^128 - 1`; and the three values `β^m c̄ mod p` sum to 0, p or 2p, so `T_0 + T_1 + T_2 ≡ -k` with k in {0,
1, 2} the sum of their fractional parts.

**The column (F1C, 27).** `G = floor(β*2^384/p)` has limbs G0..G5 (x16, x7..x11). `F1 = Σ_i [lo(c_i*G_(5-i)) +
hi(c_i*G_(4-i))] mod 2^64` is column 5 of `c*G` without its carry-in. Column 4 has the eight words `lo(c_i*G_(4-i))` and
`hi(c_i*G_(3-i))`; F1C adds `hi(c_0*floor(G3/8))` (x15 = G3 >> 3) and the other seven shifted right by 3 into s: each
term is below 2^61, so `s < 2^64`, and U_1's top word is `c6 = F1 + (s >> 61)`. **Bound.** As `c*G/2^384 = cβ/p -
cδ/2^384` with δ in [0, 1) and `cβ/p ≡ c̄β/p (mod 1)`, `floor(2^128 {c*G/2^384})` is `T_1` or `T_1 - 1` mod 2^128, and
it is columns 4 and 5 of `c*G` with the carry from columns 0 to 3. U_1 leaves out that carry (columns 0 to 3 have 1, 3,
5 and 7 words, so their carries out are at most 0, 2, 4 and 6) and the bits the shifts drop: `Σ (t mod 8)` over seven
words, at most 49, and `hi(c_0 G3) - 8 hi(c_0 floor(G3/8))`, at most 12 as `G3 mod 8 = 5`.

**Corrections and sums.** A sign's leaves are `L_t = α + β^t β' + β^(2t) γ` with `α = α' - b` and `β' = β'' - a` (3.1).
`keys` adds `-T_m(b) mod 2^128` to `U_m(α')` and `-T_m(a)` to `U_m(β'')`, m = 0, 1, before forming β''s U_2; as `T_m(c -
y) ≡ T_m(c) - T_m(y) - {0, 1}`, the corrected values estimate α's and β''s fractions. As `β^m L_t` is the sum of `β^m
α`, `β^(m+t) β'` and `β^(m+2t) γ` in F_p, taking representatives in [0, p) gives `T(v_m) ≡ T_m(α) + T_(m+t)(β') +
T_(m+2t)(γ) + κ`, κ in {0, 1, 2}; so `W_m = U_m(α) + U_(m+t)(β') + U_(m+2t)(γ) mod 2^128` errs by its three terms'
errors plus κ, and `N' = ~(W_0 + W_1)` by `-(T(v_0) - W_0) - (T(v_1) - W_1) - k + 1`. For every stored c below 2^256:

| estimate | T − estimate, mod 2^128 | from |
| --- | --- | --- |
| `U_0(c)`, `U_1(c)`, `U_2(c)` | [0, 1], [0, 68], [−71, 0] | e0 (and c ≥ p above); the column; `-(T_0 - U_0) - (T_1 - U_1) - k` |
| corrected `U_0`, `U_1` of α, β'; `U_2` of β' | [−1, 1], [−1, 68]; [−71, 2] | less {0, 1}; as `U_2`, from the corrected terms |
| `W_0` for t = 0, 1, 2 | [−2, 5], [−73, 71], [−72, 73] | `U_0(α) + U_t(β') + U_2t(γ)`, plus κ |
| `W_1` for t = 0, 1, 2 | [−2, 206], [−72, 73], [−73, 71] | `U_1(α) + U_(1+t)(β') + U_(1+2t)(γ)`, plus κ |
| `N'` for t = 0, 1, 2 | [−212, 5], [−145, 146], [−145, 146] | the two rows above, less k, plus 1 |

Over t: `T(v_0) - W_0 ∈ [-73, 73]`, `T(v_1) - W_1 ∈ [-73, 206]` and `T(v_2) - N' ∈ [-212, 146]` (indices of U mod 3).
The certificate recomputes this table from the component ranges and ties D''s 783 residues and the 359-unit widening to
it.

**The index.** KEY decides j from the top words F0 of W_0 and F1 of W_1 (3.2): with `σ = F0 + F1 mod 2^64`, j = 2 if
`min(F0, F1) + σ` carries, else 1 if F0 ≥ F1, else 0. Let `X_m = 2^64 v_m/p`; `X_0 + X_1 + X_2` is 2^64 or 2^65. The
**window set** W is where (W1) X_0 or X_1 is below 1 or above `2^64 - 1`, (W2) X_2 is below 1 or above `2^64 - 4`, (W3)
`|X_0 - X_1| < 2`, or (W4) `X_0 - X_2` or `X_1 - X_2` lies in (−1, 4). **Exactness.** Outside (W1),
`2^64 - 1 <= T(v_m) <= 2^128 - 2^64` for m = 0, 1, so `W_m = T(v_m) - e_m` without wrap, and as
`2^128 v_m/p - T(v_m) ∈ [0, 1)`, both `a' = X_0 - F0` and `b' = X_1 - F1` lie in `[-73*2^-64, 1 + 207*2^-64)`. Outside
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
across a block of 2^28: the **label-error set** D' is where `T(v_0) mod 2^28` lies in `[0, 73) ∪ [2^28 - 73, 2^28)`,
`T(v_1) mod 2^28` in `[0, 206) ∪ [2^28 - 73, 2^28)`, or `T(v_2) mod 2^28` in `[0, 146) ∪ [2^28 - 212, 2^28)` (a wrap mod
2^128 is at these residues too). Outside W and D' a leaf emits its class's label, whatever its stored components. As
each value of T is taken by at most 2^128 values of v and `x -> v_m` is a bijection, D' has at most `783*2^228` values
of x, each of probability at most `2/(r-1)`: the row `d_lab < 2^-17.38` of 6.4, a missed collision like `d_wj`; no leaf
is dropped.

**Mass bound.** secp256k1's Frobenius is `π = A + Bω` with `B = 0xe4437ed6010e88286f547fa90abfe4c3` and
`A = 0x114ca50f7a8e2f3f657c1108d9d44cfd9`: `A² - AB + B² = p`, `B*β ≡ A - B` and `(A - B)*β ≡ -A (mod p)` (checked by
the run and the certificate). So `Λ = {(v, u) ∈ Z²: u ≡ βv (mod p)}` has the basis `(B, A - B), (A - B, -A)` of
determinant -p, and as `A > B > 0` its half-open parallelogram P spans [0, A) in v and (-A, A - B) in u. For each point
λ of Λ in a box `[v0, v0 + a) × [u0, u0 + b)` of integer sides a and b, `λ + P` lies in `[v0, v0 + a - 1 + A) × (u0 - A,
u0 + b - 1 + A - B)`; these translates are disjoint, as P tiles the plane under Λ, and each has area p, so the box holds
at most `(a + A)(b + 2A - B)/p` points of Λ. Fixing f of a label's 200 bits leaves `2^(200-f)` pairs of component
values; a component value is an interval of at most 2^156 values of v (of u = βv mod p for the second), widened at each
end by its fraction's error range, at most 359 units of `p/2^128` in all (N''s [−212, 146]), for a computed fraction
(one that wraps mod 2^128 may be two end pieces of [0, p), which count as one interval as `pZ² ⊂ Λ`). So on each path at
most `2^(200-f)(2^156 + w + A)(2^156 + w + 2A - B)/p < (1 + 2^-18)2^-f p` values of `v_j`, and of x, give those bits,
with w the widening (the certificate checks the factor), whatever index KEY returns and whatever the components. With
three paths and two points per x, the emitted pair takes given values on f of its bits with probability at most `KB(f) =
6*2^(256-f)(1 + 2^-18)/(r-1)`.

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts
`B_(k+1) = SHX(B_k)`; the rows of an X constant (E, F, G's) are tabled exactly. RF and RA sum `Z = Σ a_k*B_k < 2^322`
exactly (each high word is at most 2^64 - 2), AD4 adds b⁴ or a⁴ (below p) for E and F, `Z' < 2^322 + 2^256`, and a fold
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ13 forms the square of a weak a as
`V = L + H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), and F5 folds V the same way. Only these
steps drop a carry: SHX's last ADCS and the final ADC of FDX and of F5. Without a dropped carry each returns a value
below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; of
  `Z' < 2^322 + 2^256` (E, F, G), in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`.

A fold returns a value `>= p` only when the true value is below C0; a component stored so is still its value mod p, and
3.3's bounds hold for every stored value below 2^256, so no canonical value is needed. Every other operation is exact on
its inputs: RF and RA; AD4; FS and FADDW, whose operands (`y_Y²` and `y_X²` in `fwd`, `y_Y` and `y_X` in `bwd`) are
canonical table values; the label's sums, corrections and path additions (3.3); KEY's integer arithmetic (3.3 covers
its index); the emit; and F's `fe_mul_f` and `fe_sqr_f`, which fold the 512-bit product, fold its
fifth limb again and add C0 exactly when that last sum carries past 2^256, so fe_inv (255 squarings, 15 multiplications)
is exact for every input below 2^256. Generator construction, the table entries, the header rows and values, the
correction words, the resolver and `[x]P = Q` run on `fe_mul_f`, `fe_sqr_f` and fe_inv, each coordinate brought below p
(an ADDS chain of C0 detects `t >= p`, CSELs keep `t - p`) before it is subtracted, compared or stored, with exact
arithmetic mod r for the solve. So they have no sites; the tables hold exact canonical values (no point has x = 0, as 7
is not a square mod p, so a, b ≠ 0).

**Distribution.** Given X ≠ O and every other generator, `Y_l` is uniform on the r - 1 points of G \ {O}. A value of
x_Y comes from at most 2 points, and of y_Y from at most 3. So `D_l = x_Y³ - x_X³`, `acc_l = acc_(l-1)*D_l`,
`R_l = 1/acc_l` and `inv_l = 1/D_l` take each value with probability at most `6/(r-1)` (the rest are bijections while
no `D_m = 0` for m ≤ l, the degenerate case below); `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)` (as
`D = (y_Y - y_X)(y_Y + y_X)`) at most `3/(r-1)`; u² and v² at most `6/(r-1)`; E, F and G, of degree 4, 3 and 2 in x_Y,
at most `8/(r-1)`; and each component at most `8/(r-1)`. For `c = E*u²`: `c(y_Y + y_X)² = E` and `y_Y² = x_Y³ + 7` give
`2c*y_X*y_Y = E - c(x_Y³ + 7 + y_X²)`, so for c ≠ 0 (y_X ≠ 0, as r is odd) y_Y is fixed by x_Y and squaring gives a
nonzero polynomial in x_Y of degree 8 (leading term x_Y⁸); for c = 0, `E = 0` has at most 3 roots. F and G (and v², by
the sign of y_X) give nonzero polynomials of degree at most 6 the same way (when the x_Y⁶ term of F's cancels, the x_Y³
term `-4c²y_X²` is left). A row `B_k ≡ v*2^(64k)` determines its operand v, so at most 2^97 values of v make a given row
shift carry.

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
| bwd: folds of E, F, G | 3q2 | `Y_l`: E, F, G | (2^66 + 1)*C0 | 8 |
| bwd: row shifts of E, F, G | 9q2 | `Y_l`: E, F, G | 2^97 | 8 |
| bwd: folds of the six components | 6q2 | `Y_l`: the components | 2^66*C0 | 8 |


A run fails with probability at most its bad count times its points per value over r - 1. At l = 0 three row shifts
and three folds cannot fail, and are not in the table:
- fwd's rows of `acc_(-1) = 1` are 1, 2^64 and 2^128, and its fold of `D_0*1` has `Z < 2^256`;
- bwd's fold of `inv_0 = R_0*1` has `Z = R_0 < 2^256`;
- bwd's last update `R_0*D_0 ≡ 1` folds to 1 or p + 1, as nothing else congruent to 1 lies below `2^256 + 2^66*C0`. It
  is left in x25..x28 after the last pass and never read.

The degenerate cases are X = O (probability 1/r) and `D_l = 0`, that is `Y_l = ±φ^t(X)` for some t (at most `6/(r-1)`
per incidence), which zeroes acc and spoils the batch. The first site of a batch whose predicate holds has exact
operands, so a batch has a defect with probability at most
`d_batch`: the table's sum plus `6q2/(r-1) + 1/r` (`SITES` in `certificate.py`; 6.4).

**Rare carries by branch.** In SHX the carry out of its first ADDS, in FDX and F5 the carry into limb 2, in FS the
borrow of its C0 subtraction and in FADDW the carry of its C0 addition are rare. Each macro branches over its **fixup**
(B.LO; B.HS for FS) unless that carry or borrow occurred. The fixup is the rest of the inline chain, two instructions
(ADCS and ADCS or ADC) or three for FS and FADDW, and `STADD x4, [sp]`. Where the branch is taken the skipped
instructions would add or subtract 0, so on every input every register equals the inline chain's, and so does every
deviation of the table above; only the flags differ, and every path from the join sets them before reading them or
returns (`fixups.flags_dead_at_every_join`; kernels may leave the flags changed). A fixup adds 3 instructions (SHX, FDX,
F5) or 4 (FS, FADDW) and C0 to the word at [sp], which opens a 16-byte frame pushed as zero by `fwd` (STR) and by
`bwd` beside its saved x30 (STP); each returns the word in x0, at most `37q2*C0 < 2^63` per call.
The driver adds x0/C0 to the node's fixup count N and starts a batch only while `N <= 2^45 - 37q2` (37 sites per
incidence: 18 SHX, 14 FDX, 2 F5, 2 FS and 1 FADDW). So N is at most 2^45 on every run, and section 7 charges `4*2^45`
per node.

**How often.** A fixup needs a word in a small set:
- SHX's carry needs the row B it shifts to have its low word above `2^64 - C0`, as the word it adds is below C0: at most
  `C0*2^192` rows.
- FDX's carry into limb 2 leaves its output with low 128 bits below `F*C0 <= (2^66 + 1)*C0` and the rest nonzero, where
  `F = floor(Z/2^256) <= 2^66 + 1` (`Z < 2^322`, or `Z' < 2^322 + 2^256` for E and F): at most
    `F_w = (2^66 + 1)*2^128*C0` words.
- F5's carry into limb 2 leaves its output's low 128 bits below 2^65, as its fifth word is at most C0 + 1, so that word
  times C0 has a high word of at most 1: at most 2^193 words.
- FS's borrow leaves an output whose low word is at least `2^64 - C0`, and FADDW's carry one whose low word is below C0:
  at most `C0*2^192` words each.

While no site of the batch up to and including this one has deviated, the row or output is a word below 2^256 congruent
to the site's operand value (times `2^(64k)` for the k-th row), which therefore determines that value. So, given X and
every other generator, a site takes its fixup with probability at most its set's size times the points of `Y_l` per
operand value over r - 1, with the fresh generator and multiplicities of the table above; `y_Y ± y_X` have 3 points per
value. At l = 0 the six sites the table leaves out cannot take a fixup either: rows 1, 2^64
and 2^128 have a low word of at most 1, folds of values below 2^256 add nothing, and `R_0*D_0` folds to 1 or p + 1,
neither of which is at least 2^128 with its low 128 bits below `2^66*C0`. Summed over a batch (`FIXUPS` and `E_FIX` in
`certificate.py`), the expected number of fixups a batch takes before its first deviation, and so of a batch without a
defect, is at most `E_FIX`, and over a node's batches at most `q1*E_FIX < 2^29.15`.

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
| bwd: folds of E, F, G | 3q2 | `F_w` | 8 |
| bwd: row shifts of E, F, G | 9q2 | `C0*2^192` | 8 |
| bwd: folds of the six components | 6q2 | `F_w` | 8 |

**The meter row.** A node stops early only if N exceeds `2^45 - 37q2` at the start of some batch. On the event that none
of its batches has a defect, N is that clean count, so by Markov's inequality a node stops with probability at most
`q1*E_FIX/(2^45 - 37q2) + q1*d_batch < 2^-15.85` (`d_f`, a node-local row of 6.4).

## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (7), and per block of k = 256 records its flush (8) and `dlvk` (1026) | 11.0390625 | 11.046875 |
| hop 2: `h2s` (7; loop 0.25), and per block its flush (8; 9 if dropped) and `dlvk` (1026) | 11.2890625 | 11.296875 |
| hash join `hchg` fixed (11; loop 0.25) | 11.25 | 11.25 |
| hash walks (expected; probe 6, advance 3, chain end 2) | ~5.4 | E-rows |
| K2 verifications and listings (expected; global pair counting) | ~0.3 | E-rows |
| **fixed total** | **33.578125** | **33.59375** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlvk`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. A destination's full blocks number at most its records over k, and a
dropped block costs 7 or 9 and sends nothing, so 1036 per k records bounds each hop's blocks; final blocks are `fin`'s
rows (section 7). `h2s` and `hchg` take a slot's records in groups
of 1, 2 and 4, as the low bits of its count say, then 8 at a time; the 8-record loop's CMP and B.LO, 2 per pass, are
charged as 0.25 per record (a slot of n records runs floor(n/8) passes), and the tests and set-up, 11 and 13 per
slot (`hchg` loads the slot's read cursor and stores E2 back), are per-node rows. **Addresses.** One instruction
forms each per-record table address: the emit's AND or BIC (3.2); `h2s`'s BFI of d2 (w3's low 36 bits) into bits
3..38 of x10, which holds 2^39 elsewhere (d2's row LC entry and destination word); the flushes' ORR of 2^43 (the RC
entry); `hchg`'s AND of the bucket mask with K2 and ADD of the head address (8-byte heads). **Blocks and cursors.**
The emit and `h2s` take a record's slot by one LDADD of 32 on its LC entry, which holds the absolute slot address:
block base (low 13 bits zero) plus offset 32n below 2^14, so bit 13 is set exactly after slot k (32k = 2^13). One TBZ
of bit 13 of LDADD's old value finds it (checked); there the **flush** returns
the LC entry to slot 1 (SUB, SUB, STR), forms the RC entry (ORR) and advances the RC cursor by 32k with LDADD, which returns
the old cursor c. If c is below the sender's slot end, which is the same at every destination and held in a register,
the STP of the header (destination word, c) sends the block; otherwise it is dropped (CMP and B.HS, or B.LO in `h2s`),
and `h2s` stores the end back, so a hop-2 cursor stays at most its slot end; a hop-1 cursor can pass its end, but only
by the node's own leaves, 32 bytes each, so it stays below 2^63.2 (checked). `fin` takes the L = (LC - slot 1)/32
records left in each destination's block, sends it by the header (8d, c) if L > 0 and c is below the end, and writes
the final cursor `min(c + 32L, end)` (CMP, CSEL, STR) both to the sender RC and the remote E2. **Keep lemma.** Each guarded R2 slot reserves k extra records. Group-final partial deliveries may start unaligned and a full-width delivery can reach the guard; only records below the clamped logical end are hashed. Later groups overwrite the stale suffix starting at the saved cursor after the barrier. Hop-1 slots restart aligned at each group, and cap_pg is a multiple of k. The guard proof, physical next-array mapping and edge tests are in guard-repair.md. **Hash table.** Each node
owns 2^57 head words, zeroed once by `zero4096` (2048 STPs per 4096 words, plus SUBS, B.NE and one final RET), and a packed next-array
with one 8-byte next word per R2 slot position (section 5). Bucket `K2 mod 2^57` (AND with 2^57 - 1, held in x5): same K2 always
lands in the same bucket, so same-key records meet with no ordering step. **Join.** `hchg` steps its record pointer
with a pre-indexed LDP and reads each slot from its read cursor to its E2 cursor, then stores E2 back as the new
cursor (LDR of E2, LDR of the read cursor, STR of E2: the per-slot 13); the per-group ranges partition each slot, so
every record is read exactly once across all group calls, however the slots' previous ends differ, and cross-group
state otherwise travels only through the cumulative table. Per record: LDP, AND-mask, head address/load/save,
CBZ to insert when empty, else the walk: each visit forms the next-address (LSR, ADD of C), loads next word and K2,
CMP, B.NE; a K2 hit verifies d4 (LDR, two UBFX, CMP, B.NE, 5); a full hit lists (CMP, B.HS, STP, 4) unless the list
is full (CMP, B.HS, 2); then prepends (head-address recompute, next-address pair, two STRs). A non-hit visit advances
(MOV, CBZ, B back, 9 with the probe) to the first full match or the chain end; uncapped, so every new record links to
its most recent same-key predecessor on every run. The packed next-array (one 8-byte word per R2 slot position, at
(A>>2)+C, C a call argument) fits every receive pattern up to the slot caps by construction; next words need no
zeroing (written at insert before any read); records stay intact in R2 for phase 4. A hop's
receiving node executes nothing on arrival; it reads the records in the next phase, with the LDPs of `h2s` (hop 1)
or `hchg` (hop 2). Per-node and per-group costs are rows of section 7, and every address a kernel loads or stores
lies in the regions of section 5. (`bkt`, `pass4`, `scan_ranges` stay in `route.s` for the base path evidence of
part B; the grouped run never executes them.)

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
| below 2^52 (2^51.3297 used) | code, stack, the five barrier words (one reservation granule each), header H and keys' block K, X table of the running group (512 bytes per X), X affine store (64 bytes per X, all groups), Y table (192 bytes per Y), scalars, batch entries (2^33 bytes at a multiple of 2^34), the ports' message blocks (at 2^48 - 64 and 2^48 - 128), their LC tables at 0 and 2^39 and RC tables at 2^43 and 2^43 + 2^39 (2^39 bytes each), the hash-read cursors at 2^48 (2^39 bytes), the record blocks at 2^50 and 2^51 (2^50 bytes each), E1 and E2 (2^41 each) and the phase-4 buffers and tables (`FIXED` in `certificate.py`, whose d3/d4 tables stay reserved but unused by the grouped run); all but the aligned tables lie between 2^44 and 2^48 | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | per-group hop-1 slots of cap_pg (`2^59.9986`), hash heads (`2^60`), hash next-array (`2^61`), together 2^61.9989 | per group |
| R2 after R1, ending below 2^63.998987 | cumulative hop-2 pair slots (phase 4 reads candidates' records at their addresses) | groups 1-114 |

R2 uses physical stride `32*(cap_p+k)`; its exact end and the enlarged next-array bounds are in guard-repair.md. So every node stays below 2^64 bytes, and the peak is `2^72*(fixed + |R1| + |R2|)`, at most 2^135.9989
bytes (log2 ≈ 135.9988560). A remote item is reached only by a message.

## 6. Success (no heuristic)

In `success.md`: the model (6.1), Chen–Stein on cross-node pairs (6.2), usefulness (6.3), the defect rows (6.4), the
bound (6.5) and the grouped hash-chain abort (6.6: stopping-time bound, chain detection E1-E7, uncapped walks with
history-conditional walk means and a defective-run allowance).

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 645,410 per point (645,414 for a grouped X, whose four affine-store STPs the X row charges):** a complete Jacobian
  double-and-add, exceptional-case tests, two inversions with the affine conversion, all on the exact arithmetic of 3.4:
  `fe_mul_f` and `fe_sqr_f` at 105 and 80 instructions as executed, plus 4 and 3 for the call as in fe_inv, and an
  addition or subtraction with its operands brought below p at 27; 633,952 for the 256 doublings, 256 mixed additions
  and two conversions; the exceptional-case tests (Z and H brought below p and tested for zero per iteration, 13 each,
  and 64 for the paths they select) add 6,720. The scalar's bit tests and the pointer step add at most 6 per iteration,
  and with set-up and return 1,600 in all. A Y's table entry (y², b², b⁴, b³, βb, 3b², the correction words `-T(b)` and
  `-T(βb)` at 30 each: `T(v) = (v >> 128) + e`, e = 1 iff `2^128 (v mod 2^128) + (v >> 128)*C0 >= p`, by 3.3's e0; seven
  values brought below p at 9 each; fourteen STPs: 658) or an X's mixed addition `R_g + [a_k]P` with its header (a², a⁴,
  y², a³, βa, 2a³, 2a, the nine rows `c*2^(64k) mod p`, k = 1..3, of 2a³, 2a and a², `-T(a)` and `-T(βa)`, seventeen
  values brought below p, thirty-two STPs: 3,138) brings the total to at most 645,410, the allowance (charged at its
  itemization as in aburan28's 7e8d13c7 and ecc3ec1b).
- **Node set-up, 2^20:** a program of at most 2^16 bytes (kernel code, generator construction,
  resolver, phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD of 8 to its destination word, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 1,491,011:** two requests and two replies of at most 8 records over 2 hops at 40; each record's three
  possible leaves `φ^t(X_k) + sY_l` (one inversion of `D = x_Y³ - x_X³`, the inverses `Q_t/D` of the three denominators,
  and three affine additions: 25,466); nine unit tests at 600; the solve mod r and `[x]P = Q`, on the same exact
  arithmetic; 745,505, charged at twice that plus 1 (as in 1fbd26e8).
- **Phase-4 count exchange, 2^41 per node and sub-phase:** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 88 per message (its count,
  the count message, the offset reply and the write, at most 22 each with in-flight count and delivery) for at most
  8*2^29 messages (2^29 requests of up to 8 records): below 2^40.5.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^16 per node:** 12 kernel calls per group at most 33 each (arguments and BL), the stop test and
  group loop (8), over 114 groups: below 2^16.
- **Hash-head zeroing, 2050 per 4096 words:** `zero4096` (2048 STPs, SUBS, B.NE); the row's +1 is its RET.
- **Hash-read-cursor init, 4 per slot:** `init3` (STR, ADD, SUBS, B.NE) of each slot's base into its 8-byte entry at
  2^48; the row's +1 is its RET.
- **Keys' block 29** (G0..G5: six words of MOVZ and three MOVK, LSR of G3, four STPs),
  **batch set-up 160** (X's header copied into H by 28 LDPs and 28 STPs and its correction words into K, y_X² into
  registers, acc = 1, pointers, sp for `bwd` and `keys`, the index word, the slot end for `keys`, both meters' updates
  and tests and the loop: at most 155, itemized in `certificate.py`), **inversion call 16, request served 2^10, output
  2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on
  the failure flag, plus a range test. RNDR is this machine's coin source: it returns 64 fresh uniform bits, the
  algorithm coins of the rules' probability space, and never reports failure, so the charged B.EQ is never taken and
  coins add no failure row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, the barriers and the output. A record row uses
its group's records; E-averaged rows are sums over groups weighted by PHI. With q1g the grouped q1 (145307752):

| row (P: preprocessing; E: E-averaged; W: worst-case) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1g + q2 + 1)*64 (P) | 105.4143 |
| Y generators with their table entries: q2*645410 (P) | 116.2999 |
| R_g generator | 91.2999 |
| keys' block: 29 (P) | 76.8580 |
| row LC and RC tables: (4*2^36 + 1) + (3*2^36 + 1) (P) | 110.8074 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| hash-head zeroing: zero4096 2^57/4096*2050 + 1 (P) | 128.0014 |
| hash-read-cursor init: init3 4*2^36 + 1 (P) | 110.0000 |
| leaf work, E-averaged: sum_g PHI[g] (X gens, fwd, inv, bwd, keys w/o emits, batch set-up) (E) | 134.5306 |
| arithmetic fixups: 2^45*4 (W) | 119.0000 |
| j = 1 and j = 2 key paths: cap_j*5 (W) | 128.4368 |
| route work, E-averaged: sum_g PHI[g] (emits, blocks, h2s, fin, column tables) (E) | 130.9979 |
| group work, E-averaged: sum_g PHI[g] (hchg fixed, main-pair walks, stop test) (E) | 130.5068 |
| same-X walks: global pair counting (E) | 30.0394 |
| K2 verifications: global pair counting (E) | 124.8909 |
| listings: global pair counting (E) | 75.9840 |
| defective-run hash rows: DEFV + DEFK2 + DEFL (E) | 63.8676 |
| candidates: 2^28*1491011 (W) | 120.5079 |
| requests served: 2^29*2^10 (W) | 111.0000 |
| phase-4 count exchanges: 114*4*2^41 (W) | 121.8329 |
| barriers: 2048*2^82 (W) | 93.0000 |
| phase control: 2^16 (W) | 88.0000 |
| output + stop broadcast (W) | 83.3219 |
| **expected total E** | **134.763035** |

`E` is the exact `E_total` of `certificate.py` (expected instructions), at most 267.5303 per leaf (the run's
`per_leaf`), and `W_E = E/1000` with log2 W_E ≈ 124.7972511546. **The claim, 124.797251154654, is the least multiple
of 10^-12 with 2^claim >= W_E**; the certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds
on ln 2. Rows marked P are preprocessing: their sum over 1000 has log2 ≈ 118.0360724 and is declared as 2^118.0361,
the least multiple of 10^-4 at or above it (checked). The worst case R_worst = 2^188.604980 (all groups run, uncapped
walks to HOME_MAX visits with verification and listing on every record) is disclosed alongside, not scored; it is the
least multiple of 10^-6 at or above log2 R_worst (checked).

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "seconds": 15}` (`seconds` 5 to 15). Budget-limited checks run
guaranteed minimums (2,000 random and 600 near-tie class-rule inputs, 3,000 class-label classes), then their share of
the remaining budget, up to 26,000 inputs and 12,000 classes. **Part A** runs the leaf kernels and fe_inv on production
batches of 2, 4, 6 and 10 incidences: every component against exact arithmetic (`component_mismatches`), 3.1's identity
on the stored components against exact affine addition of `φ^t(X)` and `±Y` (`identity_mismatches`), each record, slot
and LC entry against the bit-exact model of `keys` and the emit, and each leaf's index and label against the exact class
rule and its class's label outside W and D' (`wrong_index`, `keys_not_canonical`). Random leaves and near ties of the
six forms whose small values make W, each given random stored components (two random below p or a small value plus p,
the third solved) and correction words, get the exact index outside W (`class_rule`); on at least 3,000 classes, a third
with a member at an edge of W or D', every member so represented emits its class's label unless it lies in W or D', and
the computed fractions' errors stay in [−73, 73], [−73, 206] and [−212, 146] (`class_labels`). Two batches replace each
sign's components by such a representation of a leaf at an edge of W or of D', four of them inside it, and check
records, key-path count and count against the model and formula (`window_batch`, `boundary_batch`). The `fwd` and `bwd`
calls of the production batches (`fwd_bwd_formulas_hold`), the planted-point batch and the forced copy below count their
formula plus 3 or 4 per fixup they execute and return C0 times their fixups. Planted points take chosen fixups in both
passes (fwd's FS and the next incidence's row shift, bwd's FS and FADDW) with exact components and counts
(`fixup_batch`). A copy of `leaf.s` whose fixup branches are never taken runs all 37 fixups in place on 2 incidences, as
inline chains: exact components, 114 more instructions and 37 fixups per incidence (`forced_fixups`). For each site type
of 3.4 (`site_checks` in `sites.py`: row shift, fold, the fold of E and F with their constant, both square folds),
inputs built to meet its predicate deviate exactly when it holds, other boundary and random inputs are exact, and all
match its bit-level model; each fixup runs exactly when its carry or borrow occurs, adds C0 to the word and 3 or 4 to
the count, and its row or output lies in its set of 3.4 (`fixups`); FS and FADDW on canonical inputs, `fe_mul_f` and
`fe_sqr_f` are exact. A batch with `Y_l = φ(X)` and `Y_m = -X` gets acc = 0 and loses all its other components; the last
R update is ≡ 1 and unread; the SHX, FDX, F5, FS and FADDW invocations per incidence match the tables
(`sites_match_certificate`), and every path from a fixup's join sets the flags before reading them or returns
(`fixups.flags_dead_at_every_join`). **Part B** (`machine.py`) runs a 4 × 4 machine through `leaf.s` and `route.s`
(`init3`, `tinit`, `zero4`; `keys` on 4 incidences per node whose components are planted by the inverse DFT of 3.1 with
zero correction words, for leaves with labels with d1, d2, d3 < 4 and d4 < 16 at production field widths, each class's
v* the point of 3.3's lattice nearest chosen quarter points of its parts' blocks, so that its fractions lie far from
block boundaries and every record carries its class's label (checked), some classes shared; `fin`; `h2s`; `fin`; `bkt`,
`prefix4`, `pass4` and `scan_ranges` at each home), with pair slots of k records (so its blocks are final blocks) and
buckets of 8, delivering each message by executing `dlvk` or `dlv`, which must return x0 to the link's window.
`block_paths` forces the other paths: `keys` with 534, 140, 128, 2 and 0 records to five destinations, on the j = 0,
1, 2, 0 and 1 paths, from cursors 1, 2, 1, 1 and 3 blocks below the end (blocks sent and dropped, a cursor past the end,
final blocks sent with new, older and unwritten slots, held, and empty), then `fin`; `h2s` on slots of 0 to 19 records
and two of 300 and 280, with 560, 140, 50 and 20 records to four full-width columns (blocks sent, and dropped with the
end stored back). Messages, LC and RC entries, slots, pointer order and candidates equal an exact model, each count its
formula, and delivered slots the one-record rule's, each delivery store inside its physical slot including the guard. `unrolled_paths`
runs `bkt` and `pass4` (full-width d3, d4) on slots and buckets of 0 to 19 records (every group combination, up to two
8-record passes), some filling (a full bucket's LDADD undone, its cursor left at its end), and `scan_ranges` alone on
ranges of 0 to 129 records (every remainder of pairs with either parity, up to two 64-record blocks) with K2 from three
values or random, on all-equal ranges, with a full list at a pair's first and second append, and with no range:
candidates and counts equal the path model, and every instruction of `scan_ranges` runs. `grouped_paths` runs three
groups through `zero4096`, `init3` (cursors) and `hchg` on the true cumulative R2 layout (fixed slot bases,
persistent E2, per-slot cursors, 16 toy buckets, planted cross-group same-key pairs with differing K3 high bits):
the stop case halts after two groups, the full case runs all three, and a full-list case runs all three with a
1-entry list; heads, next words and candidates equal the arrival-order chain model, every count its certificate
formula, cursors advance to E2, head-following reaches every appended address exactly once, slot records are
unmodified, every listed pair is same-key, and every true-global-grouping pair is connected in the listed graph.
`grouped_paths` also plants per-group X affine entries: group-1 entries persist after group 3 and a mock request
returns them (`old_x_served`).
**Part D** runs `certificate.py`; `D_certificate.checks` lists each of its checks with its value, and
`counts_match_certificate` compares its constants with parts A and B.

The 4 × 4 run is evidence at its own scope only.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1, whose links keep no state beyond the message in flight. The sender's
  instructions (including those that form the address), each delivery (as `dlv` or `dlvk`) and the ports' in-flight
  counts are charged. Each hop is one link from the sender's port to its destination; a network that forwarded a message
  through further nodes would add their work at each extra hop, and that machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's meter tests, stop test and phase-4 exchanges), the in-flight counts, and the production attack,
  machine and links. The toy stop fires on a listed candidate; production stops only on a resolver-confirmed x, a
  stronger rule, so the E-bound (which charges confirmed stopping) is conservative over the toy.
- **Provenance**, not a premise: kaankolcu's d91b6dc5, the base (this proof's 2.1-2.2, 3.1-3.4, resolver, success rows;
  itself merging bradcush's 25e43d87 Fourier-component leaves and kaankolcu's k=256/TBZ machine with 46-bit d4, with
  e59be87d's record blocks/`dlvk`/keep lemma, bbe1cd85's two-fraction label, 10c5d938's address forms, package F's rows,
  02feb284's fe_inv, 80e9e8af/1efcef88/892a3df8/b00400fa's scan_ranges, 1fbd26e8's candidate, 7e8d13c7's generator and
  807e4bf5/a9297fb0's meter); kaankolcu's b63af001, the v1 grouped attempt whose E-ledger this package keeps
  (hash chains repair its concatenation flaw); kaankolcu's 7c8264bf,
  whose hash-chain design this package keeps (v4 repaired its cursors, conditioning and overflow tail); kaankolcu's
  03438c81, whose grouped ledger this package keeps (v6 repairs its old-X serving and global listings);
  Galbraith, Wang and Zhang for X ± Y from one inverse; Arratia,
  Goldstein and Gordon for Chen–Stein. New here: next-array hash chains with read cursors, uncapped walks,
  history-conditional KB(129) bounds and a defective-run allowance (2.3, 4, 5, 6.6, 7), per-group fresh hop-1
  slots with cumulative hop-2 slots (4, 5), retained X affine points, and hchg/zero4096; the run and `certificate.py` re-verify the code.
