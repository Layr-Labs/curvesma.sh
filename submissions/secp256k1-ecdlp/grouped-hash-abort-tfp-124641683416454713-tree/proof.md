# Stored x-only leaf birthday search on secp256k1: six leaves per inverse labelled by truncated fraction products, grouped hash-chain abort, hoisted emit base, key-path tree (124.641683416454713)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; its expected total work is at most
`time_log2 = 124.641683416454713 = 124641683416454713/10^15` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit), and on every run at most R_worst = 2^188.505443 (disclosed: all 114 groups run, every hash walk
visits its whole chain of at most HOME_MAX links). Peak memory is at most 2^135.9989 bytes; preprocessing, at most
2^118.0361 units, is included; there is no advice.
`certificate.py` checks the declared figures (`DECLARED`): work, worst case, memory, preprocessing and work per leaf
(240.1834) are the least multiples of 10^-15, 10^-6, 10^-4, 10^-4 and 10^-4 at or above their exact values; 0.390000001 is
below the success bound.

**Reviewer's guide.** **New vs parent 38df975c** (124.645118129451422): x3 << 30 hoisted per sign, path tree (3.2): emit 7 → 6,
`keys` 140.5 → 138.5, JPATH 4 → 3 + a B per sign, same records; else the parent's. Run keys: 3.1 `fraction_mismatches` (a bit-exact model of `bwd`),
`fractions_outside_bounds`, `identity_mismatches`; 3.2 the key-pass check, `keys_formula_holds`, `loop_ends_on_bit_34`;
3.3 the checks "3.3's fraction ranges" and "3.3's error table", `tfp_terms`, `class_rule`, `class_labels`,
`window_batch`, `boundary_batch`; 3.4 `site_checks` (with `t_shxt`), `sites_match_certificate`, `forced_fixups`;
`grouped_paths` at `hchg`'s new counts (4). The run prints every certificate check with its value (`D_certificate.checks`). The bound
rests on:
- **Early abort with E-priced group costs (2.3, 6.6, 7).** Group g's work runs iff T > g - 1; its walk and key-path work
  has mean bounded given the history H_g, over which the stopping indicator is fixed, so E[work] = E_pre + sum_g PHI[g] W_g
  + global rows + worst-case rest (no independence claimed); PHI[g] is the Chen-Stein prefix tail plus prefix defect rows.
  Stopping needs a resolver-confirmed x ([x]P = Q), which latches; E1-E7 give the ungrouped run's success event.
- **Order-free hash chains (4).** `hchg` links each new record to its K2 head through the packed next-array, walks to the
  first full match or the chain end (uncapped), verifies d4 and lists: every colliding pair with no intruder is listed.
- **Six leaves per inverse, never formed (3.1).** One batched inverse of `D = y_Y^2 - y_X^2` gives u, v; with w = u² or v²
  the leaves `(Ew - b) + β^t(Fw - a) + β^(2t)Gw` are a Fourier transform; `bwd` stores two 128-bit fractions per component.
- **Fractions by truncated products (3.3).** 16 single-word products on tabled rows give `T(cM mod p)` less at most 64;
  a leaf's fractions are sums of its components' within [−132, 194], [−132, 194], [−389, 264]; W and D' are charged
  missed collisions (below 2^-58.47 and 2^-16.65 per leaf); a Frobenius lattice count bounds label boxes at 1 + 2^-17.
- **Rare carries by branch, metered (3.4):** results bit-identical to the inline chains; fixups cost 3 or 4, the meter's
  cap `2^45 - 33q2` is charged as `4*2^45` per node and its failure by Markov (below 2^-16.32 per leaf).
- **No uncharged transfer or bookkeeping (2.1, 4):** emit 6, `h2s` 7.25, `hchg` 6.25 plus 5.5 per expected visit, blocks
  as executed `dlvk` plus in-flight counts, kept or dropped whole; links keep no state.
- **Success and work** are evaluated with exact Fractions (grouped bound 0.3900000031) and every kernel count is checked
  against the certificate (`counts_match_certificate`). **Not relied on:** any heuristic, prior acceptance or score, or
  extrapolation from the 4 × 4 run.

## 1. Overview

Each node computes the leaves `x(φ^t(X) ± Y)`, t = 0, 1, 2, six per pair from one batched inverse (Galbraith, Wang and
Zhang, ePrint 2015/605, share one inverse between X ± Y; φ(x, y) = (βx, y) extends it to X's three twists), held as
Fourier components and labelled by truncated fraction products; records go in two bulk-synchronous hops to a home node,
where equal labels give the logarithm. An analytic bound, not a practical break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each a scalar AArch64 processor (base set with FEAT_RNG
`MRS Xt, RNDR`, section 7, and FEAT_LSE `LDADD`, `STADD`, `SWP`, sections 3.4 and 4) with a private memory below 2^64
bytes. The ledger counts executed instructions (each one step, whatever memory it touches) and, per message, its
delivery and in-flight count; storage and wiring carry no price.
- **Ports and links.** A node has a column port and a row port, device registers in its address space. A port's
  48-byte **message block** sends one record: a destination word at +0 (bits 3..38 the destination index, a row of its
  column or a column of its row, its own included; other bits ignored), an address at +8 and four data words at
  +16..+47; the STP to +32 completes it. Its **record blocks**, one per destination d at `PB + 2^14 d` (PB = 2^50 column,
  2^51 row), hold a header (destination word, address, byte count) at +0 and k = 256 slots, slot n at +32n; the header's
  STP completes a k-record message, its count word (32 per record, 8192 full) stored just before. The link writes a
  message's data (four words, or the k records in slot order) at that address of the destination and reports the write
  to its port, except hop-1 column blocks, which the destination's `h2s` reads from the link's window (section 4). It
  keeps no flags, rings or counters, only the message it carries. A store completing a message stalls until the link has
  read it, and a stalled instruction counts once: a port holds at most one message.
- **Delivery is charged** as executed `route.s` code: `dlv`, 6 per one-record message (three LDPs, two STPs, RET),
  `dlvk`, `4k + 2 = 1026` per k-record row block (header LDP, per record two LDPs and two post-indexed STPs, RET; the
  count word unread), or hop-1 arrival `h2s` (section 4). Run by the processor as a loop with B for RET (`dlvk` returns
  x0 to its start), its first LDP waits until x0, the link's window, holds a complete message; so this is all of a
  delivery's work, with no call or poll, and the destination index selecting the link's output is wiring.
- **Completion.** Each port counts its messages in flight (up when the link takes one, down when the destination's
  write completes), charged 2 per message (sections 4, 7). DSB completes once the node's stores are written and its
  ports' counts are zero; a node executes DSB before joining a barrier, so after a barrier every earlier message is in
  memory.
- **Barriers.** A combining tree of epoch writes (children of `(i, j)`: `(2i+1, j)`, `(2i+2, j)` and, for i = 0,
  `(0, 2j+1)`, `(0, 2j+2)`): a node waits for its children's arrival words, executes DSB, writes its parent's, waits for
  its release word and writes its children's. A wait is LDAXR, CMP, B.EQ and WFE; a store into the word clears the
  exclusive monitor and wakes the WFE, each word has its own reservation granule, and nothing else wakes it.

### 2.2 Leaves and records

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145307571` and `q2
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
1-5. **Groups** g = 1..114, each a contiguous range of X's, all nodes in lockstep, run steps 1-5 below after generating
  their X's; the group loop first tests the latched stop flag. Barriers separate every phase of every group.

1. **Hop 1:** under fresh per-group caps of cap_pg = cap_p/8 records (column LC/RC reset by `init3`/`tinit`), the batch
   kernels run for each X_k; the emit writes each record into d1's column block, and after slot k the block goes to
   `(d1, j)` at its RC cursor, or is dropped if that cursor is at or past the slot end (section 4). Then `fin` sends each row d's final block by the same test
   and writes its final cursor, clamped at the slot end (drops advance cursors too), into entry i of d's table E1.
2. **Hop 2:** `h2s` on each arriving hop-1 block writes each record through the row port to `(own row, d2)` by the same
   rule into the cumulative R2 slots of cap_p (row cursors persist); `fin` sends the final blocks and writes
   each home's final cursor into its `E2[j]`.
3. **Grouping:** `hchg` reads each R2 slot from its read cursor to E2, stores E2 back, and hash-joins the new records
   into the cumulative table (section 4): per-group [cursor, E2) ranges partition each slot, so every record is read and
   inserted once (list cap 2^28); the next-array covers every receive pattern up to the slot caps; no buckets, passes or scans.
4. **Resolution,** in four sub-phases: a home `(h_r, h_c)` sends two requests per new candidate to the
  sources `(g_r, g_c)` via the relay `(g_r, h_c)`, and replies of at most 8 records (one leaf's scalars and
  points) return via `(h_r, g_c)`.
  Each sub-phase is the base count exchange (counts, bases in sender order, drop marks at 2^29 requests at a source
  and 2^34 records at a relay, three barriers). The source serves a request from its tables — X_k's affine point
  from the retained all-groups X affine store (section 5), Y_l from its batch entry's Y words, and the scalars — copied into the
  reply (under 100 instructions of address forms, loads and stores, inside 2^10); old groups serve identically.
  As w3 gives s but not t, the home forms each record's three possible
  leaves `L = φ^t(X_k) + s*Y_l = [m]P + [n]Q`, t = 0, 1, 2, tests the nine pairs in a fixed order for a unit u with
  `L' = uL`, and with the first it finds that has `n' - u*n ≠ 0` outputs `x = (u*m - m')/(n' - u*n) mod r` after
  checking `[x]P = Q`.
5. **Stop test:** if any home resolved an x this group, the stop flag (latched at a fixed address, set by
  the resolving home) broadcasts over 2 hops and every node halts after the barrier; the latched x is the output. The
  test is 5 instructions per group per node (LDR, CMP, B, group-loop ADD/STR), charged in W_g.

**Barriers:** at most 1943 (pre 4, per group at most 17, output 1), within the allowance of 2048. **Caps:** arithmetic fixups (3.4),
pair slots (cap_pg per group at hop 1, cap_p cumulative at hop 2), hash walks (uncapped, at most
HOME_MAX visits per record since a chain holds one home's records of one fold), candidates, requests, relay records
and coin retries are capped or bounded in total, so the worst-case ledger holds on every run and the E-ledger (7)
averages only group-granular executed work.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's (q2 a multiple of 8: `fwd` takes eight incidences per pass, `bwd` and `keys` two). With `D_l = y(Y_l)^2 - y(X)^2`, the
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
F5: a square and its fold); βu² and βv² (PAIRB: products on the tabled rows `β*2^(64k) mod p`, folded; limb pairs 0..3 in (x29,x30), (x27,x28), (x7,x8), (x1,x20)); E, F and G,
each a product of a Y word (b, b³ or 3b²) on the four tabled rows `c*2^(64k) mod p`, k = 0..3, of an X constant c
(2a³, 2a or a²), plus b⁴ or a⁴ for E and F (AD4), folded (FDX); then for each of E, F and G its rows (SHX, SHX, then
SHXT, which forms only the words 1..3 of its row that the products read) and four **truncated products** (TR0, TRK;
3.3) with the sixteen limbs of u², v², βu² and βv², all in registers, finalized, with the correction words for α' and β'', and stored. 66043fba's three-pair β cache (section 9) is widened to four: the Y words in the batch entry free x1, and CP and TFP4 fit 13 working registers; the run checks every stored fraction bit for bit against the model of `bwd` (8). The header
H (bwd's caller's sp) holds the rows of 2a³, 2a and a², y_X, a⁴, `-T(a)` and `-T(βa)`; the rows of β lie at H − 144,
below bwd's 16-byte frame, written once per node; batch entry l (512 bytes) holds Y_l's words `y², y, b, b³, 3b², b⁴`
at +192, +224, +352, +384, +416, +448 and its correction words `-T(b)` and `-T(βb)` at +160 (the Y table, never
overwritten), D_l (+0) and acc_(l-1) (+256) from fwd (two STPs post-indexed by 256 step it), R' staged at +0 (bwd's
last LDP steps −512), and the fractions (V0, V1) of α'−, β''−, γ−, γ+, β''+ and α'+ at +64, +96,
+128, +256, +288 and +320 (γ's V1 as (s, hi), its low word 8s).
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
LDADD advances the entry by 32 and returns the record's slot; w3 = d2' + (x3 << 30) is one ADD x18, W, LSR #28 (j = 0, 2)
or EOR x14, S, LSR #28 (j = 1, S the low word of W0 + W1), with x18 = x3 << 30 and x14 = x18 | (2^36 - 1) set once per
sign (LSL, ORR); x3 is a multiple of 64, so x18's low 36 bits are 0 and the EOR equals the old EOR-then-ADD. Two STPs
and a TBZ of bit 13 (slot k, section 4): 6 instructions. After slot k, FLUSH (9, or 7 for a dropped block) ends the block (section 4). The entry pointer x3, moved
by one post-indexed LDP per sign (−192, then −320), is entry + 128 at a sign +'s emits and entry − 192 at a sign −'s;
entries are 512 bytes at a 2^35-aligned base, and `x3 << 30` mod 2^64, which the emit adds to d2, is `2<<36 | l<<39` or
`5<<36 | (l - 1 mod 2^25)<<39`. JH and JT run on into copies of the sign's later leaves (C2, C3)
that end in one B to the sign's end. **Key paths.** The exact index is uniform on {0, 1, 2} for a leaf other
than O, and a leaf takes JH or JT only if its exact index is 1 or 2 or it lies in W (3.3): an event of its value, of
probability at most `2/3 + 2d_wj`. Their work (3 more per such leaf; the B, charged on every sign) is E-priced in 6.6 by the same history-conditional
domination as the walks; there is no key-path meter, and the worst case charges every leaf (7). Executed counts per incidence, q2 a multiple of 8 (fixed part in parentheses): `fwd` 97.25 (9) and `bwd` 1281.5 (5), each
plus 3 or 4 per fixup (3.4), and `keys` 138.5 (4) for six leaves on the j = 0 path, including their six emits of 6, plus 3
per leaf on JH or JT, 1 per sign with one, and its flushes; fe_inv (`secp_fe.s`) is 22549. Per sign `keys` runs 6 LDPs, 28 sums, LSL, ORR
and per leaf ADD, CMP, B.HS, CMN, B.HS and the emit: 69, plus a TBZ per two incidences; a leaf on JH or
JT executes 3 more (ADDS, ADC, MVN of N'). Kernels may
overwrite x0-x29 and the flags and preserve sp; `bwd` and fe_inv save and restore x30, which `keys`, `fwd` and the
transport kernels do not write; `fwd` and `bwd` return their fixup word in x0 (3.4); the driver keeps no registers
across calls and sets sp to H for `bwd`. `bwd` keeps R in x25..x28 from
pass to pass; `fwd` takes eight incidences per pass (SUBS of 8, B.HS) and any rest two per pass (9 fixed), and `bwd` and `keys` end each pass with a TBZ on
bit 34 of the entry pointer: the 2^25 entries of 512 bytes fill the lower 2^34 bytes of their 2^35-aligned block, so
bit 34 is clear at every entry and set one entry below the block (both loop ends test bit 34 and the 156 emit TBZs test
bit 13: `loop_ends_on_bit_34`). In `keys` the j = 1, 2 blocks (JHB, JTB) and copies sit after the loop's RET.

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
- SHX's (and SHXT's) carry needs the row B it shifts to have its low word above `2^64 - C0`, as the word it adds is
  below C0: at most `C0*2^192` rows.
- FDX's carry into limb 2 leaves its output with low 128 bits below `F*C0 <= (2^66 + 1)*C0` and the rest nonzero, where
  `F = floor(Z/2^256) <= 2^66 + 1`: at most `F_w = (2^66 + 1)*2^128*C0` words.
- F5's carry into limb 2 leaves its output's low 128 bits below 2^65, as its fifth word is at most C0 + 1, so that word
  times C0 has a high word of at most 1: at most 2^193 words.
- FS's borrow leaves an output whose low word is at least `2^64 - C0`, and FADDW's carry one whose low word is below C0:
  at most `C0*2^192` words each.

While no site of the batch up to and including this one has deviated, the row or output is a word below 2^256 congruent
to the site's operand value (times `2^(64k)` for the k-th row), which therefore determines that value. So, given X and
every other generator, a site takes its fixup with probability at most its set's size times the points of `Y_l` per
operand value over r - 1, with the table's fresh generator and multiplicities. Summed over a batch (`FIXUPS` and
`E_FIX` in `certificate.py`), the expected number of fixups a batch takes before its first deviation, and so of a batch
without a defect, is at most `E_FIX`, and over a node's batches at most `q1*E_FIX < 2^28.68`.

**The meter row.** A node stops early only if N exceeds `2^45 - 33q2` at the start of some batch. On the event that none
of its batches has a defect, N is that clean count, so by Markov's inequality a node stops with probability at most
`q1*E_FIX/(2^45 - 33q2) + q1*d_batch < 2^-16.32` (`d_f`, a node-local row of 6.4).


## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (6), and per column block its flush (9); no `dlvk` (fused arrival) | 6.03515625 | 6.04296875 |
| arrival and hop 2: `h2s` on the block (7; loop 0.25; 11 per block), and per row block its flush (9; 9 if dropped) and `dlvk` (1026) | 11.3359375 | 11.34375 |
| hash join `hchg` fixed (6; loop 0.25) | 6.25 | 6.25 |
| hash walks (expected; 5 per visit, one B after an even walk: 5.5 per expected visit) | ~2.7 | E-rows |
| K2 verifications and listings (expected; global pair counting) | ~0.3 | E-rows |
| **fixed total** | **23.62109375** | **23.63671875** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlvk`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. A destination's full blocks number at most its records over k, and a
dropped block costs 7 or 9 and sends nothing, so 1037 per k records bounds each hop's row blocks (column blocks: flush 9
plus in-flight 2, no `dlvk`); final blocks are `fin`'s rows (section 7). `h2s` takes an arriving block's records, `hchg`
a slot's, in groups
of 1, 2 and 4, as the low bits of its count say, then 8 at a time; the 8-record loop's CMP and B.LO, 2 per pass, are
charged as 0.25 per record (n records run floor(n/8) passes), and the tests and set-up, 11 per block for `h2s`
and 13 per slot for `hchg` (which loads the slot's read cursor and stores E2 back), are per-block and per-node rows. **Addresses.** One instruction
forms each per-record table address: the emit's AND or BIC (3.2); `h2s`'s BFI of d2 (w3's low 36 bits) into bits
3..38 of x10, which holds 2^39 elsewhere (d2's row LC entry and destination word); the flushes' ORR of 2^43 (the RC
entry); `hchg`'s BFI of K2's low 57 bits into bits 3..59 of x10, which holds 3*2^61 elsewhere (the head's address,
8-byte heads at 3*2^61), and its ADD of C to A >> 1 (a record's next entry; section 5). **Blocks and cursors.**
The emit and `h2s` take a record's slot by one LDADD of 32 on its LC entry, which holds the absolute slot address:
block base (low 13 bits zero) plus offset 32n below 2^14, so bit 13 is set exactly after slot k (32k = 2^13). One TBZ
of bit 13 of LDADD's old value finds it (checked); there the **flush** returns
the LC entry to slot 1 (SUB, SUB, STR), forms the RC entry (ORR) and advances the RC cursor by 32k with LDADD, which returns
the old cursor c. If c is below the sender's slot end, which is the same at every destination and held in a register,
the STR of the byte count (32k) then the STP of the header (destination word, c) sends the block; otherwise it is dropped (CMP and B.HS, or B.LO in `h2s`),
and `h2s` stores the end back, so a hop-2 cursor stays at most its slot end; a hop-1 cursor can pass its end, but only
by the node's own leaves, 32 bytes each, so it stays below 2^63.2 (checked). `fin` takes the L = (LC - slot 1)/32
records left in each destination's block, sends it by the header (8d, c, 32L) if L > 0 and c is below the end, and writes
the final cursor `min(c + 32L, end)` (CMP, CSEL). **Keep lemma.** Cursors start at the slot base and move by 32k per
block, and `cap_p` is a multiple of k, so c and the end are congruent mod 32k, and c below the end gives `c + 32k <=
end`: a block is kept or dropped whole, and as its records take the cursors c, c + 32, ..., the records kept, their
slot order and the final cursor equal those of sending each record alone at its own cursor (the one-record rule). A
final block's k - L slots after its L records hold older records or none: at hop 1 the arrival handler stops at the
header count L and never reads them; at hop 2 they land in `[c + 32L, c + 32k)`, inside the
sender's slot and at or past the final cursor, where every later read stops. At hop 2 the row port's LC and RC
entries persist across groups and `fin` only reads them: the next group keeps filling the block from slot L + 1,
and its first full block to d rewrites `[c, c + 32L)` with the same L records (identical bytes, below that group's
read cursor) before the new ones at `[c + 32L, c + 32k)`; cursors stay block-aligned. The keep lemma holds at cap_pg too
(a multiple of k, checked), so the per-group hop-1 sender keeps or drops whole blocks the same way. **Fused arrival.** Hop-1
records skip R1 slots: each arriving block runs `h2s` (7 per record, the H2R body: 2 LDP, BFI, LDADD, 2 STP, TBZ),
routing as the slot loop did (same row blocks byte-for-byte: part B checks `port.blocks == snd.msgs`). R1/E1/RC are
retained (same layout and `fin` writes, still charged) as a conservative control; only their reads are gone. The
sender-side cap test stays, with its drops in the per-group row (`d_p1`, 6.6). **Hash table.** Each node
owns 2^57 head words at 3*2^61, zeroed once by `zero4096` (2048 STPs per pass, 2050 per 4096 words, plus RET), and a
packed next-array at 2^61 with one 16-byte entry (next address, K2) per R2 slot position (section 5). Bucket `K2 mod
2^57`: same K2 always lands in the same bucket, so same-key records meet with no ordering step. **Join.** `hchg` sets x10 to 3*2^61
once (MOV), steps its record pointer with a pre-indexed LDP and reads each slot from its read cursor to its E2
cursor, then stores E2 back as the new cursor (LDR of E2, LDR of the read cursor, STR of E2: the per-slot 13); the
per-group ranges partition each slot, so every record is read exactly once across all group calls, however the
slots' previous ends differ, and cross-group state otherwise travels only through the cumulative table. Per record
(6): LDP; BFI of the bucket into x10 (the head's address); SWP of A into the head, returning the old head; ADD of C to
A >> 1 (the new record's entry); STP of the old head and K2 into it; CBZ to the end when the bucket was empty (the walk
reads no head, so making A the head first changes no load). The walk is unrolled by two, the prior's address alternating
between x11 and x21: a visit forms the prior's entry address (ADD of C to its address >> 1), loads its next address
into the other register and its K2 with one LDP, then CMP, B.EQ and CBZ or CBNZ (5); after an even number of visits
a chain end takes one B to the end (at most floor(v/2) per record: 5.5 per visit bounds a walk). A K2 hit verifies d4 (LDR, EOR, TST, B.NE, 4;
a mismatch rejoins the walk); a full hit lists (CMP, B.HS, STP, B: 4; 3 from the second half) unless the list is
full (CMP, B.HS, 2); the matched visit skips its CBZ or CBNZ (the run counts 4, the ledger 5). Uncapped, so every
new record links to its most recent same-key predecessor on every run. The packed next-array (one 16-byte entry per R2
slot position, at (A>>1)+C, C = 2^61 - R2/2 a call argument) fits every receive pattern up to the slot caps by
construction; entries need no zeroing (written at insert before any read); records stay intact in R2 for phase 4. A hop-1
receiving node executes `h2s` on each arriving column block (fused arrival, above); a hop-2 receiving node executes
nothing on arrival and reads the records in the next phase with the LDPs of `hchg`. Per-node and per-group costs are rows of section 7, and every address a kernel loads or stores
lies in the regions of section 5. (`bkt`, `pass4`, `scan_ranges` stay for part B's base-path evidence; the
grouped run never executes them.)

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
| below 2^52 (2^51.3297 used) | code, stack, the five barrier words (one reservation granule each), header H with the rows of β below it, X table of the running group (512 bytes per X), X affine store (64 bytes per X, all groups), scalars, batch entries with the Y words (512 bytes per Y: 2^34 bytes at a multiple of 2^35), the ports' message blocks (at 2^48 - 64 and 2^48 - 128), their LC tables at 0 and 2^39 and RC tables at 2^43 and 2^43 + 2^39 (2^39 bytes each), the hash-read cursors at 2^48 (2^39 bytes), the record blocks at 2^50 and 2^51 (2^50 bytes each), E1 and E2 (2^41 each) and the phase-4 buffers and tables (`FIXED` in `certificate.py`, whose d3/d4 tables stay reserved but unused by the grouped run); all but the aligned tables lie between 2^44 and 2^48 | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | per-group hop-1 slots of cap_pg at 2^52 (`2^59.9986`, allocated but never written under fused arrival: a conservative control, section 4), the hash next-array at 2^61 (`2^61.9986`, R2 address A's entry at (A>>1) + C), the hash heads at 3*2^61 (`2^60`: a head's address is 3*2^61 with the bucket in bits 3..59); disjoint, inside R1 (checked), together 2^62.5838 | per group |
| R2 after R1, ending at 2^63.998943 | cumulative hop-2 pair slots (phase 4 reads candidates' records at their addresses) | groups 1-32 |

So every node stays below 2^64 bytes, and the peak is `2^72*(fixed + |R1| + |R2|)`, at most 2^135.9989
bytes (log2 ≈ 135.9988120). A remote item is reached only by a message.

## 6. Success (no heuristic)

In `success.md`: the model (6.1), Chen–Stein on cross-node pairs (6.2), usefulness (6.3), the defect rows (6.4), the
bound (6.5) and the grouped hash-chain abort (6.6: stopping-time bound, chain detection E1-E7, uncapped walks with
history-conditional walk means and a defective-run allowance).

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 645,410 per point (645,414 for a grouped X, whose four affine-store STPs, 64 bytes, the X row charges):** a complete Jacobian
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
  itemization as in aburan28's 7e8d13c7 and ecc3ec1b). A Y (no mixed addition, no X header) costs 642,930
  (`GEN_MULY`, checked; the idea of pepedesigner's p256 1772672b); the base ledger R keeps 645,410.
- **Node set-up, 2^20:** a program of at most 2^16 bytes (13,022 kernel instructions expanded, generator construction,
  resolver, phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD of 8 to its destination word, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 1,557,745:** two requests and two replies of at most 8 records over 2 hops at 40; each record's three
  possible leaves `φ^t(X_k) + sY_l` (one inversion of `D = x_Y³ - x_X³`, the inverses `Q_t/D` of the three denominators,
  and three affine additions: 26,753); nine unit tests at 600; the solve mod r and `[x]P = Q`, on the same exact
  arithmetic; 778,872, charged at twice that plus 1 (as in 1fbd26e8).
- **Phase-4 count exchange, 3*2^39 per node and sub-phase (E-averaged over the groups that run; the base row keeps
  2^41):** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 88 per message (its count,
  the count message, the offset reply and the write, at most 22 each with in-flight count and delivery) for at most
  8*2^29 messages (2^29 requests of up to 8 records): below 2^40.5.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^16 per node:** 12 kernel calls per group at most 33 each (arguments and BL), the stop test and
  group loop (8), over 114 groups: below 2^16.
- **Hash-head zeroing, 2050 per 4096 words:** `zero4096` (2048 post-indexed STPs of (xzr, xzr), SUBS, B.NE; its
  macros only shorten the listing); the row's +1 is its RET.
- **Hash-read-cursor init, 4 per slot:** `init3` (STR, ADD, SUBS, B.NE) of each slot's base into its 8-byte entry at
  2^48; the row's +1 is its RET.
- **Rows of β 72** (sixteen words of MOVZ and three MOVK, eight STPs, once per node),
  **batch set-up 160** (X's header with its correction words copied into H by 30 LDPs and 30 STPs, y_X² into
  registers, acc = 1, pointers, sp for `bwd`, the index word, the slot end for `keys`, the fixup meter's update and test
  and the loop: at most 145, itemized in `certificate.py`), **inversion call 16, request served 2^10, output 2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on
  the failure flag, plus a range test. RNDR is this machine's coin source: it returns 64 fresh uniform bits, the
  algorithm coins of the rules' probability space, and never reports failure, so the charged B.EQ is never taken and
  coins add no failure row.

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
| leaf work, E-averaged: sum_g PHI[g] (X gens, fwd, inv, bwd, keys w/o emits, batch set-up) (E) | 134.4274 |
| arithmetic fixups: 2^45*4 (W) | 119.0000 |
| j = 1 and j = 2 key paths, E-averaged: sum_g PHI[g] s_g (3 PJ + 1/3) + DEFJ (E; worst 128.4365) | 127.7022 |
| route work, E-averaged: sum_g PHI[g] (emits, arrival, row blocks, fin, column tables) (E) | 130.6459 |
| group work, E-averaged: sum_g PHI[g] (hchg fixed, main-pair walks, stop test) (E) | 129.7031 |
| same-X walks: global pair counting (E) | 29.3289 |
| K2 verifications: global pair counting (E) | 124.5689 |
| listings: global pair counting (E) | 75.9840 |
| defective-run rows: DEFV + DEFK2 + DEFL + DEFJ (E) | 63.1800 |
| candidates: 2^28*1491011 (W) | 120.5079 |
| requests served: 2^29*2^10 (W) | 111.0000 |
| phase-4 count exchanges, E-averaged: sum_g PHI[g] 4 * 3*2^39 (E; worst case 114 * 4 * 3*2^39 = 2^121.4179) | 121.1981 |
| barriers: 2048*2^82 (W) | 93.0000 |
| phase control: 2^16 (W) | 88.0000 |
| output + stop broadcast (W) | 85.1548 |
| **expected total E** | **134.607468** |

`E` is the exact `E_total` of `certificate.py` (expected instructions), at most 240.1834 per leaf (`per_leaf`), and `W_E = E/1000` with log2 W_E ≈ 124.6416834165. **The claim, 124.641683416454713, is the least multiple
of 10^-15 with 2^claim >= W_E** (the base claim stays on 10^-12); the certificate encloses 2^c by an exact series for
e^(f*ln 2) with interval bounds on ln 2, both ways (2^(claim - 10^-15) < W_E), from an 80-digit decimal start.
Rows marked P are preprocessing: their sum over 1000 has log2 ≈ 118.0360707 and is declared as 2^118.0361,
the least multiple of 10^-4 at or above it (checked). The worst case R_worst = 2^188.505443 (all groups run, uncapped
walks to HOME_MAX visits with verification and listing on every record) is disclosed alongside, not scored; it is the
least multiple of 10^-6 at or above log2 R_worst (checked).

## 8. Evidence

`experiment.py` takes `{"seed": 1, "seconds": 12}` (`seconds` 5 to 15); budget-limited checks run guaranteed minimums
(2,000 random and 600 near-tie class-rule inputs, 3,000 class-label classes), then their budget share up to 26,000 and 12,000.
**Part A** runs the leaf kernels and fe_inv on production batches of 2, 4, 6 and 10 incidences: every stored fraction
against a bit-exact model of `bwd` (`fraction_mismatches`) and its bound (`fractions_outside_bounds`), each product's
terms (`tfp_terms`), 3.1's identity against exact affine addition (`identity_mismatches`), each record, slot and LC
entry against the bit-exact model of `keys`, and each index and label outside W and D' (`wrong_index`,
`keys_not_canonical`). Random and near-tie leaves in random representations (`represent`) get the exact index outside W
(`class_rule`); on at least 3,000 classes, a third at an edge of W or D', every member emits its class's label unless in
W or D' (`class_labels`); `window_batch` and `boundary_batch` plant such leaves in batches. `fwd`/`bwd` counts match
their formulas plus 3 or 4 per fixup (`fwd_bwd_formulas_hold`); planted fixups (`fixup_batch`) and an always-inline copy
(`forced_fixups`: 33 fixups, 102 more per incidence) keep exact fractions. `site_checks` (`sites.py`) tests each site
type of 3.4 against its predicate, bit-level model and fixup; a `Y_l = φ(X)`, `Y_m = -X` batch gets acc = 0 and no
other fraction near its value; site macros match the tables (`sites_match_certificate`); flags are dead at every fixup
join. **Part B** (`machine.py`) runs a 4 × 4 machine through `leaf.s` and `route.s`: `init3`, `tinit`, `zero4`, `keys`
on planted exact components (labels with small digits at production widths, some classes shared), `fin`, fused `h2s`
arrival (byte count from the header, row blocks byte-exact), `bkt`, `prefix4`, `pass4`, `scan_ranges`, delivering by
executing `dlvk` or `dlv`. `block_paths` and `unrolled_paths` force the remaining block, cursor, bucket, pass and scan
paths (`h2s` on 1..19, 24, 44, 256-record and empty blocks); `two_group_finals` checks group 2 rewriting group 1's
partial final block. `grouped_paths` runs three groups through `zero4096`, `init3` and `hchg` on the true cumulative R2
layout with planted cross-group same-key pairs (stop, full and 1-entry-list cases): heads, next entries and candidates
equal the chain model, counts equal their formulas, cursors reach E2, every appended address is reached once, every
listed pair is same-key and every global pair is connected; old X entries persist and are served (`old_x_served`).
**Part D** runs `certificate.py` (`D_certificate.checks`) and `counts_match_certificate` compares its constants with A and B.

The 4 × 4 run is evidence at its own scope only.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1 (links keep no state beyond the message in flight; the sender's
  instructions, each delivery as `dlv`, `dlvk` or fused `h2s` arrival, and the in-flight counts are charged); a forwarding network would
  add work per extra hop and is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up
  with the driver's fixup-meter test, stop test and phase-4 exchanges), the in-flight counts, and the production attack,
  machine and links. The toy stop fires on a listed candidate; production stops only on a resolver-confirmed x, a
  stronger rule, so the E-bound (which charges confirmed stopping) is conservative over the toy.
- **Provenance**, not a premise: kaankolcu's 38df975c (the declared parent; this is it plus `changes.md` 1-2) and ob2 (the grouped truncated-fraction line: 2.1-2.2, resolver, success rows,
  Fourier leaves, k=256/TBZ machine, record blocks, label, address forms, fe_inv, scan_ranges, candidate, generator,
  meter, E-ledger, hash chains, grouped ledger, old-X serving, 10^-15 grid); Galbraith, Wang and Zhang for X ± Y from
  one inverse; Arratia, Goldstein and Gordon for Chen–Stein. The parent's change vs ob2: `changes.md`. The run and
  `certificate.py` re-verify the code.
