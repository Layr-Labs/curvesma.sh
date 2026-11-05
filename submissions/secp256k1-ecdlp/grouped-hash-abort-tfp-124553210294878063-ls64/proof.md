# Stored x-only leaf birthday search on secp256k1, six leaves per inverse, truncated fraction labels, grouped hash-chain abort (124.553210294878063)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; its expected total work is at most
`time_log2 = 124.553210294878063 = 124553210294878063/10^15` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit), and on every run at most R_worst = 2^188.604975 (disclosed: all 114 groups run, every hash walk
visits its whole chain of at most HOME_MAX links). Peak memory is at most 2^135.9990 bytes; preprocessing, at most
2^117.0744 units, is included; there is no advice.
`certificate.py` checks the declared figures (`DECLARED`): work, worst case, memory, preprocessing and work per leaf
(225.8939) are the least multiples of 10^-15, 10^-6, 10^-4, 10^-4 and 10^-4 at or above their exact values; 0.390000001 is
below the success bound.

**Reviewer's guide.** Provenance and the changes against the parent 38df975c: section 9. Run keys: 3.1 `fraction_mismatches` (a bit-exact model of `bwd`),
`fractions_outside_bounds`, `identity_mismatches`; 3.2 the key-pass check, `keys_formula_holds`, `loop_ends_on_bit_34`;
3.3 the checks "3.3's fraction ranges" and "3.3's error table", `tfp_terms`, `class_rule`, `class_labels`,
`window_batch`, `boundary_batch`; 3.4 `site_checks` (with `t_shxt`), `sites_match_certificate`, `forced_fixups`;
`grouped_paths` on `dlvh` (4); `ls64` (2.1). The run prints each certificate check and value (`D_certificate.checks`). The bound rests
on:
- **Early abort with E-priced group costs (2.3, 6.6, 7).** Group g's work runs iff T > g - 1; its walk and key-path work has mean bounded given the history H_g, which fixes the stopping indicator, so E[work] = E_pre + sum_g PHI[g] W_g + global rows + worst-case rest (no independence claimed); PHI[g] is the Chen-Stein prefix tail plus prefix defect rows. Stopping needs a latched resolver-confirmed x ([x]P = Q); E1-E7 carry 6.1-6.5's success argument to the grouped run.
- **Order-free hash chains (4).** `dlvh` links each new record to its K2 head via packed next-array, walks to the first full match or chain end (uncapped), verifies d4 and lists every colliding pair with no intruder.
- **Six leaves per inverse, never formed (3.1).** One batched inverse of `D = y_Y^2 - y_X^2` gives u, v; with w = u² or v² the leaves `(Ew - b) + β^t(Fw - a) + β^(2t)Gw` are a Fourier transform; `bwd` stores two 128-bit fractions per component.
- **Fractions by truncated products (3.3).** 16 single-word products on tabled rows give `T(cM mod p)` less ≤ 64; a leaf's fractions are sums of its components' within [−132, 194], [−132, 194], [−389, 264]; W and D' are charged missed collisions (< 2^-58.47 and 2^-16.65 per leaf); a Frobenius lattice count bounds label boxes at 1 + 2^-17.
- **Rare carries by branch, metered (3.4):** results bit-identical to inline chains; fixups cost 3, 4 or 5, meter's cap `2^45 - 33q2` is charged as `5*2^45` per node and its failure by Markov (< 2^-16.32 per leaf).
- **No uncharged transfer or bookkeeping (2.1, 4):** emit 6, `h2s` 6 (LS64 pairs), join 4 plus 4.5 per expected visit, blocks as executed `dlvh` plus in-flight counts, kept or dropped whole; links keep no state.
- **Success and work** use exact Fractions (grouped bound 0.3900000046); every kernel count is checked against certificate (`counts_match_certificate`). **Not relied on:** any heuristic, prior acceptance or score, or extrapolation from the 4 × 4 run.

## 1. Overview

Each node computes leaves `x(φ^t(X) ± Y)`, t = 0, 1, 2, six per pair from one batched inverse (Galbraith, Wang and
Zhang, ePrint 2015/605, share one inverse between X ± Y; φ(x, y) = (βx, y) extends it to X's three twists), held as
Fourier components and labelled by truncated fraction products; records go in two bulk-synchronous hops to a home node,
where equal labels give the logarithm. An analytic bound, not a practical break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each with a main processor and one delivery processor per kind of message (`dlv`, `dlvk`, hop-1 arrival `h2s`, hop-2 arrival `dlvh`), all scalar AArch64 (base set with FEAT_RNG
`MRS Xt, RNDR`, section 7, and FEAT_LSE `LDADD`, `STADD`, `SWP`, sections 3.4 and 4; base-ISA `DC ZVA` with DCZID_EL0.BS = 4, a 64-byte block, DZP = 0 and SCTLR_EL1.DZE = 1 so EL0 may use it, on cacheable Normal memory, section 7; FEAT_LS64 `LD64B`/`ST64B` with SCTLR_EL1.EnALS = 1, section 4) with private memory < 2^64
bytes (a flat address space: a property of this machine, as in 38df975c's 2.1). The ledger counts every processor's executed instructions (one step each, whatever memory they touch) and each message's delivery
and in-flight count; storage and wiring are free.
- **Ports and links.** A node's column and row ports are device registers in its address space. A port's 48-byte **message block** sends one record: a destination word at +0 (bits 3..38 the destination index, a row of its column or a column of its row, own included; other bits ignored), an address at +8 and four data words at +16..+47; the STP to +32 completes it. Its **record blocks**, one per destination d at `PB + 2^14 d` (PB = 2^50 column, 2^51 row), hold a header (destination word, address, byte count; at hop 2, p at +24) at +0 and k = 256 slots, slot n at +32n; header's STP, after count word (32 per record, 8192 full), completes a k-record message. The destination's delivery code writes a message's data (four words, or the k records in slot order) at that address, completing the write reported to the sender's port; hop-1 column blocks are instead consumed by the destination's `h2s` from the link's window (section 4). It keeps no flags, rings or counters, only the message it carries. A store completing a message stalls until the link reads it (counted once): a port holds at most one message and accepts one completing store at a time from any of its node's processors, the others stalling.
- **LS64 (premise).** `LD64B`/`ST64B Xt, [Xn]` move Xt..Xt+7 in one 64-byte access: t even, ≤ 22; no offset or writeback; Xn 64-byte aligned (DDI0487; KA006230). Supported locations are IMPLEMENTATION DEFINED: here the link windows (Device; in each delivery processor's own mapping) and R2 (Normal Non-cacheable), nowhere else (`a64.py` asserts all of this). Types change no count.
- **Delivery is charged** as executed `route.s` code: `dlv`, 6 per one-record message (three LDPs, two STPs, RET), `dlvk`, `4k + 2 = 1026` per k-record row block (header LDP, per record two LDPs and two post-indexed STPs, RET; count word unread), hop-1 arrival `h2s` or hop-2 arrival `dlvh` (section 4). Run by its kind's delivery processor as a loop with B for RET (`dlvk` returns x0 to its start), its first load of the window (x0; x1 for `h2s`) waits until the window holds a complete message: this is a delivery's whole work, with no call or poll (the link presents a message in that window; only delivery code writes memory); the destination index selecting link's output is wiring. `h2s`'s processor runs `h2s0` once at node set-up (node set-up first loads x12, the node's hop-2 slot end; x10, x18, x12 and the p-word address x3 then persist; BFI rewrites only bits 3..38 of x10); each run resets x1 to the window, waits on it (the count LDR stalls the processor until the message is complete, ordering every later load after it) and only then loads p from the node word at 2^48 − 168, which the main processor writes for group g (CSEL, CSEL, ADD, STR; phase control) at set-up for g = 1 and otherwise after group g − 1's hop-1/hop-2 barrier and before its DSB at group g − 1's last barrier; a group-g block arrives only after that barrier, so its run loads p_g.
- **Completion.** Each port counts its messages in flight (up when the link takes one, down when destination's write completes), charged 2 per message (sections 4, 7). A hop-1 block's write completes when its `h2s` run returns, a row block's when its `dlvh` run returns. DSB completes once all of the node's processors' stores are written and its ports' counts are zero; a node executes DSB before joining a barrier, so afterwards every earlier message of the main processors is written and every hop-1 block consumed (hence `fin` follows the last `h2s` on the row port); the row blocks those `h2s` runs sent complete by the next barrier. **No deadlock:** a main processor's completing store waits only for a destination's delivery processor; `h2s`'s row-port stores wait only for `dlvh`, which writes memory and never waits on a port; this wait graph is acyclic.
- **Barriers.** A combining tree of epoch writes (children of `(i, j)`: `(2i+1, j)`, `(2i+2, j)` and, for i = 0, `(0, 2j+1)`, `(0, 2j+2)`): a node waits for its children's arrival words, executes DSB, writes its parent's, waits for its release word and writes its children's. A wait is LDAXR, CMP, B.EQ and WFE; a store into word clears exclusive monitor and wakes the WFE (nothing else does: a property of this machine's event register, which raises no spurious events), and each word has its own reservation granule.

### 2.2 Leaves and records

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145307202` and `q2
= 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the x-coordinates
of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and `S = 2^72*n0 = 2^126.699493`
in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and `n = λ^t*b_g`. A point Z ≠ O
has class `{±Z, ±φ(Z), ±φ²(Z)}`; its x-coordinates `{x, βx, β²x}` (mod p) have canonical member `v* = min(x, βx
mod p, β²x mod p)`, and its **label** is `(comp(v*), comp(βv* mod p))`, `comp(v) = floor(2^100 v/p)`: two 100-bit
parts shared by class (3.3). A record is `[K2, K3, w2, w3]` with K2 and K3 the top 64 bits of the first and
second part, `w2 = i | k<<36` and `w3 = d2 | c<<36 | l'<<39`, where c = 2 and l' = l for s = +1, and c = 4 and
`l' = l - 1 mod 2^25` for s = −1, bits 36..63 complemented (c = 5 or 3, `~l'`) on a j = 1 leaf (3.2): c in {2, 4, 5, 3}
decodes uniquely, and w3 gives the sign, not t. Its digits are d1 and d2, the parts' low 36 bits
(hop-1 row, hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3 mod 2^46`. At its home a record's source is `(w2 mod 2^36,
j)`, j its arrival slot's index (2.3).

### 2.3 Phases

A pair slot holds `cap_p = 2^22 - 2^12` records of 32 bytes, a multiple of k. At hop 1 a node in row i writes only into
slot i of R1, at `R1 + 32*cap_p*i`, at each node of its column; at hop 2 a node in column j writes only into slot j of
R2 at each node of its row. Each port has two tables of 2^36 entries, one per destination d: the **LC** entry is the
next slot of d's record block (`init3` sets it to slot 1), the **RC** entry holds sender's cursor in its pair slot
at d (`tinit` sets it to slot's base); the slot ends at the same address at every destination, so the sender keeps
that end in a register. The column port's LC table is at address 0 and row port's at 2^39, so d's entry, 8d or 2^39 +
8d, is also a destination word for d; the RC tables lie 2^43 above (ORR of 2^43).

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, Y generators with their table entries, rows of β, row port's LC/RC tables, hash-head zeroing (`zero4096` over 2^57 words), hash-read cursors (`init3` of each slot's base into its 8-byte entry at 2^48) and the β tables (`tbld`, 2^51 bytes at 2^64 − 2^51; its six slot deltas β*2^(43k) mod p, k = 0..5, are 192 bytes of program constants, its call inside node set-up). Then Q is broadcast in 2 hops and R_g generated.
1-5. **Groups** g = 1..114, each a contiguous range of X's, all nodes in lockstep, run steps 1-5 below after generating their X's, the group loop first testing the latched stop flag; barriers separate all phases of all groups.

1. **Hop 1:** under fresh per-group caps of cap_pg = cap_p/8 records (column LC/RC reset by `init3`/`tinit`), the batch kernels run per X_k; the emit writes each record into d1's column block, which after slot k goes to `(d1, j)` at its RC cursor, or is dropped if that cursor is at or past slot end (section 4). `fin` then sends each row d's final block by the same test and writes its final cursor, clamped at slot end (drops advance cursors too), into entry i of d's table E1.
2. **Hop 2:** `h2s` on each arriving hop-1 block writes each record through row port to `(own row, d2)` by the same rule into cumulative R2 slots of cap_p (row cursors persist); `fin` sends the final blocks and writes each home's final cursor into its `E2[j]` (two tables, by group parity).
3. **Grouping,** at arrival: `dlvh` copies each arriving row block into R2 and hash-joins its unjoined records into cumulative table (section 4): joining only records at or past slot's previous final cursor F joins each once (list cap 2^28); no sorting pass.
4. **Resolution,** in four sub-phases: a home `(h_r, h_c)` sends two requests per new candidate to sources `(g_r, g_c)` via relay `(g_r, h_c)`, and replies of ≤ 8 records (one leaf's scalars and points) return via `(h_r, g_c)`. Each sub-phase is a count exchange (counts, bases in sender order, drop marks at 2^29 requests at a source and 2^34 records at a relay, three barriers). The source serves a request from its tables (X_k's affine point from retained all-groups X affine store, section 5; Y_l from its batch entry's Y words; the scalars), copied into reply (under 100 instructions of address forms, loads and stores, inside 2^10); old groups serve identically. As w3 gives s but not t, the home forms each record's three possible leaves `L = φ^t(X_k) + s*Y_l = [m]P + [n]Q`, t = 0, 1, 2, tests the nine pairs in a fixed order for a unit u with `L' = uL`, and with the first having `n' - u*n ≠ 0` outputs `x = (u*m - m')/(n' - u*n) mod r` after checking `[x]P = Q`.
5. **Stop test:** if a home resolved an x this group, the stop flag (latched at a fixed address by resolving home) broadcasts over 2 hops and all nodes halt after barrier, outputting the latched x. The test, 5 instructions per group per node (LDR, CMP, B, group-loop ADD/STR), is charged in W_g.

**Barriers:** ≤ 1943 (pre 4, per group ≤ 17, output 1), within the 2048 allowance. **Caps:** arithmetic fixups (3.4),
pair slots (cap_pg per group at hop 1, cap_p cumulative at hop 2), hash walks (uncapped, at most
HOME_MAX visits per record since a chain holds one home's records of one fold), candidates, requests, relay records
and coin retries are capped or bounded in total, so the worst-case ledger holds on every run and the E-ledger (7)
averages only group-granular executed work.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's (q2 a multiple of 32: `fwd` takes 32 incidences per pass, `bwd` four, `keys` two). With `D_l = y(Y_l)^2 - y(X)^2`, `fwd`
stores `D_l` and `acc_(l-1)` in entry l and forms `acc_l = acc_(l-1)*D_l` (`acc_(-1) = 1`); fe_inv gives
`R = 1/acc_(q2-1)`; `bwd`, for `l = q2-1..0`, sets `inv = R*acc_(l-1)`, `R <- R*D_l`, `u = inv*(y_Y - y_X)`
and `v = inv*(y_Y + y_X)`, and stores twelve **fractions**, two of each **component** `α' = Ew`, `β'' = Fw` and `γ = Gw`
(3.3), with `a = x_X`, `b = x_Y`, `E = b^4 + 2a^3 b`, `F = a^4 + 2ab^3`, `G = 3a^2 b^2` and w = u² (s = +1) or v²
(s = −1). For t = 0, 1, 2 the leaves, like the components never formed, are

    x(φ^t(X) + sY) = (α' - b) + β^t (β'' - a) + β^(2t) γ   (mod p).

**Why.** As R = 1/acc_l before step l, inv = 1/D_l, so `u = 1/(y_Y + y_X)` and `v = 1/(y_Y - y_X)`. As
`y² = x³ + 7`, `D = b³ - a³`, and as β³ = 1, `(b - β^t a)*Q_t = D` with `Q_t = b² + β^t ab + β^(2t) a²`. The slope of
`φ^t(X) ± Y` with `φ^t(X) = (β^t a, y_X)` is `(±y_Y - y_X)/(b - β^t a) = (±y_Y - y_X)*Q_t/D`, which is `Q_t*u` or
`-Q_t*v`; its square, less the two x-coordinates, is the leaf `Q_t²w - β^t a - b`, and expanding with `β^(3t) = 1` gives
`Q_t² = E + β^t F + β^(2t) G`, hence the display. D ≠ 0 unless `b³ = a³`, i.e. `Y = ±φ^t(X)` (3.4).

**Arithmetic per pair**, after W2S (`inv`, the new R) and W2L (u and v, on shared rows of inv): u² and v² (SQ13 and
F5: a square and its fold); βu² and βv² (PAIRB: sums of six table entries `β*t*2^(43i) mod p`, t the i-th 43-bit window of w, at 2^64 - 2^51 + 256t + 32i (written by `tbld`, 20 per entry), and a fold by C0 whose rare carry adds C0 (BLK0, 52 each, 105 per pair with the base MOV; the fold's branch tests the carry out of limb 1, and the rare fixup finishes limbs 2 and 3 and adds C0 on their carry). A window's entry address is one ORR of the base 0xfff8<<48 with the window's word shifted left 8: the base's bits 51-63 are all set, so bits above the window's 43 absorb into them and no AND is needed. Window 0 is limb 0 as it stands; windows 1, 2 and 4 take one EXTR first and windows 3 and 5 one LSR; limb pairs 0..3 in (x29,x30), (x27,x28), (x7,x8), (x1,x20)); E, F and G,
each a product of a Y word (b, b³ or 3b²) on the four tabled rows `c*2^(64k) mod p`, k = 0..3, of an X constant c
(2a³, 2a or a²), plus b⁴ or a⁴ for E and F (CPK: added right after the first row), folded (FDX); then for each of E, F and G its rows (SHX, SHX, then
SHXT, forming only words 1..3, those the products read) and four **truncated products** (TR0, TRK;
3.3) with sixteen limbs of u², v², βu² and βv², all in registers, finalized with correction words for α' and β'' and stored. 66043fba's three-pair β cache (section 9) is widened to four: the Y words in batch entry free x1, and CP and TFP4 fit 13 working registers; the run checks every stored fraction bit-exactly against the `bwd` model (8). The header
H (bwd's caller's sp) holds rows of 2a³, 2a and a², y_X, a⁴, `-T(a)` and `-T(βa)`; the rows of β (kept, no longer read) lie at H − 144,
below bwd's 16-byte frame, written once per node; batch entry l (512 bytes) holds Y_l's words `y², y, b, b³, 3b², b⁴`
at +192, +224, +352, +384, +416, +448 and its correction words `-T(b)` and `-T(βb)` at +160 (the Y table, never
overwritten), D_l (+0) and acc_(l-1) (+256) from fwd (two STPs post-indexed by 256 step it), R' staged at +0 (bwd's
last LDP steps −512), and fractions (V0, V1) of α'−, β''−, γ−, γ+, β''+ and α'+ at +64, +96,
+128, +256, +288 and +320 (γ's V0 and V1 as (s, hi), low words 8s).
The rare carries of SHX, SHXT, FDX, F5, FS and FADDW are taken by branches over fixups giving inline chains'
results (3.4); the only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks entries downwards, one sign at a time, loading its six fractions (6 LDPs) and forming its three leaves'
W0 and W1 by 13 two-word sums and differences (26 in all): with a, b and g the fractions of α, β' and γ (3.3),
they are `a0 + b0 + g0` and `a1 + b1 + g1` for t = 0 (4); then `i1 = b0 - g0`, `i2 = b1 - g1`, `i3 = a0 - b0`, `i5 = a1
- b1`, `i6 = i1 + i2` and for t = 1 `i3 + i6 = a0 + b1 - g0 - g1` and `i5 - i1 = a1 + g0 - b0 - b1` (7), and for t = 2
`i3 - i2 = a0 + g1 - b0 - b1` and `i5 + i6 = a1 + b0 - g0 - g1` (2): the same values mod 2^128 as the parent's 14, γ's V0 and
V1 entering shifted (LSL #3, by ADDS or SUBS of a shifted register). Each leaf's KEY
index j comes from top words F0 of W0 and F1 of W1: with `σ = F0 + F1`, if F0 ≥ F1 (CMP, B.HS) it branches to block
JH, where j = 2 if `F1 + σ` carries (CMN, B.HS to JT), else 1; otherwise, inline, j = 2 if `F0 + σ` carries (CMN, B.HS
to JT), else 0. The j = 0 path emits (W0, W1); JH forms `N' = ~(W0 + W1)` (ADDS, ADC, MVN) and emits (W1, N'); JT forms
N' and emits (N', W0). The **emit** forms 8*d1 by one AND (j = 0, 1) or BIC (j = 2) of x1 = (2^36 - 1) << 3 with low
word of pair's first fraction shifted right by 25 bits: with column LC table at address 0 this is both d1's LC
entry and destination word. LDADD advances entry by 32 and returns record's slot; w3 = d2' + (x3 << 30) is
one ADD x18, W, LSR #28 (j = 0, 2) or EON x18, S, LSR #28 (j = 1, S the low word of W0 + W1, whose complement is N''s
low word), with x18 = x3 << 30 set once per sign (LSL); x3 is a multiple of 64, so x18's low 36 bits are 0: the ADD
puts d2 below x18, and the EON puts N''s d2 below x18's complemented bits 36..63. Two STPs and a TBZ of bit 13 (slot k, section 4): 6 instructions. After slot k, FLUSH (9,
or 7 for a dropped block) ends block (section 4). The entry pointer x3, moved by one post-indexed LDP per sign
(−128, then −384), is entry + 128 at a sign +'s emits and entry − 256 at a sign −'s; entries are 512 bytes at a
2^35-aligned base, and `x3 << 30` mod 2^64, added by emit to d2, is `2<<36 | l<<39` or
`4<<36 | (l - 1 mod 2^25)<<39` (complemented on j = 1, where the codes map 2 -> 5 and 4 -> 3, so the four codes 2, 4, 5, 3 are distinct and give s and j's complement). JH and JT end in one B back to the j = 0 path's end. **Key paths.** The exact index is uniform on {0, 1, 2} for a leaf other than O, and a leaf takes JH or JT
only if its exact index is 1 or 2 or it lies in W (3.3): an event of its value, of probability ≤ `2/3 + 2d_wj`.
Their work (4 more per such leaf) is E-priced in 6.6 by walks' history-conditional
domination; no key-path meter exists, and the worst case charges every leaf (7). Executed counts per incidence, q2 a
multiple of 32 (fixed part in parentheses): `fwd` 97.0625 (9) and `bwd` 1230.25 (5), each plus 3, 4 or 5 per fixup (3.4), and
`keys` 132.5 (4) for six leaves on j = 0 path, including their six emits of 6, plus 4 per leaf on JH or JT, and its flushes; fe_inv (`secp_fe.s`) is 22549. Per sign `keys` runs 6 LDPs, 26 sums, LSL and per
leaf ADD, CMP, B.HS, CMN, B.HS and emit: 66, plus a TBZ per two incidences; a leaf on JH or JT executes 4 more
(ADDS, ADC, MVN of N', B back). Kernels may overwrite x0-x29 and flags and preserve sp; `bwd` and fe_inv save and restore
x30, which `keys`, `fwd` and transport kernels never write; `fwd` and `bwd` return their fixup word in x0 (3.4); the
driver keeps no registers across calls and sets sp to H for `bwd`. `bwd` keeps R in x25..x28 between passes; `fwd` takes
32 incidences per pass (SUBS of 32, B.HS) and any rest two per pass (9 fixed); `bwd` (four bodies, as aburan28's c3e0530e) and `keys` end each pass with a
TBZ on bit 34 of entry pointer: the 2^25 entries of 512 bytes fill lower 2^34 bytes of their 2^35-aligned block,
so bit 34 is clear at every entry and set one entry below it (both loop ends test bit 34 and the 36 emit TBZs test bit
13: `loop_ends_on_bit_34`). In `keys` the j = 1, 2 blocks (JHB, JTB) sit after loop's RET.

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
exactly (each high word is at most 2^64 - 2), CPK adds b⁴ or a⁴ (below p) for E and F to the first row's five words (`b_0*c_0 + K <= (2^64 - 1)(2^256 - 1) + 2^256 - 1 < 2^320` for any 256-bit row, no carry out; `cpk_exact` tests it at the extremes), `Z' < 2^322 + 2^256`, and a fold
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
inputs: RF and RA; CPK; FS and FADDW, whose operands (`y_Y²` and `y_X²` in `fwd`, `y_Y` and `y_X` in `bwd`) are
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
| bwd: window-sum folds of βu² and βv² (BLK0; exact, a carry adds C0 once) | 2q2 | 0 | 30*C0*2^128 | 6 |
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
SHXT, FDX, F5), 4 (FS, FADDW) or 5 (BLK0: ADCS, ADCS of limbs 2 and 3, CSEL of C0 on their carry, ADD, STADD; the CSEL makes a
forced run add 0) and C0 to the word at [sp], which opens a 16-byte frame pushed as zero by `fwd` (STR)
and by `bwd` beside its saved x30 (STP); each returns the word in x0, at most `33q2*C0 < 2^63` per call. The driver adds
x0/C0 to the node's fixup count N and starts a batch only while `N <= 2^45 - 33q2` (33 sites per incidence: 18 SHX, 3
of them SHXT, 8 FDX, 2 BLK0, 2 F5, 2 FS and 1 FADDW). So N is at most 2^45 on every run, and section 7 charges `5*2^45`
per node.

**How often.** A fixup needs a word in a small set:
- SHX's (and SHXT's) carry needs the row B it shifts to have its low word above `2^64 - C0`, as the word it adds is
  below C0: at most `C0*2^192` rows.
- FDX's carry into limb 2 leaves its output with low 128 bits below `F*C0 <= (2^66 + 1)*C0` and the rest nonzero, where
  `F = floor(Z/2^256) <= 2^66 + 1`: at most `F_w = (2^66 + 1)*2^128*C0` words.
- F5's carry into limb 2 leaves its output's low 128 bits below 2^65, as its fifth word is at most C0 + 1, so that word
  times C0 has a high word of at most 1: at most 2^193 words.
- FS's borrow leaves an output whose low word is at least `2^64 - C0`, and FADDW's carry one whose low word is below C0:
  at most `C0*2^192` words each.
- BLK0's branch is skipped (a fixup) only when the fold's carry leaves limb 1, which needs the window sum (below 6p,
  t = sum >> 256 <= 5) to have its low 128 bits >= `2^128 - t*C0`: at most 15*C0*2^128 sums, each fixing w mod p, so
  at most `30*C0*2^128` input words w (u² or v²); with the inline fixup the result equals the full carry chain. Its result is exact (no defect site).

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
| arrival and hop 2: `h2s` on the block (6; a full block 256 unrolled, 8 per block), and per row block its flush (9; 9 if dropped) and `dlvh`'s copy (531) | 8.140625 | 8.1484375 |
| hash join in `dlvh` (4) | 4 | 4 |
| hash walks (expected; ≤ 4.5 per visit) | ~2.5 | E-rows |
| K2 verifications and listings (expected; global pair counting) | ~0.3 | E-rows |
| **fixed total** | **18.17578125** | **18.19140625** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlvh`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. A destination's full blocks number at most its records over k, cumulatively (a group's prefix sums are dominated and PHI decreases, so the PHI-weighted rows dominate too); a dropped block
costs 7 or 9 and sends nothing, so 542 per k records bounds hop 2's row blocks (hop 1's column blocks: flush 9 plus
in-flight 2, no delivery code); final blocks are `fin`'s rows (section 7). **LS64 pairs.** Window slot n (base + 32n) and its R2 copy (c + 32(n - 1), c = 32 mod 64, section 5) are 64-byte
aligned exactly for even n: slots (2m, 2m + 1) move as one access, a lone slot (slot 1, an unpaired last) by two LDPs
(and STPs). `h2s` unrolls a full block (MOVN, LDR count, LDR p, MOV, TBNZ of bit 13, RET, lone slots 1 and 256: 8),
else slot 1, pairs 1, 2, then 4 at a time, a lone last slot (≤ 16 fixed, CMP and B.LO per 8-record pass; per node and
group ≤ SIDE + 1 such blocks, at 16 + 62); `dlvh`'s 15 per fresh full block are per-block rows (`dlvh`'s other blocks: the final rows, below). **Addresses.** One
instruction forms each per-record table address: emit's AND or BIC (3.2); `h2s`'s BFI of d2 (w3's low 36 bits) into
bits 3..38 of x10, else 2^39 (d2's row LC entry and destination word); the flushes' ORR of 2^43 (RC entry); join's
BFI of K2's low 57 bits into bits 3..59 of x10, else 3*2^61 (head address, 8-byte heads at 3*2^61); the entry address x13, below. **Blocks and cursors.** The emit and `h2s` take a record's slot by one LDADD
of 32 on its LC entry, the absolute slot address: block base (low 13 bits zero) plus offset 32n < 2^14, so bit 13 is
set exactly after slot k (32k = 2^13), found by one TBZ of bit 13 of LDADD's old value (checked). There the **flush**
resets the LC entry to slot 1 (SUB, SUB, STR), forms the RC entry (ORR) and advances the RC cursor by 32k with LDADD,
returning the old cursor c. If c is below sender's slot end (2.3), the STR of byte count (32k; at hop 2 an STP
with p, below) and header's STP (destination word, c) send block; else it is dropped (CMP and B.HS, or B.LO in
`h2s`) and `h2s` stores end back, so a hop-2 cursor stays at most its slot end; a hop-1 cursor passes its end only
by node's own leaves, 32 bytes each, staying < 2^63.2 (checked). `fin` takes the L = (LC - slot 1)/32 records
left in each destination's block, sends it by header (8d, c, 32L) if L > 0 and c is below end, and writes the
final cursor `min(c + 32L, end)` (CMP, CSEL). **Keep lemma.** Cursors start at slot base and move by 32k per block,
and `cap_p` is a multiple of k, so c and end are congruent mod 32k, and c below end gives `c + 32k <=
end`: a block is kept or dropped whole, and as its records take cursors c, c + 32, ..., the records kept, their slot
order and the final cursor match sending each record alone at its own cursor (the one-record rule). A final block's k -
L slots after its L records hold older records or none: at hop 1 the arrival handler stops at header count L, never
reading them; at hop 2 they land in `[c + 32L, c + 32k)`, inside sender's slot at or past the final cursor, where
every later read stops. At hop 2 the row port's LC and RC entries persist across groups and `fin` only reads them: the
next group keeps filling block from slot L + 1, its first full block to d rewriting `[c, c + 32L)` with the same L
records (identical bytes, below F, the previous final cursor) before the new ones at `[c + 32L, c + 32k)`; cursors stay
block-aligned. The keep lemma holds at cap_pg too (a multiple of k, checked): the per-group hop-1 sender also keeps or
drops whole blocks. **Fused arrival.** Hop-1 records skip R1 slots: each arriving block runs `h2s` (6 per record: half of
LD64B and ADD, then BFI, LDADD, 2 STP, TBZ), routing as the slot loop did (same row blocks byte-for-byte: part B checks
`port.blocks == snd.msgs`). R1/E1/RC stay (same layout and `fin` writes, still charged) as a conservative control; only
their reads go. The sender-side cap test stays, its drops in per-group row (`d_p1`, 6.6). **Hash table.** Each node
owns 2^57 head words at 3*2^61, zeroed once by `zero4096` (512 DC ZVA, ADD #64 per pass, 1026 per 4096 words, plus RET), and a
packed next-array at 2^61, one 16-byte entry (next address, K2) per R2 slot position (section 5). Bucket `K2 mod
2^57`: equal K2 share a bucket, so same-key records meet without ordering. **Join at arrival (`dlvh`).** A hop-2 sender
in column j writes into header word 3 the address p of F, its slot's final cursor at destination after the previous
group: its E2 entry of group g - 1 (`fin` writes group g's into the other parity table) or, in group 1, the slot base in
the hash-read cursors at 2^48 (an STP of count and p replaces `h2s`'s and `fin`'s STR of count). `dlvh` loads c,
count and p (LDP, LDR), reads F = *p, loads list cursor, its end, C and -2C (LDP, LDP) from join state at 2^48 - 160, sets x10 =
3*2^61, copies all k slots to c, c + 32, ... (per pair LD64B into x16..x23, ADD, ST64B, ADD: 2 per slot; a lone slot 4) and joins exactly the records
in [max(F, c), c + count) while in registers. Within a group a sender's blocks to d cover disjoint ranges at or past F
in any arrival order (its cursor advances by 32k per full block; the final block is last); across groups the barrier
orders them, and the only repeats, the previous final block's L records, come first in the next block to d (identical
bytes, keep lemma) below F. So every record is joined once, in arrival order; stale slots are copied, never joined. A
fresh full block (F <= c, count 8192) unrolls slot 1, 127 pairs, slot 256 (14 fixed and the ADD of x13: 15; copies 516); any other (a final block, or a group's first
block after one) runs three segments (copied, copied and joined, copied): CBZ if empty, else a lone odd first slot (TBZ of x1 bit 5, SUB), pairs 1, 2, then 4 at a time, a lone last slot (7 fixed, 2 per pass);
the maximum over all splits, 633 with copies, is charged on each hop-2 final row with one surcharge (102) for
the next block (section 7). Heads and next fields hold **entry addresses** E = (A>>1) + C (never 0), not slot addresses. Before the joined
segment one ADD sets x13 = E of its first record (x13 = (x2>>1) + C from join state's C); consecutive slots
are 32 bytes apart, so their entries are 16 apart. Per record (4, copy's LD64B or LDP having loaded K2 and K3): BFI of bucket
into x10 (head address); SWP of x13 into head, returning the old head; STP of the old head and K2 at x13,
post-indexed by 16 (the next record's entry); CBZ to end if the bucket was empty (the walk reads no head, so making the
record the head first changes no load). The walk is unrolled by two, prior's entry address alternating between x11
and x5: a visit loads its next entry address into the other register and its K2 with one LDP at prior's entry, then
CMP, B.EQ and CBZ or CBNZ (4); a chain ending after the second half (an even visit count) adds a B to the join's end,
so a walk costs at most 4.5 per visit. A K2 hit verifies d4 (ADD of the slot address A = 2E - 2C with -2C from join
state, LDR of A + 8, EOR, TST, B.NE: 5; a mismatch rejoins the walk); a full hit lists the slot addresses (CMP, B.HS,
STP and, after the first half, B: 4 or 3; + 1 ADD or SUB of x12 = x2 ± 32 for a pair's second record or a lone slot) unless the list is full (CMP, B.HS, 2); the matched visit skips its CBZ or CBNZ (the
run counts 3, the ledger 4). The next-array (C = 2^61 - R2/2) fits every receive pattern up to slot caps by
construction; entries need no zeroing (written at insert before any read); records stay intact in R2 for phase 4. Per-node and per-group costs are rows of section 7; every address a kernel accesses lies in regions of
section 5.

**Chains, no sorting.** A home's records share d1 and d2; records with equal (K2, d4), the 182 label bits,
land in one bucket and form one arrival-order chain, each member linked on every run (walks are uncapped) to its most
recent same-key predecessor. A class's records form one connected chain unless an **intruder**, a record with the same
142 bits not carrying class's label, arrives between two members (6.4, 6.6); every listed pair agrees on d1, d2, K2
and d4, so besides a class's edges it can list pairs of different classes, which the resolver rejects (no unit relates
their points) and 6.4 charges to candidate, request and relay caps through KB(182). Writing the members `u_i*Z` with
Q coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`, an equivalence relation, so a chain with a useful pair
has a useful chain edge.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.0111 used) | code, stack, the five barrier words (one reservation granule each), header H with the rows of β below it, X table of the running group (512 bytes per X), X affine store (64 bytes per X, all groups), scalars, batch entries with the Y words (512 bytes per Y: 2^34 bytes at a multiple of 2^35), the ports' message blocks (at 2^48 - 64 and 2^48 - 128), their LC tables at 0 and 2^39 and RC tables at 2^43 and 2^43 + 2^39 (2^39 bytes each), the hash-read cursors at 2^48 (2^39 bytes), the join state at 2^48 - 160 (32 bytes) and `h2s`'s p word at 2^48 - 168 (8 bytes), the record blocks at 2^50 and 2^51 (2^50 bytes each), E1 and two E2 tables by group parity (2^41 each) and the phase-4 buffers and tables (`FIXED` in `certificate.py`); all but the aligned tables lie between 2^44 and 2^48 | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | per-group hop-1 slots of cap_pg at 2^52 (`2^59.9986`, allocated but never written under fused arrival: a conservative control, section 4), the hash next-array at 2^61 (`2^61.9986`, R2 address A's entry at (A>>1) + C), the hash heads at 3*2^61 (`2^60`: a head's address is 3*2^61 with the bucket in bits 3..59); disjoint, inside R1 (checked), together 2^62.5838 | per group |
| R2 at R1's end + 32 (32 mod 64, for LS64 pairs), ending at 2^63.998943 | cumulative hop-2 pair slots (phase 4 reads candidates' records at their addresses) | groups 1-32 |
| 2^64 - 2^51 to 2^64 | β tables: 2^43 entries of 256 bytes, slot k at +32k holds window k (k = 0..5; `tbld`) | all phases |

So every node stays < 2^64 bytes, and the peak is `2^72*(fixed + |R1| + 32 + |R2| + 2^51)`, ≤ 2^135.9990
bytes (log2 ≈ 135.9989443). A remote item is reached only by a message.

## 6. Success (no heuristic)

In `success.md`: model (6.1), Chen–Stein on cross-node pairs (6.2), usefulness (6.3), defect rows (6.4), bound (6.5) and
grouped hash-chain abort (6.6: stopping-time bound, chain detection E1-E7, uncapped walks with history-conditional walk
means and a defective-run allowance).

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances itemized in comments of
`certificate.py`, which checks generator construction, node set-up, a broadcast message, the candidate, the count
exchange and barrier against itemized work.
- **Generator construction, 645,410 per point (645,414 for a grouped X, whose four affine-store STPs, 64 bytes, the X row charges):** a complete Jacobian double-and-add, exceptional-case tests and two inversions with affine conversion, all on the exact arithmetic of 3.4: `fe_mul_f` and `fe_sqr_f` at 105 and 80 instructions as executed, plus 4 and 3 per call as in fe_inv, and an addition or subtraction with operands brought below p at 27; 633,952 for the 256 doublings, 256 mixed additions and two conversions; the exceptional-case tests (Z and H brought below p and zero-tested per iteration, 13 each, and 64 for paths they select) add 6,720. The scalar's bit tests and pointer step add ≤ 6 per iteration, 1,600 in all with set-up and return. A Y's table entry (y², b², b⁴, b³, βb, 3b², the correction words `-T(b)` and `-T(βb)` at 30 each: `T(v) = (v >> 128) + e`, e = 1 iff `2^128 (v mod 2^128) + (v >> 128)*C0 >= p`, by 3.3's e0; seven values brought below p at 9 each; fourteen STPs: 658) or an X's mixed addition `R_g + [a_k]P` with its header (a², a⁴, y², a³, βa, 2a³, 2a, the nine rows `c*2^(64k) mod p`, k = 1..3, of 2a³, 2a and a², `-T(a)` and `-T(βa)`, seventeen values brought below p, thirty-two STPs: 3,138) brings the total to ≤ 645,410, the allowance (charged at its itemization as in aburan28's 7e8d13c7 and ecc3ec1b). A Y (no mixed addition, no X header) costs 642,930 (`GEN_MULY`, checked; pepedesigner's idea, p256 1772672b).
- **Node set-up, 2^20:** a program of ≤ 2^18 bytes (26,648 kernel instructions expanded, generator construction, resolver, phase control, `h2s0` once), 2^13 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs, three STPs, ADD of 8 to its destination word, SUBS, B.NE), plus 2 for in-flight count and 6 for delivery: 16.
- **Candidate, 1,557,745:** two requests and two replies of ≤ 8 records over 2 hops at 40; each record's three possible leaves `φ^t(X_k) + sY_l` (one inversion of `D = x_Y³ - x_X³`, inverses `Q_t/D` of the three denominators, and three affine additions: 26,753); nine unit tests at 600; the solve mod r and `[x]P = Q`, on the same exact arithmetic; 778,872, charged at twice that plus 1 (as in 1fbd26e8).
- **Phase-4 count exchange, 3*2^39 per node and sub-phase (E-averaged over groups run):** ≤ 12 per entry of node's 2^36-entry table of received counts (load, epoch test, select, prefix, zero test, cap test, loop) and ≤ 88 per message (its count, count message, the offset reply and write, ≤ 22 each with in-flight count and delivery) for ≤ 8*2^29 messages (2^29 requests of up to 8 records): < 2^40.5.
- **Barrier, 2^10 per node:** ≤ 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), ≤ 5 messages of 14 (6, plus 2 for in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^16 per node:** 12 kernel calls per group ≤ 33 each (arguments and BL), stop test and group loop (8), over 114 groups: < 2^16.
- **Hash-head zeroing, 1026 per 4096 words:** `zero4096` (512 of `DC ZVA, x1` and `ADD x1, x1, #64`, one 64-byte block each, then SUBS, B.NE; DC ZVA has no offset form, hence the ADD; its macros only shorten listing); row's +1 is its RET. The heads at 3*2^61 are 64-byte aligned in base and size (`a64.py` asserts the alignment) and in Normal memory, never a port: DC ZVA on Device memory would fault.
- **Hash-read-cursor init, 4 per slot:** `init3` (STR, ADD, SUBS, B.NE) of each slot's base into its 8-byte entry at 2^48; row's +1 is its RET.
- **Rows of β 72** (sixteen words of MOVZ and three MOVK, eight STPs, once per node), **batch set-up 160** (X's header with its correction words copied into H by 30 LDPs and 30 STPs, y_X² into registers, acc = 1, pointers, sp for `bwd`, index word, the slot end for `keys`, fixup meter's update and test and loop: ≤ 145, itemized in `certificate.py`), **inversion call 16, request served 2^10, output 2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on failure flag, plus a range test. RNDR, this machine's coin source, returns 64 fresh uniform bits (the algorithm coins of rules' probability space) and never reports failure, so the charged B.EQ is never taken and coins add no failure row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, barriers and output. A record row uses
its group's records; E-averaged rows are sums over groups weighted by PHI. With q1g the grouped q1 (145307202):

| row (P: preprocessing; E: E-averaged; W: worst-case) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1g + q2 + 1)*64 (P) | 105.4143 |
| Y generators with their table entries: q2*642930 (P) | 116.2943 |
| R_g generator | 91.2999 |
| rows of β: 72 (P) | 78.1699 |
| row LC and RC tables: (4*2^36 + 1) + (3*2^36 + 1) (P) | 110.8074 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| hash-head zeroing: zero4096 2^57/4096*1026 + 1 (P) | 127.0028 |
| hash-read-cursor init: init3 4*2^36 + 1 (P) | 110.0000 |
| leaf work, E-averaged: sum_g PHI[g] (X gens, fwd, inv, bwd, keys w/o emits, batch set-up) (E) | 134.3704 |
| β tables: `tbld` 4 + 72 + 20*(5*2^43 + 2^41) (P) | 121.7142 |
| arithmetic fixups: 2^45*5 (W) | 119.3219 |
| j = 1 and j = 2 key paths, E-averaged: sum_g PHI[g] s_g 4 PJ + DEFJ (E; worst 128.6995) | 127.8948 |
| route work, E-averaged: sum_g PHI[g] (emits, arrival, row blocks, fin, column tables) (E) | 130.3576 |
| group work, E-averaged: sum_g PHI[g] (join 4 per record, main-pair walks at 4.5 per visit, stop test) (E) | 129.1857 |
| same-X walks: global pair counting (E) | 29.0394 |
| K2 verifications: global pair counting (E) | 124.8909 |
| listings: global pair counting (E) | 76.3059 |
| defective-run rows: DEFV + DEFK2 + DEFL + DEFJ (E) | 63.2308 |
| candidates: 2^28*1491011 (W) | 120.5079 |
| requests served: 2^29*2^10 (W) | 111.0000 |
| phase-4 count exchanges, E-averaged: sum_g PHI[g] 4 * 3*2^39 (E; worst case 114 * 4 * 3*2^39 = 2^121.4179) | 121.1981 |
| barriers: 2048*2^82 (W) | 93.0000 |
| phase control: 2^16 (W) | 88.0000 |
| output + stop broadcast (W) | 85.1548 |
| **expected total E** | **134.518995** |

`E` is the exact `E_total` of `certificate.py` (expected instructions), ≤ 225.8939 per leaf (`per_leaf`), and
`W_E = E/1000` with log2 W_E ≈ 124.5532102949. **The claim, 124.553210294878063, is the least multiple of 10^-15 with
2^claim >= W_E**; the certificate encloses 2^c by an exact series for e^(f*ln 2) with
interval bounds on ln 2, both ways (2^(claim - 10^-15) < W_E), from an 80-digit decimal start. Rows marked P are
preprocessing: their sum over 1000 has log2 ≈ 117.0743478, declared as 2^117.0744, the least multiple of 10^-4 at or
above it (checked). The worst case R_worst = 2^188.604975 (all groups run, uncapped walks to HOME_MAX visits with
verification and listing on every record) is disclosed, not scored: the least multiple of 10^-6 at or above log2 R_worst
(checked).
**The score is the expected cost** (total classical work in expectation over the input and the algorithm coins, early abort included): a run that never stops early (probability about 0.61) costs about 2^124.773 units, disclosed (certificate's E with all PHI[g] = 1). **DC ZVA** (change 4) and **LD64B/ST64B** (changes 17, 18) count one step like every instruction; priced as STPs (DC ZVA) or LDP/STP pairs (LS64), the claim would be ~0.0076 and 0.0169 bits higher.

## 8. Evidence

`a64_load` asserts every branch's architectural reach (TBZ/TBNZ ±32 KiB, B.cond/CBZ ±1 MiB, B ±128 MiB) on the laid-out code; `certificate.py` checks its figures against literals and submission.json's statement. `experiment.py` takes `{"seed": 1, "seconds": 12}` (`seconds` 5 to 15); budget-limited checks run guaranteed minimums
(2,000 random and 600 near-tie class-rule inputs, 3,000 class-label classes), then their budget share up to 26,000 and
12,000. **Part A** runs leaf kernels and fe_inv on production batches of 4, 8, 12, 20 and 36 incidences, checking every
stored fraction against a bit-exact model of `bwd` (`fraction_mismatches`) and its bound (`fractions_outside_bounds`),
each product's terms (`tfp_terms`), 3.1's identity against exact affine addition (`identity_mismatches`), each record,
slot and LC entry against bit-exact model of `keys`, and each index and label outside W and D' (`wrong_index`,
`keys_not_canonical`). Random and near-tie leaves in random representations (`represent`) get the exact index outside W
(`class_rule`); on ≥ 3,000 classes, a third at an edge of W or D', every member not in W or D' emits its class's
label (`class_labels`); `window_batch` and `boundary_batch` plant such leaves in batches. `fwd`/`bwd` counts match their
formulas plus 3, 4 or 5 per fixup (`fwd_bwd_formulas_hold`); planted fixups (`fixup_batch`) and an always-inline copy
(`forced_fixups`: 33 fixups, 106 more per incidence; `window_sum_fold` forces BLK0's carry with crafted entries, `beta_tables` builds short tables with `tbld`) keep exact fractions. `site_checks` (`sites.py`) tests each site
type of 3.4 against its predicate, bit-level model and fixup; a `Y_l = φ(X)`, `Y_m = -X` batch gets acc = 0 and no other
fraction near its value; site macros match tables (`sites_match_certificate`); flags are dead at every fixup join. **Part
B** (`machine.py`) runs a 4 × 4 machine through `leaf.s` and `route.s`: `init3`, `tinit`, `keys` on planted
exact components (labels with small digits at production widths, some classes shared), `fin`, fused `h2s` arrival (byte
count from header, row blocks byte-exact) and hop 2, delivering by executed `dlvk` or `dlv`: messages, R2 slots and
source tags equal the sender model. `block_paths` forces the remaining block and cursor paths (`h2s` on
1..19, 24, 44, 256-record and empty blocks); `two_group_finals` checks group 2 rewriting group 1's partial final block.
`grouped_paths` delivers three groups of row blocks from two interleaved senders (fresh full, repeated-prefix and final
blocks, stale slots: `fast_slow_blocks`; group-1 p from `h2s`/`fin`, later p from the model) through `zero4096`, `init3` and `dlvh` on true cumulative R2 layout with
planted cross-group same-key pairs (stop, full and 1-entry-list cases): copies, heads, next entries and candidates equal
the chain model, counts their formulas, cursors reach each slot's end, every record is joined once, every listed pair is
same-key and every global pair connected; old X entries persist and are served (`old_x_served`); `ls64` rejects odd or high Xt, misalignment, other regions, ST64B to the window, R2 at 0 mod 64, LD64B at slot 1. **Part D** runs
`certificate.py` (`D_certificate.checks`) and `counts_match_certificate` compares its constants with A and B.

The 4 × 4 run is evidence at its own scope only.

## 9. Limitations and provenance

- **Model.** The bound is for machine of §2.1 (links keep no state beyond the message in flight; sender's instructions, each delivery as `dlv`, `dlvk`, fused `h2s` arrival or `dlvh`, and in-flight counts are charged); a forwarding network, adding work per extra hop, is not claimed. FEAT_LS64 on the windows and R2 is a machine premise (2.1).
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control, batch set-up with driver's fixup-meter test, stop test and phase-4 exchanges), the in-flight counts, and production attack, machine and links. The toy stop fires on a listed candidate; production stops only on a resolver-confirmed x, a stronger rule, so the E-bound (charging confirmed stopping) is conservative over toy.
- **Provenance**, not a premise: kaankolcu's 38df975c (declared parent); ob2 (the grouped truncated-fraction line: machine, resolver, success rows, Fourier leaves, records, label, address forms, fe_inv, generator, meter, E-ledger, hash chains); 66043fba (three-pair β cache); Galbraith, Wang and Zhang for X ± Y from one inverse; Arratia, Goldstein and Gordon for Chen–Stein; kaankolcu's 8c3852be (declared) for DC ZVA zeroing; aburan28's c3e0530e (declared) for the four-body `bwd` loop. The run and `certificate.py` re-verify code. Changes against 38df975c (bits, each measured on the code before it): 1 hoisted emit base (−0.001716); 2 fused hop-2 join in `dlvh` (−0.005979); 3 β by 43-bit-window tables, unused base path removed (−0.030329); 4 DC ZVA head zeroing (−0.007606; priced as four STPs per block: +0.0076); 5 entry STP by writeback (−0.005305); 6 entry-address chains (−0.002642); 7 13 key sums (−0.003569); 8 one-ORR window addresses (−0.007165); 9-11 CPK after the first row, γ's V0 in (s, hi) form, BLK0 fold branching at limb 1 (together −0.007225); 12-14 fwd unrolled by 32, EON record word, h2s unrolled on full blocks (together −0.003337); 15 delivery processors stated, `h2s0`, p after the window wait, and 16 self-checks (branch reach in `a64_load`, figures against submission.json): no figure change; 17 `h2s` reads pairs by LD64B (7 to 6 per record; −0.005381501447682); 18 `dlvh` copies pairs by LD64B/ST64B (4 to 2 per slot; R2 at 32 mod 64, +1 on some listings; −0.011537090797751); 19 `bwd` four bodies per loop test (1230.5 to 1230.25; −0.000228499490883). Items 2-14 were measured in a sequence that also held three levers since removed (a key-path tree, a duplicated window-3 slot, out-of-line hit blocks); their deltas then (−0.004991) and their removal on the final code (+0.005104) net +0.000112257162957 against the rounded deltas of items 1-14 (−0.074873). Sum: −0.091907834573359 = 124.553210294878063 − 124.645118129451422.
