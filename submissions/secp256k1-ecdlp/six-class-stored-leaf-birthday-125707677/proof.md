# Stored x-only leaf birthday search on secp256k1 (125.707677)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm of this document
recovers `x` with probability at least 0.39001 over the uniform challenge and fresh algorithm coins. The success bound is
proved in section 6 and uses no heuristic. On every run, its total work is at most

```text
time_log2 = 125.707677 = 125707677/1000000
```

target-group-operation equivalents (C = 1000 scalar AArch64 instructions per unit). Peak memory is below 2^135.7542 bytes.
Preprocessing (work independent of Q) is 2^111.0349 units and is included in the total. Nonuniform advice is zero.
`certificate.py` recomputes the parameters, the success bound, the ledger, the claim, the memory sum and the preprocessing
with exact integers and Fractions, and prints PASS. This is an analytic bound, not a practical break.

## 1. What this is

Pollard rho needs y for its next step. A stored birthday search does not. It draws independent points X and Y and stores only
the x-coordinates of X + Y and X − Y (the "leaves"). Both leaves share one batched inverse, and no y is ever computed for them.
Each leaf goes into a 4-word record keyed by its six-class canonical x. The records are routed to home nodes and grouped
there by sorting. Equal keys are then resolved into the logarithm. X ± Y from one inverse is Galbraith, Wang and Zhang's
(ePrint 2015/605). The field schedules, the canonical-index row and the counting interpreter come from our rho package F
(submission 8f733826). Their counts here are re-executed, not inherited. No prior acceptance or score is a premise.

Three properties carry the claim:
- **Rigorous success.** Leaves from different nodes are independent uniform points. Their collision events are linear
  relations among independent generators, so the Chen–Stein method bounds the probability of no collision (section 6).
- **Executed per-record work.** The leaf arithmetic, the key, the emit, both handlers' event paths and bodies, the second
  hop's send, the local pass and the scan are shipped AArch64 kernels, run by the counting interpreter `kernels/a64.py`
  (section 8). The leaf kernels run on production points against exact arithmetic. The transport kernels run on a 4 × 4
  machine with windows, event queues and flags, against an exact model of the same order.
- **Priced transfers.** Each record hop is priced at 40 instructions: its executed protocol, 17.5 for hop 1 and 29 for
  hop 2, plus a margin. The local pass is priced at 20 per record (executed 9) and the scan at its executed 7 (section 4.1).

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`. A node has a private memory below 2^64 bytes (section 5) and two
scalar AArch64 processors that share it: a main processor and a handler processor. The cost model counts the
instructions every processor executes; storage and wiring carry no price.
- **Links.** Each ordered pair (s, t) of nodes in one row has one row link, and each in one column has one column link. A node's
  row and column links to itself are its two loopback links, so a record addressed to its own node takes the same path. Each
  link ends in a
  64-byte window in t's memory: a 32-byte record slot and a full flag. A node holds 2^37 windows (2^43 bytes): one for each
  node of its column and one for each node of its row. It clears their flags in preprocessing.
- **Send.** The sender writes a destination index and a record to its 40-byte port in three stores:
  - STR of the index at +0;
  - STP of K2 and K3 at +8;
  - STP of w2 and w3 at +24.

  A node has a column port and a row port, each a fixed device register. The index selects the link and is not carried.
  The third store (w2 and w3 at +24) completes the message. The link then carries the four record words into the window and
  sets the flag; it executes nothing.
- **Backpressure.** A store to a link whose window is full stalls until the flag clears. A stalled store is one executed
  instruction: the rules count instructions, "not CPU cycles".
- **Receive.** Each arrival appends the window's index to the receiver's event queue. The queue has 2^37 entries, and a
  window holds at most one pending arrival. The handler:
  - loads the entry;
  - forms the window address;
  - reads the record;
  - runs the phase's receive body;
  - clears the flag.

  That is 17.5 instructions per record for hop 1 (with the sender's emit) and 29 for hop 2 (with `h2send`), executed in
  section 8. The handler waits for arrivals with a load that stalls, and its pointer wraps at the end of the ring; this
  dispatch is itemized in 4.1 and fits the hop price.
- **Progress.** A handler never sends or stalls, and it has its own processor and registers. So every window drains, and the
  main processor's registers are untouched.
- **Two hops.** 2^37 windows fit a node's 64-bit space; windows to all 2^72 nodes would not.
- **Barriers.** A node joins a barrier when two things hold:
  - its main program has finished the phase's sends;
  - its ports' counters of uncleared windows read zero (one load, which stalls until zero).

  The barrier is a combining tree over rows and then columns. It completes when every node has joined, so no record is then
  in flight.

### 2.2 Generators and leaves

Node `g` draws fresh uniform scalars from its coins:
- `b_g` in `[0, r)`;
- `a_1..a_q1` in `[0, r)`;
- `c_1..c_q2` in `[1, r)`.

Here `q1 = 108,966,605` and `q2 = 2^27`. The node computes `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P` (section 7),
and stores the points and the scalars. The node's **leaves** are the x-coordinates of `X_k + Y_l` and `X_k − Y_l` for all `k`
and `l` (a complete grid). There are `n0 = 2·q1·q2 = 29,250,500,301,946,880` leaves per node and `S = 2^72·n0 = 2^126.699311` in
all.

Leaf `X_k + s·Y_l` equals `[m]P + [n]Q` with `m = a_k + s·c_l` and `n = b_g`.

### 2.3 Class key, record and digits

Each point Z other than O has a class `{±Z, ±φ(Z), ±φ²(Z)}` of 6 points, where φ(x, y) = (βx, y) = [λ]. Its **key** is
`K = min(x, βx mod p, β²x mod p)`, an integer in `[0, p)` with limbs `K0..K3`. K is the same for all six members and differs
between classes.

| field | bits of K | role |
| --- | --- | --- |
| `d1 = K0 mod 2^36` | 0..35 | hop-1 destination row |
| `d2 = K1 mod 2^36` | 64..99 | hop-2 destination column |
| `K2`, `K3` | 128..255 | stored key words |
| `d3 = K2 mod 2^24`, `d4 = K3 mod 2^48` | 128..151, 192..239 | home sort digits |

A record is `[K2, K3, w2, w3]`:
- `w2 = i | k<<36` (source row, X index);
- `w3 = d2 | l<<36 | s<<63` (Y index, sign).

The source column `j` replaces d2 in hop 2.
- **Sort bits.** Records at one home share d1 and d2 and are grouped by d3 and d4: 36 + 36 + 24 + 48 = 144 bits.
- **Compared bits.** The scan compares K2 and K3, so 36 + 36 + 128 = 200 bits of K.

### 2.4 Phases

0. **Preprocessing.**
   - Program delivery: 2 hops from node (0, 0).
   - Node state and coins.
   - The Y generators and the header constants.
   - The window flags cleared, the d4 counters zeroed, and the d3 cursor table.

   Q is then broadcast to every node in 2 hops, 2 records each; this is not preprocessing.
1. **Generation and hop 1.** The main processor runs the batch kernels of section 3 for each X_k. The emit sends each leaf
   through the column port to `(d1, j)`. That node's handler runs `h1recv`, which appends the record to its hop-1 receive
   buffer of `cap_c` records. A record arriving at a full buffer is dropped; its flag is still cleared.
2. **Hop 2.** The main processor runs `h2send` over its buffer. It sends every record through the row port to
   `(own row, d2)`, writing `j` over d2. The home's handler runs `h2recv`. It appends the record to its d3 bucket, of
   capacity `cap_3` (a full bucket drops it), and counts d4.
3. **Grouping.** The home node forms exclusive prefix sums of the d4 counts, then runs one local pass. That pass writes each
   record's address into a pointer array at its d4 cursor, bucket by bucket in d3 order. The order is therefore by (d4, d3),
   stable within a bucket. A scan compares K2 and K3 of adjacent records and lists equal pairs as candidates, up to 2^20. A
   full list ends the scan.
4. **Resolution**, in four bulk-synchronous sub-phases. A home `(h_r, h_c)` and a source `(g_r, g_c)` exchange messages
   through two relays: `(g_r, h_c)` for requests and `(h_r, g_c)` for replies.
   - **4a.** For each candidate the home sends two requests along its column to the relays.
   - **4b.** Each relay's handler has buffered its requests. Its main processor forwards them along its row to the sources.
   - **4c.** Each source answers at most 2^21 requests. A reply carries the scalars and points of one leaf, at most 8
     records, and goes along the source's column to its relay.
   - **4d.** The relays forward the replies along their rows to the homes.

   A relay buffer holds at most 2^22 records per sub-phase and drops the excess. The home then recomputes both leaf points
   `L = [m]P + [n]Q` and `L'` exactly and finds the unit u with `L' = uL`, if one exists. If `n' − u·n ≠ 0`, it solves
   `x = (u·m − m')/(n' − u·n) mod r` and checks `[x]P = Q` before output.

**Barriers.** Thirteen barriers are used:
- program delivery, 2;
- Q broadcast, 2;
- the ends of phases 0 to 3, 4;
- phase-4 sub-phases, 4;
- output, 1.

The allowance covers 16.

### 2.5 Caps, meters and drops

- **Hard work bound.** The ledger of section 7 holds on every run, because each data-dependent cost is capped or bounded in
  total:
  - the class-index products (section 3.4): before each batch, a node continues only if `count + 2q2 <= cap_j`;
  - handler events: at most S per hop in all, each at most its main path;
  - buffer, bucket, relay and request stores;
  - candidates;
  - coin retries: one RNDR retry per word and at most two draws per scalar (section 7).
- **Drops.** A capped or dropped item costs at most its main path. A drop can lose success only if it hits a record of the
  collision that would have succeeded. Section 6.4 bounds that event.

## 3. The leaf kernels (`kernels/leaf.s`)

### 3.1 Field arithmetic

p = 2^256 − C0, C0 = 2^32 + 977, held in x4. A field element is four 64-bit limbs. "Weak" means `[0, 2^256)`.

- **Row products** (package F's RF and RA macros, unchanged). For a weak b, its rows are `B_0 = b` and
  `B_(k+1) = [lo, b0 + hi, b1, b2]`, where `(lo, hi) = b3·C0`, so `B_(k+1) ≡ B_k·2^64`. For a weak a, `Z = Σ a_k·B_k ≡ a·b`
  and `Z < 2^322`. RF (12 instructions) forms `a_0·B_0` in five limbs. RA (18) adds `a_k·B_k`. No carry is lost: a high word is
  at most 2^64 − 2.
- **SHX (5)** is F's row shift with the carry taken inline: MUL, UMULH, ADDS, ADCS, ADCS.
- **FDX (7)** is F's fold `t = (Z mod 2^256) + (z4 + 2^64·z5)·C0` with the carry into limbs 2 and 3 taken inline. Its MADD
  forms `z5·C0 + hi <= 4·C0 − 1`.
- **No rare exits.** Neither SHX nor FDX has a branch, so the main path is the only path. The only remaining error is a carry
  out of the top limb, which drops 2^256. That is a rare-value event (section 3.5): it gives a wrong value at the same cost.
  A product on shared rows is 72 + 1 = 73, a pair on one operand's rows `2·73 + 3·5 = 161`, and a single product
  `73 + 15 = 88`.
- **FS (9)** is `a − b` for canonical a and b: SUBS, 3 SBCS, CSEL of C0 on borrow, SUBS, 2 SBCS, SBC. The result is canonical.
- **FADDW (9)** is `a + b` for canonical a and b: ADDS, 3 ADCS, CSEL of C0 on carry, ADDS, 2 ADCS, ADC. The result is weak:
  `a + b − 2^256 + C0 < 2^256 − C0` after a carry.
- **SQX (78)** is the fused square `λ^2 − x_X − x̃_Y + 2^257`, where `x̃_Y = x_Y + 2C0` is stored in the Y entry:
  - products, carry chains, doubling and diagonals: 41, as in F's fe_sqr;
  - `H·C0` and `s = L + H·C0`: 17, with `s4 <= C0`;
  - two subtrahends: 4 LDP and 2·5 instructions. Each borrow chain ends with ADC on s4 instead of SBC, so
    `s' = s − x_X − x_Y + 2p >= 0` and `s4' <= C0 + 2`;
  - fold of `s4'·C0 < 2^65`: MUL, UMULH, ADDS, ADCS, ADCS, ADC.

  That is 41 + 17 + 14 + 6 = 78, and 64 without the two subtrahends. The result is ≡ λ² − x_X − x_Y.

### 3.2 One batch: forward, inversion, backward

A batch is one X and all q2 Y's of its node.
- **Slot.** Slot l, 64 bytes, holds `dx_l` and `acc_(l−1)`.
- **Header H.**
  - `x_X` at +0 and `y_X` at +32;
  - a spill for R at +64;
  - `s−` at +96 and `s+` at +128;
  - `G1..G5` at +160;
  - the rows of β at +208 and of β² at +336.
- **Y entry**, 96 bytes: x, y, x̃.

```text
forward  (l = 0..q2-1):  dx_l = x(Y_l) - x(X);  slot_l = (dx_l, acc_(l-1));  acc_l = acc_(l-1)*dx_l      (acc_(-1) = 1)
inverse:                 R = acc_(q2-1)^(-1)                                                          (fe_inv)
backward (l = q2-1..0):  inv = R*acc_(l-1), R <- R*dx_l;  s+ = y_Y - y_X, s- = y_Y + y_X;
                         lam+ = inv*s+, lam- = inv*s-;   x+- = lam+-^2 - x_X - x_Y -> slot_l
```

Before step l, `R = 1/acc_l`, so `inv = acc_(l−1)/acc_l = 1/dx_l`. Then `λ+ = (y_Y − y_X)/(x_Y − x_X)` is the slope of
`X + Y_l`, and `λ− = (y_Y + y_X)/(x_Y − x_X)` is minus the slope of `X − Y_l`. Only `λ^2` is used, so `x+ = x(X + Y_l)` and
`x− = x(X − Y_l)`.

| kernel (executed) | per incidence | fixed | content |
| --- | ---: | ---: | --- |
| `fwd` | **105** | 3 | 2 LDP, FS 9, 4 STP, single product 88, SUBS, B.NE |
| `fe_inv` (`secp_fe.s`, F's) | — | **23,646** | addition-chain inversion, executed and checked against `pow` |
| `bwd` | **523** | 5 | W2S 165; 2 STP of R; 4 LDP of y_Y, y_X; FS 9; FADDW 9; 4 STP of s±; W2L 165; SQP 78 and 2 STP; SQM 78 and 2 STP; 2 LDP of R; SUB; SUBS; B.NE |

- **The pair macros.** W2S is F's pass-2 pair (W2P) with SHX and FDX: 4 LDP + 161 = 165. Its operands are at
  `[x3, #32]` and `[x3]`, and x3 steps back 64. W2L is the same schedule with x10..x13 and x25..x28 exchanged. Its shared
  operand is therefore inv, and its operands are s+ and s− in H.
- **Registers.** Both pairs use the 25 pool registers, x30 included, so R is spilled around W2L. The calling convention:
  - a kernel may overwrite x0–x29 and the flags, and preserves sp;
  - `bwd` and `keys` save x30 on entry and restore it before RET. Their macros use x30 as a scratch register; this is 2 of
    their 5 and 11 fixed instructions;
  - `fwd` and the transport kernels do not write x30.
- **Batch driver.** It keeps nothing in registers across calls. It keeps the node's product count in memory and loads it into
  x18 only for `keys`, so fe_inv's use of x18 is harmless.
- **Inversion call.** The inversion is called once per batch. Storing acc, setting pointers, the BL and loading R add 16.

### 3.3 The class index (macro KEY, 21 instructions)

This is F's canonical-index row. Here a window means the leaf is dropped instead of taking a fallback.
- **Registers.** `G1..G5` are limbs 1..5 of `G = floor(β·2^384/p)`. F0 is the top limb of x.
- **F1.** `F1 = Σ_i [lo(x_i·G_(5−i)) + hi(x_i·G_(4−i))] mod 2^64`: 4 UMULH, 3 ADD, 4 MADD.
- **Index.** `S = F0 + F1`, `d = F0 − F1` (C iff F0 ≥ F1), `m = min(F0, F1)`, `e = m + S`. Then j = 2 if e carries, else 1 if
  F0 ≥ F1, else 0. The row base is selected by two CSELs, 0 for j = 0.
- **Windows.** One CMN, two CCMN and one CCMP leave C set iff none of these holds:
  - (W1) `F1 >= 2^64 − 16`;
  - (W2) `S >= 2^64 − 17`;
  - (W3) `e mod 2^64 >= 2^64 − 33`;
  - (W4) `d < 16`.

  B.LO drops the leaf.

**Errors of F0 and F1.** Let `v0 = x`, `v1 = βx mod p` and `v2 = β²x mod p`. They are distinct, since 7 is a non-residue and
no point has x = 0, and `v0 + v1 + v2 ∈ {p, 2p}`. Let the reals `V_i = 2^64·v_i/p`. Then `V0 + V1 + V2 ∈ {2^64, 2^65}`.
- **F0.** `F0 = V0 − a` with `0 <= a < 2`, and `a < 1` unless `x mod 2^192 > 2^192 − C0`.
- **T.** Let T be bits 320..383 of `x·G`. Terms with `i + k >= 6` vanish mod 2^384. The terms with `i + k = 5` and the high
  halves of `i + k = 4` are F1. The rest are the low halves at 2^256 (sum < 4·2^320) and `Σ_(i+k<=3) x_i·G_k·2^(64(i+k))`
  (< 4.01·2^320). They carry at most 8 into bit 320, so `T = F1 + c (mod 2^64)` with `0 <= c <= 8`.
- **F1.** `G = β·2^384/p − θ` with `0 <= θ < 1`, so `T = floor(V1 − ε) mod 2^64` with `0 <= ε < 2^-64`. Hence
  `F1 ≡ V1 − b (mod 2^64)` with `0 <= b < 10 < 16`. Outside (W1), `F1 = V1 − b`.
- **F2\*.** With `F2* = 2^64 − S`, outside (W1) and (W2) `F2* = V2 + a + b`.

**Exactness.** Outside the windows the index is argmin_j v_j, except when j = 0 and `x mod 2^192 > 2^192 − C0`. Take each
returned value in turn.
- **j = 2.** `V0 >= F2* + a >= V2` and `V1 >= F2* + b >= V2`.
- **j = 1.** (W3) gives `F2* − F1 >= 34` and (W4) gives `F0 − F1 >= 16`. So `V0 − V1 > 16 − 16 = 0` and
  `V2 − V1 > 34 − 2 − 32 = 0`.
- **j = 0.** `V2 − V0 > 34 − 4 − 16 > 0` and `V1 − V0 >= 1 − a > 0` when `a < 1`.

In each case the values are distinct, so the inequality is strict.

**Windows are rare.** Each window confines one of the following residues to fewer than 2^199 values:
- (W1) puts v1 within 2^196 of 0 or p;
- (W2) puts `(1+β)x mod p` within 2^197 of 0;
- (W4) puts `(1−β)x mod p` within 2^197 of 0;
- (W3) puts `(2+β)x` or `(1+2β)x mod p` within 2^198 of 0.

Multiplication by these nonzero constants is a bijection, so at most `13·2^197` x-values hit a window. A uniform point's x
takes each value with probability at most `2/(r−1)`, so a leaf is dropped with probability below 2^-54.29. The wrong-index
case needs `x mod 2^192 > 2^192 − C0`, at most `2^64·C0` values of x, so its probability is below 2^-157.

### 3.4 Key pass and emit (`keys`)

**Per leaf.** The pass works one leaf at a time:
- it loads x (2 LDP), runs KEY, and branches B.LO to drop the leaf;
- it then takes CBZ on the row base, 0 for j = 0;
- for j = 1, 2 it runs KMUL, then ADD to the product counter x18. KMUL (81) computes `K = x·β^j` on the constant rows in H:
  the row limbs are loaded in pairs, 8 LDP, interleaved with RF, three row accumulations and FDX;
- the emit takes 6 instructions, 7 for the minus leaf:
  - AND of d1 and of d2;
  - ADD of the slot address shifted left by 30, which puts l in bits 36..62 because slots start at a multiple of 2^34;
  - ORR of the sign bit, for the minus leaf only;
  - STR of the destination at `[x0]`, then STP at `[x0, #8]` and `[x0, #24]`. The column port is a fixed device register,
    so x0 does not move.

**Costs.**
- A j = 0 leaf costs 31 (plus) or 32 (minus).
- A dropped leaf costs 24.
- A product adds 83: KMUL 81, the counter and one B.
- An incidence costs **66** with no products: two leaves (63, including emits 6 + 7 = 13) and the loop (3). The fixed part
  is 11, including the save and restore of x30.

For every leaf other than O, j is uniform on {0, 1, 2}, because a uniform point is a uniform member of its class. The node
caps its products at `cap_j = ceil(2·n0/3) + 2^42` (section 2.5).

### 3.5 Rare values

The kernels compute exactly unless some field value lands in one of a few tiny sets:
- **Row shift.** A carry out of a row shift needs the row's low 192 bits within 2^33 of 2^192.
- **Fold.** A carry out of a fold needs the product to be below 2^100 mod p.
- **Square.** A carry out of SQX's fold needs `x± < 2^66`.
- **Non-canonical result.** A non-canonical x± or K needs the true value below C0.
- **Index.** The wrong-index case of section 3.3.
- **Degenerate point.** `dx = 0` (Y = ±X) or X = O.

In each case the relevant value is a bijective image of a coordinate determined by one fresh generator, given the others:
- the batch's acc, inv and R are bijective in the newest dx;
- λ± is the slope through X and `Y_l`, which takes each value for at most two `Y_l`;
- x± is the x of `X ± Y_l`, which is uniform on at least r − 1 points.

So each site has probability at most `2^101/(r−1)`. A batch has fewer than `40·q2 + 64` sites, so a batch has a rare value with
probability below 2^-122.6 (section 6.4). The rare values change no instruction count.

## 4. Transport and grouping (`kernels/route.s`)

### 4.1 Prices

A record hop is priced at **40** instructions: its executed protocol plus a margin.

| step | executed per record | priced |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (6 or 7) + `h1recv` (11) | 17.5 | 40 |
| hop 2: `h2send` (9) + `h2recv` (20) | 29 | 40 |
| local pass `pass4` (LDR, AND, ADD, LDR, STR, STR, ADD, CMP, B.LO = 9; 7 per bucket) | 9 | 20 (one-word item) |
| scan (LDR, LDP, CMP, CCMP, B.EQ, SUBS, B.EQ or B.NE = 7) | 7 | 7 (moves nothing) |

- **`h1recv` (11).**
  - Event path (3): LDR of the event, ADD of the window address, and the STR that clears the flag.
  - Body (8): 2 LDP of the record, CMP and B.HS for a full buffer, 2 STP appending it, SUBS and B.NE.
  - A dropped record costs 9.
- **`h2send` (9).** 2 LDP, AND, BFI, STR, 2 STP, SUBS, B.NE.
- **`h2recv` (20).** The event path (3) and 17 more:
  - 2 LDP of the record;
  - d3 and its cursor-table entry (AND, ADD, LDP cur, end);
  - CMP and B.HS for a full bucket;
  - 2 STP of the record, STR of the cursor;
  - the d4 count (AND, ADD, LDR, ADD, STR);
  - the loop (SUBS, B.NE).

  A record dropped at a full bucket costs 12.

**Margins.** The margins are 22.5 instructions per record for hop 1, 11 for hop 2 and 11 for the local pass. The destination
word is consumed as link selection, not carried; flow control is the window flag.
- **Dispatch.** The handler runs `h1recv` or `h2recv` as an endless loop over its event queue. The queue is a ring of 2^37
  entries whose posted count is a device register. When the handler has taken every posted entry, its load of that register
  stalls until the next arrival. Per record this adds at most 4 to the executed bodies:
  - the register load, a compare and a branch;
  - the taken count;
  - a BIC that wraps the ring pointer;
  - with B in place of the counted SUBS and B.NE.

  If the kernel were instead called once per event, dispatch would cost at most 7 per record: the wait, the call, the
  kernel's CBZ and RET, and the wrap. That would make 24.5 and 36.
- **Hardware work charged in software.** Suppose the event append is charged to the sender (an atomic fetch-add and a store),
  and the uncleared-window counters to both ends (an atomic increment and decrement). With dispatch included, hop 1 then costs
  at most 25.5 and hop 2 at most 37, both below 40. `certificate.py` checks both cases.

The d3 masks and table bases are registers, so the 4 × 4 run of section 8 uses the production code with small digits.

### 4.2 Fixed per-node work

- **Window flags.** `wclear` clears 4 window flags in 7 instructions; there are 2^37 windows.
- **d3 cursor table.** `init3` writes the cursor table, 5 per bucket, 2^24 buckets.
- **d4 counters.** `zero4` clears the 2^48 counters, one instruction per counter; `prefix4` forms their exclusive prefix
  sums as pointer-array addresses, 14 per four counters.
- **Bucket loop.** `pass4` spends 7 per d3 bucket.
- **Scan.** The scan is `7n − 3` for n records, plus 4 per listed candidate. A full list ends it, so it spends at most
  `4·(2^20 + 1)` on candidates.

### 4.3 Why one local pass suffices

- **Grouping.** After hop 2 the home holds records sharing d1 and d2. The d3 buckets and the stable d4 pass leave records with
  equal (d3, d4) adjacent, so 144 key bits are grouped.
- **Runs.** Records of one class have equal K, hence equal sort bits, and form one run unless an **intruder** sits inside it:
  a record outside the class with the same 144 bits. Such a run has all its adjacent pairs listed unless the candidate list
  is full (6.4).
- **Resolution.** Write the run's members as `u_i·Z` with Q coefficients `n_i`. The pair (i, j) is useless when
  `n_j = (u_j/u_i)·n_i`, that is `n_i/u_i = n_j/u_j`. That is an equivalence relation, so if two members of a run form a
  useful pair, some adjacent pair is useful.
- **Intruders.** Section 6.4 bounds them.

## 5. Node layout and memory

All regions lie in the node's single 64-bit address space. The table gives sizes in bytes; `certificate.py` sums them.

| region | bytes |
| --- | --- |
| code, constants, stack, port and flag state | 2^20 |
| header H | 512 |
| X table (x, y) | 64·q1 |
| Y table (x, y, x + 2C0) | 96·q2 |
| scalars | 32·(q1 + q2 + 1) |
| batch slots (start at a multiple of 2^34) | 64·q2 |
| hop-1 receive buffer | 32·cap_c, cap_c = 6·n0 + 2^50 |
| d3 buckets | 2^24·cap_3·32, cap_3 = ceil(m_3) + 2^32 = 14,755,761,378 (m_3 in section 6.4) |
| d3 cursor table | 2^24·16 |
| d4 counters | 2^48·8 |
| pointer array | 2^24·cap_3·8 |
| candidate list | 2^20·16 |
| request and reply buffers (2^21 requests; 2^21 replies of up to 8 records) | 2^21·64 + 2^24·64 |
| relay buffer | 2^22·64 = 2^28 |
| windows (record slot and flag) | 2^37·64 = 2^43 |
| event queue | 2^37·8 = 2^40 |
| row and column ports (device registers) | 2·64 |

- **Per node.** The total is 2^63.7542 bytes, below 2^64, which leaves 2^61.33 bytes free.
- **Addresses.** The ports are written at fixed offsets, the event-queue pointer wraps at the ring's end (4.1), and every
  other address the kernels form lies in the regions above.
  A remote item is reached only by a message.
- **Peak.** The peak is `2^72` times the per-node total, 2^135.7542 bytes. It counts every region as live throughout,
  including code, tables, processor state and retained coins.

## 6. Success (no heuristic)

### 6.1 Model

**Generators.**
- Coins are uniform and independent, and Q is uniform. Each `Y_l = [c_l]P` is uniform on G \ {O}.
- Given R_g, each `X_k = R_g + [a_k]P` is uniform, and different k are independent. So the X's are independent uniform points,
  independent of the Y's, and nodes are independent.

**Leaves.** A leaf `L = X_k + s·Y_l` is uniform on G. Two leaves of different nodes are independent.

**b_g.** Given all generator points, each `b_g` is still uniform, and the `b_g` are independent across nodes. For each value
of b_g exactly one tuple `(a_k)` fits the points.

**Γ.** Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² − 2^72·n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`. Then `P(I_α = 1) = P1 = (6(r−1) + 1)/r²`, and

```text
mu = |Gamma| * P1 = 0.494342942...
```

Same-node pairs are left out of W; any collision they give only adds to success.

### 6.2 No twins: Chen–Stein on cross-node pairs

Let `W = Σ_(α ∈ Γ) I_α`.
- **Neighbourhoods.** The neighbourhood `B_α` is the set of pairs that share a generator point with α. I_α depends only on
  α's four generators, which are distinct because they come from two nodes. So I_α is independent of
  `{I_β : β ∉ B_α}` (b3 = 0).
- **Size.** A generator lies in `2q2` (an X) or `2q1` (a Y) leaves, so `|B_α| <= 4(q1 + q2)·S`.

For a unit u, the event `L' = u·L` with `L = X + s·Y` and `L' = X' + s'·Y'` is the linear relation
`X' + s'·Y' − u·X − u·s·Y = O` in the generators' logarithms. `I_α` is at most the sum of these 6 events. Take α ≠ β in Γ that
share a generator, and units u, v.

**Claim: the two relations are not proportional.** Proportional relations have the same support. Each node then contributes
exactly one X and one Y to it. The leaf of β in that node must use both, so it is the same incidence as α's leaf there. The
two relations can then differ only in signs, and comparing coefficients forces β = α.

So some 2 × 2 minor is invertible. Given all other generators, the two relations fix two generators, each uniform on at least
r − 1 points. So each pair of events has probability at most `1/(r−1)²`, and `E[I_α I_β] <= 36/(r−1)²`.

The Chen–Stein method (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989), Theorem 1) then gives

```text
|P(W = 0) - e^-mu| <= b1 + b2 <= |Gamma| * 4(q1 + q2)S * (P1^2 + 36/(r-1)^2) < 2^-96.87.
```

A single pool loses this. In one complete grid, the pairs `{X+Y, X'+Y'}` and `{X−Y', X'−Y}` give the same relation for u = 1.
These are the twins of section 8. Across nodes they cannot occur.

### 6.3 Usefulness and output

**Detection.** A cross-node colliding pair is **detected** if all of the following hold:
- both records are computed exactly and emitted;
- both survive the caps;
- no intruder enters their run;
- the candidate, relay and request caps are not reached first.

A detected run is resolved, and its useful adjacent pair (section 4.3) gives x. The `[x]P = Q` check prevents wrong output.

**Usefulness.** The pair is useful iff `n' ≠ u·n`, where n and n' are the b of the two different nodes. Given all points, this
fails with probability exactly 1/r.

So

```text
P(fail) <= P(W = 0) + E #{alpha : I_alpha = 1, alpha useless or not detected}.
```

### 6.4 Defect rows

**Reduction to per-leaf rates.** Take a defect of L that depends only on L's node, such as a batch rare value, a window, the
wrong index or the meter. The node of L' is independent of it and L' is uniform given it, so
`P(I_α ∧ defect) <= (6/r)·P(defect)`. Summed over Γ this is at most `μ(1 + 1/r)` times the per-leaf rate.

**Caps and intruders.** For these, condition on α's four generators.
- **Uniformity.** Every other leaf except the ± partners of α's two incidences is then still uniform on at least r − 1 points,
  and at most `6·2^(256−f)/(r−1)` of that mass falls on given values of f key bits.
- **Fixed records.** The partners and L' are fixed, and each cap leaves room for them:
  - the hop-1 and d3 caps allow 3 fixed competitors;
  - the candidate list allows 4 fixed records;
  - a source allows 8 fixed requests;
  - a relay allows 16 fixed records (α's two replies).
- **Defective batches.** The leaves of a defective batch elsewhere are not uniform. The whole-run row "a computational defect
  in any batch" excludes them.

**Second moments.** Covariances between the indicators of two other leaves are tiny:
- leaves with no common generator are independent: 0;
- the ± pair of one incidence is uniform on ordered pairs of distinct points: at most 0;
- leaves sharing a free Y are each uniform given it: 0;
- leaves sharing a free X are independent given it and each uniform on G minus X: at most `1/(4(r−1)²)`;
- the at most `2(q1 + q2)` ± pairs through a generator of α: each at most the mean of one leaf.

These add less than 1 to every count used here. So each count's variance is at most its mean plus 1, and Chebyshev applies.

| row | bound (per leaf of the pair, per pair, or for the whole run) |
| --- | --- |
| useless collision | 1/r per pair |
| batch rare values, X = O, dx = 0 (3.5) | < 2^-122.6 |
| window drop (3.3) | < 2^-54.29 |
| wrong class index (3.3) | < 2^-157 |
| product meter: `Var <= (2/9)n0 + 1`, slack `2^42 − 2q2` | < 2^-31.47 |
| hop-1 buffer: mean <= `m_c = 2^36·n0·6·2^220/(r−1) + 1`, cap_c = 6n0 + 2^50 | < 2^-42.71 |
| d3 bucket: mean <= `m_3 = S·6·2^160/(r−1) + 1`, cap_3 = ceil(m_3) + 2^32 | < 2^-30.71 |
| intruder: `S·6·2^112/(r−1)` | < 2^-14.71 per pair |
| candidate cap (Markov on other candidates at the home) | < 2^-33.43 |
| request cap (Markov on other requests at the source) | < 2^-36.01 |
| relay cap (Markov on other records at a relay; a reply counts 8) | < 2^-31.26 |
| a computational defect in any batch (whole run, `2^72·q1` batches) | < 2^-23.97 |
| coin range rejection twice, `x_Y >= p − C0` (whole run) | < 2^-120 |
| RNDR failure (whole run) | 0 under the rules' ideal coins (section 7) |

The means use the crude rule that f fixed key bits cost at most `6·2^(256−f)/(r−1)`. It needs no distribution of K beyond
"one key per class, six points per class". That is why `cap_c` is 6n0 and not n0.

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.39001.
```

That is q1 = 108,966,605. The left side is then 0.3900100037, using an upper enclosure of e^-μ and exact Fractions. The target
0.39001 leaves 10^-5 of slack above the rules' 0.39.

## 7. Ledger and claim

All rows are scalar instruction counts. The analytic allowances have explicit content:
- **Generator construction, 2^21.** It covers:
  - a complete double-and-add in Jacobian coordinates;
  - 256 doublings (2M + 5S + 10 additions);
  - at most 256 mixed additions (8M + 3S + 10 additions);
  - exceptional-case tests;
  - two inversions of 23,646 with the affine conversion;
  - the table write.

  With a product at 88 (the single product of 3.1), a square at 64 and an addition at 9, a doubling is 586 and a mixed
  addition 986. Without the tests and loop control (a few thousand) the total is 450,752 = 2^18.78; the allowance is above
  4.5 times the whole.
- **Coins, 64 per scalar.**
  - Each coin word is `MRS Xt, RNDR` (FEAT_RNG), with a B.EQ that retries once on its failure flag; one STP stores two
    words.
  - A draw is 4 words, at most 18 instructions with retries, plus a range test of at most 7.
  - A scalar takes at most two draws: at most 50.
  - In this machine RNDR is the coin source. It returns 64 fresh uniform bits, the algorithm coins of the rules' probability
    space, and never reports failure. The B.EQ is charged for hardware that can report one. In the model it is never taken,
    so coins add no failure row.
- **Node set-up, 2^20.**
  - The shipped kernels are 1,793 instructions: 892 in leaf.s, 123 in route.s and 778 in secp_fe.s.
  - With generator construction, the resolver and phase control, the program is at most 2^16 bytes, or 2^11 records. Its
    delivery over 2 hops at 40 per record is 2^17.3.
  - The rest covers port, buffer and pointer state.
- **Header constants, 148** per node.
- **Batch set-up, 64.** Copying X into H, the index word, acc = 1, pointers, the product count, the meter test and the batch
  loop.
- **Inversion call, 16.**
- **Candidate, 2^25.** It covers:
  - two requests and two replies of at most 8 records over 2 hops each, at 40;
  - two exact additions with inversion;
  - the six-unit test;
  - the solve mod r: an inversion mod r by 384 products of at most 200 instructions;
  - `[x]P = Q`.

  In total 554,530 < 2^20.
- **Request served, 2^10.**
- **Barrier, 2^82.** Up to 4 messages per node at 40, plus each row and column leader's 2^36 handler steps at 11: below
  2^79.51.
- **Output, 2^20.**

These allowances are 2^-14.720 of the total. The hop and local-pass margins of 4.1 are part of the per-record prices.

| row | log2 instructions |
| --- | ---: |
| program delivery, node and port state: 2^72·2^20 (P) | 92.0000 |
| coins: (q1 + q2 + 1)·64 per node (P) | 105.8575 |
| Y generators: q2·2^21 per node (P) | 120.0000 |
| R_g and X generators: (q1 + 1)·2^21 per node | 119.6993 |
| header constants: 148 per node (P) | 79.2095 |
| window flags cleared: 7·2^35 + 1 per node (P) | 109.8074 |
| Q broadcast: 2·(2^36 + 2^72) records at 40 | 78.3219 |
| d4 counters zeroed: 2^48 + 1 per node (P) | 120.0000 |
| d3 cursor table: 5·2^24 + 1 per node (P) | 98.3219 |
| forward: q1·(105·q2 + 3) | 132.4136 |
| inversion: q1·(23,646 + 16) | 113.2296 |
| backward: q1·(523·q2 + 5) | 134.7300 |
| key pass without emits: q1·(53·q2 + 11) | 131.4272 |
| class-index products: cap_j·83 per node | 132.4897 |
| batch set-up: q1·64 | 104.6993 |
| hop 1: S·40 | 132.0212 |
| hop 2: S·40 | 132.0212 |
| d4 prefix sums: 14·2^46 + 1 per node | 121.8074 |
| local pass: S·20 + 7·2^24 + 1 per node | 131.0212 |
| scan: S·7 + 4 + 4·(2^20 + 1) per node | 129.5067 |
| candidates: 2^20·2^25 per node | 117.0000 |
| requests served: 2^21·2^10 per node | 103.0000 |
| barriers: 16·2^82 | 86.0000 |
| output | 20.0000 |
| **total R** | **135.673461** |

**Totals.**
- `R = 69,467,457,275,035,229,547,858,820,770,057,608,495,104` instructions, 502.9079 per leaf.
- `W = R/1000 = 2^125.7076768`.
- **The claim, 125.707677, is the least multiple of 10^-6 with 2^claim >= W.** The certificate encloses 2^c by an exact series
  for e^(f·ln 2) with interval bounds on ln 2.

**Other reported numbers.**
- **Executed counts.** If the hops and the local pass were priced at their executed counts (17.5, 29 and 9, the event path
  included) instead of 40, 40 and 20, the work would be 458.4079 per leaf and the bound 125.574015. The claim uses the prices.
- **Preprocessing.** The rows marked P come to 2^111.0349 units, and are included in R.
- **Memory.** 2^135.7542 bytes. No advice.

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "trials": 6000, "control_trials": 1000, "curve_trials": 100,
"seconds": 15}`. The schema allows `seconds` from 5 to 15.

**Budget.** The class-rule checks and the toy trials are limited by the budget. The other checks have fixed sizes.
- Each budget-limited count first runs a guaranteed minimum, or its requested count if that is smaller.
- It then continues while its share of the remaining budget lasts, so the requested counts are caps.
- The minimums are:
  - 2,000 random and 600 near-tie class-rule inputs (the near ties come from a second random stream if the budget stops the
    main loop first);
  - 400 multi-node trials;
  - 100 trials for each control;
  - 10 curve checks.
- A run therefore takes about the larger of `seconds` and the fixed work plus the minimums.

**Run times.** In our sandbox the declared input ran in 14.6 s in each of two runs. A second CPU-bound process then shared
its core, which made a fixed workload 2.2 times slower. Under that load it ran in 14.6 s, and with `seconds: 5` in 7.9 s.
Peak memory was 12.5 MiB.
The counts a run reaches vary with the runner, so the larger counts are quoted only from the census below.

**A. Production identities and leaf kernels.**
- The identities: β³ = 1, λ³ = 1, [λ]G = (βx, y), 7 a non-residue, and the five window forms nonzero.
- The KEY model against the exact class rule on at least 2,000 random and 600 near-tie inputs (up to 20,000 and 6,000 as the
  budget allows): no mismatch outside windows.
- `fwd`, `bwd` and `keys` on production batches of 1, 2, 5, 9 and 16 incidences.
  - Every slot matches affine addition.
  - A device hook captures the messages written to the fixed port, and each matches the exact key, digits and index words.
  - Each call returns through an unchanged sentinel x30, which `a64.py` asserts at the outermost RET.
- Counts:
  - 105 and 523 per incidence, with fixed parts 3 and 5;
  - `keys = 11 + 66·q2 + 83·(products) − (7 or 8 per dropped leaf)`.
- A batch with constructed window inputs drops exactly those leaves.
- fe_inv: 23,646 instructions, correct.
- Program size: 892, 123 and 778 instructions.

**B. Transport on a 4 × 4 machine.**
- **Records.** 640 records with planted equal keys run through the routing kernels with small digit masks:
  - each node runs `wclear`, `init3` and `zero4` before hop 1;
  - hop 1 arrives in rounds, at most one record per window per round, as backpressure allows. Each round's events go to
    `h1recv`;
  - hop 2 leaves through the row port (device hook, `h2send`) and arrives the same way into `h2recv`;
  - `prefix4`, `pass4` and `scan` run after hop 2.
- **Results.**
  - Every flag is clear after every round.
  - The candidates equal a Python model of the same stable order exactly.
  - Every record reaches its home tagged with its source row and column.
  - The executed totals equal the per-call formulas, evaluated with the certificate's constants.
- **Full paths.**
  - `h1recv` with a full buffer: 62 = 2 + 3·11 + 3·9.
  - `h2recv` with a full bucket: 46 = 2 + 20 + 2·12, with d4 counted once.
  - The scan with a one-entry list, which ends at the second equal pair: 36, as modelled.
- **Intruders.** The run reports groups minus candidates. It is nonzero because 8 toy sort bits let intruders split groups,
  the row of 6.4.

**C. Toy pipeline at q25** (q = 33556867, r = 33563419, j = 0).
- **Multi-node grid.** 74 nodes of 4 × 4 grids give 2,368 leaves and μ_cross = 0.494435, near the production 0.494343.
- **Single-pool controls.** One complete 28 × 42 grid, either with b per X or with production's shared b.
- **Log form.** Trials use the log form of the same search: the class of `[t]G` is `t^6 mod r`, since the six units are the
  sixth roots of unity mod r.
- **Curve check.** A curve pipeline runs x-only leaves by batched inversion, `min(x, βx, β²x)` keys, resolution with points
  and `[x]G = Q`. On the same coins it gives identical outcomes.

The run reports measured success with its standard error next to:
- the Poisson value `1 − e^-μ`;
- the twin model `1 − e^(−5μ/6)` (u = 1 and −1 come in twins);
- the shared-b model `1 − e^(−3μ/4)` (u = 1 is also useless with one b);
- the Chen–Stein error bound of section 6.2 evaluated at the toy's own size, 0.0134. Production's is below 2^-96.87.

**Census.** `census-toy.json` holds a longer off-runner run of the same code (section 9). Its results:

| run | trials | success (± se) | reference |
| --- | ---: | ---: | --- |
| multi-node grid, all pairs | 100000 | 0.39182 ± 0.00154 | cross-node μ plus 3/4 of the same-node μ (0.500997 − 0.494435): 0.39308 |
| multi-node grid, cross-node pairs only | 100000 | 0.38849 ± 0.00154 | Poisson 0.39008; toy Chen–Stein bound 0.0134 |
| one pool, b per X | 20000 | 0.34230 ± 0.00336 | twin model 0.33759 |
| one pool, shared b | 20000 | 0.31400 ± 0.00328 | shared-b model 0.30974 |
| curve pipeline against log form | 1000 | 1000 identical outcomes | — |

The toy's cross-node success neither establishes nor contradicts production success. In the census it lies 1.03 standard
errors below the Poisson value and well inside the toy's Chen–Stein bound. A declared run has far fewer trials, so its
figure should be read with its own standard error.

**D. Certificate.** The run imports `certificate.py`, reports PASS, and checks these executed counts against parts A and B:
- `FWD_INC`, `FWD_FIX`, `BWD_INC`, `BWD_FIX`;
- `KEY_INC`, `KEY_FIX`, `KMUL_PATH` (through the keys formula);
- `INV`;
- `H1RECV`, `H1RECV_DROP`, `H2SEND`, `H2RECV`, `H2RECV_DROP` (the event path inside);
- `ZERO4_PER4`, `PREFIX4_PER4`, `INIT3_BKT`, `WCLR_PER4`;
- `PASS4_REC`, `PASS4_BKT`, `SCAN_REC`.

The analytic allowances of section 7 are not executed.

The toy and the 4 × 4 run are evidence at their own scope only. The production success is the proof of section 6, and the
production cost is the ledger of section 7.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1. Only the two endpoints of a hop execute instructions for it, and every one
  of them is charged. On a machine whose links pass through intermediate processors, their instructions would be added; that
  machine is not claimed.
- **Analytic allowances.** These are bounded analytically and not executed:
  - generator construction and coins;
  - node set-up, header constants, batch set-up and the inversion call;
  - candidates, requests, barriers and output.

  Together they are 2^-14.720 of the total. Each is at least the itemized work of section 7 it covers: generator
  construction by a factor above 4.5, candidates above 60, barriers above 5.6.

  The handler's dispatch loop (4.1) is itemized, not executed. The hop prices cover it.
- **Executed per-record work.** It is shipped assembly run in `kernels/a64.py`. The interpreter is package F's, extended by:
  - W registers;
  - BL/RET nesting, with a check that the outermost RET returns through an unchanged x30;
  - CINC and CSET;
  - register ORR, AND and EOR;
  - LSL;
  - a device-register hook for stores.
- **Not executed at production scale.** The production attack, the machine and its links are not executed.
- **Provenance.** It is not a premise:
  - package F (8f733826): RF, RA, W2P, FS's CSEL form, fe_inv, the F1 canonical row and its error analysis, a64.py;
  - Galbraith, Wang and Zhang (ePrint 2015/605): X + Y and X − Y from one inverse, x-coordinates only.

  Arratia–Goldstein–Gordon (1989) is cited for the Chen–Stein bound.
- **Census file.** `census-toy.json` comes from `experiment.py` with input
  `{"seed": 20261030, "trials": 100000, "control_trials": 20000, "curve_trials": 1000, "seconds": 1500, "kernels": false}`.
  This run is off-runner: its `seconds` is above the hosted schema's maximum of 15. With that budget every count reaches its
  requested value.
