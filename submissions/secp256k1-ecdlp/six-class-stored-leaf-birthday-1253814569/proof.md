# Stored x-only leaf birthday search on secp256k1, six leaves per inverse (125.3814569)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.3900001 over the challenge and fresh algorithm coins; on every run its total work is at most
`time_log2 = 125.3814569 = 1253814569/10000000` target-group-operation equivalents (C = 1000 scalar AArch64
instructions per unit). Peak memory is at most 2^135.9988 bytes; preprocessing, at most 2^110.2061 units, is included;
there is no advice. `certificate.py` checks these declared figures (`DECLARED`): the work, memory and preprocessing
figures are the least multiples of 10^-7, 10^-4 and 10^-4 at or above their exact values, the work per leaf 401.1366 the
least multiple of 10^-4 at or above R/S, and 0.3900001 is below the success bound.

**Reviewer's guide.** The bound rests on five things.
- **Six leaves per inverse (3.1).** On secp256k1 (j = 0), `D = y_Y^2 - y_X^2 = x_Y^3 - x_X^3` is the product of the
  three denominators `x_Y - β^t x_X`, t = 0, 1, 2. So one batched inverse of D gives the slopes of all six points
  `φ^t(X) ± Y`, and the backward pass computes their x-coordinates, the leaves, at 1453 instructions per pair (X, Y).
  The run checks every leaf of its production batches of 1, 2, 5 and 9 incidences against exact affine arithmetic.
- **No uncharged transfer or bookkeeping (2.1).** A node's per-record steps are instructions of its own processor,
  charged at their executed counts: the emit inside `keys` (12 per record) at hop 1, `h2s` (11.25) at hop 2, and `bkt`,
  `pass4` and `scan` (21.5), the unrolled loops' control included (section 4). Each message also pays 8 at each hop for
  its link: its delivery as `dlv` in `route.s` (6), which the run executes for every message of its 4 × 4 machine, and
  its port's in-flight count, which DSB waits on (2). So the hops are charged 20 and 19.25. Links keep no flags, rings or
  counters. Every address, cursor, slot count and cap test is a counted instruction of the sender or of the node that
  later reads the data.
- **Success (section 6):** Chen–Stein on cross-node leaf pairs with b3 = 0 (6.2; its one non-routine step is the
  non-proportionality claim, for six leaf types), usefulness (6.3), the defect table (6.4), whose batch
  row comes from the arithmetic-site table of 3.4 and whose cap rows bound every pair of leaves, dependent ones included,
  and the final inequality (6.5). `certificate.py` evaluates the terms with exact Fractions and prints the bound
  0.3900001034 (the run's `success_lower_bound`).
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
extends this to the three twists of X. Keyed by a six-class canonical x, each leaf goes in two bulk-synchronous hops to
a home node, where equal keys give the logarithm. This is an analytic bound, not a practical break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each with one scalar AArch64 processor and a private memory
below 2^64 bytes. The ledger counts the instructions the processors execute and, per message, its delivery and its
in-flight count (below); storage and wiring carry no price. The processor implements the base scalar AArch64
instruction set with FEAT_RNG (`MRS Xt, RNDR`, section 7) and FEAT_LSE (`STADD`, section 4); every counted step is one
such instruction, whatever memory accesses it performs.
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

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 145286863` and
`q2 = 2^25`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are the
x-coordinates of `φ^t(X_k) + s*Y_l` for all k, l, t in {0, 1, 2} and s = ±1: `n0 = 6*q1*q2` per node and
`S = 2^72*n0 = 2^126.699292` in all. As φ = [λ], leaf `φ^t(X_k) + s*Y_l` is `[m]P + [n]Q` with `m = λ^t*a_k + s*c_l` and
`n = λ^t*b_g`. A point Z ≠ O has the class `{±Z, ±φ(Z), ±φ²(Z)}`, and the **key** `K = min(x, βx mod p, β²x mod p)`
(limbs K0..K3) is shared by the class and differs between classes. A record is `[K2, K3, w2, w3]` with `w2 = i | k<<36`
and `w3 = d2 | c<<36 | l<<39`, where `c = 1 + 2t + [s = -1]` in 1..6 is the leaf type; its digits are `d1 = K0 mod 2^36`
(hop-1 row), `d2 = K1 mod 2^36` (hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3 mod 2^48`. At its home a record's source
is `(w2 mod 2^36, j)`, where j is the index of the slot it arrived in (2.3).

### 2.3 Phases

A pair slot holds `cap_p = 2^22 - 2^12` records of 32 bytes. A node in row i writes at hop 1 only into slot i of R1, at
`R1 + 32*cap_p*i`, at each node of its column, and a node in column j writes at hop 2 only into slot j of R2 at each
node of its row. Each port's cursor table holds, per destination, the sender's cursor in its slot there and the slot's
end; `tinit` sets each of its 2^36 entries to the slot's bounds.

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators with their table
   entries, header constants, the two cursor tables, the d4 counters (`zero4`) and the d3 cursor table (`init3`). Then Q
   is broadcast in 2 hops.
1. **Hop 1:** the batch kernels (section 3) run for each X_k, and the emit writes each leaf through the column port to
   `(d1, j)` at the cursor of d1, then advances it; a leaf whose slot is full is dropped. Then `fin` writes the final
   cursor of each row d into entry i of d's table E1, and a barrier follows.
2. **Hop 2:** `h2s` reads slot i of R1 up to `E1[i]`, for every i, and writes each record through the row port to
   `(own row, d2)` by the same cursor rule; `fin` writes each home's final cursor into its `E2[j]`; a barrier follows.
3. **Grouping:** `bkt` reads each slot of R2 up to its E2 entry and appends the record's address and d4 to its d3 bucket
   of `cap_3` entries, or drops the record if the bucket is full; a stored record counts its d4. `prefix4` forms the d4
   cursors, `pass4` writes each entry's address into the pointer array at its d4 cursor, bucket by bucket in d3 order,
   and `scan` lists adjacent records with equal K2 and K3 as candidates, up to 2^20; a full list ends it.
4. **Resolution,** in four sub-phases: a home `(h_r, h_c)` sends two requests per candidate to the sources `(g_r, g_c)`
   via the relay `(g_r, h_c)`, and replies of at most 8 records (one leaf's scalars and points) return via `(h_r, g_c)`.
   Each sub-phase is a count exchange: senders write each receiver their count and the epoch; after a barrier each
   receiver returns each sender's base, in sender order, or a drop mark if its records would pass the cap (2^21 requests
   at a source, 2^22 records at a relay); after a barrier senders write the records at their bases; a barrier ends it.
   The home recomputes `L = [m]P + [n]Q` and `L'` from the records' leaf types, finds the unit u with `L' = uL` if any,
   and, if `n' - u*n ≠ 0`, outputs `x = (u*m - m')/(n' - u*n) mod r` after checking `[x]P = Q`.

**Barriers:** 21 (program delivery 2, Q broadcast 2, the ends of phases 0 to 3 4, phase 4 12, output 1), within the
allowance of 32. **Caps:** products (3.2), pair slots, d3 buckets, candidates, requests, relay records and coin retries
are capped or bounded in total, so the ledger holds on every run.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's. With `D_l = y(Y_l)^2 - y(X)^2`, the forward pass stores `D_l` and `acc_(l-1)` in entry
l and forms `acc_l = acc_(l-1)*D_l` (`acc_(-1) = 1`); fe_inv gives `R = 1/acc_(q2-1)`; the backward pass, for
`l = q2-1..0`, sets `inv = R*acc_(l-1)`, `R <- R*D_l`, `u = inv*(y_Y - y_X)` and `v = inv*(y_Y + y_X)`, and for t = 0, 1,
2 writes

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
the six leaves overwrite the 256-byte entry in the order `c - 1 = 2t + [s = -1]`. The arithmetic is F's row schedule, with
the row-shift and fold carries taken inline rather than by branches; its only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` walks the entries downwards and takes each entry's six leaves: KEY (3.3) gives the class index j or drops the
leaf; for j = 1, 2, KMUL computes `K = x*β^j` and counts one product; the **emit** forms d1, loads the cursor and slot
end of d1 from the table at port + 64, and drops the leaf if the slot is full (5 instructions); otherwise it forms w3,
writes the message with three STPs and stores the advanced cursor (7). The leaf pointer, advanced 32 bytes per leaf by a
post-indexed LDP, gives w3's high bits: entries are 256 bytes at a 2^34-aligned base, so `pointer << 31` is
`c<<36 | l<<39`. As j is uniform on {0, 1, 2} for a leaf other than O, a node starts a batch only while its product count is at
most `cap_j - 6q2`, `cap_j = ceil(2*n0/3) + 2^42`. Two leaves of one incidence are an invertible linear image of (X, Y),
and leaves sharing one generator are independent given it up to one excluded point, so a node's product count has
variance at most `(2/9)n0 + 1` (6.4's meter row, while no batch of the node is defective). Executed counts per incidence
(fixed part in parentheses): `fwd` 105 (3), `bwd` 1453 (5), `keys` 224 for six leaves, including their six emits of 12
(11), plus 83 per product; fe_inv (`secp_fe.s`) is 23646. Kernels may overwrite x0-x29 and the flags and preserve sp;
`bwd`, `keys` and fe_inv save and restore x30, which `fwd` and the transport kernels do not write; the driver keeps no
registers across calls and loads the product count into x18 for `keys`. `bwd` keeps R in x25..x28 from pass to pass;
`bwd` and `keys` end with a TBZ on bit 33 of the entry pointer: the 2^25 entries of 256 bytes fill the lower 2^33 bytes
of their 2^34-aligned block, so bit 33 is clear at every entry and set one entry below the block.

### 3.3 The class index (macro KEY, 21 instructions)

This is F's canonical-index row, except that a window drops the leaf. `G1..G5` are limbs 1..5 of `G = floor(β*2^384/p)`,
F0 is the top limb of x, `F1 = Σ_i [lo(x_i*G_(5-i)) + hi(x_i*G_(4-i))] mod 2^64`, `σ = F0 + F1`, `d = F0 - F1` and
`e = min(F0, F1) + σ`. Then j = 2 if e carries, else 1 if F0 ≥ F1, else 0. The windows are (W1) `F1 >= 2^64 - 16`, (W2)
`σ >= 2^64 - 17`, (W3) `e mod 2^64 >= 2^64 - 33` and (W4) `d < 16`.

**Exactness.** Let `V_i = 2^64*v_i/p` for `v0 = x`, `v1 = βx mod p` and `v2 = β²x mod p`, which are distinct and sum to
p or 2p. Then `F0 = V0 - a` with `0 <= a < 2`, and `a < 1` unless `x mod 2^192 > 2^192 - C0`. Bits 320..383 of `x*G`
equal `F1 + c` and `floor(V1 - ε)` mod 2^64, with `0 <= c <= 8` (the terms left out of F1 sum below 8.01*2^320 mod
2^384) and `0 <= ε < 2^-64` (as `0 <= β*2^384/p - G < 1`). So outside (W1) `F1 = V1 - b` with `0 <= b < 10`, and outside
(W1) and (W2) `2^64 - σ = V2 + a + b`. Hence j = 2 gives `V0, V1 >= V2`; j = 1 gives `V0 - V1 > 16 - 10` and
`V2 - V1 > 34 - 2 - 20` by (W4) and (W3); j = 0 gives `V2 - V0 > 34 - 4 - 10` and `V1 - V0 >= 1 - a > 0` when `a < 1`.
Outside the windows j is therefore `argmin v_j`, except when j = 0 and `x mod 2^192 > 2^192 - C0`.

**Rarity.** (W1) puts v1 within 2^196 of 0 or p; (W2) and (W4) put `(1±β)x mod p` within 2^197 of 0; (W3) puts `(2+β)x`
or `(1+2β)x mod p` within 2^198 of 0. As these maps are bijections, at most `13*2^197` values of x hit a window and at
most `2^64*C0` give the wrong index; a uniform point's x takes each value with probability at most `2/(r-1)`.

### 3.4 Arithmetic sites

**Steps that can lose a carry.** A product on the rows of a weak b (below 2^256) uses `B_0 = b` and row shifts
`B_(k+1) = SHX(B_k)`; RF and RA sum `Z = Σ a_k*B_k < 2^322` exactly (each high word is at most 2^64 - 2), and a fold
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. SQ5 forms the square of a weak a as
`V = L + H*C0` exactly (`a² = L + 2^256*H`; V's top limb is at most C0 + 1), AC4 may add a constant below p (top limb at
most C0 + 2), and F5 folds V the same way. SB6 forms `Z' = Z - s'' + 2^256` exactly (top word at most 4). Only these
steps drop a carry or a borrow: SHX's last ADCS, the final ADC of FDX and of F5, and FS when its subtrahend is a fold's
output at or above p. Without a dropped carry each returns a value below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`; the fused fold of
  `Z' < 2^322 + 2^256` likewise needs the leaf in `[C0, (2^66 + 2)*C0)`.
- **square fold**: lost iff its sum is at least 2^256, which needs the true square in `[C0, C0² + 2C0)`, or Q_t in
  `[C0, C0² + 3C0)` with the constant added.
- **FS** `a - b` with `a < p` and `b < 2^256`: wrong iff `b - a > p`, which needs `b >= p`; the subtrahends Q_0 and Q_1
  are fold outputs, at or above p only when the true value is below C0.

A fold returns a value `>= p` only when the true value is below C0. That matters only where a canonical value is
needed: the leaf entering KEY, the emitted key K, and FS's subtrahend (above). Every other operation is exact on its
inputs: RF and RA; SB6; KMUL's sum over the constant rows of β^j (below `4*2^64*p < 2^322`; only its FDX can drop a
carry); FS and FADDW on canonical inputs; FADDW of a canonical value and `x_Y + C0` (the sum is below `p + 2^256`, so
one conditional subtraction of p leaves it below 2^256); KEY's integer arithmetic (3.3 covers its index); the emit; and
F's `fe_mul_f` and `fe_sqr_f`, which fold the 512-bit product, fold its fifth limb again and add C0 exactly when that
last sum carries past 2^256, so fe_inv (255 squarings, 15 multiplications) is exact for every input below 2^256.
Generator construction, the table entries, the header values, the resolver and `[x]P = Q` run on `fe_mul_f`,
`fe_sqr_f` and fe_inv, each coordinate brought below p (an ADDS chain of C0 detects `t >= p`, CSELs keep `t - p`) before
it is subtracted, compared or stored, with exact arithmetic mod r for the solve. So they have no sites; the tables hold
exact canonical values, and `x_Y + C0 < 2^256` as `x_Y < p`.

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
| keys: KMUL folds of the six leaves | 6q2 | `Y_l`: the leaves | 2*2^66*C0 | 2 |

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

## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (12; 5 if the slot is full), and `dlv` (6) | 18 | 20 |
| hop 2: `h2s` (11; 7 if the slot is full; loop 0.25), and `dlv` (6) | 17.25 | 19.25 |
| bucket pass `bkt` (11; 6 if the bucket is full; loop 0.25) | 11.25 | 11.25 |
| local pass `pass4` (4; loop 0.25) | 4.25 | 4.25 |
| `scan`: 6 per record (paired comparisons), plus 4 per listed candidate | 6 | 6 |
| **total** | **56.75** | **60.75** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`, `dlv`'s included; the charged count adds
the in-flight count of 2.1, 2 per message per hop. `h2s`, `bkt` and `pass4` take a slot's or bucket's records in
blocks of 1, 2 and 4, as the low bits of its count say, then 8 at a time; the 8-record loop's CMP and B.LO, 2 per pass,
are charged as 0.25 per record (a slot of n records runs floor(n/8) passes), and the tests and set-up, 11, 12 and 11
per slot or bucket, are per-node rows. `h2s` adds 64 to its port argument and addresses the register block at -64, -48
and -32, so its cursor entry needs no offset and its advance is one STADD of 32. `bkt` steps its record pointer with a
pre-indexed LDP and counts d4 with one STADD at the counter's address; `pass4` moves its d4 cursors with register-offset
LDR and STR. `scan` runs its comparisons in pairs: a prologue handles an odd first comparison, then each pair of
comparisons closes with one CMP and B.LO, so a record costs its pointer's LDR, its key's LDP, CMP, CCMP, B.EQ and half
the end test. A hop's receiving node executes nothing on arrival; it reads the records in the next phase, with the
LDPs of `h2s` (hop 1) or `bkt` (hop 2). Per-node costs are rows of section 7, and every address a kernel loads or
stores lies in the regions of section 5.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 key bits. A class forms one run unless an **intruder**, a record outside it with
the same 144 bits, sits inside (6.4); a run's adjacent pairs are all listed unless the candidate list is full. Writing
the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`, an equivalence relation, so a
run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums the regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.0085 used) | code, stack, the five barrier words (one reservation granule each), header, X table (eight header values and y² per X, 288 bytes), Y table (160 bytes per Y), scalars, batch entries (2^33 bytes at a multiple of 2^34), the two ports with cursor tables of 2^40 bytes each, E1 and E2 (2^41 each), the d3 cursor table, the d4 counters (2^51) and the phase-4 buffers and tables (`FIXED` in `certificate.py`) | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | hop-1 pair slots | phases 1-2 |
| the same range, from phase 3 | d3 buckets `16*2^24*cap_3` and pointer array `8*2^24*cap_3`, together 2^62.3655 | phase 3 |
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
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494329715` (the run's `mu`).
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

A cross-node colliding pair is **detected** if both records are computed exactly and emitted, both survive the caps, no
intruder enters their run, and the candidate, request and relay caps are not reached first. Its run is then resolved, a
useful adjacent pair (section 4) gives x, and `[x]P = Q` is checked before output. The pair is useful iff `n' ≠ u*n`,
where `n = λ^t*b_g` and `n' = λ^(t')*b_(g')` for its two nodes; given all points, this fails with probability exactly
1/r. So `P(fail) <= P(W = 0) + E #{α : I_α = 1, α useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (a batch defect of 3.4, a window, the wrong index, the meter) is a function of node g's
coins, and L' is uniform given them, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs at most `μ(1 + 1/r)`
times its per-leaf rate.

**Caps and intruders.** Condition on α's four generators. A leaf with a generator outside them, a **free** generator,
is then uniform on at least r - 1 points, so at most `6*2^(256-f)/(r-1)` (`KB(f)`) of its mass falls on given values of
f key bits ("one key per class, six points per class"). The 12 leaves of α's two incidences are fixed; the six leaves of
an incidence through one of α's generators and a free one are functions of that free generator.
- **Pair slots and d3 buckets** count single leaves. A pair slot has mean at most `m_p = n0*KB(36) + 1` from the free
  leaves (at hop 1 it holds one node's leaves with a given d1, at hop 2 one column's leaves with a given (d1, d2), the
  same bound), and a d3 bucket `m_3 = S*KB(96) + 1`; each cap leaves room for the at most 12 fixed leaves. Two leaves'
  indicators have covariance 0 without a common free generator or with a common free Y; at most `1/(4(r-1)²)` with a
  common free X (at most `KB(f)/(r-1)` when one of them also uses α's Y); at most `1/(r²(r-1))` times the two cells'
  sizes for two leaves of one incidence with both generators free, whose pair is uniform on r(r - 1) pairs (an
  invertible linear image of (X, Y), with Y ≠ O); and at most one leaf's mean for two leaves of one incidence through
  α's generators, at most 60(q1 + q2) ordered pairs. So a count's variance is at most its mean plus
  `1 + 60(q1 + q2)*KB(f)`, and Chebyshev bounds the drop of a given record by `(m + 1 + 60(q1 + q2)KB(f))/(c - 12 - m)²` for a cap c
  and mean bound m (`d_p`, `d_3`).
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
  `other_h <= ρ(S²*KB(72)*KB(200) + (2S + 20)*KB(200))` (pairs with equal 200 key bits at the home; L and L' are
  there); a source receives at most `other_g <= 2ρ*n0*S*KB(200)` requests; a relay holds at most
  `other_v <= 18ρ((S/2^36)*S*KB(36)*KB(200) + (2S + 20)*KB(200))` records (a request is one record and a reply at
  most 8).
  Markov with the caps less α's own entries (4 in the candidate list, 8 requests, 16 records) gives `d_cl`, `d_rq` and
  `d_rl`. Leaves of a defective batch elsewhere are not uniform; the whole-run row "a computational defect in any batch"
  excludes them.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, D = 0, X = O | `d_batch` | < 2^-126 |
| window drop; wrong class index (3.3) | `d_win`; `d_wj` | < 2^-54.29; in the total |
| product meter, Chebyshev with slack `2^42 - 6q2 - 1` | `d_j` | < 2^-31 |
| pair slot at hop 1, and again at hop 2 (Chebyshev, `cap_p`) | `d_p`, counted twice | < 2^-19.9 each |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `ρ(S + 12)KB(144)` | `lam_int` | < 2^-14.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-33.4; < 2^-35.9; < 2^-31.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26.9 |
| two coin range rejections (whole run) | in `glob` | all of `glob` < 2^-26.9 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.3900001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 145286863, where the left side is at least 0.3900001034, and
the target is 10^-7 above the rules' 0.39.

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 2^20 per point:** a complete Jacobian double-and-add, exceptional-case tests, two inversions
  with the affine conversion, all on the exact arithmetic of 3.4: `fe_mul_f` and `fe_sqr_f` at 115 and 82 instructions
  as executed, plus 4 and 3 for the call as in fe_inv, and an addition or subtraction with its operands brought below p
  at 27; 665,910 for the 256 doublings, 256 mixed additions and two conversions. The scalar's bit tests, the pointer step
  and the return add at most 6 per iteration, 2,048 in all. A Y's table entry (y², 3x², x + C0 and ten STPs, 239) or an
  X's mixed addition `R_g + [a_k]P` with its eight header values, y² and eighteen STPs (2,379) brings the total to at
  most 670,337, within the allowance (1.5642x).
- **Node set-up, 2^20:** a program of at most 2^16 bytes (3,670 kernel instructions, generator construction, resolver,
  phase control), 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD, SUBS, B.NE), plus 2 for its in-flight count and 6 for its delivery: 16.
- **Candidate, 2^21:** two requests and two replies of at most 8 records over 2 hops at 40, two exact additions with
  inversion, the six-unit test, the solve mod r and `[x]P = Q`, on the same exact arithmetic; 769,058 (2.726x).
- **Phase-4 count exchange, 2^41 per node and sub-phase:** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 88 per message (its count,
  the count message, the offset reply and the write, at most 22 each with in-flight count and delivery) for at most
  2^24 messages: below 2^39.6.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 14 (6, plus 2 for the in-flight count and 6 for delivery), DSB and 20 of set-up: 131.
- **Phase control, 2^10 per node:** 11 kernel calls between barriers, each at most 33 (eight argument registers of at
  most four MOVZ or MOVK each, and BL), the E-entry addresses and `scan`'s n: below 2^9.
- **d4 prefix sums, 10 per four counters:** two LDPs, two STPs, four ADDs, SUBS and B.NE; the RET is the row's +1.
- **Header 148, batch set-up 128** (X's eight header values copied into H by 16 LDPs and 16 STPs, y_X² into registers,
  acc = 1, pointers, the index word, the meter test and the loop), **inversion call 16, request served 2^10,
  output 2^20.**
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
| Y generators with their table entries: q2*2^20 (P) | 117.0000 |
| R_g and X generators with their header values: (q1 + 1)*2^20 | 119.1143 |
| header: 148 (P) | 79.2095 |
| cursor tables: 2*(3*2^36 + 1) (P) | 110.5850 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| d4 counters: 2^48 + 1 (P) | 120.0000 |
| d3 cursor table: 5*2^24 + 1 (P) | 98.3219 |
| forward: q1*(105*q2 + 3) | 130.8286 |
| inversion: q1*(23646 + 16) | 113.6446 |
| backward: q1*(1453*q2 + 5) | 134.6191 |
| key pass without emits: q1*(152*q2 + 11) | 131.3623 |
| class-index products: cap_j*83 | 132.4897 |
| batch set-up: q1*128 | 106.1143 |
| hop 1, the emits, in-flight counts and deliveries: q1*120*q2 | 131.0212 |
| hop 2, `h2s`, in-flight counts and deliveries: S*19.25, plus 11*2^36 + 3 per node | 130.9661 |
| final cursors, in-flight counts and deliveries, both hops: 2*(15*2^36 + 3) | 112.9069 |
| bucket pass `bkt`: S*11.25, plus 12*2^36 + 2 per node | 130.1911 |
| d4 prefix sums: 2^48/4*10 + 1 | 121.3219 |
| local pass `pass4`: S*4.25, plus 11*2^24 + 1 per node | 128.7868 |
| scan: S*6, plus 24 + 4*(2^20 + 1) per node | 129.2843 |
| candidates: 2^20*2^21 | 113.0000 |
| requests served: 2^21*2^10 | 103.0000 |
| phase-4 count exchanges: 4*2^41 | 115.0000 |
| barriers: 32*2^82 | 87.0000 |
| phase control: 2^10 | 82.0000 |
| output | 20.0000 |
| **total R** | **135.347241** |

`R = 55408883461773899612282413748807787347968` instructions, at most 401.1366 per leaf (the run's `per_leaf`), and
`W = R/1000` with log2 W ≈ 125.38145681. **The claim, 125.3814569, is the least multiple of 10^-7 with 2^claim >= W**;
the certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds on ln 2. Rows marked P are
preprocessing: their sum over 1000 has log2 ≈ 110.20607 and is declared as 2^110.2061, the least multiple of 10^-4
at or above it (checked).

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "trials": 6000, "curve_trials": 100, "seconds": 15}` (`seconds` 5
to 15); it finds the package files beside itself, under files/ or under /snapshot/files. Budget-limited checks run
guaranteed minimums (2,000 random and 600 near-tie class-rule inputs, 400 toy trials, 10 curve checks), then their share
of the remaining budget. **Part A** runs the leaf kernels on production batches of 1, 2, 5 and 9 incidences, and fe_inv,
against exact arithmetic: every leaf against exact affine addition of `φ^t(X)` and `±Y`, each port message and the
cursor table against an exact model of the emit (production 36-bit digits, sentinel x30), and the KEY model against the
exact class rule. For each site type of 3.4 (`site_checks`: row shift, fold, fused fold, square fold with and without
the constant, FS on a fold output, KMUL), inputs built to meet its predicate deviate exactly when it holds, other
boundary and random inputs are exact, and all match its bit-level model; FS and FADDW on canonical inputs, FADDW of a
canonical and a weak value, `fe_mul_f` and `fe_sqr_f` are exact. A batch with `Y_l = φ(X)` and `Y_m = -X` gets acc = 0
and loses all its other leaves; the last R update is ≡ 1 and unread; the SHX, FDX and F5 invocations per pass match the
table (`sites_match_certificate`). **Part B** runs a 4 × 4 machine through `leaf.s` and `route.s` (`tinit`, `init3`,
`zero4`; `keys` on 4 planted incidences of six leaves per node, d1, d2 < 4, some keys shared; `fin`; `h2s`; `fin`; `bkt`,
`prefix4`, `pass4` and `scan` at each home), with pair slots of 7 and buckets of 8, delivering each message by executing
`dlv` on it (`counts.dlv` is 6 × `deliveries`). Every call's count equals its formula from the certificate's constants;
messages, slots, bucket entries, pointer order and candidates equal an exact model; full slots, a full bucket and a full
list are forced; `unrolled_paths` runs `h2s`, `bkt` and `pass4` on slots and buckets of 0 to 19 records (every block
combination, up to two 8-record passes). **Part C** runs the q25 toy in log form (the class of `[t]G` is `t^6 mod r`,
φ is `[λ]`): 25 nodes of 4 × 4 incidences with six leaves each, and a curve pipeline that computes each incidence's six
leaves from one batched inverse of D by the formulas of 3.1 and must match the log form. **Part D** runs
`certificate.py`, and `counts_match_certificate` compares its constants with parts A and B.

**Census.** `census-toy.json` is an off-runner run of the same code (`seconds` 1500, above the hosted maximum of 15;
each toy count reaches its request). Over 100000 toy trials, success is 0.3953 with all pairs and 0.39027 with
cross-node pairs only, against Poisson 0.38997 at the cross-node μ = 0.494252 and a toy Chen–Stein bound of 0.0407,
and the curve pipeline matches the log form in 1000 of 1000 trials. The toy and the 4 × 4 run are evidence at their
own scope only; a declared run's toy figure has its own standard error.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1, whose links keep no state beyond the message in flight. The sender's
  instructions (including those that form the address), each delivery (as `dlv`) and the ports' in-flight counts are
  charged. Each hop is one link from the sender's port to its destination; a network that forwarded a message through
  further nodes would add their work at each extra hop, and that machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control and phase-4
  exchanges), the in-flight counts, and the production attack, machine and links.
- **Provenance**, not a premise: package F (submission 8f733826) for the row products (RF, RA, and the pair schedule of
  W2S and W2L), FS's CSEL form, fe_inv, the F1 row and its error analysis, and `kernels/a64.py`, whose extensions its
  docstring lists; submissions 3828873f, 6ef3e7b3, d637fb0a, 9798f5f4 and d9e3dbf4 for `h2s`'s and `bkt`'s STADD, the
  paired `scan`, `prefix4`, the 2^20 generator and 2^21 candidate allowances and the 7-decimal claim; Galbraith, Wang
  and Zhang for X + Y and X - Y from one inverse; Arratia, Goldstein and Gordon for Chen–Stein. The six-leaf sharing
  of 3.1 is derived here.
