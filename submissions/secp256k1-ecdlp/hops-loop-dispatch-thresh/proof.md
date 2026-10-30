# Stored x-only leaf birthday search on secp256k1 (125.629556)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.390000001 over the challenge and fresh algorithm coins; on every run its total work is at most
`time_log2 = 125.629556 = 125629556/1000000` target-group-operation equivalents (C = 1000 scalar AArch64 instructions
per unit). Peak memory is 2^135.7542 bytes; preprocessing, 2^111.0349 units, is included; there is no advice.

**Reviewer's guide.** The bound rests on three things.
- **Success (section 6):** Chen–Stein on cross-node leaf pairs with b3 = 0 (6.2; its one non-routine step is the
  non-proportionality claim), usefulness (6.3), the defect table (6.4), whose batch row comes from the arithmetic-site
  table of 3.4, and the final inequality (6.5). `certificate.py` evaluates the terms with exact Fractions (e^-μ by an
  upper enclosure), checks each row and the total, and prints the bound 0.3900000040 (the run's `success_lower_bound`).
- **Work (section 7):** `certificate.py` recomputes the ledger exactly and checks that the claim is the least multiple
  of 10^-6 at or above log2 of the total; the run executes the kernels and checks their counts against the certificate
  (`counts_match_certificate`). The two data hops are priced at their itemized 26.5 and 38 (executed endpoint,
  implemented loop dispatch, link work and flag store), the local pass at its executed 9. The 40 ceiling is retained
  only for the analytic messages (Q broadcast, candidates, barriers, program delivery).
- **Machine (2.1; scope in section 9):** it charges every processor's instructions and nothing else, per the accounting
  convention classical-group-work-v2: "instruction counts, not CPU cycles" (a stalled store counts once); "Sum all
  processors' work, including coordination" (handlers, dispatch, barriers); "Memory has no scalar contribution" (windows
  and queues are memory in section 5; clearing the window flags is charged).
- **Not relied on:** any heuristic, any prior acceptance or score, or extrapolation from the toy or the 4 × 4 run.

## 1. Overview

Each node computes x(X ± Y), the **leaves**, for independent points X and Y, with one batched inverse per pair
(Galbraith, Wang and Zhang, ePrint 2015/605). Keyed by a six-class canonical x, each leaf goes in two hops to a home
node, where equal keys give the logarithm. This is an analytic bound, not a practical break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each with a private memory below 2^64 bytes and two scalar
AArch64 processors sharing it, a main processor and a handler processor. The cost model counts the instructions every
processor executes; storage and wiring carry no price.
- **Links.** Each ordered pair (s, t) of nodes in a row has a row link, and in a column a column link; a node's links to
  itself are its two loopback links, so a record for its own node takes the same path. A link ends in a 64-byte window
  in t's memory, a 32-byte record slot and a full flag. A node has 2^37 windows, one per node of its row and of its
  column, with flags cleared in preprocessing; windows to all 2^72 nodes would not fit its 64-bit space, hence two hops.
- **Send.** The column and row ports are fixed 40-byte device registers. The sender stores the destination index at +0,
  K2 and K3 at +8, and w2 and w3 at +24; the index selects the link and is not carried. The third store completes the
  message: the link carries the four words into the window and sets the flag, executing nothing.
- **Backpressure.** A store to a link whose window is full stalls until the flag clears; a stalled store is one executed
  instruction.
- **Receive.** Each arrival appends the window's index to the receiver's event queue, a ring of 2^37 entries (a window
  holds at most one pending arrival). The handler loads the entry, forms the window address, reads the record, runs the
  phase's receive body and clears the flag; it waits for arrivals with a stalling load (section 4).
- **Progress.** A handler never sends, so it never stalls on a link, and it has its own processor and registers: every
  window drains, and the main processor's registers are untouched.
- **Barriers.** A node joins once its main program has sent the phase's records and its ports' counters of uncleared
  windows read zero (one stalling load); a combining tree over rows, then columns, completes when all have joined, so no
  record is then in flight.

### 2.2 Leaves and records

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 108964788` and
`q2 = 2^27`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are x(X_k ± Y_l)
for all k and l: `n0 = 2*q1*q2` per node and `S = 2^72*n0 = 2^126.699287` in all. Leaf `X_k + s*Y_l` is `[m]P + [n]Q`
with `m = a_k + s*c_l` and `n = b_g`. A point Z ≠ O has the class `{±Z, ±φ(Z), ±φ²(Z)}`, where φ(x, y) = (βx, y) = [λ],
and the **key** `K = min(x, βx mod p, β²x mod p)` (limbs K0..K3) is shared by the class and differs between classes. A
record is `[K2, K3, w2, w3]` with `w2 = i | k<<36` and `w3 = d2 | l<<36 | s<<63`; its digits are `d1 = K0 mod 2^36`
(hop-1 row), `d2 = K1 mod 2^36` (hop-2 column, then replaced by the source column), `d3 = K2 mod 2^24` and
`d4 = K3 mod 2^48`.

### 2.3 Phases

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators, header constants,
   window flags, the d4 counters and the d3 cursor table. Then Q is broadcast in 2 hops.
1. **Hop 1:** the batch kernels (section 3) run for each X_k, and the emit sends each leaf through the column port to
   `(d1, j)`, whose `h1recv` appends it to a buffer of `cap_c` records or drops it when the buffer is full (the flag is
   still cleared).
2. **Hop 2:** `h2send` sends each buffered record through the row port to `(own row, d2)`, whose `h2recv` appends it to
   a d3 bucket of capacity `cap_3`, or drops it when the bucket is full, and counts d4.
3. **Grouping:** prefix sums of the d4 counts; one local pass writes each record's address into a pointer array at its
   d4 cursor, bucket by bucket in d3 order; a scan lists adjacent records with equal K2 and K3 as candidates, up to
   2^20, and a full list ends it.
4. **Resolution,** in four bulk-synchronous sub-phases: a home `(h_r, h_c)` sends two requests per candidate to the
   sources `(g_r, g_c)` via the relay `(g_r, h_c)`, and replies of at most 8 records (one leaf's scalars and points)
   return via `(h_r, g_c)`. A source answers at most 2^21 requests; a relay holds at most 2^22 records per sub-phase and
   drops the excess. The home recomputes `L = [m]P + [n]Q` and `L'`, finds the unit u with `L' = uL` if any, and, if
   `n' - u*n ≠ 0`, outputs `x = (u*m - m')/(n' - u*n) mod r` after checking `[x]P = Q`.

**Barriers:** 13 (program delivery 2, Q broadcast 2, the ends of phases 0 to 3 4, the phase-4 sub-phases 4, output 1),
within the allowance of 16. **Caps:** products (3.2), handler events (at most S per hop, each at most its main path),
stores, candidates and coin retries are capped or bounded in total, so the ledger holds on every run.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's. The forward pass stores `dx_l = x(Y_l) - x(X)` and `acc_(l-1)` in slot l and forms
`acc_l = acc_(l-1)*dx_l` (`acc_(-1) = 1`); fe_inv gives `R = 1/acc_(q2-1)`; the backward pass, for `l = q2-1..0`, sets
`inv = R*acc_(l-1)`, `R <- R*dx_l`, `λ± = inv*(y_Y ∓ y_X)` and `x± = λ±² - x_X - x_Y`. As `R = 1/acc_l` before step l,
`inv = 1/dx_l` and λ± is plus or minus the slope of `X ± Y_l`, so `x± = x(X ± Y_l)`. The arithmetic is F's row schedule,
with the row-shift and fold carries taken inline rather than by branches; its only errors are the sites of 3.4.

### 3.2 Key pass and counts

`keys` takes each leaf: KEY (3.3) gives the class index j or drops the leaf; for j = 1, 2, KMUL computes `K = x*β^j` and
counts one product; the emit writes the record to the column port. As j is uniform on {0, 1, 2} for a leaf other than O,
a node starts a batch only while its product count is at most `cap_j - 2q2`, `cap_j = ceil(2*n0/3) + 2^42`. Executed
counts per incidence (fixed part in parentheses): `fwd` 105 (3), `bwd` 523 (5), `keys` 66 for two leaves, including
their 13 emit instructions (11), plus 83 per product; fe_inv (`secp_fe.s`) is 23646. Kernels may overwrite x0-x29 and
the flags and preserve sp; `bwd`, `keys` and fe_inv save and restore x30, and `fwd` and the transport kernels do not
write it. The driver keeps no registers across calls and loads the product count into x18 only for `keys` (fe_inv uses
x18).

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
(FDX) computes `t = (Z mod 2^256) + floor(Z/2^256)*C0` and keeps t mod 2^256. A fused square (SQP, SQM) forms
`λ² = L + 2^256*H` and `s' = L + H*C0 - x_X - x̃_Y + 2^257` exactly, with `x̃_Y = x_Y + 2C0` (below 2^256 unless
`x_Y >= p - C0`, a whole-run row of 6.4) and top limb at most `C0 + 2`, and folds `s'` the same way. Only these three
steps drop a carry: SHX's last ADCS, FDX's ADC and the fused square's final ADC. Without a dropped carry each returns a
value below 2^256 congruent to the exact one:
- **row shift** of `B = (b0, b1, b2, b3)`: the carry is lost iff `(B mod 2^192) + floor(b3*C0/2^64) >= 2^192`, which
  needs `B mod 2^192 > 2^192 - 2^33`: at most 2^97 rows. Otherwise `SHX(B) = B*2^64 - b3*p`.
- **fold** of `Z < 2^322`: lost iff `t >= 2^256`, which needs `Z mod p = t - p` in `[C0, 2^66*C0)`.
- **fused square**: lost iff its fold sum is at least 2^256, which needs the true `x±` in `[C0, C0² + 3C0)`.

A fold or fused square returns a value `>= p` only when the true value is below C0. That matters only where a canonical
value is needed: `x±` entering KEY, and the emitted key K. Every other operation is exact on its inputs: RF and RA;
KMUL's sum over the constant rows of β^j (below `4*2^64*p < 2^322`; only its FDX can drop a carry); FS and FADDW on
canonical inputs; KEY's integer arithmetic (3.3 covers its index); the emit; and F's `fe_mul_f` and `fe_sqr_f`, which
fold the 512-bit product, fold its fifth limb again and add C0 exactly when that last sum carries past 2^256, so fe_inv
(255 squarings, 15 multiplications) is exact for every input below 2^256. Generator construction, the resolver and the
check `[x]P = Q` run outside the batch kernels on `fe_mul_f`, `fe_sqr_f` and fe_inv, with exact arithmetic mod r for the
solve. Each coordinate is brought below p (an ADDS chain of C0 detects `t >= p`, and CSELs keep `t - p`) before it is
subtracted, compared or stored. So these steps are exact and have no sites, and the X and Y tables hold exact canonical
coordinates.

**Distribution.** Given X ≠ O and every other generator, `Y_l` is uniform on the r - 1 points of G \ {O}. So `x(Y_l)`,
`dx_l = x(Y_l) - x(X)`, `acc_l = acc_(l-1)*dx_l`, `R_l = 1/acc_l` and `inv_l = 1/dx_l` take each value with probability
at most `2/(r-1)` (two points per x; the rest are bijections while no `Y_m = ±X` for m ≤ l, the degenerate case below).
So do the slopes λ±, as a line through X meets E in at most two more points, and `x(X ± Y_l)`, as `X ± Y_l` is uniform
on r - 1 points. A row `B_k ≡ v*2^(64k)` determines its operand v, so at most 2^97 values of v make a given row shift
carry.

| site, at incidence l | runs per batch | fresh generator: operand | bad operand values |
| --- | --- | --- | --- |
| fwd: row shifts of `acc_(l-1)`, l ≥ 1 | 3(q2 - 1) | `Y_(l-1)`: `acc_(l-1)` | 2^97 |
| fwd: fold of `acc_l`, l ≥ 1 | q2 - 1 | `Y_l`: `acc_l` | 2^66*C0 |
| bwd: row shifts of `R_l` | 3q2 | `Y_l`: `1/acc_l` | 2^97 |
| bwd: fold of `inv_l = R_l*acc_(l-1)`, l ≥ 1 | q2 - 1 | `Y_l`: `1/dx_l` | 2^66*C0 |
| bwd: fold of `R_(l-1) = R_l*dx_l`, l ≥ 1 | q2 - 1 | `Y_(l-1)`: `1/acc_(l-1)` | 2^66*C0 |
| bwd: row shifts of `inv_l` | 3q2 | `Y_l`: `1/dx_l` | 2^97 |
| bwd: folds of `λ± = inv_l*(y_Y ∓ y_X)` | 2q2 | `Y_l`: slopes | 2^66*C0 |
| bwd: fused squares | 2q2 | `Y_l`: `x(X ± Y_l)` | C0² + 3C0 |
| keys: KMUL fold of `K = β^j*x±`, j = 1, 2 | 2q2 | `Y_l`: `x(X ± Y_l)` | 2*2^66*C0 |

A run fails with probability at most `2/(r-1)` times its bad count. At l = 0 three row shifts and three folds cannot
fail, and are not in the table:
- fwd's rows of `acc_(-1) = 1` are 1, 2^64 and 2^128, and its fold of `dx_0*1` has `Z < 2^256`;
- bwd's fold of `inv_0 = R_0*1` has `Z = R_0 < 2^256`;
- bwd's last update `R_0*dx_0 ≡ 1` folds to 1 or p + 1, as nothing else congruent to 1 lies below `2^256 + 2^66*C0`. It
  is stored at H + 64 and reloaded once at the end of the last pass, then unused.

The degenerate cases are X = O (probability 1/r) and `dx_l = 0`, that is `Y_l = ±X` (at most `2/(r-1)` per incidence),
which zeroes acc and spoils the batch. The first site of a batch whose predicate holds has exact operands, so a batch
has a defect with probability at most `d_batch`: the table's sum plus `2q2/(r-1) + 1/r` (`SITES` in `certificate.py`;
6.4).

## 4. Transport and grouping (`kernels/route.s`)

| step | executed per record | priced |
| --- | ---: | ---: |
| hop 1: emits (13 per two leaves) and `h1recv` (11; 9 if dropped) | 17.5 | 26.5 |
| hop 2: `h2send` (9) and `h2recv` (20; 12 if dropped) | 29 | 38 |
| local pass `pass4`, plus 7 per d3 bucket | 9 | 9 |
| scan: `7n - 3` for n records, plus 4 per listed candidate | 7 | 7 |

The handler bodies include the event path (LDR of the event, ADD of the window address, STR clearing the flag);
`wclear`, `init3`, `zero4` and `prefix4` are priced as executed. **Dispatch:** a handler runs its kernel as an endless
loop over its event queue, whose posted count is a device register that it reads with a stalling load. This adds at most
4 per record (load, compare, branch, taken count and a BIC wrapping the ring pointer, with B in place of the counted
SUBS and B.NE), or 7 with one kernel call per event (wait, call, CBZ, RET, wrap). Charging the hardware's work in
software, the event append to the sender (atomic fetch-add and store) and the uncleared-window counters to both ends
(atomic increment and decrement), adds 4 more. The link's store that sets the window's full flag adds 1 more. The data
hops are therefore priced at the implemented loop-dispatch total plus this flag store: `17.5 + 4 + 4 + 1 = 26.5` and
`29 + 4 + 4 + 1 = 38` per record. The same itemization under per-call dispatch would be 41 on hop 2, above the old 40
ceiling, so this package uses the endless-loop dispatch that the machine of 2.1 actually runs. The 40 ceiling is kept
for the analytic messages only, so no planned transfer is priced below its itemized cost.

**Local pass.** The grouping pass runs entirely on the home node's main processor. `pass4` reads each record's stored
d4 digit, forms the d4 cursor address, loads the cursor, writes the record's address into the pointer array and stores
the advanced cursor: 9 executed instructions per record, plus 7 per d3 bucket for the bucket bounds and pointer
advance. That is the whole pass; the record transfer is separately priced in the hop rows. The ledger therefore prices
`pass4` at its executed 9 rather than at the earlier conservative multiple of 20, and `certificate.py` checks the
constant against the kernel's executed count.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 key bits. A class forms one run unless an **intruder**, a record outside it with
the same 144 bits, sits inside (6.4); a run's adjacent pairs are all listed unless the candidate list is full. Writing
the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`, an equivalence relation, so a
run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

All regions lie in the node's single 64-bit address space. `certificate.py` sums them, checks that a node is below 2^64
bytes, and takes the peak as 2^72 node totals, every region live throughout.

| region | bytes |
| --- | --- |
| code, constants, stack, port and flag state; header | 2^20; 512 |
| X table (x, y); Y table (x, y, x + 2C0); scalars | 64*q1; 96*q2; 32*(q1 + q2 + 1) |
| batch slots (starting at a multiple of 2^34) | 64*q2 |
| hop-1 buffer, `cap_c = 6*n0 + 2^50` | 32*cap_c |
| d3 buckets, `cap_3 = ceil(m_3) + 2^32` (6.4); pointer array | 2^24*cap_3*32; 2^24*cap_3*8 |
| d3 cursor table; d4 counters | 2^24*16; 2^48*8 |
| candidates; requests; replies (up to 8 records each) | 2^20*16; 2^21*64; 2^24*64 |
| relay buffer; windows; event queue | 2^22*64; 2^37*64; 2^37*8 |
| row and column ports (device registers) | 2*64 |

The event-queue pointer wraps at the ring's end, and every other address the kernels form, the fixed ports included,
lies in these regions; a remote item is reached only by a message.

## 6. Success (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494342851` (the run's `mu`).
Same-node pairs are left out; their collisions only add to success.

### 6.2 Chen–Stein on cross-node pairs

Let `W = Σ_(α∈Γ) I_α`, and let `B_α` be the pairs that share a generator point with α. I_α depends only on α's four
generators, distinct because they come from two nodes, so I_α is independent of `{I_β : β ∉ B_α}` (b3 = 0). A generator
lies in 2q2 (an X) or 2q1 (a Y) leaves, so `|B_α| <= 4(q1 + q2)*S`.

For a unit u, the event `L' = u*L` with `L = X + s*Y` and `L' = X' + s'*Y'` is the linear relation
`X' + s'*Y' - u*X - u*s*Y = O` in the generators' logarithms; I_α is at most the sum of these 6 events. Take α ≠ β in Γ
sharing a generator, and units u and v. **The two relations are not proportional:** proportional relations have the same
support, so each node contributes exactly one X and one Y to both; the leaf of β in that node must use both, so it is
α's incidence there; the relations then differ only in signs, and comparing coefficients forces β = α. So some 2 × 2
minor is invertible: given all other generators, the relations fix two generators, each uniform on at least r - 1
points, and `E[I_α I_β] <= 36/(r-1)²`. The Chen–Stein method (Arratia, Goldstein and Gordon, Ann. Probab. 17 (1989),
Theorem 1) gives `|P(W = 0) - e^-μ| <= b1 + b2 <= |Γ|*4(q1 + q2)S*(P1² + 36/(r-1)²)`, which the certificate checks is
below 2^-96. In one pool, by contrast, `{X+Y, X'+Y'}` and `{X-Y', X'-Y}` give the same relation for u = 1 (the twins of
section 8).

### 6.3 Usefulness

A cross-node colliding pair is **detected** if both records are computed exactly and emitted, both survive the caps, no
intruder enters their run, and the candidate, request and relay caps are not reached first. Its run is then resolved, a
useful adjacent pair (section 4) gives x, and `[x]P = Q` is checked before output. The pair is useful iff `n' ≠ u*n`,
where n and n' are the b of its two nodes; given all points, this fails with probability exactly 1/r. So
`P(fail) <= P(W = 0) + E #{α : I_α = 1, α useless or not detected}`.

### 6.4 Defect rows

A **node-local defect** of L (a batch defect of 3.4, a window, the wrong index, the meter) is independent of L', which
is uniform given it, so `P(I_α ∧ defect) <= (6/r)*P(defect)`; over Γ a row costs at most `μ(1 + 1/r)` times its per-leaf
rate. **Caps and intruders** are bounded conditional on α's four generators. Every other leaf except the ± partners of
α's two incidences is then uniform on at least r - 1 points, so at most `6*2^(256-f)/(r-1)` of its mass falls on given
values of f key bits ("one key per class, six points per class"; hence cap_c is 6n0, not n0). A hop-1 buffer and a d3
bucket thus have means at most `m_c = 2^36*n0*6*2^220/(r-1) + 1` and `m_3 = S*6*2^160/(r-1) + 1`. The partners and L'
are fixed, and each cap leaves room for them: 3 fixed competitors at the hop-1 and d3 caps, 4 records in the candidate
list, 8 requests at a source, and 16 records (α's two replies) at a relay. Leaves of a defective batch elsewhere are not
uniform; the whole-run row "a computational defect in any batch" excludes them. Two other leaves' indicators have
covariance 0 without a common generator, at most 0 for one incidence's ± pair, 0 when they share a free Y, and at most
`1/(4(r-1)²)` when they share a free X; each of the at most `2(q1 + q2)` ± pairs through a generator of α adds at most
one leaf's mean. So each count's variance is at most its mean plus 1, and Chebyshev applies.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, dx = 0, X = O | `d_batch` | < 2^-125 |
| window drop; wrong class index (3.3) | `d_win`; `d_wj` | < 2^-54.29; in the total |
| product meter, Chebyshev with slack `2^42 - 2q2 - 1` | `d_j` | < 2^-31 |
| hop-1 buffer; d3 bucket (Chebyshev) | `d_c`; `d_3` | < 2^-42; < 2^-30 |
| intruder (per pair), `S*6*2^112/(r-1)` | `lam_int` | < 2^-14.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-33.4; < 2^-35.9; < 2^-31.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26 |
| two coin range rejections; `x_Y >= p - C0` (whole run) | in `glob` | all of `glob` < 2^-26 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.390000001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 108964788, where the left side is at least 0.3900000040, and
the target is 10^-5 above the rules' 0.39.

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, the candidate and the barrier
against itemized work. **Generator construction, 2^21 per point:** a complete Jacobian double-and-add, exceptional-case
tests, two inversions with the affine conversion and the table write, all on the exact arithmetic of 3.4: `fe_mul_f` and
`fe_sqr_f` at 115 and 82 instructions as executed, plus 4 and 3 for the call as in fe_inv, and an addition or
subtraction with its operands brought below p at 27. That is 665,910 without the tests and loop control. **Node set-up,
2^20:** a program of at most 2^16 bytes (the 1,793 kernel instructions, generator construction, resolver and phase
control), or 2^11 records over 2 hops at 40, plus port, buffer and pointer state. **Candidate, 2^25:** two requests and
two replies of at most 8 records over 2 hops at 40, two exact additions with inversion, the six-unit test, the solve mod
r and `[x]P = Q`, on the same exact arithmetic; 769,058. **Barrier, 2^82:** up to 4 messages per node at 40, plus each
row and column leader's 2^36 handler steps at 11. **Header 148, batch set-up 64, inversion call 16, request served 2^10,
output 2^20.** **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ
retrying once on the failure flag, plus a range test. In this machine RNDR is the coin source: it returns 64 fresh
uniform bits, the algorithm coins of the rules' probability space, and never reports failure. The B.EQ is charged for
hardware that can report one; in the model it is never taken, so coins add no failure row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, the barriers and the output.

| row (P: preprocessing) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1 + q2 + 1)*64 (P) | 105.8575 |
| Y generators: q2*2^21 (P) | 120.0000 |
| R_g and X generators: (q1 + 1)*2^21 | 119.6993 |
| header: 148 (P) | 79.2095 |
| window flags: 7*2^35 + 1 (P) | 109.8074 |
| Q broadcast: 2*(2^36 + 2^72) records at 40 | 78.3219 |
| d4 counters: 2^48 + 1 (P) | 120.0000 |
| d3 cursor table: 5*2^24 + 1 (P) | 98.3219 |
| forward: q1*(105*q2 + 3) | 132.4136 |
| inversion: q1*(23646 + 16) | 113.2296 |
| backward: q1*(523*q2 + 5) | 134.7300 |
| key pass without emits: q1*(53*q2 + 11) | 131.4272 |
| class-index products: cap_j*83 | 132.4897 |
| batch set-up: q1*64 | 104.6993 |
| hop 1: S*53/2 | 131.4272 |
| hop 2: S*38 | 131.9472 |
| d4 prefix sums: 14*2^46 + 1 | 121.8074 |
| local pass: S*9, plus 7*2^24 + 1 per node | 129.8692 |
| scan: S*7, plus 4 + 4*(2^20 + 1) per node | 129.5067 |
| candidates: 2^20*2^25 | 117.0000 |
| requests served: 2^21*2^10 | 103.0000 |
| barriers: 16*2^82 | 86.0000 |
| output | 20.0000 |
| **total R** | **135.595340** |

`R = 65805873179022909518403237047120574808064` instructions, 476.4079 per leaf (the run's `per_leaf`), and
`W = R/1000 = 2^125.62955586`. **The claim, 125.629556, is the least multiple of 10^-6 with 2^claim >= W**; the
certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds on ln 2. Rows marked P are
preprocessing.

## 8. Evidence

`experiment.py` has the declared input
`{"seed": 1, "trials": 6000, "control_trials": 1000, "curve_trials": 100, "seconds": 15}` (`seconds` 5 to 15); it finds
the package files beside itself, under files/, or under /snapshot/files, as the hosted runner lays them out.
Budget-limited checks run a guaranteed minimum (2,000 random and 600 near-tie class-rule inputs, 400 multi-node trials,
100 per control, 10 curve checks), then continue on their share of the remaining budget; requested counts are caps, and
larger counts are quoted only from the census. Part A runs the leaf kernels on production batches, and fe_inv, against
exact arithmetic, checks each port message (via a device hook) and a sentinel x30, and checks the KEY model against the
exact class rule. It also checks the sites of 3.4 (`site_checks`). For the row shift, the fold, both fused squares (SQP
and SQM) and KMUL, inputs built to meet the predicate deviate exactly when it holds, and boundary and random inputs are
exact; each also equals the bit-level model of its arithmetic. FS, FADDW, `fe_mul_f` and `fe_sqr_f` are exact on edge,
random and small-residue inputs, including square roots of small residues; the run reports how many take the routines'
final +C0. A batch with `Y_l = ±X` gets acc = 0 and loses all its other leaves. The last R update is ≡ 1, and
overwriting it changes no message. The SHX, FDX, SQP and SQM invocations per loop pass of fwd, bwd and keys match the
table (`sites_match_certificate`). Part B runs `route.s` on a 4 × 4 machine with windows, event queues and flags
(`wclear`, `init3` and `zero4` before hop 1; `prefix4`, `pass4` and `scan` after hop 2) against a Python model of the
same stable order, including the full-buffer, full-bucket and full-list paths; groups minus candidates is nonzero
because the toy's short digits let intruders split groups. Part C runs the q25 toy in log form (the class of `[t]G` is
`t^6 mod r`): 74 nodes of 4 × 4 grids, one-pool 28 × 42 controls with b per X or a shared b, and a curve pipeline that
must match the log form. Part D runs `certificate.py`, and `counts_match_certificate` compares its count constants with
parts A and B.

**Census.** `census-toy.json` is an off-runner run of the same code (`seconds` 1500, above the hosted maximum of 15;
each toy count reaches its request). Over 100000 multi-node trials, success is 0.39182 with all pairs and 0.38849 with
cross-node pairs only, against Poisson 0.39008 and a toy Chen–Stein bound of 0.0134. The one-pool controls give 0.3423
and 0.314 against their twin and shared-b models, 0.33759 and 0.30974, and the curve pipeline matches in 1000 of 1000
trials. The toy neither establishes nor contradicts production success; read a declared run's figure with its own
standard error. The toy and the 4 × 4 run are evidence at their own scope only.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1. Only the two endpoints of a hop execute instructions for it, and
  every one of them is charged. On a machine whose links pass through intermediate processors, their instructions would
  be added; that machine is not claimed.
- **Not executed:** the allowances of section 7, the dispatch loop (section 4; the hop prices cover it), and the
  production attack, machine and links.
- **Provenance**, not a premise: package F (submission 8f733826) for the row products (RF, RA, and the pair schedule of
  W2S and W2L), FS's CSEL form, fe_inv, the F1 row and its error analysis, and `kernels/a64.py`, whose extensions its
  docstring lists; Galbraith, Wang and Zhang for X + Y and X - Y from one inverse; Arratia, Goldstein and Gordon for
  Chen–Stein.
