# Stored x-only leaf birthday search on secp256k1 (125.567689)

**Claim.** For a uniform full-production secp256k1 input `Q = [x]P`, the classical cold algorithm below recovers x with
probability at least 0.3900001 over the challenge and fresh algorithm coins; on every run its total work is at most
`time_log2 = 125.567689 = 125567689/1000000` target-group-operation equivalents (C = 1000 scalar AArch64 instructions
per unit). Peak memory is 2^135.9988 bytes; preprocessing, 2^111.0353 units, is included; there is no advice.

**Reviewer's guide.** The bound rests on four things.
- **No uncharged bookkeeping (2.1).** Every per-record step is an instruction of the node's own processor, and the
  ledger charges each at its executed count: hop 1 is the emit inside `keys` (12.5 per record), hop 2 is `h2s` (14), and
  grouping is `bkt`, `pass4` and `scan` (30). Each port's in-flight count, which DSB waits on, is charged as 2
  instruction-equivalents per message, so the hops are charged 14.5 and 16. The only uncharged action is a link's
  delivery: it writes a message's four data words at the address its sender stored and reports the write back to the
  port, with no flags, rings or counters. Every address, cursor, count and cap test is a counted instruction of the
  sender or of the node that later reads the data. Barriers, phase control and the phase-4 count exchanges are software
  with itemized allowances (section 7).
- **Success (section 6):** Chen–Stein on cross-node leaf pairs with b3 = 0 (6.2; its one non-routine step is the
  non-proportionality claim), usefulness (6.3), the defect table (6.4), whose batch row comes from the arithmetic-site
  table of 3.4, and the final inequality (6.5). `certificate.py` evaluates the terms with exact Fractions (e^-μ by an
  upper enclosure), checks each row and the total, and prints the bound 0.3900001047 (the run's `success_lower_bound`).
- **Work (section 7):** `certificate.py` recomputes the ledger exactly and checks that the claim is the least multiple
  of 10^-6 at or above log2 of the total. The run executes the leaf kernels and fe_inv on production points, and the
  emit and each kernel of `route.s` on a 4 × 4 machine, and checks each call's count against the certificate's constants
  (`counts_match_certificate`).
- **Machine (2.1; scope in section 9):** per the accounting convention classical-group-work-v2, "loads/stores,
  addressing and control instructions count" (all are the nodes'), "instruction counts, not CPU cycles" (a stalled store
  or WFE counts once), "Sum all processors' work, including coordination" (barriers and count exchanges are charged) and
  "Memory has no scalar contribution" (slots and tables, section 5).
- **Not relied on:** any heuristic, any prior acceptance or score, or extrapolation from the toy or the 4 × 4 run.

## 1. Overview

Each node computes x(X ± Y), the **leaves**, for independent points X and Y, with one batched inverse per pair
(Galbraith, Wang and Zhang, ePrint 2015/605). Keyed by a six-class canonical x, each leaf goes in two bulk-synchronous
hops to a home node, where equal keys give the logarithm. This is an analytic bound, not a practical break.

## 2. The algorithm

### 2.1 Machine

The machine has 2^72 nodes `(i, j)`, `0 <= i, j < 2^36`, each with one scalar AArch64 processor and a private memory
below 2^64 bytes. The cost model counts the instructions the processors execute; storage and wiring carry no price.
- **Ports and links.** A node has a column port and a row port, 48-byte blocks of device registers in its address space.
  To send, it stores a destination index (a row of its column or a column of its row, its own included) at +0, a
  destination address at +8 and four data words at +16..+47; the STP to +32 completes the message. The link writes the
  four words at that address of the destination's memory, reports the write back to its port, and does nothing else: no
  flags, rings or counters, and no instructions; its only state is the message it carries. A store that finds the link
  busy stalls until the link takes it; a stalled instruction counts once.
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

Node g draws fresh uniform `b_g` and `a_1..a_q1` in `[0, r)` and `c_1..c_q2` in `[1, r)`, with `q1 = 108965147` and
`q2 = 2^27`, and stores them with `R_g = [b_g]Q`, `X_k = R_g + [a_k]P` and `Y_l = [c_l]P`. Its leaves are x(X_k ± Y_l)
for all k and l: `n0 = 2*q1*q2` per node and `S = 2^72*n0 = 2^126.699292` in all. Leaf `X_k + s*Y_l` is `[m]P + [n]Q`
with `m = a_k + s*c_l` and `n = b_g`. A point Z ≠ O has the class `{±Z, ±φ(Z), ±φ²(Z)}`, where φ(x, y) = (βx, y) = [λ],
and the **key** `K = min(x, βx mod p, β²x mod p)` (limbs K0..K3) is shared by the class and differs between classes. A
record is `[K2, K3, w2, w3]` with `w2 = i | k<<36` and `w3 = d2 | l<<36 | s<<63`; its digits are `d1 = K0 mod 2^36`
(hop-1 row), `d2 = K1 mod 2^36` (hop-2 column), `d3 = K2 mod 2^24` and `d4 = K3 mod 2^48`. At its home a record's source
is `(w2 mod 2^36, j)`, where j is the index of the slot it arrived in (2.3).

### 2.3 Phases

A pair slot holds `cap_p = 2^22 - 2^12` records of 32 bytes. A node in row i writes at hop 1 only into slot i of R1, at
`R1 + 32*cap_p*i`, at each node of its column, and a node in column j writes at hop 2 only into slot j of R2 at each
node of its row. Each port's cursor table holds, per destination, the sender's cursor in its slot there and the slot's
end; `tinit` sets each of its 2^36 entries to the slot's bounds.

0. **Preprocessing:** program delivery (2 hops from node (0, 0)), node state, coins, the Y generators, header constants,
   the two cursor tables, the d4 counters (`zero4`) and the d3 cursor table (`init3`). Then Q is broadcast in 2 hops.
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
   The home recomputes `L = [m]P + [n]Q` and `L'`, finds the unit u with `L' = uL` if any, and, if `n' - u*n ≠ 0`,
   outputs `x = (u*m - m')/(n' - u*n) mod r` after checking `[x]P = Q`.

**Barriers:** 21 (program delivery 2, Q broadcast 2, the ends of phases 0 to 3 4, phase 4 12, output 1), within the
allowance of 32. **Caps:** products (3.2), pair slots, d3 buckets, candidates, requests, relay records and coin retries
are capped or bounded in total, so the ledger holds on every run.

## 3. Leaf kernels (`kernels/leaf.s`)

### 3.1 One batch

A batch is one X and all q2 Y's. The forward pass stores `dx_l = x(Y_l) - x(X)` and `acc_(l-1)` in slot l and forms
`acc_l = acc_(l-1)*dx_l` (`acc_(-1) = 1`); fe_inv gives `R = 1/acc_(q2-1)`; the backward pass, for `l = q2-1..0`, sets
`inv = R*acc_(l-1)`, `R <- R*dx_l`, `λ± = inv*(y_Y ∓ y_X)` and `x± = λ±² - x_X - x_Y`. As `R = 1/acc_l` before step l,
`inv = 1/dx_l` and λ± is plus or minus the slope of `X ± Y_l`, so `x± = x(X ± Y_l)`. The arithmetic is F's row schedule,
with the row-shift and fold carries taken inline rather than by branches; its only errors are the sites of 3.4.

### 3.2 Key pass, emit and counts

`keys` takes each leaf: KEY (3.3) gives the class index j or drops the leaf; for j = 1, 2, KMUL computes `K = x*β^j` and
counts one product; the **emit** forms d1, loads the cursor and slot end of d1 from the table at port + 64, and drops
the leaf if the slot is full (5 instructions); otherwise it forms w3, writes the message with three STPs and stores the
advanced cursor (7, or 8 with the sign bit). As j is uniform on {0, 1, 2} for a leaf other than O, a node starts a batch
only while its product count is at most `cap_j - 2q2`, `cap_j = ceil(2*n0/3) + 2^42`. Executed counts per incidence
(fixed part in parentheses): `fwd` 105 (3), `bwd` 523 (5), `keys` 78 for two leaves, including their emits of 12 and 13
(11), plus 83 per product; fe_inv (`secp_fe.s`) is 23646. Kernels may overwrite x0-x29 and the flags and preserve sp;
`bwd`, `keys` and fe_inv save and restore x30, and `fwd` and the transport kernels do not write it. The driver keeps no
registers across calls and loads the product count into x18 only for `keys` (fe_inv uses x18).

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

## 4. Transport and grouping (`kernels/route.s` and the emit)

| step, per record | executed | charged |
| --- | ---: | ---: |
| hop 1: the emit in `keys` (12 for the + leaf, 13 for the - leaf; 5 if the slot is full) | 12.5 | 14.5 |
| hop 2: `h2s` (9 if the slot is full) | 14 | 16 |
| bucket pass `bkt` (9 if the bucket is full) | 16 | 16 |
| local pass `pass4` | 7 | 7 |
| `scan`: `7n - 3` for n records, plus 4 per listed candidate | 7 | 7 |
| **total** | **56.5** | **60.5** |

Each executed count is its loop's instruction list in `route.s` or `leaf.s`; the charged count adds the in-flight
count of 2.1, 2 per message, at both hops. A hop's receiving node executes nothing on arrival;
it reads the records in the next phase, with the LDPs of `h2s` (hop 1) or `bkt` (hop 2). Per-node costs are rows of
section 7, and every address a kernel forms lies in the regions of section 5.

**One local pass suffices.** A home's records share d1 and d2; the d3 buckets and the stable d4 pass make records with
equal (d3, d4) adjacent, grouping 144 key bits. A class forms one run unless an **intruder**, a record outside it with
the same 144 bits, sits inside (6.4); a run's adjacent pairs are all listed unless the candidate list is full. Writing
the members `u_i*Z` with Q coefficients `n_i`, a pair is useless iff `n_i/u_i = n_j/u_j`, an equivalence relation, so a
run with a useful pair has a useful adjacent pair.

## 5. Node layout and memory

Every node uses the same addresses. `certificate.py` sums the regions and checks each bound below.

| range | regions | live |
| --- | --- | --- |
| below 2^52 (2^51.0085 used) | code, stack, the five barrier words (one reservation granule each), header, X and Y tables, scalars, batch slots (at a multiple of 2^34), the two ports with their cursor tables (2^40 bytes each), E1 and E2 (2^41 each), the d3 cursor table, the d4 counters (2^51) and the phase-4 buffers and tables (`FIXED` in `certificate.py`) | all phases |
| R1 at 2^52, `32*cap_p*2^36 = 2^62.9986` bytes | hop-1 pair slots | phases 1-2 |
| the same range, from phase 3 | d3 buckets `16*2^24*cap_3` and pointer array `8*2^24*cap_3`, together 2^62.3655 | phase 3 |
| R2 after R1, ending at 2^63.998943 | hop-2 pair slots (phase 4 reads candidates' records at their addresses) | phases 2-4 |

So every node stays below 2^64 bytes, and the peak, in phase 2, is `2^72*(fixed + |R1| + |R2|) = 2^135.9988` bytes. A
remote item is reached only by a message.

## 6. Success (no heuristic)

### 6.1 Model

Coins are uniform and independent, and Q is uniform. Each `Y_l` is uniform on G \ {O}; given R_g, the
`X_k = R_g + [a_k]P` are independent and uniform. So the X's are independent uniform points, independent of the Y's, and
nodes are independent; a leaf is uniform on G, and leaves of different nodes are independent. Given all generator
points, the `b_g` are still independent and uniform: for each value of b_g exactly one tuple `(a_k)` fits the points.

Γ is the set of unordered pairs of leaves from different nodes, `|Γ| = (S² - 2^72*n0²)/2`. For α = {L, L'} ∈ Γ let
`I_α = 1[L' ∈ class(L)]`; then `P(I_α = 1) = P1 = (6(r-1) + 1)/r²` and `μ = |Γ|*P1 = 0.494329713` (the run's `mu`).
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
values of f key bits ("one key per class, six points per class"; hence the factor 6 in the means). A pair slot has mean
at most `m_p = n0*6*2^220/(r-1) + 1`: at hop 1 it holds one node's leaves with a given d1 (36 bits), at hop 2 one
column's 2^36*n0 leaves with a given (d1, d2) (72 bits), the same bound. A d3 bucket has mean at most
`m_3 = S*6*2^160/(r-1) + 1`. The partners and L' are fixed, and each cap leaves room for them: 3 fixed competitors at a
pair slot and at a d3 bucket, 4 records in the candidate list, 8 requests at a source, and 16 records (α's two replies)
at a relay. Leaves of a defective batch elsewhere are not uniform; the whole-run row "a computational defect in any
batch" excludes them. Two other leaves' indicators have covariance 0 without a common generator, at most 0 for one
incidence's ± pair, 0 when they share a free Y, and at most `1/(4(r-1)²)` when they share a free X; each of the at most
`2(q1 + q2)` ± pairs through a generator of α adds at most one leaf's mean. So each count's variance is at most its mean
plus 1, and Chebyshev applies: a pair slot or d3 bucket with cap c and mean bound m drops a given record with
probability at most `(m + 1)/(c - 2 - m)²`.

| row (per leaf of the pair unless stated) | term in `certificate.py` | checked |
| --- | --- | --- |
| useless collision (per pair) | `1/r` | - |
| batch defect: the sites of 3.4, dx = 0, X = O | `d_batch` | < 2^-125 |
| window drop; wrong class index (3.3) | `d_win`; `d_wj` | < 2^-54.29; in the total |
| product meter, Chebyshev with slack `2^42 - 2q2 - 1` | `d_j` | < 2^-31 |
| pair slot at hop 1, and again at hop 2 (Chebyshev, `cap_p`) | `d_p`, counted twice | < 2^-19.9 each |
| d3 bucket (Chebyshev, `cap_3 = ceil(m_3) + 2^32`) | `d_3` | < 2^-30 |
| intruder (per pair), `S*6*2^112/(r-1)` | `lam_int` | < 2^-14.6 |
| candidate, request and relay caps (Markov; a reply counts 8 records) | `d_cl`; `d_rq`; `d_rl` | < 2^-33.4; < 2^-35.9; < 2^-31.2 |
| a computational defect in any batch (whole run) | `g_def = 2^72*q1*d_batch` | < 2^-26 |
| two coin range rejections; `x_Y >= p - C0` (whole run) | in `glob` | all of `glob` < 2^-26 |
| RNDR failure (whole run) | `g_rndr` = 0 under the rules' ideal coins (section 7) | - |

### 6.5 The bound

`certificate.py` picks the least q1 for which

```text
1 - e^-mu - (b1 + b2) - mu(1 + 1/r)(1/r + intruder + 2 * sum of per-leaf rows) - whole-run rows >= 0.3900001,
```

with an upper enclosure of e^-μ and exact Fractions: q1 = 108965147, where the left side is at least 0.3900001047, and
the target is 10^-7 above the rules' 0.39.

## 7. Ledger and claim

Rows are scalar instruction counts: the executed counts of sections 3 and 4, and allowances whose content
`certificate.py` lists in its comments; it checks generator construction, node set-up, a broadcast message, the
candidate, the count exchange and the barrier against itemized work.
- **Generator construction, 2^21 per point:** a complete Jacobian double-and-add, exceptional-case tests, two inversions
  with the affine conversion and the table write, all on the exact arithmetic of 3.4: `fe_mul_f` and `fe_sqr_f` at 115
  and 82 instructions as executed, plus 4 and 3 for the call as in fe_inv, and an addition or subtraction with its
  operands brought below p at 27; 665,910 without the tests and loop control.
- **Node set-up, 2^20:** a program of at most 2^16 bytes (the 1,821 kernel instructions, generator construction,
  resolver and phase control), or 2^11 messages over 2 hops at 40, plus node state.
- **A message outside the hops, 40 per hop** (program, Q broadcast, phase 4): a broadcast message executes 8 (two LDPs,
  three STPs, ADD, SUBS, B.NE), plus 2 for its in-flight count.
- **Candidate, 2^25:** two requests and two replies of at most 8 records over 2 hops at 40, two exact additions with
  inversion, the six-unit test, the solve mod r and `[x]P = Q`, on the same exact arithmetic; 769,058.
- **Phase-4 count exchange, 2^41 per node and sub-phase:** at most 12 per entry of the node's 2^36-entry table of
  received counts (load, epoch test, select, prefix, zero test, cap test, loop) and at most 64 per message (its count,
  the count message, the offset reply and the write, at most 16 each with in-flight counts) for at most 2^24 messages:
  below 2^39.6.
- **Barrier, 2^10 per node:** at most 5 waits of LDAXR, CMP, B.EQ and at most one WFE, B and re-check (8 each), at most
  5 messages of 8 (6, plus 2 for the in-flight count), DSB and 20 of set-up: 101.
- **Phase control, 2^10 per node:** the kernel calls between barriers (`tinit` twice, `init3`, `zero4`, `fin` twice,
  `h2s`, `bkt`, `prefix4`, `pass4` and `scan`, with register set-up and BL, at most 8 each), the E-entry addresses, and
  n for `scan` from the last d4 cursor: below 2^7.
- **Header 148, batch set-up 64, inversion call 16, request served 2^10, output 2^20.**
- **Coins, 64 per scalar:** at most two draws of 4 words `MRS Xt, RNDR` (FEAT_RNG), each with a B.EQ retrying once on
  the failure flag, plus a range test. RNDR is this machine's coin source: it returns 64 fresh uniform bits, the
  algorithm coins of the rules' probability space, and never reports failure, so the charged B.EQ is never taken and
  coins add no failure row.

Formulas are per node, times 2^72, except those with S, the Q broadcast, the barriers and the output. A record row uses
S, as each leaf gives at most one record.

| row (P: preprocessing) | log2 instructions |
| --- | ---: |
| node set-up: 2^20 (P) | 92.0000 |
| coins: (q1 + q2 + 1)*64 (P) | 105.8575 |
| Y generators: q2*2^21 (P) | 120.0000 |
| R_g and X generators: (q1 + 1)*2^21 | 119.6993 |
| header: 148 (P) | 79.2095 |
| cursor tables: 2*(3*2^36 + 1) (P) | 110.5850 |
| Q broadcast: 2*(2^36 + 2^72) messages at 40 | 78.3219 |
| d4 counters: 2^48 + 1 (P) | 120.0000 |
| d3 cursor table: 5*2^24 + 1 (P) | 98.3219 |
| forward: q1*(105*q2 + 3) | 132.4135 |
| inversion: q1*(23646 + 16) | 113.2296 |
| backward: q1*(523*q2 + 5) | 134.7300 |
| key pass without emits: q1*(53*q2 + 11) | 131.4272 |
| class-index products: cap_j*83 | 132.4897 |
| batch set-up: q1*64 | 104.6993 |
| hop 1, the emits and in-flight counts: q1*29*q2 | 130.5573 |
| hop 2, `h2s` and in-flight counts: S*16, plus 7*2^36 + 1 per node | 130.6993 |
| final cursors and in-flight counts, both hops: 2*(9*2^36 + 3) | 112.1699 |
| bucket pass `bkt`: S*16, plus 7*2^36 + 1 per node | 130.6993 |
| d4 prefix sums: 2^48/4*14 + 1 | 121.8074 |
| local pass `pass4`: S*7, plus 7*2^24 + 1 per node | 129.5066 |
| scan: S*7, plus 4 + 4*(2^20 + 1) per node | 129.5066 |
| candidates: 2^20*2^25 | 117.0000 |
| requests served: 2^21*2^10 | 103.0000 |
| phase-4 count exchanges: 4*2^41 | 115.0000 |
| barriers: 32*2^82 | 87.0000 |
| phase control: 2^10 | 82.0000 |
| output | 20.0000 |
| **total R** | **135.533473** |

`R = 63043548009814235640392594672626495913984` instructions, 456.4082 per leaf (the run's `per_leaf`), and
`W = R/1000 = 2^125.5676882`. **The claim, 125.567689, is the least multiple of 10^-6 with 2^claim >= W**; the
certificate encloses 2^c by an exact series for e^(f*ln 2) with interval bounds on ln 2. Rows marked P are
preprocessing.

## 8. Evidence

`experiment.py` has the declared input `{"seed": 1, "trials": 6000, "control_trials": 1000, "curve_trials": 100,
"seconds": 15}` (`seconds` 5 to 15); it finds the package files beside itself, under files/ or under /snapshot/files.
Budget-limited checks run guaranteed minimums (2,000 random and 600 near-tie class-rule inputs, 400 multi-node trials,
100 per control, 10 curve checks), then use their share of the remaining budget; larger counts are quoted only from the
census. **Part A** runs the leaf kernels on production batches, and fe_inv, against exact arithmetic. It checks each
port message and the cursor table against an exact model of the emit with the production 36-bit digits, a sentinel x30,
and the KEY model against the exact class rule. For the sites of 3.4 (`site_checks`), inputs built to meet a row-shift,
fold, SQP, SQM or KMUL predicate deviate exactly when it holds, other boundary and random inputs are exact, and each
matches the bit-level model of its arithmetic; FS, FADDW, `fe_mul_f` and `fe_sqr_f` are exact on edge, random and
small-residue inputs. A batch with `Y_l = ±X` gets acc = 0 and loses all its other leaves; the last R update is ≡ 1, and
overwriting it changes no message; and the SHX, FDX, SQP and SQM invocations per loop pass of fwd, bwd and keys match
the table (`sites_match_certificate`). **Part B** runs a 4 × 4 machine through `leaf.s` and `route.s`: `tinit`, `init3`,
`zero4`; `keys` on 12 planted incidences per node (d1, d2 < 4, some keys shared between nodes), whose emits fill pair
slots of 7 records; `fin`; `h2s` into pair slots of 7; `fin`; then `bkt` (four buckets of 8), `prefix4`, `pass4` and
`scan` at each home. Its links are Python code that writes each message's four words at its address, nothing else. Every
kernel call's count equals its formula from the certificate's constants; messages, slots, bucket entries, pointer order
and candidates equal an exact Python model of the same order; full slots (`keys`, `h2s`), a full bucket and a full list
are forced too. **Part C** runs the q25 toy in log form (the class of `[t]G` is `t^6 mod r`): 74 nodes of 4 × 4 grids,
one-pool 28 × 42 controls with b per X or a shared b, and a curve pipeline that must match the log form. **Part D** runs
`certificate.py`, and `counts_match_certificate` compares its count constants with parts A and B.

**Census.** `census-toy.json` is an off-runner run of the same code (`seconds` 1500, above the hosted maximum of 15;
each toy count reaches its request). Over 100000 multi-node trials, success is 0.39182 with all pairs and 0.38849 with
cross-node pairs only, against Poisson 0.39008 and a toy Chen–Stein bound of 0.0134. The one-pool controls give 0.3423
and 0.314 against their twin and shared-b models, 0.33759 and 0.30974, and the curve pipeline matches in 1000 of 1000
trials. The toy and the 4 × 4 run are evidence at their own scope only; a declared run's toy figure has its own standard
error.

## 9. Limitations and provenance

- **Model.** The bound is for the machine of §2.1. Its links write each message's four words at the address the sender
  stored and keep no state beyond the message in flight; the sender's instructions, including those that form the
  address, and the ports' in-flight counts are charged. On a machine
  whose links pass through intermediate processors, or whose delivery executes instructions, those would be added; that
  machine is not claimed.
- **Not executed:** the allowances of section 7 (generators, resolver, broadcasts, barriers, phase control and phase-4
  exchanges), the in-flight counts, and the production attack, machine and links.
- **Provenance**, not a premise: package F (submission 8f733826) for the row products (RF, RA, and the pair schedule of
  W2S and W2L), FS's CSEL form, fe_inv, the F1 row and its error analysis, and `kernels/a64.py`, whose extensions its
  docstring lists; Galbraith, Wang and Zhang for X + Y and X - Y from one inverse; Arratia, Goldstein and Gordon for
  Chen–Stein.
