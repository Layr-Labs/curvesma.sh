# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks and sort-based matching

Exploratory claim: a classical deterministic attack with a total-work bound at the scoring value `n = 2^20`, success
probability at least 0.39 under disclosed heuristics, and no advice. Not a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.7839` (exact `1167839/10000`); success `>= 0.39` (H1 predicts `0.39347`; the per-node caps and
node overflow cost at most `2^-35` over all `2^84` nodes and `2^-20` is budgeted for modelling error, so `P >= 0.3934`);
peak memory `<= 2^124.1` bytes (every node counted at its full `2^40` bytes); preprocessing `<= 2^115.80` units (included);
no advice. Unit = one batched incomplete affine
G1 addition = 2000 AArch64 instructions (`C = 2000`). Only `T` and `A_d`, `d = 2^20`, are read, so the bound holds for every
`n >= 2^20`.

## 1. Setting

`r` is the BLS12-381 scalar prime, `G1 = <P0>`. `r - 1 = 2^32 * 3 * 11 * 19 * 10177 * ...`; `d = 2^20`, `m = (r-1)/d`
(`log2 m = 234.857`), `H = {x : x^m = 1}` the subgroup of `d`-th powers. `mu6 = {u : u^6 = 1}`; since `6d | r-1`, `mu6 subset H`.
The GLV map `phi(x,y) = (beta x, y)` acts on `G1` as `[omega]`, `omega` a primitive cube root of unity; with negation the
orbit of `[v]P0` under `{+-1, +-phi, +-phi^2}` is `{[u v]P0 : u in mu6}`. `N` below is the number of elements per side and
`n_b` the number of base primes per side (not the input parameter `n`, not the input point `T`).

## 2. Reduction (Cheon, d | r-1)

For `s` uniform in `{1..r-1}`, `z := s^d` is uniform in `H` and `A_d = [z]P0`.

*Lemma 1.* Given `z`, the secret `s` is recovered with at most `2^21` units.
Proof. Take 20 successive square roots `z_{i+1}^2 = z_i`, `z_0 = z = w^{2^20}`. Each `z_i` (`i <= 19`) is `eps_i w^{2^(20-i)}` with
`eps_i` of order dividing `2^i`; both factors are squares (`2^19 | (r-1)/2`; `20 - i >= 1`), so all roots exist and `s_0 = z_20` has
`s_0^d = z`. The `d`-th roots are `s_0 omega_d^j`, `j < 2^20`; exactly one has `[s_j]P0 = T`. With `R = [s_0]P0` and
`j = j_1 + 2^10 j_2` we need `[omega_d^{j_1}]R = [omega_d^{-2^10 j_2}]T`; compute both sides for all `j_1, j_2 < 2^10` (`2^11` scalar
multiplications of at most 383 reference operations each) and match in a table: `< 2^20` units. The square roots cost `< 2^21`
instructions. QED

## 3. Baby and giant sets

**H-primes.** A prime `q` is an *H-prime* if `q^m = 1 mod r`. **H2:** the classes `q mod r` of primes `q << r` equidistribute in
`F_r^*/H` (index `d`): a fraction `1/d` of primes are H-primes, i.e. density `1/(d ln q)` among integers, about `X/(d ln X)` below
`X`, with gaps like independent exponentials of mean `d ln q`. Check on this `r`: of the 6542 primes below `2^16`, 402 (0.0614) are
16th powers (expected 0.0625) and 27 (0.0041) are 256th powers (expected 0.0039).

**Walk list `W`:** the first `L = 2^50` H-primes above `2^39`, `w_0 < ... < w_{L-1}`, even gaps `g_j`. **H2'** (checked by the
computation that builds `W`): `w_{L-1} < 2^77` and all `g_j <= G = 2^33` (under H2 a gap exceeds `2^33` with probability
`< exp(-2^7.3) < 2^-228`, over `2^50` gaps `< 2^-178`).
**Base primes:** H-primes in `(2^77,2^97)` (about `2^70` exist by H2; `2 n_b ~ 2^66.6` are needed) dealt alternately into a baby list
`a_i` and a disjoint giant list `a'_i`, `n_b = ceil(N/L)` each.
**Elements:** baby `(a,w)`: `x = a w`, point `[x]P0`; giant `(a',w')`: `x' = a' w'`, point `[x']A_d = [x' z]P0`.

*Lemma 2.* (i) All `x, x'` lie in `H`. (ii) `x` determines `(a, w)` (`a` is the only prime factor above `2^77`), so the `N` scalars of a side
are distinct. (iii) For pairs `(x,x')`, `(y,y')`, `x = aw, x' = a' w', y = bu, y' = b'u'`: `x/x' = y/y'` as rationals implies the pairs
coincide or `a = b, a' = b', w = w', u = u'` ("diagonal", ratio `a/a'`); diagonal pairs are a `1/L = 2^-50` fraction.
Proof. (i) closure; (ii) unique factorisation; (iii) `a w b' u' = b u a' w'`: the primes above `2^77` give `{a,b'} = {b,a'}`, and
baby/giant bases are disjoint, so `a = b`, `a' = b'`; then `w u' = u w'` gives `{w,u'} = {u,w'}`. QED

## 4. Walk, key, match

**Chain.** For a base prime (baby `Q = [a]P0`, giant `Q = [a']A_d`) build `tau[g] = [g]Q` for even `g <= G` (`tau[2]` by a doubling,
`tau[4] = 2tau[2]`, `tau[g] = tau[g-2] + tau[2]`) and walk `X_{j+1} = X_j + tau[g_j]`, `X_j = [w_j]Q`. For batching on one node the
walk is cut into `2^16` equal segments started by scalar multiplications `[w_s]Q` (the table likewise into `2^16` blocks); the segments
advance in lockstep and share one Montgomery inversion per round: the batched incomplete affine reference operation (6M).

*Lemma 3.* Every addition has distinct x-coordinates and no infinity operand.
Proof. `0 < g <= 2^33 < 2^39 < w_0 <= w_j < r/2` gives `w_j != +-g mod r`, so `X_j != +-tau[g]` as `Q` has order `r`; in the table
`g - 2 != +-2` for `g >= 6`, and likewise for block starts `g_0 >= 6`. QED

**Key.** `y*` = the even one of `y, p - y` (`p` odd); `kappa(Y)` = bits 1 to 224 of `y*`.
*Lemma 4.* Points in one `phi`/negation orbit have equal `kappa`; equal `y*` implies one orbit (`y' = +-y` gives `x'^3 = y'^2 - 4 = x^3`).
Different orbits agree on `kappa` only by chance (H3).

**Match.** A baby `[x]P0` and giant `[x'z]P0` in one orbit give `x = u x' z`, `u in mu6`, so `z = x (u x')^{-1}`. Try the six `u`; accept
the one with `[z]P0 = A_d` (at most 6 scalar multiplications of 383 reference operations: `<= 2^12` units). An orbit match always
yields the true `z`; an equal-`kappa` match of different orbits fails the test.

## 5. Machine model and matching

**Machine.** A shared-nothing cluster: `S = 2^84` owner nodes, `2^20` generating nodes, `2^31` storage nodes (for `W` and the base
lists). A node is a scalar AArch64 core with `2^40` bytes of local RAM, addressed by ordinary 64-bit loads and stores; **no
instruction addresses another node's memory**. Owner identifiers are 84-bit strings read as digits of 28, 28, 28 bits. Nodes
exchange blocks of `b = 25` records (1200 bytes) with any node differing from them in one digit; hardware moves the bytes,
software charges at most 60 instructions to post and 40 to poll a block. Posting is counted per destination and per block, not
per record: a node with `n_d` records for destination `d` posts `ceil(n_d/b)` blocks, hence at most `mu/b + 2^28` blocks in a
pass (Section 5.2). The `mu/b` part is the `100/25 = 4` per record of the routing row; the second part is the final blocks,
charged once per node per pass as fixed work (Section 7.2): at the end of a pass every destination receives exactly one *final*
block holding its `0..24` remaining records, which is also the end-of-pass marker (a node knows a pass is complete when its
64-bit final-block counter reaches the number of nodes that send to it in that pass, at most `2^28`). **Transport is priced as if the CPU copied the bytes:** one extra scalar copy of the
record at every hop (3 `ldp` + 3 `stp` = 6 instructions per 48-byte record per hop), on top of the sender-side copy into the
outgoing block. Arriving blocks land in the receiver's receive ring (`2^20` slots of 3 KiB); a forwarding node reads them there
(the 3 `ldp` of the next pass), and an owner copies its final-pass records out of the ring into its source buffer (the owner
drain, charged per record, Section 7.1). `W` (5 bytes per gap, `2^52.3` bytes) and the base lists stream to a generating node
in 1 MiB blocks: posting and polling under 0.001 instruction per step plus 1 instruction per 8 bytes of transport, i.e. 1 per
step; a label `(base index, step)` is resolved by two block reads (under 200 instructions, only inside the verification cap).

**Where the owner count, the digit widths and the block size come from.** A node can hold `C` records, so the owner count must
satisfy `S >= 2N/C`; the smallest power of two that fits is `2^84` (Section 5.1), and the owner id is then exactly 84 bits, read
as three 28-bit routing digits. The pass count `p` and digit width `w` (`p w >= 84`) are fixed by exact integers from two
opposing terms: every pass costs `17 + ceil(100/b) + 6` instructions per record, and every pass costs `2^w` final blocks of
107 instructions plus the table work of Section 7.2 per node, `110.875 * 2^w` instructions in all, which over `2^84` nodes is
`110.875 * 2^w * 2^84 / 2000` units.
The alternatives, each at its cheapest block size, with every other row of Sections 7.1 and 7.2 unchanged:

| passes x width | per-record row | per-node fixed row | `time_log2` |
| --- | ---: | ---: | ---: |
| 3 x 28, 25-record blocks (chosen) | 202 | `2^109.44` | `116.7839` |
| 4 x 21, 100-record blocks | 217 | `2^104.20` | `116.7851` |
| 5 x 17, 100-record blocks | 241 | `2^103.57` | `116.8005` |

so three passes of 28 bits are the cheapest. `w = 28` is the only width with `p = 3`: `w >= 29` needs more than `2^29` open
blocks of 25 records (`2^39.23` bytes) and leaves no room for the record buffers, and `w <= 27` needs four passes. The block
size is fixed the other way, by memory: the open-block area `2^28 * b * 48` bytes grows with `b`, `C/mu = 1.0858` at `b = 25`
and falls below 1 at `b = 30`, while the per-record charge `ceil(100/b)` is 4 for every `b` in `25..33`, so `b = 25` is the
cheapest block that fits. The local sort must raise the sorted prefix `P` from 84 bits to the value Section 6 needs; two
28-bit LSD passes give `P = 140`, the cheapest `P` whose run-pair cap is still affordable, while a 35-bit digit would not
fit: its two `2^35`-entry count and offset arrays are `2^39` bytes, which with the `2^38.229`-byte open-block area leaves less
than `96 mu` bytes for the record buffers, and its clearing and prefix sums (2.75 instructions per entry per pass, Section 7.2)
would add `2 x 2^35 x 2.75 = 2^37.46` instructions per node.

**Record**: 48 bytes: `kappa` 28 B, label 16 B (67-bit base index over the `2 n_b` base primes, 50-bit step, type bit), 4 B padding.
**Routing pass** (per record): the node keeps one open 25-record block per destination and appends to it, so no per-digit
counting loop and no counted scatter are charged. Digit `extr`+`ubfx` (2), `ldr` of the open block's offset in the block area
(1), 3 `ldp` + 3 `stp` (6), fill-count `ldr/add/str` (3), block-full test (2), loop (3) = 17 control, plus 4 for posting and
polling the full 25-record block (`100/25`), plus 6 for the transport copy: charged **27**. At the end of the pass every
destination gets its final block, `0..24` records (Sections 5.2 and 7.2).
**Owner drain** (per record, after the third pass): the poll of every block arriving in the third pass (inside its 40) returns the ring slot
and the record count `n`; the loop then runs `n` times: capacity test of the cursor against the last record slot of the
`C`-record buffer, before any write (`cmp`, `b.hi` to *overflow*: 2), 3 `ldp` from the slot and 3 `stp` to the source buffer at
the cursor, both post-indexed so the slot pointer and the cursor advance by 48 (6), loop (`subs`, `b.ne`: 2): charged **10**. The local sort reads only the source buffer.
**Local LSD pass**: counting loop (`ldr`, `extr`, `ubfx`, `ldr/add/str`, loop: 9) + scatter (`extr`, `ubfx`, `ldr` offset, `madd`,
`add`, 6 copies, `add/str`, loop: 16) = 25: charged **25**.

**Matching.** Each record is routed by 3 routing passes on the first 84 bits of `kappa` to its owner (about `2^32.64` records, one
84-bit prefix), sorted by 2 LSD passes on `kappa[85..140]` (28-bit digits) and scanned once. Records with equal 140-bit prefix are
adjacent (a run); inside a run all opposite-type pairs are compared on the full 224-bit `kappa`, and pairs with equal `kappa` are
verified by the Match step. **Caps, per node:** every owner keeps two 64-bit counters, run-pair comparisons and verifications; a
node that would make its `2^20 + 1`-th comparison or its `129`-th verification stops and signals *fail* (the counters never exceed
`2^20 + 1` and `129`, far below `2^64`). No global counter and no combining step exist: each cap is local, and the failure of any
node is a failure of the run. They are charged at the caps on every owner: `2^84 * 2^20` comparisons of 32 instructions
(`2^98.03` units) and `2^84 * 2^7` verifications of `2^12` units (`2^103` units). A comparison, itemised: pair loop (index
`add`, `cmp`, `b`: 3), partner record address (1), partner type bit (`ldrb`, `tbz`: 2), both 28-byte keys (`ldp`, `ldr`, `ldr w`
each: 6), 224-bit compare (4 `eor`, 3 `orr`, `cbnz`: 8), counter (`add`, `cmp`, `b.hi`: 3): 23, charged 32. All other work is
fixed: the number of records and passes does not depend on the data.

### 5.1 Placement, capacity and overflow

`2N = 2^116.6361` records, `mu = 2N/S = 2^32.636 = 6.675 * 10^9` per owner. **Every node that forwards records**, owners included,
holds the open-block area of a 28-bit pass: `2^28` blocks of `b = 25` records, `2^28 * 1200 B = 2^38.229` bytes. Per-node
allocation (bytes), all concurrently resident: the open-block area `2^38.229`; two record buffers `C * 96` (48 bytes source +
48 bytes destination); LSD count and offset arrays `2 * 2^28 * 8 = 2^32` (8-byte entries: the cumulative offset reaches `C`);
routing offset and cursor tables `2 * 2^28 * 8 = 2^32`; receive ring `2^20 * 3 KiB = 2^31.6`; code and stack `<= 2^30`; reserved
`2^36`. The auxiliary total is `403726925824 = 2^38.5546` bytes, so `C = (2^40 - 4.037e11)/96 = 7247757312 = 2^32.7549` records,
`96 C + aux = 2^40` exactly, and `C/mu = 1.0858`. A node that would exceed `C` signals *overflow* and the attack fails (the
drain's capacity test; the cursor is a byte address below `2^40`, the record count at most `C + 1 < 2^33`, both 64-bit); by the
Chernoff bound of Section 6 this costs far less than the `2^-20` budgeted for modelling error.
A generating node holds a chain table `2^32 * 96 B = 2^38.585`, the same open-block area `2^38.229`, stream buffers `<= 2^30` and
code `<= 2^30`: `<= 2^39.4 < 2^40` bytes. A storage node holds a share of `W` and the base lists.

### 5.2 Blocks sent per pass

A pass routes on one 28-bit digit, so a node has `2^28` destinations. With `n_d` records for destination `d` and
`sum_d n_d <= mu`, the blocks it posts are `sum_d ceil(n_d/b) <= mu/b + #{d : n_d > 0} <= mu/b + 2^28`, for every possible
distribution of the records over the destinations: `ceil` wastes at most one block per non-empty destination and there are at
most `2^28` of them. The `mu/b` term is the `100/25 = 4` of each routing pass; the second is the final block, sent to every
destination (empty or not) once per pass, charged as `2^28` blocks of 107 instructions of fixed work per node per pass
(Section 7.2). No concentration assumption enters the posting bound: there is exactly one final block per destination whatever
the counts are, and the record counts per node enter only through the overflow rule of Section 5.1.

## 6. Algorithm and success probability

```text
1. Build W and the base lists; read A_d, T.
2. Generate all N baby and N giant records (each generating node walks one chain as 2^16 segments); route each record as produced.
3. Match as in Section 5 with the caps; on the first verified z run Lemma 1 and output s.
```

`K = 0.50`, `N = n_b L` with `n_b = ceil(sqrt(K m/6)/L)`: `N = 2^115.6361`, `6N^2/m = 0.50`. Without caps and overflow the attack
succeeds exactly on `S = {x/(u x')}`, a union of `6N^2` values of `H`, so `P(success) = |S|/m`.
**H1:** these values behave like independent uniform elements of `H`: `|S|/m = 1 - e^{-0.50} = 0.39347`. Support: by Lemma 2(iii) they
are distinct as rationals except a `2^-50` diagonal; within a row/column they are distinct (negative correlation); the other
coincidences are modular events among integers below `2^174` reduced modulo a 255-bit prime with no known structure.
**H3 (random-function form, on orbits):** `kappa(record) = F(orbit of the record)` for a uniformly random function `F` from
orbits to 224-bit strings, independent of which orbits occur (a heuristic about coordinates of structured points, as H1; the baby
keys are fixed, not randomised by the challenge). Records of one orbit share their key (Lemma 4); keys of different orbits are
independent and uniform. An orbit carries at most 12 records: its 6 elements of `H`, each at most once per side (Lemma 2(ii),
and distinct `x'` give distinct `x' z`), so at most 6 baby and 6 giant records, and two orbits give at most `6*6 + 6*6 = 72`
opposite-type record pairs.
**Overflow.** Orbits fall independently and uniformly on the `2^84` owners (key bits 1..84), each adding 1 to 12 records, so the
Chernoff bound for independent `[0,12]`-bounded summands with mean `mu` gives `P(load >= C) <= exp(-mu (a ln a - a + 1)/12)`,
`a = C/mu = 1.0858`, i.e. `2^-2.876e6` per node and still `2^-2.876e6` after the union over `2^84` nodes (`evidence.py` prints the
exponent).
**Per-node caps.** Fix an owner holding `M <= C` records (else it has already failed by overflow), hence `M' <= C` orbits. Given
which orbits it holds (key bits 1..84), their key bits 85..140 are independent and uniform on `2^56` values. A compared pair is
an opposite-type pair in one run; a same-orbit one has equal full keys, is verified, yields `z` and ends the attack, so at most one is
compared, and every other compared pair joins two different orbits with equal bits 85..140.
(A) Some value is held by at least 7 orbits: probability `<= C(M',7) * 2^(-56*6) <= C^7/(7! * 2^336) = 2^-119.01` per
node, `2^-35.01` over `2^84` nodes.
(B) Let `D'` be the number of orbits whose value repeats that of an earlier orbit (any fixed order). Each orbit repeats with
probability at most `M'/2^56` whatever the earlier ones did, so `D'` is dominated by `Bin(M', M'/2^56)`, mean
`<= C^2/2^56 = 729.0`. Outside (A) every value holds `j <= 6` orbits and `C(j,2) <= 3(j-1)`, so the orbit pairs with equal value
are at most `3 D'` and the compared pairs at most `1 + 72 * 3 D' = 1 + 216 D'`. Exceeding `2^20` therefore needs `D' >= 4855`, and
`P(Bin >= k) <= (e * 729.0/k)^k` gives `2^-6276` per node, `2^-6192` over all nodes.
(V) Verifications at a node are at most `1 + 72 * Q`, `Q` the number of orbit pairs with equal 224-bit keys. The cap 128 is
exceeded only if `Q >= 2`, i.e. two distinct orbit pairs each agreeing on the 140 bits 85..224 (disjoint or sharing one orbit,
probability `2^-280` either way): `<= C(C(M',2),2) * 2^-280 <= C^4/8 * 2^-280 = 2^-151.98` per node, `2^-67.98` over all nodes.
**Success.** `P(success) >= 0.393469 - 2^-35.01 - 2^-6192 - 2^-67.98 - 2^-2.876e6 - 2^-20 > 0.393468`; `2^-20` is budgeted for
modelling error. **0.39** is claimed, `0.0034` above the requirement.

**Toy check** (`evidence.py coverage`, scalar level, three toy primes with `6 * 2^k | r-1`, chain length `L = 60`, `K = 0.50`, seed 2026,
`8000/5000/3000` trials, the same trials the manifest's declared `mode: all` run executes). At toy scale the diagonal pairs are a `1/L` fraction, so the Poisson
prediction uses `K' = K(1 - (L-1)/L^2)`. Measured success for `r ~ 2^27.6, 2^29.7, 2^31.7`: `0.4255, 0.4242, 0.3850` against
Poisson(`K'`) `0.4249, 0.4192, 0.3894` (|z| <= 0.72; Poisson(`K`) is `0.4302, 0.4245, 0.3944`). Evidence at toy scale only: it
supports the Poisson form and the duplicate bookkeeping, not the production value of `m`.

## 7. Cost ledger (unit = 2000 instructions)

### 7.1 Per record, on top of the reference addition

| item | instructions |
| --- | ---: |
| walk bookkeeping, 28 itemised: gap load and copy 3, gap-stream transport 1, table address `base + (gap/2)*96` 3, segment loop 3, label 3, batch pointer/index updates outside the Montgomery-trick multiplications 12, stream-refill check 3 | 28 |
| key, 23 itemised: `ldp` of 6 `y` limbs 3, parity test and select mask 2, `p - y` by `subs`, 5 `sbcs`, 6 `csel` 12, `extr` to bits 1..224 4, `stp` of the key 2 | 23 |
| 3 routing passes x 27 (17 control + 4 for posting a full 25-record block + 6 transport copy) | 81 |
| owner drain, receive ring to source buffer (3 `ldp` + 3 `stp` 6, capacity test 2, loop 2) | 10 |
| 2 local LSD passes x 25 (counting 9 + scatter 16) | 50 |
| run scan (load the 140-bit sort prefix, compare with the previous record, branch) | 10 |
| **itemised; charged** | **202** |

Element cost `c = 1 + 202/2000 = 1.101`.

**Posting is not an average.** The `4` in each routing pass is the `mu/b` part of `sum_d ceil(n_d/b) <= mu/b + 2^28`
(Section 5.2). The `2^28` part is not per record but per node and per pass, and it is the final-block row of Section 7.2. Charging `100/b = 4` per record and one final block per destination bounds every distribution of the
records over the destinations, so this row carries no expectation and no concentration assumption.

**Why no percentage allowance.** The earlier form of this row listed 236 instructions and charged 260, about 10% for "unlisted loop
and spill instructions". Those instructions are named here instead of covered by a percentage: the loop counters are the "loop 3" in
each pass row and the "segment loop 3" of the walk row; the spill and reload traffic is the "batch pointer/index updates 12" of the
walk row and the 3 `ldp` + 3 `stp` copies inside each pass row; the block-boundary work is the "block-full test 2" and the
amortised posting. Everything a record does not touch per unit of work is fixed work and is in the rows of Section 7.2: clearing and
prefix-summing the histograms and the routing tables, the final blocks of Section 5.2, the receive ring and node setup (the
itemised per-node rows), and the chain tables. For orientation, `+10` unlisted instructions per record would move the bound by
0.0066 bits.

### 7.2 Fixed terms

| term | units |
| --- | ---: |
| elements `2Nc` | `2^116.774878` |
| tables: `2 n_b` chains x `2^32` entries x (1 addition + 12 stores/2000) | `2^98.645` |
| segment and table-block starts, `[a]P0` (`<= 2^25` per chain) | `2^91.636` |
| per-node fixed work, `(2^84 + 2^20 + 2^31) nodes x 90768936960 = 2^36.4015` instructions, itemised below | `2^109.436` |
| H-prime search in `(2^77,2^97)`: `2^86.6` primes x `353 * 383/2000 = 67.6`, `2^88.3` Miller-Rabin candidates x 27.8, sieve, 4x margin | `2^97` |
| H-primes of `W` below `2^77`: `2^70` primes, same pricing, 4x margin | `2^84` |
| verification cap, per node: `2^84` owners x `2^7` verifications x `2^12` units | `2^103` |
| run-pair cap, per node: `2^84` owners x `2^20` comparisons x 32/2000 | `2^98.034` |
| input access `2^20`; Lemma 1 `2^21` | `<= 2^21` each |

**Per-node fixed rows** (instructions per node; every node, generating and storage nodes included, is charged the owner's rows).
Every array is 8-byte entries, `2^28` of them; each loop handles 16 entries per iteration and pays `subs` + `b.ne` (2) per
iteration, and `2^28` is a multiple of 16, so there is no remainder.

| array and pass | per 16 entries | entries x passes | instructions |
| --- | ---: | ---: | ---: |
| routing fill-count table, clear (8 `stp xzr, xzr` + 2) | 10 | `2^28 x 3` | 503316480 |
| routing offset table, `offset[v] = 1200 v` (per 2 entries `stp`, 2 `add`: 24 + 2) | 26 | `2^28 x 3` | 1308622848 |
| final-block scan (per 2 entries `ldp`, 2 `cbnz`: 24 + 2) | 26 | `2^28 x 3` | 1308622848 |
| final block per destination: call and return 2, block address 1, count header 1, post 60, poll 40, final-block counter `add`, `cmp`, `b.eq` 3 | `107 x 16` | `2^28 x 3` | 86167781376 |
| LSD count array, clear (8 `stp xzr, xzr` + 2) | 10 | `2^28 x 2` | 335544320 |
| LSD exclusive prefix sum into the offset array (per 2 entries `ldp`, `add`, `stp`, `add`: 32 + 2) | 34 | `2^28 x 2` | 1140850688 |
| receive ring, `2^20` slot descriptors (`str`, `add`, `subs`, `b.ne` per slot) | - | `2^20` slots | 4194304 |
| node setup: code, stack, cursors and the three 64-bit counters | - | once | 4096 |
| **per node** | | | **90768936960** |

The 8-byte entries are wide enough everywhere: a fill count is at most 25, an offset below `2^39`, an LSD count or offset at most
`C < 2^33`, a final-block counter at most the number of senders of a pass, at most `2^28`.

All F_r multiplications in prime searches are priced at the 383 instructions measured for `blst_fp_mul` (381-bit), an upper bound. The
rows other than `2Nc` sum to `2^109.4539`; `W = 2^116.774878 + 2^109.4539 = 2^116.783872 <= 2^116.7839`, score `1167839/10000`
(`evidence.py ledger` recomputes every row with exact integers and checks `W <= 2^116.7839` by integer exponentiation).
Preprocessing (lists, tables, baby generation) `<= Nc + 2^109.454 = 2^115.7928 <= 2^115.80`.

Memory counted at full allocation: `(2^84 + 2^20 + 2^31) * 2^40 = 2^124.0`, claimed `<= 2^124.1` (records alone: `2N * 48 B * 2 =
2^123.24`; `W` `2^52.3` and base lists `2^70.2` bytes fit the `2^31` storage nodes). Memory has no scalar weight and is exponential,
as for any BSGS-type Cheon attack.

## 8. Heuristics, evidence, provenance

H1 decides success; H3 decides the failure probability of the caps and of node overflow (supported only by the toy occupancy test and
the absence of known structure); H2/H2' decide preprocessing and table size (checked by building `W`). The cluster model is idealised
and the instruction counts are listed, not measured on hardware. The cost is data independent except for the two caps, charged at
their limits.

`evidence.py` (one initial run `mode: all`, about 6 s; each mode can be rerun with more trials): `ledger` (exact integers:
Sections 5.1-7, the digit and block widths, the blocks-sent bound `mu/b + 2^w` and the final blocks, every per-node row, the
per-node caps and their bounds (A), (B), (V), capacity, Chernoff exponent, the success line, memory, `6d | r-1`,
`mu6 subset H`, H2 check); `coverage` (Section 6); `ec`: Algorithm A end to end on
genuine `j = 0` curves of prime order `r ~ 2^28..2^32`: affine
incomplete additions with Montgomery batch inversion (asserting no zero denominator), the parity key truncated to 18/20/22 bits so
that false key matches occur (0.5-1.8 per run) and are verified and rejected under a cap, sort, equal-key runs, `z`, `d`-th root;
every recovered `z`, `s` is checked, success agrees with the scalar model in every trial, and keys binned into 32 pseudo-nodes show
chi-square z-scores in `[-0.8, +1.0]`. Limits: tables are built only for occurring gaps (so `tau[g] = tau[g-2] + tau[2]` and the
segments rest on Lemma 3), and the toy sorts locally without distributed routing.

Provenance: the idea of walking small integers of `H` with one addition per step, the GLV six-fold key and Cheon's `d | r-1`
reduction come from the published note of submission `cbadde77`; the semiprime construction, Lemmas 2-3, matching, caps, machine
model and overhead listing were derived here.
