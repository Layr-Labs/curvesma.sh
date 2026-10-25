# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks and sort-based matching

Exploratory claim: a classical deterministic attack with a total-work bound at the scoring value `n = 2^20`, success
probability at least 0.39 under disclosed heuristics, and no advice. Not a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.7764` (exact `291941/2500`); success `>= 0.39` (H1 predicts `0.39347`; the two caps cost at
most `0.00060` and node overflow is far below `2^-20`, so `0.3928` is claimed); peak memory `<= 2^124.1` bytes (every node
counted at its full `2^40` bytes); preprocessing `<= 2^115.79` units (included); no advice. Unit = one batched incomplete affine
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
pass (Section 5.2). The `mu/b` part is the `100/25 = 4` per record of the routing row; the second part is a partial-block flush
charged once per node per pass as fixed work (Section 7.2). **Transport is priced as if the CPU copied the bytes:** one extra
scalar copy of the record at every hop (3 `ldp` + 3 `stp` = 6 instructions per 48-byte record per hop), on top of the
sender-side copy into the outgoing block. `W` (5 bytes per gap, `2^52.3` bytes) and the base lists stream to a generating node
in 1 MiB blocks: posting and polling under 0.001 instruction per step plus 1 instruction per 8 bytes of transport, i.e. 1 per
step; a label `(base index, step)` is resolved by two block reads (under 200 instructions, only inside the verification cap).

**Where the owner count, the digit widths and the block size come from.** A node can hold `C` records, so the owner count must
satisfy `S >= 2N/C`; the smallest power of two that fits is `2^84` (Section 5.1), and the owner id is then exactly 84 bits, read
as three 28-bit routing digits. The pass count `p` and digit width `w` (`p w >= 84`) are fixed by exact integers from two
opposing terms: every pass costs `17 + 100/b + 6` instructions per record, and every pass costs a partial-block flush of at most
`100 * 2^w` instructions per node, which over `2^84` nodes is `100 * 2^w * 2^84 / 2000` units. The alternatives, each at its
cheapest block size:

| passes x width | per-record row | flush | `time_log2` |
| --- | ---: | ---: | ---: |
| 3 x 28 (chosen) | 192 | `2^109.26` | `116.7764` |
| 4 x 21, 64-record blocks | 211 | `2^102.68` | `116.7810` |
| 5 x 17 | 246 | `2^99.00` | `116.8035` |

so three passes of 28 bits are the cheapest. `w = 28` is the only width with `p = 3`: `w >= 29` needs more than `2^29` open
blocks of 25 records (`2^39.23` bytes) and leaves no room for the record buffers, and `w <= 27` needs four passes. The block
size is fixed the other way, by memory: the open-block area `2^28 * b * 48` bytes grows with `b`, `C/mu = 1.0858` at `b = 25`
and falls below 1 at `b = 30`, while the per-record charge `ceil(100/b)` is 4 for every `b` in `25..33`, so `b = 25` is the
cheapest block that fits. The local sort must raise the sorted prefix `P` from 84 bits to the value Section 6 needs; two
28-bit LSD passes give `P = 140`, the cheapest `P` whose run-pair cap is still affordable, while a 35-bit digit would not
fit: its `2^35`-entry count and offset arrays are `2^38` bytes, which with the `2^38.229`-byte open-block area leaves less
than `96 mu` bytes for the record buffers, and its clearing and prefix sums would add `2^38` instructions per node.

**Record**: 48 bytes: `kappa` 28 B, label 16 B (66-bit base index, 50-bit step, type bit), 4 B padding.
**Routing pass** (per record): the node keeps one open 25-record block per destination and appends to it, so no per-digit
counting loop and no counted scatter are charged. Digit `extr`+`ubfx` (2), `ldr` of the open block's offset in the block area
(1), 3 `ldp` + 3 `stp` (6), fill-count `ldr/add/str` (3), block-full test (2), loop (3) = 17 control, plus 4 for posting and
polling the full 25-record block (`100/25`), plus 6 for the transport copy: charged **27**. At the end of the pass every
destination still holding `1..24` records posts one partial block (Section 7.2).
**Local LSD pass**: counting loop (`ldr`, `extr`, `ubfx`, `ldr/add/str`, loop: 9) + scatter (`extr`, `ubfx`, `ldr` offset, `madd`,
`add`, 6 copies, `add/str`, loop: 16) = 25: charged **25**.

**Matching.** Each record is routed by 3 routing passes on the first 84 bits of `kappa` to its owner (about `2^32.64` records, one
84-bit prefix), sorted by 2 LSD passes on `kappa[85..140]` (28-bit digits) and scanned once. Records with equal 140-bit prefix are
adjacent (a run); inside a run all opposite-type pairs are compared on the full 224-bit `kappa`, and pairs with equal `kappa` are
verified by the Match step. **Caps:** after more than `2^104` run-pair comparisons or `2^20` verifications the attack stops with
*fail*; they are charged at the caps (`2^98.03` and `2^32` units). All other work is fixed: the number of records and passes does
not depend on the data.

### 5.1 Placement, capacity and overflow

`2N = 2^116.6361` records, `mu = 2N/S = 2^32.636 = 6.675 * 10^9` per owner. **Every node that forwards records**, owners included,
holds the open-block area of a 28-bit pass: `2^28` blocks of `b = 25` records, `2^28 * 1200 B = 2^38.229` bytes. Per-node
allocation (bytes), all concurrently resident: the open-block area `2^38.229`; two record buffers `C * 96` (48 bytes source +
48 bytes destination); LSD count and offset arrays `2 * 2^28 * 8 = 2^32` (8-byte entries: the cumulative offset reaches `C`);
routing offset and cursor tables `2 * 2^28 * 8 = 2^32`; receive pool `2^20 * 3 KiB = 2^31.6`; code and stack `<= 2^30`; reserved
`2^36`. The auxiliary total is `403726925824 = 2^38.5546` bytes, so `C = (2^40 - 4.037e11)/96 = 7247757312 = 2^32.7549` records,
`96 C + aux = 2^40` exactly, and `C/mu = 1.0858`. A node that would exceed `C` signals *overflow* and the attack fails; by the
Chernoff bound of Section 6 this costs far less than the `2^-20` budgeted for modelling error.
A generating node holds a chain table `2^32 * 96 B = 2^38.585`, the same open-block area `2^38.229`, stream buffers `<= 2^30` and
code `<= 2^30`: `<= 2^39.4 < 2^40` bytes. A storage node holds a share of `W` and the base lists.

### 5.2 Blocks sent per pass

A pass routes on one 28-bit digit, so a node has `2^28` destinations. With `n_d` records for destination `d` and
`sum_d n_d <= mu`, the blocks it posts are `sum_d ceil(n_d/b) <= mu/b + #{d : n_d > 0} <= mu/b + 2^28`, for every possible
distribution of the records over the destinations: `ceil` wastes at most one block per non-empty destination and there are at
most `2^28` of them. The `mu/b` term is the `100/25 = 4` of each routing pass; the second is charged once per node per pass as
`100 * 2^28 = 2^36.229` instructions of fixed work (Section 7.2). No concentration assumption enters the posting bound: the
flush is at most one block per destination whatever the counts are, and the record counts per node enter only through the
overflow rule of Section 5.1.

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
**H3 (random-function form):** the 140-bit key prefixes of the `2N` records behave like independent uniform strings (a heuristic
about coordinates of structured points, as H1; the baby keys are fixed, not randomised by the challenge). Used: (a) records of
different orbits agree on `b` prescribed bits with probability `2^-b`; (b) node occupancies are binomial with mean `mu`.
**Caps, overflow, success.** By (a) the expected number of record pairs with equal 140-bit prefix is `2^92.3`, so the `2^104` cap is
exceeded with probability `< 2^-11.7` (Markov); equal full `kappa` from different orbits: `2^8.3` pairs, so the `2^20` cap fails
with probability `< 2^-11.7`. By (b), Chernoff for `C = 1.0858 mu` gives `2^-3.45e7` per node, union over `2^84` nodes still
`2^-3.45e7` (`evidence.py` prints the exponent); `2^-20` is budgeted for modelling error. So
`P(success) >= 0.39347 - 2^-10.7 - 2^-20 > 0.3928`; **0.39** is claimed, `0.0028` above the requirement.

**Toy check** (`evidence.py coverage`, scalar level, three toy primes with `6 * 2^k | r-1`, chain length `L = 60`, `K = 0.50`, seed 2026,
8k/5k/3k trials run locally; the manifest run uses fewer). At toy scale the diagonal pairs are a `1/L` fraction, so the Poisson
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
| 2 local LSD passes x 25 (counting 9 + scatter 16) | 50 |
| run scan (load the 140-bit sort prefix, compare with the previous record, branch) | 10 |
| **itemised; charged** | **192** |

Element cost `c = 1 + 192/2000 = 1.096`.

**Posting is not an average.** The `4` in each routing pass is the `mu/b` part of `sum_d ceil(n_d/b) <= mu/b + 2^28`
(Section 5.2). The `2^28` part is not per record but per node and per pass, and it is the `2^36.229` instruction flush row of
Section 7.2. Charging `100/b = 4` per record and one flush block per destination bounds every distribution of the
records over the destinations, so this row carries no expectation and no concentration assumption.

**Why no percentage allowance.** The earlier form of this row listed 236 instructions and charged 260, about 10% for "unlisted loop
and spill instructions". Those instructions are named here instead of covered by a percentage: the loop counters are the "loop 3" in
each pass row and the "segment loop 3" of the walk row; the spill and reload traffic is the "batch pointer/index updates 12" of the
walk row and the 3 `ldp` + 3 `stp` copies inside each pass row; the block-boundary work is the "block-full test 2" and the
amortised posting. Everything a record does not touch per unit of work is fixed work and is in the rows of Section 7.2: clearing and
prefix-summing the histograms and the routing tables (`per_node_fixed`), the partial-block flush of Section 5.2, the receive pool,
node setup, and the chain tables. For orientation, `+40` unlisted instructions per record would move the
bound by 0.026 bits, so an unlisted-item error of up to 20 instructions would still leave the claim below `2^116.79`.

### 7.2 Fixed terms

| term | units |
| --- | ---: |
| elements `2Nc` | `2^116.7683` |
| tables: `2 n_b` chains x `2^32` entries x (1 addition + 12 stores/2000) | `2^98.645` |
| segment and table-block starts, `[a]P0` (`<= 2^25` per chain) | `2^91.636` |
| per-node fixed work, `(2^84 + 2^20 + 2^31) nodes x 2^36.248` instr = partial-block flush `3 x 2^28` blocks x 100 instr plus clearing and prefix-summing the routing and LSD tables, receive pool, node setup | `2^109.282` |
| H-prime search in `(2^77,2^97)`: `2^86.6` primes x `353 * 383/2000 = 67.6`, `2^88.3` Miller-Rabin candidates x 27.8, sieve, 4x margin | `2^97` |
| H-primes of `W` below `2^77`: `2^70` primes, same pricing, 4x margin | `2^84` |
| input access `2^20`; Lemma 1 `2^21`; verification cap `2^20` verifications x `2^12` units `= 2^32`; run-pair cap `2^104` pairs x 32/2000 `= 2^98.034` | `<= 2^98.034` each |

All F_r multiplications in prime searches are priced at the 383 instructions measured for `blst_fp_mul` (381-bit), an upper bound. The
rows other than `2Nc` sum to `2^109.2840`; `W = 2^116.7683 + 2^109.2840 = 2^116.776346 <= 2^116.7764`, score `291941/2500`
(`evidence.py ledger` recomputes every row with exact integers and checks `W <= 2^116.7764` by integer exponentiation).
Preprocessing (lists, tables, baby generation) `<= Nc + 2^109.284 = 2^115.7843 <= 2^115.79`.

Memory counted at full allocation: `(2^84 + 2^20 + 2^31) * 2^40 = 2^124.0`, claimed `<= 2^124.1` (records alone: `2N * 48 B * 2 =
2^123.24`; `W` `2^52.3` and base lists `2^70.2` bytes fit the `2^31` storage nodes). Memory has no scalar weight and is exponential,
as for any BSGS-type Cheon attack.

## 8. Heuristics, evidence, provenance

H1 decides success; H3 decides the failure probability of the caps and of node overflow (supported only by the toy occupancy test and
the absence of known structure); H2/H2' decide preprocessing and table size (checked by building `W`). The cluster model is idealised
and the instruction counts are listed, not measured on hardware. The cost is data independent except for the two caps, charged at
their limits.

`evidence.py` (one initial run `mode: all`, about 6 s; each mode can be rerun with more trials): `ledger` (exact integers:
Sections 5.1-7, the digit and block widths, the blocks-sent bound `mu/b + 2^w` and the flush, cap probabilities, capacity,
Chernoff exponent, memory, `6d | r-1`, `mu6 subset H`, H2 check); `coverage` (Section 6); `ec`: Algorithm A end to end on
genuine `j = 0` curves of prime order `r ~ 2^28..2^32`: affine
incomplete additions with Montgomery batch inversion (asserting no zero denominator), the parity key truncated to 18/20/22 bits so
that false key matches occur (0.5-1.8 per run) and are verified and rejected under a cap, sort, equal-key runs, `z`, `d`-th root;
every recovered `z`, `s` is checked, success agrees with the scalar model in every trial, and keys binned into 32 pseudo-nodes show
chi-square z-scores in `[-0.8, +1.0]`. Limits: tables are built only for occurring gaps (so `tau[g] = tau[g-2] + tau[2]` and the
segments rest on Lemma 3), and the toy sorts locally without distributed routing.

Provenance: the idea of walking small integers of `H` with one addition per step, the GLV six-fold key and Cheon's `d | r-1`
reduction come from the published note of submission `cbadde77`; the semiprime construction, Lemmas 2-3, matching, caps, machine
model and overhead listing were derived here.
