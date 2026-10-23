# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks and sort-based matching

Exploratory claim under `rules.md`: a classical deterministic attack with a stated total-work bound at the
scoring value `n = 2^20`, success probability at least 0.40 under the disclosed heuristics, and no advice.
Nothing here is a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.8661` (exact rational `1168661/10000`); success `>= 0.40` (H1 predicts
`0.4055`, the caps and the overflow event cost at most `0.0007`, rules require `0.39`); peak memory `<= 2^123.3` bytes; preprocessing
`<= 2^115.87` units (included in the total); no advice. Unit = one batched incomplete affine G1 addition =
2000 AArch64 instructions (`C = 2000`). The bound holds for every `n >= 2^20` and does not depend on `n`.

## 1. Setting and notation

* `r` is the BLS12-381 scalar prime, `p` the base-field prime, `G1 = <P0>` has order `r`.
* `r - 1 = 2^32 * 3 * 11 * 19 * 10177 * ...`. Put `d = 2^20` and `m = (r-1)/d` (`log2 m = 234.857`).
  `H = {w^d : w in F_r^*} = {x : x^m = 1}` is the subgroup of order `m`.
* `mu6 = {u : u^6 = 1} = {+-1, +-omega, +-omega^2}`, `omega` a primitive cube root of unity mod `r`. Since
  `6d | r-1`, `mu6 subset H` (checked by `evidence.py ledger`).
* The GLV map `phi(x,y) = (beta x, y)` acts on `G1` as `[omega]` for one primitive cube root `omega`.
  With negation, the orbit of `[v]P0` under `{+-1, +-phi, +-phi^2}` is `{[u v]P0 : u in mu6}`.
* The input gives `T = A_1 = [s]P0` and `A_i = [s^i]P0`, `i <= n`. The attack reads only `T` and `A_d`,
  `d = 2^20`. (`N` below is the number of elements per side, `n_b` the number of base primes per side; neither
  is the input parameter `n` or the input point `T`.)

## 2. Reduction (Cheon, variant d | r-1)

For `s` uniform in `{1..r-1}`, `z := s^d` is uniform in `H` (`s -> s^d` is `d`-to-1 onto `H`) and
`A_d = [z]P0`. The rest of the attack searches for `z`.

*Lemma 1.* Given `z in H`, the secret `s` is recovered with at most `2^21` units.

Proof. Let `z_0 = z` and let `z_{i+1}` be any square root of `z_i`. Writing `z = w^{2^20}`, each `z_i`
(`i <= 19`) is `eps_i * w^{2^(20-i)}` with `eps_i` a root of unity of order dividing `2^i`; both factors
are squares (`eps_i` because its order divides `2^19`, which divides `(r-1)/2`; the other because `20-i >= 1`),
so all 20 square roots exist and `s_0 := z_20` has `s_0^d = z`. All `d`-th roots are `s_0 omega_d^j`,
`j < 2^20`, for a primitive `d`-th root of unity `omega_d`, and exactly one satisfies `[s_j]P0 = T`
because `P0` has order `r`. Put `R = [s_0]P0`; with `j = j_1 + 2^10 j_2` we need
`[omega_d^{j_1}]R = [omega_d^{-2^10 j_2}]T`. Compute both sides for all `j_1, j_2 < 2^10`
(`2^11` scalar multiplications by double-and-add, each at most 383 reference operations because the scalars
lie below `2^255`) and match them in a table. Cost `< 2^11 * 383 < 2^20` units; 20 Tonelli-Shanks square
roots cost `< 2^21` instructions. QED

## 3. Baby and giant sets

**H-primes.** A prime `q` is an *H-prime* if `q^m = 1 mod r`. Heuristic **H2**: for primes `q << r` the classes
`q mod r` equidistribute in `F_r^*/H` (index `d`), so H-primes have density `1/(d ln q)` among primes, and the
gaps between consecutive H-primes behave like independent exponentials of mean `d ln q`. Sanity check on this
`r` (`evidence.py ledger`): of the 6542 primes below `2^16`, 402 (0.0614) are 16th powers (expected 0.0625)
and 27 (0.0041) are 256th powers (expected 0.0039).

**Walk list `W`.** The first `L = 2^50` H-primes above `2^39`: `w_0 < w_1 < ... < w_{L-1}`. Gaps
`g_j = w_{j+1} - w_j` are even. The attack computes `W` in preprocessing. **H2'** (a property of the fixed
constants `r`, checked by that very computation): `w_{L-1} < 2^77` and every gap is at most `G = 2^33`.
Under H2 the probability that one gap exceeds `2^33` is about `exp(-2^7.3) < 2^-228`, so over `2^50` gaps H2'
fails with probability below `2^-178`.

**Base primes.** H-primes in `(2^77, 2^97)` in increasing order, dealt alternately into a baby list
`a_1, a_2, ...` and a disjoint giant list `a'_1, a'_2, ...`, `n_b = ceil(N/L)` of each, `N` as in Section 6.
H2 gives about `2^70` H-primes in this range against `2 n_b ~ 2^66.7` needed.

**Elements.** baby `(a, w)`: scalar `x = a w`, point `X = [x]P0`. giant `(a', w')`: scalar `x' = a' w'`, point
`X' = [x']A_d = [x' z]P0`.

*Lemma 2.* (i) Every `x`, `x'` lies in `H`. (ii) `x` determines `(a, w)` (`a` is its only prime factor above
`2^77`), so the `N` baby scalars are pairwise distinct, likewise the giants. (iii) For baby/giant pairs
`(x,x')`, `(y,y')` with `x = a w`, `x' = a' w'`, `y = b u`, `y' = b' u'`: if `x/x' = y/y'` as rationals then
either the pairs coincide or `a = b`, `a' = b'`, `w = w'`, `u = u'` (a "diagonal" pair, whose ratio is
`a/a'` for every `w`). The fraction of diagonal pairs among all pairs is `1/L = 2^-50`.

Proof. (i) is closure of `H` under products; (ii) is unique factorisation. (iii) `x/x' = y/y'` means
`a w b' u' = b u a' w'`. The prime factors above `2^77` are `{a, b'}` on the left and `{b, a'}` on the right;
the baby and giant bases are disjoint, so `a = b` and `b' = a'`. Then `w u' = u w'` gives `{w, u'} = {u, w'}`:
either `w = u` and `u' = w'` (same pair) or `w = w'` and `u = u'`. QED

## 4. The walk, the key and the match

**Chain.** For a base prime (baby `Q = [a]P0`; giant `Q = [a']A_d`) build the table `tau[g] = [g]Q` for the even
`g <= G` (`tau[2]` by one doubling, `tau[4] = 2 tau[2]`, then `tau[g] = tau[g-2] + tau[2]`), and walk
`X_{j+1} = X_j + tau[g_j]` with `X_j = [w_j]Q` (baby: `[a w_j]P0`; giant: `[a' w_j]A_d`). To batch inversions on
one node, the walk of a chain is cut into `2^16` segments of equal length, each started by a scalar
multiplication `[w_s]Q` (and the table into `2^16` blocks, each started by a scalar multiplication of `Q`);
the `2^16` segments advance in lockstep and share one Montgomery inversion per round. This is the batched
incomplete affine reference operation (6M with S = M, A = 0; the 3M batch amortisation is inside the reference
price).

*Lemma 3 (no degenerate addition).* Every addition adds two affine points with different x-coordinates, neither
being infinity.

Proof. In the walk, `X_j = [w_j]Q` and `tau[g] = [g]Q` with `0 < g <= G = 2^33 < 2^39 < w_0 <= w_j < r/2`, so
`w_j != +-g mod r`; `Q` has order `r`, hence `X_j != +-tau[g]`. In the table, `tau[g-2] + tau[2]` with `g >= 6`
has `g-2 != +-2 mod r`; a block start `[g_0]Q` with `g_0 >= 6` plus `tau[2]` likewise. QED

**Key.** For `Y = (x, y)` let `y*` be the even one of `y`, `p - y` (`p` is odd, exactly one is even) and let
`kappa(Y)` be bits 1 to 224 of `y*` (bit 0 is always 0).

*Lemma 4.* Two points in one orbit of `{+-1, +-phi, +-phi^2}` have equal `kappa`. If two points have equal `y*`
they are in one orbit. Proof. The orbit points are `(beta^i x, +-y)` and `y*` is defined by the unordered pair
`{y, p-y}`. If `y' = +-y` then `x'^3 = y'^2 - 4 = x^3`, so `x' = beta^i x`. QED
For points in different orbits, `kappa` agrees only by chance (hypothesis **H3**, Section 6).

**Match.** If a baby `X = [x]P0` and a giant `X' = [x' z]P0` are in one orbit, then `x = u x' z mod r` for some
`u in mu6`, so `z = x (u x')^{-1}`. The attack tries the six `u` and accepts the one with `[z]P0 = A_d` (one
scalar multiplication, at most 383 reference operations; at most 6 per verification, `<= 2^12` units).
Every orbit match yields the true `z`; a match of equal `kappa` between different orbits fails this test.

## 5. Machine model, records and matching

**Machine.** A shared-nothing cluster of `S <= 2^85` identical nodes. Each node is a scalar AArch64 core with at
most `2^40` bytes of local RAM, addressed by ordinary 64-bit loads and stores; **no instruction addresses another
node's memory**. Node identifiers are 85-bit strings read as four digits of 22, 22, 22 and 19 bits. Data moves between
nodes only in blocks of 1024 records (48 KiB), over an interconnect in which a node can send a block to any node
whose identifier differs from its own in one digit. The interconnect hardware moves the bytes; the software work
for each block is charged: at most 60 instructions to post a send and 40 to poll its completion, i.e. under 0.1
instruction per record per hop, plus the record-by-record copying counted below.

**Gap and base streams.** `W` (the gaps, 5 bytes each, `2^52.3` bytes) and the base-prime lists are stored on dedicated nodes. A
generating node walks one chain as `2^16` segments; each segment needs a contiguous piece of the gap list (`2^34` gaps), which
reaches the node in blocks of 1 MiB (posting and polling at most 100 instructions per block, under 0.001 per step). The
per-step copy into the local stream buffer and the load of the gap are part of the 40-instruction walk bookkeeping (about 3
instructions). A label `(base index, step)` is resolved to `x = a w` by two block reads from these nodes (under 200
instructions per resolution; resolutions happen only inside the verification cap, row `2^32`).

**Record.** 48 bytes: `kappa` (28 B), label (66-bit base index, 50-bit step, type bit: 16 B), 4 B padding.

**Scatter pass.** One pass over a stream of records: per record `ubfx` of the digit (1), `ldr` of the bucket's
write pointer from a table of `2^22` entries (1), copy of 6 words by 3 `ldp` and 3 `stp` (6), increment and store
of the bucket fill count (3), block-full test (2), loop control (3): 16 instructions, plus under 0.1 for block
posting. The *local LSD pass* adds a counting loop (load digit, load/increment/store counter, loop control: 8)
and has the same scatter; both are charged **23** instructions per record per pass.

**Matching.** After generation, every record is routed by four scatter hops on the first 85 bits of `kappa`
(digits 22, 22, 22, 19) to its owner node, which then holds about `2^31.7` records, all with the same 85-bit
prefix. The owner sorts them with three LSD passes on `kappa[86..154]` (digits of 23 bits) and scans once. Records
with equal 154-bit prefix are adjacent (a *run*); within each run the scan compares full `kappa` of all pairs of
opposite type and verifies each equal-`kappa` pair by the Match step of Section 4. Total per record: 4 + 3 = 7
passes x 23 = 161, scan 10.

**Caps.** The attack stops with *fail* after more than `2^90` run-pair comparisons or more than `2^20`
verifications. Their costs are charged at the caps (`2^90 * 20 instr / 2000 < 2^84` and
`2^20 * 6 * 383 < 2^32` units). All other work is data independent: the number of records, scatter passes and sort passes
is fixed.

**Placement, capacity and overflow.** There are `2N = 2^116.66` records and `S = 2^85` owner nodes, so the mean
occupancy is `mu = 2N/S = 2^31.66` records per node. A node can hold `C = 2^40 / (2 * 48) = 2^33.42` records
(two buffers of 48-byte records for the LSD passes), i.e. `C/mu = 3.37` times the mean. Each node also keeps, per routing hop, a
table of `2^22` write pointers and one partially filled 48 KiB block per bucket (`2^22 * 48 KiB = 2^37.6` bytes), which with the
records stays below `2^40` bytes; a generating node additionally holds one chain table of `2^32 * 96 B = 2^38.6` bytes
(records of a generating node are streamed out in blocks and not retained). A node that is about to exceed `C` records
signals *overflow* and the whole attack stops with *fail* (no spill mechanism is used and none is needed under H3
below). The overflow event is the only storage failure: the buffers of the routing hops are bounded by construction and
a block is sent when full.

## 6. Algorithm, success probability and the caps

```text
Algorithm A (deterministic).
 1. Build W and the base-prime lists; read A_d and T.
 2. Generate all N baby and all N giant records (waves of nodes, each node walks one chain as 2^16 lockstep
    segments); each record enters the routing hops as it is produced.
 3. Matching as in Section 5, with the two caps. On the first verified z run Lemma 1 and output s.
```

Take `K = 0.52` and `N = n_b L` with `n_b = ceil(sqrt(K m / 6) / L)`: `N = 2^115.6644`, `6 N^2/m = 0.52`
(`evidence.py ledger` recomputes this with integers).

Without caps, A succeeds exactly for `z in S = { x/(u x') : baby x, giant x', u in mu6 }`, a union of `6 N^2`
values of `H`, so `P(success) = |S|/m` for uniform `z`.

**H1 (coverage).** The `6 N^2` values behave like independent uniform elements of `H` (Poisson coverage):
`|S|/m = 1 - exp(-6 N^2/m) = 1 - exp(-0.52) = 0.4055`. Support. By Lemma 2(iii) the values are pairwise distinct as
rationals except for a `2^-50` fraction of diagonal pairs, so there is no systematic coincidence; within a row or a
column of the `N x N` array the values are distinct (a negative correlation that can only raise coverage); the
remaining coincidences are modular events `x1 x2' = x2 x1' mod r` among integers below `2^174`, products reduced
modulo a 255-bit prime, for which no algebraic structure is known.

**H3 (keys, random-function form).** The map from the `2N` scalars `x` of the attack to the key bits `kappa(+-[x]P0 or
+-[x]A_d)` behaves like a uniformly random function: the 154-bit prefixes of the records are indistinguishable from
independent uniform strings (a heuristic about coordinates of structured elliptic-curve points, in the same spirit as H1;
the keys of the baby records are fixed and not randomised by the challenge, so this is a statement about those fixed
keys). Two consequences are used: (a) pairs of records from different orbits agree on a prescribed `b` bits of `kappa`
with probability `2^-b`; (b) the numbers of records per owner node are binomial with mean `mu`, hence Chernoff-concentrated.

**Caps, overflow and success.** Under H3(a), with `2N` records, the expected number of record pairs with equal 154-bit
prefix is `C(2N,2) 2^-154 < 2^78.4`, so by Markov the cap `2^90` is exceeded with probability `< 2^-11.6`; the expected number
of pairs with equal full `kappa` from different orbits is `C(2N,2) 2^-224 < 2^8.4`, so the cap `2^20` on
verifications is exceeded with probability `< 2^-11.6` (true orbit matches stop the run at the first verified one). Under H3(b)
the Chernoff bound for a node with `C = a mu`, `a = 3.37`, is `exp(-mu (a ln a - a + 1)) = 2^-(8.4 * 10^9)`, and the union
bound over `S = 2^85` nodes is `2^(85 - 8.4 * 10^9)` (`evidence.py ledger` prints the exponent); we budget `2^-20` for the
overflow event to absorb modelling error. Hence `P(success) >= 0.4055 - 2^-10.6 - 2^-20 > 0.4048`, and **0.40** is claimed,
`0.01` above the requirement.

**Toy check** (`evidence.py coverage`, scalar level, same prime structure; three toy primes `r` with `6 * 2^k | r-1`,
chain length `L = 60`, `K ~ 0.55`, 40k trials each with seed 2026 run locally as `{"mode":"coverage","toy":"t28",
"trials":40000,"seed":2026,"K_permille":520}`; the manifest's run uses fewer). At toy scale the diagonal pairs are a
visible `1/L` fraction, so the Poisson prediction is applied with `K' = K (1 - (L-1)/L^2)`:

| toy `r` | measured success | Poisson(`K`) | Poisson(`K'`) | z vs `K'` |
| --- | ---: | ---: | ---: | ---: |
| `2^27.6` (d=32) | 0.4243 | 0.4302 | 0.4249 | -0.25 |
| `2^29.7` (d=64) | 0.4209 | 0.4245 | 0.4192 | +0.67 |
| `2^31.7` (d=128) | 0.4140 | 0.4187 | 0.4135 | +0.19 |

The corrected model fits within one standard error; at production size `1/L = 2^-50`. This is evidence at toy
scale only: it supports the Poisson form and the bookkeeping of duplicates, not the production value of `m`.

## 7. Cost ledger (unit = 2000 instructions)

Per record (baby or giant), beyond the reference addition, in AArch64 instructions:

| item | instructions |
| --- | ---: |
| walk bookkeeping (table address from the gap, segment counter, label, batch-pointer updates; the multiplications of the batched inversion are inside the 6M) | 40 |
| key `kappa` (parity test of `y`, `p - y` over 6 limbs, select, extract 224 bits, write the record) | 30 |
| four routing scatter hops, 23 each (Section 5) | 92 |
| three local LSD passes, 23 each | 69 |
| run scan | 10 |
| **itemised** | **241** |
| **budget used (about 24% slack)** | **300** |

The element cost is `c = 1 + 300/2000 = 1.15` units.

| term | bound (units) |
| --- | ---: |
| baby + giant elements `2 N c` | `2^116.8660` |
| tables: `2 n_b` chains x `2^32` entries x (1 addition + 12 stores/2000) | `2^98.67` |
| segment and table-block starts, `[a]P0` (<= `2^25` per chain) | `2^91.66` |
| per-node fixed work (histograms of `2^23` counters, bucket tables, setup): `S x 2^28` instr | `2^102.03` |
| H-primes in `(2^77, 2^97)`: about `2^86.7` primes at `353 x 383 / 2000 = 67.6` units per modexp, `2^88.4` Miller-Rabin candidates at 27.8, sieve; 4x margin | `2^97` |
| H-primes of `W` below `2^77`: about `2^70` primes tested, same pricing, 4x margin | `2^84` |
| locate/decode `A_d`, `T` in the input | `2^20` |
| Lemma 1 | `2^21` |
| false-match cap `2^20` verifications | `2^32` |
| run-pair comparison cap `2^90` | `2^84` |

All F_r multiplications in the prime search are priced at the 383 instructions measured for `blst_fp_mul`
(381-bit), an upper bound for a 255-bit product. All rows other than `2Nc` sum to `2^102.2087`, and
`W = 2^116.8660 + 2^102.2087 = 2^116.86604 <= 2^116.8661` (`evidence.py ledger` recomputes every row with exact
integers). The score is the exact rational `1168661/10000 = 116.8661`.

*Preprocessing* (target independent: lists, tables, baby generation): `<= N c + 2^102.21 = 2^115.8661 <= 2^115.87`.
*Advice*: none; every searched object is rebuilt and charged. *Memory*: records `2 N x 48 B x 2` (two buffers for the
LSD passes) `= 2^123.2493`; tables of `2^20` active chain nodes `2^20 * 2^32 * 96 B = 2^58.6`; `W`
`2^50 * 5 B = 2^52.3`; base-prime lists `2 n_b * 12 B = 2^70`; total `<= 2^123.3` bytes. With `S = 2^85` nodes this is
about `2^38.3` bytes per node. Memory is exponential, as for any BSGS-type Cheon attack.

## 8. Heuristics and limitations

* **H1** (coverage) decides the success probability; **H3** (key bits, random-function form) decides the failure probability of
  the two caps and of node overflow. If H3(b) failed, some nodes would exceed `C` and the attack would fail more often;
  H3 is supported only by the toy occupancy test of Section 9 and by the absence of any known structure in these coordinates.
* **H2/H2'** (H-prime density, gap bound, `w_{L-1} < 2^77`) decide the preprocessing cost and the table size; the
  computation of `W` checks them, and the table row is `2^98.67` of `2^116.87`.
* The cluster model is idealised: no instruction accesses remote memory; block transfers are charged by posting and
  polling work and by the record copies of each pass. The instruction counts are itemised, not measured on hardware.
* The cost bound is data independent except for the two caps, which are charged at their limits.
* Full-size runs are impossible; the evidence is analytic plus toy experiments.

## 9. Evidence (`evidence.py`)

The manifest declares one initial run (`mode: all`, at most one is allowed): `ledger`, `coverage` and `ec` with small
fixed trial counts (about 5 s); each can be rerun on its own with more trials as a follow-up.

* `ledger`: exact integer recomputation of Sections 6-7 (including the cap probabilities and the memory sum), the facts
  `d | r-1`, `6d | r-1`, `mu6 subset H`, `v_2(r-1) = 32`, the H2 sanity check and the rounded score.
* `coverage`: scalar-level toy of the semiprime baby/giant sets with the six-fold orbit table.
* `ec`: end-to-end toy of Algorithm A on genuine `j = 0` curves of prime order `r ~ 2^28..2^32`: lockstep walks of affine
  incomplete additions with Montgomery batch inversion (asserting no zero denominator), the parity key **truncated** to
  18/20/22 bits so that false key matches occur (about 0.5-1.8 per run) and are verified and rejected under a cap, sort,
  scan of equal-key runs, candidate `z`, `d`-th root and `mu_d` scan. Every recovered `z` and `s` is checked and success
  agrees with the scalar model (`z in S`) in every trial. It also bins the truncated keys by their top 5 bits (32 "nodes") and reports the
  chi-square z-score of the bin counts of the (challenge independent) baby records and, summed over trials, of the giants, with the
  largest bin over the mean against the Poisson median-maximum (locally: z in `[-0.9, +1.0]` for all three toys and maxima `1.3-1.6`
  versus Poisson `1.23-1.33`), i.e. no sign of non-uniform node occupancy at toy scale. Limits: the toy builds tables only for the gaps that occur, so
  `tau[g] = tau[g-2] + tau[2]` and the segment structure are covered by Lemma 3 only.

## 10. Provenance

No earlier accepted result is used. The idea of walking small integers of `H` with one addition per step, the GLV
six-fold key and Cheon's `d | r-1` reduction are the starting point of the published note of submission `cbadde77`; the
semiprime H-prime construction, the distinctness lemma, the non-degeneracy argument by chain start above the gap bound,
the sort-based matching with caps, the machine model and the itemised overhead were derived here and are checked
independently by this folder's experiments. A previous version of this package (identical algorithm, hash-table matching)
was rejected for a memory figure rounded below its computed value and for an unbounded false-match budget; a second version
(sort-based matching with caps) was rejected because the per-node occupancy had no maximum or overflow rule. Sections 5-6 now
state the capacity, the overflow event and the concentration assumption H3, and Section 9 tests node occupancy at toy scale.
