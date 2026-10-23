# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks (maximum over challenges)

Exploratory claim under `rules.md`: a classical deterministic attack with a stated maximum-over-challenges
work bound `W` (under H2 and expected-`O(1)`-probe hashing) at the scoring value `n = 2^20`, success probability at least 0.40 under the disclosed
heuristic H1, and no advice. Nothing here is a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.8665` (exact rational `233733/2000`), success `>= 0.40` (H1 predicts
`0.4055`, the rules require `0.39`), peak memory `<= 2^121.4` bytes, preprocessing `<= 2^115.87` units
(included in the total), no advice. Unit = one batched incomplete affine G1 addition = 2000 AArch64
instructions (`C = 2000`).

## 1. Setting and notation

* `r` is the BLS12-381 scalar prime, `p` the base-field prime, `G1 = <P0>` has order `r`.
* `r - 1 = 2^32 * 3 * 11 * 19 * 10177 * ...`. Put `d = 2^20` and `m = (r-1)/d` (`log2 m = 234.857`).
  `H = {w^d : w in F_r^*} = {x : x^m = 1}` is the subgroup of order `m`.
* `mu6 = {u : u^6 = 1} = {+-1, +-omega, +-omega^2}` with `omega` a primitive cube root of unity mod `r`.
  Because `6d | r-1`, `mu6 subset H` (checked by `evidence.py ledger`).
* The GLV map `phi(x,y) = (beta x, y)` acts on `G1` as `[omega]` for one primitive cube root `omega`.
  Together with negation, the orbit of a point `[v]P0` under `{+-1, +-phi, +-phi^2}` is
  `{[u v]P0 : u in mu6}`.
* The input supplies `T = A_1 = [s]P0` and `A_i = [s^i]P0`, `i <= n`. The attack reads only `T` and
  `A_d` with `d = 2^20`, so it applies for every `n >= 2^20` and `W` does not depend on `n`.

## 2. Reduction (Cheon, variant d | r-1)

For `s` uniform in `{1..r-1}`, `z := s^d` is uniform in `H` (the map `s -> s^d` is `d`-to-1 onto `H`) and
`A_d = [z]P0`. The remainder of the attack searches for `z`.

*Lemma 1.* Given `z in H`, the secret `s` is recovered with at most `2^23` units.

Proof. `s` is a `d`-th root of `z`. Let `z_0 = z` and let `z_{i+1}` be any square root of `z_i`.
Writing `z = w^{2^20}`, each `z_i` (i <= 19) has the form `eps_i * w^{2^(20-i)}` with `eps_i` of 2-power
order at most `2^i`; both factors are squares (the root of unity because its order divides `2^19`,
which divides `(r-1)/2`; the second factor because `20 - i >= 1`), so all 20 square roots exist, and
`s_0 := z_20` satisfies `s_0^d = z`. All `d`-th roots are `s_0 * omega_d^j`, `j < 2^20`, for a primitive
`d`-th root of unity `omega_d`. Exactly one satisfies `[s_j]P0 = T` since `P0` has order `r`. Finding it:
let `R = [s_0]P0` (one scalar multiplication). With `j = j_1 + 2^10 j_2` we need `[omega_d^{j_1}]R =
[omega_d^{-2^10 j_2}]T`; compute both sides for all `j_1, j_2 < 2^10` by `2^11` scalar multiplications
(each at most 383 reference additions by double-and-add, as scalars are in `F_r`) and match them in a
table. Cost `<= 2^11 * 383 < 2^20` units. Twenty Tonelli-Shanks square roots cost `< 2^21`
instructions. QED

## 3. Baby and giant sets

**H-primes.** A prime `q` is an *H-prime* if `q^m = 1 mod r`. Heuristic **H2**: for primes `q << r`
the classes `q mod r` equidistribute in `F_r^*/H` (index `d`), so H-primes have density `1/(d ln q)`
among primes. Sanity check on this `r` (`evidence.py ledger`): of the 6542 primes below `2^16`,
402 (0.0614) are 16th powers (expected 0.0625) and 27 (0.0041) are 256th powers (expected 0.0039).
H2 only affects the cost of *finding* H-primes; Section 6 charges that with a 4x margin.

**Walk list `W`.** The first `L = 2^50` H-primes above `2^39`: `w_0 < w_1 < ... < w_{L-1}`. By H2,
`w_{L-1} ~ 2^50 * 2^20 * ln(2^76) ~ 2^75.7 < 2^76`. The gaps `g_j = w_{j+1} - w_j` are even. Put
`G = 2^38`. A gap above `G` (about `exp(-2^12)` per gap) is cut into chunks of length at most `G`; this adds
at most `2^76/G = 2^38` additions per chain, a deterministic bound.

**Base primes.** H-primes in `(2^76, 2^96)` taken in increasing order and dealt alternately into a baby
list `a_1, a_2, ...` and a disjoint giant list `a'_1, a'_2, ...`, `n_b` of each with `n_b = ceil(N/L)` and the set size `N` as in Section 5 (`N` is not the input point `T`, and
`n_b` is not the power count `n`). H2 gives about `2^69.9` H-primes in this range against `2 n_b ~ 2^66.7` needed.

**Elements.**
* baby `(a, w)`: scalar `x = a w`, point `X = [x]P0`;
* giant `(a', w')`: scalar `x' = a' w'`, point `X' = [x']A_d = [x' z]P0`.

*Lemma 2.* (i) Every `x`, `x'` lies in `H`. (ii) `x` determines `(a, w)` (`a` is its only prime factor
above `2^76`), so the `N` baby scalars are pairwise distinct, and likewise the giants. (iii) Let
`(x, x')` and `(y, y')` be two baby/giant pairs with `x = a w`, `x' = a' w'`, `y = b u`, `y' = b' u'`.
If `x/x' = y/y'` as rationals, then `a = b`, `a' = b'` and `w = w'`, `u = u'` (the "diagonal" pairs), or the
pairs coincide. The fraction of diagonal pairs among all pairs is `1/L = 2^-50`.

Proof. (i) is closure of `H` under products; (ii) is unique factorisation. (iii) `x/x' = y/y'` means
`a w b' u' = b u a' w'`. The prime factors above `2^76` are `{a, b'}` on the left and `{b, a'}` on the
right; baby and giant bases are disjoint, so `a = b` and `b' = a'`. The remaining equality
`w u' = u w'` of products of two primes of `W` gives `{w, u'} = {u, w'}` as multisets: either `w = u,
u' = w'` (same pair) or `w = w'` and `u = u'`. QED

## 4. The walk, the key and the match

**Chain.** For a base prime (baby: `Q = [a]P0`; giant: `Q = [a']A_d`) the attack builds
`tau[g] = [g]Q` for every even `g <= G` (`tau[2]` is one doubling, `tau[4] = 2 tau[2]`, then
`tau[g] = tau[g-2] + tau[2]`), computes `X_0 = [w_0]Q` and then `X_{j+1} = X_j + tau[g_j]`. Hence
`X_j = [w_j]Q`, i.e. `[a w_j]P0` for a baby and `[a' w_j]A_d` for a giant.

*Lemma 3 (no degenerate addition).* Every addition above adds two affine points with different
x-coordinates, neither being infinity.

Proof. `X_j = [w_j]Q` and `tau[g] = [g]Q` with `0 < g <= G = 2^38 < 2^39 < w_0 <= w_j < r/2`, so
`w_j != +-g mod r`, and `Q` has order `r`, hence `X_j != +-tau[g]`. A chunked gap adds chunks `<= G`
to points `[w_j + t]Q` with `w_j + t > G`: same argument. In the table, `g - 2 != +-2 mod r` for
`g >= 6`. QED

So every step is one incomplete affine addition.

**Batching.** `2^16` chains of one wave run in lockstep. One round is `2^16` independent additions sharing
one Montgomery inversion, the batched incomplete affine reference operation (6M with S = M, A = 0; the
3M batch amortisation is inside the reference price). All chains of a wave use the same gap `g_j` in
round `j`.

**Key.** For a point `Y = (x, y)` let `y*` be the even one of `y`, `p - y` (`p` is odd, so exactly one is
even) and let `kappa(Y)` be bits 1 to 224 of `y*` (bit 0 is always 0), a 224-bit value.

*Lemma 4.* Two points of `G1` in one orbit of `{+-1, +-phi, +-phi^2}` have equal `kappa`. Conversely
equal `kappa` of two points in different orbits has probability `2^-224` per pair.

Proof. The six orbit points are `(beta^i x, +-y)`, and the parity rule picks the same representative
of `{y, p-y}`. If two points have `y' = +-y` then `x'^3 = y'^2 - 4 = x^3`, so `x' = beta^i x` and they are
in one orbit; otherwise their `y*` values differ, and bits 1 to 224 of two different even representatives of
381-bit field elements agree with probability `2^-224` under the usual uniformity of `y` bits. QED

**Match.** If a baby `X = [x]P0` and a giant `X' = [x' z]P0` are in one orbit, then `x = u x' z mod r`
for some `u in mu6`, i.e. `z = x * (u x')^{-1}`. The attack tries the six `u` and accepts the one with
`[z]P0 = A_d` (one scalar multiplication, at most 383 reference additions). False matches (about
`N^2 2^-224 = 2^7.3` over the whole run) merely fail this test.

## 5. Algorithm and success probability

```text
Algorithm A (deterministic).
 1. Build W and the base-prime lists. Read A_d and T from the input.
 2. Baby phase: for every baby base, in waves of 2^16 chains and rounds j = 0..L-1, compute the
    elements, kappa(X), and insert (kappa, label, type) into the sharded hash table.
 3. Giant phase: the same for giant bases; every element is looked up. On a hit run the Section 4 test.
 4. On the first verified z run Lemma 1 and output s. If all T giants are processed, output fail.
```

Take `K = 0.52` and `N = n_b L` with `n_b = ceil(sqrt(K m / 6) / L)`: `N = 2^115.6644`, `6 N^2/m = 0.52`
(`evidence.py ledger` recomputes this with integers).

Algorithm A succeeds exactly for `z in S = { x / (u x') : baby x, giant x', u in mu6 }`. `S` is the union
of `6 N^2` values of `H`, so `P(success) = |S|/m` for uniform `z`.

**H1 (coverage).** The `6 N^2` values behave like independent uniform elements of `H` (Poisson coverage):
`|S|/m = 1 - exp(-6 N^2/m) = 1 - exp(-0.52) = 0.4055`.

Support. By Lemma 2(iii) the values are pairwise distinct as rationals except for a `2^-50` fraction of
diagonal pairs, so there is no systematic coincidence; within a row or column of the `N x N` array the values
are distinct (a negative correlation, which can only raise coverage); the remaining coincidences are modular
events `x1 x2' = x2 x1' mod r` among integers below `2^171`, products reduced modulo a 255-bit prime, for
which no algebraic structure is known. We claim `0.40`, which leaves `0.0055` below H1 and `0.01`
above the requirement.

**Toy check** (`evidence.py coverage`, scalar level, same prime structure; three toy primes `r` with
`6 * 2^k | r-1`, chain length `L = 60`, `K ~ 0.55`, 40k trials each with seed 2026, run locally as
`{"mode":"coverage","toy":"t28","trials":40000,"seed":2026,"K_permille":520}`; the manifest's run uses fewer).
At toy scale the diagonal pairs of Lemma 2(iii) are a visible `1/L` fraction and repeat the ratio `a/a'`, so
the Poisson prediction is applied with `K' = K (1 - (L-1)/L^2)`. Measured against that prediction:

| toy `r` | measured success | Poisson(`K`) | Poisson(`K'`) | z vs `K'` | mean stop / `N` measured | Poisson(`K'`) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `2^27.6` (d=32) | 0.4243 | 0.4302 | 0.4249 | -0.25 | 0.8430 | 0.8426 |
| `2^29.7` (d=64) | 0.4209 | 0.4245 | 0.4192 | +0.67 | 0.8453 | 0.8449 |
| `2^31.7` (d=128) | 0.4140 | 0.4187 | 0.4135 | +0.19 | 0.8471 | 0.8473 |

The corrected model fits within one standard error, including the early-stopping statistic; at production
size `1/L = 2^-50`. This is evidence at toy scale only: it supports the Poisson form and the bookkeeping of
duplicates, not the production value of `m`.

## 6. Cost ledger (unit = 2000 instructions)

Per element (baby or giant), beyond the reference addition, in AArch64 instructions:

| item | instructions |
| --- | ---: |
| walk bookkeeping (table address, step counter, label, batch-pointer updates; the multiplications of the batched inversion are inside the 6M) | 40 |
| canonical `kappa` (parity test, `p - y` over 6 limbs, select, store the 28-byte key) | 30 |
| local probe/insert in the owning shard (index, tag compare, expected 1.5 probes in a bucketed table at load 0.9, 6 stores; table initialisation is one store per slot, under one instruction per record) | 61 |
| routing: 6 radix hops of `2^16` buckets, 18 instructions per hop for the 48-byte record (key 28 B, label 16 B, type bit) | 108 |
| **itemised** | **239** |
| **budget used (about 25% slack)** | **300** |

The table is sharded by key over `S <= 2^96` shards (at `2^40` bytes per shard `S ~ 2^82`);
`ceil(82/16) = 6` scatter passes suffice. Records move in blocks of `2^10`, so message posting costs
well under one instruction per record. The element cost is `c = 1 + 300/2000 = 1.15` units.

| term | bound (units) |
| --- | ---: |
| baby + giant elements `2 N c` | `2^116.8660` |
| tables: `2 n_b` chains x `2^37` entries x (1 addition + 12 stores/2000) | `2^103.67` |
| chain starts `[a]P0`, `[w_0]Q` (<= 512 per chain) | `2^75.7` |
| gap-chunk bound `2 n_b * 2^38` | `2^104.66` |
| H-prime search in `(2^76, 2^96)`: about `2^86.7` primes at about 21 units (255-bit modexp), plus sieve/Miller-Rabin, with 4x margin | `2^96` |
| H-primes of `W` below `2^76`: about `2^70.3` primes tested | `2^84` |
| locate/decode `A_d`, `T` in the input | `2^20` |
| false matches (`2^7.3` x 6 x 383), `z` candidates, Lemma 1 | `2^23` |

All rows other than `2Nc` sum to `2^105.2546`, and `W = 2^116.8660 + 2^105.2546 = 2^116.86645
<= 2^116.8665` (the printed ledger of `evidence.py` recomputes every row with exact integers). The score
is the exact rational `233733/2000 = 116.8665`.

*Preprocessing* (target independent: lists, tables, baby phase): `<= N c + 2^105.26 = 2^115.87`.
*Advice*: none; every searched object is rebuilt and charged. *Memory*: baby records `N * 48 B / 0.9 = 2^121.40` (bucketed table at load 0.9);
per-wave tables `2^16 * 2^37 * 96 B = 2^59.6`; the walk list `2^50 * 5 B = 2^52.3`; base-prime lists
`2 n_b * 12 B = 2^70`; total `<= 2^121.4` bytes. Memory is exponential, as for any BSGS-type Cheon attack.

## 7. Heuristics and limitations

* **H1** (Section 5) decides the success probability; `0.40` is claimed.
* **H2** H-prime density decides only the preprocessing cost, which is below `2^-20` of the total even
  with the stated margin.
* **Hashing**: expected `O(1)` probes; the sharded table with 6 routing hops is an idealised cluster
  model charged only through the instruction counts above; a physical network is not modelled.
* Instruction counts are itemised, not measured on hardware.
* Full-size runs are impossible; the evidence is analytic plus toy experiments.
* The bound is the maximum over challenges *under H2 and the expected-`O(1)`-probe hash assumption*;
  Algorithm A also stops at the first verified match.

## 8. Evidence (`evidence.py`)

The manifest declares one initial run (`mode: all`, at most one is allowed), which executes the three modes below
with small fixed trial counts (about 5 s); each mode can be rerun on its own with more trials as a follow-up.

* `ledger`: exact integer recomputation of Sections 5-6, the facts `d | r-1`, `6d | r-1`,
  `mu6 subset H`, `v_2(r-1) = 32`, the H2 sanity check and the rounded score.
* `coverage`: scalar-level toy of the semiprime baby/giant sets with the six-fold orbit table.
* `ec`: end-to-end toy on genuine `j = 0` curves of prime order `r ~ 2^28..2^32`: lockstep waves of affine
  incomplete additions with Montgomery batch inversion (asserting no zero denominator), the parity key,
  invariance of the key under `phi` and `-1`, the unified table of the interleaved schedule (the success set `S` is the
  same as for the two-phase Algorithm A), candidate `z`, `d`-th root and `mu_d` scan;
  every recovered `z` and `s` is checked and the stopping round equals the scalar model's. Limits: the toy
  builds tables only for the gaps that occur and never needs chunking (its chain start exceeds every gap), so
  `tau[g] = tau[g-2] + tau[2]` and chunking are covered by Lemma 3 only.

## 9. Provenance

No earlier accepted result is used. The idea of walking small integers of `H` with one addition per step,
the GLV six-fold key and Cheon's `d | r-1` reduction are the starting point of the published note of
submission `cbadde77`; the semiprime H-prime construction, the distinctness lemma, the non-degeneracy
argument by chain start above the gap bound, the itemised overhead and the parity key were re-derived
here and are checked independently by this folder's experiments.
