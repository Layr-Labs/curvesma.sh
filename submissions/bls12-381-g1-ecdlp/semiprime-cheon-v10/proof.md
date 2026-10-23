# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks and sort-based matching

Exploratory claim: a classical deterministic attack with a total-work bound at the scoring value `n = 2^20`, success
probability at least 0.395 under disclosed heuristics, and no advice. Not a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.8271` (exact `1168271/10000`); success `>= 0.395` (H1 predicts `0.3995`; caps and node overflow cost
at most `0.0007`; required `0.39`); peak memory `<= 2^125.1` bytes (every node counted at its full `2^40` bytes);
preprocessing `<= 2^115.83` units (included); no advice. Unit = one batched incomplete affine G1 addition = 2000 AArch64
instructions (`C = 2000`). Only `T` and `A_d`, `d = 2^20`, are read, so the bound holds for every `n >= 2^20`.

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
**Base primes:** H-primes in `(2^77, 2^97)` (about `2^70` exist by H2; `2 n_b ~ 2^66.7` are needed) dealt alternately into a baby list
`a_i` and a disjoint giant list `a'_i`, `n_b = ceil(N/L)` each.
**Elements:** baby `(a,w)`: `x = a w`, point `[x]P0`; giant `(a',w')`: `x' = a' w'`, point `[x']A_d = [x' z]P0`.

*Lemma 2.* (i) All `x, x'` lie in `H`. (ii) `x` determines `(a, w)` (`a` is the only prime factor above `2^77`), so the `N` scalars of a side
are distinct. (iii) For pairs `(x,x')`, `(y,y')`, `x = aw, x' = a'w', y = bu, y' = b'u'`: `x/x' = y/y'` as rationals implies the pairs
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

**Machine.** A shared-nothing cluster: `S = 2^85` owner nodes, `2^20` generating nodes, `2^31` storage nodes (for `W` and the base
lists). A node is a scalar AArch64 core with `2^40` bytes of local RAM, addressed by ordinary 64-bit loads and stores; **no
instruction addresses another node's memory**. Owner identifiers are 85-bit strings read as digits of 22, 22, 22, 19 bits. Nodes
exchange blocks of 64 records (3 KiB) with any node differing from them in one digit; hardware moves the bytes, software charges
at most 60 instructions to post and 40 to poll a block (1.6 per record per hop) plus the copies counted below. **Transport is priced as
if the CPU copied the bytes:** one extra scalar copy of the record at every hop (3 `ldp` + 3 `stp` = 6 instructions per 48-byte record per hop), on top
of the sender-side copy into the bucket. `W` (5 bytes per gap, `2^52.3` bytes) and the base lists stream to a generating node in 1 MiB blocks:
posting and polling under 0.001 instruction per step plus 1 instruction per 8 bytes of transport, i.e. 1 per step; a label
`(base index, step)` is resolved by two block reads (under 200 instructions, only inside the verification cap).

**Record**: 48 bytes: `kappa` 28 B, label 16 B (66-bit base index, 50-bit step, type bit), 4 B padding.
**Routing pass** (per record): `extr`+`ubfx` of the digit (2), `ldr` of the bucket pointer (1), 3 `ldp` + 3 `stp` (6), fill count
`ldr/add/str` (3), block-full test (2), loop (3) = 17, plus 1.6 posting and 6 for the transport copy = 24.6: charged **25**.
**Local LSD pass**: counting loop (`ldr`, `extr`, `ubfx`, `ldr/add/str`, loop: 9) + scatter (`extr`, `ubfx`, `ldr` offset, `madd`, `add`,
6 copies, `add/str`, loop: 16) = 25: charged **25**.

**Matching.** Each record is routed by 4 routing passes on the first 85 bits of `kappa` to its owner (about `2^31.7` records, one 85-bit
prefix), sorted by 3 LSD passes on `kappa[86..154]` (23-bit digits) and scanned once. Records with equal 154-bit prefix are adjacent
(a run); inside a run all opposite-type pairs with equal `kappa` are verified by the Match step. **Caps:** after more than `2^90`
run-pair comparisons or `2^20` verifications the attack stops with *fail*; they are charged at the caps (`2^84` and `2^32` units).
All other work is fixed: the number of records and passes does not depend on the data.

**Placement and overflow.** `2N = 2^116.66` records, `mu = 2N/S = 2^31.65` per owner. Per-node allocation (bytes): two record buffers
`C * 96 <= 2^40 - 2^36`; routing blocks `2^22 * 3 KiB = 2^33.6`; pointer tables `4 * 2^22 * 8 = 2^29`; receive pool `2^20 * 3 KiB = 2^31.6`;
LSD histograms `2 * 2^23 * 4 = 2^26`; code, stack `<= 2^30`; auxiliary total `2^34.0`, reserved `2^36`, concurrently, no phased reuse.
So `C = (2^40 - 2^36)/96 = 2^33.32` records, `C/mu = 3.19`. A node that would exceed `C` signals *overflow* and the attack fails.
A generating node holds a chain table `2^32 * 96 B = 2^38.6`, routing blocks `2^33.6` and stream buffers `<= 2^30` (`<= 2^38.7`).

## 6. Algorithm and success probability

```text
1. Build W and the base lists; read A_d, T.
2. Generate all N baby and N giant records (each generating node walks one chain as 2^16 segments); route each record as produced.
3. Match as in Section 5 with the caps; on the first verified z run Lemma 1 and output s.
```

`K = 0.51`, `N = n_b L` with `n_b = ceil(sqrt(K m/6)/L)`: `N = 2^115.6503`, `6N^2/m = 0.51`. Without caps and overflow the attack succeeds
exactly on `S = {x/(u x')}`, a union of `6N^2` values of `H`, so `P(success) = |S|/m`.
**H1:** these values behave like independent uniform elements of `H`: `|S|/m = 1 - e^{-0.51} = 0.3995`. Support: by Lemma 2(iii) they are
distinct as rationals except a `2^-50` diagonal; within a row/column they are distinct (negative correlation); the other coincidences are
modular events among integers below `2^174` reduced modulo a 255-bit prime with no known structure.
**H3 (random-function form):** the 154-bit key prefixes of the `2N` records behave like independent uniform strings (a heuristic about
coordinates of structured points, as H1; the baby keys are fixed, not randomised by the challenge). Used: (a) records of different orbits
agree on `b` prescribed bits with probability `2^-b`; (b) node occupancies are binomial with mean `mu`.
**Caps, overflow, success.** By (a) the expected number of record pairs with equal 154-bit prefix is `< 2^78.4`, so the `2^90` cap is exceeded
with probability `< 2^-11.6` (Markov); equal full `kappa` from different orbits: `< 2^8.4` pairs, so the `2^20` cap fails with probability `< 2^-11.6`.
By (b), Chernoff for `C = 3.19 mu` gives `2^-(7.3 * 10^9)` per node, union over `2^85` nodes (`evidence.py` prints the exponent); `2^-20` is
budgeted for modelling error. So `P(success) >= 0.3995 - 2^-10.6 - 2^-20 > 0.3988`; **0.395** is claimed, `0.005` above the requirement.

**Toy check** (`evidence.py coverage`, scalar level, three toy primes with `6 * 2^k | r-1`, chain length `L = 60`, `K ~ 0.55`, seed 2026, 40k trials run locally;
the manifest run uses fewer). At toy scale the diagonal pairs are a `1/L` fraction, so the Poisson prediction uses `K' = K(1 - (L-1)/L^2)`. Measured
success for `r ~ 2^27.6, 2^29.7, 2^31.7`: `0.4243, 0.4209, 0.4140` against Poisson(`K'`) `0.4249, 0.4192, 0.4135` (|z| < 0.7; Poisson(`K`) is 0.4302,
0.4245, 0.4187). Evidence at toy scale only: it supports the Poisson form and the duplicate bookkeeping, not the production value of `m`.

## 7. Cost ledger (unit = 2000 instructions)

| per record (instructions beyond the reference addition) | |
| --- | ---: |
| walk bookkeeping (listed 28: gap load and copy 3, gap-stream transport 1, table address `base + (gap/2)*96` 3, segment loop 3, label 3, batch pointer/index updates outside the Montgomery-trick multiplications 12, stream-refill check 3) | 28 |
| key (listed 23: `ldp` of 6 `y` limbs 3, parity test and select mask 2, `p - y` by `subs`, 5 `sbcs`, 6 `csel` 12, `extr` to bits 1..224 4, `stp` of the key 2) | 23 |
| 4 routing passes x 25 (17 + 1.6 posting + 6 transport) | 100 |
| 3 local LSD passes x 25 | 75 |
| run scan (load the 69-bit sort prefix, compare with previous, branch) | 10 |
| **listed 236; charged (about 10% slack for unlisted loop and spill instructions)** | **260** |

Element cost `c = 1 + 260/2000 = 1.13`.

| term | units |
| --- | ---: |
| elements `2Nc` | `2^116.8267` |
| tables: `2 n_b` chains x `2^32` entries x (1 addition + 12 stores/2000) | `2^98.67` |
| segment and table-block starts, `[a]P0` (`<= 2^25` per chain) | `2^91.66` |
| per-node fixed work, `S x 2^31` instr (histograms, bucket and pointer tables, receive pool, flush of up to `4 x 2^22` partial blocks at 100 instr, setup) | `2^105.03` |
| H-prime search in `(2^77,2^97)`: `2^86.7` primes x `353 * 383/2000 = 67.6`, `2^88.4` Miller-Rabin candidates x 27.8, sieve, 4x margin | `2^97` |
| H-primes of `W` below `2^77`: `2^70` primes, same pricing, 4x margin | `2^84` |
| input access `2^20`; Lemma 1 `2^21`; verification cap `2^32`; run-pair cap `2^84` | `<= 2^84` each |

All F_r multiplications in prime searches are priced at the 383 instructions measured for `blst_fp_mul` (381-bit), an upper bound. The rows
other than `2Nc` sum to `2^105.0571`; `W = 2^116.8267 + 2^105.0571 = 2^116.82708 <= 2^116.8271`, score `1168271/10000` (`evidence.py ledger`
recomputes every row with exact integers). Preprocessing (lists, tables, baby generation) `<= Nc + 2^105.06 = 2^115.8275 <= 2^115.83`.
Memory counted at full allocation: `(2^85 + 2^20 + 2^31) * 2^40 < 2^125.0001`, claimed `<= 2^125.1` (records alone: `2N * 48 B * 2 = 2^123.24`;
`W` `2^52.3` and base lists `2^70.2` bytes fit the `2^31` storage nodes). Memory has no scalar weight and is exponential, as for any BSGS-type Cheon attack.

## 8. Heuristics, evidence, provenance

H1 decides success; H3 decides the failure probability of the caps and of node overflow (supported only by the toy occupancy test and the absence of
known structure); H2/H2' decide preprocessing and table size (checked by building `W`). The cluster model is idealised and the instruction counts are
listed, not measured on hardware. The cost is data independent except for the two caps, charged at their limits.

`evidence.py` (one initial run `mode: all`, about 5 s; each mode can be rerun with more trials): `ledger` (exact integers: Sections 6-7, cap
probabilities, capacity, memory, `6d | r-1`, `mu6 subset H`, H2 check); `coverage` (Section 6); `ec`: Algorithm A end to end on genuine `j = 0` curves
of prime order `r ~ 2^28..2^32`: affine incomplete additions with Montgomery batch inversion (asserting no zero denominator), the parity key truncated to
18/20/22 bits so that false key matches occur (0.5-1.8 per run) and are verified and rejected under a cap, sort, equal-key runs, `z`, `d`-th root;
every recovered `z`, `s` is checked, success agrees with the scalar model in every trial, and keys binned into 32 pseudo-nodes show chi-square z-scores
in `[-0.8, +1.0]`. Limits: tables are built only for occurring gaps (so `tau[g] = tau[g-2] + tau[2]` and the segments rest on Lemma 3), and the toy
sorts locally without distributed routing.

Provenance: the idea of walking small integers of `H` with one addition per step, the GLV six-fold key and Cheon's `d | r-1` reduction come from the
published note of submission `cbadde77`; the semiprime construction, Lemmas 2-3, matching, caps, machine model and overhead listing were derived here.
