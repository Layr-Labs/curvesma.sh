# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks and sort-based matching

Exploratory claim: a classical deterministic attack with a total-work bound at the scoring value `n = 2^20`, success
probability at least 0.39 under disclosed heuristics, and no advice. Not a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.8172` (exact `292043/2500`); success `>= 0.39` (H1's birthday model gives `0.4054795`, the
model allowance is `0.01`, the per-node caps, node overflow, the pass loads, the prime-search caps and the sums charged
in Section 6 cost at most `3.0e-11` in total, so `P >= 0.3954794`); peak memory `<= 2^124.1` bytes (every node counted
at its full `2^40` bytes); preprocessing `<= 2^115.83` units (included); no advice. Unit = one batched incomplete affine
G1 addition = 2000 AArch64 instructions (`C = 2000`). Only `T` and `A_d`, `d = 2^20`, are read, so the bound holds for
every `n >= 2^20`.

## 1. Setting

`r` is the BLS12-381 scalar prime, `G1 = <P0>`. `r - 1 = 2^32 * 3 * 11 * 19 * 10177 * ...`; `d = 2^20`, `m = (r-1)/d`
(`log2 m = 234.857`), `H = {x : x^m = 1}` the subgroup of `d`-th powers. `mu6 = {u : u^6 = 1}`; since `6d | r-1`, `mu6 subset H`.
The GLV map `phi(x,y) = (beta x, y)` acts on `G1` as `[omega]`, `omega` a primitive cube root of unity; with negation the
orbit of `[v]P0` under `{+-1, +-phi, +-phi^2}` is `{[u v]P0 : u in mu6}`. `N` below is the number of elements per side and
`n_b` the number of base primes per side (not the input parameter `n`, not the input point `T`).

## 2. Reduction (Cheon, d | r-1)

For `s` uniform in `{1..r-1}`, `z := s^d` is uniform in `H` and `A_d = [z]P0`.

*Lemma 1.* Given `z`, the secret `s` is recovered with at most `2^21` units.
Proof. Take 20 successive square roots `z_{i+1}^2 = z_i`, `z_0 = z = w^{2^20}`. Each `z_i` (`i <= 19`) is
`eps_i w^{2^(20-i)}` with `eps_i` of order dividing `2^i`; both factors are squares (`2^19 | (r-1)/2`; `20 - i >= 1`),
so all roots exist and `s_0 = z_20` has `s_0^d = z`. The `d`-th roots are `s_0 omega_d^j`, `j < 2^20`; exactly one has
`[s_j]P0 = T`. With `R = [s_0]P0` and `j = j_1 + 2^10 j_2` we need `[omega_d^{j_1}]R = [omega_d^{-2^10 j_2}]T`; compute
both sides for all `j_1, j_2 < 2^10` (`2^11` scalar multiplications of at most 383 reference operations each) and match
in a table: `< 2^20` units. The square roots cost `< 2^21` instructions. QED

*Lemma 1' (input randomization; makes the secret uniform over the group).* The attack draws `t` uniform in `F_r^*` and
replaces `T` by `T' = [t]T = [ts]P0` and `A_d` by `A'_d = [t^d]A_d = [(ts)^d]P0`; it recovers `s' = ts` by Lemma 1 and
outputs `s = s' t^{-1} mod r`. The map `x -> x^d` on `F_r^*` has image `H` and kernel of size `d`, so `t^d` is uniform on
`H`; therefore `z' = (s')^d = t^d z` is uniform on `H` for every fixed input. Cost: 2 scalar multiplications and one
`F_r` inversion, inside the `coins_input_lemma1` row of Section 7 (charged `2^22` units, about `2^21` more than Lemma 1
alone). Lemma 1 applies to `(T', s')` unchanged.

## 3. Baby and giant sets, and the capped H-prime search

**H-primes.** A prime `q` is an *H-prime* if `q^m = 1 mod r`. **H2:** the classes `q mod r` of primes `q << r`
equidistribute in `F_r^*/H` (index `d`): a fraction `1/d` of primes are H-primes, i.e. density `1/(d ln q)` among
integers, about `X/(d ln X)` below `X`, with gaps like independent exponentials of mean `d ln q`. Checked on this `r`
(`evidence.py ledger`): of the 6542 primes below `2^16`, 402 (0.06145) are 16th powers (expected 0.0625) and 27
(0.00413) are 256th powers (expected 0.0039); of the 5709 primes in `(2^16, 2^17)`, 0.06516 and 0.00473.

**Capped search, and what it is a bound on.** Both lists are built by one procedure, and its cost is charged at the
cap, not at the expected count:
1. Scan the odd integers of the interval in order, in segments of `2^42`, and **sieve each segment by every prime below
   the square root of the interval's top** — `2^48` for the base lists, whose `2^94` odd candidates above `2^77` span
   `2^95` and never reach `2^96`, and `2^40` for `W`, whose `2^78` odd candidates above `2^39` span `2^79`. A survivor
   of this sieve has no prime factor `<= sqrt(q)`, so **every survivor is prime by complete trial division**; no
   primality assumption, no Miller-Rabin, no base set. (Cost in Section 7.2: `2^96.93` instructions of crossing out
   and `2^96.58` for the per-segment first multiple and 8-byte table read for the base lists.)
2. Every surviving prime `q` is tested for membership in `H` by one square-and-multiply of the 234-bit exponent `m` mod
   `r`: at most 468 `F_r` multiplications of 128 instructions each (Section 7.2), `59904` instructions `= 29.95` units.
   The test is charged once per survivor, i.e. once per prime of the interval — `pi(2^95) <= 2^88.99` for the base
   lists and `pi(2^79) <= 2^73.26` for `W` (Dusart's bound `pi(x) <= (x/ln x)(1 + 1.2762/ln x)`, and the interval lies
   below `2^95` and `2^79` respectively) — not once per *accepted* H-prime.
3. **Why not Miller-Rabin.** A published deterministic base set reaches `n < 3,317,044,064,679,887,385,961,981 =
   2^81.46` with the first 13 primes `{2,3,5,7,11,13,17,19,23,29,31,37,41}` (Sorenson and Webster, *Math. Comp.* 86
   (2017) 985-1003); no published set covers `2^96`, and 12 unspecified bases cover nothing below `2^96` at all. Cutting
   the base primes down to `2^81.46` cannot be done: `H` has index `d = 2^20` in `F_r^*`, so only `pi(2^81.46)/2^20 =
   2^55.67` H-primes exist below that bound against the `2 n_b = 2^66.664` the construction needs, `11` bits short.
   Complete trial division needs no published bound, holds for every candidate the cap can reach, and is *cheaper*
   here than the sieve-to-`2^16` plus Miller-Rabin it replaces (the price of that pair is in the paragraph below), so
   the primality of every accepted base value is established outright and no false-accept term exists.
4. **Candidate cap and abort rule.** `W` (the `2^50` walk H-primes above `2^39`) scans at most `2^78` odd integers; the
   base lists (`2 n_b = 2^66.664 H-primes in (2^77, 2^96)`) scan at most `2^94`. If a cap is reached with too few
   H-primes the run **fails** (there is no retry, no second interval and no hidden advice); that failure is charged
   `2^-40` in the success term of Section 6. The caps exceed the expected scan by a factor of about `5` (`2^94` odd
   candidates against the `2^91.7` needed for `2^66.664` H-primes at density `1/(d ln q)`) and about `7` for `W`
   (`2^78` against `2^75.1`), so an abort needs the H-prime density to be more than 5x below the H2 prediction, and
   the two measured density ratios above are 0.983 and 1.043. The H-primes the caps can deliver are `2^68.99` for the
   base lists (from `2^88.99` primes) and `2^53.26` for `W`, both above what is needed.
   Charged: `2^79` units for `W` (the search itemises at `2^78.17`) and `2^95` for the base lists (itemised at
   `2^93.90`), both including the sieve, the H-test on every prime of the capped range and the scan of the interval.

**What the earlier form of this row got wrong, and the numbers it should have carried.** The row used to price the
H-test as if only *primes* survived a sieve to `2^16` (`A_ODDS/31` and `W_ODDS/27`). A sieve to `2^16` leaves
`2 prod_{p<2^16}(1 - 1/p) = 2 e^{-gamma}/ln(2^16) = 0.101` of the odd candidates, not `0.0304` (the prime density),
so at the cap that is `0.101 * 2^94 = 2^90.69` H-tests `= 2^95.59` units for the base lists, above the whole `2^95`
allocation, and another `2^95.55` units of Miller-Rabin — even before any primality argument. `W` would carry
`0.101 * 2^78 = 2^74.69` H-tests `= 2^79.59` units plus `2^79.55` of Miller-Rabin, against its `2^79` charge. The
repaired procedure above (sieve to the square root, so survivors *are* primes, one H-test per prime) costs
`2^93.90` and `2^78.17` units for the two lists, which is `1.7` and `1.4` bits *below* the rejected figures and is
what the ledger now charges; `prime_search` in `evidence.py` recomputes both from the sieve bound, the Mertens sum
`sum_{p<X} 1/p <= ln ln X + 0.30` and the `pi` bound.

**Walk list `W`:** the first `L = 2^50` H-primes above `2^39`, `w_0 < ... < w_{L-1}`, even gaps `g_j`. **H2'** (checked
by the computation that builds `W`): `w_{L-1} < 2^77` and all `g_j <= G = 2^33` (under H2 a gap exceeds `2^33` with
probability `< exp(-2^7.3) < 2^-228`, over `2^50` gaps `< 2^-178`).
**Base primes:** H-primes in `(2^77,2^96)` (about `2^69.5` exist by H2, and the cap's `2^88.99` primes deliver
`2^68.99` of them; `2 n_b ~ 2^66.664` are needed) dealt alternately into a baby list `a_i` and a disjoint giant list
`a'_i`, `n_b = ceil(N/L)` each.
**Elements:** baby `(a,w)`: `x = a w`, point `[x]P0`; giant `(a',w')`: `x' = a' w'`, point `[x']A_d = [x' z]P0`.

*Lemma 2.* (i) All `x, x'` lie in `H`. (ii) `x` determines `(a, w)` (`a` is the only prime factor above `2^77`), so the
`N` scalars of a side are distinct. (iii) For pairs `(x,x')`, `(y,y')`, `x = aw, x' = a' w', y = bu, y' = b'u'`:
`x/x' = y/y'` as rationals implies the pairs coincide or `a = b, a' = b', w = w', u = u'` ("diagonal", ratio `a/a'`);
diagonal pairs are a `1/L = 2^-50` fraction.
Proof. (i) closure; (ii) unique factorisation; (iii) `a w b' u' = b u a' w'`: the primes above `2^77` give
`{a,b'} = {b,a'}`, and baby/giant bases are disjoint, so `a = b`, `a' = b'`; then `w u' = u w'` gives `{w,u'} = {u,w'}`.
QED Note that (iii) is a statement about integers, so it is unchanged by the uniform base scale `tau` of Section 6,
which multiplies both `x` and `x'` by the same value.

## 4. Walk, key, match

**Chain.** For a base prime (baby `Q = [a]P0`, giant `Q = [a']A_d`) build `tau[g] = [g]Q` for even `g <= G` (`tau[2]` by a
doubling, `tau[4] = 2tau[2]`, `tau[g] = tau[g-2] + tau[2]`) and walk `X_{j+1} = X_j + tau[g_j]`, `X_j = [w_j]Q`. For
batching on one node the walk is cut into `2^16` equal segments started by scalar multiplications `[w_s]Q` (the table
likewise into `2^16` blocks); the segments advance in lockstep and share one Montgomery inversion per round: the batched
incomplete affine reference operation (6M).

*Lemma 3.* Every addition has distinct x-coordinates and no infinity operand.
Proof. `0 < g <= 2^33 < 2^39 < w_0 <= w_j < r/2` gives `w_j != +-g mod r`, so `X_j != +-tau[g]` as `Q` has order `r`; in
the table `g - 2 != +-2` for `g >= 6`, and likewise for block starts `g_0 >= 6`. QED

**Key.** `y*` = the even one of `y, p - y` (`p` odd); `kappa(Y)` = bits 1 to 224 of `y*`.
*Lemma 4.* Points in one `phi`/negation orbit have equal `kappa`; equal `y*` implies one orbit (`y' = +-y` gives
`x'^3 = y'^2 - 4 = x^3`). Different orbits agree on `kappa` only by chance (H3').

**Match.** A baby `[x]P0` and giant `[x'z]P0` in one orbit give `x = u x' z`, `u in mu6`, so `z = x (u x')^{-1}`. Try the six
`u`; accept the one with `[z]P0 = A_d` (at most 6 scalar multiplications of 383 reference operations: `<= 2^12` units). An
orbit match always yields the true `z`; an equal-`kappa` match of different orbits fails the test.

## 5. Machine model and matching

**Machine.** A shared-nothing cluster with three kinds of node, every count derived from the record capacity `C` of
Section 5.1: `S = 2^84` forwarding and owner nodes, `2^20` generating nodes and `2^31` storage nodes (for `W` and the
base lists). A node is a scalar AArch64 core with `2^40` bytes of local RAM, addressed by ordinary 64-bit loads and
stores; **no instruction addresses another node's memory**, so every byte that crosses a node boundary is a message
and is priced below. Nodes exchange blocks of `b = 25` records (1200 bytes); hardware moves the bytes, software charges
at most 60 instructions to post and 40 to poll a block.

**Addresses, and what a pass does.** A record's 84-bit owner id is its key bits 1..84, read as three 28-bit digits
`(a_1, a_2, a_3)`; the forwarding nodes are the same `2^84` triples. Pass `k` is an exchange along axis `k`: a node
holding records of digit `d_k = v` sends them to the node that differs from it only in coordinate `k`, which is set to
`v`. After pass 3 the record with owner id `(a_1,a_2,a_3)` is on node `(a_1,a_2,a_3)`: pass 1 fixes `a_1` over the whole
node set, and every node of the axis-2 line `(a_1, . , v_3)` already carries `a_1` in coordinate 1, so pass 2 can fix
`a_2` on it and pass 3 fixes `a_3`. Each pass has **`2^28` destinations per node**, not `2^84`: a node's open-block area
holds one 25-record block per destination digit, `2^28 * 25` records `= 2^38.229` bytes, and that area is the largest
that fits a `2^40`-byte node beside two full record buffers (Section 5.1). Fan-out, not arithmetic, is why three passes
are needed, and it is why a pass is a digit exchange rather than a direct scatter to the owner.

**Where a pass's records are before it runs.** A generating node holds no records. It walks chains and posts each
record as it is produced: it holds `2^46.664` chains of `2^50` records and streams `2N/2^20 = 2^96.6644` records in
total, so its store is the open-block area (`2^38.229` bytes) and nothing else. Its posting chooses the receivers of
pass 1: for a record of digit `a_1` the block goes to the receiver `(a_1, r_2, r_3)` whose `r_2, r_3` are the high and
low 28 bits of a per-digit posting counter, so successive blocks for one digit walk the `2^56` nodes of the slab
`(a_1, . , .)` in a fixed round robin. The sender never reads digits 2 and 3, so the assignment is independent of the
key, and each generator's block count per receiver is `floor` or `ceil` of its own count for that digit over the
`2^56` nodes of the slab.

**Per-node inbound load at every stage, from `C`.** `2N = 2^116.6644` records. Section 5.1 fixes `C = 2^32.7549`
records in a `2^40`-byte node; hence a pass needs at least `2N/C = 2^83.9095` receiving nodes, and `S = 2^84` (a
product of three 28-bit digits) is the smallest admissible count. `mu = 2N/S = 2^32.6644` is the mean inbound load of
every routing pass and also the mean final bucket, and `C/mu = 1.0648`:

| stage | nodes | inbound records per node | where they are held while the next stage runs |
| --- | ---: | ---: | --- |
| generate (pass-1 senders) | `2^20` | `2^96.6644` streamed, none held | open-block area, `2^38.229` B |
| pass 1, digit 1 | `2^84` | `<= mu + 25*2^20 = 2^32.6699` given its digit slab's count; `> C` with probability `<= 2^-2^76` | pass buffer A |
| pass 2, digit 2 | `2^84` | `<=` its digit slab's count / `2^56`, `> C` with probability `<= 2^-2^76` | pass buffer B, while A is routed out |
| pass 3, digit 3 | `2^84` | the same bound | pass buffer A, while B is routed out |
| owner bucket, sort and scan | `2^84` | `mu`, `> C` with probability `<= 2^-1.68e6` | sort pair A, B (count into B, scatter back) |
| storage of `W` and base lists | `2^31` | `2^39.664` B of 16-byte fields | local RAM, `<= 2^40` B |

**The inbound of pass 1, stated as what it is.** Round robin runs *inside one first digit*. A generator `g` produces
`n(g,a)` records of digit `a` and spreads them over the `2^56` nodes of the slab `(a, ., .)`, so one receiver `r` of
that slab gets `floor` or `ceil` of `n(g,a)/2^56` records from `g`, `<= n(g,a)/2^56 + 1`. Summed over the `2^20`
generators and adding the final block each of them sends the receiver (at most `25` record slots), the inbound is

```text
inbound(r in slab a)  <=  sum_g n(g,a)/2^56 + 2^20 + 25*2^20  =  N_a/2^56 + 26*2^20 ,
```

where `N_a = sum_g n(g,a)` is the number of records whose first digit is `a`. This is **not** the unconditional
`mu + 25*2^20`: `N_a` is a random count over the digit slabs, and a slab that is handed more records than `C*2^56`
overloads its nodes. `N_a` is the same quantity for every pass — the records a receiver of a digit exchange must hold
are exactly those whose routed digit equals its own value — so passes 1, 2 and 3 get one bound, derived from `N_a`
alone (Section 6): with `a = (C - 26*2^20)/mu = 1.0608`, `N_a` exceeds the threshold with probability at most
`2^-2^76` per slab and, after the union over the `3 * 2^28` (pass, digit) pairs, at most `2^-2^76` in all. The
`26*2^20` of rounding is what the reviewer's `sum_g n(g,a)/2^56 + 25*2^20` becomes with the `+1` per generator of the
`floor`/`ceil` split made explicit; the ledger's `pass1_end_tokens` row (one marker block per (generator, owner) pair)
prices the marker traffic and is unchanged. Section 6 charges this failure probability in the success line, so pass 1
is no longer claimed to be deterministic. Because the threshold factor `1.0608` still exceeds 1 by a wide margin in
the exponent, `S = 2^84` does not change and neither does any node count, row or the `S x 2^31`-style fixed work.
Both record buffers are `C` records; a node reads one while it fills the other, and at the end the pair is the sort
buffer. No stage of the run asks a node to hold more than `C` records except on the `2^-1.68e6` (bucket) and
`2^-2^76` (pass) events charged in Section 6, and the bound is recomputed from `C` at each stage instead of being read
off the final owner load.

**Record.** 48 bytes: `kappa` 28 B, label 16 B (67-bit base index over the `2 n_b` base primes, 50-bit step, type
bit), 4 B padding. The label holds **indexes, not values**, so the record size does not depend on the prime width.

**Stored primes: one 16-byte field, and the read path.** `W` and the base lists are arrays of **16-byte fields**: bits
0..96 hold the prime (a base prime is below `2^96`, the capped scan of Section 3 ending below `2^95`; a `W` value below
`2^77` uses the same field), bits 112..127 hold a 16-bit tag (kind: walk, baby base, giant base, with the parity bit
used as a check), bits 97..111 are zero. A 96-bit value therefore needs no truncation and loses no tag bits, which a
10-byte field could not do. The arrays live on the `2^31` storage nodes (`2^54` bytes for `W`, `2^70.664` bytes for the
base lists, `2^39.664` bytes per storage node), together with the `2^43`-entry sieve table of Section 3 (`2^46.4` bytes
in all, 41 bytes per storage node).
A label is resolved by a message, not by a remote load: the owner sends a request naming the index and receives the
16-byte field back. That read happens **only inside the verification cap** (two fields per verified pair, two
round trips of about 250 instructions each, against the `2^12` units = `8.192*10^6` instructions charged per
verification), and once per chain for the base prime `a`. The walk stream reads `W` sequentially, one field per
record, and needs no prefix sum: the chain accumulates the gap `w_{j+1} - w_j` from the two fields it already holds.

**Routing pass** (per record, on one 28-bit digit). Every destination has its own open block at a fixed address, so the
loop keeps the record's destination computed and written back:

```text
d  = (kappa >> 28*k) & (2^28 - 1)      extr, ubfx                         2
o  = off[d]                            ldr   static offset, 1200*d         1
c  = fill[d]                           ldr   this destination's fill count 1
b0 = area + o                          add   block base in the open area   1
a  = b0 + 48*c                         madd  x48 held in a register        1
ldp x3 from src; stp x3 to [a]         copy                                6
fill[d] = c + 1                        add, str   fill-count writeback    2
cmp c+1, #25; b.ne                     block-full test                     2
src += 48; subs; b.ne                  loop, source pointer advance        3
if the block is full: post(b0), fill[d] = 0, wait for a credit   100/25 = 4
```

= **19 control + 4 posting + 6 transport = 29** per pass. The store address `b0 + 48*c` is formed by two instructions,
both listed: the `add` of the block base `area + off[d]` and the `madd` that adds `48 * fill[d]` to it; `48` is loaded
once per pass (node setup row), and the address needs no carry instruction because `b0 <= 2^38.23` and
`a = b0 + 48*24 < 2^39 < 2^64`. The post is counted per block (Section 5.2), the transport copy is the receiver's copy
out of the ring into its pass buffer (Section 5.3), and the credit wait is inside the 60 of the post. The fill-count
writeback is the `str` above. Pass 1 uses the same loop; its extra work is the receiver choice of the last paragraph,
two instructions per block inside the 60 of the post.
**Owner sort drain** (per record, after the third pass): the poll of every block arriving in the third pass returns the
ring slot and the record count `n`; the records have been copied into the pass buffer by the third pass's transport
copy, and the loop then runs `n` times: capacity test of the cursor against the last record slot of the `C`-record sort
buffer, before any write (`cmp`, `b.hi` to *overflow*: 2), 3 `ldp` from the pass buffer and 3 `stp` to the sort buffer
at the cursor, both post-indexed so the source pointer and the cursor advance by 48 (6), loop (`subs`, `b.ne`: 2):
charged **10**. The local sort reads only the sort buffer.
**Local LSD pass**: counting loop (`ldr`, `extr`, `ubfx`, `ldr/add/str`, loop: 9) + scatter (`extr`, `ubfx`, `ldr`
offset, `madd`, `add`, 6 copies, `add/str`, loop: 16) = 25: charged **25**.

**Matching.** Each record is routed by 3 routing passes on the first 84 bits of `kappa` to its owner (the `mu` of the
Section 5.1 table, one 84-bit prefix), sorted by 2 LSD passes on `kappa[85..140]` (28-bit digits) and scanned once.
Records with equal 140-bit prefix are adjacent (a run); inside a run all opposite-type pairs are compared on the full
224-bit `kappa`, and pairs with equal `kappa` are verified by the Match step. **Caps, per node:** every owner keeps two
64-bit counters, run-pair comparisons and verifications; a node that would make its `2^20 + 1`-th comparison or its
`129`-th verification stops and signals *fail* (the counters never exceed `2^20 + 1` and `129`, far below `2^64`). No
global counter and no combining step exist: each cap is local, and the failure of any node is a failure of the run.
They are charged at the caps on every owner: `2^84 * 2^20` comparisons of 32 instructions (`2^98.03` units) and
`2^84 * 2^7` verifications of `2^12` units (`2^103` units). A comparison, itemised: pair loop (index `add`, `cmp`,
`b`: 3), partner record address (1), partner type bit (`ldrb`, `tbz`: 2), both 28-byte keys (`ldp`, `ldr`, `ldr w`
each: 6), 224-bit compare (4 `eor`, 3 `orr`, `cbnz`: 8), counter (`add`, `cmp`, `b.hi`: 3): 23, charged 32. A
verification resolves its two labels to 16-byte fields by message (Section 5, read path), so the two fetches are inside
its `2^12` units. All other work is fixed: the number of records and passes does not depend on the data.

### 5.1 Placement, capacity and buffers

`2N = 2^116.6644` records. Every node that forwards records, owners, generating nodes and storage nodes alike, holds the
same four record stores, all concurrently resident, each with one lifetime:
1. the **open-block area**, `2^28` blocks of `b = 25` records = `2^38.229` bytes, one block per destination at the
   fixed address `1200 d`; it holds this node's outgoing partial blocks and is reused as each block is posted;
2. the **pass buffer**, `C` records (48 bytes each), holding the inbound of the pass being received;
3. the **second record buffer**, `C` records, holding the inbound of the pass being routed out (the pair alternates
   roles per pass and, at the end, is the sort buffer: counting pass into one, scatter back into the other);
4. the **receive ring**, `2^20` slots of one block (`2^31.6` bytes), holding blocks in flight.
Rounding out the allocation: LSD count and offset arrays `2 * 2^28 * 8 = 2^32`; routing offset and cursor tables
`2 * 2^28 * 8 = 2^32`; code and stack `<= 2^30`; reserved `2^36`. The auxiliary total is
`403726925824 = 2^38.5546` bytes, so `C = (2^40 - 4.037e11)/96 = 7247757312 = 2^32.7549` records and
`96 C + aux = 2^40` exactly (the two record buffers are `96 C`). **`C` is what fixes the node count**: a pass must
receive `2N` records somewhere, so it needs at least `2N/C = 2^83.9095` receiving nodes, and `S = 2^84` — a product of
three 28-bit digits, the smallest such count above that bound — gives `mu = 2N/S = 2^32.6644 = 6.81 * 10^9` records per
node and `C/mu = 1.0648`. See the stage table of Section 5 for the load at every stage; a node whose incoming pass
would exceed `C` records signals *overflow* and the attack fails (the drain's capacity test; the cursor is a byte
address below `2^40`, the record count at most `C + 1 < 2^33`, both 64-bit). By the Chernoff bounds of Section 6 that
costs at most `2^-1.68e6` over the final bucket and `2^-2^76` over the three routing passes (the digit-slab bound, which
covers pass 1 as well).

A generating node holds a chain table `2^32 * 96 B = 2^38.585`, the same open-block area `2^38.229`, stream buffers
`<= 2^30` and code `<= 2^30`: `<= 2^39.4 < 2^40` bytes, and no record store at all, because it posts every record as
it is produced. A storage node holds `2^39.664` bytes of 16-byte fields for its share of `W` and the base lists.

### 5.2 Blocks sent per pass, summed over the actual nodes

**The posting charge is one global total over the run, not a per-node bound.** A node routes on one 28-bit digit, so
it has `2^28` destinations, while a generating node's pass-1 receivers are the `2^84` owner nodes (Section 5), and
nothing below assumes otherwise. Let `n_{v,d}` be the record count of node `v` for destination `d` in one pass and
`M_v = sum_d n_{v,d}`, which is `2^96.6644` at a generating node and `mu` at an owner. Each destination receives
`ceil(n_{v,d}/25)` blocks: its full blocks plus the one *final* block holding its `0..24` remaining records (sent
whether or not the destination is empty, since it is also the end-of-pass marker of Section 5.3). Summed over every
node of one pass, with `sum_v M_v = 2N` because every record of the attack is in exactly one node during a pass (a
generating node's records are produced into its own open blocks and counted in its own `M_v`),

```text
full blocks      = sum_v sum_d floor(n_{v,d}/25) <= sum_v M_v/25 = 2N/25
final blocks     = sum over forwarding nodes of 2^28 per pass, plus 2^20 x 2^84 for pass 1
                   (its sender set is the generating nodes, whose marker goes to every owner)
total posting    = 3 x (2N/25) x 100 + 321 x (2^84 + 2^20 + 2^31) x 2^28 + 107 x 2^20 x 2^84
                 = 2^120.2493 + 2^120.3264 + 2^110.7415 = 2^121.2894 instructions  (= 2^110.3236 units) .
```

`#nodes = 2^84 + 2^20 + 2^31` is counted at every node that can hold a store, generating and storage nodes included,
each charged a full owner's rows, so the `2^28`-per-pass final-block term is an over-charge for the nodes that send no
records; the third term is pass 1's own marker, one block per (generating node, owner) pair. The ledger carries
exactly these three terms and nothing else, under the names `100/25 = 4` instructions of full-block posting per record
per pass (60 to post plus 40 to poll a 25-record block, inside the routing row of Section 7.1, `3 x 4 = 12` per
record), the final-block row `3 x 2^28 x 107` instructions per node and the `pass1_end_tokens` row
`107 x 2^20 x 2^84` (Section 7.2); `evidence.py` asserts that the three rows equal this total. Since `48 x 25 = 1200`
bytes and a final block is at most 24 records, no block is ever double counted. No statement of the form "a node
handles at most `mu` records" is used: the per-node inequality `sum_d ceil(n/25) <= M_v/25 + 2^28` holds for every
`M_v`, and the global rows above depend only on `sum_v M_v = 2N`, on the node count and on pass 1's sender set. `mu`
enters only the capacity of Section 5.1 and the per-stage load table of Section 5.

### 5.3 Communication schedule: flow control, buffers, termination

**Flow control.** A ring slot is a credit. A sender posts a block only against a credit it has been granted for that
destination, and the receiver's poll returns the slot as a credit once the block's records have been copied out. Post
and poll are charged 60 and 40 instructions per block (Section 5.2), which cover the credit test, the header write, the
slot release and the wait.

**Why no sender waits forever.** Take one pass and a node `w` receiving it. The buffer `w` fills for that pass is empty
when the pass starts, its capacity is `C` records, and the whole inbound of the pass is at most `C` records with
probability at least `1 - 2^-2^76` (the digit-slab bound of Section 5.1 for passes 1, 2 and 3 alike, pass 1's
`26*2^20` of block rounding included; a node whose inbound would exceed `C` signals *overflow* and the run fails, and
that event is one of the terms of the success line of Section 6). So `w`'s copy-out of a ring slot (6 instructions per
record, the transport copy of the routing row) is never blocked by storage while the pass runs, and a slot is returned
after at most 25 records of work, whatever `w`'s own sending state is. Sending and receiving use disjoint stores: the open blocks and the buffer being routed out
are `w`'s, the inbound lands in the other record buffer. A node that has filled a destination's block posts it and, if
no credit is free, waits; meanwhile its receiving continues. So blocked senders wait only on receivers, no receiver's
copy-out ever waits on a sender, and the wait-for graph has no cycle: the schedule is deadlock free without any global
order, barrier or priority rule. A node has at most one block outstanding at a time (the routing loop posts a block
before refilling it), so pass 1's `2^20` senders fit the ring exactly; in passes 2 and 3 the sender set is `2^28`, so
at most `2^28 - 2^20` senders wait for a credit at a node, each granted one within one block's copy-out.

**Termination and detection.** Each node keeps three 64-bit counters (charged in `node_setup`). The pass counter counts
the final blocks it has received in the current pass. The senders of a pass are fixed by the digit route: the `2^28`
nodes that agree with the receiver in the other two digits, except in pass 1, whose senders are the `2^20` generating
nodes, each of which sends one final block to every one of the `2^84` owners (this is the `pass1_end_tokens` row). The
count is exact because every sender sends exactly one final block per destination, so a node knows the pass is complete
when the count reaches its pass's sender count, `2^20` then `2^28` twice. A node then clears its counters and starts
the next pass; the completion of pass 3 is the end of the routing, after which each owner sorts and scans. `fill[d]`
and `off[d]` are reset per pass by the routing offset init row of Section 7.2, and the open blocks are reused in place,
so there is no allocation during a pass.

**Peak memory** is the four stores above at one node, `96 C + aux = 2^40` bytes exactly (Section 5.1); the credits are
the ring slots themselves and need no extra store.

## 6. Algorithm and success probability

```text
1. Build W and the base lists (capped search of Section 3); read A_d, T.
2. Draw tau, t uniform in F_r^*; set Q = [tau a]P0, Q' = [tau a']A_d for every base prime; T' = [t]T, A'_d = [t^d]A_d.
3. Generate all N baby and N giant records (each generating node walks one chain as 2^16 segments); route each record as
   produced. Run the attack on (T', A'_d) and recover s' = ts; output s = s' t^{-1} mod r.
```

**The two coins, and what they buy.** `t` makes the secret uniform over `H` for every input (Lemma 1'), so the success
event is not a statement about a random `s`: for a fixed input the algorithm succeeds on exactly the coin values with
`z' = t^d z` in `S`, the set of values `x (u x')^{-1}`, and `z'` is uniform on `H`, so

```text
P(success) = |S| / m      exactly, for every input;  |S| <= 6N^2,  and 6N^2/m = K = 0.52 .
```

`tau` multiplies every base point and therefore every record key by the same scalar, so it changes no ratio (Lemma
2(iii)) and it makes the multiset of keys a function of the algorithm's own coins rather than of a fixed key set; this
is the ensemble H3' below is stated on.

**H1 (success set).** The `6N^2 = 0.52 m` values behave like independent uniform elements of `H`, so
`|S| = m (1 - e^{-K})`: support. (a) Lemma 2(iii) proves the values are distinct as rationals except a `2^-50`
diagonal. (b) The toy measurement of `coverage` matches the birthday prediction at three primes (below). (c) The
`distinct` measurement counts the occupied residues directly: at `t28` (`N = 780`, `m = 2^22.63`) 571652 of the 608400
`u = 1` ratios are distinct, a loss of 6.04% against a predicted 4.5% birthday loss plus a `1.65%` diagonal loss
(`1/L = 1.67%`), a deviation of 0.11 percentage points. At production the diagonal term is `2^-50`, so the birthday
model is the whole loss. (d) Within the model, the occupancy concentrates: `|S|` changes by at most 1 when one value
moves, so McDiarmid gives `P(|S| <= m(1-e^{-K}) - t) <= exp(-2 t^2 / 6N^2)`, which for `t = 0.01 m` is
`exp(-4*10^-4 m) ≈ 2^-2.9e67`.
**Model allowance.** The claim charges `0.01` as a flat allowance for the model itself: `|S| >= (1 - e^{-K})m - 0.01 m`.
That is 10x the largest deviation measured at toy scale (0.11 percentage points) and `2^117` times the model's own
fluctuation (`O(sqrt(m)) = 2^117.4`, i.e. `2^-117.4` of `m`). There is no unexplained constant: the earlier `2^-20`
"modelling error" is replaced by this measured allowance plus the computed tail terms below.

**H3' (keys, on the coins).** For fixed inputs and the coins `(tau, z)`, the records' 224-bit keys and their 56-bit
windows (bits 85..140) behave like independent uniform draws, jointly over the records of a node. Support: Lemma 4 is
exact on orbits, so only different orbits can collide; the toy `ec` run bins keys into 32 pseudo-nodes and measures
chi-square `z`-scores in `[-0.79, +0.98]`; and the `coverage`/`ec` agreement shows no sign of structure. Because every
key is `kappa` of `[tau x]P0` or `[tau x' z]P0` with `tau` uniform, the ensemble is over the algorithm's coins, as the
review standard requires, and the baby records are not held fixed.
**Runs and orbits.** An orbit carries at most 12 records: its 6 elements of `H`, each at most once per side (Lemma
2(ii), and distinct `x'` give distinct `x' z`), so at most 6 baby and 6 giant records, and two orbits give at most
`6*6 + 6*6 = 72` opposite-type record pairs.
**Overflow of a bucket or of a pass.** Orbits fall independently and uniformly on the `2^84` owners (key bits 1..84),
each adding 1 to 12 records, so the Chernoff bound for independent `[0,12]`-bounded summands with mean `mu` gives
`P(load >= C) <= exp(-mu (a ln a - a + 1)/12)`, `a = C/mu = 1.0648`, i.e. `2^-1.68e6` per node and still `2^-1.68e6`
after the union over `2^84` nodes (`evidence.py` prints the exponent); this covers the final bucket. A routing pass is a
digit exchange on uniform keys (H3': the routed digit is uniform given the coins, and a pass-1 sender's choice of the
other two coordinates never reads the key), so the records a receiver must hold are those of its own **digit slab**:
`N_a` for the slab `a`, whose mean is `2N/2^28 = mu * 2^56 = 2^88.6644` and which is a sum of the same
`[0,12]`-bounded orbit contributions. With the `26*2^20` of round-robin and marker rounding of Section 5.1 the receiver
overflows only if `N_a > (C - 26*2^20) * 2^56`, a factor `a = 1.0608` above the mean, so

```text
P(overflow in one slab)  <= exp(-2^88.6644 (a ln a - a + 1)/12) <= 2^-2^76 ,   a = 1.0608 ,
```

and the union over the three passes and the `2^28` digit values (`2^29.6` slabs) leaves it at `2^-2^76`. This one
bound covers pass 1, pass 2 and pass 3 — pass 1 has no deterministic bound and is not claimed to have one.
**Per-node caps.** Fix an owner holding `M <= C` records (else it has already failed by overflow), hence `M' <= C`
orbits. Given which orbits it holds (key bits 1..84), their key bits 85..140 are independent and uniform on `2^56`
values. A compared pair is an opposite-type pair in one run; a same-orbit one has equal full keys, is verified, yields
`z` and ends the attack, so at most one is compared, and every other compared pair joins two different orbits with
equal bits 85..140.
(A) Some value is held by at least 7 orbits: probability `<= C(M',7) * 2^(-56*6) <= C^7/(7! * 2^336) = 2^-119.01` per
node, `2^-35.01` over `2^84` nodes.
(B) Let `D'` be the number of orbits whose value repeats that of an earlier orbit (any fixed order). Each orbit repeats
with probability at most `M'/2^56` whatever the earlier ones did, so `D'` is dominated by `Bin(M', M'/2^56)`, mean
`<= C^2/2^56 = 729.0`. Outside (A) every value holds `j <= 6` orbits and `C(j,2) <= 3(j-1)`, so the orbit pairs with
equal value are at most `3 D'` and the compared pairs at most `1 + 72 * 3 D' = 1 + 216 D'`. Exceeding `2^20` therefore
needs `D' >= 4855`, and `P(Bin >= k) <= (e * 729.0/k)^k` gives `2^-6276` per node, `2^-6192` over all nodes.
(V) Verifications at a node are at most `1 + 72 * Q`, `Q` the number of orbit pairs with equal 224-bit keys. The cap
128 is exceeded only if `Q >= 2`, i.e. two distinct orbit pairs each agreeing on the 140 bits 85..224 (disjoint or
sharing one orbit, probability `2^-280` either way): `<= C(C(M',2),2) * 2^-280 <= C^4/8 * 2^-280 = 2^-151.98` per node,
`2^-67.98` over all nodes.

**Success.** `P(success) >= 0.4054795 - 0.01 - 2^-35.01 - 2^-6192 - 2^-67.98 - 2^-1.68e6 - 2^-2^76 - 2^-40 = 0.3954794`;
`0.39` is claimed, `0.0055` below the bound, and the tails total `3.0e-11` (`2^-2^76` is far below the rest). Primality
contributes no term: every accepted base value is proven prime by complete trial division (Section 3).

**Toy checks** (`evidence.py`, `mode all`, seed 2026; `8000/5000/3000` trials, the same trials the manifest's declared
initial input runs). Scalar model, three toy primes with `6 * 2^k | r-1`, chain length `L = 60`: measured success for
`r ~ 2^27.6, 2^29.7, 2^31.7` is `0.4255, 0.4242, 0.3850` against the diagonal-corrected Poisson `0.42493, 0.41923,
0.38940` (`|z| <= 0.71`); at toy scale the diagonal pairs are a `1/L` fraction, so the prediction uses
`K' = K(1 - (L-1)/L^2)`. The `ec` mode runs Algorithm A end to end on genuine `j = 0` curves of prime order
`r ~ 2^28..2^32`: affine incomplete additions with Montgomery batch inversion (asserting no zero denominator), the
parity key truncated to 18/20/22 bits so that false key matches occur (0.5-1.8 per run) and are rejected under a cap,
sort, equal-key runs, `z`, and the `d`-th root **selected by matching `T`** (the `d` roots are scanned with one scalar
multiplication each and the one whose point is `T` is taken; the recovered scalar is then compared with the secret):
`28/14/8` recoveries in `80/30/20` trials, all selected from `T`, success agreeing with the scalar model in every trial,
chi-square `z`-scores in `[-0.79, +0.98]`. Limits: tables are built only for occurring gaps (so `tau[g] = tau[g-2] +
tau[2]` and the segments rest on Lemma 3), and the toy sorts locally without distributed routing. Evidence at toy scale
only: it supports the Poisson form, the duplicate bookkeeping and the root selection, not the production value of `m`.

## 7. Cost ledger (unit = 2000 instructions)

### 7.1 Per record, on top of the reference addition

| item | instructions |
| --- | ---: |
| walk bookkeeping, 30 itemised: gap load and arithmetic 3, absolute-prime stream transport (16 B per 8 = 2 loads) 2, stream block post/poll amortised over a 128-field block 1, table address `base + (gap/2)*96` 3, segment loop 3, label 3, batch pointer/index updates outside the Montgomery-trick multiplications 12, stream-refill check 3 | 30 |
| key, 23 itemised: `ldp` of 6 `y` limbs 3, parity test and select mask 2, `p - y` by `subs`, 5 `sbcs`, 6 `csel` 12, `extr` to bits 1..224 4, `stp` of the key 2 | 23 |
| 3 routing passes x 29 (19 control: destination `extr`/`ubfx`, offset `ldr`, fill-count `ldr`, block-base `add`, address `madd`, six loads/stores, fill-count `str`, full test, loop; 4 for posting a full 25-record block, 6 transport copy) | 87 |
| owner sort drain, pass buffer to the sort buffer (3 `ldp` + 3 `stp` 6, capacity test 2, loop 2) | 10 |
| 2 local LSD passes x 25 (counting 9 + scatter 16) | 50 |
| run scan (load the 140-bit sort prefix, compare with the previous record, branch) | 10 |
| **itemised; charged** | **210** |

Element cost `c = 1 + 210/2000 = 1.105`.

**Posting is one global total.** The `4` in each routing pass is the full-block share `100/25` of the global count
`sum over nodes of sum_d floor(n_{v,d}/25) <= 2N/25` per pass, and the `2^28` final blocks per node per pass are not per
record but the final-block row of Section 7.2. Summed over the whole run and over every node, the posting charge is
`12 x 2^116.6644 + 321 x (2^84 + 2^20 + 2^31) x 2^28 + 107 x 2^20 x 2^84 = 2^121.2894` instructions, exactly the
sum of the share inside this row, the final-block row and the `pass1_end_tokens` row (Section 5.2 carries the
derivation and the identity). No per-node premise of the form "a node handles at most `mu` records" is used anywhere:
the per-node inequality used is `sum_d ceil(n_{v,d}/25) <= M_v/25 + 2^28`, which holds for every `M_v` (a generating
node's `2^96.6644` included), and the global rows use only `sum_v M_v = 2N` per pass, the node count and pass 1's
sender set.

**Why no percentage allowance.** The earlier form of this row listed 236 instructions and charged 260, about 10% for
"unlisted loop and spill instructions". Those instructions are named here instead of covered by a percentage: the loop
counters are the "loop 3" in each pass row and the "segment loop 3" of the walk row; the spill and reload traffic is the
"batch pointer/index updates 12" of the walk row and the 3 `ldp` + 3 `stp` copies inside each pass row; the
block-boundary work is the "block-full test 2" and the amortised posting. Everything a record does not touch per unit of
work is fixed work and is in the rows of Section 7.2: clearing and prefix-summing the histograms and the routing
tables, the final blocks of Section 5.2, the receive ring and node setup (the itemised per-node rows), and the chain
tables. For orientation, `+10` unlisted instructions per record would move the bound by 0.0065 bits.

### 7.2 Fixed terms

| term | units |
| --- | ---: |
| elements `2Nc` | `2^116.8084` |
| tables: `2 n_b` chains x `2^32` entries x (1 addition + 12 stores/2000) | `2^98.673` |
| segment and table-block starts, `[tau a]P0` (`<= 2^25` per chain, `2^12` for the two base scales) | `2^91.665` |
| per-node fixed work, `(2^84 + 2^20 + 2^31) nodes x 90768936960 = 2^36.4015` instructions, itemised below | `2^109.4357` |
| capped H-prime search in `(2^77,2^96)`: `2^94` odd candidates, sieve by every prime `< 2^48` plus one H-test per prime, itemised `2^93.90`, charged | `2^95` |
| capped H-prime search for `W` above `2^39`: `2^78` odd candidates, sieve by every prime `< 2^40` plus one H-test per prime, itemised `2^78.17`, charged | `2^79` |
| coins, input access, Lemma 1 (square roots and the `2^20`-unit root match), `t` setup | `2^22` |
| verification cap, per node: `2^84` owners x `2^7` verifications x `2^12` units (the 16-byte label fetches are inside it) | `2^103` |
| run-pair cap, per node: `2^84` owners x `2^20` comparisons x 32/2000 | `2^98.034` |
| pass-1 end-of-pass markers, one per (generating node, owner): `107 x 2^20 x 2^84` instructions | `2^99.7757` |
| walk stream set-up and base-prime fetches: `2^66.664` chains x (`2^16` segment streams + 1) x ~100 instructions | `2^79` |

**Per-node fixed rows** (instructions per node; every node, generating and storage nodes included, is charged the
owner's rows). Every array is 8-byte entries, `2^28` of them; each loop handles 16 entries per iteration and pays
`subs` + `b.ne` (2) per iteration, and `2^28` is a multiple of 16, so there is no remainder.

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

The 8-byte entries are wide enough everywhere: a fill count is at most 25, an offset below `2^39`, an LSD count or
offset at most `C < 2^33`, a final-block counter at most the number of senders of a pass, at most `2^28`.

**The two search rows, itemised.** A segmented sieve crosses out every odd multiple of every sieving prime `p`, at 2
instructions per multiple (advance by `p`, store 0), so the number of instructions is
`2 * (candidates) * sum_{p < X} 1/p`, and `sum_{p<X} 1/p <= ln ln X + 0.30` (the true constant is `0.2615`; this is the
Mertens/Rosser-Schoenfeld form). Per segment of `2^42` odd candidates, each sieving prime also costs a first-multiple
division and one 8-byte read of its entry in the sieve table (3 instructions). Base lists, `X = 2^48`,
`sum <= 3.805`:

```text
crossing out            2 x 2^94 x 3.805              = 2^96.93 instructions
first multiple + read   3 x pi(2^48) x 2^52 segments  = 2^96.58
sieve-table build       2 x 2^48 x 3.805              = 2^51.9   (negligible)
H-tests                 2^88.99 primes x 59904        = 2^104.86
total                                                 = 2^104.87 instructions = 2^93.90 units
```

`W`, `X = 2^40`, `sum <= 3.622`, `pi(2^79) <= 2^73.26` H-tests: `2 x 2^78 x 3.622 = 2^80.86`,
`3 x pi(2^40) x 2^36 = 2^72.86`, H-tests `2^89.13`, total `= 2^78.17` units. The H-test dominates both rows, as it
should: it is `29.95` units per prime of the scanned interval.

**`F_r` arithmetic (an implementation bound, not a measurement).** All arithmetic mod `r` is 4x64-bit CIOS Montgomery
on a scalar core: the product is 16 `mul` + 16 `umulh` with 24 `adds/adc`, the reduction 4 `mul` + 4 `umulh` with 20
`adds/adc`, and the final conditional subtraction 7 `subs/sbc` + 3 `csel`: 94 instructions, charged **128**. An
inversion is binary extended Euclid, `<= 1.4 * 255` iterations of 8 instructions, charged `2^16`. The one expensive
`F_r` operation in the attack is the H-test of Section 3, a 234-bit square-and-multiply: at most 468 multiplications,
`59904` instructions `= 29.95` units; it is inside the two search rows. The Match step's inversions (6 per verified
pair) and Lemma 1's `2^11` scalar multiplications are inside the verification cap and the `coins_input_lemma1` row.
Nothing in the ledger is priced from a measured `F_p` figure.

**Gap and base storage.** `W` and the base lists are arrays of **16-byte fields** (Section 5): the prime in bits
0..96 (every accepted base value is below `2^96`), a 16-bit tag in bits 112..127, the rest zero. `W` is `2^50` fields
(`2^54` bytes) and the base lists `2 n_b` fields (`2^70.664` bytes), so a record label `(base index, step)` is resolved
by a request/reply to the storage node that owns the field, which returns 16 bytes; the walk stream reads `W`
sequentially at one field (16 B = 2 loads) per record, plus one amortised block post/poll per 128 fields, which is the
`2 + 1` of the walk row, and it needs no prefix sum. Construction is the capped search of Section 3 plus writing the
arrays at 1 instruction per 8 bytes (`2^40.4` and `2^56.4` units, inside the search rows, the wider field included).
Both fit the `2^31` storage nodes: `2^39.664` bytes each of `2^40`, with the `2^46.4`-byte sieve table of Section 3
spread over the same nodes at 41 bytes each.

All rows other than `2Nc` sum to `2^109.4554`; `W = 2^116.8084 + 2^109.4554 = 2^116.817200 <= 2^116.8172`, score
`292043/2500` (`evidence.py ledger` recomputes every row with exact integers and finds the smallest integer `k` with
`W^10000 <= 2^k` by integer exponentiation, `k = 1168172`; the rounding slack is `3.1e-7` bits).
Preprocessing (lists, tables, baby generation) `<= 2^115.826 <= 2^115.83`.

**Layout choice.** The alternatives, each at its cheapest block size, with the per-record row, the per-node rows, the
node count its own capacity forces and the charge rows recomputed the same way (`evidence.py ledger` recomputes them
from the same formulas and this table is its output):

| passes x width | per-record row | owner bits / `S` | `time_log2` |
| --- | ---: | ---: | ---: |
| 3 x 28, 25-record blocks (chosen) | 210 | 84 / `2^84` | `116.8172` |
| 4 x 21, 100-record blocks | 227 | 84 / `2^84` | `116.8198` |
| 5 x 17, 100-record blocks | 253 | 85 / `2^85` | `116.8367` |

`w = 28` is forced for three passes: `w >= 29` needs `2^29` open blocks of 25 records (`2^39.23` bytes), which leaves
`C < mu` and breaks the capacity rule, and `w <= 27` needs four passes. The block size `b = 25` is the other way:
`mu/b = 2^28` exactly at `b = 25`, so the open-block area holds one block per destination at the mean load; a larger
`b` makes the area exceed `mu` records and leaves no room for `C`.

Memory counted at full allocation: `(2^84 + 2^20 + 2^31) * 2^40 = 2^124.0`, claimed `<= 2^124.1` (records alone:
`2N * 48 B * 2 = 2^123.25`; `W` `2^54` bytes and base lists `2^70.664` bytes fit the `2^31` storage nodes at
`2^39.664` bytes each). Memory has no scalar weight and is exponential, as for any BSGS-type Cheon attack.

## 8. Heuristics, evidence, provenance

H1 decides the size of the success set (with the flat model allowance of Section 6); H3' decides the failure
probability of the caps and of node overflow, on the algorithm's own coins; H2/H2' decide preprocessing, table size and
the search caps (checked by building `W` and by the two measured density ratios). The cluster model is idealised and
the instruction counts are listed, not measured on hardware. The cost is data independent except for the two caps,
charged at their limits, and the capacity of a bucket or of a pass, charged by the Chernoff bound.

`evidence.py` (one initial run `mode: all`, about 8 s; each mode can be rerun with more trials): `ledger` (exact
integers: Sections 5.1-7, the digit and block widths, the record and 16-byte prime-field layouts, the derived node
count `2N/C` and the per-stage load table (`stages`), the global posting total `2^121.2894` instructions, its identity
with the three rows that carry it and pass 1's `2^20 x 2^84` markers, every per-node row, the layout alternatives,
the capped searches with their charged and cap-priced costs, the `F_r` multiply bound, the per-node caps and their
bounds (A), (B), (V), capacity and the two Chernoff exponents (bucket and pass), the storage-node field budget, the
success line, memory, `6d | r-1`, `mu6 subset H`, the H2 density checks); `coverage` (Section 6, plus the distinct-ratio count); `distinct` (the occupied-residue measurement on
one toy); `ec` (Algorithm A end to end, with the `d`-th root selected by matching `T`). Limits: the toy sorts locally
without distributed routing, and the ledger is arithmetic on declared costs.

Provenance: the idea of walking small integers of `H` with one addition per step, the GLV six-fold key and Cheon's
`d | r-1` reduction come from the published note of submission `cbadde77`; the semiprime construction, Lemmas 2-3,
matching, caps, machine model, the two coins, the communication schedule and the overhead listing were derived here.
