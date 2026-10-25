# bls12-381-g1-ecdlp: Cheon reduction with H-prime semiprime walks and sort-based matching

Exploratory claim: a classical deterministic attack with a total-work bound at the scoring value `n = 2^20`, success
probability at least 0.39 under disclosed heuristics, and no advice. Not a practical break; the cost is exponential.

**Result.** `time_log2 <= 116.7923`; success `>= 0.39` (H1's birthday model gives `0.4054794`, the model allowance is
`0.01`, and the per-node caps, node overflow and the pass loads charged in Section 6 cost `1.2e-11` in total,
so `P >= 0.3954794`); peak memory `<= 2^124.1` bytes (each forwarding node counted at its
full `2^50` bytes and each storage node at its `2^40`); preprocessing `<= 2^115.80` units (included); no advice. Unit =
one batched incomplete affine G1 addition = 2000 AArch64 instructions (`C = 2000`). Only `T` and `A_d`, `d = 2^20`, are
read, so the bound holds for every `n >= 2^20`.

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
`F_r` inversion, inside the `coins_input_lemma1` row of Section 7 (charged `2^22` units). Lemma 1 applies to `(T', s')` unchanged.

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
   lists and `pi(2^79) <= 2^73.26` for `W` (Dusart's bound `pi(x) <= (x/ln x)(1 + 1.2762/ln x)` for the two intervals) — not
   once per *accepted* H-prime.
3. **Why not Miller-Rabin.** A published deterministic base set reaches `n < 3,317,044,064,679,887,385,961,981 =
   2^81.46` with the first 13 primes `{2,3,5,7,11,13,17,19,23,29,31,37,41}` (Sorenson and Webster, *Math. Comp.* 86
   (2017) 985-1003); no published set covers `2^96`, and 12 unspecified bases cover nothing below `2^96` at all. Cutting
   the base primes down to `2^81.46` cannot be done: `H` has index `d = 2^20` in `F_r^*`, and H2 (a density
   prediction, not a bound) puts only `pi(2^81.46)/2^20 = 2^55.67` H-primes below that bound against the
   `2 n_b = 2^66.664` the construction needs, `11` bits short.
   Complete trial division needs no published bound, holds for every candidate the cap can reach, and is *cheaper*
   here than the sieve-to-`2^16` plus Miller-Rabin it replaces (the price of that pair is in the paragraph below), so
   the primality of every accepted base value is established outright and no false-accept term exists.
4. **Candidate cap, unconditional.** `W` (the `2^50` walk H-primes above `2^39`) scans at most `2^78` odd integers; the
   base lists (`2 n_b = 2^66.664 H-primes in (2^77, 2^96)`) scan at most `2^94`. Both scans are deterministic functions
   of the fixed modulus `r` and are built before `t` and `tau` are drawn: a cap is exhausted either for this `r` for
   every input and every coin choice, or for none. That is a property of `r`, not an event over the algorithm's coins,
   so **no probability is charged for it**, and the run executes both scans in full and pays the full cap on every run
   whatever they find. The two rows are `2^79` and `2^95` units against a `2^116.7923` total, so charging them
   unconditionally costs the score nothing. The sieves **count the H-primes they find** and the construction proceeds
   on those: `W` is the first `2^50` of them above `2^39`, and the base lists are the next `2 n_b` in `(2^77, 2^96)`,
   read off the sieve's own output. Whether either list is long enough is a property of the fixed modulus that the
   search *decides*, and it is carried as premise **H2''** below rather than inferred from a prime-counting bound: no
   inequality of the form "`pi(x)/d` H-primes are available" is used, and the search is charged in full whether it
   finds the lists or not. H2
   sizes the caps, predicting `2^68.99` H-primes for the base lists (from the `2^88.99` primes a bound allows in the
   capped range) and `2^53.26` for `W`, against the `2^66.664` and `2^50` needed, about 5x and 7x headroom over the
   predicted scan; the two measured density ratios above (0.983 and 1.043) are the toy-scale check of that prediction,
   and H2 is used to size intervals, never as a probability in the success line or as a lower bound on what exists.
   Charged: `2^79` units for `W` (the search itemises at `2^78.17`) and `2^95` for the base lists (itemised at
   `2^93.90`), both including the sieve, the H-test on every prime of the capped range and the scan of the interval.

**What the earlier form of this row got wrong.** It priced the H-test as if only *primes* survived a sieve to `2^16`: that
sieve leaves `2 prod_{p<2^16}(1 - 1/p) = 2 e^{-gamma}/ln(2^16) = 0.101` of the odd candidates, not `0.0304`, so at the
cap it charges `0.101 * 2^94 = 2^90.69` H-tests `= 2^95.59` units for the base lists plus `2^95.55` of Miller-Rabin, and
`2^74.69` H-tests `= 2^79.59` units plus `2^79.55` for `W`, above the `2^95` and `2^79` allocations. The repaired
procedure (sieve to the square root, so survivors *are* primes, one H-test per prime) costs `2^93.90` and `2^78.17`
units, `1.7` and `1.4` bits *below* the rejected figures, which is what the ledger charges; `prime_search` in
`evidence.py` recomputes both from the sieve bound, the Mertens sum `sum_{p<X} 1/p <= ln ln X + 0.30` and the `pi`
bound.

**Walk list `W`:** the first `L = 2^50` H-primes above `2^39`, `w_0 < ... < w_{L-1}`, even gaps `g_j`. **H2'** (checked
by the computation that builds `W`): `w_{L-1} < 2^77` and all `g_j <= G = 2^33` (under H2 one gap exceeds `2^33` with
probability `exp(-2^7.3) = 2^-227.35`, `G` being `2^7.3` mean gaps `d ln 2^77 = 2^25.7`; over `2^50` gaps `2^-177.35`,
both far below the smallest tail of the success line, `2^-36.27`).

**H2'' (availability): a premise of the success claim, not a consequence of charging the search.** Three facts about the
fixed modulus `r`, all decided before `t` and `tau` are drawn: (i) the declared walk range `(2^39, 2^77)` contains at
least `2^50` H-primes, so `W`, the first `2^50` of them, exists with `w_{L-1} < 2^77`; (ii) the declared base range
`(2^77, 2^96)` contains at least `2 n_b = 2^66.664` H-primes, so both base lists exist; (iii) every consecutive walk gap
satisfies `g_j <= G = 2^33`. Both searches are deterministic functions of `r`, so H2'' holds for this `r` for every input
and every coin choice, or fails for all of them: a 0/1 property of the modulus, not an event over the coins, which no
averaging removes. Its failure enters the success line of Section 6 as the factor `[H2'']`, not as a probability: if
H2'' fails the run stops at the search with no scalar and the success probability is 0.

*Support, strongest first.* (1) **The run decides it**: the charged scan of steps 1-2 enumerates every odd candidate of
both intervals, proves each survivor prime by complete trial division and applies one H-test to it, and the construction
proceeds on the H-primes it counts, so (i) and (ii) are tested at the price already charged (`2^79`, `2^95` units), and a
short list ends the run before a record is generated. (2) **No lower bound of the form "`2^50` H-primes exist below
`2^77`" is proved here, and a prime-counting argument cannot give one**: with `q < 2^77 < r` an H-prime is a prime value of
`g^(d k) mod r` over an interval of length `2^77`, where even the integer Polya-Vinogradov bound for the indicator of `H`
(`sqrt(r) ln r = 2^135`) is vacuous, and no prime-value analogue exists; an upper bound on `pi` divided by `d` (`2^53.26`
below `2^77`, `2^68.99` in the base range) predicts the count, it does not bound it. (3) **What is
measured**: the H2 density prediction is calibrated on this `r` by the two ratios above (0.983, 1.043 of `1/d`), and
`evidence.py availability` runs the same procedure at toy scale on real moduli — both lists, the counts against
`L + 2 n_b`, and the gap condition — finding `10.5x` to `56.7x` the needed count with the largest gap `2.06` to `6.48`
mean gaps (table in Section 6). Evidence at toy scope, which is what this track asks of an experiment: it supports the
procedure, the counting and the gap rule, not the production counts.

**Base primes:** H-primes in `(2^77,2^96)`: H2 predicts about `2^69.5` in the whole interval and `2^68.99` from the
`2^88.99` primes the cap holds; the search reports the ones it finds and the construction uses those, the parameter
needing `2 n_b ~ 2^66.664` of them, dealt alternately into a baby list `a_i` and a disjoint giant list
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
`x'^3 = y'^2 - 4 = x^3`). Different orbits agree on `kappa` only by chance (H3).

**Match.** A baby `[x]P0` and giant `[x'z]P0` in one orbit give `x = u x' z`, `u in mu6`, so `z = x (u x')^{-1}`. Try the six
`u`; accept the one with `[z]P0 = A_d` (at most 6 scalar multiplications of 383 reference operations: `<= 2^12` units). An
orbit match always yields the true `z`; an equal-`kappa` match of different orbits fails the test.

## 5. Machine model and matching

**Machine.** A shared-nothing cluster with three kinds of node, every count derived from the record capacity `C` of
Section 5.1: `S = 2^74` forwarding and owner nodes, `2^20` generating nodes and `2^31` storage nodes (for `W` and the
base lists). A forwarding node is a scalar AArch64 core with `2^50` bytes of local RAM, addressed by ordinary 64-bit
loads and stores; **no instruction addresses another node's memory**, so every byte that crosses a node boundary is a
message and is priced below. Nodes exchange blocks of `b = 25` records (1200 bytes); hardware moves the bytes, software
charges at most 60 instructions to post and 40 to poll a block. The storage nodes keep the smaller `2^40`-byte budget of
Section 5.1; only the forwarding nodes need the larger one, and memory is counted at full allocation, not scored.

**Addresses, and what a pass does.** A record's 74-bit owner id is its key bits 1..74, read as two 37-bit digits
`(a_1, a_2)`; the forwarding nodes are the same `2^74` pairs. Pass `k` is an exchange along axis `k`: a node holding
records of digit `d_k = v` sends them to the node that differs from it only in coordinate `k`, which is set to `v`.
After pass 2 the record with owner id `(a_1,a_2)` is on node `(a_1,a_2)`: pass 1 fixes `a_1` over the whole node set,
and every node of the axis-2 line `(a_1, . )` already carries `a_1` in coordinate 1, so pass 2 fixes `a_2` on it. Each
pass has **`2^37` destinations per node**, not `2^74`: a node's open-block area holds one 25-record block per
destination digit, `2^37 * 25` records `= 2^47.229` bytes, and that area together with the two record buffers is what
sets the node's `2^50`-byte RAM (Section 5.1). Fan-out, not arithmetic, is why two passes are needed, and it is why a
pass is a digit exchange rather than a direct scatter to the owner.

**Where a pass's records are before it runs.** A generating node holds no records. It walks chains and posts each
record as it is produced: it holds `2^46.664` chains of `2^50` records and streams `2N/2^20 = 2^96.6644` records in
total, so its store is the open-block area (`2^47.229` bytes) and nothing else. Its posting chooses the receivers of
pass 1: for a record of digit `a_1` the block goes to the receiver `(a_1, r_2)` whose `r_2` is the 37-bit value of a
per-digit posting counter, so successive blocks for one digit walk the `2^37` nodes of the slab `(a_1, . )` in a fixed
round robin. The sender never reads digit 2, so the assignment is independent of the key, and each generator's block
count per receiver is `floor` or `ceil` of its own count for that digit over the `2^37` nodes of the slab.

**Per-node inbound load at every stage, from `C`.** `2N = 2^116.6644` records. Section 5.1 fixes `C = 2^43.1827`
records in a `2^50`-byte node; hence a pass needs at least `2N/C = 2^73.4817` receiving nodes, and `S = 2^74` (a
product of two 37-bit digits) is the smallest admissible count. `mu = 2N/S = 2^42.6644` is the mean inbound load of
every routing pass and also the mean final bucket, and `C/mu = 1.4323`:

| stage | nodes | inbound records per node | where they are held while the next stage runs |
| --- | ---: | ---: | --- |
| generate (pass-1 senders) | `2^20` | `2^96.6644` streamed, none held | open-block area, `2^47.229` B |
| pass 1, digit 1 | `2^74` | round-robin share of its first-digit slab; `> C` with probability `<= 2^-2^73` (union bound over the `2^37` slabs) | pass buffer A |
| pass 2, digit 2, owner bucket | `2^74` | its owner bucket, mean `mu`; `> C` with probability `<= 2^-6.9e10` (union bound over the `2^74` nodes) | pass buffer B, while A is routed out |
| local sort and scan | `2^74` | the same bucket, `> C` with probability `<= 2^-6.9e10` (the same event) | sort pair A, B (count into B, scatter back) |
| storage of `W` and base lists | `2^31` | `2^39.664` B of 16-byte fields | local RAM, `<= 2^40` B |

**The stage loads are probabilistic, not deterministic; no all-pass inequality is claimed.** The capacity test of
Section 5.1 happens at run time, and the statement "every node in every pass is within `C`" is **not** made. What is
proved, under H3 of Section 6 (the orbit keys behave as independent uniform strings at the relevant width), is a union
bound for each pass separately.

*Pass 1.* A sender spreads a digit's records over the `2^37` nodes of the first-digit slab by a posting counter that
never reads the key, so a receiver `r` in slab `a` gets `floor` or `ceil` of each generator's `n(g,a)/2^37` records,
plus that generator's final marker block (at most `25` record slots, `2^20` generators):

```text
inbound(r in slab a)  <=  sum_g n(g,a)/2^37 + 2^20 + 25*2^20  =  N_a/2^37 + 26*2^20 ,
```

with `N_a = sum_g n(g,a)` the records whose first digit is `a`. This inequality is conditional on `N_a`; `N_a` itself
is a sum of independent `[0,12]`-bounded orbit contributions of mean `mu*2^37 = 2^79.6644`, so by Chernoff with
threshold factor `a_1 = (C - 26*2^20)/mu = 1.4323`

```text
P(overflow in one slab)  <= exp(-2^79.6644 (a_1 ln a_1 - a_1 + 1)/12) <= 2^-2^73.005 ,
```

and the union over the `2^37` first-digit slabs gives `2^37 * 2^-2^73.005 <= 2^-2^73`. That union also covers all
`2^74` receiving nodes, because the event that any receiver of a slab exceeds `C` is contained in the slab event (the
round-robin ceiling is exact given `N_a`). Nothing here is a deterministic bound on the pass: `N_a` is a random count
and the probability is over the keys, i.e. over the drawn coins.

*Pass 2.* A record's second digit *is* the coordinate of its destination, so pass 2 has no round-robin split and a
final marker block carries no records. A receiver's inbound is exactly its owner bucket, mean `mu = 2^42.6644`, a sum
of independent `[0,12]`-bounded orbit contributions, and with the same factor `a = C/mu = 1.4323`

```text
P(overflow in one bucket)  <= exp(-mu (a ln a - a + 1)/12) <= 2^-6.897e10 ,
```

the union over the `2^74` nodes giving `2^74 * 2^-6.897e10 = 2^-6.897e10` (the exponent loses only `74`). The local
sort and scan read the same bucket, so they add no new event.

*Deduction.* The two events are charged in the success line of Section 6 as `2^-6.9e10 + 2^-2^73`, both far below the
`2^-36.27` term, so the tail total stays `1.2e-11` and success `0.3954794`; the `0.0055` margin is untouched. **The
deterministic reading of the stage table is withdrawn: no pass's inbound is bounded by a slab count divided by a node
count, and no statement of the form "for all passes, inbound `<= C`" appears in this proof.** The stage table's
inbound column is what the capacity test compares against, with the event's bound stated beside it, not a proved
per-pass bound.

Both record buffers are `C` records; a node reads one while it fills the other, and at the end the pair is the sort
buffer. No stage of the run asks a node to hold more than `C` records except on the `2^-6.9e10` (bucket) and
`2^-2^73` (pass 1 slab) events charged in Section 6, and the capacity test is recomputed from `C` at each stage
instead of being read off the final owner load.

**Record.** 48 bytes: `kappa` 28 B, label 16 B (67-bit base index over the `2 n_b` base primes, 50-bit step, type
bit), 4 B padding. The label holds **indexes, not values**, so the record size does not depend on the prime width.

**Stored primes: one 16-byte field, and the read path.** `W` and the base lists are arrays of **16-byte fields**: bits
0..96 hold the prime (a base prime is below `2^96`, the capped scan of Section 3 ending below `2^95`; a `W` value below
`2^77` uses the same field), bits 112..127 hold a 16-bit tag (kind: walk, baby base, giant base, with the parity bit
used as a check), bits 97..111 are zero. A 96-bit value therefore needs no truncation and loses no tag bits, which a
10-byte field could not do. The arrays live on the `2^31` storage nodes (`2^54` bytes for `W`, `2^70.664` bytes for the
base lists, `2^39.664` bytes per storage node), together with the `2^43`-entry sieve table of Section 3 (`2^46.4` bytes
in all, `2^15.4` bytes each).
A label is resolved by a message, not by a remote load: the owner sends a request naming the index and receives the
16-byte field back. That read happens **only inside the verification cap** (two fields per verified pair, two
round trips of about 250 instructions each, against the `2^12` units = `8.192*10^6` instructions charged per
verification), and once per chain for the base prime `a`. The walk stream reads `W` sequentially, one field per
record, and needs no prefix sum: the chain accumulates the gap `w_{j+1} - w_j` from the two fields it already holds.

**Routing pass** (per record, on one 37-bit digit). Every destination has its own open block at a fixed address, so the
loop keeps the record's destination computed and written back:

```text
d  = (kappa >> 37*k) & (2^37 - 1)      extr, ubfx                         2
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
once per pass (node setup row), and the address needs no carry instruction because `b0 <= 2^47.23` and
`a = b0 + 48*24 < 2^47.3 < 2^64`. The post is counted per block (Section 5.2), the transport copy is the receiver's copy
out of the ring into its pass buffer (Section 5.3), and the credit wait is inside the 60 of the post. The fill-count
writeback is the `str` above. Pass 1 uses the same loop; its extra work is the receiver choice of the last paragraph,
two instructions per block inside the 60 of the post.
**Owner sort drain** (per record, after the second pass): the poll of every block arriving in the second pass returns the
ring slot and the record count `n`; the records have been copied into the pass buffer by the second pass's transport
copy, and the loop then runs `n` times: capacity test of the cursor against the last record slot of the `C`-record sort
buffer, before any write (`cmp`, `b.hi` to *overflow*: 2), 3 `ldp` from the pass buffer and 3 `stp` to the sort buffer
at the cursor, both post-indexed so the source pointer and the cursor advance by 48 (6), loop (`subs`, `b.ne`: 2):
charged **10**. The local sort reads only the sort buffer.
**Local LSD pass**: counting loop (`ldr`, `extr`, `ubfx`, `ldr/add/str`, loop: 9) + scatter (`extr`, `ubfx`, `ldr`
offset, `madd`, `add`, 6 copies, `add/str`, loop: 16) = 25: charged **25**.

**Matching.** Each record is routed by 2 routing passes on the first 74 bits of `kappa` to its owner (the `mu` of the
Section 5.1 table, one 74-bit prefix), sorted by 2 LSD passes on `kappa[75..142]` (34-bit digits) and scanned once.
Records with equal 142-bit prefix are adjacent (a run); inside a run all opposite-type pairs are compared on the full
224-bit `kappa`, and pairs with equal `kappa` are verified by the Match step. **Caps, per node:** every owner keeps two
64-bit counters, run-pair comparisons and verifications; a node that would make its `2^28 + 1`-th comparison or its
`129`-th verification stops and signals *fail* (the counters never exceed `2^28 + 1` and `129`, far below `2^64`). No
global counter and no combining step exist: each cap is local, and the failure of any node is a failure of the run.
They are charged at the caps on every owner: `2^74 * 2^28` comparisons of 32 instructions (`2^96.03` units) and
`2^74 * 2^7` verifications of `2^12` units (`2^93` units). A comparison, itemised: pair loop (index `add`, `cmp`,
`b`: 3), partner record address (1), partner type bit (`ldrb`, `tbz`: 2), both 28-byte keys (`ldp`, `ldr`, `ldr w`
each: 6), 224-bit compare (4 `eor`, 3 `orr`, `cbnz`: 8), counter (`add`, `cmp`, `b.hi`: 3): 23, charged 32. A
verification resolves its two labels to 16-byte fields by message (Section 5, read path), so the two fetches are inside
its `2^12` units. All other work is fixed: the number of records and passes does not depend on the data.

### 5.1 Placement, capacity and buffers

`2N = 2^116.6644` records. Every node that forwards records, owners and generating nodes alike, holds the
same four record stores, all concurrently resident, each with one lifetime:
1. the **open-block area**, `2^37` blocks of `b = 25` records = `2^47.229` bytes, one block per destination at the
   fixed address `1200 d`; it holds this node's outgoing partial blocks and is reused as each block is posted;
2. the **pass buffer**, `C` records (48 bytes each), holding the inbound of the pass being received;
3. the **second record buffer**, `C` records, holding the inbound of the pass being routed out (the pair alternates
   roles per pass and, at the end, is the sort buffer: counting pass into one, scatter back into the other);
4. the **receive ring**, `2^20` slots of one block (`2^31.6` bytes), holding blocks in flight.
Rounding out the allocation: LSD count and offset arrays `2 * 2^34 * 8 = 2^38`; routing offset and cursor tables
`2 * 2^37 * 8 = 2^41`; code and stack `<= 2^30`; reserved `2^36`. The auxiliary total is
`167473659772928 = 2^47.251` bytes, so `C = (2^50 - 1.6747e14)/96 = 9983606740309 = 2^43.1827` records and
`96 C + aux = 2^50 - 32` (the floor in `C` leaves 32 bytes unused; the two record buffers are `96 C = 2^49.767`).
**`C` is what fixes the node count**: a
pass must receive `2N` records somewhere, so it needs at least `2N/C = 2^73.4817` receiving nodes, and `S = 2^74` — a
product of two 37-bit digits, the smallest such count above that bound — gives `mu = 2N/S = 2^42.6644 = 6.97 * 10^12`
records per node and `C/mu = 1.4323`. See the stage table of Section 5 for the load at every stage; a node whose
incoming pass would exceed `C` signals *overflow* and the attack fails (the drain's capacity test; the cursor is a byte
address below `2^50`, the count at most `C + 1 < 2^44`, both 64-bit). By the per-pass union bounds of
Section 6 that costs at most `2^-2^73` for pass 1 (union over the `2^37` first-digit slabs) and `2^-6.9e10` for pass 2
and the final bucket (union over the `2^74` nodes). This is a probabilistic statement under H3, not a deterministic
bound on any pass.

A generating node holds a chain table `2^32 * 96 B = 2^38.585`, the same open-block area `2^47.229`, stream buffers
`<= 2^30` and code `<= 2^30`: `<= 2^47.3 < 2^50` bytes, and no record store at all, because it posts every record as
it is produced. A storage node is unchanged: `2^40` bytes of local RAM holding `2^39.664` bytes of 16-byte fields for
its share of `W` and the base lists.

### 5.2 Blocks sent per pass, summed over the actual nodes

**The posting charge is one global total over the run, not a per-node bound.** A node routes on one 37-bit digit, so
it has `2^37` destinations, while a generating node's pass-1 receivers are the `2^74` owner nodes (Section 5), and
nothing below assumes otherwise. Let `n_{v,d}` be the record count of node `v` for destination `d` in one pass and
`M_v = sum_d n_{v,d}`, which is `2^96.6644` at a generating node and `mu` at an owner. Each destination receives
`ceil(n_{v,d}/25)` blocks: its full blocks plus the one *final* block holding its `0..24` remaining records (sent
whether or not the destination is empty, since it is also the end-of-pass marker of Section 5.3). Summed over every
node of one pass, with `sum_v M_v = 2N` because every record of the attack is in exactly one node during a pass (a
generating node's records are produced into its own open blocks and counted in its own `M_v`),

```text
full blocks      = sum_v sum_d floor(n_{v,d}/25) <= sum_v M_v/25 = 2N/25
final blocks     = sum over forwarding nodes of 2^37 per pass, plus 2^20 x 2^74 for pass 1
                   (its sender set is the generating nodes, whose marker goes to every owner)
total posting    = 2 x (2N/25) x 100 + 214 x (2^74 + 2^20 + 2^31) x 2^37 + 107 x 2^20 x 2^74
                 = 2^119.6644 + 2^118.7415 + 2^100.7415 = 2^120.2755 instructions  (= 2^109.3097 units) .
```

`#nodes = 2^74 + 2^20 + 2^31` is counted at every node that can hold a store, generating and storage nodes included,
each charged a full owner's rows, so the `2^37`-per-pass final-block term is an over-charge for the nodes that send no
records; the third term is pass 1's own marker, one block per (generating node, owner) pair. The ledger carries
exactly these three terms and nothing else, under the names `100/25 = 4` instructions of full-block posting per record
per pass (60 to post plus 40 to poll a 25-record block, inside the routing row of Section 7.1, `2 x 4 = 8` per
record), the final-block row `2 x 2^37 x 107` instructions per node and the `pass1_end_tokens` row
`107 x 2^20 x 2^74` (Section 7.2); `evidence.py` asserts that the three rows equal this total. Since `48 x 25 = 1200`
bytes and a final block is at most 24 records, no block is ever double counted. No statement of the form "a node
handles at most `mu` records" is used: the per-node inequality `sum_d ceil(n/25) <= M_v/25 + 2^37` holds for every
`M_v`, and the global rows above depend only on `sum_v M_v = 2N`, on the node count and on pass 1's sender set. `mu`
enters only the capacity of Section 5.1 and the per-stage load table of Section 5.

### 5.3 Communication schedule: flow control, buffers, termination

**Flow control.** A ring slot is a credit. A sender posts a block only against a credit it has been granted for that
destination, and the receiver's poll returns the slot as a credit once the block's records have been copied out. Post
and poll are charged 60 and 40 instructions per block (Section 5.2), which cover the credit test, the header write, the
slot release and the wait.

**Why no sender waits forever.** Take one pass and a node `w` receiving it. The buffer `w` fills for that pass is empty
when the pass starts, its capacity is `C` records, and the whole inbound of the pass is at most `C` records with
probability at least `1 - 2^-2^73` for pass 1 (its union bound over the `2^37` first-digit slabs in Section 5, its
`26*2^20` of block rounding included) and at least `1 - 2^-6.9e10` for pass 2 (the bucket union bound over the `2^74`
nodes; a node whose inbound would exceed `C` signals *overflow* and the
run fails, and that event is one of the terms of the success line of Section 6). So `w`'s copy-out of a ring slot
(6 instructions per record, the transport copy of the routing row) is never blocked by storage while the pass runs,
and a slot is returned after at most 25 records of work, whatever `w`'s own sending state is. Sending and receiving use
disjoint stores: the open blocks and the buffer being routed out
are `w`'s, the inbound lands in the other record buffer. A node that has filled a destination's block posts it and, if
no credit is free, waits; meanwhile its receiving continues. So blocked senders wait only on receivers, no receiver's
copy-out ever waits on a sender, and the wait-for graph has no cycle: the schedule is deadlock free without any global
order, barrier or priority rule. A node has at most one block outstanding at a time (the routing loop posts a block
before refilling it), so pass 1's `2^20` senders fit the ring exactly; in pass 2 the sender set is `2^37`, so at most
`2^37 - 2^20` senders wait for a credit at a node, each granted one within one block's copy-out.

**Termination and detection.** Each node keeps three 64-bit counters (charged in `node_setup`). The pass counter counts
the final blocks it has received in the current pass. The senders of a pass are fixed by the digit route: the `2^37`
nodes that agree with the receiver in the other digit, except in pass 1, whose senders are the `2^20` generating
nodes, each of which sends one final block to every one of the `2^74` owners (this is the `pass1_end_tokens` row). The
count is exact because every sender sends exactly one final block per destination, so a node knows the pass is complete
when the count reaches its pass's sender count, `2^20` then `2^37`. A node then clears its counters and starts
the next pass; the completion of pass 2 is the end of the routing, after which each owner sorts and scans. `fill[d]`
and `off[d]` are reset per pass by the routing offset init row of Section 7.2, and the open blocks are reused in place,
so there is no allocation during a pass.

**Peak memory** is the four stores above at one forwarding node, `96 C + aux = 2^50` bytes exactly (Section 5.1); the
credits are the ring slots themselves and need no extra store.

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
is the ensemble H3 below is stated on.

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

**H3 (keys, on the coins).** For fixed inputs and the coins `(tau, z)`, the records' 224-bit keys and their 68-bit
windows (bits 75..142) behave like independent uniform strings at the relevant width, jointly over the records of a
node. This is the single assumption behind both the per-node caps below and the per-pass routing loads of Section 5:
it makes each record's routed digits uniform and independent across orbits, so a receiver's load is a sum of
independent `[0,12]`-bounded orbit contributions. Support: Lemma 4 is
exact on orbits, so only different orbits can collide; the toy `ec` run bins keys into 32 pseudo-nodes and measures
chi-square `z`-scores in `[-0.79, +0.98]`; and the `coverage`/`ec` agreement shows no sign of structure. Because every
key is `kappa` of `[tau x]P0` or `[tau x' z]P0` with `tau` uniform, the ensemble is over the algorithm's coins, as the
review standard requires, and the baby records are not held fixed.
**Runs and orbits.** An orbit carries at most 12 records: its 6 elements of `H`, each at most once per side (Lemma
2(ii), and distinct `x'` give distinct `x' z`), so at most 6 baby and 6 giant records, and two orbits give at most
`6*6 + 6*6 = 72` opposite-type record pairs.
**Overflow of a bucket or of a pass: a union bound under H3, never a deterministic statement.** Orbits fall
independently and uniformly on the `2^74` owners (key bits 1..74), each adding 1 to 12 records, so at any receiver the
load is a sum of independent `[0,12]`-bounded contributions of mean `mu = 2^42.6644`, and the Chernoff bound for such a
sum gives `P(load >= C) <= exp(-mu (a ln a - a + 1)/12)`, `a = C/mu = 1.4323`, i.e. `2^-6.897e10` per node. This is
pass 2's inbound — its destination is fixed by the second digit, so its inbound is exactly its owner bucket, with no
round-robin split and no marker record slot — and also the final bucket the sort reads. The union over the `2^74`
nodes is `2^74 * 2^-6.897e10 = 2^-6.897e10` (the exponent loses only `74`; `evidence.py` prints all three numbers).
A **pass-1** receiver is a digit exchange on uniform keys (H3), so the records it must hold are those of its own
**first-digit slab**: `N_a` for the slab `a`, whose mean is `2N/2^37 = mu*2^37 = 2^79.6644` and which is a sum of the
same `[0,12]`-bounded orbit contributions. With the `26*2^20` of round-robin and marker rounding of Section 5, the
receiver overflows only if `N_a > (C - 26*2^20) * 2^37`, a factor `a = 1.4323` above the mean, so

```text
P(overflow in one slab)  <= exp(-2^79.6644 (a ln a - a + 1)/12) <= 2^-2^73.005 ,   a = 1.4323 ,
```

and the union over the `2^37` first-digit slabs — which contains the union over the `2^74` receiving nodes — is
`2^-2^73`. This is the only pass whose inbound is argued through a slab count, and the argument is probabilistic in
`N_a`. **No pass's inbound is deterministically bounded by a slab count divided by a node count, and no "for all
passes" overflow statement is made:** pass 2 is an owner bucket and needs no slab argument, and the final bucket is the
same event as pass 2. The success line below charges `2^-6.9e10 + 2^-2^73` for these events.
**Per-node caps.** Fix an owner holding `M <= C` records (else it has already failed by overflow), hence `M' <= C`
orbits. Given which orbits it holds (key bits 1..74), their key bits 75..142 are independent and uniform on `2^68`
values. A compared pair is an opposite-type pair in one run; a same-orbit one has equal full keys, is verified, yields
`z` and ends the attack, so at most one is compared, and every other compared pair joins two different orbits with
equal bits 75..142.
(A) Some value is held by at least 7 orbits: probability `<= C(M',7) * 2^(-68*6) <= C^7/(7! * 2^408) = 2^-118.02` per
node, `2^-44.02` over `2^74` nodes.
(B) Let `D'` be the number of orbits whose value repeats that of an earlier orbit (any fixed order). Each orbit repeats
with probability at most `M'/2^68` whatever the earlier ones did, so `D'` is dominated by `Bin(M', M'/2^68)`, mean
`<= C^2/2^68 = 337703`. Outside (A) every value holds `j <= 6` orbits and `C(j,2) <= 3(j-1)`, so the orbit pairs with
equal value are at most `3 D'` and the compared pairs at most `1 + 72 * 3 D' = 1 + 216 D'`. The cap is `3.68` times the
mean of that quantity, `1 + 216 * 337703 = 2^26.12`, i.e. `2^28`, so exceeding it needs `D' >= 1242757`, and `P(Bin >= k) <= (e * 337703/k)^k` gives
`2^-543111` per node, `2^-543037` over all nodes.
(V) Verifications at a node are at most `1 + 72 * Q`, `Q` the number of orbit pairs with equal 224-bit keys. The cap
128 is exceeded only if `Q >= 2`, i.e. two distinct orbit pairs each agreeing on the 150 bits 75..224 (disjoint or
sharing one orbit, probability `2^-280` either way): `<= C(C(M',2),2) * 2^-280 <= C^4/8 * 2^-280 = 2^-110.27` per node,
`2^-36.27` over all nodes.

**Success.** With `[H2'']` the indicator of the availability premise of Section 3 (1 for this modulus or 0, the same for
every input and every coin choice),

```text
P(success) >= [H2''] * 0.4054794 - 0.01 - 2^-44.02 - 2^-543037 - 2^-36.27 - 2^-6.9e10 - 2^-2^73 = 0.3954794 ,
```

so `0.39` holds for every input **when H2'' holds for this `r`** (`0.0055` below the bound), and a failure of H2'' gives
0: the run ends at the search before any record is generated. The tails total `1.2e-11` (`2^-36.27` dominates; `2^-6.9e10`
and `2^-2^73` are far below it). Primality contributes no term: every accepted base value is proven prime by complete
trial division, and both scans are deterministic in `r` and charged in full. Availability is the one premise of the line
that is not an event over the coins; it is neither averaged over (no `0.39 - P(short list)` term) nor retried into. The
model allowance is the one flat term and H2'' the one modulus-dependent premise, named and tested.

**Toy checks** (`evidence.py`, `mode all`, seed 2026; `8000/5000/3000` trials). Scalar model, three toy primes with `6 * 2^k | r-1`, chain length `L = 60`: measured success for
`r ~ 2^27.6, 2^29.7, 2^31.7` is `0.4255, 0.4242, 0.3850` against the diagonal-corrected Poisson `0.42493, 0.41923,
0.38940` (`|z| <= 0.712246`); at toy scale the diagonal pairs are a `1/L` fraction, so the prediction uses
`K' = K(1 - (L-1)/L^2)`. The `ec` mode runs Algorithm A end to end on genuine `j = 0` curves of prime order
`r ~ 2^28..2^32`: affine incomplete additions with Montgomery batch inversion (asserting no zero denominator), the
parity key truncated to 18/20/22 bits so that false key matches occur (0.5-1.8 per run) and are rejected under a cap,
sort, equal-key runs, `z`, and the `d`-th root **selected by matching `T`** (the `d` roots are scanned with one scalar
multiplication each and the one whose point is `T` is taken; the recovered scalar is then compared with the secret):
`28/14/8` recoveries in `80/30/20` trials, all selected from `T`, success agreeing with the scalar model in every trial,
chi-square `z`-scores in `[-0.79, +0.98]`. Limits: tables are built only for occurring gaps (`tau[g] = tau[g-2] +
tau[2]`, the segments resting on Lemma 3), and the toy sorts locally without distributed routing. Evidence at toy scale
only: it supports the Poisson form, the duplicate bookkeeping and the root selection, not the production value of `m`.

**Toy availability (H2'' of Section 3).** `evidence.py availability` (also in `mode all`) builds the analogue of both
charged lists with the toy parameters and the same procedure — one segmented scan crossing out the multiples of every
prime below the square root of the range top (so each survivor is prime by complete trial division), then one H-test per
survivor — takes the first `L` H-primes as the walk list and the next `2 n_b` as the base lists, re-checks every accepted
prime by trial division, asserts the lists disjoint, and checks every consecutive walk gap against `G` scaled to the toy
by production's own multiple (`G = 2^33` is `2^7.3` mean gaps, a mean gap being `d ln q`).

| toy | declared range | odd candidates | primes | H-primes | needed | headroom | max gap | max gap / mean |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| t28 | `2^11`..`2^21` | 1047552 | 155302 | 4875 | 86 | 56.7x | 958 | 2.06 |
| t30 | `2^12`..`2^21` | 1046528 | 155047 | 2395 | 98 | 24.4x | 3732 | 4.01 |
| t32 | `2^14`..`2^21` | 1040384 | 153711 | 1160 | 110 | 10.5x | 12070 | 6.48 |

Evidence at toy scope only, a real modulus and procedure at toy parameters. The toy's window is wider in relative terms
than production's, so its headroom exceeds the `7x` and `5x` predicted at production size; it establishes that the
procedure returns the needed counts with gaps inside the multiple of the mean gap that `G` allows, not the production
values of `2^50` and `2^66.664`.

## 7. Cost ledger (unit = 2000 instructions)

### 7.1 Per record, on top of the reference addition

| item | instructions |
| --- | ---: |
| walk bookkeeping, 30 itemised: gap load and arithmetic 3, absolute-prime stream transport (16 B per 8 = 2 loads) 2, stream block post/poll amortised over a 128-field block 1, table address `base + (gap/2)*96` 3, segment loop 3, label 3, batch pointer/index updates outside the Montgomery-trick multiplications 12, stream-refill check 3 | 30 |
| key, 23 itemised: `ldp` of 6 `y` limbs 3, parity test and select mask 2, `p - y` by `subs`, 5 `sbcs`, 6 `csel` 12, `extr` to bits 1..224 4, `stp` of the key 2 | 23 |
| 2 routing passes x 29 (19 control: destination `extr`/`ubfx`, offset `ldr`, fill-count `ldr`, block-base `add`, address `madd`, six loads/stores, fill-count `str`, full test, loop; 4 for posting a full 25-record block, 6 transport copy) | 58 |
| owner sort drain, pass buffer to the sort buffer (3 `ldp` + 3 `stp` 6, capacity test 2, loop 2) | 10 |
| 2 local LSD passes x 25 (counting 9 + scatter 16) | 50 |
| run scan (load the 142-bit sort prefix, compare with the previous record, branch) | 10 |
| **itemised; charged** | **181** |

Element cost `c = 1 + 181/2000 = 1.0905`.

**Posting is one global total.** The `4` in each routing pass is the full-block share `100/25` of the global count
`sum over nodes of sum_d floor(n_{v,d}/25) <= 2N/25` per pass, and the `2^37` final blocks per node per pass are not per
record but the final-block row of Section 7.2. Summed over the whole run and over every node, the posting charge is
`8 x 2^116.6644 + 214 x (2^74 + 2^20 + 2^31) x 2^37 + 107 x 2^20 x 2^74 = 2^120.2755` instructions, exactly the
sum of the share inside this row, the final-block row and the `pass1_end_tokens` row (Section 5.2 carries the
derivation and the identity). No per-node premise of the form "a node handles at most `mu` records" is used anywhere:
the per-node inequality used is `sum_d ceil(n_{v,d}/25) <= M_v/25 + 2^37`, which holds for every `M_v` (a generating
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
| elements `2Nc` | `2^116.7893` |
| tables: `2 n_b` chains x `2^32` entries x (1 addition + 12 stores/2000) | `2^98.673` |
| segment and table-block starts, `[tau a]P0` (`<= 2^25` per chain, `2^12` for the two base scales) | `2^91.665` |
| per-node fixed work, `(2^74 + 2^20 + 2^31) nodes x 30571581411328 = 2^44.797` instructions, itemised below | `2^107.8315` |
| capped H-prime search in `(2^77,2^96)`: `2^94` odd candidates, sieve by every prime `< 2^48` plus one H-test per prime, itemised `2^93.90`, charged | `2^95` |
| capped H-prime search for `W` above `2^39`: `2^78` odd candidates, sieve by every prime `< 2^40` plus one H-test per prime, itemised `2^78.17`, charged | `2^79` |
| coins, input access, Lemma 1 (square roots and the `2^20`-unit root match), `t` setup | `2^22` |
| verification cap, per node: `2^74` owners x `2^7` verifications x `2^12` units (the 16-byte label fetches are inside it) | `2^93` |
| run-pair cap, per node: `2^74` owners x `2^28` comparisons x 32/2000 | `2^96.034` |
| pass-1 end-of-pass markers, one per (generating node, owner): `107 x 2^20 x 2^74` instructions | `2^89.7757` |
| walk stream set-up and base-prime fetches: `2^66.664` chains x (`2^16` segment streams + 1) x ~100 instructions | `2^79` |

**Per-node fixed rows** (instructions per node; every node, generating and storage nodes included, is charged the
owner's rows). Every array is 8-byte entries, `2^37` for the routing tables and `2^34` for the local sort; each loop
handles 16 entries per iteration and pays `subs` + `b.ne` (2) per iteration, and both sizes are multiples of 16, so
there is no remainder.

| array and pass | per 16 entries | entries x passes | instructions |
| --- | ---: | ---: | ---: |
| routing fill-count table, clear (8 `stp xzr, xzr` + 2) | 10 | `2^37 x 2` | 171798691840 |
| routing offset table, `offset[v] = 1200 v` (per 2 entries `stp`, 2 `add`: 24 + 2) | 26 | `2^37 x 2` | 446676598784 |
| final-block scan (per 2 entries `ldp`, 2 `cbnz`: 24 + 2) | 26 | `2^37 x 2` | 446676598784 |
| final block per destination: call and return 2, block address 1, count header 1, post 60, poll 40, final-block counter `add`, `cmp`, `b.eq` 3 | `107 x 16` | `2^37 x 2` | 29411936043008 |
| LSD count array, clear (8 `stp xzr, xzr` + 2) | 10 | `2^34 x 2` | 21474836480 |
| LSD exclusive prefix sum into the offset array (per 2 entries `ldp`, `add`, `stp`, `add`: 32 + 2) | 34 | `2^34 x 2` | 73014444032 |
| receive ring, `2^20` slot descriptors (`str`, `add`, `subs`, `b.ne` per slot) | - | `2^20` slots | 4194304 |
| node setup: code, stack, cursors and the three 64-bit counters | - | once | 4096 |
| **per node** | | | **30571581411328** |

The 8-byte entries are wide enough everywhere: a fill count is at most 25, an offset below `2^48`, an LSD count or
offset at most `C < 2^44`, a final-block counter at most the number of senders of a pass, at most `2^37`.

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
spread over the same nodes at `2^15.4` bytes each.

All rows other than `2Nc` sum to `2^107.8347`; `W = 2^116.7893 + 2^107.8347 = 2^116.792250 <= 2^116.7923`, score
`1167923/10000` (`evidence.py ledger` recomputes every row with exact integers and finds the smallest integer `k` with
`W^10000 <= 2^k` by integer exponentiation, `k = 1167923`; the rounding slack is `5.0e-5` bits).
Preprocessing (lists, tables, baby generation) `<= 2^115.796 <= 2^115.80`.

**The scanned layout table.** `evidence.py ledger` recomputes each row below from the same formulas the ledger uses:
per-record row from the pass widths, auxiliary reservation from the digit widths, `C` from that and the declared node
RAM, node count from its own digit split, per-node and posting rows from the same counts, total. `C` is the capacity the
node RAM leaves, `floor((RAM - aux)/96)` (a record occupies one 48-byte slot in each of the node's two `C`-record
buffers); `records/node` is the row's mean load `2N/S`; `C/mu` is the load factor. Every figure of the form
`2^x records` in this file, in the note and in the claim is one of these two numbers *of its own row*, and
`evidence.py quotes` fails on any other.
Each total is its exact value rounded to 4 decimals: the shipped row's exact `W = 2^116.792250` shows as
`2^116.7922` (row 3), and the claim `1167923/10000 = 2^116.7923` is its ceiling at that precision, the inequality above.

| passes | route x width | local LSD | per-record row | node RAM | `S` | `C` | records/node | `C/mu` | window | fixed rows | `time_log2` |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 3 x 28 | 2 x 28 | 210 | `2^40` | `2^84` | `2^32.7549` | `2^32.6644` | 1.06 | 56 | `2^109.436` | `116.8172` |
| 4 | 2 x 37 | 2 x 28 | 181 | `2^50` | `2^74` | `2^43.1831` | `2^42.6644` | 1.43 | 56 | `2^107.827` | `116.7956` |
| **4** | **2 x 37** | **2 x 34** | **181** | **`2^50`** | **`2^74`** | **`2^43.1827`** | **`2^42.6644`** | **1.43** | **68** | **`2^107.831`** | **`116.7922`** |
| 4 | 2 x 37 | 2 x 38 | 181 | `2^50` | `2^74` | `2^43.1765` | `2^42.6644` | 1.43 | 76 | `2^107.897` | `116.7924` |
| 4 | 2 x 38 | 2 x 34 | 181 | `2^50` | `2^76` | `2^42.9063` | `2^40.6644` | 4.730200 | 68 | `2^110.829` | `116.8123` |
| 4 | 2 x 36 | 2 x 34 | 181 | `2^50` | `2^72` | `2^43.3033` | `2^44.6644` | < 1 | 68 | - | infeasible: `2N/C = 2^73.36` |
| 4 | 2 x 30 | 2 x 43 | 181 | `2^64` | `2^60` | `2^57.4150` | `2^56.6644` | 1.68 | 86 | `2^94.501` | `116.7894` |
| 4 | 3 x 28 | 1 x 56 | 185 | `2^59` | `2^84` | - | `2^32.6644` | < 1 | 56 | - | infeasible: area `2^59` |
| 4 | 1 x 29 | 3 x 47 | 177 | `2^95` | `2^29` | `2^88.4150` | `2^87.6644` | 1.68 | 141 | `2^70.401` | `116.7867` |
| 6 | 4 x 21 | 2 x 28 | 239 | `2^40` | `2^84` | `2^33.3063` | `2^32.6644` | 1.56 | 56 | `2^104.201` | `116.8276` |

Four passes are cheaper than five by construction: one routing pass is 29 instructions and one local LSD pass 25, so a
4-pass row costs 181 against 210 (0.019 bits of `2Nc`). Three routing
passes plus **one** local pass fails on memory, not on instructions: a single local pass must cover the whole window,
so its histogram has `2^window` entries, `2^59` bytes at a 56-bit window, and the open-block area of a 28-bit pass
plus that histogram leaves no `C` at all (row 8). Two routing passes plus two local passes keeps the `2^34` histogram
affordable, and it is the shape that wins.

Two 37-bit digits are the largest that fit: the final-block row is `passes x 2^w x 107 x #nodes`, so `k` passes
on `w` bits cost about `2^(w(k+1))` instructions in total; `w = 38` with two passes is `2^114` against `2^111` at
`w = 37` (`+0.020` bits, row 5). `w = 36` is cheaper in fixed rows but its 72 owner bits do
not cover `2N/C = 2^73.36`, so the capacity rule refuses it (row 6). So `2 x 37` is the cheapest admissible pair of routing digits, and the node RAM that
admits it is the smallest power of two that holds the `2^47.229`-byte open-block area and the two `C`-record buffers:
`2^50`.

The local width is the other trade: the window `m*l` keeps `216 C^2 / 2^(m*l)` under the comparison cap, while the
histogram costs `m * 2^l` entries per node. `2 x 28` needs a `2^40` cap and pays `0.0033` bits for it (`2^116.7956`
against `2^116.7922`, row 2). At `2 x 34` the window is 68 bits, the expected comparisons per node
`216 * C^2/2^68 = 2^26.1`, so a `2^28` cap is `3.68` times that mean, the margin its Chernoff bound needs, and the
histogram pair is `2^38` bytes per node (two `2^34`-entry arrays of 8-byte counters; row 3, chosen). Wider local digits
(`2 x 38`, row 4) cost more in the histogram than they save in the cap.

**A bigger node does not buy a better total.** Rows 7 and 9 show the flat part of the surface: at `2^64`-byte nodes
(row 7) the fixed rows fall to `2^94.501` and the total reaches the floor `2^116.7894` set by `2Nc` alone; the `2^60`
nodes of that row carry `2^56.6644 records on average, 2^57.4150 of capacity`, and a 48-byte record leaves at most
`2^58.4150` in a `2^64`-byte node, so no node of row 7 holds the `2^87` records of row 9. At `2^95`-byte nodes with one
routing pass (row 9) the total is `2^116.7867`, and its `2^29` nodes carry
`2^87.6644 records on average, 2^88.4150 of capacity`. Both are worse as machines than the `2^50` node: a `2^95`-byte node sorts its records
in one place, the flat-memory model, and the `0.0029` bits row 7 gains are below the model allowance the claim already
carries; the chosen layout keeps a routed cluster of `2^74` nodes of `2^50` bytes and pays `0.0029` bits for it.

The block size stays `b = 25`: `mu/b = 2^37.664`, the open-block area is `2^47.229` bytes, and the two `C`-record
buffers fit `2^50`.

Memory counted at full allocation: `(2^74 + 2^20) * 2^50 + 2^31 * 2^40 = 2^124.0`, claimed `<= 2^124.1` (records alone
`2N * 48 B * 2 = 2^123.25`, which fixes the floor; `W` `2^54` and base lists `2^70.664` bytes fit the `2^31` storage
nodes at `2^39.664` each in `2^40`). Memory has no scalar weight and is exponential,
as for any BSGS-type Cheon attack.

## 8. Heuristics, evidence, provenance

H1 decides the size of the success set (with the flat model allowance of Section 6); H3 decides the failure
probability of the caps and of the per-pass routing loads, on the algorithm's own coins; H2/H2' size preprocessing, the
table and the search intervals (checked by building `W` and by the two measured density ratios) and are predictions of a
count, never lower bounds on it; H2'' is the availability premise the success line multiplies by `[H2'']`, a 0/1 property
of the fixed modulus, decided by the run's own charged scan and measured at toy scale by `availability`, and H2 is not
used in the success line. The cluster model is idealised and
the instruction counts are listed, not measured on hardware. The cost is data independent: the two prime caps are
deterministic functions of `r`, charged at their limits in every run, and a pass's capacity enters only the success
line, through the union bounds of Section 6.

`evidence.py` (one initial run `mode: all`, about 10 s, which also runs the `quotes` check; each mode can be rerun with
more trials): `ledger` (exact integers: Sections 5.1-7, the digit and block widths, the record and 16-byte prime-field
layouts, the derived node count `2N/C` and the per-stage load table (`stages`), the global posting total `2^120.2755`
instructions, its identity with the three rows that carry it and pass 1's `2^20 x 2^74` markers, every per-node row, the
scanned layout table (`layout_scan`: passes, digit widths, node RAM, node count, `C`, mean records per node, window, fixed
rows and total for each shape, including the infeasible ones and why) with `assert_layouts`, which checks for every layout
that `C = floor((RAM - aux)/96)`, that `C <= floor(RAM/48)`, and that a feasible row's mean load fits its `C`; the capped
searches with their charged and cap-priced costs, the `F_r` multiply bound, the per-node caps and their bounds (A), (B),
(V), capacity and the two Chernoff exponents (bucket and pass), the storage-node field budget, the success line, memory,
`6d | r-1`, `mu6 subset H`, the H2 density checks); `coverage` (Section 6, plus the distinct-ratio count); `availability`
(the toy analogue of both charged lists by the same segmented sieve and H-test, counts against need, gap condition,
Section 6); `distinct` (the occupied-residue measurement on one toy); `ec` (Algorithm A end to end, with the `d`-th root
selected by matching `T`); `quotes` (every derived quantity the text quotes at the display precision, with the check that
each printed string appears verbatim in this file, in the claim, in `README.md` and in the note, including every layout's
table line and its per-node record counts, so a text edit that drifts from the ledger fails the run; every file is
resolved from the directory `evidence.py` lives in, never from the working directory). Limits: the toy sorts locally
without distributed routing, and the ledger is arithmetic on declared costs.

Provenance: the idea of walking small integers of `H` with one addition per step, the GLV six-fold key and Cheon's
`d | r-1` reduction come from the published note of submission `cbadde77`; the semiprime construction, Lemmas 2-3,
matching, caps, machine model, the two coins, the communication schedule and the overhead listing were derived here.
