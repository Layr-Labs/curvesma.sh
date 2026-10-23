# Cheon p-1 rectangle with mu_6 orbit keys on BLS12-381 G1

Status: ready for independent review. This is a concrete classical attack with a
complete cost ledger. It is not a best-known-attack claim and it is not
qualified by any local run. Everything below is stated at the scoring length
`n = 2^20`; no other `n` is given a score.

## 1. Target, inputs and claim

The target group is G1 = <P0> of prime order

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

on `y^2 = x^3 + 4` over the protected Fp. The input is `T = A_1 = [s]P0` and the
powers `A_i = [s^i mod r]P0` for `i = 1..n`, with `n = 2^20 = 1048576` and `s`
uniform in `[1, r-1)`. The attack reads `A_1` and `A_n` only; it never scans the
rest of the array. It outputs a scalar `a` and accepts only when `[a]P0 = T`.

The claim is

```text
For uniform s in [1, r-1) and n = 2^20, the attack in Sections 3-4 recovers s
with probability at least 0.395 and with total work at most
W = 3.001960433063365e35 = 2^117.853388 target-group-operation equivalents,
including table construction, every walk, sorting, lookup, lifting and
verification. Peak storage is at most 2^123.466026 bytes.
The submitted score is log2(W) <= 117.86.
```

The attack is deterministic. It uses no randomness, so the success probability
is a property of the challenge alone and no algorithm coins need to be averaged.

## 2. Parameters and their verification

`d = 2^20` and `d | r-1`; this is the scoring length, and `A_d = A_n` is
supplied. The factorization used for the primitive-root check is

```text
r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2 * 2508409
      * 2529403 * 52437899 * 254760293^2
```

`7` is a primitive root of `Z_r^*`: `7^((r-1)/q) != 1 (mod r)` for each of the
twelve distinct prime factors `q` above. Set

```text
N = (r-1)/d = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff000
```

so `N` has 235 bits, `6 | N` (`N = 2^12 * 3 * ...`), and `d*N = r-1`. Define

```text
gamma = 7,  eta = gamma^d mod r,  zeta = gamma^N mod r,  lambda = gamma^((r-1)/3) mod r.
eta  = 0x628fc184374a1219e85d9e943d57d0c20bad4d8f813859f6121b51c927eaa4dc
zeta = 0x3e1c54bcb947035a57a6e07cb98de4a2f69e02d265e09d9fece7e0e39898d4b
lambda = 0xac45a4010001a40200000000ffffffff
```

`eta` has order exactly `N`, `zeta` has order exactly `d`, and `lambda^2 +
lambda + 1 = 0 (mod r)`; all three are checked by direct exponentiation.

The curve has `j = 0`, so the order-six automorphism group acts on G1:
`phi(x, y) = (lambda*x, y)` is the order-three endomorphism (its eigenvalue on
G1 is `lambda`), and negation is `[k]P0 -> [-k]P0`. The six scalars

```text
mu_6 = {1, -1, lambda, -lambda, lambda^2, -lambda^2}  (mod r)
```

are distinct. Since each has order dividing 6 and `gamma` generates `Z_r^*`,
each equals `gamma^((r-1)t/6)` for some `t`, hence

```text
log_eta(u) = t*N/6  (mod N),  u in mu_6,  t in {0..5},
```

which is integral because `6 | N`. Direct evaluation gives the exact pairing

```text
t=0 -> 1,  t=1 -> -lambda^2,  t=2 -> lambda,
t=3 -> -1, t=4 -> lambda^2,   t=5 -> -lambda.
```

The six values `t*N/6` are distinct modulo `N`. This identity is the only place
the endomorphism enters, and it is checked numerically rather than assumed; the
toy experiment in `experiment.py` re-derives it in a small group with the same
structure.

## 3. Construction

Write `a = log_gamma(s)`, which is uniform in `[0, r-1)` because `s` is uniform
in `[1, r-1)` and `gamma` is primitive. Set `k = a mod N`.

1. **Giant base.** `A_d = [s^d]P0 = [gamma^(a*d)]P0 = [eta^a]P0 = [eta^k]P0`,
   because `eta` has order `N`. The problem instance therefore *is* a discrete
   logarithm of `k` to the known base `eta` in the cyclic group `G1`, with
   `0 <= k < N`.
2. **Digit tables.** Split a scalar `e < r` into three 85-bit digits
   `e = v_0 + v_1*2^85 + v_2*2^170` with `0 <= v_i < 2^85`; this is exact
   because `3*85 = 255` and `r < 2^255`. Build six affine tables
   `T_i^P[v] = [v*2^(85i)]P0` and `T_i^A[v] = [v*2^(85i)]A_d`, for `i = 0,1,2`
   and `v = 0..2^85-1`. Every entry of `T^P` is built by one repeated addition
   from `[2^(85i)]P0`, and every entry of `T^A` the same way from
   `[2^(85i)]A_d`. Then for any scalar `e`,
   `[e]P0 = T_0^P[v_0] + T_1^P[v_1] + T_2^P[v_2]` (two additions) and
   `[e]A_d = T_0^A[v_0] + T_1^A[v_1] + T_2^A[v_2]` (two additions).
   The giant tables are built from `A_d`, not from `P0`: this is what makes
   `[eta^(-G i)]A_d` available without knowing `s^d` as a scalar.
3. **Baby side.** For `j = 0..G-1`: maintain `e_j = eta^j mod r` by one
   scalar-field product per step, split it into digits, and materialise
   `B_j = [e_j]P0` with two additions. Store the record
   `(min(y(B_j), p-y(B_j)), j)`.
4. **Giant side.** For `i = 0..G-1`: maintain `e_i = eta^(-G i) mod r` by one
   scalar-field product per step, split it into digits, and materialise
   `C_i = [e_i]A_d = [eta^(k-G i)]P0` with two additions. Store the record
   `(min(y(C_i), p-y(C_i)), i)`.
5. **Dictionary.** Sort each of the two record arrays by its 381-bit key with a
   deterministic six-pass LSD radix sort on 64-bit digits and 128-bit counters,
   then merge-join the two sorted arrays on key equality. Every equality found
   is one candidate match `(i, j)`.
6. **Lift.** For each match `(i, j)`, for each `t` in `{0..5}`, form
   `k_t = (G*i + j + t*N/6) mod N`; for each `t'` in `[0, d)`, form
   `a' = k_t + t'*N` and test `[a']P0 = A_1`. Return the first `a'` that
   passes. If no match exists, or no candidate passes, report failure.

Both sides are streamed in blocks: a block of `2^20` points is materialised,
its two addition layers are batched, records are emitted, and the block is
discarded. The scalar updates are sequential (they are cheap integer products);
the *point* additions are independent inside a layer, which is what allows one
batch inversion per layer.

## 4. Correctness and success probability

**Key lemma.** For a point `Q` in G1, `min(y(Q), p-y(Q))` determines the
mu_6 orbit of `Q`. If `min(y(Q1), p-y(Q1)) = min(y(Q2), p-y(Q2))` then
`y(Q2) = +-y(Q1)`, hence `x(Q2)^3 = y(Q1)^2 - 4 = x(Q1)^3`, hence
`x(Q2) = lambda^c x(Q1)` for some `c`, hence `Q2` is one of the six points
`(lambda^c x(Q1), +-y(Q1))`. Conversely all six share the same `y` up to sign.
All points here have `y != 0` because G1 has odd prime order, so no orbit
degenerates.

**Match characterisation.** `C_i` has key equal to `B_j`'s key exactly when
`eta^(k-G i) = u*eta^j` for some `u` in `mu_6`, i.e. exactly when

```text
k = G*i + j + t*N/6  (mod N)  for some t in {0..5}.
```

No probabilistic hashing is involved: the radix sort and the merge-join compare
the full keys, so this is an exact equivalence, not a heuristic.

**Covering.** Let `S = { (G*i + j + t*N/6) mod N : 0 <= i,j < G, 0 <= t < 6 }`.
The map `(i,j) -> G*i + j` is injective into `[0, G^2)`, and
`G^2 < N/6` for the chosen `G` (checked: `G^2 = 3.292e69` and
`N/6 = 8.334e69`). The six translates of `[0, G^2)` by multiples of `N/6` are
therefore pairwise disjoint, so `|S| = 6*G^2` exactly.

**Uniformity and probability.** `s` uniform on `[1, r-1)` makes `a` uniform on
`[0, dN)` and `k = a mod N` exactly uniform on `[0, N)`, since `N | r-1`. Hence

```text
Pr[success] = |S| / N = 6*G^2/N = 0.395  >= 0.39
```

with `G = 57376919326004192012165322720755522` (116 bits). This is a worst-case
figure: the algorithm either covers `k` or does not, and it covers with
probability exactly `6G^2/N`.

**Lift.** If `k` is covered then some match `(i,j)` satisfies
`k = G*i + j + t*N/6` for the true `t`. Write `a = k + t_a*N` with
`t_a = floor(a/N) in [0, d)`. The enumeration in step 6 contains
`t = t_true` and `t' = t_a`, hence contains `a`, and the verification
`[a]P0 = A_1` accepts it. At most six matches exist (for each `t`, the value
`(k - t*N/6) mod N` lies in `[0, G^2)` for at most one pair `(i,j)`), so the
enumeration is bounded by `36*d` tests. No false accept is possible: the final
test is an exact group equality and `s -> [s]P0` is injective on `[0, r)`.

Edge cases: `s = 1` gives `a = 0`, `k = 0`, which lies in `S` through
`(i,j,t) = (0,0,0)`; the baby record `j = 0` is `[eta^0]P0 = P0`, so the
previous off-by-one failure mode does not arise. Digits `v_i = 0` produce the
identity in a table slot; the block layer drops identity operands from the
batch and the screening step below charges that branch. Equal-x and
inverse-point pairs inside a layer are handled by the same screening.

## 5. Cost ledger

The reference unit is one batched incomplete affine BLS12-381 G1 addition
(6M with `S=M`, `A=0`) at `C = 2000`. Non-reference work is charged as its
scalar AArch64 instruction count divided by 2000. `H` counts reference
additions, `R` counts instructions.

### 5.1 Per point (each side)

| Item | Bound |
| --- | ---: |
| Two batched incomplete affine additions | 2 reference units |
| Scalar-field product `e <- e*rho mod r`, 4x4 CIOS + Montgomery reduction + final subtraction | 512 instructions |
| Digit split (three 85-bit shifts/masks, three stores) | 48 |
| Exceptional screening, 2 additions x 128 (x-equality test, identity test, select) | 256 |
| Record build (`min(y, p-y)` select and index store) | 64 |
| Radix-sort share, 6 passes x 48 instructions | 288 |
| Merge-join share | 64 |

`1232/2000 = 0.616` reference units of non-reference work plus 2 reference
additions, so `2.616` reference units per point. The scalar product is charged
at 512 instructions although a four-limb CIOS product with a four-limb
Montgomery fold is about 208; the excess is deliberate headroom for loads,
stores, addressing and loop control. Likewise the screening charge is a full
128 instructions per addition for a check that is a handful of comparisons.

### 5.2 Main loop

`2G = 1.1475e35` points at `2.616` reference units each:
`2G*2.616 = 2^117.853388` (this term is the whole bound to within `2^-27`).

### 5.3 Preprocessing

Six tables of `2^85` affine entries. Each entry is one sequential addition in
Jacobian coordinates (no per-step inversion) plus its share of a batch
normalisation to affine; charged at 6 reference units per entry:

```text
6 * 2^85 * 6 = 1.392e27 = 2^90.170
```

Setup (primitive-root check, `eta`, `zeta`, `lambda`, the six shift constants,
block scratch) is charged `2^30`. Preprocessing total `2^90.170`, all inside
`W`.

### 5.4 Dictionary

Records are 64 bytes: a 48-byte key (381 bits, zero-extended) and a 16-byte
index. Two sort buffers hold `4G*64 = 2^123.466` bytes. Radix counters are
128-bit values, so a bucket count of at most `G < 2^116` cannot overflow; the
previous four-byte-counter failure mode is excluded. Counter clearing and
prefix sums over `2^64` buckets for six passes cost `6*2^64*8 = 2^69.6`
instructions, i.e. `2^66.6` reference units, charged explicitly. Record moves
are already inside the per-point sort share; the merge-join is a single linear
pass over `2G` records and is charged 64 instructions per record.

### 5.5 Recovery

At most 6 matches, each expanded to 6 shifts and `d` lifts, so at most `36*d`
candidates. Each is tested by a full binary scalar multiplication, charged 512
reference units:

```text
36 * 2^20 * 512 = 2^34.170
```

Candidate scalar construction (`36*d` scalar-field products) is far below this
and is covered by the same term.

### 5.6 Total

```text
W = 2G*2.616 + 36*2^85 + 36*2^20*512 + 2^66.6 + 2^30
  = 3.001960433063365e35
log2(W) = 117.85338828289322
```

The submitted score is `117860000/1000000 = 117.86`, a ceiling of `log2(W)`
with a `2^0.0066 = 1.0046` allowance above the ledger, in addition to the
per-step headroom described in 5.1.

## 6. Memory and advice

```text
peak bytes = 4G*64 (two sort buffers) + 6*2^85*96 (six affine tables)
             + 2*2^64*16 (radix counters) + block scratch
           = 2^123.466026 + 2^94.170 + 2^69 + 2^26
peak     <= 2^123.466026 bytes
```

Advice embedded in code is `0` bytes. The tables are built at run time from the
public inputs `P0` and `A_d`; their construction is the `2^90.170` term of `W`
and their bytes are counted in peak storage, so nothing is precomputed for free
and nothing is amortised over future targets.

## 7. Heuristics, limitations and uncertainty

- No heuristic is used for the matching step. The dictionary is an exact sort
  and merge; there is no occupancy model, no overflow probability and no
  randomised hash family.
- The covering argument is exact, not a Poisson or birthday approximation:
  `|S| = 6G^2` by disjointness, and `k` is exactly uniform because `N | r-1`.
- The `0.395` figure is a worst-case success probability, not an average over
  algorithm coins; the algorithm is deterministic.
- The only analytical quantities are instruction ceilings for the scalar-field
  product, the exceptional screening, the sort moves and the Jacobian table
  build. Each is padded well above its expected value and the total claim
  carries a further multiplicative allowance.
- No production-size execution is claimed. The toy experiment is evidence only
  at its own scope: it checks the orbit-shift identity, the exact covered-set
  size `6G^2`, and that every covered scalar lifts to the correct one, in a
  small multiplicative group with `d = 8`. It does not establish the production
  instruction counts.
- The attack is a cold attack on one target; it does not amortise anything over
  other instances and uses no quantum or uncharged oracle step.
