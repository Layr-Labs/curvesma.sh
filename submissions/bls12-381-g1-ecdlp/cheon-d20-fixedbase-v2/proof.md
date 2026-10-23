# bls12-381-g1-ecdlp: Cheon d=2^20 with fixed-base windowed BSGS

Status: ready for exploratory review. This is an analytic attack and accounting
claim, not a full-size execution, a novelty claim or a practical break.

## Target and notation

The protected profile fixes the BLS12-381 G1 group of prime order

```text
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
```

and a generator `P0`. At the scoring length `n = 2^20 = 1048576` the challenge
supplies `A_i = [s^i mod r] P0` for `i = 1..n` and `T = A_1 = [s]P0`, where `s` is
uniform in `[1, r-1]`. Write `G = <P0>`, so every point is `[x]P0` for a unique
`x mod r`, and write `F_r^*` for the multiplicative scalar field.

The attack uses only the single supplied power `A_d` with `d = 2^20`. Since
`2^20 | r-1`, this is the largest admissible `d <= n`; no larger divisor of `r-1`
is at most `n`. The exponent `2^20` is fixed, so the same work function applies at
`n = 2^20` and at every larger `n`; the bound is not claimed for smaller `n`.

Set

```text
M1 = (r-1)/d = 50006747412801924209068050869165387952509453297164571592906626415194112
m  = ceil(sqrt(M1)) = 223621884914696853578958571250384896
g1 = ceil(M1/m)     = 223621884914696853578958571250384896
N  = m + g1         = 447243769829393707157917142500769792
```

One checks `m^2 - M1 = m*g1 - M1 = 218380746987008646073201729736704 > 0` and
`(m-1)^2 < M1`, so `m` and `g1` are the exact BSGS parameters. Also `N =
2*sqrt(M1)` up to the same small ceiling, and `log2 N = 118.42854470652...`.

## The two Cheon phases

Because `d | r-1`, let `zeta` be an element of order `r-1` in `F_r^*` (a primitive
root modulo `r`), and write `s = zeta^a` with `0 <= a < r-1`. Define

```text
zeta_d = zeta^d        (order M1 = (r-1)/d)
omega  = zeta^M1       (order d)
```

**Phase 1.** Since `zeta_d` has order `M1`, the point `A_d = [s^d]P0 =
[zeta_d^a]P0` determines `a mod M1`. Write `k0 = a mod M1`. The standard
representation is `k0 = u + m*v` with `u = k0 mod m` and `v = floor(k0/m)`, so
`0 <= u < m` and `0 <= v < g1`; the pair exists. Because `m*g1 > M1`, the
scanned rectangle can contain one more pair with `u + m*v = k0 + M1` when
`k0 < m*g1 - M1 = 218380746987008646073201729736704`; reducing modulo `M1`
returns the same `k0`. The identity used is

```text
[zeta_d^{-u}] A_d = [zeta_d^{m*v}] P0.
```

The left side is computed for all `u` in `[0,m)` (baby steps, base `A_d`); the
right side for all `v` in `[0,g1)` (giant steps, base `P0`). A match returns
`k0 = (u + m*v) mod M1`. This is equation (1) of Cheon's Lemma 1 adapted to the
supplied `A_d`.

**Phase 2.** The value `a` is `k0 + M1*k1` for a unique `k1` in `[0,d)` because
`0 <= a < M1*d = r-1`. Define

```text
T' = [zeta^{-k0}] T = [omega^{k1}] P0.
```

Write `k1 = u' + m2*v'` with `m2 = ceil(sqrt(d)) = 1024`, `0 <= u' < m2`, and
`0 <= v' < g2 = ceil(d/m2) = 1024`. The tested identity is

```text
[omega^{-u'}] T' = [omega^{m2*v'}] P0.
```

Baby steps over `u'` use the point `T'` as a variable base (one variable-base
scalar multiplication each, charged at 384 additions); giant steps over `v'`
reuse the `P0` tables. A match returns `k1 = u' + m2*v'`. The output is

```text
s = zeta^{k0 + M1*k1} mod r,
```

followed by the mandatory attack-side check `[s]P0 = T` (one more fixed-base
multiplication with the `P0` tables). Returning a point, a relation or a next
power is excluded by this check.

The only primitive-root machinery is: factor `r-1`, test candidate `zeta` in
`F_r^*` until `zeta^{(r-1)/q} != 1` for every prime `q | r-1`. The factorization
is found by trial division up to `2^28` (the largest prime factor of `r-1` is
`254760293 < 2^28`), at a charged cost far below `2^40` instructions; the test is
charged in the same budget. The smallest primitive root modulo `r` is `7`
(candidates `2,3,4,5,6,7`), so the search is only a few tens of thousands of
255-bit modular multiplications; trial division to `2^28` and sequentially
reading the input sequence up to `A_n` also cost far below the `2^40` budget.
The `2^40` setup charge therefore covers the input access with a wide margin. No
factorization or primitive root is embedded as uncharged advice. For reference,
the factorization checked by the experiment is

```text
r-1 = 2^32 * 3 * 11 * 19 * 10177 * 125527 * 859267 * 906349^2
      * 2508409 * 2529403 * 52437899 * 254760293^2.
```

## Fixed-base windowed multiplication

All scalars in the two phases are known field elements: `zeta_d^{-u}` for the
baby steps, `zeta_d^{m*v}` for the phase-1 giants, and `omega^{...}` for phase 2. Only
two point bases occur in phase 1: `A_d` and `P0`. Because `r < 2^255`, three
85-bit windows cover every scalar:

```text
x = x0 + 2^85 x1 + 2^170 x2,   0 <= xj < 2^85.
```

For each base `B` in `{A_d, P0}` and each `j` in `{0,1,2}` build the table

```text
T[B,j][c] = [c * 2^(85*j)] B,   c = 0..2^85-1.
```

Then `[x]B = T[B,2][x2] + T[B,1][x1] + T[B,0][x0]`: three lookups and exactly two
point additions (the first table entry initializes the accumulator). The tables
are built by repeated addition from `B`, `[2^85]B`, and `[2^170]B`; the latter two
cost 85 and 170 doublings. Table construction costs

```text
6*(2^85 - 1) additions + 2*(85+170) doublings
```

for the two bases, and is charged conservatively at three reference operations
per entry (a projective addition plus batched affine conversion), giving
`H_pre = 3*6*(2^85-1) + 510 = 696341272098026404630757868`.

The main loop therefore performs

```text
H_main = 2*(m + g1) = 2*N = 894487539658787414315834285001539584
```

reference additions: two per scalar multiplication, for `N = 447243769829393707157917142500769792`
scalar multiplications. This is the dominant term, `log2 H_main = 119.42854470652...`.

The `6M` convention is the track's batched incomplete affine addition. In the
main loop the two additions per scalar multiplication are processed in blocks
(block size `2^18`): the first addition for every element of a block is batched,
then the second, so the `3M` batch-amortisation part of the `6M` convention is
realized with a batched inversion, and a block's inverse and output buffers stay
below `2^26` bytes in total. The table-construction additions are charged
at three reference operations each, which covers the less convenient batched
conversion there.

Identity table entries (`c=0`) are skipped by a branch whose instruction cost is
included in the lookup and digit charges. For every nonzero scalar the two
additions need no exceptional correction, and this is proved directly rather
than bounded statistically. The base `B` is `A_d` or `P0`, both of order `r`,
so `[u]B` and `[v]B` have equal x-coordinates exactly when `u = v` or
`u + v = r`. For the first addition, `T[B,2][x2] = [2^170 x2]B` and
`T[B,1][x1] = [2^85 x1]B` can have equal x-coordinates only if
`2^85 x2 = +/- x1`; the two sides differ by less than `r` because
`2^170 + 2^85 < r`, so the congruence is an integer equality, and `0 <= x1`
then forces `x2 = x1 = 0` (a skipped addition). For the second addition let
`S = [x - x0]B` be the first sum and `T[B,0][x0] = [x0]B` the third entry.
Equal x-coordinates give `x = 2x0` (equivalently `x - x0 = x0`), `x = 0` or
`x = r` in the integers; for `0 < x < r` only `x = 2x0` survives, and then
`x0 = x - x0 = 2^85*(2^85*x2 + x1)`, impossible since `0 <= x0 < 2^85`.
Hence no incomplete-addition exception occurs in either Cheon phase or in the
final verification, and the digit-zero skip branches are already charged in the
lookup and digit terms. `H_exc = 1024` remains as rounding and branch slack; the
`2*(85+170)` base doublings and the entry `c = 2` of each table are already
charged in `H_pre`.

## Matching

The baby points `B_u = [zeta_d^{-u}]A_d` are distinct: their exponents
`zeta_d^{-u}` are distinct for `u` in `[0,m)`. Build an open-addressing hash
table with `2^118` slots. Each slot stores the full affine point `(x,y)` (96
bytes) and the index `u` (15 bytes; `m < 2^118`), padded to 112 bytes. The
initial slot is selected by a 118-bit index derived from a 128-bit hash of the
full point; the hash is computed per operation and charged in `R_hash`, and the
128-bit output covers all `2^118` slots uniformly. A probe compares the full
96-byte point, so a hit is an exact point equality, not an x-only or
sign-ambiguous match. The table is zero-initialized; the all-zero encoding is
not an affine point, so it serves as the empty sentinel. The load factor is
`m/2^118 = 2^-0.571... = 0.673... < 0.71`. The probe count is the only
average-case quantity: at this load factor the expected probes are about 2 for
a successful lookup and about 5 for an unsuccessful one. It is not charged at
its expectation: the itemised insert path (104 instructions) plus 36 per extra
probe is charged at 500 per operation, which is exactly the 12-probe budget of
P2 below, and the algorithm aborts (returns FAILURE) if its global probe budget
is exhausted, so the probe work of every attempt is bounded by construction.
Correctness never depends on the hash; any salt is charged.

Each giant step `G_v = [zeta_d^{m*v}]P0` looks up its full affine point. On a
hit, the attack computes `k0 = (u + m*v) mod M1`, runs phase 2, and accepts only
if `[s]P0 = T`. Every full-point match satisfies `u + m*v ≡ k0 (mod M1)` and
therefore returns the correct `k0`; the possible second rectangle pair (only for
`k0 < 218380746987008646073201729736704`) returns the same `k0`. Hence every
accepted hit is correct, and hash collisions only cost probes and cannot create
a false full-point match. A hash miss or a failed check continues the fixed
scan. The true match is present, so the last loop always returns the correct
`s`.

## Itemised instruction accounting for the non-addition terms

`R` is charged in scalar AArch64 retirement counts; loads, stores, addressing,
control, hashing and memory transfers are included. The model used here: 64-bit
little-endian limbs, an `Fp` element = 6 limbs = 48 bytes, an affine G1 point
`(x, y)` = 96 bytes = 12 limbs, one open-addressing slot = 112 bytes = 14 limbs
(96-byte point + 16-byte index/pad), one comb entry = the same 96-byte affine
point. Each table below prices **an explicit instruction sequence**, and the
charged figure is the itemised sum times the stated factor. Nothing here is a
timing measurement, and nothing is asserted without its itemisation.

### P1 — one comb lookup (3 per fixed-base scalar multiplication), charged 200

The comb tables are batch-normalised to affine at construction, so a lookup is a
read, a conditional negation and an install; the group addition itself is priced
separately in `H` (two reference additions per scalar multiplication).

| step | instructions |
| --- | ---: |
| digit extract (`LSR` + `AND` for the 85-bit digit) | 2 |
| sign test and conditional-negate mask (`TBZ`, mask build) | 3 |
| byte offset `digit*96` (`LSL`, `ADD`, `ADD` base) | 3 |
| 12 aligned 8-byte loads (`LDR`) | 12 |
| conditional negation, 12 limbs (`EOR`+`SUB`, 1.5 avg) | 18 |
| install into accumulator registers (12 x `MOV`) | 12 |
| **digit-zero exceptional case** (identity entry) test and branch | 3 |
| window loop control, amortised over 3 digits | 2 |
| **itemised sum** | **55** |
| **charged** | **200** (3.6x) |

The `c = 0` entry is the identity and is not an affine point; the branch is not
optional and is priced above, and the no-incomplete-addition proof of the
previous section shows no other exceptional path occurs.

### P2 — one hash-table operation (1 per fixed-base scalar multiplication), charged 500

Open addressing over `2^118` slots of 112 bytes, load factor `m/2^118 = 0.673`,
128-bit hash of the full 96-byte affine point, full 96-byte point compare on
every probe (never x-only), insert on the first empty slot.

| step | instructions |
| --- | ---: |
| 128-bit hash of the 12-word point (12 x 3-step mix + 2 final mixes) | 40 |
| slot index (`AND` 118-bit mask, fold) | 3 |
| slot address `base + idx*112` (`LSL`, `MADD`) | 3 |
| load 14-word slot (14 x `LDR`) | 14 |
| empty-sentinel test and branch | 3 |
| 96-byte full-point compare (12 `CMP`, 11 `ORR`, `CBZ`) | 15 |
| store on insert (14 x `STR`) | 14 |
| next-probe index (odd-constant multiply, mask) | 4 |
| probe-loop control, amortised | 8 |
| **insert path** | **104** |
| each additional probe (load, sentinel, compare, step, control) | 36 |
| **12-probe bound: 104 + 11 x 36** | **500** |

The probe count is the **only average-case quantity in this candidate**, and it
is handled deterministically rather than at its expectation. The algorithm
carries a global probe budget of `12*N` probe steps over its `N` operations and
returns FAILURE if the budget is exhausted, so every attempt costs at most
`N*(104 + 11*36) = 500*N` probe instructions by construction. At load factor
0.673 the expected probe count is `~2` for a successful and `~5.2` for an
unsuccessful lookup, and unsuccessful lookups dominate (`N` giants, at most a
few hits), so the expected total is `<= 6*N`, comfortably inside the `12*N`
budget. The probe tail is geometric with ratio `m/2^118` (the standard
linear-probing law, disclosed here as the only probabilistic input); a Chernoff
bound over the `N ~ 2^58.6` near-independent lookups gives

```text
Pr[total probes > 12*N] <= exp(-N * D(12/6 || 1)) <= exp(-0.35*N) <= 2^-1.6e17,
```

so the budget is a bound and not an event, in the same sense as the hard
evaluation cap of the two-phase loop. The four Cheon rejections on this board
named exactly this: an allowance asserted at its expectation. It is now derived
line by line and capped.

### P3 — one scalar chain (1 per `zeta` power), charged 300

| step | instructions |
| --- | ---: |
| load 8 operand limbs | 8 |
| 4x4 schoolbook (16 x `UMULH`+`MADD`+carry) | 40 |
| Montgomery reduction, 4 rounds x 6 | 24 |
| conditional subtraction of `r` (4 x `CSEL`+`SUB`, borrow) | 8 |
| store 4 limbs | 4 |
| loop control, amortised | 4 |
| **itemised sum** | **88** |
| **charged** | **300** (3.4x) |

This is also below the package's earlier comparison point, the measured
383-instruction 6-limb base-field multiply with full reduction: a 4-limb
Montgomery multiply is a strictly smaller operation, and 300 exceeds it by 3.4x.

### P4 — digit extraction and loop control (1 per scalar multiplication), charged 100

| step | instructions |
| --- | ---: |
| 3 x (digit extract, offset, control) | 21 |
| accumulator loop control and address updates | 10 |
| **itemised sum** | **31** |
| **charged** | **100** (3.2x) |

### The charged numbers do not carry the claim

Doubling all three contested allowances moves the claim by +0.33 bits and 8x
moves it by +1.47 bits (the ledger arithmetic is in `experiment.py` and in the
recomputation below), so this derivation buys *support*, not score:

| comb | hash | chain | claim log2 |
| ---: | ---: | ---: | ---: |
| 200 | 500 | 300 | 119.8911 (submitted) |
| 400 | 1000 | 600 | 120.2176 |
| 1600 | 4000 | 2400 | 121.3651 |

### Layout consistency

The comb entries are 96-byte affine points, produced by batch normalisation from
a 144-byte Jacobian staging buffer over one block of `2^18` entries = 2^25.17
bytes, with a prefix/inverse array of `2^18` field elements = 2^25.17 bytes and
an affine output buffer of `2^18` points = 2^25.58 bytes. The three together are
2^26.2 bytes, i.e. below the `2^27`-byte block/code/scratch term of the memory
ledger, and they are transient inside one block: only the normalised 96-byte
entries persist in the `6*(2^85-1)*96`-byte tables. So the 96-byte comb stride
and the 144-byte Jacobian stride do not contradict each other, and the block size
appears in no counted quantity except the (negligible) per-block inversion.

## Complete cost ledger

The track charges `W = H + R/C` with `C = 2000` for BLS12-381 G1 and
`H` counting reference additions. `R` is in scalar AArch64 instruction
equivalents; loads, stores, addressing, control, hashing and memory transfers are
included. No reference operation and its implementation are double counted.

| Component | Count | Reason |
|---|---:|---|
| `H_pre` | 696341272098026404630757868 | 6*(2^85-1) table additions + 510 doublings, at 3 ref ops per table entry |
| `H_main` | 894487539658787414315834285001539584 | 2 additions x N fixed-base scalar multiplications |
| `H_p2` | 395650 | 1025 variable-base phase-2 scalar multiplications (`T'` shift plus 1024 baby steps) at 384 additions, plus 1025 fixed-base giant/verification multiplications at 2 additions |
| `H_exc` | 1024 | rounding and branch slack; table doublings are in `H_pre` |
| **H** | **894487540355128686413860689632694126** | |

| Instruction component | Count `R` | Reason |
|---|---:|---|
| comb lookups | `N*3*200 = 268346261897636224294750285500461875200` | P1: itemised 55, charged 200 (3.6x); 3 lookups per scalar multiplication |
| hash operations | `N*500 = 223621884914696853578958571250384896000` | P2: itemised 104 insert + 36 per extra probe; 500 covers the 12-probe global budget (104+11*36), abort on exhaustion |
| scalar chains | `N*300 = 134173130948818112147375142750230937600` | P3: itemised 88, charged 300 (3.4x), for one 4-limb Montgomery multiplication plus canonical reduction per `zeta` power |
| digit extraction/loop | `N*100 = 44724376982939370715791714250076979200` | P4: itemised 31, charged 100 (3.2x); three 85-bit digits, accumulator addressing |
| table stores | `6*(2^85-1)*100 = 23211375736600880154358578600` | writing 96-byte affine entries during precomputation |
| hash-table zeroing | `2^122 = 5316911983139663491615228241121378304` | paired 16-byte zero stores over `2^124.807` bytes plus loop control, rounded up |
| factorization/primitive root/setup | `2^40 = 1099511627776` | trial division to `2^28`, primitive-root tests (smallest root `7`), parameter arithmetic, sequential input access |
| phase-2 support | `2^22 = 4194304` | modular exponentiation of the `zeta^{...}`/`omega^{...}` scalars, dict/hash and loop overhead beyond counted additions |
| output/verification support | `2^12 = 4096` | encoding `s`, final comparisons |
| **R** | **676182566750441599965091823246150471080** | |

Therefore

```text
R/C = 338091283375220799982545911623075235.54
W   = H + R/C = 1232578823730349486396406601255769361.54
log2 W = 119.8910913262144667779170456173938538953...
```

The submitted bound is `2^119.891092`, which is at least `2^(log2 W)`; the
rounding margin is under `10^-6` bits. All preprocessing, table construction,
the 2^122-instruction zeroing, the setup and the final verification are included
in `W`; none is amortized over another target.

## Memory, preprocessing and advice

Peak memory is dominated by the open-addressing table: `2^118` slots of 112
bytes = `2^124.8073549...` bytes. The six comb tables add `6*(2^85-1)*96` bytes
(`2^94.17`), and the block buffers, the two retained input points (`A_d`, `P0`),
code and scratch add less than `2^27` bytes. The exact declared accounting value
is

```text
(2^118)*112 + 6*(2^85-1)*96 + 2^27 = 37218383904260565148443442636168101312 bytes
log2 = 124.807354922921356...  <= 124.807355.
```

This memory is astronomical and is required by the BSGS table; it is reported as
a bound, not a practical resource budget. Memory has no score contribution.

Preprocessing is declared as the fixed-base comb construction, the factorization
and primitive-root search, and the hash-table initialization:

```text
H_pre + (R_stores + R_setup + R_init)/2000
= 2658456699516791712134459378253660.34
log2 = 111.03421609952...  <= 111.034217.
```

It is included in `W`. No table or search result is embedded in the candidate
files or supplied as nonuniform advice: `nonuniform_advice_log2_bytes = 0`. The
input powers are public challenge inputs; their generation is organizer work and
is not charged, but every attack-side access to `A_d` and every runtime table is
charged and included in peak memory. All tables depend on the target (`A_d`) or
are rebuilt from public constants; nothing is reused across targets.

## Success probability and applicability

For every valid challenge the representation `a = k0 + M1*k1` is unique, the
phase-1 pair `(u,v)` exists, and the phase-2 pair `(u',v')` exists. The scans are
fixed and complete, the hash table stores every baby point, and the final check
rejects any false candidate. The algorithm is deterministic apart from the
deterministic small-integer search for `zeta`, which succeeds because `F_r^*` is
cyclic. Hence the success probability is 1, above the required 0.39, over the
uniform choice of `s` and any fresh coins. No heuristic assumption is needed for
correctness; the open-addressing probe count is the only average-case quantity
and is charged well above its expectation. The bound holds
for `n = 2^20` and any `n > 2^20`; the attack reads only `A_d` with `d = 2^20`.

## Why this bound, and rejected alternatives

The phase-1 search is the BSGS of Cheon's Lemma 1 with the supplied `d = 2^20`.
Cheon's `O(sqrt((r-1)/d) + sqrt(d))` counts group exponentiations; replacing each
known-scalar exponentiation by a three-window fixed-base lookup+addition chain
reduces the per-step cost from a full variable-base scalar multiplication to two
reference additions. The `sqrt(d)` phase is `2^10` and negligible. The choice
`d = 2^20` minimizes the phase-1 modulus even without the divisor restriction:
for any exponent `d` the power `A_d` determines `a` modulo
`(r-1)/gcd(d, r-1)`, and `gcd(d, r-1) <= d <= n = 2^20`, so the phase-1 modulus
is smallest at `d = 2^20`; a smaller `d` increases `M1` and the work, while phase
2 stays at `sqrt(gcd(d, r-1)) <= 2^10`. The `p+1` variant applies to divisors of
`r+1` instead; any divisor `d <= n = 2^20` of `r+1` gives
`(r+1)/d >= (r+1)/2^20 > M1`, so its phase-1 search is larger and it is not used.
A second supplied power `A_e` with `e <= n` yields only a second residue modulo
`(r-1)/gcd(e, r-1) >= M1`; recovering it independently costs another BSGS of at
least the same size (another `N` scalar multiplications), and combining the two
residues leaves an ambiguity of at most `2^20`, no better than the `2^10` phase 2.
So no multi-power combination is used, and the other supplied powers are not
needed. The map `a -> [zeta_d^a]P0` used in phase 1 is not a group
homomorphism on the index, so Pollard-rho-style collision search does not replace
the BSGS table here; no memoryless shortcut is claimed or used.

This is not claimed to be the best possible attack, a security lower bound, or a
practical break. The bound is analytic and could be tightened by a more detailed
instruction model; the declared charges are deliberately conservative in the
lookup, hash, scalar-chain and initialization terms, which together contribute
about 0.46 bits to `log2 W`.

## Evidence and reproducibility

`experiment.py` (stdlib only) is declared in `submission.json`. On a toy
prime-order group of about 31 bits it runs the same two-phase Cheon algebra with
the same fixed-window lookup/addition structure, full-point matching, both
phases, and the final check, for seeded random and boundary scalars; it also
checks the production 3x85-bit digit structure on boundary and seeded sample
scalars (the no-exception proof itself is analytic and covers every
`0 < x < r`), verifies the factorization of `r-1` and that `7` is
the smallest primitive root modulo `r`, recomputes the production `M1`, `m`,
`g1`, exact integer ledger, score bound, preprocessing and memory bound, and
asserts the claimed bound is not exceeded. The
experiment is evidence only at toy scope; it does not execute a production-size
attack or validate production memory. The repository's unit tests check the
protected BLS12-381 profile constants, the subgroup-order factor `2^20`, the
supplied auxiliary powers and the pinned pairing conventions.

No production-size run is claimed. The ledger is analytic; all integer
components are stated above so a reviewer can recompute the bound. The main
uncertainty is the instruction accounting for the large-table lookups and the
hash operations, charged far above their minimal instruction counts; the dominant
`2N` additions are independent of those charges.

## References

- J. H. Cheon, "Discrete Logarithm Problems with Auxiliary Inputs", Journal of
  Cryptology 23 (2010), Lemma 1, Theorem 1 and the KKM precomputed-table
  speed-up cited there. The identities used above are equations (1)-(3) of the
  p-1 case, adapted to the supplied `A_d`.
- Track rules, target profile and `classical-group-work-v2` cost model in this
  repository; the BLS12-381 G1 `6M` addition convention with `C = 2000`.
