# Truncated baby-step/giant-step with a seeded bounded-probe baby table on BLS12-381 G1

Status: complete claim for independent review under `exploratory-v1` and
`classical-group-work-v2`. Every table, record field, probe, sign branch and
failure path is stated; every probe window is bounded by construction, so the
cost bound is worst case and does not depend on hash behaviour. The submitted
score is `time_log2 = 127.929648` in BLS12-381 G1 reference operations.

## 1. Target

- `E: y^2 = x^3 + 4` over `F_p`,
  `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`.
- `G1 = E(F_p)[r]`, cyclic of prime order
  `r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`
  (`log2 r = 254.857089`).
- `P0` is the profile generator; the challenge is `T = [x]P0` with `x` uniform
  in `[0, r)`, including `x = 0` (`T` the identity). There are no auxiliary
  powers and no `n` parameter. Output is the scalar `x`.

## 2. Parameters

```
m = h = 0x6e4ed412a3c552fd2d010d8328cef768 = 146624380038933968876996171018221320040
P = 8                      # probe cap, identical for insertion and lookup
alpha = m / 2^128 = 0.43089032607146766
alpha^P = 0.0011883214889591687
```

`m = h` because the charged cost `m + h` (plus the per-probe table term) is
minimised at `m = h` for fixed `m*h`; the coverage is fixed first and the sizes
follow from it. Both fit in 127 bits, so every index and every product `j*m + i`
fits one 256-bit register pair, and every slot index fits 128 bits.

## 3. Algorithm

### 3.1 Record and slot layout

- Record, 64 bytes = 8 little-endian 64-bit words: words 0..5 are the 384-bit
  compressed `x`-coordinate of the baby point (only 381 bits are used, the top
  three are zero); words 6..7 are one 128-bit integer holding `i` in bits 0..126
  and the sign bit of `y([i]P0)` in bit 127. The sign bit is stored because the
  compressed form omits `y`.
- Slot array: `2^128` entries of 16 bytes each, holding `i+1` as a 128-bit
  little-endian integer, so `0` means "empty". `i + 1 <= m < 2^127`, so the
  encoding is injective and the `0` sentinel is never a valid record index.
- The record array holds exactly `m` records, `i = 0 .. m-1`, including the
  identity at `i = 0`. A record is either *indexed* (its `i+1` sits in one slot)
  or *unindexed*; an unindexed record stays in the record array but is never
  found by a lookup.

### 3.2 Key

For a point `Q` in `G1` the key is a 6-word (384-bit) vector: the six words of
the compressed `x`-coordinate, with `key = 0` for the identity `O`.

The key is unambiguous. `x(Q) = x(-Q)`, so two points share a key iff they are
equal or opposite, and the identity is the only self-opposite point. No finite
point of `G1` has `x = 0`: `(0, ±2)` lies on the curve with `2*(0,2) = (0,-2)`,
so it has order 3, while `G1 = E[r]` with `r` prime and `r != 3`. The all-zero
key therefore identifies the identity exactly.

### 3.3 Setup (charged once, no reuse credit)

1. Draw the probe family from fresh coins: `P = 8` independent pairs
   `(a_t, b_t)` with `a_t` uniform odd and `b_t` uniform in `Z/2^128`, i.e. 2048
   bits of seed, charged in the `2^24` setup term. The family is the seeded
   multiply-shift probe `s_t(k) = (a_t * (k mod 2^128) + b_t) mod 2^128`: for a
   fixed key the eight positions are independent uniform, and two records share
   a position only when their 384-bit keys agree modulo `2^128`. It is public
   and key-independent once the seed is fixed.
2. `T_0 = O` (identity). For `i = 0 .. m-1`, append the record for `T_i` and set
   `T_{i+1} = T_i + P0` by one mixed affine addition, retaining `T_m = [m]P0`.
3. Index each record: for `t = 0 .. P-1` compute `s_t = (a_t * key + b_t) mod
   2^128`; if slot `s_t` is empty, store `i+1` there and stop. If all `P` slots
   are occupied, record `i` is left unindexed. Keys are distinct (Lemma 1), so
   insertion never needs to compare keys; a collision only costs a probe.
4. Set `D = -[m]P0` (one point negation, coordinate-wise).

### 3.4 Per challenge

```
input T; if T is the identity, return 0
gamma = T
for j = 0 .. h-1:
    key = key(gamma)
    for t = 0 .. P-1:
        s = (a_t * key + b_t) mod 2^128
        e = slot[s]
        if e == 0: break                      # this j misses
        i = e - 1
        if record[i].key == key:
            # gamma = [i]P0 or gamma = -[i]P0; for i = 0 both give j*m
            if record[i].sign == sign_of_y(gamma):
                a = (j*m + i) mod r
            else:
                a = (j*m - i) mod r
            verify [a]P0 == T; if it holds, return a
            break                             # unreachable, defensive
    gamma = gamma + D
return FAIL
```

The candidate check `[a]P0 == T` is a double-and-add over 255 bits, charged in
the ledger. It is not needed for correctness (Theorem 1 shows every hit is a
true hit) but it is executed, and the cost is charged.

### 3.5 Probe-window cost

One probe is at most 23 scalar instructions: 3 for the multiply-shift
(`multiply`, `add`, `mask`), 2 for the 16-byte slot load, 2 for the empty test
on the 128-bit slot, 8 record-word loads, 6 key-word comparisons, 2 for branch
and index arithmetic. Eight probes plus loop control are at most 184 + 8 = 192;
an insertion additionally stores the 8-word record and the 2-word slot, so the
window charge is 208 instructions. Every insertion and every lookup therefore
costs at most 208 instructions regardless of the table state: the window is
capped at `P = 8` probes by construction, and no probe sequence can run longer.

## 4. Correctness

**Lemma 1 (keys are distinct).** For `0 <= i < i' < m`, `[i]P0 != +/-[i']P0`.

*Proof.* `[i]P0 = [i']P0` implies `r | (i'-i)`, impossible since
`0 < i'-i < m < r`. `[i]P0 = -[i']P0` implies `r | (i+i')`, impossible since
`0 < i+i' < 2m < r` (`m < r/2`). Distinct points therefore have distinct
compressed keys, so no two records share a key and no slot ever holds two
records. QED

**Lemma 2 (interval coverage).** Write `x = j*m + i` with `0 <= i < m`. If
`x < m*h` then `0 <= j < h`.

*Proof.* `j = floor(x/m) <= (m*h-1)/m < h`. QED

**Lemma 3 (a placed record is always found).** If record `i` was stored at probe
position `t < P` of its key, every lookup of that key probes the same positions
in the same order and finds `i+1` at position `t`, with a full key comparison
that matches.

*Proof.* The probe family and the table are fixed after setup; the lookup
recomputes the identical sequence `s_0, s_1, ...` from the same key. No slot is
ever overwritten (keys distinct, one store per record), so the word `i+1` is
still there at position `t`. QED

**Theorem 1 (exact success conditional on the indexed set).** Let `A` be the
number of indexed records. For uniform `x` in `[0, r)` the algorithm returns `x`
on at least `A*h + A - 1` residues and never returns a wrong scalar.

*Proof.* If `x < m*h`, Lemma 2 gives `j < h` and the giant step
`gamma = T - [j*m]P0 = [x - j*m]P0 = [i]P0` with `i = x - j*m`. If `i` is
indexed, Lemma 3 finds it, and the sign branch returns `j*m + i = x` (the
opposite sign branch returns `j*m - i`, which equals `x` only if `i = 0`, where
both branches agree). If `m*h <= x <= r - m`, no `j < h` has `x - j*m` in
`[0, m)` and no `j < h` has `j*m - x` in `(0, m)`: for `j = 0` that would need
`x > r - m`, and for `j >= 1`, `j*m - x < j*m <= (h-1)m < m*h`. If
`r - m < x < r`, write `x = r - i` with `0 < i < m`; at `j = 0`,
`gamma = T = -[i]P0`, whose key is that of record `i`, and if `i` is indexed the
opposite-sign branch returns `-i mod r = x`. The identity input is answered by
the pre-check. Counting: the indexed `i` give `A*h` residues below `m*h`, and
the indexed `i > 0` give `A - 1` residues above `r - m`; the residue `x = r - m`
is the single uncovered top residue (`(0 - (r-m)) mod r = m`, not `< m`), so the
success count is at least `A*h + A - 1`. A hit always means
`gamma = +/-[i]P0`, hence `T = [j*m +/- i]P0`; with `r` prime and `P0` of order
`r` the scalar in `[0, r)` is unique, so no wrong answer is possible. QED

**Lemma 4 (indexing loss).** Averaged over the seed, the expected number of
unindexed records is at most `m * alpha^P = 1.7424e35`, i.e. `0.001188` of the
table.

*Proof.* For a fixed record the `P` probe positions are independent uniform in
`[0, 2^128)` (each `b_t` is uniform and independent, and `a_t` odd makes
`s_t` a bijection of the reduced key). At the moment record `i` is inserted, at
most `m` slots are occupied. Under the random-hash model for the seeded family,
the probability that all `P` positions are occupied is at most
`(m/2^128)^P = alpha^P`. Linearity of expectation over the `m` records gives
the bound. QED

**Corollary (success).** Take the tolerance `D <= 0.005*m`. Then
`A >= 0.995*m` and Theorem 1 gives success at least
`((0.995*m)*h + 0.995*m - 1)/r = 0.4079500000000000000000000000000000000000...`,
above the required `0.39` with margin `0.018`. The model of Lemma 4 predicts a
loss of `0.001188*m`, so the tolerance is `4.21` times the expected loss: even
a fourfold failure of the random-hash model leaves the success bound intact.
Under the model, `Var[D] <= m*alpha^P` and Chebyshev gives
`Pr[D > 0.005*m] <= 5.6e-37`, below `2^-120`; the bound is therefore not
sensitive to the tail.

The single modelling premise in this claim is the random-hash model of Lemma 4,
used only to count unindexed records. It is disclosed, it is not used in the
cost bound, and the `4.21x` margin absorbs a large violation of it.

## 5. Ledger

Reference operation: one batched incomplete affine BLS12-381 G1 addition
(6M convention, `S = M`, `A = 0`), price `C = 2000` scalar AArch64 instructions.

### 5.1 Reference operations `H`

| Component | Count |
| --- | ---: |
| baby table: `m` mixed additions (`T_{i+1} = T_i + P0`) | `m` |
| giant steps: `h` mixed additions (`gamma += D`) | `h` |
| verification: two 255-bit double-and-add scalar multiplications, `2*(255+255)` | `1020` |
| **`H`** | **293248760077867937753992342036442641100** |

`H = m + h + 1020`. The doubling and addition cases of the affine formulas are
priced as the same reference operation, which is conservative for doublings.

### 5.2 Other instructions `R`

| Component | Instructions |
| --- | ---: |
| setup: seed draw, decoding, constants, table metadata | `2^24` |
| slot array initialisation, two 64-bit stores per 16-byte slot | `2^129` |
| table build: `m` bounded probe windows, 208 per window | `208*m` |
| giant search: `h` bounded probe windows, 208 per window | `208*h` |
| verification/output: `[a]P0` check, encoding, serialisation | `2^24` |
| **`R`** | **61676306830038407979757156358443639113984** |

`R = 2^25 + 2^129 + 208*(m + h)`. The 208-instruction window is the worst case
of section 3.5 and is a hard bound, not an average: each window runs at most
eight probes because the loop is capped at `P`. No constant-time table oracle,
cache model or free hashing is assumed.

### 5.3 Total

```
W = H + R/C = 293248760077867937753992342036442641100 + 61676306830038407979757156358443639113984/2000
  = 648173826985774283487741840431328921313984 / 2000
  = 3.2408691349288714e38
time_log2 = log2(W) = 127.929648372071...  -> submitted 127.929648
```

The submitted rational is `127929648/1000000`, rounded upward from
`127.92964837207101` at six decimals. The exact integer numerator above is the
audit authority.

### 5.4 Preprocessing, memory, advice

- Preprocessing (the whole baby table, its slot array, the seed and setup; no
  reuse credit, fully included in `W`):
  `P_pre = m + (2^24 + 2^129 + 208*m)/2000 = 1.6221359792990404e38`, so
  `preprocessing_log2 = 126.931162`.
- Peak memory: slot array `2^128 * 16 = 2^132` bytes, records `64*m` bytes,
  `2^20` bytes for code, constants, scratch and stacks:
  `14828478193226789423541748664074456914432` bytes, `memory_log2_bytes =
  133.445494`. This is astronomical and mandatory for the construction; it is
  not a practical attack, and memory carries no scalar contribution.
- Nonuniform advice: `0` bytes; the algorithm embeds no searched results. Every
  table entry and the seed are produced and charged inside `W`.

## 6. Scope, evidence and limitations

The correctness argument and the cost calculation are analytic; no
production-size execution is claimed and none is needed. The only supporting
experiment is `reference.py`, a scope-limited analogue over a generic cyclic
group with the same quotient-by-negation key, the same bounded probe window and
the same drop rule; it exercises table construction, the capped probe, the sign
branch, the coverage lower bound and the failure path on small orders, and it
is evidence only at that scope.

The cost bound is worst case: the probe cap makes every insertion and lookup a
fixed 208-instruction window, so no hash behaviour can enlarge the charge. The
success bound is exact conditional on the indexed set (Theorem 1) and averaged
over the seed with a `4.21x` tolerance margin (Lemma 4 and its corollary). The
construction is a truncated, memory-heavy baby-step/giant-step; the claim is
only that this is a valid accounted attack at the stated bound, not that it is
the best possible method.
