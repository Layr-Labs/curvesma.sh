# Deterministic truncated baby-step/giant-step with a hashed baby table on BLS12-381 G1

Status: complete claim for independent review under `exploratory-v1` and
`classical-group-work-v2`. The algorithm below is fully specified: every table,
record field, hash step, probe, sign branch and failure path is stated, and the
success fraction is exact rather than heuristic. The submitted score is
`time_log2 = 127.875678` in BLS12-381 G1 reference operations.

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
```

`m = h` because the charged cost `m + h` (plus a small per-element table term)
is minimised at `m = h` for fixed `m*h`. The coverage is fixed first; the sizes
follow from it. Both fit in 127 bits, so every index and every product `i*m + j`
fits one 256-bit register pair, and every slot index fits 128 bits.

## 3. Algorithm

### 3.1 Record and slot layout

- Record, 64 bytes = 8 little-endian 64-bit words: words 0..5 are the 384-bit
  compressed `x`-coordinate of the baby point (only 381 bits are used, the top
  three are zero); word 6 is the tag, `0` for a finite point and `1` for the
  identity; word 7 holds `i` in bits 0..31 and the sign bit of `y(i*P)` in bit
  32.
- Slot array: `2^128` entries of 8 bytes each, holding `i+1` (so `0` means
  "empty"). `i < 2^127`, so the encoding is injective.
- The record array holds exactly `m` records, `i = 0 .. m-1`, and includes the
  identity at `i = 0`.

### 3.2 Setup (charged once, no reuse credit)

1. `T_0 = O` (identity). For `i = 0 .. m-1`, append the record for `T_i` and
   set `T_{i+1} = T_i + P0` by one mixed affine addition, retaining
   `T_m = [m]P0` after the final update.
2. Hash each record into the slot array: `s = hash(record) mod 2^128`; while
   slot `s` is occupied, `s = (s + 1) mod 2^128`; store `i+1` in slot `s`.
   Keys are distinct (Lemma 1), so no record is ever displaced and the load
   factor is `m / 2^128 < 0.5`.
3. Set `D = -[m]P0` (one point negation, coordinate-wise).

The hash is a fixed 64-bit multiplicative-xor mixing of the six key words and
the tag word: `acc = 0xcbf29ce484222325`; for each of the seven words
`acc = (acc xor w) * 0x100000001b3 mod 2^64`; then
`acc = acc xor (acc >> 29)`, `acc = acc * 0xbf58476d1ce4e5b9 mod 2^64`,
`acc = acc xor (acc >> 32)`. This is a public, key-independent function; it is
not a random oracle, and nothing below assumes collision resistance. A collision
only causes an extra probe.

### 3.3 Per challenge

```
input T; if T is the identity, return 0
gamma = T
for j = 0 .. h-1:
    if gamma is the identity: (tag, key) = (1, 0,0,0,0,0,0)
    else:                     (tag, key) = (0, compressed_x(gamma))   # 6 words
    s = hash(key, tag) mod 2^128
    while slot s is occupied:
        i = slot[s] - 1
        if record[i].key == key and record[i].tag == tag:
            # gamma = [i]P0 or gamma = -[i]P0; for i = 0 both branches give j*m
            if record[i].sign == sign_of_y(gamma):
                a = (j*m + i) mod r
            else:
                a = (j*m - i) mod r
            verify [a]P0 == T; if it holds, return a   # it always holds
            break                                   # unreachable, defensive
        s = (s + 1) mod 2^128
    gamma = gamma + D
return FAIL
```

The candidate check `[a]P0 == T` is a double-and-add over 255 bits, charged in
the ledger. It is not needed for correctness (Theorem 1 shows every hit is a
true hit) but it is executed, and the cost is charged.

## 4. Correctness

**Lemma 1 (keys are distinct).** For `0 <= i < i' < m`, `[i]P0 != +/-[i']P0`.

*Proof.* `[i]P0 = [i']P0` implies `r | (i'-i)`, impossible since
`0 < i'-i < m < r`. `[i]P0 = -[i']P0` implies `r | (i+i')`, impossible since
`0 < i+i' < 2m < r` (`m < r/2`). Distinct points therefore have distinct
compressed keys, so no two records share a slot key and the probe sequence
always terminates on the correct record or on an empty slot.

**Lemma 2 (interval coverage).** Write `x = j*m + i` with `0 <= i < m`. If
`x < m*h` then `0 <= j < h`.

*Proof.* `j = floor(x/m) <= (m*h-1)/m < h`. QED

**Theorem 1 (exact success).** For uniform `x` in `[0, r)` the algorithm returns
`x` with probability exactly `(m*h + m - 1)/r`, and it never returns a wrong scalar.

*Proof.* If `x < m*h`, Lemma 2 gives `j < h` and the giant step
`gamma = T - [j*m]P0 = [x - j*m]P0 = [i]P0`, which is record `i`. The lookup
finds it, and the sign branch returns `j*m + i = x` (the opposite sign branch
would return `j*m - i`, which is `x` only if `i = 0`; in that case both
branches agree). The verification accepts.
If `m*h <= x <= r - m`, no `j < h` has `x - j*m` in `[0, m)`, and no `j` has
`j*m - x` in `[0, m)` either, so no giant step hits and the answer is `FAIL`.
If `r - m < x < r`, write `x = r - i` with `0 < i < m`; at `j = 0`,
`gamma = T = -[i]P0`, whose key is record `i`, and the opposite-sign branch
returns `-i mod r = x`. The three sets `[0, m*h)`, `[m*h, r-m]` and `(r-m, r)`
partition `[0, r)`, and `x = r-m` is the single uncovered residue of the top
block (`(0 - (r-m)) mod r = m`, which is not `< m`), so the success count is
`m*h + (m-1)` out of `r` residues, exactly.
A hit always means `gamma = +/-[i]P0`, hence `T = [j*m +/- i]P0`, so every
returned `a` satisfies `[a]P0 = T`; with `r` prime and `P0` of order `r` the
scalar in `[0, r)` is unique, so no wrong answer is possible. QED

**Corollary.** `(m*h + m - 1)/r > 0.41`, above the required `0.39`. With
`m*h = 21498708821801738096573573608356245993507346963331534251623150184000105601600`
and `r` as above, `m*h/r = 0.4100000000000000000000000000000000000010...`, so
the success bound has margin `0.02` and does not rely on any heuristic.

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
| setup: decoding, constants, table metadata | `2^24` |
| slot array initialisation, one 64-bit store per slot | `2^128` |
| table build: `m` probes, 128 per insert (hash 7 words, probe, 8-word record store, slot store) | `128*m` |
| giant search: `h` probes, 128 per lookup (hash, probe, 8-word record load, 6-word key compare, branch) | `128*h` |
| verification/output: `[a]P0` check, encoding, serialisation | `2^24` |
| **`R`** | **37876123656888034495974394388096459696128** |

`R = 2^25 + 2^128 + 128*(m + h)`. The 128-instruction probe bound is roughly
twice the 60–70 instructions a straightforward 64-bit implementation of the
hash and the open-addressing probe needs; it also covers loop control and
index arithmetic. No constant-time table oracle, cache model or free hashing is
assumed.

### 5.3 Total

```
W = H + R/C = 293248760077867937753992342036442641100 + 37876123656888034495974394388096459696128/2000
  = 624373643812623910003959078460981741896128 / 2000
  = 3.1218682190631197e38
time_log2 = log2(W) = 127.875677244984...  -> submitted 127.875678
```

The submitted rational is `127875678/1000000`, rounded upward from
`127.8756772449841` at six decimals. The exact integer `W` numerator above is
the audit authority.

### 5.4 Preprocessing, memory, advice

- Preprocessing (the whole baby table, its slot array and setup; no reuse
  credit, fully included in `W`):
  `P = m + (2^24 + 2^128 + 128*m)/2000 = 1.5617848154488621e38`, so
  `preprocessing_log2 = 126.876464`.
- Peak memory: slot array `2^128 * 8 = 2^131` bytes, records `64*m` bytes,
  `2^20` bytes for code, constants, scratch and stacks:
  `12106219257859281715834751804620311222784` bytes, `memory_log2_bytes =
  133.152873`. This is astronomical and mandatory for the construction; it is
  not a practical attack, and memory carries no scalar contribution.
- Nonuniform advice: `0` bytes; the algorithm embeds no searched results.
  Every table entry is built and charged inside `W`.

## 6. Scope, evidence and limitations

The correctness argument and the cost calculation are analytic; no
production-size execution is claimed and none is needed. The only supporting
experiment is `reference.py`, a scope-limited analogue over a generic cyclic
group with the same quotient-by-negation key structure; it exercises the exact
table construction, open-addressing probe, sign branch, interval coverage and
failure path on small orders, and it is evidence only at that scope.

There are no heuristic premises: the success fraction is an exact count over
the uniform challenge distribution, the algorithm uses no randomness, and the
failure set is stated explicitly (`[m*h, r-m]`). The construction is a
deterministic, truncated, memory-heavy baby-step/giant-step, and the claim is
only that this is a valid accounted attack at the stated bound, not that it is
the best possible method.
