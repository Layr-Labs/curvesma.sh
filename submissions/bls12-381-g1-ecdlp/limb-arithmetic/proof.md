# bls12-381-g1-ecdlp: deterministic generic accounting control

Status: ready for independent exploratory review; no qualification or emitted score
is claimed beyond the submitted rational bound. This is a conservative accounting
control, not a best-known-attack claim and not a practical execution.

The normative problem, encoding and acceptance rule are those of the protected
`target-profile.json`; the accounting convention is `classical-group-work-v2` and
`tracks/docs/COST_MODEL.md`. No restricted challenge distribution, uncharged
advice, random-oracle premise, quantum operation or amortisation over future
targets is used.

## Target and parameters

- `p` = `1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab` (381 bits).
- `r` = `73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001` (255 bits).
- `G1 = <P0>` has prime order `r`; the target is `T = [s]P0` and the powers are
  `A_i = [s^i]P0` for `i = 1..n`.
- Scoring value `n = 2^20 = 1048576`.

Set `b = bit_length(r) = 255`, `L = floor(b/2) = 127`, `m = 2^L = 2^127`,
`h = ceil(r/m) = 308190375243917312857644309053803048971`. Then
`h*m >= r`, so every `x` in `[0,r)` has a unique representation
`x = i*m + j` with `0 <= j < m` and `0 <= i < h`.

## Construction and success

The construction reads only `A_1 = T` (and, when `n >= 2`, `A_2` for one
independent post-recovery check). This is a permitted subset of the supplied
powers: `target-profile.json` states that the attack may use a subset and need
not read or validate the whole sequence. The auxiliary sequence is therefore not
required by this construction, and `W(n)` is independent of `n`.

```text
Input: P0 (public generator), T = A_1, and (optionally) A_2.
Use P = P0; input Q = T.

Baby table. Start at infinity and, for j = 0..m-1, store the record
(x(P_j), y(P_j), j) where P_j = [j]P0, producing P_j from P_{j-1} by one
affine addition; retain P_m = [m]P0 after the final update.
Set D = -P_m and gamma = Q.
Sort the m records by the x-coordinate with a deterministic radix sort
(24 passes of 16-bit digits over the 381-bit x-coordinate; x-collisions are
resolved by comparing y, and an exact (x,y) match yields j).

Giant search. For i = 0..h-1, binary-search gamma in the sorted table. On an
exact match with index j, set a = (i*m + j) mod r, verify [a]P0 = Q (and, when
n >= 2, [a^2]P0 = A_2), and return a. Otherwise set gamma = gamma + D.
Reject malformed inputs; exhaust the fixed loops with failure.
```

For every valid input there is a unique `s` in `[0,r)` with `Q = [s]P0`. Write
`s = i*m + j`. At giant step `i`, `gamma = Q - [i*m]P0 = [s - i*m]P0 = [j]P0`,
so the table contains an exact match and the search returns `a = s` after the
reduction modulo `r`. The independent equality checks accept it. Infinity gives
`s = 0` and is represented explicitly. The algorithm succeeds for every valid
input, so its success probability is 1, which is at least the required 0.39.
There is no randomness, retry or amplification term.

## Reference operations and instruction-equivalent ledger

One reference operation is one batched incomplete affine BLS12-381 G1 addition
under the 6M convention (`S = M`, `A = 0`), fixed at `C = 2000` scalar AArch64
instruction equivalents. Reference additions are charged as `H`; all other work
is charged as instructions `R` and divided by `C`. An addition's internal
instructions are not charged twice.

All non-reference work is bounded on a scalar 64-bit basis. A stored record is
`14` words (six words for `x`, six for `y`, two for the index). A word
load/store pass is bounded by 16 instructions, covering addressing, alignment,
bounds and control, which is the same conservative rate used by the repository's
secp256k1 accounting control.

| Kernel | Upper bound in scalar instructions |
| --- | ---: |
| Record write (14 words) | `16*14 + 64 = 288` |
| Record read (14 words) | `16*14 + 64 = 288` |
| Key comparison (x, 6 words) | `16*6 + 64 = 160` |
| Radix pass per record (read + write + digit/bucket/index) | `288 + 288 + 64 = 640` |
| Radix sort per record (24 passes) | `24*640 = 15360` |
| Binary search per query (128 comparisons + record access) | `128*160 + 288 + 256 = 21024` |
| Setup and output scratch | `2^17` |

Setup computes `[m]P0` by 127 doublings and decodes `P0` and `T`. Verification
performs `[a]P0` and `[a^2]P0`, each bounded by `2*b = 510` affine additions, and
serialises the 32-byte scalar. Output construction, the attack-side equality
checks, the malformed-input rejection and the loop exhaustion are covered by the
scratch reserve.

Charged reference operations:

```text
H = m (baby additions) + h (giant additions) + 127 (setup) + 510 + 510 (verification)
  = m + h + 1147
  = 478331558704386544589331612769687155846.
```

Charged other work:

```text
R = m*(288 + 15360) + h*21024 + 2^17
  = 15648*m + 21024*h + 131072
  = 9141763687917540123656556882093309788129120 instructions.
```

Thus the total is

```text
W = H + R/C
  = 478331558704386544589331612769687155846
    + 9141763687917540123656556882093309788129120/2000
  = 126230335066578915160440251345408551247764/25
  = 5.0492134026631566... * 10^39
  target-group-operation equivalents.
```

`log2(W) = 131.89125435...`, rounded upward to six decimals: **131.891255**.
Because the construction ignores `n`, `W(n) = W` for every `n >= 1`; the bound
applies in particular at the scoring value `n = 2^20`.

Preprocessing (the entire baby table plus its sort, with no reuse credit) is

```text
m + m*(288 + 15360)/2000 = 187665725356897562600051095998620168617984/125
                         = 1.5013... * 10^39,
log2 = 130.141433 (rounded upward).
```

It is included in `W`; no preprocessing is excluded or amortised.

## Memory

The construction stores two record arrays of `m` records for the radix sort,
`2*m*112` bytes, plus `24*65536*4` bytes of counting-sort buckets and a `2^24`
byte reserve for code, public constants, field/point scratch, loop state and
stacks:

```text
peak bytes = 38111625095145107907897956032358062751744,
log2(peak bytes) = 134.807355 (rounded upward).
```

This memory is astronomical and mandatory for this construction; it is not a
practical attack. Memory has no scalar contribution or tie-break, and is reported
for completeness. Nonuniform advice is `0` bytes: no table is precomputed from
the secret or embedded in code, and the baby table is built at attack time and
fully charged.

## Heuristics, uncertainty and evidence

Heuristic premises: none. The construction is deterministic, exhaustive and
success probability 1; the radix sort and binary search are exact, and the
instruction bounds are analytic. The only extrapolation is the stated
instruction-per-word and per-comparison allowance, which is deliberately
conservative relative to the pinned measurement basis in `COST_MODEL.md`.

The correctness argument and cost calculation are analytic; no full-size
execution is claimed. Small experiments would be evidence only at their actual
scope, never a production complexity proof. The exponential memory is a property
of the accounting control, not a claim of feasibility.

## Submitted computation bound

- Total: `W = 126230335066578915160440251345408551247764/25` target-group-operation
  equivalents; `log2(W(2^20)) = 131.891255` (rounded upward).
- Preprocessing: included, `log2 = 130.141433`.
- Peak memory: `log2 = 134.807355` bytes.
- Nonuniform advice: `log2 = 0` bytes.
- Success probability: `1 >= 0.39`.
