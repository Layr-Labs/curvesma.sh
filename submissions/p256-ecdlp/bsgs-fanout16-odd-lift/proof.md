# Odd-lift half-table BSGS with a 16-ary local merge

Status: analytic accounting for exploratory review. No production execution and
no qualification is claimed. Normative inputs are the protected
`target-profile.json` (SHA-256
`802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605`) and
`classical-group-work-v2` (SHA-256
`e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`).
`C = 1000`. Work is scalar AArch64 instruction equivalents. The score is
`log2(T/C)`.

## Claim

`T = 34594432522413595855694169752482540321299522` satisfies both

```
T < 1000 * 2^135
(6/5) * T < 1000 * 2^135
```

so `T/1000 < 2^(1351/10)`. The submitted score is the rational `1351/10`.
Preprocessing `P_pre = 11451624578700418392471984676283525182691180` satisfies
`P_pre < 1000 * 2^134` and is included in `T`. Peak storage
`8920298079412249256614287309059344603463680` bytes is `< 2^143`. Advice is
0 bytes. Success probability is exactly 1. No heuristic premise is used.
Comparisons below are integer comparisons; a floating logarithm is not an
acceptance input.

## Odd lift and coverage

Let `L = 2^128`, `m = L/2 = 2^127`, and `h = ceil(r/L)`. The protected order is

```
r = 115792089210356248762697446949407573529996955224135760342422259061068512044369
h = 340282366841710300967557013911933812736
```

Integer checks: `(h-1)*L < r <= h*L`, `h < L`, `2*(L-1) < r`, and
`m`, `L` and `h` are all divisible by 1024.

The order `r` is odd. For canonical `x` in `[0, r)` the odd lift is
`z = x` when `x` is odd and `z = x+r` when `x` is even (including `x = 0`,
where `z = r`). Then `z` is odd and `1 <= z <= 2r-1`, and `[z]P = [x]P = Q`.

The odd integers partition into blocks indexed by `i >= 0`: block `i` is the
odd integers in the open-closed interval `(2*i*L, 2*(i+1)*L)`. Equivalently,
with center `c_i = (2*i+1)*L`,

```
z = c_i + delta,    delta odd,    -L < delta < L.
```

Because `z <= 2r-1 <= 2*h*L - 1`, the block index satisfies `0 <= i < h`.
Write `delta = s*(2*j+1)` with `s` in `{+1,-1}` and `0 <= j < m`. This
representation is unique. Therefore

```
z = (2*i+1)*L + s*(2*j+1).
```

Store the positive odd babies `B_j = [2*j+1]P` for `0 <= j < m`, and stream
`G_i = Q - [c_i]P` for `0 <= i < h`. The witness gives `G_i = s*B_j`.
An equal finite `x`-coordinate means `G_i = B_j` or `G_i = -B_j`, because a
second solution would force `r` to divide a nonzero even integer strictly
smaller than `2*(L-1)`, contradicting `2*(L-1) < r`. The stored `y` selects
the sign. Both signed candidates are still given a full verification budget.
The raw integer lies in `[1, 2*h*L)` and `2*h*L < 3*r`, so at most two
subtractions of `r` produce a canonical residue in `[0, r)`. Output only after
`[c]P = Q`. The point at infinity is `x = 0`, `z = r`, and is inside the same
partition; it is not a side path.

The baby chain starts at `P` and adds the fixed affine `2P`. The giant chain
starts at `Q - [L]P` and adds the fixed affine `-2L*P`, which is computed
once. Numeric order of the odd positives is exactly the positive half of the
length-`L` odd lattice. No per-record 128-bit scalar multiplication is used
to build the table.

## Local 16-ary merge, not a wide array

A flat array of `2^127` records is not a 64-bit address. This schedule never
forms one. Baby records and giant records are ordered by two separate trees.
Every internal node has at most 16 children. Leaf counts are `m` and `L`,
both `<= 16^32 = 2^128`, so every leaf-to-root path has at most 32 edges.
One central engine edge sits above each root. Every record is therefore
charged for `33` hops, including records whose real path is shorter and
including the `L - h` unused giant leaves, which emit END.

A node owns `2^14` bytes and nothing else:

| Region | Bytes |
| --- | ---: |
| code | 8192 |
| 16 cached heads, 96 bytes each | 1536 |
| 16 one-record queues, 64 bytes each | 1024 |
| 16 heap slots, 8 bytes each | 128 |
| counters, links, tags, stack | 5504 |
| total | 16384 |

Port indices and heap indices are 16-bit local offsets. A load inside a node
is one `LDR` from a pointer the node already holds. No 128-bit pointer is
materialized. Child count is the constant 16, not a function of `m`.

The merge at a node is a 16-way binary heap of child heads. Extract-min is
four sift steps. The node emits the smaller head, stable on ties by child
index, and requests a replacement only from that child. By induction the
stream leaving the root is sorted by `(tag, x)` and preserves multiplicity.
END travels the same edges once the child is exhausted. The engine intersects
the two sorted streams. It may stop at the first equal key; both streams are
charged in full anyway.

Internal-node count of a 16-ary tree on `n` leaves is at most
`(n+14)//15`. Charged nodes are

```
n_b = m + (m+14)//15 = 181483929024500513847133123963609712777
n_g = L + (L+14)//15 = 362967858049001027694266247927219425554
n   = n_b + n_g     = 544451787073501541541399371890829138331
```

## Generation

The fixed addend is converted to Montgomery form once, during setup. The
chain state stays Jacobian. A step copies the state into the batch and then
advances the projective state. The next denominator does not depend on the
inversion of the previous point: every `Z` in the batch already exists before
the product is inverted.

Mixed unequal addition is the `a = -3` Jacobian schedule `7M+4S` plus adds.
The equal-point branch is the `a = -3` double, `8M` plus adds. Four
multiplications are paid before the branch, so either branch is inside
`12M + 32A`. Identity and `Y = 0` return a tagged state inside the same
bound. Every step is charged

```
J = 16*M + 32*A + 512 = 14288
```

with `M = 765` and `A = 48`. The `512` covers copies, full-limb tests and
branches. Squaring uses the multiplication kernel.

`M = 765` is `ceil(6/5 * 637)` for one four-limb CIOS multiply modulo this
`p`. The raw 637 counts 16 schoolbook `MUL`/`UMULH` pairs and their carries
(32 + 144), four REDC passes of four limb-pairs because `p' = 1` (32 + 144),
quotient copies (4), operand traffic (24), scratch zeroing (9), the final
five-word conditional subtraction (16), callee-save allowance (128), round
control (64) and entry/exit (40). `A = 48` is `ceil(6/5 * 40)` for a
four-limb add or subtract with conditional correction. Carries are ordinary
`ADDS`/`ADCS` on 64-bit registers. No sparse-modulus shortcut is taken beyond
`p' = 1`, which removes a multiply by the Montgomery constant and nothing else.

For each exact block of 1024 states: replace a zero `Z` by Montgomery one;
prefix and backward product, `3M` per point; one inversion of the product by
the addition chain for `p-2`, charged

```
I = 512*(M+64) + 256 = 424704
```

and `ceil(I/1024) = 415` per point. Then `2M` to scale `X`, `3M` to scale
`Y` (the second `Z^-2` is recomputed, not aliased), `2M` to leave Montgomery
form on both coordinates, `256` for masks and `512` for record assembly and
local send. Per record,

```
K = J + 10*M + 415 + 256 + 512 = 23121.
```

Both charged chain lengths `m` and `L` are divisible by 1024, as is the real
giant length `h`. There is no partial-block inversion left uncharged. A full
scalar verification is

```
V = 512*J + I + 7*M + 512 = 7746027.
```

## Instruction ledger

Hop price `H = 576`, both endpoints included:

| Row | Instructions |
| --- | ---: |
| four heap sifts, 32 each (tag word, four `x` limbs, select, index swap) | 128 |
| 12-word record transfer, both endpoints | 48 |
| request and response | 52 |
| state and END checks | 96 |
| scheduling and calls | 128 |
| local pointer arithmetic | 64 |
| flags, digit, counters | 60 |
| total | 576 |

Placement on an edge is an additional `128` (transfer plus route-bit and
queue control), the same shape as a one-edge handoff, charged again so the
merge price is not asked to cover movement. Node initialisation is `18432`:
the `16384` local-memory schedule (code distribution at both ends of the
link, zeroing, metadata, control) plus `2048` for sixteen port descriptors.

```
E_b = 33*m
E_g = 33*L
E   = E_b + E_g = 16843977162586453941437043067872526467072

R_init  = 18432 * n
        = 10035335339338780413691073222691762677716992
R_gen   = 23121 * (m+L)
        = 11801502908368527320605026447644869225611264
R_baby  = 64 * m
        = 10889035741470030830827987437816582766592
R_place = 128 * E
        = 2156029076811066104503941512687683387785216
R_sort  = 576 * E + 2 * 576 * n
        = 10329339304358471246123428883512810412390784
R_merge = 512 * (m+L)
        = 261336857795280739939871698507597986398208
R_setup = 4*V + 4*J + 2^20
        = 32089836
R_verify = 2*V + 2^20
        = 16540630
```

`R_baby` is index and odd-label traffic on the baby chain. `R_setup` pays
four full verifications and four Jacobian updates for `P`, `2P`, `3P`,
`[m]P` and the two giant seeds, which dominates `3V` plus those updates, plus
`2^20` words of central scratch. `R_verify` pays two attack-side scalar
multiplications. `R_sort` includes an END request and response at every node.
`R_gen` charges `L` giant records even though only `h < L` are real.

```
T = R_init + R_gen + R_baby + R_place + R_sort + R_merge + R_setup + R_verify
  = 34594432522413595855694169752482540321299522
```

Preprocessing is setup, the baby chain, baby labels, baby nodes, baby
placement and the baby sort, and not the giant chain:

```
P_pre = R_setup + 23121*m + R_baby + 18432*n_b + 128*E_b + 576*E_b + 2*576*n_b
      = 11451624578700418392471984676283525182691180
```

Peak bytes are the charged nodes at `2^14` bytes plus the `2^20` central
reserve:

```
n * 2^14 + 2^20 = 8920298079412249256614287309059344603463680 < 2^143.
```

The inequalities `T < 1000*2^135`, `(6/5)*T < 1000*2^135` and
`P_pre < 1000*2^134` are exact integer statements. The factor `6/5` reprices
every instruction, including those already inside `M` and `A`. It still lands
strictly under the submitted `1351/10`.

## Heuristics

None. Coverage is an integer partition of the odd residues. The trees have a
fixed shape. Batch inversions run on Jacobian `Z` coordinates that the chain
already produced. Success is 1 on every valid input, which is at least 0.39.
The toy experiment checks the odd-lift identity on a small additive group
only. It is not a measurement of `T` and it is not used as a premise for the
production bound.

## Limitations

The ledger is an analytical instruction schedule. It is not a Callgrind count
and not a claim that the trees fit in any real machine. Memory on this scale
is mandatory for a generic half-table of `2^127` points. Beating the integer
bound is not a practical break and does not change the hardness assumption
used by P-256. The prior odd-lift half-table with a binary merge is earlier
work on the same coverage identity; the routing degree, hop price and totals
in this file are self-contained and do not depend on that package.
