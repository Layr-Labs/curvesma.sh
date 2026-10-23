# Odd-lift half-table BSGS with byte-radix local pages

Status: analytic accounting for exploratory review. No production execution and
no qualification is claimed. Normative inputs are the protected
`target-profile.json` (SHA-256
`802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605`) and
`classical-group-work-v2` (SHA-256
`e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294`).
`C = 1000`. Work is scalar AArch64 instruction equivalents. The score is
`log2(T/C)`.

## Claim

`T = 28068769515497374754164342230842034976561602` satisfies both

```
T < 1000 * 2^(1348/10)
(6/5) * T < 1000 * 2^(1348/10)
```

The submitted score is the rational `1348/10`. Preprocessing
`P_pre = 9276403576395011358628708835736690067785068` satisfies
`P_pre < 1000 * 2^134` and `(6/5)*P_pre < 1000 * 2^134`, and is included in
`T`. Peak storage `278759314981632789269196478408104519892992` bytes is
`< 2^143`. Advice is 0 bytes. Success probability is exactly 1. No heuristic
premise is used. Comparisons below are integer comparisons; a floating
logarithm is not an acceptance input.

## Odd lift and coverage

Let `L = 2^128`, `m = L/2 = 2^127`, and `h = ceil(r/L)`. The protected order is

```
r = 115792089210356248762697446949407573529996955224135760342422259061068512044369
h = 340282366841710300967557013911933812736
```

Integer checks: `(h-1)*L < r <= h*L`, `h < L`, `2*(L-1) < r`, and
`m`, `L` and `h` are all divisible by 1024 and by 32.

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

## Pages, byte radix, bounded ports

A flat array of `2^127` records is not addressed by a scaled 128-bit index.
This schedule never forms one. Baby records and giant records live in two
separate trees of 16384-byte pages. A page holds `32` records and nothing
addressed from outside except through a child slot whose base pointer was
planted at initialisation. The planted pointer is the register base the node
already holds. Port indices and record slots are 16-bit offsets from that
base. No pointer is materialised as `global_base + ordinal * stride`.

Each record is 80 bytes: the affine `x` (32), the affine `y` (32) and the
chain index (16). Two banks of 32 records occupy `2 * 32 * 80 = 5120` bytes.
The rest of the page is the radix state:

| Region | Bytes |
| --- | ---: |
| code | 4096 |
| 256 histogram counters, 2 bytes each | 512 |
| 256 cursors, 2 bytes each | 512 |
| 256 ports, capacity 4, 16-bit offsets | 2048 |
| links, pass index, path digits, spare | 4096 |
| two record banks | 5120 |
| total | 16384 |

The in-page order is a stable LSD byte radix, 32 passes, one pass per byte of
`x`, low byte first. A pass is a counting partition, not a heap:

1. Clear the 256 counters. Every counter is cleared on every pass, including
   digits that do not occur.
2. Histogram the current key byte of each of the 32 records.
3. Turn the histogram into 256 exclusive prefix cursors. Cursor `b` is the
   16-bit start offset of a destination run whose length is `hist[b]` and is
   therefore at most 32.
4. Walk the source bank in order. Push the source slot onto port `b` until
   the port holds 4 offsets, then spill those four records into the
   destination bank at the cursor and advance the cursor. A port that still
   holds a short tail at the end of the pass is spilled once.

Stability is the arrival order inside each port. The destination bank becomes
the source bank on the next pass. After 32 passes the page is fully ordered
by `x`, because every byte has been a radix digit and the partition is
stable. Ties on `x` keep chain order. Inside one tree an equal finite `x`
does not occur: the same `2*(L-1) < r` bound used above makes the baby
`x`-coordinates distinct and the giant `x`-coordinates distinct. The radix
does not rely on that distinctness; it still terminates on duplicates.

The pass price is `128` instructions per record, both banks included:

| Row | Instructions |
| --- | ---: |
| load the local slot | 8 |
| extract the current key byte | 8 |
| update the 16-bit histogram or cursor | 12 |
| push a 16-bit offset into a port of capacity 4 | 16 |
| spill predicate | 8 |
| amortized spill of an 80-byte record, both endpoints | 48 |
| bank base maintenance | 16 |
| pass index, flags, counters | 12 |
| total | 128 |

The amortized spill is `floor(32/4) = 8` full spills per page plus at most
one short spill per port. The short spills are not left inside the `128`.
They are charged again, for every port, on every pass, whether or not the
port is occupied:

```
R_tail = 16 * 256 * 32 * P
```

with `16` instructions of residual bookkeeping per port.

Page counts, exact because `m` and `L` are divisible by 32:

```
P_b = m/32 = 5316911983139663491615228241121378304
P_g = L/32 = 10633823966279326983230456482242756608
P   = P_b + P_g = 15950735949418990474845684723364134912
```

`R_radix = 128 * 32 * (m+L)` charges every baby record and every one of the
`L` giant slots, including the `L - h` unused giant slots.

## Local 16-ary merge of sorted pages

Sorted pages are the leaves of two complete 16-ary merge trees, one per
chain. Every internal node has at most 16 children and the same local-memory
shape already used for a 16-way head merge: 8192 bytes of code, 16 cached
heads of 96 bytes, 16 one-record queues of 64 bytes, 16 heap slots, and the
remaining bytes for counters and links, totalling 16384. A head is the page's
current `(x, y, index)` triple plus the local child slot. Child count is the
constant 16. Port and heap indices are 16-bit local offsets. A load inside a
node is one `LDR` from a pointer the node already holds.

Internal-node count of a 16-ary tree on `n` leaves is at most `(n+14)//15`:

```
I_b = (P_b+14)//15 = 354460798875977566107681882741425221
I_g = (P_g+14)//15 = 708921597751955132215363765482850441
I   = I_b + I_g   = 1063382396627932698323045648224275662
```

Both `P_b = 2^122` and `P_g = 2^123` are strictly below `16^32 = 2^128`, so
every leaf-to-root path has at most 32 edges. The schedule charges 32 hops
per record anyway, including shorter paths and including unused giant leaves,
which emit END. The engine that intersects the two roots is not one of those
hops; it is `R_merge` below.

The merge at an internal node is a 16-way binary heap of child heads.
Extract-min is four sift steps. The node emits the smaller head, stable on
ties by child index, and requests a replacement only from that child. By
induction the stream leaving the root is sorted by `x` and preserves
multiplicity. END travels the same edges once the child is exhausted. The
engine intersects the two sorted streams. It may stop at the first equal key;
both streams are charged in full anyway.

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

The stored record is 10 words. The transfer row charges 12 words. Placement
on an edge is an additional `128`, charged separately from the hop so the
merge price is not asked to cover movement.

Node initialisation is `18432` for a page and for an internal node: the
`16384` local-memory schedule (code distribution at both ends of the link,
zeroing, metadata, control) plus `2048` for the port-descriptor block. On a
page that block is the 256 offset queues. On an internal node it is the
sixteen port descriptors. Charged nodes are the pages plus the internal nodes,

```
N_nodes = P + I = 17014118346046923173168730371588410574.
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
the local send into the page slot. Per record,

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

```
E_b = 32*m = 5444517870735015415413993718908291383296
E_g = 32*L = 10889035741470030830827987437816582766592
E   = E_b + E_g = 16333553612205046246241981156724874149888

R_init  = 18432 * N_nodes
        = 313604229354336887927846038209117583699968
R_gen   = 23121 * (m+L)
        = 11801502908368527320605026447644869225611264
R_baby  = 64 * m
        = 10889035741470030830827987437816582766592
R_place = 128 * E
        = 2090694862362245919518973588060783891185664
R_sort  = 576 * E + 2 * 576 * I
        = 9409351897151022016303849294860281875898112
R_radix = 128 * 32 * (m+L)
        = 2090694862362245919518973588060783891185664
R_tail  = 16 * 256 * 32 * P
        = 2090694862362245919518973588060783891185664
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
multiplications. `R_sort` includes an END request and response at every
internal node. `R_gen` charges `L` giant records even though only `h < L`
are real. `R_radix` and `R_tail` charge `L` giant slots the same way.

```
T = R_init + R_gen + R_baby + R_place + R_sort + R_radix + R_tail
    + R_merge + R_setup + R_verify
  = 28068769515497374754164342230842034976561602
```

Preprocessing is setup, the baby chain, baby labels, baby pages, the baby
internal nodes, baby placement, the baby merge, and the baby radix. It is
not the giant chain:

```
P_pre = R_setup + 23121*m + R_baby + 18432*(P_b+I_b) + 128*E_b
        + 576*E_b + 2*576*I_b + 128*32*m + 16*256*32*P_b
      = 9276403576395011358628708835736690067785068
```

Peak bytes are the charged pages and internal nodes at `2^14` bytes plus the
`2^20` central reserve. The radix banks sit inside the page, so they are not
a second copy.

```
N_nodes * 2^14 + 2^20
  = 278759314981632789269196478408104519892992 < 2^143.
```

The inequalities `T < 1000*2^(1348/10)`, `(6/5)*T < 1000*2^(1348/10)`,
`P_pre < 1000*2^134` and `(6/5)*P_pre < 1000*2^134` are exact integer
statements. Equivalently, with denominator 10,

```
T^10 < 1000^10 * 2^1348
(6*T)^10 < 5^10 * 1000^10 * 2^1348
```

The factor `6/5` reprices every instruction, including those already inside
`M` and `A`. It still lands strictly under the submitted `1348/10`.

## Heuristics

None. Coverage is an integer partition of the odd residues. Page geometry and
both trees have a fixed shape. The radix is 32 exact passes over a fixed
80-byte record; it does not stop early and it does not assume a distribution
of `x`. Batch inversions run on Jacobian `Z` coordinates that the chain
already produced. Success is 1 on every valid input, which is at least 0.39.
The toy experiment checks the odd-lift identity on a small additive group and
checks that the capacity-4 port flush matches a stable sort on one 32-record
page. It is not a measurement of `T` and it is not used as a premise for the
production bound.

## Limitations

The ledger is an analytical instruction schedule. It is not a Callgrind count
and not a claim that the pages fit in any real machine. Memory on this scale
is mandatory for a generic half-table of `2^127` points together with the
page tree that addresses them. Beating the integer bound is not a practical
break and does not change the hardness assumption used by P-256. Prices `M`,
`A`, `J`, `K`, `V` and the 16-ary hop `H` are derived in this file from the
protected prime and from the local-node rule. They do not depend on another
package's totals.
