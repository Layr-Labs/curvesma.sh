# Signed BSGS with an itemised distributed merge network

This candidate claims the conservative rational upper bound **time_log2 = 1305909/10000 = 130.5909** in protected target-group-operation equivalents, with **success_probability = 1/2**, **memory_log2_bytes = 137**, **preprocessing_log2 = 130.5909** (the same bound, applied to the complete attack) and **nonuniform_advice_log2_bytes = 0**.

These are analytic scalar-instruction bounds, not a production execution, not a timing measurement and not a qualified result. The score is `log2(T/1000)` where `T` is the exact integer instruction total derived below and `C = 1000` is the protected fixed divisor for this track. Reference group additions are charged at the measured six-multiplication calibration `6*179 = 1074` instructions each, which is more than their protected price of 1000; every other instruction is charged separately and explicitly. There is no free address space, no free code distribution, no amortisation over other targets, no randomness other than the fixed public constants, and no heuristic.

## Target

The protected secp256k1 target is used unchanged:

```
p = fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
r = fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
curve: y^2 = x^3 + 7 over F_p, cofactor 1
P = (79be667e..., 483ada77...)   the standard generator
challenge: Q = [x]P with x uniform in [0, r), only Q given
```

The attack is one cold classical attack. It performs no precomputation for any other target and receives no advice beyond the public profile constants and Q itself.

Parameters:

```
m = h = 2^127 = 170141183460469231731687303715884105728
d = 2m - 1   = 340282366920938463463374607431768211455 = 2^128 - 1
D = -[d]P
B = 2^16     chains per list
F = 64       fan-out of each merge tree
L = 22       levels; 21*6 = 126 < 127 <= 132 = 22*6
```

## Construction

A **record** occupies 8 words (64 bytes) in storage: one tag word (1 for the point at infinity, 0 otherwise), four 64-bit limbs holding the canonical x-coordinate in [0, p), two words holding the index (j for a baby record, i for a giant record), and one padding word that keeps the storage stride at 64 bytes. Seven words carry information; the padding word is never transmitted on a link and is charged only in storage.

The two lists are

- baby records `(x([j]P), j)` for `0 <= j < m`;
- giant records `(x(Q + [i]D), i)` for `0 <= i < h`.

A finite point's merge key is `(0, x)` compared lexicographically: tag first, then the four x limbs starting at the most significant limb. Infinity has tag 1 and x = 0. The baby record at j = 0 is the point at infinity and is the only tag-1 baby record, since [j]P is finite for every 1 <= j < m = 2^127 < r. The giant record at i is infinity for at most one i, because i*d mod r is injective on 0 <= i < h as h*d < r; when it occurs it is the only tag-1 giant record.

Both lists are written into two complete 64-ary merge trees, one per list, with m leaves each and one record per leaf. Each tree's root produces a nondecreasing key stream. A single central engine intersects the two streams: it caches one head from each root, discards only the smaller head, and on equal keys performs the match work. On an equal key the giant point `Q + [i]D` equals `[j]P` or `-[j]P`; the engine forms `a = i*d + j mod r` or `a = i*d - j mod r`, verifies `[a]P = Q` inside the fixed reserve, and outputs the canonical 32-byte scalar. Exhausting either stream returns failure.

Duplicate giant keys are possible and are not dropped. Internal nodes emit the left child's head first on equal keys and mark only that head consumed; the other head stays cached and is emitted later. Every record of both lists is emitted exactly once in nondecreasing key order, so an equal-key run is never skipped and never requires more than one final scalar verification.

## Exact coverage

Work modulo r. For a fixed giant index i, the giant key equals a baby key exactly when

```
x = i*d + j  or  x = i*d - j   (mod r),   0 <= j < m.
```

As j runs over [0, m), the residues `i*d + j` form the consecutive integers `[id, id + m - 1]` and the residues `i*d - j` form `[id - m + 1, id]`. Their union is the integer interval `[id - (m-1), id + (m-1)]` of `2m - 1 = d` consecutive integers.

Consecutive giant intervals are adjacent and disjoint because the next one starts at `(i+1)d - (m-1) = id + (m-1) + 1`. The union over `0 <= i < h` is therefore the single integer interval

```
[-(m-1), (h-1)d + (m-1)]
```

with exactly `h*d` members. Since `h*d < r`, reduction modulo r is injective on that interval, including across zero. Exact integers:

```
COV = h*d = 57896044618658097711785492504343953926464851149359812787997104700240680714240
2*COV - r = 92138019644721193389046258963199934143 > 0
```

so `COV/r > 1/2`. The construction is deterministic: no coins, no retries, no amplification. It deliberately fails for the uncovered scalars, which are charged in full because the complete drain is charged whether or not a match occurs.

## Generation: the protected reference operation, charged at 1074

Each list is generated in `B = 2^16` contiguous chains. Chain c covers the contiguous index range `[c*m/B, (c+1)*m/B)`, of length `m/B = 2^111`, and starts from the affine point `[c*m/B]P` (baby list) or `Q + [c*m/B]D` (giant list). All chains advance in lockstep, one affine addition per chain per round.

Each round's `B` independent affine additions are performed with one Montgomery batch inversion: form the `B` prefix products of the denominators, invert the final product once with a single field inversion, and recover each inverse in a backward pass by multiplying the running inverse by the preceding prefix. The protected reference operation for this track is exactly this object — "batched incomplete affine addition, 6M convention, with S = M and A = 0" — priced at `C = 1000` target-group-operation equivalents, with the batch work (3M of the 6M) inside that price and **no standalone inversion charged per addition** (tracks/docs/COST_MODEL.md).

The measured six-multiplication calibration for this track is `6*179 = 1074` instructions, which the organizer rounds down to 1000. This submission charges the un-rounded **1074 instructions per addition**, so it is strictly more expensive than the protected price and no reference operation is double counted: the 1074 covers the formula's three multiplications, the batch-amortised three, the additions, the inversion amortised over `B = 2^16` (one inversion of 2^19 instructions is 8 instructions per addition, inside the 3M batch allowance), the chain-state updates and the loop work.

Exceptional steps are handled by a charged complete fallback, never by applying the incomplete formula to a doubling or to infinity. A running point equals `+` or `-` its fixed increment at most twice per chain: both are integer ranges modulo the prime r, so `[t]P = [u]P` forces `t = u` and `[t]P = -[u]P` forces `t + u = r`, each with at most one solution per chain per increment because the increment is a fixed nonzero multiple of P and the chain covers a contiguous index interval. The bound is `2` exceptional steps per chain, `2^18` in total for both lists, each charged a complete Jacobian point operation of `JAC = 16*512 + 32*64 + 4096 = 14336` instructions (the field-operation allowances are itemised: `M = 512` for a four-limb modular multiplication, `A = 64` for a modular addition, and 4096 for the branch, zero-test and copy work of a complete operation).

Chain start points: `2^16` per list, `2^17` in total, each a scalar multiplication charged at `512` complete Jacobian additions, i.e. `512*14336 = 7340032` instructions. This is charged outside the fixed reserve, so no start-up work is hidden.

Total reference additions: `m + h = 2^128`, each charged 1074 instructions.

## Machine model

The machine is the protected scalar AArch64 model: 64-bit general-purpose registers; every arithmetic, logical, load, store, address and control instruction counts; no vector units, no free oracles, and no flat global address space. The storage is the distributed merge network described below: leaves, internal nodes and the central engine are separate scalar processors with local memory and explicit ports. All work of all processors is summed. A record is never accessed by a wide address in one instruction; it is moved word by word along links, and every word transfer, address update and control instruction is charged.

The two trees are complete `F = 64`-ary trees with `m` leaves each. A leaf is a scalar processor holding one record, its address latch, two ports and flags; it answers a NEXT request with its record and then with END. An internal node holds `F = 64` cached head slots of 8 words each, `F` input ports, one output port, a tournament over its `F` heads, and routing metadata. Nodes hold no curve arithmetic; all curve work happens in the generating engine and the central engine.

The engine addresses a leaf by 127 bits, `log2 F = 6` bits per level, so a route is `L = 22` hops (21 levels of 6 bits fully used, one level carrying the remaining bit). No intermediate node performs a routing-table lookup: it forwards the record on the next 6-bit field of the destination it received.

A node holds `576` words of state (512 words of cached heads, 64 words of port state, flags and routing metadata) and a `512`-word read-only code image of 1024 AArch64 instructions, so 8704 bytes per node. Node initialisation is charged at 4 instructions per state word, 4 instructions per code word for the image transfer and store at both endpoints, and 512 fixed control instructions: `576*4 + 512*4 + 512 = 4864` instructions per node, covering the clears, the initial slot flags, the END markers, the routing metadata, the code delivery and the start-up acknowledgements. Code distribution is therefore not free and no node is assumed to receive a program without paying for the transfer. The number of nodes is charged as `ceil(2m/63)` per list, which is twice the exact count `(m-1)/63` of internal nodes of a complete 64-ary tree with m leaves; both lists together are charged `10802614822886935348043638331167244810` nodes.

### Record write, itemised per hop

The generating engine pushes each record down its fixed route, one hop per level. Per hop:

| Item | Instructions |
| --- | ---: |
| sender loads the 7 information words | 7 |
| sender stores them into the link port | 7 |
| receiver loads them from the port | 7 |
| receiver stores them into its buffer | 7 |
| address update at each end, 2 per word | 14 |
| control: 6-bit field extract (2), child port address (2), handshake (1), credit write (1), loop counter and branch (2) | 8 |
| **hop total** | **50** |

A write is `L*50 = 1100` instructions. Every one of the `2m` records is written exactly once. The padding word is not transmitted; it exists only so leaf storage has a 64-byte stride, and the storage bytes are charged in the memory bound.

### Record read for the merge, itemised per hop

A record moves up one hop when its node emits it, and the receiving node performs its tournament update. Per hop:

| Item | Instructions |
| --- | ---: |
| the transfer above | 50 |
| six tournament key comparisons (`log2 64 = 6`), 16 instructions each: 2 tag-word loads (2), CMP (1), conditional branch on the tag (1), 8 x-limb loads (8), 4 SUBS/SBCS forming the borrow chain (4) | 96 |
| six sift pointer steps, 2 instructions each: shift the 6-bit slot index (1), select the child/sibling slot (1) | 12 |
| slot state: mark the consumed slot empty, mark the refilled slot full, decrement the filled count, set the END flag | 4 |
| ACK/credit: write the credit word, set the ready flag, address update, branch | 4 |
| **hop total** | **166** |

A read is `L*166 = 3652` instructions. Every one of the `2m` records is read exactly once. The head copy into the node's cached slot is not charged twice: the receiving buffer of the transfer *is* the cached head slot, so the store in the transfer table is the head copy. Early exit from the key comparison is permitted only inside this fully charged comparison; the bound assumes the worst case of all four limbs differing at the last limb.

A complete drain of both streams upper-bounds the actual work, because the engine may abort at the first match. Cached lookahead on an aborted run therefore cannot exceed the budget. END markers are counted as records' worth of control in the per-hop charges and in the per-node initialisation; there is no zero-cost sentinel.

### Root intersection

The two root streams are intersected inside one charged engine holding one head from each root. Before each advance it compares the two cached heads (16 instructions, the comparison itemised above), moves the winning or discarded 8-word record between two local buffers (8 loads and 8 stores, 16 instructions) and updates its own state (4 instructions): `INTER = 36` instructions per head consumed. At most `2m` heads are consumed, because the drain is complete; the charge is `2m*36`.

### Generation outer work

Per generated record, the engine's block counter, leaf-address computation, chain-advance control, block-normalisation bookkeeping and loop branch are charged `OUTER = 8` instructions, `2m*8` in total.

## Fixed reserve

`2^40 = 1099511627776` instructions cover, all with margin: input parsing and canonical checks (curve equation, on-curve, `[r]Q = O`, duplicate-key rejection); one scalar multiplication for `[r]Q` and one for `D = -[d]P`, each bounded by 512 complete Jacobian additions at 14336 instructions; preparation of the two chain increments and of the batch-inversion constants; the reconstruction of `a` and one independent `[a]P = Q` verification, including the multi-limb arithmetic modulo r; canonical output; and fixed cleanup. The chain start points are **not** in this reserve: they are charged as their own ledger line, so no work is counted twice and none is hidden.

## Exact ledger

Every term is an exact integer. `T` is the sum of all of them; the reference part is `1074*(m+h)`.

| Component | Instructions |
| --- | ---: |
| reference additions, 1074*(m+h) | 365463262073087909759664328381719059103744 |
| record writes, 2m*1100 | 374310603613032309809712068174945032601600 |
| record reads, 2m*3652 | 1242711203995267268568244066340817508237312 |
| node initialisation, 10802614822886935348043638331167244810*4864 | 52543918498522053532884256842797478755840 |
| root intersection, 2m*36 | 12250165209153784684681485867543655612416 |
| generation outer work, 2m*8 | 2722258935367507707706996859454145691648 |
| exceptional fallbacks, 2^18*14336 | 3758096384 |
| chain start points, 2^17*512*14336 | 962072674304 |
| fixed reserve | 1099511627776 |
| **total T** | **2050001412324430834062893202469342222401024** |

The diagnostic logarithm is `log2(T/1000) = 130.5908206043`. The submitted bound is the four-decimal rational `1305909/10000 = 130.5909`, justified without floating-point rounding by the exact integer inequality

```
T^10000 < 1000^10000 * 2^1305909
```

which holds, and which fails one step of 1e-4 lower:

```
NOT ( T^10000 < 1000^10000 * 2^1305908 ).
```

Both sides are exact Python big integers. `T < 1000 * 2^130.5909` may be read directly from the same inequality.

## Memory, preprocessing and advice

Peak storage is charged in full, including code, node state, records and scratch:

```
leaf records        2m * 64 = 21778071482940061661655974875633165533184 bytes
node state + code   10802614822886935348043638331167244810 * 8704
                            = 94025959418407885269371828034479698826240 bytes
engine, batch scratch, output buffers
                            =    16777216 bytes
peak                        = 115804030901347946931027802910112881136640 bytes
```

`log2(peak) = 136.4107`, so `memory_log2_bytes = 137`. No second copy of either list is retained; records are written into their leaves and streamed out.

Preprocessing — storage construction, node initialisation, chain start points, and the complete baby list — is a subset of the charged attack, so the conservative `preprocessing_log2 = 130.5909` bound is justified by the same `T`. No work is amortised over another target. Nonuniform advice is zero bytes: only the public profile constants and Q are read, and the reported logarithm is the required 0.

## No heuristic, no expectation used as a bound

There is no random-walk mixing assumption, no collision-rate or birthday estimate, no random oracle, no ideal hash table, no empirical success extrapolation and no uncharged precomputed table. Success is a **counted set**, not a probability model: the covered residues are exactly `h*(2m-1)` and `2*COV - r > 0`, so success is at least `1/2` for every uniform challenge scalar, with no retry and no amplification. Failure outside the covered set is stated and its full cost charged. Every bound above is a worst-case instruction count: the comparison costs assume the worst case with no early exit, the read and write traversals charge the complete drain rather than the work up to a hoped-for match, and the exceptional-step count is a deterministic bound on the number of solutions of `[t]P = +-[u]P` in a contiguous index interval. No expected value is used in place of a maximum anywhere in the ledger. The only probabilistic object in the document is the challenge scalar's own stated uniform distribution.

## Limitations

The field-operation prices, node capacity, hop charges and routing model are analytic scalar-instruction bounds, itemised above, not compiled measurements; the only measured figure used is the protected 179-instruction field multiplication, from which the 1074 reference price comes. The construction is an accounting object: the peak storage of about `2^136.4` bytes is not physically available, and no production attack was run. The companion experiment checks bounded instances of the algebra, the F-ary merge with duplicate keys and left-first ties, Montgomery batch inversion against direct inversion, the exact coverage count on the production integers and the exact score inequality; it is evidence at its executed scope only, and it is not a production complexity proof or an instruction-count measurement.
