# p256-ecdlp: deterministic BSGS with projective mixed additions and block batch inversion

Status: ready for independent review; no qualification or emitted score.
Parent tree: `6423c2728dc03d5f2d647c7f22adf1dac0643a66`. This is a conservative
generic accounting control for the P-256 target, not a best-known-attack claim.
`claim.json` binds `tracks/p256-ecdlp/target-profile.json` and
`tracks/classical-group-work-v2.json` by SHA-256. No advice, no auxiliary
input, no `n` parameter, no quantum step, no reuse across targets.

## 1. Target and inputs

The protected profile fixes SEC 2 Section 2.4.2 (secp256r1): p = 2^256 -
2^224 + 2^192 + 2^96 - 1, a = -3, b = 5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b,
r = ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551 (prime,
cofactor 1), and the standard generator P. Input Q = [x]P with x uniform in
[0, r); Q = infinity is included and means x = 0. Output is x. Nothing else is
given: no powers, no interval, no distinguished-point advice.

All field elements are canonical integers in [0, p), carried as four 64-bit
limbs. Coordinates are Jacobian (X, Y, Z) with affine (X/Z^2, Y/Z^3); the point
at infinity is Z = 0.

## 2. Construction

b = 256, L = 128, m = 2^128, h = ceil(r/m) = 2^128 - 2^96 + 2^64
(340282366841710300967557013911933812736). Both m and h are fixed constants of
the profile; m*h >= r, so every x in [0, r) has x = i*m + j with 0 <= i < h and
0 <= j < m.

```text
Build:  B_0 = infinity.  For j = 1..m-1:  B_j = B_{j-1} + P  (mixed Jacobian
        addition, addend P affine and fixed).  Every BLK = 1024 consecutive
        results are buffered in Jacobian form; when a block is complete, all
        its Z values are batch-inverted, each point is converted to affine
        (x_j, y_j), and (x_j, y_j, j) is inserted into the hash table.
        B_0 = infinity is stored under an explicit infinity tag.
Search: G = [m]P (affine, computed once at setup), D = -G.  S_0 = Q.
        For i = 0..h-1: normalise S_i in blocks of 1024 exactly as above and
        look up its affine coordinates in the table.  On a hit with index j,
        return a = (i*m + j) mod r.  Otherwise S_{i+1} = S_i + D.
Final:  decode checks, [r]Q = infinity, [a]P = Q, serialise a.
```

## 3. Correctness: coverage, identity, no false stop

Coverage. Let i = x >> 128 and j = x mod 2^128. Then i < h and j < m. At giant
step i the walk holds S_i = Q - [i*m]P = [x - i*m]P = [j]P, because
x - i*m = j >= 0. That is exactly baby entry j, so the table contains S_i and
the lookup at step i hits with index j. The returned value is
(i*m + j) mod r = x mod r = x since x < r. The search never runs past step i.

No false stop. Suppose the walk hits baby entry j' at step i'. Then
S_{i'} = [j']P, and by the walk's definition Q = [i'*m + j']P. Since P has
prime order r, i'*m + j' is congruent to x modulo r, hence
a = (i'*m + j') mod r satisfies [a]P = Q, and by uniqueness of the discrete
logarithm in a group of prime order a = x. Every hit, not only the intended
one, therefore returns the correct scalar; the independent [a]P = Q check can
only confirm. No restart or amplification is needed.

Identity. x = 0 gives Q = infinity. S_0 = infinity is found under the infinity
tag at j = 0, so a = 0. Zero coordinates are never treated as a field element:
the tag is checked before the coordinate comparison.

Failure modes outside the coverage proof: only the hash table, discussed in
section 5. Its failure probability is at most 2^-128, so algorithmic success is
at least 1 - 2^-128, which exceeds the required 0.39 over uniform x and fresh
coins.

## 4. Batching schedule: where the denominators come from

The earlier affine-BSGS obstruction does not apply here, and the schedule is
explicit. The recurrence is carried in Jacobian coordinates, where an addition
never needs an inversion:

- Each walk step is a *mixed* addition (add-2007-bl): the running point is
  Jacobian and the addend is affine and fixed (P for the baby walk, D = -[m]P
  for the giant walk). The output is (X3, Y3, Z3) with Z3 = (Z1+H)^2 - Z1Z1 - HH.
  No inversion is performed inside the loop, so no step waits on another.
- After BLK = 1024 sequential additions the block buffer holds 1024 Jacobian
  points and therefore 1024 denominators Z_1..Z_1024, all simultaneously
  available. They are inverted with one Montgomery batch inversion: 1023
  prefix products, one field inversion of the total product, and two
  multiplications per element for back substitution.
- Each element is then converted once: Zinv = Z^-1, Zinv2 = Zinv^2 (1S),
  Zinv3 = Zinv2*Zinv (1M), x = X*Zinv2 (1M), y = Y*Zinv3 (1M).
- Z = 0 occurs only for the point at infinity (baby j = 0, and the giant step
  that would produce infinity, which is the answer step). Such an entry is
  excluded from the batch (its Z is replaced by 1), keeps its infinity tag and
  is charged the same conversion cost.

Initialisation, conversion, buffer traffic and the single inversion are all
charged in section 5. There is no uncharged projective-to-affine conversion and
no per-step inversion anywhere in the claim.

## 5. Instruction ledger (scalar AArch64, 64-bit limbs)

Every constant is an upper bound on executed scalar instructions, including
loads, stores, addressing and control. No primitive is unit-cost.

| Kernel | Bound | Derivation |
| --- | ---: | --- |
| A, modular add/sub | 32 | 4 ADD + 4 ADC + 4 SUB + 4 SBC, 8 compare/select, 8 load/store |
| M, field multiplication | 256 | 16 MUL/UMULH pairs (32), 16 partial accumulations (64), carry propagation (32), P-256 fast reduction: nine shifts and sixteen add/sub with carry chains (96), operand traffic (32) |
| S, field squaring | 256 | charged as M |
| I, field inversion | 86016 | fixed 4-bit window exponentiation by p-2: 256 squarings + 64 multiplications = 320 field multiplications, plus 4096 for the 16-entry window table, exponent extraction and loop control |
| C_mix, mixed addition | 3584 | 8M + 3S + 16A + 256; the 256 covers infinity/equal/opposite branch tests, case selection and loop bookkeeping |
| C_dbl, doubling (a = -3) | 2816 | 3M + 5S + 16A + 256; used only in the final scalar-multiplication check |
| N, per-point normalisation | 2304 | block batch: 3M per element (prefix product, amortised back substitution, and the I/1024 = 84-instruction share of the single inversion) + 3M + 2S for Zinv2, Zinv3, x, y; 256 for buffer traffic, the zero-Z tag and control |
| H, table insert or lookup | 2144 | 32 for the multiplicative hash (4 MUL/UMULH, carry, shift, address), 64 probes x 32 (8 key words loaded, 8 words compared, tag check, index advance, branch), 64 control |
| step | 8288 | C_mix + N + H + 256 (walk control, offset arithmetic, record store) |

Both walks use the same per-step charge: one mixed addition, one block
normalisation, one table operation, one step of control. The giant step's
addend D is affine and fixed, so it is the same 8M + 3S kernel.

Hash table. The table has 2^132 slots for 2^128 records: load factor 2^-4.
A record is 80 bytes (32-byte x, 32-byte y, 16-byte index j). Probing is linear
from h(x) = (x * mult) >> 192 with a 64-bit multiplier `mult` drawn once at
setup (charged). For a random multiplier, an insertion at load 2^-4 exceeds 64
probes with probability at most (2^-4)^64 = 2^-256; a union bound over 2^128
insertions and 2^128 lookups gives total failure at most 2^-128. Every insert
and every lookup is nevertheless charged the full 64-probe cap, so the work
bound is hard and does not depend on the multiplier. Keys are the full affine
pair (x, y), so the two points with equal x-coordinate are distinguished; a
spurious x-only match is therefore impossible.

Setup (2^24 = 16777216): profile decoding and canonicality checks, computing
[m]P by 256 doublings and 256 mixed additions, its negation, one inversion for
the affine form of D, the multiplier draw, header initialisation.

Table initialisation (20 * 2^132 = 108890357414700308308279874378165827665920):
zeroing 2^132 slots of 80 bytes at 16 instructions per 64-byte line. This is
memory traffic and is charged, not assumed free.

Verification and output (20447232 = 4*256*C_mix + 2^24): decoding Q, the
[r]Q = infinity check, [a]P = Q by double-and-add, scalar reconstruction and
serialisation, with slack for the curve and encoding checks.

## 6. Exact totals

| Component | Instructions |
| --- | ---: |
| baby walk, m = 2^128 steps | 2820260257040737985184448746394494936547328 |
| giant walk, h steps | 2820260256384094974419112531302107439955968 |
| table initialisation | 108890357414700308308279874378165827665920 |
| setup | 16777216 |
| verification and output | 20447232 |

Total work T = **5749410870839533267911841152074768241393664** instructions.
Preprocessing (setup, table initialisation and the whole baby walk, no reuse
credit) = **2929150614455438293492728620772660780990464** instructions.
Peak memory = 80 * 2^132 + 2^24 = **435561429658801233233119497512663327440896**
bytes; nonuniform advice = 0 bytes.

Dividing by the fixed reference price C = 1000 instructions per
target-group-operation equivalent:

- time_log2 = log2(T/1000) = **132.078610** (132078610/1000000, rounded up)
- preprocessing_log2 = log2(prep/1000) = **131.105679** (131105679/1000000)
- memory_log2_bytes = **138.321929** (138321929/1000000)

All table access, initialisation traffic, failure paths, verification and
output are inside T. Memory has no scalar contribution and is reported, not
scored.

## 7. Heuristics and limitations

Heuristic premises: none for correctness, coverage or work. The only
probabilistic element is the hash multiplier, and it is bounded above by
2^-128 and folded into the success probability, so no expected-value premise
enters the work bound. There is no birthday argument, no cycle-escape
assumption, no distinguished-point density assumption and no interval
assumption.

Limitations, stated plainly. The 2^138.3-byte table is astronomically
infeasible; it is mandatory for this construction and is not a practical
attack. The bound is a conservative generic control: Pollard rho with
distinguished points, automorphism or negation-map quotients and kangaroo
variants, once each is given its own complete success and cost argument under
this same C = 1000 normalisation, will price below it. Beating this control
alone does not establish novel cryptanalysis or a weakened hardness assumption.

## 8. Evidence

`experiment.py` runs the *same* update on the real P-256 curve at toy size: the
same Jacobian mixed-addition kernel, the same block batch inversion with the
same prefix-product/back-substitution schedule, the same open-addressed lookup
with the same probe cap, and the same coverage loop, with m = ceil(sqrt(n))
instead of 2^128. It recovers every scalar in [0, n) exhaustively, checks each
recovered value against [a]P = Q, and includes x = 0. That is toy-scope
evidence for the update and the coverage argument; it is not evidence for the
production instruction counts, which are analytic bounds.

References: HAC section 3.6 (baby-step giant-step and generic DLP bounds);
add-2007-bl and dbl-2001-b formulas as catalogued by the Explicit-Formulas
Database. The construction and ledger above are self-contained; review does not
need to fetch them.
