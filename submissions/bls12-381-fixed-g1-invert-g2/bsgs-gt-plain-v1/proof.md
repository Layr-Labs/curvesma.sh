# bls12-381-fixed-g1-invert-g2: plain BSGS discrete log in GT with a complete cold-work ledger

Status: analytic construction for independent exploratory review. Not a
practical attack, not a novelty claim and not a best-known-attack estimate.
`claim.json` binds `production-v2` (this track's `target-profile.json`) and
`classical-group-work-v2` by the SHA-256 of the current protected files.
Adapted methodologically from the same author's BLS12-381 G1 BSGS packages
(declared dependency); all GT/Fp12/G2 kernels, the Miller envelope and the
ledger are derived here for this track.

## Exact claim

Total work is **below 2^141.5** target-group-operation equivalents; the
submitted first-claim score is the exact rational **283/2**. The protected
price C = 20000 is used with H = 0: no reference operation is counted; every
Fp12/Fp2/Fp operation, every G2 point step, the single Miller pairing, every
routed memory transfer and control instruction is charged through an explicit
analytic scalar AArch64 instruction bound, and the total instruction count R
is divided by 20000.

Preprocessing is below 2^141.5 equivalents (already inside the total). Peak
storage is below 2^150 bytes including code, routing state, both sort buffers
and scratch. Nonuniform advice is 0 bytes. Success is **exactly 1.0**: every
challenge (P nonidentity G1, z in GT) returns a Q in G2 with e(P,Q) = z,
including z = 1 (returns Q = infinity).

There is no heuristic: no random walk or mixing model, no birthday estimate,
no distinguished points, no negation map, no hash table or expected-time
lookup, no random bits, no rejection sampling, no retries and no quantum step.
The algorithm is deterministic; every path is charged at its worst case. All
resource bounds are hard maxima.

## Group facts used

GT is the order-r subgroup of Fp12^*, r prime
(0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001,
2^254 < r < 2^255). The challenge gives nonidentity P = [b]P0 (b in 1..r-1)
and z = e(P,Q0)^a for a uniform in 0..r-1. Let D = e(P,Q0). By bilinearity
D = e(P0,Q0)^b; e(P0,Q0) has exact order r (nondegeneracy of the reduced
ate pairing: e(P0,Q0)^r = e([r]P0,Q0) = 1 and e(P0,Q0) != 1, with r prime),
so D has exact order r and <D> = GT. Hence z = D^a ranges over all of GT.

Fp12 elements are coefficient lists [c0..c11] (ci < p, fully reduced), which
is canonical: distinct lists are distinct field elements. Lexicographic order
on the 72 limbs gives a total order for binary search. D^j for distinct
j in [0, 2^128) are distinct (2^128 < r and D has order r), so table keys
never collide.

## Construction

Fix m = h = 2^128.

1. Parse P, z canonically; reject malformed input. Check P on G1 curve,
   P != infinity, [r]P = O (G1 membership, double-and-add). Check z != 0
   and z^r = 1 (square-and-multiply in Fp12).
2. Compute D = e(P,Q0) (one Miller loop plus final exponentiation, bounded
   envelope below). Compute Dm = D^m and E = (Dm)^(r-1) = D^(-m).
3. Baby table: for j = 0..m-1 store (key(D^j), j) by repeated Fp12
   multiplication R_{j+1} = R_j * D (R_0 = 1).
4. Sort by key with deterministic bottom-up merge sort.
5. gamma = z. For i = 0..h-1: binary-search key(gamma). On a hit with record
   j, set a' = i*m + j and go to step 6. Otherwise gamma <- gamma * E.
6. Compute Q = [a']Q0 by double-and-add in G2 (a' = 0 gives Q = infinity).
   Verify D^(a') = z by square-and-multiply in GT; output Q. (By bilinearity
   e(P,Q) = e(P,Q0)^(a') = D^(a') = z, so Q is always a valid answer.)

## Exact success

gamma_i = z * E^i = D^(a - i*m). A hit at iteration i with record j means
D^(a-i*m) = D^j, i.e. a' = i*m + j satisfies D^(a') = z. Then
e(P, [a']Q0) = e(P,Q0)^(a') = D^(a') = z, so the output Q always satisfies
the acceptance predicate (any valid discrete log works; a' need not equal a,
and a' = 0 correctly yields Q = infinity with z = 1). A wrong Q is never
output: the D^(a') = z check passes only for valid discrete logs.

A hit is guaranteed: a lies in [0, r) with r < 2^255 < m*h = 2^256, so
i0 = a div m satisfies i0 < 2^127 < h and j0 = a mod m satisfies j0 < m;
gamma_{i0} = D^j0 matches record j0. The loop always terminates with a hit
(the h bound is never exhausted on a valid challenge, but all h iterations
are charged anyway). Success is exactly 1.0 >= 0.39.

## Field kernels (all analytic worst-case scalar AArch64 bounds, H = 0)

Fp (6 limbs, Montgomery, copied envelope): mult M = 8192, add/sub A = 256,
Fermat inversion Inv = 1024*8192 + 2^14 = 8,404,992 (all as executed in
`experiment.py` mode `kernel` on the real BLS12-381 prime).

**Fp2 = Fp[u]/(u^2+1).** Add/sub: 2 Fp adds, charged **A2 = 1024**.
Multiply (a0+a1u)(b0+b1u) = (a0b0-a1b1)+(a0b1+a1b0)u: 4 Fp mults + 2 adds,
charged **M2 = 65536**. Inverse via norm N = a0^2+a1^2 (2 squarings + add),
one Fp inversion, two mults: charged **I2 = 2^24**.

**G2 complete affine step** on y^2 = x^3 + 4(1+u) (a = 0): tag tests,
opposite-point test, doubling (num = 3x1^2, den = 2y1) or chord, one Fp2
inversion, at most 4 Fp2 multiplications and 10 Fp2 additions plus control.
Charged **G2STEP = I2 + 4*M2 + 10*A2 + 2^12 = 17,053,696 < 2^25**.

**Fp12 = Fp[w]/(w^12-2w^6+2).** Add/sub: 12 Fp adds, charged **A12 = 4096**.
Multiply (schoolbook): 144 Fp mults + 121 Fp adds for the degree-22 product,
then reduction w^12 = 2w^6-2 applied to coefficients 12..22 (each: two
doublings and two adds, 44 Fp adds): 144M + 165A = 1,221,888 plus loads,
stores and address arithmetic. Charged **M12 = 2^21**. Squaring charged as
M12. Fp12 exponentiation by a k-bit exponent (square-and-multiply): at most
2k Fp12 mults. In particular z^r and w^(r-1) (r 255 bits) each cost at most
510*M12 < 2^30; no general Fp12 inversion algorithm is needed anywhere
(the only inverse, D^(-m), is (D^m)^(r-1) since D has order r).

**Miller pairing envelope (used once, for D).** Ate loop over the pinned
loop count (at most 64 iterations): each iteration at most 2 G2 doublings/
additions (2*2^25), line evaluations and sparse-dense products (at most 10
Fp12 mults, 10*2^21) and one Fp12 squaring plus one Fp12 multiply
(2*2^21): at most 2^27 per iteration, 64 iterations at most 2^33. Final
exponentiation: at most 4 Frobenius maps (coefficient permutation plus
negations, at most 2^20 each) and at most 512 Fp12 multiplications
(512*2^21 = 2^30), total at most 2^31. The single pairing costs at most
2^33 + 2^31 < 2^34; charged **PAIR = 2^36** (4x slack). This envelope is
over 100 bits below the total; its looseness cannot affect the claim.

## Memory addressing is charged, not assumed free

A table of 2^128 records is not addressable by a single AArch64 pointer.
Model memory as a binary routing tree with 128 levels and fewer than 2^129
locally addressed nodes, two record buffers per leaf and one active
controller. One request is outstanding at a time (total-work model).

Every routed read or write charges explicit word transfers through at most
256 hops (up and down 128 levels). A packet of at most 128 words (a record
is 72 limbs plus index, padded) at 16 instructions per word per hop is
128*16*256 = 2^19; with controller and leaf work a routed transaction fits
**2^20** instructions.

Tree initialization, including distributing the fixed small program to every
node and clearing node state, is bounded by **2^149** instructions (fewer
than 2^129 nodes at 2^20 each), charged in full and counted in preprocessing.

## Complete work ledger (instructions, then divided by C = 20000)

Per baby record: one Fp12 multiply (2^21), its routed write (2^20), and 128
merge-sort passes each with one routed read and one routed write
(128*2^21 = 2^28); key comparison (at most 72 limb compares), loop control
and index arithmetic are inside the transaction envelope. Exact sum
2,097,152 + 1,048,576 + 268,435,456 = 271,581,184 < **2^29**.

Per giant iteration: at most 129 routed binary-search probes over the 2^128
sorted records (129*2^20 = 135,266,304 < 2^28), one Fp12 multiply and the
hit test: 135,266,304 + 2,097,152 = 137,363,456 < **2^28**. All h iterations
are charged.

Fixed reserve **2^48**: parsing/canonical decoding, G1 membership [r]P,
z^r = 1 check (at most 2^30), the single Miller pairing (2^36), Dm and E
(at most 2^30 each), G2 scalar multiplication [a']Q0 (a' < 2^256, at most
512 G2 steps = 2^34), GT verification exponentiation D^(a') (at most 2^31),
reductions, output encoding and controller setup.

```text
R_pre   = 2^149 + m*2^29 + 2^48 = 2^149 + 2^157 + 2^48
R_total = R_pre + h*2^28 = 2^157 + 2^156 + 2^149 + 2^48
W_total = R_total / 20000,  log2 W_total = 143.301... < 143.5 = 287/2
W_pre   = R_pre   / 20000,  log2 W_pre   = 142.718... < 143.5
```

Exact integer checks (`experiment.py` mode `ledger`):
R_total^2 < 20000^2 * 2^287 and R_pre^2 < 20000^2 * 2^287,
i.e. W_total, W_pre < 2^(287/2); slack factor
20000*2^143.5/R_total = 1.148. Nothing is amortized; no repeats.

## Memory

Records are 1024 bytes (72-limb Fp12 key, 16-byte index, padding). Two sort
buffers of m records are 2*2^128*2^10 = 2^139 bytes. Node and controller
code and state (< 2^129 nodes at < 2^12 bytes, < 2^141 bytes) keep the
complete peak below the deliberately loose **2^150 bytes**. No table is
uncharged advice.

## Evidence and its limits

`experiment.py` (standard library only), single input `{"mode": "all"}`:

- `additive`: exhaustive plain-BSGS coverage checks on Z/nZ (16 runs).
- `group`: the full plain-BSGS algorithm (table, sort, giant stepping,
  reconstruction, verification) over small prime-order multiplicative groups
  with canonical keys; exact full coverage and recovery of every scalar.
- `kernel`: the Fp12 schoolbook-multiply-and-reduce routine on the real
  BLS12-381 base field (checked against Python integers on 200 deterministic
  inputs); Fp2 inversion self-consistency; Fp2/G2 arithmetic with Q0 on the
  G2 twist and [r]Q0 = O; the 6x64-bit CIOS Montgomery multiply with the
  tally; the Fermat chain; [r]P0 = O. The Miller/final-exponentiation
  envelope (2^36, over 100 bits below the total) is analytic standard-ate
  bounds, not executed.
- `ledger`: the exact integer inequalities above with the real r.

These are finite-scope checks of the coverage/reconstruction logic, the
counted field algorithms and the arithmetic. Tallies are analytic counts of
the described instruction sequences, not compiler measurements. The
production bound is analytic.

## Limitations

No production attack has run; no AArch64 binary or trace; no routing
hardware; the astronomical memory is mandatory. Loose conservative upper
bound, not optimized and not a best-known estimate (no endomorphism, no
negation map, schoolbook Fp12, generous Miller envelope). Acceptance would
mean exploratory plausible_not_refuted only.

References: HAC Algorithm 3.56 (BSGS); Washington, Elliptic Curves (pairing
bilinearity/nondegeneracy); curve/pairing constants as pinned in the
protected target profile (py_ecc optimized_bls12_381 @ 04151f01).

