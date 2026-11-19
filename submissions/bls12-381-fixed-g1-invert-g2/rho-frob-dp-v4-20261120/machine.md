# Finite scalar implementation and work bounds

These are **analytical scalar AArch64 instruction-equivalent bounds**, not
timings, on one processor with 64-bit scalar registers. Python is an
executable mathematical reference, not the priced implementation.

## The reference operation (one unit per walk step)

The protected cost model (tracks/docs/COST_MODEL.md) defines this track's unit
as the pinned blst dense Fp12 multiplication `blst_fp12_mul`, in blst's tower
representation (Fp2/Fp6/Fp12). Its count includes reduction, internal
additions and memory work: measured at 22,996 instructions, rounded to
C = 20,000. "A reference operation costs one group-equivalent unit", and its
implementation is not double-counted.

Each walk step, and each fruitless-cycle escape (a squaring), is exactly one
call `blst_fp12_mul(out, X, T)` on canonical GT elements. That call is
charged 1 unit; everything else is charged by instruction count divided
by C.

**Fp multiplication price.** The same document records the same pinned blst
build at **383 instructions per `blst_fp_mul`** (six-limb Montgomery
multiplication). For the few per-step Fp products in the frame computation,
we charge **512 instructions each**, which includes call overhead. We also
charge an explicit final conditional subtraction (48 instructions), so we do
not rely on blst's output range. Our own looser analytical bound (5,728,
arithmetic.c) is used for all one-target fixed work below.

The tower is Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-(1+u)),
Fp12 = Fp6[w]/(w^2-v). Its w is the protected w (u = w^6 - 1).
`poly_to_tower` maps the protected basis to it, and the experiment checks
products agree.

## Per-step non-reference work: E <= 12,288 instructions

| Work | Bound |
|---|---:|
| Canonicalize all 12 Fp coefficients of the product: at most 9 conditional subtractions of p each (2^384 < 10p), <= 48 instructions each | 5,760 |
| Key t0 = coefficient 0, which is already canonical: index = word0 & 255, DP test on bits 64..173 (two masked words), branch | 32 |
| Frame: V_0 = copy of Re c_{1,0}, V_3 = Im c_{1,0} (copies). V_1, V_2, V_4, V_5 are one Fp multiplication by a fixed constant each, because gamma[i][2] is real or purely imaginary (printed in development; follows from v^(p^i-1) = xi^((p^i-1)/3)). That is 4 x (512 + 48) | 2,240 |
| Frame: argmin of six 6-word values (5 comparisons, <= 32 each) | 160 |
| Frame: W = Re(sigma^i0(c_{0,1}) gamma[i0][1]), at most 2 Fp multiplications, 1 subtraction mod p (<= 96) and canonicalization; compare W with (p-1)/2 (<= 32); select i*; record pointer i*·256+j | 2*(512+48) + 96 + 32 + 32 = 1,280 |
| Labels: a += T_c[i*][j], b += T_d[i*][j] mod r (four-limb add_mod, <= 256 each, loads included) | 512 |
| Ring buffer (64 entries): compare key word0 with 64 stored word0s (3 instructions each); a full 6-word compare only on word0 equality; store t0 (6 words), X (72 words), a, b (8 words), i* (1 word) into the next slot (2 instructions per word) | 192 + 200 = 392 |
| Step counter: three-word increment and compare with S | 16 |
| Pointer swap, blst call/return, spills, loop control | 256 |
| **Total** | **10,648 < E = 12,288** |

The sign handling needs no extra multiplication: Re(sigma(c)*(f*u)) = b*f
when c = a + b*u, and Re(c*(e)) = a*e. The table of 3,072 twisted records,
T[i][j] = pi^{-i}(M_j) with 576 bytes of element plus 64 bytes of scalars
each, occupies 1,966,080 bytes.

**Escape** (rate about 1/3072 per step; also triggered by a DP revisit with equal canonical labels): select the minimal key among at most
65 ring entries (65 x 32 instructions), double two scalars mod r (512), make
one blst_fp12_mul (1 unit, counted as a step) and clear the ring (16).
Fewer than 4,096 instructions beyond a normal step. Every escape is counted
as one step of the cap S. Escapes beyond S/256 count as failure; the expected
number is about S/3072, so this limit is exceeded with negligible probability.
The escapes' non-reference work is therefore at most S/256 * 4,096 = 16*S
instructions, i.e. **64 instructions per step amortized**. This fits in E's
slack: 10,648 + 16 = 10,664 < 12,288. The ring is not cleared on escape.

**A genuine relation at a ring match:** two four-limb multiplications by
lam^{i} mod r on each side, a comparison, then one inverse mod r (< 2^22)
and verification g^x = z. This happens at most 64 times (otherwise failure),
covered in F.

## Distinguished points (inside F)

On a DP, which happens with probability 2^-110 per step:
- compute the canonical labels (2 four-limb Montgomery multiplications by
  lam^{i*}, < 2^14);
- probe an open-addressing table of 2^21 slots of 128 bytes (full 48-byte
  key t0, 64 bytes of labels, occupancy), with at most 2^10 probes. Each
  probe compares up to 6 key words, <= 64 instructions.

That is < 2^17 instructions per DP event. DP events of every kind (new,
fruitless revisit or genuine match) are counted, and the 2^22 + 1st event
ends the run as a failure. The total is therefore < 2^22 * 2^17 = 2^39,
whatever the walk does. At most 64 verifications g^x = z are made, at
<= 510 dense products each with our 2^21 bound: < 2^36.

## Fixed work F = 2^42 instructions (one target, charged in full)

* **Input.** Decoding and curve checks, z^r = 1: < 2^31.
* **Pairing g = e(P,Q0)**, following the protected py_ecc algorithm with
  our dense product G <= 2^21:
  - twisting and casting: <= 64 product-equivalents;
  - Miller loop: <= 64 iterations, each with tangent and chord lines,
    numerator/denominator updates, and projective doubling and addition,
    < 200 Fp12 operations: < 12,800 products;
  - Fermat inversion of f_den: <= 9,144 products;
  - final exponent (< 4,320 bits): <= 8,640 products.

  That is < 2^15 products, < 2^36 instructions. Only the function's value
  is used, and it equals the protected value on this exact P.
* **Basis change and Montgomery conversion of g and z:** < 2^20.
* **Randomness:** 514 scalars, at most 1024 draws of 255 bits each, at 64
  instructions per bit: < 2^34.
* **Records:** 256 records and the start, each two exponentiations
  (<= 510 products) and one product, with our 2^21 product bound:
  < 257 * 1021 * 2^21 < 2^39.1.
* **Twisted table.** 11 Frobenius images per record, each 6 Fp2
  multiplications by constants plus conjugations (< 2^17 instructions per
  image with our bounds): 2,816 images, < 2^29. Label scaling: 6,144
  multiplications mod r, < 2^26. The Frobenius constants
  gamma[i][n] = xi^(n(p^i-1)/6) (reduced mod p^2-1) can be computed by 72
  Fp2 exponentiations of <= 762 bits, < 2^36, or loaded as public constants.
* **DP handling and verifications:** < 2^39 + 2^36.
* **Genuine ring relations:** <= 64 x (2^22 + 510 x 2^21) < 2^37.
* **Output.** Affine G2 double-and-add for [x]Q0 and the optional
  [r]Q = infinity check: < 2^35. The optional [r]P check: < 2^34.
* **Fixed overhead:** parsing, constants, allocation, code loading and
  shutdown: < 2^32.

The sum is < 2^41 < F = 2^42.

## Peak storage

| Item | Size |
|---|---:|
| DP table, 2^21 x 128 bytes | 256 MiB |
| Twisted table | 1.9 MiB |
| Ring buffer, 64 x ~720 bytes | 46 KiB |
| Walk buffers, frame scratch, pairing and curve state | < 256 KiB |
| Code, constants, stacks and output | < 2 MiB |

The total is < 2^29 bytes. Advice is zero; no trace is retained.

## Total and score

    S = floor(31*(isqrt(floor(r/12))+1)/20) + 2^115
      = 102501770541727644560508313235544969372 steps (escapes included)
    W = S*(1 + 12288/20000) + 2^42/20000
      < 165478858362565109378484620887683700880 units
    log2(W) <= 126.959915   (rounded up at the sixth decimal; experiment.py recomputes it)

For comparison: without Frobenius classes, the plain distinguished-point walk
(1.5*sqrt(r) steps at 1 + 8192/20000 units each) gives 128.508968.
