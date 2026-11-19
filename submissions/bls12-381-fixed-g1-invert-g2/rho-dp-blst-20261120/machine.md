# Finite scalar implementation and work bounds

These are **analytical scalar AArch64 instruction-equivalent bounds**, not
timings, on one processor with 64-bit scalar registers. Python is an
executable mathematical reference, not the priced implementation.

## The reference operation (one unit per walk step)

The protected cost model (tracks/docs/COST_MODEL.md) defines the unit for this
track as the pinned blst dense Fp12 multiplication `blst_fp12_mul`
(commit f7a88067...). It runs in blst's tower representation
Fp2/Fp6/Fp12, with reduction, internal additions and memory work included,
measured at 22,996 instructions and rounded to C = 20,000. The rules state
that "a reference operation costs one group-equivalent unit" and that its
implementation must not be double-counted.

Each walk step makes exactly one call `blst_fp12_mul(out, Y, M_j)` on two
canonical GT elements in that representation. That call is charged 1 unit;
everything else is charged by instruction count / C. We assume nothing about
blst's output beyond what is canonicalized below: blst keeps Fp residues as
six 64-bit words, i.e. values < 2^384. We re-establish full canonical form
ourselves.

The tower is Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-(1+u)),
Fp12 = Fp6[w]/(w^2-v). Its w equals the protected basis' w, since w^6 = 1+u
gives u = w^6-1. `rho.py` gives the exact basis map and checks that the tower
product agrees with the protected product (experiment section 2).

## Per-step non-reference work: E <= 8,192 instructions

| Work | Bound |
|---|---:|
| Canonicalize 12 coefficients: at most 9 conditional subtractions of p each (2^384 < 10p), each a 6-word subtract-with-borrow plus masked select, <= 48 instructions | 12*10*48 = 5,760 |
| Index: load word 0 of t_0, AND 31, multiply by the 576-byte record stride and add the base | 8 |
| DP test: load two words of t_1, mask the low 110 bits, OR, compare and branch | 8 |
| Labels: a += c_j, b += d_j mod r (four-limb add_mod from arithmetic.c, <= 256 each, including loads of c_j and d_j from a separate 64-byte record) | 512 |
| Step counter: three-word increment and compare against S, branch | 16 |
| Swap Y/out pointers, call/return overhead for blst_fp12_mul, spills, loop control | 256 |
| Total | 6,560 < **E = 8,192** |

The walk state is a 576-byte GT element plus 64 bytes of labels, kept in two
buffers alternated by pointer swap; no per-step copy is needed. The table has
32 records of 576 + 64 bytes. Canonicalization makes the coefficient words a
function of the group element, so the index, DP flag, fingerprint and
equality are well defined.

## Distinguished points: charged inside F

DP handling runs only when the 110-bit test passes. Per DP:
- Load the 128-bit fingerprint.
- Probe an open-addressing table of 2^21 slots of 96 bytes (fingerprint,
  a, b, occupancy), stopping at an empty slot or a match, with probes capped
  at 2^10 (a full probe run counts as failure).
- Store, or on a match, compute and test x.

The cost is at most 2^10 probes * 32 instructions + 2^10 = 2^16 instructions
per DP. At most 2^21 DPs can be stored (overflow is a failure), and at most
2^21 + 64 DP events are handled before any terminal outcome (each event
either stores or is one of at most 64 verifications), so the total is
< 2^38 instructions. Each of at most 64 match verifications costs:
- one modular inverse mod r (< 2^22),
- a subtraction and a multiplication mod r,
- g^x via at most 510 Fp12 products, each computed either with blst_fp12_mul
  (1 unit = 20,000 instructions-equivalent) or with our own bounded 2^21
  routine.

The total is < 64 * 2^30 = 2^36. Probe-length failure: with load at most
1/8 and a random fingerprint hash, a probe run of 2^10 has probability
< 2^-1000.

## Fixed work F = 2^40 instructions (all one-target, charged in full)

The bounds below are those of the six-limb Montgomery arithmetic (8,192 per Fp
multiplication, 512 per add/sub) and of the schoolbook protected-basis Fp12
product G <= 2^21 (144*8192 + 400*512 + 2^16). The arithmetic and its
derivations are in `arithmetic.c`:

* **Fp multiplication cost.** 36+36 product cells at <= 48 instructions, 54
  carry cells at <= 16, 12 q-passes at <= 32, plus 1,024 for init,
  subtraction and calls: 5,728 < 8,192. REDC output is < 2p < 2^384, so one
  conditional subtraction makes it canonical.
* **Input.** Bounded canonical decoding of P and z, curve and z != 0 checks,
  and the z^r = 1 check (<= 510 products): < 2^31.
* **Pairing g = e(P,Q0)**, following the protected py_ecc algorithm with
  dense products:
  - twisting Q0 and casting P to Fp12: <= 64 product-equivalents;
  - the Miller loop over the 64-bit protected ate loop count
    d201000000010000: <= 64 iterations, each with one tangent line, the two
    numerator/denominator square-and-multiply updates, one projective
    doubling and, on the 5 set bits after the leading bit, one chord line,
    two multiplications and one addition. Each line, doubling or addition is
    < 32 Fp12 operations, so < 200 per iteration and < 12,800 products;
  - the quotient f_num/f_den by Fermat inversion f_den^(p^12-2): <= 9,144
    products;
  - the final exponent (p^12-1)/r (< 4,320 bits): <= 8,640 products.

  That is < 2^15 products, < 2^36 instructions. Only the function's value is
  used, and it equals the protected value on this exact P because it is the
  same deterministic computation.
* **Basis change and Montgomery conversion.** The basis change for g and z is
  24 additions per element. Conversion is 12 Fp multiplications per element.
  Both are < 2^20.
* **Randomness.** 66 scalars: 255-bit candidates with at most 1024 draws
  each, at 64 instructions per fresh bit: < 2^31.
* **Records.** 32 records and the start: two exponentiations (<= 510
  products each) and one product per record, <= 33*1021 products, < 2^36.
* **DP handling and verifications** (above): < 2^38 + 2^36.
* **Output.** Inverse mod r; affine double-and-add for [x]Q0 with optional
  [r]Q = infinity check, < 2^35. Fp2 multiplication is 4 Fp
  multiplications plus adds, < 2^16. An Fp2 inverse uses an Fp inverse by
  p-2 (<= 762 multiplications), < 2^24. A point operation is < 2^25, and
  1,020 point operations are < 2^35. The optional [r]P check on G1 is
  < 2^34.
* **Fixed overhead.** Constants, parsing, encoding, allocation, code
  loading and shutdown: < 2^32.

The sum is < 2^39, so F = 2^40 leaves slack.

## Peak storage

- DP table: 2^21 * 96 bytes = 192 MiB.
- Walk records: 32 * 640 bytes = 20 KiB.
- Walk buffers and scratch: < 64 KiB.
- Pairing and curve state: < 64 KiB.
- Code, constants (including the final-exponent bit string), stacks and
  output: < 1 MiB.

The total is < 2^28 bytes. Advice is zero; no trace or random tape is
retained.

## Total and score

    S = floor(3*(isqrt(r)+1)/2) + 2^115 = 343524753603842645718308609411224961024
    W = S * (1 + 8192/20000) + 2^40/20000
      < 484232492679976593404527815826117680641 units
    log2(W) <= 128.508968   (rounded up at the sixth decimal)

For comparison, the same walk with our own fully priced schoolbook product
(<= 2^21 instructions instead of 1 unit) would be about 6.4 bits higher. The
reduction comes only from using the protected reference operation as the
protected rules allow, not from an unpriced oracle.
