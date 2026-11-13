# bn254-gt-dlp: gated segmented distinguished-point Pollard rho

**What this is.** A generic (square-root) attack on the BN254 target-group discrete
logarithm, priced end to end with an analytic AArch64 instruction bound: one dense
Fp12 multiplication per walk step, a bounded distinguished-point table, a
fresh-coin run gate, every one-time cost charged to this target. Submitted score:
`log2 E[W] = 128.296605831` reference units, success at least 0.39 over uniform
challenges and fresh coins.

**What this is not.** Not the best known attack and not near this track's research
frontier: published STNFS submissions on the same benchmark reach about `2^88.4`.
Nothing here reproduces, checks or competes with that line; no optimality is
asserted. It is a generic baseline whose walk, arithmetic, checks and accounting
are executable and reproducible.

## 0. What changed from the two rejected revisions

Two earlier candidates on this track were rejected, and this revision removes the
specific defects rather than arguing around them.

1. `5d74a156` charged a borrowed BLS12-381 reference price for the Fp12
   multiplication. It is gone: the price below is a cell-level AArch64 bound for
   this arithmetic.
2. `d24ab0e7` had three defects, all fixed here:
   - the 2,000,000-step C walk started from an arbitrary Fp2 element instead of
     `g^a z^b`, so it never entered GT. `walk_stats` now starts at `g^a z^b` for
     fixed public labels, and the partition and distinguished-point windows are
     read from the **canonical** form of the same coefficient (one Montgomery
     reduction per step, charged below) instead of Montgomery residues.
   - the tower-addition count was an estimate. It is now the count taken from the
     shipped kernel: **287** Fp add/sub/neg operations (90 in the 18 Fp2 products,
     156 in the three Fp6 layers, 41 at the Fp12 layer including multiplication by
     xi = 9+u).
   - the segment loop tested its budget only between segments and could overrun.
     It is now clamped to the remaining budget, and the toy campaign under the
     clamp reports `max_u = 2.0` rather than 2.1991.

## 1. Target, tower, encoding

Protected profile (`target-profile.json`, copied here): Ethereum alt_bn128, GT ⊂
Fp12 of prime order r, fixed generator g = e(P0,Q0); input z = g^x with x uniform
in [0,r); output x; identity answered as 0. No auxiliary powers, no source-group
preimage, no n parameter.

From the profile's `Fp[w]/(w^12-18w^6+82)` with `u = w^6-9`:

    Fp2 = Fp[u]/(u^2+1),  Fp6 = Fp2[v]/(v^3-(9+u)),  Fp12 = Fp6[w]/(w^2-v)
    a_m = (c_m + 9 c_{m+6}) + c_{m+6} u,  m = 0..5,  m = j + 2k

`evidence.py` converts both ways, checks the Karatsuba tower product against an
independent schoolbook product in the w-basis, and the C kernel re-checks
`g^r = 1` and `g != 1` for the profile generator.

## 2. Algorithm

Public increments `c_j = SHA256("bn254gt-c"||j) mod r`, `d_j` likewise,
`M_j = g^{c_j} z^{d_j}`. State (Y, a, b) with Y = g^a z^b. One step:

    t   <- canonical value of Y's first Fp coefficient        (one Montgomery reduction)
    j   <- (t >> 120) & 31                                    (five middle bits)
    Y   <- Y · M_j                                            (one dense Fp12 multiplication)
    a   <- a + c_j mod r ;  b <- b + d_j mod r
    DP  iff (low 128 bits of the new canonical t) < D

Segments: budget `T = ceil(1.6·floor(sqrt(r)))` =
236,714,811,010,862,895,392,714,499,903,555,895,296 steps, at most 12,700
segments (mean 2^113.843877). A segment starts at a fresh uniform state
`(a, b)` — built as `g^a z^b`, whose exponentiation pair is charged — and stops at
the first distinguished point, `D = floor(12,700·2^128/T) = 18,256`
(probability 2^-113.843917). Each DP is stored as (128-bit SHA-256 fingerprint, a,
b) in 80 bytes; on a fingerprint match `x = (a−a')/(b'−b) mod r`, verified by
`g^x == z` over all twelve coefficients. Zero denominator (probability exactly
1/r) or a failed check is a failed trial. Gate: a fresh 64-bit coin opens the run
with probability `p_open = 0.39/0.66 = 0.5909091`; a closed run costs one unit.

## 3. Cost

| item | count | each | bound |
|---|---:|---:|---:|
| Fp multiplication (32 products of MUL+UMULH, accumulator and carry cells, conditional subtract/select, operand and result traffic) | 54 | 640 | 34,560 |
| Fp add / sub / neg (carry fix-up or conditional add of p), counted in the shipped kernel | 287 | 56 | 16,072 |
| loads, stores, spills, call/return, loop control of the multiplication | — | — | 8,192 |
| **one dense Fp12 multiplication** | | | **58,824** |
| one Montgomery reduction to canonical form (windows are read canonically) | 1 | 736 | 736 |
| partition, two four-limb label updates mod r, DP test, loop | — | — | 200 |
| **per step** | | | **59,760 = 2.9880 units** |

Totals:

    R_run  = T · 59,760 + 2^40                     = 2^143.343 instructions
    W_run  = R_run / 20,000                        = 2^129.055598 units
    E[W]   = p_open · W_run + (1 - p_open) · 1     = 2^128.296605831 units
    memory = 12,700·80 + 32·448 + 4,096            = 1,034,432 bytes < 2^20
    success= p_open · 0.66                         = 0.39

One-time term, 2^40, charged once for this target: 32 multipliers (64
exponentiations, 2^30); **12,700 fresh-segment exponentiations pairs** (25,400
exponentiations of ≤ 508 multiplications ≈ 2^39.4); DP bookkeeping (≤ 12,700
fingerprint/table operations, 2^32); recovery, inversion mod r, `g^x == z`,
output (2^22); randomness and seed expansion (2^20). Nothing is amortised over
future targets; no advice is assumed. Public profile constants and the algorithm
design are the only excluded items. The unit is the protected one, R/20,000,
with the BN254 arithmetic charged through its own instruction count.

The shipped portable C kernel measures **52,313** instructions per multiplication
(callgrind `Ir`, x86-64, gcc 14.2.0 `-O3 -fno-tree-vectorize -fno-builtin`), i.e.
below the analytic bound, as an upper bound should be; it is a consistency check,
never the price, and no cross-architecture transfer enters the charge.

## 4. Correctness

**Invariant.** Y_0 = g^{a_0} z^{b_0}; each step multiplies by `g^{c_j} z^{d_j}` and
adds `(c_j, d_j)`, so Y = g^a z^b throughout. `evidence.py` re-derives
`g^{a+bx}` on the real GT every 250 steps of 1,500 (0 bad steps); `kernel.c`
mode `s` runs 2,000,000 production steps from `g^a z^b` in the counted
implementation.

**Recovery.** Equal states give `a − a' = x(b' − b) mod r`; r is prime so a
nonzero difference has a unique inverse. `evidence.py` builds four genuine
colliding pairs on the real GT with known logarithms and runs the production
recovery: all four recover x and verify `g^x = z`.

**Degeneracy.** Conditional on the true x and the visited group elements, any
unconditioned `d_j` stays uniform, so a zero denominator is a nonzero combination
of independent uniform scalars: probability exactly 1/r.

## 5. Success and its support

Random-map value at u = 1.6: `1 − exp(−1.28) = 0.721963`. Charged **0.66**, below
both the model and both measurements; the gate turns it into an overall 0.39.

Toy campaign with the production walk class (same SHA-256 increments, same
middle-bit partition, same trailing-window DP, same segment bookkeeping, same
recovery), on prime-order subgroups of F_q^* with q = 2n+1 prime, budget
2·sqrt(n) **clamped hard at the budget**, 256 slots, 100 trials per size:

| n | successes | success by u = 1.2 / 1.4 / 1.6 / 1.8 / 2.0 |
|---|---|---|
| ≈2^27 | 82/100 | 0.46 / 0.61 / 0.69 / 0.79 / 0.82 |
| ≈2^31 | 83/100 | 0.45 / 0.56 / 0.71 / 0.80 / 0.83 |
| model | — | 0.513 / 0.625 / 0.722 / 0.802 / 0.865 |

Production-encoding statistics on the real GT (`gt-stats.json`, kernel mode `s`,
2,000,000 steps from `g^a z^b`, DP threshold scaled to 2^-10): partition
chi-square **18.51** on 31 degrees of freedom (95% point 50) and **1,995 DP hits
against 1,953 expected**. The 1,500-step Python run reports chi-square 27.5 and
1 hit against 1.46.

Disclosed heuristics: **H1** the 32-adding walk with this partition and these
increments behaves as a random map on the 2^253.6 group, so collisions appear with
the birthday constant; toys and encoding statistics support this at their size,
they do not prove it at production size. **H2** segment starts and increments are
independent across fresh coins. **H3** the gate coin is fresh and independent.
Only the success claim depends on H1; the work bound does not.

## 6. Position, limitations, provenance

- `2^128.30` here; `2^135.216266` for this benchmark's promoted rho baseline;
  about `2^88.4` for the published STNFS line, which this entry neither reproduces
  nor assesses.
- No amortisation, no advice, no hidden oracle: the partition and DP test are bit
  extractions from a computed, canonically reduced coefficient; hashing, table
  access, recovery and verification are charged.
- Provenance: Pollard (1978); Teske, *Math. Comp.* 70 (2001) 809-825;
  van Oorschot and Wiener (1999); Brent (1980) for the cycle-detection alternative
  not used. Field constants and generator are the protected profile's. No prior
  submission's code, data or measured constant is used; no dependencies declared.

## 7. Reproduction

```
echo '{}'                                        | python3 -I -B evidence.py
echo '{"part":"toys","trials":100,"bits":[26,30]}' | python3 -I -B evidence.py
echo '{"part":"rebuild"}'                         | python3 -I -B evidence.py
gcc -O3 -fno-tree-vectorize -fno-builtin -o kernel kernel.c
./kernel g                 # g^r = 1, g != 1
./kernel s 2000000 10      # production walk statistics on the real GT
valgrind --tool=callgrind ./kernel 4096 ; valgrind --tool=callgrind ./kernel 8192
```
