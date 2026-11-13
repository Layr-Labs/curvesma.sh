# bn254-gt-dlp: gated segmented distinguished-point Pollard rho

**What this is.** A generic (square-root) attack on the BN254 target-group
discrete logarithm, priced end to end with an analytic AArch64 instruction bound:
one dense Fp12 multiplication per walk step, a bounded distinguished-point table,
a fresh-coin run gate, and every one-time cost charged to this target. The
submitted score is `log2 E[W] = 128.169282977` reference units, success at least
0.39 over uniform challenges and fresh coins.

**What this is not.** Not the best known attack on this target and not near this
track's research frontier: published STNFS submissions on the same benchmark reach
about `2^88.4`. Nothing here reproduces, checks or competes with that line, and no
optimality is asserted. It is a generic baseline whose walk, arithmetic, tamper
checks and cost accounting are executable and reproducible.

---

## 1. Target, tower, encoding

Protected profile only (`target-profile.json`, copied into this folder): Ethereum
alt_bn128, GT ⊂ Fp12 of prime order r, fixed generator g = e(P0,Q0). Input z = g^x
with x uniform in [0,r); output x; the identity target is answered as x = 0. No
auxiliary powers, no source-group preimage, no n parameter.

The profile's `Fp[w]/(w^12-18w^6+82)` with `u = w^6-9` gives the tower

    Fp2  = Fp[u]/(u^2+1)
    Fp6  = Fp2[v]/(v^3 - xi),   xi = 9+u
    Fp12 = Fp6[w]/(w^2 - v)

with the w-basis to tower map

    a_m = (c_m + 9 c_{m+6}) + c_{m+6} u,    m = 0..5,   m = j + 2k,

j the Fp6 half and k the v^k slot. The coupling through `9 c_{m+6}` is not
cosmetic: `evidence.py` converts both ways and checks the Karatsuba tower product
against an independent schoolbook product in the w-basis on deterministic test
elements, and the shipped C kernel re-checks `g^r = 1` and `g != 1` in Fp12 for the
profile generator.

## 2. Algorithm

**Fixed public increments.** For j = 0..31,

    c_j = int(SHA256("bn254gt-c" || j)) mod r,     d_j likewise for "bn254gt-d".

No advice, no secret input; these are public constants of the algorithm. (An
earlier revision used `int("...c"||j) * K mod r`, which is an arithmetic
progression in j; that is a defect and was replaced by the hash above.)

**Step.** State (Y, a, b) with Y = g^a z^b. One step:

    j  <- five bits from the middle of Y's first Fp coefficient, (coef >> 120) & 31
    Y  <- Y · M_j            one dense Fp12 multiplication, M_j = g^{c_j} z^{d_j}
    a  <- a + c_j mod r ;   b <- b + d_j mod r

**Segments and distinguished points.** W = 12,700 table slots, total budget
T = ceil(1.6·floor(sqrt(r))) = 236,714,811,010,862,895,392,714,499,903,555,895,296
steps, mean segment length L = T/W = 2^113.843877. A segment starts at a fresh
uniform state (a, b uniform mod r) and runs to the first distinguished point:

    DP(Y)  iff  (low 128 bits of the imaginary part of Y's c1[2]) < D,
    D = floor(W · 2^128 / T) = 18,256,    P(DP) = 2^-113.843917 per step.

The pair (128-bit SHA-256 fingerprint of the DP, a, b) is stored in 80 bytes. A
match means two segments ended at the same group element; recovery is

    x = (a - a')/(b' - b) mod r,

verified by g^x == z over all twelve coefficients before output. A zero
denominator (probability exactly 1/r) or a failed verification is a failed trial.
The identity input is answered as x = 0.

**Gate.** A fresh 64-bit coin opens the run with probability
p_open = 0.39/0.66 = 0.5909091; a closed run costs one unit and performs no group
arithmetic.

## 3. Cost: analytic AArch64 instruction bound

One dense Fp12 multiplication in this tower is 54 Fp multiplications (18 Fp2
Karatsuba products of 3 Fp multiplications each) plus ~210 Fp add/sub/neg
operations (Fp2/Fp6/Fp12 Karatsuba combinations and the two multiplications by
xi = 9+u, which are additions only). Cell-level AArch64 bounds:

| item | count | each | bound |
|---|---:|---:|---:|
| Fp multiplication (32 products: MUL+UMULH, accumulator and carry cells, four-word conditional subtract/select, operand and result traffic) | 54 | 640 | 34,560 |
| Fp add / sub / neg (four limbs with carry fix-up or the conditional add of p) | 210 | 56 | 11,760 |
| loads, stores, spills, call/return, loop control of the whole multiplication | — | — | 8,192 |
| **one dense Fp12 multiplication** | | | **54,512** |
| partition, two four-limb label updates mod r, DP test, loop | — | — | 200 |
| **per step** | | | **54,712 = 2.7356 units** |

The 640 per Fp multiplication is deliberately loose (32 multiplies + 32 carried
adds + 24 accumulator cells + a 64-instruction select + 96 memory instructions is
350 on a cell-by-cell count); it is charged loose so that no optimism is hidden in
the cells. The shipped portable C kernel executes the same tower in **52,313**
instructions per multiplication (callgrind `Ir`, x86-64, gcc 14.2.0 `-O3
-fno-tree-vectorize -fno-builtin`), i.e. below the analytic bound, which is the
expected direction for an upper bound; no cross-architecture transfer is used for
the charge. Earlier revisions of this entry charged a *borrowed BLS12-381
reference price* (22,996) instead; that was not an AArch64 bound for this
arithmetic and it has been withdrawn.

One-time work, charged once for this target:

| item | bound |
|---|---:|
| 32 multipliers M_j = g^{c_j} z^{d_j}: 64 exponentiations of ~508 multiplications each | 2^30 |
| DP handling: ≤ 12,700 fingerprint/serialise/table operations | 2^32 |
| recovery (one inversion mod r, one Fp12 exponentiation for g^x == z), output | 2^22 |
| **total** | **2^36** |

Totals:

    R_run  = T · 54,712 + 2^36                    = 2^143.216 instructions
    W_run  = R_run / 20,000                       = 2^128.928275 units
    E[W]   = p_open · W_run + (1 - p_open) · 1    = 2^128.169282977 units
    memory = 12,700·80 + 32·448 + 4,096 = 1,034,432 bytes < 2^20
    success = p_open · 0.66 = 0.39

The unit is the protected one: R/20,000 for R scalar AArch64 instruction
equivalents, with the BN254 arithmetic charged through its own instruction count
rather than given a unit price. Nothing is amortised over other targets; the
public profile constants, the algorithm design and the constant table derivation
are the only excluded items.

## 4. Correctness

**Invariant.** Y_0 = g^{a_0} z^{b_0}; each step multiplies Y by g^{c_j} z^{d_j} and
adds (c_j, d_j) to the labels, so Y = g^a z^b at every step. `evidence.py` runs the
production walk on the real GT and re-derives g^{a+bx} independently every 250
steps of 1,500 (0 bad steps); `kernel.c` mode `s` does the same over 2,000,000
production steps in the counted implementation (histogram and DP rate recorded in
`gt-stats.json`), and mode `g` re-checks g^r = 1.

**Recovery.** Equal states give a - a' = x(b' - b) mod r; r is prime, so a nonzero
b' - b has a unique inverse. `evidence.py` builds four *genuine* colliding pairs on
the real GT with known logarithms (choosing b' ≠ b and solving for a') and runs the
production recovery on them: all four recover x and verify g^x = z.

**Degeneracy.** For a fixed true x, any d_j not conditioned by the visited group
elements remains uniform, so a zero denominator is a nonzero linear combination of
independent uniform scalars and has probability exactly 1/r.

## 5. Success probability and its support

Random-map model: colliding segment pairs are approximately Poisson with mean
u²/2 = 1.28 at u = 1.6, so P(success at T) ≈ 1 - exp(-1.28) = 0.721963. We charge
**0.66**, below both the model and the measurements, and the gate turns it into an
overall 0.39.

**Toy campaign with the production walk.** `evidence.py part=toys` runs *the same
code* — same SHA-256 increment derivation, same middle-bit partition placement,
same trailing-window DP rule, same segment bookkeeping, same recovery — on
prime-order subgroups of F_q^* with q = 2n+1 prime (the subgroup of squares, dense
enough that both bit rules have their nominal probabilities). 100 trials per size,
budget 2·sqrt(n), 256 slots:

| n | q | successes | success by u = 1.2 / 1.4 / 1.6 / 1.8 / 2.0 |
|---|---|---|---|
| ≈2^27 | 28 bits | 82/100 | 0.46 / 0.61 / 0.69 / 0.79 / 0.82 |
| ≈2^31 | 32 bits | 83/100 | 0.45 / 0.56 / 0.71 / 0.80 / 0.83 |
| model | — | — | 0.513 / 0.625 / 0.722 / 0.802 / 0.865 |

The measured rate at u = 1.6 is 0.69 and 0.71 against the model's 0.722; charging
0.66 leaves a margin below both. These runs are evidence at their actual size:
they do not prove collision behaviour at 2^253.6 and are not a production DLP.

**Production-encoding statistics.** On the real GT, `gt-stats.json` records
2,000,000 production steps at a DP threshold scaled to 2^-10: partition histogram
chi-square 19.45 on 31 degrees of freedom (95% point 50), and 1,898 DP hits against
1,953 expected. The 1,500-step Python run reports chi-square 27.5 and 1 hit against
1.46. Both the partition and the DP rule therefore behave as specified on the
actual encoding, at production arithmetic, over a large number of steps.

**Disclosed heuristics.** H1: the 32-adding walk with this partition and these
increments behaves like a random map on the 2^253.6-element group, so a collision
appears with the birthday constant above; the toy campaign and the encoding
statistics support this but do not prove it at production size. H2: segment starts
and increments are independent across fresh coins. H3: the gate coin is fresh and
independent. Only the *success* claim depends on H1; the work bound T does not.

## 6. Position, limitations, provenance

- Generic square-root attack: `2^128.17` here, `2^135.216266` for this
  benchmark's promoted rho baseline, and about `2^88.4` for the published STNFS
  line, which this entry neither reproduces nor assesses.
- No amortisation, no advice, no hidden oracle: the partition and DP test are bit
  extractions from computed elements, and all hashing, table access, recovery and
  verification are charged.
- No optimality claim: finite-field index calculus and NFS remain eligible and
  stronger on this target.
- Provenance: Pollard (1978) for rho; Teske, *Math. Comp.* 70 (2001) 809-825 for
  the r-adding walk; van Oorschot and Wiener (1999) for distinguished points;
  Brent (1980) for the cycle-detection alternative not used here. Field constants
  and generator are the protected profile's. No prior submission's code, data or
  measured constant is used, so no dependencies are declared.

## 7. Reproduction

```
echo '{}'                                        | python3 -I -B evidence.py
echo '{"part":"toys","trials":100,"bits":[26,30]}' | python3 -I -B evidence.py
echo '{"part":"ledger"}'                          | python3 -I -B evidence.py
echo '{"part":"rebuild"}'                         | python3 -I -B evidence.py
gcc -O3 -fno-tree-vectorize -fno-builtin -o kernel kernel.c
./kernel g                 # g^r = 1 and g != 1 in Fp12
./kernel s 2000000 10      # production walk statistics on the real GT
valgrind --tool=callgrind ./kernel 4096 ; valgrind --tool=callgrind ./kernel 8192
```
