# bls12-381-fixed-g1-invert-g2: STNFS discrete logarithm in GT

Status: ready for independent review under the exploratory policy. The attack is
the standard Special-Tower Number Field Sieve discrete logarithm in
`F_{p^12}`, applied to the GT side of the fixed-argument pairing inversion. It
is not claimed to be novel or to have been executed at production size. The
bound is analytic, reproduces the published Guillevic-Singh STNFS cost for
BLS12-381 from its published parameters, and then applies a deliberately
conservative instruction-equivalent and overhead conversion.

## 1. Claim

For every valid challenge instance (nonidentity `P` in G1, `z` in GT, including
the identity target `z=1`), the following algorithm returns `Q` in G2 with
`e(P,Q)=z` in the pinned py_ecc v8.0.0 convention `pairing(Q,P)`:

1. Fix a public nonidentity `Q0` in G2 (for example the standard G2
   generator). Compute `g = e(P,Q0) = pairing(Q0,P,final_exponentiate=True)`.
   Since the pinned pairing is nondegenerate and both inputs have prime order
   `r`, `g != 1`, so `g` generates the prime-order group GT.
2. Solve `a = log_g z` by STNFS-DL in `F_{p^12}`, modulo `r`.
3. Return `Q = [a]Q0` and verify `e(P,Q) = z` and `[r]Q = O`.

The discrete logarithm exists and is unique modulo `r` for every `z` in GT, so
the Las Vegas algorithm succeeds with probability 1 after retries. The complete
per-instance work is at most

```text
R  <= 2^133.146602 instruction equivalents
W  =  R / C = 2^133.146602 / 20000 = 2^118.858891
```

target-group-operation equivalents, where `C=20000` is the fixed dense-Fp12
multiplication price of the cost model. The submitted score is
`log2(W) = 118.858891`. Peak storage is at most `2^72.389801` bytes;
preprocessing is at most `2^118.858891` target-operation equivalents and is
included in the total; nonuniform advice is 0. No auxiliary powers are used.

The number of target-group operations `H` is zero for this attack: STNFS uses
small-field arithmetic, sieving and sparse linear algebra, not dense Fp12
multiplications. The cost model charges this work through `R/C`, exactly as
directed by "Sieving, smoothness testing, linear algebra and descent count
through their computation; they do not need to use field multiplication."
The final pairing check and the `[a]Q0` scalar multiplication do use reference
operations, but fewer than `2^16` in total; both are already included in the
128x conversion margin and are below `2^-100` of the total, so stating `H=0`
changes nothing at the claimed six decimals.

## 2. Reduction and correctness

Let `r` be the 255-bit prime subgroup order from the protected target profile,
and let `e: G1 x G2 -> GT` be the pinned pairing. GT is cyclic of prime order
`r`. For uniform nonidentity `P` in G1 and any nonidentity `Q0` in G2, the
pairing value `g = e(P,Q0)` is nonidentity by nondegeneracy, hence is a
generator of GT. Every target `z` in GT therefore has a unique
`a = log_g z` in `[0,r)`, including `z=1` (where `a=0`). Bilinearity gives

```text
e(P, [a]Q0) = e(P,Q0)^a = g^a = z .
```

So returning `Q=[a]Q0` is correct for every valid input, and `z=1` correctly
returns the identity of G2. The output is checked with the exact pinned pairing
and a subgroup check. `P` is used only through `g`; it is never replaced, and
no auxiliary powers or structured advice are assumed.

The STNFS implementation returns logarithms relative to an internal generator
`g0` of the order-`r` subgroup. One descent is performed for `z` and one for
the base element `g` (or the implementation is normalized directly to base
`g`); since `r` is prime and `g` is a generator, `log_{g0} g != 0 mod r` and

```text
a = log_g z = (log_{g0} z) * (log_{g0} g)^{-1}  mod r .
```

Both descents, the final modular inverse, the base conversion and the
verification are covered by the conversion margin of Section 4. This is the
standard way index calculus answers a subgroup DLP; it does not require
computing logarithms of all of GT.

## 3. The STNFS attack and its published cost

### 3.1 Parameters

The field is `F_{p^12}` for the BLS12-381 prime `p` (381 bits) and subgroup
order `r` (255 bits), with `p^12` of 4569 bits. The published STNFS polynomial
selection for this curve is Guillevic-Singh (GS21, Table 11):

```text
h(Y)   = Y^6 - Y^2 + 1
f(Y,X) = X^12 - 2Y X^10 + 2Y^3 X^6 + Y^5 X^2 + Y^2 - 1
g(Y,X) = X^2 + 15132376222941642752 Y
```

with `deg h = 6`, `deg_X f = 12`, `deg_X g = 2`, `k=12`. The factor base is
the set of prime ideals of norm below `B` in the two number fields, and
relations are collected from tuples

```text
(a_0,...,a_5, b_0,...,b_5) in [-A,A]^12,  lc(b)>0,
```

whose two norms are B-smooth. GS21's simulator chose `A = 686` and
`B = 2^65.316` (window `2^65.316`, average norms `2^724.89` and `2^497.04`,
B-smooth probability `2^-62.5660`). The resulting data are (GS21, Table 10,
curve BLS-381):

```text
sieving space       2^124.08
factor base         2^60.8480
relations obtained  2^61.3898
total cost          2^126
```

The simulator also reports, for these parameters, `wt = 200` nonzeros per
matrix row, filtering factor `flt = 20`, and linear algebra over the 255-bit
prime `r`, i.e. `ceil(255/64) = 4` 64-bit words per modular multiplication.

### 3.2 Reproduction of the published arithmetic

GS21 Section 6.1 models relation collection as the sieving volume times
`log(log B)` (natural logarithms) and block-Wiedemann linear algebra as
`cnst * wt * (#B/flt)^2`. Recomputing from the published rounded parameters:

```text
space      = (2*686+1)^12 / 2                       = 2^124.077390921765
sieve cost = 2^124.077390921765 * ln(ln(2^65.316)) = 2^126.008213062200
LA cost    = 4 * 200 * (2^60.848 / 20)^2            = 2^122.696000000000
total      = 2^126.008213062200 + 2^122.696        = 2^126.146602479037
```

This reproduces GS21's reported total `2^126` (their table rounds to integer
bits); we use the larger recomputed value `2^126.146602479037`, not the rounded
`2^126`, so the ledger is conservative. The same 126-bit figure appears in the
IETF pairing-friendly-curves draft, which states that the security of
BLS12-381 "is expected to be 126 rather than 128 bits" (citing GMT19), and the
126-bit figure is tabulated via GS21 in subsequent work. A 2023
survey applies the Barbulescu-Duquesne (BD19) SexTNFS model to BLS12-381 and
reports the *lower* figure `2^122.54`, but that figure is not reproducible
from the survey's own stated parameters: its sieving space
`(2*639+1)^12/2 = 2^122.85` already exceeds the reported total, which
contradicts BD19's convention that one sieving-domain element costs one
operation, and its factor-base count `2^65.20` does not match
`2B/(A log B) = 2^63.32` for its stated `A=2`. GS21's STNFS simulation is
self-consistent and is the basis of the IETF 126-bit statement, so it is used
here; the larger published cost is the conservative choice, and a valid lower
estimate could only reduce the submitted bound.

Relations are converted to linear equations modulo `r` with the standard
Schirokauer maps; the sparse system is filtered by `flt=20` and solved with
block Wiedemann. The individual logarithm of the target is obtained by the
standard special-q descent. GS21 states explicitly that "the cost of individual
discrete logarithm phase is very small in comparison to the other two steps",
and BP14 gives the asymptotic analysis; Section 4 nonetheless charges the worst
measured sieving constant plus a further allowance 1.25 times that sieving work
for the descent and the other omitted phases.

## 4. Instruction-equivalent ledger and normalization

The cost model defines `cost = H + R/C`, with `C=20000` scalar AArch64
instructions per dense Fp12 multiplication, and `R` the work in scalar 64-bit
instruction equivalents. Here `H=0`.

### 4.1 Conversion from the published model unit

BD19, whose cost model GS21 uses, states: "We use the same convention as in
[Len01a, Section 2.4.6] and count a clock cycle as one operation." GS21's
sieving term is therefore a cycle-equivalent count. BD19 calibrate this model
on records: their Table 6 reports `log2(c_sieve)` between -2 and 3.5, i.e.
measured/model ratios between 0.25 and 11.3, and the GGMT17 dimension-three
STNFS record is the worst case at about 27 cycles per tuple. GS21's model
charges `ln ln B = 3.8129` per tuple, so the GGMT17 calibration is about
`27 / 3.8129 = 7.08` measured cycles per GS21 model operation.

The widest current AArch64 cores retire at most 8 scalar instructions per
cycle (Apple M-series, Cortex-X and Neoverse V2 are 8-wide); narrower cores
retire fewer, which cannot increase the instruction count of the algorithm.
Vector work is counted as multiple scalar equivalents, so vectorization cannot
raise the scalar count beyond this bound either. The worst measured sieving
calibration therefore needs at most

```text
57 <= 8 * 7.08 instructions per GS21 model operation.
```

We charge `128` instructions per model operation: 57 for the worst measured
sieving constant at the maximum retirement width, plus 71 (1.25x the sieving
requirement) for the descent, ECM, retries, filtering, verification, I/O and
output that GS21's total omits. Hence

```text
R <= 128 * 2^126.146602479037 = 2^133.146602479037 instructions.
```

### 4.2 Unmodeled phases and success amplification

GS21's total contains relation collection and linear algebra and explicitly
excludes the individual logarithm (calling it very small) and the ECM
cofactorization (calling it negligible). The cost model requires charging all
preprocessing, advice construction, random sampling, trials and failures,
sorting, lookup, transfers, verification, output and amplification.

The 128-instruction charge per model operation consists of the 57-instruction
worst measured sieving calibration and a remaining 71 instructions per model
operation (1.25x the sieving requirement) for those unmodeled phases:

- individual-logarithm descent and base conversion: GS21 states that this
  phase "is very small in comparison to the other two steps", and BP14
  analyses it as subdominant; an allowance of 1.25 times the entire worst-case
  sieving work is far above "very small";
- ECM cofactorization: recomputes to
  `L_B(1/2,sqrt(2)) * 2B/ln B ~= 2^87.6` operations, 38 bits below the core;
- failed relations, rank deficiencies, descent failures and retries: the
  attack is Las Vegas and repeats with fresh coins; the simulated relation
  count exceeds the factor base by 1.45x and full-rank systems are the normal
  case, so the expected retry work is a small fraction of the 1.25x allowance;
- filtering, relation storage and I/O, the field-representation isomorphism,
  the scalar multiplication `[a]Q0`, the final pairing verification, and
  serialization: all are below 2^70 instructions individually.

The total instruction bound is

```text
R <= 128 * 2^126.146602479037 = 2^133.146602479037.
```

This factor of 128 is not a claim that every listed phase actually costs that
much; it is an upper-bound allowance calibrated on the worst measured sieving
constant at the maximum retirement width, with the remainder reserved for the
phases GS21 omits. Section 5 records the margin explicitly.

### 4.3 Final normalized score

```text
W = 2^133.146602479037 / 20000
  = 2^133.146602479037 / 2^14.287712379549...
  = 2^118.858890099487...
score = log2(W) = 118.858891   (rounded upward to six decimals)
```

The exact rational claim is `118858891/1000000`. Preprocessing is the
field-only relation collection and linear algebra, which dominates the attack
and which the attacker could in principle reuse across instances but does
**not** get credit for here; because the descent is very small and no reuse
credit is taken, we report the full per-instance charged total
`2^118.858891` as the preprocessing upper bound. It is included in the total.
Every instruction is counted exactly once: the `8x` is the maximum instruction
retirement width, the `16x` covers the worst measured sieving calibration plus
the omitted-phase allowance, and the division by 20000 converts instructions to
the fixed target unit.

## 5. Why the 128x conversion is adequate

- **Published model spread.** BD19's less conservative BLS12-381 estimate is
  `2^122.54`; GS21's is `2^126.15`. The difference is 3.6 bits; the lower
  figure is not reproducible from its own stated parameters (Section 3.2), so
  the larger self-consistent total is used.
- **Measured sieving constant.** BD19's Table 6 reports `log2(c_sieve)` from
  -2 to 3.5, i.e. measured/model ratios from 0.25 to 11.3, and the GGMT17
  dimension-three STNFS record is the worst case at about 27 cycles per tuple.
  GS21's model charges `ln ln B = 3.8129` per tuple, so GGMT17 corresponds to
  about `27 / 3.8129 = 7.08` measured cycles per GS21 model operation. At the
  8-wide retirement limit this is at most 57 instructions per operation. The
  128-instruction charge is 2.25x that worst measured requirement and leaves
  71 instructions per operation, 1.25x the sieving requirement, for the phases
  GS21 omits.
- **Descent and ECM.** GS21 states that the descent "is very small in
  comparison to the other two steps" and calls the ECM cost negligible; the
  ECM term recomputes to `2^87.6`, 38 bits below the core. An allowance of 1.25
  times the entire worst-case sieving work is far above GS21's description.
- **Retries.** The relation count `2^61.39` exceeds the factor base `2^60.85`
  by 1.45x, which is the simulated oversieving target; rank and descent
  failures are handled by fresh-coins repetition, and under the simulated
  oversieving they are rare relative to the 1.25x allowance.

The conversion is therefore an upper bound even under the worst measured
sieving constant and the maximum retirement width, and it still budgets
1.25 times that worst-case sieving for the phases that GS21's total omits. A
discrepancy large enough to invalidate it would require the dimension-twelve
STNFS implementation to exceed the worst published record constant by more
than 2.25x at full retirement width, or to have a descent far larger than
GS21's explicit "very small" statement; both are far outside the published
record data and the published analysis.

## 6. Success probability and edge cases

- `P` is nonidentity by the challenge definition. If `P` were the identity the
  pairing input would be invalid; the implementation rejects it.
- `z` is in GT by the challenge definition; `z=1` yields `a=0` and `Q=O`, and
  is verified by `e(P,O)=1`. Every other `z` has a unique discrete logarithm
  modulo the prime `r`.
- STNFS relation collection and descent are Las Vegas procedures: if a run
  fails to produce a full-rank system or a smooth enough descent, the attack
  restarts with fresh randomness inside the charged expected work. Success
  probability is therefore 1, above the required 0.39.
- The challenge requires correctness for the exact supplied `P` and `z`; the
  algorithm uses `P` through `g=e(P,Q0)` and its output is checked against
  `e(P,Q)=z`. No assumption is made that `P` is a fixed generator.

## 7. Memory, preprocessing and advice

- Filtered matrix: `2^60.848/20 = 2^56.526072` rows, 200 nonzeros per row,
  40 bytes per nonzero (8-byte index plus four 64-bit value words), giving
  `2^69.491856` bytes. Block-Wiedemann vectors add about `2^61.5` bytes.
- Relation file: at most 256 entries per relation (a generous bound for two
  degree-6 sides at these norm sizes), `2^61.3898 * 256 * 8 = 2^72.3898`
  bytes. The sieve itself is streamed in lines/chunks and is smaller.
- Reported peak storage is the larger figure, `2^72.389801` bytes. Memory has
  no scalar contribution in this track; for scale, the organizer's archived
  generic control for the sibling G1 orientation
  (`tracks/docs/archive/BASELINE_QUALIFICATION_V1.md`) reports
  `2^137.643857` bytes, about 65 bits more.
- Preprocessing is the field-only relation collection and linear algebra; it
  dominates the attack, and the full per-instance charged total
  (`2^118.858891` target-operation equivalents) is reported as its upper bound.
  It is charged in full with no reuse credit. The polynomials and simulator
  constants are public profile constants from GS21 Table 11; algorithm design
  effort and public constants are excluded by the rules.
- Nonuniform advice is 0. The factor base is constructed by the attack and is
  charged through the relation-collection and linear-algebra models; it is not
  supplied.

## 8. Heuristics and limitations

Heuristics, all disclosed:

1. Smoothness and relation counts use Dickman's rho with the GS21 alpha
   correction, estimated by the GS21 simulator over published samples. This is
   heuristic, as in all NFS analyses.
2. The sieving cost per tuple is modeled as `ln ln B` cycles, calibrated on
   record computations; the worst measured STNFS value (GGMT17) is about 7.08
   cycles per model operation. The 128-instruction charge covers this worst
   value at the 8-wide retirement limit and still leaves 1.25x the sieving
   requirement for the omitted phases.
3. The matrix weight (200), filtering factor (20) and block-Wiedemann cost are
   GS21's modeling choices. The 1.25x omitted-phase allowance covers moderate
   deviations.
4. The descent is assumed subdominant, as stated by GS21 and analyzed by BP14;
   the omitted-phase allowance is 1.25 times the entire worst-case sieving
   work, far above a "very small" descent.
5. Floating-point reproduction of the published costs is done in binary
   floating point; all quoted bounds are rounded upward at the sixth decimal.

Limitations:

- No production-size run is claimed or possible here. The evidence is the
  analytic model, its reproduction from published parameters, and the general
  correctness reduction.
- The result is an upper bound in the published cost model, not a theorem that
  no better attack exists.
- Peak storage is astronomically large (exabytes), as in the organizer's
  generic control; memory is unscored.
- The field-only preprocessing is chargeable to this instance under the
  no-amortization rule; no benefit from hypothetical future instances is taken.
- The 8x retirement-width and 16x calibration/omitted-phase factors are
  conservative bookkeeping conventions, not measurements of a specific binary.

## 9. Verification evidence

- The repository's 15 scientific checks pass (`python3 -m unittest discover -s
  tests`), including the fixed-G1 pairing semantics and the pinned pairing
  equations used in Section 2.
- `experiment.py` recomputes the GS21 BLS-381 cost from the published Table 10
  parameters and reproduces the submitted ledger:
  `space 2^124.077390921765`, `sieve 2^126.008213062200`,
  `LA 2^122.696`, `total 2^126.146602479037`, `score 118.858890099487`.
  It asserts the reproduction is within 0.25 bits of the published total and
  that the score matches the claim. It is arithmetic evidence, not an attack
  run.
- The initial input is the Table 10 parameter set; the script uses only the
  Python standard library and runs in well under a second.

## 10. References

- A. Guillevic and S. Singh. "On the alpha value of polynomials in the Tower
  Number Field Sieve algorithm." Mathematical Cryptology 1(1), 2021, pp. 1-39.
  Tables 10-11 (BLS-381 parameters, total cost `2^126`), Section 6.1
  (relation-collection and line algebra models, descent negligible).
  https://journals.flvc.org/mathcryptology/article/view/125142 and
  https://eprint.iacr.org/2019/885
- R. Barbulescu and S. Duquesne. "Updating key size estimations for pairings."
  J. Cryptology 32, 2019; https://eprint.iacr.org/2017/334. Equation (2), the
  clock-cycle convention, Table 6 record calibration (`log2(c_sieve)` between
  -2 and 3.5), and Section 4.2 of the model.
- X. Wang, P. Zheng and Q. Xing. "Security Analysis of Pairing-based
  Cryptography." arXiv:2309.04693, 2023. Section 5.2 applies the BD19 model to
  BLS12-381 and reports `2^122.54`; the figure and its inconsistent
  factor-base arithmetic are discussed in Section 3.2 above.
- IETF draft-irtf-cfrg-pairing-friendly-curves-11, Section 4.2.1: BLS12-381
  security "expected to be 126 rather than 128 bits".
- L. Gremy, A. Guillevic, F. Morain and E. Thome. "Computing discrete
  logarithms in GF(p^6)." SAC 2017 (GGMT17), the dimension-three STNFS record
  with `c_sieve ~= 27` (about 7.08 cycles per GS21 model operation) cited in
  BD19.
- R. Barbulescu and C. Pierrot. "The multiple number field sieve for medium-
  and high-characteristic finite fields." LMS J. Comput. Math., 2014 (BP14),
  individual-logarithm analysis.
- Curvesma.sh `tracks/docs/COST_MODEL.md` and
  `tracks/classical-group-work-v2.json`: `H + R/C`, `C=20000`, instruction
  basis, and the requirement to charge sieving, linear algebra and descent
  through their computation.
