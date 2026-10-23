# bls12-381-fixed-g1-invert-g2: hard-capped STNFS discrete log in GT

Status: ready for independent review under the exploratory policy. The attack is
the standard Special-Tower Number Field Sieve discrete logarithm in
`F_{p^12}`, applied to the GT side of the fixed-argument pairing inversion. It
is not claimed to be novel or to have been executed at production size.

This revision replaces every "expected work" statement of the earlier drafts by
an **explicit, enforced work cap**: the algorithm below is truncated at a fixed
instruction budget and returns FAILURE if the budget is exhausted. Nothing in
this document claims an expected value as an upper bound. The submitted score is
the truncated algorithm's budget divided by the reference price, and the
algorithm's success probability is at least `1/2` over uniform inputs and fresh
coins, by Markov's inequality applied to its total work.

## 1. Claim

Algorithm `A`, run once per challenge instance on fresh coins:

1. Let `Q0` be the public G2 generator. Compute
   `g = e(P,Q0) = pairing(Q0,P,final_exponentiate=True)`. The pinned pairing is
   nondegenerate and `P`, `Q0` both have prime order `r`, so `g != 1` and `g`
   generates the prime-order group GT.
2. Compute `a = log_g z` by STNFS-DL in `F_{p^12}`, using the published
   Guillevic-Singh BLS12-381 polynomial and Table 10 parameters. The phase is
   run under the hard instruction cap `R_cap` of Section 4.5; it is a sequence
   of complete attempts, each of which is itself capped, and it stops (returning
   FAILURE) as soon as the cumulative instruction count reaches `R_cap`.
3. On success return `Q = [a]Q0` and verify `e(P,Q) = z` and `[r]Q = O`.

For every valid instance (`P` a nonidentity point of G1, `z` in GT, including
the identity target `z = 1`) the algorithm is correct whenever it returns: GT
is cyclic of prime order `r`, so `log_g z` is unique modulo `r`, and
`e(P,[a]Q0) = e(P,Q0)^a = g^a = z`. The identity target gives `a = 0` and
`Q = O`, verified by `e(P,O) = 1`.

The complete per-instance work of `A` is bounded by the hard cap

```text
R_cap = 2^135.146602479037 scalar instruction equivalents
W_cap = R_cap / C + H = 2^120.858890099488 + 2^14
      = 2^120.858890099488   target-group-operation equivalents
score = log2(W_cap) = 120.858891
```

where `C = 20000` is the fixed dense-Fp12 multiplication price of the cost
model and `H = 2^14` reference target operations are charged explicitly for the
final pairings, the scalar multiplication `[a]Q0` and verification. The exact
rational claim is `120858891/1000000`.

Success probability is at least `1/2` (Section 6), over both the uniform
nonidentity `P` and the uniform target `z` and fresh algorithm coins, which is
above the required `0.39`. Peak storage is at most `2^72.389801` bytes.
Preprocessing is the field-only relation collection and linear algebra; the full
per-instance charged total (`2^120.858891` target-operation equivalents) is
reported as its bound, with no reuse credit. Nonuniform advice is 0. No
auxiliary powers are used.

## 2. Reduction and correctness

Let `r` be the 255-bit prime subgroup order from the protected target profile
and `e: G1 x G2 -> GT` the pinned pairing in the library's `pairing(Q,P)`
argument order. GT is cyclic of prime order `r`. For nonidentity `P` in G1 the
map `Q -> e(P,Q)` is an isomorphism `G2 -> GT`; hence for every `z` in GT there
is a unique `a` mod `r` with `g^a = z` where `g = e(P,Q0)`, and `Q = [a]Q0` is
the unique preimage of `z` under that isomorphism.

Correctness therefore follows from computing the discrete logarithm
`a = log_g z` in the multiplicative group of `F_{p^12}` (which contains GT as
the order-`r` subgroup) and reducing it modulo `r`. The algorithm verifies its
own output with one pairing check, so a wrong `a` cannot be returned as success.
No assumption is made that `P` is a fixed generator across instances: `P` enters
only through `g = e(P,Q0)`, which is recomputed per instance and charged.

The GT target is a random element of the order-`r` subgroup. STNFS computes
logarithms in `F_{p^12}^*`; to obtain `a mod r` the standard reduction is used
(sparse linear algebra over `Z/(r)` on the relation matrix, with Schirokauer
maps for the virtual logarithms). The individual logarithm of the target is
obtained by the standard special-q descent after multiplying the target by a
small random power of `g` to obtain a descent-friendly representative; that
randomization is charged in Section 4.2.

## 3. The STNFS attack and its published cost

### 3.1 Parameters

The attack uses exactly the published Guillevic-Singh (GS21) simulation for the
BLS12-381 curve in `F_{p^12}`, Table 10 (`BLS-381` row): degree-12 special-tower
NFS with the BLS12-381 polynomial, sieving bound `A = 686`, `eta = 6`,
smoothness bound `log2 B = 65.316`, factor base `log2(#B) = 60.848`, simulated
relations `log2(#rel) = 61.3898`, matrix weight `wt = 200`, filtering factor
`flt = 20`, and block-Wiedemann with 4 words per entry.

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
126-bit figure is tabulated via GS21 in subsequent work. A 2023 survey applies
the Barbulescu-Duquesne (BD19) SexTNFS model to BLS12-381 and reports the lower
figure `2^122.54`, but that figure is not reproducible from the survey's own
stated parameters: its sieving space `(2*639+1)^12/2 = 2^122.85` already exceeds
the reported total, which contradicts BD19's convention that one sieving-domain
element costs one operation, and its factor-base count `2^65.20` does not match
`2B/(A log B) = 2^63.32` for its stated `A = 2`. GS21's STNFS simulation is
self-consistent and is the basis of the IETF 126-bit statement, so it is used
here; a valid lower estimate could only reduce the submitted bound.

### 3.3 The GT discrete logarithm is the only hard step

Everything outside STNFS is cheap and is charged explicitly or bounded:
`g = e(P,Q0)` and the final verification are single pairings (`2^14` reference
operations each, `H = 2^14` charged, see below), `[a]Q0` is a single G2 scalar
multiplication (`< 2^13` reference operations, inside `H`), serialization is
`O(2^8)` bytes, and the descent randomization costs `< 2^20` GT multiplications
(`2^-100` of the total).

## 4. Instruction-equivalent ledger and normalization

The cost model defines `cost = H + R/C`, with `C = 20000` scalar AArch64
instructions per dense Fp12 multiplication and `R` the work in scalar 64-bit
instruction equivalents. Each of the four factors of `R_cap` below is stated
with its source; none of them is an expected value.

### 4.1 Conversion from the published model unit: 128 instructions per model operation

BD19 (whose cost model GS21 uses) states: "We use the same convention as in
[Len01a, Section 2.4.6] and count a clock cycle as one operation." GS21's
sieving term is therefore a cycle-equivalent count, and one GS21 "model
operation" (the evaluation of one sieving-domain element) is one cycle in the
model.

Two published, measured quantities convert cycles to scalar instructions:

1. **Model-to-hardware calibration.** BD19's Table 6 reports measured/model
   ratios `log2(c_sieve)` from `-2` to `3.5` on the record NFS computations they
   calibrate against, i.e. at most `2^3.5 = 11.31` measured cycles per model
   operation on the worst record (the GGMT17 dimension-three STNFS computation,
   about 27 cycles per tuple against GS21's model price
   `ln ln(2^65.316) = 3.8129`). We charge **`16 = 2^4` cycle-equivalents per
   model operation**, which is `1.41x` the worst measured ratio ever reported
   for this model family, so the conversion does not extrapolate beyond
   published measurements; the `1.41x` surplus is the budget for the phases
   GS21's total omits (Section 4.2).
2. **Instructions retired per cycle.** The widest current AArch64 cores retire
   at most 8 scalar instructions per cycle (Apple M-series and Neoverse V2 are
   8-wide; narrower cores retire fewer, which cannot increase the instruction
   count of the algorithm). Vector or SIMD work is counted as the corresponding
   number of scalar instruction equivalents, so vectorization cannot raise the
   scalar count past this bound either.

Hence the single-attempt instruction budget per model operation is

```text
16 cycles/op * 8 instructions/cycle = 128 instructions per GS21 model operation
B_single = 128 * 2^126.146602479037 = 2^133.146602479037 instructions.
```

### 4.2 What the 1.41x surplus covers (all subdominant, all bounded)

- **Individual-logarithm descent.** BP14 analyse the individual logarithm for
  the medium/high-characteristic NFS at exponent `L_Q[1/3, 1.44]` against
  `L_Q[1/3, 1.923]` for relation collection. At `Q = p^12 = 2^4572`
  (`(ln Q)^{1/3}(ln ln Q)^{2/3} = 59.0`) that is `2^85` against `2^113`
  model operations, i.e. about 28 bits below the sieving term; GS21 states the
  same ("the cost of individual discrete logarithm phase is very small in
  comparison to the other two steps"). Budgeted: `<= 2^-28` of the core.
- **ECM cofactorization.** `L_B(1/2, sqrt 2) * 2B/ln B ~= 2^87.6` operations,
  about 38 bits below the core.
- **Filtering, relation storage and I/O, the field-representation isomorphism,
  descent randomization, output.** Each is below `2^70` instructions, i.e. below
  `2^-63` of the core.
- **Rank deficiency and descent failure at fixed volume.** Handled by the retry
  factor of Section 4.3, not by this surplus.

Every term above is at least 28 bits below the sieving term, so the 1.41x
surplus (0.5 bits) is not merely adequate: it is about `2^27` times larger than
the largest omitted phase's actual estimated cost.

### 4.3 Trials and failures: retry factor 2

The published simulation targets a relation set that overshoots the factor base
by `2^61.3898 / 2^60.848 = 1.45x`, which is what makes full rank the normal
case. We nevertheless charge the retry multiplicity pessimistically as
`1/p = 2`, i.e. we assume a single capped attempt succeeds with probability
only `1/2`. (Under the smoothness model the relation count is concentrated
`2^29` standard deviations above the requirement, and rank deficiency is far
rarer, so `p >= 1/2` is very loose.) Each attempt is truncated at `B_single`
instructions, so

```text
E[T]  <= 2 * B_single = 2^134.146602479037 instructions,
```

where `T` is the total instruction work of `A` until it succeeds. This is the
only place an expectation appears, and it appears only on the left of an
inequality that is then converted to a hard cap by Markov.

### 4.4 Hard truncation: Markov factor 2

`A` halts and returns FAILURE when its cumulative work reaches
`R_cap = 2 * E[T]`. By Markov's inequality applied to the nonnegative `T`,

```text
Pr[T > R_cap] <= E[T] / R_cap = 1/2,
```

so `A` returns a correct answer within its cap on at least half of the
(in-instance `P`, target `z`, fresh coins) triples. The cap is a property of the
algorithm, not an estimate: the implementation counts retired instruction
equivalents and stops.

### 4.5 Final normalized score

```text
R_cap = 2 * 2 * 128 * 2^126.146602479037 = 2^135.146602479037 instructions
W_cap = R_cap / 20000 + 2^14
      = 2^135.146602479037 / 2^14.287712379549 + 2^14
      = 2^120.858890099488 + 2^14
      = 2^120.858890099488   (the `H` term is below `2^-106` of the total)
score = log2(W_cap) = 120.858891   (rounded upward to six decimals)
```

Every instruction is counted exactly once: `2^126.146602479037` is the published
model cost, `16 x 8 = 128` converts model cycles to scalar instructions, `2x`
covers retries, `2x` is the Markov truncation, and division by `20000` converts
instructions to the fixed target unit. The exact rational claim is
`120858891/1000000`.

## 5. Why the conversion is defensible

- **No measured quantity is extrapolated.** The only calibration data used are
  the published measured/model ratios of BD19 Table 6 (`log2 c_sieve <= 3.5`)
  and the maximum AArch64 scalar retirement width (8). The charge is `2^4 * 8`,
  a `1.41x` margin over the worst measured record, and the omitted phases are
  `28+` bits below the core.
- **No "at most" from an expected value.** The only expectation, `E[T]`, is used
  in the direction Markov's inequality requires (upper bound on a tail
  probability), and the reported work is a truncation that the algorithm
  enforces, not an expectation.
- **No amortization.** The full field-only computation is charged per instance
  with no credit for reuse; the target-specific descent is charged inside the
  same budget.
- **Consequences of error.** If the reviewer deems the retry or truncation
  factors unnecessary, the *expectation* for one capped attempt is
  `2^118.858891`; that smaller figure is not claimed here, precisely because it
  is only an expectation. The submitted claim is the capped figure.
- **Direction of the comparison to rho.** The submitted number is not obtained
  by any assumption about the difficulty of rho; it is the model cost of index
  calculus converted through the track's own instruction accounting, which is
  the accounting the cost model prescribes for sieving ("Sieving, smoothness
  testing, linear algebra and descent count through their computation; they do
  not need to use field multiplication").

## 6. Success probability and edge cases

- The algorithm is truncated, so its success event is "returns a verified `Q`
  before the cap". By Section 4.4 that event has probability at least `1/2` over
  uniform nonidentity `P`, uniform `z` in GT and fresh coins, which exceeds the
  required `0.39`. This is a lower bound on success, not a confidence score.
- `P` is nonidentity by the challenge definition, so `g != 1`. `z = 1` gives
  `a = 0` and `Q = O`, accepted and verified (`e(P,O) = 1`); this is also the
  cheapest instance.
- Every other `z` has a unique logarithm modulo the prime `r`; the relation
  matrix is solved modulo `r`, so no non-`r` factor of `p^12-1` can corrupt the
  answer.
- The returned point is verified by a pairing; a failed check is reported as
  failure, never as a wrong answer. GT has prime order, so `[a]Q0` is in G2 and
  has the required order.
- The challenge requires correctness for the exact supplied `P` and `z`; the
  algorithm touches `P` only through `g = e(P,Q0)`, recomputed and charged per
  instance.

## 7. Memory, preprocessing and advice

- **Memory.** The block-Wiedemann matrix has
  `2^60.848/20 = 2^56.528` rows of weight 200 with 40 bytes per nonzero:
  `2^69.17` bytes; the relation file is `2^61.3898 * 256` entries of 8 bytes,
  `2^72.39` bytes. Total peak `<= 2^72.389801` bytes, dominated by the relation
  file. Memory has no scalar contribution in this track; for scale, the
  organizer's archived generic control for the sibling G1 orientation reports
  `2^137.643857` bytes.
- **Preprocessing.** The field-only relation collection and linear algebra are
  the attack's dominant cost and are charged in full, per instance, with no
  reuse credit; the reported preprocessing bound is the same
  `2^120.858891` target-operation equivalents, included in the total.
- **Advice.** 0. The factor base is constructed by the attack and charged
  through the relation-collection and linear-algebra models; no nonuniform
  advice is supplied.

## 8. Heuristics and limitations

Heuristics, all disclosed, and all in the direction of the published model:

1. Smoothness and relation counts use Dickman's rho with GS21's alpha
   correction, as estimated by the GS21 simulator. This is the standard
   heuristic basis of every NFS cost estimate; it is not re-proved here.
2. The model-to-hardware calibration (`16` cycles per model operation) is taken
   from BD19 Table 6's measured range, whose worst record is used; the
   conversion is a factor `1.41` above that worst measurement.
3. Subdominant phases (descent, ECM, filtering) are bounded by their published
   asymptotic estimates (BP14, GS21), which are `28+` bits below the core; the
   `1.41x` surplus absorbs them with a `2^27` margin.
4. The retry factor `2` is a deliberately loose bound; the relation
   overshoot of the published parameters makes the true per-attempt failure
   probability negligible.
5. The result is an analytic reproduction of a published cost model at
   production parameters, not an executed 2^120 attack, as the rules allow
   ("An analytic argument need not run an exponentially large attack").

## 9. Verification evidence

`experiment.py` (declared in `submission.json`, runnable with
`scripts/run_experiments.sh`) recomputes, from the published Table 10
parameters: the sieving space, the sieving and linear-algebra costs, their sum,
the instruction conversion (`16 x 8`), the retry and truncation factors, the
normalized target-operation count including the `H = 2^14` verification charge,
and the resulting score. It asserts that the recomputed sieving space matches
the published `2^124.08`, that the total is within `0.25` bits of the published
`2^126`, and that the score matches the claim to `1e-6`. It is arithmetic
evidence for the ledger, not an attack run.

`claim.json` carries the same numbers (score `120.858891`, success `0.5`,
memory `2^72.389801`, preprocessing `2^120.858891`, advice 0), and
`submission.json`'s first claim states the same bound in words.

## 10. References

- A. Guillevic and S. Singh. "The number field sieve in the medium
  characteristic setting." Mathematical Cryptology 1(1), 2021 (GS21). Table 10
  (`BLS-381`), Section 6.1 model, Eq. 6.5, and the statement that the
  individual-logarithm phase is very small.
- R. Barbulescu and S. Duquesne. "Updating key size estimations for pairings."
  J. Cryptology 32, 2019 (BD19). Clock-cycle convention (Section 2.4.6 of
  Len01a), Table 6 record calibration (`log2(c_sieve)` between `-2` and `3.5`).
- X. Wang, P. Zheng and Q. Xing. "Security Analysis of Pairing-based
  Cryptography." arXiv:2309.04693, 2023. Section 5.2 applies the BD19 model to
  BLS12-381 and reports `2^122.54`; the figure and its inconsistent
  factor-base arithmetic are discussed in Section 3.2 above.
- IETF draft-irtf-cfrg-pairing-friendly-curves-11, Section 4.2.1: BLS12-381
  security "expected to be 126 rather than 128 bits".
- L. Gremy, A. Guillevic, F. Morain and E. Thome. "Computing discrete
  logarithms in GF(p^6)." SAC 2017 (GGMT17), the dimension-three STNFS record
  with `c_sieve ~= 27`, cited in BD19 and used for the calibration cross-check.
- R. Barbulescu and C. Pierrot. "The multiple number field sieve for medium-
  and high-characteristic finite fields." LMS J. Comput. Math., 2014 (BP14),
  individual-logarithm analysis (`L[1/3, 1.44]`).
- Curvesma.sh `tracks/docs/COST_MODEL.md` and
  `tracks/classical-group-work-v2.json`: `H + R/C`, `C = 20000`, instruction
  basis, and the requirement to charge sieving, linear algebra and descent
  through their computation.
