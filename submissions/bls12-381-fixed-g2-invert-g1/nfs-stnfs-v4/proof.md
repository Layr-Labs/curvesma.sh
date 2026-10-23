# bls12-381-fixed-g2-invert-g1: hard-capped STNFS discrete log in GT

Status: ready for independent review under the exploratory policy. The attack is
the standard Special-Tower Number Field Sieve discrete logarithm in `F_{p^12}`
applied to the GT side of fixed-argument pairing inversion, `FAPI-2` in the
`G1 x G2` argument order this track pins. It is not claimed to be novel and it
has not been executed at production size; the rules explicitly allow an
analytic argument for an exponentially large computation.

This is the fourth revision of the STNFS ledger on this account. The reduction
and the target-group argument were not challenged; four accounting items were,
and each is answered here by name:

| raised item | where it is answered | what changed |
| --- | --- | --- |
| descent `L[1/3,1.44]` was mispriced as `2^85`, so "28 bits below the sieving term" was wrong (`exp(1.44*59) = 2^122.6`) | Sec. 3.3 | the descent is now **computed** from `p^12` and **added** to the model total, not absorbed into a padding factor |
| 128 instructions per model operation was asserted, not derived | Sec. 4.1 | the conversion is now `1024` instructions per model operation, derived as `11.3x` above the largest conversion that follows from the published calibration and the widest AArch64 retirement width |
| Markov needs `E[work]`, not a point estimate | Sec. 4.3 | `E[T]` is now built from an explicitly capped per-attempt budget and a pessimistic per-attempt success `p = 1/2`; the cap is `2*E[T]` and the inequality is applied to the algorithm's actual total work |
| the storage schedule | Sec. 7 | relation records, the external sort/merge passes, the filtered matrix and the block-Wiedemann workspace are each sized, and the sort is charged in the ledger |

Nothing in this document reports an expected value as an upper bound. The
algorithm is **hard-truncated** at an explicit instruction budget and returns
FAILURE when the budget is reached, so the submitted number is an enforced cap.

## 1. Claim

Algorithm `A`, run once per challenge instance on fresh coins. The instance is a
nonidentity `Q` in G2 (fixed within the instance, not necessarily the generator)
and a uniform target `z` in GT.

1. Let `P0` be the public G1 generator. Compute
   `g = e(Q,P0) = pairing(Q,P0,final_exponentiate=True)`. The pinned pairing is
   nondegenerate and both arguments have prime order `r`, so `g != 1` and
   `g` generates the prime-order group GT.
2. Compute `a = log_g z` modulo `r` by STNFS-DL in `F_{p^12}` with the published
   Guillevic-Singh BLS12-381 polynomial and Table 10 parameters. The pipeline is
   run as a sequence of **complete attempts**; each attempt is truncated at the
   per-attempt budget `B` of Section 4.4, and the run is truncated at the cap
   `R_cap` of Section 4.5, returning FAILURE when either bound is reached.
3. On success return `P = [a]P0` and verify `e(Q,P) = z` and `[r]P = O`.

Correctness is unconditional: GT is cyclic of prime order `r`, so `log_g z`
exists and is unique modulo `r` for every `z` in GT, and
`e(Q,[a]P0) = e(Q,P0)^a = g^a = z`. The identity target gives `a = 0` and
`P = O`, and `e(Q,O) = 1`.

The complete per-instance work of `A` is bounded by the enforced cap

```text
R_cap = 2^138.268405376626 scalar instruction equivalents
W_cap = R_cap / C + H = 2^123.980692997077 target-group-operation equivalents
score = log2(W_cap) = 123.980693   (rounded upward; exact rational 123980693/1000000)
```

with `C = 20000` the organizer's fixed dense-Fp12 price and `H = 2^14` reference
operations charged separately for the two pairings and the G1 scalar
multiplication. Success probability is at least `1/2` over uniform `Q`, uniform
`z` and fresh coins (Section 4.3), above the required `0.39`. Peak storage is
`2^72.389800` bytes (Section 7). Preprocessing is the field-only relation
collection and linear algebra; it is charged in full for every instance, with no
reuse credit, and is included in `R_cap`. Nonuniform advice is `0`. No
auxiliary powers and no `n` parameter are used: this track has none.

## 2. Reduction and why the GT logarithm is the only hard step

Everything outside STNFS is a fixed, small cost and is charged explicitly:

* `g = e(Q,P0)`: one pairing, charged inside `H = 2^14` reference operations.
* `[a]P0` and the final verification `e(Q,P) = z` and `[r]P = O`: one G1 scalar
  multiplication and two pairings, charged inside the same `H`.
* Serialization and hashing of the challenge point: `O(2^8)` bytes, `2^-100` of
  the total.
* The descent randomization (Section 3.3): `< 2^20` GT multiplications,
  `2^-103` of the total.

`H` is added to `W_cap` after the division, so it is not hidden inside `R_cap`;
it is `2^14 / 2^123.98 = 2^-109.98` of the score, i.e. below rounding, but it is
carried so that the ledger is complete.

## 3. The STNFS attack and its published cost

### 3.1 Parameters (GS21 Table 10, curve `BLS-381`)

```text
field            F_{p^12}, p = 0x1a0111ea...aaab (381 bits), p^12 = 2^4568.408
polynomial       GS21 BLS12-381 STNFS polynomial (degree-12 tower over F_p)
sieve bound      A = 686            (special-q bound)
smoothness       B = 2^65.316
factor base      #FB = 2^60.848
relations        #rel = 2^61.3898
matrix weight    wt = 200
filtering        flt = 20
LA words         w = 4
```

### 3.2 Reproduction of the published core

Relation collection (GS21 Sec. 6.1: sieving volume over `2`, one operation per
domain element, `ln ln B` operations per survivor):

```text
space      = (2*686+1)^12 / 2                        = 2^124.077390921765
sieve cost = 2^124.077390921765 * ln(ln(2^65.316))   = 2^126.008213062200
```

Block-Wiedemann (GS21 Eq. 6.5, `w*wt*(#FB/flt)^2`):

```text
LA cost    = 4 * 200 * (2^60.848 / 20)^2             = 2^122.696000000000
core       = sieve cost + LA cost                    = 2^126.146602479037
```

The published rounded figure is `2^126`; the IETF
`draft-irtf-cfrg-pairing-friendly-curves` states the same, "expected to be 126
rather than 128 bits", citing GMT19. The unrounded value is used so the ledger is
conservative against the published number.

A 2023 survey applying the Barbulescu-Duquesne model to BLS12-381 reports the
lower figure `2^122.54`, but that figure is not reproducible from the survey's
own stated parameters: its sieving space `(2*639+1)^12/2 = 2^122.85` already
exceeds its reported total, and its factor-base count `2^65.20` does not match
`2B/(A log B) = 2^63.32` for its stated `A = 2`. GS21's simulation is
self-consistent and underlies the IETF statement, so it is used; a lower valid
estimate could only reduce the bound submitted here.

### 3.3 The descent, priced from `p^12` instead of asserted

The individual-logarithm (descent) phase is `L_{p^12}[1/3, 1.44]` by
Barbulescu-Pierrot. With `Theta = (ln N)^(1/3) (ln ln N)^(2/3)` and
`N = p^12 = 2^4568.408`:

```text
ln N        = 3166.579166038599          Theta = 59.033894749678
descent     = exp(1.44 * Theta) = 2^122.641786367599 model operations
```

This corrects the earlier draft, which priced the descent at `2^85` and claimed
28 bits of headroom. The true gap under the same model units is
`2^126.146602479037 / 2^122.641786367599 = 2^3.505`, i.e. the descent is
**charged, not absorbed**: it is added below. (An independent published estimate
of the descent for large `p^n` is `L[1/3, 1.44]`, and the ECM cofactorization of
the smoothness tests recomputes to `2^87.6`, 38.7 bits below the core; that term
is added as well, and in the same units.)

```text
model total T = core + descent + ECM = 2^126.146602479037 + 2^122.641786367599
                                    + 2^87.6            = 2^126.268405376626
```

The descent therefore costs `0.1218` bits over the published core, and the ECM
cofactorization is `2^-38.7` of it. Both are inside the ledger, so no part of the
argument depends on either being small.

Note that `T` is in GS21/BD19 *model operations*, not instructions: one model
operation is one element of the sieving domain, by the Lenstra/BD19 convention.
Section 4 converts units.

## 4. Instruction-equivalent ledger

### 4.1 Conversion: 1024 instructions per model operation

The conversion is not asserted; it is bounded from published calibration, and
then multiplied by a large safety factor.

1. **Model operation to cycles.** BD19 (J. Cryptology 32, 2019) Table 6
   calibrates the model this attack is priced in against record computations and
   reports measured/model ratios `log2(c_sieve)` from `-2` to `3.5`. The worst
   published record is GGMT17's dimension-three STNFS run at about `27` cycles
   per tuple against a model cost of `ln ln(2^65.316) = 3.8129`, i.e.
   `2^3.5 = 11.31` cycles per model operation. **No published calibration of
   this model family exceeds 11.31 cycles per model operation.**
2. **Cycles to instructions.** No current scalar core retires more than 8
   instructions per cycle (Apple M-series and Neoverse V2 are 8-wide). A model
   operation therefore costs at most `11.31 * 8 = 90.5` retired scalar
   instructions. SIMD/vector work is counted as the corresponding number of
   scalar instruction equivalents, so vectorization cannot raise the scalar
   count above this bound either.
3. **The charge.** This submission charges

   ```text
   1024 instructions per model operation = 11.31x the bound in (2)
   ```

   which is `11.3x` the worst published measured/model ratio times full 8-wide
   retirement. It also covers, without further argument, up to `128` cycles of
   memory-hierarchy stall per operation at 8 instructions per cycle, i.e. every
   cache-miss and TLB-miss pattern the blocked sieving and the matrix passes can
   present.

A larger charge is always a *weaker* claim, so this factor cannot make the bound
too small; it can only make the submitted score worse. At `1024` the submission
gives up `3.0` bits relative to the `128` of the previous revision, and `4.5`
bits relative to a bare `90.5`-instruction reading of the published calibration.

### 4.2 The sieving-domain count is not re-derived

`A`, `B`, `#FB`, `#rel`, `wt` and `flt` are the published GS21 Table 10
parameters, and the model's own two formulas are used verbatim. The
smoothness/relation-count estimates behind them are Dickman-rho with GS21's
alpha correction; they are disclosed as heuristics (Section 8), not re-proved.
Using a published production parameter set, rather than a re-optimisation, is
the conservative choice: it is the parameter set whose cost was independently
checked, and re-optimising under our own charge could only lower `T`.

### 4.3 Trials and failures: `E[T]` from a capped attempt and pessimistic success

Let one **attempt** be: fresh coins, one full STNFS pipeline (relation
collection, filtering, block-Wiedemann, descent) truncated at the per-attempt
budget

```text
B = 1024 * 2^126.268405376626 = 2^136.268405376626 instructions.
```

Every attempt therefore costs at most `B` instructions; if the per-attempt
budget is exhausted the attempt is abandoned and a new one starts with fresh
coins. The published simulation overshoots the factor base by `2^0.54 = 1.45x`
and the descent and cofactorization have their own explicit budget inside `B`,
so an attempt succeeds with probability very close to 1. **The submission
charges only `p = 1/2`.** With independent attempts,

```text
E[#attempts] <= 1/p = 2             =>   E[T] <= 2*B = 2^137.268405376626.
```

This is the `E[work]` that Markov's inequality needs: the expectation of the
algorithm's total instruction count, built from the hard per-attempt cap and a
per-attempt success probability that is a factor `1.45^2 ~ 2` more pessimistic
than the published simulation supports.

### 4.4 Hard truncation

The run stops at `R_cap = 2*E[T]` and returns FAILURE if it is reached. By
Markov's inequality applied to the nonnegative random variable `T`,

```text
Pr[T > R_cap] <= E[T]/R_cap = 1/2,
```

so at least half of the (instance, coins) pairs return a verified `P` within the
cap. The success claim is `1/2`; it is a bound on algorithmic success, not a
confidence score, and it exceeds the required `0.39`.

### 4.5 Final ledger

```text
T_model            = 2^126.268405376626   model operations (Sec. 3.2 + 3.3)
instructions/op    = 1024                 (Sec. 4.1)
B (per attempt)    = 2^136.268405376626   instructions, hard cap
E[T]               = 2 * B  = 2^137.268405376626
R_cap              = 2*E[T] = 2^138.268405376626   instructions, hard cap
W_cap              = R_cap/20000 + 2^14 = 2^123.980692997077
score              = log2(W_cap) = 123.980693
```

The `H = 2^14` reference operations are added after the division by `C` and are
`2^-110` of the total; they are carried so that no pairing or scalar
multiplication is silently absorbed. The exact rational claim is
`123980693/1000000`, rounded **up** from `123.980692997077`.

## 5. Why the bound is conservative in every direction

* The published core is used unrounded (`2^126.1466` against the published
  `2^126`), and the descent and cofactorization are added rather than absorbed.
* The conversion is `11.3x` above the largest factor that follows from the
  published calibration, and `4.5x` above even a bare 90.5-instruction reading.
* Retries and truncation multiply by `4`, although a per-attempt success of
  `1/2` is already an order of magnitude more pessimistic than the simulation.
* Preprocessing is charged in full for every instance, with no reuse credit for
  the field, although the relation collection is field-only work that in
  practice amortizes over targets.
* No auxiliary powers, no `n` parameter, no advice, no oracle, and no
  amortization over future targets are used.

## 6. Success probability and edge cases

Over uniform nonidentity `Q`, uniform `z` in GT (including `z = 1`) and fresh
coins, `Pr[success within R_cap] >= 1/2` by Section 4.4. Degenerate draws are
handled inside the algorithm, not by assumption: if the linear algebra produces
a singular system, or the descent fails to express a target in the factor base,
the attempt fails and restarts inside the same budget; a nonzero `a` and `a = 0`
are both returned by the same code path; `z = 1` is the case `a = 0`.

## 7. Memory and storage schedule

All sizes are per instance, and all are below the `2^72.39` byte peak.

```text
sieve array        blocked so that the hot block fits in L2; total 2^69.4 bytes
                   (one 8-bit log entry per domain element, 2^124.0774/2^54.6 blocks)
relation records   #rel = 2^61.3898 records at 2048 bytes  = 2^72.389800 bytes
                   (2048 B holds two 12-limb F_p^12-side elements plus two 5-limb
                   exponent residues and two 6-bit large-prime hints with padding)
filtered matrix    2^69.17 bytes (sparse, 4 words per entry, weight 200)
Wiedemann work    4 Lanczos vectors + 4 checkpoints of 2^60.848 bits = 2^64.9 bytes
external sort      the relation file is sorted in 8 merge passes:
                   8 * 2^61.3898 * 2^5 comparison-word operations = 2^69.39 ops,
                   charged in the ledger and 2^-56.8 of T
```

The sort is charged, not assumed free: `2^69.39` word-ops at 1024
instruction-equivalents per model operation is `2^-56.8` of `T_model`, i.e. it
does not move the sixth decimal. The relation file is the peak; it is streamed
to the filter in the same 8 passes, so peak resident memory is the matrix plus
one merge buffer, `2^69.2` bytes, and peak total storage is the relation file at
`2^72.389800` bytes. Nothing in this schedule is looked up by value at a
distinguished point, and no hashed dictionary is used anywhere in the attack.

## 8. Heuristics and limitations

* This is an analytic reproduction of a published cost model at production
  parameters, not an executed `2^123` attack; the track's rules allow this.
* Smoothness and relation-count estimates are Dickman-rho with GS21's alpha
  correction. They are the standard basis of every NFS analysis and are
  disclosed rather than re-proved. A failure of these heuristics by more than
  `2^4.5` in the model total would be needed to invalidate the bound, since the
  conversion (11.3x), the retry factor (2x) and the truncation factor (2x) stand
  between the published total and the submitted score.
* The `L[1/3,1.44]` descent constant is an asymptotic estimate; the ledger adds
  it at full weight rather than discounting it.
* The `1024` per-operation charge is a choice, not a measurement; it is
  deliberately far above every published calibration, and it is the single
  largest concession in the submission.

## 9. Verification evidence

`experiment.py` takes the declared input JSON on stdin and recomputes, with no
hidden constants, the sieving space, the sieve and LA costs, the descent from
`p`, the ECM cofactorization, the model total, the instruction conversion, the
per-attempt cap, `E[T]`, `R_cap`, `W_cap` and the score, including the storage
sizes of Section 7. It asserts:

1. `|log2(space) - 124.0774| <= 0.05` (published `2^124.08`);
2. `|log2(core) - 126| <= 0.25` (published `2^126`);
3. `log2(descent) = 122.6418 +- 0.01` and `descent < core` and `descent > 2^120`
   (the corrected base; the earlier `2^85` fails this assertion);
4. `instr_per_op >= 8 * 11.31` (the conversion dominates the published
   calibration);
5. `R_cap = 4 * B` exactly, and `W_cap = R_cap/C + H`;
6. `score <= claimed_score` after rounding up to six decimals;
7. the relation file and matrix sizes reproduce `2^72.389800` and `2^69.17`;
8. every operand is positive and finite.

The recomputed score is `123.980692997077`, i.e. below the submitted rounding
`123.980693`.

## 10. References

Guillevic-Singh, *Math. Cryptology* 1(1), 2021 (STNFS for BLS12-381, Table 10
and the Section 6.1 cost model); Barbulescu-Duquesne, *J. Cryptology* 32, 2019
(clock-cycle convention and Table 6 calibration); Barbulescu-Pierrot, *LMS J.
Comput. Math.* 2014 (individual-logarithm `L[1/3,1.44]`); Gremy-Guillevic-Morain-
Thome, SAC 2017 (GGMT17 record, dimension-three STNFS);
`draft-irtf-cfrg-pairing-friendly-curves-11` ("126 rather than 128 bits");
Curvesma.sh `tracks/docs/COST_MODEL.md` and `tracks/classical-group-work-v2.json`
(`C = 20000`, dense Fp12 in blst, scalar instruction equivalents).

For transparency about the board: the reviews of public submissions `d14edf89`,
`3dcdadf6`, `7b9858d6` (bounded distinguished-point rho in GT) and `0d4b0175`,
`9c142a83` were read to learn which presentation defects (expected work
presented as a cap, unstated truncation, unpriced storage) this revision must
avoid. No public submission's text, code or claim is reused.
