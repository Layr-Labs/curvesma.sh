# bls12-381-g1-ecdlp: Cheon subgroup walk, batched affine additions, window 20 (v6)

Claim: for s uniform in Z_r^* and inputs T = [s]P0, A_i = [s^i]P0 (i = 1..n, n = 2^20), the algorithm below recovers
s with probability at least **0.40**, and it runs under a hard work meter that halts it when its total reaches
**W_max = 2^120.90418** target-group-operation equivalents at n = 2^20 (score 6045209/50000). W_max includes table
construction, the walk, distinguished-point storage and lookup, restarts, the square-root chain and the output, at
worst-case charges. Peak memory is at most 2^34.265982 bytes (20,657,995,776 bytes). Preprocessing is at most
2^24.577492 units. Nonuniform advice is 0 bytes.

This is a classical cold attack on one target. It is not a practical break. The success bound rests on the random-mapping
heuristic H1 (section 4), the fresh toy runs (section 7), and the losses of section 4, which are bounds.

Files: `accounting.py` computes every number here with exact integers; `evidence/accounting.json` is its output.
`evidence/measure-v6exact/` holds the kernel (`cheon_step_w20.c`, with exact call counters), its replica (`replica_w20.py`),
the EC module, the input generator, the summarising script (`summarize_w20b.py`), the exact summary, and the run logs. The
input file (272,315 bytes) is not in the package; `gen_inputs_v4.py <target-profile.json> 7` regenerates it byte for byte.

## 1. Reduction

Let r be the group order, d = 2^20 (which divides r−1), and H = {x^d : x in Z_r^*} the subgroup of d-th powers,
|H| = N = (r−1)/d. Since s is uniform in Z_r^*, h = s^d is uniform in H, and A_d = [h]P0. Given h, recover s by twenty
square roots in Z_r, each sign fixed by the public power A_(2^j), and the final value checked against T = A_1. The chain
costs 41 comb evaluations and twenty Tonelli–Shanks roots in Fr (2-adicity 32, at most 2·32² multiplications each).

## 2. The walk on H modulo mu6, with batched additions and window 20

mu6 = {1, −1, λ, −λ, λ², −λ²} with λ matched to β (so (βx, y) = [λ](x, y) on G1); 6 | N, so mu6 lies in H, and the walk
runs on M = N/6 classes.

**Canonical point and key.** The canonical point has the least Montgomery x-coordinate among the six images, with y of even
parity. Its key is its Montgomery x-coordinate.

**Multipliers and comb.** 64 multipliers g_j in H. From the canonical key k, j = bits 90–95 of k; c2 = c·g_j mod r. The point
[c2]B is a comb of **13 signed 20-bit windows** with carry: window k of magnitude m ∈ [0, 2^19] adds the table entry
m·2^(20k)·B (negated for a negative digit); a zero digit adds nothing. Each base has 13 × 2^19 affine table entries (1.31 GB
for both bases). The point is canonicalised, and c ← μ·c2.

**Batched additions.** Each comb addition is a batched incomplete affine G1 addition (the 6M reference operation): three field
multiplications for the formula and three for Montgomery's simultaneous-inversion trick. In each comb stage k, every record with
a nonzero digit and an accumulator contributes one addition, and the stage shares one inversion. A record's first nonzero digit
is copied, so a point needs at most **12 additions**. A batch of n additions makes exactly **5n−2 `core_mul` and n `core_sqr`**
calls (the reference core).

**Lemma (incomplete additions).** A batched addition needs x1 ≠ x2 and no infinity operand. For windows k ≤ 11, a partial sum s
and an entry e = ±m·2^(20k)·B coincide only if s = ±m·2^(20k) as integers, which forces all lower digits to vanish and then
m = 0, a contradiction. Sums never vanish. For window 12, the congruence s ≡ ±m·2^240 (mod r) has at most about 2^20 bad values
of s, so the heuristic failure probability over at most 2^117 points is at most 2^−100, under the assumption that the low digits
are uniform. A violation aborts the procedure; none occurred in any measured run. Table construction uses stages j ≥ 2, and
j·b ≠ ±b for j ∈ [2, 2^19), so it has no coincidences.

**Records and rounds.** RECORDS = 4096 records advance in rounds, each taking one walk point; even records are tame (base P0)
and odd records wild (base A_d). A trail is capped at 2^100 steps, then restarted from a fresh random start exponent (an
unbatched comb evaluation, which makes no reference-core calls).

**Distinguished points.** A key whose low 90 bits are zero is distinguished (probability 2^−90), stored as a 72-byte record (the
exact high 291 bits of the key, the exponent, the family and padding). Lookup matches the 37-byte prefix exactly.

**Table and meter.** The table has 2^20 buckets of 256 slots, a capacity of 2^27 stored points, and aborts on overflow. The
procedure halts when its worst-case-charged total reaches W_max, so the claimed total bounds the work by construction.

## 3. Points guaranteed

B_min = S(13/20) + 2^112 points, S(τ) = ceil(2·sqrt(τM)). The detection span 2^112 is RECORDS × 2^90 × 2^10. Restarts never
exceed floor(B_min/2^100) = 120,221, and the trails number 124,317.

## 4. Success probability and losses

S is the number of walk points (both families) before a collision is detected. Under H1 (random mapping on M classes),
F(S) = 1 − exp(−S²/(4M)); at S(13/20), F = 0.477954.

- **Round-robin detection.** A merged pair's detection costs up to 4096 points per round. The budget 2^112 covers it except with
  probability 2·e^{−1024}.
- **Cap truncation.** A merge near a trail's 2^100-step cap loses at most 2^−10 per trail, 2^−9 for the pair (charged).
- **Trails without a distinguished point.** Negligible (e^{−1024} per trail over 124,317 trails).
- **In-trail cycles.** At most 1.2 × 10^−5 (charged).
- **Table overflow.** At most 1.8 × 10^−16 (charged; Chernoff bound on bucket loads of mean 128).
- **Window-12 coincidence.** At most 2^−100 (heuristic; section 2).

Total loss 0.0019651. The success probability is at least 0.4227 × (1 − 0.0019651) = 0.42187 ≥ 0.40, with a margin of 0.0219.
The 0.4227 is the Wilson 95% lower bound of 4000 fresh toy trials (section 7).

## 5. Cost

Work is W_max = H + R/C with C = 2000 (the G1 price in COST_MODEL.md). H counts reference units: one per batched addition (12 per
walk point, one per table addition, 12 per initial record). R counts all other scalar instructions.

**Reference core, exact.** The reference core is measured per call type with exact call counts (`core_mul` and `core_sqr` are
noinline wrappers that count their own calls). Each call costs a fixed number of instructions (constant-time field arithmetic):

| Call | Instructions per call |
| --- | ---: |
| `core_mul` | 382 |
| `core_sqr` | 324 |

A batch of n additions therefore costs exactly (5n−2)·382 + n·324 reference instructions. For the walk's stages (n = 4096) that
is 9,149,700, or 2,233.81 per addition. The organizer's convention prices one addition at 1 unit against C = 2000. The rounding is
stated openly: at this kernel's measured core, the convention's unit is 11.7% below the measured instructions per addition. Charging
each addition at its measured core instead would raise the score by about 0.12 bits. We have not changed the charge, since the
convention fixes the unit. The counter increments inside the wrappers are part of the measured core and are not charged to R. The
earlier average of 2,204.82 per addition came from an earlier kernel without the counters and is withdrawn.

**Walk other work per point, exact.** Measured over the forced round (every digit nonzero, so each point performs 12 additions and
one copy, and all twelve stages are occupied). The round made 245,736 `core_mul` and 49,152 `core_sqr` calls. Its exact reference core
is 109,796,400 instructions, and the other work is the rest of the round:

| Quantity | Value |
| --- | ---: |
| additions per point (forced path) | 12 |
| reference core per point (exact) | 26,806.74 instructions |
| other work per point, forced maximum (exact) | **35,415,493 / 4096 = 8,646.36** instructions (4.32 units) |
| other work per point, natural mean (sanity: below the maximum) | 8,610.36 |

The forced maximum is charged for every point. The forced run includes its forcing code, a conservative overcharge. The batch
inversions are covered: a round has at most one inversion per stage, and the forced run has all twelve. The zero-digit skip and
the first-digit copy are covered, since removing a nonzero digit only removes work. The distinguished-point and restart tests are
per-point constant tests; their work is charged as global reserves below.

**Table size.** The forced round was also measured with the full 2^19-entry tables and matches the reduced-table round in every
counted instruction (difference 0).

**Restarts.** Each restart is an unbatched comb evaluation of 12 additions, with no reference-core calls (asserted in the summary,
since the forced-restart round makes the same `core_mul`/`core_sqr` calls as the forced round). Measured other work: 415,422 per
restart. At most 120,221 restarts, plus a fresh start exponent for each (sampling and 20 squarings, bounded).

**Distinguished points.** Each stored point is charged at an analytic bound of 16,000 instructions (hashing, a scan of at most 256
slots, a 72-byte store). The procedure stores at most 2^27 points, so the reserve is 2^27 × 16,000.

**Table build.** Each stage makes 26 batched additions (2 bases × 13 windows), with **128 `core_mul` and 26 `core_sqr`** calls, one
inversion, and the other work measured exactly by stage differencing: 38,444 instructions per stage, plus a fixed part of 3,077,891
(doublings, the first two entries, and the set-up). Over 2^19 − 2 stages the other work is 3,077,891 + (2^19 − 2) × 38,444 =
20,158,728,875 instructions. The earlier figure of 37,733 per stage came from an average and is withdrawn.

**Initial records.** 4096 batched comb evaluations of 12 additions (H = 49,152), other work charged at one point's exact worst case
per record (4096 × 8,646.36 = 35,415,493).

**Totals.**

| Item | value |
| --- | ---: |
| H (reference units) | 2^120.460292 |
| R/C | 2^118.987423 |
| preprocessing (units) | 2^24.577492 |
| **W_max** | **120.90418** (= 6045209/50000) |

**Range of n.** The bound uses A_1, A_(2^j) for j ≤ 20 and A_d, so it holds for every n with 2^20 ≤ n < r at the same W_max.

## 6. Memory

| Item | bytes |
| --- | ---: |
| distinguished-point table: 2^20 buckets × (4 + 256 × 72) | 19,331,547,136 |
| comb tables: 2 × 13 × 2^19 affine points × 96 bytes | 1,308,622,848 |
| records (4096 × 256) | 1,048,576 |
| code, stacks and scratch | 16,777,216 |
| **total** | **20,657,995,776 (2^34.265982)** |

## 7. Evidence

- **Fresh toy runs** (unchanged; `evidence/toy-fresh/`): 4000 trials on 32-bit and 36-bit toy groups with the μ6 action; all
  recovered, no wrong collision. At τ = 13/20: 1752/4000, Wilson 95% [0.4227, 0.4534]. The toy has none of the production losses
  of section 4; they are charged on top.
- **Kernel and replica.** `cheon_step_w20.c` and `replica_w20.py` agree on the final state of all 4096 records after one round for
  the production round (full tables), the restart path (4096 restarts) and the DP path (250 distinguished points). The exact call
  counters do not change the arithmetic.
- **Measurement.** `evidence/measure-v6exact/summary.json` and `run.log`: the `core_mul` and `core_sqr` runs (per-call costs and
  counts), the natural production round, the forced round and forced restart round, the full-table cross-check, and the table build
  at 1024 and 2048 stages.

## 8. Objections and how they are answered

- **Round 5, core subtraction.** The per-batch reference core is no longer an average. Each batch's core is the exact call count
  times the measured per-call cost; the walk, the restart and the table build are each recomputed this way (sections 5 and 7).
- **Round 5, worst-case charge.** The per-point other work is the exact forced maximum (section 5).
- **Round 5, rounding.** The reference unit is the convention's 2000 instructions; the exact measured core at this kernel is
  2,234 per addition for large batches (11.7% above), stated in section 5. The score sensitivity is about 0.12 bits.

## 9. Window and memory

The window is 20, as the user chose: 13 windows, 12 additions per point, 1.31 GB of comb tables. The window-16 package scored
121.196226 with its earlier charges; that comparison used the average core and is not revisited here.

## 10. Limitations

- H1 is an unproved random-mapping heuristic. The toy groups have 32 and 36 bits.
- The reference unit is the convention's 2000 instructions; the exact measured core is 2,234 per addition for large batches.
- The forced run includes its forcing code (a conservative overcharge).
- The restart path uses 16 precomputed exponents in the measurement; the charge covers fresh sampling.
- The window-12 coincidence and the bucket-overflow bound are heuristic or bounded as stated.
- Measurement: one compiler (GCC 14.2.0), one pinned blst, Valgrind 3.24.0, one input set.

## 11. Reproduction

```sh
python3 -B -I accounting.py
echo '{"bits":32,"k":6,"trials":60,"seed":2601}' | python3 -B -I experiment.py
```

The measurement: `gen_inputs_v4.py <target-profile.json> 7` regenerates the input file; the organizer harness runs
`cheon_step_w20.c` with `REF_FUNCS="core_mul"` and `MODES="forced forcedcap1 forcedfull"` (4096 and 8192 points), with
`REF_FUNCS="core_sqr"` for the squaring cost, and with `MODES="build"` at 1024 and 2048 stages; `summarize_w20b.py` computes
the exact summary from the three run directories; `replica_w20.py` checks a full-table dump.
