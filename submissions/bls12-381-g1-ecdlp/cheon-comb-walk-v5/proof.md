# bls12-381-g1-ecdlp: Cheon subgroup walk, batched affine additions, window 20 (v5)

Claim: for s uniform in Z_r^* and inputs T = [s]P0, A_i = [s^i]P0 (i = 1..n, n = 2^20), the algorithm below recovers
s with probability at least **0.40**, and it runs under a hard work meter that halts it when its total reaches
**W_max = 2^120.904710** target-group-operation equivalents at n = 2^20 (score 12090471/100000). W_max includes table
construction, the walk, distinguished-point storage and lookup, restarts, the square-root chain and the output, at
worst-case charges. Peak memory is at most 2^34.265982 bytes (20,657,995,776 bytes). Preprocessing is at most
2^24.566709 units. Nonuniform advice is 0 bytes.

This is a classical cold attack on one target. It is not a practical break. The success bound rests on the random-mapping
heuristic H1 (section 4), the fresh toy runs (section 7), and the losses of section 4, which are bounds.

Files: `accounting.py` computes every number here with exact integers; `evidence/accounting.json` is its output.
`evidence/measure-v5w20/` holds the kernel (`cheon_step_w20.c`), its Python replica (`replica_w20.py`), the EC module, the
input generator, the summarising script, the measurement summary and log. The input file (272,315 bytes) is not in the
package; `gen_inputs_v4.py <target-profile.json> 7` regenerates it byte for byte.

## 1. Reduction

Let r be the group order, d = 2^20 (which divides r−1), and H = {x^d : x in Z_r^*} the subgroup of d-th powers,
|H| = N = (r−1)/d. Since s is uniform in Z_r^*, h = s^d is uniform in H, and A_d = [h]P0. Given h, recover s by
twenty square roots in Z_r, each sign fixed by the public power A_(2^j), and the final value checked against T = A_1.
The chain costs 41 comb evaluations and twenty Tonelli–Shanks roots in Fr (2-adicity 32, at most 2·32² multiplications).

## 2. The walk on H modulo mu6, with batched additions and window 20

mu6 = {1, −1, λ, −λ, λ², −λ²} with λ matched to β (so (βx, y) = [λ](x, y) on G1); 6 | N, so mu6 lies in H, and the walk
runs on M = N/6 classes.

**Canonical point and key.** The canonical point has the least Montgomery x-coordinate among the six images, with y of
even parity. Its key is its Montgomery x-coordinate.

**Multipliers and comb.** 64 multipliers g_j in H. From the canonical key k, j = bits 90–95 of k; c2 = c·g_j mod r.
The point [c2]B is a comb of **13 signed 20-bit windows** with carry: window k of magnitude m ∈ [0, 2^19] adds the
table entry m·2^(20k)·B (negated for a negative digit); a zero digit adds nothing. Each base has 13 × 2^19 affine table
entries (about 1.31 GB for both bases). The point is canonicalised, and c ← μ·c2.

**Batched additions.** Each comb addition is a batched incomplete affine G1 addition (the 6M reference operation): three
field multiplications for the formula and three for Montgomery's simultaneous-inversion trick, namely the six field
multiplications/squarings of the reference core. In each comb stage k, every record with a nonzero digit and an accumulator
contributes one addition, and the stage shares one inversion. A record's first nonzero digit is copied, so a point needs
at most **12 additions**.

**Lemma (incomplete additions).** A batched addition needs x1 ≠ x2 and no infinity operand. For windows k ≤ 11, a partial
sum s (|s| < 2^(20k)) and an entry e = ±m·2^(20k)·B (|m 2^(20k)| ≤ 2^(19+20k) < r/2) coincide only if s = ±m·2^(20k) as
integers, which forces all lower digits to vanish and then m = 0, a contradiction. Sums never vanish. For window 12, the
congruence s ≡ ±m 2^240 (mod r) has at most about 2^20 bad values of s, so the heuristic failure probability over at most
2^117 points is at most 2^−100 under the assumption that the low digits are uniform. A violation aborts the procedure; none
occurred in any measured run. Table construction uses stages j ≥ 2, and j·b ≠ ±b for j ∈ [2, 2^19), so it has no
coincidences.

**Records and rounds.** RECORDS = 4096 records advance in rounds, each taking one walk point; even records are tame (base
P0) and odd records wild (base A_d). A trail is capped at 2^100 steps, then restarted from a fresh random start exponent
(an unbatched comb evaluation).

**Distinguished points.** A key whose low 90 bits are zero is distinguished (probability 2^−90), stored as a 72-byte record
(the exact high 291 bits of the key, the exponent, the family and padding). Lookup matches the 37-byte prefix exactly.

**Table and meter.** The table has 2^20 buckets of 256 slots, a capacity of 2^27 stored points, and aborts on overflow.
The procedure halts when its worst-case-charged total reaches W_max, so the claimed total bounds the work by construction.

## 3. Points guaranteed

B_min = S(13/20) + 2^112 points, S(τ) = ceil(2·sqrt(τM)). The detection span 2^112 is RECORDS × 2^90 × 2^10. Restarts never
exceed floor(B_min/2^100) = 120,221, and the trails number 124,317.

## 4. Success probability and losses

S is the number of walk points (both families) before a collision is detected. Under H1 (random mapping on M classes),
F(S) = 1 − exp(−S²/(4M)); at S(13/20), F = 0.477954.

- **Round-robin detection.** A merged pair's detection costs up to 4096 points per round. The budget 2^112 covers it except
  with probability 2·e^{−1024}.
- **Cap truncation.** A merge near a trail's 2^100-step cap loses at most 2^−10 per trail, 2^−9 for the pair (charged).
- **Trails without a distinguished point.** Negligible (e^{−1024} per trail over 124,317 trails).
- **In-trail cycles.** At most 1.2 × 10^−5 (charged).
- **Table overflow.** At most 1.8 × 10^−16 (charged; Chernoff bound on bucket loads of mean 128).
- **Window-12 coincidence.** At most 2^−100 (heuristic; section 2).

Total loss 0.0019651. The success probability is at least 0.4227 × (1 − 0.0019651) = 0.42187 ≥ 0.40, with a margin of 0.0219.
The 0.4227 is the Wilson 95% lower bound of 4000 fresh toy trials (section 7).

## 5. Cost

Work is W_max = H + R/C with C = 2000 (the G1 price in COST_MODEL.md). H counts reference units: one per batched addition
(12 per walk point, one per table addition, 12 per initial record). R counts all other scalar instructions.

**Reference core.** Each batched addition costs 2,204.82 instructions for its six field multiplications/squarings (the
wrappers `core_mul`, `core_sqr`), measured on this build. The organizer's convention prices it at 1 unit, against C = 2000:
this is a rounding of the convention, stated here openly. Charging each addition at its measured core would raise the
score by about 0.10 bits; the charge follows the convention.

**Worst-case other work per point.** Measured on the forced path: every digit nonzero (so each point performs its 12
additions, one copy, and all 12 stages occupied, which bounds every batch-boundary and inversion case), with the
canonicalisation, the digit extraction, the negation and the edge check all branch-free:

| Quantity | Value |
| --- | ---: |
| additions per point (forced) | 12 |
| other work per point, forced maximum | **8,658.35** instructions (4.33 units) |
| other work per point, natural mean (sanity: must not exceed the maximum) | 8,622.35 |

The forced run includes the forcing code, so it overcharges the production path slightly. The remaining data-dependent
paths are covered: the zero-digit skip and the copy (the forced run is the maximum, since skipping never adds work), the
batch boundary (at most 12 inversions per round, all present in the forced run), the DP test and the restart test (global
reserves, below), and the table select (pointer arithmetic, no branch).

**Cross-check of the tables.** The forced round was also measured with the full 2^19-entry tables: 143,835,734 instructions
in both cases (difference 0). The 4096-entry tables used in the measurement therefore give the same cost as the full tables.

**Restarts (global reserve).** Each restart is an unbatched comb evaluation of 12 additions, measured on the forced path at
415,422 instructions per restart. At most 120,221 restarts, charged at 415,422 each, plus a fresh start exponent (sampling
and 20 squarings, bounded at 4,000 + 20 × 1,000 instructions).

**Distinguished points (global reserve).** Each stored point costs at most 16,000 instructions (analytic bound: hashing, a scan of
at most 256 slots, a 72-byte store). The procedure stores at most 2^27 points, so the reserve is 2^27 × 16,000.

**Table build.** 26 batched additions per stage (2 bases × 13 windows), 2^19 − 2 stages, each with one inversion. The
per-stage other work is measured at 37,733 instructions and the fixed part at 3,077,891 instructions (doublings and the
first two entries); the total other work is therefore 3,077,891 + (2^19 − 2) × 37,733 = 19,785,849,401 instructions. These
are measured by differencing 1,024 and 2,048 stages.

**Initial records.** 4096 batched comb evaluations of 12 additions (H = 49,152), other work bounded at one point's worst-case
charge per record.

**Totals.**

| Item | value |
| --- | ---: |
| H (reference units) | 2^120.460292 |
| R/C | 2^118.989423 |
| preprocessing (units) | 2^24.566709 |
| **W_max** | **120.904710** (= 12090471/100000) |

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
  the production round (full tables), the restart path (4096 restarts) and the DP path (250 distinguished points). The replica
  computes table entries by real scalar multiplication.
- **Measurement.** `evidence/measure-v5w20/summary.json` and `run.log`: the reference core, the natural production round (sanity),
  the forced round and the forced restart round (4096 and 8192 points), the full-table cross-check, and the table build.

## 8. Round-5 objection and how it is answered

The worst-case charge for the per-point other work is the measured maximum over the forced path, with every data-dependent
branch made branch-free or covered by the forced run (section 5). The batch-inversion overhead per point is bounded by the
forced run's 12 inversions per round, whatever the batch sizes, since a round has at most one inversion per stage. The
reference rounding is stated openly (section 5).

## 9. Window and memory

The window is 20, as the user chose: 13 windows, 12 additions per point, 1.31 GB of comb tables. The earlier window-16
package (121.196226 with its worst-case charge) is 0.29 bits worse than this one.

## 10. Limitations

- H1 is an unproved random-mapping heuristic. The toy groups have 32 and 36 bits.
- The reference unit is the convention's 2000 instructions; the measured core is 2,204.82 per addition.
- The forced run includes its forcing code, a conservative overcharge.
- The restart path uses 16 precomputed exponents in the measurement; the charge covers fresh sampling.
- The window-12 coincidence and the bucket-overflow bound are heuristic or bounded as stated.
- Measurement: one compiler (GCC 14.2.0), one pinned blst, Valgrind 3.24.0, one input set.

## 11. Reproduction

```sh
python3 -B -I accounting.py
echo '{"bits":32,"k":6,"trials":60,"seed":2601}' | python3 -B -I experiment.py
```

The measurement: `gen_inputs_v4.py <target-profile.json> 7` regenerates the input file; the organizer harness runs
`cheon_step_w20.c` with `REF_FUNCS="core_mul core_sqr" N1=4096 N2=8192 MODES="forced forcedcap1 forcedfull"` and
`MODES="build" MODE_N1=1024 MODE_N2=2048`; `summarize_w20.py <run1> <run2>` computes the summary;
`replica_w20.py <target-profile.json> inputs.txt <n> <cfg> [dump]` checks a full-table dump.
