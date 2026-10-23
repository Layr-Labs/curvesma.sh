# bls12-381-g1-ecdlp: Cheon subgroup walk with fixed-base comb tables (v3)

Claim: for s uniform in Z_r^* and inputs T = [s]P0, A_i = [s^i]P0 (i = 1..n, n = 2^20), the
algorithm below recovers s with probability at least **0.40**, and it is run under a hard
work meter that stops it when its total reaches **W_max = 2^120.911517** target-group-operation
equivalents at n = 2^20 (score 120911517/1000000). W_max includes preprocessing, the walk, all
distinguished-point storage and lookups, restarts, the square-root chain and the output, at
their worst-case bounds. Peak memory is at most 2^34.178981 bytes. Preprocessing is at most
2^22.840069 units. Nonuniform advice is 0 bytes.

This is a classical cold attack on one target. It is not a practical break. The success bound
rests on the random-mapping heuristic H1 (section 4), on the fresh toy runs (section 7), and on
the analytic loss terms of section 4, which are stated as bounds and not measured.

Files: `accounting.py` computes every number here with exact integers (vendored `exact.py`);
`evidence/accounting.json` is its output. `evidence/measure-v3/` holds the kernel
(`cheon_step_v3.c`), its Python replica (`replica_v3.py`), the EC module (`ec.py`), the input
generator, the input file, the measurement log and the summary. `experiment.py` is the bundled
toy experiment. `evidence/toy-fresh/` holds the fresh toy runs and their driver.

## 1. Reduction

Let r be the group order, d = 2^20 (which divides r−1 since r−1 = 2^32·odd), and H = {x^d : x in
Z_r^*} the subgroup of d-th powers, |H| = N = (r−1)/d. Since s is uniform in Z_r^*, h = s^d is
uniform in H, and A_d = [h]P0.

Given h, recover s by twenty square roots in Z_r. With t_20 = h, for j = 19 down to 0 take a root
t_j of t_(j+1) with [t_j]P0 = A_(2^j). The true s^(2^j) is one of the two roots, and the other
fails the check. Finally s = t_0, checked against T = A_1. The square-root chain costs 41 comb
evaluations and twenty Tonelli–Shanks roots in Fr; 2-adicity 32 bounds a root by 2·32² modular
multiplications (section 5).

## 2. The walk on H modulo mu6

mu6 = {1, −1, λ, −λ, λ², −λ²}, with λ a cube root of unity mod r matched to β (so (βx, y) = [λ](x, y)
on G1). Since 6 | N, mu6 lies in H, and the walk runs on M = N/6 = 2^232.272 classes.

**Canonical point and key.** Among the six images of P under mu6, the canonical point is the one
whose x-coordinate is least, comparing Montgomery representations, with y of even parity (parity
of the standard representation). A canonical point is determined by its x-coordinate, so its
**key** is its Montgomery x-coordinate, a 381-bit integer.

**Multipliers.** 64 elements g_1..g_64 of H, each x^d for a random x in Z_r^* (20 squarings).

**Step.** For a canonical point with key k and tracked exponent c (so the point is [c]B):
- j = bits 90..95 of k, i.e. (k ≫ 90) mod 64, selects the multiplier;
- c2 = c·g_j mod r;
- the point [c2]B is computed by a comb: 16 signed 16-bit digits d_k with carry, magnitude m_k;
  a zero digit contributes nothing; magnitude m ≥ 1 adds the table entry m·2^(16k)·B (negated for a
  negative digit). The table of each base holds, for each window k, the entries m·2^(16k)·B for
  m = 1..2^15, so the table has 16 × 2^15 affine points per base;
- the result is canonicalised, the exponent becomes c ← μ·c2 for the mu6 element μ used.

**Bases.** Tame trails use B = P0 (public). Wild trails use B = A_d. The measurement uses
B = Q = [e]P0 in its place; the instruction path is the same.

**Records and rounds.** The walk keeps 256 records in round-robin order (even records tame, odd
wild). One round advances each record by one walk point. Each record is one trail.

**Caps and restarts.** A trail is capped at 2^100 steps. When a record's trail reaches the cap, it
is restarted from a fresh random start exponent c_0 in H (sampled as x^d) and the point [c_0]B,
canonicalised. Restarts are charged in section 5.

**Distinguished points.** A canonical point is distinguished when its key has its low 90 bits equal
to zero, which happens with probability 2^−90. Each distinguished point is stored as a 72-byte record:
the high 291 bits of the key (37 bytes; the other 90 bits are zero, so the key is recovered
exactly), the exponent (32 bytes), the family (1 byte), and 2 bytes of padding. Storage and lookup
use the bounded table of section 3. A lookup matches a stored record only when the 37-byte prefixes
are equal, which is an exact key match: two distinguished points with equal prefixes have equal keys.
No fingerprint is used, and the proof never infers point equality from a shorter match.

**Collision.** A tame point (exponent c_t) and a wild point (exponent c_w) with equal keys satisfy
[c_t]P0 = [c_w]A_d = [c_w·h]P0, so h = c_t·c_w^{−1} mod r exactly. The candidate is verified by
[h]P0 = A_d before the chain runs. Tame–tame and wild–wild matches carry no information and are
ignored.

## 3. The bounded distinguished-point table and the work meter

**Table.** 2^20 buckets of 256 slots, each slot a 72-byte record. A distinguished point is hashed
(a keyed hash of its 37-byte prefix in the production procedure; the measured kernel uses a multiplicative hash) into one bucket and scanned
against the occupied slots of that bucket, at most 256 comparisons. If the bucket has a free slot,
the record is stored; if the bucket is full, the procedure aborts and reports failure. The total
number of stored points is capped at DP_CAP = 2^27 (half the 2^28 slots); exceeding it also aborts.
There is no probing beyond the bucket.

**Failure of the table.** For a uniform hash, the number of points in a bucket, given at most DP_CAP
stored points, is binomial with mean 2^7 = 128 (DP_CAP/2^20). By the Chernoff bound
P(X ≥ 257) ≤ e^{−λ}(eλ/257)^257 with λ = 128, and a union bound over the 2^20 buckets, the
probability of a bucket overflow is at most 2^−52 (computed in `accounting.py`). The probability
that more than 2^27 points are distinguished is far smaller, since the expected number of points is
about 2^26.83. So the abort probability is at most 2^−52 and is charged in section 4.

**Work per distinguished point.** A bound of 16,000 instructions covers hashing, a scan of at most
256 slots at 60 instructions each, and the 72-byte store. The worst case charges DP_CAP points at
this bound.

**Meter.** The procedure keeps a running total of its work, using each item's bound: 16 units per
walk point and the measured per-point instructions, each restart's cost, each distinguished point's
bound, and preprocessing and the chain. It halts when the total reaches W_max minus the chain's
reserve. The total is therefore an upper bound on the work by construction, and every "at most"
claim holds for the procedure as run.

**Points guaranteed.** Let B_min = S(τ) + 2^108 with τ = 13/20, S(τ) = ceil(2·sqrt(τM)). W_max is
defined so that the meter allows at least B_min walk points even when every bound is at its worst
case. The number of restarts never exceeds floor(B_min/2^100) = 116,381, and the number of trails
is 256 + 116,381 = 116,637.

## 4. Success probability and losses

Let S be the number of walk points, both families, before a collision is detected.

**H1 (random mapping, disclosed).** The class walk behaves like a random mapping on the M classes,
so S points contain a tame–wild collision with probability F(S) = 1 − exp(−S²/(4M)). At
S(13/20), F = 1 − e^{−13/20} = 0.477954.

**Round-robin detection (H2).** Two trails that merge follow the same class sequence. The first
distinguished point on the merged path, at distance G steps from the merge, is reached by both.
The walk advances all 256 records in rounds, so the two trails' G steps cost about 256·G total
points, not 2·G. G is geometric with mean 2^90, and the bound used is G_max = max over the two
trails. The collision is detected within the budget B_min if 256·G_max ≤ B_min − S, which fails with
probability at most 2·exp(−(B_min − S)/(256·2^90)) = 2·exp(−2^108/2^98) = 2·e^{−1024}, far below 2^−1400.
The derivation uses 256 records, the value fixed by the procedure.

**Cap truncation (H2).** A merge at position p of a trail (p steps into its cap of 2^100) is detected
only if the next distinguished point on the merged path is reached before the cap. That needs
G ≤ 2^100 − p for each of the two trails. Averaging over p uniform in [0, 2^100), the probability that a
trail reaches its cap first is (μ/CAP)(1 − e^{−CAP/μ}) with μ = 2^90, so it is at most 2^−10 per trail.
For the pair it is at most 2·2^−10 = 2^−9. This loss is charged. The earlier analysis omitted it.

**Trails with no distinguished point (H2).** A trail of 2^100 steps with no distinguished point has
probability e^{−1024}; over 116,637 trails this is negligible.

**In-trail cycles (H2).** A trail of length 2^100 repeats a class with probability at most
(2^100)²/(2M) ≈ 2^−33.27. Over the 116,637 trails the total is at most 1.1 × 10^−5 ≈ 2^−16.4.

**Table overflow (section 3).** At most 2^−52.

**Statistics.** The fresh toy runs (section 7) measured success at S(13/20) with Wilson 95% lower bound
0.4227. The toy model has no cap truncation, no round-robin factor and no bounded table, so these
losses are charged on top of the toy rate as multiplicative factors.

**Claim.** The success probability is at least 0.4227 × (1 − 0.0019644) = 0.42187 ≥ 0.40, under H1 and
the losses above. The total loss is 0.0019644, dominated by cap truncation (0.0019531). The claimed
bound 0.40 has a margin of 0.0219.

## 5. Cost

Work is W_max = H + R/C with C = 2000 (G1 price, `tracks/docs/COST_MODEL.md`). H counts reference
units: one per batched incomplete affine addition or doubling, and each walk point and restart is
charged 16. R counts all other scalar AArch64 instructions, including loads, stores and control.

**Measured per walk point** (`evidence/measure-v3/summary.json`, run with the organizer's
Callgrind harness, 512 and 256 points over 256 records, production tables):

| Quantity | Value |
| --- | ---: |
| total instructions per point (differenced) | 518,340.24 |
| of which 15 reference additions (`ref_add`, 34,429.65 each) | 516,444.79 |
| other instructions per point (R) | **1,895.45** (0.94773 units) |

The 15 additions per point are executed as real affine additions, each with its own inversion. Each is
charged as one reference unit, so 16 units per point (one more than the 15 executed; a conservative
allowance for the first window). The other instructions cover the multiplier exponent product, the
digit extraction and table reads, the canonicalisation, the exponent update and the walk-state store.

**Measured DP and restart paths** (configurations with DP bits 4 and with cap 1, differenced):
- a distinguished point's store costs 762.5 instructions (measured; the charge is the bound 16,000);
- a restart's comb evaluation and canonicalisation cost 1,692.0 instructions of R plus 15 reference
  additions (charged as 16 units).

**Preprocessing.**

| Item | H (units) | R (instructions) |
| --- | ---: | ---: |
| table additions, measured count 1,048,544 | 1,048,544 | — |
| base shifts, 2 × 16 × 16 doublings | 512 | — |
| batch normalisation of every table entry (bound 10,000 each, 2^21 entries) | — | 1.049 × 10^10 |
| initial records: comb evaluations, 256 × 16 units; one point's other work each | 4,096 | 4.85 × 10^5 |
| fresh start exponents for the 256 records (sampling and 20 squarings, bound) | — | 6.1 × 10^6 |
| multipliers, 64 (sampling and 20 squarings, bound) | — | 1.5 × 10^6 |
| DP table zeroing, one store per 8-byte word of 2^28 slots | — | 2.42 × 10^9 |
| constants and input parse (bound) | — | 6.6 × 10^4 |

The table build's non-addition work is charged at the analytic bound of 10,000 instructions per entry.
The measured build total (8,429,431,517 instructions for 1,048,544 additions) and the reference
toggle agree to within 233 instructions, which is not plausible for the batch normalisation, so we do
not use that difference. The raw figures are kept in `evidence/measure-v3/summary.json`.

**Walk and finish.**

| Item | H (units) | R (instructions) |
| --- | ---: | ---: |
| walk points: 16 units each, B_min points | 16·B_min ≈ 2^120.8285 | — |
| walk points: other instructions, 1,895.45 each | — | 2.80 × 10^38 |
| restarts: 16 units each, at most 116,381 | 1,862,096 | — |
| restarts: other instructions, 1,692.0 each | — | 1.97 × 10^8 |
| restarts: fresh start exponent (bound) | — | 2.79 × 10^9 |
| distinguished points: worst case, 2^27 × 16,000 | — | 2.15 × 10^12 |
| chain: 41 comb checks | 656 | — |
| chain: one point's other work per check | — | 7.8 × 10^4 |
| chain: 20 Tonelli–Shanks roots (bound, 2·32² multiplications at 1000) | — | 4.10 × 10^7 |
| collision solve (bound) | — | 2.6 × 10^5 |

Totals: H = 2^120.828496; R/C = 2^116.751039; preprocessing 2^22.840069 units; and
W_max = H + R/C, so **time_log2 = 120.911517** (120911517/1000000).

**Range of n.** The bound uses A_1, A_(2^j) for j ≤ 20 and A_d = A_(2^20), so it holds for every n with
2^20 ≤ n < r at the same W_max. The scored value is W_max at n = 2^20; for n < 2^20 nothing is claimed.

**Accounting convention.** Every item is charged once. Reference units are priced by the organizer's
convention and their instructions are not charged again in R. The walk's additions are priced as
batched affine additions under the 6M convention. **The kernel does not batch them:** each addition
computes its own inversion, and that cost lies inside the reference unit. This is a reliance on the
convention, not a measurement of a batched implementation.

## 6. Memory

| Item | bytes |
| --- | ---: |
| DP table: 2^20 buckets × (4 + 256 × 72) bytes | 19,331,547,136 |
| comb tables: 2 × 16 × 2^15 affine points × 96 bytes | 100,663,296 |
| records (256 × 256 bytes) | 65,536 |
| code, stacks and scratch | 16,777,216 |
| **total** | **19,449,053,184 (2^34.178981)** |

Nonuniform advice is 0 bytes: the comb tables are built from the public base points at run time and
charged in section 5.

## 7. Evidence

**Fresh toy runs (unchanged from v2).** `evidence/toy-fresh/`: 4000 trials, 32-bit (k = 6, 3000 trials) and
36-bit (k = 10, 1000 trials) groups, driven by `toy_stats.py`. All 4000 secrets were recovered and no
collision gave a wrong h. At τ = 13/20: 1752/4000 = 0.438, Wilson 95% [0.4227, 0.4534]. The toy runs
do not include the production losses of section 4.

**Kernel and replica.** `evidence/measure-v3/cheon_step_v3.c` and `replica_v3.py` agree exactly on the
final state of all 256 records for these configurations: production 64 and 512 points, the
distinguished-point configuration (31 distinguished points in 512 points), the restart configuration
(64 restarts in 64 points), and the dump mode. The replica computes the table entries on demand by real
scalar multiplication, so it does not share a table with the kernel. A separate check confirmed that
each record's point equals [c]B for its exponent c, with B the record's base point, for all 256
records after 40 walk points.

**Measurement.** `evidence/measure-v3/run.log` and `summary.json`. Counts: 256 and 512 points, 7680
and 11520 reference-addition calls (the difference is the 256 × 15 steps; the 3840 extra calls are
the initial records); the restart configuration has 19,200 calls for 512 points, which is the 3840
initial calls plus 15 per restart.

## 8. Round-3 objections and how they are answered

1. **Kernel surrogate.** The kernel now uses real affine accumulation through `ref_add`, a real group
   addition; the comb tables are the production structure, 16 windows of 2^15 entries per base, built in
   the kernel from the base points by real projective additions and batch normalisation (preprocessing,
   so no smaller-table argument is needed); the key bits are 90–95; the trail counter is 128 bits with
   cap 2^100; the restart and distinguished-point paths are real code, measured in separate configurations
   and differenced. Two differences remain, stated here: the restart start exponents in the measurement
   are 16 precomputed values, while the production algorithm samples a fresh one per restart (the same
   instruction path; the sampling cost is charged as a bound); and the kernel's additions are unbatched,
   each with its own inversion, charged as reference units (section 5). **Status: the measurement now
   matches the specified step except these two points, which are handled as stated.**
2. **H2 losses.** The cap truncation (2^−9) is derived and charged. The detection bound uses the
   round-robin factor, 256 points per round, so the detection span is 2^108 (not 2^101) and the loss is
   2·e^{−1024}. The earlier e^{−512} figure is withdrawn.
3. **DP table.** The table is bounded: 2^20 buckets of 256 slots, a capacity of 2^27 stored points, abort on
   overflow, no probing beyond the bucket, a worst-case work bound of 16,000 instructions per point, and the
   failure probability 2^−52 charged.
4. **Hard cutoff.** The procedure runs under a meter that halts at W_max; W_max is the claimed total,
   computed from worst-case bounds, and the stopping loss is the detection and cap losses above.

**Not resolved by this round:** the batched-addition question (section 5), the sampled-versus-precomputed
restart exponents (objection 1), and the unmeasured analytic bounds (normalisation, sampling, Tonelli,
DP work). They are reported as bounds, not measurements.

## 9. Limitations

- H1 is the random-mapping heuristic for the class walk. The toy groups have 32 and 36 bits; nothing is
  tested at production size.
- The per-point count is one compiler (GCC 14.2.0 in the organizer's image), one pinned blst, 512 points
  over 256 records; the replica agrees exactly.
- The bucket-overflow bound assumes that the bucket index of distinct prefixes is uniform (a keyed hash
  in production), which is an additional assumption beyond H1.
- The reference-unit price is the organizer's; the kernel's own additions are unbatched (section 5).

## 10. Reproduction

```sh
python3 -B -I accounting.py
echo '{"bits":32,"k":6,"trials":60,"seed":2601}' | python3 -B -I experiment.py
```

The measurement: `evidence/measure-v3/gen_inputs_v3.py <target-profile.json> 7` regenerates the input
file (`inputs.txt`, 6,075 bytes, included). The organizer harness runs `cheon_step_v3.c` with
`REF_FUNCS="ref_add ref_proj_add" N1=256 N2=512 MODES="dp4 cap1 build" MODE_N1=256 MODE_N2=512`.
`replica_v3.py <target-profile.json> inputs.txt <n> <prod|dp4|cap1> [dump]` checks a dump.
