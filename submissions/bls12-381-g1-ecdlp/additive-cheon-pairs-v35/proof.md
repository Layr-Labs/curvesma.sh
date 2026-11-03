# 971-instruction side-output kernel and coverage retune on the bc91c5ff paired walk

## Claim

For the protected powers-of-tau target at n=2^20, the paired class-walk Cheon algorithm of source-proof.md, as amended by base-proof.md (submission bc91c5ff, kept verbatim), has total work at most 2^116.277604948446796 G1-reference units (C=2000). This holds under the disclosed hypotheses H1-H6 (source-proof.md section 7, base-proof.md section 6).
- **Score:** 116277604948446796/10^15.
- **Success:** at least 0.39.
- **Memory:** below 2^121.5854 bytes, with every node below 2^63.
- **Preprocessing:** below 2^115.2868 units, included.
- **Advice:** none.

There are three changes against bc91c5ff; everything else in base-proof.md is unchanged.
1. **Side-output kernel** (section 1). The function that computes P-T and its key is replaced. Its description replaces base-proof.md section 1 and the earlier kernel descriptions in source-proof.md sections 5.1 and 8.
2. **Coverage target**, as in submission 86d147e9. KT = 495119536/10^9 gives the exact floor 1 - exp(-K) - 0.0005 - 2^-20 - 4·2^-40 - 2^-111 = 0.390001000 ≥ 0.390001, computed by evidence.py in double precision with margin 3.4·10^-10 above 0.390001. N is the least integer with 6N^2 ≥ KT·m.

3. **Fixed-call row** (as in fd385c98, itemized): 40 instructions per kernel call instead of 512 = prologue and epilogue 32 (asserted by kernel_check.py) + terminator store, argument move, bl 3 + last exit branch 1 + margin 4 (the args cursor ends at args[n]; no caller-saved register is live across the call). Invocations (BASES+2)·LANE_OPS unchanged.

The algorithm and its outputs are unchanged: points, keys, routing and success events. Every number shown in base-proof.md (all sections) is superseded by section 2 and evidence.py; only formulas and hypotheses are inherited. Files cited but not shipped (witness scripts, generator, harness, H6 offline specs) are ancestor records, not relied on: the kernel charge rests on the static count of the shipped .S (kernel_check.py), Callgrind is corroboration, H6 rests on bc91c5ff. Submission ids are provenance only.

## 1. The kernel (measure/side_lanes.S): 971 instructions per lane, worst path

**Interface.** `side_lanes(args, n)` takes lane records {out, P, T, inv}. It writes the fully reduced x' to out+0 and y' to out+48. It writes the key to out+96: limbs 1..3 of whichever of y' and p-y' is even (the parity key of a38cb931). The loop ends at args[n] = {NULL, NULL}.

**Inputs.**
- x1, y1, x2, y2 are canonical (all stored points and table entries are fully reduced), and y != 0 on G1.
- inv = 1/(x2-x1) < p.
- Every gap-table entry also stores v2 = p-x2 and n2 = p-y2.

**Routines.** The kernel inlines blst f7a88067's `__mul_mont_384` (twice), `__sqr_384` and `__mul_by_1_mont_384`, with these edits:
- n0 is pinned in x18.
- The ret is dropped.
- The round-1 zero and the dead round-6 top-carry chain are deleted.
- Neither product does its final conditional subtraction.
- The square stays in registers, and p stays in x5..x10.

**Steps.** S(a,b) < ab/2^384 + p, with p/2^384 < 0.1016.
- **mu.** t = y1+y2 < 2p, then mu = S(t, inv) < 1.21p.
- **V.** V = S(mu²) < 1.15p.
- **x'.** X = V + v2 - x1 ∈ (-p, 2.15p). Add p once, or subtract p once or twice, using branches. Store x'.
- **w.** w = x' + v2 ∈ (0, 2p).
- **z.** z = S(w, mu) < 1.25p.
- **y'.** This is the line through P and -T: y' = mu(x'-x2) + y2, which equals mu(x'-x1) - y1. So y' ≡ D = z - n2 ∈ (-p, 1.25p). Add p, or subtract p once. Store y' and the key.

Accumulators stay below a + p < 3p and raw products below 2^384, so deleted carry words are zero.

**Duplicated continuation.** The code after x' (second product, y', key, loop tail) is emitted twice, so no x' path jumps to shared code.

**Instruction counts.** Per lane:

| region | instructions |
| --- | --- |
| first product | 316 |
| square | 92 |
| square reduction | 150 |
| second product | 316 |
| driver, worst path | 97 |
| total, worst path | 971 |

There are 24 paths, worst 971 and best 960. The prologue and epilogue total 32 per call; the last lane's exit branch is the +1 in 971n+1. Both are covered by the fixed-call row of 40 per call (change 3).

**Measurement.** measure/measurement.log (Callgrind, organizers' image): (Ir(8192) - Ir(4096))/4096 = 964.920 per lane, the mean; calibration 383. Its harness self-checks (8192 lanes vs blst's API, 64 genuine pairs, 8192 lanes biased to X ≥ 2p) passed; the harness is corroboration only.

**Hosted check.** kernel_check.py interprets the shipped .S (with a64.py) on 469 lanes against exact integer arithmetic, including 24 genuine curve pairs; rare-branch lanes are synthetic witnesses. It enumerates the 24 paths, and asserts:
- both outcomes of all 11 data-dependent branches occur;
- the lane-body count is ≤ 971n+1 and the fixed part is constant and ≤ 32;
- inputs and callee-saved registers are unchanged.

key_check.py checks the key blocks.

## 2. Ledger (evidence.py; other rows keep their formula, values recomputed at the new N)

| item | bc91c5ff | here |
| --- | --- | --- |
| kernel per lane | 1000 | 971 |
| pair (two points) | 3053 | 3024 |
| per point (F=197/200) | 1574.8525 | 1560.57 |
| KT | 0.49514 | 0.495119536 |
| N | 64239578809815195260251725009983698 | 64238251294037793068639116880053465 |
| table step | 27 | 48 (+3 ldp of p, 12 subtractions for v and n, 6 stp) |
| gap-table entry | 96 B | 192 B (gen node 2^61.70 B) |
| re-walk table rebuild | 2000+32 | 2000+50 |
| fixed-call row, per invocation | 512 | 40 |

Changed rows, in instructions:

| row | instructions |
| --- | --- |
| records | 200496575643873117458252293259010071751 |
| padding | 50390014992125786787272127997280256 |
| gap tables | 128476502588068283214693060591820800 |
| re-walks | 1269036660377395806558870306816 |
| side-output fixed calls | 83319366358042498007225220464640 |

Totals:
- W = 201408089141898897965986818227785118509 instructions.
- log2(W/2000) = 116.2776049484467956..., with ceiling 116277604948446796/10^15 (certified by log_bound.py).
- Preprocessing is 101344467732796792599227774576020440410 instructions.

Base-point entries get v and n (21) inside the 8192000-instruction scalar-multiplication charge (about 4x a double-and-add). The 192-byte stride indexes like 96 (3·2^k).

## 3. Limitations

Analytic bound; H1-H6 inherited extrapolations; no full-size run.
