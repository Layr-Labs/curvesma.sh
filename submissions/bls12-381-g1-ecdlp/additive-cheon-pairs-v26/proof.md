# Measured 972-instruction side-output kernel on the powers-of-tau paired walk

## Complete claim and scope

For the protected BLS12-381 powers-of-tau target at n=2^20, the paired class-walk Cheon construction of source-proof.md, as amended by base-proof.md, has total work at most 2^116.278065276150107 G1-reference units (C=2000), conditional on the disclosed hypotheses H1-H6 (source-proof.md section 7, base-proof.md section 6).
- The exact score is 116278065276150107/1000000000000000.
- Success is at least 0.39 under those hypotheses (floor 0.39000100, after the coverage retune of section 2).
- Memory is below 2^121.5854 bytes, with every node below 2^63 bytes. Preprocessing is below 2^115.2873 units and included. There is no advice.
- Only A_1 and A_(2^20) are accessed, and every access is charged. The same W(n) bound holds for 2^20<=n<r.

File base: base-proof.md is the complete argument of submission bc91c5ff (116.290720263097717), kept verbatim. This package changes two things: the function that computes P-T and its orbit key, and the coverage target K (section 2).
- Section 1 here replaces base-proof.md section 1 (kernel body and witness files). It also replaces the kernel descriptions in source-proof.md sections 5.1 and 8 and the kernel-audit sentence of base-proof.md section 5. Those passages describe earlier kernels (including base-proof.md's "36 fixed function instructions", "1000" and "458 lanes"): in this package, the files in measure/ are the new kernel, its Callgrind record, its generator and its harness.
- The kernel's two extra gap-table fields are charged in the table rows (section 2). All other rows, hypotheses and parameters of base-proof.md are unchanged except those listed in section 2.
- The algorithm, its points, keys, routing, candidates and success events are unchanged. The output is the same: (x', y') fully reduced, and the key is limbs 1..3 of whichever of y', p-y' is even.
- The base's files side_k2.S, kern2_witness.py, lazy_witness.py and native_check.py audited the replaced kernel and are not shipped. kernel_check.py and key_check.py audit the new one.

## 1. The kernel: measure/side_lanes.S (972 instructions per lane, worst path)

**Interface.** `side_lanes(args, n)` takes the same lane records as the base: args[i] = {out, P, T, inv}, with out receiving x' at +0, y' at +48 and the key at +96. One call processes the n lanes of a round after its batched inversion. args[n] = {NULL, NULL} ends the loop, and that store is inside the fixed-call row.

**Input preconditions.** The function relies on four facts about its inputs:
- x1, y1, x2, y2 are canonical (< p): they are walk points and gap-table entries, all stored fully reduced.
- y2 != 0, because G1 has no 2-torsion.
- inv < p: it is the reduced output of the batched inversion.
- The gap-table entry also stores the integers v2 = p - x2 (in [1, p]) and n2 = p - y2 (in [1, p)). Every T that reaches side_lanes is such an entry: pair operations, padding operations and candidate re-walks all read T[g] from tables built by the 48-instruction table step.

**Inlined blst routines.** The kernel inlines blst f7a88067's `__mul_mont_384` twice, `__sqr_384` and `__mul_by_1_mont_384`. measure/gen_side_lanes.py applies the following mechanical edits and asserts each one; measure/side_lanes.S.in is the driver.
- The ret is dropped.
- n0 is pinned in x18, so the per-round reloads are deleted.
- The round-1 `mov x17,xzr` is dropped; its consumer becomes `adc x4,xzr,xzr`.
- The dead round-6 top-carry chain is deleted, since both raw products are < 2^384.
- Both products lose their final conditional subtraction.
- The square's six result stores and the reduction's three input loads are dropped, and registers are renamed so the 768-bit square stays in registers.
- The square's temporaries are renamed, so p stays in x5..x10 for the whole lane.
- The second product reads mu sp-relative and does not touch x2 (out) or x29 (T).

**Static counts per lane.**

| region | instructions |
| --- | --- |
| first product (mu) | 316 |
| square | 92 |
| reduction of the square | 150 |
| second product (z) | 316 |
| driver, worst path | 98 |
| total, worst path | 972 |

**Bounds.** Let S(a,b) = (ab+qp)/2^384 be the raw word-serial Montgomery output, with q < 2^384. For a < 2^384 and b < p·1.21, the accumulator stays below a + p < 2^385 within the carry word, and S < ab/2^384 + p, where p/2^384 < 0.1016. The steps are:

1. t = y1+y2 < 2p, then mu = S(t, inv) < 1.21p.
2. V = S(mu^2) < 1.15p.
3. X = V + v2 - x1 ≡ x' (mod p), with X ∈ (-p, 2.15p). One sign branch adds p once, or subtracts p once or twice, each time with a borrow branch. x' is stored fully reduced.
4. w = x' + v2 ∈ (0, 2p), so w ≡ x' - x2.
5. z = S(w, mu) < 1.25p, unreduced.
6. y' = mu(x'-x2) + y2: this is the line through P and -T. It equals the base's mu(x'-x1) - y1, because mu(x2-x1) = y1+y2. So y' ≡ D = z - n2 ∈ (-p, 1.25p) (n2 <= p-1). One branch adds p, or subtracts p at most once. y' is stored fully reduced.
7. The key is limbs 1..3 of y' if y' is even, else of p - y': one tbz, plus 4 instructions on the odd path.

Because every raw output is < 2^384, the deleted top-carry words are identically zero.

**Price.** The lane body has 24 paths. kernel_check.py enumerates them from the shipped file: worst 972, best 961. The ledger charges 972 per lane, plus the base's separate fixed-call row of 512 per invocation, which also covers the one exit branch (b .Lexit) after the final lane on the D ≥ p path. kernel_check measures the actual fixed part (prologue, epilogue and terminator) at 32 instructions per call.

Callgrind (measure/measurement.log):
- Setup: the organizers' image (GCC 14.2.0, Valgrind 3.24.0), blst built as in scripts/reference_costs/run.sh, Ir toggled around one side_lanes(args,n) call.
- Result: (Ir(8192)-Ir(4096))/4096 = 965.513 per lane. That is the mean over the harness lanes, below the charged worst path. The bls_field calibration gives 383, as published.

**Correctness evidence.**
- (a) **Callgrind harness.** measure/side_kernel.c self-checks under Callgrind and natively. It compares all 120 output bytes of every lane with blst's exported reduced API on the same formula:
  - 8192 lanes: all {0, 1, p-1, p-2, random} combinations, then randoms;
  - 10^6 lazy-operand products with a < 4.15p and b < 1.21p;
  - 64 genuine G1 pairs, where x' and y' equal blst's affine P + (-T);
  - 8192 lanes with x1=0 and x2 in {0, p-1}. x2=0 gives X = V + p, so X ≥ 2p whenever V ≥ p.
- (b) **Instrumented native build.** Over the 8192-lane run plus (a)'s other lanes, an instrumented build counted 6956 entries to the X ∈ [p, 2p) block, 219 to X ≥ 2p, 6400 to D ≥ 0 and 59 to D ≥ p. All outputs matched.
- (c) **kernel_check.py (hosted).** It interprets the shipped .S with the base's a64.py on 469 lanes against exact integer field and curve arithmetic. These include 24 genuine curve pairs and 120 lanes biased to X ≥ 2p and D ≥ p. It counts every lane-body instruction and asserts four things:
  - both successors of all 9 conditional branches are reached;
  - each call's lane-body count is ≤ 972n + 1 (the final lane's exit branch);
  - a constant fixed part (32 instructions) per call;
  - unchanged inputs and callee-saved registers.
- (d) key_check.py checks the key's limb semantics and both key blocks of the shipped file.

## 2. Ledger changes (evidence.py; all other rows identical to the base)

**Coverage retune.** The base sets the coverage target at K = 6N^2/m >= 24757/50000 =
0.49514 and keeps a success floor of 0.39001347, i.e. an arithmetic margin of 1.3e-5
above the required 0.39. That margin is not required by the track. The floor is exact
arithmetic on K, floor(K) = 1 - exp(-K) - 0.0005 - 2^-20 - 4*2^-40 - 2^-111, so the
least coverage that still meets the requirement is obtained at
K_min(0.390001) = 0.4951195354440521..., and this package sets

    KT = 495119536/10^9 = 0.495119536

with the success assertion `succ >= 0.390001`. N is recomputed by evidence.py as the
least integer with 6N^2 >= KT*m (asserted in both directions), and every row that
depends on N is recomputed by the same code. No row formula, price, kernel, transport
charge, cap, allowance or heuristic changes. N goes from
64,239,578,809,815,195,260,251,725,009,983,698 to
64,238,251,294,037,793,068,639,116,880,053,465, and the total from
201,476,529,596,183,341,906,903,258,936,228,796,747 to
201,472,363,651,311,296,166,978,104,573,756,839,268 instructions.

| item | base | here |
| --- | --- | --- |
| KERNEL (per lane, worst path) | 1000 | 972 |
| pair, two points | 3053 | 3025 |
| per point (F=197/200) | 1574.8525 | 1561.0625 |
| table step listing TABLE_BK | 27 | 48: +3 ldp of p, +2 subs and 10 sbcs for v and n, +6 stp |
| gap-table entry / gen-node table memory | 96 B | 192 B (gen node 2^61.70 B, < 2^63) |
| re-walk table rebuild | 2000+32 | 2000+50 (asserted >= 2000+TABLE_BK) |

The changed rows, in instructions:

| row | instructions |
| --- | --- |
| records | 200559850321397744684424902789136924414 |
| padding at the pair price | 50406678356871860129463686240665600 |
| gap tables | 128476502588068283214693060591820800 |
| candidate re-walks | 1269036660377395806558870306816 |

Totals:
- W = 201472363651311296166978104573756839268 instructions, and log2(W/2000) = 116.2780652761501064...
- The 10^-15 ceiling, 116278065276150107/10^15, is certified by log_bound.py.
- Preprocessing is 101376283796246521510758074553421900522 instructions (< 2^115.2873 units).

The T[g] load in the pair and root listings is unchanged: the stride is 192 instead of 96 bytes (both 3·2^k), so the index takes the same instructions.

## 3. Limitations

This is an analytic bound with a measured, path-enumerated kernel. H1-H6, the toy and planner evidence, and their stated scopes are inherited unchanged from base-proof.md and source-proof.md. They remain extrapolations, and no full-size run has been done.
