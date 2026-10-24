# Kernel measurements behind C_HIT and C_ENTRY

The track metric counts executed scalar AArch64 instructions. The two dominant constants of the
cost model are **retired-instruction counts** of the kernels in this package. They were measured
with the OS counter (`/usr/bin/time -l`, "instructions retired") on Apple-silicon AArch64 with
clang -O2 and no fast-math. These are C sources: the verifier's Python sandbox does not compile
them, so they are evidence to read and re-run on an AArch64 host. `evidence.py` part [8]
re-counts the branch-free LA loop from the shipped assembly.

## C_HIT = 150: large-ideal hit (`archd_hit.c`, `archd_lattice_d12.h`)

**Benchmark lattice.** A generic index-2^40 sublattice of Z^12 (the kernel of a random form mod
2^40, LLL-reduced, isotropic |b*| ≈ 10). It is enumerated inside the fitted 12-dim ellipsoid
(`LAT_Q` is the submission's matrix) with a Schnorr–Euchner DFS: 2,097,108 points, matching the
Gaussian heuristic to 0.03%.

**Earlier validation of the enumeration** (unchanged code paths):
- exact PARI `qfminim` agreement at three radii;
- brute-force agreement in dimensions 4 and 6;
- no duplicates;
- all points inside the ellipsoid by the direct form.

**What a real hit must do, and what the kernel does.** The original benchmark (`hit9`, kept as
mode 1) built its "slab id" from the walk's lattice-basis coordinates z. A sieve must instead
address the hit by its ambient (a, b) coefficients. The realistic kernel (`enum_fast_real`,
`enum_fast_send`, `HIT_REAL_REGS`, `HIT_SEND_REGS`):

1. Carries 7 ambient state words per level, all affine in v = Σ z_i b_i:
   - the dense slab index, split into its a-part and b-part;
   - the sign keys a5 and b5;
   - the 3 in-slab coordinates.

   Each is updated with one add per node. The last level is a register loop over its exact
   interval: every z0 in [lo0, hi0] is a hit, exactly as in the reference `enumerate()`.
2. Canonicalises the point under (a, b) → (−a, −b) and the automorphism (a, b) → (a, −b): it
   flips the a-part if a5 < 0 and the b-part if b5 < 0, using X(−a) = C − X(a).
3. Loads the canonical slab's descriptor {offset, w0, w01} and computes the cell.
4. Either updates the byte bucket directly (single memory), or appends {slab, cell, logp} to an
   outgoing buffer, which `apply_buf` applies at the owner (distributed).

**Correctness:**
- `./archd_hit stcheck` compares the carried state with the directly computed v at every hit:
  0 mismatches over 2,097,108 hits.
- The naive reference (state in memory, mode 3), the optimised single-memory kernel (mode 23)
  and the optimised distributed kernel with apply (mode 24 + apply) leave a bit-identical bucket
  array (FNV digest `f704108504039503` after 2 repetitions).

**Retired instructions per enumerated point** = (count at 8 repetitions − count at 0) /
(8 × 2,097,108):

| mode | what | instr/point |
|---|---|---|
| original `enumerate`, no hit | reference walk | 92.4 |
| original `enumerate` + `hit9` | the v1 claim's kernel (claimed 103 static) | 122.4 |
| 20 | optimised walk, state carried, no hit work | 74.5 |
| 23 | realistic hit, single shared memory | 131.1 |
| 24 | realistic hit, send to slab owner | 133.3 |
| 24 + apply | send + apply at owner | **140.1** |

Repeat runs agree to ±0.01. 2^121.6 bytes of buckets are not 64-bit addressable, so the
distributed figure is the one charged. One intermediate routing hop (store, then reload and
forward, ≈ the 6.8-instruction apply delta) brings it to ≈ 147. **Claimed: C_HIT = 150.**

Production ideal lattices hold 2^53–2^92 points each, versus ~2^21 here. The walk is then even
more leaf-dominated (1.75 node visits per point here), so per-point internal-node overhead is
lower in production.

Reproduce:

    clang -O2 -o archd_hit archd_hit.c -lm
    ./archd_hit stcheck
    ./archd_hit icount 2 23 ; ./archd_hit icount 2 24 1          # same digest
    for m in 20 23 24; do /usr/bin/time -l ./archd_hit icount 0 $m; /usr/bin/time -l ./archd_hit icount 8 $m; done
    /usr/bin/time -l ./archd_hit icount 0 24 1; /usr/bin/time -l ./archd_hit icount 8 24 1

(In zsh pass the mode and the apply flag as separate words.)

## C_ENTRY = 35: block-Wiedemann entry (`la_mac.c`, `la_mac.s`, `la_finish.c`, `la_test.c`, `la_consts.h`)

**Storage.** Matrix entries after filtering are small nonzero integers (|c| < 256). Each row is
stored sign-split: its positive entries, then its negative entries, as uint8 magnitudes.

**Product.** `mac_run` accumulates Σ c·x[col] into a lazy 5-limb accumulator, run once for each
sign. `row_finish` reduces both accumulators into [0, r) by a Barrett-type step (t = v >> 255,
q̂ = (t·⌊2^275/r⌋) >> 20, ≤ 3 conditional subtractions) and returns (pos − neg) mod r.

**Costs:**
- The `mac_run` inner loop (`LBB0_2` in `la_mac.s`) is 33 instructions and branch-free apart
  from its back-edge, so it retires exactly 33 per entry. `evidence.py` part [8] re-counts this
  from the file.
- Per row: 2 loop prologues/epilogues (11 each) + `row_finish` (68 static, one branch) + calls
  ≈ 100 instructions, i.e. 0.5 per entry at weight 200.
- Total 33.5 per entry; **claimed: C_ENTRY = 35**.

**Correctness.** `la_test.c` runs 2,000 random rows: weights 1–400, including 100 rows of
all-255 coefficients, over 1,000 vectors in [0, r) including r − 1. All results match Python
big integers.

The v1 kernel this replaces was not sign-correct: its compiled sign mask was constant 0, so
negative entries were added.

    clang -O2 -o la_test la_test.c la_mac.c la_finish.c
    # x.txt: 1000 lines of 4 hex limbs (little-endian) of values in [0, r)
    ./la_test x.txt > out.txt      # each line: w npos idx:val ... | result (hex); check vs bigints

## Constants not re-measured here

C_UPD = 14, C_SET = 141, C_SCAN = 10 and C_REP = 6 are static counts of the small-ideal slab
kernels. They are charged at **× 1.25** (17.5, 176, 12.5, 7.5), which pads the static-to-retired
gap measured above (122.4 / 103 = 1.19). Together they are 31% of relation collection; a further
1.5× on all four (26.25, 264, 18.75, 11.25) moves the model total from 114.482 to 114.612
(+0.13 bit).
