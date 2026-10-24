# Kernel measurements behind C_HIT, C_SMALL and C_BIG

The track metric counts executed scalar AArch64 instructions. The dominant constants of the v3
cost model are **retired-instruction counts** of the kernels in this package, measured with
`/usr/bin/time -l` ("instructions retired") on Apple-silicon AArch64, clang -O2, no fast-math,
0-vs-8 repetition deltas over the same enumerated point set. These are C sources: the verifier's
Python sandbox does not compile them; they are evidence to read and re-run on an AArch64 host.
`evidence.py` part [8] re-counts the la_mac inner loop from the shipped assembly and takes a
census of la_big.s.

## C_HIT = 170: large-ideal hit with production addressing (`archd_hit2.c`)

Replaces v1's `archd_hit.c`, whose benchmark shortcuts the review rejected (slab table masked to
2^16 entries, apply ignoring the destination). Same walk/enumeration math as v1 (same benchmark
lattice `archd_lattice_d12.h`, same 2,097,108-point coset, 3,671,024 node visits), new hit path:

1. Carried ambient state (checked against direct v = sum z_i b_i at every hit: **0 mismatches**).
2. Canonicalisation under (a,b) -> (-a,-b) and sigma: (a,b) -> (a,-b), before the id is formed.
3. **Full-width slab id (v3 repair)**: the 9 fixed coordinates as signed balanced-radix
   digits (base **4097**, digits in [-2048, 2048]), packed into two uint64 words
   (60.0 + 48.0 bits; **108.0 bits total**). No truncation before hashing or key
   comparison. The production coordinate bound is **|v_i| <= 1699.1 at A = 2^10.00**
   (computed from the embedded covariance in `evidence.py` part [2c]); radix 4097
   covers +/-2048 with ~20% headroom. The v2 base 1497 (A = 2^9.55 = 748) aliased at
   the declared operating point: (749,0,0,0,0,1) and (-748,1,0,0,0,1) are both inside
   the ellipsoid and packed identically. The hot path does not multiply by the radix:
   the packed words are carried by the walk (precomputed `DLT` deltas) and a hit only
   negates them for canonicalisation, so changing 1497 -> 4097 changes setup constants
   and the key range, not the per-hit instruction sequence (retired count 161.8 still
   applies). The three in-slab coordinates use **exact strides** in production
   (`wa = 2B+1`, `wab = wa^2`, with the slab's bbox offsets), so id + cell is injective
   on all 12 coordinates.
4. **Descriptor lookup**: open-addressing hash table keyed by the full 128-bit (k0, k1) pair,
   exact comparison of both words, 2^22 slots at load 0.46 (1.93M distinct ids, prepopulated by
   an untimed pre-pass so measured probes hit; a miss is fatal in both the single-memory path
   and apply). Present-key probes: mean 1.42 slots; absent-key probes (control): mean 2.22.
5. **Routing**: node index in [0, 2^20) from a multiplicative hash of both id words; a real
   2^20-entry route table maps node -> one of 8 destination ring buffers (all 8 used, 262-263k
   messages each). The update message is 48 bytes {k0, k1, t0, t1, t2, logp}.
6. **Apply at the owner**: probes the descriptor with the full key, recomputes the cell
   `off + t0 + wa*t1 + wab*t2` in 64-bit arithmetic (values up to ~2^50; only the final
   byte-store index is masked into the benchmark's bucket array), adds logp.

**Correctness:** the reference (recomputing v from the (a,b) coefficients directly), the
optimised single-memory kernel, and the optimised distributed kernel + apply leave **bit-identical**
bucket arrays (FNV digest `d9f629fbb23665e3` over 2x2,097,108 hits). The two measured paths agree
(digest `9d9107728f9cb30f`).

**Retired instructions per enumerated point** ((count at 8 reps - count at 0 reps)/(8 x 2,097,108);
five independent process pairs per mode, spread < 0.15; a 16-rep cross-check agrees to < 0.05):

| mode | what | instr/point |
|---|---|---|
| 20 | walk, ambient state forced live (asm barrier), no hit work | 85.0 |
| 23 | full hit, single shared memory | 140.4 |
| 24 | full hit, send to per-destination buffers | 127.9 |
| 24+apply | send + owner-side apply (the charged path) | **161.8** |

Note on walk-only modes: v1's published 74.5 "state-carried walk" was dead-code-eliminated by
clang (removing its sink entirely gives 74.45); with a liveness barrier the walk is 85-91. The
full-hit modes are the robust figures — the state is genuinely consumed there.

**Claimed: C_HIT = 170** = 161.8 (send + apply, production addressing) + one intermediate routing
hop (store/reload/forward, ~ the 7-instr apply delta), rounded up.

Reproduce:

    clang -O2 -o archd_hit2 archd_hit2.c -lm
    ./archd_hit2 packcheck        # v3: v2 alias reproduced, 0 collisions at radix 4097
    ./archd_hit2 stcheck          # 0 mismatches, max |v_i| = 38
    ./archd_hit2 refcheck 2       # ref == single == distributed, bit-identical
    ./archd_hit2 stats            # table load, probe depths, destination histogram
    for m in 20 23 24; do /usr/bin/time -l ./archd_hit2 icount 0 $m; /usr/bin/time -l ./archd_hit2 icount 8 $m; done
    /usr/bin/time -l ./archd_hit2 icount 0 24 1; /usr/bin/time -l ./archd_hit2 icount 8 24 1

(First-touch fix: message buffers are pre-touched in every invocation so the 0-vs-8 delta is
linear; without it the first repetition pays ~43M instructions of page faults.)

## C_SMALL = 35: small-valuation entries (`la_mac.c/.s`, `la_finish.c`, `la_test.c` — v1, unchanged)

Rows stored sign-split as uint8 magnitudes (< 256); branch-free 33-instruction inner loop
(re-counted from `la_mac.s` in evidence.py part [8]); ~100 instructions of per-row overhead;
`la_test.c` checks 2,000 random rows against Python big integers: 0 mismatches.

## C_BIG = 147: full-size Schirokauer entries (`la_big.c`, `la_big.s`, `la_big_check.py`)

Each relation row carries 41 Schirokauer-map values in F_r (255-bit, 4 limbs). `la_big.c`
computes acc += a*b with a column-ordered 4x4 schoolbook product into a lazy 10-limb
accumulator and **one Barrett finish per row** (mu = floor(2^516/r), qhat = (acc>>258)*mu>>258,
<= 10 conditional subtractions). `la_big_check.py` verifies all 64 check rows **exactly against
Python big integers: 64/64 ok** (a carry-propagation bug in an early draft was caught by exactly
this check). Measured on 4096x48 entries:

- 0-vs-8 reps: 146.9 instr/entry; 0-vs-16: 146.88; 8-vs-16: 146.86 (all agree).

**Claimed: C_BIG = 147.** A fused Montgomery MAC would be somewhat cheaper; the schoolbook+
Barrett figure is charged (conservative).

Reproduce:

    clang -O2 -o la_big la_big.c && python3 la_big_check.py
    /usr/bin/time -l ./la_big bench 4096 48 0 ; /usr/bin/time -l ./la_big bench 4096 48 8
    # instr/entry = (count8 - count0) / (7 * 4096 * 48)   (the 0-rep baseline does one full pass)

## Constants not re-measured here

C_UPD = 14, C_SET = 141, C_SCAN = 10, C_REP = 6 (static, v1) charged at x 1.25 (17.5, 176, 12.5,
7.5) — together 31% of relation collection; a further 1.5x on all four moves the v3 total by
+0.1 bit (see proof.md §11).
