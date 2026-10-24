# Kernel measurements behind C_HIT, C_SMALL and C_BIG

The track metric counts executed scalar AArch64 instructions. The dominant constants of the v4
cost model are **retired-instruction counts** of the kernels in this package, measured with
`/usr/bin/time -l` ("instructions retired") on Apple-silicon AArch64, clang -O2, no fast-math,
0-vs-8 repetition deltas over the same enumerated point set. These are C sources: the verifier's
Python sandbox does not compile them; they are evidence to read and re-run on an AArch64 host.
`evidence.py` part [8] re-counts the la_mac inner loop from the shipped assembly and takes a
census of la_big.s.

## C_HIT: measured 161.8, charged 256 as a production envelope (`archd_hit2.c`)

Same walk/enumeration math as v1-v3 (benchmark lattice `archd_lattice_d12.h`, 2,097,108-point
coset, 3,671,024 node visits), production-shaped hit path:

1. Carried ambient state (checked against direct v = sum z_i b_i at every hit: 0 mismatches).
2. Canonicalisation under (a,b) -> (-a,-b) and sigma: (a,b) -> (a,-b), before the id.
3. **Full-width slab id (v3/v4):** 9 fixed coordinates in balanced radix **4097** (digits
   [-2048,2048], 60.0+48.0 = 108.0 bits in two uint64 words). No truncation before hashing or
   comparison. Production bound |v_i| <= 1699.1 at A=2^10.00 (evidence [2c]); radix 4097 covers
   +/-2048. The v2 radix 1497 (A=2^9.55=748) aliased at the operating point; `packcheck`
   reproduces that alias and finds 0 collisions at 4097 over 2^20 random+boundary tuples.
   The hot path does not multiply by the radix: packed words are carried by the walk
   (precomputed `DLT` deltas) and a hit only negates them for canonicalisation.
4. Descriptor lookup: open addressing keyed by the full 128-bit (k0,k1), exact comparison,
   2^22 slots at load 0.46; present-key probes mean 1.42 slots.
5. Routing: node index in [0,2^20) from a multiplicative hash of both id words; a 2^20-entry
   per-node route table maps the node to one of 8 per-destination ring buffers.
6. Apply at the owner: full-key probe, exact-stride cell `off + (v4-lo4) + wa*(v8-lo8) +
   wab*(v10-lo10)`, add logp.

**Correctness:** reference, single-memory and distributed+apply paths leave bit-identical bucket
arrays (FNV digest `d9f629fbb23665e3` over 2x2,097,108 hits).

**Measured retired instructions per enumerated point**
((count at 8 reps - count at 0 reps)/(8 x 2,097,108); five process pairs, spread < 0.15):

| mode | what | instr/point |
|---|---|---|
| 20 | walk, ambient state live (asm barrier), no hit work | 85.0 |
| 23 | full hit, single shared memory | 140.4 |
| 24 | full hit, send to per-destination buffers | 127.9 |
| 24+apply | send + owner-side apply (the measured path) | **161.8** |

**Why 256 is charged (v4, review objection 4).** C_HIT is a per-hit constant: the per-node
walk/control work and the per-hit canonicalisation/hash/compare/route/apply work are fixed code,
independent of how many nodes or ideals the production sieve has. The production node count
(~2^62.97) multiplies the number of hits; the sieve model charges that hit count separately. The
route table and queues are per-node fixed sizes; the bucket is per-node. The benchmark's final
bucket index is masked into 2^25 bytes while production uses the unmasked 64-bit exact cell, so
the benchmark is one instruction *more* expensive on that step; cache/latency effects are outside
the retired-instruction metric. The benchmark's 1.75 nodes/hit is representative of the
12-dimensional Schnorr-Euchner walk on an LLL-reduced basis, and every production ideal lattice
has the same dimension and walk code. To cover lattice- and node-density variation the claim
charges **C_HIT = 256 = 1.58x (161.8 + routing hop)**. Sensitivity: 170 -> 119.984,
200 -> 120.081, 250 -> 120.229, 256 -> 120.247.

Reproduce:

    clang -O2 -o archd_hit2 archd_hit2.c -lm
    ./archd_hit2 packcheck        # v2 alias reproduced, 0 collisions at radix 4097
    ./archd_hit2 stcheck          # 0 mismatches, max |v_i| = 38
    ./archd_hit2 refcheck 2       # ref == single == distributed, bit-identical
    ./archd_hit2 stats            # table load, probe depths, destination histogram
    for m in 20 23 24; do /usr/bin/time -l ./archd_hit2 icount 0 $m; /usr/bin/time -l ./archd_hit2 icount 8 $m; done
    /usr/bin/time -l ./archd_hit2 icount 0 24 1; /usr/bin/time -l ./archd_hit2 icount 8 24 1

## C_SMALL = 35: small-valuation entries (`la_mac.c/.s`, `la_finish.c`, `la_test.c` — unchanged)

Rows stored sign-split as uint8 magnitudes; branch-free 33-instruction inner loop (re-counted
from `la_mac.s` in evidence.py part [8]); ~100 instructions of per-row overhead; `la_test.c`
checks 2,000 random rows against Python big integers: 0 mismatches.

## C_BIG = 147: full-size Schirokauer entries (`la_big.c`, `la_big.s`, `la_big_check.py`)

Each relation row carries 41 Schirokauer-map values in F_r (255-bit, 4 limbs). `la_big.c` does a
column-ordered 4x4 schoolbook product into a lazy 10-limb accumulator with one Barrett finish per
row. `la_big_check.py` verifies all 64 check rows exactly against Python big integers: 64/64 ok.
Measured 146.9 instr/entry; claimed 147.

Reproduce:

    clang -O2 -o la_big la_big.c && python3 la_big_check.py
    /usr/bin/time -l ./la_big bench 4096 48 0 ; /usr/bin/time -l ./la_big bench 4096 48 8

## Memory (reported, not scored)

The bucket arrays keep one byte per bounding-box point of each live 3-dim slice. The bbox/volume
ratio is the **exact** value for the (2,3,4) slice block, not the generic 6/pi lower bound:
(6/pi)*sqrt(cof00*cof11*cof22)/det(M_TT) = **1.951093168** (recomputed from ELL_M in evidence.py
part [3b]). Peak storage is 1.951093*E with E=(2A+1)^12/64 = 2^126.008, i.e. **2^126.973 bytes**.
