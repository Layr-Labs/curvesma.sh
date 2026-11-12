# Implementation, worst-case instruction measurement and reproduction

## Files

- `fp12.h`, `fp12.c`: scalar 64-bit Montgomery `Fp`, tower `Fp2/Fp6/Fp12`
  (Karatsuba; `xi = 9+u`; `w_t^2 = v`).
- `rho.c`: profile-generator conversion, 32-adding walk, mixer, scalar
  arithmetic mod r, `selftest`, an audit mode `worstcase K j pat` (pins the
  table index and the operand pattern while running the shipped step) and
  `hashdump n` (prints the mixer index of n fixed encoded elements).
- `experiment.py`: standard-library recomputation of the ledger, the capped
  sampler and its failure bound, the mixer cross-check and the production-walk
  evidence in a small-order subgroup of the real BN254 Fp12.
- `LOOP_BOUNDS.txt`: every loop of the algorithm with its hard bound.
- `evidence/worstcase-sweep.txt`: all 704 measured configurations.
- `evidence/production-walk.txt`: frozen 64-trial output of `experiment.py`.
- `evidence/mixer-crosscheck.txt`: C `hashdump` output and the Python match.

## Build and self-test

```
gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o rho fp12.c rho.c
./rho selftest
```

Output on the build machine:

```
g_is_one=0 g^r_is_one=1
assoc=1
```

`g^r = 1` uses the profile generator converted from the profile basis to the
tower by the fixed linear map in `proof.md` section 1. This checks the
conversion and the tower multiplication together.

## Worst-case instruction count (not an average)

Machine: AArch64, 8-core Cortex-X925, Linux 6.8.0. Toolchain: gcc 13.3.0,
Valgrind 3.22.0 (`--tool=callgrind`, event `Ir`).

Why a sweep is a worst-case bound: every loop inside a step has a
data-independent trip count (`fp_mul` 4x4 limbs, the mixer 24 rounds over 48
words, two 4-limb scalar additions), so the executed instruction count depends
on the inputs only through the arms of the data-dependent decisions -- the
carry/reduction arms of `fp_add`, the borrow arm of `fp_sub`, the final
comparison in `fp_mul`, the select in `scalar_add_r`. Arms change *which*
instructions run, never how many iterations run. The sweep below therefore
covers the whole input dependence:

* table index `j` = 0..31: all 32 production table entries;
* operand pattern: 0 = production start `g`, 1 = all limbs `0xffff...`
  (carry arm everywhere), 2 = all limbs 1 (borrow arm everywhere),
  3 = all limbs `FP_ONE_MONT`, 4..19 = deterministic pseudo-random 384-byte
  elements, 20 = all limbs `p_limb - 1` (carry plus the full comparison chain),
  21 = all limbs `p_limb` (every addition overflows the modulus);
* loop lengths K = 1,000 and K = 3,000, differenced so that table build and
  conversion cancel.

Reproduce (the loop that produced `evidence/worstcase-sweep.txt`):

```
for j in $(seq 0 31); do
  for pat in $(seq 0 21); do
    for K in 1000 3000; do
      valgrind --tool=callgrind --callgrind-out-file=cg.$K.$j.$pat.out \
        ./rho worstcase $K $j $pat
    done
  done
done
```

Result over the 704 configurations (1,408 runs):

| statistic | Ir per step |
|---|---:|
| minimum | 46,319.00 |
| median | 47,248.14 |
| **maximum** (j = 13, pat = 12) | **47,264.24** |

**Charged: 47,265 instructions per step** -- the maximum rounded up to the next
instruction, with no percentage padding. The step is exactly the loop body of
`rho_step`: one `fp12_mul`, one 48-word mixer, two 4-word additions mod r and
the table load. Floyd comparison and loop control are outside the measured loop
and are charged separately at 4,096 instructions per iteration, an analytic
over-charge (48 word compares at ~3 instructions each, plus a counter and a
branch).

## Capped sampler cost

`proof.md` section 2 defines the sampler: one SHA-256 call per 256-bit word.
Worst case 66 scalars x K = 64 attempts = 4,224 attempts, charged at 2,048
instructions each = 8,650,752 instruction equivalents, present in the ledger as
its own row. The charge is the cap, not the observed attempt count.

## Memory

32 records of 448 bytes (384-byte element + 32-byte `c_j` + 32-byte `d_j`),
two walk states (about 800 bytes), 48-word hash scratch, the 4,224-word sampler
buffer is not stored (words are consumed as they are produced) and code: below
2^20 bytes. No table outlives the run, no advice is embedded, and the salt,
seed and increments are fresh per run and charged in the ledger.
