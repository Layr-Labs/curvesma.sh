# Kernel measurements

The score metric is *executed scalar AArch64 instructions*, so the two dominant ledger constants
(`C_HIT_*`, `C_ENTRY`) are measured by differential callgrind `Ir` counts on the actual pinned
toolchain, exactly as `tracks/docs/COST_MODEL.md` prices its reference kernels. `kern.c` is shipped.

Reference image and method: `scripts/reference_costs/reproduce.sh` builds
`curvesma-reference-count` from `scripts/reference_costs/Dockerfile` (Debian12, GCC, Valgrind 3.24.0,
`linux/arm64`) and the reference harness runs

```
valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/$kind-$n.cg /work/kern $kind $n
```

with `CALLGRIND_ZERO_STATS; CALLGRIND_TOGGLE_COLLECT;` around the timed region (the same toggles as
`scripts/reference_costs/fp12.c`). `kern.c` is compiled with
`cc -O3 -fno-tree-vectorize -mgeneral-regs-only -fno-builtin`, the `COST_MODEL.md` flags. Two counts,
2048 and 4096 iterations, are differenced so fixed process overhead cancels.

Reproduce (from the repository root, image already built):

```
mkdir -p /tmp/kres
docker run --rm --platform linux/arm64 --network none \
  -v "$PWD/submissions/bls12-381-fixed-g1-invert-g2/stnfs-gs12-repaired:/harness:ro" \
  -v /tmp/kres:/work curvesma-reference-count sh -c '
  cc -O3 -fno-tree-vectorize -mgeneral-regs-only -fno-builtin -o /work/kern /harness/kern.c
  for k in sieve scan smvp; do for n in 2048 4096; do
    valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/$k-$n.cg /work/kern $k $n 2>/dev/null
  done; done'
```

Measured `Collected` totals and the differentials actually used:

| kernel | Ir @2048 | Ir @4096 | \|diff\|/2048 | per unit | charged |
|---|---:|---:|---:|---|---:|
| `sieve` (log-array apply over a stride-`p` progression) | 2,637,057,710 | 5,274,114,544 | 1,287,814 | 214,400 updates/rep → **6.006 instr/update** | 16 |
| `scan` (threshold sweep, cold branch, 2^16 cells/rep) | 940,591,854 | 1,881,182,958 | 459,273 | 65,536 cells/rep → 7.008 instr/cell | 5 |
| `smvp` (sign-split 4-limb lazy entry accumulate, 4096 entries/rep) | 134,829,544 | 269,657,576 | 65,835 | 4096 entries/rep → **16.07 instr/entry** | 48 |

Notes on the pads, all in the safe direction:

* **`C_HIT_SMALL = 16`** — the measured 6.006 is the naked `ldrb; add; strb; add; cmp; b.lt` inner
  loop, which is the correct shape for ideals of norm `<= 2A+1 = 1301` (a candidate array cell per
  hit, stride `<= 1301`, cache-resident lines). The pad covers line start/stop arithmetic, the region
  predicate on both ends, and the record hand-off that moves survivors to the cofactorisation stage.
* **`C_HIT_LARGE = 64`** — for ideals of norm `> 2A+1` a line holds at most one hit, so each hit needs
  a lattice point enumeration (12 coefficient multiplications and 12 additions in the ambient `(a,b)`
  coordinates, `~30-36` instructions) plus the same apply and a survivor record. 64 is roughly 2x the
  uncounted-but-analytic cost.
* **`C_SCAN = 5`** — the measured 7.008 is a deliberately branch-heavy variant (the test's threshold
  passes on ~50% of cells, so the record store is hot). At the production survivor rate (~`2^-22`) the
  branch is cold and the loop is a load, pointer bump, compare and taken-branch; 5 is conservative.
* **`C_ENTRY = 48`** — the measured 16.07 covers index load, sign test and the 4-limb accumulate. The
  pad covers a signed two-accumulator layout (two adds or a add/sub pair per limb), the row-finish
  Barrett reduction amortised over a row, and the loop/back-edge. Compare `submissions/
  bls12-381-fixed-g2-invert-g1` which measured 33 for a heavier 5-limb signed kernel.

**Addressing (objection b).** Both kernels address a flat array; production addressing is
`addr = slab_base + i` with `i in [0, 2^32)` a dense local coordinate. The sieve is run per slab of
`2^32` bytes (4 GiB, one per node), and the global candidate is recovered from `(node, slab, i)` by the
mixed-radix decomposition of the 12 box coordinates: the `2A+1 = 1301` radix is exact
(`1301^12 = 2^124.08`), so the representation is injective with no collisions, and the per-candidate
reconstruction of the 12 coordinates costs 11 divisions -- charged inside the 1.5x pad on
`C_HIT_SMALL`/`C_SCAN` (11 divisions amortised over 6.3 hits is under 2 instructions per hit). The
survivor record stores `(slab, i, side)` in 16 bytes, and the survivor list is bucketed by slab so the
cofactorisation pass is a pure sequential scan. Nothing is masked or recycled through a small table --
unlike the sibling track's 2^16-entry descriptor table, which was its second objection.

**What is not measured.** Cache misses, TLB behaviour, network hops and the distributed bucket routing
are not instructions and are excluded by the metric. The kernels are single-node benchmarks; the
production split over `2^92`-ish slabs does not change the per-hit instruction count of a local walk.
