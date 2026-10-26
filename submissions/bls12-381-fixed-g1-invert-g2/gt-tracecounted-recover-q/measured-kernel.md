# Squaring price: analytic itemisation, compiled kernel and Callgrind record

This note derives the squaring price of `proof.md` section 5 and records the measurement that
checks it. The priced object is the straight-line wide-lazy Granger-Scott cyclotomic squaring used
by the walk: six-limb Montgomery representation, blst f7a88067 exported arm64 primitives, C glue
compiled with `cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only`. The price is
the kernel's fixed operation counts times the static instruction length of each pinned primitive
body, plus the static length of the compiled glue body and the loop. Every primitive and the glue
are straight-line (no data-dependent branch), so a static length is also the executed count.
Callgrind on the same compiled kernel gives the same figure; its sources, commands and counts are
reproduced below so that every number can be checked without leaving the package.

## Kernel operation counts

The kernel is three wide Fp4 squarings followed by the final combination. Each Fp4 squaring
contains three wide Fp2 squarings, one narrow Fp2 addition `s = a0 + a1`, eight wide 768-bit Fp
add/subs (three additions and five subtractions on 12-limb values) and four Montgomery reductions.

| Item | Count |
| --- | ---: |
| wide Fp2 squarings (`sqr_382x`) | 9 |
| wide 768-bit Fp additions (`add_mod_384x384`), 3 per Fp4 squaring | 9 |
| wide 768-bit Fp subtractions (`sub_mod_384x384`), 5 per Fp4 squaring | 15 |
| Montgomery reductions (`redc_mont_384`) | 12 |
| narrow Fp2 additions in the Fp4 bodies (`add_mod_384x`) | 3 |
| narrow Fp2 additions in the final chains (`add_mod_384x`) | 15 |
| narrow Fp2 subtractions in the final chains (`sub_mod_384x`) | 3 |
| `xi` multiplications (`mul_by_1_plus_i_mod_384x`) | 1 |

These are the 67 call sites of the compiled body: 9 + 9 + 15 + 12 + 18 + 3 + 1 = 67 `bl`
instructions, one per relocation `R_AARCH64_CALL26` of `kernel.o` (targets: 9 `sqr_382x`, 9
`add_mod_384x384`, 15 `sub_mod_384x384`, 12 `redc_mont_384`, 18 `add_mod_384x`, 3 `sub_mod_384x`,
1 `mul_by_1_plus_i_mod_384x`). The final combination has six outputs of the form `3T -/+ 2a`;
together they contribute the 15 narrow additions and 3 narrow subtractions above.

## Static instruction lengths

`llvm-objdump -d --no-show-raw-insn assembly.o` of the pinned `build/assembly.S` (assembled by the
`run.sh` below); one line per encoded instruction from the symbol to its last `ret`, the alignment
`nop`s after the `ret` excluded. Each length includes the body's prologue, epilogue, call
instructions and return. A local body that falls through into the next (`__add_mod_384` loads
its operands and falls into `__add_mod_384_ab_are_loaded`) is counted with the body it falls into.

| Primitive | Symbol | Instructions |
| --- | --- | ---: |
| 6x6 -> 12-limb multiply | `__mul_384` | 162 |
| wide Fp2 squaring body (two calls of `__mul_384`) | `sqr_382x` | 86 |
| wide 768-bit Fp add wrapper | `add_mod_384x384` | 17 |
| wide 768-bit Fp add body | `__add_mod_384x384` | 45 |
| wide 768-bit Fp sub wrapper | `sub_mod_384x384` | 17 |
| wide 768-bit Fp sub body | `__sub_mod_384x384` | 44 |
| Montgomery reduction wrapper | `redc_mont_384` | 23 |
| Montgomery reduction body | `__mul_by_1_mont_384` | 154 |
| Montgomery reduction tail | `__redc_tail_mont_384` | 27 |
| narrow Fp2 add wrapper (two calls of `__add_mod_384`) | `add_mod_384x` | 24 |
| narrow Fp2 sub wrapper (two calls of `__sub_mod_384`) | `sub_mod_384x` | 24 |
| narrow six-limb add, operands loaded (6 `ldp` + 21) | `__add_mod_384` | 27 |
| narrow six-limb add, operands in registers | `__add_mod_384_ab_are_loaded` | 21 |
| narrow six-limb sub | `__sub_mod_384` | 26 |
| `xi` wrapper (calls `__sub_mod_384` and `__add_mod_384_ab_are_loaded`) | `mul_by_1_plus_i_mod_384x` | 29 |

## Compiled glue body

`llvm-objdump -d --no-show-raw-insn kernel.o` of the compiled `kernel.c` below (built by `run.sh`
in the image below; `kernel.o` sha256 `6b1e2b15db89a6aa68d6fa0761b5acff77ff68563e8ac2815f9d27077bacdf8b`)
shows one symbol, `sqmeasure_kernel`, straight-line, 379 instructions: 236 `mov` and 34 `add`
(270 argument, pointer and address setup), 67 `bl`, 17 `ldr` and 11 `str` (28), 6 `stp` and 6
`ldp` (12), one `sub` and one `ret`. It has no branch other than the 67 calls and the return.

## Itemisation

Each line is a count times a static length from the tables above.

| Part | Instructions |
| --- | ---: |
| 9 wide Fp2 squarings: `9 x (86 + 2 x 162)` | 3,690 |
| 24 wide 768-bit Fp add/sub, 8 per Fp4 squaring: `9 x (17 + 45) + 15 x (17 + 44)` | 1,473 |
| 12 Montgomery reductions: `12 x (23 + 154 + 27)` | 2,448 |
| 18 narrow Fp2 adds: `18 x (24 + 2 x 27)` | 1,404 |
| 3 narrow Fp2 subs: `3 x (24 + 2 x 26)` | 228 |
| one `xi` multiply: `29 + 26 + 21` | 76 |
| `sqmeasure_kernel` body, straight-line: 270 `mov`/`add` setup, 67 `bl`, 28 `ldr`/`str`, 12 `ldp`/`stp`, `sub`, `ret` | 379 |
| loop control | 10 |
| **total** | **9,708** |

The six primitive lines sum to 9,319; with the 379-instruction body and the loop the total is
9,708. Callgrind measures 79,527,954 `Ir` over 8,192 calls of the same kernel (below):
8,192 x 9,708 = 79,527,936, so the measurement is the itemised total plus 18 instructions of
one-time harness entry, 9,708.0022 per call. The package charges **9,709 instructions**, the
ceiling of the larger of the two figures; the charge is at or above both.

## Library and harness

- blst at commit `f7a88067ed98385ebc6cf41ccecf21ead706a5b7` (supranational/blst, `f7a8806`).
- Reference harness sources, reproduced in full below: `cycsqr.c` (cyclotomic squaring),
  `bls_field.c` (`blst_fp_mul`), `fp12sqr.c` (`blst_fp12_sqr`) and `fp12.c` (dense
  `blst_fp12_mul`). Each binary makes one dense multiplication before `CALLGRIND_ZERO_STATS`, so
  the measured accumulator is a GT element; only the loop between `CALLGRIND_TOGGLE_COLLECT`
  calls is counted.
- Compiled-kernel harness sources, reproduced in full below: `sqmeasure.h`, `kernel.c` (the
  straight-line wide-lazy kernel written over blst's exported arm64 primitives), `check.c` (the
  bit-exactness gate against `blst_fp12_cyclotomic_sqr`), `count.c` (the Callgrind loop),
  `run.sh` (build and run) and `perpart.py` (the self-cost split).
- Build: `src/server.c` and `build/assembly.S` of the pinned tree compiled with
  `cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only` for `linux/arm64`, inside
  the arm64 `python:3.12-slim` image whose digest is in the Dockerfile below (gcc, libc6-dev,
  valgrind, perl, make). The compiled kernel is built with the same flags and linked against the
  same two objects.
- Count: each binary runs `n = 8192` calls of the op after `CALLGRIND_ZERO_STATS`, under
  `valgrind --tool=callgrind --collect-atstart=no`. The per-op figure is the Callgrind `summary:`
  (event `Ir`) divided by 8,192.

## Container image

```dockerfile
FROM python:3.12-slim@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1
RUN apt-get update && apt-get install -y --no-install-recommends gcc libc6-dev valgrind perl make && rm -rf /var/lib/apt/lists/*
WORKDIR /work
```

## Reference harness sources

### `cycsqr.c` (cyclotomic squaring)

```c
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):4096;
    blst_fp12 a,b;
    blst_p1_affine p; blst_p2_affine q;
    blst_p1_to_affine(&p,blst_p1_generator());
    blst_p2_to_affine(&q,blst_p2_generator());
    blst_miller_loop(&a,&q,&p);
    blst_final_exp(&a,&a);
    b=a;
    blst_fp12_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) blst_fp12_cyclotomic_sqr(&a,&a);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%d\n",blst_fp12_in_group(&a));
    return 0;
}
```

### `bls_field.c` (`blst_fp_mul`)

```c
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):4096;
    blst_fp a,b;
    uint64_t av[6]={17},bv[6]={42};
    blst_fp_from_uint64(&a,av);
    blst_fp_from_uint64(&b,bv);
    blst_fp_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) blst_fp_mul(&a,&a,&b);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%lu\n",(unsigned long)a.l[0]);
    return 0;
}
```

### `fp12sqr.c` (`blst_fp12_sqr`)

```c
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):4096;
    blst_fp12 a,b;
    blst_p1_affine p; blst_p2_affine q;
    blst_p1_to_affine(&p,blst_p1_generator());
    blst_p2_to_affine(&q,blst_p2_generator());
    blst_miller_loop(&a,&q,&p);
    blst_final_exp(&a,&a);
    b=a;
    blst_fp12_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) blst_fp12_sqr(&a,&a);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%d\n",blst_fp12_in_group(&a));
    return 0;
}
```

### `fp12.c` (dense `blst_fp12_mul`)

```c
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):4096;
    blst_fp12 a,b;
    blst_p1_affine p; blst_p2_affine q;
    blst_p1_to_affine(&p,blst_p1_generator());
    blst_p2_to_affine(&q,blst_p2_generator());
    blst_miller_loop(&a,&q,&p);
    blst_final_exp(&a,&a);
    b=a;
    blst_fp12_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) blst_fp12_mul(&a,&a,&b);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%d\n",blst_fp12_in_group(&a));
    return 0;
}
```

## Compiled-kernel harness sources

### `sqmeasure.h`

```c
#ifndef SQMEASURE_H
#define SQMEASURE_H

#include <stdint.h>

typedef uint64_t sm_limb;
typedef sm_limb sm_vec384[6];
typedef sm_vec384 sm_vec384x[2];
typedef sm_vec384x sm_vec384fp12[6];

/*
 * Granger-Scott cyclotomic squaring written as one straight-line function
 * over the blst f7a88067 arm64 assembly primitives.  ret and a may alias.
 */
void sqmeasure_kernel(sm_vec384fp12 ret, const sm_vec384fp12 a,
                      const sm_vec384 p, uint64_t p0);

#endif
```

### `kernel.c`

```c
/*
 * sqmeasure: the sqkernel straight-line wide-lazy Granger-Scott squaring,
 * compiled over blst f7a88067's own arm64 assembly primitives.
 *
 * Layout is blst's vec384fp12: six Fp2 coefficients c0..c5 = a[0..5],
 * each Fp2 = re + im*u.  (blst's a[0][0..2], a[1][0..2] are the same bytes.)
 *
 * The wide-lazy structure is the one in blst/src/fp12_tower.c
 * (__FP2x2__ sqr_fp4 / cyclotomic_sqr_fp12): three Fp4 squarings, each
 *    t0 = a0^2, t1 = a1^2, s = a0 + a1,
 *    out0 = redc(t0 + (u+1)t1), out1 = redc(s^2 - t0 - t1),
 * followed by the final 3T -/+ 2a chains.
 *
 * Every Fp-level operation is a call to the primitive blst exports; the
 * kernel itself is straight-line C with no loops.
 */

#include "sqmeasure.h"

typedef sm_limb vec384[6];
typedef vec384 vec384x[2];
typedef vec384x vec384fp4[2];
typedef vec384x vec384fp12[6];
typedef sm_limb vec768[12];
typedef vec768 vec768x[2];

extern void sqr_382x(vec768 ret[2], const vec384x a, const vec384 p);
extern void add_mod_384(vec384 ret, const vec384 a, const vec384 b,
                        const vec384 p);
extern void sub_mod_384(vec384 ret, const vec384 a, const vec384 b,
                        const vec384 p);
extern void add_mod_384x(vec384x ret, const vec384x a, const vec384x b,
                         const vec384 p);
extern void sub_mod_384x(vec384x ret, const vec384x a, const vec384x b,
                         const vec384 p);
extern void add_mod_384x384(vec768 ret, const vec768 a, const vec768 b,
                            const vec384 p);
extern void sub_mod_384x384(vec768 ret, const vec768 a, const vec768 b,
                            const vec384 p);
extern void redc_mont_384(vec384 ret, const vec768 a, const vec384 p,
                          uint64_t n0);
extern void mul_by_1_plus_i_mod_384x(vec384x ret, const vec384x a,
                                     const vec384 p);

/*
 * Narrow Fp2 add/sub: the reference fp12_tower.c add_fp2/sub_fp2 call the
 * exported add_mod_384x/sub_mod_384x.  -DSQMEASURE_LIMB builds the stricter
 * 6-limb variant of the ticket's primitive list (two add_mod_384 calls).
 */
#ifdef SQMEASURE_LIMB
static inline void fp2_add(vec384x ret, const vec384x a, const vec384x b,
                           const vec384 p)
{
    add_mod_384(ret[0], a[0], b[0], p);
    add_mod_384(ret[1], a[1], b[1], p);
}

static inline void fp2_sub(vec384x ret, const vec384x a, const vec384x b,
                           const vec384 p)
{
    sub_mod_384(ret[0], a[0], b[0], p);
    sub_mod_384(ret[1], a[1], b[1], p);
}
#else
static inline void fp2_add(vec384x ret, const vec384x a, const vec384x b,
                           const vec384 p)
{
    add_mod_384x(ret, a, b, p);
}

static inline void fp2_sub(vec384x ret, const vec384x a, const vec384x b,
                           const vec384 p)
{
    sub_mod_384x(ret, a, b, p);
}
#endif

/* One wide Fp4 squaring: 3 wide Fp2 squarings, 3 ADDW + 5 SUBW, 4 REDC. */
static inline void sqr_fp4_wide(vec384fp4 ret, const vec384x a0,
                                const vec384x a1, const vec384 p,
                                const uint64_t p0)
{
    vec768x t0, t1, t2;
    vec384x s;

    sqr_382x(t0, a0, p);                    /* a0^2                     */
    sqr_382x(t1, a1, p);                    /* a1^2                     */
    fp2_add(s, a0, a1, p);                  /* s = a0 + a1 (narrow)     */

    sub_mod_384x384(t2[0], t1[0], t1[1], p);/* t2 = (u+1)*t1            */
    add_mod_384x384(t2[1], t1[0], t1[1], p);
    add_mod_384x384(t2[0], t2[0], t0[0], p);/* t2 += t0                 */
    add_mod_384x384(t2[1], t2[1], t0[1], p);
    redc_mont_384(ret[0][0], t2[0], p, p0);
    redc_mont_384(ret[0][1], t2[1], p, p0);

    sqr_382x(t2, s, p);                     /* s^2                      */
    sub_mod_384x384(t2[0], t2[0], t0[0], p);
    sub_mod_384x384(t2[1], t2[1], t0[1], p);
    sub_mod_384x384(t2[0], t2[0], t1[0], p);
    sub_mod_384x384(t2[1], t2[1], t1[1], p);
    redc_mont_384(ret[1][0], t2[0], p, p0);
    redc_mont_384(ret[1][1], t2[1], p, p0);
}

/* 3*t - 2*c mod p, the negative half of the final combination. */
static inline void chain_neg(vec384x ret, const vec384x t, const vec384x c,
                             const vec384 p)
{
    fp2_sub(ret, t, c, p);
    fp2_add(ret, ret, ret, p);
    fp2_add(ret, ret, t, p);
}

/* 3*t + 2*c mod p, the positive half of the final combination. */
static inline void chain_pos(vec384x ret, const vec384x t, const vec384x c,
                             const vec384 p)
{
    fp2_add(ret, t, c, p);
    fp2_add(ret, ret, ret, p);
    fp2_add(ret, ret, t, p);
}

void sqmeasure_kernel(vec384fp12 ret, const vec384fp12 a, const vec384 p,
                      const uint64_t p0)
{
    vec384fp4 t0, t1, t2;
    vec384x xi;

    sqr_fp4_wide(t0, a[0], a[4], p, p0);
    sqr_fp4_wide(t1, a[3], a[2], p, p0);
    sqr_fp4_wide(t2, a[1], a[5], p, p0);

    /* c0..c2 = 3T - 2a */
    chain_neg(ret[0], t0[0], a[0], p);
    chain_neg(ret[1], t1[0], a[1], p);
    chain_neg(ret[2], t2[0], a[2], p);

    /* c3 = 3*(u+1)*t2[1] + 2*a[3] */
    mul_by_1_plus_i_mod_384x(xi, t2[1], p);
    chain_pos(ret[3], xi, a[3], p);

    /* c4..c5 = 3T + 2a */
    chain_pos(ret[4], t0[1], a[4], p);
    chain_pos(ret[5], t1[1], a[5], p);
}
```

### `check.c`

```c
/*
 * Bit-exactness gate: 1000 random cyclotomic-subgroup elements (final-exp
 * outputs and their powers), kernel output vs blst_fp12_cyclotomic_sqr,
 * byte for byte.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "/blst/bindings/blst.h"
#include "sqmeasure.h"

static const sm_vec384 P = {
    0xb9feffffffffaaabULL, 0x1eabfffeb153ffffULL, 0x6730d2a0f6b0f624ULL,
    0x64774b84f38512bfULL, 0x4b1ba7b6434bacd7ULL, 0x1a0111ea397fe69aULL
};

static uint64_t p0_inv(void)
{
    uint64_t x = 1;
    for (int i = 0; i < 6; i++)
        x *= 2 - P[0] * x;
    return -x;
}

static uint64_t rng_state = 0x9e3779b97f4a7c15ULL;

static uint64_t xs(void)
{
    uint64_t x = rng_state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng_state = x;
    return x;
}

int main(int argc, char **argv)
{
    int n = argc > 1 ? atoi(argv[1]) : 1000;
    uint64_t p0 = p0_inv();
    blst_p1_affine p;
    blst_p2_affine q;
    blst_fp12 g, x, y, z;
    int bad = 0;

    blst_p1_to_affine(&p, blst_p1_generator());
    blst_p2_to_affine(&q, blst_p2_generator());
    blst_miller_loop(&g, &q, &p);
    blst_final_exp(&g, &g);                 /* g in GT, cyclotomic subgroup */

    for (int i = 0; i < n; i++) {
        uint64_t e = xs() | 0x8000000000000000ULL;  /* 64-bit, nonzero */
        x = g;
        for (int b = 62; b >= 0; b--) {
            blst_fp12_cyclotomic_sqr(&x, &x);
            if ((e >> b) & 1)
                blst_fp12_mul(&x, &x, &g);
        }
        blst_fp12_cyclotomic_sqr(&y, &x);
        memset(&z, 0, sizeof(z));
        sqmeasure_kernel((sm_vec384x *)&z, (const sm_vec384x *)&x, P, p0);
        if (memcmp(&y, &z, sizeof(y)) != 0) {
            bad++;
            if (bad <= 3)
                printf("mismatch at %d exponent %016llx\n", i,
                       (unsigned long long)e);
        }
    }
    printf("checked %d random cyclotomic elements, mismatches %d, "
           "in_group %d\n", n, bad, blst_fp12_in_group(&z));
    return bad == 0 ? 0 : 1;
}
```

### `count.c`

```c
/*
 * Callgrind harness: the same loop shape as the reference cycsqr.c above.
 * "ctl" calls blst_fp12_cyclotomic_sqr n times, "kernel" calls
 * sqmeasure_kernel n times, both in place from the same GT element.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
#include "sqmeasure.h"

static const sm_vec384 P = {
    0xb9feffffffffaaabULL, 0x1eabfffeb153ffffULL, 0x6730d2a0f6b0f624ULL,
    0x64774b84f38512bfULL, 0x4b1ba7b6434bacd7ULL, 0x1a0111ea397fe69aULL
};

static uint64_t p0_inv(void)
{
    uint64_t x = 1;
    for (int i = 0; i < 6; i++)
        x *= 2 - P[0] * x;
    return -x;
}

int main(int argc, char **argv)
{
    int n = argc > 1 ? atoi(argv[1]) : 8192;
    const char *which = argc > 2 ? argv[2] : "kernel";
    uint64_t p0 = p0_inv();
    blst_fp12 a, b;
    blst_p1_affine p;
    blst_p2_affine q;

    blst_p1_to_affine(&p, blst_p1_generator());
    blst_p2_to_affine(&q, blst_p2_generator());
    blst_miller_loop(&a, &q, &p);
    blst_final_exp(&a, &a);
    b = a;
    blst_fp12_mul(&a, &a, &b);

    CALLGRIND_ZERO_STATS;
    if (strcmp(which, "ctl") == 0) {
        CALLGRIND_TOGGLE_COLLECT;
        for (int i = 0; i < n; i++)
            blst_fp12_cyclotomic_sqr(&a, &a);
        CALLGRIND_TOGGLE_COLLECT;
    } else {
        CALLGRIND_TOGGLE_COLLECT;
        for (int i = 0; i < n; i++)
            sqmeasure_kernel((sm_vec384x *)&a, (const sm_vec384x *)&a,
                             P, p0);
        CALLGRIND_TOGGLE_COLLECT;
    }
    printf("%s %d %d\n", which, n, blst_fp12_in_group(&a));
    return 0;
}
```

### `run.sh`

```sh
#!/bin/sh
set -eu
cc --version | head -1
valgrind --version
CFLAGS="-O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only"
cc $CFLAGS -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
cc $CFLAGS -I/blst/bindings -I/harness -c /harness/kernel.c -o /work/kernel.o
cc $CFLAGS -I/blst/bindings -I/harness -DSQMEASURE_LIMB -c /harness/kernel.c -o /work/kernel_limb.o
cc $CFLAGS -I/blst/bindings -I/harness /harness/check.c /work/kernel.o /work/blst.o /work/assembly.o -o /work/check
cc $CFLAGS -I/blst/bindings -I/harness /harness/check.c /work/kernel_limb.o /work/blst.o /work/assembly.o -o /work/check_limb
cc $CFLAGS -I/blst/bindings -I/harness /harness/count.c /work/kernel.o /work/blst.o /work/assembly.o -o /work/count
cc $CFLAGS -I/blst/bindings -I/harness /harness/count.c /work/kernel_limb.o /work/blst.o /work/assembly.o -o /work/count_limb

/work/check 1000
/work/check_limb 1000

for kind in ctl kernel kernel_limb; do
  case $kind in
    ctl)    bin=/work/count;      arg=ctl;;
    kernel) bin=/work/count;      arg=kernel;;
    *)      bin=/work/count_limb; arg=kernel;;
  esac
  valgrind --tool=callgrind --collect-atstart=no \
    --callgrind-out-file=/work/$kind.callgrind $bin 8192 $arg \
    > /work/$kind.stdout 2> /work/$kind.stderr
  echo "$kind $(grep -E '^summary:' /work/$kind.callgrind | head -1)"
  cat /work/$kind.stdout
done
```

### `perpart.py`

```python
#!/usr/bin/env python3
"""Print the per-part Callgrind self-cost split of out/<name>.callgrind."""
import collections
import re
import sys

BUCKETS = {
    "wide Fp2 squarings": ["__mul_384", "sqr_382x"],
    "wide add/sub": ["add_mod_384x384", "__add_mod_384x384",
                     "sub_mod_384x384", "__sub_mod_384x384"],
    "Montgomery reductions": ["redc_mont_384", "__mul_by_1_mont_384",
                              "__redc_tail_mont_384"],
    "narrow adds/chains/xi": ["add_mod_384x", "sub_mod_384x",
                              "add_mod_384", "sub_mod_384",
                              "__add_mod_384", "__sub_mod_384",
                              "mul_by_1_plus_i_mod_384x"],
    "kernel body": ["sqmeasure_kernel'2"],
}


def self_costs(path):
    names = {}
    cost = collections.Counter()
    cur = None
    pending_call = False
    for line in open(path):
        line = line.rstrip("\n")
        m = re.match(r"^(fn|cfn)=\((\d+)\)(?:\s+(.*))?$", line)
        if m:
            kind, i, name = m.groups()
            if name is not None:
                names[i] = name
            if kind == "fn":
                cur = i
                pending_call = False
            else:
                pending_call = True
            continue
        if line.startswith("calls="):
            pending_call = True
            continue
        m = re.match(r"^[0-9+\-*]+\s+([0-9]+)$", line)
        if m:
            if pending_call:
                pending_call = False
            elif cur is not None:
                cost[names.get(cur, cur)] += int(m.group(1))
    return cost


def main():
    for path in sys.argv[1:]:
        cost = self_costs(path)
        total = sum(cost.values())
        print("== %s  total %d  per call %.4f" % (path, total, total / 8192))
        seen = set()
        for title, names in BUCKETS.items():
            v = sum(cost[n] for n in names)
            seen.update(names)
            print("  %-24s %10.2f/call" % (title, v / 8192))
        loop = cost["main'2"] / 8192
        seen.add("main'2")
        print("  %-24s %10.2f/call" % ("loop/main", loop))
        rest = total / 8192 - sum(cost[n] for n in seen) / 8192
        if abs(rest) > 0.01:
            print("  %-24s %10.2f/call" % ("other", rest))


if __name__ == "__main__":
    main()
```

## Build and Callgrind commands

In the commands below, `/harness` is the directory holding the sources reproduced above, `/blst`
is the blst tree at the commit named above, and `/work` is the container's output directory set by
the Dockerfile above. None of these is a path in a solver workspace: every source needed to check
a quoted count is either reproduced in this file or fixed by the named commit.

The reference image above is built for `linux/arm64` and run with the harness and the pinned blst
tree mounted read-only, `/harness/run2.sh` as its command (job.sh):

```sh
#!/bin/bash
set -euo pipefail
RC=$(cd "$(dirname "$0")" && pwd)
docker build -q --platform linux/arm64 -t curvesma-reference-count "$RC/harness"
docker run --rm --platform linux/arm64 --network none --name refcost -v "$RC/blst:/blst:ro" -v "$RC/harness:/harness:ro" -v "$RC/results:/work" curvesma-reference-count sh /harness/run2.sh
```

Inside the container, run2.sh compiles the pinned blst objects once, then builds and counts each
reference binary; the 8,192-count rows are the ones quoted (run2.sh):

```sh
#!/bin/sh
set -eu
cc --version | head -1
valgrind --version
cc -O3 -fno-builtin -fPIC -fno-tree-vectorize -mgeneral-regs-only -I/blst/bindings -c /blst/src/server.c -o /work/blst.o
cc -c /blst/build/assembly.S -o /work/assembly.o
for kind in fp12 bls_field cycsqr fp12sqr; do
  cc -O3 -fno-tree-vectorize -mgeneral-regs-only /harness/$kind.c /work/blst.o /work/assembly.o -o /work/$kind
  for n in 4096 8192; do
    valgrind --tool=callgrind --collect-atstart=no --callgrind-out-file=/work/$kind-$n.callgrind /work/$kind "$n" > /work/$kind-$n.stdout 2> /work/$kind-$n.stderr
    echo "$kind $n $(grep -E '^summary:|^totals:' /work/$kind-$n.callgrind | head -1)"
  done
done
```

The compiled kernel is built and run by `run.sh` in the same image, with the harness at `/harness`
and the objects and output at `/work`:

```sh
docker build -q --platform linux/arm64 -t curvesma-sqmeasure .
docker run --rm --platform linux/arm64 --network none --user $(id -u):$(id -g) \
  -v "<the pinned blst checkout>:/blst:ro" -v "<the harness directory>:/harness:ro" \
  -v "<an output directory>:/work" curvesma-sqmeasure sh /harness/run.sh
python3 perpart.py out/ctl.callgrind out/kernel.callgrind out/kernel_limb.callgrind
```

## Callgrind counts (grader's arm64 image, n = 8,192)

| operation | total `Ir` | per call |
| --- | ---: | ---: |
| `blst_fp12_mul` (dense, `fp12.c`) | 188,383,252 | 22,996 |
| `blst_fp12_sqr` (`fp12sqr.c`) | 134,004,756 | 16,358 |
| `blst_fp12_cyclotomic_sqr` (control, `cycsqr.c`) | 80,306,196 | 9,803 |
| `blst_fp_mul` (`bls_field.c`) | 3,137,556 | 383 |
| compiled kernel, reference `add_fp2` form (`kernel.c`) | 79,527,954 | **9,708.0022** |
| compiled kernel, 6-limb `add_mod_384` form (`-DSQMEASURE_LIMB`) | 82,771,986 | 10,104.0022 |

The control row is the reference run's 80,306,196 (9,803.0024). In the compiled-kernel run the same
harness reads 80,306,191 = **9,803.0018**, five instructions fewer from main/loop layout; the 95.0
Ir/call difference below is to that same-run figure.

The package uses 383 for one Fp multiplication and charges the squaring branch at **9,709
instructions**, the ceiling of 79,527,954/8,192 = 9,708.002197265625. The fractional remainder is
not discarded; the charge includes the harness loop and fixed collection overhead.
The same run's control, `blst_fp12_cyclotomic_sqr`, is 9,803.0018; the reference kernel is 95.0
Ir/call below it. The 6-limb form (two `add_mod_384` calls per narrow Fp2 op) is recorded, not
used. The protected cost convention fixes one reference dense multiplication at C = 20,000 instruction
equivalents, regardless of this harness's observed 22,996 instructions; this measurement does not
redefine the unit or introduce another dense-multiplication price inside the squaring step.

## Bit-exactness gate

`check.c` builds GT elements as `blst_final_exp(blst_miller_loop(...))` and random nonidentity
powers `g^e` of them (`e` a 64-bit xorshift value with bit 63 set), and compares the kernel's 576
output bytes with those of `blst_fp12_cyclotomic_sqr` with `memcmp`. Both compiled forms report
`checked 1000 random cyclotomic elements, mismatches 0, in_group 1`; no count would have been
reported on a mismatch. The gate's accumulator is built with the library squaring, so the kernel is
tested only on cyclotomic elements, which is where it is used.

## Per-part split of the compiled kernel (`perpart.py` self cost, Ir/call)

Reference `add_fp2` form:

| part | measured |
| --- | ---: |
| 9 wide Fp2 squarings (18 `__mul_384` at 2,916 + 9 `sqr_382x` bodies at 774) | 3,690.00 |
| 24 wide 768-bit Fp add/sub, 8 per Fp4 squaring (`add_mod_384x384` 153 + `__add_mod_384x384` 405 + `sub_mod_384x384` 255 + `__sub_mod_384x384` 660) | 1,473.00 |
| 12 Montgomery reductions (`redc_mont_384` 276 + `__mul_by_1_mont_384` 1,848 + `__redc_tail_mont_384` 324) | 2,448.00 |
| narrow sums, chains and `xi` (`add_mod_384x` 432 + `sub_mod_384x` 72 + `__add_mod_384` 993 + `__sub_mod_384` 182 + `mul_by_1_plus_i_mod_384x` 29) | 1,708.00 |
| C glue body (`sqmeasure_kernel`: argument/address setup, loads/stores, calls) | 378.95 |
| other (body instructions Callgrind attributes outside the symbol, and harness entry) | 0.05 |
| loop | 10.00 |
| **total** | **9,708.0022** |

Against the itemisation above: the narrow line is 1,404 + 228 + 76 = 1,708, where `__add_mod_384`
993 = 36 x 27 + 21 (the `xi` wrapper enters at `__add_mod_384_ab_are_loaded`, which the profile
books under `__add_mod_384`) and `__sub_mod_384` 182 = 7 x 26; the body 378.95 plus the 0.05 other
is the 379 static instructions of `sqmeasure_kernel` and the 18-instruction entry.

The profile is identical to the control's for every shared primitive: 9 `sqr_382x` at 774, 18
`__mul_384` at 2,916, 24 wide 768-bit Fp add/subs, 12 Montgomery reductions, and the narrow chains. The 95.0
Ir/call saving over the control is entirely the C glue: the control pays `sqr_fp4` 336 plus the
`cyclotomic_sqr_fp12` / `blst_fp12_cyclotomic_sqr` wrappers 142 plus the loop 6 = 484, while the
flattened kernel pays 378.95 plus the loop 10 = 389.
