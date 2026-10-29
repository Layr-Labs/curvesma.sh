/* Measurement and self-check harness for side_fused (see side_fused.S).
 * Built and counted exactly as the organizers' scripts/reference_costs pipeline:
 * GCC -O3 -fno-tree-vectorize -mgeneral-regs-only, Valgrind Callgrind Ir,
 * price = (Ir(8192) - Ir(4096)) / 4096, collection toggled around the lane loop.
 * Everything after the counted region is a self-check; exit status non-zero on
 * any mismatch. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"

void side_fused(uint64_t out[15], const uint64_t P[12], const uint64_t T[12],
                const uint64_t inv[6]);

static const uint64_t PM[6] = {0xb9feffffffffaaabULL, 0x1eabfffeb153ffffULL,
                               0x6730d2a0f6b0f624ULL, 0x64774b84f38512bfULL,
                               0x4b1ba7b6434bacd7ULL, 0x1a0111ea397fe69aULL};
#define MAXL 8192
static uint64_t Pp[MAXL][12], Tt[MAXL][12], Inv[MAXL][6], Out[MAXL][15];

static uint64_t rs = 0x9e3779b97f4a7c15ULL;
static uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static void rfp(uint64_t a[6]) { /* uniform-ish residue below p (top limb < 2^60 < p's) */
    for (int i = 0; i < 6; i++) a[i] = rnd();
    a[5] &= 0x0fffffffffffffffULL;
}
static void special(uint64_t a[6], int k) { /* 0, 1, p-1, p-2, random */
    memset(a, 0, 48);
    if (k == 1) a[0] = 1;
    else if (k == 2) { memcpy(a, PM, 48); a[0] -= 1; }
    else if (k == 3) { memcpy(a, PM, 48); a[0] -= 2; }
    else if (k == 4) rfp(a);
}
/* integer helpers on 6 limbs */
static int sub6(uint64_t r[6], const uint64_t a[6], const uint64_t b[6]) {
    unsigned __int128 br = 0;
    for (int i = 0; i < 6; i++) {
        unsigned __int128 t = (unsigned __int128)a[i] - b[i] - br;
        r[i] = (uint64_t)t; br = (t >> 64) & 1;
    }
    return (int)br;
}
static void add6(uint64_t r[6], const uint64_t a[6], const uint64_t b[6]) {
    unsigned __int128 c = 0;
    for (int i = 0; i < 6; i++) { c += (unsigned __int128)a[i] + b[i]; r[i] = (uint64_t)c; c >>= 64; }
}
static int less6(const uint64_t a[6], const uint64_t b[6]) {
    for (int i = 5; i >= 0; i--) if (a[i] != b[i]) return a[i] < b[i];
    return 0;
}
/* reference side output through blst's exported, fully reduced field API */
static void side_ref(uint64_t out[15], const uint64_t P[12], const uint64_t T[12],
                     const uint64_t inv[6]) {
    blst_fp t, mu, mu2, xp, w, yp, n;
    blst_fp_add(&t, (const blst_fp *)(P + 6), (const blst_fp *)(T + 6));
    blst_fp_mul(&mu, &t, (const blst_fp *)inv);
    blst_fp_sqr(&mu2, &mu);
    blst_fp_sub(&xp, &mu2, (const blst_fp *)P);
    blst_fp_sub(&xp, &xp, (const blst_fp *)T);
    blst_fp_sub(&w, &xp, (const blst_fp *)P);
    blst_fp_mul(&yp, &mu, &w);
    blst_fp_sub(&yp, &yp, (const blst_fp *)(P + 6));
    memcpy(out, xp.l, 48); memcpy(out + 6, yp.l, 48);
    sub6(n.l, PM, yp.l);
    const uint64_t *k = less6(yp.l, n.l) ? yp.l : n.l;
    memcpy(out + 12, k, 24);
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 4096, fails = 0;
    if (n > MAXL) n = MAXL;
    for (int i = 0; i < n; i++) {
        if (i < 625) { int k = i; special(Pp[i], k % 5); k /= 5; special(Pp[i] + 6, k % 5); k /= 5;
                       special(Tt[i], k % 5); k /= 5; special(Tt[i] + 6, k % 5); }
        else { rfp(Pp[i]); rfp(Pp[i] + 6); rfp(Tt[i]); rfp(Tt[i] + 6); }
        rfp(Inv[i]);
    }
    side_fused(Out[0], Pp[0], Tt[0], Inv[0]);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for (int i = 0; i < n; i++) side_fused(Out[i], Pp[i], Tt[i], Inv[i]);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%lu\n", (unsigned long)Out[n - 1][14]);

    /* check 1: every lane against the exported-API computation */
    for (int i = 0; i < n; i++) {
        uint64_t ref[15];
        side_ref(ref, Pp[i], Tt[i], Inv[i]);
        if (memcmp(ref, Out[i], 120)) { fails++; if (fails < 5) printf("lane %d mismatch\n", i); }
    }
    printf("check1 lanes=%d mismatches=%d\n", n, fails);
    /* check 2: Lemma 2a, mul_mont_384 with a = r + k p (k = 0..3, a < 4p), b < p */
    int f2 = 0;
    for (int it = 0; it < 200000; it++) {
        uint64_t r[6], b[6], a[6]; blst_fp ref, got;
        if (it == 0) { memcpy(r, PM, 48); r[0]--; memcpy(b, r, 48); } else { rfp(r); rfp(b); }
        blst_fp_mul(&ref, (const blst_fp *)r, (const blst_fp *)b);
        memcpy(a, r, 48);
        for (int k = 1; k <= 3; k++) {
            add6(a, a, PM);
            blst_fp_mul(&got, (const blst_fp *)a, (const blst_fp *)b);
            if (memcmp(ref.l, got.l, 48)) f2++;
        }
    }
    printf("check2 lazy-operand products=600000 mismatches=%d\n", f2); fails += f2;
    /* check 3: genuine G1 points: (x', y') == affine(P - T) */
    int f3 = 0;
    for (int it = 0; it < 64; it++) {
        byte sa[32], sb[32];
        for (int j = 0; j < 32; j++) { sa[j] = (byte)rnd(); sb[j] = (byte)rnd(); }
        blst_p1 A, B, D; blst_p1_affine a, b, nb, d;
        blst_p1_mult(&A, blst_p1_generator(), sa, 256);
        blst_p1_mult(&B, blst_p1_generator(), sb, 256);
        blst_p1_to_affine(&a, &A); blst_p1_to_affine(&b, &B);
        nb = b; blst_fp_cneg(&nb.y, &nb.y, 1);
        blst_p1_add_or_double_affine(&D, &A, &nb);
        blst_p1_to_affine(&d, &D);
        uint64_t P[12], T[12], inv[6], out[15]; blst_fp den, iv;
        memcpy(P, a.x.l, 48); memcpy(P + 6, a.y.l, 48);
        memcpy(T, b.x.l, 48); memcpy(T + 6, b.y.l, 48);
        blst_fp_sub(&den, &b.x, &a.x); blst_fp_inverse(&iv, &den); memcpy(inv, iv.l, 48);
        side_fused(out, P, T, inv);
        if (memcmp(out, d.x.l, 48) || memcmp(out + 6, d.y.l, 48)) f3++;
    }
    printf("check3 curve-point pairs=64 mismatches=%d\n", f3); fails += f3;
    printf(fails ? "selfcheck FAILED\n" : "selfcheck ok\n");
    return fails != 0;
}
