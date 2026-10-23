/* Non-reference work of one walk step of the Cheon comb walk (work-cheon/RESEARCH.md).
 *
 * Counted per step (the R term of the accounting):
 *   - exponent multiply c*g_j mod r, multiplier index j from the walk key
 *   - 16 signed 16-bit digits of the exponent, one table read per window with
 *     conditional negation (the 15 batched additions are reference units and are
 *     not executed here; the last read stands in for the walk point)
 *   - mu6 canonicalisation: one Fp multiplication by beta, one addition and one
 *     negation, least-of-three selection, parity of y, conditional negation
 *   - exponent update by mu6 (one mod-r multiplication), distinguished-point test
 *
 * The table holds TBL entries per window: an index costs the same instructions at any
 * table size, so the modelled table is small and its values come from the input file.
 *
 * Usage: kernel <inputs> <n> [dump]
 *   inputs: see work-cheon/measure/gen_inputs.py (hex, one value per line)
 * Functions ref_fp_mul and ref_fr_mul are noinline leaves for REF_FUNCS attribution.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef NATIVE_CHECK
#define CALLGRIND_ZERO_STATS
#define CALLGRIND_TOGGLE_COLLECT
#else
#include <valgrind/callgrind.h>
#endif
#ifndef BLST_H
#define BLST_H "/blst/bindings/blst.h"
#endif
#include BLST_H

#define WINDOWS 16
#define TBL 64
#define NMULT 64
#define NMU 6
#define DP_LOW_BITS 90

static blst_p1_affine comb[WINDOWS][TBL];
static blst_fr mult[NMULT], mu[NMU];
static blst_fp beta;
static unsigned long dps;

typedef struct { blst_p1_affine pt; blst_fr c; } walk_t;

__attribute__((noinline)) static void ref_fp_mul(blst_fp *r, const blst_fp *a, const blst_fp *b) {
    blst_fp_mul(r, a, b);
}

__attribute__((noinline)) static void ref_fr_mul(blst_fr *r, const blst_fr *a, const blst_fr *b) {
    blst_fr_mul(r, a, b);
}

static int hexval(int ch) {
    if (ch >= '0' && ch <= '9') return ch - '0';
    if (ch >= 'a' && ch <= 'f') return ch - 'a' + 10;
    if (ch >= 'A' && ch <= 'F') return ch - 'A' + 10;
    return 0;
}

static void hex2bytes(const char *hex, unsigned char *out, size_t n) {
    for (size_t i = 0; i < n; i++) out[i] = (unsigned char)((hexval(hex[2 * i]) << 4) | hexval(hex[2 * i + 1]));
}

static void read_fp(FILE *f, blst_fp *out) {
    char buf[97];
    unsigned char b[48];
    if (fscanf(f, "%96s", buf) != 1) exit(2);
    hex2bytes(buf, b, 48);
    blst_fp_from_bendian(out, b);
}

static void read_fr(FILE *f, blst_fr *out) {
    char buf[65];
    unsigned char b[32];
    blst_scalar s;
    if (fscanf(f, "%64s", buf) != 1) exit(2);
    hex2bytes(buf, b, 32);
    blst_scalar_from_bendian(&s, b);
    blst_fr_from_scalar(out, &s);
}

static int fp_less(const blst_fp *a, const blst_fp *b) {
    for (int i = 5; i >= 0; i--) {
        if (a->l[i] != b->l[i]) return a->l[i] < b->l[i];
    }
    return 0;
}

static void step(walk_t *w) {
    blst_fr c2, c3;
    unsigned j = (unsigned)(w->pt.x.l[1] & (NMULT - 1));
    ref_fr_mul(&c2, &w->c, &mult[j]);

    blst_scalar s;
    blst_scalar_from_fr(&s, &c2);
    unsigned carry = 0;
    blst_p1_affine acc;
    for (int k = 0; k < WINDOWS; k++) {
        unsigned v = s.b[2 * k] | ((unsigned)s.b[2 * k + 1] << 8);
        v += carry;
        carry = v >= 0x8000;
        unsigned neg = carry;
        if (neg) v = 0x10000 - v;
        acc = comb[k][v & (TBL - 1)];
        if (neg) blst_fp_cneg(&acc.y, &acc.y, 1);
    }

    blst_fp bx, nx, best, ynew;
    ref_fp_mul(&bx, &acc.x, &beta);
    blst_fp_add(&nx, &acc.x, &bx);
    blst_fp_cneg(&nx, &nx, 1);
    unsigned sel = 0;
    best = acc.x;
    if (fp_less(&bx, &best)) { best = bx; sel = 1; }
    if (fp_less(&nx, &best)) { best = nx; sel = 2; }

    ynew = acc.y;
    unsigned char yb[48];
    blst_bendian_from_fp(yb, &ynew);
    unsigned sign = yb[47] & 1;
    if (sign) blst_fp_cneg(&ynew, &ynew, 1);

    ref_fr_mul(&c3, &c2, &mu[sel * 2 + sign]);
    w->pt.x = best;
    w->pt.y = ynew;
    w->c = c3;
    if (((best.l[0]) | (best.l[1] & ((1ULL << (DP_LOW_BITS - 64)) - 1))) == 0) dps++;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: kernel <inputs> <n> [dump]\n");
        return 1;
    }
    int n = atoi(argv[2]);
    int dump = argc > 3;
    FILE *f = fopen(argv[1], "r");
    if (!f) return 1;
    read_fp(f, &beta);
    walk_t w;
    read_fp(f, &w.pt.x);
    read_fp(f, &w.pt.y);
    read_fr(f, &w.c);
    for (int i = 0; i < NMULT; i++) read_fr(f, &mult[i]);
    for (int i = 0; i < NMU; i++) read_fr(f, &mu[i]);
    for (int k = 0; k < WINDOWS; k++) {
        for (int i = 0; i < TBL; i++) {
            read_fp(f, &comb[k][i].x);
            read_fp(f, &comb[k][i].y);
        }
    }
    fclose(f);

#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    for (int i = 0; i < n; i++) step(&w);
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif

    fprintf(stderr, "steps=%d dps=%lu\n", n, dps);
    if (dump) {
        unsigned char xb[48], yb[48], cb[32];
        blst_scalar cs;
        blst_bendian_from_fp(xb, &w.pt.x);
        blst_bendian_from_fp(yb, &w.pt.y);
        blst_scalar_from_fr(&cs, &w.c);
        blst_bendian_from_scalar(cb, &cs);
        for (int i = 0; i < 48; i++) printf("%02x", xb[i]);
        printf("\n");
        for (int i = 0; i < 48; i++) printf("%02x", yb[i]);
        printf("\n");
        for (int i = 0; i < 32; i++) printf("%02x", cb[i]);
        printf("\n");
    }
    return 0;
}
