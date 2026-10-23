/* One walk step of the Cheon comb walk, as the walk runs it (proof.md section 4).
 *
 * Measured per step (all executed loads, stores and control are counted):
 *   - family selection of the comb table (P0 for tame, A_d for wild)
 *   - multiplier index from the key, exponent multiply c*g_j mod r
 *   - 16 signed 16-bit digits, one table entry copied per window with conditional y negation
 *   - the 16 batched affine additions are reference units (6M each, charged as H): they are
 *     represented by ref_add, a noinline copy that is counted here as an upper bound on the
 *     bookkeeping around each addition
 *   - canonicalisation: beta*x (ref_fp_mul), x + beta*x, negation, least-of-three, parity
 *   - exponent update by the mu6 element (ref_fr_mul)
 *   - walk record store (point, exponent, steps since trail start, family, DP flag)
 *   - distinguished-point test; a distinguished point calls dp_store, which is not executed
 *     in the measured range and is charged analytically (proof.md section 5)
 *   - trail-cap test; a restart is not executed in the measured range and is charged
 *     analytically (proof.md section 6)
 *
 * The loop runs n steps over WALKS records in round-robin order: even records are tame
 * (P0 table), odd records wild (A_d table). Each record is one trail's state.
 *
 * Usage: kernel <inputs> <n> [dump]
 *   dump prints every record's final state.
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
#define WALKS 256                 /* walk records cycled by the loop (half tame, half wild) */
#define RESTART_AT (1ULL << 62)   /* trail cap is 2^100 steps; the measured range never reaches it */

typedef struct {
    blst_p1_affine pt;
    blst_fr c;
    uint64_t since;   /* steps since the trail started */
    uint64_t fam;     /* 0 tame (P0 table), 1 wild (A_d table) */
    uint64_t dp;      /* 1 if the last point was distinguished */
} walk_t;

static blst_p1_affine combP[WINDOWS][TBL], combA[WINDOWS][TBL];
static blst_fr mult[NMULT], mu[NMU];
static blst_fp beta;
static walk_t recs[WALKS];
static unsigned long dps, restarts;

__attribute__((noinline, noipa)) static void ref_fp_mul(blst_fp *r, const blst_fp *a, const blst_fp *b) {
    blst_fp_mul(r, a, b);
}

__attribute__((noinline, noipa)) static void ref_fr_mul(blst_fr *r, const blst_fr *a, const blst_fr *b) {
    blst_fr_mul(r, a, b);
}

/* Stands for one batched affine addition (a reference unit). Here it copies the entry. */
__attribute__((noinline, noipa)) static void ref_add(blst_p1_affine *acc, const blst_p1_affine *e) {
    *acc = *e;
}

/* Not executed in the measured range: a distinguished point has probability 2^-90. */
__attribute__((noinline, noipa)) static void dp_store(const walk_t *w) {
    (void)w;
    dps++;
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
    const blst_p1_affine(*tbl)[TBL] = w->fam ? combA : combP;
    blst_fr c2, c3;
    unsigned j = (unsigned)(w->pt.x.l[1] & (NMULT - 1));
    ref_fr_mul(&c2, &w->c, &mult[j]);

    blst_scalar s;
    blst_scalar_from_fr(&s, &c2);
    unsigned carry = 0;
    blst_p1_affine acc, e;
    for (int k = 0; k < WINDOWS; k++) {
        unsigned v = s.b[2 * k] | ((unsigned)s.b[2 * k + 1] << 8);
        v += carry;
        carry = v >= 0x8000;
        unsigned neg = carry;
        if (neg) v = 0x10000 - v;
        e = tbl[k][v & (TBL - 1)];
        if (neg) blst_fp_cneg(&e.y, &e.y, 1);
        ref_add(&acc, &e);
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
    w->since += 1;
    w->dp = ((best.l[0] | (best.l[1] & ((1ULL << (DP_LOW_BITS - 64)) - 1))) == 0);
    if (w->dp) dp_store(w);
    if (w->since >= RESTART_AT) restarts++;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: kernel <inputs> <n> [dump]\n");
        return 1;
    }
    int n = atoi(argv[2]);
    int dump = argc > 3 && strcmp(argv[3], "dump") == 0;
    FILE *f = fopen(argv[1], "r");
    if (!f) return 1;
    read_fp(f, &beta);
    blst_fp x0, y0;
    blst_fr c0;
    read_fp(f, &x0);
    read_fp(f, &y0);
    read_fr(f, &c0);
    for (int i = 0; i < NMULT; i++) read_fr(f, &mult[i]);
    for (int i = 0; i < NMU; i++) read_fr(f, &mu[i]);
    for (int k = 0; k < WINDOWS; k++)
        for (int i = 0; i < TBL; i++) {
            read_fp(f, &combP[k][i].x);
            read_fp(f, &combP[k][i].y);
        }
    for (int k = 0; k < WINDOWS; k++)
        for (int i = 0; i < TBL; i++) {
            read_fp(f, &combA[k][i].x);
            read_fp(f, &combA[k][i].y);
        }
    fclose(f);
    for (int w = 0; w < WALKS; w++) {
        recs[w].pt.x = x0;
        recs[w].pt.y = y0;
        recs[w].c = c0;
        recs[w].since = 0;
        recs[w].fam = (uint64_t)(w & 1);
        recs[w].dp = 0;
    }

#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    for (int i = 0; i < n; i++) step(&recs[i & (WALKS - 1)]);
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif

    fprintf(stderr, "steps=%d dps=%lu restarts=%lu\n", n, dps, restarts);
    if (dump) {
        for (int w = 0; w < WALKS; w++) {
            unsigned char xb[48], yb[48], cb[32];
            blst_scalar cs;
            blst_bendian_from_fp(xb, &recs[w].pt.x);
            blst_bendian_from_fp(yb, &recs[w].pt.y);
            blst_scalar_from_fr(&cs, &recs[w].c);
            blst_bendian_from_scalar(cb, &cs);
            for (int i = 0; i < 48; i++) printf("%02x", xb[i]);
            printf("\n");
            for (int i = 0; i < 48; i++) printf("%02x", yb[i]);
            printf("\n");
            for (int i = 0; i < 32; i++) printf("%02x", cb[i]);
            printf("\n");
            printf("%016llx\n", (unsigned long long)recs[w].since);
        }
    }
    return 0;
}
