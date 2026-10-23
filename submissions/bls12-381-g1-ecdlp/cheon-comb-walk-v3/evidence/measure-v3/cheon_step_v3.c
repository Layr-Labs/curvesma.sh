/* The Cheon comb walk as specified in proof.md (v3), measured with Callgrind.
 *
 * Production structure:
 *   - comb tables of the two bases P0 (tame) and Q = [e]P0 (wild): 16 windows, each with
 *     2^15 affine entries, entry j of window k is (j+1) * 2^(16k) * base. The tables are built
 *     in the kernel from the base points by real projective additions and batch normalisation
 *     (preprocessing; build cost measured with config "build").
 *   - walk records: 256 round-robin records, even tame, odd wild.
 *   - digits: signed 16-bit windows with carry, magnitude m in [1, 2^15] added as entry m-1
 *     (negated for negative digits), zero digits skipped.
 *   - the accumulation is real: ref_add is a real affine group addition (reference unit).
 *
 * Per walk point (measured, everything except ref_* functions is counted):
 *   key bits 90..95 select the multiplier; c2 = c * g_j; comb evaluation (16 windows, up to
 *   15 ref_add); canonicalisation (beta*x, x + beta*x, negation, least-of-three in Montgomery
 *   order, parity, conditional negation); c3 = c2 * mu6; walk-state store; 128-bit trail
 *   counter, and at CAP a restart (comb evaluation of a start exponent, canonicalisation);
 *   distinguished-point test; a distinguished point is written into a DP bucket (capacity 256).
 *
 * Configurations (argv[3]): prod (default: DP bits 90, CAP 2^100), dp4 (DP bits 4), cap1
 * (CAP 1), build (build the tables only, toggled). Optional fourth argument dump prints the
 * final records. Usage: kernel <inputs> <n> [prod|dp4|cap1|build|dump]
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

typedef unsigned __int128 u128;

#define WINDOWS 16
#define TBLN (1u << 15)              /* entries per window: magnitudes 1..2^15 */
#define NMULT 64
#define NMU 6
#define NSTART 16
#define WALKS 256
#define DP_BUCKETS 1024              /* measurement table; production: 2^20 buckets */
#define BUCKET_SLOTS 256
#define KEY_BYTES 37                 /* high 291 bits of the 381-bit Montgomery key */

typedef struct {
    blst_p1_affine pt;
    blst_fr c;
    u128 since;
    uint64_t fam;
} walk_t;

typedef struct {
    uint8_t key[KEY_BYTES];
    uint8_t c[32];
    uint8_t fam;
    uint8_t pad[2];
} slot_t;                            /* 72 bytes */

typedef struct {
    uint32_t n;
    slot_t slot[BUCKET_SLOTS];
} bucket_t;

static blst_p1_affine combP[WINDOWS][TBLN], combA[WINDOWS][TBLN];
static blst_p1 work[TBLN];
static blst_fp prefix[TBLN];
static blst_fr mult[NMULT], mu[NMU], starts[NSTART];
static blst_fp beta;
static walk_t recs[WALKS];
static bucket_t buckets[DP_BUCKETS];
static unsigned long dps, matches, overflow, restarts;
static unsigned long n_add, n_proj;   /* call counters for the reference functions */
static unsigned dp_bits = 90;
static u128 cap;

__attribute__((noinline, noipa)) static void ref_fp_mul(blst_fp *r, const blst_fp *a, const blst_fp *b) {
    blst_fp_mul(r, a, b);
}

__attribute__((noinline, noipa)) static void ref_fr_mul(blst_fr *r, const blst_fr *a, const blst_fr *b) {
    blst_fr_mul(r, a, b);
}

/* One table addition in projective coordinates (table build; reference unit). */
__attribute__((noinline, noipa)) static void ref_proj_add(blst_p1 *out, const blst_p1 *a,
                                                          const blst_p1_affine *b) {
    n_proj++;
    blst_p1_add_or_double_affine(out, a, b);
}

/* One real batched affine G1 addition (walk accumulation; reference unit). */
__attribute__((noinline, noipa)) static void ref_add(blst_p1_affine *out, const blst_p1_affine *a,
                                                     const blst_p1_affine *b) {
    blst_p1 pa, pr;
    n_add++;
    blst_p1_from_affine(&pa, a);
    blst_p1_add_or_double_affine(&pr, &pa, b);
    blst_p1_to_affine(out, &pr);
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

/* Affine normalisation of n projective points with one inversion (Montgomery's trick). */
static void batch_to_affine(const blst_p1 *in, blst_p1_affine *out, unsigned n) {
    prefix[0] = in[0].z;
    for (unsigned i = 1; i < n; i++) blst_fp_mul(&prefix[i], &prefix[i - 1], &in[i].z);
    blst_fp inv;
    blst_fp_inverse(&inv, &prefix[n - 1]);
    for (unsigned i = n; i-- > 0;) {
        blst_fp zinv;
        if (i > 0) blst_fp_mul(&zinv, &inv, &prefix[i - 1]);
        else zinv = inv;
        blst_fp_mul(&inv, &inv, &in[i].z);
        blst_fp zinv2, zinv3;
        blst_fp_sqr(&zinv2, &zinv);
        blst_fp_mul(&zinv3, &zinv2, &zinv);
        blst_fp_mul(&out[i].x, &in[i].x, &zinv2);
        blst_fp_mul(&out[i].y, &in[i].y, &zinv3);
    }
}

/* tbl[k][j] = (j+1) * 2^(16k) * base, built by projective additions, then normalised. */
static void build_tables(const blst_p1_affine *base, blst_p1_affine (*tbl)[TBLN]) {
    blst_p1 cur;
    blst_p1_from_affine(&cur, base);
    for (int k = 0; k < WINDOWS; k++) {
        blst_p1_affine b;
        blst_p1_to_affine(&b, &cur);
        work[0] = cur;
        for (unsigned j = 1; j < TBLN; j++) ref_proj_add(&work[j], &work[j - 1], &b);
        batch_to_affine(work, tbl[k], TBLN);
        for (int i = 0; i < 16; i++) blst_p1_double(&cur, &cur);
    }
}

/* [c]B for the table: signed 16-bit digits with carry; magnitude m adds entry m-1. */
static void comb_eval(blst_p1_affine *out, const blst_fr *c, blst_p1_affine (*tbl)[TBLN]) {
    blst_scalar s;
    blst_scalar_from_fr(&s, c);
    unsigned carry = 0;
    int has = 0;
    blst_p1_affine acc, e;
    for (int k = 0; k < WINDOWS; k++) {
        unsigned v = s.b[2 * k] | ((unsigned)s.b[2 * k + 1] << 8);
        v += carry;
        carry = v >= 0x8000;
        unsigned neg = carry;
        if (neg) v = 0x10000 - v;
        if (v == 0) continue;
        e = tbl[k][v - 1];
        if (neg) blst_fp_cneg(&e.y, &e.y, 1);
        if (!has) {
            acc = e;
            has = 1;
        } else {
            ref_add(&acc, &acc, &e);
        }
    }
    *out = acc;
}

/* Canonical point of pt; returns the mu6 index sel*2+sign. */
static unsigned canon(const blst_p1_affine *pt, blst_p1_affine *out) {
    blst_fp bx, nx, best, ynew;
    ref_fp_mul(&bx, &pt->x, &beta);
    blst_fp_add(&nx, &pt->x, &bx);
    blst_fp_cneg(&nx, &nx, 1);
    unsigned sel = 0;
    best = pt->x;
    if (fp_less(&bx, &best)) { best = bx; sel = 1; }
    if (fp_less(&nx, &best)) { best = nx; sel = 2; }
    ynew = pt->y;
    unsigned char yb[48];
    blst_bendian_from_fp(yb, &ynew);
    unsigned sign = yb[47] & 1;
    if (sign) blst_fp_cneg(&ynew, &ynew, 1);
    out->x = best;
    out->y = ynew;
    return sel * 2 + sign;
}

static blst_p1_affine (*table_of(uint64_t fam))[TBLN] {
    return fam ? combA : combP;
}

static void dp_store(const blst_p1_affine *pt, const blst_fr *c, uint64_t fam) {
    const blst_fp *l = &pt->x;
    uint64_t w[5];
    w[0] = (l->l[1] >> 26) | (l->l[2] << 38);
    w[1] = (l->l[2] >> 26) | (l->l[3] << 38);
    w[2] = (l->l[3] >> 26) | (l->l[4] << 38);
    w[3] = (l->l[4] >> 26) | (l->l[5] << 38);
    w[4] = l->l[5] >> 26;
    unsigned b = (unsigned)((w[0] * 0x9E3779B97F4A7C15ULL) >> (64 - 10));
    bucket_t *bk = &buckets[b];
    for (uint32_t i = 0; i < bk->n; i++) {
        if (memcmp(bk->slot[i].key, w, KEY_BYTES) == 0) matches++;
    }
    if (bk->n < BUCKET_SLOTS) {
        slot_t *s = &bk->slot[bk->n];
        memcpy(s->key, w, KEY_BYTES);
        blst_scalar cs;
        blst_scalar_from_fr(&cs, c);
        blst_bendian_from_scalar(s->c, &cs);
        s->fam = (uint8_t)fam;
        bk->n++;
    } else {
        overflow++;
    }
    dps++;
}

static void step(walk_t *w) {
    blst_fr c2, c3;
    unsigned j = (unsigned)((w->pt.x.l[1] >> 26) & (NMULT - 1));
    ref_fr_mul(&c2, &w->c, &mult[j]);
    blst_p1_affine pt, cpt;
    comb_eval(&pt, &c2, table_of(w->fam));
    unsigned idx = canon(&pt, &cpt);
    ref_fr_mul(&c3, &c2, &mu[idx]);
    w->pt = cpt;
    w->c = c3;
    w->since += 1;
    if (w->since >= cap) {
        blst_fr cs = starts[restarts % NSTART];
        restarts++;
        comb_eval(&pt, &cs, table_of(w->fam));
        idx = canon(&pt, &cpt);
        ref_fr_mul(&c3, &cs, &mu[idx]);
        w->pt = cpt;
        w->c = c3;
        w->since = 0;
    }
    const blst_fp *l = &w->pt.x;
    int dp;
    if (dp_bits <= 64) dp = (l->l[0] & ((1ULL << dp_bits) - 1)) == 0;
    else dp = (l->l[0] == 0) && ((l->l[1] & ((1ULL << (dp_bits - 64)) - 1)) == 0);
    if (dp) dp_store(&w->pt, &w->c, w->fam);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: kernel <inputs> <n> [prod|dp4|cap1|build|dump]\n");
        return 1;
    }
    int n = atoi(argv[2]);
    const char *cfg = argc > 3 ? argv[3] : "prod";
    int dump = argc > 4 || strcmp(cfg, "dump") == 0;   /* dump: production, print records */
    if (strcmp(cfg, "dp4") == 0) dp_bits = 4;
    cap = strcmp(cfg, "cap1") == 0 ? (u128)1 : ((u128)1 << 100);
    int build_only = strcmp(cfg, "build") == 0;
    FILE *f = fopen(argv[1], "r");
    if (!f) return 1;
    blst_fp x0, y0, x1, y1;
    read_fp(f, &beta);
    read_fp(f, &x0);
    read_fp(f, &y0);
    read_fp(f, &x1);
    read_fp(f, &y1);
    for (int i = 0; i < NMULT; i++) read_fr(f, &mult[i]);
    for (int i = 0; i < NMU; i++) read_fr(f, &mu[i]);
    for (int i = 0; i < NSTART; i++) read_fr(f, &starts[i]);
    fclose(f);

    blst_p1_affine P0 = {x0, y0}, Q = {x1, y1};

#ifndef REF_ONLY
    if (build_only) {
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
    }
#endif
    build_tables(&P0, combP);
    build_tables(&Q, combA);
#ifndef REF_ONLY
    if (build_only) {
        CALLGRIND_TOGGLE_COLLECT;
    }
#endif
    if (build_only) {
        fprintf(stderr, "build done\n");
        return 0;
    }

    /* Initial records: c_w = starts[w & 15] * mult[(w >> 4) & 63], pt = canon([c_w]B). */
    for (int w = 0; w < WALKS; w++) {
        recs[w].fam = (uint64_t)(w & 1);
        blst_fr c;
        blst_fr_mul(&c, &starts[w & (NSTART - 1)], &mult[(w >> 4) & (NMULT - 1)]);
        blst_p1_affine pt, cpt;
        comb_eval(&pt, &c, table_of(recs[w].fam));
        unsigned idx = canon(&pt, &cpt);
        blst_fr c3;
        blst_fr_mul(&c3, &c, &mu[idx]);
        recs[w].pt = cpt;
        recs[w].c = c3;
        recs[w].since = 0;
    }

#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    for (int i = 0; i < n; i++) step(&recs[i & (WALKS - 1)]);
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif

    fprintf(stderr, "steps=%d dps=%lu matches=%lu overflow=%lu restarts=%lu adds=%lu proj=%lu\n", n, dps,
            matches, overflow, restarts, n_add, n_proj);
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
            printf("%016llx%016llx\n", (unsigned long long)(recs[w].since >> 64),
                   (unsigned long long)recs[w].since);
        }
    }
    return 0;
}
