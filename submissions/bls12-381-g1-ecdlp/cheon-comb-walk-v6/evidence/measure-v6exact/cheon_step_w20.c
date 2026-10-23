/* The Cheon comb walk with batched affine additions, window 20 (v5 with w = 20), measured with Callgrind.
 *
 * Reference operation (6M convention, COST_MODEL.md): one batched incomplete affine G1 addition, whose
 * six field multiplications/squarings are the noinline wrappers core_mul and core_sqr. Everything else
 * is other work: subtractions, copies, loads and stores, the batch inversion (one per stage), the comb
 * digits, the canonicalisation (branch-free), the exponent products, the distinguished-point test and
 * the restarts.
 *
 * Comb: 13 signed 20-bit windows with carry; a window of magnitude m >= 1 adds the table entry
 * m * 2^(20k) * B, stored at index m - 1 (2^19 entries per window and base, full table).
 *
 * Modes (argv[3]):
 *   full        full tables (2^19 entries per window and base); used for the native replica check
 *   (default)   lite: tables of LITE_N entries per window; entry indices are masked to LITE_N (measurement
 *               only: the instruction count does not depend on the table length or values, proof.md section 5)
 *   forced      lite + every digit nonzero, fixed sign pattern (maximal point path; measurement only)
 *   forcedfull  full tables + forced digits (cross-check of the lite measurement)
 *   forcedcap1  lite + forced digits + cap 1 (restart path, measurement only)
 *   dp4         lite, distinguished points with 4 bits (DP path)
 *   cap1        lite, cap 1 (restart path)
 *   build       lite table build for n stages (differenced over stage counts)
 * Optional fourth argument "dump" prints every record.
 * Usage: kernel <inputs> <n> [mode] [dump]; n is a multiple of RECS for walk modes.
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

#define WBITS 20
#define WINDOWS 13                   /* ceil(256 / 20) */
#define TBLN_FULL (1u << 19)         /* magnitudes 1 .. 2^19 */
#define LITE_N 4096u
#define HALF (1u << 19)
#define WMASK ((1u << 20) - 1)
#define RECS 4096
#define NMULT 64
#define NMU 6
#define NSTART 16
#define DP_BUCKETS 1024
#define BUCKET_SLOTS 256
#define KEY_BYTES 37

typedef struct {
    blst_p1_affine pt;
    blst_fr c;
    u128 since;
    uint64_t fam;
} rec_t;

typedef struct {
    uint8_t key[KEY_BYTES];
    uint8_t c[32];
    uint8_t fam;
    uint8_t pad[2];
} slot_t;

typedef struct {
    uint32_t n;
    slot_t slot[BUCKET_SLOTS];
} bucket_t;

static blst_p1_affine (*combP)[TBLN_FULL], (*combA)[TBLN_FULL];   /* [WINDOWS][TBLN_FULL], allocated lazily */
static unsigned tbl_n = TBLN_FULL;     /* entries built and used per window */
static unsigned tbl_mask = TBLN_FULL - 1;
static blst_fr mult[NMULT], mu[NMU], starts[NSTART], initc[RECS];
static blst_fp beta;
static rec_t recs[RECS];
static blst_fr c2s[RECS];
static blst_p1_affine ptsbuf[RECS];
static uint64_t fams[RECS];
static uint32_t mag[RECS][WINDOWS];
static uint8_t negs[RECS][WINDOWS];
static blst_p1_affine acc[RECS];
static uint8_t started[RECS];
static blst_p1_affine bP[RECS], bQ[RECS], bres[RECS];
static uint32_t bidx[RECS];
static blst_fp dxs[RECS], pre[RECS];
static bucket_t buckets[DP_BUCKETS];
static unsigned long dps, matches, overflow, restarts, n_adds, n_lemma, n_rounds;
static unsigned long n_mul, n_sqr;   /* exact call counts of the reference core wrappers */
static unsigned dp_bits = 90;
static unsigned forced = 0;
static u128 cap;

/* Reference core: the six field multiplications of one batched affine addition. */
__attribute__((noinline, noipa)) static void core_mul(blst_fp *r, const blst_fp *a, const blst_fp *b) {
    n_mul++;
    blst_fp_mul(r, a, b);
}

__attribute__((noinline, noipa)) static void core_sqr(blst_fp *r, const blst_fp *a) {
    n_sqr++;
    blst_fp_sqr(r, a);
}

static void dp_store(const blst_p1_affine *pt, const blst_fr *c, uint64_t fam);

static unsigned fp_lt_ct(const blst_fp *a, const blst_fp *b) {
    uint64_t borrow = 0;
    for (int i = 0; i < 6; i++) {
        unsigned __int128 t = (unsigned __int128)a->l[i] - b->l[i] - borrow;
        borrow = (uint64_t)(t >> 64) & 1;
    }
    return (unsigned)borrow;
}

static void fp_select(blst_fp *out, unsigned flag, const blst_fp *a, const blst_fp *b) {
    uint64_t m = 0 - (uint64_t)flag;
    for (int i = 0; i < 6; i++) out->l[i] = (a->l[i] & m) | (b->l[i] & ~m);
}

/* Branch-free table select (fam is 0 or 1). */
static blst_p1_affine (*table_of(uint64_t fam))[TBLN_FULL] {
    uintptr_t a = (uintptr_t)combP, b = (uintptr_t)combA;
    uintptr_t m = (uintptr_t)0 - (uintptr_t)(fam & 1u);
    return (blst_p1_affine (*)[TBLN_FULL])(a ^ ((a ^ b) & m));
}

/* One batch of n affine additions bP[j] + bQ[j] -> bres[j], one inversion. */
static void batch_add(unsigned n) {
    for (unsigned j = 0; j < n; j++) {
        blst_fp_sub(&dxs[j], &bQ[j].x, &bP[j].x);
        n_lemma += (unsigned long)(dxs[j].l[0] | dxs[j].l[1] | dxs[j].l[2] | dxs[j].l[3] | dxs[j].l[4] | dxs[j].l[5]) == 0;
    }
    pre[0] = dxs[0];
    for (unsigned j = 1; j < n; j++) core_mul(&pre[j], &pre[j - 1], &dxs[j]);
    blst_fp inv;
    blst_fp_inverse(&inv, &pre[n - 1]);
    for (unsigned j = n; j-- > 0;) {
        blst_fp invj;
        if (j > 0) core_mul(&invj, &inv, &pre[j - 1]);
        else invj = inv;
        core_mul(&inv, &inv, &dxs[j]);
        blst_fp dy, lam, lam2, x3, t, y3;
        blst_fp_sub(&dy, &bQ[j].y, &bP[j].y);
        core_mul(&lam, &dy, &invj);
        core_sqr(&lam2, &lam);
        blst_fp_sub(&x3, &lam2, &bP[j].x);
        blst_fp_sub(&x3, &x3, &bQ[j].x);
        blst_fp_sub(&t, &bP[j].x, &x3);
        core_mul(&y3, &lam, &t);
        blst_fp_sub(&bres[j].y, &y3, &bP[j].y);
        bres[j].x = x3;
    }
    n_adds += n;
}

/* Unbatched affine addition (other work): restarts. */
static void add_unbatched(blst_p1_affine *out, const blst_p1_affine *a, const blst_p1_affine *b) {
    blst_p1 pa, pr;
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

/* Signed window k (bits [20k, 20k+20)) of a scalar held in the padded little-endian buffer buf[40]. */
static void window_digit(const uint8_t buf[40], unsigned k, unsigned *carry, unsigned *m, unsigned *neg) {
    unsigned q = (WBITS * k) >> 3, r = (WBITS * k) & 7u;
    uint32_t w = (uint32_t)buf[q] | ((uint32_t)buf[q + 1] << 8) | ((uint32_t)buf[q + 2] << 16) |
                 ((uint32_t)buf[q + 3] << 24);
    unsigned v = (unsigned)((w >> r) & WMASK);
    v += *carry;
    *carry = v >= HALF;
    *neg = *carry;
    unsigned mk = 0u - *neg;
    *m = ((0x100000u - v) & mk) | (v & ~mk);
}

static void scalar_buf(const blst_fr *c, uint8_t buf[40]) {
    blst_scalar s;
    blst_scalar_from_fr(&s, c);
    memset(buf, 0, 40);
    memcpy(buf, s.b, 32);
}

/* Measurement only: every magnitude nonzero and a fixed sign pattern (maximal point path). */
static void force_digits(unsigned i) {
    for (int k = 0; k < WINDOWS; k++) {
        mag[i][k] = 1 + ((i * 7u + (unsigned)k * 13u) % (HALF - 1u));
        negs[i][k] = (uint8_t)((i + (unsigned)k) & 1u);
    }
}


/* Comb evaluation of count exponents cs[i] with bases fam[i], batched per stage. */
static void comb_batch(const blst_fr *cs, const uint64_t *fam, unsigned count, blst_p1_affine *out) {
    for (unsigned i = 0; i < count; i++) {
        uint8_t buf[40];
        scalar_buf(&cs[i], buf);
        unsigned carry = 0;
        for (int k = 0; k < WINDOWS; k++) {
            unsigned m, neg;
            window_digit(buf, (unsigned)k, &carry, &m, &neg);
            mag[i][k] = m;
            negs[i][k] = (uint8_t)neg;
        }
        if (forced) force_digits(i);
        started[i] = 0;
    }
    for (int k = 0; k < WINDOWS; k++) {
        unsigned n = 0;
        for (unsigned i = 0; i < count; i++) {
            unsigned m = mag[i][k];
            if (m == 0) continue;
            blst_p1_affine e = table_of(fam[i])[k][(m - 1) & tbl_mask];
            blst_fp_cneg(&e.y, &e.y, negs[i][k]);
            if (!started[i]) {
                acc[i] = e;
                started[i] = 1;
            } else {
                bP[n] = acc[i];
                bQ[n] = e;
                bidx[n] = i;
                n++;
            }
        }
        if (n) {
            batch_add(n);
            for (unsigned j = 0; j < n; j++) acc[bidx[j]] = bres[j];
        }
    }
    for (unsigned i = 0; i < count; i++) out[i] = acc[i];
}

/* Unbatched comb evaluation of one exponent (restarts). */
static void comb_single(blst_p1_affine *out, const blst_fr *c, uint64_t fam) {
    uint8_t buf[40];
    scalar_buf(c, buf);
    unsigned carry = 0;
    int has = 0;
    blst_p1_affine a;
    for (int k = 0; k < WINDOWS; k++) {
        unsigned v, neg;
        window_digit(buf, (unsigned)k, &carry, &v, &neg);
        if (forced) {
            v = 1 + (((unsigned)k * 13u) % (HALF - 1u));
            neg = (unsigned)(k & 1);
        }
        if (v == 0) continue;
        blst_p1_affine e = table_of(fam)[k][(v - 1) & tbl_mask];
        blst_fp_cneg(&e.y, &e.y, neg);
        if (!has) {
            a = e;
            has = 1;
        } else {
            add_unbatched(&a, &a, &e);
        }
    }
    *out = a;
}

/* Canonical point; returns the mu6 index sel*2+sign. Branch-free. */
static unsigned canon(const blst_p1_affine *pt, blst_p1_affine *out) {
    blst_fp bx, nx, best, ynew;
    blst_fp_mul(&bx, &pt->x, &beta);
    blst_fp_add(&nx, &pt->x, &bx);
    blst_fp_cneg(&nx, &nx, 1);
    unsigned lt1 = fp_lt_ct(&bx, &pt->x);
    fp_select(&best, lt1, &bx, &pt->x);
    unsigned lt2 = fp_lt_ct(&nx, &best);
    fp_select(&best, lt2, &nx, &best);
    unsigned sel = (lt2 << 1) | (lt1 & (lt2 ^ 1u));
    ynew = pt->y;
    unsigned char yb[48];
    blst_bendian_from_fp(yb, &ynew);
    unsigned sign = yb[47] & 1;
    blst_fp_cneg(&ynew, &ynew, sign);
    out->x = best;
    out->y = ynew;
    return sel * 2 + sign;
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

static void restart(unsigned i) {
    blst_fr cs = starts[restarts % NSTART];
    restarts++;
    blst_p1_affine pt, cpt;
    comb_single(&pt, &cs, recs[i].fam);
    unsigned idx = canon(&pt, &cpt);
    blst_fr c3;
    blst_fr_mul(&c3, &cs, &mu[idx]);
    recs[i].pt = cpt;
    recs[i].c = c3;
    recs[i].since = 0;
}

static void round_points(void) {
    for (unsigned i = 0; i < RECS; i++) {
        unsigned j = (unsigned)((recs[i].pt.x.l[1] >> 26) & (NMULT - 1));
        blst_fr_mul(&c2s[i], &recs[i].c, &mult[j]);
        fams[i] = recs[i].fam;
    }
    comb_batch(c2s, fams, RECS, ptsbuf);
    for (unsigned i = 0; i < RECS; i++) {
        blst_p1_affine cpt;
        unsigned idx = canon(&ptsbuf[i], &cpt);
        blst_fr c3;
        blst_fr_mul(&c3, &c2s[i], &mu[idx]);
        recs[i].pt = cpt;
        recs[i].c = c3;
        recs[i].since += 1;
        if (recs[i].since >= cap) restart(i);
        const blst_fp *l = &recs[i].pt.x;
        int dp;
        if (dp_bits <= 64) dp = (l->l[0] & ((1ULL << dp_bits) - 1)) == 0;
        else dp = (l->l[0] == 0) && ((l->l[1] & ((1ULL << (dp_bits - 64)) - 1)) == 0);
        if (dp) dp_store(&recs[i].pt, &recs[i].c, recs[i].fam);
    }
    n_rounds++;
}

/* Tables: 2 bases x WINDOWS chains, batched additions per stage j = 2 .. tbl_n - 1 (at most `stages` stages). */
static void build_tables(const blst_p1_affine *P0, const blst_p1_affine *Q, unsigned stages) {
    static blst_p1_affine basek[2][WINDOWS];
    const blst_p1_affine *base[2] = {P0, Q};
    for (int b = 0; b < 2; b++) {
        blst_p1 cur;
        blst_p1_from_affine(&cur, base[b]);
        for (int k = 0; k < WINDOWS; k++) {
            blst_p1_to_affine(&basek[b][k], &cur);
            for (int t = 0; t < WBITS; t++) blst_p1_double(&cur, &cur);
        }
    }
    for (int b = 0; b < 2; b++) {
        for (int k = 0; k < WINDOWS; k++) {
            blst_p1_affine *tb = b ? combA[k] : combP[k];
            tb[0] = basek[b][k];
            blst_p1 pb, pd;
            blst_p1_from_affine(&pb, &basek[b][k]);
            blst_p1_double(&pd, &pb);
            blst_p1_to_affine(&tb[1], &pd);
        }
    }
    unsigned last = 2 + stages;
    if (last > tbl_n) last = tbl_n;
    for (unsigned j = 2; j < last; j++) {
        unsigned n = 0;
        for (int b = 0; b < 2; b++) {
            for (int k = 0; k < WINDOWS; k++) {
                blst_p1_affine *tb = b ? combA[k] : combP[k];
                bP[n] = tb[j - 1];
                bQ[n] = basek[b][k];
                n++;
            }
        }
        batch_add(n);
        n = 0;
        for (int b = 0; b < 2; b++) {
            for (int k = 0; k < WINDOWS; k++) {
                blst_p1_affine *tb = b ? combA[k] : combP[k];
                tb[j] = bres[n];
                n++;
            }
        }
    }
}

static void dump_records(void) {
    for (int w = 0; w < RECS; w++) {
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
        printf("%016llx%016llx\n", (unsigned long long)(recs[w].since >> 64), (unsigned long long)recs[w].since);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: kernel <inputs> <n> [full|forced|forcedfull|forcedcap1|dp4|cap1|build|lite] [dump]\n");
        return 1;
    }
    int n = atoi(argv[2]);
    const char *mode = argc > 3 ? argv[3] : "lite";
    int dump = argc > 4 || strcmp(mode, "dump") == 0;
    int build_only = strcmp(mode, "build") == 0;
    int full = strncmp(mode, "full", 4) == 0 || strcmp(mode, "forcedfull") == 0;
    forced = strcmp(mode, "forced") == 0 || strcmp(mode, "forcedfull") == 0 || strcmp(mode, "forcedcap1") == 0;
    if (strcmp(mode, "dp4") == 0 || strcmp(mode, "fulldp4") == 0) dp_bits = 4;
    cap = (strcmp(mode, "cap1") == 0 || strcmp(mode, "forcedcap1") == 0 || strcmp(mode, "fullcap1") == 0)
              ? (u128)1 : ((u128)1 << 100);
    tbl_n = full ? TBLN_FULL : LITE_N;
    tbl_mask = tbl_n - 1;
    combP = calloc((size_t)WINDOWS * TBLN_FULL, sizeof(blst_p1_affine));
    combA = calloc((size_t)WINDOWS * TBLN_FULL, sizeof(blst_p1_affine));
    if (!combP || !combA) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
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
    for (int i = 0; i < RECS; i++) read_fr(f, &initc[i]);
    fclose(f);
    blst_p1_affine P0 = {x0, y0}, Q = {x1, y1};

#ifndef REF_ONLY
    if (build_only) {
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
    }
#endif
    build_tables(&P0, &Q, build_only ? (unsigned)n : (tbl_n - 2));
#ifndef REF_ONLY
    if (build_only) CALLGRIND_TOGGLE_COLLECT;
#endif
    if (strcmp(mode, "tbl") == 0) {
        unsigned char xb[48];
        blst_bendian_from_fp(xb, &combP[0][0].x); for (int i = 0; i < 4; i++) printf("%02x", xb[i]); printf(" P0.x\n");
        blst_bendian_from_fp(xb, &combP[0][1].x); for (int i = 0; i < 4; i++) printf("%02x", xb[i]); printf(" 2P0.x\n");
        blst_bendian_from_fp(xb, &combP[1][0].x); for (int i = 0; i < 4; i++) printf("%02x", xb[i]); printf(" 2^20P0.x\n");
        blst_bendian_from_fp(xb, &combA[0][0].x); for (int i = 0; i < 4; i++) printf("%02x", xb[i]); printf(" Q.x\n");
        blst_bendian_from_fp(xb, &combP[0][2].x); for (int i = 0; i < 4; i++) printf("%02x", xb[i]); printf(" 3P0.x\n");
        return 0;
    }
    if (build_only) {
        fprintf(stderr, "build done adds=%lu lemma=%lu stages=%d mul=%lu sqr=%lu\n", n_adds, n_lemma, n, n_mul, n_sqr);
        return 0;
    }

    for (unsigned w = 0; w < RECS; w++) {
        recs[w].fam = (uint64_t)(w & 1);
        fams[w] = recs[w].fam;
    }
    static blst_fr cinit[RECS];
    memcpy(cinit, initc, sizeof(cinit));
    comb_batch(cinit, fams, RECS, ptsbuf);
    for (unsigned w = 0; w < RECS; w++) {
        blst_p1_affine cpt;
        unsigned idx = canon(&ptsbuf[w], &cpt);
        blst_fr c3;
        blst_fr_mul(&c3, &cinit[w], &mu[idx]);
        recs[w].pt = cpt;
        recs[w].c = c3;
        recs[w].since = 0;
    }

    int rounds = n / RECS;
    if (rounds * RECS != n) {
        fprintf(stderr, "n must be a multiple of %d\n", RECS);
        return 1;
    }
#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    for (int r = 0; r < rounds; r++) round_points();
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif

    fprintf(stderr, "steps=%d dps=%lu matches=%lu overflow=%lu restarts=%lu adds=%lu lemma=%lu rounds=%lu mul=%lu sqr=%lu\n", n,
            dps, matches, overflow, restarts, n_adds, n_lemma, n_rounds, n_mul, n_sqr);
    if (dump) dump_records();
    return 0;
}
