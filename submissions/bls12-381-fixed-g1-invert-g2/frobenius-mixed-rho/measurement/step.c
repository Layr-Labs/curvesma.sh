/* blst GT Frobenius-class walk: step kernel (work-gt/SPEC-gt-walk.md sections 2 and 4-6).
 *
 * Usage: kernel <inputs.bin> <rounds> [dump|sqr]
 *   rounds: each round steps every walk in the input file once (walks = 64).
 *   dump:   print the final walk states to stdout; replica.py must reproduce it exactly.
 *   sqr:    <rounds> cyclotomic squarings of a GT element, to measure kappa.
 * Event counts go to stderr as a final key=value line. Inputs are produced by gen_inputs.py.
 *
 * Differences from the spec, all of them ordering or bookkeeping choices that keep the
 * classes the same:
 *   - the minimum over the 12 values uses the Montgomery integer of each Fp value (a fixed
 *     injection of Fp, so the minimum is still a class invariant);
 *   - the anchor period counts main steps only, not the rare cycle-check steps.
 * Compile with the organizer's flags (measure/container_run.sh): no floating point.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"

#define ANCHOR_PERIOD 64
#define MAX_CYCLE 64
#define DP_RING 65536
#define MAGIC 0x3130304f48525447ULL /* "GTRHO001" little-endian */

typedef struct { uint64_t v[4]; } u256;
typedef struct { blst_fp12 x; u256 a, b; blst_fp K; int kn; } walk_t;
typedef struct { blst_fp K; u256 a, b; int kn; } dp_t;

static const u256 ORD = {{0xffffffff00000001ULL, 0x53bda402fffe5bfeULL,
                          0x3339d80809a1d805ULL, 0x73eda753299d7d48ULL}};

static const uint8_t *cur;
static uint32_t RA, ADDER_BITS, DP_BITS, S_SQ, NW;
static uint64_t DP_MASK;
static blst_fp12 *T;
static u256 *EC, *ED;
static blst_fp ZETA, G0, G1;
static int FAM[12], JJ[12], EPS[12];
static walk_t *W;
static blst_fp *prevK, *anchorK;
static int *prevValid, *anchorValid;
static dp_t *DPR;
static uint64_t dp_head;
static uint64_t c_steps, c_adds, c_sqrs, c_dps, c_cyc2, c_cycl, gstep;

static uint64_t rd64(void) {
    uint64_t x;
    memcpy(&x, cur, 8);
    cur += 8;
    return x;
}

static u256 rd256(void) {
    u256 x;
    for (int i = 0; i < 4; i++) x.v[i] = rd64();
    return x;
}

static void rd_fp(blst_fp *o) {
    blst_fp_from_bendian(o, cur);
    cur += 48;
}

static void rd_fp12(blst_fp12 *o) {
    for (int i = 0; i < 12; i++) rd_fp(&o->fp6[i / 6].fp2[(i % 6) / 2].fp[i % 2]);
}

static int fp_lt(const blst_fp *a, const blst_fp *b) {
    for (int i = 5; i >= 0; i--)
        if (a->l[i] != b->l[i]) return a->l[i] < b->l[i];
    return 0;
}

static int fp_eq(const blst_fp *a, const blst_fp *b) {
    for (int i = 0; i < 6; i++)
        if (a->l[i] != b->l[i]) return 0;
    return 1;
}

/* (a + b) mod r, with a, b < r. */
static u256 addmod(u256 a, u256 b) {
    u256 s, t;
    unsigned __int128 c = 0;
    for (int i = 0; i < 4; i++) {
        c += (unsigned __int128)a.v[i] + b.v[i];
        s.v[i] = (uint64_t)c;
        c >>= 64;
    }
    uint64_t borrow = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t x = s.v[i], y = ORD.v[i];
        t.v[i] = x - y - borrow;
        borrow = (x < y) || (x == y && borrow);
    }
    return (c || !borrow) ? t : s;
}

/* Class key (spec section 2): min over the 12 values, and the index of the minimiser. */
static void key_of(const blst_fp12 *x, blst_fp *K, int *kn) {
    const blst_fp *a0 = &x->fp6[1].fp2[0].fp[0];
    const blst_fp *a1 = &x->fp6[1].fp2[0].fp[1];
    blst_fp B, t, az, bz, fam[2][3], val, best;
    blst_fp_mul(&B, a0, &G0);
    blst_fp_mul(&t, a1, &G1);
    blst_fp_add(&B, &B, &t);
    blst_fp_mul(&az, a0, &ZETA);
    blst_fp_mul(&bz, &B, &ZETA);
    fam[0][0] = *a0;
    fam[0][1] = az;
    blst_fp_sub(&fam[0][2], &az, a0);
    fam[1][0] = B;
    fam[1][1] = bz;
    blst_fp_sub(&fam[1][2], &bz, &B);
    int best_n = -1;
    for (int n = 0; n < 12; n++) {
        val = fam[FAM[n]][JJ[n]];
        if (EPS[n]) {
            blst_fp neg;
            blst_fp_cneg(&neg, &val, 1);
            val = neg;
        }
        if (best_n < 0 || fp_lt(&val, &best)) {
            best = val;
            best_n = n;
        }
    }
    *K = best;
    *kn = best_n;
}

/* The reference operation (1 unit): a dense Fp12 multiplication. Kept as a noinline wrapper
 * so that --toggle-collect=ref_mul isolates it in the REF_FUNCS pass. */
__attribute__((noinline)) static void ref_mul(blst_fp12 *r, const blst_fp12 *a, const blst_fp12 *b) {
    blst_fp12_mul(r, a, b);
}

static uint64_t dp_bits_of(const blst_fp *K) {
    unsigned pos = ADDER_BITS, li = pos / 64, off = pos % 64;
    uint64_t v = K->l[li] >> off;
    if (off && li + 1 < 6) v |= K->l[li + 1] << (64 - off);
    return v & DP_MASK;
}

/* One walk step (spec section 5). count = 0 for the rare cycle-check steps. */
static void step(const walk_t *s, walk_t *o, int count) {
    uint32_t h = (uint32_t)(s->K.l[0] & (RA - 1));
    if (h < S_SQ) {
        blst_fp12_cyclotomic_sqr(&o->x, &s->x);
        o->a = addmod(s->a, s->a);
        o->b = addmod(s->b, s->b);
        if (count) c_sqrs++;
    } else {
        size_t idx = (size_t)h * 12 + (size_t)s->kn;
        ref_mul(&o->x, &s->x, &T[idx]);
        o->a = addmod(s->a, EC[idx]);
        o->b = addmod(s->b, ED[idx]);
        if (count) c_adds++;
    }
    key_of(&o->x, &o->K, &o->kn);
}

/* Escape (spec section 6): square the least-key state of a cycle. */
static void escape(const walk_t *m, walk_t *o) {
    blst_fp12_cyclotomic_sqr(&o->x, &m->x);
    o->a = addmod(m->a, m->a);
    o->b = addmod(m->b, m->b);
    key_of(&o->x, &o->K, &o->kn);
}

static void run_rounds(uint64_t rounds) {
    for (uint64_t r = 0; r < rounds; r++) {
        for (uint32_t w = 0; w < NW; w++) {
            walk_t *c = &W[w];
            walk_t nxt, nn;
            step(c, &nxt, 1);
            c_steps++;
            gstep++;
            int two = 0;
            if (prevValid[w] && fp_eq(&prevK[w], &nxt.K)) {
                step(&nxt, &nn, 0);
                if (fp_eq(&nn.K, &c->K)) two = 1;
            }
            if (two) {
                c_cyc2++;
                escape(fp_lt(&nxt.K, &c->K) ? &nxt : c, &nxt);
            } else if (anchorValid[w] && fp_eq(&anchorK[w], &nxt.K)) {
                walk_t mn = nxt, t, d;
                step(&nxt, &t, 0);
                int len = 1;
                while (!fp_eq(&t.K, &nxt.K) && len <= MAX_CYCLE) {
                    if (fp_lt(&t.K, &mn.K)) mn = t;
                    step(&t, &d, 0);
                    t = d;
                    len++;
                }
                if (fp_eq(&t.K, &nxt.K)) {
                    c_cycl++;
                    escape(&mn, &nxt);
                }
            }
            prevK[w] = c->K;
            prevValid[w] = 1;
            if (gstep % ANCHOR_PERIOD == 0) {
                anchorK[w] = nxt.K;
                anchorValid[w] = 1;
            }
            *c = nxt;
            if (dp_bits_of(&c->K) == 0) {
                dp_t *rec = &DPR[dp_head % DP_RING];
                rec->K = c->K;
                rec->a = c->a;
                rec->b = c->b;
                rec->kn = c->kn;
                dp_head++;
                c_dps++;
            }
        }
    }
}

static void load(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", path);
        exit(2);
    }
    fseek(f, 0, SEEK_END);
    long len = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t *buf = malloc((size_t)len);
    if (fread(buf, 1, (size_t)len, f) != (size_t)len) {
        fprintf(stderr, "short read\n");
        exit(2);
    }
    fclose(f);
    cur = buf;
    if (rd64() != MAGIC || rd64() != 1) {
        fprintf(stderr, "bad magic or version\n");
        exit(2);
    }
    RA = (uint32_t)rd64();
    ADDER_BITS = (uint32_t)rd64();
    DP_BITS = (uint32_t)rd64();
    S_SQ = (uint32_t)rd64();
    NW = (uint32_t)rd64();
    u256 ord = rd256();
    if (memcmp(&ord, &ORD, sizeof ord) != 0) {
        fprintf(stderr, "group order mismatch\n");
        exit(2);
    }
    DP_MASK = (1ULL << DP_BITS) - 1;
    rd_fp(&ZETA);
    rd_fp(&G0);
    rd_fp(&G1);
    for (int n = 0; n < 12; n++) {
        FAM[n] = (int)rd64();
        JJ[n] = (int)rd64();
        EPS[n] = (int)rd64();
    }
    T = malloc(sizeof(blst_fp12) * RA * 12);
    EC = malloc(sizeof(u256) * RA * 12);
    ED = malloc(sizeof(u256) * RA * 12);
    for (size_t i = 0; i < (size_t)RA * 12; i++) {
        rd_fp12(&T[i]);
        EC[i] = rd256();
        ED[i] = rd256();
    }
    W = malloc(sizeof(walk_t) * NW);
    prevK = malloc(sizeof(blst_fp) * NW);
    anchorK = malloc(sizeof(blst_fp) * NW);
    prevValid = calloc(NW, sizeof(int));
    anchorValid = calloc(NW, sizeof(int));
    DPR = malloc(sizeof(dp_t) * DP_RING);
    for (uint32_t w = 0; w < NW; w++) {
        rd_fp12(&W[w].x);
        W[w].a = rd256();
        W[w].b = rd256();
        key_of(&W[w].x, &W[w].K, &W[w].kn);
    }
}

static void dump_states(void) {
    uint8_t buf[48];
    for (uint32_t w = 0; w < NW; w++) {
        const walk_t *s = &W[w];
        printf("%u %016lx%016lx%016lx%016lx %016lx%016lx%016lx%016lx %d ", w,
               (unsigned long)s->a.v[3], (unsigned long)s->a.v[2], (unsigned long)s->a.v[1],
               (unsigned long)s->a.v[0], (unsigned long)s->b.v[3], (unsigned long)s->b.v[2],
               (unsigned long)s->b.v[1], (unsigned long)s->b.v[0], s->kn);
        blst_bendian_from_fp(buf, &s->K);
        for (int i = 0; i < 48; i++) printf("%02x", buf[i]);
        for (int i = 0; i < 12; i++) {
            const blst_fp *f = &s->x.fp6[i / 6].fp2[(i % 6) / 2].fp[i % 2];
            blst_bendian_from_fp(buf, f);
            printf(" ");
            for (int j = 0; j < 48; j++) printf("%02x", buf[j]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: kernel <inputs.bin> <rounds> [dump|sqr]\n");
        return 2;
    }
    uint64_t rounds = strtoull(argv[2], NULL, 10);
    const char *mode = argc > 3 ? argv[3] : "";
    load(argv[1]);
    if (strcmp(mode, "sqr") == 0 || strcmp(mode, "mul") == 0) {
        /* Price loops for the reference multiplication and for kappa. The loop is the
         * organizer's method (cf. scripts/reference_costs/fp12.c): one call per iteration. */
        int is_mul = strcmp(mode, "mul") == 0;
        blst_fp12 buf = T[0];
#ifndef REF_ONLY
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
#endif
        if (is_mul) {
            for (uint64_t i = 0; i < rounds; i++) blst_fp12_mul(&buf, &buf, &T[1]);
        } else {
            for (uint64_t i = 0; i < rounds; i++) blst_fp12_cyclotomic_sqr(&buf, &buf);
        }
#ifndef REF_ONLY
        CALLGRIND_TOGGLE_COLLECT;
#endif
        fprintf(stderr, "%s=%llu\n", is_mul ? "muls" : "sqrs", (unsigned long long)rounds);
        return 0;
    }
#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    run_rounds(rounds);
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif
    if (strcmp(mode, "dump") == 0) dump_states();
    fprintf(stderr, "steps=%llu adds=%llu sqrs=%llu dps=%llu cyc2=%llu cycl=%llu\n",
            (unsigned long long)c_steps, (unsigned long long)c_adds, (unsigned long long)c_sqrs,
            (unsigned long long)c_dps, (unsigned long long)c_cyc2, (unsigned long long)c_cycl);
    return 0;
}
