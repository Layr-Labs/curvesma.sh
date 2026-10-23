/* blst GT Frobenius-class walk step kernel (walk-spec v11, sections 2 and 4-6).
 * Usage: kernel <inputs.bin> <rounds> [mode]; modes: dump, sqr, mul, restartbench, prfbench, storebench,
 * vectors. Event counts go to stderr. Inputs: gen_inputs.py (header version 6).
 * Test builds only: -DCAP_STEPS_OVERRIDE=n (trail cap), -DTEST_ESCAPE_E0 (escape exponent 0).
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
#define MAGIC 0x3130304f48525447ULL
#define VERSION 6
#define BRANCH_KEY 0x5851F42D4C957F2DULL /* kappa_b (spec section 5, measurement value) */
#define PRF_KEY_TEXT "gt-v11-measurement" /* PRF key KP = SHA-256 of this string (spec section 6, measurement. */
#define COMB_BYTES 32
#define COMB_VALUES 255
#define COMB_ENTRIES (COMB_BYTES * COMB_VALUES)
#ifndef CAP_STEPS_OVERRIDE
#define CAP_STEPS ((unsigned __int128)1 << 66) /* trail cap: CAP_STEPS main steps of one trail (spec section 6) */
#else
#define CAP_STEPS ((unsigned __int128)(CAP_STEPS_OVERRIDE))
#endif

typedef struct { uint64_t v[4]; } u256;
typedef struct { blst_fp12 x; u256 a, b; blst_fp K; int kn; } walk_t;
typedef struct { uint64_t fp[2]; u256 a, b; } dp_t;

static const u256 ORD = {{0xffffffff00000001ULL, 0x53bda402fffe5bfeULL,
                          0x3339d80809a1d805ULL, 0x73eda753299d7d48ULL}};

static const u256 TWO256_MOD_R = {{0x00000001fffffffeULL, 0x5884b7fa00034802ULL,
                                   0x998c4fefecbc4ff5ULL, 0x1824b159acc5056fULL}};
static const uint32_t SHA256_K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
static const uint32_t SHA256_H0[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                      0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};

static const uint8_t *cur;
static uint32_t RA, ADDER_BITS, DP_BITS, S_SQ, NW;
static uint64_t DP_MASK;
static blst_fp12 *T;
static blst_fp12 *COMB;     /* COMB[i * COMB_VALUES + v - 1] = g^(v 2^(8i)) (spec section 6) */
static blst_fr FR_TWO256;
static uint8_t PRF_KEY[32];
static u256 *EC, *ED;
static blst_fp ZETA, G0, G1;
static blst_fr LAMBDA[12];
static int FAM[12], JJ[12], EPS[12];
static walk_t *W;
static uint64_t *RC; /* per-walk restart counter t (spec section 6): incremented before each. */
static blst_fp *prevK, *anchorK;
static int *prevValid, *anchorValid;
static dp_t *DPR;
static uint64_t dp_head;
static unsigned __int128 *AGE;
static uint64_t c_steps, c_adds, c_sqrs, c_dps, c_cyc2, c_cycl, c_restarts, c_cap_restarts, c_chained;
static uint64_t c_escapes, c_escape_dps;
static uint64_t c_comb_mults;

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

static uint64_t fp_lt_mask(const blst_fp *a, const blst_fp *b) {
    uint64_t lt = 0;
    for (int i = 0; i < 6; i++) {
        uint64_t x = a->l[i], y = b->l[i];
        uint64_t ltw = (uint64_t)(x < y);
        uint64_t ne = ltw | (uint64_t)(x > y);
        uint64_t m = 0 - ne;
        lt = (ltw & m) | (lt & ~m);
    }
    return 0 - lt;
}

static int fp_lt(const blst_fp *a, const blst_fp *b) {
    return (int)(fp_lt_mask(a, b) & 1);
}

static int fp_eq(const blst_fp *a, const blst_fp *b) {
    uint64_t diff = 0;
    for (int i = 0; i < 6; i++) diff |= a->l[i] ^ b->l[i];
    return diff == 0;
}

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
        borrow = (uint64_t)(x < y) | ((uint64_t)(x == y) & borrow);
    }

    uint64_t m = 0 - (uint64_t)((c != 0) | (borrow == 0));
    u256 r;
    for (int i = 0; i < 4; i++) r.v[i] = (t.v[i] & m) | (s.v[i] & ~m);
    return r;
}

/* Class key (spec section 2): min over the 12 values, and the index of. */
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
        {
            blst_fp neg;
            blst_fp_cneg(&neg, &val, (bool)EPS[n]);
            val = neg;
        }
        if (n == 0) {
            best = val;
            best_n = 0;
        } else {
            uint64_t m = fp_lt_mask(&val, &best);
            for (int i = 0; i < 6; i++) best.l[i] = (val.l[i] & m) | (best.l[i] & ~m);
            best_n = (int)(((uint64_t)n & m) | ((uint64_t)best_n & ~m));
        }
    }
    *K = best;
    *kn = best_n;
}

__attribute__((noinline)) static void ref_mul(blst_fp12 *r, const blst_fp12 *a, const blst_fp12 *b) {
    blst_fp12_mul(r, a, b);
}

/* Canonical exponent e * lambda^kn mod r (spec section 3): the exponent. */
static u256 canon_exp(const u256 *e, int kn) {
    blst_scalar sc, out;
    blst_fr f;
    u256 res;
    memcpy(sc.b, e->v, 32);
    blst_fr_from_scalar(&f, &sc);
    blst_fr_mul(&f, &f, &LAMBDA[kn]);
    blst_scalar_from_fr(&out, &f);
    memcpy(res.v, out.b, 32);
    return res;
}

static void fingerprint(const blst_fp *K, uint64_t fp[2]) {
    fp[0] = (K->l[1] >> 6) | (K->l[2] << 58);
    fp[1] = (K->l[2] >> 6) | (K->l[3] << 58);
}

static uint64_t dp_bits_of(const blst_fp *K) {
    unsigned pos = ADDER_BITS, li = pos / 64, off = pos % 64;
    uint64_t v = K->l[li] >> off;
    if (off && li + 1 < 6) v |= K->l[li + 1] << (64 - off);
    return v & DP_MASK;
}

/* mix (spec section 5): splitmix64 finaliser, all arithmetic mod 2^64. */
static uint64_t mix64(uint64_t z) {
    z ^= z >> 30;
    z *= 0xBF58476D1CE4E5B9ULL;
    z ^= z >> 27;
    z *= 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return z;
}

/* One group operation of the walk (spec section 5): u = mix(K_lo ^. */
static void step(const walk_t *s, walk_t *o) {
    uint64_t u = mix64(s->K.l[0] ^ BRANCH_KEY);
    uint32_t h = (uint32_t)(u & (RA - 1));
    if ((u >> (64 - ADDER_BITS)) < S_SQ) {
        blst_fp12_cyclotomic_sqr(&o->x, &s->x);
        o->a = addmod(s->a, s->a);
        o->b = addmod(s->b, s->b);
        c_sqrs++;
    } else {
        size_t idx = (size_t)h * 12 + (size_t)s->kn;
        ref_mul(&o->x, &s->x, &T[idx]);
        o->a = addmod(s->a, EC[idx]);
        o->b = addmod(s->b, ED[idx]);
        c_adds++;
    }
    key_of(&o->x, &o->K, &o->kn);
}

static uint32_t rotr32(uint32_t x, unsigned n) {
    return (x >> n) | (x << (32 - n));
}

static void sha256_block(uint32_t h[8], const uint8_t *p) {
    uint32_t w[64];
    for (int i = 0; i < 16; i++)
        w[i] = ((uint32_t)p[4 * i] << 24) | ((uint32_t)p[4 * i + 1] << 16) | ((uint32_t)p[4 * i + 2] << 8) |
               (uint32_t)p[4 * i + 3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = rotr32(w[i - 15], 7) ^ rotr32(w[i - 15], 18) ^ (w[i - 15] >> 3);
        uint32_t s1 = rotr32(w[i - 2], 17) ^ rotr32(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int i = 0; i < 64; i++) {
        uint32_t S1 = rotr32(e, 6) ^ rotr32(e, 11) ^ rotr32(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = hh + S1 + ch + SHA256_K[i] + w[i];
        uint32_t S0 = rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + maj;
        hh = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }
    h[0] += a;
    h[1] += b;
    h[2] += c;
    h[3] += d;
    h[4] += e;
    h[5] += f;
    h[6] += g;
    h[7] += hh;
}

static void sha256(const uint8_t *msg, size_t len, uint8_t out[32]) {
    uint32_t h[8];
    memcpy(h, SHA256_H0, sizeof h);
    size_t i = 0;
    for (; len - i >= 64; i += 64) sha256_block(h, msg + i);
    uint8_t tail[128];
    memset(tail, 0, sizeof tail);
    size_t rem = len - i;
    memcpy(tail, msg + i, rem);
    tail[rem] = 0x80;
    size_t tl = rem < 56 ? 64 : 128;
    uint64_t bits = (uint64_t)len * 8;
    for (int k = 0; k < 8; k++) tail[tl - 1 - k] = (uint8_t)(bits >> (8 * k));
    sha256_block(h, tail);
    if (tl == 128) sha256_block(h, tail + 64);
    for (int k = 0; k < 8; k++) {
        out[4 * k] = (uint8_t)(h[k] >> 24);
        out[4 * k + 1] = (uint8_t)(h[k] >> 16);
        out[4 * k + 2] = (uint8_t)(h[k] >> 8);
        out[4 * k + 3] = (uint8_t)h[k];
    }
}

static void be64(uint8_t *o, uint64_t v) {
    for (int i = 0; i < 8; i++) o[i] = (uint8_t)(v >> (56 - 8 * i));
}

static void fp48_be(uint8_t out[48], const blst_fp *K) {
    for (int i = 0; i < 6; i++) {
        uint64_t v = K->l[i];
        for (int j = 0; j < 8; j++) out[40 - 8 * i + 7 - j] = (uint8_t)(v >> (8 * j));
    }
}

/* PRF(m) (spec section 6): the 64 bytes SHA-256(KP || m || 0x00) ||. */
static void prf_hash(const uint8_t *m, size_t mlen, uint8_t out[64]) {
    uint8_t buf[32 + 64 + 1];
    memcpy(buf, PRF_KEY, 32);
    memcpy(buf + 32, m, mlen);
    buf[32 + mlen] = 0x00;
    sha256(buf, 33 + mlen, out);
    buf[32 + mlen] = 0x01;
    sha256(buf, 33 + mlen, out + 32);
}

static u256 prf_reduce(const uint8_t h[64]) {
    blst_scalar s0, s1, o;
    blst_fr f0, f1, p, sum;
    u256 res;
    blst_scalar_from_bendian(&s0, h);
    blst_scalar_from_bendian(&s1, h + 32);
    blst_fr_from_scalar(&f0, &s0);
    blst_fr_from_scalar(&f1, &s1);
    blst_fr_mul(&p, &f0, &FR_TWO256);
    blst_fr_add(&sum, &p, &f1);
    blst_scalar_from_fr(&o, &sum);
    memcpy(res.v, o.b, 32);
    return res;
}

static u256 restart_prf(uint32_t w, uint64_t t) {
    uint8_t m[17], h[64];
    m[0] = 0x01;
    be64(m + 1, w);
    be64(m + 9, t);
    prf_hash(m, sizeof m, h);
    return prf_reduce(h);
}

static u256 escape_prf(const blst_fp *K) {
    uint8_t m[49], h[64];
    m[0] = 0x02;
    fp48_be(m + 1, K);
    prf_hash(m, sizeof m, h);
    return prf_reduce(h);
}

/* Comb (spec section 6): COMB[i][v] = g^(v 2^(8i)) for i < 32, v =. */
static void comb_build(const blst_fp12 *g) {
    COMB = malloc(sizeof(blst_fp12) * COMB_ENTRIES);
    blst_fp12 base = *g;
    for (int i = 0; i < COMB_BYTES; i++) {
        blst_fp12 *row = &COMB[(size_t)i * COMB_VALUES];
        row[0] = base;
        for (int v = 1; v < COMB_VALUES; v++) blst_fp12_mul(&row[v], &row[v - 1], &base);
        if (i + 1 < COMB_BYTES)
            for (int k = 0; k < 8; k++) blst_fp12_sqr(&base, &base);
    }
}

static uint64_t comb_mul(blst_fp12 *x, const u256 *e) {
    const uint8_t *eb = (const uint8_t *)e->v;
    uint64_t n = 0;
    for (int i = 0; i < COMB_BYTES; i++) {
        uint8_t b = eb[i];
        if (b) {
            blst_fp12_mul(x, x, &COMB[(size_t)i * COMB_VALUES + (b - 1)]);
            n++;
        }
    }
    return n;
}

/* Restart of walk w (spec section 6): the counter t is incremented. */
static void restart_walk(walk_t *c, uint32_t w) {
    uint64_t t = ++RC[w];
    u256 e = restart_prf(w, t);
    c_comb_mults += comb_mul(&c->x, &e);
    c->a = addmod(c->a, e);
    key_of(&c->x, &c->K, &c->kn);
    AGE[w] = 0;
    prevValid[w] = 0;
    anchorValid[w] = 0;
}

static void sigma_pow(blst_fp12 *out, const blst_fp12 *in, uint32_t n) {
    blst_fp12 u = *in, v;
    for (uint32_t i = 0; i < n; i++) {
        blst_fp12_frobenius_map(&v, &u, 1);
        u = v;
    }
    *out = u;
}

/* Escape (spec section 6) from the least-key state m of a detected. */
static void escape(const walk_t *m, walk_t *o) {
    blst_fp12 cs;
    sigma_pow(&cs, &m->x, m->kn);
#ifdef TEST_ESCAPE_E0
    u256 e = {{0, 0, 0, 0}};
#else
    u256 e = escape_prf(&m->K);
#endif
    o->x = cs;
    c_comb_mults += comb_mul(&o->x, &e);
    o->a = addmod(canon_exp(&m->a, m->kn), e);
    o->b = canon_exp(&m->b, m->kn);
    key_of(&o->x, &o->K, &o->kn);
    c_escapes++;
}

/* A distinguished point (spec section 6): the state c is a DP (its test. */
static void chain_dp(walk_t *c, uint32_t w) {
    for (;;) {
        dp_t *rec = &DPR[dp_head % DP_RING];
        fingerprint(&c->K, rec->fp);
        rec->a = canon_exp(&c->a, c->kn);
        rec->b = canon_exp(&c->b, c->kn);
        dp_head++;
        c_dps++;
        restart_walk(c, w);
        c_restarts++;
        if (dp_bits_of(&c->K) != 0) break;
        c_chained++;
    }
}

/* One round of walk w (spec sections 5-6): produce the new point (one. */
static void run_rounds(uint64_t rounds) {
    for (uint64_t r = 0; r < rounds; r++) {
        for (uint32_t w = 0; w < NW; w++) {
            walk_t *c = &W[w];
            walk_t nxt, nn, esc;
            step(c, &nxt);
            c_steps++;
            if (dp_bits_of(&nxt.K) == 0) {
                chain_dp(&nxt, w);
                *c = nxt;
                continue;
            }
            if (++AGE[w] >= CAP_STEPS) {
                restart_walk(&nxt, w);
                c_cap_restarts++;
                if (dp_bits_of(&nxt.K) == 0) chain_dp(&nxt, w);
                *c = nxt;
                continue;
            }
            int two = 0, detected = 0, restarted = 0;
            if (prevValid[w] && fp_eq(&prevK[w], &nxt.K)) {
                step(&nxt, &nn);
                if (fp_eq(&nn.K, &c->K)) two = 1;
            }
            if (two) {
                c_cyc2++;
                escape(fp_lt(&nxt.K, &c->K) ? &nxt : c, &esc);
                detected = 1;
            } else if (anchorValid[w] && fp_eq(&anchorK[w], &nxt.K)) {
                walk_t mn = nxt, t, d;
                step(&nxt, &t);
                int len = 1;
                while (!fp_eq(&t.K, &nxt.K) && len <= MAX_CYCLE) {
                    if (fp_lt(&t.K, &mn.K)) mn = t;
                    step(&t, &d);
                    t = d;
                    len++;
                }
                if (fp_eq(&t.K, &nxt.K)) {
                    c_cycl++;
                    escape(&mn, &esc);
                    detected = 1;
                }
            }
            if (detected) {

                prevValid[w] = 0;
                anchorValid[w] = 0;
                nxt = esc;
                if (dp_bits_of(&nxt.K) == 0) {
                    c_escape_dps++;
                    chain_dp(&nxt, w);
                    restarted = 1;
                }
            } else {
                prevK[w] = c->K;
                prevValid[w] = 1;
            }
            if (!restarted && (r + 1) % ANCHOR_PERIOD == 0) {
                anchorK[w] = nxt.K;
                anchorValid[w] = 1;
            }
            *c = nxt;
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
    if (rd64() != MAGIC || rd64() != VERSION) {
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
    for (int k = 0; k < 12; k++) {
        u256 lp = rd256();
        blst_scalar sc;
        memcpy(sc.b, lp.v, 32);
        blst_fr_from_scalar(&LAMBDA[k], &sc);
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
    RC = calloc(NW, sizeof(uint64_t));
    prevK = malloc(sizeof(blst_fp) * NW);
    anchorK = malloc(sizeof(blst_fp) * NW);
    prevValid = calloc(NW, sizeof(int));
    AGE = calloc(NW, sizeof(unsigned __int128));
    anchorValid = calloc(NW, sizeof(int));
    DPR = malloc(sizeof(dp_t) * DP_RING);
    for (uint32_t w = 0; w < NW; w++) {
        rd_fp12(&W[w].x);
        W[w].a = rd256();
        W[w].b = rd256();
        key_of(&W[w].x, &W[w].K, &W[w].kn);
    }
    /* the generator g (spec section 6, the comb base); the comb table is. */
    blst_fp12 gen;
    rd_fp12(&gen);
    comb_build(&gen);
    if (cur != buf + len) {
        fprintf(stderr, "trailing bytes in the input\n");
        exit(2);
    }
    free(buf);
    sha256((const uint8_t *)PRF_KEY_TEXT, sizeof PRF_KEY_TEXT - 1, PRF_KEY);
    blst_scalar two256;
    memcpy(two256.b, TWO256_MOD_R.v, 32);
    blst_fr_from_scalar(&FR_TWO256, &two256);
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

static void print_hex(const uint8_t *p, size_t n) {
    for (size_t i = 0; i < n; i++) printf("%02x", p[i]);
}

static void print_u256_be(const u256 *x) {
    for (int j = 3; j >= 0; j--) printf("%016lx", (unsigned long)x->v[j]);
}

static uint64_t prng_state;
static uint64_t prng_next(void) {
    prng_state += 0x9E3779B97F4A7C15ULL;
    return mix64(prng_state);
}

static void vectors(void) {
    uint8_t h[32], h64[64], k48[48], in64[64];
    printf("kp ");
    print_hex(PRF_KEY, 32);
    printf("\n");
    sha256((const uint8_t *)"abc", 3, h);
    printf("sha abc ");
    print_hex(h, 32);
    printf("\n");
    sha256((const uint8_t *)"", 0, h);
    printf("sha empty ");
    print_hex(h, 32);
    printf("\n");
    const char *m448 = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    sha256((const uint8_t *)m448, strlen(m448), h);
    printf("sha m448 ");
    print_hex(h, 32);
    printf("\n");
    size_t big = 1000000;
    uint8_t *mil = malloc(big);
    memset(mil, 'a', big);
    sha256(mil, big, h);
    free(mil);
    printf("sha million_a ");
    print_hex(h, 32);
    printf("\n");
    uint64_t ws[3] = {0, 1, 63}, ts[3] = {1, 2, 1000};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            uint8_t m[17];
            m[0] = 0x01;
            be64(m + 1, ws[i]);
            be64(m + 9, ts[j]);
            prf_hash(m, sizeof m, h64);
            printf("prf01 %llu %llu ", (unsigned long long)ws[i], (unsigned long long)ts[j]);
            print_hex(h64, 64);
            printf(" ");
            u256 r = prf_reduce(h64);
            print_u256_be(&r);
            printf("\n");
        }
    prng_state = 0x1234567ULL;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 48; j += 8) {
            uint64_t v = prng_next();
            for (int q = 0; q < 8; q++) k48[j + q] = (uint8_t)(v >> (8 * q));
        }
        uint8_t m[49];
        m[0] = 0x02;
        memcpy(m + 1, k48, 48);
        prf_hash(m, sizeof m, h64);
        u256 r = prf_reduce(h64);
        printf("prf02 ");
        print_hex(k48, 48);
        printf(" ");
        print_hex(h64, 64);
        printf(" ");
        print_u256_be(&r);
        printf("\n");
    }
    for (int i = 0; i < 64; i++) {
        if (i == 0) memset(in64, 0x00, 64);
        else if (i == 1) memset(in64, 0xff, 64);
        else
            for (int j = 0; j < 64; j += 8) {
                uint64_t v = prng_next();
                for (int q = 0; q < 8; q++) in64[j + q] = (uint8_t)(v >> (8 * q));
            }
        u256 r = prf_reduce(in64);
        printf("reduce ");
        print_hex(in64, 64);
        printf(" ");
        print_u256_be(&r);
        printf("\n");
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: kernel <inputs.bin> <rounds> [mode]\n");
        return 2;
    }
    uint64_t rounds = strtoull(argv[2], NULL, 10);
    const char *mode = argc > 3 ? argv[3] : "";
    load(argv[1]);
    if (strcmp(mode, "vectors") == 0) {
        vectors();
        return 0;
    }
    if (strcmp(mode, "restartbench") == 0) {
        /* Cost of one restart (spec section 6): restart_walk on walk 0, as the. */
        walk_t s = W[0];
#ifndef REF_ONLY
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
#endif
        for (uint64_t i = 0; i < rounds; i++) restart_walk(&s, 0);
#ifndef REF_ONLY
        CALLGRIND_TOGGLE_COLLECT;
#endif
        fprintf(stderr, "restartbench=%llu comb_mults=%llu\n", (unsigned long long)rounds,
                (unsigned long long)c_comb_mults);
        return 0;
    }
    if (strcmp(mode, "prfbench") == 0) {

        volatile uint64_t sink = 0;
#ifndef REF_ONLY
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
#endif
        for (uint64_t i = 0; i < rounds; i++) {
            u256 e = restart_prf(0, i + 1);
            sink += e.v[0];
        }
#ifndef REF_ONLY
        CALLGRIND_TOGGLE_COLLECT;
#endif
        fprintf(stderr, "prfbench=%llu\n", (unsigned long long)rounds);
        return 0;
    }
    if (strcmp(mode, "storebench") == 0) {
        /* Cost of one distinguished-point record (spec section 6): the. */
        dp_t rec;
        volatile uint64_t sink = 0;
        walk_t s = W[0];
#ifndef REF_ONLY
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
#endif
        for (uint64_t i = 0; i < rounds; i++) {
            fingerprint(&s.K, rec.fp);
            rec.a = canon_exp(&s.a, s.kn);
            rec.b = canon_exp(&s.b, s.kn);
            sink += rec.fp[0] ^ rec.a.v[0] ^ rec.b.v[0];
        }
#ifndef REF_ONLY
        CALLGRIND_TOGGLE_COLLECT;
#endif
        fprintf(stderr, "storebench=%llu\n", (unsigned long long)rounds);
        return 0;
    }
    if (strcmp(mode, "sqr") == 0 || strcmp(mode, "mul") == 0) {

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
    if (strcmp(mode, "dump") == 0) {
        dump_states();
        for (uint64_t i = 0; i < dp_head && i < DP_RING; i++) {
            const dp_t *rec = &DPR[i];
            printf("dp %016lx %016lx ", (unsigned long)rec->fp[0], (unsigned long)rec->fp[1]);
            for (int j = 3; j >= 0; j--) printf("%016lx", (unsigned long)rec->a.v[j]);
            printf(" ");
            for (int j = 3; j >= 0; j--) printf("%016lx", (unsigned long)rec->b.v[j]);
            printf("\n");
        }
    }
    fprintf(stderr, "steps=%llu adds=%llu sqrs=%llu dps=%llu cyc2=%llu cycl=%llu restarts=%llu cap_restarts=%llu chained=%llu escapes=%llu escape_dps=%llu comb_mults=%llu\n",
            (unsigned long long)c_steps, (unsigned long long)c_adds, (unsigned long long)c_sqrs,
            (unsigned long long)c_dps, (unsigned long long)c_cyc2, (unsigned long long)c_cycl,
            (unsigned long long)c_restarts, (unsigned long long)c_cap_restarts, (unsigned long long)c_chained,
            (unsigned long long)c_escapes, (unsigned long long)c_escape_dps, (unsigned long long)c_comb_mults);
    return 0;
}
