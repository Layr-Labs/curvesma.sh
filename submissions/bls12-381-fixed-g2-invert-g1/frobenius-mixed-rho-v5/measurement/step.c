/* blst GT Frobenius-class walk: step kernel (work-gt/SPEC-gt-walk.md sections 2 and 4-6).
 *
 * Usage: kernel <inputs.bin> <rounds> [dump|sqr|mul|restartbench]
 *   rounds: each round steps every walk in the input file once (walks = 64).
 *   dump:   print the final walk states to stdout; replica.py must reproduce it exactly.
 *   sqr|mul: <rounds> reference operations, to measure kappa and the reference multiplication.
 *   restartbench: <rounds> restarts of walk 0 (the restart alone, spec section 6).
 * Event counts go to stderr as a final key=value line. Inputs are produced by gen_inputs.py.
 *
 * Differences from the spec, all of them ordering or bookkeeping choices that keep the
 * classes the same:
 *   - the minimum over the 12 values uses the Montgomery integer of each Fp value (a fixed
 *     injection of Fp, so the minimum is still a class invariant);
 *   - cycle detection is per walk: a walk's anchor is refreshed after every ANCHOR_PERIOD of
 *     its own steps (each round steps every walk once, so this is the round counter);
 *   - a distinguished point stores a 16-byte fingerprint of the key (the low 128 bits of its
 *     Montgomery integer) and the two exponents: 80 bytes per record;
 *   - a distinguished point ends the trail: after the record the walk restarts with the restart element
 *     Y_j of the restart table, j = mix((w + 1) * GAMMA + t) mod TR, t the walk's restart counter
 *     (incremented before use) and w its global index (spec section 6, v10); a trail that reaches the
 *     trail cap restarts the same way, without a record.
 * Event counts: steps counts the main walk steps; adds and sqrs count every group operation,
 * including the extra steps of cycle checks, and escapes (squarings); restarts counts the
 * restarts at a distinguished point (including chained ones) and cap_restarts those at the trail cap;
 * chained counts the restarts whose state is itself a distinguished point (spec section 6).
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
#define VERSION 5
#define GAMMA 0x9E3779B97F4A7C15ULL /* restart index multiplier (spec section 6) */
#define TR 1024 /* restart table size T_R (spec section 6, v10) */
#ifndef CAP_STEPS_OVERRIDE
#define CAP_STEPS ((unsigned __int128)1 << 66) /* trail cap: 64 x 2^DP_BITS steps without a distinguished point (spec section 6) */
#else
#define CAP_STEPS ((unsigned __int128)(CAP_STEPS_OVERRIDE)) /* test builds only: a low cap exercises the restart at the cap */
#endif

typedef struct { uint64_t v[4]; } u256;
typedef struct { blst_fp12 x; u256 a, b; blst_fp K; int kn; } walk_t;
typedef struct { uint64_t fp[2]; u256 a, b; } dp_t;

static const u256 ORD = {{0xffffffff00000001ULL, 0x53bda402fffe5bfeULL,
                          0x3339d80809a1d805ULL, 0x73eda753299d7d48ULL}};

static const uint8_t *cur;
static uint32_t RA, ADDER_BITS, DP_BITS, S_SQ, NW;
static uint64_t DP_MASK;
static blst_fp12 *T;
static blst_fp12 *RY; /* restart table: Y_j = g^(f_j), j < TR (spec section 6, v10) */
static u256 *RF;      /* restart table: f_j mod r */
static u256 *EC, *ED;
static blst_fp ZETA, G0, G1;
static blst_fr LAMBDA[12]; /* lambda^k in the Montgomery form of Fr */
static int FAM[12], JJ[12], EPS[12];
static walk_t *W;
static uint64_t *RC; /* per-walk restart counter t (spec section 6, v9): incremented before each restart */
static blst_fp *prevK, *anchorK;
static int *prevValid, *anchorValid;
static dp_t *DPR;
static uint64_t dp_head;
static unsigned __int128 *AGE; /* main steps since the walk's last distinguished point or restart */
static uint64_t c_steps, c_adds, c_sqrs, c_dps, c_cyc2, c_cycl, c_restarts, c_cap_restarts, c_chained;
static uint64_t c_escapes, c_escape_dps; /* escapes (one per detected cycle); escape targets that are DPs */

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

/* Constant-instruction comparisons: every limb is visited, with no data-dependent branch or exit.
 * fp_lt_mask returns all ones when a < b (as integers), else zero; the limb loop runs from the least
 * significant limb up, so the most significant differing limb decides. */
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
        borrow = (uint64_t)(x < y) | ((uint64_t)(x == y) & borrow);
    }
    /* Select t when the sum overflowed or did not borrow, without a branch. */
    uint64_t m = 0 - (uint64_t)((c != 0) | (borrow == 0));
    u256 r;
    for (int i = 0; i < 4; i++) r.v[i] = (t.v[i] & m) | (s.v[i] & ~m);
    return r;
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

/* The reference operation (1 unit): a dense Fp12 multiplication. Kept as a noinline wrapper
 * so that --toggle-collect=ref_mul isolates it in the REF_FUNCS pass. */
__attribute__((noinline)) static void ref_mul(blst_fp12 *r, const blst_fp12 *a, const blst_fp12 *b) {
    blst_fp12_mul(r, a, b);
}

/* Canonical exponent e * lambda^kn mod r (spec section 6): the exponent of sigma^kn(x). */
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

/* Fingerprint of a distinguished point: key bits 70..197 of the Montgomery integer. Bits 10..69 are
 * zero for every distinguished point (the test), so bits 70..197 are the 128 that vary. */
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

/* One group operation of the walk (spec section 5): one squaring or one adder multiplication. */
static void step(const walk_t *s, walk_t *o) {
    uint32_t h = (uint32_t)(s->K.l[0] & (RA - 1));
    if (h < S_SQ) {
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

/* mix (spec section 6, v9): splitmix64 finaliser, all arithmetic mod 2^64. */
static uint64_t mix64(uint64_t z) {
    z ^= z >> 30;
    z *= 0xBF58476D1CE4E5B9ULL;
    z ^= z >> 27;
    z *= 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return z;
}

/* Restart of walk w (spec section 6, v10): the walk's counter t is incremented before use, the restart
 * index is j = mix((w + 1) * GAMMA + t) mod TR, and the walk moves by the restart element Y_j = g^(f_j) of the
 * restart table (not an adder): x <- x * Y_j, a <- a + f_j, b unchanged. The detector state (previous key,
 * anchor, trail age) is cleared. */
static void restart_walk(walk_t *c, uint32_t w) {
    uint64_t t = ++RC[w];
    uint32_t j = (uint32_t)(mix64((uint64_t)(w + 1) * GAMMA + t) & (TR - 1));
    ref_mul(&c->x, &c->x, &RY[j]);
    c->a = addmod(c->a, RF[j]);
    key_of(&c->x, &c->K, &c->kn);
    AGE[w] = 0;
    prevValid[w] = 0;
    anchorValid[w] = 0;
}

/* sigma^n(x) = x^(p^n) for n < 12: n applications of the p-power map (blst's Frobenius map is exact for n = 1). */
static void sigma_pow(blst_fp12 *out, const blst_fp12 *in, uint32_t n) {
    blst_fp12 u = *in, v;
    for (uint32_t i = 0; i < n; i++) {
        blst_fp12_frobenius_map(&v, &u, 1);
        u = v;
    }
    *out = u;
}

/* Escape (spec section 6, v10): from the least-key state s of a detected cycle, restart at c_s * Y_j with the canonical
 * representative c_s = sigma^kn(s), the canonical exponents (lambda^kn a_s + f_j, lambda^kn b_s) and
 * j = mix(low64(K_s)) mod TR. It uses no restart counter. The caller clears the detector. */
static void escape(const walk_t *m, walk_t *o) {
    blst_fp12 cs;
    sigma_pow(&cs, &m->x, m->kn);
    uint32_t j = (uint32_t)(mix64(m->K.l[0]) & (TR - 1));
    ref_mul(&o->x, &cs, &RY[j]);
    o->a = addmod(canon_exp(&m->a, m->kn), RF[j]);
    o->b = canon_exp(&m->b, m->kn);
    key_of(&o->x, &o->K, &o->kn);
    c_escapes++;
}

/* A distinguished point (spec section 6): the state c is a DP (its test passed); its record is stored before the
 * restart. The restarted state may itself be a distinguished point: it is stored the same way and the walk restarts
 * again, until the state is not a distinguished point (chained restarts). */
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

/* One round of walk w (spec sections 5-6, v10): produce the new point; the DP test comes first (record and restart);
 * then the detector (2-cycle, anchor): a detected cycle is escaped; the trail cap restarts a trail that runs CAP_STEPS
 * steps. The anchor refreshes on the global round, unless the walk restarted in this round. */
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
                /* the escape clears the detector; its target is a DP candidate (chained, as a restart) */
                prevValid[w] = 0;
                anchorValid[w] = 0;
                AGE[w] = 0;
                nxt = esc;
                if (dp_bits_of(&nxt.K) == 0) {
                    c_escape_dps++;
                    chain_dp(&nxt, w);
                    restarted = 1;
                }
            } else {
                prevK[w] = c->K;
                prevValid[w] = 1;
                /* Trail cap (spec section 6): a trail that runs CAP_STEPS steps without a DP restarts as at a DP,
                 * storing no record (chained if the restarted state is a DP). Unreachable at 2^66 steps. */
                if (++AGE[w] >= CAP_STEPS) {
                    restart_walk(&nxt, w);
                    c_cap_restarts++;
                    if (dp_bits_of(&nxt.K) == 0) chain_dp(&nxt, w);
                    restarted = 1;
                }
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
    /* restart table block (appended after the walks): count TR, then TR elements Y_j, then TR exponents f_j */
    if (rd64() != TR) {
        fprintf(stderr, "restart table size mismatch\n");
        exit(2);
    }
    RY = malloc(sizeof(blst_fp12) * TR);
    RF = malloc(sizeof(u256) * TR);
    for (uint32_t j = 0; j < TR; j++) rd_fp12(&RY[j]);
    for (uint32_t j = 0; j < TR; j++) RF[j] = rd256();
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
        fprintf(stderr, "usage: kernel <inputs.bin> <rounds> [dump|sqr|mul]\n");
        return 2;
    }
    uint64_t rounds = strtoull(argv[2], NULL, 10);
    const char *mode = argc > 3 ? argv[3] : "";
    load(argv[1]);
    if (strcmp(mode, "sortbench") == 0) {
        /* Cost of one comparison and one record exchange of the sort in lookup (spec section 9): a
         * constant-instruction 128-bit comparison and an 80-byte swap, per iteration. */
        volatile uint64_t sink = 0;
#ifndef REF_ONLY
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
#endif
        for (uint64_t i = 0; i < rounds; i++) {
            dp_t *a = &DPR[i & (DP_RING - 1)], *b = &DPR[(i + 1) & (DP_RING - 1)];
            uint64_t hi = (uint64_t)(b->fp[1] < a->fp[1]);
            uint64_t eq = (uint64_t)(b->fp[1] == a->fp[1]);
            uint64_t lo = (uint64_t)(b->fp[0] < a->fp[0]);
            uint64_t m = 0 - (hi | (eq & lo));
            dp_t tmp = *a;
            *a = *b;
            *b = tmp;
            sink += m;
        }
#ifndef REF_ONLY
        CALLGRIND_TOGGLE_COLLECT;
#endif
        fprintf(stderr, "sortbench=%llu\n", (unsigned long long)rounds);
        return 0;
    }
    if (strcmp(mode, "restartbench") == 0) {
        /* Cost of one restart (spec section 6): restart_walk on walk 0, as the walk does it. */
        walk_t s = W[0];
#ifndef REF_ONLY
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
#endif
        for (uint64_t i = 0; i < rounds; i++) restart_walk(&s, 0);
#ifndef REF_ONLY
        CALLGRIND_TOGGLE_COLLECT;
#endif
        fprintf(stderr, "restartbench=%llu\n", (unsigned long long)rounds);
        return 0;
    }
    if (strcmp(mode, "storebench") == 0) {
        /* Cost of one distinguished-point record (spec section 6): the fingerprint and the two canonical
         * exponents of walk 0's state. The volatile sink keeps the record computed. */
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
    fprintf(stderr, "steps=%llu adds=%llu sqrs=%llu dps=%llu cyc2=%llu cycl=%llu restarts=%llu cap_restarts=%llu chained=%llu escapes=%llu escape_dps=%llu\n",
            (unsigned long long)c_steps, (unsigned long long)c_adds, (unsigned long long)c_sqrs,
            (unsigned long long)c_dps, (unsigned long long)c_cyc2, (unsigned long long)c_cycl,
            (unsigned long long)c_restarts, (unsigned long long)c_cap_restarts, (unsigned long long)c_chained,
            (unsigned long long)c_escapes, (unsigned long long)c_escape_dps);
    return 0;
}
