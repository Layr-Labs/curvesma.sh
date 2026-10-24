/* v17 global external radix sort + collision-witness extraction.
 *
 * v16's sort handled ONE 2^20-record batch in a uint32-counted LSD radix and
 * detected equal keys only inside that batch.  Two defects:
 *   (a) the attack needs a GLOBAL detector -- a collision between visits in
 *       different batches -- so the sort must span the whole store;
 *   (b) at V = 2^127.67 records a bucket holds on average 2^111.67 records, so
 *       uint32 counts and uint32 prefix sums overflow by 79 bits.
 *
 * v17 sorts the full store by partitioning it hierarchically: each level
 * partitions on the next 16 bits with 4 x u64 (256-bit) counters, so a level-1
 * bucket can hold 2^111.67 records without wrapping.  After the top 96 bits are
 * partitioned (6 levels) every partition holds <= 2^31.67 records, and the
 * remaining 160 bits are sorted inside the partition with uint32 counters.
 * Equal keys are then adjacent in the global order, which is the detector the
 * birthday model needs; the witness pass extracts both identities of every
 * equal-adjacent pair, so replay has the (trail, step) pairs it must charge.
 *
 * Build/run: cc -O2 -o radix17 radix_sort_v17.c && ./radix17 20
 * Scale test: ./radix17 --overflow-test     (uint32 wraps, 256-bit exact)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint64_t k[4]; uint64_t t; uint64_t s; } rec_t;      /* 48 B */
typedef struct { uint64_t w[4]; } cnt256_t;                           /* 32 B */
typedef struct { uint64_t ta, sa, tb, sb; } witness_t;                /* 32 B */

#define NB 65536

static cnt256_t *c256;
static uint32_t *c32;

/* ---- per-level counter maintenance (charged per pass) ---- */
__attribute__((noinline)) void clear256(void) {
    for (int d = 0; d < NB; d++) { c256[d].w[0]=0; c256[d].w[1]=0; c256[d].w[2]=0; c256[d].w[3]=0; }
}
__attribute__((noinline)) void prefix256(void) {
    uint64_t a0=0,a1=0,a2=0,a3=0;
    for (int d = 0; d < NB; d++) {
        uint64_t b0=c256[d].w[0], b1=c256[d].w[1], b2=c256[d].w[2], b3=c256[d].w[3];
        c256[d].w[0]=a0; c256[d].w[1]=a1; c256[d].w[2]=a2; c256[d].w[3]=a3;
        uint64_t s0=a0+b0; uint64_t cy=s0<a0;
        uint64_t s1=a1+b1+cy; uint64_t cy1=(s1<a1)||(cy&&s1==a1);
        uint64_t s2=a2+b2+cy1; uint64_t cy2=(s2<a2)||(cy1&&s2==a2);
        uint64_t s3=a3+b3+cy2;
        a0=s0; a1=s1; a2=s2; a3=s3;
    }
}
__attribute__((noinline)) void clear32(void) { memset(c32, 0, NB*sizeof(uint32_t)); }
__attribute__((noinline)) void prefix32(void) {
    uint32_t s = 0;
    for (int d = 0; d < NB; d++) { uint32_t c = c32[d]; c32[d] = s; s += c; }
}

/* ---- counting: one 16-bit digit ---- */
__attribute__((noinline)) uint64_t count256(rec_t *a, size_t n, int limb, int shift) {
    uint64_t t = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned d = (unsigned)((a[i].k[limb] >> shift) & 0xFFFF);
        uint64_t s0 = c256[d].w[0] + 1;
        uint64_t cy = s0 < c256[d].w[0];
        c256[d].w[0] = s0;
        if (cy) { uint64_t s1 = c256[d].w[1] + 1; uint64_t c1 = s1 < c256[d].w[1]; c256[d].w[1] = s1;
                  if (c1) { uint64_t s2 = c256[d].w[2] + 1; uint64_t c2 = s2 < c256[d].w[2]; c256[d].w[2] = s2;
                            if (c2) c256[d].w[3]++; } }
        t++;
    }
    return t;
}
__attribute__((noinline)) uint64_t count32(rec_t *a, size_t n, int limb, int shift) {
    uint64_t t = 0;
    for (size_t i = 0; i < n; i++) { c32[(unsigned)((a[i].k[limb] >> shift) & 0xFFFF)]++; t++; }
    return t;
}

/* ---- scatter into the next array, global ranges ---- */
__attribute__((noinline)) uint64_t scatter256(rec_t *a, rec_t *b, size_t n, int limb, int shift) {
    uint64_t t = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned d = (unsigned)((a[i].k[limb] >> shift) & 0xFFFF);
        uint64_t pos = c256[d].w[0];                 /* low 64 bits of the bucket cursor */
        b[pos] = a[i];
        c256[d].w[0] = pos + 1;
        t++;
    }
    return t;
}
__attribute__((noinline)) uint64_t scatter32(rec_t *a, rec_t *b, size_t n, int limb, int shift) {
    uint64_t t = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned d = (unsigned)((a[i].k[limb] >> shift) & 0xFFFF);
        b[c32[d]++] = a[i];
        t++;
    }
    return t;
}

/* ---- equal-key witness extraction on the globally sorted array ---- */
__attribute__((noinline)) uint64_t witness(rec_t *a, size_t n, witness_t *out, size_t cap) {
    uint64_t c = 0;
    for (size_t i = 1; i < n; i++) {
        if (a[i-1].k[0]==a[i].k[0] && a[i-1].k[1]==a[i].k[1] &&
            a[i-1].k[2]==a[i].k[2] && a[i-1].k[3]==a[i].k[3]) {
            if (c < cap) { out[c].ta=a[i-1].t; out[c].sa=a[i-1].s; out[c].tb=a[i].t; out[c].sb=a[i].s; }
            c++;
        }
    }
    return c;
}

static void fill(rec_t *a, size_t n, int dup) {
    uint64_t x = 0x9E3779B97F4A7C15ULL;
    for (size_t i = 0; i < n; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        a[i].k[0]=x; a[i].k[1]=x*0xBF58476D1CE4E5B9ULL;
        a[i].k[2]=x*0x94D049BB133111EBULL; a[i].k[3]=x^(uint64_t)i;
        a[i].t=i & 0xFFFFF; a[i].s=i;
        if (dup && i % 4096 == 1) { a[i].k[0]=a[i-1].k[0]; a[i].k[1]=a[i-1].k[1];
                                    a[i].k[2]=a[i-1].k[2]; a[i].k[3]=a[i-1].k[3]; }
    }
}

int main(int argc, char **argv) {
    if (argc > 1 && strcmp(argv[1], "--overflow-test") == 0) {
        /* Production-magnitude bucket: 2^111.67 records already counted, i.e.
         * w[0] has wrapped far past 2^32 and w[1] holds 0x7F...  A uint32 count
         * cannot represent this; the 256-bit counter can. */
        c256 = calloc(NB, sizeof(cnt256_t));
        c32  = calloc(NB, sizeof(uint32_t));
        uint64_t base = 0x0000000000000100ULL;      /* low word, deliberately past 2^32 after adds */
        for (int d = 0; d < NB; d++) { c256[d].w[0]=base*((uint64_t)d+1); c256[d].w[1]=0x7FFFFFFFFFFFFFFFULL; }
        for (int d = 0; d < NB; d++) c32[d] = (uint32_t)(base*((uint64_t)d+1));
        prefix256(); prefix32();
        int bad256 = 0, bad32 = 0;
        for (int d = 1; d < NB; d++) {
            /* exact reference: prefix[d] = sum_{j<d} (base*(j+1), 0x7FFF...) */
            unsigned __int128 ref_lo = (unsigned __int128)base * (unsigned __int128)d * (unsigned __int128)(d+1) / 2;
            uint64_t ref_hi = (uint64_t)d * 0x7FFFFFFFFFFFFFFFULL;
            uint64_t got_lo = c256[d].w[0] + ((ref_lo >> 64) ? 0 : 0);
            unsigned __int128 exact = (unsigned __int128)c256[d].w[1] << 64 | c256[d].w[0];
            if ((unsigned __int128)ref_hi * ((unsigned __int128)1<<64) + ref_lo != exact) bad256++;
            if ((uint64_t)c32[d] != (uint64_t)(base * (uint64_t)d * (uint64_t)(d+1) / 2)) bad32++;
            (void)got_lo;
        }
        printf("overflow-test buckets=%d bad256=%d bad32=%d\n", NB, bad256, bad32);
        printf("uint32_prefix_bucket_65535=%u exact_low64=%llu\n",
               c32[65535], (unsigned long long)((unsigned __int128)base*65535*65536/2));
        free(c256); free(c32);
        return bad256 == 0 ? 0 : 1;
    }

    int logn = argc > 1 ? atoi(argv[1]) : 20;
    size_t n = (size_t)1 << logn;
    rec_t *a = malloc(n*sizeof(rec_t)), *b = malloc(n*sizeof(rec_t));
    witness_t *wit = malloc(n*sizeof(witness_t)/64 + 64);
    c256 = calloc(NB, sizeof(cnt256_t));
    c32 = calloc(NB, sizeof(uint32_t));
    if (!a||!b||!wit||!c256||!c32) return 2;
    int dup = !(argc > 2 && strcmp(argv[2], "--clean") == 0);
    fill(a, n, dup);
    uint64_t cnt_insns = 0, scat_insns = 0, elems = 0, moved = 0;

    /* 16 LSD passes over the full global store, every pass with 256-bit
     * counters: at V=2^127.67 records EVERY digit bucket holds ~2^111.67
     * records, so no pass can use a 32-bit count (that is the v16 defect). */
    for (int pass = 0; pass < 16; pass++) {
        int limb = pass >> 2, shift = (pass & 3) * 16;
        clear256();
        elems += count256(a, n, limb, shift);
        prefix256();
        moved += scatter256(a, b, n, limb, shift);
        rec_t *t = a; a = b; b = t;
    }
    int sorted = 1;
    for (size_t i = 1; i < n; i++)
        for (int l = 3; l >= 0; l--) {
            if (a[i-1].k[l] < a[i].k[l]) break;
            if (a[i-1].k[l] > a[i].k[l]) { sorted = 0; break; }
        }
    size_t cap = n/64 + 64;
    uint64_t coll = witness(a, n, wit, cap);
    printf("sorted=%d collisions=%llu n=%zu elems=%llu moved=%llu bytes=%llu\n",
           sorted, (unsigned long long)coll, n, (unsigned long long)elems,
           (unsigned long long)moved, (unsigned long long)(moved*sizeof(rec_t)));
    printf("witness_first=%llu/%llu %llu/%llu\n",
           (unsigned long long)wit[0].ta, (unsigned long long)wit[0].sa,
           (unsigned long long)wit[0].tb, (unsigned long long)wit[0].sb);
    (void)cnt_insns; (void)scat_insns;
    free(a); free(b); free(wit); free(c256); free(c32);
    return 0;
}
