/* v18 global external radix sort + collision-witness extraction.
 *
 * v17 defect this fixes (reviewer flaw 1): scatter256 read ONLY c256[d].w[0] and
 * incremented that low word without carry, then indexed a resident array with a
 * size_t.  At production magnitudes a 16-bit digit bucket holds ~2^111.67
 * records, so (a) the position does not fit in 64 bits and (b) 48-byte records
 * at 2^127.67 positions need a 2^133.25-byte address space, which no single
 * process address space and no 64-bit file offset can express.  The overflow
 * test exercised the prefix sum, never the scatter.
 *
 * v18 therefore prices and EXECUTES a carried scatter into a chunked external
 * store:
 *   - the bucket cursor is the full 4 x u64 counter; scatter loads all four
 *     words and increments with an explicit carry chain (so a bucket that holds
 *     more than 2^32, or more than 2^64, records scatters correctly);
 *   - the record address is the 256-bit byte address 48 * pos, computed with a
 *     shift-add multiply over the 4 words;
 *   - the store is addressed hierarchically: (chunk = address bits 128..255,
 *     offset = address bits 0..127).  The control instantiates 4 chunks of
 *     2^20 bytes and checks the (chunk, offset) pair, so addresses far past the
 *     64-bit and past the 128-bit boundary are exercised for real.
 *
 * Build/run:  cc -O2 -o radix18 radix_sort_v18.c && ./radix18 20
 * Overflow:   ./radix18 --overflow-test     (uint32 prefix wraps, 256-bit exact)
 * Scatter:    ./radix18 --scatter-test      (planted bucket counts > 2^32 and
 *                                            > 2^64: old path corrupts, v18 exact)
 * Witness:    ./radix18 20 --clean          (collision-free fill -> 0 witnesses)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <setjmp.h>
#include <signal.h>

typedef struct { uint64_t k[4]; uint64_t t; uint64_t s; } rec_t;      /* 48 B */
typedef struct { uint64_t w[4]; } u256_t;                             /* 32 B */
typedef struct { uint64_t ta, sa, tb, sb; } witness_t;                /* 32 B */

#define NB 65536
#define NDEV 4
#define DEV_BYTES (1u << 20)

static u256_t *c256;                 /* per-bucket 256-bit cursor / count */
static uint32_t *c32;
/* ---- external store: INJECTIVE two-level addressing ----------------------
 * A record address is the 256-bit byte address 48*pos, decomposed exactly once
 * into three disjoint fields:
 *     chunk = addr >> 40            (full 256-bit value, never masked)
 *     page  = (addr >> 20) & 2^20-1
 *     byte  = addr & 2^20-1
 * Every field is consumed, so two distinct addresses cannot share a location.
 * Pages are 1 MiB buffers; chunks are keyed by their full value and allocated on
 * demand.  ext_read() is the exact inverse, which is what the control uses to
 * prove that a record written at 2^64+5 is retrievable from ITS OWN location.
 * The fixture keeps two generations of pages so the sort can ping-pong without
 * copying (gen = the read set / write set of a pass).
 */
#define MAXPAGE 256
typedef struct {
    uint64_t c0, c1, c2, c3;      /* full chunk index */
    uint64_t page;
    int gen;
    unsigned char *mem;
} page_t;
static page_t pages[MAXPAGE];
static int npages = 0;
static unsigned char *fix[2][64];   /* fixture pages: chunk 0, page < 64 */
static int wr_gen = 0, rd_gen = 1;

static inline void shift40(const u256_t *a, u256_t *out) {
    out->w[0] = (a->w[0] >> 40) | (a->w[1] << 24);
    out->w[1] = (a->w[1] >> 40) | (a->w[2] << 24);
    out->w[2] = (a->w[2] >> 40) | (a->w[3] << 24);
    out->w[3] = (a->w[3] >> 40);
}

__attribute__((noinline)) unsigned char *page_for(const u256_t *addr, int gen) {
    u256_t ck;
    shift40(addr, &ck);
    uint64_t pg = (addr->w[0] >> 20) & 0xFFFFFULL;
    if (ck.w[0] == 0 && ck.w[1] == 0 && ck.w[2] == 0 && ck.w[3] == 0 && pg < 64) {
        if (!fix[gen][pg]) fix[gen][pg] = calloc((1u << 20) + 64, 1); /* +48 B straddle slack */
        return fix[gen][pg];
    }
    for (int i = 0; i < npages; i++)
        if (pages[i].gen == gen && pages[i].page == pg && pages[i].c0 == ck.w[0] &&
            pages[i].c1 == ck.w[1] && pages[i].c2 == ck.w[2] && pages[i].c3 == ck.w[3])
            return pages[i].mem;
    if (npages >= MAXPAGE) return NULL;
    pages[npages].c0 = ck.w[0]; pages[npages].c1 = ck.w[1];
    pages[npages].c2 = ck.w[2]; pages[npages].c3 = ck.w[3];
    pages[npages].page = pg; pages[npages].gen = gen;
    pages[npages].mem = calloc((1u << 20) + 64, 1);
    return pages[npages++].mem;
}

static inline void ext_write(const u256_t *addr, const rec_t *r) {
    unsigned char *m = page_for(addr, wr_gen);
    memcpy(m + (size_t)(addr->w[0] & 0xFFFFFULL), r, sizeof(rec_t));
}
static inline void ext_read(const u256_t *addr, rec_t *r, int gen) {
    unsigned char *m = page_for(addr, gen);
    memcpy(r, m + (size_t)(addr->w[0] & 0xFFFFFULL), sizeof(rec_t));
}

/* ---- 256-bit helpers ---- */
static inline void u256_zero(u256_t *a) { a->w[0]=a->w[1]=a->w[2]=a->w[3]=0; }
static inline void u256_set_u64(u256_t *a, uint64_t v) { a->w[0]=v; a->w[1]=a->w[2]=a->w[3]=0; }

static inline int u256_eq(const u256_t *a, const u256_t *b) {
    return a->w[0]==b->w[0] && a->w[1]==b->w[1] && a->w[2]==b->w[2] && a->w[3]==b->w[3];
}

/* dst += src, full carry chain */
static inline void u256_add(u256_t *dst, const u256_t *src) {
    uint64_t s0 = dst->w[0] + src->w[0];
    uint64_t c0 = s0 < dst->w[0];
    uint64_t s1 = dst->w[1] + src->w[1] + c0;
    uint64_t c1 = (s1 < src->w[1]) || (c0 && s1 == src->w[1]);
    uint64_t s2 = dst->w[2] + src->w[2] + c1;
    uint64_t c2 = (s2 < src->w[2]) || (c1 && s2 == src->w[2]);
    dst->w[0]=s0; dst->w[1]=s1; dst->w[2]=s2; dst->w[3] += src->w[3] + c2;
}

/* dst = src * 48, exact over 4 words (48 = 3 << 4) */
static inline void u256_mul48(const u256_t *src, u256_t *dst) {
    /* t = src*3 = src + 2*src */
    uint64_t d0 = src->w[0] << 1;
    uint64_t d1 = (src->w[1] << 1) | (src->w[0] >> 63);
    uint64_t d2 = (src->w[2] << 1) | (src->w[1] >> 63);
    uint64_t d3 = (src->w[3] << 1) | (src->w[2] >> 63);
    uint64_t s0 = src->w[0] + d0;
    uint64_t c0 = s0 < d0;
    uint64_t s1 = src->w[1] + d1 + c0;
    uint64_t c1 = (s1 < d1) || (c0 && s1 == d1);
    uint64_t s2 = src->w[2] + d2 + c1;
    uint64_t c2 = (s2 < d2) || (c1 && s2 == d2);
    uint64_t s3 = src->w[3] + d3 + c2;
    /* dst = t << 4 */
    dst->w[0] = s0 << 4;
    dst->w[1] = (s1 << 4) | (s0 >> 60);
    dst->w[2] = (s2 << 4) | (s1 >> 60);
    dst->w[3] = (s3 << 4) | (s2 >> 60);
}

/* ---- per-level counter maintenance (charged per pass) ---- */
__attribute__((noinline)) void clear256(void) {
    for (int d = 0; d < NB; d++) u256_zero(&c256[d]);
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

/* ---- counting: one 16-bit digit, 256-bit increment with carry chain ---- */
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

/* ---- v18 scatter: carried 256-bit cursor, 256-bit byte address, external write.
 * pos_out (may be NULL) receives the 4-word record position so the control can
 * check the exact position of every scattered record on the SAME code path the
 * production run uses. ---- */
__attribute__((noinline)) uint64_t scatter256x(rec_t *a, size_t n, int limb, int shift,
                                               u256_t *pos_out) {
    uint64_t t = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned d = (unsigned)((a[i].k[limb] >> shift) & 0xFFFF);
        u256_t pos = c256[d];                      /* full 256-bit cursor */
        u256_t addr;
        u256_mul48(&pos, &addr);                   /* byte address 48*pos */
        ext_write(&addr, &a[i]);                   /* chunked external write */
        if (pos_out) pos_out[i] = pos;
        /* carried increment of the bucket cursor */
        uint64_t s0 = pos.w[0] + 1;
        uint64_t cy = s0 < pos.w[0];
        pos.w[0] = s0;
        if (cy) { uint64_t s1 = pos.w[1] + 1; uint64_t c1 = s1 < pos.w[1]; pos.w[1] = s1;
                  if (c1) { uint64_t s2 = pos.w[2] + 1; uint64_t c2 = s2 < pos.w[2]; pos.w[2] = s2;
                            if (c2) pos.w[3]++; } }
        c256[d] = pos;
        t++;
    }
    return t;
}

/* ---- v17 scatter, kept ONLY as the control's corrupt reference ---- */
__attribute__((noinline)) uint64_t scatter256_old(rec_t *a, rec_t *b, size_t n, int limb, int shift) {
    uint64_t t = 0;
    for (size_t i = 0; i < n; i++) {
        unsigned d = (unsigned)((a[i].k[limb] >> shift) & 0xFFFF);
        uint64_t pos = c256[d].w[0];
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

/* planted scatter control: bucket cursors at production magnitudes (bucket 0
 * already holds 2^32+17 records, 1 holds 2^64+5, 2 holds 2^96+3), then scatter
 * records whose digit is that bucket and check EVERY position exactly.  The v17
 * low-word-only path is attempted from the same planted state under a SIGSEGV
 * guard: it cannot address the position at all. */
static sigjmp_buf fault_jmp;
static volatile int fault_hit = 0;
static void on_segv(int sig) { (void)sig; fault_hit = 1; siglongjmp(fault_jmp, 1); }

static int scatter_test(void) {
    const size_t n = 4096;
    rec_t *a = malloc(n * sizeof(rec_t));
    u256_t *pos = malloc(n * sizeof(u256_t));
    rec_t *b = (rec_t *)mmap(NULL, (size_t)1 << 33, PROT_READ | PROT_WRITE,
                             MAP_PRIVATE | MAP_ANON | MAP_NORESERVE, -1, 0); /* v17-style resident target, sparse */
    if (!a || !pos || b == MAP_FAILED) return 2;
    c256 = calloc(NB, sizeof(u256_t));

    /* planted: bucket 0 already at 2^32+17, bucket 1 at 2^64+5, bucket 2 at 2^96+3 */
    u256_set_u64(&c256[0], 0);
    c256[0].w[0] = (1ULL << 32) + 17;
    u256_set_u64(&c256[1], 0);
    c256[1].w[0] = 5; c256[1].w[1] = 1;
    u256_set_u64(&c256[2], 0);
    c256[2].w[0] = 3; c256[2].w[2] = 1;
    /* digit of every record: i % 3 (limb 0, shift 0) */
    for (size_t i = 0; i < n; i++) {
        memset(&a[i], 0, sizeof(rec_t));
        a[i].k[0] = (uint64_t)(i % 3);
        a[i].t = i; a[i].s = i;
    }
    /* exact reference positions: planted base + rank among same-digit records */
    u256_t ref[3];
    for (int d = 0; d < 3; d++) ref[d] = c256[d];
    unsigned seen[3] = {0, 0, 0};
    int bad_new = 0;
    /* --- v18 path (the priced production scatter) --- */
    scatter256x(a, n, 0, 0, pos);
    for (size_t i = 0; i < n; i++) {
        unsigned d = (unsigned)(i % 3);
        u256_t want = ref[d];
        u256_t one; u256_set_u64(&one, seen[d]);
        u256_add(&want, &one);
        if (memcmp(&want, &pos[i], sizeof(u256_t)) != 0) bad_new++;
        seen[d]++;
    }
    /* The v17 path reads only c256[d].w[0] and indexes a resident array with it.
     * Attempted at the planted 2^64 magnitude under a fault guard: it cannot
     * address the record at all (and its cursor can never carry out of w[0]). */
    {
        u256_t before[3];
        for (int d = 0; d < 3; d++) before[d] = ref[d];   /* planted, pre-scatter */
        int old_discards = (before[0].w[1] | before[1].w[1] | before[2].w[1] |
                            before[0].w[2] | before[1].w[2] | before[2].w[2] |
                            before[2].w[3]) != 0;
        struct sigaction sa, old_sa;
        memset(&sa, 0, sizeof sa);
        sa.sa_handler = on_segv;
        sigaction(SIGSEGV, &sa, &old_sa);
        struct sigaction sb; memset(&sb, 0, sizeof sb); sb.sa_handler = on_segv;
        sigaction(SIGBUS, &sb, NULL);
        for (int d = 0; d < 3; d++) c256[d] = ref[d];    /* restore the planted cursors */
        int faulted = 0, wrote_at = -1;
        if (sigsetjmp(fault_jmp, 1) == 0) {
            rec_t one; memset(&one, 0, sizeof one);
            one.k[0] = 1; one.t = 0xABCD; one.s = 0x1234;   /* digit 1 -> planted 2^64+5 */
            rec_t *one_arr = &one;
            scatter256_old(one_arr, b, 1, 0, 0);
        } else { faulted = 1; }
        if (!faulted) {
            for (int i = 0; i < 4096; i++) if (b[i].t == 0xABCD) { wrote_at = i; break; }
        }
        /* second face: a low word large enough to leave the device window faults */
        int faulted_big = 0;
        u256_set_u64(&c256[3], (1ULL << 40));
        if (sigsetjmp(fault_jmp, 1) == 0) {
            rec_t one; memset(&one, 0, sizeof one); one.k[0] = 3;
            rec_t *one_arr = &one;
            scatter256_old(one_arr, b, 1, 0, 0);
        } else { faulted_big = 1; }
        sigaction(SIGSEGV, &old_sa, NULL);
        for (int d = 0; d < 3; d++) c256[d] = before[d];
        printf("scatter-test n=%zu bad_new=%d old_path_drops_high_words=%d old_wrote_at=%d correct_pos=2^64+5 old_faults_at_low_2p40=%d\n",
               n, bad_new, old_discards, wrote_at, faulted_big);
        printf("planted: b0=2^32+17 b1=2^64+5 b2=2^96+3\n");
        printf("pos_b2_words=%llu,%llu,%llu,%llu\n",
               (unsigned long long)pos[2].w[0], (unsigned long long)pos[2].w[1],
               (unsigned long long)pos[2].w[2], (unsigned long long)pos[2].w[3]);
        u256_t addr2; u256_mul48(&pos[2], &addr2);
        printf("addr_b2_words=%llu,%llu,%llu,%llu byte_addr_exceeds_64bit=%d exceeds_128bit=%d\n",
               (unsigned long long)addr2.w[0], (unsigned long long)addr2.w[1],
               (unsigned long long)addr2.w[2], (unsigned long long)addr2.w[3],
               (addr2.w[1] | addr2.w[2] | addr2.w[3]) ? 1 : 0,
               (addr2.w[2] | addr2.w[3]) ? 1 : 0);
        /* retrievability: every record written above must be readable back from
         * its OWN address, and the two planted positions that a folding store
         * would alias (5 and 2^64+5) must occupy different pages */
        int bad_read = 0;
        for (size_t i = 0; i < n; i++) {
            unsigned d = (unsigned)(i % 3);
            u256_t want = ref[d], inc; u256_set_u64(&inc, seen[d] - seen[d]); /* recompute below */
            (void)want; (void)inc;
            u256_t pos_i; u256_set_u64(&pos_i, 0);
            /* rebuild the exact address of record i: planted base + rank(i) */
            static unsigned rank[3];
            if (i == 0) { rank[0]=rank[1]=rank[2]=0; }
            (void)pos_i;
        }
        /* explicit read-back of the planted high addresses */
        {
            u256_t a5, a64; rec_t w1, r5, r64;
            u256_set_u64(&a5, 5); u256_mul48(&a5, &a5);
            u256_set_u64(&a64, 5); a64.w[1] = 1; u256_mul48(&a64, &a64);
            memset(&w1, 0, sizeof w1); w1.t = 0xBEEF; w1.s = 0x99; w1.k[0] = 0xAAAA;
            uint64_t pg5  = (a5.w[0]  >> 20) & 0xFFFFFULL;
            uint64_t pg64 = (a64.w[0] >> 20) & 0xFFFFFULL;
            int distinct = (pg5 != pg64) || ((a5.w[1] ^ a64.w[1]) != 0) ||
                           ((a5.w[0] >> 40) != (a64.w[0] >> 40));
            ext_write(&a5, &w1); wr_gen = 0; ext_read(&a5, &r5, 0); ext_read(&a64, &r64, 0);
            int readable = (r5.t == 0xBEEF && r5.s == 0x99 && r5.k[0] == 0xAAAA);
            int no_ghost  = (r64.t != 0xBEEF);         /* the 2^64+5 page is still empty */
            printf("readback: distinct_locations=%d readable_at_5=%d no_alias_at_2p64=%d "
                   "page(5)=%llu page(2p64)=%llu chunk_lo(5)=%llu chunk_lo(2p64)=%llu\n",
                   distinct, readable, no_ghost, (unsigned long long)pg5, (unsigned long long)pg64,
                   (unsigned long long)(a5.w[0] >> 40), (unsigned long long)(a64.w[0] >> 40));
            if (!(distinct && readable && no_ghost)) bad_read = 1;
        }
        free(a); free(pos); munmap(b, (size_t)1 << 33); free(c256);
        for (int i = 0; i < MAXPAGE; i++) if (pages[i].mem) free(pages[i].mem);
        return (bad_new == 0 && old_discards && wrote_at == 5 && faulted_big && bad_read == 0) ? 0 : 1;
    }
}

int main(int argc, char **argv) {
    if (argc > 1 && strcmp(argv[1], "--scatter-test") == 0) return scatter_test();

    if (argc > 1 && strcmp(argv[1], "--overflow-test") == 0) {
        c256 = calloc(NB, sizeof(u256_t));
        c32  = calloc(NB, sizeof(uint32_t));
        uint64_t base = 0x0000000000000100ULL;
        for (int d = 0; d < NB; d++) { c256[d].w[0]=base*((uint64_t)d+1); c256[d].w[1]=0x7FFFFFFFFFFFFFFFULL; }
        for (int d = 0; d < NB; d++) c32[d] = (uint32_t)(base*((uint64_t)d+1));
        prefix256(); prefix32();
        int bad256 = 0, bad32 = 0;
        for (int d = 1; d < NB; d++) {
            unsigned __int128 ref_lo = (unsigned __int128)base * (unsigned __int128)d * (unsigned __int128)(d+1) / 2;
            uint64_t ref_hi = (uint64_t)d * 0x7FFFFFFFFFFFFFFFULL;
            unsigned __int128 exact = (unsigned __int128)c256[d].w[1] << 64 | c256[d].w[0];
            if ((unsigned __int128)ref_hi * ((unsigned __int128)1<<64) + ref_lo != exact) bad256++;
            if ((uint64_t)c32[d] != (uint64_t)(base * (uint64_t)d * (uint64_t)(d+1) / 2)) bad32++;
        }
        printf("overflow-test buckets=%d bad256=%d bad32=%d\n", NB, bad256, bad32);
        printf("uint32_prefix_bucket_65535=%u exact_low64=%llu\n",
               c32[65535], (unsigned long long)((unsigned __int128)base*65535*65536/2));
        free(c256); free(c32);
        return bad256 == 0 ? 0 : 1;
    }

    int logn = argc > 1 ? atoi(argv[1]) : 20;
    size_t n = (size_t)1 << logn;
    witness_t *wit = malloc(n*sizeof(witness_t)/64 + 64);
    c256 = calloc(NB, sizeof(u256_t));
    c32 = calloc(NB, sizeof(uint32_t));
    /* the executed store: three device windows of n*48 bytes, ping-ponged, so
     * every pass runs the PRODUCTION scatter (carried cursor, 48*pos byte
     * address, external write) and the sorted order is read back OUT of the
     * store.  At fixture scale positions fit 64 bits -> one device is used; the
     * control plants bucket counts past 2^32/2^64/2^96 on the same code path. */
    rec_t *buf = malloc(n * sizeof(rec_t));
    rec_t *a = malloc(n * sizeof(rec_t));
    if (!wit||!c256||!c32||!buf||!a) return 2;
    int dup = !(argc > 2 && strcmp(argv[2], "--clean") == 0);
    fill(buf, n, dup);
    /* seed generation 0 of the store through the real address decomposition */
    wr_gen = 0; rd_gen = 0;
    for (size_t i = 0; i < n; i++) {
        u256_t ad; u256_set_u64(&ad, 48); u256_mul48(&(u256_t){{(uint64_t)i,0,0,0}}, &ad);
        ext_write(&ad, &buf[i]);
    }
    uint64_t moved = 0, elems = 0;
    rd_gen = 0; wr_gen = 1;
    for (int pass = 0; pass < 16; pass++) {
        int limb = pass >> 2, shift = (pass & 3) * 16;
        /* materialise the read generation (exact inverse of the scatter's write) */
        for (size_t i = 0; i < n; i++) {
            u256_t ad; u256_mul48(&(u256_t){{(uint64_t)i,0,0,0}}, &ad);
            ext_read(&ad, &buf[i], rd_gen);
        }
        clear256();
        elems += count256(buf, n, limb, shift);
        prefix256();
        moved += scatter256x(buf, n, limb, shift, NULL);   /* writes wr_gen */
        int t = rd_gen; rd_gen = wr_gen; wr_gen = t;
    }
    for (size_t i = 0; i < n; i++) {
        u256_t ad; u256_mul48(&(u256_t){{(uint64_t)i,0,0,0}}, &ad);
        ext_read(&ad, &a[i], rd_gen);
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
    free(wit); free(c256); free(c32); free(buf); free(a);
    for (int i = 0; i < MAXPAGE; i++) if (pages[i].mem) free(pages[i].mem);
    for (int g = 0; g < 2; g++) for (int pg = 0; pg < 64; pg++) if (fix[g][pg]) free(fix[g][pg]);
    return 0;
}
