/* dict_test.c -- 20000-operation differential test of the shipped dictionary
 * against a Python model.  The driver prints every operation, then a digest of
 * the whole table every 1000 ops, plus the hit and eviction counters; the model
 * replays the ops exactly and must reproduce every digest and counter.
 *
 * Compile: clang -O2 -DDICT_TEST dict26.c dict_test.c -o dict_test
 * Mutation control: add -DMUTATE, which flips one bit of the hash inside
 * dict_insert only; the digests must then differ from the model.
 */
#include <stdio.h>
#include <stdint.h>

typedef uint64_t u64;
void state_key48(const u64 v[12], u64 out[6]);
int dict_lookup(const u64 key[6], u64 *tab, u64 mask);
int dict_insert(const u64 key[6], u64 *tab, u64 mask, u64 *evicted);
int dict_lf_ok(u64 n, u64 cap);

#define CAP  (1u << 12)
#define MASK (CAP - 1)
#define OPS  20000
#define SPACE 8192

static u64 rng_state = 0x243F6A8885A308D3ULL;
static u64 xorshift(void) {
    u64 x = rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rng_state = x;
    return x;
}

static u64 tab[CAP * 7];
static u64 keys[SPACE][6];
static int nkeys = 0;

static u64 fnv(const unsigned char *p, size_t n) {
    u64 h = 0xcbf29ce484222325ULL;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 0x100000001b3ULL; }
    return h;
}

int main(void) {
    u64 evicted = 0, hits = 0, miss = 0;
    for (int op = 0; op < OPS; op++) {
        u64 r = xorshift();
        int kind;                       /* 0,1 insert-new; 2 insert; 3 lookup */
        int idx;
        u64 key[6];
        if (nkeys < 8 || ((r % 4) == 0 && nkeys < 2560)) {
            for (int w = 0; w < 6; w++) key[w] = xorshift();
            for (int w = 0; w < 6; w++) keys[nkeys][w] = key[w];
            idx = nkeys++;
            kind = (nkeys <= 8) ? 1 : 0;
        } else {
            idx = (int)((r >> 8) % (u64)nkeys);
            for (int w = 0; w < 6; w++) key[w] = keys[idx][w];
            kind = ((r >> 3) & 1) ? 2 : 3;
        }
        int res;
        if (kind == 3) {
            res = dict_lookup(key, tab, MASK);
            if (res >= 0) hits++; else miss++;
        } else {
            res = dict_insert(key, tab, MASK, &evicted);
        }
        printf("%d %d %d %016llx %016llx %016llx %016llx %016llx %016llx\n",
               op, kind, res, (unsigned long long)key[0], (unsigned long long)key[1],
               (unsigned long long)key[2], (unsigned long long)key[3],
               (unsigned long long)key[4], (unsigned long long)key[5]);
        if ((op + 1) % 1000 == 0)
            printf("D %d %016llx\n", op + 1,
                   (unsigned long long)fnv((const unsigned char *)tab, sizeof(tab)));
    }
    printf("HITS %llu MISS %llu EVICTED %llu LFOK %d\n",
           (unsigned long long)hits, (unsigned long long)miss,
           (unsigned long long)evicted, dict_lf_ok(nkeys, CAP));
    return 0;
}
