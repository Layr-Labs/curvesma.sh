/* R3: the whole relation-collection record route, counted in executed AArch64 instructions
   (STORAGE.md, measurements12.txt "R3" lines). Emitter: direct chunk scatter at emission (the
   STORAGE.md hot path); a full page goes to the chunk's page directory, a fresh page is taken.
   Transfer: the owner copies every page (16 B/record) into its receive area, kept for recovery.
   Owner: walks its directory page by page into apply. modes: 0 source, 2 +append, 3 +apply, 1 route */
#include <stdint.h>
#include <stdlib.h>
typedef struct { uint64_t cell, ideal; } rec_t;
#ifndef PAGE_BITS
#define PAGE_BITS 16                     /* 2^16 records = 1 MiB pages */
#endif
#define PAGE_RECS (1ull << PAGE_BITS)
#define PAGE_MASK ((PAGE_RECS << 4) - 1)
#define CHUNK_BITS 30
#define SUM_MASK ((1ull << CHUNK_BITS) - 1)   /* production mask; stand-in offsets stay below 2^16 */
typedef struct { rec_t **page; uint64_t n; } dir_t;
static rec_t *pool; static uint64_t pool_used; static dir_t *dir;
__attribute__((noinline)) rec_t *emit_page(uint64_t k, rec_t *end) {
  dir_t *d = &dir[k]; d->page[d->n++] = end - PAGE_RECS;
  rec_t *p = pool + pool_used; pool_used += PAGE_RECS; return p;
}
static inline void emit(rec_t **tail, uint64_t cell, uint64_t ideal) {
  uint64_t k = cell >> CHUNK_BITS; rec_t *t = tail[k];
  t->cell = cell; t->ideal = ideal; t++;
  if (__builtin_expect(((uintptr_t)t & PAGE_MASK) == 0, 0)) t = emit_page(k, t);
  tail[k] = t;
}
#define LCG(x) ((x) * 6364136223846793005ULL + 1442695040888963407ULL)
#define CELL(x) ((((x) >> 62) << CHUNK_BITS) | (((x) >> 20) & 0xffff))   /* 4 chunks */
__attribute__((noinline)) void walk_scatter(rec_t **tail, uint64_t n, uint64_t x) {
  for (uint64_t i = 0; i < n; i++) { x = LCG(x); emit(tail, CELL(x), i); }
}
__attribute__((noinline)) rec_t *walk_linear(rec_t *t, uint64_t n, uint64_t x) {
  for (uint64_t i = 0; i < n; i++) { x = LCG(x); t->cell = CELL(x); t->ideal = i; t++; }
  return t;
}
__attribute__((noinline)) void walk_nostore(volatile uint64_t *sink, uint64_t n, uint64_t x) {
  uint64_t acc = 0;
  for (uint64_t i = 0; i < n; i++) { x = LCG(x); acc += CELL(x) ^ i; }
  *sink = acc;
}
typedef unsigned __int128 u128;
__attribute__((noinline)) void recv_page(rec_t *restrict dst, const rec_t *restrict src, uint64_t n) {
  u128 *d = (u128 *)dst; const u128 *s = (const u128 *)src;
  for (uint64_t i = 0; i < n; i += 4) { u128 a = s[i], b = s[i+1], c = s[i+2], e = s[i+3]; d[i] = a; d[i+1] = b; d[i+2] = c; d[i+3] = e; }
}
__attribute__((noinline)) void apply(uint32_t *sums, const rec_t *r, uint64_t n) {
  for (uint64_t i = 0; i < n; i++) sums[r[i].cell & SUM_MASK] += (uint32_t)r[i].ideal & 0xff;
}
__attribute__((noinline)) void owner_run(dir_t *d, rec_t *open, rec_t *tail, rec_t *rbuf, uint32_t *sums) {
  rec_t *r = rbuf;
  for (uint64_t j = 0; j < d->n; j++, r += PAGE_RECS) { recv_page(r, d->page[j], PAGE_RECS); apply(sums, r, PAGE_RECS); }
  uint64_t m = (uint64_t)(tail - open); recv_page(r, open, (m + 3) & ~3ull); apply(sums, r, m);
}
int main(int argc, char **argv) {
  uint64_t n = strtoull(argv[1], 0, 10), mode = strtoull(argv[2], 0, 10), nch = 4;
  rec_t **tail = malloc(nch * sizeof *tail);
  dir = calloc(nch, sizeof *dir);
  for (uint64_t k = 0; k < nch; k++) dir[k].page = malloc((n / PAGE_RECS + 4) * sizeof(rec_t *));
  pool = aligned_alloc(PAGE_RECS << 4, ((n / PAGE_RECS) + 2 * nch + 4) * (PAGE_RECS << 4)); pool_used = 0;
  rec_t *rbuf = aligned_alloc(64, (n + 2 * PAGE_RECS) * sizeof(rec_t));
  uint32_t *sums = calloc(1 << 16, 4);
  volatile uint64_t sink;
  if (mode == 1) {
    for (uint64_t k = 0; k < nch; k++) { tail[k] = pool + pool_used; pool_used += PAGE_RECS; }
    rec_t *open[4]; walk_scatter(tail, n, 7);
    for (uint64_t k = 0; k < nch; k++) open[k] = (rec_t *)((uintptr_t)tail[k] & ~(uintptr_t)PAGE_MASK);
    for (uint64_t k = 0; k < nch; k++) owner_run(&dir[k], open[k], tail[k], rbuf, sums);
  } else if (mode == 2) walk_linear(pool, n, 7);
  else if (mode == 3) apply(sums, pool, walk_linear(pool, n, 7) - pool);
  else walk_nostore(&sink, n, 7);
  sink = sums[5];
  return (int)(sink & 1);
}
