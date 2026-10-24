/* lookup_count.c -- the per-step table lookup the ledger cites, counted by
 * disassembly.  Open addressing over the Rg-slot pool: hash the six key words,
 * mask to the table, probe two slots, compare six 64-bit words on a hit.
 * Returns the slot index or -1.  This is the symbol priced as LOOKUP in
 * ledger25.json (scripts/ir25/lookup_count.c + evidence/ir25.json).
 */
#include <stdint.h>
typedef uint64_t u64;

static inline u64 mix64(u64 x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

int lookup_probe(const u64 key[6], const u64 *tab, u64 mask) {
    u64 h = 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < 6; i++)
        h = mix64(h ^ key[i]);
    u64 s0 = h & mask;
    const u64 *e0 = tab + (s0 << 3);
    if (e0[0] == key[0] && e0[1] == key[1] && e0[2] == key[2] &&
        e0[3] == key[3] && e0[4] == key[4] && e0[5] == key[5])
        return (int)s0;
    u64 s1 = (s0 + 1) & mask;
    const u64 *e1 = tab + (s1 << 3);
    if (e1[0] == key[0] && e1[1] == key[1] && e1[2] == key[2] &&
        e1[3] == key[3] && e1[4] == key[4] && e1[5] == key[5])
        return (int)s1;
    return -1;
}
