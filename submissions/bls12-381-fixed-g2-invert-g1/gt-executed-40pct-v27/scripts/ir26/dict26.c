/* dict26.c -- THE SHIPPED LOOKUP DICTIONARY, implemented for real (v26).
 *
 * v25 priced `lookup_probe`: a 64-bit hash, two probes, six-word compare, no
 * insertion, no overflow policy, no load factor, and a key builder that was not
 * the one Group.state uses.  The v25 verdict (d2be26ee flaw 2) rejected exactly
 * that gap.  This file is the structure the ledger now cites, and every symbol
 * below is disassembly-priced (scripts/ir26/count_ir26.py -> evidence/ir26.json):
 *
 *   state_key48   sort the twelve field values the shipped Group.state sorts,
 *                 then fold them into the 48-byte (six-word) pool key.
 *   dict_lookup   open addressing, FOUR probes at a fixed stride, six-word
 *                 compare on an occupied slot, stop at the first empty slot.
 *   dict_insert   first empty of the four probes; if all four are occupied,
 *                 EVICT the occupant with the smallest leading key word (ties
 *                 break to the lowest probe index) and count the eviction.
 *   dict_lf_ok    the load-factor guarantee the caller enforces: 8n <= 5*cap.
 *
 * Layout: one entry is EIGHT 64-bit words -- six key words, one record
 * identifier (the index of the A,D,y record the hit resolves to, which is what
 * v25's int-returning probe never carried), and a tag word (0 = empty).  A slot
 * is 64 bytes; the index holds S slots at the 8n <= 5cap load factor.
 */
#include <stdint.h>

typedef uint64_t u64;

#define TAG_EMPTY 0ULL
#define TAG_FULL  1ULL
#define PROBES    4
#define ENTRY     8
#define STRIDE    0x9E3779B97F4A7C15ULL

static inline u64 mix64(u64 x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

/* state_key48: the twelve field values, insertion-sorted, folded to six words.
 * Sorting is what makes the key invariant to the rotation kappa; the fold keeps
 * the record at 48 B (the y itself is stored in the record, not the key). */
void state_key48(const u64 v[12], u64 out[6]) {
    u64 a[12];
    for (int i = 0; i < 12; i++) a[i] = v[i];
    for (int i = 1; i < 12; i++) {
        u64 t = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > t) { a[j + 1] = a[j]; j--; }
        a[j + 1] = t;
    }
    for (int j = 0; j < 6; j++)
        out[j] = mix64(a[j] + STRIDE * a[j + 6]);
}

static inline u64 dict_hash(const u64 key[6]) {
    u64 h = STRIDE;
    for (int i = 0; i < 6; i++) h = mix64(h ^ key[i]);
    return h;
}

static inline int dict_cmp(const u64 *e, const u64 key[6]) {
    for (int i = 0; i < 6; i++)
        if (e[i] != key[i]) return 0;
    return 1;
}

/* returns the slot index of key, or -1 if absent (stop at first empty slot). */
int dict_lookup(const u64 key[6], u64 *tab, u64 mask) {
    u64 s = dict_hash(key) & mask;
    for (int i = 0; i < PROBES; i++) {
        u64 *e = tab + (((s + (u64)i * STRIDE) & mask) * ENTRY);
        if (e[7] == TAG_EMPTY) return -1;
        if (dict_cmp(e, key)) return (int)(((s + (u64)i * STRIDE) & mask));
    }
    return -1;
}

/* insert key; returns slot index written, *evicted bumped when all four probes
 * were occupied and the smallest-leading-word occupant was displaced. */
int dict_insert(const u64 key[6], u64 *tab, u64 mask, u64 rec, u64 *evicted) {
    u64 s = dict_hash(key) & mask;
    u64 slot[PROBES];
    u64 *ent[PROBES];
    for (int i = 0; i < PROBES; i++) {
        slot[i] = (s + (u64)i * STRIDE) & mask;
        ent[i] = tab + slot[i] * ENTRY;
        if (ent[i][7] == TAG_EMPTY) {
            for (int w = 0; w < 6; w++) ent[i][w] = key[w];
            ent[i][6] = rec;                 /* resolves the hit to its record */
            ent[i][7] = TAG_FULL;
            return (int)slot[i];
        }
        if (dict_cmp(ent[i], key)) { ent[i][6] = rec; return (int)slot[i]; }
    }
    int v = 0;
    for (int i = 1; i < PROBES; i++)
        if (ent[i][0] < ent[v][0]) v = i;                  /* ties -> lowest i */
    for (int w = 0; w < 6; w++) ent[v][w] = key[w];
    ent[v][6] = rec;
    ent[v][7] = TAG_FULL;
    *evicted += 1;
    return (int)slot[v];
}

/* the record a hit resolves to; only meaningful for an occupied slot. */
u64 dict_record(const u64 *tab, int slot) { return tab[(u64)slot * ENTRY + 6]; }

/* the load-factor guarantee: 8n <= 5cap, i.e. at most 62.5% occupied. */
int dict_lf_ok(u64 n, u64 cap) {
    return 8ULL * n <= 5ULL * cap;
}
