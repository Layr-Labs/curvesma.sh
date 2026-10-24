/* v16 executed radix-sort price: LSD radix sort of 2^20 records keyed on the full
 * 256-bit class x, 4 x 64-bit limbs, 16 passes of 16 bits (16*16 = 256).
 * Each record is 48 bytes (32B key + 16B identity payload).  The program is
 * compiled from this shipped source and its inner loops are disassembled; the
 * ledger row is the static instruction count of the count and scatter loops times
 * the trip counts of the executed run, so it is an implementation bound, not a
 * constant.
 *
 * Build/run: cc -O2 -o radix_sort_v16 radix_sort_v16.c && ./radix_sort_v16 20
 * Disassemble: objdump -d radix_sort_v16
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint64_t k[4]; uint64_t t; uint64_t s; } rec_t;   /* 48 B */

static uint32_t *cnt;

int main(int argc, char **argv) {
    int logn = argc > 1 ? atoi(argv[1]) : 20;
    size_t n = (size_t)1 << logn;
    rec_t *a = malloc(n * sizeof(rec_t));
    rec_t *b = malloc(n * sizeof(rec_t));
    cnt = calloc(65536, sizeof(uint32_t));
    if (!a || !b || !cnt) return 2;
    uint64_t x = 0x9E3779B97F4A7C15ULL;
    for (size_t i = 0; i < n; i++) {                     /* deterministic fill */
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        a[i].k[0] = x; a[i].k[1] = x * 0xBF58476D1CE4E5B9ULL;
        a[i].k[2] = x * 0x94D049BB133111EBULL; a[i].k[3] = x ^ (uint64_t)i;
        a[i].t = i & 0xFFFFF; a[i].s = i;
    }
    uint64_t elems = 0, moved = 0;
    for (int pass = 0; pass < 16; pass++) {              /* LSD, 16 bits per pass */
        memset(cnt, 0, 65536 * sizeof(uint32_t));
        int limb = pass >> 2, shift = (pass & 3) * 16;
        for (size_t i = 0; i < n; i++) {                 /* count loop */
            cnt[(a[i].k[limb] >> shift) & 0xFFFF]++;
            elems++;
        }
        uint32_t sum = 0;
        for (int d = 0; d < 65536; d++) { uint32_t c = cnt[d]; cnt[d] = sum; sum += c; }
        for (size_t i = 0; i < n; i++) {                 /* scatter loop */
            unsigned d = (a[i].k[limb] >> shift) & 0xFFFF;
            b[cnt[d]++] = a[i];
            moved++;
        }
        rec_t *tmp = a; a = b; b = tmp;
    }
    int sorted = 1;
    for (size_t i = 1; i < n; i++) {
        for (int l = 3; l >= 0; l--) {                   /* k[3] is most significant */
            if (a[i-1].k[l] < a[i].k[l]) break;
            if (a[i-1].k[l] > a[i].k[l]) { sorted = 0; break; }
        }
        if (!sorted) break;
    }
    printf("sorted=%d\n", sorted);
    /* equal-key collision detection pass: emit adjacent equal-key witnesses */
    uint64_t collisions = 0;
    for (size_t i = 1; i < n; i++) {
        if (a[i-1].k[0] == a[i].k[0] && a[i-1].k[1] == a[i].k[1] &&
            a[i-1].k[2] == a[i].k[2] && a[i-1].k[3] == a[i].k[3]) collisions++;
    }
    printf("collisions=%llu\n", (unsigned long long)collisions);
    printf("n=%zu elems_counted=%llu records_scattered=%llu bytes_moved=%llu\n",
           n, (unsigned long long)elems, (unsigned long long)moved,
           (unsigned long long)(moved * sizeof(rec_t)));
    printf("sample=%016llx %016llx\n", (unsigned long long)a[0].k[0],
           (unsigned long long)a[n - 1].k[0]);
    return 0;
}
