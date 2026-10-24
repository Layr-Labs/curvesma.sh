
#include <stdint.h>
int gap_scan(const uint64_t *idx, uint64_t n, uint64_t limit, uint64_t *out_max)
{
    uint64_t prev = idx[0], mx = 0;
    for (uint64_t i = 1; i < n; i++) {
        uint64_t g = idx[i] - prev;
        prev = idx[i];
        if (g > mx) mx = g;
        if (mx > limit) { *out_max = mx; return 0; }
    }
    *out_max = mx;
    return 1;
}
