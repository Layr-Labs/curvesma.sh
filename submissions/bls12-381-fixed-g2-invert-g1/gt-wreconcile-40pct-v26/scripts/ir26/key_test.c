/* key_test.c -- prints state_key48(v) for fixed pseudo-random twelve-tuples so
 * the Python model can check the sort+fold bit for bit. */
#include <stdio.h>
#include <stdint.h>
typedef uint64_t u64;
void state_key48(const u64 v[12], u64 out[6]);

int main(void) {
    u64 s = 0x9E3779B97F4A7C15ULL;
    for (int t = 0; t < 32; t++) {
        u64 v[12], o[6];
        for (int i = 0; i < 12; i++) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            v[i] = s;
        }
        state_key48(v, o);
        printf("%d", t);
        for (int i = 0; i < 6; i++) printf(" %016llx", (unsigned long long)o[i]);
        printf("\n");
    }
    return 0;
}
