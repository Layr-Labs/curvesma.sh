/* driver_align_v15n.c -- extern definitions (separate TU so clang cannot fold) */
typedef unsigned long long u64;
typedef struct { u64 x[4], y[4], slot, len, esc, pad; } rec_t;
#define NREC_MAX 4096
rec_t A[NREC_MAX];
u64 NREC = 0;
u64 NODE_HEAD[65536];
u64 REPORT_STORE[8];
u64 REPORT_REC[4];
const u64 BETA[4] = {0xC1396C28719501EEULL, 0x9CF0497512F58995ULL,
                     0x6E64479EAC3434E9ULL, 0x7AE96A2B657C0710ULL};
