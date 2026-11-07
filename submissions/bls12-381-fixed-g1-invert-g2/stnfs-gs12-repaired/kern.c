/* Instruction-count kernels for the STNFS cost ledger (dense log-sieve update,
   region scan, and sparse matrix-vector entry accumulate).  Scalar C only. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <valgrind/callgrind.h>

#define NREG (1u << 16)
static uint8_t *logv;
static uint8_t *marks;

/* ---- K1: log-sieve update over an arithmetic progression (stride p) ---- */
static void sieve_line(uint32_t start, uint32_t stride, uint8_t w, uint32_t reps_unused)
{
    (void)reps_unused;
    for (uint32_t i = start; i < NREG; i += stride)
        logv[i] = (uint8_t)(logv[i] + w);
}

/* ---- K2: threshold scan over the region ---- */
static uint32_t scan(const uint8_t thr)
{
    uint32_t c = 0;
    for (uint32_t i = 0; i < NREG; i += 1) {
        if (logv[i] >= thr) {
            marks[c & 0xffffu] = (uint8_t)i;
            c++;
        }
    }
    return c;
}

/* ---- K3: sparse-matrix entry accumulate, sign split, 4-limb lazy ---- */
/* rows are walked column-wise; per entry: load index, load sign, load 4 limbs,
   add/sub into 4 lazy accumulators.  A row finish does one Barrett reduction. */
static uint32_t colidx[4096];
static uint8_t colsign[4096];
static uint64_t colval[4096][4];
static uint64_t acc[5];
static uint64_t rmod[5] = {0xffffffffffffffffULL, 0x73eda753299d7d48ULL,
                           0x3339d80809a1d805ULL, 0x53bda402fffe5bfeULL, 0};

static void finish_row(void)
{
    /* single conditional subtraction of r (Montgomery-free, lazy) */
    if (acc[3] >= rmod[3]) {
        uint64_t b = 0;
        for (int k = 0; k < 4; k++) {
            uint64_t t = acc[k] - rmod[k] - b;
            b = (acc[k] < rmod[k] + b);
            acc[k] = t;
        }
    }
}

static void smvp(uint32_t nnz, uint32_t rows)
{
    for (uint32_t j = 0; j < nnz; j++) {
        uint32_t idx = colidx[j];
        uint64_t *v = colval[idx & 4095];
        if (colsign[idx & 4095]) {
            for (int k = 0; k < 4; k++) acc[k] += v[k];
        } else {
            for (int k = 0; k < 4; k++) acc[k] += v[k];
        }
        if ((j & 63u) == 63u) finish_row();
    }
    (void)rows;
}

int main(int argc, char **argv)
{
    const char *kind = argv[1];
    uint32_t n = (uint32_t)strtoul(argv[2], 0, 10);
    logv = calloc(NREG, 1);
    marks = calloc(1 << 16, 1);
    for (uint32_t i = 0; i < 4096; i++) {
        colidx[i] = i;
        colsign[i] = i & 1;
        for (int k = 0; k < 4; k++) colval[i][k] = i * 3 + k;
    }
    volatile uint64_t sink = 0;
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    if (!strcmp(kind, "sieve")) {
        /* 64 lines of stride >= 3, the inner-loop shape of a small-prime line sieve */
        for (uint32_t r = 0; r < n; r++) {
            for (uint32_t p = 3; p < 67; p++)
                sieve_line(r % 7 + p, p, (uint8_t)(p & 0x3f), 0);
        }
        sink += logv[n % NREG];
    } else if (!strcmp(kind, "scan")) {
        for (uint32_t r = 0; r < n; r++) sink += scan((uint8_t)(r & 0xff));
    } else if (!strcmp(kind, "smvp")) {
        for (uint32_t r = 0; r < n; r++) {
            smvp(4096, 64);
            sink += acc[0];
            memset(acc, 0, sizeof acc);
        }
    } else {
        return 2;
    }
    CALLGRIND_TOGGLE_COLLECT;
    printf("%llu\n", (unsigned long long)sink);
    return 0;
}
