/* Batched sigma-walk step for Callgrind counts; every instruction is charged. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "p256_field.h"
#include "p256_sub.h"
#include "p256_fused.h"

/* the hot backward loop's field work, 576 instructions */
#define BLOOP_FIELD \
"ldp x7,x8,[x23,#32]\nmul x0,x19,x7\numulh x15,x19,x7\nmul x2,x20,x7\numulh x16,x20,x7\nmul x18,x21,x7\numulh x17,x21,x7\nmul x30,x22,x7\numulh x1,x22,x7\nadds x2,x2,x15\nadcs x18,x18,x16\nlsr x16,x0,#32\n" \
"adcs x30,x30,x17\nadc x13,xzr,x1\nsubs x17,x0,x0,lsl #32\nsbc x1,x0,x16\nadds x0,x2,x0,lsl #32\nmul x15,x19,x8\nadcs x2,x18,x16\nmul x16,x20,x8\nadcs x18,x30,x17\nmul x17,x21,x8\nadcs x30,x13,x1\nmul x1,x22,x8\n" \
"adc x13,xzr,xzr\nadds x0,x0,x15\numulh x15,x19,x8\nadcs x2,x2,x16\numulh x16,x20,x8\nadcs x18,x18,x17\numulh x17,x21,x8\nadcs x30,x30,x1\numulh x1,x22,x8\nadc x13,x13,xzr\nldp x7,x8,[x23,#48]\nadds x2,x2,x15\n" \
"adcs x18,x18,x16\nlsr x16,x0,#32\nadcs x30,x30,x17\nadcs x13,x13,x1\nadc x14,xzr,xzr\nsubs x17,x0,x0,lsl #32\nsbc x1,x0,x16\nadds x0,x2,x0,lsl #32\nmul x15,x19,x7\nadcs x2,x18,x16\nmul x16,x20,x7\nadcs x18,x30,x17\n" \
"mul x17,x21,x7\nadcs x30,x13,x1\nmul x1,x22,x7\nadc x13,x14,xzr\nadds x0,x0,x15\numulh x15,x19,x7\nadcs x2,x2,x16\numulh x16,x20,x7\nadcs x18,x18,x17\numulh x17,x21,x7\nadcs x30,x30,x1\numulh x1,x22,x7\n" \
"adc x13,x13,xzr\nadds x2,x2,x15\nadcs x18,x18,x16\nlsr x16,x0,#32\nadcs x30,x30,x17\nadcs x13,x13,x1\nadc x14,xzr,xzr\nsubs x17,x0,x0,lsl #32\nsbc x1,x0,x16\nadds x0,x2,x0,lsl #32\nmul x15,x19,x8\nadcs x2,x18,x16\n" \
"mul x16,x20,x8\nadcs x18,x30,x17\nmul x17,x21,x8\nadcs x30,x13,x1\nmul x1,x22,x8\nadc x13,x14,xzr\nadds x0,x0,x15\numulh x15,x19,x8\nadcs x2,x2,x16\numulh x16,x20,x8\nadcs x18,x18,x17\numulh x17,x21,x8\n" \
"adcs x30,x30,x1\numulh x1,x22,x8\nadc x13,x13,xzr\nadds x2,x2,x15\nadcs x18,x18,x16\nlsr x16,x0,#32\nadcs x30,x30,x17\nadcs x13,x13,x1\nadc x14,xzr,xzr\nsubs x17,x0,x0,lsl #32\nsbc x1,x0,x16\nadds x0,x2,x0,lsl #32\n" \
"adcs x2,x18,x16\nadcs x18,x30,x17\nadcs x30,x13,x1\nadc x13,x14,xzr\nneg x7,x13\nand x8,x7,#-0x100000000\nadds x0,x0,x13\nand x13,x7,#0xfffffffe\nadcs x2,x2,x8\nadcs x18,x18,x7\nadc x30,x30,x13\nldp x5,x6,[x23,#16]\n" \
"ldp x3,x4,[x23],#-64\nmul x9,x3,x19\numulh x15,x3,x19\nmul x10,x4,x19\numulh x16,x4,x19\nmul x11,x5,x19\numulh x17,x5,x19\nmul x12,x6,x19\numulh x1,x6,x19\nadds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\n" \
"adcs x12,x12,x17\nadc x13,xzr,x1\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nmul x15,x3,x20\nadcs x10,x11,x16\nmul x16,x4,x20\nadcs x11,x12,x17\nmul x17,x5,x20\nadcs x12,x13,x1\nmul x1,x6,x20\n" \
"umulh x14,x6,x20\nadc x13,x14,xzr\nadds x9,x9,x15\numulh x15,x3,x20\nadcs x10,x10,x16\numulh x16,x4,x20\nadcs x11,x11,x17\numulh x17,x5,x20\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x10,x10,x15\nadcs x11,x11,x16\n" \
"lsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nmul x15,x3,x21\nadcs x10,x11,x16\nmul x16,x4,x21\nadcs x11,x12,x17\nmul x17,x5,x21\nadcs x12,x13,x1\n" \
"mul x1,x6,x21\numulh x14,x6,x21\nadc x13,x14,xzr\nadds x9,x9,x15\numulh x15,x3,x21\nadcs x10,x10,x16\numulh x16,x4,x21\nadcs x11,x11,x17\numulh x17,x5,x21\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x10,x10,x15\n" \
"adcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nmul x15,x3,x22\nadcs x10,x11,x16\nmul x16,x4,x22\nadcs x11,x12,x17\nmul x17,x5,x22\n" \
"adcs x12,x13,x1\nmul x1,x6,x22\numulh x14,x6,x22\nadc x13,x14,xzr\nadds x9,x9,x15\numulh x15,x3,x22\nadcs x10,x10,x16\numulh x16,x4,x22\nadcs x11,x11,x17\numulh x17,x5,x22\nadcs x12,x12,x1\nadc x13,x13,xzr\n" \
"adds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nadcs x10,x11,x16\nadcs x11,x12,x17\nadcs x12,x13,x1\nadc x13,xzr,xzr\n" \
"neg x7,x13\nand x8,x7,#-0x100000000\nadds x19,x9,x13\nand x13,x7,#0xfffffffe\nadcs x20,x10,x8\nadcs x21,x11,x7\nadc x22,x12,x13\nldr x7,[x26,#8]\nldp x3,x4,[x7,#32]\nldp x5,x6,[x7,#48]\nldp x9,x10,[x24]\nldp x11,x12,[x24,#16]\n" \
"subs x3,x3,x9\nsbcs x4,x4,x10\nsbcs x5,x5,x11\nsbcs x6,x6,x12\nsbc x13,xzr,xzr\nlsr x14,x13,#32\nand x15,x13,#-0xffffffff\nadds x3,x3,x13\nadcs x4,x4,x14\nadcs x5,x5,xzr\nadc x6,x6,x15\nmul x9,x3,x0\n" \
"umulh x15,x3,x0\nmul x10,x4,x0\numulh x16,x4,x0\nmul x11,x5,x0\numulh x17,x5,x0\nmul x12,x6,x0\numulh x1,x6,x0\nadds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,xzr,x1\n" \
"subs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nmul x15,x3,x2\nadcs x10,x11,x16\nmul x16,x4,x2\nadcs x11,x12,x17\nmul x17,x5,x2\nadcs x12,x13,x1\nmul x1,x6,x2\numulh x14,x6,x2\nadc x13,x14,xzr\n" \
"adds x9,x9,x15\numulh x15,x3,x2\nadcs x10,x10,x16\numulh x16,x4,x2\nadcs x11,x11,x17\numulh x17,x5,x2\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\nadcs x12,x12,x17\n" \
"adc x13,x13,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nmul x15,x3,x18\nadcs x10,x11,x16\nmul x16,x4,x18\nadcs x11,x12,x17\nmul x17,x5,x18\nadcs x12,x13,x1\nmul x1,x6,x18\numulh x14,x6,x18\n" \
"adc x13,x14,xzr\nadds x9,x9,x15\numulh x15,x3,x18\nadcs x10,x10,x16\numulh x16,x4,x18\nadcs x11,x11,x17\numulh x17,x5,x18\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x10,x10,x15\nadcs x11,x11,x16\nlsr x16,x9,#32\n" \
"adcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nmul x15,x3,x30\nadcs x10,x11,x16\nmul x16,x4,x30\nadcs x11,x12,x17\nmul x17,x5,x30\nadcs x12,x13,x1\nmul x1,x6,x30\n" \
"umulh x14,x6,x30\nadc x13,x14,xzr\nadds x9,x9,x15\numulh x15,x3,x30\nadcs x10,x10,x16\numulh x16,x4,x30\nadcs x11,x11,x17\numulh x17,x5,x30\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x10,x10,x15\nadcs x11,x11,x16\n" \
"lsr x16,x9,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x9,x9,lsl #32\nsbc x1,x9,x16\nadds x9,x10,x9,lsl #32\nadcs x10,x11,x16\nadcs x11,x12,x17\nadcs x12,x13,x1\nadc x13,xzr,xzr\nadds x15,x9,#1\nsbcs x16,x10,x28\n" \
"sbcs x17,x11,xzr\nadcs x1,x12,x29\nsbcs xzr,x13,xzr\ncsel x0,x9,x15,lo\ncsel x2,x10,x16,lo\ncsel x18,x11,x17,lo\ncsel x30,x12,x1,lo\nmul x4,x2,x0\numulh x16,x2,x0\nmul x5,x18,x0\numulh x17,x18,x0\nmul x6,x30,x0\n" \
"umulh x11,x30,x0\nadds x5,x5,x16\nmul x15,x18,x2\numulh x16,x18,x2\nadcs x6,x6,x17\nmul x17,x30,x2\numulh x1,x30,x2\nadc x11,x11,xzr\nmul x12,x30,x18\numulh x13,x30,x18\nadds x16,x16,x17\nmul x3,x0,x0\n" \
"adc x17,x1,xzr\nadds x6,x6,x15\numulh x10,x0,x0\nadcs x11,x11,x16\nmul x16,x2,x2\nadcs x12,x12,x17\numulh x15,x2,x2\nadc x13,x13,xzr\nadds x4,x4,x4\nmul x17,x18,x18\nadcs x5,x5,x5\numulh x8,x18,x18\n" \
"adcs x6,x6,x6\nmul x1,x30,x30\nadcs x11,x11,x11\numulh x9,x30,x30\nadcs x12,x12,x12\nadcs x13,x13,x13\nadc x14,xzr,xzr\nadds x4,x4,x10\nadcs x5,x5,x16\nadcs x6,x6,x15\nadcs x11,x11,x17\nadcs x12,x12,x8\n" \
"adcs x13,x13,x1\nlsr x16,x3,#32\nadc x14,x14,x9\nsubs x17,x3,x3,lsl #32\nsbc x1,x3,x16\nadds x3,x4,x3,lsl #32\nadcs x4,x5,x16\nadcs x5,x6,x17\nlsr x16,x3,#32\nadc x6,x1,xzr\nsubs x17,x3,x3,lsl #32\nsbc x1,x3,x16\n" \
"adds x3,x4,x3,lsl #32\nadcs x4,x5,x16\nadcs x5,x6,x17\nlsr x16,x3,#32\nadc x6,x1,xzr\nsubs x17,x3,x3,lsl #32\nsbc x1,x3,x16\nadds x3,x4,x3,lsl #32\nadcs x4,x5,x16\nadcs x5,x6,x17\nlsr x16,x3,#32\nadc x6,x1,xzr\n" \
"subs x17,x3,x3,lsl #32\nsbc x1,x3,x16\nadds x3,x4,x3,lsl #32\nadcs x4,x5,x16\nadcs x5,x6,x17\nadc x6,x1,xzr\nadds x3,x3,x11\nadcs x4,x4,x12\nadcs x5,x5,x13\nadcs x6,x6,x14\nadc x11,xzr,xzr\n" \
"adds x15,x3,#1\nsbcs x16,x4,x28\nsbcs x17,x5,xzr\nadcs x1,x6,x29\nsbcs xzr,x11,xzr\ncsel x3,x3,x15,lo\ncsel x4,x4,x16,lo\ncsel x5,x5,x17,lo\ncsel x6,x6,x1,lo\nldp x8,x9,[x24,#-32]\nldp x10,x11,[x24,#-16]\n" \
"subs x3,x3,x8\nsbcs x4,x4,x9\nsbcs x5,x5,x10\nsbcs x6,x6,x11\nsbc x15,xzr,xzr\nlsr x16,x15,#32\nand x17,x15,#-0xffffffff\nadds x3,x3,x15\nadcs x4,x4,x16\nadcs x5,x5,xzr\nadc x6,x6,x17\nldp x12,x13,[x7]\n" \
"ldp x14,x1,[x7,#16]\nsubs x3,x3,x12\nsbcs x4,x4,x13\nsbcs x5,x5,x14\nsbcs x6,x6,x1\nsbc x15,xzr,xzr\nlsr x16,x15,#32\nand x17,x15,#-0xffffffff\nadds x3,x3,x15\nadcs x4,x4,x16\nadcs x5,x5,xzr\nadc x6,x6,x17\n" \
"stp x3,x4,[x25,#-48]\nstp x5,x6,[x25,#-32]\nsubs x8,x8,x3\nsbcs x9,x9,x4\nsbcs x10,x10,x5\nsbcs x11,x11,x6\nsbc x12,xzr,xzr\nlsr x13,x12,#32\nand x14,x12,#-0xffffffff\nadds x8,x8,x12\nadcs x9,x9,x13\nadcs x10,x10,xzr\n" \
"adc x11,x11,x14\nmul x5,x8,x0\numulh x15,x8,x0\nmul x6,x9,x0\numulh x16,x9,x0\nmul x7,x10,x0\numulh x17,x10,x0\nmul x12,x11,x0\numulh x1,x11,x0\nadds x6,x6,x15\nadcs x7,x7,x16\nlsr x16,x5,#32\n" \
"adcs x12,x12,x17\nadc x13,xzr,x1\nsubs x17,x5,x5,lsl #32\nsbc x1,x5,x16\nadds x5,x6,x5,lsl #32\nmul x15,x8,x2\nadcs x6,x7,x16\nmul x16,x9,x2\nadcs x7,x12,x17\nmul x17,x10,x2\nadcs x12,x13,x1\nmul x1,x11,x2\n" \
"umulh x14,x11,x2\nadc x13,x14,xzr\nadds x5,x5,x15\numulh x15,x8,x2\nadcs x6,x6,x16\numulh x16,x9,x2\nadcs x7,x7,x17\numulh x17,x10,x2\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x6,x6,x15\nadcs x7,x7,x16\n" \
"lsr x16,x5,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x5,x5,lsl #32\nsbc x1,x5,x16\nadds x5,x6,x5,lsl #32\nmul x15,x8,x18\nadcs x6,x7,x16\nmul x16,x9,x18\nadcs x7,x12,x17\nmul x17,x10,x18\nadcs x12,x13,x1\n" \
"mul x1,x11,x18\numulh x14,x11,x18\nadc x13,x14,xzr\nadds x5,x5,x15\numulh x15,x8,x18\nadcs x6,x6,x16\numulh x16,x9,x18\nadcs x7,x7,x17\numulh x17,x10,x18\nadcs x12,x12,x1\nadc x13,x13,xzr\nadds x6,x6,x15\n" \
"adcs x7,x7,x16\nlsr x16,x5,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x5,x5,lsl #32\nsbc x1,x5,x16\nadds x5,x6,x5,lsl #32\nmul x15,x8,x30\nadcs x6,x7,x16\nmul x16,x9,x30\nadcs x7,x12,x17\nmul x17,x10,x30\n" \
"adcs x12,x13,x1\nmul x1,x11,x30\numulh x14,x11,x30\nadc x13,x14,xzr\nadds x5,x5,x15\numulh x15,x8,x30\nadcs x6,x6,x16\numulh x16,x9,x30\nadcs x7,x7,x17\numulh x17,x10,x30\nadcs x12,x12,x1\nadc x13,x13,xzr\n" \
"adds x6,x6,x15\nadcs x7,x7,x16\nlsr x16,x5,#32\nadcs x12,x12,x17\nadc x13,x13,xzr\nsubs x17,x5,x5,lsl #32\nsbc x1,x5,x16\nadds x5,x6,x5,lsl #32\nadcs x6,x7,x16\nadcs x7,x12,x17\nadcs x12,x13,x1\nadc x13,xzr,xzr\n" \
"adds x15,x5,#1\nsbcs x16,x6,x28\nsbcs x17,x7,xzr\nadcs x1,x12,x29\nsbcs xzr,x13,xzr\ncsel x0,x5,x15,lo\ncsel x2,x6,x16,lo\ncsel x18,x7,x17,lo\ncsel x30,x12,x1,lo\nldp x10,x11,[x24,#16]\nldp x8,x9,[x24],#-64\nsubs x0,x0,x8\n" \
"sbcs x2,x2,x9\nsbcs x18,x18,x10\nsbcs x30,x30,x11\nsbc x12,xzr,xzr\nlsr x13,x12,#32\nand x14,x12,#-0xffffffff\nadds x0,x0,x12\nadcs x2,x2,x13\nadcs x18,x18,xzr\nadc x30,x30,x14\nstp x0,x2,[x25,#-16]\nstp x18,x30,[x25],#-64\n" \
""

#ifndef ADDER_BITS
#define ADDER_BITS 16
#endif
#define NA (1 << ADDER_BITS)
#ifndef B
#define B 65536
#endif
#define ANCHOR_PERIOD 64
#define MAX_CYCLE 64
#ifndef DP_BITS
#define DP_BITS 60
#endif
#define SPARES 1024
/* A trail is abandoned after more than 2^DP_BITS anchor boundaries. */
#ifndef CAP_PERIODS
#define CAP_PERIODS (1ULL << DP_BITS)
#endif
/* The filters compare the low key bits; a test build can narrow them. */
#ifndef FILTER_MASK
#define FILTER_MASK (~0ULL)
#endif
#if DP_BITS + ADDER_BITS <= 64
#define IS_DP(kl, p) ((((kl) >> ADDER_BITS) & ((1ULL << DP_BITS) - 1)) == 0)
#define IS_DP1 IS_DP
#else
#define IS_DP(kl, p) (((kl) >> ADDER_BITS) == 0 && \
                      ((p)->x.v[1] & ((1ULL << (DP_BITS + ADDER_BITS - 64)) - 1)) == 0)
/* IS_DP as one compare. */
#define IS_DP1(kl, p) ((((kl) >> ADDER_BITS) | \
                        ((p)->x.v[1] & ((1ULL << (DP_BITS + ADDER_BITS - 64)) - 1))) == 0)
#endif

typedef struct { fe x, y; } pt; /* affine, Montgomery form */

static const fe P = {{0xffffffffffffffffULL, 0x00000000ffffffffULL, 0, 0xffffffff00000001ULL}};
static const fe ZERO = {{0}};
/* p - 2 */
static const uint64_t PM2[4] = {0xfffffffffffffffdULL, 0x00000000ffffffffULL, 0, 0xffffffff00000001ULL};

static pt T[NA]; /* adders R_i as read */
static pt T2[NA][2]; /* T2[i][s] = (X_i, (-1)^s Y_i) */
static pt WA[B], WB[B];
static pt *W = WA, *WN = WB; /* current and next states, swapped every batch step */
/* DP[w] = {D[w], PRE[w]}, PRE[w] = D[0]*...*D[w-1] mod p, below 2^256; DP is aligned to its size, so the hot
 * loop ends on its pointer's low bits */
typedef struct { fe d, pre; } dpre;
static char dp_mem[2 * B * sizeof(dpre)];
#define DP_SIZE ((uintptr_t)B * sizeof(dpre))
typedef struct { uint64_t prev, anchor, cur; const pt *sel; } lane; /* two steps back, anchor, current, adder */
static lane L[B];
static uint64_t tper[B];
static pt SPARE[SPARES];
#define RING_K 4
static uint64_t ring[B][RING_K];
static int ring_n[B], ring_pos[B];
static uint64_t loops;
static uint64_t two_cycles, long_cycles, dps, exceptional, unbatched;
static uint64_t abandoned, next_spare, ahead_dps;

static inline int feq(const fe *a, const fe *b) { return memcmp(a->v, b->v, sizeof a->v) == 0; }
static inline int fis_zero(const fe *a) { return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0; }
static inline int fcmp(const fe *a, const fe *b) {
    for (int i = 3; i >= 0; i--) {
        if (a->v[i] != b->v[i]) return a->v[i] < b->v[i] ? -1 : 1;
    }
    return 0;
}

/* r = a^(p-2) = a^-1 (Fermat), left-to-right square and multiply. */
static __attribute__((noinline)) void finv(fe *r, const fe *a) {
    fe acc = *a;
    for (int i = 254; i >= 0; i--) {
        k_sqr(&acc, &acc);
        if ((PM2[i >> 6] >> (i & 63)) & 1) k_mul(&acc, &acc, a);
    }
    *r = acc;
}

static inline uint64_t key_lo(const pt *p) { return p->x.v[0]; }

/* sigma(R_i) from W. */
static inline const pt *select_adder(const pt *p) {
    return &T2[p->x.v[0] & (NA - 1)][p->y.v[0] & 1];
}

/* select_adder; the empty asm keeps the mask one instruction at every A. */
static inline const pt *select_adder1(const pt *p) {
    uint64_t i = p->x.v[0] & (NA - 1);
    __asm__("" : "+r"(i));
    return &T2[i][p->y.v[0] & 1];
}

/* Distinguished-point test for a look-ahead point. */
static int point_is_dp(const pt *p) {
    uint64_t k = key_lo(p);
    return IS_DP(k, p);
}

/* Folds each distinguished point's class key into a digest the replica checks: it shows which point every trail records. */
static uint64_t dp_digest;
static void record_dp(const pt *p) {
    for (int i = 0; i < 4; i++) dp_digest = dp_digest * 1099511628211ULL + p->x.v[i];
}

/* Start a new trail in walk w: in production the finished trail's key and coins go to the
 * table and fresh coins give the new start point (charged separately); this harness takes
 * the spare points of the inputs. */
static __attribute__((noinline)) void new_trail(int w, pt *nw, uint64_t *kl, const pt **sel) {
    *nw = SPARE[next_spare++ % SPARES];
    *kl = key_lo(nw);
    *sel = select_adder(nw);
    L[w].prev = L[w].anchor = ~0ULL;
    tper[w] = 0;
    ring_n[w] = ring_pos[w] = 0;
}

/* Diagnostic (loops); no control flow depends on it. */
static inline void esc_seen(int w, uint64_t kl) {
    int n = ring_n[w] < RING_K ? ring_n[w] : RING_K;
    for (int i = 0; i < n; i++)
        if (ring[w][i] == kl) { loops++; return; }
    ring[w][ring_pos[w]] = kl;
    ring_pos[w] = (ring_pos[w] + 1) % RING_K;
    ring_n[w]++;
}

/* Record the distinguished point nw and start a new trail; out of line. */
static __attribute__((noinline)) void dp_restart(int w, pt *nw, uint64_t *kl, const pt **sel) {
    record_dp(nw);
    new_trail(w, nw, kl, sel);
}

static int pt_eq(const pt *a, const pt *b) { return feq(&a->x, &b->x) && feq(&a->y, &b->y); }

/* Unbatched affine addition/doubling with its own inversion (rare paths). */
static int add_single(pt *r, const pt *a, const pt *b) {
    fe num, den, lam, t, x3, y3;
    unbatched++;
    if (feq(&a->x, &b->x)) {
        if (!feq(&a->y, &b->y) || fis_zero(&a->y)) return 0;
        fe three = {{0}};
        k_sqr(&num, &a->x);
        k_add(&t, &num, &num);
        k_add(&num, &num, &t);
        /* a = -3: num = 3x^2 - 3 */
        three.v[0] = 3;
        k_to(&three, &three);
        k_sub(&num, &num, &three);
        k_add(&den, &a->y, &a->y);
    } else {
        k_sub(&num, &b->y, &a->y);
        k_sub(&den, &b->x, &a->x);
    }
    finv(&den, &den);
    k_mul(&lam, &num, &den);
    k_sqr(&x3, &lam);
    k_sub(&x3, &x3, &a->x);
    k_sub(&x3, &x3, &b->x);
    k_sub(&t, &a->x, &x3);
    k_mul(&y3, &lam, &t);
    k_sub(&y3, &y3, &a->y);
    r->x = x3;
    r->y = y3;
    return 1;
}

/* One walk step from the actual point a, unbatched. */
static int f_single(pt *r, const pt *a) { return add_single(r, a, select_adder(a)); }

/* Deterministic exit from a fruitless cycle. */
static int escape(pt *r, const pt *cyc, int len) {
    int m = 0;
    for (int i = 1; i < len; i++)
        if (fcmp(&cyc[i].x, &cyc[m].x) < 0) m = i;
    pt c = cyc[m];
    if (c.y.v[0] & 1) k_sub(&c.y, &ZERO, &c.y);
    return add_single(r, &c, &c);
}

static inline int same_class(const pt *a, const pt *b) { return feq(&a->x, &b->x); }

/* Rare paths, out of line. */
static __attribute__((noinline)) void rare_2cycle(int w, pt *nw, const pt *Wc, uint64_t *kl, const pt **sel,
                                                  int *dp, int *esc) {
    pt g, cyc[2];
    if (f_single(&g, nw) && same_class(&g, &Wc[w])) {
        cyc[0] = Wc[w];
        cyc[1] = *nw;
        if (!escape(nw, cyc, 2)) { exceptional++; abort(); }
        two_cycles++;
        *esc = 1;
        *kl = key_lo(nw);
        *sel = select_adder(nw);
        *dp = IS_DP(*kl, nw);
        if (!*dp) esc_seen(w, *kl);
    }
}

static __attribute__((noinline)) void rare_anchor(int w, pt *nw, uint64_t *kl, const pt **sel, int *dp) {
    pt cyc[MAX_CYCLE + 1], s;
    int len = 1;
    cyc[0] = *nw;
    if (!f_single(&s, nw)) { exceptional++; abort(); }
    while (!same_class(&s, nw) && len <= MAX_CYCLE) {
        if (point_is_dp(&s)) {
            /* s is the trail's distinguished point and its record. */
            *nw = s;
            *kl = key_lo(nw);
            *dp = 1;
            ahead_dps++;
            break;
        }
        cyc[len++] = s;
        if (!f_single(&s, &s)) { exceptional++; abort(); }
    }
    if (!*dp && same_class(&s, nw)) {
        if (!escape(nw, cyc, len)) { exceptional++; abort(); }
        long_cycles++;
        *kl = key_lo(nw);
        *sel = select_adder(nw);
        *dp = IS_DP(*kl, nw);
        if (!*dp) esc_seen(w, *kl);
    }
}

/* A step's tail on a distinguished, matching or over-cap point; out of line. */
static __attribute__((noinline)) void step_rare(int w, pt *nw, const pt *Wc, int anchor_now) {
    uint64_t kl = key_lo(nw);
    const pt *sel = select_adder(nw);
    /* The distinguished-point test comes first, on every point the walk produces
     * (new, look-ahead and escape points), so an escape never skips one. */
    int dp = IS_DP1(kl, nw), esc = 0;
    if (!dp && (kl & FILTER_MASK) == L[w].prev)
        rare_2cycle(w, nw, Wc, &kl, &sel, &dp, &esc);
    /* Not "else": a false predecessor match must not suppress the anchor check, or a walk
     * whose anchor matches the class two steps back could never escape. */
    if (!dp && !esc && (kl & FILTER_MASK) == L[w].anchor)
        rare_anchor(w, nw, &kl, &sel, &dp);
    if (dp) {
        dps++;
        dp_restart(w, nw, &kl, &sel);
    } else if (anchor_now) {
        L[w].anchor = kl & FILTER_MASK;
        if (++tper[w] > CAP_PERIODS) {
            abandoned++;
            new_trail(w, nw, &kl, &sel);
        } else {
            L[w].prev = L[w].cur;
        }
    } else {
        L[w].prev = L[w].cur;
    }
    L[w].cur = kl & FILTER_MASK;
    L[w].sel = sel;
}

/* The hot loop's step tail after BLOOP_FIELD: x3, x4 = x3[0], x3[1], x0 = y3[0], x26 = &L[w].cur. Z is set when the
 * point is distinguished or f matches prev or anchor; then the block exits with x0 = 1 before the lane update. */
#if DP_BITS + ADDER_BITS <= 64
#define BL_DP "tst x3,#%c[dpm]\n"
#define BL_DPM (((1ULL << DP_BITS) - 1) << ADDER_BITS)
#else
#define BL_DP "and x8,x4,#%c[dpm]\norr x8,x8,x3,lsr #%c[ab]\ncmp x8,#0\n"
#define BL_DPM ((1ULL << (DP_BITS + ADDER_BITS - 64)) - 1)
#endif
#if FILTER_MASK == ~0ULL
#define BL_F "x3"
#define BL_FM ""
#else
#define BL_F "x10"
#define BL_FM "and x10,x3,#%c[fm]\n"
#endif
#define BLOOP "str x29,[sp,#-16]!\nmov w28,#-1\nmov w29,#-2\n1:\n" BLOOP_FIELD \
    "ldr x5,[x26]\nldp x6,x7,[x26,#-16]\n" BL_FM BL_DP "ccmp " BL_F ",x6,#4,ne\nccmp " BL_F ",x7,#4,ne\nb.eq 2f\n" \
    "ubfiz x9,x3,#7,#%c[ab]\nbfi x9,x0,#6,#1\nadd x9,x27,x9\nstr x5,[x26,#-16]\nstp " BL_F ",x9,[x26],#-32\n" \
    "tst x23,#%c[dm]\nb.ne 1b\nmov x0,#0\nb 3f\n2:\nmov x0,#1\n3:\nldr x29,[sp],#16\n"

static __attribute__((noinline)) void batch_step(uint64_t t) {
    fe acc, S[4]; /* inv, iw, lam, FIN constants */
    S[3].v[0] = 0xffffffff, S[3].v[1] = 0xfffffffe;
    dpre *DP = (dpre *)(((uintptr_t)dp_mem + DP_SIZE - 1) & ~(DP_SIZE - 1));
    pt *Wc = W, *Wn = WN;
    const int anchor_now = t % ANCHOR_PERIOD == 0;
    k_subb_bl(&DP[0].d, &L[0].sel->x, &Wc[0].x);
    DP[1].pre = DP[0].d;
    uintptr_t d = (uintptr_t)&DP[1].d, r = (uintptr_t)&DP[2].pre, c = (uintptr_t)&Wc[1];
    for (int w = 1; w < B - 1; w++)
        k_submul(&d, &r, &L[w].sel->x, &c);
    k_subb_bl(&DP[B - 1].d, &L[B - 1].sel->x, &Wc[B - 1].x);
    k_subm(&acc, &L[B - 1].sel->x, &Wc[B - 1].x, &DP[B - 1].pre);
    if (fis_zero(&acc)) { exceptional++; abort(); }
    finv(&S[0], &acc);
    int w0 = B - 1;
    if (!anchor_now) {
        /* Steps B-1..1 in one block: inv, the pointers and the reduction constants stay in registers */
        register uint64_t i0 __asm__("x19") = S[0].v[0], i1 __asm__("x20") = S[0].v[1];
        register uint64_t i2 __asm__("x21") = S[0].v[2], i3 __asm__("x22") = S[0].v[3];
        register uintptr_t pd __asm__("x23") = (uintptr_t)&DP[B - 1], pc __asm__("x24") = (uintptr_t)&Wc[B - 1].y;
        register uintptr_t pn __asm__("x25") = (uintptr_t)&Wn[B - 1].y + 16, pl __asm__("x26") = (uintptr_t)&L[B - 1].cur;
        register uintptr_t t2 __asm__("x27") = (uintptr_t)T2;
        register uint64_t rare __asm__("x0");
        for (;;) {
            __asm__ volatile(BLOOP : "+r"(i0), "+r"(i1), "+r"(i2), "+r"(i3), "+r"(pd), "+r"(pc), "+r"(pn), "+r"(pl), "=&r"(rare)
                : "r"(t2), [ab] "i"(ADDER_BITS), [dm] "i"(DP_SIZE - 1), [dpm] "i"(BL_DPM), [fm] "i"(FILTER_MASK)
                : "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16",
                  "x17", "x18", "x28", "x30", "cc", "memory");
            if (!rare) break;
            int w = (int)((pl - (uintptr_t)&L[0].cur) / sizeof(lane) & (B - 1)); /* the mask keeps w's range [0, B) */
            step_rare(w, &Wn[w], Wc, 0);
            pl -= sizeof(lane);
            if (w == 1) break;
        }
        S[0].v[0] = i0;
        S[0].v[1] = i1;
        S[0].v[2] = i2;
        S[0].v[3] = i3;
        w0 = 0;
    }
    uintptr_t m = (uintptr_t)&DP[B].pre, n = (uintptr_t)&DP[B].d;
    for (int w = w0; w >= 0; w--) {
        const pt *R = L[w].sel;
        pt *nw = &Wn[w];
        if (w > 0)
            k_mul2(S, &m, &n); /* iw = inv*PRE[w], inv = inv*D[w]: below 2^256, not canonical */
        else
            S[1] = S[0];
        k_lxy(S, &R->x, &Wc[w].x, &nw->x);
        uint64_t kl = key_lo(nw), f = kl & FILTER_MASK;
        const pt *sel = select_adder1(nw);
        /* step_rare runs the whole tail; on any other step none of its rare branches applies */
        if (__builtin_expect(IS_DP1(kl, nw) | (f == L[w].prev) | (f == L[w].anchor) |
                             (anchor_now && tper[w] >= CAP_PERIODS), 0)) {
            step_rare(w, nw, Wc, anchor_now);
            continue;
        }
        if (anchor_now) {
            L[w].anchor = f;
            tper[w]++;
        }
        L[w].prev = L[w].cur;
        L[w].cur = f;
        L[w].sel = sel;
    }
    W = Wn;
    WN = Wc;
}

static void read_fe(FILE *f, fe *r) {
    char hex[65];
    if (fscanf(f, "%64s", hex) != 1 || strlen(hex) != 64) { fprintf(stderr, "bad input\n"); exit(2); }
    fe a;
    for (int i = 0; i < 4; i++) {
        char limb[17];
        memcpy(limb, hex + 16 * (3 - i), 16);
        limb[16] = 0;
        a.v[i] = strtoull(limb, NULL, 16);
    }
    k_to(r, &a);
}

static void print_fe(const fe *a) {
    fe t;
    k_from(&t, a);
    for (int i = 3; i >= 0; i--) printf("%016llx", (unsigned long long)t.v[i]);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: step <inputs> <batch-steps> [dump|single]\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) return 2;
    for (int i = 0; i < NA; i++) { read_fe(f, &T[i].x); read_fe(f, &T[i].y); }
    for (int w = 0; w < B; w++) { read_fe(f, &W[w].x); read_fe(f, &W[w].y); }
    for (int i = 0; i < SPARES; i++) { read_fe(f, &SPARE[i].x); read_fe(f, &SPARE[i].y); }
    fclose(f);
    for (int i = 0; i < NA; i++) {
        T2[i][0] = T2[i][1] = T[i];
        k_sub(&T2[i][1].y, &ZERO, &T[i].y);
    }
    for (int w = 0; w < B; w++) {
        L[w].sel = select_adder(&W[w]);
        L[w].cur = key_lo(&W[w]) & FILTER_MASK;
        L[w].prev = L[w].anchor = ~0ULL;
        tper[w] = 0;
    }
    int n = atoi(argv[2]);
    if (argc > 3 && strcmp(argv[3], "single") == 0) {
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
        for (int i = 0; i < n; i++) {
            if (!f_single(&W[i % B], &W[i % B])) { exceptional++; abort(); }
        }
        CALLGRIND_TOGGLE_COLLECT;
        fprintf(stderr, "unbatched=%llu\n", (unsigned long long)unbatched);
        return 0;
    }
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for (int t = 1; t <= n; t++) batch_step((uint64_t)t);
    CALLGRIND_TOGGLE_COLLECT;
    fprintf(stderr, "walk_steps=%llu two_cycles=%llu long_cycles=%llu dps=%llu ahead_dps=%llu abandoned=%llu exceptional=%llu unbatched=%llu dp_digest=%016llx loops=%llu\n",
            (unsigned long long)n * B, (unsigned long long)two_cycles, (unsigned long long)long_cycles,
            (unsigned long long)dps, (unsigned long long)ahead_dps, (unsigned long long)abandoned, (unsigned long long)exceptional,
            (unsigned long long)unbatched, (unsigned long long)dp_digest, (unsigned long long)loops);
    if (argc > 3) {
        for (int w = 0; w < B; w++) { print_fe(&W[w].x); printf(" "); print_fe(&W[w].y); printf("\n"); }
    }
    return 0;
}
