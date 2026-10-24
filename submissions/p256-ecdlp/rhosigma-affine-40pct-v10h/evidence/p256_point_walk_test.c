/* point_walk_test.c — REAL P-256 POINT WALK through the SHIPPED counted step.
 *
 * v10h (verdict b8b00688 item 2): the v10g kernel test called only field
 * primitives and canon_y, so it never exercised the counted point-walk routine.
 * This driver #includes evidence/step_v10f_o3.c VERBATIM -- the exact source
 * that stepcount.json instruction-counts -- and runs it on REAL P-256 points
 * (y^2 = x^3 - 3x + b over p = 2^256 - 2^224 + 2^192 + 2^96 - 1) with the real
 * 2048-row two-table schedule, so every fe_mul/fe_add/fe_sub/fe_inv on the
 * counted path executes on the actual P-256 field.
 *
 * Representation: Montgomery (fe_mul(a,b) = a*b*R^-1). Points are held as Mont
 * pairs; canon_y and the index rule therefore act on the Mont-form limbs exactly
 * as they do in production. Canonical output is obtained with ONE fe_mul by 1.
 *
 * Emits (stdout, one record per line):
 *   TBL/TB2 <j> <xcanon> <ycanon>       all 2048 rows of each table
 *   DIGEST <h1> <h2>                    FNV-1a over the canonical row bytes
 *   WALK <block> <slot> <s_hex> <x> <y> start scalar and final point of a 512-step walk
 *   CHK <block> <slot> <x> <y>          mid-walk (step 256) checkpoint, block 0
 * evidence/point_walk_test.py recomputes every one of them in Python big
 * integers at the real P-256 modulus and writes point_walk.json.
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "step_v10f_o3.c"     /* the shipped counted step: includes p256_fe_cc.h */

/* the step's extern storage (defined here, exactly as the counted TU expects) */
apt TBL[2048], TBL2[2048];
awalker A[16];
uint64_t PROD[16][4], SUF[16][4];

/* P-256 constants, canonical little-endian limbs */
static const uint64_t CX[4] = {0xf4a13945d898c296ULL, 0x77037d812deb33a0ULL,
                               0xf8bce6e563a440f2ULL, 0x6b17d1f2e12c4247ULL};   /* Gx */
static const uint64_t CY[4] = {0xcbb6406837bf51f5ULL, 0x2bce33576b315eceULL,
                               0x8ee7eb4a7c0f9e16ULL, 0x4fe342e2fe1a7f9bULL};   /* Gy */
static const uint64_t R2[4] = {0x0000000000000003ULL, 0xfffffffbffffffffULL,
                               0xfffffffffffffffeULL, 0x00000004fffffffdULL};   /* R^2 mod p */
static const uint64_t ONE[4] = {1, 0, 0, 0};
static const uint64_t THREE[4] = {3, 0, 0, 0};
/* a = -3 in MONTGOMERY form: 3R mod p (fe_sub/fe_add are representation-preserving,
 * so the curve constant must be carried in the same representation as the point) */
static const uint64_t AMONT[4] = {0x0000000000000003ULL, 0xfffffffd00000000ULL,
                                  0xffffffffffffffffULL, 0x00000002fffffffcULL};

static inline void to_mont(uint64_t r[4], const uint64_t a[4]) { fe_mul(r, a, R2); }
static inline void from_mont(uint64_t r[4], const uint64_t a[4]) { fe_mul(r, a, ONE); }

static inline uint64_t splitmix64(uint64_t *s)
{ uint64_t z = (*s += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }

/* affine P+Q in Montgomery form, using ONLY shipped header primitives */
static void pt_add(uint64_t x3[4], uint64_t y3[4], const uint64_t x1[4], const uint64_t y1[4],
                   const uint64_t x2[4], const uint64_t y2[4])
{
    uint64_t d[4], num[4], inv[4], lam[4], t[4];
    fe_sub(d, x2, x1);
    fe_sub(num, y2, y1);
    fe_inv(inv, d);                 /* R^2/d : Montgomery representation of 1/d */
    fe_mul(lam, num, inv);          /* (num/d) * R */
    fe_mul(t, lam, lam);
    fe_sub(x3, t, x1); fe_sub(x3, x3, x2);
    fe_sub(t, x1, x3);
    fe_mul(y3, lam, t); fe_sub(y3, y3, y1);
}

/* affine doubling, a = -3 : lam = 3(x^2 - 1) / 2y */
static void pt_dbl(uint64_t x3[4], uint64_t y3[4], const uint64_t x1[4], const uint64_t y1[4])
{
    uint64_t num[4], den[4], inv[4], lam[4], t[4];
    fe_mul(t, x1, x1);
    fe_add(num, t, t); fe_add(num, num, t);      /* 3x^2 */
    fe_sub(num, num, AMONT);                     /* + a, a = -3 in Montgomery form */
    fe_add(den, y1, y1);
    fe_inv(inv, den);
    fe_mul(lam, num, inv);
    fe_mul(t, lam, lam);
    fe_sub(x3, t, x1); fe_sub(x3, x3, x1);
    fe_sub(t, x1, x3);
    fe_mul(y3, lam, t); fe_sub(y3, y3, y1);
}

static uint64_t fnv(uint64_t h, const uint64_t a[4])
{ for (int i = 0; i < 4; i++) { uint64_t v = a[i]; for (int k = 0; k < 8; k++) { h ^= (v & 0xff); h *= 0x100000001b3ULL; v >>= 8; } } return h; }

static void emit(const uint64_t x[4], const uint64_t y[4])
{
    for (int i = 3; i >= 0; i--) printf("%016llx", (unsigned long long)x[i]);
    for (int i = 3; i >= 0; i--) printf("%016llx", (unsigned long long)y[i]);
}

/* [s]G on Montgomery forms, double-and-add, shipped primitives only */
static void mul_g(uint64_t rx[4], uint64_t ry[4], uint64_t s)
{
    uint64_t gx[4], gy[4], ax[4], ay[4];
    to_mont(gx, CX); to_mont(gy, CY);
    memcpy(ax, gx, sizeof ax); memcpy(ay, gy, sizeof ay);
    int started = 0;
    for (int b = 63; b >= 0; b--) {
        if (started) { uint64_t tx[4], ty[4]; pt_dbl(tx, ty, ax, ay); memcpy(ax, tx, sizeof tx); memcpy(ay, ty, sizeof ty); }
        if ((s >> b) & 1) {
            if (!started) { memcpy(ax, gx, sizeof ax); memcpy(ay, gy, sizeof ay); started = 1; }
            else { uint64_t tx[4], ty[4]; pt_add(tx, ty, ax, ay, gx, gy); memcpy(ax, tx, sizeof tx); memcpy(ay, ty, sizeof ty); }
        }
    }
    memcpy(rx, ax, sizeof ax); memcpy(ry, ay, sizeof ay);
}

int main(void)
{
    /* ---- tables: TBL[j] = [j+1]G, TBL2[j] = [j+1]*K2 G, both canon (even Mont y) ---- */
    const uint64_t K2 = 0x9e3779b97f4a7c15ULL;
    {
        uint64_t gx[4], gy[4], px[4], py[4];
        to_mont(gx, CX); to_mont(gy, CY);
        memcpy(px, gx, sizeof px); memcpy(py, gy, sizeof py);
        for (int j = 0; j < 2048; j++) {
            canon_y(py);
            memcpy(TBL[j].x, px, sizeof px); memcpy(TBL[j].y, py, sizeof py);
            uint64_t nx[4], ny[4];
            if (j == 0) pt_dbl(nx, ny, px, py);      /* [2]G: the chord formula needs x1!=x2 */
            else pt_add(nx, ny, px, py, gx, gy);
            memcpy(px, nx, sizeof nx); memcpy(py, ny, sizeof ny);
        }
        uint64_t k2x[4], k2y[4];
        mul_g(k2x, k2y, K2);
        memcpy(px, k2x, sizeof px); memcpy(py, k2y, sizeof py);
        for (int j = 0; j < 2048; j++) {
            canon_y(py);
            memcpy(TBL2[j].x, px, sizeof px); memcpy(TBL2[j].y, py, sizeof py);
            uint64_t nx[4], ny[4];
            if (j == 0) pt_dbl(nx, ny, px, py);
            else pt_add(nx, ny, px, py, k2x, k2y);
            memcpy(px, nx, sizeof nx); memcpy(py, ny, sizeof ny);
        }
    }

    uint64_t h1 = 0xcbf29ce484222325ULL, h2 = 0xcbf29ce484222325ULL;
    for (int j = 0; j < 2048; j++) {
        uint64_t cx[4], cy[4];
        from_mont(cx, TBL[j].x); from_mont(cy, TBL[j].y);
        h1 = fnv(fnv(h1, cx), cy);
        printf("TBL %d ", j); emit(cx, cy); printf("\n");
        from_mont(cx, TBL2[j].x); from_mont(cy, TBL2[j].y);
        h2 = fnv(fnv(h2, cx), cy);
        printf("TB2 %d ", j); emit(cx, cy); printf("\n");
    }
    printf("DIGEST %016llx %016llx\n", (unsigned long long)h1, (unsigned long long)h2);

    /* ---- 13 blocks x 16 walkers x 512 steps = 208 walks, 106,496 counted steps ---- */
    uint64_t rs = 0x1234abcd5678ef90ULL;
    const int NBLK = 13, NSTEP = 512, MID = 256;
    for (int blk = 0; blk < NBLK; blk++) {
        uint64_t ss[16];
        for (int i = 0; i < 16; i++) {
            uint64_t s = splitmix64(&rs) | 1ULL;
            ss[i] = s;
            mul_g(A[i].x, A[i].y, s);
            for (int k = 0; k < 4; k++) { PROD[i][k] = (uint64_t)(k == 0); SUF[i][k] = (uint64_t)(k == 0); }
            A[i].len = 0; A[i].esc = 0;
        }
        for (int st = 0; st < NSTEP; st++) {
            step_v10f_o3(16);          /* ONE step for each of the 16 walkers */
            if (st + 1 == MID && blk == 0) {
                for (int i = 0; i < 16; i++) {
                    uint64_t cx[4], cy[4];
                    from_mont(cx, A[i].x); from_mont(cy, A[i].y);
                    printf("CHK %d %d ", blk, i); emit(cx, cy); printf("\n");
                }
            }
        }
        for (int i = 0; i < 16; i++) {
            uint64_t cx[4], cy[4];
            from_mont(cx, A[i].x); from_mont(cy, A[i].y);
            printf("WALK %d %d %016llx ", blk, i, (unsigned long long)ss[i]);
            emit(cx, cy); printf("\n");
        }
    }
    return 0;
}
