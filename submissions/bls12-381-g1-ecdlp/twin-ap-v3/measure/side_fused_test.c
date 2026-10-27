/* twin-ap-v3: count and check the fused side-output kernel tw_side_fused (side_fused.S), exactly like
   scripts/reference_costs/bls_field.c: Callgrind Ir over n lanes, price = (Ir(8192) - Ir(4096)) / 4096, in the
   organizers' pinned AArch64 image and flags. Counted region: one C call of tw_side_fused per lane.
   With "check" as the second argument, after the counted region the program verifies (exit code != 0 on failure):
   (1) every lane against the same side output computed with blst's own field layer (mul_fp, sqr_fp, add_fp,
       sub_fp, all fully reduced), with the first 625 lanes built from the extreme operands {0, 1, p-1, p-2, random};
   (2) the lazy-input Montgomery multiplication used by the kernel: tw__mul_mont_384(r + k p, b) = mul_fp(r, b) for
       k = 0..3, 200000 pairs (r, b) including r = b = p - 1;
   (3) on 64 pairs of genuine G1 points ([3k+5]P0, [7k+11]P0), y' equals the y of blst's affine P + (-G). */
#include "/blst/src/server.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#define LANES 8192
typedef struct { vec384 x, y; } aff;
static aff P[LANES], G[LANES], OUT[LANES];
static vec384 INV[LANES];
void tw_side_fused(aff *out, const aff *p, const aff *g, const limb_t *inv, const limb_t *p_and_3p, limb_t n0);
void tw_mul_mont_384(limb_t *ret, const limb_t *a, const limb_t *b, const limb_t *mod, limb_t n0);
static limb_t P3[12];                                   /* p, then 3p as a plain 384-bit integer */
static void add_int(limb_t *r, const limb_t *a, const limb_t *b)   /* r = a + b, 6 limbs, no reduction */
{
    unsigned __int128 c = 0;
    for (int k = 0; k < 6; k++) { c += (unsigned __int128)a[k] + b[k]; r[k] = (limb_t)c; c >>= 64; }
}
static inline void side_ref(aff *out, const aff *p, const aff *g, const vec384 i)
{
    vec384 t, mu, w;
    add_fp(t, p->y, g->y); mul_fp(mu, t, i); sqr_fp(w, mu);
    sub_fp(w, w, p->x); sub_fp(out->x, w, g->x); sub_fp(w, out->x, p->x);
    mul_fp(w, mu, w); sub_fp(out->y, w, p->y);
}
__attribute__((noinline)) static void fused_lanes(int n)
{
    for (int j = 0; j < n; j++) tw_side_fused(&OUT[j], &P[j], &G[j], INV[j], P3, p0);
}
int main(int argc, char **argv)
{
    if (argc < 2) return 1;
    int n = atoi(argv[1]);
    if (n < 1 || n > LANES) return 1;
    vec384 a, b;
    unsigned long long av[6] = {17}, bv[6] = {42};
    blst_fp_from_uint64(a, av);
    blst_fp_from_uint64(b, bv);
    for (int j = 0; j < LANES; j++) {
        mul_fp(a, a, b); vec_copy(P[j].x, a, 48);
        mul_fp(a, a, b); vec_copy(P[j].y, a, 48);
        mul_fp(a, a, b); vec_copy(G[j].x, a, 48);
        mul_fp(a, a, b); vec_copy(G[j].y, a, 48);
        mul_fp(a, a, b); vec_copy(INV[j], a, 48);
    }
    vec_copy(P3, BLS12_381_P, 48);
    { limb_t t2[6]; add_int(t2, BLS12_381_P, BLS12_381_P); add_int(P3 + 6, t2, BLS12_381_P); }
    /* extreme operands in the first 625 lanes: every field of a lane from {0, 1, p-1, p-2, random} */
    {
        vec384 ext[5]; vec_zero(ext[0], 48); vec_zero(ext[1], 48); ext[1][0] = 1;
        vec_copy(ext[2], BLS12_381_P, 48); ext[2][0] -= 1; vec_copy(ext[3], BLS12_381_P, 48); ext[3][0] -= 2;
        for (int j = 0; j < 625; j++) {
            int c = j;
            vec_copy(ext[4], P[j].x, 48);
            vec_copy(P[j].x, ext[c % 5], 48); c /= 5;
            vec_copy(P[j].y, ext[c % 5], 48); c /= 5;
            vec_copy(G[j].x, ext[c % 5], 48); c /= 5;
            vec_copy(G[j].y, ext[c % 5], 48);
        }
    }
    fused_lanes(1);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    fused_lanes(n);
    CALLGRIND_TOGGLE_COLLECT;
    if (argc < 3 || strcmp(argv[2], "check")) { printf("%lu\n", (unsigned long)OUT[n - 1].y[0]); return 0; }
    for (int j = 0; j < n; j++) {
        aff ref; side_ref(&ref, &P[j], &G[j], INV[j]);
        if (!vec_is_equal(ref.y, OUT[j].y, 48)) { printf("lane mismatch %d\n", j); return 2; }
    }
    /* lazy-input Montgomery multiplication: a = r + k p for k = 0..3 (so a < 4p), b < p */
    {
        vec384 r, bb, ref, got, lz; unsigned long long rv[6] = {5}, bv2[6] = {9};
        blst_fp_from_uint64(r, rv); blst_fp_from_uint64(bb, bv2);
        for (int t = 0; t < 200000; t++) {
            mul_fp(r, r, bb); add_fp(bb, bb, r);
            if (t == 0) { vec_copy(r, BLS12_381_P, 48); r[0] -= 1; vec_copy(bb, BLS12_381_P, 48); bb[0] -= 1; }
            mul_fp(ref, r, bb);
            for (int k = 0; k < 4; k++) {
                vec_copy(lz, r, 48);
                for (int q = 0; q < k; q++) add_int(lz, lz, BLS12_381_P);
                tw_mul_mont_384(got, lz, bb, BLS12_381_P, p0);
                if (!vec_is_equal(ref, got, 48)) { printf("lazy mul mismatch t=%d k=%d\n", t, k); return 4; }
            }
        }
    }
    for (int k = 1; k <= 64; k++) {
        POINTonE1 Pj, Gj, R; POINTonE1_affine Pa, Ga, nGa, Ra; aff pt, gt, ot; vec384 dx, iv;
        unsigned v1 = 3 * k + 5, v2 = 7 * k + 11;
        byte e1[2] = {(byte)(v1 & 255), (byte)(v1 >> 8)}, e2[2] = {(byte)(v2 & 255), (byte)(v2 >> 8)};
        blst_p1_mult(&Pj, blst_p1_generator(), e1, 16);
        blst_p1_mult(&Gj, blst_p1_generator(), e2, 16);
        blst_p1_to_affine(&Pa, &Pj); blst_p1_to_affine(&Ga, &Gj);
        vec_copy(pt.x, Pa.X, 48); vec_copy(pt.y, Pa.Y, 48); vec_copy(gt.x, Ga.X, 48); vec_copy(gt.y, Ga.Y, 48);
        sub_fp(dx, gt.x, pt.x); blst_fp_inverse(iv, dx);
        tw_side_fused(&ot, &pt, &gt, iv, P3, p0);
        nGa = Ga; cneg_fp(nGa.Y, Ga.Y, 1);
        blst_p1_add_or_double_affine(&R, &Pj, &nGa);
        blst_p1_to_affine(&Ra, &R);
        if (!vec_is_equal(Ra.Y, ot.y, 48)) { printf("curve mismatch %d\n", k); return 3; }
    }
    printf("selfcheck ok: %d lanes (625 with extreme operands) match the blst field-layer kernel, 800000 lazy-input Montgomery products match, 64 curve pairs match y(P + (-G))\n", n);
    return 0;
}
