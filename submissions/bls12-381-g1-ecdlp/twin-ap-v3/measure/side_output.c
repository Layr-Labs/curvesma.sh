/* twin-ap-v2: instruction count of the twin step's side output P - G (proof.md Section 2).
   Counted exactly like scripts/reference_costs/bls_field.c: Callgrind Ir over n iterations, with
   price = (Ir(8192) - Ir(4096)) / 4096, inside the organizers' pinned AArch64 image and flags.

   One lane: P = (x1, y1) is the chain point, G = (x2, y2) the gap-table entry and i = 1/(x2 - x1) the lane inverse
   from the Montgomery batch of the reference addition P + G. The side output is
       mu = (y1 + y2) i,   x' = mu^2 - x1 - x2,   y' = mu (x' - x1) - y1
   i.e. 2 multiplications, 1 squaring, 1 addition and 4 subtractions mod p. Every operand is read from memory
   and both x' and y' are stored (x' is not needed by the record, so storing it is conservative).

   mode "lanes":    the charged kernel. blst's own field layer (src/fields.h: mul_fp, sqr_fp, add_fp, sub_fp, which
                    call the AArch64 assembly mul_mont_384, sqr_mont_384, add_mod_384, sub_mod_384), compiled with
                    the library, inside the lane loop of one batch round.
   mode "call":     the same kernel as a separate (noinline) function called once per lane.
   mode "exported": the same kernel through the exported API blst_fp_add/mul/sqr/sub, a noinline function per lane.
   The two extra modes only show how much call overhead the measurement contains. After the counted region the
   program checks its outputs (see the end of main) and exits non-zero on any mismatch. */
#include "/blst/src/server.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#define LANES 8192
typedef struct { vec384 x, y; } aff;
static aff P[LANES], G[LANES], OUT[LANES];
static vec384 INV[LANES];

static inline void side_output(aff *out, const aff *p, const aff *g, const vec384 i)
{
    vec384 t, mu, w;
    add_fp(t, p->y, g->y);                  /* y1 + y2          */
    mul_fp(mu, t, i);                       /* mu               */
    sqr_fp(w, mu);                          /* mu^2             */
    sub_fp(w, w, p->x);                     /* mu^2 - x1        */
    sub_fp(out->x, w, g->x);                /* x'               */
    sub_fp(w, out->x, p->x);                /* x' - x1          */
    mul_fp(w, mu, w);                       /* mu (x' - x1)     */
    sub_fp(out->y, w, p->y);                /* y'               */
}
__attribute__((noinline)) static void side_lanes(int n)
{
    for (int j = 0; j < n; j++) side_output(&OUT[j], &P[j], &G[j], INV[j]);
}
__attribute__((noinline)) static void side_call(aff *out, const aff *p, const aff *g, const vec384 i)
{
    side_output(out, p, g, i);
}
__attribute__((noinline)) static void side_exported(aff *out, const aff *p, const aff *g, const vec384 i)
{   /* exported entry points of exports.c (the signatures there take vec384 inside the library build) */
    vec384 t, mu, w;
    blst_fp_add(t, p->y, g->y);
    blst_fp_mul(mu, t, i);
    blst_fp_sqr(w, mu);
    blst_fp_sub(w, w, p->x);
    blst_fp_sub(out->x, w, g->x);
    blst_fp_sub(w, out->x, p->x);
    blst_fp_mul(w, mu, w);
    blst_fp_sub(out->y, w, p->y);
}
int main(int argc, char **argv)
{
    if (argc < 3) return 1;
    const char *mode = argv[1];
    int n = atoi(argv[2]);
    if (n < 1 || n > LANES) return 1;
    vec384 a, b;
    unsigned long long av[6] = {17}, bv[6] = {42};
    blst_fp_from_uint64(a, av);
    blst_fp_from_uint64(b, bv);
    for (int j = 0; j < LANES; j++) {       /* distinct valid Montgomery residues */
        mul_fp(a, a, b); vec_copy(P[j].x, a, 48);
        mul_fp(a, a, b); vec_copy(P[j].y, a, 48);
        mul_fp(a, a, b); vec_copy(G[j].x, a, 48);
        mul_fp(a, a, b); vec_copy(G[j].y, a, 48);
        mul_fp(a, a, b); vec_copy(INV[j], a, 48);
    }
    int lanes = !strcmp(mode, "lanes"), call = !strcmp(mode, "call"), exported = !strcmp(mode, "exported");
    if (!(lanes || call || exported)) return 1;
    if (lanes) side_lanes(1); else if (call) side_call(&OUT[0], &P[0], &G[0], INV[0]);
    else side_exported(&OUT[0], &P[0], &G[0], INV[0]);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    if (lanes) side_lanes(n);
    else if (call) for (int j = 0; j < n; j++) side_call(&OUT[j], &P[j], &G[j], INV[j]);
    else for (int j = 0; j < n; j++) side_exported(&OUT[j], &P[j], &G[j], INV[j]);
    CALLGRIND_TOGGLE_COLLECT;
    /* self-check, outside the counted region: (1) every counted lane equals the exported-API computation;
       (2) on 16 pairs of genuine G1 points, the kernel's (x', y') equals blst's own affine P + (-G). */
    for (int j = 0; j < n; j++) {
        aff ref;
        side_exported(&ref, &P[j], &G[j], INV[j]);
        if (!vec_is_equal(ref.y, OUT[j].y, 48) || !vec_is_equal(ref.x, OUT[j].x, 48)) { printf("lane mismatch\n"); return 2; }
    }
    for (int k = 1; k <= 16; k++) {
        POINTonE1 Pj, Gj, R; POINTonE1_affine Pa, Ga, nGa, Ra; aff pt, gt, ot; vec384 dx, iv;
        byte e1[2] = {(byte)(3 * k + 5), 0}, e2[2] = {(byte)(7 * k + 11), 0};
        blst_p1_mult(&Pj, blst_p1_generator(), e1, 8);
        blst_p1_mult(&Gj, blst_p1_generator(), e2, 8);
        blst_p1_to_affine(&Pa, &Pj); blst_p1_to_affine(&Ga, &Gj);
        vec_copy(pt.x, Pa.X, 48); vec_copy(pt.y, Pa.Y, 48); vec_copy(gt.x, Ga.X, 48); vec_copy(gt.y, Ga.Y, 48);
        sub_fp(dx, gt.x, pt.x); blst_fp_inverse(iv, dx);
        side_output(&ot, &pt, &gt, iv);
        nGa = Ga; cneg_fp(nGa.Y, Ga.Y, 1);
        blst_p1_add_or_double_affine(&R, &Pj, &nGa);
        blst_p1_to_affine(&Ra, &R);
        if (!vec_is_equal(Ra.X, ot.x, 48) || !vec_is_equal(Ra.Y, ot.y, 48)) { printf("curve mismatch\n"); return 3; }
    }
    printf("selfcheck ok: %d lanes match the exported API, 16 curve pairs match P + (-G)\n", n);
    printf("%lu\n", (unsigned long)OUT[n - 1].y[0]);
    return 0;
}
