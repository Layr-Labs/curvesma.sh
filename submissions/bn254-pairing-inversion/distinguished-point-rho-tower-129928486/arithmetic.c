/* Scalar-word reference for the analytical backend in machine.md.
 * No SIMD, wide address space, arbitrary-size machine operations or tables.
 * __uint128_t is a pair of scalar 64-bit registers in the priced expansion.
 * add_mod/mont_mul follow the rho-baseline-v2-20261110 loop structure; the
 * tower routines implement the Karatsuba formulas priced in machine.md and
 * checked by rho.py. The priced bounds live in machine.md; this file pins
 * the algorithms those bounds count. All loops have fixed public bounds.
 * Tower outputs must not alias inputs; temporaries keep routines safe.
 */
#include <stdint.h>
typedef uint64_t word;
typedef __uint128_t wide;
typedef word fe[4];
typedef fe fp2[2];
typedef fp2 fp6[3];
typedef fp6 fp12[2];

/* m is p or r, both < 2^254. Aliasing out with either operand is permitted. */
static void add_mod(fe out, const fe a, const fe b, const fe m) {
    word s[4], d[4], carry = 0, borrow = 0;
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)a[i]+b[i]+carry;
        s[i]=(word)v; carry=(word)(v>>64);
    }
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)m[i]+borrow;
        d[i]=(word)((wide)s[i]-v);
        borrow=(wide)s[i]<v;
    }
    word mask=0-borrow;
    for (unsigned i=0; i<4; ++i)
        out[i]=(s[i]&mask)|(d[i]&~mask);
}

static void sub_mod(fe out, const fe a, const fe b, const fe m) {
    word s[4], t[4], borrow = 0, carry = 0;
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)a[i]-b[i]-borrow;
        s[i]=(word)v; borrow=(wide)a[i]<(wide)b[i]+borrow;
    }
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)s[i]+m[i]+carry;
        t[i]=(word)v; carry=(word)(v>>64);
    }
    word mask=0-borrow;
    for (unsigned i=0; i<4; ++i)
        out[i]=(t[i]&mask)|(s[i]&~mask);
}

static void copy_fe(fe out, const fe a) {
    for (unsigned i=0; i<4; ++i) out[i]=a[i];
}

/* n0 = -m[0]^{-1} mod 2^64; a,b are canonical Montgomery residues.
 * Output is a*b*2^(-256) mod m. All propagation loops have public bounds.
 */
static void mont_mul(fe out, const fe a, const fe b, const fe m, word n0) {
    word t[9]={0};
    for (unsigned i=0; i<4; ++i) {
        word carry=0;
        for (unsigned j=0; j<4; ++j) {
            wide v=(wide)a[i]*b[j]+t[i+j]+carry;
            t[i+j]=(word)v; carry=(word)(v>>64);
        }
        for (unsigned k=i+4; k<9; ++k) {
            wide v=(wide)t[k]+carry;
            t[k]=(word)v; carry=(word)(v>>64);
        }
    }
    for (unsigned i=0; i<4; ++i) {
        word q=t[i]*n0, carry=0;
        for (unsigned j=0; j<4; ++j) {
            wide v=(wide)q*m[j]+t[i+j]+carry;
            t[i+j]=(word)v; carry=(word)(v>>64);
        }
        for (unsigned k=i+4; k<9; ++k) {
            wide v=(wide)t[k]+carry;
            t[k]=(word)v; carry=(word)(v>>64);
        }
    }
    /* REDC is <2m<2^255, hence t[8]=0 and one subtraction suffices. */
    fe d;
    word borrow=0;
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)m[i]+borrow;
        d[i]=(word)((wide)t[i+4]-v);
        borrow=(wide)t[i+4]<v;
    }
    word mask=0-borrow;
    for (unsigned i=0; i<4; ++i)
        out[i]=(t[i+4]&mask)|(d[i]&~mask);
}

/* Fp2 = Fp[u]/(u^2+1). Elements are (c0, c1) = c0 + c1*u. */
static void fp2_add(fp2 out, const fp2 a, const fp2 b, const fe m) {
    add_mod(out[0], a[0], b[0], m);
    add_mod(out[1], a[1], b[1], m);
}

static void fp2_sub(fp2 out, const fp2 a, const fp2 b, const fe m) {
    sub_mod(out[0], a[0], b[0], m);
    sub_mod(out[1], a[1], b[1], m);
}

/* Complex Karatsuba: 3 Fp muls + 5 Fp add/subs. */
static void fp2_mul(fp2 out, const fp2 a, const fp2 b, const fe m, word n0) {
    fe ac, bd, s, t0, t1;
    mont_mul(ac, a[0], b[0], m, n0);
    mont_mul(bd, a[1], b[1], m, n0);
    add_mod(t0, a[0], a[1], m);
    add_mod(t1, b[0], b[1], m);
    mont_mul(s, t0, t1, m, n0);
    sub_mod(out[0], ac, bd, m);
    sub_mod(t0, s, ac, m);
    sub_mod(out[1], t0, bd, m);
}

static void mul9(fe out, const fe a, const fe m) {
    fe d;
    add_mod(d, a, a, m);
    add_mod(d, d, d, m);
    add_mod(d, d, d, m);
    add_mod(out, d, a, m);
}

/* Multiply by xi = 9+u: (9a0-a1) + (a0+9a1)u. */
static void fp2_mul_xi(fp2 out, const fp2 a, const fe m) {
    fe n0_, n1_;
    mul9(n0_, a[0], m);
    mul9(n1_, a[1], m);
    sub_mod(out[0], n0_, a[1], m);
    add_mod(out[1], a[0], n1_, m);
}

/* Fp6 = Fp2[v]/(v^3-xi). Elements are (b0, b1, b2). */
static void fp6_add(fp6 out, const fp6 a, const fp6 b, const fe m) {
    for (unsigned i=0; i<3; ++i) fp2_add(out[i], a[i], b[i], m);
}

static void fp6_sub(fp6 out, const fp6 a, const fp6 b, const fe m) {
    for (unsigned i=0; i<3; ++i) fp2_sub(out[i], a[i], b[i], m);
}

/* 6 Fp2 muls; c0 = v0+xi*v5, c1 = v3+xi*v2, c2 = v4+v1. */
static void fp6_mul(fp6 out, const fp6 a, const fp6 b, const fe m, word n0) {
    fp2 v0, v1, v2, v3, v4, v5, s0, s1, t;
    fp2_mul(v0, a[0], b[0], m, n0);
    fp2_mul(v1, a[1], b[1], m, n0);
    fp2_mul(v2, a[2], b[2], m, n0);
    fp2_add(s0, a[0], a[1], m);
    fp2_add(s1, b[0], b[1], m);
    fp2_mul(t, s0, s1, m, n0);
    fp2_sub(t, t, v0, m);
    fp2_sub(v3, t, v1, m);
    fp2_add(s0, a[0], a[2], m);
    fp2_add(s1, b[0], b[2], m);
    fp2_mul(t, s0, s1, m, n0);
    fp2_sub(t, t, v0, m);
    fp2_sub(v4, t, v2, m);
    fp2_add(s0, a[1], a[2], m);
    fp2_add(s1, b[1], b[2], m);
    fp2_mul(t, s0, s1, m, n0);
    fp2_sub(t, t, v1, m);
    fp2_sub(v5, t, v2, m);
    fp2_mul_xi(t, v5, m);
    fp2_add(out[0], v0, t, m);
    fp2_mul_xi(t, v2, m);
    fp2_add(out[1], v3, t, m);
    fp2_add(out[2], v4, v1, m);
}

/* Multiply by v: (xi*a2) + a0*v + a1*v^2. Alias-safe via temporaries. */
static void fp6_mul_v(fp6 out, const fp6 a, const fe m) {
    fp2 t, c0, c1;
    fp2_mul_xi(t, a[2], m);
    copy_fe(c0[0], a[0][0]); copy_fe(c0[1], a[0][1]);
    copy_fe(c1[0], a[1][0]); copy_fe(c1[1], a[1][1]);
    copy_fe(out[0][0], t[0]); copy_fe(out[0][1], t[1]);
    copy_fe(out[1][0], c0[0]); copy_fe(out[1][1], c0[1]);
    copy_fe(out[2][0], c1[0]); copy_fe(out[2][1], c1[1]);
}

/* Fp12 = Fp6[t]/(t^2-v). Elements are (a0, a1). Karatsuba: 3 Fp6 muls. */
static void fp12_mul(fp12 out, const fp12 a, const fp12 b,
                     const fe m, word n0) {
    fp6 v0, v1, v2, s0, s1, t;
    fp6_mul(v0, a[0], b[0], m, n0);
    fp6_mul(v1, a[1], b[1], m, n0);
    fp6_add(s0, a[0], a[1], m);
    fp6_add(s1, b[0], b[1], m);
    fp6_mul(t, s0, s1, m, n0);
    fp6_sub(t, t, v0, m);
    fp6_sub(v2, t, v1, m);
    fp6_mul_v(t, v1, m);
    fp6_add(out[0], v0, t, m);
    for (unsigned i=0; i<3; ++i)
        for (unsigned j=0; j<2; ++j) copy_fe(out[1][i][j], v2[i][j]);
}
