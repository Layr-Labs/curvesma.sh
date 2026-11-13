/* Scalar-word reference for the analytical backend in machine.md.
 * No SIMD, wide address space, arbitrary-size machine operations or tables.
 * __uint128_t is a pair of scalar 64-bit registers in the priced expansion.
 *
 * The priced Fp12 arithmetic uses the tower basis of machine.md:
 *   Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-xi), xi = 9+u,
 *   Fp12 = Fp6[w]/(w^2-v).
 * All functions below follow the Karatsuba formulas in machine.md exactly.
 * A residue is four little-endian Montgomery-form words; an Fp2 element is
 * two residues, an Fp6 element six, an Fp12 element twelve.
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
    word d[4], s[4], borrow = 0, carry = 0;
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)b[i]+borrow;
        d[i]=(word)((wide)a[i]-v);
        borrow=(wide)a[i]<v;
    }
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)d[i]+m[i]+carry;
        s[i]=(word)v; carry=(word)(v>>64);
    }
    word mask=0-borrow;
    for (unsigned i=0; i<4; ++i)
        out[i]=(d[i]&~mask)|(s[i]&mask);
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

/* Fp2 add/sub: componentwise. */
static void fp2_add(fp2 out, const fp2 a, const fp2 b, const fe m) {
    add_mod(out[0], a[0], b[0], m);
    add_mod(out[1], a[1], b[1], m);
}
static void fp2_sub(fp2 out, const fp2 a, const fp2 b, const fe m) {
    sub_mod(out[0], a[0], b[0], m);
    sub_mod(out[1], a[1], b[1], m);
}

/* Fp2 multiply, Karatsuba: 3 Fp mults, 5 Fp adds. u^2 = -1. */
static void fp2_mul(fp2 out, const fp2 a, const fp2 b, const fe m, word n0) {
    fe t0, t1, t2, s0, s1;
    mont_mul(t0, a[0], b[0], m, n0);
    mont_mul(t1, a[1], b[1], m, n0);
    add_mod(s0, a[0], a[1], m);
    add_mod(s1, b[0], b[1], m);
    mont_mul(t2, s0, s1, m, n0);
    sub_mod(out[0], t0, t1, m);
    sub_mod(t2, t2, t0, m);
    sub_mod(out[1], t2, t1, m);
}

/* Multiply an Fp2 element by xi = 9+u: (9*x0-x1) + (x0+9*x1)*u.
 * 9*x = 8*x+x costs 3 doubles + 1 add; 10 Fp additions total. */
static void fp2_mul_xi(fp2 out, const fp2 a, const fe m) {
    fe x9_0, x9_1, t;
    add_mod(t, a[0], a[0], m);
    add_mod(t, t, t, m);
    add_mod(x9_0, t, t, m);
    add_mod(x9_0, x9_0, a[0], m);
    add_mod(t, a[1], a[1], m);
    add_mod(t, t, t, m);
    add_mod(x9_1, t, t, m);
    add_mod(x9_1, x9_1, a[1], m);
    sub_mod(out[0], x9_0, a[1], m);
    add_mod(out[1], a[0], x9_1, m);
}

/* Fp6 add/sub: componentwise over Fp2. */
static void fp6_add(fp6 out, const fp6 a, const fp6 b, const fe m) {
    for (unsigned i=0; i<3; ++i) fp2_add(out[i], a[i], b[i], m);
}
static void fp6_sub(fp6 out, const fp6 a, const fp6 b, const fe m) {
    for (unsigned i=0; i<3; ++i) fp2_sub(out[i], a[i], b[i], m);
}

/* Fp6 multiply, Karatsuba: 6 Fp2 mults, 15 Fp2 adds, 2 xi-mults. v^3 = xi. */
static void fp6_mul(fp6 out, const fp6 a, const fp6 b, const fe m, word n0) {
    fp2 v0, v1, v2, v3, v4, v5, s0, s1, t;
    fp2_mul(v0, a[0], b[0], m, n0);
    fp2_mul(v1, a[1], b[1], m, n0);
    fp2_mul(v2, a[2], b[2], m, n0);
    fp2_add(s0, a[0], a[1], m); fp2_add(s1, b[0], b[1], m);
    fp2_mul(v3, s0, s1, m, n0);
    fp2_add(s0, a[1], a[2], m); fp2_add(s1, b[1], b[2], m);
    fp2_mul(v4, s0, s1, m, n0);
    fp2_add(s0, a[0], a[2], m); fp2_add(s1, b[0], b[2], m);
    fp2_mul(v5, s0, s1, m, n0);
    /* c0 = v0 + xi*(v4-v1-v2) */
    fp2_sub(t, v4, v1, m); fp2_sub(t, t, v2, m);
    fp2_mul_xi(t, t, m);
    fp2_add(out[0], v0, t, m);
    /* c1 = v3-v0-v1 + xi*v2 */
    fp2_sub(t, v3, v0, m); fp2_sub(t, t, v1, m);
    fp2_mul_xi(s0, v2, m);
    fp2_add(out[1], t, s0, m);
    /* c2 = v5-v0-v2+v1 */
    fp2_sub(t, v5, v0, m); fp2_sub(t, t, v2, m);
    fp2_add(out[2], t, v1, m);
}

/* Multiply an Fp6 element by v: v*(d0+d1*v+d2*v^2) = xi*d2 + d0*v + d1*v^2.
 * Alias-safe: the input is copied to scratch before any output is written. */
static void fp6_mul_v(fp6 out, const fp6 a, const fe m) {
    fp6 tmp;
    for (unsigned i=0; i<3; ++i)
        for (unsigned j=0; j<2; ++j)
            for (unsigned k=0; k<4; ++k)
                tmp[i][j][k]=a[i][j][k];
    fp2_mul_xi(out[0], tmp[2], m);
    for (unsigned j=0; j<2; ++j)
        for (unsigned k=0; k<4; ++k) {
            out[1][j][k]=tmp[0][j][k];
            out[2][j][k]=tmp[1][j][k];
        }
}

/* Fp12 multiply, Karatsuba: 3 Fp6 mults, 5 Fp6 adds, 1 v-mult. w^2 = v. */
static void fp12_mul(fp12 out, const fp12 a, const fp12 b, const fe m, word n0) {
    fp6 t0, t1, t2, s0, s1;
    fp6_mul(t0, a[0], b[0], m, n0);
    fp6_mul(t1, a[1], b[1], m, n0);
    fp6_add(s0, a[0], a[1], m); fp6_add(s1, b[0], b[1], m);
    fp6_mul(t2, s0, s1, m, n0);
    fp6_mul_v(s0, t1, m);
    fp6_add(out[0], t0, s0, m);
    fp6_sub(t2, t2, t0, m);
    fp6_sub(out[1], t2, t1, m);
}

/* Subtraction and negation use subtraction/borrow, conditional addition of m,
 * and masking of the zero case. Their fixed four-limb bounds are in machine.md.
 * Higher operations follow the explicit bounded algorithms in that file.
 */
