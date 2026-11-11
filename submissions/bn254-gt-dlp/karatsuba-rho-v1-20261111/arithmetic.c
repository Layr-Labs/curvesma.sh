/* Scalar-word reference for the analytical backend in machine.md.
 * No SIMD, wide address space, arbitrary-size machine operations or tables.
 * __uint128_t is a pair of scalar 64-bit registers in the priced expansion.
 */
#include <stdint.h>
typedef uint64_t word;
typedef __uint128_t wide;
typedef word fe[4];

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

/* Canonical a-b mod m. The first subtraction is modulo 2^256; if it
 * underflows, adding m modulo 2^256 gives the canonical nonnegative answer.
 * Inputs are copied into d before writing out, so aliasing is safe.
 */
static void sub_mod(fe out, const fe a, const fe b, const fe m) {
    fe d;
    word borrow=0;
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)b[i]+borrow;
        d[i]=(word)((wide)a[i]-v);
        borrow=(wide)a[i]<v;
    }
    word mask=0-borrow, carry=0;
    for (unsigned i=0; i<4; ++i) {
        wide v=(wide)d[i]+(m[i]&mask)+carry;
        out[i]=(word)v; carry=(word)(v>>64);
    }
}

static void zero_poly(fe *out, unsigned n) {
    for (unsigned i=0; i<n; ++i)
        for (unsigned j=0; j<4; ++j) out[i][j]=0;
}

/* Fixed degree-2 schoolbook base case: nine coefficient multiplications.
 * Polynomial helper outputs do not alias their inputs.
 */
static void poly3(fe out[5], const fe a[3], const fe b[3],
                  const fe p, word n0) {
    zero_poly(out,5);
    for (unsigned i=0; i<3; ++i)
        for (unsigned j=0; j<3; ++j) {
            fe t;
            mont_mul(t,a[i],b[j],p,n0);
            add_mod(out[i+j],out[i+j],t,p);
        }
}

/* a=a0+X^3*a1; (a0+a1)(b0+b1)-a0*b0-a1*b1 is the cross term.
 * Three poly3 calls, six input additions, ten cross subtractions,
 * fifteen output additions. All coefficients remain canonical.
 */
static void poly6(fe out[11], const fe a[6], const fe b[6],
                  const fe p, word n0) {
    fe sa[3],sb[3],z0[5],z1[5],z2[5];
    for (unsigned i=0; i<3; ++i) {
        add_mod(sa[i],a[i],a[i+3],p);
        add_mod(sb[i],b[i],b[i+3],p);
    }
    poly3(z0,a,b,p,n0);
    poly3(z2,a+3,b+3,p,n0);
    poly3(z1,sa,sb,p,n0);
    zero_poly(out,11);
    for (unsigned i=0; i<5; ++i) {
        sub_mod(z1[i],z1[i],z0[i],p);
        sub_mod(z1[i],z1[i],z2[i],p);
        add_mod(out[i],out[i],z0[i],p);
        add_mod(out[i+3],out[i+3],z1[i],p);
        add_mod(out[i+6],out[i+6],z2[i],p);
    }
}

/* The second Karatsuba level: three poly6 calls, hence 81 mont_mul calls.
 * Twelve input additions, 22 cross subtractions, 33 output additions.
 */
static void poly12(fe out[23], const fe a[12], const fe b[12],
                   const fe p, word n0) {
    fe sa[6],sb[6],z0[11],z1[11],z2[11];
    for (unsigned i=0; i<6; ++i) {
        add_mod(sa[i],a[i],a[i+6],p);
        add_mod(sb[i],b[i],b[i+6],p);
    }
    poly6(z0,a,b,p,n0);
    poly6(z2,a+6,b+6,p,n0);
    poly6(z1,sa,sb,p,n0);
    zero_poly(out,23);
    for (unsigned i=0; i<11; ++i) {
        sub_mod(z1[i],z1[i],z0[i],p);
        sub_mod(z1[i],z1[i],z2[i],p);
        add_mod(out[i],out[i],z0[i],p);
        add_mod(out[i+6],out[i+6],z1[i],p);
        add_mod(out[i+12],out[i+12],z2[i],p);
    }
}

/* Protected basis: X^12=18*X^6-82. Descending reduction includes the
 * coefficient generated at degree 16 by reducing degree 22.
 * Caller output may alias a or b: convolution finishes before copying out.
 */
static void fp12_mul_karatsuba(fe out[12], const fe a[12], const fe b[12],
                              const fe p, word n0) {
    fe c[23];
    poly12(c,a,b,p,n0);
    for (int k=22; k>=12; --k) {
        fe x2,x4,x8,x16,x18,x32,x64,x80,x82;
        add_mod(x2,c[k],c[k],p);
        add_mod(x4,x2,x2,p);
        add_mod(x8,x4,x4,p);
        add_mod(x16,x8,x8,p);
        add_mod(x18,x16,x2,p);
        add_mod(x32,x16,x16,p);
        add_mod(x64,x32,x32,p);
        add_mod(x80,x64,x16,p);
        add_mod(x82,x80,x2,p);
        add_mod(c[k-6],c[k-6],x18,p);
        sub_mod(c[k-12],c[k-12],x82,p);
    }
    for (unsigned i=0; i<12; ++i)
        for (unsigned j=0; j<4; ++j) out[i][j]=c[i][j];
}
