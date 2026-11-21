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

/* Subtraction and negation use subtraction/borrow, conditional addition of m,
 * and masking of the zero case. Their fixed four-limb bounds are in machine.md.
 * Higher operations follow the explicit bounded algorithms in that file.
 */
