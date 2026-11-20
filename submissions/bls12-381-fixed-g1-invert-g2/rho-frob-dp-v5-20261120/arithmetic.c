/* Scalar-word reference for the analytical backend in machine.md.
 * No SIMD, wide address space, arbitrary-size machine operations or tables.
 * __uint128_t is a pair of scalar 64-bit registers in the priced expansion.
 * NL is the public limb count: NL=6 for the BLS12-381 field prime p
 * (p < 2^381) and NL=4 for the subgroup order r (r < 2^255). Two separate
 * instantiations are emitted; every loop bound below is a public constant.
 */
#include <stdint.h>
typedef uint64_t word;
typedef __uint128_t wide;

/* m < 2^(64*NL-1), so a+b < 2m cannot overflow NL words.
 * Aliasing out with either operand is permitted. */
#define DEFINE_ADD_MOD(NL, name)                                         \
static void name(word out[NL], const word a[NL], const word b[NL],       \
                 const word m[NL]) {                                     \
    word s[NL], d[NL], carry = 0, borrow = 0;                            \
    for (unsigned i=0; i<NL; ++i) {                                      \
        wide v=(wide)a[i]+b[i]+carry;                                    \
        s[i]=(word)v; carry=(word)(v>>64);                               \
    }                                                                    \
    for (unsigned i=0; i<NL; ++i) {                                      \
        wide v=(wide)m[i]+borrow;                                        \
        d[i]=(word)((wide)s[i]-v);                                       \
        borrow=(wide)s[i]<v;                                             \
    }                                                                    \
    word mask=0-borrow;                                                  \
    for (unsigned i=0; i<NL; ++i)                                        \
        out[i]=(s[i]&mask)|(d[i]&~mask);                                 \
}

/* n0 = -m[0]^{-1} mod 2^64; a,b are canonical Montgomery residues.
 * Output is a*b*2^(-64*NL) mod m. All propagation loops have public bounds.
 * REDC result is < 2m < 2^(64*NL), hence t[2*NL]=0 and one subtraction
 * of m makes it canonical.
 */
#define DEFINE_MONT_MUL(NL, name)                                        \
static void name(word out[NL], const word a[NL], const word b[NL],       \
                 const word m[NL], word n0) {                            \
    word t[2*NL+1]={0};                                                  \
    for (unsigned i=0; i<NL; ++i) {                                      \
        word carry=0;                                                    \
        for (unsigned j=0; j<NL; ++j) {                                  \
            wide v=(wide)a[i]*b[j]+t[i+j]+carry;                         \
            t[i+j]=(word)v; carry=(word)(v>>64);                         \
        }                                                                \
        for (unsigned k=i+NL; k<2*NL+1; ++k) {                           \
            wide v=(wide)t[k]+carry;                                     \
            t[k]=(word)v; carry=(word)(v>>64);                           \
        }                                                                \
    }                                                                    \
    for (unsigned i=0; i<NL; ++i) {                                      \
        word q=t[i]*n0, carry=0;                                         \
        for (unsigned j=0; j<NL; ++j) {                                  \
            wide v=(wide)q*m[j]+t[i+j]+carry;                            \
            t[i+j]=(word)v; carry=(word)(v>>64);                         \
        }                                                                \
        for (unsigned k=i+NL; k<2*NL+1; ++k) {                           \
            wide v=(wide)t[k]+carry;                                     \
            t[k]=(word)v; carry=(word)(v>>64);                           \
        }                                                                \
    }                                                                    \
    word d[NL], borrow=0;                                                \
    for (unsigned i=0; i<NL; ++i) {                                      \
        wide v=(wide)m[i]+borrow;                                        \
        d[i]=(word)((wide)t[i+NL]-v);                                    \
        borrow=(wide)t[i+NL]<v;                                          \
    }                                                                    \
    word mask=0-borrow;                                                  \
    for (unsigned i=0; i<NL; ++i)                                        \
        out[i]=(t[i+NL]&mask)|(d[i]&~mask);                              \
}

DEFINE_ADD_MOD(6, add_mod_p)
DEFINE_MONT_MUL(6, mont_mul_p)
DEFINE_ADD_MOD(4, add_mod_r)
DEFINE_MONT_MUL(4, mont_mul_r)

/* Subtraction and negation use subtraction/borrow, conditional addition of m,
 * and masking of the zero case. Their fixed bounds are in machine.md.
 * Higher operations follow the explicit bounded algorithms in that file.
 */
