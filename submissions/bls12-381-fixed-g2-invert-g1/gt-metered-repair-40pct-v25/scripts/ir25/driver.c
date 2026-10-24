/* driver25.c -- vector driver for the counted kernels (v25).
 *
 * Each mode reads 6-word little-endian limbs on stdin, one vector per line
 * (whitespace separated), applies one kernel, and prints the result limbs.
 * All inputs are Montgomery residues (a*2^384 mod p), so the Python bigint
 * reference is exact integer arithmetic with one conversion per element.
 *   f2mul/f2sqr/f2inv : 2 Fp2 operands (24 limbs) / 1 operand, 12 out
 *   f6mul/f6sqr/f6inv : 2 Fp6 operands (36 limbs) / 1 operand, 36 out
 *   f12mul/f12sqr/f12inv/csqr : 2 Fp12 operands (72) / 1 operand, 72 out
 */
#include <stdio.h>
#include <stdint.h>
#include "kernels.c"

static int rd(u64 *x, int n) {
    for (int i = 0; i < n; i++)
        if (scanf("%llx", (unsigned long long *)&x[i]) != 1) return 0;
    return 1;
}
static void pr(const u64 *x, int n) {
    for (int i = 0; i < n; i++) printf("%016llx%c", (unsigned long long)x[i],
                                       i == n - 1 ? '\n' : ' ');
}

int main(int argc, char **argv) {
    if (argc < 2) return 2;
    const char *m = argv[1];
    u64 in[144], out[144];
    fp2 *a2 = (fp2 *)in, *b2 = (fp2 *)(in + 12), *r2 = (fp2 *)out;
    fp6 *a6 = (fp6 *)in, *b6 = (fp6 *)(in + 36), *r6 = (fp6 *)out;
    fp12 *a12 = (fp12 *)in, *b12 = (fp12 *)(in + 72), *r12 = (fp12 *)out;
    if (m[0] == 'f' && m[1] == '2') {
        int nin = (m[2] == 'm') ? 24 : 12;
        int sq = (m[2] == 's'), inv = (m[2] == 'i' && m[3] == 'n');
        while (rd(in, nin)) {
            if (nin == 24) {
                if (sq) f2sqr(r2, a2); else if (inv) f2inv(r2, a2); else f2mul(r2, a2, b2);
            } else {
                if (sq) f2sqr(r2, a2); else f2inv(r2, a2);
            }
            pr(out, 12);
        }
        return 0;
    }
    if (m[0] == 'f' && m[1] == '6') {
        int nin = (m[2] == 'm') ? 72 : 36;
        int sq = (m[2] == 's'), inv = (m[2] == 'i');
        while (rd(in, nin)) {
            if (nin == 72) {
                if (sq) f6sqr(r6, a6); else if (inv) f6inv(r6, a6); else f6mul(r6, a6, b6);
            } else {
                if (sq) f6sqr(r6, a6); else f6inv(r6, a6);
            }
            pr(out, 36);
        }
        return 0;
    }
    /* f12* and csqr */
    int nin = (m[0] == 'c') ? 72 : ((m[3] == 'm') ? 144 : 72);
    while (rd(in, nin)) {
        if (m[0] == 'c') cyclotomic_sqr(r12, a12);
        else if (m[3] == 'm') f12mul(r12, a12, b12);
        else if (m[3] == 's') f12sqr(r12, a12);
        else f12inv(r12, a12);
        pr(out, 72);
    }
    return 0;
}
