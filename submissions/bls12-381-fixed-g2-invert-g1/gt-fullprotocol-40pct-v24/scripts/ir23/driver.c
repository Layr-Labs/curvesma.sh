/* driver.c -- vector driver for the counted kernels (v23).
 *   mode "mul": reads 24 limbs (a0,a1,b0,b1 x 6), prints f2mul in Montgomery form
 *   mode "sqr": reads 72 limbs (12 Fp2 components x 6), prints cyclotomic_sqr
 * All values are 6x64-bit little-endian limbs, all inputs are Montgomery residues
 * (a*2^384 mod p) so the exact bigint check is a plain integer computation.
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
    if (argv[1][0] == 'm') {
        fp2 a, b, r;
        while (rd(&a.c[0][0], 6) && rd(&a.c[1][0], 6) && rd(&b.c[0][0], 6) &&
               rd(&b.c[1][0], 6)) {
            f2mul(&r, &a, &b);
            pr(&r.c[0][0], 6); pr(&r.c[1][0], 6);
        }
        return 0;
    }
    fp12 a, r;
    u64 *w = &a.c[0].c[0].c[0][0];
    while (rd(w, 72)) {
        cyclotomic_sqr(&r, &a);
        pr(&r.c[0].c[0].c[0][0], 72);
    }
    return 0;
}
