/* canon: units of Z[y]/(h) with |coeff|<=H (mod sign), then for each input point (12 ints a0..a5 b0..b5)
 * report: a-canonical flag (Q_a(a) <= Q_a(eps a) for all units), number of units with (eps a, eps b) in E,
 * Q-min flag.  Usage: canon h0 h1 h2 h3 h4 h5 H A0 matrixfile < points  (h monic deg 6, M 12x12 rows) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
static int hc[6];
static double M[12][12];
typedef struct { long c[6]; } el;
static el *units; static int nunits;
static void mulmod(const long *a, const long *u, long *out) {
  long t[11] = {0};
  for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) t[i + j] += a[i] * u[j];
  for (int d = 10; d >= 6; d--) { long c = t[d]; if (!c) continue; t[d] = 0; for (int k = 0; k < 6; k++) t[d - 6 + k] -= c * hc[k]; }
  for (int k = 0; k < 6; k++) out[k] = t[k];
}
static double qf(const double *v, int off) {
  double s = 0; for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) s += M[off + i][off + j] * v[i] * v[j]; return s;
}
int main(int argc, char **argv) {
  for (int i = 0; i < 6; i++) hc[i] = atoi(argv[1 + i]);
  int H = atoi(argv[7]); double A0 = atof(argv[8]);
  FILE *fm = fopen(argv[9], "r");
  for (int i = 0; i < 12; i++) for (int j = 0; j < 12; j++) if (fscanf(fm, "%lf", &M[i][j]) != 1) { fprintf(stderr, "bad matrix\n"); return 1; }
  fclose(fm);
  /* roots for norm */
  double complex r[6]; for (int i = 0; i < 6; i++) r[i] = cpow(0.4 + 0.9 * I, i + 1);
  for (int it = 0; it < 2000; it++) for (int i = 0; i < 6; i++) {
    double complex num = 1; for (int k = 5; k >= 0; k--) num = num * r[i] + hc[k];
    double complex den = 1; for (int j = 0; j < 6; j++) if (j != i) den *= (r[i] - r[j]);
    r[i] -= num / den;
  }
  int cap = 1 << 16; units = malloc(sizeof(el) * cap); nunits = 0;
  long c[6];
  for (c[0] = -H; c[0] <= H; c[0]++) for (c[1] = -H; c[1] <= H; c[1]++) for (c[2] = -H; c[2] <= H; c[2]++)
  for (c[3] = -H; c[3] <= H; c[3]++) for (c[4] = -H; c[4] <= H; c[4]++) for (c[5] = -H; c[5] <= H; c[5]++) {
    double complex nrm = 1;
    for (int i = 0; i < 6; i++) { double complex v = 0; for (int k = 5; k >= 0; k--) v = v * r[i] + c[k]; nrm *= v; }
    double n = creal(nrm);
    if (fabs(n - 1) < 1e-6 || fabs(n + 1) < 1e-6) {
      /* keep one of +-c: first nonzero coefficient positive; skip 1 */
      int k = 0; while (k < 6 && c[k] == 0) k++;
      if (c[k] < 0) continue;
      if (k == 0 && c[0] == 1 && !c[1] && !c[2] && !c[3] && !c[4] && !c[5]) continue;
      memcpy(units[nunits++].c, c, sizeof(c));
    }
  }
  fprintf(stderr, "units(mod sign, excl 1): %d\n", nunits);
  long a[6], b[6], ea[6], eb[6]; double va[6], vb[6];
  while (scanf("%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld", &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &b[0], &b[1], &b[2], &b[3], &b[4], &b[5]) == 12) {
    for (int i = 0; i < 6; i++) { va[i] = a[i] / A0; vb[i] = b[i] / A0; }
    double qa = qf(va, 0), qb = qf(vb, 6);
    int canon_a = 1, canon_q = 1, m = 0;
    for (int u = 0; u < nunits; u++) {
      mulmod(a, units[u].c, ea); mulmod(b, units[u].c, eb);
      for (int i = 0; i < 6; i++) { va[i] = ea[i] / A0; vb[i] = eb[i] / A0; }
      double qea = qf(va, 0), qeb = qf(vb, 6);
      if (qea < qa - 1e-12) canon_a = 0;
      if (qea + qeb < qa + qb - 1e-12) canon_q = 0;
      if (qea + qeb <= 1.0) m++;
    }
    printf("%d %d %d\n", canon_a, canon_q, m);
  }
  return 0;
}
