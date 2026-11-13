/* Per-call instruction microbenchmark for the audit: worst and best arm. */
#include "fp12.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static fp A, B, Z;
static void setp(fp *z, uint64_t w0, uint64_t w1, uint64_t w2, uint64_t w3) {
  z->v[0]=w0; z->v[1]=w1; z->v[2]=w2; z->v[3]=w3;
}
int main(int argc, char **argv) {
  long n = (argc >= 3) ? atol(argv[2]) : 100000;
  int mode = atoi(argv[1]);
  const uint64_t P0=0x3c208c16d87cfd47ULL,P1=0x97816a916871ca8dULL,
                 P2=0xb85045b68181585dULL,P3=0x30644e72e131a029ULL;
  uint64_t acc = 0;
  if (mode == 0) {                 /* fp_add: s = p exactly -> reduction + full compare chain */
    setp(&A, P0-1, P1, P2, P3); setp(&B, 1, 0, 0, 0);
    for (long i=0;i<n;i++) { fp_add(&Z,&A,&B); acc ^= Z.v[0]; }
  } else if (mode == 1) {          /* fp_add: no carry, short chain */
    setp(&A, 1, 0, 0, 0); setp(&B, 2, 0, 0, 0);
    for (long i=0;i<n;i++) { fp_add(&Z,&A,&B); acc ^= Z.v[0]; }
  } else if (mode == 2) {          /* fp_sub: borrow -> add p (long arm) */
    setp(&A, 0, 0, 0, 0); setp(&B, 1, 0, 0, 0);
    for (long i=0;i<n;i++) { fp_sub(&Z,&A,&B); acc ^= Z.v[0]; }
  } else if (mode == 3) {          /* fp_sub: no borrow (short arm) */
    setp(&A, 5, 0, 0, 0); setp(&B, 2, 0, 0, 0);
    for (long i=0;i<n;i++) { fp_sub(&Z,&A,&B); acc ^= Z.v[0]; }
  } else if (mode == 4) {          /* fp_mul: dense operands near p */
    setp(&A, P0-1, P1-1, P2-1, P3-1); setp(&B, P0-3, P1-3, P2-3, P3-3);
    for (long i=0;i<n;i++) { fp_mul(&Z,&A,&B); acc ^= Z.v[0]; }
  } else if (mode == 5) {          /* fp_mul: small operands */
    setp(&A, 3, 0, 0, 0); setp(&B, 7, 0, 0, 0);
    for (long i=0;i<n;i++) { fp_mul(&Z,&A,&B); acc ^= Z.v[0]; }
  } else if (mode == 6) {          /* empty loop baseline */
    for (long i=0;i<n;i++) { acc ^= (uint64_t)i; }
  } else if (mode == 7) {          /* fp2_mul */
    fp2 x = {A,B}, y = {B,A}, z;
    setp(&A, P0-1, P1-1, P2-1, P3-1); setp(&B, P0-3, P1-3, P2-3, P3-3);
    x.a=A; x.b=B; y.a=B; y.b=A;
    for (long i=0;i<n;i++) { fp2_mul(&z,&x,&y); acc ^= z.a.v[0]; }
  }
  printf("%llu\n", (unsigned long long)acc);
  return 0;
}
