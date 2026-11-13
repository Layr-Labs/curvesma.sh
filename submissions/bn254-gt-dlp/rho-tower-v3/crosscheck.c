/* Local correctness cross-check for the optimised tower, built separately from
 * the attack binary:
 *   gcc -O3 -fno-tree-vectorize -mgeneral-regs-only -o crosscheck crosscheck.c fp12.c
 *   ./crosscheck
 * It compares fp6_mul/fp12_mul against a schoolbook reference using the same
 * fp2_mul primitive (9 fp2 products in Fp6), and checks fp_add/fp_sub
 * round-trips, on 2,000 deterministic pseudo-random pairs each.
 */
#include "fp12.h"
#include <stdio.h>
#include <string.h>

static uint64_t st = 0x123456789abcdefULL;
static uint64_t rnd(void){ st^=st<<13; st^=st>>7; st^=st<<17; return st; }
static void red(fp *z){
  for (int k=0;k<4;k++){ uint64_t d[4]; uint64_t br=0;
    for(int i=0;i<4;i++){ __uint128_t u=(__uint128_t)z->v[i]-FP_P.v[i]-br; d[i]=(uint64_t)u; br=(uint64_t)(u>>64)&1ULL; }
    if(!br) memcpy(z->v,d,sizeof d); else break; }
}
static void rndfp(fp*z){ for(int i=0;i<4;i++) z->v[i]=rnd(); z->v[3]&=0x3fffffffffffffffULL; red(z); }
static void rnd2(fp2*z){ rndfp(&z->a); rndfp(&z->b); }
static void mxi(fp2 *z, const fp2 *m){
  fp a2,a8,a9,b2,b8,b9;
  fp_add(&a2,&m->a,&m->a); fp_add(&a8,&a2,&a2); fp_add(&a8,&a8,&a8); fp_add(&a9,&a8,&m->a);
  fp_add(&b2,&m->b,&m->b); fp_add(&b8,&b2,&b2); fp_add(&b8,&b8,&b8); fp_add(&b9,&b8,&m->b);
  fp_sub(&z->a,&a9,&m->b); fp_add(&z->b,&m->a,&b9);
}
static void ref6(fp6 *z,const fp6*x,const fp6*y){
  fp2 t0,t1,t2,s0,s1,s2,w;
  fp2_mul(&t0,&x->c[0],&y->c[0]); fp2_mul(&t1,&x->c[1],&y->c[1]); fp2_mul(&t2,&x->c[2],&y->c[2]);
  fp2_mul(&s0,&x->c[1],&y->c[2]); fp2_mul(&s1,&x->c[2],&y->c[1]); fp2_add(&s0,&s0,&s1);
  mxi(&w,&s0); fp2_add(&s0,&t0,&w);
  fp2_mul(&s1,&x->c[0],&y->c[1]); fp2_mul(&s2,&x->c[1],&y->c[0]); fp2_add(&s1,&s1,&s2);
  mxi(&w,&t2); fp2_add(&s1,&s1,&w);
  fp2_mul(&s2,&x->c[0],&y->c[2]); fp2_add(&s2,&s2,&t1); fp2_mul(&w,&x->c[2],&y->c[0]); fp2_add(&s2,&s2,&w);
  z->c[0]=s0; z->c[1]=s1; z->c[2]=s2;
}
static void ref12(fp12*z,const fp12*x,const fp12*y){
  fp6 t0,t1,t2,a,b,v1,w;
  ref6(&t0,&x->c[0],&y->c[0]); ref6(&t1,&x->c[1],&y->c[1]);
  fp6_add(&a,&x->c[0],&x->c[1]); fp6_add(&b,&y->c[0],&y->c[1]); ref6(&t2,&a,&b);
  mxi(&w.c[0],&t1.c[2]); v1.c[2]=t1.c[1]; v1.c[1]=t1.c[0]; v1.c[0]=w.c[0];
  fp6_add(&z->c[0],&t0,&v1);
  fp6_sub(&t2,&t2,&t0); fp6_sub(&z->c[1],&t2,&t1);
}
int main(void){
  for(int it=0; it<2000; it++){
    fp6 x,y,o,r; for(int i=0;i<3;i++){ rnd2(&x.c[i]); rnd2(&y.c[i]); }
    fp6_mul(&o,&x,&y); ref6(&r,&x,&y);
    for(int i=0;i<3;i++) if(!fp_eq(&o.c[i].a,&r.c[i].a)||!fp_eq(&o.c[i].b,&r.c[i].b)){ printf("fp6 mismatch %d\n",it); return 1; }
    fp12 X,Y,O,R; for(int i=0;i<2;i++) for(int j=0;j<3;j++){ rnd2(&X.c[i].c[j]); rnd2(&Y.c[i].c[j]); }
    fp12_mul(&O,&X,&Y); ref12(&R,&X,&Y);
    for(int i=0;i<2;i++) for(int j=0;j<3;j++) if(!fp_eq(&O.c[i].c[j].a,&R.c[i].c[j].a)||!fp_eq(&O.c[i].c[j].b,&R.c[i].c[j].b)){ printf("fp12 mismatch %d\n",it); return 1; }
  }
  for(int it=0; it<2000; it++){ fp a,b,c,s; rndfp(&a); rndfp(&b); fp_add(&s,&a,&b); fp_sub(&c,&s,&b); if(!fp_eq(&c,&a)){ printf("add/sub mismatch\n"); return 1; } }
  printf("fp6/fp12 vs schoolbook: 2000 pairs ok\nadd/sub round-trip: 2000 ok\n");
  return 0;
}
