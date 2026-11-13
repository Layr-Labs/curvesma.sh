#include "fp12.h"
#include <stdio.h>
#include <stdlib.h>
int main(){
  srand(7);
  int bad=0;
  for (int t=0;t<200;t++){
    fp2 x,y,z,r;
    for(int i=0;i<4;i++){ x.a.v[i]=(((uint64_t)rand()<<40)^(uint64_t)rand())&0x0fffffffffffffffULL;
      x.b.v[i]=(((uint64_t)rand()<<40)^(uint64_t)rand())&0x0fffffffffffffffULL;
      y.a.v[i]=(((uint64_t)rand()<<40)^(uint64_t)rand())&0x0fffffffffffffffULL;
      y.b.v[i]=(((uint64_t)rand()<<40)^(uint64_t)rand())&0x0fffffffffffffffULL; }
    fp2_mul(&z,&x,&y);
    fp t0,t1,t2,sa,sb;
    fp_mul(&t0,&x.a,&y.a); fp_mul(&t1,&x.b,&y.b);
    fp_add(&sa,&x.a,&x.b); fp_add(&sb,&y.a,&y.b); fp_mul(&t2,&sa,&sb);
    fp_sub(&r.a,&t0,&t1);
    fp_sub(&t2,&t2,&t0); fp_sub(&r.b,&t2,&t1);
    if(!(fp_eq(&z.a,&r.a)&&fp_eq(&z.b,&r.b))){
      bad++;
      if(bad<3){
        printf("t%d a:",t); for(int i=0;i<4;i++)printf(" %016lx",z.a.v[i]);
        printf(" ref:"); for(int i=0;i<4;i++)printf(" %016lx",r.a.v[i]);
        printf("  b:"); for(int i=0;i<4;i++)printf(" %016lx",z.b.v[i]);
        printf(" ref:"); for(int i=0;i<4;i++)printf(" %016lx",r.b.v[i]); printf("\n");
      }
    }
  }
  printf("bad=%d/200\n",bad);
  return bad!=0;
}
