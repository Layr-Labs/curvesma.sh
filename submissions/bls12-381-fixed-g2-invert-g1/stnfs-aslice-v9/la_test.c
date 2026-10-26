#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef struct { uint64_t v[4]; } fe_t;
typedef struct { uint64_t v[5]; } acc_t;
void row_dot(acc_t*,acc_t*,const fe_t*,const uint64_t*,const uint8_t*,long,long);
void row_finish(fe_t*,const acc_t*,const acc_t*);
static uint64_t s=88172645463325252ULL; static uint64_t rnd(void){s^=s<<13;s^=s>>7;s^=s<<17;return s;}
int main(int argc,char**argv){
  int NX=1000; fe_t *x=malloc(NX*sizeof(fe_t));
  FILE*f=fopen(argv[1],"r"); for(int i=0;i<NX;i++) for(int k=0;k<4;k++) fscanf(f,"%llx",(unsigned long long*)&x[i].v[k]); fclose(f);
  for(int row=0;row<2000;row++){
    long w = 1 + rnd()%400; long npos = rnd()%(w+1);
    uint64_t idx[400]; uint8_t val[400];
    for(long j=0;j<w;j++){ idx[j]=rnd()%NX; val[j]= (row<100)?255:(1+rnd()%255);} 
    acc_t P={{0}},N={{0}}; fe_t y; row_dot(&P,&N,x,idx,val,npos,w); row_finish(&y,&P,&N);
    printf("%ld %ld",w,npos); for(long j=0;j<w;j++) printf(" %u:%u",idx[j],val[j]);
    printf(" | %016llx%016llx%016llx%016llx\n",(unsigned long long)y.v[3],(unsigned long long)y.v[2],(unsigned long long)y.v[1],(unsigned long long)y.v[0]);
  }
}
