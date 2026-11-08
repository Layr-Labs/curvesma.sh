/* Projection through X = U*D: U is a fixed sparse unit matrix (3 rows per column, random sign),
   D a random diagonal.  Per column element: one 64-bit multiply by D[j], then one signed add
   into each of the 3 rows of that column.  acc = (U D) v. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define M 28
static uint64_t V[1<<20][4], ACC[1<<12][4]; static uint32_t ROW[1<<20][3]; static uint32_t R1[1<<20]; static uint64_t D[1<<20]; static int8_t SG[1<<20];
int main(int argc,char**argv){
  int reps=(argc>1&&argv[1][0]>='0'&&argv[1][0]<='9')?atoi(argv[1]):(argc>2?atoi(argv[2]):0); size_t N=1<<20;
  if(argc>1&&!strcmp(argv[1],"selftest")){ uint64_t a=0,v=5,d=3,mask=~0ULL; /* v-d with d*v scaled */
    uint64_t s=(v*d)&0xffffffffffffffULL; (void)mask; printf("scale then add: s=%llu acc=%llu\n",(unsigned long long)s,(unsigned long long)(a+s)); return 0; }
  for(size_t i=0;i<N;i++){for(int q=0;q<4;q++)V[i][q]=i*2654435761u+q*7; D[i]=i*40503u+7;
    for(int k=0;k<3;k++)ROW[i][k]=(uint32_t)((i*2654435761u+k*97)&(M-1)); R1[i]=(uint32_t)((i*2654435761u)&(M-1)); SG[i]=(int8_t)((i&1)?1:-1);}
  int unit=(argc>1&&!strcmp(argv[1],"unit")); int k1=(argc>1&&!strcmp(argv[1],"k1"));
  for(int r=0;r<reps;r++){
    for(size_t j=0;j<N;j++){
      uint64_t mask=(uint64_t)((int64_t)SG[j]>>63);
      if(k1){ uint64_t *a=ACC[R1[j]]; unsigned char c=(unsigned char)(mask&1);
        for(int q=0;q<4;q++){ uint64_t t=(unit?V[j][q]:V[j][q]*D[j]); unsigned __int128 z=(unsigned __int128)a[q]+(t^mask)+c; a[q]=(uint64_t)z; c=(unsigned char)(z>>64);} }
      else for(int k=0;k<3;k++){ uint64_t *a=ACC[ROW[j][k]]; unsigned char c=(unsigned char)(mask&1);
        for(int q=0;q<4;q++){ uint64_t t=(unit?V[j][q]:V[j][q]*D[j]); unsigned __int128 z=(unsigned __int128)a[q]+(t^mask)+c; a[q]=(uint64_t)z; c=(unsigned char)(z>>64);} }
    }
  }
  printf("%llu\n",(unsigned long long)(ACC[0][0]+ACC[5][1])); return 0;
}
