/* unit-X projection: acc[j] += s * v[idx] for a 4-limb v and +-1 signs; the dense
   variant's bw_dot2 is 59 per element.  Sparse unit entries: 4 index loads, 4 addr, 4 addc. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
static uint64_t V[1<<20][4], ACC[8][4]; static uint32_t IDX[1<<20]; static int8_t SG[1<<20];
int main(int argc,char**argv){
  int reps=argc>1?atoi(argv[1]):0; size_t N=1<<20;
  for(size_t i=0;i<N;i++){V[i][0]=i*3; V[i][1]=i*5; V[i][2]=i*7; V[i][3]=i*9; IDX[i]=(uint32_t)(i*2654435761u)&(N-1); SG[i]=(int8_t)((i&1)?1:-1);}
  for(int r=0;r<reps;r++){
    for(size_t i=0;i<N;i++){
      const uint64_t *v=V[IDX[i]]; uint64_t mask= (uint64_t)((int64_t)SG[i]>>63); uint64_t *a=ACC[i&7];
      unsigned char c=0;
      for(int j=0;j<4;j++){ unsigned __int128 t=(unsigned __int128)a[j]+(v[j]^mask)+c; a[j]=(uint64_t)t; c=(unsigned char)(t>>64);}
      { unsigned __int128 t=(unsigned __int128)a[3]+mask; a[3]=(uint64_t)t; }
    }
  }
  printf("%llu\n",(unsigned long long)ACC[0][0]+ACC[3][1]); return 0;
}
