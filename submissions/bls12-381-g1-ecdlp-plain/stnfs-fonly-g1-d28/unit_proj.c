/* unit-X projection: acc += s * v[idx], 4-limb, s = +-1 two's complement (v^mask, carry-in 1
   when negative); the dense variant's bw_dot2 is 59 per element. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
static uint64_t V[1<<20][4], ACC[8][4]; static uint32_t IDX[1<<20]; static int8_t SG[1<<20];
int main(int argc,char**argv){
  int reps=argc>1?atoi(argv[1]):0; size_t N=1<<20;
  if(argc>1&&!strcmp(argv[1],"selftest")){ uint64_t a[4],v[4];
    for(int s=-1;s<=1;s+=2){ for(int j=0;j<4;j++){a[j]=0;v[j]=(uint64_t)(j*7+3);}
      uint64_t mask=(uint64_t)((int64_t)s>>63); unsigned char c=(unsigned char)(mask&1);
      for(int j=0;j<4;j++){ unsigned __int128 t=(unsigned __int128)a[j]+(v[j]^mask)+c; a[j]=(uint64_t)t; c=(unsigned char)(t>>64);}
      printf("s=%d -> %llx %llx %llx %llx\n",s,(unsigned long long)a[0],(unsigned long long)a[1],(unsigned long long)a[2],(unsigned long long)a[3]); }
    for(int j=0;j<4;j++){a[j]=0;v[j]=0;} { uint64_t mask=~0ULL; unsigned char c=1;
      for(int j=0;j<4;j++){ unsigned __int128 t=(unsigned __int128)a[j]+(v[j]^mask)+c; a[j]=(uint64_t)t; c=(unsigned char)(t>>64);}
      printf("zero - zero -> %llx %llx %llx %llx\n",(unsigned long long)a[0],(unsigned long long)a[1],(unsigned long long)a[2],(unsigned long long)a[3]); }
    return 0; }
  for(size_t i=0;i<N;i++){V[i][0]=i*3; V[i][1]=i*5; V[i][2]=i*7; V[i][3]=i*9; IDX[i]=(uint32_t)(i*2654435761u)&(N-1); SG[i]=(int8_t)((i&1)?1:-1);}
  for(int r=0;r<reps;r++){
    for(size_t i=0;i<N;i++){
      const uint64_t *v=V[IDX[i]]; uint64_t mask= (uint64_t)((int64_t)SG[i]>>63); uint64_t *a=ACC[i&7];
      unsigned char c=(unsigned char)(mask & 1);           /* two's complement: carry-in 1 when negative */
      for(int j=0;j<4;j++){ unsigned __int128 t=(unsigned __int128)a[j]+(v[j]^mask)+c; a[j]=(uint64_t)t; c=(unsigned char)(t>>64);}
    }
  }
  printf("%llu\n",(unsigned long long)ACC[0][0]+ACC[3][1]); return 0;
}
