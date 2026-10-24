#include <stdio.h>
#include <stdint.h>
void key_mixed_final(const uint64_t*, uint64_t*);
int main(void){uint64_t l[6],o[3];
 while(fread(l,sizeof(uint64_t),6,stdin)==6){key_mixed_final(l,o);
  fwrite(o,sizeof(uint64_t),3,stdout);} return 0;}
