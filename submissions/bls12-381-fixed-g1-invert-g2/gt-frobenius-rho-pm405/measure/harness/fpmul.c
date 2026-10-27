#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):1;
    long off=argc>2?atol(argv[2]):0;
    blst_fp a,b;
    uint64_t av[6]={17},bv[6]={42};
    blst_fp_from_uint64(&a,av);
    blst_fp_from_uint64(&b,bv);
    for(long i=0;i<off;i++) blst_fp_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) blst_fp_mul(&a,&a,&b);
    CALLGRIND_TOGGLE_COLLECT;
    printf("fp %lu\n",(unsigned long)a.l[0]);
    return 0;
}
