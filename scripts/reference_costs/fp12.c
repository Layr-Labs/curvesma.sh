#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):4096;
    blst_fp12 a,b;
    blst_p1_affine p; blst_p2_affine q;
    blst_p1_to_affine(&p,blst_p1_generator());
    blst_p2_to_affine(&q,blst_p2_generator());
    blst_miller_loop(&a,&q,&p);
    blst_final_exp(&a,&a);
    b=a;
    blst_fp12_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) blst_fp12_mul(&a,&a,&b);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%d\n",blst_fp12_in_group(&a));
    return 0;
}
