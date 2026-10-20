#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#define SECP256K1_WIDEMUL_INT128 1
#define SECP256K1_INT128_NATIVE 1
#include "/secp/src/int128_impl.h"
#include "/secp/src/field_impl.h"
__attribute__((noinline)) static void reference_field_mul(secp256k1_fe *r, const secp256k1_fe *a, const secp256k1_fe *b) {
    secp256k1_fe_mul(r,a,b);
}
int main(int argc,char **argv) {
    int n=argc>1?atoi(argv[1]):4096;
    secp256k1_fe a,b;
    secp256k1_fe_set_int(&a,17);
    secp256k1_fe_set_int(&b,42);
    reference_field_mul(&a,&a,&b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for(int i=0;i<n;i++) reference_field_mul(&a,&a,&b);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%lu\n",(unsigned long)a.n[0]);
    return 0;
}
