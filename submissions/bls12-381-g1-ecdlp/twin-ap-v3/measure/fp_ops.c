/* twin-ap-v2: single exported blst Fp calls, counted exactly like scripts/reference_costs/bls_field.c
   (a = op(a, b) n times; price = (Ir(8192) - Ir(4096)) / 4096). Mode "mul" is the organizers' harness. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
int main(int argc, char **argv)
{
    if (argc < 3) return 1;
    const char *op = argv[1];
    int n = atoi(argv[2]);
    blst_fp a, b;
    uint64_t av[6] = {17}, bv[6] = {42};
    blst_fp_from_uint64(&a, av);
    blst_fp_from_uint64(&b, bv);
    blst_fp_mul(&a, &a, &b);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    if (!strcmp(op, "mul")) for (int i = 0; i < n; i++) blst_fp_mul(&a, &a, &b);
    if (!strcmp(op, "sqr")) for (int i = 0; i < n; i++) blst_fp_sqr(&a, &a);
    if (!strcmp(op, "add")) for (int i = 0; i < n; i++) blst_fp_add(&a, &a, &b);
    if (!strcmp(op, "sub")) for (int i = 0; i < n; i++) blst_fp_sub(&a, &a, &b);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%lu\n", (unsigned long)a.l[0]);
    return 0;
}
