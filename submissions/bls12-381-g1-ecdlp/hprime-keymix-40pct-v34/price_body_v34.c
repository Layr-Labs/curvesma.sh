#include <stdint.h>
typedef struct { uint64_t l[6]; } fe381;          
typedef struct { uint64_t k[3]; uint64_t gap[2]; uint32_t cursor, bound; } rec_t;
__attribute__((noinline)) uint32_t op_walk_cursor_and_bound(uint32_t *cursor, uint32_t bound,
                                                            uint32_t step) {
    uint32_t c = *cursor + step;                  
    if (c >= bound) c -= bound;                   
    *cursor = c;
    return (c < bound) ? 1u : 0u;                 
}
__attribute__((noinline)) uint64_t op_gap_load_2limb(const uint64_t *gap) {
    return gap[0] ^ (gap[1] << 1);                
}
__attribute__((noinline)) uint64_t op_gap_index(uint64_t gap, uint64_t mask) {
    return (gap * 0x9E3779B97F4A7C15ULL) & mask;  
}
__attribute__((noinline)) void op_table_point_load(const fe381 *tab, uint64_t i, fe381 *out) {
    out->l[0] = tab[i].l[0]; out->l[1] = tab[i].l[1]; out->l[2] = tab[i].l[2];
    out->l[3] = tab[i].l[3]; out->l[4] = tab[i].l[4]; out->l[5] = tab[i].l[5];
}
__attribute__((noinline)) void op_point_store_into_record(rec_t *r, const fe381 *p) {
    r->k[0] = p->l[0]; r->k[1] = p->l[1]; r->k[2] = p->l[2];
    r->gap[0] = p->l[3]; r->gap[1] = p->l[4] + p->l[5];
}
__attribute__((noinline)) uint32_t op_walk_loop_control(uint32_t i, uint32_t n) {
    return (i + 1u < n) ? (i + 1u) : 0u;          
}
__attribute__((noinline)) void op_p_minus_y(const fe381 *y, const fe381 *p, fe381 *out) {
    uint64_t borrow = 1;                          
    for (int i = 0; i < 6; i++) {
        uint64_t t = (~y->l[i]) + (i == 0 ? 1u : 0u);   
        uint64_t s = t + (borrow & (i == 0 ? 1u : 0u));
        out->l[i] = p->l[i] + s;
    }
}
__attribute__((noinline)) void op_even_rep_select(const fe381 *y, const fe381 *pmy, fe381 *out) {
    uint64_t odd = y->l[0] & 1u;                   
    for (int i = 0; i < 6; i++) out->l[i] = odd ? pmy->l[i] : y->l[i];
}
__attribute__((noinline)) uint64_t op_route_4hops_x_21bits(uint64_t key, const uint32_t *tab) {
    uint64_t h = 0;                                
    for (int j = 0; j < 4; j++) {
        uint32_t d = (uint32_t)((key >> (21 * j)) & 0x1FFFFFu);
        h += tab[d & 0x3FFu];
    }
    return h;
}
__attribute__((noinline)) void op_sort_3_passes_x_25(uint64_t *k, int n, uint64_t *scratch) {
    for (int pass = 0; pass < 3; pass++) {         
        for (int i = 0; i < n; i++) scratch[i] = k[i];
        for (int i = 1; i < n; i++) {              
            uint64_t v = scratch[i]; int j = i - 1;
            while (j >= 0 && ((scratch[j] >> (25 * pass)) & 0x1FFFFFFULL) >
                             ((v >> (25 * pass)) & 0x1FFFFFFULL)) {
                scratch[j + 1] = scratch[j]; j--;
            }
            scratch[j + 1] = v;
        }
        for (int i = 0; i < n; i++) k[i] = scratch[i];
    }
}
__attribute__((noinline)) uint32_t op_group_scan(const uint64_t *k, int n, uint64_t target) {
    uint32_t hits = 0;
    for (int i = 0; i < n; i++) hits += (k[i] == target);   
    return hits;
}
#define M381_LIMBS 6
__attribute__((noinline)) void op_fp381_mulmod(const fe381 *a, const fe381 *b, fe381 *out) {
    unsigned __int128 acc[M381_LIMBS];
    for (int i = 0; i < M381_LIMBS; i++) acc[i] = 0;
    for (int i = 0; i < M381_LIMBS; i++) {                 
        unsigned __int128 carry = 0;
        for (int j = 0; j < M381_LIMBS; j++) {
            unsigned __int128 t = (unsigned __int128)a->l[i] * b->l[j] + acc[j] + carry;
            acc[j] = t & 0xFFFFFFFFFFFFFFFFull;
            carry = t >> 64;
        }
        acc[i] += carry;
    }
    for (int i = M381_LIMBS - 1; i >= 0; i--) out->l[i] = (uint64_t)acc[i];  
}
__attribute__((noinline)) void op_subgroup_membership(const fe381 *x, uint64_t e, fe381 *out) {
    fe381 base = *x, acc = *x;
    for (int bit = 62; bit >= 0; bit--) {
        op_fp381_mulmod(&acc, &acc, &acc);
        if ((e >> bit) & 1u) op_fp381_mulmod(&acc, &base, &acc);
    }
    *out = acc;
}
__attribute__((noinline)) uint32_t op_prime_screen(uint64_t cand) {
    static const uint32_t small[16] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59};
    uint32_t rem = 0;
    for (int i = 0; i < 16; i++) rem += (uint32_t)(cand % small[i]);
    return rem;
}
__attribute__((noinline)) void op_point_construct(const fe381 *p, const fe381 *q, fe381 *out) {
    fe381 t, u;
    op_fp381_mulmod(p, q, &t);
    op_fp381_mulmod(&t, p, &u);
    op_fp381_mulmod(&u, q, &t);
    *out = t;
    for (int i = 0; i < 8; i++) op_fp381_mulmod(out, out, out);   
}
#ifdef PRICE_COUNTER_MAIN
#include <stdio.h>
int main(void) {
    uint32_t cursor = 0, bound = 1u << 18; uint64_t gap[2] = {12345, 6789};
    fe381 tab[8] = {{0}}, p = {{~0ull, ~0ull, ~0ull, ~0ull, ~0ull, 0x1aull}}, y, out;
    rec_t r; uint64_t key = 0xDEADBEEFCAFEBABEULL; uint32_t rtab[1024] = {0};
    uint64_t ks[16] = {0}, scratch[16] = {0};
    long n = 1000000, acc = 0;
    for (long i = 0; i < n; i++) {
        acc += op_walk_cursor_and_bound(&cursor, bound, 3);
        acc += (long)op_gap_load_2limb(gap);
        acc += (long)op_gap_index(key + i, 0xFFFFF);
        op_table_point_load(tab, i & 7, &out);
        op_point_store_into_record(&r, &out);
        acc += op_walk_loop_control((uint32_t)i & 15u, 16u);
        y = out; op_p_minus_y(&y, &p, &out);
        op_even_rep_select(&y, &out, &out);
        acc += (long)op_route_4hops_x_21bits(key ^ (uint64_t)i, rtab);
        ks[i & 15] = key + i; op_sort_3_passes_x_25(ks, 16, scratch);
        acc += op_group_scan(ks, 16, key + i);
    }
    printf("calls_per_op=%ld checksum=%ld cursor=%u\n", n, acc, cursor);
    return 0;
}
#endif
