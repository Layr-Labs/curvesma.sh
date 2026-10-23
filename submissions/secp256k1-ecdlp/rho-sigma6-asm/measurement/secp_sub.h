/* Branchy secp256k1 subtractions for reduced inputs (free register allocation).
 * On borrow, a - b + 2^256 is corrected by subtracting c = 2^32 + 977 (adds p mod 2^256);
 * the result is then in (0, p) with no further borrow. */
#define SECP_SUB_C_ON_BORROW(lbl)                                \
    "b.cs " lbl "f\n\t"                                          \
    "mov x9, #0x3d1\n\t"                                         \
    "movk x9, #1, lsl #32\n\t"                                   \
    "subs %[t0], %[t0], x9\n\t"                                  \
    "sbcs %[t1], %[t1], xzr\n\t"                                 \
    "sbcs %[t2], %[t2], xzr\n\t"                                 \
    "sbc %[t3], %[t3], xzr\n"                                    \
    lbl ":\n\t"

/* r = a - b mod p */
static inline void s_subb(fe *r, const fe *a, const fe *b) {
    uint64_t t0, t1, t2, t3, u0, u1, u2, u3;
    __asm__ volatile(
        "ldp %[t0], %[t1], [%[a]]\n\t"
        "ldp %[u0], %[u1], [%[b]]\n\t"
        "ldp %[t2], %[t3], [%[a], #16]\n\t"
        "ldp %[u2], %[u3], [%[b], #16]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        SECP_SUB_C_ON_BORROW("1")
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3)
        : [r] "r"(r->v), [a] "r"(a->v), [b] "r"(b->v)
        : "x9", "cc", "memory");
}

/* r = a - b - c mod p */
static inline void s_sub2(fe *r, const fe *a, const fe *b, const fe *c) {
    uint64_t t0, t1, t2, t3, u0, u1, u2, u3;
    __asm__ volatile(
        "ldp %[t0], %[t1], [%[a]]\n\t"
        "ldp %[u0], %[u1], [%[b]]\n\t"
        "ldp %[t2], %[t3], [%[a], #16]\n\t"
        "ldp %[u2], %[u3], [%[b], #16]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        SECP_SUB_C_ON_BORROW("1")
        "ldp %[u0], %[u1], [%[c]]\n\t"
        "ldp %[u2], %[u3], [%[c], #16]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        SECP_SUB_C_ON_BORROW("2")
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3)
        : [r] "r"(r->v), [a] "r"(a->v), [b] "r"(b->v), [c] "r"(c->v)
        : "x9", "cc", "memory");
}
