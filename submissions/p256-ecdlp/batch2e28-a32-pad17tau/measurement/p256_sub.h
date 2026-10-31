#define P256_ADD_P_ON_BORROW(lbl)                                \
    "b.cs " lbl "f\n\t"                                          \
    "mov w9, #-1\n\t"                                            \
    "neg x10, x9\n\t"                                            \
    "subs %[t0], %[t0], #1\n\t"                                  \
    "adcs %[t1], %[t1], x9\n\t"                                  \
    "adcs %[t2], %[t2], xzr\n\t"                                 \
    "adc %[t3], %[t3], x10\n"                                    \
    lbl ":\n\t"

static inline void k_subb(fe *r, const fe *a, const fe *b) {
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
        P256_ADD_P_ON_BORROW("1")
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3)
        : [r] "r"(r->v), [a] "r"(a->v), [b] "r"(b->v)
        : "x9", "x10", "cc", "memory");
}

static inline void k_sub2(fe *r, const fe *a, const fe *b, const fe *c) {
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
        P256_ADD_P_ON_BORROW("1")
        "ldp %[u0], %[u1], [%[c]]\n\t"
        "ldp %[u2], %[u3], [%[c], #16]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        P256_ADD_P_ON_BORROW("2")
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3)
        : [r] "r"(r->v), [a] "r"(a->v), [b] "r"(b->v), [c] "r"(c->v)
        : "x9", "x10", "cc", "memory");
}

#define P256_ADD_P_MASKED                                        \
    "sbc x9, xzr, xzr\n\t"                                       \
    "lsr x10, x9, #32\n\t"                                       \
    "lsl x11, x9, #32\n\t"                                       \
    "sub x11, x11, x9\n\t"                                       \
    "adds %[t0], %[t0], x9\n\t"                                  \
    "adcs %[t1], %[t1], x10\n\t"                                 \
    "adcs %[t2], %[t2], xzr\n\t"                                 \
    "adc %[t3], %[t3], x11\n\t"

static inline void k_subb_bl(fe *r, const fe *a, const fe *b) {
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
        P256_ADD_P_MASKED
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3)
        : [r] "r"(r->v), [a] "r"(a->v), [b] "r"(b->v)
        : "x9", "x10", "x11", "cc", "memory");
}

static inline void k_sub2_bl(fe *r, const fe *a, const fe *b, const fe *c) {
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
        P256_ADD_P_MASKED
        "ldp %[u0], %[u1], [%[c]]\n\t"
        "ldp %[u2], %[u3], [%[c], #16]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        P256_ADD_P_MASKED
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3)
        : [r] "r"(r->v), [a] "r"(a->v), [b] "r"(b->v), [c] "r"(c->v)
        : "x9", "x10", "x11", "cc", "memory");
}
