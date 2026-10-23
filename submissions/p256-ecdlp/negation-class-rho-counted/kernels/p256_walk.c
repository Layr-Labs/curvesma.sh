/* P-256 batch round of the class walk, for the p256-ecdlp cost ledger.
 *
 * 4 x 64-bit limbs, Montgomery form (R = 2^256), freestanding C: no libc, no
 * vector code.  The file has two uses and nothing else:
 *   1. correctness of the arithmetic and of the round against Python big
 *      integers and affine P-256 arithmetic (kernels/test_native.py loads it
 *      through ctypes), and
 *   2. instruction counts.  The field kernels are straight-line code (the
 *      reduction rounds are written out); pass1, pass2, pass1_dbl, pass2_dbl and
 *      fe_inv are loops, and their cost is the longest path of one iteration
 *      through the compiled control-flow graph (computed by experiment.py).
 *
 * Build for the count (flags of the reference measurements in
 * tracks/docs/COST_MODEL.md, except that the compiler is clang, see proof.md):
 *   clang -target aarch64-linux-gnu -O3 -fno-tree-vectorize -mgeneral-regs-only \
 *         -ffreestanding -fno-builtin -fPIC -S -o p256_walk.s p256_walk.c
 *
 * p = 2^256 - 2^224 + 2^192 + 2^96 - 1 = -1 mod 2^64, so -p^-1 = 1 mod 2^64 and
 * the Montgomery quotient digit of each reduction round is the limb itself. */
typedef unsigned long long u64;
typedef unsigned __int128 u128;
typedef struct { u64 l[4]; } fe;

#define P0 0xFFFFFFFFFFFFFFFFull
#define P1 0x00000000FFFFFFFFull
#define P3 0xFFFFFFFF00000001ull
#define AI static inline __attribute__((always_inline))

#define ADC(x, y, c, co) __builtin_addcll((x), (y), (c), &(co))
#define SBB(x, y, c, bo) __builtin_subcll((x), (y), (c), &(bo))
AI u64 mulhi(u64 a, u64 b) { return (u64)(((u128)a * b) >> 64); }

/* column accumulate: (r2:r1:r0) += a * b */
#define MA(a, b) { u64 lo_ = (a) * (b), hi_ = mulhi((a), (b)), c_;  \
                   r0 = ADC(r0, lo_, 0, c_); r1 = ADC(r1, hi_, c_, c_); r2 = ADC(r2, 0, c_, c_); }
#define NEXT(t) { t = r0; r0 = r1; r1 = r2; r2 = 0; }

/* r = a * b * 2^-256 mod p, fully reduced to [0, p). r may alias a, b.
 * Product by columns (the 16 limb products, carries in adc chains), then four
 * Montgomery rounds.  In round i, m = t[i] and t += m * p * 2^(64 i): the
 * limbs of p are (2^64 - 1, 2^32 - 1, 0, 2^64 - 2^32 + 1), so limb i becomes 0,
 * and with the carry m of the first limb the next limb gets m << 32 (a 96-bit
 * quantity, two limbs), limb i+3 gets m * (2^64 - 2^32 + 1) whose low half is
 * m - (m << 32) and whose high half is one `umulh`. */
AI void mulm(fe *r, const fe *a, const fe *b)
{
    u64 a0 = a->l[0], a1 = a->l[1], a2 = a->l[2], a3 = a->l[3];
    u64 b0 = b->l[0], b1 = b->l[1], b2 = b->l[2], b3 = b->l[3];
    u64 t[8], r0 = 0, r1 = 0, r2 = 0;

    MA(a0, b0);                                     NEXT(t[0])
    MA(a0, b1); MA(a1, b0);                         NEXT(t[1])
    MA(a0, b2); MA(a1, b1); MA(a2, b0);             NEXT(t[2])
    MA(a0, b3); MA(a1, b2); MA(a2, b1); MA(a3, b0); NEXT(t[3])
    MA(a1, b3); MA(a2, b2); MA(a3, b1);             NEXT(t[4])
    MA(a2, b3); MA(a3, b2);                         NEXT(t[5])
    MA(a3, b3);                                     t[6] = r0; t[7] = r1;

    u64 top = 0;
#define ROUND(i, PROP)                                                         \
    {                                                                          \
        u64 c, m = t[i], lo = m << 32, hi = m >> 32;                           \
        u64 lo3 = m - lo, hi3 = mulhi(m, P3);                                  \
        t[(i) + 1] = ADC(t[(i) + 1], lo, 0, c);                                \
        t[(i) + 2] = ADC(t[(i) + 2], hi, c, c);                                \
        t[(i) + 3] = ADC(t[(i) + 3], lo3, c, c);                               \
        t[(i) + 4] = ADC(t[(i) + 4], hi3, c, c);                               \
        PROP                                                                   \
    }
#define PC(j) t[j] = ADC(t[j], 0, c, c);
    ROUND(0, PC(5) PC(6) PC(7) top += c;)
    ROUND(1, PC(6) PC(7) top += c;)
    ROUND(2, PC(7) top += c;)
    ROUND(3, top += c;)
#undef ROUND
#undef PC

    /* t[4..7] + top * 2^256 < 2p: subtract p once if needed */
    u64 bw, d0, d1, d2, d3;
    d0 = SBB(t[4], P0, 0, bw);
    d1 = SBB(t[5], P1, bw, bw);
    d2 = SBB(t[6], 0, bw, bw);
    d3 = SBB(t[7], P3, bw, bw);
    u64 keep = bw & ~top;
    r->l[0] = keep ? t[4] : d0;
    r->l[1] = keep ? t[5] : d1;
    r->l[2] = keep ? t[6] : d2;
    r->l[3] = keep ? t[7] : d3;
}

/* r = a - b mod p, inputs in [0, p) */
AI void subm(fe *r, const fe *a, const fe *b)
{
    u64 bw, c, d0, d1, d2, d3;
    d0 = SBB(a->l[0], b->l[0], 0, bw);
    d1 = SBB(a->l[1], b->l[1], bw, bw);
    d2 = SBB(a->l[2], b->l[2], bw, bw);
    d3 = SBB(a->l[3], b->l[3], bw, bw);
    u64 mk = 0 - bw;                                /* all ones if a < b */
    r->l[0] = ADC(d0, P0 & mk, 0, c);
    r->l[1] = ADC(d1, P1 & mk, c, c);
    r->l[2] = ADC(d2, 0, c, c);
    r->l[3] = ADC(d3, P3 & mk, c, c);
}

/* r = a + b mod p, inputs in [0, p) */
AI void addm(fe *r, const fe *a, const fe *b)
{
    u64 c, bw, s0, s1, s2, s3, d0, d1, d2, d3;
    s0 = ADC(a->l[0], b->l[0], 0, c);
    s1 = ADC(a->l[1], b->l[1], c, c);
    s2 = ADC(a->l[2], b->l[2], c, c);
    s3 = ADC(a->l[3], b->l[3], c, c);
    d0 = SBB(s0, P0, 0, bw);
    d1 = SBB(s1, P1, bw, bw);
    d2 = SBB(s2, 0, bw, bw);
    d3 = SBB(s3, P3, bw, bw);
    u64 keep = bw & ~c;
    r->l[0] = keep ? s0 : d0;
    r->l[1] = keep ? s1 : d1;
    r->l[2] = keep ? s2 : d2;
    r->l[3] = keep ? s3 : d3;
}

/* exported single-operation entry points (standalone counts, ctypes tests) */
void fe_mul(fe *r, const fe *a, const fe *b) { mulm(r, a, b); }
void fe_add(fe *r, const fe *a, const fe *b) { addm(r, a, b); }
void fe_sub(fe *r, const fe *a, const fe *b) { subm(r, a, b); }

/* r = a^(p-2) (Fermat), left-to-right binary: 255 squarings, 127 multiplications */
void fe_inv(fe *r, const fe *a)
{
    /* p - 2 = 2^256 - 2^224 + 2^192 + 2^96 - 3, limbs little endian */
    static const u64 E[4] = { 0xFFFFFFFFFFFFFFFDull, 0x00000000FFFFFFFFull,
                              0x0000000000000000ull, 0xFFFFFFFF00000001ull };
    fe acc = *a, base = *a;
    for (int bit = 254; bit >= 0; bit--) {
        mulm(&acc, &acc, &acc);
        if ((E[bit >> 6] >> (bit & 63)) & 1) mulm(&acc, &acc, &base);
    }
    *r = acc;
}

/* ---- walk state and batched step ------------------------------------- */

#define RB 11                       /* R = 2^11 table partitions */
typedef struct { fe x, y; } pt;     /* affine point, Montgomery coordinates */

typedef struct {
    pt P;                           /* current member of the class */
    u64 km1[3], km2[3], kcp[3];     /* fp history: limbs 0, 2, 3 of x */
    u64 len_lo, len_hi;             /* operations so far, a 128-bit counter */
    u64 flags;                      /* bit 0: pending escape or traversal: every arrival is an event */
    u64 idx;                        /* table entry 2h+t chosen in pass 1 */
} walker;

/* One round of the batch schedule (proof.md section 5.1).  A walker is either
 * in the addition list W[0..n) or, when its last arrival ended in a detected
 * cycle, in the doubling queue Q[0..nd) (the slow path puts it there).  Both
 * lists share one running product and one inversion:
 *
 *   acc = 1;  pass1(W) ; pass1_dbl(Q) ;  inv = acc^-1 ;
 *   pass2_dbl(Q) ; pass2(W)
 *
 * pass1 / pass2 / pass1_dbl / pass2_dbl are loops over the bodies below, and
 * they are what is counted (iteration cost = longest path through the loop).
 *
 * body_a: sign test t = [p - y < y], table address 2h+t, dx = x - xT,
 *         zero-denominator test, running product.
 * body_b: peel the inverse (two products), lambda, lambda^2, y3, then the
 *         arrival checks of the new class (distinguished, length cap,
 *         pending escape, 2-cycle, checkpoint) and the history update.  Any
 *         event goes through one combined branch to the separately priced
 *         slow path, which owns the walker's history from then on. */
extern void slow_path(walker *w);   /* rare events: priced separately in proof.md */
#ifdef WITH_SLOW_PATH_STUB
void slow_path(walker *w) { w->flags |= 2; }  /* test builds only */
#endif

/* arrival processing of the new class with key x3; returns 1 if an event was
 * handed to the slow path.  The operation count is 128 bits wide because the
 * production cap is L_max = 20 * 2^64 (69 bits): lmax_hi = 20 means L_max =
 * lmax_hi * 2^64, and len >= L_max is len_hi >= lmax_hi. */
AI u64 arrival(walker *w, const fe *x3, u64 lmax_hi)
{
    u64 f0 = x3->l[0], f2 = x3->l[2], f3 = x3->l[3];
    u64 lo = w->len_lo + 1;
    u64 hi = w->len_hi + (lo == 0);
    w->len_lo = lo; w->len_hi = hi;
    u64 ev = (x3->l[1] == 0);                       /* distinguished */
    ev |= (hi >= lmax_hi);                          /* trail cap L_max = lmax_hi * 2^64 */
    ev |= (w->flags & 1);                           /* pending escape */
    ev |= (f0 == w->km2[0]) & (f2 == w->km2[1]) & (f3 == w->km2[2]);
    ev |= (f0 == w->kcp[0]) & (f2 == w->kcp[1]) & (f3 == w->kcp[2]);
    if (__builtin_expect(ev, 0)) { slow_path(w); return 1; }
    if ((lo & 63) == 0) { w->kcp[0] = f0; w->kcp[1] = f2; w->kcp[2] = f3; }
    w->km2[0] = w->km1[0]; w->km2[1] = w->km1[1]; w->km2[2] = w->km1[2];
    w->km1[0] = f0;        w->km1[1] = f2;        w->km1[2] = f3;
    return 0;
}

AI void body_a(walker *w, const pt *T, fe *dx_out, fe *pre_out, fe *acc, u64 *bad)
{
    u64 h = w->P.x.l[0] & ((1u << RB) - 1);
    u128 d;
    u64 y0 = w->P.y.l[0], y1 = w->P.y.l[1], y2 = w->P.y.l[2], y3 = w->P.y.l[3];
    u64 n0, n1, n2, n3, b;
    d = (u128)P0 - y0;                     n0 = (u64)d; b = (u64)(d >> 64) & 1;
    d = (u128)P1 - y1 - b;                 n1 = (u64)d; b = (u64)(d >> 64) & 1;
    d = (u128)0  - y2 - b;                 n2 = (u64)d; b = (u64)(d >> 64) & 1;
    d = (u128)P3 - y3 - b;                 n3 = (u64)d;
    /* (p - y) < y  <=>  borrow out of (p - y) - y */
    d = (u128)n0 - y0;                     b = (u64)(d >> 64) & 1;
    d = (u128)n1 - y1 - b;                 b = (u64)(d >> 64) & 1;
    d = (u128)n2 - y2 - b;                 b = (u64)(d >> 64) & 1;
    d = (u128)n3 - y3 - b;                 b = (u64)(d >> 64) & 1;
    u64 idx = (h << 1) | b;
    w->idx = idx;
    fe dx;
    subm(&dx, &w->P.x, &T[idx].x);
    *bad += ((dx.l[0] | dx.l[1] | dx.l[2] | dx.l[3]) == 0);
    *dx_out = dx;
    mulm(acc, acc, &dx);
    *pre_out = *acc;
}

AI u64 body_b(walker *w, const pt *T, const fe *dx, const fe *pre_prev, fe *inv_acc, u64 lmax_hi)
{
    const pt *t = &T[w->idx];
    fe inv, dy, lam, l2, x3, d, y3;
    mulm(&inv, inv_acc, pre_prev);
    mulm(inv_acc, inv_acc, dx);
    subm(&dy, &w->P.y, &t->y);
    mulm(&lam, &dy, &inv);
    mulm(&l2, &lam, &lam);
    subm(&x3, &l2, &w->P.x);
    subm(&x3, &x3, &t->x);
    subm(&d, &w->P.x, &x3);
    mulm(&y3, &lam, &d);
    subm(&y3, &y3, &w->P.y);
    w->P.x = x3; w->P.y = y3;
    return arrival(w, &x3, lmax_hi);
}

/* Doubling escape (a = -3), same two-pass shape: denominator 2y, then
 * lambda = 3(x - 1)(x + 1) / (2y), x3 = lambda^2 - 2x, y3 = lambda (x - x3) - y.
 * ONE is the Montgomery form of 1.  The escape clears the history and the
 * pending flag, and the new class then arrives like any other. */
AI void body_a_dbl(walker *w, fe *dx_out, fe *pre_out, fe *acc, u64 *bad)
{
    fe dx;
    addm(&dx, &w->P.y, &w->P.y);
    *bad += ((dx.l[0] | dx.l[1] | dx.l[2] | dx.l[3]) == 0);
    *dx_out = dx;
    mulm(acc, acc, &dx);
    *pre_out = *acc;
}

AI u64 body_b_dbl(walker *w, const fe *ONE, const fe *dx, const fe *pre_prev, fe *inv_acc, u64 lmax_hi)
{
    fe inv, xm1, xp1, num, lam, l2, x3, d, y3;
    mulm(&inv, inv_acc, pre_prev);
    mulm(inv_acc, inv_acc, dx);
    subm(&xm1, &w->P.x, ONE);
    addm(&xp1, &w->P.x, ONE);
    mulm(&num, &xm1, &xp1);
    addm(&xm1, &num, &num);
    addm(&num, &xm1, &num);
    mulm(&lam, &num, &inv);
    mulm(&l2, &lam, &lam);
    subm(&x3, &l2, &w->P.x);
    subm(&x3, &x3, &w->P.x);
    subm(&d, &w->P.x, &x3);
    mulm(&y3, &lam, &d);
    subm(&y3, &y3, &w->P.y);
    w->P.x = x3; w->P.y = y3;
    w->km1[0] = w->km1[1] = w->km1[2] = 0;
    w->km2[0] = w->km2[1] = w->km2[2] = 0;
    w->kcp[0] = w->kcp[1] = w->kcp[2] = 0;
    w->flags = 0;
    return arrival(w, &x3, lmax_hi);
}

/* The batch loops.  ONE is the Montgomery form of 1 (R mod p).  pass2_dbl
 * takes PRE0, the product before its first walker: the last prefix product of
 * the addition list (ONE if that list is empty); pass2 starts from ONE. */
u64 pass1(walker *W, const pt *T, fe *DX, fe *PRE, fe *acc, u64 n)
{
    u64 bad = 0;
    for (u64 i = 0; i < n; i++) body_a(&W[i], T, &DX[i], &PRE[i], acc, &bad);
    return bad;
}
u64 pass2(walker *W, const pt *T, const fe *DX, const fe *PRE, const fe *ONE, fe *inv_acc, u64 n, u64 lmax_hi)
{
    u64 events = 0;
    for (u64 i = n; i-- > 0; ) events += body_b(&W[i], T, &DX[i], i ? &PRE[i - 1] : ONE, inv_acc, lmax_hi);
    return events;
}
u64 pass1_dbl(walker **Q, fe *DX, fe *PRE, fe *acc, u64 n)
{
    u64 bad = 0;
    for (u64 i = 0; i < n; i++) body_a_dbl(Q[i], &DX[i], &PRE[i], acc, &bad);
    return bad;
}
u64 pass2_dbl(walker **Q, const fe *DX, const fe *PRE, const fe *PRE0, const fe *ONE, fe *inv_acc, u64 n, u64 lmax_hi)
{
    u64 events = 0;
    for (u64 i = n; i-- > 0; ) events += body_b_dbl(Q[i], ONE, &DX[i], i ? &PRE[i - 1] : PRE0, inv_acc, lmax_hi);
    return events;
}

/* Slow-path pieces that touch the batch lists (counted; proof.md 5.1 step 4).
 *
 * move_to_dbl: a walker due its doubling is copied into the doubling array, its
 *   pointer is appended to the queue, and the hole in the addition array is
 *   filled by the last walker (swap-remove).
 * move_back: after its doubling the walker returns to the addition array.
 * trav_update: one step of a cycle traversal; keeps the member with the minimum
 *   192-bit fingerprint (limbs 3, 2, 0 compared as one integer) and counts the
 *   steps (more than 64 abandons the trail). */
typedef struct { u64 minfp[3]; pt member; u64 steps; } trav;

/* word-wise copy of a walker record (21 words), written out so that the
 * compiler does not turn it into a call to memcpy */
#define CP(i) d[i] = a[i];
AI void copy_walker(walker *dst, const walker *src)
{
    u64 *d = (u64 *)dst;
    const u64 *a = (const u64 *)src;
    CP(0) CP(1) CP(2) CP(3) CP(4) CP(5) CP(6) CP(7) CP(8) CP(9) CP(10)
    CP(11) CP(12) CP(13) CP(14) CP(15) CP(16) CP(17) CP(18) CP(19) CP(20)
}
#undef CP

void move_to_dbl(walker *A, u64 *na, walker *D, walker **Q, u64 *nd, u64 i)
{
    u64 k = *nd;
    walker *d = &D[k];
    copy_walker(d, &A[i]);
    Q[k] = d;
    *nd = k + 1;
    u64 last = *na - 1;
    copy_walker(&A[i], &A[last]);
    *na = last;
}

void move_back(walker *A, u64 *na, const walker *d)
{
    u64 k = *na;
    copy_walker(&A[k], d);
    *na = k + 1;
}

u64 trav_update(trav *t, const walker *w)
{
    u64 f0 = w->P.x.l[0], f2 = w->P.x.l[2], f3 = w->P.x.l[3];
    u64 less = (f3 < t->minfp[2]) | ((f3 == t->minfp[2]) & ((f2 < t->minfp[1]) |
               ((f2 == t->minfp[1]) & (f0 < t->minfp[0]))));
    if (less) {
        t->minfp[0] = f0; t->minfp[1] = f2; t->minfp[2] = f3;
        t->member = w->P;
    }
    u64 n = t->steps + 1;
    t->steps = n;
    return n > 64;
}
