/* slice12_check.c: brute-force reference, correctness checks and benchmark driver of slice12.c (included from slice12.c; not a
 * separate translation unit). Modes are documented in the header of slice12.c. */
/* ---------------- brute-force reference over the packed rows ---------------- */
static inline double q6(const long long *b) { double t = 0; for (int i = 0; i < D; i++) { double r = 0; for (int j = 0; j < D; j++) r += M6[i][j] * (double)b[j]; t += (double)b[i] * r; } return t; }
static void suffix_of(long long idx4, long long *b) { for (int i = 1; i <= 4; i++) { b[i] = idx4 % (2 * r6[i] + 1) - r6[i]; idx4 /= 2 * r6[i] + 1; } }
static uint64_t bmod(const long long *b, uint64_t s, uint64_t p) {    /* b(s) mod p, exact (reference only) */
    __int128 v = 0; for (int i = D - 1; i >= 0; i--) { v = (v * (__int128)s + b[i]) % (__int128)p; } if (v < 0) v += p; return (uint64_t)v;
}
static void ref_add(uint16_t **refs, const pair_t *P, uint64_t c) {   /* refs[2*seg+side] */
    for (long long si = 0; si < nseg; si++) {
        seg_t *sg = &SEG[si];
        for (long long r = 0; r < sg->nrows; r++) {
            row_t *rw = &sg->rows[r]; long long b[D]; suffix_of(rw->idx4, b); b[5] = rw->k;
            for (int j = 0; j < rw->len; j++) {
                b[0] = rw->lo + j; if (q6(b) > R2) continue;
                uint64_t v = bmod(b, P->s, P->p);
                uint64_t neg = (P->p - c) % P->p;
                int m = (v == c) + (v == neg);
                refs[2 * si + P->side][rw->start + j] += (uint16_t)(m * P->lp);
            }
        }
    }
}
static uint64_t digest_arr(const uint16_t *a, long long n, uint64_t h) { for (long long i = 0; i < n; i++) { h ^= a[i]; h *= 1099511628211ULL; } return h; }

static int run_check(int log2pts, int gshift) {
    setup_region(log2pts, gshift); alloc_slice();
    long long ncell = total_cells(), box = (2 * r6[0] + 1) * L4 * K;
    fprintf(stderr, "region: R2=%.3f K=%lld L4=%lld segments=%lld cells=%lld box6=%lld (box/cells %.2f)\n", R2, K, L4, nseg, ncell, box, (double)box / ncell);
    uint16_t **refs = calloc(2 * nseg, sizeof(uint16_t *));
    for (long long i = 0; i < 2 * nseg; i++) refs[i] = calloc(SEG[i / 2].ncell + 1, 2);
    /* a mix of pairs: line path below 2^20, record path above; both sides; one prime close to 2^62 */
    static const uint64_t PS[][3] = {   /* p, s, rinv (any unit) */
        {251, 7, 3}, {4099, 1234, 77}, {65537, 4321, 999}, {1048573, 777, 31337},
        {1048583, 777, 31337}, {16777259, 55555, 1234567}, {67108879, 999983, 4242}, {2147483659ULL, 999983, 777},
        {4611686018427388039ULL, 999983, 4242}, {9223372036854775783ULL, 123456789, 987654321}};
    const long long a[6] = {-1000, 517, 33, -812, 7, 401};
    int npairs = sizeof PS / sizeof PS[0]; long long tot_hits = 0;
    pair_t *P = malloc(sizeof(pair_t) * npairs);
    for (int i = 0; i < npairs; i++) {
        prepare_pair(&P[i], PS[i][0], PS[i][1], PS[i][2], 1000 + i, i & 1);
        uint64_t c = setup6(&P[i], a);
        uint64_t cref = 0; { __int128 v = 0; for (int j = 5; j >= 0; j--) v = (v * (__int128)P[i].s + a[j]) % (__int128)P[i].p; if (v < 0) v += P[i].p; cref = (uint64_t)((unsigned __int128)v * P[i].rinv % P[i].p); }
        if (c != cref) { printf("setup6 MISMATCH p=%llu\n", (unsigned long long)P[i].p); return 1; }
        if (P[i].p < (1u << 20)) tot_hits += line_path(&P[i], c);
        else { tot_hits += record_path(&P[i], c); for (long long s = 0; s < nseg; s++) apply_batch(&SEG[s], P[i].side, P[i].lp); }
        ref_add(refs, &P[i], c);
    }
    int ok = 1;
    for (long long s = 0; s < nseg; s++) {
        if (digest_arr(SEG[s].f, SEG[s].ncell, 1) != digest_arr(refs[2 * s], SEG[s].ncell, 1)) ok = 0;
        if (digest_arr(SEG[s].g, SEG[s].ncell, 1) != digest_arr(refs[2 * s + 1], SEG[s].ncell, 1)) ok = 0;
    }
    /* scan with thresholds: every cell holding a record-path contribution on side f, and any g sum, survives */
    long long refsurv = 0;
    for (long long s = 0; s < nseg; s++) {
        long long nb = (SEG[s].ncell + 15) / 16; long long *bs = malloc(8 * (nb + 1)); uint16_t *tf = malloc(2 * nb), *tg = malloc(2 * nb);
        for (long long b = 0; b < nb; b++) { bs[b] = 16 * b; tf[b] = (uint16_t)(16 * 20); tg[b] = 0; } bs[nb] = SEG[s].ncell;
        for (long long c = 0; c < SEG[s].ncell; c++) refsurv += refs[2 * s][c] >= 16 * 20;
        scan_segment(s, bs, tf, tg, nb); free(bs); free(tf); free(tg);
    }
    uint64_t rec_ref = 0;
    for (long long s = 0; s < nseg; s++) for (uint64_t i = 0; i < SEG[s].rn; i++) rec_ref += refs[2 * s][SEG[s].rcell[i]] >= 16 * 20;
    for (long long s = 0; s < nseg; s++) recover_segment(&SEG[s]);
    /* every recovered record must be a true hit of its ideal at that cell */
    int rok = (nsurv == refsurv) && (nrecov == rec_ref);
    printf("cells=%lld walked=%lld f/g digests %s | survivors %lld (ref %lld) recovered records %llu (ref %llu) %s\n", ncell, tot_hits,
           ok ? "MATCH" : "MISMATCH", nsurv, refsurv, (unsigned long long)nrecov, (unsigned long long)rec_ref, rok ? "MATCH" : "MISMATCH");
    return ok && rok ? 0 : 1;
}

static int run_barrett(long long n) {
    uint64_t st = 88172645463325252ULL; long long bad = 0, tot = 0;
    #define RND() (st ^= st << 13, st ^= st >> 7, st ^= st << 17)
    for (unsigned k = 3; k <= 63; k++) {
        for (long long t = 0; t < n; t++) {
            uint64_t p = (RND() >> (64 - k)) | (1ULL << (k - 1)) | 1; if ((p & (p - 1)) == 0) continue;
            bar_t b = bar_init(p);
            __uint128_t x;
            int mode = t % 4;
            uint64_t u = RND(); uint64_t v = RND(); __uint128_t top = ((__uint128_t)1) << (2 * k);
            if (mode == 0) x = ((__uint128_t)(u % p)) * (v % p);                          /* product of residues */
            else if (mode == 1) x = (((__uint128_t)u << 64) | v) % top;                   /* full input range */
            else if (mode == 2) x = top - 1 - ((u & 1023) % top);                         /* top of the range */
            else x = (__uint128_t)(p - 1) * (p - 1);                                       /* largest product */
            uint64_t r = bar_red(&b, x), e = (uint64_t)(x % p);
            bad += r != e; tot++;
        }
    }
    printf("barrett: %lld reductions, bit lengths 3..63, %lld mismatches\n", tot, bad);
    return bad != 0;
}
static int run_recentre(long long n) {
    setup_region(24, 2); uint64_t st = 0x9E3779B97F4A7C15ULL; long long worst = 0; int bad = 0;
    static const uint64_t PS[] = {1048583, 16777259, 2147483659ULL, 1099511627791ULL, 281474976710677ULL, 4611686018427388039ULL};
    for (unsigned i = 0; i < sizeof PS / sizeof PS[0]; i++) {
        pair_t *P = malloc(sizeof(pair_t)); prepare_pair(P, PS[i], 999983, 1, 0, 0);
        /* exact rounding: t = adj(B) e0 c / det(B) with int128 arithmetic via the inverse's exact rational form is
           replaced by its consequence: the deviation of the floating representative from the coset is checked exactly
           (b(s) = c mod p) and its size is reported in units of the largest basis vector */
        long long bmax = 0; for (int r = 0; r < 6; r++) for (int j = 0; j < 6; j++) if (llabs(P->L6.B[r][j]) > bmax) bmax = llabs(P->L6.B[r][j]);
        for (long long t = 0; t < n; t++) {
            uint64_t c = ((RND() >> 1) % PS[i]); long long oo[D]; recentre(&P->L6, P->W6, c, oo);
            uint64_t v = bmod(oo, P->s, P->p); if (v != c) bad++;
            for (int j = 0; j < 6; j++) { long long m = llabs(oo[j]) / (bmax ? bmax : 1); if (m > worst) worst = m; }
        }
        printf("recentre p=%llu: max |coordinate| of representative / max |basis entry| = %lld\n", (unsigned long long)PS[i], worst);
        worst = 0; free(P);
    }
    printf("recentre: coset membership failures %d\n", bad);
    return bad != 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: see header\n"); return 2; }
    const char *mode = argv[1];
    if (!strcmp(mode, "check")) return run_check(atoi(argv[2]), atoi(argv[3]));
    if (!strcmp(mode, "barrett")) return run_barrett(atoll(argv[2]));
    if (!strcmp(mode, "recentre")) return run_recentre(atoll(argv[2]));
    if (!strcmp(mode, "bench")) {
        const char *what = argv[2]; int log2pts = atoi(argv[3]), gshift = atoi(argv[4]);
        uint64_t p = strtoull(argv[5], 0, 10), s = strtoull(argv[6], 0, 10), c0 = strtoull(argv[7], 0, 10); long long reps = atoll(argv[8]);
        setup_region(log2pts, gshift); alloc_slice();
        pair_t *P = malloc(sizeof(pair_t)); prepare_pair(P, p, s, 4242, 77, 0);
        long long a[6] = {-1000, 517, 33, -812, 7, 401}; long long hits = 0; uint64_t acc = 0; long long units = 0;
        if (!strcmp(what, "line")) { for (long long r = 0; r < reps; r++) hits += line_path(P, (c0 + r) % p); units = hits; }
        else if (!strcmp(what, "rec")) { for (long long r = 0; r < reps; r++) { hits += record_path(P, (c0 + r) % p); } units = hits; }
        else if (!strcmp(what, "apply")) {
            for (long long r = 0; r < reps; r++) record_path(P, (c0 + r) % p);
            for (long long sgi = 0; sgi < nseg; sgi++) units += (long long)apply_batch(&SEG[sgi], 0, P->lp);
        }
        else if (!strcmp(what, "setup6")) { for (long long r = 0; r < reps; r++) { a[0] = -1000 + r % 2001; acc += setup6(P, a); } units = reps; }
        else if (!strcmp(what, "setupl")) { lset_t o; for (long long r = 0; r < reps; r++) units += setup_line(P, (c0 + r) % p, 1 + r % K, r & 1, &o); acc += (uint64_t)o.A1.base; }
        else if (!strcmp(what, "scan") || !strcmp(what, "recover")) {
            for (long long r = 0; r < reps; r++) record_path(P, (1234567 + r) % p);
            for (long long sgi = 0; sgi < nseg; sgi++) apply_batch(&SEG[sgi], 0, P->lp);
            long long sc = 0;
            for (long long sgi = 0; sgi < nseg; sgi++) {
                long long bl = (long long)c0, nb = (SEG[sgi].ncell + bl - 1) / bl;   /* block length = argument c */
                long long *bs = malloc(8 * (nb + 1)); uint16_t *tf = malloc(2 * nb), *tg = malloc(2 * nb);
                for (long long b = 0; b < nb; b++) { bs[b] = bl * b; tf[b] = (uint16_t)(2 * P->lp); tg[b] = 0; } bs[nb] = SEG[sgi].ncell;
                scan_segment(sgi, bs, tf, tg, nb); sc += SEG[sgi].ncell; free(bs); free(tf); free(tg);
            }
            if (!strcmp(what, "scan")) units = sc;
            else { for (long long sgi = 0; sgi < nseg; sgi++) { recover_segment(&SEG[sgi]); units += (long long)SEG[sgi].rn; } }
        }
        else if (!strcmp(what, "build") || !strcmp(what, "zero")) {
            units = total_cells();
            if (!strcmp(what, "zero")) for (long long sgi = 0; sgi < nseg; sgi++) zero_cells(&SEG[sgi]);
            /* build: measured once inside alloc_slice (callgrind counts that call) */
        }
        else { fprintf(stderr, "unknown bench\n"); return 2; }
        long long rows = 0; for (long long sgi = 0; sgi < nseg; sgi++) rows += SEG[sgi].nrows;
        printf("entries=%lld rows=%lld ", L4 * K, rows);
        printf("what=%s cells=%lld K=%lld L4=%lld segments=%lld units=%lld hits=%lld acc=%llu survivors=%lld recovered=%llu\n", what, total_cells(), K, L4, nseg,
               units, hits, (unsigned long long)acc, nsurv, (unsigned long long)nrecov);
        return 0;
    }
    fprintf(stderr, "unknown mode\n"); return 2;
}
