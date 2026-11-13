#include "fp12.h"
#include "sha256.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const uint64_t R_ORDER[4] = {
  0x43e1f593f0000001ULL, 0x2833e84879b97091ULL,
  0xb85045b68181585dULL, 0x30644e72e131a029ULL };

/* Target-profile generator g = e(P0,Q0), coefficients in the profile basis. */
static const char *G_HEX[12] = {
 "28c6e04df059260df7d2d2a1f9b5f77676d1939847852c4ed50d2318744c1d5f",
 "17bb74adab1705c26133af1dac87044a3833ac011018e8158da48382bbd2dcd6",
 "0d3bd72f54d742f78ea9e6015c8ea2f2e7fbb728c9c905ec531dcf7de5b246f0",
 "090cb8ee97e091a667af03882b06c3ecb4e437993cbd1b05b98c7f9dfcfe9c40",
 "016b6d855b5cbf76f9829a309db52f5c442f65ae29f996af59d65f85f4afe78a",
 "0a0272204db51dadc0342bd318b9302a44faec12ff500bdd4d4b012ffe45f36f",
 "084f330485b09e866bc2f2ea2b897394deaf3f12aa31f28cb0552990967d4704",
 "27ed208e7a0b55ae6e710bbfbd2fd922669c026360e37cc5b2ab862411536104",
 "2067586885c3318eeffa1938c754fe3c60224ee5ae15e66af6b5104c47c8c5d8",
 "279db296f9d479292532c7c493d8e0722b6efae42158387564889c79fc038ee3",
 "2b03614464f04dd772d86df88674c270ffc8747ea13e72da95e3594468f222c4",
 "108c19d15f9446f744d0f110405d3856d6cc3bda6c4d537663729f5257628417" };

static void from_u64(fp *z, uint64_t x) {
  fp r; memset(&r,0,sizeof r); r.v[0]=x;
  memset(z,0,sizeof *z);
  /* multiply by R2 -> x*R mod p */
  fp_mul(z,&r,(const fp*)&(const uint64_t[4]){0xf32cfc5b538afa89ULL,0xb5e71911d44501fbULL,0x47ab1eff0a417ff6ULL,0x06d89f71cab8351fULL});
}

static void parse_hex_fp(fp *z, const char *h) {
  uint64_t w[4] = {0,0,0,0};
  for (int i = 0; i < 64; i++) {
    char c = h[i]; unsigned d = (c<='9')?(c-'0'):(c-'a'+10);
    int bitpos = 4*(63-i);
    w[bitpos>>6] |= (uint64_t)d << (bitpos & 63);
  }
  fp c; memcpy(c.v,w,sizeof w);
  extern const fp FP_ONE_MONT;
  fp_mul(z,&c,&(const fp){{0xf32cfc5b538afa89ULL,0xb5e71911d44501fbULL,0x47ab1eff0a417ff6ULL,0x06d89f71cab8351fULL}}); /* canonical -> Montgomery */
}

static void target_to_tower(fp12 *z, const char *hx[12]) {
  fp t[12]; for (int i=0;i<12;i++) parse_hex_fp(&t[i], hx[i]);
  fp nine; from_u64(&nine, 9);
  fp A,B,C,D,E,F,G,H,I,J,K,L, tmp;
  B=t[6]; fp_mul(&tmp,&nine,&B); fp_add(&A,&t[0],&tmp);
  H=t[7]; fp_mul(&tmp,&nine,&H); fp_add(&G,&t[1],&tmp);
  D=t[8]; fp_mul(&tmp,&nine,&D); fp_add(&C,&t[2],&tmp);
  J=t[9]; fp_mul(&tmp,&nine,&J); fp_add(&I,&t[3],&tmp);
  F=t[10]; fp_mul(&tmp,&nine,&F); fp_add(&E,&t[4],&tmp);
  L=t[11]; fp_mul(&tmp,&nine,&L); fp_add(&K,&t[5],&tmp);
  z->c[0].c[0].a=A; z->c[0].c[0].b=B;
  z->c[0].c[1].a=C; z->c[0].c[1].b=D;
  z->c[0].c[2].a=E; z->c[0].c[2].b=F;
  z->c[1].c[0].a=G; z->c[1].c[0].b=H;
  z->c[1].c[1].a=I; z->c[1].c[1].b=J;
  z->c[1].c[2].a=K; z->c[1].c[2].b=L;
}

/* 64-bit mixing hash over the 12 Montgomery coefficients (48 words). */
static inline uint32_t hash_index(const fp12 *y, uint64_t salt) {
  const uint64_t *w = (const uint64_t*)y;
  uint64_t h1 = 0x9e3779b97f4a7c15ULL ^ salt;
  uint64_t h2 = 0xbf58476d1ce4e5b9ULL + salt;
  for (int i = 0; i < 24; i++) {
    uint64_t a = w[2*i], b = w[2*i+1];
    h1 = (h1 ^ a) * 0x100000001b3ULL;
    h2 = (h2 ^ b) * 0x9e3779b97f4a7c15ULL;
    h2 ^= h2 >> 29;
  }
  h1 ^= h1 >> 33; h2 ^= h2 >> 27;
  return (uint32_t)((h1 ^ h2) & 31);
}

static inline void scalar_add_r(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]) {
  uint64_t s[5]; uint64_t carry=0;
  for (int i=0;i<4;i++){ __uint128_t u=(__uint128_t)x[i]+y[i]+carry; s[i]=(uint64_t)u; carry=(uint64_t)(u>>64); }
  s[4]=carry;
  uint64_t d[4]; uint64_t borrow=0;
  for (int i=0;i<4;i++){ __uint128_t u=(__uint128_t)s[i]-R_ORDER[i]-borrow; d[i]=(uint64_t)u; borrow=(uint64_t)(u>>64)&1; }
  if (s[4] || !borrow) memcpy(z,d,sizeof d); else memcpy(z,s,sizeof d);
}


/* ---- capped SHA-256 sampler and production table construction ---- */
#define K_MAX 64
static void stream_word(const uint8_t seed[32], uint64_t ctr, uint8_t out[32]) {
  uint8_t buf[40];
  memcpy(buf, seed, 32);
  for (int i = 0; i < 8; i++) buf[32 + i] = (uint8_t)(ctr >> (8 * i));
  sha256(buf, 40, out);
}
static int sample_scalar(const uint8_t seed[32], uint64_t *ctr, uint64_t x[4]) {
  for (int k = 0; k < K_MAX; k++) {
    uint8_t w[32];
    stream_word(seed, (*ctr)++, w);
    uint64_t v[4];
    for (int i = 0; i < 4; i++) { v[i] = 0; for (int b = 0; b < 8; b++) v[i] |= (uint64_t)w[8*i+b] << (8*b); }
    /* accept iff v < r */
    int lt = 0, eq = 1;
    for (int i = 3; i >= 0; i--) { if (v[i] != R_ORDER[i]) { lt = v[i] < R_ORDER[i]; eq = 0; break; } }
    if (eq || lt) { memcpy(x, v, sizeof v); return k + 1; }
  }
  return -1;
}
static void fp12_pow_u64(fp12 *z, const fp12 *x, const uint64_t e[4]) { fp12_pow(z, x, e); }

typedef struct { fp12 m; uint64_t c[4], d[4]; } rec;

static uint64_t step_checksum;
/* Worst-case audit: PIN_J < 32 pins the table index, 32 = production hash. */
static uint32_t PIN_J = 32;
static inline void rho_step(fp12 *y, uint64_t a[4], uint64_t b[4], const rec *tab, uint64_t salt) {
  uint32_t j = hash_index(y, salt);
  if (PIN_J < 32) j = PIN_J;
  fp12_mul(y, y, &tab[j].m);
  scalar_add_r(a, a, tab[j].c);
  scalar_add_r(b, b, tab[j].d);
  step_checksum ^= y->c[0].c[0].a.v[0];
}

static void build_table(rec tab[32], const fp12 *g) {
  fp12 cur = *g;
  for (int j=0;j<32;j++) {
    tab[j].m = cur;
    fp12_mul(&cur,&cur,g);
    for (int i=0;i<4;i++){ tab[j].c[i]= (uint64_t)(j*0x9e3779b97f4a7c15ULL) ^ (uint64_t)i; tab[j].d[i]= (uint64_t)(j*0xbf58476d1ce4e5b9ULL) + i; }
  }
}

int main(int argc, char **argv) {
  fp12 g; target_to_tower(&g, G_HEX);
  if (argc >= 2 && !strcmp(argv[1],"selftest")) {
    fp12 one, gr; fp12_set_one(&one);
    fp12_pow(&gr,&g,R_ORDER);
    printf("g_is_one=%d g^r_is_one=%d\n", fp12_eq(&g,&one), fp12_eq(&gr,&one));
    fp12 gg,x,y;
    fp12_mul(&gg,&g,&g);
    fp12_mul(&x,&gg,&g);   /* (g*g)*g */
    fp12_mul(&y,&g,&gg);   /* g*(g*g)  */
    printf("assoc=%d\n", fp12_eq(&x,&y));
    return 0;
  }
  if (argc >= 4 && !strcmp(argv[1],"mktable")) {
    /* mktable <seed-int> <out>: production table M_j = g^{c_j} z^{d_j} with
       scalars from the capped SHA-256 sampler of proof.md section 2. */
    uint8_t sigma[32], le[8];
    uint64_t sd = (uint64_t)strtoull(argv[2], 0, 10);
    for (int i = 0; i < 8; i++) le[i] = (uint8_t)(sd >> (8 * i));
    sha256(le, 8, sigma);
    uint64_t ctr = 0, att[66], c[32][4], d[32][4], x[4], a0[4], b0[4];
    for (int j = 0; j < 32; j++) {
      if (sample_scalar(sigma, &ctr, c[j]) < 0) { fprintf(stderr, "abort c%d\n", j); return 3; }
      if (sample_scalar(sigma, &ctr, d[j]) < 0) { fprintf(stderr, "abort d%d\n", j); return 3; }
    }
    if (sample_scalar(sigma, &ctr, a0) < 0) return 3;
    if (sample_scalar(sigma, &ctr, b0) < 0) return 3;
    if (sample_scalar(sigma, &ctr, x) < 0) return 3;
    fp12 z, t1, t2;
    fp12_pow(&z, &g, x);
    rec tab[32];
    for (int j = 0; j < 32; j++) {
      fp12_pow(&t1, &g, c[j]);
      fp12_pow(&t2, &z, d[j]);
      fp12_mul(&tab[j].m, &t1, &t2);
      memcpy(tab[j].c, c[j], sizeof c[j]);
      memcpy(tab[j].d, d[j], sizeof d[j]);
    }
    FILE *f = fopen(argv[3], "wb");
    if (!f) { perror("fopen"); return 2; }
    fwrite(x, 8, 4, f); fwrite(a0, 8, 4, f); fwrite(b0, 8, 4, f);
    fwrite(tab, sizeof(rec), 32, f);
    fwrite(&ctr, 8, 1, f);
    fclose(f);
    printf("mktable seed=%llu words=%llu\n", (unsigned long long)sd, (unsigned long long)ctr);
    return 0;
  }
  if (argc >= 5 && !strcmp(argv[1],"stepwc")) {
    /* stepwc K j file: K steps of the shipped step on the production table.
       j in 0..31 pins that table entry; j = 32 uses the production mixer. */
    long k = atol(argv[2]);
    int pj = atoi(argv[3]);
    FILE *f = fopen(argv[4], "rb");
    if (!f) { perror("fopen"); return 2; }
    uint64_t x[4], a0[4], b0[4];
    rec tab[32];
    if (fread(x, 8, 4, f) != 4 || fread(a0, 8, 4, f) != 4 || fread(b0, 8, 4, f) != 4 ||
        fread(tab, sizeof(rec), 32, f) != 32) { fprintf(stderr, "short table\n"); return 2; }
    fclose(f);
    PIN_J = (pj < 32) ? (uint32_t)pj : 32u;
    fp12 z, y, t1, t2;
    fp12_pow(&z, &g, x);
    fp12_pow(&t1, &g, a0);
    fp12_pow(&t2, &z, b0);
    fp12_mul(&y, &t1, &t2);
    uint64_t a[4], b[4];
    memcpy(a, a0, sizeof a); memcpy(b, b0, sizeof b);
    for (long i = 0; i < k; i++) rho_step(&y, a, b, tab, 0x123456789abcdefULL);
    printf("stepwc k=%ld j=%d checksum=%016llx\n", k, pj, (unsigned long long)step_checksum);
    return 0;
  }
  if (argc >= 3 && !strcmp(argv[1],"hashdump")) {
    /* Print the production mixer index for n fixed encoded elements. */
    int n = atoi(argv[2]);
    for (int e = 0; e < n; e++) {
      fp12 y; memset(&y, 0, sizeof y);
      fp *slot = (fp *)&y;
      for (int i = 0; i < 12; i++) from_u64(&slot[i], (uint64_t)(e * 97 + i + 1));
      uint64_t wchk = 0;
      for (int i = 0; i < 48; i++) wchk = wchk * 0x100000001b3ULL + ((uint64_t *)&y)[i];
      printf("e=%d j=%u wchk=%016llx\n", e, hash_index(&y, 0x123456789abcdefULL),
             (unsigned long long)wchk);
    }
    return 0;
  }
  if (argc >= 5 && !strcmp(argv[1],"worstcase")) {
    /* worstcase K j pattern: pinned table index j, operand pattern. */
    long k = atol(argv[2]);
    int pj = atoi(argv[3]);
    int pat = atoi(argv[4]);
    if (pj < 0 || pj > 31) { fprintf(stderr, "bad index\n"); return 2; }
    PIN_J = (uint32_t)pj;
    rec tab[32]; build_table(tab,&g);
    fp12 y = g;
    if (pat == 1) { for (int i=0;i<48;i++) ((uint64_t*)&y)[i] = 0xffffffffffffffffULL; }
    else if (pat == 2) { for (int i=0;i<48;i++) ((uint64_t*)&y)[i] = 1ULL; }
    else if (pat == 3) { for (int i=0;i<48;i++) ((uint64_t*)&y)[i] = FP_ONE_MONT.v[i&3]; }
    else if (pat == 20) { /* limbs = p-1: forces carry/borrow and full compare chains */
      const uint64_t pl[4] = {0x3c208c16d87cfd47ULL,0x97816a916871ca8dULL,
                              0xb85045b68181585dULL,0x30644e72e131a029ULL};
      for (int i=0;i<48;i++) ((uint64_t*)&y)[i] = pl[i&3] - 1;
    }
    else if (pat == 21) { /* limbs = p: every addition overflows the modulus */
      const uint64_t pl[4] = {0x3c208c16d87cfd47ULL,0x97816a916871ca8dULL,
                              0xb85045b68181585dULL,0x30644e72e131a029ULL};
      for (int i=0;i<48;i++) ((uint64_t*)&y)[i] = pl[i&3];
    }
    else if (pat >= 4) { /* deterministic pseudo-random operand pattern */
      uint64_t s = (uint64_t)pat * 0x9e3779b97f4a7c15ULL + 0x1234567ULL;
      for (int i=0;i<48;i++) {
        s += 0x9e3779b97f4a7c15ULL;
        uint64_t z = s;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z ^= z >> 31;
        ((uint64_t*)&y)[i] = z % 0x30644e72e131a029ULL;
      }
    }
    uint64_t a[4]={1,0,0,0}, b[4]={2,0,0,0};
    for (long i=0;i<k;i++) rho_step(&y,a,b,tab,0x123456789abcdefULL);
    printf("worstcase k=%ld j=%d pat=%d checksum=%016llx\n", k, pj, pat,
           (unsigned long long)step_checksum);
    return 0;
  }
  long k = (argc>=2)? atol(argv[1]) : 10000;
  rec tab[32]; build_table(tab,&g);
  fp12 y=g; uint64_t a[4]={1,0,0,0}, b[4]={2,0,0,0};
  for (long i=0;i<k;i++) rho_step(&y,a,b,tab,0x123456789abcdefULL);
  printf("checksum=%016llx y0=%016llx\n",(unsigned long long)step_checksum,(unsigned long long)y.c[0].c[0].a.v[0]);
  return 0;
}
