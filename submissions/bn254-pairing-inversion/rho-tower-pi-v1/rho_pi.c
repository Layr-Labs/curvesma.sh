/* bn254-pairing-inversion: GT rho-tower kernel plus G1 scalar multiplication.
 *
 * Two independent parts, both shipped as evidence:
 *   1. the measured GT walk step (same scalar tower arithmetic as the
 *      bn254-gt-dlp rho-tower-v1 package);
 *   2. G1 arithmetic (Jacobian doubling, mixed addition by the affine fixed
 *      base P0=(1,2), double-and-add scalar multiplication, affine output).
 *
 * The attack described in proof.md runs part 1 to recover x with g^x = z,
 * then part 2 to emit P = [x]P0.  Q is the profile generator Q0, so
 * e(P,Q0) = e(P0,Q0)^x = g^x = z by bilinearity of the reduced pairing.
 */
#include "fp12.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const uint64_t R_ORDER[4] = {
  0x43e1f593f0000001ULL, 0x2833e84879b97091ULL,
  0xb85045b68181585dULL, 0x30644e72e131a029ULL };

/* p - 2, little-endian 64-bit limbs, for Fermat inversion in Fp. */
static const uint64_t P_MINUS_2[4] = {
  0x3c208c16d87cfd45ULL, 0x97816a916871ca8dULL,
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

/* ---------------- part 1: GT walk step (measured) ---------------- */

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

/* ---------------- part 2: G1 arithmetic ---------------- */

typedef struct { fp X, Y, Z; } g1;   /* Jacobian; Z == 0 is infinity */

static int fp_is_zero(const fp *x) {
  return (x->v[0]|x->v[1]|x->v[2]|x->v[3]) == 0;
}

static void fp_inv(fp *z, const fp *x) {
  fp r = FP_ONE_MONT;
  for (int i = 254; i >= 0; i--) {
    fp_mul(&r,&r,&r);
    if ((P_MINUS_2[i>>6] >> (i&63)) & 1) fp_mul(&r,&r,x);
  }
  *z = r;
}

static void g1_set_inf(g1 *p) { memset(p,0,sizeof *p); }

static int g1_is_inf(const g1 *p) { return fp_is_zero(&p->Z); }

static int g1_on_curve_affine(const fp *x, const fp *y) {
  fp x2,x3,rhs,three,t;
  fp_mul(&x2,x,x);
  fp_mul(&x3,&x2,x);
  from_u64(&three,3);
  fp_add(&rhs,&x3,&three);
  fp_mul(&t,y,y);
  return fp_eq(&t,&rhs);
}

/* d = 2*p, a = 0 */
static void g1_double(g1 *d, const g1 *p) {
  if (g1_is_inf(p) || fp_is_zero(&p->Y)) { g1_set_inf(d); return; }
  fp A,B,C,D,E,F,t,Z3;
  fp_mul(&Z3,&p->Y,&p->Z); fp_add(&Z3,&Z3,&Z3);   /* before X/Y are overwritten */
  fp_mul(&A,&p->X,&p->X);
  fp_mul(&B,&p->Y,&p->Y);
  fp_mul(&C,&B,&B);
  fp_add(&t,&p->X,&B); fp_mul(&t,&t,&t); fp_sub(&t,&t,&A); fp_sub(&t,&t,&C);
  fp_add(&D,&t,&t);
  fp_add(&E,&A,&A); fp_add(&E,&E,&A);            /* 3A */
  fp_mul(&F,&E,&E);
  fp_add(&t,&D,&D);
  fp_sub(&d->X,&F,&t);
  fp_sub(&t,&D,&d->X); fp_mul(&t,&E,&t);
  fp_add(&C,&C,&C); fp_add(&C,&C,&C); fp_add(&C,&C,&C);  /* 8C */
  fp_sub(&d->Y,&t,&C);
  d->Z = Z3;
}

/* d = p + q with q affine (Z=1): mixed addition */
static void g1_add_mixed(g1 *d, const g1 *p, const fp *qx, const fp *qy) {
  if (g1_is_inf(p)) {
    d->X = *qx; d->Y = *qy;
    memset(&d->Z,0,sizeof d->Z); d->Z.v[0] = 0xd35d438dc58f0d9dULL;
    d->Z.v[1] = 0x0a78eb28f5c70b3dULL; d->Z.v[2] = 0x666ea36f7879462cULL; d->Z.v[3] = 0x0e0a77c19a07df2fULL;
    return;
  }
  fp Z1Z1,U2,S2,H,HH,I,J,r,V,t;
  fp_mul(&Z1Z1,&p->Z,&p->Z);
  fp_mul(&U2,qx,&Z1Z1);
  fp_mul(&S2,&p->Z,&Z1Z1); fp_mul(&S2,qy,&S2);
  fp_sub(&H,&U2,&p->X);
  if (fp_is_zero(&H)) {
    if (fp_eq(&S2,&p->Y)) { g1_double(d,p); return; }
    g1_set_inf(d); return;
  }
  fp_mul(&HH,&H,&H);
  fp_add(&I,&HH,&HH); fp_add(&I,&I,&I);           /* 4HH */
  fp_mul(&J,&H,&I);
  fp_sub(&r,&S2,&p->Y); fp_add(&r,&r,&r);         /* 2(S2-Y1) */
  fp_mul(&V,&p->X,&I);
  fp_mul(&t,&r,&r);
  fp_sub(&d->X,&t,&J); fp_sub(&d->X,&d->X,&V); fp_sub(&d->X,&d->X,&V);
  fp_sub(&t,&V,&d->X); fp_mul(&t,&r,&t);        /* r*(V-X3) */
  fp_mul(&J,&p->Y,&J); fp_add(&J,&J,&J);        /* 2*Y1*J */
  fp_sub(&d->Y,&t,&J);
  fp_add(&t,&p->Z,&H); fp_mul(&t,&t,&t); fp_sub(&t,&t,&Z1Z1); fp_sub(&t,&t,&HH);
  d->Z = t;
}

/* P = [k]P0 with fixed affine base (1,2), double-and-add from the top bit. */
static void g1_mul_base(g1 *out, const uint64_t k[4]) {
  g1 r; g1_set_inf(&r);
  fp one, two;
  from_u64(&one,1); from_u64(&two,2);
  for (int i = 255; i >= 0; i--) {
    g1_double(&r,&r);
    if ((k[i>>6] >> (i&63)) & 1) g1_add_mixed(&r,&r,&one,&two);
  }
  *out = r;
}

static void g1_to_affine(fp *x, fp *y, const g1 *p) {
  if (g1_is_inf(p)) { memset(x,0,sizeof *x); memset(y,0,sizeof *y); return; }
  fp zi, zi2, zi3;
  fp_inv(&zi,&p->Z);
  fp_mul(&zi2,&zi,&zi);
  fp_mul(&zi3,&zi2,&zi);
  fp_mul(x,&p->X,&zi2);
  fp_mul(y,&p->Y,&zi3);
}

static void print_fp_hex(const char *label, const fp *x) {
  /* convert out of Montgomery form, then print 64 lowercase hex digits */
  fp one, c;
  memset(&one,0,sizeof one); one.v[0]=1;   /* plain 1: fp_mul(x,1) divides by R */
  fp_mul(&c,x,&one);
  printf("%s=", label);
  for (int i = 3; i >= 0; i--) printf("%016llx", (unsigned long long)c.v[i]);
  printf("\n");
}

static void parse_scalar(uint64_t k[4], const char *h) {
  memset(k,0,sizeof(uint64_t)*4);
  size_t n = strlen(h);
  for (size_t i = 0; i < n && i < 64; i++) {
    char c = h[n-1-i]; unsigned d = (c<='9')?(c-'0'):(c-'a'+10);
    k[i>>4] |= (uint64_t)d << (4*(i&15));
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
    fp12_mul(&x,&gg,&g);
    fp12_mul(&y,&g,&gg);
    printf("assoc=%d\n", fp12_eq(&x,&y));
    fp onef, twof; from_u64(&onef,1); from_u64(&twof,2);
    printf("P0_on_curve=%d\n", g1_on_curve_affine(&onef,&twof));
    g1 rp; g1_mul_base(&rp, R_ORDER);
    printf("r_times_P0_is_inf=%d\n", g1_is_inf(&rp));
    g1 two_p; g1_mul_base(&two_p, (const uint64_t[4]){2,0,0,0});
    fp ax, ay; g1_to_affine(&ax,&ay,&two_p);
    print_fp_hex("2P0.x", &ax); print_fp_hex("2P0.y", &ay);
    return 0;
  }
  if (argc >= 3 && !strcmp(argv[1],"point")) {
    uint64_t k[4]; parse_scalar(k, argv[2]);
    g1 p; g1_mul_base(&p, k);
    if (g1_is_inf(&p)) { printf("P=null\n"); return 0; }
    fp ax, ay; g1_to_affine(&ax,&ay,&p);
    print_fp_hex("P.x", &ax); print_fp_hex("P.y", &ay);
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
