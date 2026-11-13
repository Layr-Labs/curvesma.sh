/* BN254 GT multiplication kernel (CIOS), scalar 64-bit limbs only.
 * Counted like the challenge reference costs: gcc -O3 -fno-tree-vectorize
 * -fno-builtin, valgrind callgrind event Ir, loops of 4096 and 8192.
 * "t" mode prints the 96 words of a fixed product for external checking.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
typedef struct { uint64_t v[4]; } fp;

static const uint64_t P[4] = {0x3c208c16d87cfd47ULL, 0x97816a916871ca8dULL,
                              0xb85045b68181585dULL, 0x30644e72e131a029ULL};
static uint64_t N0;

/* Derived from tracks/bn254-gt-dlp/target-profile.json:
   GCOEF[i] = coefficient i of generator_g, in Montgomery form (x*2^256 mod p);
   REXP = subgroup order r; RMONT = 2^256 mod p; NINEMONT = 9*2^256 mod p. */
static const uint64_t GCOEF[12][4] = {
  {0x53e1d9fc3a8329efULL,0x9254a1949ff5465fULL,0x3d01af561fad5084ULL,0x2ae045056b3b7c1eULL},
  {0xd514f44343aae165ULL,0x6d330380ed14f311ULL,0xe6490e6130d72457ULL,0x0e220885bfce1f4dULL},
  {0x7f138852111c85f6ULL,0x6e42d4f25e64fff5ULL,0x1e2c29f438c5faf0ULL,0x05a226123d0a333dULL},
  {0x0221eab5ce4d20adULL,0x69b3e161979120e2ULL,0xea0115f3d02e5cd4ULL,0x18d86b9b7ecf1950ULL},
  {0x5869976aa182adffULL,0xb87c7bb0c0a87668ULL,0xf7d1d1a0ddced101ULL,0x07918bfca0015959ULL},
  {0x0be63e0e6309bae7ULL,0x2efbc923df47f05fULL,0xf5f5b66cc5d7e7f8ULL,0x253dac9cb214fa73ULL},
  {0x2e02a64acbd60549ULL,0xd618018ea58e4addULL,0x14d585f1a45ba647ULL,0x1832226987c434fcULL},
  {0x9458abcb56d24998ULL,0xb17540bd2a9e5adbULL,0x9a9983c82e401a9fULL,0x1614817a84c16291ULL},
  {0x172d1f257a4d598eULL,0xddf5bc7b7ffb5ac0ULL,0xae0b22c0bbb0f602ULL,0x1b158f3c2fae9b18ULL},
  {0x7002907c28ebfe11ULL,0x7b0591d3d080da67ULL,0xde7e5aa2181f138eULL,0x210e437dfc43d951ULL},
  {0xf16c96d081754cdbULL,0xce0394312bceeb55ULL,0x644e4dcf1f01ff0aULL,0x0cbea85ee0b236ccULL},
  {0xd34bab373157aa84ULL,0x3511ed44fd0d8598ULL,0x67e42a0bc2ced972ULL,0x2b8f1d5dfd20c55bULL},
};
static const uint64_t REXP[4] = {0x43e1f593f0000001ULL,0x2833e84879b97091ULL,0xb85045b68181585dULL,0x30644e72e131a029ULL };
static const uint64_t RMONT[4] = {0xd35d438dc58f0d9dULL,0x0a78eb28f5c70b3dULL,0x666ea36f7879462cULL,0x0e0a77c19a07df2fULL };
static const uint64_t NINEMONT[4] = {0xf60647ce410d7ff7ULL,0x2f3d6f4dd31bd011ULL,0x2943337e3940c6d1ULL,0x1d9598e8a7e39857ULL };


/* Public increment table: c_j = SHA256("bn254gt-c"||j) mod r, d_j likewise. */
static const uint64_t INC_C[32][4] = {
  {0x5dee7fd75d73f0e0ULL,0xe733562e34a9f9c4ULL,0x1f3466ca35d44d2dULL,0x1f9fda5d5acfcbd3ULL},
  {0xf967d6ea583bc8a5ULL,0x2207e5653d72609cULL,0x5afa9fbed2be15e1ULL,0x074ae5d48b757b1dULL},
  {0x794de3ecdecb90dcULL,0xeea89ef80384faa0ULL,0xd309fe692c8fd03fULL,0x1bfafe7c69c58864ULL},
  {0x5419d149722f6683ULL,0xce28cf4d76592fddULL,0x7e63803ae896cc79ULL,0x248878499db3ff80ULL},
  {0x923ec3584bac33f3ULL,0xb2517e61369f83e4ULL,0x7f5ca5cae4cd3ad4ULL,0x2c85255159ef4d4dULL},
  {0x651cc30c04ff53c8ULL,0x9788d66fb183a508ULL,0xce0c9ad728b200eaULL,0x142bdbefd127fdceULL},
  {0x778855c0eed6820fULL,0xd84e8f738350c8abULL,0xab21d24ddfd69533ULL,0x1d322cf62a4111adULL},
  {0x218e7bdcb125c10cULL,0x652fdca7fb0ffda4ULL,0xc52eb8d62155c301ULL,0x105f4c0c363e8224ULL},
  {0x8dff525a0abf70e6ULL,0x7e618c5f7699bb7dULL,0xced0944388692bccULL,0x00c7b85d67517d4aULL},
  {0x12c57ef0b1d5fa74ULL,0x94cd61b169990ac9ULL,0x73fa336dc16e8124ULL,0x001c7154a3aa82b7ULL},
  {0x2a9e04020642e410ULL,0x0742c9e77d0fc3c4ULL,0x3da29ba3a4ceb88dULL,0x1a6ba85b0dbedd4eULL},
  {0x0762ae536fc8c313ULL,0x0efb5a67e944a83aULL,0x5cccf0b119abedcaULL,0x0f665e42f9882c60ULL},
  {0x2ecdcdc32af97aaeULL,0xdac64f2d946035f8ULL,0x3a8f2f233044310aULL,0x0878f8b321c75729ULL},
  {0x92a6019fb49728a8ULL,0x970436f123f6d99fULL,0x563f60f9bfc2d5d2ULL,0x1cd4244675bc9aabULL},
  {0x133541d40f3808d9ULL,0xd87614fe6e7abc82ULL,0xcd14ede6b5580f18ULL,0x2b1e98bb9edf6a32ULL},
  {0xd74fa2db62473ab2ULL,0xa9b4c6598fb6f222ULL,0xce7a9fdb63a35785ULL,0x095ee9caa3ccd7a5ULL},
  {0x7eb1e9cd6a3aba14ULL,0x69f5aaa3c6b65fd6ULL,0x4f402613be8f3569ULL,0x27c645ba1a07f867ULL},
  {0x070b54ad5a2fa098ULL,0x7313027291943ddfULL,0xd00c1d144fb02075ULL,0x0e7540eba6a1909fULL},
  {0xa8fbf7023c19f104ULL,0x0efe92fab787a0c2ULL,0x73668b8392dc6b3fULL,0x2a984a0abb4b564fULL},
  {0xb8539094ec23d12cULL,0x539357b3b2466692ULL,0x146a0af2fc9695b1ULL,0x2e29cd9e15a0454aULL},
  {0x1b67215768e1fd62ULL,0xb7c4db13201e5041ULL,0x184b7f8d01548a61ULL,0x1d037d20185b6187ULL},
  {0x452ea17ca2319d12ULL,0xd7d8619a61b72afcULL,0x4ef0091550018d95ULL,0x07163d9440f6e475ULL},
  {0x593dbe0ee5f9060bULL,0x8a86db697c1cbc3aULL,0x4751e77ee08e9a29ULL,0x1c212168c100577fULL},
  {0x547b4b16c9c7ad99ULL,0x084bb0d034817ab8ULL,0xa30333548ff21440ULL,0x150c735bb90b7c64ULL},
  {0xd1bb658351dbe723ULL,0x37d6a5411557bf1bULL,0xcb6fc5f69022cc14ULL,0x30010232e153bfdeULL},
  {0x78fb5b41255c4d08ULL,0xf768d021f31db2e8ULL,0xb91e996f7b6ad74dULL,0x0894a31b118d9c1dULL},
  {0xa45c8c99f1014b91ULL,0xf751aff4bcab62e4ULL,0x9d67a9b60e1a30eaULL,0x129ea2ec3483f000ULL},
  {0x5615098ed98dadc7ULL,0x0ba70799c008399eULL,0xefea21c26ce4f10dULL,0x295dbd46053acbdcULL},
  {0x84fb5d093c98b815ULL,0xf9bf8f2525c6db06ULL,0x689b6c1610d1c99eULL,0x0533e1404e169c3cULL},
  {0x8f9fe9edcb8c0b69ULL,0xf17feaffab029232ULL,0x562ba88e4221eb6eULL,0x13254a0c3733ea6dULL},
  {0x08ea844ce63bf1bcULL,0x0fdb816ea3a17dc1ULL,0x9504807eec162241ULL,0x1061c5782f94239bULL},
  {0xed5fb8362ff878d2ULL,0xb58c397940de9ac7ULL,0xc46fa2a01d5c114fULL,0x2e6ffa09d8b0a58aULL},
};
static const uint64_t INC_D[32][4] = {
  {0x0e40a9eb347a30beULL,0xecdc7eec2bff97a4ULL,0xddebda55e150e2d1ULL,0x00bf385ce4f5c0fcULL},
  {0x7706f342687ffcb9ULL,0x2db2fab842c6ecd1ULL,0x2d949a286a5e2d04ULL,0x0fc43f8da0d4455eULL},
  {0xa3694644392876c1ULL,0xe7408e66d355733eULL,0x361ca1bbdbd3e84aULL,0x11cba121ff7c940eULL},
  {0x366629d867c9b28cULL,0x1d4411990d9d6863ULL,0x8ea4a3cc3ea436baULL,0x1d3e9ce13df29dcaULL},
  {0xf5fb258c0236a77dULL,0xc301c86b6b874160ULL,0x4a4ff2c3e029c2bbULL,0x1e71e605ba4dfa0fULL},
  {0xcd855769ee8b834dULL,0x332e39fa656c65c3ULL,0xb8256d0ab590e20fULL,0x115a9553638eb15eULL},
  {0xf906cda49ce8a622ULL,0x768a9baf95d23b64ULL,0x229641e12da93d64ULL,0x05323cb36a0e7087ULL},
  {0x7bce9dcf29ee5831ULL,0x4bfaede46c8007c9ULL,0x24a081053fc5d73dULL,0x2ae5667a7108b678ULL},
  {0xbad50ed2a2e152ddULL,0x19923192b7c4bcebULL,0xfe37d1639c7ca595ULL,0x163c2365708b059cULL},
  {0xc30e7abdce56a8c9ULL,0x2bc51979bd3ac007ULL,0x78ce7b4132c7d3b1ULL,0x21a2b14cb0000a97ULL},
  {0x8f3757e0535f4ae9ULL,0xdc66a2e303f17580ULL,0x6395fc7a4ebf9685ULL,0x258e7c861e83520bULL},
  {0x5541cdb2bb2f7184ULL,0x7e08d19882e34420ULL,0x3bdd107d91803e39ULL,0x1323b29cdbb042c4ULL},
  {0xa9202572e14a9f80ULL,0xd1b0dd7a57487e62ULL,0xa788fa6772b163bbULL,0x0daba608e4dffc76ULL},
  {0xae035b010f43a24cULL,0x3cbfc727fb832041ULL,0xa3564aee6738a40eULL,0x08c5c888d81780fdULL},
  {0x56c7cba703d484adULL,0xa1e9d623b5b6e89eULL,0xee16e8938645bf73ULL,0x228d4998b51b2d25ULL},
  {0x7d7bf099c11bf094ULL,0x97f77f917aaf8819ULL,0xd0faf285ba904214ULL,0x244caf9df81b40e2ULL},
  {0x00745b7e92ed6dfbULL,0x03eae4b5e0e5c161ULL,0xf442b4629b301a41ULL,0x2056e9e1dabb0057ULL},
  {0xed4dc0103e2a8dbdULL,0x0e81de3dc70f36c9ULL,0x00d1f5263897fa1dULL,0x224289e74d1d9097ULL},
  {0x20a03e99dac1c9a2ULL,0xf973d8ea8b4b2dfbULL,0xf6a63df0459808c5ULL,0x2744f1334bdaeea5ULL},
  {0xd488144c6ccaecb5ULL,0x420bc20e1f64ca46ULL,0x4f172e4ad3c9a55fULL,0x28072e0d80d7a04aULL},
  {0x6717a1d84652e41aULL,0x4e11e5285375493fULL,0x095f7d96bac5548eULL,0x01f37361a85ed7d3ULL},
  {0x10fd2223d6539e9eULL,0x0ef3bc4dcde8a176ULL,0xfb62e371be19c80fULL,0x030288d574de786dULL},
  {0xe63fedc6e58f51d5ULL,0xe2bf4d8083a289fdULL,0x8790c333bf2a113dULL,0x20284f25f8df8293ULL},
  {0x91889b8d1c7f8430ULL,0x137ae9ba7a527edeULL,0x83d92973653ed12eULL,0x2035c96232143c26ULL},
  {0x3f610f21c7317ec0ULL,0xa0b51136642c7914ULL,0x4f0a73fe01c0c2caULL,0x083f1e3da05a4cf8ULL},
  {0x4b002086eae5a763ULL,0x22ddae1482503a9dULL,0xd29abb549d469f98ULL,0x0c80b0bf2b281ebfULL},
  {0x64610e281c313261ULL,0xff12cc4342634784ULL,0x7d0943cb3af52dfdULL,0x12d55905d5a85093ULL},
  {0x7609f5129fb8637cULL,0xa37dfccdccd5e821ULL,0xe1f103ca0112c62eULL,0x102473e173d7152aULL},
  {0x86f7c20c31d333a9ULL,0x8b72bb5932c916bbULL,0xc76065116b8e0908ULL,0x26087539a6f521c5ULL},
  {0x49b3f4119aae8f35ULL,0x3a9b6b5ee49c52a6ULL,0xa493d04bcdf0a48eULL,0x1407de323f607a9aULL},
  {0x0887a8a4b267f364ULL,0xa7347cafb201d257ULL,0x5597f2c5fe32a41eULL,0x105e9dcb8b8e9b96ULL},
  {0x9de7732a5dbe271dULL,0xa7155662ba434f0bULL,0x27d77c9aa32aaf07ULL,0x18323e7649d958e1ULL},
};


static inline void fadd(fp *r, const fp *a, const fp *b) {
    uint64_t t0, t1, t2, t3, c;
    u128 w;
    w = (u128)a->v[0] + b->v[0]; t0 = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)a->v[1] + b->v[1] + c; t1 = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)a->v[2] + b->v[2] + c; t2 = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)a->v[3] + b->v[3] + c; t3 = (uint64_t)w;
    uint64_t d = 0, s0, s1, s2, s3;
    w = (u128)t0 - P[0]; s0 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t1 - P[1] - d; s1 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t2 - P[2] - d; s2 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t3 - P[3] - d; s3 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    uint64_t m = (uint64_t)0 - (uint64_t)(d == 0);
    r->v[0] = (s0 & m) | (t0 & ~m); r->v[1] = (s1 & m) | (t1 & ~m);
    r->v[2] = (s2 & m) | (t2 & ~m); r->v[3] = (s3 & m) | (t3 & ~m);
}

static inline void fsub(fp *r, const fp *a, const fp *b) {
    uint64_t t0, t1, t2, t3, d;
    u128 w;
    w = (u128)a->v[0] - b->v[0]; t0 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)a->v[1] - b->v[1] - d; t1 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)a->v[2] - b->v[2] - d; t2 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)a->v[3] - b->v[3] - d; t3 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    uint64_t m = (uint64_t)0 - d, c;
    w = (u128)t0 + (P[0] & m); r->v[0] = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)t1 + (P[1] & m) + c; r->v[1] = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)t2 + (P[2] & m) + c; r->v[2] = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)t3 + (P[3] & m) + c; r->v[3] = (uint64_t)w;
}

/* Montgomery CIOS: r = a*b*R^{-1} mod P, a,b < P, fully register-resident */
static inline void fmul(fp *r, const fp *a, const fp *b) {
    const uint64_t *A = a->v;
    uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t bi = b->v[i], c, m;
        u128 w;
        w = (u128)A[0] * bi + t0;        t0 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)A[1] * bi + t1 + c;    t1 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)A[2] * bi + t2 + c;    t2 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)A[3] * bi + t3 + c;    t3 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)t4 + c;                t4 = (uint64_t)w;
        m = t0 * N0;
        w = (u128)m * P[0] + t0;         c = (uint64_t)(w >> 64);
        w = (u128)m * P[1] + t1 + c;     t0 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)m * P[2] + t2 + c;     t1 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)m * P[3] + t3 + c;     t2 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)t4 + c;                t3 = (uint64_t)w; t4 = (uint64_t)(w >> 64);
    }
    uint64_t d, s0, s1, s2, s3;
    u128 w;
    w = (u128)t0 - P[0];         s0 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t1 - P[1] - d;     s1 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t2 - P[2] - d;     s2 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t3 - P[3] - d;     s3 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    uint64_t m = (uint64_t)0 - (uint64_t)(d == 0);
    r->v[0] = (s0 & m) | (t0 & ~m); r->v[1] = (s1 & m) | (t1 & ~m);
    r->v[2] = (s2 & m) | (t2 & ~m); r->v[3] = (s3 & m) | (t3 & ~m);
}

typedef struct { fp c0, c1; } fp2;
typedef struct { fp2 c[3]; } fp6;
typedef struct { fp6 c0, c1; } fp12;

static inline void f2add(fp2 *r, const fp2 *a, const fp2 *b) { fadd(&r->c0, &a->c0, &b->c0); fadd(&r->c1, &a->c1, &b->c1); }
static inline void f2sub(fp2 *r, const fp2 *a, const fp2 *b) { fsub(&r->c0, &a->c0, &b->c0); fsub(&r->c1, &a->c1, &b->c1); }

static inline void f2mul(fp2 *r, const fp2 *a, const fp2 *b) {
    fp t0, t1, s0, s1;
    fmul(&t0, &a->c0, &b->c0);
    fmul(&t1, &a->c1, &b->c1);
    fadd(&s0, &a->c0, &a->c1);
    fadd(&s1, &b->c0, &b->c1);
    fmul(&s0, &s0, &s1);
    fsub(&s0, &s0, &t0);
    fsub(&s0, &s0, &t1);
    fsub(&r->c0, &t0, &t1);
    r->c1 = s0;
}

static inline void mul_by_xi(fp2 *r, const fp2 *a) { /* (a0+a1u)(9+u) = (9a0-a1) + (9a1+a0)u */
    fp2 d, t, au;
    f2add(&d, a, a); f2add(&d, &d, &d); f2add(&d, &d, &d);   /* 8a */
    f2add(&t, &d, a);                                        /* 9a */
    au.c0.v[0] = au.c0.v[1] = au.c0.v[2] = au.c0.v[3] = 0;
    au.c0.v[0] = 0;
    fsub(&au.c0, &au.c0, &a->c1);
    au.c1 = a->c0;
    f2add(r, &t, &au);
}

static inline void f6add(fp6 *r, const fp6 *a, const fp6 *b) { for (int i = 0; i < 3; i++) f2add(&r->c[i], &a->c[i], &b->c[i]); }
static inline void f6sub(fp6 *r, const fp6 *a, const fp6 *b) { for (int i = 0; i < 3; i++) f2sub(&r->c[i], &a->c[i], &b->c[i]); }

static inline void f6mul(fp6 *r, const fp6 *a, const fp6 *b) {
    fp2 v0, v1, v2, t0, t1, t2, s0, s1, x;
    f2mul(&v0, &a->c[0], &b->c[0]);
    f2mul(&v1, &a->c[1], &b->c[1]);
    f2mul(&v2, &a->c[2], &b->c[2]);
    f2add(&s0, &a->c[1], &a->c[2]); f2add(&s1, &b->c[1], &b->c[2]);
    f2mul(&t0, &s0, &s1); f2sub(&t0, &t0, &v1); f2sub(&t0, &t0, &v2);
    f2add(&s0, &a->c[0], &a->c[1]); f2add(&s1, &b->c[0], &b->c[1]);
    f2mul(&t1, &s0, &s1); f2sub(&t1, &t1, &v0); f2sub(&t1, &t1, &v1);
    f2add(&s0, &a->c[0], &a->c[2]); f2add(&s1, &b->c[0], &b->c[2]);
    f2mul(&t2, &s0, &s1); f2sub(&t2, &t2, &v0); f2sub(&t2, &t2, &v2);
    mul_by_xi(&x, &t0); f2add(&r->c[0], &v0, &x);
    mul_by_xi(&x, &v2); f2add(&r->c[1], &t1, &x);
    f2add(&r->c[2], &t2, &v1);
}

static inline void mul_by_v(fp6 *r, const fp6 *a) {
    fp2 x;
    mul_by_xi(&x, &a->c[2]);
    r->c[0] = x; r->c[1] = a->c[0]; r->c[2] = a->c[1];
}

static inline void f12mul(fp12 *r, const fp12 *a, const fp12 *b) {
    fp6 ac, bd, s0, s1, t, vbd;
    f6mul(&ac, &a->c0, &b->c0);
    f6mul(&bd, &a->c1, &b->c1);
    f6add(&s0, &a->c0, &a->c1);
    f6add(&s1, &b->c0, &b->c1);
    f6mul(&t, &s0, &s1);
    f6sub(&t, &t, &ac);
    f6sub(&t, &t, &bd);
    mul_by_v(&vbd, &bd);
    f6add(&r->c0, &ac, &vbd);
    r->c1 = t;
}

static uint64_t seed_state = 0x243f6a8885a308d3ULL;
static inline uint64_t xr(void) { seed_state ^= seed_state << 13; seed_state ^= seed_state >> 7; seed_state ^= seed_state << 17; return seed_state; }

/* w-basis (profile) -> tower: for m=0..5, j=m&1, k=m>>1,
   a = (c_m + 9 c_{m+6}) + c_{m+6} u  (u = w^6-9). */
static void from_profile(fp12 *z, const uint64_t c[12][4]) {
    for (int i = 0; i < 48; i++) ((uint64_t *)z)[i] = 0;
    for (int m = 0; m < 6; m++) {
        int j = m & 1, k = m >> 1;
        fp nine, cbeta, cm, ninec, re;
        for (int t = 0; t < 4; t++) nine.v[t] = NINEMONT[t];
        for (int t = 0; t < 4; t++) cbeta.v[t] = c[m + 6][t];
        for (int t = 0; t < 4; t++) cm.v[t] = c[m][t];
        fmul(&ninec, &nine, &cbeta);
        fadd(&re, &cm, &ninec);
        fp2 a; a.c0 = re; a.c1 = cbeta;
        if (j == 0) z->c0.c[k] = a; else z->c1.c[k] = a;
    }
}
static int check_g(void) {
    fp12 base, one, acc;
    from_profile(&base, GCOEF);
    for (int i = 0; i < 48; i++) ((uint64_t *)&one)[i] = 0;
    for (int t = 0; t < 4; t++) ((uint64_t *)&one)[t] = RMONT[t];
    acc = one;
    int top = 0;
    for (int i = 255; i >= 0; i--) if ((REXP[i / 64] >> (i % 64)) & 1) { top = i; break; }
    for (int i = top; i >= 0; i--) {
        fp12 t; f12mul(&t, &acc, &acc); acc = t;
        if ((REXP[i / 64] >> (i % 64)) & 1) { fp12 u; f12mul(&u, &acc, &base); acc = u; }
    }
    int eq = 1, nz = 0;
    for (int i = 0; i < 48; i++) {
        if (((uint64_t *)&acc)[i] != ((uint64_t *)&one)[i]) eq = 0;
        if (((uint64_t *)&base)[i] != ((uint64_t *)&one)[i]) nz = 1;
    }
    printf("g^r == 1: %s ;  g != 1: %s\n", eq ? "YES" : "NO", nz ? "YES" : "NO");
    return eq && nz;
}


/* ---- production walk statistics on the real BN254 GT (mode 's') ---- */
static void f12pow_scalar(fp12 *out, const fp12 *base, const uint64_t *e4) {
    fp12 acc; for (int i = 0; i < 48; i++) ((uint64_t *)&acc)[i] = 0;
    for (int t = 0; t < 4; t++) ((uint64_t *)&acc)[t] = RMONT[t];
    int top = 0;
    for (int i = 255; i >= 0; i--) if ((e4[i / 64] >> (i % 64)) & 1) { top = i; break; }
    for (int i = top; i >= 0; i--) {
        fp12 t; f12mul(&t, &acc, &acc); acc = t;
        if ((e4[i / 64] >> (i % 64)) & 1) { fp12 u; f12mul(&u, &acc, base); acc = u; }
    }
    *out = acc;
}
static int walk_stats(long steps, int dp_shift) {
    fp12 g, z, x;                 /* g from the profile, z = g^12345 */
    from_profile(&g, GCOEF);
    uint64_t xsmall[4] = {12345ULL, 0, 0, 0};
    f12pow_scalar(&z, &g, xsmall);
    /* multiplier table M_j = g^{c_j} z^{d_j} */
    static fp12 M[32];
    for (int j = 0; j < 32; j++) {
        fp12 a, b;
        f12pow_scalar(&a, &g, INC_C[j]);
        f12pow_scalar(&b, &z, INC_D[j]);
        f12mul(&M[j], &a, &b);
    }
    /* The walk state must be a GT element: start at g^a z^b for fixed public labels. */
    uint64_t alab[4] = {0x1234567ULL, 0, 0, 0}, blab[4] = {0x89ABULL, 0, 0, 0};
    fp12 ga, zb, y;
    f12pow_scalar(&ga, &g, alab);
    f12pow_scalar(&zb, &z, blab);
    f12mul(&y, &ga, &zb);
    fp one_raw; one_raw.v[0] = 1; one_raw.v[1] = one_raw.v[2] = one_raw.v[3] = 0;
    fp canon; fmul(&canon, &y.c0.c[0].c0, &one_raw);   /* Montgomery residue -> canonical */
    long hist[32]; for (int i = 0; i < 32; i++) hist[i] = 0;
    long dps = 0; uint64_t dt_hi = 1ULL << (64 - dp_shift);   /* 128-bit window threshold = 2^(128-dp_shift) */
    for (long i = 0; i < steps; i++) {
        int j = (int)((canon.v[1] >> 56) & 31);        /* five middle bits of the canonical coefficient */
        fp12 t; f12mul(&t, &y, &M[j]); y = t;
        hist[j]++;
        fmul(&canon, &y.c0.c[0].c0, &one_raw);         /* one conversion per step, shared by both windows */
        if (canon.v[1] < dt_hi) dps++;                 /* low 128 bits of the same canonical coefficient */
    }
    double exp = (double)steps / 32.0, chi = 0.0;
    for (int i = 0; i < 32; i++) { double d = hist[i] - exp; chi += d * d / exp; }
    printf("{\"steps\": %ld, \"chi2_31\": %.2f, \"dp_hits\": %ld, \"dp_expected\": %.3f, \"dp_shift\": %d, \"histogram\": [",
           steps, chi, dps, (double)steps / (double)(1L << dp_shift), dp_shift);
    for (int i = 0; i < 32; i++) printf("%s%ld", i ? "," : "", hist[i]);
    printf("]}\n");
    return 0;
}

int main(int argc, char **argv) {
    uint64_t n0 = 1;
    for (int i = 0; i < 6; i++) n0 *= 2 - P[0] * n0;
    N0 = (uint64_t)0 - n0;


    if (argc > 1 && argv[1][0] == 'g') return check_g() ? 0 : 1;
    if (argc > 1 && argv[1][0] == 's') return walk_stats(atol(argv[2]), atoi(argv[3]));
    if (argc > 1 && argv[1][0] == 't') {
        fp12 x, y, z;
        uint64_t *px = (uint64_t *)&x, *py = (uint64_t *)&y;
        for (int i = 0; i < 48; i++) {
            px[i] = (uint64_t)((i * 0x9e3779b97f4a7c15ULL + 12345) & 0x0fffffffffffffffULL);
            py[i] = (uint64_t)((i * 0xc2b2ae3d27d4eb4fULL + 7) & 0x0fffffffffffffffULL);
        }
        f12mul(&z, &x, &y);
        uint64_t *pz = (uint64_t *)&z;
        for (int i = 0; i < 48; i++) printf("%016llx\n", (unsigned long long)pz[i]);
        return 0;
    }
    long iters = argc > 1 ? atol(argv[1]) : 4096;
    fp12 a, b;
    uint64_t *pa = (uint64_t *)&a;
    for (int i = 0; i < 48; i++) pa[i] = xr() & 0x0fffffffffffffffULL;
    b = a; b.c1.c[2].c0.v[0] ^= 1ULL;
    volatile uint64_t sink = 0;
    for (long i = 0; i < iters; i++) {
        fp12 c; f12mul(&c, &a, &b);
        b = c;
        sink ^= c.c0.c[0].c0.v[0];
    }
    printf("%llu\n", (unsigned long long)sink);
    return 0;
}
