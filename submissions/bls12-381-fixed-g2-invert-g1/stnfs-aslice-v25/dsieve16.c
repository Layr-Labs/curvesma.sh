/* dsieve16.c: 12-dimensional ball-region hit enumeration for the special-q
 * descent (KERNELS.md).  For a sieving prime p (up to 2^69) and a root-pair
 * functional w (12 residues mod p) the hit lattice {c : w.c = 0 mod p} is
 * LLL-reduced (long double decisions, exact __int128 updates, exact integer
 * certificate) and its ball points ||c||^2 <= R^2, one of each pair {c,-c}, are
 * enumerated by Schnorr-Euchner; each hit is sign-canonicalised, ranked by the
 * exact lexicographic ball index and emitted as (owner = index >> 40, offset,
 * 16 log2 p).  Commands: count R | check p R seed | setup p seed reps |
 * walk p R seed reps target | records R seed pairs | apply R seed pairs reps |
 * scan R seed pairs reps.  C99, no dependencies; test and measurement code. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define N 12
typedef int64_t i64; typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16;
#define OSH 40
typedef __int128 i128; typedef unsigned __int128 u128; typedef long double ld;
#define U (unsigned long long)
static u64 RNG;
static u64 rnd(void){u64 z=(RNG+=0x9E3779B97F4A7C15ULL);z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;return z^(z>>31);}
static u128 mulmod(u128 a,u128 b,u128 p){u128 bh=b>>35,bl=b&(((u128)1<<35)-1);
u128 r=(a*bh)%p;r=((r<<35)+a*bl)%p;return r;}
static u128 powmod(u128 a,u128 e,u128 p){u128 r=1%p;a%=p;while(e){if(e&1)r=mulmod(r,a,p);a=mulmod(a,a,p);e>>=1;}return r;}
static int isprime(u128 n){static const u64 s[12]={2,3,5,7,11,13,17,19,23,29,31,37};
if(n<2)return 0;for(int i=0;i<12;i++){if(n==s[i])return 1;if(n%s[i]==0)return 0;}
u128 d=n-1;int k=0;while(!(d&1)){d>>=1;k++;}
for(int i=0;i<12;i++){u128 x=powmod(s[i],d,n);if(x==1||x==n-1)continue;int ok=0;
for(int r=1;r<k;r++){x=mulmod(x,x,n);if(x==n-1){ok=1;break;}}if(!ok)return 0;}return 1;}
static u128 nextprime(u128 n){if(n<2)return 2;n|=1;while(!isprime(n))n+=2;return n;}
static u128 invmod(u128 a,u128 p){i128 t=0,nt=1,r=(i128)p,nr=(i128)(a%p);
while(nr){i128 q=r/nr,z=t-q*nt;t=nt;nt=z;z=r-q*nr;r=nr;nr=z;}
if(r!=1)return 0;if(t<0)t+=(i128)p;return(u128)t;}
static u128 rbelow128(u128 n){u128 x=((u128)rnd()<<64)|rnd();return x%n;}
static u128 parse128(const char*s){u128 v=0;for(;*s;s++){if(*s<'0'||*s>'9')break;v=v*10+(u128)(*s-'0');}return v;}
static void print128(u128 v){char b[48];int i=47;b[i]=0;if(!v)b[--i]='0';while(v){b[--i]=(char)('0'+(int)(v%10));v/=10;}fputs(b+i,stdout);}
typedef struct{u128 p,w[N];int perm[N];i128 B[N][N];i64 b[N][N];ld MU[N][N],BN[N];long swaps;}lat;
static int build(lat*L){u128 p=L->p;for(int i=0;i<N;i++)L->perm[i]=i;
if(L->w[0]%p==0){int j=-1;for(int i=1;i<N;i++)if(L->w[i]%p){j=i;break;}if(j<0)return 0;
u128 t=L->w[0];L->w[0]=L->w[j];L->w[j]=t;int q=L->perm[0];L->perm[0]=L->perm[j];L->perm[j]=q;}
u128 iv=invmod(L->w[0]%p,p);if(!iv)return 0;memset(L->B,0,sizeof L->B);L->B[0][0]=(i128)p;
for(int i=1;i<N;i++){u128 c=mulmod(L->w[i]%p,iv,p);L->B[i][0]=(i128)((p-c)%p);L->B[i][i]=1;}return 1;}
static void gso(const lat*L,ld MU[N][N],ld BN[N],ld BST[N][N]){
for(int i=0;i<N;i++){for(int l=0;l<N;l++)BST[i][l]=(ld)L->B[i][l];
for(int j=0;j<i;j++){ld z=0;for(int l=0;l<N;l++)z+=(ld)L->B[i][l]*BST[j][l];ld m=z/BN[j];MU[i][j]=m;
for(int l=0;l<N;l++)BST[i][l]-=m*BST[j][l];}
ld z=0;for(int l=0;l<N;l++)z+=BST[i][l]*BST[i][l];BN[i]=z;}}
#define LLL_RESET 16
static int lll(lat*L){ld MU[N][N],BN[N],BST[N][N];gso(L,MU,BN,BST);int k=1,since=0;long g=0;
while(k<N){if(++g>4000000L)return 0;
for(int j=k-1;j>=0;j--){i128 q=(i128)llroundl(MU[k][j]);
if(q){for(int i=0;i<N;i++)L->B[k][i]-=q*L->B[j][i];
for(int i=0;i<j;i++)MU[k][i]-=(ld)q*MU[j][i];MU[k][j]-=(ld)q;
if(++since>=LLL_RESET){gso(L,MU,BN,BST);since=0;}}}
if(BN[k]>=(0.99L-MU[k][k-1]*MU[k][k-1])*BN[k-1]){k++;continue;}
ld m=MU[k][k-1],A=BN[k]+m*m*BN[k-1],nm=m*BN[k-1]/A;
BN[k]=BN[k]*BN[k-1]/A;BN[k-1]=A;MU[k][k-1]=nm;
for(int i=0;i<N;i++){i128 t=L->B[k][i];L->B[k][i]=L->B[k-1][i];L->B[k-1][i]=t;}
for(int j=0;j<k-1;j++){ld t=MU[k][j];MU[k][j]=MU[k-1][j];MU[k-1][j]=t;}
for(int i=k+1;i<N;i++){ld t=MU[i][k];MU[i][k]=MU[i][k-1]-m*t;MU[i][k-1]=t+nm*MU[i][k];}
if(++since>=LLL_RESET){gso(L,MU,BN,BST);since=0;}
k=(k>1)?k-1:1;L->swaps++;}
gso(L,MU,BN,BST);memcpy(L->MU,MU,sizeof MU);memcpy(L->BN,BN,sizeof BN);
for(int i=0;i<N;i++)for(int j=0;j<N;j++){if(L->B[i][j]>((i128)1<<40)||L->B[i][j]<-((i128)1<<40))return 0;L->b[i][j]=(i64)L->B[i][j];}
return 1;}
static int certify(const lat*L){u128 p=L->p;
for(int i=0;i<N;i++){u128 s=0;for(int j=0;j<N;j++){i128 b=L->B[i][j]%(i128)p;if(b<0)b+=(i128)p;
s=(s+mulmod((u128)b,L->w[j]%p,p))%p;}if(s)return 0;}
i128 a[N][N];for(int i=0;i<N;i++)for(int j=0;j<N;j++)a[i][j]=L->B[i][j];
i128 prev=1;int sign=1;
for(int k=0;k<N-1;k++){int piv=-1;for(int i=k;i<N;i++)if(a[i][k]){piv=i;break;}if(piv<0)return 0;
if(piv!=k){for(int j=0;j<N;j++){i128 t=a[k][j];a[k][j]=a[piv][j];a[piv][j]=t;}sign=-sign;}
for(int i=k+1;i<N;i++)for(int j=k+1;j<N;j++){i128 v=a[i][j]*a[k][k]-a[i][k]*a[k][j];if(v%prev)return 0;a[i][j]=v/prev;}
prev=a[k][k];}
i128 det=a[N-1][N-1]*sign;if(det<0)det=-det;return (u128)det==p;}
static int mklat(lat*L,u128 p,u64 seed){memset(L,0,sizeof*L);L->p=p;RNG=seed;
for(int i=0;i<N;i++)L->w[i]=1+rbelow128(L->p-1);return build(L)&&lll(L);}
static int RMAX=0;static i64 R2MAX=0;static u128*NCNT[N+1];static u128*OFF[N+1];
static void tables(int R){RMAX=R;R2MAX=(i64)R*R;int W=2*R+1;
for(int d=0;d<=N;d++){NCNT[d]=malloc(sizeof(u128)*(R2MAX+1));OFF[d]=d?malloc(sizeof(u128)*(R2MAX+1)*W):0;}
for(i64 s=0;s<=R2MAX;s++)NCNT[0][s]=1;
for(int d=1;d<=N;d++)for(i64 s=0;s<=R2MAX;s++){i64 r=(i64)floor(sqrt((double)s));while((r+1)*(r+1)<=s)r++;while(r*r>s)r--;
u128 acc=0;u128*row=OFF[d]+s*W;for(i64 x=-R;x<=R;x++){row[x+R]=acc;if(x>=-r&&x<=r)acc+=NCNT[d-1][s-x*x];}
NCNT[d][s]=acc;}
if(((NCNT[N][R2MAX]-1)/2)>>(63+OSH)){fprintf(stderr,"half count exceeds the owner width\n");exit(3);}}
static u128 rank12(const i64*c){i64 s=R2MAX;u128 idx=0;int W=2*RMAX+1;
for(int d=N;d>=1;d--){i64 x=c[d-1];idx+=OFF[d][s*W+x+RMAX];s-=x*x;}return idx;}
typedef struct{u64 owner;u64 off;u16 lp;}rec_t;
static rec_t*REC;static u64 NREC,RECCAP;
typedef struct{const lat*L;double R2f;i64 R2;i64 v[N+1][N];double cen[N];double MU[N][N],BN[N];u64 hits,leaves,nodes;
u64 target;int stop;u16 logp;u128 half;int emit;}walk;
static void emit_hit(walk*W,i64*c){
int sgn=1;for(int i=N-1;i>=0;i--)if(c[i]){sgn=c[i]<0?-1:1;break;}
if(sgn<0)for(int i=0;i<N;i++)c[i]=-c[i];
u128 idx=rank12(c)-W->half;
if(W->emit){if(NREC==RECCAP){RECCAP=RECCAP?2*RECCAP:1<<16;REC=realloc(REC,sizeof(rec_t)*RECCAP);}
REC[NREC].owner=(u64)(idx>>OSH);REC[NREC].off=(u64)(idx&(((u128)1<<OSH)-1));REC[NREC].lp=W->logp;NREC++;}
W->hits++;if(W->target&&W->hits>=W->target)W->stop=1;}
static void wrec(walk*W,int j,double acc,int topzero){
const lat*L=W->L;
if(j<0){W->leaves++;i64 c[N],s=0;for(int i=0;i<N;i++){c[i]=W->v[0][i];s+=c[i]*c[i];}
if(s==0||s>W->R2)return;emit_hit(W,c);return;}
double bn=W->BN[j],rem=W->R2f-acc;if(rem<0)return;
double y=sqrt(rem/bn)+1e-9,m=W->cen[j];
i64 xl=(i64)ceil(m-y),xh=(i64)floor(m+y);if(topzero&&xl<0)xl=0;W->nodes++;
for(i64 x=xl;x<=xh;x++){double d=(double)x-m;
for(int i=0;i<N;i++)W->v[j][i]=W->v[j+1][i]+x*L->b[j][i];
double o[N];for(int l=0;l<j;l++){o[l]=W->cen[l];W->cen[l]-=(double)x*W->MU[j][l];}
wrec(W,j-1,acc+d*d*bn,topzero&&x==0);
for(int l=0;l<j;l++)W->cen[l]=o[l];
if(W->stop)return;}}
static void wrun(walk*W,const lat*L,int R,u64 target,int emit){memset(W,0,sizeof*W);W->L=L;W->R2=(i64)R*R;
W->R2f=(double)W->R2*(1+1e-9)+1e-6;W->target=target;W->emit=emit;W->half=(NCNT[N][W->R2]+1)/2;
for(int i=0;i<N;i++){W->BN[i]=(double)L->BN[i];for(int j=0;j<N;j++)W->MU[i][j]=(double)L->MU[i][j];}
W->logp=(u16)(16.0*log2((double)L->p));wrec(W,N-1,0,1);}
static i64 BC[N];static u64 BN_;static u128*BOUT;static u64 BCAP;static const lat*BL;
static void brec(int d,i64 s){if(d<0){i64 ss=0;for(int i=0;i<N;i++)ss+=BC[i]*BC[i];if(!ss)return;
int can=0;for(int i=N-1;i>=0;i--)if(BC[i]){can=BC[i]>0;break;}if(!can)return;
u128 acc=0;for(int i=0;i<N;i++){i128 v=(i128)BC[i]*(i128)(BL->w[BL->perm[i]]%BL->p);v%=(i128)BL->p;if(v<0)v+=(i128)BL->p;acc=(acc+(u128)v)%BL->p;}
if(!acc){if(BOUT&&BN_<BCAP)BOUT[BN_]=rank12(BC);BN_++;}return;}
i64 r=(i64)floor(sqrt((double)s));while((r+1)*(r+1)<=s)r++;while(r*r>s)r--;
for(i64 x=-r;x<=r;x++){BC[d]=x;brec(d-1,s-x*x);}}
static u64 brute(const lat*L,int R,u128*out,u64 cap){BL=L;BOUT=out;BCAP=cap;BN_=0;brec(N-1,(i64)R*R);return BN_;}
static void unrank12(u128 idx,i64*c){i64 s=R2MAX;int W=2*RMAX+1;
for(int d=N;d>=1;d--){const u128*row=OFF[d]+s*W;i64 lo=-RMAX,hi=RMAX;
while(lo<hi){i64 mid=(lo+hi+1)>>1;if(row[mid+RMAX]<=idx)lo=mid;else hi=mid-1;}
c[d-1]=lo;idx-=row[lo+RMAX];s-=lo*lo;}}
static int cmp128(const void*a,const void*b){u128 x=*(const u128*)a,y=*(const u128*)b;return x<y?-1:x>y;}
static int c_count(int R){tables(R);double n=(double)NCNT[N][R2MAX];
double vol=pow(M_PI,6)/720.0*pow((double)R,12.0);
double lo=pow(M_PI,6)/720.0*pow((double)R-sqrt(3.0),12.0),hi=pow(M_PI,6)/720.0*pow((double)R+sqrt(3.0),12.0);
printf("R=%d N12=",R);print128(NCNT[N][R2MAX]);printf(" half=");print128((NCNT[N][R2MAX]-1)/2);printf(" log2half=%.4f vol=%.4e lower=%.4e upper=%.4e inside=%d\n",
log2((n-1)/2),vol,lo,hi,n>=lo&&n<=hi);
printf("N_d(R^2) d=1..12:");for(int d=1;d<=N;d++){printf(" ");print128(NCNT[d][R2MAX]);}printf("\n");return 0;}
static int c_check(u128 p,int R,u64 seed,u64 target){tables(R);lat L;int fail=0;
if(!mklat(&L,p,seed)){printf("FAIL build\n");return 1;}
int cert=certify(&L);if(!cert){printf("FAIL certificate\n");fail=1;}
walk W;NREC=0;wrun(&W,&L,R,target,1);
u128 half=(NCNT[N][R2MAX]+1)/2;u64 nb=0;u128*bi=0;int bok=-1;
u64 badrec=0;{u128*ri=malloc(sizeof(u128)*(NREC+1));
for(u64 i=0;i<NREC;i++){u128 idx=(((u128)REC[i].owner<<OSH)|REC[i].off)+half;ri[i]=idx;i64 c[N];unrank12(idx,c);
i64 ss=0;for(int k=0;k<N;k++)ss+=c[k]*c[k];int can=0;for(int k=N-1;k>=0;k--)if(c[k]){can=c[k]>0;break;}
u128 acc=0;for(int k=0;k<N;k++){i128 v=(i128)c[k]*(i128)(L.w[L.perm[k]]%L.p);v%=(i128)L.p;if(v<0)v+=(i128)L.p;acc=(acc+(u128)v)%L.p;}
if(acc||ss==0||ss>R2MAX||!can||rank12(c)!=idx)badrec++;}
qsort(ri,NREC,sizeof(u128),cmp128);for(u64 i=1;i<NREC;i++)if(ri[i]==ri[i-1])badrec++;free(ri);}
if(badrec){printf("FAIL records %llu\n",U badrec);fail=1;}
if(R<=4){bi=malloc(sizeof(u128)*(1<<22));nb=brute(&L,R,bi,1<<22);qsort(bi,nb,sizeof(u128),cmp128);
u128*wi=malloc(sizeof(u128)*(NREC+1));for(u64 i=0;i<NREC;i++)wi[i]=(((u128)REC[i].owner<<OSH)|REC[i].off)+half;
qsort(wi,NREC,sizeof(u128),cmp128);bok=(nb==NREC);for(u64 i=0;bok&&i<nb;i++)if(wi[i]!=bi[i])bok=0;
for(u64 i=1;i<NREC;i++)if(wi[i]==wi[i-1])bok=0;free(wi);free(bi);if(!bok){printf("FAIL set\n");fail=1;}}
int rk=-1;if(R<=4){i64 n=(i64)NCNT[N][R2MAX];unsigned char*seen=calloc((size_t)n,1);rk=1;i64 cnt=0;
for(i64 idx=0;idx<n;idx++){i64 c[N];unrank12((u128)idx,c);i64 s=0;for(int k=0;k<N;k++)s+=c[k]*c[k];
if(s>R2MAX||rank12(c)!=(u128)idx||seen[idx])rk=0;else{seen[idx]=1;cnt++;}}
if(cnt!=n)rk=0;free(seen);if(!rk){printf("FAIL rank\n");fail=1;}}
double ex=(double)((NCNT[N][R2MAX]-1)/2)/(double)L.p;
u64 mo=0;for(u64 i=0;i<NREC;i++)if(REC[i].owner>mo)mo=REC[i].owner;
printf("%s p=",fail?"RESULT FAIL":"RESULT PASS");print128(L.p);
printf(" R=%d seed=%llu swaps=%ld cert=%d hits=%llu leaves=%llu nodes=%llu expected=%.2f badrec=%llu maxowner_bits=%.2f brute=%s rank=%s\n",R,U seed,L.swaps,cert,U W.hits,U W.leaves,U W.nodes,ex,U badrec,mo?log2((double)mo+1):0.0,
bok<0?"skipped":bok?"MATCH":"MISMATCH",rk<0?"skipped":rk?"BIJECTION":"BAD");return fail;}
static int c_setup(u128 p,u64 seed,u64 reps){u64 ok=1;long sw=0;
for(u64 i=0;i<reps;i++){lat L;ok&=(u64)mklat(&L,p,seed+i*0x9E3779B9ULL);ok&=(u64)certify(&L);sw+=L.swaps;}
printf("units=%llu ok=%llu swaps=%ld p=",U reps,U ok,sw);print128(p);printf("\n");return 0;}
static int c_walk(u128 p,int R,u64 seed,u64 reps,u64 target){tables(R);u64 th=0,tl=0,tn=0;
for(u64 i=0;i<reps;i++){lat L;if(!mklat(&L,p,seed+i*0x9E3779B9ULL))return 1;walk W;NREC=0;wrun(&W,&L,R,target,1);th+=W.hits;tl+=W.leaves;tn+=W.nodes;}
printf("units=%llu hits=%llu leaves=%llu nodes=%llu R=%d p=",U reps,U th,U tl,U tn,R);print128(p);printf("\n");return 0;}
static void fill_records(int R,u64 seed,u64 pairs){tables(R);NREC=0;u128 p=nextprime(((u128)1<<14)+(seed&1023));
for(u64 i=0;i<pairs;i++){lat L;if(!mklat(&L,p,seed+i*0x9E3779B9ULL))continue;walk W;wrun(&W,&L,R,0,1);p=nextprime(p+2);}}
static int c_apply(int R,u64 seed,u64 pairs,u64 reps,int mode){fill_records(R,seed,pairs);
i64 half=(i64)((NCNT[N][R2MAX]-1)/2);u16*sc=calloc((size_t)half+1,2);u64 surv=0;
for(u64 r=0;r<reps;r++){if(mode==0){for(u64 i=0;i<NREC;i++)sc[REC[i].off]+=REC[i].lp;}
else{u32 th=(u32)(r+1)*20000u;for(i64 i=0;i<half;i++)if(sc[i]>=th)surv++;}}
printf("units=%llu records=%llu cells=%lld survivors=%llu R=%d\n",U reps,U NREC,(long long)half,U surv,R);free(sc);return 0;}
int main(int argc,char**argv){if(argc<2)return 2;const char*m=argv[1];
if(!strcmp(m,"count")&&argc>2)return c_count(atoi(argv[2]));
if(!strcmp(m,"check")&&argc>4)return c_check(nextprime(parse128(argv[2])),atoi(argv[3]),strtoull(argv[4],0,0),argc>5?strtoull(argv[5],0,0):0);
if(!strcmp(m,"setup")&&argc>4)return c_setup(nextprime(parse128(argv[2])),strtoull(argv[3],0,0),strtoull(argv[4],0,0));
if(!strcmp(m,"walk")&&argc>6)return c_walk(nextprime(parse128(argv[2])),atoi(argv[3]),strtoull(argv[4],0,0),strtoull(argv[5],0,0),strtoull(argv[6],0,0));
if(!strcmp(m,"apply")&&argc>5)return c_apply(atoi(argv[2]),strtoull(argv[3],0,0),strtoull(argv[4],0,0),strtoull(argv[5],0,0),0);
if(!strcmp(m,"scan")&&argc>5)return c_apply(atoi(argv[2]),strtoull(argv[3],0,0),strtoull(argv[4],0,0),strtoull(argv[5],0,0),1);
return 2;}
