/* 12-dim special-q descent lattice sieve: LLL(delta=.99) + Schnorr-Euchner
* walk of L_p={c:sum w_i c_i=0 mod p} in the box |c_i|<=k, bucketed records,
* per-segment apply/scan.  C99, no deps.  Every command prints key=value. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define N 12
typedef int64_t i64; typedef uint64_t u64; typedef uint32_t u32;
typedef uint16_t u16; typedef __int128 i128; typedef unsigned __int128 u128;
typedef long double ld;
#define LLL_RESET 16
#define U (unsigned long long)
static u64 RNG;
static u64 rnd(void){u64 z=(RNG+=0x9E3779B97F4A7C15ULL);z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;z=(z^(z>>27))*0x94D049BB133111EBULL;return z^(z>>31);}
static u64 rbelow(u64 n){return rnd()%n;}
static u64 mm(u64 a,u64 b,u64 m){return(u64)((u128)a*b%m);}
static u64 pm(u64 a,u64 e,u64 m){u64 r=1%m;a%=m;while(e){if(e&1)r=mm(r,a,m);a=mm(a,a,m);e>>=1;}return r;}
static int isp(u64 n){static const u64 s[12]={2,3,5,7,11,13,17,19,23,29,31,37};
if(n<2)return 0;for(int i=0;i<12;i++){if(n==s[i])return 1;if(n%s[i]==0)return 0;}
u64 d=n-1;int k=0;while(!(d&1)){d>>=1;k++;}
for(int i=0;i<12;i++){u64 a=s[i]%n,x=pm(a,d,n);if(!a||x==1||x==n-1)continue;int ok=0;
for(int r=1;r<k;r++){x=mm(x,x,n);if(x==n-1){ok=1;break;}}if(!ok)return 0;}return 1;}
static u64 nprime(u64 n){if(n<2)return 2;n|=1;while(n<((u64)1<<63)&&!isp(n))n+=2;return n;}
static u64 invm(u64 a,u64 p){i128 t=0,nt=1,r=p,nr=a%p;while(nr){i128 q=r/nr,z=t-q*nt;t=nt;nt=z;z=r-q*nr;r=nr;nr=z;}
if(r!=1)return 0;t%=p;if(t<0)t+=p;return(u64)t;}
static u64 mod128(i128 x,u64 p){i128 r=x%p;if(r<0)r+=p;return(u64)r;}
typedef struct{u64 p,w[N];int perm[N];i128 B[N][N];ld MU[N][N],BN[N],BST[N][N];long swaps;}lat;
static int build(lat*L){u64 p=L->p;for(int i=0;i<N;i++)L->perm[i]=i;
if(L->w[0]%p==0){int j=-1;for(int i=1;i<N;i++)if(L->w[i]%p){j=i;break;}if(j<0)return 0;
u64 t=L->w[0];L->w[0]=L->w[j];L->w[j]=t;int q=L->perm[0];L->perm[0]=L->perm[j];L->perm[j]=q;}
u64 iv=invm(L->w[0]%p,p);if(!iv)return 0;memset(L->B,0,sizeof L->B);L->B[0][0]=p;
for(int i=1;i<N;i++){u64 c=(u64)((u128)(L->w[i]%p)*iv%p);L->B[i][0]=(i128)((p-c)%p);L->B[i][i]=1;}return 1;}
static void gso(const lat*L,ld MU[N][N],ld BN[N],ld BST[N][N]){
for(int i=0;i<N;i++){for(int l=0;l<N;l++)BST[i][l]=L->B[i][l];
for(int j=0;j<i;j++){ld z=0;for(int l=0;l<N;l++)z+=(ld)L->B[i][l]*BST[j][l];ld m=z/BN[j];MU[i][j]=m;
for(int l=0;l<N;l++)BST[i][l]-=m*BST[j][l];}
ld z=0;for(int l=0;l<N;l++)z+=BST[i][l]*BST[i][l];BN[i]=z;}}
static int lll(lat*L){ld MU[N][N],BN[N],BST[N][N];gso(L,MU,BN,BST);int k=1,since=0;long g=0;
while(k<N){if(++g>4000000L)return 0;
for(int j=k-1;j>=0;j--){i128 q=llroundl(MU[k][j]);
if(q){for(int i=0;i<N;i++)L->B[k][i]-=q*L->B[j][i];
for(int i=0;i<j;i++)MU[k][i]-=q*MU[j][i];MU[k][j]-=q;
if(++since>=LLL_RESET){gso(L,MU,BN,BST);since=0;}}}
if(BN[k]>=(0.99L-MU[k][k-1]*MU[k][k-1])*BN[k-1]){k++;continue;}
ld m=MU[k][k-1],A=BN[k]+m*m*BN[k-1],nm=m*BN[k-1]/A;
BN[k]=BN[k]*BN[k-1]/A;BN[k-1]=A;MU[k][k-1]=nm;
for(int i=0;i<N;i++){i128 t=L->B[k][i];L->B[k][i]=L->B[k-1][i];L->B[k-1][i]=t;}
for(int j=0;j<k-1;j++){ld t=MU[k][j];MU[k][j]=MU[k-1][j];MU[k-1][j]=t;}
for(int i=k+1;i<N;i++){ld t=MU[i][k];MU[i][k]=MU[i][k-1]-m*t;MU[i][k-1]=t+nm*MU[i][k];}
if(++since>=LLL_RESET){gso(L,MU,BN,BST);since=0;}
k=(k>1)?k-1:1;L->swaps++;}
gso(L,MU,BN,BST);memcpy(L->MU,MU,sizeof MU);memcpy(L->BN,BN,sizeof BN);memcpy(L->BST,BST,sizeof BST);return 1;}
static int verify(const lat*L){u64 p=L->p;
for(int i=0;i<N;i++){u64 s=0;for(int j=0;j<N;j++){u64 b=mod128(L->B[i][j],p);s=(u64)(((u128)s+(u128)b*(L->w[j]%p))%p);}if(s)return 0;}
ld q=1;for(int i=0;i<N;i++){if(!(L->BN[i]>0))return 0;q*=L->BN[i];}
ld p2=(ld)p*p;return fabsl(q-p2)/p2<=1e-9L;}
static u64 mix(u128 x){u64 a=(u64)x,b=(u64)(x>>64),h=a*0x9E3779B97F4A7C15ULL^(b+0x9E3779B97F4A7C15ULL+(a<<6)+(a>>2));
h^=h>>33;h*=0xff51afd7ed558ccdULL;h^=h>>33;return h;}
typedef struct{u128 seg;u32*off;u16*lp;u32 n,cap;}bk_t;
typedef struct{lat*L;u64 k,base,win;u128 pp[N];i64 lo[N],hi[N];ld ctr[N],G[N],L1[N],R2,half;
i128 pv[N];ld c[N];u64 hits,leaves,nodes;int stop;u64 target;u16 logp;
bk_t*bk;u32 nbk,bkcap;u32*map;u128*ms;u64 mcap,mused;}walk;
static void mrebuild(walk*W){memset(W->map,0,4*W->mcap);memset(W->ms,0,sizeof(u128)*W->mcap);
u64 mk=W->mcap-1;for(u32 b=0;b<W->nbk;b++){u64 h=mix(W->bk[b].seg)&mk;while(W->map[h])h=(h+1)&mk;W->map[h]=b+1;W->ms[h]=W->bk[b].seg;}}
static void badd(walk*W,u128 seg,u32 off){u64 mk=W->mcap-1,h=mix(seg)&mk;
while(W->map[h]){if(W->ms[h]==seg)break;h=(h+1)&mk;}
u32 bi;if(W->map[h])bi=W->map[h]-1;
else{if(W->nbk==W->bkcap){W->bkcap=W->bkcap?2*W->bkcap:1024;W->bk=realloc(W->bk,sizeof(bk_t)*W->bkcap);}
bi=W->nbk++;W->bk[bi]=(bk_t){seg,0,0,0,0};W->map[h]=bi+1;W->ms[h]=seg;W->mused++;
if(W->mused*10>W->mcap*7){W->mcap<<=1;W->map=realloc(W->map,4*W->mcap);W->ms=realloc(W->ms,sizeof(u128)*W->mcap);mrebuild(W);}}
bk_t*b=&W->bk[bi];if(b->n==b->cap){b->cap=b->cap?2*b->cap:4;b->off=realloc(b->off,4*b->cap);b->lp=realloc(b->lp,2*b->cap);}
b->off[b->n]=off;b->lp[b->n]=W->logp;b->n++;}
static int winit(walk*W,lat*L,u64 k,u64 win,u64 target){memset(W,0,sizeof*W);W->L=L;W->k=k;W->base=2*k+1;W->win=win;W->target=target;
W->logp=(u16)(log2((double)L->p)*16.0);
if(win==0){for(int i=0;i<N;i++){W->lo[i]=-(i64)k;W->hi[i]=(i64)k;}W->half=k;W->R2=(ld)N*k*k;for(int i=0;i<N;i++)W->ctr[i]=0;}
else{for(int i=0;i<N;i++){W->lo[i]=-(i64)k;W->hi[i]=-(i64)k+(i64)win-1;}W->half=((ld)win-1)/2;W->R2=(ld)N*W->half*W->half;
for(int i=0;i<N;i++)W->ctr[i]=(ld)(-(i64)k)+W->half;}
u128 pw[N],b=1;for(int i=0;i<N;i++){pw[i]=b;b*=(u128)W->base;}for(int i=0;i<N;i++)W->pp[i]=pw[L->perm[i]];
for(int j=0;j<N;j++){ld z=0,l=0;for(int i=0;i<N;i++){z+=W->ctr[i]*L->BST[j][i];l+=fabsl(L->BST[j][i]);}W->G[j]=z/L->BN[j];W->L1[j]=l;}
W->mcap=1;u64 want=(target?target:1)*4+16;while(W->mcap<want)W->mcap<<=1;
W->map=calloc(W->mcap,4);W->ms=calloc(W->mcap,sizeof(u128));W->bkcap=target?(u32)(target+16):1024;W->bk=calloc(W->bkcap,sizeof(bk_t));
return W->map&&W->ms&&W->bk;}
static void wfree(walk*W){for(u32 b=0;b<W->nbk;b++){free(W->bk[b].off);free(W->bk[b].lp);}free(W->bk);free(W->map);free(W->ms);}
static void emit(walk*W){u128 id=0;for(int i=0;i<N;i++)id+=(u128)(W->pv[i]-W->lo[i])*W->pp[i];badd(W,id>>16,(u32)(id&0xffff));}
static void wrec(walk*W,int j,ld acc){if(W->stop)return;
if(j<0){W->leaves++;for(int i=0;i<N;i++)if(W->pv[i]<W->lo[i]||W->pv[i]>W->hi[i])return;
W->hits++;emit(W);if(W->target&&W->hits>=W->target)W->stop=1;return;}
ld bn=W->L->BN[j],r2=W->R2-acc;if(r2<0)return;
ld ym=sqrtl(r2/bn),sl=W->half*W->L1[j]/bn;if(sl<ym)ym=sl;
ld m=W->c[j]-W->G[j];
long long xl=(long long)ceill(-ym-m)-1,xh=(long long)floorl(ym-m)+1;
if(xl<-100000000LL)xl=-100000000LL;if(xh>100000000LL)xh=100000000LL;W->nodes++;
for(long long x=xl;x<=xh;x++){ld y=x+m;
for(int i=0;i<N;i++)W->pv[i]+=(i128)x*W->L->B[j][i];
ld o[N];for(int l=0;l<j;l++){o[l]=W->c[l];W->c[l]+=x*W->L->MU[j][l];}
wrec(W,j-1,acc+y*y*bn);
for(int l=0;l<j;l++)W->c[l]=o[l];
for(int i=0;i<N;i++)W->pv[i]-=(i128)x*W->L->B[j][i];
if(W->stop)return;}}
static void wrun(walk*W){for(int i=0;i<N;i++){W->pv[i]=0;W->c[i]=0;}wrec(W,N-1,0);}
static u64 prun(walk*W,u32 th,int ap,int sc){u32*s=malloc(262144);if(!s)return 0;u64 sv=0;
for(u32 b=0;b<W->nbk;b++){memset(s,0,262144);bk_t*k=&W->bk[b];
if(ap)for(u32 t=0;t<k->n;t++)s[k->off[t]]+=k->lp[t];
if(sc)for(int i=0;i<65536;i++)if(s[i]>=th)sv++;}
free(s);return sv;}
static u64 brute(lat*L,const i64 lo[N],const i64 hi[N],u64 base,u128*out,u64 cap){u64 p=L->p,cnt=0;
i128 pw[N],b=1;for(int i=0;i<N;i++){pw[i]=b;b*=(i128)base;}
i64 c[N];for(int i=0;i<N;i++)c[i]=lo[i];
for(;;){u64 s=0;for(int i=0;i<N;i++)s=(u64)(((u128)s+(u128)(((i128)c[i]*L->w[i])%p+p)%p)%p);
if(!s){if(out&&cnt<cap){u128 id=0;for(int i=0;i<N;i++)id+=(u128)(c[i]-lo[i])*(u128)pw[L->perm[i]];out[cnt]=id;}cnt++;}
int i=0;for(;i<N;i++){if(c[i]<hi[i]){c[i]++;break;}c[i]=lo[i];}if(i==N)break;}
return cnt;}
static int cmpu(const void*a,const void*b){u128 x=*(const u128*)a,y=*(const u128*)b;return x<y?-1:x>y;}
static int check_one(u64 p,u64 k,u64 seed,u64 win){int fail=0,full=!win;RNG=seed;lat L;memset(&L,0,sizeof L);L.p=nprime(p);
for(int i=0;i<N;i++)L.w[i]=1+rbelow(L.p-1);
if(!build(&L)||!lll(&L)){printf("FAIL build p=%llu\n",U L.p);return 1;}
int v=verify(&L);if(!v){printf("FAIL latt p=%llu k=%llu s=%llu\n",U L.p,U k,U seed);fail=1;}
walk W;if(!winit(&W,&L,k,win,0)){printf("FAIL alloc\n");return 1;}wrun(&W);
u64 base=2*k+1,bad=0;long long wd=full?(long long)(2*k+1):(long long)win;
for(u32 b=0;b<W.nbk;b++)for(u32 t=0;t<W.bk[b].n;t++){u128 id=((u128)W.bk[b].seg<<16)|W.bk[b].off[t];
i64 D[N],cp[N];for(int e=0;e<N;e++){D[e]=(i64)(id%base);id/=base;}
for(int i=0;i<N;i++)cp[i]=D[L.perm[i]]-(i64)k;
u64 s=0;for(int i=0;i<N;i++)s=(u64)(((u128)s+(u128)(((i128)cp[i]*L.w[i])%L.p+L.p)%L.p)%L.p);
if(s)bad++;for(int i=0;i<N;i++)if(cp[i]<-(i64)k||cp[i]>-(i64)k+wd-1)bad++;}
if(bad){printf("FAIL memb p=%llu k=%llu s=%llu b=%llu\n",U L.p,U k,U seed,U bad);fail=1;}
u64 nr=0;for(u32 b=0;b<W.nbk;b++)nr+=W.bk[b].n;
u128*ids=malloc(sizeof(u128)*(nr?nr:1));u64 t2=0;
for(u32 b=0;b<W.nbk;b++)for(u32 t=0;t<W.bk[b].n;t++)ids[t2++]=((u128)W.bk[b].seg<<16)|W.bk[b].off[t];
qsort(ids,nr,sizeof(u128),cmpu);u64 dup=0;for(u64 i=1;i<nr;i++)if(ids[i]==ids[i-1])dup++;
if(dup){printf("FAIL dup p=%llu k=%llu s=%llu d=%llu\n",U L.p,U k,U seed,U dup);fail=1;}
int bok=1;u64 bc=0;
if(full?(2*k+1)<=5:win<=5){u128*bi=malloc(sizeof(u128)*(nr+8));bc=brute(&L,W.lo,W.hi,base,bi,nr+8);
if(bc!=nr){printf("FAIL cnt p=%llu k=%llu s=%llu b=%llu w=%llu\n",U L.p,U k,U seed,U bc,U nr);fail=1;bok=0;}
else{qsort(bi,bc,sizeof(u128),cmpu);for(u64 i=0;i<nr;i++)if(bi[i]!=ids[i]){bok=0;break;}
if(!bok){printf("FAIL set p=%llu k=%llu s=%llu\n",U L.p,U k,U seed);fail=1;}}
free(bi);}
double ex=pow((double)(full?2*k+1:win),12.0)/(double)L.p;double ra=ex>0?(double)nr/ex:0;
printf(" L p=%llu k=%llu s=%llu sw=%ld v=%s\n",U L.p,U k,U seed,L.swaps,v?"P":"F");
printf(" W p=%llu k=%llu w=%llu hits=%llu leaves=%llu exp=%.1f ratio=%.4f brute=%s\n",
U L.p,U k,U win,U nr,U W.leaves,ex,ra,(full?(2*k+1)<=5:win<=5)?(bok?"MATCH":"MISMATCH"):"skipped");
free(ids);wfree(&W);return fail;}
static u64 mklat(lat*L,u64 p,u64 seed){memset(L,0,sizeof*L);L->p=nprime(p);RNG=seed;
for(int i=0;i<N;i++)L->w[i]=1+rbelow(L->p-1);if(!build(L)||!lll(L))return 0;return L->p;}
static int c_setup(u64 p,u64 seed,u64 r){u64 ok=1,lp=0;
for(u64 i=0;i<r;i++){lat L;lp=mklat(&L,p,seed+i*0x9E3779B9ULL);ok&=verify(&L);}
printf("units=%llu hits=0 p=%llu ok=%llu\n",U r,U lp,U ok);return 0;}
static int c_walk(u64 p,u64 k,u64 seed,u64 r,u64 tg,u64 win){u64 th=0,tl=0,lp=0;
for(u64 i=0;i<r;i++){lat L;lp=mklat(&L,p,seed+i*0x9E3779B9ULL);walk W;if(!winit(&W,&L,k,win,tg))return 1;
wrun(&W);th+=W.hits;tl+=W.leaves;wfree(&W);}
printf("units=%llu hits=%llu leaves=%llu p=%llu k=%llu\n",U r,U th,U tl,U lp,U k);return 0;}
static int c_pass(u64 p,u64 k,u64 seed,u32 th,u64 r,u64 tg,u64 win,int mode){lat L;mklat(&L,p,seed);walk W;
if(!winit(&W,&L,k,win,tg))return 1;wrun(&W);u64 nr=0;for(u32 b=0;b<W.nbk;b++)nr+=W.bk[b].n;u64 sv=0;
for(u64 i=0;i<r;i++)sv=prun(&W,th,mode==0,mode<=1);
u64 cl=(u64)W.nbk*65536ULL;
printf("units=%llu records=%llu cells=%llu survivors=%llu\n",U r,U nr,U cl,U sv);
wfree(&W);return 0;}
static int c_wc(u64 p,u64 k,u64 sd,u64 r,u64 win){u64 P=nprime(p),th=0,tp=0;
for(u64 z=0;z<r;z++){RNG=sd+z*0x9E3779B9ULL;u64 w[N];for(int i=0;i<N;i++)w[i]=1+rbelow(P-1);
u64 iv=invm(w[11],P);i64 lo=-(i64)k,hi=win?-(i64)k+(i64)win-1:(i64)k,c[N];
for(int i=0;i<N;i++)c[i]=lo;u64 h=0,q=0;
for(;;){q++;u64 s=0;for(int i=0;i<11;i++)s=(u64)(((u128)s+mod128((i128)c[i]*(i128)w[i],P))%P);
i64 rr=(i64)((u64)((u128)((P-s)%P)*iv%P));if(rr>=lo&&rr<=hi)h++;if(rr-(i64)P>=lo&&rr-(i64)P<=hi)h++;
int i=0;for(;i<11;i++){if(c[i]<hi){c[i]++;break;}c[i]=lo;}if(i==11)break;}
th+=h;tp+=q;}
printf("units=%llu hits=%llu prefixes=%llu p=%llu k=%llu win=%llu\n",U r,U th,U tp,U P,U k,U win);return 0;}
static u64 au(const char*s){return strtoull(s,0,0);}
int main(int c,char**v){int a=1;if(c>1&&!strcmp(v[1],"bench"))a=2;if(a>=c)return 2;
const char*m=v[a];
if(!strcmp(m,"setup")&&a+3<c)return c_setup(au(v[a+1]),au(v[a+2]),au(v[a+3]));
if(!strcmp(m,"walk")&&a+4<c){u64 tg=a+5<c?au(v[a+5]):0,win=a+6<c?au(v[a+6]):0;
return c_walk(au(v[a+1]),au(v[a+2]),au(v[a+3]),au(v[a+4]),tg,win);}
if((!strcmp(m,"apply")||!strcmp(m,"scan")||!strcmp(m,"zero"))&&a+5<c){int md=!strcmp(m,"apply")?0:!strcmp(m,"scan")?1:2;
u64 tg=a+6<c?au(v[a+6]):0,win=a+7<c?au(v[a+7]):0;
return c_pass(au(v[a+1]),au(v[a+2]),au(v[a+3]),(u32)au(v[a+4]),au(v[a+5]),tg,win,md);}
if(!strcmp(m,"walkc")&&a+4<c){u64 win=a+5<c?au(v[a+5]):0;
return c_wc(au(v[a+1]),au(v[a+2]),au(v[a+3]),au(v[a+4]),win);}
if(!strcmp(m,"check")&&a+3<c){u64 p=au(v[a+1]),k=au(v[a+2]),sd=au(v[a+3]),win=a+4<c?au(v[a+4]):0;
int fb=a+5<c?atoi(v[a+5]):1;
int f=check_one(p,k,sd,fb?0:win);printf("%s p=%llu k=%llu seed=%llu win=%llu\n",f?"RESULT FAIL":"RESULT PASS",
U nprime(p),U k,U sd,U win);return f;}
return 2;}
