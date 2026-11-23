"""boot_dp.py: production sample for the boot rule of descent_bn254_v6.py (distinct-norm O_K-primitive candidates).

For lattice j the target is T = g^e with e = random.Random(SEED0 + j).randrange(1, r) (uniform in [1, r-1], so T is
distributed as z*g^e for any z in the order-r subgroup); the boot lattice {(a, b) : b = R a mod p}, R = T_1/T_0, given by the
six rows (0, p e_i) followed by the six rows (e_k, R Z_K-coordinates mod p, centred) with weights (S_A |u|, U) in lowest
terms ((S_A, 1) at production u; the basis and weights of the input-uniform bound of KERNELS.md section 0.2), is
reduced by the same 2-stage Bareiss LLL as the attack, and the 144 fixed-rule vectors {b_k} U {b_i +- b_j} are listed
in index order (k, then i<j with + before -). A candidate is KEPT when it passes the pre-norm gate (all a-coordinates
< 2^128, all b-coordinates < 2^160; KERNELS.md section 4), is O_K-primitive and its exact
|N_g| = |N_K(u w a^2 - b^2)| differs from the |N_g| of every earlier candidate. sample(j) returns
[#distinct |N_g|, #primitive, kept indices, off, [number of kept pairs whose norms share a prime > 2^16, largest
bit length of such a shared part], s16] with off = log2|N_g| - 6 log2 p - 6 log2 u and s16 = log2 of the
part of |N_g| made of primes <= 2^16, both in units of 0.01 b, rounded to nearest (so the true value is within
0.005 b; descent_bn254_v6.py rounds success sizes up and survival sizes down by 0.005 b). boot_dp.b85 stores items
0-4 and boot_s16.b85 item 5 for lattices j = 0..n-1; replay(js) recomputes lattices js from scratch and compares."""
import base64,json,lzma,math,random
from math import log2
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC
SEED0=20261123
ZK3I=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]]
# g = ((3+y) - (1+2y)x)/((3+y) + (1+2y)x) ^ ((p^6+1)/r) for production, stored to save one 1,270-bit exponentiation;
# ctx(full=True) recomputes it and asserts equality (part "boot_dp" of the audit).
GG_PROD=[[8848624762422320342039834634586426602462006419422836673594466152942042939263,3953904918844897013355600426584044535472843156280848009193875963389324194241,9000031084490751966903281724280738519946796249428522219577568268430274370667,4158084932877046750858551790800463158328858469633640003738661023425764580467,6881999913871205692009270439592732937424986232541888714958876276586897236271,19804309580111041286957140466698868208415105794595284631199539753957958835504],[21498944321226087805386443261305559411103068419820109153752133157703483520502,18891827141438104355613965194300477026139509703679662316728538796643886137461,5604065960675545813525672771897121050395761811287632441961451180583116809804,16654477346818037882472273354750524721413278497121867099610410633311788391791,21740983272139554571848662512071575188866064072045298809389188042016135931170,9026018834305333961717328041149559038729054309629997075117026754871650427504]]
_CTX={}
def ctx(u=None,full=False):
 u=M.U_BN254 if u is None else u
 if u in _CTX and not full:return _CTX[u]
 p=36*u**4+36*u**3+24*u**2+6*u+1;r=36*u**4+36*u**3+18*u**2+6*u+1
 pm,pinv=M._pmul_h,M._pinv_h;W=M.W_COEFFS
 iuw=pinv([(u*c)%p for c in W],p)
 def mul(x,y):return([(a+b)%p for a,b in zip(pm(x[0],y[0],p),pm(pm(x[1],y[1],p),iuw,p))],[(a+b)%p for a,b in zip(pm(x[0],y[1],p),pm(x[1],y[0],p))])
 def inv(x):
  d=pinv([(a-b)%p for a,b in zip(pm(x[0],x[0],p),pm(pm(x[1],x[1],p),iuw,p))],p);return pm(x[0],d,p),pm([-c%p for c in x[1]],d,p)
 def pw(x,e):
  y=([1,0,0,0,0,0],[0]*6)
  for bt in bin(e)[2:]:
   y=mul(y,y)
   if bt=="1":y=mul(y,x)
  return y
 one=([1,0,0,0,0,0],[0]*6)
 canon=lambda:pw(mul(([3,1,0,0,0,0],[p-1,p-2,0,0,0,0]),inv(([3,1,0,0,0,0],[1,2,0,0,0,0]))),(p**6+1)//r)
 gg=(list(GG_PROD[0]),list(GG_PROD[1]))if u==M.U_BN254 and not full else canon()
 assert gg!=one and pw(gg,r)==one   # r is prime, so gg generates the order-r subgroup
 if full and u==M.U_BN254:assert gg==(list(GG_PROD[0]),list(GG_PROD[1]))
 G6=M.build_exact_int_geom(M.load_data()["t4_geom"])[0]
 from fractions import Fraction as Fr
 wf=Fr(12163341763596768109*abs(u),M.U_BN254)   # weights in lowest terms: (S_A, 1) at production u
 c=dict(u=u,p=p,r=r,pw=pw,gg=gg,G6=G6,i3=pow(3,-1,p),ref=6*log2(p)+6*log2(abs(u)),wa=wf.numerator,wb=wf.denominator)
 if not full:_CTX[u]=c
 return c
def lattice(c,j):
 """Reduced boot basis (12 rows) of stored production lattice j."""
 return lattice_e(c,random.Random(SEED0+j).randrange(1,c["r"]))
def lattice_e(c,e):
 """Reduced boot basis (12 rows) of the target T = g^e."""
 p=c["p"];pm,pinv=M._pmul_h,M._pinv_h;ZKN=M.ZK_NUM
 T=c["pw"](c["gg"],e)
 R=pm(T[1],pinv(T[0],p),p);B=[[p if t==6+k else 0 for t in range(12)]for k in range(6)]
 for k in range(6):
  bp=pm(R,[(x*c["i3"])%p for x in ZKN[k]],p);row=[1 if t==k else 0 for t in range(6)]
  for i in range(6):vv=sum(ZK3I[i][t]*bp[t]for t in range(6))%p;row.append(vv-p if vv>p//2 else vv)
  B.append(row)
 DESC.bareiss_lll_12d(B,c["G6"],c["wa"],c["wb"],False,11,40,(9,10))
 for row in B:assert M._pmul_h(R,_zkp(row[:6],p),p)==_zkp(row[6:],p)
 return B
def _zkp(v,p):
 return [(sum(v[k]*M.ZK_NUM[k][t]for k in range(6))*pow(3,-1,p))%p for t in range(6)]
def cands(B):
 return [B[i]for i in range(12)]+[[x+s*y for x,y in zip(B[i],B[j])]for i in range(12)for j in range(i+1,12)for s in(1,-1)]
def ng(c,v):
 ZKN,W,mh=M.ZK_NUM,M.W_COEFFS,M._mh_int
 an=[sum(v[k]*ZKN[k][t]for k in range(6))for t in range(6)];bn=[sum(v[6+k]*ZKN[k][t]for k in range(6))for t in range(6)]
 w=mh(W,mh(an,an));b2=mh(bn,bn);N=abs(M._norm_Kh_exact([c["u"]*w[t]-b2[t]for t in range(6)]))
 assert N%9**6==0;return N//9**6
_PP=[None]
def s16(N):
 if _PP[0] is None:
  a=list(M._primes_upto(1<<16))
  while len(a)>1:a=[a[i]*a[i+1]if i+1<len(a)else a[i]for i in range(0,len(a),2)]
  _PP[0]=a[0]
 x=N;g=math.gcd(_PP[0],x);sp=1
 while g>1:
  while x%g==0:x//=g;sp*=g
  g=math.gcd(g,x)
 return log2(sp)
def sample(j,u=None):
 c=ctx(u);seen=set();keep=[];nprim=0
 for k,v in enumerate(cands(lattice(c,j))):
  if max(abs(z)for z in v[:6]).bit_length()>=128 or max(abs(z)for z in v[6:]).bit_length()>=160:continue
  N=ng(c,v)
  if N==0:continue
  new=N not in seen;seen.add(N)
  if not M.is_ideal_primitive_oh(v[:6],v[6:]):continue
  nprim+=1
  if new:keep.append((k,round(100*(log2(N)-c["ref"])),round(100*s16(N)),N))
 sh=0;mx=0
 for i in range(len(keep)):
  for t in range(i+1,len(keep)):
   g=math.gcd(keep[i][3],keep[t][3])
   if g>1:
    g=_rough16(g)
    if g>1:sh+=1;mx=max(mx,g.bit_length())
 return [len(seen),nprim,[x[0]for x in keep],[x[1]for x in keep],[sh,mx],[x[2]for x in keep]]
def _rough16(g):
 """g with every prime <= 2^16 removed."""
 s16(1)
 t=math.gcd(_PP[0],g)
 while t>1:
  while g%t==0:g//=t
  t=math.gcd(t,g)
 return g
def gate_norm_bits():
 """Upper bound on log2|N_g| for every vector that passes the pre-norm gate (|a_k| < 2^128, |b_k| < 2^160): the same
 derivation as bt_ng_hi in verify_bn254_v6.py (Z_K multiplication tables, column norms, Hadamard over 6 embeddings)."""
 from fractions import Fraction as Fr
 zk=[[Fr(c,3)for c in row]for row in M.ZK_NUM]
 def pmf(a,b):
  r=[Fr(0)]*11
  for i in range(6):
   for t in range(6):r[i+t]+=a[i]*b[t]
  for d in range(10,5,-1):
   cc=r[d]
   if cc:
    for t in range(6):r[d-6+t]-=cc*M.H_COEFFS[t]
    r[d]=Fr(0)
  return r[:6]
 p2zk=lambda q:[int(sum(ZK3I[r][cc]*(3*q[cc])for cc in range(6)))//3 for r in range(6)]
 T=[[p2zk(pmf(zk[i],zk[t]))for t in range(6)]for i in range(6)]
 C_mul=max(sum(abs(T[i][t][k])for i in range(6)for t in range(6))for k in range(6))
 C_col2=math.sqrt(max(sum((sum(abs(T[i][t][k])for i in range(6)))**2 for k in range(6))for t in range(6)))
 w=p2zk([Fr(cc)for cc in M.W_COEFFS])
 op_w=max(sum(abs(sum(w[i]*T[i][t][k]for i in range(6)))for t in range(6))for k in range(6))
 return 6.0*(log2(C_col2)+log2(M.U_BN254*op_w*(2.0**256)+(2.0**320))+log2(C_mul))
def _ld(name):
 return json.loads(lzma.decompress(base64.b85decode((Path(__file__).resolve().parent/name).read_bytes())).decode())
_S=[None]
def load():
 if _S[0] is None:
  S=_ld("boot_dp.b85");S16=_ld("boot_s16.b85")["s16"]
  assert len(S16)==len(S["rec"])==S["n"]
  for r,s in zip(S["rec"],S16):assert len(s)==len(r[2])==len(r[3]);r.append(s)
  # record: [#distinct, #primitive, indices, off, [pairs sharing a prime > 2^16, max bits of that shared part], s16]
  _S[0]=S
 return _S[0]
def replay(js,full=False):
 if full:ctx(full=True)
 S=load();bad=[j for j in js if sample(j)!=S["rec"][j]]
 return dict(replayed=list(js),mismatches=bad,generator_recomputed=full,ok=not bad)
def alpha16(S=None):
 """Mean log2 content of primes <= 2^16 per kept vector vs a uniform integer (sum log2 q/(q-1)); the lattice is the
 sampling unit (ratio of per-lattice sums, linearised SE across lattices)."""
 S=S or load();sl=[sum(rec[5])/100.0 for rec in S["rec"]];kl=[len(rec[5])for rec in S["rec"]]
 n=len(sl);m=sum(sl)/sum(kl);kb=sum(kl)/n
 se=math.sqrt(sum((x-m*k)**2 for x,k in zip(sl,kl))/(n-1)/n)/kb
 er=sum(log2(q)/(q-1)for q in M._primes_upto(1<<16))
 # pair moment of the small-prime parts at the production slope k = 0.08 (descent): kept pairs of one lattice vs
 # pairs from different lattices, both exact over the stored records, ratios to the squared mean
 k=0.08;X=[[2.0**(k*x/100.0)for x in rec[5]]for rec in S["rec"]];N=sum(kl);mx=sum(map(sum,X))/N
 pw=[(sum(l)**2-sum(x*x for x in l))/mx**2 for l in X];dw=[len(l)*(len(l)-1)for l in X];rw=sum(pw)/sum(dw)
 sw=math.sqrt(sum((a-rw*b)**2 for a,b in zip(pw,dw))/(n-1)/n)/(sum(dw)/n)
 rb=((N*mx)**2-sum(sum(l)**2 for l in X))/mx**2/(N*N-sum(len(l)**2 for l in X))
 shp=sum(rec[4][0]for rec in S["rec"]);npairs=sum(len(rec[5])*(len(rec[5])-1)//2 for rec in S["rec"])
 return dict(kept_pairs=npairs,kept_pairs_sharing_a_prime_above_2_16=shp,max_shared_rough_bits=max(rec[4][1]for rec in S["rec"]),kept=sum(kl),lattices=n,mean_kept=round(kb,3),mean_s16=round(m,4),se_lattice=round(se,4),e_uniform=round(er,4),alpha_emp=round(er-m,4),alpha_emp_plus_3se=round(er-m+3*se,4),
 pair_moment_k008_within_lattice=round(rw,4),pair_moment_k008_within_se=round(sw,4),pair_moment_k008_between_lattices=round(rb,4))
