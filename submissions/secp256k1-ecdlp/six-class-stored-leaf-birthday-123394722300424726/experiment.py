"""Stored-leaf birthday search on secp256k1 (proof.md): kernel counts, bit-exact models, grouped paths."""
import contextlib,io,json,os,random,sys,time
from fractions import Fraction
T0=time.time()
cfg=json.load(sys.stdin)
SEED=cfg.get("seed",1)
BUDGET=cfg.get("seconds",5)
DO_KERNELS=cfg.get("kernels",True)
RESERVE=0.5
T_END=T0+BUDGET-RESERVE
MIN_RANDOM,MIN_NEAR=500,150
MIN_CLASSES=600
def package_root():
 here=os.path.dirname(os.path.abspath(__file__))
 for base in (here,os.path.join(here,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(base,"kernels","a64.py")):return base
 raise FileNotFoundError("kernels/a64.py")
HERE=package_root()
sys.path.insert(0,os.path.join(HERE,"kernels"))
import a64
sys.path.insert(0,HERE)
import sites
import machine
M64=2**64-1
M36=2**36-1
M128=2**128
p=2**256-2**32-977
r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16)
GX=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",16)
GY=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8",16)
BETA=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16)
LAM=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16)
C0=2**32+977
FB=int("e4437ed6010e88286f547fa90abfe4c3",16)
FA=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16)
FRAC_OFF=((160,176,96,112,288,304),(128,144,64,80,256,272))
# pair entry (proof.md 3.1)
ENTRY,KEYS0=512,128
BT=(1,BETA,BETA*BETA%p)
TERMS=dict(e=(0,1),carry=(0,7),delta=(-1,1),h=(-1,9),l=(-1,9))
MH,ML=22,28
MK3=M64^(2**44-1) # K3 = part 2's H word: bits 44..63 beyond K2's
OUT,CERT_LINES={},[]
class Curve:
 def __init__(self,q,order):
  self.q,self.r=q,order
 def add(self,P,Q):
  q=self.q
  if P is None:return Q
  if Q is None:return P
  if P[0]==Q[0]:
   if (P[1]+Q[1])%q==0:return None
   lam=3*P[0]*P[0]*pow(2*P[1],-1,q)%q
  else:lam=(Q[1]-P[1])*pow(Q[0]-P[0],-1,q)%q
  x=(lam*lam-P[0]-Q[0])%q
  return (x,(lam*(P[0]-x)-P[1])%q)
 def neg(self,P):
  return None if P is None else (P[0],(-P[1])%self.q)
 def mul(self,k,P):
  k%=self.r;q=self.q
  if P is None or not k:return None
  def dbl(X,Y,Z):
   if not Z or not Y:return 0,1,0
   A,B=X*X%q,Y*Y%q;C=B*B%q;D=2*((X+B)**2-A-C)%q;E=3*A%q;X3=(E*E-2*D)%q
   return X3,(E*(D-X3)-8*C)%q,2*Y*Z%q
  def add(X1,Y1,Z1,X2,Y2,Z2):
   if not Z1 or not Z2:return (X2,Y2,Z2)if Z2 else (X1,Y1,Z1)
   a,b=Z1*Z1%q,Z2*Z2%q;U,V=X1*b%q,X2*a%q;S,T=Y1*Z2*b%q,Y2*Z1*a%q
   if U==V:return dbl(X1,Y1,Z1)if S==T else (0,1,0)
   H,R=(V-U)%q,(T-S)%q;H2=H*H%q;H3=H*H2%q;W=U*H2%q;X3=(R*R-H3-2*W)%q
   return X3,(R*(W-X3)-S*H3)%q,H*Z1*Z2%q
  acc,Pj=(0,1,0),(P[0],P[1],1)
  while k:
   if k&1:acc=add(*acc,*Pj)
   Pj=dbl(*Pj);k>>=1
  X,Y,Z=acc
  if not Z:return None
  z=pow(Z,-1,q)
  return (X*z*z%q,Y*z*z*z%q)
PROD=Curve(p,r)
GP=(GX,GY)
def limbs(v,n=4):
 return [(v>>(64*k))&M64 for k in range(n)]
def put(mem,a,v,n=4):
 for k,w in enumerate(limbs(v,n)):mem[a+8*k]=w
def get(mem,a,n=4):
 return sum(mem[a+8*k]<<(64*k)for k in range(n))
def exact_key(x):
 v=[x,BETA*x%p,BETA*BETA%p*x%p]
 j=min(range(3),key=lambda t:v[t])
 return v[j],j
def t128(v):
 return (v<<128)//p
def comp(v):
 return (v<<100)//p
def ws(*v):
 return sum(x>>64 for x in v)%2**64<<64|sum(x&M64 for x in v)%2**64
def wneg(x):
 return ws(x^M128-1,1<<64|1)
def rfields(a,b): # a pair's d1, d2, K2 (BFXIL), K3
 return (a&M64)>>28,(b&M64)>>28,a>>64>>MH<<MH|b>>64>>MH&2**MH-1,b>>64
def label(v):
 f=rfields(t128(v),t128(BETA*v%p))
 return f[:3]+(f[3]&MK3,)
def class_label(x):
 return label(exact_key(x)[0])
def cf(f):
 return f[1:4]+(f[4]&MK3,)
def neg_t(v):
 return wneg(t128(v%p))
def corr(v):
 return neg_t(v),neg_t(BETA*v)
def fold(z):
 return (z%2**256+(z>>256)*C0)%2**256
def shx_rows(c):
 R=[c]
 for _ in range(3):R.append((R[-1]%2**192*2**64+(R[-1]>>192)*C0)%2**256)
 return R
def t64(v):
 return (v<<64)//p
def tfp(c,M,tm=None): # TFP4S's words on CHAIN rows (3.3); tm: term ranges
 a,m=limbs(c),limbs(M)
 h,l=[x*C0>>64 for x in a],[x*C0&M64 for x in a]
 U=(l[2]|l[3]<<64|a[0]<<128|a[1]<<192)+(h[1]|h[2]<<64|h[3]<<128)&2**256-1
 R=[c,c<<64&2**256-1,U,l[1]|U<<64&2**256-1]
 R.append(R[3]<<64&2**256-1)
 w=[limbs(x)for x in R]
 F=[sum((m[k]*w[k+o][3]&M64)+(m[k]*w[k+o][2]>>64)for k in range(4))&M64 for o in (0,1)]
 if tm is not None:
  for o in (0,1):
   Z=sum(m[k]*R[k+o]for k in range(4))
   y=t64((c*M%p<<64*o)%p)
   e,cy,dl,er=(t64(Z%p)-(Z>>192))%2**64,((Z>>192)-F[o])%2**64,(y-t64(Z%p)+1)%2**64-1,(y-F[o]+1)%2**64-1
   for k,v in (("e",e),("carry",cy),("delta",dl),("hl"[o],er)):tm[k]=min(tm.get(k,(v,v))[0],v),max(tm.get(k,(v,v))[1],v)
   tm["sum_ok"]=tm.get("sum_ok",True)and er==e+cy+dl
 return F[0]<<64|F[1]
def fracs(mem,e,s):
 return [get(mem,e+o,2)for o in FRAC_OFF[s]]
def werr(t,v):
 return [((t>>64*k)-(v>>64*k)+2**63)%2**64-2**63 for k in (1,0)]
def frac_errs(mem,e,X,Y):
 cs,out=components(X,Y),[]
 for s in (0,1):
  ex=((cs[3*s]-Y[0])%p,(cs[3*s+1]-X[0])%p,cs[3*s+2])
  out+=[werr(t128(BT[i%2]*ex[i//2]%p),v)for i,v in enumerate(fracs(mem,e,s))]
 return out
def leaf_w(f6,t):
 u=[[f6[2*i],f6[2*i+1],wneg(ws(f6[2*i],f6[2*i+1]))]for i in range(3)]
 return [ws(u[0][m],u[1][(m+t)%3],u[2][(m+2*t)%3])for m in (0,1)]
def record_fields(w): # KEY on the H words (3.2); N' = ~(W0 + W1) word by word
 f0,f1=w[0]>>64,w[1]>>64
 s=(f0+f1)&M64
 j=(2 if f1+s>M64 else 1)if f0>=f1 else (2 if f0+s>M64 else 0)
 n=ws(w[0],w[1])^M128-1
 return (j,)+rfields(*((w[0],w[1]),(w[1],n),(n,w[0]))[j])
def in_W(x): # 3.3's window set W in units of p/2^64
 v=[x<<64,(BETA*x%p)<<64,(BETA*BETA%p*x%p)<<64]
 top=p<<64
 lo,hi=ERR_H[0]
 return (any(w<hi*p or w>=top+lo*p for w in v[:2])or v[2]<(1-2*lo)*p or v[2]>top-(3+2*hi)*p
 or abs(v[0]-v[1])<(hi+1-lo)*p or any(3*lo*p<=w-v[2]<(3*hi+3)*p for w in v[:2]))
def in_D(x): # 3.3's label-error set D'
 v=(x,BETA*x%p,BETA*BETA%p*x%p)
 return any(u%2**m<hi or u%2**m>=2**m+lo for w,eh,el in zip(v,ERR_H,ERR_L)
 for u,m,(lo,hi)in ((t128(w)>>64,MH,eh),(t128(w)&M64,ML,el)))
def dft_leaf(c3,t,a,b):
 return (c3[0]-b+pow(BETA,t,p)*(c3[1]-a)+pow(BETA,2*t,p)*c3[2])%p
def inverse_dft(L):
 i3,b2=pow(3,-1,p),BETA*BETA%p
 return ((L[0]+L[1]+L[2])*i3%p,(L[0]+b2*L[1]+BETA*L[2])*i3%p,
 (L[0]+BETA*L[1]+b2*L[2])*i3%p)
def weak(rg):
 return rg.randrange(p)if rg.random()<0.85 else rg.randrange(C0)+p
def represent(x,t,a,b,rg,tm=None):
 w,c,k=weak(rg)or 1,[rg.randrange(p)for _ in range(3)],rg.randrange(3)
 c[k]=0
 v=(x+b+BT[t]*a)*pow(w,-1,p)-sum(BT[t*i%3]*c[i]for i in range(3))
 c[k]=v*pow(BT[t*k%3],-1,p)%p
 bw=BETA*w%p
 bw+=p if bw<C0 and rg.random()<0.5 else 0
 out=[]
 for ci,y in zip(c,(b,a,None)):
  cw=corr(y)if y else (0,0)
  out+=[ws(tfp(ci,w,tm),cw[0]),ws(tfp(ci,bw,tm),cw[1])]
 return out,((c[0]*w-b)%p,(c[1]*w-a)%p,c[2]*w%p)
def edge_leaf(rg,kind):
 m=rg.randrange(3)
 if kind=="D":
  b=rg.randrange(2)
  mb=(MH,ML)[b]
  u=rg.getrandbits(64-mb)<<mb|rg.choice((rg.randrange(70),2**mb-1-rg.randrange(70),rg.randrange(2**mb)))
  v=-(-(u<<64|rg.getrandbits(64)if b==0 else rg.getrandbits(64)<<64|u)*p//M128)
  return v*pow(BETA,-m,p)%p or 1
 off=rg.randrange(160<<192)
 w=off if rg.random()<0.5 else p-1-off
 if kind=="Xm":return w*pow(BETA,-m,p)%p or 1
 return w*pow((1-BETA,1-BETA*BETA,BETA-BETA*BETA)[m]%p,-1,p)%p or 1
def deadline(frac):
 now=time.time()
 return now+frac*max(0.0,T_END-now)
def near_tie_x(rg):
 mlt=rg.choice([1,BETA,BETA*BETA,1-BETA,1-BETA*BETA,BETA-BETA*BETA])
 dd=rg.randrange(-(1<<rg.randrange(150,250)),1<<rg.randrange(150,250))
 return dd*pow(mlt%p,-1,p)%p or 1
def index_check(x,rg):
 t,a,b=rg.randrange(3),rg.randrange(1,p),rg.randrange(1,p)
 j=record_fields(leaf_w(represent(x,t,a,b,rg)[0],t))[0]
 return in_W(x),j!=exact_key(x)[1]
def part_a_identities(rng):
 res={}
 res["beta_cube_root_of_1_mod_p"]=pow(BETA,3,p)==1 and BETA!=1
 res["lambda_cube_root_of_1_mod_r"]=pow(LAM,3,r)==1 and LAM!=1
 res["lambda_G_equals_beta_x"]=PROD.mul(LAM,GP)==(BETA*GX%p,GY)
 res["seven_nonresidue_no_x0_point"]=pow(7,(p-1)//2,p)==p-1
 res["window_forms_nonzero"]=all(c%p for c in (1-BETA,1-BETA*BETA,BETA-BETA*BETA,3))
 res["label_constants"]=(FA*FA-FA*FB+FB*FB==p and FB*BETA%p==FA-FB
 and (FA-FB)*BETA%p==p-FA)
 ok=0 # P1 (3.1): M -+ N = 1/(y_Y +- y_X)^2, D = y_Y^2 - y_X^2, w = 1/D^2
 for _ in range(400):
  yx,yy=rng.randrange(1,p),rng.randrange(1,p)
  D=(yy*yy-yx*yx)%p
  if D:
   w=pow(D,-2,p)
   M,N=(pow(D,-1,p)+2*yx*yx*w)%p,2*yx*yy*w%p
   ok+=(M-N)%p==pow(yy+yx,-2,p)and (M+N)%p==pow(yy-yx,-2,p)
 res["p1_identity"]=dict(cases=400,ok=ok==400)
 bad=drops=wrong=done=0
 kinds={"random":0,"near_tie":0}
 dl=deadline(0.10)
 while done<26000 and (done<MIN_RANDOM or time.time()<dl):
  if done<20000:x,kind=rng.randrange(1,p),"random"
  else:x,kind=near_tie_x(rng),"near_tie"
  kinds[kind]+=1
  done+=1
  win,err=index_check(x,rng)
  drops+=win
  wrong+=win and err
  bad+=err and not win
 extra=max(0,MIN_NEAR-kinds["near_tie"])
 rng2=random.Random(SEED+7919)
 for _ in range(extra):
  win,err=index_check(near_tie_x(rng2),rng2)
  kinds["near_tie"]+=1
  drops+=win
  wrong+=win and err
  bad+=err and not win
 res["class_rule"]=dict(checked=kinds,in_window_set=drops,wrong_index_in_window_set=wrong,
 mismatches_outside_window_set=bad,near_ties_from_second_stream=extra)
 res["class_labels"]=class_labels(random.Random(SEED+104729))
 return res
def class_labels(rg):
 out=dict(classes=0,edge_classes=0,members=[0,0,0],in_W=0,in_D=0,wrong_index=0,mismatches=0,
 errors=[[99,-99]for _ in range(6)],tfp_terms={})
 dl=deadline(0.05)
 while out["classes"]<MIN_CLASSES or (out["classes"]<4*MIN_CLASSES and time.time()<dl):
  edge=out["classes"]%3==2
  x0=edge_leaf(rg,("D","Xm","diff")[out["classes"]//3%3])if edge else rg.randrange(1,p)
  v=exact_key(x0)[0]
  out["classes"]+=1
  out["edge_classes"]+=edge
  lab=label(v)
  for m in range(3):
   x=v*pow(BETA,(3-m)%3,p)%p
   t,a,b=rg.randrange(3),rg.randrange(1,p),rg.randrange(1,p)
   w=leaf_w(represent(x,t,a,b,rg,out["tfp_terms"])[0],t)
   for k,c in enumerate((w[0],w[1],ws(w[0],w[1])^M128-1)):
    for i,e in zip((k,k+3),werr(t128(pow(BETA,k,p)*x%p),c)):out["errors"][i]=[min(out["errors"][i][0],e),max(out["errors"][i][1],e)]
   f=record_fields(w)
   if in_W(x):
    out["in_W"]+=1
    continue
   if f[0]!=m:
    out["wrong_index"]+=1
    continue
   out["members"][m]+=1
   if in_D(x):out["in_D"]+=1
   elif cf(f)!=lab:out["mismatches"]+=1
 tm=out["tfp_terms"]
 out["ok"]=(out["mismatches"]==0 and out["wrong_index"]==0 and out["in_D"]>0 and out["in_W"]>0
 and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors"],ERR_H+ERR_L))
 and tm["sum_ok"]and all(lo<=tm[k][0]and tm[k][1]<=hi for k,(lo,hi)in TERMS.items()))
 return out
SENTINEL=0x5E5E5E5E5E5E5E5E
A_SB,A_YR,A_SP=1<<42,(1<<41)+0x10000000,(1<<41)+0x7000000
TB0,MXX=2**64-3*2**57+2**49,1024 # 3.1
W7W=(37,)*4+(36,)*3
W7O=[sum(W7W[:k])for k in range(7)]
PW=[pow(2,o,p)for o in W7O]
def win(v,k):
 return v>>W7O[k]&(1<<W7W[k])-1
class LazyTable(dict):
 def __init__(self,tables,xt=None):
  super().__init__()
  self.tables,self.xt=tables,xt
 def __missing__(self,a):
  if a>=BT0:
   o=a-BT0
   return BETA*(o>>8)*pow(2,43*BT_WIN[o>>5&7],p)%p>>64*(o>>3&3)&M64 if o>>5&7<7 else SENTINEL
  if self.xt and TB0<=a<TB0+(1<<58):
   o=a-TB0
   pr,sl=o>>49,o>>46&7
   if sl<6 and not (sl in (0,3)and pr%2==(sl==3)):
    h,t,k=o>>5&1,o>>9&(1<<37)-1,o>>6&7
    c,add=self.xt[(2*pr-2*h,0)if sl==0 else (2*pr+h,sl)if sl<3 else (2*(pr+h),4)if sl==3 else (2*pr+sl-4,6+h)]
   else:raise KeyError(a)
   if k>6 or k>3 and t>>36: # never built
    raise KeyError(a)
   return (c*t*PW[k]+(0 if k else add))%p>>64*(o>>3&3)&M64
  for lo,hi,first,stride in self.tables:
   if lo<=a<hi:return first+(a-lo)//8*stride
  raise KeyError(a)
BT0,BT_WIN=0xfff8<<48,(0,1,3,2,4,5,3)
def bt_addr(w,i):
 t=w>>43*i&(1<<43)-1
 return BT0+(t<<8)+32*((0,1,3,2,4,5)[i]if i!=3 else 2+4*(w>>128&1))
def blook(w):
 z=sum(BETA*(w>>43*i&(1<<43)-1)*pow(2,43*i,p)%p for i in range(6))
 z=z%2**256+(z>>256)*C0
 return z-2**256+C0 if z>>256 else z
def rows(c):
 return [c*pow(2,64*k,p)%p for k in range(4)]
def xtab(X): # (j, i) -> (c, add) (proof.md 3.1)
 t,A={},[y*y%p for _,y in X]
 for j,(a,y)in enumerate(X):
  t[j,0],t[j,1],t[j,2],t[j,6],t[j,7]=(2*a**3%p,0),(2*a%p,pow(a,4,p)),(a*a%p,0),(2*A[j]%p,0),(2*y%p,0)
  if j%2==0 and j+1<len(X):
   g=(A[j]+A[j+1])%p
   t[j,4]=(-g%p,(A[j]*A[j+1]-7*g)%p)
 return t
def look(t,j,i,v):
 c,add=t[j,i]
 return sum((c*win(v,k)*PW[k]+(0 if k else add))%p for k in range(7))
def yrecord(mem,u,Y): # 3.1 (ygen's in production)
 b,y=Y
 a0,B=A_YR+512*u,y*y%p
 mem[a0]=mem[a0+8]=0
 for o,v in ((16,B),(48,B*B%p),(80,pow(b,4,p))):put(mem,a0+o,v)
 for k,v in enumerate(rows(y)):put(mem,a0+128+32*k,v)
 for k,w in enumerate(corr(b)):put(mem,a0+(112,256)[k],w,2)
 put(mem,a0+464,b)
 for n,v in enumerate((b,pow(b,3,p),3*b*b%p)):
  for k in range(8):mem[a0+272+64*n+8*k]=win(v,k)*512+64*k+(n<<46)if k<7 else 0
def xwords(mem,l,X): # 3.1
 j=l%MXX
 a,y=X[j]
 e,pb=A_SB+512*l,TB0+(j>>1<<49)
 for k,w in enumerate(corr(a)):put(mem,e+320+16*k,w,2)
 put(mem,e+192,y*y%p)
 mem[e+224]=pb
 for h in range(7):mem[e+384+8*h]=pb+(4+(j&1)<<46)+64*h
 mem[e+440]=0 if j%2 else delta_of(X,j)&M64
 if j%2==0:
  mem[e+352]=2*delta_of(X,j)&M64
  mem[e+232],mem[e+240]=TB0+(j>>2<<50)+(2<<46),0
CUBE_E,SQRT_E=(p+2)//9,(p+1)//4
PC,PD,PH,HCAP,LCAP=102,63,150,48,64 # proof.md 2.3: cell, close, head bits; caps
def real_A(A):
 if not 0<A<p or A==7:return None
 y=pow(A,SQRT_E,p)
 if y*y%p!=A:return None
 x=pow(A-7,CUBE_E,p)
 return (x,y)if x**3%p==(A-7)%p else None
def near_A(A0,lo,hi,rg):
 for _ in range(4000):
  P=real_A(A0+lo+rg.randrange(hi-lo))
  if P:return P
def opairs(cands,cw=PC,dw=PD): # reference, by sorting
 l=sorted(cands)
 return [(l[i],l[k])for i in range(len(l))for k in range(i+1,len(l))if l[k][0]-l[i][0]<1<<dw and l[k][0]>>cw==l[i][0]>>cw]
def mine(cands,cw=PC,dw=PD,hs=PH,hb=6,hc=HCAP,lc=LCAP): # 2.3's dictionary: 2^hb heads a home (2^34 in production) on A's
 homes,flag,most={},0,0 # bits from hs, links, <= hc visits an insert, <= lc pairs
 for A,t in cands:
  hd,nx,rec,ls=homes.setdefault(A>>hs+hb,([-1]*(1<<hb),[],[],[]))
  h=A>>hs&(1<<hb)-1;j=hd[h];n=0
  while j>=0:
   if n==hc:flag=1;break
   n+=1;B=rec[j][0]
   if B==A or B>>cw==A>>cw and abs(A-B)<1<<dw and len(ls)==lc:flag=1
   elif B>>cw==A>>cw and abs(A-B)<1<<dw:ls.append(tuple(sorted((rec[j],(A,t)))))
   j=nx[j]
  most=max(most,n);nx.append(hd[h]);hd[h]=len(rec);rec.append((A,t))
 return homes,flag,most
CMPS=[0,0]
def msort(v): # the collector's merge sort by lower A; comparisons against n ceil(log2 n)
 w,n=1,len(v)
 CMPS[1]+=n*max(n-1,0).bit_length()
 while w<n:
  o=[]
  for i in range(0,n,2*w):
   a,b=v[i:i+w],v[i+w:i+2*w]
   while a and b:
    CMPS[0]+=1
    o.append((b if b[0][0][0]<a[0][0][0]else a).pop(0))
   o+=a+b
  v,w=o,2*w
 return v
def psel(cands,need,pc=1<<20,**k): # 2.3's rule on mine's lists; None = abort, no padding
 homes,flag,_=mine(cands,**k)
 pr=[x for v in homes.values()for x in v[3]]
 used=[t for a,b in pr for t in (a[1],b[1])]
 if flag or len(pr)>pc or len(set(used))<len(used)or len(pr)<need:return None
 return msort(pr)[:need]
def phase0_gen(rg,n=80,cw=250,dw=246): # real [s]P, toy widths: 4 homes of 4 heads of 4 cells
 c=[(P[1]*P[1]%p,i)for i in range(n)for P in [PROD.mul(rg.randrange(1,r),GP)]]
 k=dict(cw=cw,dw=dw,hs=252,hb=2,hc=n,lc=n)
 pr=opairs(c,cw,dw);f=len(pr)
 o=c[:];rg.shuffle(o)
 homes,flag,most=mine(c,**k)
 s1=psel(c,f,**k)
 ts=[t for a,b in pr for t in (a[1],b[1])]
 sh=len(set(ts))<len(ts)
 hh=[A>>252 for A,_ in c];full=max(map(hh.count,hh));lm=max(len(v[3])for v in homes.values())
 res=dict(candidates=n,found=f,shared_candidate=int(sh),
 lists_equal_reference=not flag and sorted(x for v in homes.values()for x in v[3])==pr,
 heads_with_several_cells=sum(len({A>>cw for A,_ in c if A>>252==h})>1 for h in set(hh)),fullest_head=full,
 visits_at_most_fullest_less_one=most==full-1,head_cap_aborts=mine(c,**dict(k,hc=full-2))[1]==1 and mine(c,**dict(k,hc=full-1))[1]==0,
 list_cap_aborts=lm==0 or mine(c,**dict(k,lc=lm-1))[1]==1,collector_cap_aborts=sh or f==0 or psel(c,f,pc=f-1,**k)is None,
 sort_is_reference_order=sh or s1==pr,sort_comparisons_within_n_log_n=CMPS[0]<=CMPS[1],
 arrival_order_invariant=s1==psel(o,f,**k),abort_if_one_more_needed=psel(c,f+1,**k)is None,selected_iff_no_shared=(s1 is None)==sh)
 res["ok"]=all(v for q,v in res.items()if q not in ("found","shared_candidate"))
 return res
def phase0_toy(rg): # production widths, 64 heads a home
 cd={}
 def add(P,k):
  t=(k,len(cd));cd[t]=P[1]*P[1]%p;return t
 while True: # a pair across a 2^63 boundary inside one cell
  A0=rg.randrange(1,p>>PC)<<PC|rg.randrange(1<<PC)
  B=((A0>>PD)+1)<<PD
  Pl,Pu=near_A(B,-(1<<40),0,rg),near_A(B,0,1<<40,rg)
  if Pl and Pu and B>>PC==A0>>PC:
   st=(add(Pl,'S'),add(Pu,'S'))
   break
 for _ in range(5):
  P=near_A(rg.randrange(1,p>>PC)<<PC|rg.randrange(1<<PC)&~((1<<PD)-1),0,1<<40,rg)
  Q_=P and near_A(P[1]*P[1]%p,1,1<<PD,rg)
  if Q_:add(P,'R');add(Q_,'R')
 cl=[(A,t)for t,A in cd.items()]
 n=len(opairs(cl));sel=psel(cl,n)
 sh=cl[:];rg.shuffle(sh)
 lo,hi=cd[st[0]],cd[st[1]]
 Pw=near_A(hi,1,1<<40,rg)
 H0=rg.randrange(1,p>>PH)<<PH
 H=[(P[1]**2%p,('H',i))for i in range(HCAP+2)for P in [near_A(H0+(i<<120),0,1<<40,rg)]]
 Hc=(near_A(H[0][0],1,1<<62,rg)[1]**2%p,('H',99))
 res=dict(pairs=n,straddle_crosses_bucket=lo>>PD!=hi>>PD,
 straddle_selected=sel is not None and any({a[1],b[1]}==set(st)for a,b in sel),
 arrival_order_invariant=sel==psel(sh,n),abort_short=psel(cl,n+1)is None,
 abort_shared_candidate=psel(cl+[(Pw[1]*Pw[1]%p,('W',0))],1)is None,abort_equal_A=psel(cl+[(lo,('E',0))],1)is None,
 one_head_50_cells=len({A>>PH for A,_ in H})==1 and len({A>>PC for A,_ in H})==len(H),
 head_of_49_passes=mine(cl+H[:-1])[1:]==(0,HCAP)and psel(cl+H[:-1],n)==sel,
 head_of_50_aborts=mine(cl+H)[1]==1 and psel(cl+H,n)is None,
 close_pair_found_in_full_head=(H[0],Hc)in(psel(cl+H[:HCAP]+[Hc],n+1)or()),
 selected_real_oriented=sel is not None and all(real_A(a[0])and real_A(b[0])and 0<b[0]-a[0]<2**63 and a[0]>>PC==b[0]>>PC for a,b in sel))
 res["ok"]=all(v for k,v in res.items()if k!="pairs")
 return res
def gf128(a,b): # GF(2)[x]/(x^7 + x + 1)
 z=0
 for i in range(7):
  if b>>i&1:z^=a<<i
 for i in range(13,6,-1):
  if z>>i&1:z^=0b10000011<<(i-7)
 return z
def class_pool(c): # class c = (u, w): 4 lines b = u a + v, v in 4w..4w+3, as ranks 128 a + b
 u,w=c%128,c//128
 return [128*a+(gf128(u,a)^v)for v in range(4*w,4*w+4)for a in range(128)]
def design_check(rg):
 got={}
 for _ in range(12):
  c,d=rg.randrange(4096),rg.randrange(4096)
  got.setdefault("same" if c==d else "slope" if c%128==d%128 else "diff",set()).add(len(set(class_pool(c))&set(class_pool(d))))
 c=rg.randrange(4096)
 got.setdefault("slope",set()).add(len(set(class_pool(c))&set(class_pool((c+128)%4096))))
 got.setdefault("same",set()).add(len(set(class_pool(c))))
 return dict(shared_pairs=str({k:sorted(v)for k,v in got.items()}),ok=got.get("diff")=={16}and got["slope"]=={0}and got["same"]=={512})
def resolver_swap(rg): # X public, Y = [b_g]Q + [c]P: x = (u mu - mu')/(nu' - u nu)
 ok=0
 for _ in range(4):
  x=rg.randrange(1,r);Q=PROD.mul(x,GP)
  sj,sj2,bg,bg2,c=(rg.randrange(1,r)for _ in range(5))
  t,t2,e,e2=rg.randrange(3),rg.randrange(3),rg.choice((1,-1)),rg.choice((1,-1))
  u=rg.choice((1,-1))*pow(LAM,rg.randrange(3),r)%r
  mu,nu=(pow(LAM,t,r)*sj+e*c)%r,e*bg%r
  c2=(e2*(u*(mu+nu*x)-pow(LAM,t2,r)*sj2)-bg2*x)%r
  mu2,nu2=(pow(LAM,t2,r)*sj2+e2*c2)%r,e2*bg2%r
  L,L2=(PROD.add(PROD.mul(m,GP),PROD.mul(n_,Q))for m,n_ in ((mu,nu),(mu2,nu2)))
  xr=(u*mu-mu2)*pow(nu2-u*nu,-1,r)%r
  ok+=L2==PROD.mul(u,L)and xr==x and PROD.mul(xr,GP)==Q
 x=rg.randrange(1,r);Q=PROD.mul(x,GP);bg,sj,sj2,c=(rg.randrange(1,r)for _ in range(4));e,e2=1,-1;u=e*e2%r
 mu,nu=(sj+e*c)%r,e*bg%r
 c2=(e2*(u*(mu+nu*x)-sj2)-bg*x)%r
 mu2,nu2=(sj2+e2*c2)%r,e2*bg%r
 same=PROD.add(PROD.mul(mu2,GP),PROD.mul(nu2,Q))==PROD.mul(u,PROD.add(PROD.mul(mu,GP),PROD.mul(nu,Q)))and (nu2-u*nu)%r==0
 return dict(cases=4,recovered=ok,same_node_collision_useless=same,ok=ok==4 and same)
def delta_of(X,j): # A_(j+1) - A_j (A = y^2 mod p); in (0, 2^63) on mined pairs
 return (X[j+1][1]**2%p-X[j][1]**2%p)%p
def close_point(A0,rg,dmax=2**63):
 while True:
  d=rg.choice((1+rg.randrange(dmax-1),dmax-1-rg.randrange(2**20),1+rg.randrange(2**20)))
  P=real_A(A0+d)
  if P:return (P[0]*BT[rg.randrange(3)]%p,P[1]if rg.randrange(2)else p-P[1])
PSEL=[]
def mined_xs(rg,n,first): # psel's output on a shuffled pool of close pairs, decoys, a cell-split pair
 out,P=[],first
 while len(out)<n:
  if len(out)//2%5==3:
   B=rg.randrange(1,p>>63)<<63
   a,b=near_A(B,-(1<<20),0,rg),near_A(B,0,1<<20,rg)
   if a and b and B%2**PC:out+=[a,b];continue
  out+=[P,close_point(P[1]**2%p,rg)]
  P=PROD.add(P,GP)
 B=rg.randrange(1,p>>PC)<<PC
 dec=[PROD.mul(rg.randrange(1,r),GP)for _ in range(8)]+[near_A(B,-(1<<40),0,rg),near_A(B,0,1<<40,rg)]
 pool=[(Z[1]**2%p,Z)for Z in out+dec];rg.shuffle(pool)
 X=[Z for a,b in psel(pool,n//2)for Z in (a[1],b[1])]
 PSEL.append(sorted(X)==sorted(out)and psel(pool,n//2+1)is None)
 return X
def components(X,Y):
 a,b=X[0],Y[0]
 E,F,G=(b**4+2*a**3*b)%p,(a**4+2*a*b**3)%p,3*a*a*b*b%p
 out=[]
 for w in (pow(Y[1]+X[1],-2,p),pow(Y[1]-X[1],-2,p)):out+=[E*w%p,F*w%p,G*w%p]
 return out
def lfold(z): # LF: z < 2^257 folded by C0, the carry adds C0 once more
 z=z%2**256+(z>>256)*C0
 return z-2**256+C0 if z>>256 else z
def fsub(a,b): # FS, bit for bit (wrong mod p iff b > a + p: 3.4)
 d=a-b
 return d if d>=0 else (d+2**256-C0)%2**256
def fadd(a,b): # FADDW, bit for bit (wrong mod p iff a + b >= 2^257 - C0: 3.4)
 z=a+b
 return z if z<2**256 else (z-2**256+C0)%2**256
def prod4(m,RR):
 return sum((m>>64*k&M64)*RR[k]for k in range(4))
def fwd_model(X,Ys,t,q2,rho):
 acc,out=rho,[]
 for l in range(0,q2,2):
  b,y=Ys[l//MXX]
  e=lfold(pow(y,4,p)+look(t,l%MXX,4,pow(b,3,p)))
  out.append((e,acc))
  acc=acc*e%p
 return out,acc
def bwd_model(X,Ys,t,q2,ea,R,tm=None):
 """bwd's twelve stored words per incidence, bit for bit (proof.md 3.1, 3.3, 3.4)."""
 out={}
 for pr in range(q2//2-1,-1,-1):
  l1=2*pr+1
  j1=l1%MXX
  b,yy=Ys[l1//MXX]
  RR=shx_rows(R)
  ie,R=fold(prod4(ea[pr][1],RR)),fold(prod4(ea[pr][0],RR))
  r1=fold(prod4(ie,shx_rows(fsub(yy*yy%p,X[j1][1]**2%p))))
  for l,inv in ((l1,lfold(r1+delta_of(X,j1-1)*ie)),(l1-1,r1)): # DT: 1/D_j + delta inv_e, folded as LF
   j=l%MXX
   a,yx=X[j]
   ww=fold(inv*inv%2**256+(inv*inv>>256)*C0)
   Mv,Lv=lfold(inv+look(t,j,6,ww)),lfold(look(t,j,7,ww))
   Nv=fold(prod4(Lv,rows(yy)))
   w=[fsub(Mv,Nv),fadd(Mv,Nv)]
   bw=[blook(v)for v in w]
   ee=lfold(pow(b,4,p)+look(t,j-j%2,0,b)) # EQUAD
   cs=(lfold(ee+2*delta_of(X,j-1)*b)if j%2 else ee,lfold(look(t,j,1,pow(b,3,p))),lfold(look(t,j,2,3*b*b%p)))
   fr={}
   for i,(c,y)in enumerate(zip(cs,(b,a,None))):
    cw=corr(y)if y else (0,0)
    for sg in (0,1):
     fr[sg,2*i]=ws(tfp(c,w[sg],tm),cw[0])
     fr[sg,2*i+1]=ws(tfp(c,bw[sg],tm),cw[1])
   out[l]=[fr[sg,i]for sg in (0,1)for i in range(6)]
 return out,R
def six_leaves(X,Y):
 out=[]
 for t in range(3):
  PX=(pow(BETA,t,p)*X[0]%p,X[1])
  out+=[PROD.add(PX,Y)[0],PROD.add(PX,PROD.neg(Y))[0]]
 return out
def fixup_extra(prog):
 out,last={},0
 for i,(mn,_)in enumerate(prog):
  if mn=="stadd":out[i]=i-last
  elif mn in ("b.lo","b.hs","cbnz"):last=i
 return out
def flags_dead_after_fixups(prog,labels):
 reads=("adc","adcs","sbc","sbcs","csel","cinc","cset","ccmp","ccmn","b.eq","b.ne","b.lo","b.hs")
 def target(t,pc):
  if t[:-1].isdigit():
   c_=labels[t[:-1]]
   return [q for q in c_ if q<pc][-1]if t[-1]=="b" else [q for q in c_ if q>=pc][0]
  return labels[t][0]
 for i in fixup_extra(prog):
  todo,seen=[i+1],set()
  while todo:
   j=todo.pop()
   if j in seen:continue
   seen.add(j)
   mn,o=prog[j]
   if mn in reads:return False
   if mn in ("adds","subs","cmp","cmn","tst","ret"):continue
   if mn in ("b","tbz","tbnz","cbz","cbnz"):
    todo.append(target(o[-1],j+1))
    if mn=="b":continue
   todo.append(j+1)
 return True
LS_F,LS_B=([(TB0,TB0+(1<<58))],[(A_SB,A_SB+(1<<34))]),([(A_SB,A_SB+(1<<34)),(A_YR,A_YR+(1<<24)),(TB0,TB0+(1<<58)),(BT0,1<<64)],[(A_SB,A_SB+(1<<34))])
def batch_mem(X,Ys,q2,tables):
 t,L=xtab(X[:min(q2,MXX)]),(q2+MXX-1)//MXX
 mem=LazyTable(tables,t)
 for u in range(L):yrecord(mem,u,Ys[u])
 for l in range(q2):xwords(mem,l,X)
 return mem,t,L
def flip_batch(prog,labels,X,Ys,q2,cert,W2=0x123456789,override=None,rho=None):
 """One batch on fwd, bwd, keys: X_(l mod 1024) against Y_(l >> 10), l < q2 (proof.md 3.1)."""
 fx=fixup_extra(prog)
 vis={k:dict.fromkeys(fx,0)for k in ("fwd","bwd")}
 BASE,END=1<<62,(1<<62)+(1<<40)
 mem,t,L=batch_mem(X,Ys,q2,[(0,8<<36,machine.PB1+32,machine.BLK),(machine.RC,machine.RC+(8<<36),BASE,0)])
 rho=rho or random.Random(q2).randrange(1,p)
 regs=dict(x0=q2,x1=A_YR,x3=A_SB+512,sp=A_SP,x30=SENTINEL)
 regs.update(as_regs(rho,("x24","x25","x26","x27")))
 out={}
 nf=a64.a64_run(prog,labels,"fwd",regs,mem,out_regs=out,visits=vis["fwd"],ls64=LS_F)
 words,acc,sp_kept=[out["x0"]],reg_int(out,("x24","x25","x26","x27")),out["sp"]==A_SP
 ea=[(get(mem,A_SB+1024*i+512),get(mem,A_SB+1024*i+544))for i in range(q2//2)]
 mf,macc=fwd_model(X,Ys,t,q2,rho)
 R0=pow(acc,p-2,p)
 regs=dict(x3=A_SB+ENTRY*(q2-1),sp=A_YR+512*(L-1),x30=SENTINEL)
 regs.update(as_regs(R0,("x25","x26","x27","x28")))
 out={}
 nb=a64.a64_run(prog,labels,"bwd",regs,mem,out_regs=out,visits=vis["bwd"],ls64=LS_B)
 words.append(out["x0"])
 sp_kept&=out["sp"]==A_YR
 final_R=reg_int(out,("x25","x26","x27","x28"))
 tm={}
 model,mR=bwd_model(X,Ys,t,q2,ea,R0,tm)
 bad_model,bad_bound,bad_leaf,err,leaves=int(mR!=final_R),0,0,[[99,-99]for _ in range(4)],[]
 for l in range(q2):
  Xj,Y,e=X[l%MXX],Ys[l//MXX],A_SB+ENTRY*l
  cs=components(Xj,Y)
  bad_model+=sum(get(mem,e+o,2)!=w for o,w in zip(FRAC_OFF[0]+FRAC_OFF[1],model[l]))
  for i,ds in enumerate(frac_errs(mem,e,Xj,Y)):
   g=i%6//4
   for k,d in enumerate(ds):
    err[2*g+k]=[min(err[2*g+k][0],d),max(err[2*g+k][1],d)]
    bad_bound+=not TERMS["hl"[k]][0]+g-1<=d<=TERMS["hl"[k]][1]
  six=six_leaves(Xj,Y)
  bad_leaf+=sum(dft_leaf(cs[3*sg:3*sg+3],tt,Xj[0],Y[0])!=six[2*tt+sg]for sg in (0,1)for tt in range(3))
  leaves.append(six)
 for (l,sg),(x,tt,rg)in (override or {}).items():
  f6,c3=represent(x,tt,X[l%MXX][0],Ys[l//MXX][0],rg)
  for o,v in zip(FRAC_OFF[sg],f6):put(mem,A_SB+ENTRY*l+o,v,2)
  for t_ in range(3):leaves[l][2*t_+sg]=(c3[0]+BT[t_]*c3[1]+BT[2*t_%3]*c3[2])%p
 port=machine.Port(mem,machine.MB1,machine.PB1)
 regs=dict(x3=A_SB+ENTRY*(q2-1)+KEYS0,x12=END,x17=W2,x18=0,sp=A_SP,x30=SENTINEL)
 out={}
 nk=a64.a64_run(prog,labels,"keys",regs,mem,out_regs=out,dev=port.dev,ls64=([(A_SB,A_SB+(1<<34))],[]))
 sp_kept&=out["sp"]==A_SP
 snd=machine.Sender(lambda d:BASE,END,False,cert)
 meter,fields=machine.keys_reference(sys.modules[__name__],mem,A_SB,q2,W2,snd,cert)
 flat=[leaves[l][2*tt+sg]for l in range(q2-1,-1,-1)for sg in (0,1)for tt in range(3)]
 win=[in_W(x)for x in flat]
 dset=[in_D(x)for x in flat]
 wrong=sum(not w and f[0]!=exact_key(x)[1]for w,f,x in zip(win,fields,flat))
 canon=sum(cf(f)!=class_label(x)for w,d,f,x in zip(win,dset,fields,flat)if not w and not d)
 nfix={k:[sum(v for i,v in vis[k].items()if fx[i]==e)for e in (3,4,6)]for k in vis}
 fok=all(any(a<=(x-A)%512<a+b for a,b in F)for A,F in ((A_SB,cert.ENTRY_F),(A_YR,cert.YREC_F))for x in mem if A<=x<A+(1<<34))
 ext={k:3*nfix[k][0]+4*nfix[k][1]+6*nfix[k][2]for k in vis}
 return dict(fwd=nf,bwd=nb,keys=nk,acc_ok=acc%p==macc%p,fwd_stores_ok=ea==mf,fraction_mismatches=bad_model,
 fractions_outside_bounds=bad_bound,fraction_errors=err,tfp_terms=tm,identity_mismatches=bad_leaf,
 fixups_3_4=nfix,fixup_extra=ext,
 fixup_words_ok=words==[C0*sum(nfix[k])for k in ("fwd","bwd")],
 fwd_bwd_formulas=nf==cert.fwd_count(q2)+ext["fwd"]and nb==cert.bwd_count(q2)+ext["bwd"],
 records_ok=port.blocks==snd.msgs and not port.msgs and machine.state_ok(mem,snd,0,machine.PB1),
 records=len(snd.seq),key_paths=meter,in_W=sum(win),
 in_D=sum(d and not w for w,d in zip(win,dset)),wrong_index=wrong,keys_not_canonical=canon,
 keys_formula=nk==cert.KEY_FIX+cert.KEY_INC*q2+meter+snd.extra,
 final_R_is_inverse_rho=final_R*rho%p==1,sp_kept=sp_kept,fields_ok=fok)
def reg_int(regs,names):
 return sum(regs[n]<<(64*k)for k,n in enumerate(names))
def as_regs(v,names):
 return {n:w for n,w in zip(names,limbs(v,len(names)))}
def macro_counts(src,entry,names):
 macros,cur,lines={},None,[ln.split("//")[0].strip()for ln in src.split("\n")]
 for s in lines:
  if s.startswith(".macro"):
   cur=s.split()[1].rstrip(",")
   macros[cur]=[]
  elif s.startswith(".endm"):cur=None
  elif cur is not None:macros[cur].append(s)
 out=dict.fromkeys(names,0)
 def walk(seq):
  for s in seq:
   w=s.split()[0]if s else ""
   if w in out:out[w]+=1
   if w in macros:walk(macros[w])
 start=lines.index(entry+":")
 walk(lines[start+1:lines.index("ret",start)])
 return out
def sqrt_mod_p(v):
 s_=pow(v,(p+1)//4,p)
 return s_ if s_*s_%p==v%p else None
def degenerate_batch(prog,labels,X,Ys,q2):
 mem=batch_mem(X,Ys,q2,[])[0]
 regs=dict(x0=q2,x1=A_YR,x3=A_SB+512,sp=A_SP,x30=SENTINEL)
 regs.update(as_regs(1,("x24","x25","x26","x27")))
 out={}
 a64.a64_run(prog,labels,"fwd",regs,mem,out_regs=out,ls64=LS_F)
 acc=reg_int(out,("x24","x25","x26","x27"))
 regs=dict(x3=A_SB+ENTRY*(q2-1),sp=A_YR,x30=SENTINEL)
 regs.update(as_regs(pow(acc,p-2,p),("x25","x26","x27","x28")))
 a64.a64_run(prog,labels,"bwd",regs,mem,ls64=LS_B)
 good=0
 for l in range(q2):
  Xj,Y=X[l],Ys[0]
  if (Y[1]*Y[1]-Xj[1]*Xj[1])%p:good+=sum(TERMS["h"][0]-1<=h<=TERMS["h"][1]and TERMS["l"][0]-1<=l_<=TERMS["l"][1]for h,l_ in frac_errs(mem,A_SB+ENTRY*l,Xj,Y))
 return acc,good
def edge_batch(prog,labels,X,Ys,q2,cert,rg,kinds,need):
 over={}
 for l in range(q2):
  for sg in (0,1):
   while True:
    x=edge_leaf(rg,kinds[(2*l+sg)%len(kinds)])
    if len(over)>=4 or (in_W(x)if need=="W" else in_D(x)and not in_W(x)):break
   over[(l,sg)]=(x,rg.randrange(3),rg)
 return flip_batch(prog,labels,X,Ys,q2,cert,override=over)
def part_a_kernels(rng,cert):
 res={}
 src=open(os.path.join(HERE,"kernels","leaf.s")).read()
 prog,labels=a64.a64_load(src)
 X=mined_xs(rng,MXX,PROD.mul(rng.randrange(1,r),GP))
 res["x_pairs_one_word"]=all(0<delta_of(X,j)<2**63 for j in range(0,MXX,2))
 Ys=[PROD.mul(rng.randrange(1,r),GP)for _ in range(2)]
 runs={q2:flip_batch(prog,labels,X,Ys,q2,cert)for q2 in (32,64,MXX+32)} # MXX + 32 crosses the Y switch at MXX
 res["runs"]=runs
 net={q:{k:d[k]-d["fixup_extra"][k]for k in ("fwd","bwd")}for q,d in runs.items()}
 res["fwd_bwd_formulas_hold"]=all(d["fwd_bwd_formulas"]and d["fixup_words_ok"]for d in runs.values())
 res["counts_ok"]=all(net[q]["fwd"]==cert.fwd_count(q)and net[q]["bwd"]==cert.bwd_count(q)for q in runs)
 res["per_incidence"]={k:float(Fraction(net[MXX+32][k]-net[32][k],MXX))for k in ("fwd","bwd")}
 res["fixups"]=dict(sites=len(fixup_extra(prog)),flags_dead_at_every_join=flags_dead_after_fixups(prog,labels))
 tbzs=[o for mn,o in prog if mn in ("tbz","tbnz")]
 res["loop_ends_on_bit_34"]=[o for o in tbzs if o[0]=="x3"]==[["x3","#34","2f"],["x3","#34","1b"]] \
  and all(o[1]=="#13" for o in tbzs if o[0]=="x16")
 res["keys_formula"]="KEY_FIX + KEY_INC q2 + JPATH, SJOIN + KFLUSH or KEY_BLOCK_DROP per block"
 res["keys_formula_holds"]=all(d["keys_formula"]for d in runs.values())
 res["tfp_terms"]=dict(ranges={k:[min(d["tfp_terms"][k][0]for d in runs.values()),max(d["tfp_terms"][k][1]for d in runs.values())]
 for k in TERMS},bounds=TERMS,
 sums_exact=all(d["tfp_terms"]["sum_ok"]for d in runs.values()))
 res["fraction_errors"]=[[min(d["fraction_errors"][g][0]for d in runs.values()),max(d["fraction_errors"][g][1]
 for d in runs.values())]for g in range(4)]
 res["all_correct"]=res["x_pairs_one_word"]and all(d["acc_ok"]and d["fwd_stores_ok"]and d["fraction_mismatches"]==0
 and d["fractions_outside_bounds"]==0 and d["identity_mismatches"]==0 and d["records_ok"]
 and d["keys_not_canonical"]==0 and d["wrong_index"]==0 and d["sp_kept"]and d["final_R_is_inverse_rho"]and d["fields_ok"]
 for d in runs.values())
 erng=random.Random(SEED+49979687)
 for name,kinds,need in (("window_batch",("Xm","diff"),"W"),("boundary_batch",("D",),"D")):
  w=edge_batch(prog,labels,X,Ys,32,cert,erng,kinds,need)
  res[name]=dict(leaves=192,in_W=w["in_W"],in_D_outside_W=w["in_D"],wrong_index=w["wrong_index"],
 keys_not_canonical=w["keys_not_canonical"],records=w["records"],records_ok=w["records_ok"],
 key_paths=w["key_paths"],keys=w["keys"],keys_formula=w["keys_formula"],
 ok=(w[{"W":"in_W","D":"in_D"}[need]]>=4 and w["wrong_index"]==0
 and w["keys_not_canonical"]==0 and w["records"]==192 and w["records_ok"]
 and w["keys_formula"]and w["fwd_bwd_formulas"]and w["fixup_words_ok"]and w["fields_ok"]))
 gp,gl=a64.a64_load(src.replace("\nb.lo 9f\n","\ncbnz xzr,9f\n").replace("\nb.hs 9f\n","\ncbnz xzr,9f\n"))
 g=flip_batch(gp,gl,X,Ys,32,cert)
 per_pair=sum(cert.FIX_SITES.values())
 ext=sum(cert.FIX_EXTRA[m]*cert.FIX_SITES[m]for m in cert.FIX_SITES)
 res["forced_fixups"]=dict(incidences=32,fixups_3_4=g["fixups_3_4"],fixup_extra=g["fixup_extra"],
 fraction_mismatches=g["fraction_mismatches"],formulas=g["fwd_bwd_formulas"],
 ok=g["fraction_mismatches"]==0 and g["fractions_outside_bounds"]==0 and g["fwd_stores_ok"]
 and g["identity_mismatches"]==0 and g["fwd_bwd_formulas"]and g["fixup_words_ok"]and g["records_ok"]and g["fields_ok"]
 and sum(map(sum,g["fixups_3_4"].values()))==16*per_pair and sum(g["fixup_extra"].values())==16*ext)
 fprog,flabels=a64.a64_load(open(os.path.join(HERE,"kernels","secp_fe.s")).read())
 cnts,ok=set(),True
 for _ in range(3):
  a=rng.randrange(1,p)
  mem={}
  put(mem,0x1000,a)
  cnts.add(a64.a64_run(fprog,flabels,"fe_inv",dict(x0=0x2000,x1=0x1000,sp=0x100000,x30=SENTINEL),mem,
 limit=10**6))
  ok&=get(mem,0x2000)%p==pow(a,p-2,p)
 res["fe_inv"]=dict(counts=sorted(cnts),correct=ok)
 C_,T_,S_=(PROD.mul(rng.randrange(1,r),GP)for _ in range(3))
 yg=dict(runs=[],words_bad=0,counts_ok=True,clean_flags_zero=True)
 for M,N,CC in ((2,3,C_),(4,2,C_),(1,2,S_)):
  mem={}
  for a_,P in ((0x1000,CC),(0x1040,T_),(0x1080,S_)):
   put(mem,a_,P[0])
   put(mem,a_+32,P[1])
  out={}
  n=a64.a64_run(fprog,flabels,"ygen",dict(x0=A_YR,x1=M,x2=N,x3=0x1000,x4=0x1040,x5=0x1080,x6=0x100000,
 sp=0x8000000,x30=SENTINEL),mem,limit=10**9,out_regs=out)
  yg["runs"].append(dict(m=M,n=N,count=n,flag=out["x0"]))
  yg["counts_ok"]&=n==cert.ygen_count(M,N)
  if CC is S_:
   yg["defect_flagged"]=out["x0"]!=0
   continue
  yg["clean_flags_zero"]&=out["x0"]==0
  ref={}
  for nn in range(N):
   Ym=CC
   for m in range(M):
    Y=Ym
    for _ in range(nn):Y=PROD.add(Y,S_)
    yrecord(ref,nn*M+m,Y)
    Ym=PROD.add(Ym,T_)
  yg["words_bad"]+=sum(mem.get(a_)!=v for a_,v in ref.items())
 yg["ok"]=yg["words_bad"]==0 and yg["counts_ok"]and yg["clean_flags_zero"]and yg["defect_flagged"]
 res["ygen"]=yg
 names=("SHX","FDX","F5","FS","FADDW","BLK6","LF","DF")
 counts={k:macro_counts(src,k,names+("PFWD2","BWDPA","BWDPB"))for k in ("fwd","bwd","keys")}
 copies={"fwd":2*counts["fwd"]["PFWD2"],"bwd":counts["bwd"]["BWDPA"]+counts["bwd"]["BWDPB"]}
 per={n:sum(counts[k][n]//copies[k]for k in copies)for n in names}
 per["BLK"]=per.pop("BLK6")+per.pop("LF")+per.pop("DF")
 res["site_counts"]=dict(counts,per_pair=per,pair_copies_in_listing=copies)
 res["sites_match_certificate"]=(all(counts[k][n]%copies[k]==0 for k in copies for n in names)
 and {n:per[n]for n in cert.FIX_SITES}==cert.FIX_SITES and sum(counts["keys"].values())==0
 and cert.FIX_LF==sum(range(9))*C0<<192 and cert.FIX_BLK==2*sum(range(6))*C0<<192 and cert.FIX_EXTRA["BLK"]==sites.EXTRA["t_lf"]==sites.EXTRA["t_blk"])
 srng=random.Random(SEED+15485863)
 res["site_checks"]=sites.site_checks(sys.modules[__name__],srng,fprog,flabels)
 Yd=PROD.mul(srng.randrange(1,r),GP)
 Xd=mined_xs(srng,32,PROD.mul(srng.randrange(1,r),GP))
 Xd[5]=(BETA*Yd[0]%p,p-Yd[1]) # Y = -phi^2(X_5): x(Y) = beta^2 x(X_5), so D = 0
 Xd[4]=near_A(Yd[1]**2%p,-(1<<40),0,srng)
 acc0,good=degenerate_batch(prog,labels,Xd,[Yd],32)
 res["d_zero_batch"]=dict(acc_zero=acc0%p==0,fractions_within_bounds=good,all_other_fractions_wrong=good==0)
 res["kernel_xs_from_psel"]=dict(calls=len(PSEL),ok=len(PSEL)==2 and all(PSEL))
 res["final_R_inverse_rho_and_unused"]=all(d["final_R_is_inverse_rho"]and d["records_ok"]for d in runs.values())
 prog_r,lab_r=a64.a64_load(open(os.path.join(HERE,"kernels","route.s")).read())
 res["program_instructions"]=dict(leaf=len(prog),route=len(prog_r),secp_fe=len(fprog),
 total=len(prog)+len(prog_r)+len(fprog))
 return res
def lattice_point(tv,tu):
 c1=(2*(FA*tv+(FA-FB)*tu)+p)//(2*p)
 c2=(2*((FA-FB)*tv-FB*tu)+p)//(2*p)
 return c1*FB+c2*(FA-FB),c1*(FA-FB)-c2*FA
def planted_x(rg,d1,d2,v=None,j=None):
 while True:
  vc=v
  if vc is None:
   c1,c2=rg.getrandbits(64)<<36|d1,rg.getrandbits(64)<<36|d2
   vv,uu=lattice_point((c1<<28|0x8200000)*p>>128,(c2<<28|0x8200000)*p>>128)
   if not (0<vv<p and uu==BETA*vv%p and comp(vv)==c1 and comp(uu)==c2 and exact_key(vv)[0]==vv):continue
   vc=vv
  m=rg.randrange(3)if j is None else j
  x=vc*pow(BETA,(3-m)%3,p)%p
  if exact_key(x)[1]==m and not in_W(x)and not in_D(x):return x,vc
def load_certificate():
 buf=io.StringIO()
 with contextlib.redirect_stdout(buf):import certificate as cert
 lines=[" ".join(ln.split())for ln in buf.getvalue().splitlines()]
 return cert,[ln for ln in lines if " ok" in ln or "FAIL" in ln or ln.startswith(("q1 ","E_total","declared:"))]
def part_d(cert):
 return dict(passed=not cert.FAILED,failed=cert.FAILED,checks=CERT_LINES,claim="%d.%015d"%divmod(int(cert.CLAIM_E*10**15),10**15),q1=cert.Q1GD,
 per_leaf_declared=format(float(cert.DECLARED["per_leaf"]),".4f"),
 mu=float(cert.TFULL["mu"]),mu_good=float(cert.TFULL["mu_g"]),success_lower_bound="0.39",success_floor="0.%025d"%int(cert.P_GATE*cert.TFULL["lo"]*10**25),run_success=float(cert.TFULL["lo"]),per_leaf=float(cert.PER_LEAF_E),
 groups=cert.G_GROUPS,phi_16=float(cert.PHI[16]),
 d_batch=float(cert.D_BATCH),g_def=float(cert.TFULL["g_def"]),d_fixup_meter=float(cert.TFULL["d_f"]),
 d_label_error=float(cert.D_LAB),gate_K=cert.GATE_K)
def canon_key(G,D,WT,p=32779,K=12): # success.md 6.4 at small scale
 b=next(pow(g,(p-1)//3,p)for g in range(2,99)if pow(g,(p-1)//3,p)!=1)
 M,(lo,hi),n,bad,gm=1<<K,ERR_H[0],0,0,0
 def key(f0,f1):
  s=(f0+f1)%M
  return (2 if f1+s>=M else 1)if f0>=f1 else (2 if f0+s>=M else 0)
 F=((1,-1,0),(1,1,M),(1,2,M),(1,2,2*M),(2,1,M),(2,1,2*M))
 for x in range(1,p):
  v=(x,b*x%p,b*b*x%p)
  T=[(u<<K)//p for u in v]
  if min(T)<WT or max(T)>=M-WT:continue
  n+=1
  a0,z0,a1,z1=T[0]-hi,T[0]-lo,T[1]-hi,T[1]-lo
  if all(len({u*f0+w*f1>=c for f0 in (a0,z0)for f1 in (a1,z1)})==1 for u,w,c in F):js={key(a0,a1)}
  else:js={key(f0,f1)for f0 in range(a0,z0+1)for f1 in {a1,z1}|{q for q in (f0,f0+1,M-f0-1,M-f0)if a1<=q<=z1}}
  for j in js:
   i,l,g=T[j],T[(j+1)%3],M*(v[j]-min(v))
   gm=max(gm,g)
   bad+=g>=G*p or not (l>=i-D and (2*i+l<=M+D or (i+l>=M-D and 2*i+l<=2*M+D)))
 return dict(p=p,K=K,x_checked=n,violations=bad,max_gap=gm/p,ok=bad==0 and n>p//2)
def main():
 global CERT_LINES,ERR_H,ERR_L
 rng=random.Random(SEED)
 cert,CERT_LINES=load_certificate()
 ERR_H,ERR_L=cert.ERRH,cert.ERRL # 3.3's error tables (W0, W1, N') for H, L words
 OUT["input"]=cfg
 OUT["budget"]=dict(seconds=BUDGET,reserve=RESERVE,minimums=dict(class_rule_random=MIN_RANDOM,
 class_rule_near_tie=MIN_NEAR,class_labels=MIN_CLASSES))
 OUT["A_identities"]=part_a_identities(rng)
 if DO_KERNELS:
  OUT["A_kernels"]=a=part_a_kernels(random.Random(SEED+1299709),cert)
  OUT["B_routing"]=b=machine.part_b(sys.modules[__name__],cert)
  OUT["counts_match_certificate"]=(
 a["counts_ok"]and a["keys_formula_holds"]and a["fwd_bwd_formulas_hold"]and a["all_correct"]
 and a["window_batch"]["ok"]and a["boundary_batch"]["ok"]
 and a["forced_fixups"]["ok"]and OUT["A_identities"]["class_labels"]["ok"]
 and OUT["A_identities"]["class_rule"]["mismatches_outside_window_set"]==0
 and a["site_checks"]["fixups_ok"]and a["site_checks"]["all_ok"]and a["fixups"]["flags_dead_at_every_join"]
 and a["fe_inv"]["counts"]==[cert.INV]and a["fe_inv"]["correct"]and a["ygen"]["ok"]and a["program_instructions"]["total"]==cert.KERNEL_INSTR and b["formulas_ok"]
 and b["block_paths"]["ok"]and b["grouped_paths"]["ok"]and b["ls64"]["ok"]and b["isa"]["ok"]
 and a["sites_match_certificate"]
 and a["tfp_terms"]["sums_exact"]
 and all(lo<=a["tfp_terms"]["ranges"][k][0]and a["tfp_terms"]["ranges"][k][1]<=hi for k,(lo,hi)in TERMS.items())
 and a["d_zero_batch"]["acc_zero"]and a["d_zero_batch"]["all_other_fractions_wrong"]and OUT["A_identities"]["p1_identity"]["ok"]
 and a["final_R_inverse_rho_and_unused"]and a["loop_ends_on_bit_34"]and a["kernel_xs_from_psel"]["ok"])
 OUT["processor_views"]=machine.view_checks(sys.modules[__name__])
 OUT["counts_match_certificate"]=OUT.get("counts_match_certificate")and OUT["processor_views"]["ok"]
 OUT["phase0"]=p0=dict(constructed=phase0_toy(random.Random(SEED+611953)),generator=phase0_gen(random.Random(SEED+7)),
 design=design_check(random.Random(SEED+31)),resolver=resolver_swap(random.Random(SEED+97)))
 OUT["canon_key"]=ck=[canon_key(cert.GAPE,cert.PAD,cert.WRAPT)]
 OUT["counts_match_certificate"]=OUT["counts_match_certificate"]and all(c["ok"]for c in ck)and all(v["ok"]for v in p0.values())
 OUT["D_certificate"]=part_d(cert)
 OUT["seconds"]=round(time.time()-T0,2)
 print(json.dumps(OUT,indent=1))
if __name__=="__main__":main()
