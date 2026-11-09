import contextlib,functools,io,json,os,random,sys,time
random.Random.rr=random.Random.randrange
from collections import defaultdict
from fractions import Fraction
T0=time.time()
cfg=json.load(sys.stdin)
SEED=cfg.get("seed",1)
BUDGET=cfg.get("seconds",12)
DO_KERNELS=cfg.get("kernels",True)
RESERVE=min(14.5,BUDGET-0.3)
T_END=T0+BUDGET-RESERVE
MIN_RANDOM,MIN_NEAR,MIN_CLASSES=300,100,300
def package_root():
 here=os.path.dirname(os.path.abspath(__file__))
 for base in(here,os.path.join(here,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(base,"kernels","a64.py")):return base
 raise FileNotFoundError("kernels/a64.py not found")
HERE=package_root()
sys.path.insert(0,os.path.join(HERE,"kernels"))
import a64
a64.a64_load=functools.lru_cache(maxsize=None)(a64.a64_load)
sys.path.insert(0,HERE)
import sites,machine
M64=2**64-1
M41,M42=2**40-1,2**41-1
PART=105
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
FRAC_OFF=((64,80,128,144,192,208),(96,112,160,176,224,240))
ENTRY,ACC,KEYS0,MX,Q2M=512,32,64,8,2**24-1
YR_B4,YR_B,YR_B2,YR_YROW,YR_NT=256,288,64,128,320
W6S=(0,43,86,128,171,214)
W6M=tuple((1<<b)-1 for b in(43,43,42,43,43,42))
KO=(1,2,3,4,6,7)
BT=(1,BETA,BETA*BETA%p)
MH_A,MH_B,ML=20,21,23
TERMS=dict(e=1,carry=7,delta=1,h=8,l=9)
ERR_H,ERR_L=((-20,26),(-20,26),(-53,40)),((-22,29),(-22,29),(-59,44))
OUT,CERT_LINES={},[]
class Curve:
 def __init__(self,q,order):self.q,self.r=q,order
 def add(self,P,Q):
  q=self.q
  if P is None:return Q
  if Q is None:return P
  if P[0]==Q[0]:
   if(P[1]+Q[1])%q==0:return None
   lam=3*P[0]*P[0]*pow(2*P[1],-1,q)%q
  else:lam=(Q[1]-P[1])*pow(Q[0]-P[0],-1,q)%q
  x=(lam*lam-P[0]-Q[0])%q
  return(x,(lam*(P[0]-x)-P[1])%q)
 def neg(self,P):return None if P is None else(P[0],(-P[1])%self.q)
 def mul(self,k,P):
  k%=self.r
  if not k or P is None:return None
  q,X,Y,Z=self.q,P[0],P[1],1
  for b in bin(k)[3:]:
   if X is not None:
    if Y==0:X=None
    else:S,M=4*X*Y*Y%q,3*X*X%q;X2=(M*M-2*S)%q;Y2=(M*(S-X2)-8*pow(Y,4,q))%q;Z=2*Y*Z%q;X,Y=X2,Y2
   if b=="1":
    if X is None:X,Y,Z=P[0],P[1],1
    else:
     Z2=Z*Z%q;U2,S2=P[0]*Z2%q,P[1]*Z2*Z%q
     if X==U2:
      if(Y+S2)%q==0 or Y==0:X=None
      else:S,M=4*X*Y*Y%q,3*X*X%q;X2=(M*M-2*S)%q;Y2=(M*(S-X2)-8*pow(Y,4,q))%q;Z=2*Y*Z%q;X,Y=X2,Y2
     else:
      H,R=(U2-X)%q,(S2-Y)%q;H2=H*H%q;H3=H*H2%q;V=X*H2%q
      X=(R*R-H3-2*V)%q;Y=(R*(V-X)-Y*H3)%q;Z=Z*H%q
  if X is None or Z%q==0:return None
  iz=pow(Z,-1,q);iz2=iz*iz%q
  return(X*iz2%q,Y*iz2*iz%q)
PROD=Curve(p,r);GP=(GX,GY)
def limbs(v,n=4):return[(v>>(64*k))&M64 for k in range(n)]
def put(mem,a,v,n=4):
 for k,w in enumerate(limbs(v,n)):mem[a+8*k]=w
def get(mem,a,n=4):return sum(mem[a+8*k]<<(64*k)for k in range(n))
def rgi(regs,names):return sum(regs[n]<<(64*k)for k,n in enumerate(names))
def ars(v,names):return{n:w for n,w in zip(names,limbs(v,len(names)))}
def exk(x):
 v=[x,BETA*x%p,BETA*BETA%p*x%p]
 j=min(range(3),key=lambda t:v[t])
 return v[j],j
def t128(v):return(v<<128)//p
def t64(v):return(v<<64)//p
def ws(*v):return sum(x>>64 for x in v)%2**64<<64|sum(x&M64 for x in v)%2**64
def wneg(x):return ws(x^(M128-1),1<<64|1)
def rfields(a,b):
 ha,la,hb,lb=a>>64,a&M64,b>>64,b&M64
 d1,d2=la>>23,lb>>23
 return d1,d2,(ha>>41<<41)|d2,(hb>>21<<21)|((ha>>20)&((1<<21)-1))
def label(v):return rfields(t128(v),t128(BETA*v%p))
def class_label(x):return label(exk(x)[0])
def neg_t(v):return wneg(t128(v%p))
def corr(v):return neg_t(v),neg_t(BETA*v)
def fold(z):return(z%2**256+(z>>256)*C0)%2**256
def lkfold(z):
 z=z%2**256+(z>>256)*C0
 return z-2**256+C0 if z>>256 else z
def fs_op(a,b):return(a-b if a>=b else a-b+p)%2**256
def faddw_op(a,b):return(a+b if a+b<2**256 else a+b-p)%2**256
def shr(c):
 R=[c]
 for _ in range(3):R.append((R[-1]%2**192*2**64+(R[-1]>>192)*C0)%2**256)
 return R
def tfp(R,M,tm=None):
 m,w=limbs(M),[limbs(x)for x in R]
 w.append([0,0,w[3][1],w[3][2]])
 H,L=[sum((m[k]*w[k+o][3]&M64)+(m[k]*w[k+o][2]>>64)for k in range(4))&M64 for o in(0,1)]
 if tm is not None:
  Z=sum(m[k]*R[k]for k in range(4))
  Z1=sum(m[k]*R[k+1]for k in range(3))+m[3]*(R[3]<<64&(2**256-1))
  z=Z%p
  e,c,e1,c1,dl=(t64(z)-(Z>>192))%2**64,((Z>>192)-H)%2**64,(t64(Z1%p)-(Z1>>192))%2**64,((Z1>>192)-L)%2**64,(t64(z<<64)%2**64-t64(Z1%p))%2**64
  th,tl=(t64(z)-H)%2**64,(t64(z<<64)-L)%2**64
  for k,v in(("e",max(e,e1)),("carry",max(c,c1)),("delta",dl),("h",th),("l",tl)):tm[k]=max(tm.get(k,0),v)
  tm["sum_ok"]=tm.get("sum_ok",True)and th==e+c and tl==dl+e1+c1
 return H<<64|L
def fracs(mem,e,s):return[get(mem,e+o,2)for o in FRAC_OFF[s]]
def components(X,Y):
 a,b=X[0],Y[0]
 E,F,G=(b**4+2*a**3*b)%p,(a**4+2*a*b**3)%p,3*a*a*b*b%p
 out=[]
 for w in(pow(Y[1]+X[1],-2,p),pow(Y[1]-X[1],-2,p)):out+=[E*w%p,F*w%p,G*w%p]
 return out
def werr(t,v):return[((t>>64*k)-(v>>64*k)+2**63)%2**64-2**63 for k in(1,0)]
def frac_errs(mem,e,X,Y):
 cs,out=components(X,Y),[]
 for s in(0,1):
  ex=((cs[3*s]-Y[0])%p,(cs[3*s+1]-X[0])%p,cs[3*s+2])
  out+=[werr(t128(BT[i%2]*ex[i//2]%p),v)for i,v in enumerate(fracs(mem,e,s))]
 return out
def leaf_w(f6,t):
 u=[[f6[2*i],f6[2*i+1],wneg(ws(f6[2*i],f6[2*i+1]))]for i in range(3)]
 return[ws(u[0][m],u[1][(m+t)%3],u[2][(m+2*t)%3])for m in(0,1)]
def record_fields(w):
 f0,f1=w[0]>>64,w[1]>>64
 s=(f0+f1)&M64
 j=(2 if f1+s>M64 else 1)if f0>=f1 else(2 if f0+s>M64 else 0)
 n=ws(w[0],w[1])^(M128-1)
 return(j,)+rfields(*((w[0],w[1]),(w[1],n),(n,w[0]))[j])
def in_W(x):
 v=[x<<64,(BETA*x%p)<<64,(BETA*BETA%p*x%p)<<64];top=p<<64;lo,hi=ERR_H[0]
 return(any(w<hi*p or w>=top+lo*p for w in v[:2])or v[2]<(1-2*lo)*p or v[2]>top-(3+2*hi)*p or abs(v[0]-v[1])<(hi+1-lo)*p or any(3*lo*p<=w-v[2]<(3*hi+3)*p for w in v[:2]))
def in_D(x):
 v=(x,BETA*x%p,BETA*BETA%p*x%p)
 return any(u%2**mb<hi or u%2**mb>=2**mb+lo for w,eh,el in zip(v,ERR_H,ERR_L)for u,mb,(lo,hi)in((t128(w)>>64,MH_A,eh),(t128(w)&M64,ML,el)))
def dft_leaf(c3,t,a,b):
 return(c3[0]-b+pow(BETA,t,p)*(c3[1]-a)+pow(BETA,2*t,p)*c3[2])%p
def inverse_dft(L):
 i3,b2=pow(3,-1,p),BETA*BETA%p
 return((L[0]+L[1]+L[2])*i3%p,(L[0]+b2*L[1]+BETA*L[2])*i3%p,(L[0]+BETA*L[1]+b2*L[2])*i3%p)
def weak(rg):return rg.rr(p)if rg.random()<0.85 else rg.rr(C0)+p
def represent(x,t,a,b,rg,tm=None):
 w,c,k=weak(rg)or 1,[rg.rr(p)for _ in range(3)],rg.rr(3);c[k]=0
 v=(x+b+BT[t]*a)*pow(w,-1,p)-sum(BT[t*i%3]*c[i]for i in range(3))
 c[k]=v*pow(BT[t*k%3],-1,p)%p
 bw=BETA*w%p;bw+=p if bw<C0 and rg.random()<0.5 else 0
 out=[]
 for ci,y in zip(c,(b,a,None)):
  R,cw=shr(ci),corr(y)if y else(0,0)
  out+=[ws(tfp(R,w,tm),cw[0]),ws(tfp(R,bw,tm),cw[1])]
 return out,((c[0]*w-b)%p,(c[1]*w-a)%p,c[2]*w%p)
def edge_leaf(rg,kind):
 m=rg.rr(3)
 if kind=="D":
  b=rg.rr(2);mb=(MH_A,ML)[b]
  u=rg.getrandbits(64-mb)<<mb|rg.choice((rg.rr(70),2**mb-1-rg.rr(70),rg.rr(2**mb)))
  v= -(-(u<<64|rg.getrandbits(64)if b==0 else rg.getrandbits(64)<<64|u)*p//M128)
  return v*pow(BETA,-m,p)%p or 1
 off=rg.rr(160<<192);w=off if rg.random()<0.5 else p-1-off
 if kind=="Xm":return w*pow(BETA,-m,p)%p or 1
 return w*pow((1-BETA,1-BETA*BETA,BETA-BETA*BETA)[m]%p,-1,p)%p or 1
def deadline(frac):now=time.time();return now+frac*max(0.0,T_END-now)
def near_tie_x(rg):
 mlt=rg.choice([1,BETA,BETA*BETA,1-BETA,1-BETA*BETA,BETA-BETA*BETA])
 dd=rg.rr(-(1<<rg.rr(150,250)),1<<rg.rr(150,250))
 return dd*pow(mlt%p,-1,p)%p or 1
def index_check(x,rg):
 t,a,b=rg.rr(3),rg.rr(1,p),rg.rr(1,p)
 j=record_fields(leaf_w(represent(x,t,a,b,rg)[0],t))[0]
 return in_W(x),j!=exk(x)[1]
def class_labels(rg):
 out=dict(classes=0,edge_classes=0,members=[0,0,0],in_W=0,in_D=0,wrong_index=0,mismatches=0,errors_h=[[99,-99]for _ in range(3)],errors_l=[[99,-99]for _ in range(3)],tfp_terms={})
 dl=deadline(0.05)
 while out["classes"]<MIN_CLASSES or(out["classes"]<4*MIN_CLASSES and time.time()<dl):
  edge=out["classes"]%3==2
  x0=edge_leaf(rg,("D","Xm","diff")[out["classes"]//3%3])if edge else rg.rr(1,p)
  v=exk(x0)[0];out["classes"]+=1;out["edge_classes"]+=edge;lab=label(v)
  for m in range(3):
   x=v*pow(BETA,(3-m)%3,p)%p;t,a,b=rg.rr(3),rg.rr(1,p),rg.rr(1,p)
   w=leaf_w(represent(x,t,a,b,rg,out["tfp_terms"])[0],t)
   for k,c in enumerate((w[0],w[1],ws(w[0],w[1])^(M128-1))):
    eh,el=werr(t128(pow(BETA,k,p)*x%p),c)
    out["errors_h"][k]=[min(out["errors_h"][k][0],eh),max(out["errors_h"][k][1],eh)]
    out["errors_l"][k]=[min(out["errors_l"][k][0],el),max(out["errors_l"][k][1],el)]
   f=record_fields(w)
   if in_W(x):out["in_W"]+=1;continue
   if f[0]!=m:out["wrong_index"]+=1;continue
   out["members"][m]+=1
   if in_D(x):out["in_D"]+=1
   elif f[1:]!=lab:out["mismatches"]+=1
 tm=out["tfp_terms"]
 out["ok"]=(out["mismatches"]==0 and out["wrong_index"]==0 and out["in_D"]>0 and out["in_W"]>0 and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors_h"],ERR_H))and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors_l"],ERR_L))and tm["sum_ok"]and all(tm[k]<=v for k,v in TERMS.items()))
 return out
def part_a_identities(rng):
 res={}
 res["beta_cube_root_of_1_mod_p"]=pow(BETA,3,p)==1 and BETA!=1
 res["lambda_cube_root_of_1_mod_r"]=pow(LAM,3,r)==1 and LAM!=1
 res["lambda_G_equals_beta_x"]=PROD.mul(LAM,GP)==(BETA*GX%p,GY)
 res["seven_nonresidue_no_x0_point"]=pow(7,(p-1)//2,p)==p-1
 res["window_forms_nonzero"]=all(c%p for c in(1-BETA,1-BETA*BETA,BETA-BETA*BETA,3))
 res["label_constants"]=(FA*FA-FA*FB+FB*FB==p and FB*BETA%p==FA-FB and(FA-FB)*BETA%p==p-FA)
 bad=drops=wrong=done=0;kinds={"random":0,"near_tie":0};dl=deadline(0.10)
 while done<26000 and(done<MIN_RANDOM or time.time()<dl):
  x,kind=(rng.rr(1,p),"random")if done<20000 else(near_tie_x(rng),"near_tie")
  kinds[kind]+=1;done+=1;win,err=index_check(x,rng)
  drops+=win;wrong+=win and err;bad+=err and not win
 extra=max(0,MIN_NEAR-kinds["near_tie"]);rng2=random.Random(SEED+7919)
 for _ in range(extra):
  win,err=index_check(near_tie_x(rng2),rng2);kinds["near_tie"]+=1;drops+=win;wrong+=win and err;bad+=err and not win
 res["class_rule"]=dict(checked=kinds,in_window_set=drops,wrong_index_in_window_set=wrong,mismatches_outside_window_set=bad,near_ties_from_second_stream=extra)
 res["class_labels"]=class_labels(random.Random(SEED+104729))
 return res
SNT=0x5E5E5E5E5E5E5E5E
TIB=1<<40
A_SP=(34<<40)+0x7000000
EB_BASE=(1<<45)+(1<<40)
EB_END=EB_BASE+(1<<34)
YR0,YR_STRIDE,XT0,BT0=(1<<45)+(1<<31),1<<31,24576*TIB,0xfff8<<48
XT_END=(24576+12*2560)*TIB
AXYX0=4096*TIB
AXYX_HI=20480*TIB
D6O=tuple(c<<48 for c in(0,2,8,4,6,9))
def ax_bases(j):
 return(AXYX0+((4*j+0)<<49),AXYX0+((4*j+1)<<49),AXYX_HI+((2*j+0)<<48),AXYX0+((4*j+2)<<49),AXYX0+((4*j+3)<<49),AXYX_HI+((2*j+1)<<48))
SLOTW=(0,1,3,2,4,5,3)
def offs6(v):return[D6O[k]+(((v>>sh)&m)<<6)for k,(sh,m)in enumerate(zip(W6S,W6M))]
def xcs(a):return 2*a**3%p,2*a%p,a*a%p,pow(a,4,p)
def y_words(b):return b,pow(b,3,p),3*b*b%p
def bt_addr(w,i):
 t=w>>43*i&(1<<43)-1
 return BT0+(t<<8)+32*SLOTW.index(i)+(128*(w>>128&1)if i==3 else 0)
def blook(w):return lkfold(sum(BETA*(w>>43*i&(1<<43)-1)*pow(2,43*i,p)%p for i in range(6)))
def lookup6(cval,w,add0=0):
 return lkfold(sum(((cval*((w>>sh)&m)*pow(2,sh,p)+(add0 if k==0 else 0))%p)for k,(sh,m)in enumerate(zip(W6S,W6M))))
class _TraceLC1Dict(dict):
 def __init__(self,rprog,rlab):
  super().__init__()
  self._rp,self._rl,self.touched=rprog,rlab,set()
 def __getitem__(self,a):
  self.touched.add(a)
  return super().__getitem__(a)
 def __setitem__(self,a,v):
  self.touched.add(a)
  return super().__setitem__(a,v)
 def get(self,a,d=None):
  if a in self:self.touched.add(a)
  return super().get(a,d)
 def __missing__(self,a):
  if 0<=a<16*machine.H1 and a%8==0:
   d=a//8;b=machine.blk1(d)
   a64.a64_run(self._rp,self._rl,"init3",dict(x0=1,x1=a,x2=b+32,x3=machine.B1,x30=SNT),self)
   self.touched.add(a);return super().__getitem__(a)
  for lo,hi,first,stride in machine.LC1_TABLES[2:]:
   if lo<=a<hi:
    self.touched.add(a);super().__setitem__(a,0);return 0
  raise KeyError(hex(a))
def pop60(rprog,rlab,mem,tb,c1,ad1,c2,ad2,w):
 DSC=(34<<40)+0x180000
 for k,(sh,m)in enumerate(zip(W6S,W6M)):
  tk=(w>>sh)&m;ad=tb+D6O[k]+(tk<<6)
  if ad in mem:continue
  put(mem,DSC,0);put(mem,DSC+32,(c1*tk*pow(2,sh,p)+(ad1 if k==0 else 0))%p)
  put(mem,DSC+64,0);put(mem,DSC+96,(c2*tk*pow(2,sh,p)+(ad2 if k==0 else 0))%p)
  a64.a64_run(rprog,rlab,"tbld",dict(x0=ad,x1=DSC,x2=1,x3=1,x4=2,x5=64,x6=0,x30=SNT),mem)
def cube_root(v):
 c_=pow(v,(p+2)//9,p)
 return c_ if c_*c_*c_%p==v%p else None
def make_bday_xs(rg):
 xs=[]
 for _ in range(MX//2):
  while True:
   X0=PROD.mul(rg.rr(1,r),GP);A0=X0[1]*X0[1]%p;d=rg.rr(1,1<<63);A1=A0+d
   if A1>=p or(A0>>63)!=(A1>>63):continue
   y1,a1=sqrt_mod_p(A1),cube_root((A1-7)%p)
   if y1 is not None and a1 is not None:xs+=[X0,(a1,y1)];break
 return xs
def popax(rprog,rlab,mem,j,X,w):
 DSC=(34<<40)+0x182000;cax,cyx,bs=2*(X[1]*X[1]%p)%p,2*X[1]%p,ax_bases(j)
 for k,(sh,m)in enumerate(zip(W6S,W6M)):
  tk=(w>>sh)&m;ad=bs[k]+(tk<<6)
  if ad in mem:continue
  put(mem,DSC,0);put(mem,DSC+32,cax*tk*pow(2,sh,p)%p)
  put(mem,DSC+64,0);put(mem,DSC+96,cyx*tk*pow(2,sh,p)%p)
  a64.a64_run(rprog,rlab,"tbld",dict(x0=ad,x1=DSC,x2=1,x3=1,x4=2,x5=64,x6=0,x30=SNT),mem)
def popbt(rprog,rlab,mem,w):
 DSC=(34<<40)+0x183000
 for t_val in{((w>>(43*i))&((1<<43)-1))for i in range(6)}:
  ba=BT0+(t_val<<8)
  if ba in mem:continue
  for idx,wi in enumerate(SLOTW):
   tk=t_val&((1<<41)-1)if idx==5 else t_val
   put(mem,DSC+64*idx,0);put(mem,DSC+64*idx+32,BETA*tk*pow(2,43*wi,p)%p)
  a64.a64_run(rprog,rlab,"tbld",dict(x0=ba,x1=DSC,x2=1,x3=1,x4=7,x5=256,x6=0,x30=SNT),mem)
  put(mem,ba+224,0)
def popbat(rprog,rlab,mem,Xs,Ys,rho):
 acc=rho;acc_pairs,e_pairs=[],[]
 for yi,Y in enumerate(Ys):
  b,yy=Y;B=yy*yy%p;b3,tb2=pow(b,3,p),3*b*b%p
  for q in range(len(Xs)//4):
   a0,a2=Xs[4*q][0],Xs[4*q+2][0]
   pop60(rprog,rlab,mem,XT0+q*2560*TIB,2*pow(a0,3,p)%p,0,2*pow(a2,3,p)%p,0,b)
   A1,A2,A3,A4=Xs[4*q][1]**2%p,Xs[4*q+1][1]**2%p,Xs[4*q+2][1]**2%p,Xs[4*q+3][1]**2%p
   sA,pA,sB,pB=(A1+A2)%p,(A1*A2)%p,(A3+A4)%p,(A3*A4)%p
   pop60(rprog,rlab,mem,XT0+(10+q)*2560*TIB,(-sB)%p,(pB-7*sB)%p,(-sA)%p,(pA-7*sA)%p,b3)
  for pr in range(len(Xs)//2):
   (a1,y1),(a2,y2)=Xs[2*pr],Xs[2*pr+1]
   c10,c11,c12,a14=xcs(a1);c20,c21,c22,a24=xcs(a2)
   pop60(rprog,rlab,mem,XT0+(2+pr)*2560*TIB,c11,a14,c21,a24,b3)
   pop60(rprog,rlab,mem,XT0+(6+pr)*2560*TIB,c12,0,c22,0,tb2)
   A1p,A2p=y1*y1%p,y2*y2%p
   sAp,pAp=(A1p+A2p)%p,(A1p*A2p)%p
   ep_v=lkfold(sum((((-sAp)%p*((b3>>sh)&m)*pow(2,sh,p)+((pAp-7*sAp)%p if k==0 else 0))%p)for k,(sh,m)in enumerate(zip(W6S,W6M)))+(B*B%p))
   acc_pairs.append(acc);e_pairs.append(ep_v)
   acc=fold(sum(((ep_v>>64*k)&M64)*shr(acc)[k]for k in range(4)))
 R=pow(acc%p,p-2,p)if acc%p else 0
 prod=lambda m_,rr:sum(((m_>>64*k)&M64)*rr[k]for k in range(4))
 for jp in range(len(acc_pairs)-1,-1,-1):
  yi=(2*jp)//len(Xs);b,yy=Ys[yi];B=yy*yy%p
  pr=jp%(len(Xs)//2);(a1,y1),(a2,y2)=Xs[2*pr],Xs[2*pr+1]
  A1,A2=y1*y1%p,y2*y2%p;dA=(A2-A1)%p
  RR=shr(R);inv_e,R=fold(prod(acc_pairs[jp],RR)),fold(prod(e_pairs[jp],RR))
  inv1=fold(prod(inv_e,shr(fs_op(B,A2))))
  trf=dA*inv_e;inv2=faddw_op(inv1,fold(trf%2**256+(trf>>256)*C0))
  for j,X,inv in((2*pr+1,(a2,y2),inv2),(2*pr,(a1,y1),inv1)):
   isq=fold(inv*inv%2**256+(inv*inv>>256)*C0)
   popax(rprog,rlab,mem,j,X,isq)
   Mv=lkfold(inv+lookup6(2*(X[1]*X[1]%p)%p,isq))
   Pv=lookup6(2*X[1]%p,isq)
   Nv=fold(sum(((Pv>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in range(4)))
   popbt(rprog,rlab,mem,fs_op(Mv,Nv))
   popbt(rprog,rlab,mem,faddw_op(Mv,Nv))
def x_words(mem,e,Xs,l):
 j=l%MX;pr=j>>1;X=Xs[j]
 for k,w in enumerate(corr(X[0])):put(mem,e+256+16*k,w,2)
 put(mem,e+288,X[1]*X[1]%p)
 for k,b_ in enumerate(ax_bases(j)):mem[e+320+8*k]=b_
 if j&1:mem[e+368]=(Xs[2*pr+1][1]**2-Xs[2*pr][1]**2)%p
 if j%4==0:
  q=j//4;d0=(Xs[4*q+1][1]**2-Xs[4*q][1]**2)%p;d1=(Xs[4*q+3][1]**2-Xs[4*q+2][1]**2)%p
  mem[e+376],mem[e+384],mem[e+392]=XT0+q*2560*TIB,2*d0,2*d1
 if(j&1)==0:mem[e+400],mem[e+408]=XT0+(2+pr)*2560*TIB,XT0+(6+pr)*2560*TIB
 if j==0:mem[e+416],mem[e+424]=XT0+10*2560*TIB,XT0+11*2560*TIB
def put_yrec(mem,e,Y):
 b,y=Y;B=y*y%p
 mem[e]=0;put(mem,e+YR_B4,pow(b,4,p));put(mem,e+YR_B,B);put(mem,e+YR_B2,B*B%p)
 for k in range(4):put(mem,e+YR_YROW+32*k,y*pow(2,64*k,p)%p)
 put(mem,e+YR_NT,neg_t(b),2);put(mem,e+YR_NT+16,neg_t(BETA*b),2)
 ob,ob3,o3b2=offs6(b),offs6(pow(b,3,p)),offs6(3*b*b%p)
 for k in range(4):mem[e+96+8*k]=ob3[k]
 for k in range(2):mem[e+352+8*k]=ob3[4+k]
 bl=limbs(b,4);mem[e+368],mem[e+376]=bl[0],bl[1]
 for k in range(6):mem[e+384+8*k]=ob[k]
 for k in range(6):mem[e+440+8*k]=o3b2[k]
 mem[e+496],mem[e+504]=bl[2],bl[3]
def bwd_model(Xs,Ys,mem,acc_pairs,e_pairs,R,tm=None,ybase=None):
 if ybase is None:ybase=YR0
 out,prod={},lambda m,rr:sum((m>>64*k&M64)*rr[k]for k in range(4))
 for jp in range(len(acc_pairs)-1,-1,-1):
  l1,l2=2*jp,2*jp+1
  (a1,yx1),(a2,yx2)=Xs[l1%MX],Xs[l2%MX]
  b,yy=Ys[l2//MX]
  B,A1,A2=yy*yy%p,yx1*yx1%p,yx2*yx2%p
  RR=shr(R)
  inv_e,R=fold(prod(acc_pairs[jp],RR)),fold(prod(e_pairs[jp],RR))
  D2=fs_op(B,A2);dA=(A2-A1)%p
  inv1=fold(prod(inv_e,shr(D2)))
  trf=dA*inv_e;inv2=faddw_op(inv1,fold(trf%2**256+(trf>>256)*C0))
  for l,(a,yx),inv in((l2,(a2,yx2),inv2),(l1,(a1,yx1),inv1)):
   b_l,yy_l=Ys[l//MX];AX=yx*yx%p
   yr=ybase+512*(l//MX)
   yrows=[get(mem,yr+YR_YROW+32*k)for k in range(4)]
   inv_sq=fold(inv*inv%2**256+(inv*inv>>256)*C0)
   Mval=lkfold(inv+lookup6(2*AX%p,inv_sq))
   Pval=lookup6(2*yx%p,inv_sq)
   Nval=fold(prod(Pval,yrows))
   w=[fs_op(Mval,Nval),faddw_op(Mval,Nval)]
   bw=[blook(v)for v in w]
   c0,c1,c2,a4=xcs(a);b4=get(mem,yr+YR_B4);yw=y_words(b_l)
   if(l%2)==0:e_val=lkfold(sum(((c0*((b_l>>sh)&m)*pow(2,sh,p))%p)for sh,m in enumerate(W6M)if False)if False else sum(((c0*((b_l>>sh)&m)*pow(2,sh,p))%p)for sh,m in zip(W6S,W6M))+b4)
   else:
    c0e=xcs(a1)[0];eev=lkfold(sum(((c0e*((b_l>>sh)&m)*pow(2,sh,p))%p)for sh,m in zip(W6S,W6M))+b4)
    te=(2*dA)*b_l;e_val=faddw_op(eev,fold(te%2**256+(te>>256)*C0))
   cs=[e_val]+[lkfold(sum(((ci_v*((wv>>sh)&m)*pow(2,sh,p)+(a4 if ci==1 and k==0 else 0))%p)for k,(sh,m)in enumerate(zip(W6S,W6M))))for ci,(ci_v,wv)in enumerate(zip((c0,c1,c2),yw))if ci>0]
   fr={}
   for i,(c,yc)in enumerate(zip(cs,(b_l,a,None))):
    Rc,cw=shr(c),corr(yc)if yc else(0,0)
    for s in(0,1):
     fr[s,2*i]=ws(tfp(Rc,w[s],tm),cw[0])
     fr[s,2*i+1]=ws(tfp(Rc,bw[s],tm),cw[1])
   out[l]=[fr[s,i]for s in(0,1)for i in range(6)]
 return out,R
def six_leaves(X,Y):
 out=[]
 for t in range(3):
  PX=(pow(BETA,t,p)*X[0]%p,X[1])
  out+=[PROD.add(PX,Y)[0],PROD.add(PX,PROD.neg(Y))[0]]
 return out
def fixup_extra(prog):
 out,last={},-1
 for i,(mn,_)in enumerate(prog):
  if mn in("b.lo","b.hs","cbnz"):last=i
  elif mn=="stadd":out[i]=i-last
 return out
def flags_dead_after_fixups(prog,labels):
 reads=("adc","adcs","sbc","sbcs","csel","cinc","cset","ccmp","ccmn","b.eq","b.ne","b.lo","b.hs")
 def target(t,pc):
  if t[: -1].isdigit()and t[-1:]in("f","b"):
   c_=labels[t[: -1]]
   return[q for q in c_ if q<pc][-1]if t[-1]=="b"else[q for q in c_ if q>=pc][0]
  return labels[t][0]
 for i in fixup_extra(prog):
  todo,seen=[i+1],set()
  while todo:
   j=todo.pop()
   if j in seen:continue
   seen.add(j);mn,o=prog[j]
   if mn in reads:return False
   if mn in("adds","subs","cmp","cmn","tst","ret"):continue
   if mn in("b","tbz","tbnz","cbz","cbnz"):
    todo.append(target(o[-1],j+1))
    if mn=="b":continue
   todo.append(j+1)
 return True
CBS,CB0,YSP=2304,(34<<40)+0x2000000,(34<<40)+0x4000000
def chain_n(i):
 rho,q=divmod(i,33)
 return 33*rho+(0 if q==32 else(16-q//2)*(1-2*(q%2)))
@functools.lru_cache(maxsize=None)
def chain_y(ch,i):return PROD.add(ch[0],PROD.mul(chain_n(i)%r,ch[1]))
def y_generator(fprog,flabels,rng,cert):
 mem,chains={},[]
 for c in range(2):
  C,S=PROD.mul(rng.rr(1,r),GP),PROD.mul(rng.rr(1,r),GP);chains.append((C,S));cb=CB0+CBS*c
  for i,P in[(0,C)]+[(i,PROD.mul(i,S))for i in list(range(1,17))+[33]]:
   put(mem,cb+(1088 if i==33 else 64*i),P[0]);put(mem,cb+(1088 if i==33 else 64*i)+32,P[1])
  mem[cb+2240]=YR0+c*YR_STRIDE
 put(mem,CB0-128,1)
 for a in range(YSP-4096,YSP+1024,8):mem[a]=0
 put(mem,YSP+352,BETA);mem[YSP+384],mem[YSP+424]=C0,CB0
 counts=[]
 for nch in(2,1):
  mem[YSP+416]=CB0+CBS*nch
  counts.append(a64.a64_run(fprog,flabels,"yround",dict(sp=YSP,x30=SNT),mem,limit=10**7))
 bad=nrec=0
 for c,rd in((0,0),(1,0),(0,1)):
  for q in range(33):
   i=33*rd+q;(b,y),e=chain_y(chains[c],i),YR0+c*YR_STRIDE+512*i;B=y*y%p
   ok=(get(mem,e+YR_B4)%p==pow(b,4,p)and get(mem,e+YR_B)==B and get(mem,e+YR_B2)==B*B%p)
   ok&=all(get(mem,e+YR_YROW+32*k)==y*pow(2,64*k,p)%p for k in range(4))
   ok&=get(mem,e+YR_NT,2)==neg_t(b)and get(mem,e+YR_NT+16,2)==neg_t(BETA*b)
   ob,ob3,o3b2=offs6(b),offs6(pow(b,3,p)),offs6(3*b*b%p)
   ok&=([mem[e+96+8*k]for k in range(4)]+[mem[e+352+8*k]for k in range(2)])==ob3
   ok&=[mem[e+384+8*k]for k in range(6)]==ob
   bl_=limbs(b,4);ok&=[mem[e+440+8*k]for k in range(6)]==o3b2 and(mem[e+368],mem[e+376],mem[e+496],mem[e+504])==(bl_[0],bl_[1],bl_[2],bl_[3])
   bad+=not ok;nrec+=1
 cen=[(get(mem,CB0+CBS*c),get(mem,CB0+CBS*c+32))for c in range(2)]
 ok=bad==0 and cen==[chain_y(chains[0],33*2+32),chain_y(chains[1],33+32)]
 cok=counts==[cert.YR_FIX+2*cert.YR_BLK,cert.YR_FIX+cert.YR_BLK]
 return dict(records=nrec,mismatches=bad,centres_ok=ok,counts=counts,counts_ok=cok),mem,chains[0],chains[1]
def leaf_batch(prog,labels,Xs,Ys,q2,ymem,cert,rprog=None,rlab=None,W2=0x123456789,KN=5<<26,override=None,chain=None,rho=None,ybase=YR0):
 assert q2%32==0 and ybase%YR_STRIDE==0
 if rprog is None:rprog,rlab=a64.a64_load(open(os.path.join(HERE,"kernels","route.s")).read())
 fx=fixup_extra(prog);vis={k:dict.fromkeys(fx,0)for k in("fwd","bwd")}
 nY,SP=q2//MX,A_SP
 sb=EB_END-ENTRY*q2
 LS=([(AXYX0,XT_END),(BT0,1<<64),(ybase,ybase+512*nY),(sb,EB_END)],[(sb,EB_END)])
 mem=_TraceLC1Dict(rprog,rlab)
 rho=rho or random.Random(q2).rr(1,p)
 popbat(rprog,rlab,mem,Xs,Ys[:nY],rho)
 if ymem is not None:
  for a in range(ybase-64,ybase+512*nY,8):mem[a]=ymem.get(a,0)
 else:
  for a in range(ybase-64,ybase,8):mem[a]=0
  for yi in range(nY):put_yrec(mem,ybase+512*yi,Ys[yi])
 for l in range(q2):x_words(mem,sb+ENTRY*l,Xs,l)
 regs=dict(x0=q2,x1=ybase+64,x3=sb,sp=SP,x30=SNT);regs.update(ars(rho,("x24","x25","x26","x27")))
 out={};nf=a64.a64_run(prog,labels,"fwd",regs,mem,out_regs=out,visits=vis["fwd"],ls64=LS)
 words,acc=[out["x0"]],rgi(out,("x24","x25","x26","x27"))
 accs=[get(mem,sb+ENTRY*(2*j)+ACC)for j in range(q2//2)]
 es=[get(mem,sb+ENTRY*(2*j))for j in range(q2//2)]
 R0=pow(acc,p-2,p);top=ybase+512*(nY-1)
 regs=dict(x3=EB_END,x1=ybase-32,sp=top,x30=SNT);regs.update(ars(R0,("x0","x2","x14","x15")))
 out={};nb=a64.a64_run(prog,labels,"bwd",regs,mem,out_regs=out,visits=vis["bwd"],ls64=LS)
 words.append(out["x5"])
 want=rho
 for l in range(q2):want=want*(Ys[l//MX][1]**2-Xs[l%MX][1]**2)%p
 final_R=rgi(out,("x0","x2","x14","x15"));sp_kept=out["sp"]==top
 tm={};model,mR=bwd_model(Xs,Ys,mem,accs,es,R0,tm,ybase)
 bad_model,bad_bound,bad_leaf,err,leaves=int(mR!=final_R),0,0,[[99,-99],[99,-99]],[]
 for l in range(q2):
  X,Y=Xs[l%MX],Ys[l//MX];e,cs=sb+ENTRY*l,components(X,Y)
  bad_model+=sum(get(mem,e+o,2)!=w for o,w in zip(FRAC_OFF[0]+FRAC_OFF[1],model[l]))
  for i,(dh,dl)in enumerate(frac_errs(mem,e,X,Y)):
   g=i%6//4;err[0]=[min(err[0][0],dh),max(err[0][1],dh)];err[1]=[min(err[1][0],dl),max(err[1][1],dl)]
   bad_bound+=not(g-1<=dh<=8 and g-1<=dl<=9)
  six=six_leaves(X,Y)
  bad_leaf+=sum(dft_leaf(cs[3*s:3*s+3],t,X[0],Y[0])!=six[2*t+s]for s in(0,1)for t in range(3))
  leaves.append(six)
 for(l,s),(x,t,rg)in(override or{}).items():
  f6,c3=represent(x,t,Xs[l%MX][0],Ys[l//MX][0],rg)
  for o,v in zip(FRAC_OFF[s],f6):put(mem,sb+ENTRY*l+o,v,2)
  for t_ in range(3):leaves[l][2*t_+s]=(c3[0]+BT[t_]*c3[1]+BT[2*t_%3]*c3[2])%p
 port=machine.Port1(mem)
 regs=dict(x2=W2,x3=sb+KEYS0,x26=KN,sp=SP,x30=SNT)
 out={};nk=a64.a64_run(prog,labels,"keys",regs,mem,out_regs=out,dev=port.dev,ls64=([(sb,EB_END)],[]))
 sp_kept&=out["sp"]==SP
 snd=machine.Sender1(cert)
 meter,fields=machine.keys_reference(sys.modules[__name__],mem,sb,q2,W2,KN,snd,cert)
 flat=[leaves[l][2*t+s]for l in range(q2)for s in(0,1)for t in range(3)]
 win=[in_W(x)for x in flat];dset=[in_D(x)for x in flat]
 wrong=sum(not w and f[0]!=exk(x)[1]for w,f,x in zip(win,fields,flat))
 canon=sum(f[1:]!=class_label(x)for w,d,f,x in zip(win,dset,fields,flat)if not w and not d)
 dec=0
 for k,(d,rec)in enumerate(snd.seq if chain else()):
  l=(rec[3]>>2)&((1<<24)-1);s=k//3%2
  Y,X=chain_y(chain,l//MX),Xs[l%MX]
  cand=[class_label(PROD.add((pow(BETA,t,p)*X[0]%p,X[1]),PROD.neg(Y)if s else Y)[0])for t in range(3)]
  dec+=l!=k//6 or Y!=Ys[l//MX]or(not win[k]and not dset[k]and(d,rec[0]&M42,rec[0],rec[1])not in cand)
 nfix={k:[sum(v for i,v in vis[k].items()if fx[i]==e)for e in(3,4,6)]for k in vis}
 ext={k:3*nfix[k][0]+4*nfix[k][1]+6*nfix[k][2]for k in vis}
 return dict(fwd=nf,bwd=nb,keys=nk,acc_ok=acc%p==want,fraction_mismatches=bad_model,fractions_outside_bounds=bad_bound,fraction_errors=err,tfp_terms=tm,identity_mismatches=bad_leaf,fixups_3_4=nfix,fixup_extra=ext,record_decode_mismatches=dec,fixup_words_ok=words==[C0*sum(nfix[k])for k in("fwd","bwd")]and all(sum(nfix[k])<=w>>32<=sum(nfix[k])*C0>>32 for w,k in zip(words,("fwd","bwd"))),fwd_bwd_formulas=nf==cert.fwd_count(q2)+ext["fwd"]and nb==cert.bwd_count(q2)+ext["bwd"],records_ok=port.blocks==snd.msgs and machine.state1_ok(mem,snd),records=len(snd.seq),key_paths=meter,in_W=sum(win),in_D=sum(d and not w for w,d in zip(win,dset)),wrong_index=wrong,keys_not_canonical=canon,keys_formula=nk==cert.KEY_FIX+cert.KEY_INC*q2+meter+snd.extra,final_R_times_rho_is_one=final_R*rho%p==1,sp_kept=sp_kept,reg_trace_ok=all(any(lo<=a<hi for rl,_,lo,hi,_ in cert.REGISTRY if rl=='G')for a in mem.touched))
def macro_counts(src,entry,names):
 macros,cur,lines={},None,[ln.split("//")[0].strip()for ln in src.split("\n")]
 for s in lines:
  if s.startswith(".macro"):cur=s.split()[1].rstrip(",");macros[cur]=[]
  elif s.startswith(".endm"):cur=None
  elif cur is not None:macros[cur].append(s)
 out=dict.fromkeys(names,0)
 def walk(seq):
  for s in seq:
   if not s:continue
   w=s.split()[0].rstrip(",")
   if w in out:out[w]+=1
   if w in macros:walk(macros[w])
 start=lines.index(entry+":")
 walk(lines[start+1:lines.index("ret",start)])
 return out
def sqrt_mod_p(v):
 s_=pow(v,(p+1)//4,p)
 return s_ if s_*s_%p==v%p else None
def degenerate_batch(prog,labels,Xs,Ys):
 rprog,rlab=a64.a64_load(open(os.path.join(HERE,"kernels","route.s")).read())
 q2,nY=32,4;sb=EB_END-ENTRY*q2;Ys4=Ys*nY
 LS=([(AXYX0,XT_END),(BT0,1<<64),(YR0,YR0+512*nY),(sb,EB_END)],[(sb,EB_END)])
 mem=_TraceLC1Dict(rprog,rlab)
 popbat(rprog,rlab,mem,Xs,Ys4,7)
 for a in range(YR0-64,YR0,8):mem[a]=0
 for yi in range(nY):put_yrec(mem,YR0+512*yi,Ys4[yi])
 for l in range(q2):x_words(mem,sb+ENTRY*l,Xs,l)
 regs=dict(x0=q2,x1=YR0+64,x3=sb,x24=7,x25=0,x26=0,x27=0,sp=A_SP,x30=SNT)
 out={};a64.a64_run(prog,labels,"fwd",regs,mem,out_regs=out,ls64=LS)
 acc=rgi(out,("x24","x25","x26","x27"))
 regs=dict(x3=EB_END,x1=YR0-32,sp=YR0+512*(nY-1),x30=SNT)
 regs.update(ars(pow(acc,p-2,p),("x0","x2","x14","x15")))
 a64.a64_run(prog,labels,"bwd",regs,mem,ls64=LS)
 good=0
 for l in range(q2):
  X=Xs[l%MX]
  if(Ys[0][1]**2-X[1]**2)%p:
   good+=sum(-1<=dh<=8 and-1<=dl<=9 for dh,dl in frac_errs(mem,sb+ENTRY*l,X,Ys[0]))
 return acc,good
def edge_batch(prog,labels,Xs,Ys,ymem,cert,rg,kinds,need):
 over={}
 for l in range(32):
  for s in(0,1):
   while True:
    x=edge_leaf(rg,kinds[(2*l+s)%len(kinds)])
    if len(over)>=4 or(in_W(x)if need=="W"else in_D(x)and not in_W(x)):break
   over[(l,s)]=(x,rg.rr(3),rg)
 return leaf_batch(prog,labels,Xs,Ys,32,ymem,cert,override=over)
def tbld_production_tables(leaf_src,rng):
 rsrc=open(os.path.join(HERE,"kernels","route.s")).read()
 rprog,rlab=a64.a64_load(rsrc)
 mrp,mrl=a64.a64_load(rsrc.replace("\ttst x6,#1\n\tcsel x18,x3,x2,ne\n\tlsr x6,x6,#1","\tcmp x8,#2\n\tcsel x18,x3,x2,eq"))
 hprog,hlab=a64.a64_load(leaf_src+sites.SITE_KERNELS+sites.T_KERNELS)
 PSN,DSC,sb=0xDEADBEEFCAFEF00D,0x300000,EB_END-ENTRY*MX
 LS_H=([(AXYX0,XT_END),(BT0,1<<64),(YR0,YR0+512),(sb,EB_END)],[(sb,EB_END)])
 wval6=lambda ts:sum(t<<sh for t,sh in zip(ts,W6S))
 def psn(m,b,nw):
  for i in range(nw):m[b+8*i]=PSN
 def bp60(m,tb,c1,a1,c2,a2,rp_=rprog,rl_=rlab):
  tot=0
  for k,sh in enumerate(W6S):
   psn(m,tb+D6O[k],64);dk=DSC+128*k
   put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,a1 if k==0 else 0)
   put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,a2 if k==0 else 0)
   tot+=a64.a64_run(rp_,rl_,"tbld",dict(x0=tb+D6O[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k in(2,5)),x30=SNT),m)
  return tot
 def run_h(ent,regs,m):
  o={};a64.a64_run(hprog,hlab,ent,dict(regs,sp=YR0,x30=SNT),m,out_regs=o,ls64=LS_H);return o
 A1,Bv,cd=rng.rr(1,p//2),rng.rr(1,p),rng.rr(1,1<<63);A2=A1+cd;inv_e,ep=rng.rr(1,p),sb+576
 mc={};put(mc,YR0+288,Bv);put(mc,ep+224,A2);mc[ep+304]=cd
 oc=run_h("t_c1b",dict(x3=ep,**ars(inv_e,("x25","x26","x27","x28"))),mc)
 w1=inv_e*(Bv-A2)%p;w2=(w1+cd*inv_e)%p
 bday_c1b_ok=(get(mc,ep-64)%p,rgi(oc,("x8","x9","x10","x11"))%p)==(w1,w2)
 ma,cax,cyx,yYv,bs0={},rng.rr(1,p),rng.rr(1,p),rng.rr(1,p),ax_bases(0)
 inv_p1=3+(1<<43)+(1<<86);isq=inv_p1*inv_p1;mam={};na=0;q=bs0[5]+544;ma[q]=mam[q]=PSN
 for k,sh in enumerate(W6S):
  dk=DSC+128*k;put(ma,dk,cax*pow(2,sh,p)%p);put(ma,dk+32,0);put(ma,dk+64,cyx*pow(2,sh,p)%p);put(ma,dk+96,0)
  for a_ in range(dk,dk+128,8):mam[a_]=ma[a_]
  for rp_,rl_,m_ in((rprog,rlab,ma),(mrp,mrl,mam)):
   n_=a64.a64_run(rp_,rl_,"tbld",dict(x0=bs0[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5)),x30=SNT),m_)
   if m_ is ma:na+=n_
 for m_ in(ma,mam):
  for k,b_ in enumerate(bs0):m_[ep+256+8*k]=b_
  for k in range(4):put(m_,YR0+YR_YROW+32*k,yYv*pow(2,64*k,p)%p)
 oa=run_h("t_p1",dict(x3=ep,**ars(inv_p1,("x8","x9","x10","x11"))),ma)
 Mv,Nv=(inv_p1+cax*isq)%p,cyx*isq%p*yYv%p
 u2_got=oa["x28"]|oa["x29"]<<64|ma[YR0+496]<<128|oa["x1"]<<192
 axyx_ok=(na==6*26+20*(4*128+2*16)and u2_got%p==(Mv-Nv)%p and rgi(oa,("x0","x2","x5","x6"))%p==(Mv+Nv)%p and ma[q]==PSN and mam[q]!=PSN)
 mp,c1s,ad1,c2s,ad2,B2v,wb3={},rng.rr(1,p),rng.rr(1,p),rng.rr(1,p),rng.rr(1,p),rng.rr(1,p),wval6([5,3,1,6,7,1])
 np49,ob3=bp60(mp,XT0,c1s,ad1,c2s,ad2),offs6(wb3)
 put(mp,YR0+64,B2v);put(mp,YR0+96,B2v);mp[sb+416]=XT0
 for k in range(4):mp[YR0+96+8*k]=ob3[k]
 for k in range(2):mp[YR0+352+8*k]=ob3[4+k]
 ol=run_h("t_lkp",dict(x1=YR0+64,x3=sb),mp)
 mpm=dict(mp);bp60(mpm,XT0,c1s,ad1,c2s,ad2,mrp,mrl);olm=run_h("t_lkp",dict(x1=YR0+64,x3=sb),mpm)
 p49_ok=(np49==6*26+20*(4*16+2*4)and rgi(ol,("x16","x17","x18","x19"))%p==(c1s*wb3+ad1+B2v)%p and rgi(ol,("x20","x21","x22","x23"))%p==(c2s*wb3+ad2+B2v)%p and rgi(olm,("x16","x17","x18","x19"))%p!=(c1s*wb3+ad1+B2v)%p)
 me,wb,w3b2,b4v={},wval6([2,6,1,5,3,1]),wval6([3,4,1,2,6,1]),rng.rr(1,p)
 put(me,YR0+256,b4v);put(me,YR0+288,0);bl_w=limbs(wb,4);me[YR0+368],me[YR0+376],me[YR0+496],me[YR0+504]=bl_w[0],bl_w[1],bl_w[2],bl_w[3]
 for k,v in enumerate(offs6(wb3)[:4]):me[YR0+96+8*k]=v
 for k,v in enumerate(offs6(wb3)[4:]):me[YR0+352+8*k]=v
 for k,v in enumerate(offs6(wb)):me[YR0+384+8*k]=v
 for k,v in enumerate(offs6(w3b2)):me[YR0+440+8*k]=v
 eq_c=[(rng.rr(1,p),rng.rr(1,p),rng.rr(1,1<<63),rng.rr(1,1<<63))for _ in range(2)]
 fg_c=[[(rng.rr(1,p),rng.rr(1,p)if mi==0 else 0,rng.rr(1,p),rng.rr(1,p)if mi==0 else 0)for mi in range(2)]for _ in range(4)]
 for q in range(2):
  tb=XT0+q*2560*TIB;bp60(me,tb,eq_c[q][0],0,eq_c[q][1],0)
  e_q=sb+ENTRY*(4*q);me[e_q+376],me[e_q+384],me[e_q+392]=tb,2*eq_c[q][2],2*eq_c[q][3]
 for pr in range(4):
  for mi in range(2):
   tb=XT0+(2+4*mi+pr)*2560*TIB;bp60(me,tb,*fg_c[pr][mi]);me[sb+ENTRY*(2*pr)+400+8*mi]=tb
 run_h("t_efg",dict(x3=EB_END-448,x0=1,x2=2,x14=3,x15=4),me)
 efg_ok=all(get(me,sb+ENTRY*(4*q+2*s_)+64)%p==(eq_c[q][s_]*wb+b4v)%p and get(me,sb+ENTRY*(4*q+2*s_)+96)%p==(eq_c[q][s_]*wb+b4v+2*eq_c[q][2+s_]*wb)%p for q in range(2)for s_ in(0,1))
 efg_ok&=all(get(me,sb+ENTRY*(2*pr)+128+64*mi)%p==(fg_c[pr][mi][0]*wi+fg_c[pr][mi][1])%p and get(me,sb+ENTRY*(2*pr)+160+64*mi)%p==(fg_c[pr][mi][2]*wi+fg_c[pr][mi][3])%p for pr in range(4)for mi,wi in enumerate((wb3,w3b2)))
 mb,mbm={},{};psn(mb,BT0,256);psn(mbm,BT0,256)
 for idx,wi in enumerate(SLOTW):
  for m_ in(mb,mbm):put(m_,DSC+64*idx,BETA*pow(2,43*wi,p)%p);put(m_,DSC+64*idx+32,0)
 nbt=a64.a64_run(rprog,rlab,"tbld",dict(x0=BT0,x1=DSC,x2=8,x3=2,x4=7,x5=256,x6=1<<5,x30=SNT),mb)
 a64.a64_run(rprog,rlab,"tbld",dict(x0=BT0,x1=DSC,x2=8,x3=2,x4=7,x5=256,x6=1<<4,x30=SNT),mbm)
 u2=sum(t<<(43*i)for i,t in enumerate([6,5,7,5,4,1]))
 v2=sum(t<<(43*i)for i,t in enumerate([4,7,6,6,7,1]))
 ob=run_h("t_blk",dict(**ars(u2,("x28","x29","x11","x1")),**ars(v2,("x0","x2","x5","x6"))),mb)
 obm=run_h("t_blk",dict(**ars(u2,("x28","x29","x11","x1")),**ars(v2,("x0","x2","x5","x6"))),mbm)
 bt_ok=(nbt==4+7*11+20*50 and rgi(ob,("x24","x25","x26","x27"))%p==BETA*u2%p and rgi(ob,("x20","x21","x22","x23"))%p==BETA*v2%p and rgi(obm,("x24","x25","x26","x27"))%p!=BETA*u2%p)
 ok=bday_c1b_ok and axyx_ok and p49_ok and efg_ok and bt_ok
 return dict(ok=ok,bday_c1b_ok=bday_c1b_ok,axyx_ok=axyx_ok,pair49_ok=p49_ok and efg_ok,bt_ok=bt_ok)
def part_a_kernels(rng,cert):
 res={}
 src=open(os.path.join(HERE,"kernels","leaf.s")).read()
 prog,labels=a64.a64_load(src)
 fprog,flabels=a64.a64_load(open(os.path.join(HERE,"kernels","secp_fe.s")).read())
 rprog,rlab=a64.a64_load(open(os.path.join(HERE,"kernels","route.s")).read())
 yg,ymem,ch0,ch1=y_generator(fprog,flabels,rng,cert)
 res["y_generator"]=yg
 Xs=make_bday_xs(rng)
 Ys=[chain_y(ch0,i)for i in range(12)]
 runs={}
 for q2 in(32,64):
  runs[q2]=leaf_batch(prog,labels,Xs,Ys,q2,ymem,cert,rprog=rprog,rlab=rlab,chain=ch0)
 ov_top={(0,0):(planted_x(rng,(1<<41)-4,0)[0],0,rng),(0,1):(planted_x(rng,(1<<41)-3,1)[0],1,rng)}
 c1_run=leaf_batch(prog,labels,Xs,[chain_y(ch1,i)for i in range(4)],32,ymem,cert,rprog=rprog,rlab=rlab,override=ov_top,chain=ch1,ybase=YR0+YR_STRIDE)
 res["yslot1_and_lc1_top_ok"]=(c1_run["records_ok"]and c1_run["keys_formula"]and c1_run["fwd_bwd_formulas"]and c1_run["sp_kept"]and c1_run["reg_trace_ok"])
 res["registry_trace_ok"]=all(d["reg_trace_ok"]for d in runs.values())and c1_run["reg_trace_ok"]
 res["runs"]={q:dict(fwd=d["fwd"],bwd=d["bwd"],keys=d["keys"],fixup_extra=d["fixup_extra"],records_ok=d["records_ok"])for q,d in runs.items()}
 net={q:{k:d[k]-d["fixup_extra"][k]for k in("fwd","bwd")}for q,d in runs.items()}
 res["fwd_bwd_formulas_hold"]=all(d["fwd_bwd_formulas"]and d["fixup_words_ok"]for d in runs.values())
 res["counts_ok"]=all(net[q]["fwd"]==cert.fwd_count(q)and net[q]["bwd"]==cert.bwd_count(q)for q in runs)
 res["fixups"]=dict(sites=len(fixup_extra(prog)),flags_dead_at_every_join=flags_dead_after_fixups(prog,labels))
 tbzs=[o for mn,o in prog if mn in("tbz","tbnz")]
 res["loop_ends_on_bit_34"]=[o for o in tbzs if o[0]=="x3"]==[["x3","#34","1b"]]and all(o[1]=="#11"for o in tbzs if o[0]=="x5")
 res["keys_formula_holds"]=all(d["keys_formula"]for d in runs.values())
 res["tfp_terms"]=dict(largest={k:max(d["tfp_terms"][k]for d in runs.values())for k in TERMS},bounds=TERMS,sums_exact=all(d["tfp_terms"]["sum_ok"]for d in runs.values()))
 res["fraction_errors"]=[[min(d["fraction_errors"][g][0]for d in runs.values()),max(d["fraction_errors"][g][1]for d in runs.values())]for g in(0,1)]
 res["all_correct"]=(yg["mismatches"]==0 and yg["centres_ok"]and yg["counts_ok"]and res["yslot1_and_lc1_top_ok"]and all(d["acc_ok"]and d["fraction_mismatches"]==0 and d["fractions_outside_bounds"]==0 and d["identity_mismatches"]==0 and d["record_decode_mismatches"]==0 and d["records_ok"]and d["keys_not_canonical"]==0 and d["wrong_index"]==0 and d["final_R_times_rho_is_one"]and d["sp_kept"]and d["reg_trace_ok"]for d in runs.values()))
 erng=random.Random(SEED+49979687)
 for name,kinds,need in(("window_batch",("Xm","diff"),"W"),("boundary_batch",("D",),"D")):
  w=edge_batch(prog,labels,Xs,Ys,ymem,cert,erng,kinds,need)
  res[name]=dict(leaves=6*32,in_W=w["in_W"],in_D_outside_W=w["in_D"],wrong_index=w["wrong_index"],keys_not_canonical=w["keys_not_canonical"],records=w["records"],records_ok=w["records_ok"],ok=(w[{"W":"in_W","D":"in_D"}[need]]>=4 and w["wrong_index"]==0 and w["keys_not_canonical"]==0 and w["records"]==6*32 and w["records_ok"]and w["keys_formula"]and w["fwd_bwd_formulas"]and w["fixup_words_ok"]))
 gp,gl=a64.a64_load(src.replace("\tb.lo 9f\n","\tcbnz xzr,9f\n").replace("\tb.hs 9f\n","\tcbnz xzr,9f\n"))
 g=leaf_batch(gp,gl,Xs,Ys,32,ymem,cert,rprog=rprog,rlab=rlab,chain=ch0)
 tot_quad=sum(cert.FIX_SITES_QUAD.values())
 ext_quad=sum(cert.FIX_EXTRA[m]*cert.FIX_SITES_QUAD[m]for m in cert.FIX_SITES_QUAD)
 res["forced_fixups"]=dict(incidences=32,fixups_3_4=g["fixups_3_4"],fixup_extra=g["fixup_extra"],fraction_mismatches=g["fraction_mismatches"],ok=(g["fraction_mismatches"]==0 and g["fractions_outside_bounds"]==0 and g["identity_mismatches"]==0 and g["fwd_bwd_formulas"]and g["fixup_words_ok"]and g["records_ok"]and sum(map(sum,g["fixups_3_4"].values()))==8*tot_quad and sum(g["fixup_extra"].values())==8*ext_quad))
 cnts,ok=set(),True
 for _ in range(3):
  a=rng.rr(1,p);mem={};put(mem,0x1000,a)
  cnts.add(a64.a64_run(fprog,flabels,"fe_inv",dict(x0=0x2000,x1=0x1000,sp=0x100000,x30=SNT),mem,limit=10**6))
  ok&=get(mem,0x2000)%p==pow(a,p-2,p)
 res["fe_inv"]=dict(counts=sorted(cnts),correct=ok and cnts=={cert.INV})
 names=("SHX","SHXT","FDX","F5","FS","FADDW","LKFOLD4")
 counts={k:macro_counts(src,k,names+("LKP6_PAIR","BWDPAIR"))for k in("fwd","bwd","keys")}
 copies_f,copies_b=2*counts["fwd"]["LKP6_PAIR"],counts["bwd"]["BWDPAIR"]
 pquad={}
 for n in names:pquad[n]=(counts["fwd"][n]//2+counts["bwd"][n])//4
 pquad["SHX"]+=pquad.pop("SHXT")
 pquad["BLK"]=pquad.pop("LKFOLD4")
 res["site_counts"]=dict(per_x_quad=pquad,pair_copies=dict(fwd=copies_f,bwd=copies_b))
 res["sites_match_certificate"]=(pquad==cert.FIX_SITES_QUAD and copies_f==16 and copies_b==8 and sum(counts["keys"].values())==0)
 srng=random.Random(SEED+15485863)
 res["site_checks"]=sites.site_checks(sys.modules[__name__],srng,fprog,flabels)
 res["tbld_production"]=tbld_production_tables(src,random.Random(SEED+27182818))
 phY=(BETA*Ys[0][0]%p,Ys[0][1])
 acc0,good=degenerate_batch(prog,labels,[Xs[0],phY,Xs[1],PROD.neg(Ys[0])]+Xs[2:6],[Ys[0]])
 res["d_zero_batch"]=dict(acc_zero=acc0%p==0,fractions_within_bounds=good,all_other_fractions_wrong=good==0)
 prog_r,_=a64.a64_load(open(os.path.join(HERE,"kernels","route.s")).read())
 tot_ins=len(prog)+len(prog_r)+len(fprog)
 res["program_instructions"]=dict(leaf=len(prog),route=len(prog_r),secp_fe=len(fprog),total=tot_ins,fits_2_18_bytes=4*tot_ins<2**18)
 return res
def lattice_point(tv,tu):
 c1=(2*(FA*tv+(FA-FB)*tu)+p)//(2*p)
 c2=(2*((FA-FB)*tv-FB*tu)+p)//(2*p)
 return c1*FB+c2*(FA-FB),c1*(FA-FB)-c2*FA
def planted_x(rg,d1,d2,v=None,j=None):
 while True:
  vc=v
  if vc is None:
   ha=(rg.getrandbits(23)<<41)|(rg.rr(16)<<20)|rg.rr(200,(1<<20)-200)
   la=(d1<<23)|rg.rr(200,(1<<23)-200)
   hb=(rg.getrandbits(14)<<50)|(rg.rr(2)<<20)|rg.rr(200,(1<<20)-200)
   lb=(d2<<23)|rg.rr(200,(1<<23)-200)
   vv,uu=lattice_point(((ha<<64|la)*p)>>128,((hb<<64|lb)*p)>>128)
   if not(0<vv<p and uu==BETA*vv%p and exk(vv)[0]==vv):continue
   rf=label(vv)
   if not(rf[0]==d1 and rf[1]==d2 and(rf[3]&((1<<50)-1))<16):continue
   vc=vv
  m=rg.rr(3)if j is None else j
  x=vc*pow(BETA,(3-m)%3,p)%p
  if exk(x)[1]==m and not in_W(x)and not in_D(x):return x,vc
def load_certificate():
 buf=io.StringIO()
 with contextlib.redirect_stdout(buf):import certificate as cert
 lines=[" ".join(ln.split())for ln in buf.getvalue().splitlines()]
 return cert,[ln for ln in lines if" ok"in ln or"FAIL"in ln or ln.startswith(("q1 ","E_total","declared:"))]
def part_c_dependent_pairs():
 from collections import Counter
 pc=43
 ptc=[(x,y)for x in range(pc)for y in range(pc)if(y*y-x**3-7)%pc==0]
 rc=len(ptc)+1
 C43=Curve(pc,rc)
 bc=next(b for b in range(2,pc)if pow(b,3,pc)==1)
 phc=lambda P,t:None if P is None else(P[0]*pow(bc,t,pc)%pc,P[1])
 G0=ptc[0]
 pta=[C43 .mul(k,G0)for k in range(rc)]
 same_orb=lambda P,Q:(P is None and Q is None)or(P is not None and Q is not None and pow(P[0],3,pc)==pow(Q[0],3,pc))
 X1,X2=pta[3],pta[7]
 c1_tot=sum(1 for Y1_ in pta for Y2_ in pta if same_orb(C43 .add(X1,Y1_),C43 .add(X2,Y2_))for X3_ in pta if same_orb(C43 .add(X3_,Y1_),C43 .add(X2,Y2_)))
 c1_ok=(c1_tot==36*(rc-1)+1 and c1_tot*((rc-1)**2)<36*(rc**3))
 c2b_max=0
 for t3 in range(3):
  for s3 in(1,-1):
   for t4 in range(3):
    for s4 in(1,-1):
     hits=0
     for A in pta:
      L1=C43 .add(X1,A)
      for Y2_ in pta:
       if not same_orb(L1,C43 .add(X2,Y2_)):continue
       L4=C43 .add(phc(X2,t4),Y2_ if s4>0 else C43 .neg(Y2_))
       for B in pta:
        if B==A:continue
        L3=C43 .add(phc(X1,t3),B if s3>0 else C43 .neg(B))
        if same_orb(L3,L4):hits+=1
     if hits>c2b_max:c2b_max=hits
 lam=next(l for l in range(2,rc)if pow(l,3,rc)==1 and phc(G0,1)==C43 .mul(l,G0))
 phi_ok=all(phc(pta[k],t)==pta[k*pow(lam,t,rc)%rc]for k in range(rc)for t in range(3))
 cls=[min(j for j in range(rc)if same_orb(pta[j],pta[k]))for k in range(rc)]
 def hist(t,s_):
  h=Counter()
  for x in range(rc):
   for y in range(rc):h[(cls[(x+y)%rc],cls[(pow(lam,t,rc)*x+s_*y)%rc])]+=1
  return h
 c2a_max=0
 for t3 in range(3):
  for s3 in(1,-1):
   for t4 in range(3):
    for s4 in(1,-1):
     if(t3,s3,t4,s4)==(0,1,0,1):continue
     h1,h2=hist(t3,s3),hist(t4,s4)
     c2a_max=max(c2a_max,sum(v*h2 .get(k_,0)for k_,v in h1 .items()))
 c2_ok=(phi_ok and c2b_max<=36*rc and c2b_max*((rc-1)**2)<36*(rc**2*(rc-1))and c2a_max*((rc-1)**2)<36*rc**4)
 deg4_max=max(max(Counter((D[0]-pow(bc,m,pc)*A[0])%pc for A in ptc if(D:=C43 .add(A,T))is not None).values())for T in ptc for m in range(3))
 ok=c1_ok and c2_ok and deg4_max<=4
 return dict(ok=ok,curve_p=pc,curve_r=rc,case1_hits=c1_tot,case2a_max=c2a_max,case2b_max=c2b_max,lemma_l_deg4_max_fiber=deg4_max)
def part_c_vandermonde_cs(cert):
 from collections import Counter
 def g16(a,b):
  r_=0
  for i in range(4):
   if b>>i&1:r_^=a<<i
  for i in range(6,3,-1):
   if r_>>i&1:r_^=0x13<<(i-4)
  return r_
 mg,rp,nc=16,4,64
 H=[[((c&3)^g16(c>>2&15,m))&3 for m in range(mg)]for c in range(nc)]
 pcnt=Counter(rp*m+H[c][m]for c in range(nc)for m in range(mg))
 sm=sum(1 for a in range(nc)for b in range(a+1,nc)for m in range(mg)if H[a][m]!=H[b][m])
 sb=(nc*nc//(rp*rp)-nc)//2
 want_g=sm+sb*(mg*(mg-1)//2)+(nc*(nc-1)//2-sb)*mg*(mg-1)
 got_g=sum((mg-s)+mg*(mg-1)-s*(s-1)//2 for a in range(nc)for b in range(a+1,nc)for s in[sum(H[a][m]==H[b][m]for m in range(mg))])
 pc=61;ptc=[(x,y)for x in range(pc)for y in range(pc)if(y*y-x**3-7)%pc==0]
 rc=len(ptc)+1;C61=Curve(pc,rc);bc=next(b for b in range(2,pc)if pow(b,3,pc)==1)
 phc=lambda P,t:(P[0]*pow(bc,t,pc)%pc,P[1])
 G0=ptc[0];lc=next(k for k in range(1,rc)if C61.mul(k,G0)==phc(G0,1))
 un=[1,rc-1,lc,rc-lc,lc*lc%rc,rc-lc*lc%rc]
 et={(a-b)%rc for a in un for b in un if(a-b)%rc}
 br={z*a*pow(b,-1,rc)%rc for z in un for a in et for b in et}|set(un)
 X0,X1=G0,C61.mul(next(k for k in range(2,rc)if k not in br),G0)
 ds=[(t1,s1,t2,s2)for t1 in range(3)for s1 in(1,-1)for t2 in range(3)for s2 in(1,-1)]
 def bvec(A,B,t1,s1,t2,s2):
  v=i=0
  for Y1 in ptc:
   L1=C61.add(phc(A,t1),Y1 if s1>0 else C61.neg(Y1))
   for Y2 in ptc:
    L2=C61.add(phc(B,t2),Y2 if s2>0 else C61.neg(Y2))
    if L1 and L2 and pow(L1[0],3,pc)==pow(L2[0],3,pc):v|=1<<i
    i+=1
  return v
 vs=[bvec(X0,X1,*d)for d in ds]
 mxg=max((vs[i]&vs[j]).bit_count()for i in range(len(ds))for j in range(i+1,len(ds)))
 swh=(vs[0]&bvec(X1,X0,0,1,0,1)).bit_count()
 mxf=max(max(Counter((pow(D[0],3,pc)-pow(Q[0],3,pc))%pc for Q in ptc if(D:=C61.add(T,C61.neg(Q)))).values())for T in ptc[:12])
 p0=[cert.pool_of_gen(m)for m in range(1024)]
 sh32=sum(cert.pool_of_gen(m)==cert.pool_of_gen((37<<15)+(cert.gf1024_mul(37,5)&31)*1024+m)for m in range(1024))
 ok=set(pcnt.values())=={nc//rp}and got_g==want_g and swh>36 and mxg<=30 and mxf<=12 and len(set(p0))==1024 and sh32==32
 return dict(ok=ok,vandermonde_exact_pairs=got_g,prod_cluster0_pools=len(set(p0)),prod_dv37_shared=sh32,old_swapped_hits=swh,swap_free_gamma_max_hits=mxg,deg12_morphism_max_fiber=mxf)
def part_d(cert):
 return dict(passed=not cert.FAILED,failed=cert.FAILED,checks_passed=len(CERT_LINES),tail_summary=CERT_LINES[-3:],claim="%d.%015d"%divmod(int(cert.CLAIM_E*10**15),10**15),q1=cert.Q1GD,per_leaf_declared=format(float(cert.DECLARED["per_leaf"]),".4f"),mu=float(cert.TFULL["mu"]),opened_success_lower_bound=float(cert.LO_RUN),gate_prob=float(cert.P_GATE),success_lower_bound=float(cert.GATED_LO),per_leaf=float(cert.PER_LEAF_OPEN),groups=cert.G_GROUPS)
def main():
 global CERT_LINES
 rng=random.Random(SEED)
 cert,CERT_LINES=load_certificate()
 OUT["input"]=cfg
 OUT["budget"]=dict(seconds=BUDGET,reserve=RESERVE,minimums=dict(class_rule_random=MIN_RANDOM,class_rule_near_tie=MIN_NEAR,class_labels=MIN_CLASSES))
 A_id=OUT["A_identities"]=part_a_identities(rng)
 id_ok=(A_id["beta_cube_root_of_1_mod_p"]and A_id["lambda_cube_root_of_1_mod_r"]and A_id["lambda_G_equals_beta_x"]and A_id["seven_nonresidue_no_x0_point"]and A_id["window_forms_nonzero"]and A_id["label_constants"]and A_id["class_rule"]["mismatches_outside_window_set"]==0 and A_id["class_labels"]["ok"])
 if DO_KERNELS:
  A_k=OUT["A_kernels"]=part_a_kernels(rng,cert)
  B_m=OUT["B_machine"]=machine.part_b(sys.modules[__name__],cert)
  OUT["counts_match_certificate"]=cm=(A_k["fwd_bwd_formulas_hold"]and A_k["counts_ok"]and A_k["fixups"]["flags_dead_at_every_join"]and A_k["loop_ends_on_bit_34"]and A_k["keys_formula_holds"]and A_k["tfp_terms"]["sums_exact"]and A_k["all_correct"]and A_k["window_batch"]["ok"]and A_k["boundary_batch"]["ok"]and A_k["forced_fixups"]["ok"]and A_k["fe_inv"]["correct"]and A_k["sites_match_certificate"]and A_k["site_checks"]["all_true"]and A_k["tbld_production"]["ok"]and A_k["d_zero_batch"]["acc_zero"]and A_k["d_zero_batch"]["all_other_fractions_wrong"]and A_k["program_instructions"]["fits_2_18_bytes"]and B_m["phase0"]["ok"]and B_m["formulas_ok"]and B_m["block_paths"]["ok"]and B_m["grouped_paths"]["ok"]and B_m["ls64"]["ok"]and B_m["va56"]["ok"]and B_m["encoding"]["ok"]and B_m["phase4_round"]["ok"]and B_m["y_ring"]["ok"]and B_m["gate"]["ok"]and B_m["wait_dispatch"]["ok"]and B_m["registry_trace"]["ok"])
 else:
  cm=True
 C_dep=OUT["C_dependent_pairs"]=part_c_dependent_pairs()
 C_van=OUT["C_vandermonde_cs"]=part_c_vandermonde_cs(cert)
 D_c=OUT["D_certificate"]=part_d(cert)
 OUT["all_checks_pass"]=bool(id_ok and cm and(not DO_KERNELS or A_k["registry_trace_ok"])and C_dep["ok"]and C_van["ok"]and D_c["passed"])
 OUT["seconds"]=round(time.time()-T0,2)
 json.dump(OUT,sys.stdout,separators=(", ",": "))
 print()
 if not OUT["all_checks_pass"]:
  sys.exit(1)
if __name__=="__main__":
 main()
