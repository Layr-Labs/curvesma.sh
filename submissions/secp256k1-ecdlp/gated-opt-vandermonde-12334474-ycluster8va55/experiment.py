import contextlib,functools,io,json,os,random,sys,time;random.Random.rr=random.Random.randrange;from collections import Counter,defaultdict;from fractions import Fraction;T0=time.time();cfg=json.load(sys.stdin);SEED=cfg.get("seed",1);BUDGET=cfg.get("seconds",12);DO_KERNELS=cfg.get("kernels",True);RESERVE=min(14.5,BUDGET-0.3);T_END=T0+BUDGET-RESERVE;MIN_RANDOM,MIN_NEAR,MIN_CLASSES=300,100,300
def package_root():
 here=os.path.dirname(os.path.abspath(__file__))
 for base in(here,os.path.join(here,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(base,"kernels","a64.py")):return base
 raise FileNotFoundError("kernels/a64.py not found")
HERE=package_root();sys.path.insert(0,os.path.join(HERE,"kernels"));import a64;a64 .a64_load=functools.lru_cache(maxsize=None)(a64 .a64_load);sys.path.insert(0,HERE);import sites,machine;M64=2**64-1;M41,M42=2**40-1,2**41-1;PART=105;M128=2**128;p=2**256-2**32-977;r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16);GX=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",16);GY=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8",16);BETA=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16);LAM=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16);C0=2**32+977;FB=int("e4437ed6010e88286f547fa90abfe4c3",16);FA=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16);FRAC_OFF=((64,80,128,144,192,208),(96,112,160,176,224,240));ENTRY,ACC,KEYS0,MX,Q2M=512,32,64,8,2**24-1;YR_B4,YR_B,YR_B2,YR_YROW,YR_NT=256,288,64,128,320;W6S=(0,43,86,128,171,214);W6M=tuple((1<<b)-1 for b in(43,43,42,43,43,42));KO=(1,2,3,4,6,7);BT=(1,BETA,BETA*BETA%p);MH_A,MH_B,ML=20,21,23;TERMS=dict(h_min= -1,h_max=9,l_min= -1,l_max=9);ERR_H=ERR_L=((-24,29),(-24,29),(-59,48));OUT,CERT_LINES={},[]
class Curve:
 def __init__(self,q,order):self.q,self.r=q,order
 def add(self,P,Q):
  q=self.q
  if P is None:return Q
  if Q is None:return P
  if P[0]==Q[0]:
   if(P[1]+Q[1])%q==0:return None
   o=3*P[0]*P[0]*pow(2*P[1],-1,q)%q
  else:o=(Q[1]-P[1])*pow(Q[0]-P[0],-1,q)%q
  x=(o*o-P[0]-Q[0])%q;return(x,(o*(P[0]-x)-P[1])%q)
 def neg(self,P):return None if P is None else(P[0],(-P[1])%self.q)
 def mul(self,k,P):
  k%=self.r
  if not k or P is None:return None
  q,X,Y,Z=self.q,P[0],P[1],1
  for b in bin(k)[3:]:
   if X is not None:
    if Y==0:X=None
    else:S,M=4*X*Y*Y%q,3*X*X%q;o=(M*M-2*S)%q;Y2=(M*(S-o)-8*pow(Y,4,q))%q;Z=2*Y*Z%q;X,Y=o,Y2
   if b=="1":
    if X is None:X,Y,Z=P[0],P[1],1
    else:
     Z2=Z*Z%q;U2,S2=P[0]*Z2%q,P[1]*Z2*Z%q
     if X==U2:
      if(Y+S2)%q==0 or Y==0:X=None
      else:S,M=4*X*Y*Y%q,3*X*X%q;o=(M*M-2*S)%q;Y2=(M*(S-o)-8*pow(Y,4,q))%q;Z=2*Y*Z%q;X,Y=o,Y2
     else:
      H,R=(U2-X)%q,(S2-Y)%q;H2=H*H%q;H3=H*H2%q;V=X*H2%q;X=(R*R-H3-2*V)%q;Y=(R*(V-X)-Y*H3)%q;Z=Z*H%q
  if X is None or Z%q==0:return None
  J=pow(Z,-1,q);I=J*J%q;return(X*I%q,Y*I*J%q)
PROD=Curve(p,r);GP=(GX,GY)
def limbs(v,n=4):return[(v>>(64*k))&M64 for k in range(n)]
def put(mem,a,v,n=4):
 for k,w in enumerate(limbs(v,n)):mem[a+8*k]=w
def get(mem,a,n=4):return sum(mem[a+8*k]<<(64*k)for k in range(n))
def rgi(regs,o):return sum(regs[n]<<(64*k)for k,n in enumerate(o))
def ars(v,o):return{n:w for n,w in zip(o,limbs(v,len(o)))}
def exk(x):v=[x,BETA*x%p,BETA*BETA%p*x%p];j=min(range(3),key=lambda t:v[t]);return v[j],j
def t128(v):return(v<<128)//p
def t64(v):return(v<<64)//p
def ws(*v):return sum(x>>64 for x in v)%2**64<<64|sum(x&M64 for x in v)%2**64
def wneg(x):return ws(x^(M128-1),1<<64|1)
def rfields(a,b):o,I,J,lb=a>>64,a&M64,b>>64,b&M64;d1,d2=I>>23,lb>>23;return d1,d2,(o>>41<<41)|d2,(J>>21<<21)|((o>>20)&((1<<21)-1))
def label(v):return rfields(t128(v),t128(BETA*v%p))
def class_label(x):return label(exk(x)[0])
def neg_t(v):return wneg(t128(v%p))
def corr(v):return neg_t(v),neg_t(BETA*v)
def fold(z):return(z%2**256+(z>>256)*C0)%2**256
def lkfold(z):z=z%2**256+(z>>256)*C0;return z-2**256+C0 if z>>256 else z
def fs_op(a,b):return(a-b if a>=b else a-b+p)%2**256
def faddw_op(a,b):return(a+b if a+b<2**256 else a+b-p)%2**256
def shr(c):
 R=[c]
 for _ in range(3):R.append((R[-1]%2**192*2**64+(R[-1]>>192)*C0)%2**256)
 return R
def shr_chain(c):a=limbs(c);U=((a[2]*C0&M64)+(a[1]*C0>>64))+(((a[3]*C0&M64)+(a[2]*C0>>64))<<64)+((a[0]+(a[3]*C0>>64))<<128)+(a[1]<<192);u=limbs(U%(1<<256));l1=(a[1]*C0)&M64;return[c,(a[1]<<128)|(a[2]<<192),U%(1<<256),l1|(u[0]<<64)|(u[1]<<128)|(u[2]<<192),(u[0]<<128)|(u[1]<<192)]
def tfp(R,M,I=None):
 m,w=limbs(M),[limbs(x)for x in R];H,L=[sum((m[k]*w[k+o][3]&M64)+(m[k]*w[k+o][2]>>64)for k in range(4))&M64 for o in(0,1)]
 if I is not None:z=(R[0]*M)%p;dh,J=(t64(z)-H+(1<<63))%(1<<64)-(1<<63),(t64(z<<64)%2**64-L+(1<<63))%(1<<64)-(1<<63);I["h_min"],I["h_max"]=min(I.get("h_min",99),dh),max(I.get("h_max",-99),dh);I["l_min"],I["l_max"]=min(I.get("l_min",99),J),max(I.get("l_max",-99),J);I["sum_ok"]=I.get("sum_ok",True)and-1<=dh<=9 and-1<=J<=9
 return H<<64|L
def fracs(mem,e,s):return[get(mem,e+o,2)for o in FRAC_OFF[s]]
def components(X,Y):
 a,b=X[0],Y[0];E,F,G=(b**4+2*a**3*b)%p,(a**4+2*a*b**3)%p,3*a*a*b*b%p;out=[]
 for w in(pow(Y[1]+X[1],-2,p),pow(Y[1]-X[1],-2,p)):out+=[E*w%p,F*w%p,G*w%p]
 return out
def werr(t,v):return[((t>>64*k)-(v>>64*k)+2**63)%2**64-2**63 for k in(1,0)]
def frac_errs(mem,e,X,Y):
 cs,out=components(X,Y),[]
 for s in(0,1):o=((cs[3*s]-Y[0])%p,(cs[3*s+1]-X[0])%p,cs[3*s+2]);out+=[werr(t128(BT[i%2]*o[i//2]%p),v)for i,v in enumerate(fracs(mem,e,s))]
 return out
def leaf_w(o,t):u=[[o[2*i],o[2*i+1],wneg(ws(o[2*i],o[2*i+1]))]for i in range(3)];return[ws(u[0][m],u[1][(m+t)%3],u[2][(m+2*t)%3])for m in(0,1)]
def record_fields(w):I,o=w[0]>>64,w[1]>>64;s=(I+o)&M64;j=(2 if o+s>M64 else 1)if I>=o else(2 if I+s>M64 else 0);n=ws(w[0],w[1])^(M128-1);return(j,)+rfields(*((w[0],w[1]),(w[1],n),(n,w[0]))[j])
def in_W(x):v=[x<<64,(BETA*x%p)<<64,(BETA*BETA%p*x%p)<<64];o=p<<64;lo,hi=ERR_H[0];return(any(w<hi*p or w>=o+lo*p for w in v[:2])or v[2]<(1-2*lo)*p or v[2]>o-(3+2*hi)*p or abs(v[0]-v[1])<(hi+1-lo)*p or any(3*lo*p<=w-v[2]<(3*hi+3)*p for w in v[:2]))
def in_D(x):v=(x,BETA*x%p,BETA*BETA%p*x%p);return any(u%2**mb<hi or u%2**mb>=2**mb+lo for w,I,o in zip(v,ERR_H,ERR_L)for u,mb,(lo,hi)in((t128(w)>>64,MH_A,I),(t128(w)&M64,ML,o)))
def dft_leaf(o,t,a,b):return(o[0]-b+pow(BETA,t,p)*(o[1]-a)+pow(BETA,2*t,p)*o[2])%p
def inverse_dft(L):o,I=pow(3,-1,p),BETA*BETA%p;return((L[0]+L[1]+L[2])*o%p,(L[0]+I*L[1]+BETA*L[2])*o%p,(L[0]+BETA*L[1]+I*L[2])*o%p)
def weak(o):return o.rr(p)if o.random()<0.85 else o.rr(C0)+p
def represent(x,t,a,b,o,J=None):
 w,c,k=weak(o)or 1,[o.rr(p)for _ in range(3)],o.rr(3);c[k]=0;v=(x+b+BT[t]*a)*pow(w,-1,p)-sum(BT[t*i%3]*c[i]for i in range(3));c[k]=v*pow(BT[t*k%3],-1,p)%p;I=BETA*w%p;I+=p if I<C0 and o.random()<0.5 else 0;out=[]
 for ci,y in zip(c,(b,a,None)):R,cw=shr_chain(ci),corr(y)if y else(0,0);out+=[ws(tfp(R,w,J),cw[0]),ws(tfp(R,I,J),cw[1])]
 return out,((c[0]*w-b)%p,(c[1]*w-a)%p,c[2]*w%p)
def edge_leaf(o,kind):
 m=o.rr(3)
 if kind=="D":b=o.rr(2);mb=(MH_A,ML)[b];u=o.getrandbits(64-mb)<<mb|o.choice((o.rr(70),2**mb-1-o.rr(70),o.rr(2**mb)));v= -(-(u<<64|o.getrandbits(64)if b==0 else o.getrandbits(64)<<64|u)*p//M128);return v*pow(BETA,-m,p)%p or 1
 off=o.rr(160<<192);w=off if o.random()<0.5 else p-1-off
 if kind=="Xm":return w*pow(BETA,-m,p)%p or 1
 return w*pow((1-BETA,1-BETA*BETA,BETA-BETA*BETA)[m]%p,-1,p)%p or 1
def deadline(I):o=time.time();return o+I*max(0.0,T_END-o)
def near_tie_x(o):I=o.choice([1,BETA,BETA*BETA,1-BETA,1-BETA*BETA,BETA-BETA*BETA]);dd=o.rr(-(1<<o.rr(150,250)),1<<o.rr(150,250));return dd*pow(I%p,-1,p)%p or 1
def index_check(x,o):t,a,b=o.rr(3),o.rr(1,p),o.rr(1,p);j=record_fields(leaf_w(represent(x,t,a,b,o)[0],t))[0];return in_W(x),j!=exk(x)[1]
def class_labels(o):
 out=dict(classes=0,edge_classes=0,members=[0,0,0],in_W=0,in_D=0,wrong_index=0,mismatches=0,errors_h=[[99,-99]for _ in range(3)],errors_l=[[99,-99]for _ in range(3)],tfp_terms={});dl=deadline(0.05)
 while out["classes"]<MIN_CLASSES or(out["classes"]<4*MIN_CLASSES and time.time()<dl):
  I=out["classes"]%3==2;x0=edge_leaf(o,("D","Xm","diff")[out["classes"]//3%3])if I else o.rr(1,p);v=exk(x0)[0];out["classes"]+=1;out["edge_classes"]+=I;a4=label(v)
  for m in range(3):
   x=v*pow(BETA,(3-m)%3,p)%p;t,a,b=o.rr(3),o.rr(1,p),o.rr(1,p);w=leaf_w(represent(x,t,a,b,o,out["tfp_terms"])[0],t)
   for k,c in enumerate((w[0],w[1],ws(w[0],w[1])^(M128-1))):eh,el=werr(t128(pow(BETA,k,p)*x%p),c);out["errors_h"][k]=[min(out["errors_h"][k][0],eh),max(out["errors_h"][k][1],eh)];out["errors_l"][k]=[min(out["errors_l"][k][0],el),max(out["errors_l"][k][1],el)]
   f=record_fields(w)
   if in_W(x):out["in_W"]+=1;continue
   if f[0]!=m:out["wrong_index"]+=1;continue
   out["members"][m]+=1
   if in_D(x):out["in_D"]+=1
   elif f[1:]!=a4:out["mismatches"]+=1
 tm=out["tfp_terms"];J=(out["mismatches"]==0 and out["wrong_index"]==0 and out["in_D"]>0 and out["in_W"]>0 and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors_h"],ERR_H))and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors_l"],ERR_L))and tm["sum_ok"]);return dict(ok=J,classes=out["classes"],in_W=out["in_W"],in_D=out["in_D"])
def part_a_identities(ac):
 I={};I["beta3"]=pow(BETA,3,p)==1 and BETA!=1;I["lam3"]=pow(LAM,3,r)==1 and LAM!=1;I["lamG_betax"]=PROD.mul(LAM,GP)==(BETA*GX%p,GY);I["no_x0"]=pow(7,(p-1)//2,p)==p-1;I["win_nz"]=all(c%p for c in(1-BETA,1-BETA*BETA,BETA-BETA*BETA,3));I["label_constants"]=(FA*FA-FA*FB+FB*FB==p and FB*BETA%p==FA-FB and(FA-FB)*BETA%p==p-FA);bad=drops=aa=a8=0;a6={"random":0,"near_tie":0};dl=deadline(0.10)
 while a8<26000 and(a8<MIN_RANDOM or time.time()<dl):x,kind=(ac.rr(1,p),"random")if a8<20000 else(near_tie_x(ac),"near_tie");a6[kind]+=1;a8+=1;a7,a9=index_check(x,ac);drops+=a7;aa+=a7 and a9;bad+=a9 and not a7
 extra=max(0,MIN_NEAR-a6["near_tie"]);ab=random.Random(SEED+7919)
 for _ in range(extra):a7,a9=index_check(near_tie_x(ab),ab);a6["near_tie"]+=1;drops+=a7;aa+=a7 and a9;bad+=a9 and not a7
 I["class_rule"]=dict(checked=a6,in_W=drops,mismatches_outside_window_set=bad);I["class_labels"]=class_labels(random.Random(SEED+104729));J=random.Random(SEED+918273);o,c132_ok=1<<64,True
 def _in_r(i_,j_,d_):ai=max(0,i_-d_);return 0<=i_<o and((ai<=j_<=min(o-1,o+d_-2*i_))or(max(ai,o-d_-i_)<=j_<=min(o-1,2*o+d_-2*i_)))
 for ak in range(600):ae,ad=J.choice([0,1,o-1,o//2,J.getrandbits(64)]),J.choice([0,1,o-1,J.getrandbits(64)]);aj=record_fields((ae<<64,ad<<64))[0];ag=(ae,ad,(~((ae+ad)%o))&M64);c132_ok&=_in_r(ag[aj],ag[(aj+1)%3],0);af=near_tie_x(J)if ak%3==0 else edge_leaf(J,("Xm","diff","D")[ak%3]);t_=J.rr(3);al=record_fields(leaf_w(represent(af,t_,J.rr(1,p),J.rr(1,p),J)[0],t_))[0];ah=pow(BETA,al,p)*af%p;a5,a4=t64(ah),t64(BETA*ah%p);c132_ok&=(_in_r(a5,a4,90)or a5<64 or a5>=o-64 or a4<64 or a4>=o-64 or min(abs(a5+a4-o),abs(a5+a4-2*o),abs(2*a5+a4-o),abs(2*a5+a4-2*o),abs(a5-a4))<=64)
 I["c132_ok"]=c132_ok;return I
SNT=0x5E5E5E5E5E5E5E5E;TIB=1<<40;A_SP=(34<<40)+0x7000000;EB_BASE=(1<<45)+(1<<40);EB_END=EB_BASE+(1<<34);YR0,YR_STRIDE,XT0,BT0=(1<<45)+(1<<31),1<<31,(1<<64)-32768*TIB,0xfff8<<48;XT_END=(1<<64)-(32768-12*2560)*TIB;AXYX0=4096*TIB;AXY_END=24576*TIB;AXYX_HI=20480*TIB;D6O=tuple(c<<48 for c in(0,2,8,4,6,9))
def ax_bases(j):return(AXYX0+((4*j+0)<<49),AXYX0+((4*j+1)<<49),AXYX_HI+((2*j+0)<<48),AXYX0+((4*j+2)<<49),AXYX0+((4*j+3)<<49),AXYX_HI+((2*j+1)<<48))
SLOTW=(0,1,3,2,4,5,3)
def offs6(v):return[D6O[k]+(((v>>o)&m)<<6)for k,(o,m)in enumerate(zip(W6S,W6M))]
def xcs(a):return 2*a**3%p,2*a%p,a*a%p,pow(a,4,p)
def y_words(b):return b,pow(b,3,p),3*b*b%p
def bt_addr(w,i):t=w>>43*i&(1<<43)-1;return BT0+(t<<8)+32*SLOTW.index(i)+(128*(w>>128&1)if i==3 else 0)
def blook(w):return lkfold(sum(BETA*(w>>43*i&(1<<43)-1)*pow(2,43*i,p)%p for i in range(6)))
def lookup6(o,w,I=0):return lkfold(sum(((o*((w>>J)&m)*pow(2,J,p)+(I if k==0 else 0))%p)for k,(J,m)in enumerate(zip(W6S,W6M))))
class _TraceLC1Dict(dict):
 def __init__(self,rprog,rlab):super().__init__();self._rp,self._rl,self.touched=rprog,rlab,set()
 def __getitem__(self,a):self.touched.add(a);return super().__getitem__(a)
 def __setitem__(self,a,v):self.touched.add(a);return super().__setitem__(a,v)
 def get(self,a,d=None):
  if a in self:self.touched.add(a)
  return super().get(a,d)
 def __missing__(self,a):
  if 0<=a<16*machine.H1 and a%8==0:d=a//8;b=machine.blk1(d);a64 .a64_run(self._rp,self._rl,"init3",dict(x0=1,x1=a,x2=b+32,x3=machine.B1,x30=SNT),self);self.touched.add(a);return super().__getitem__(a)
  for lo,hi,first,o in machine.LC1_TABLES[2:]:
   if lo<=a<hi:self.touched.add(a);super().__setitem__(a,0);return 0
  raise KeyError(hex(a))
def pop60(rprog,rlab,mem,tb,c1,a5,c2,a4,w,bs=None):
 o=(34<<40)+0x180000
 for k,(J,m)in enumerate(zip(W6S,W6M)):
  I=(w>>J)&m;ad=(bs[k]if bs else tb+D6O[k])+(I<<6)
  if ad in mem:continue
  put(mem,o,0);put(mem,o+32,(c1*I*pow(2,J,p)+(a5 if k==0 else 0))%p);put(mem,o+64,0);put(mem,o+96,(c2*I*pow(2,J,p)+(a4 if k==0 else 0))%p);a64 .a64_run(rprog,rlab,"tbld",dict(x0=ad,x1=o,x2=1,x3=1,x4=2,x5=64,x6=0,x30=SNT),mem)
def cube_root(v):c_=pow(v,(p+2)//9,p);return c_ if c_*c_*c_%p==v%p else None
def make_bday_xs(o):
 xs=[]
 for q in range(MX//2):
  while True:
   I=PROD.mul(o.rr(1,r),GP);A0=I[1]*I[1]%p;lo=A0&((1<<63)-1)
   if q==0 and lo<(1<<62):continue
   d=o.rr((1<<63)-lo,1<<63)if q==0 else o.rr(1,1<<63);J=A0+d
   if J>=p:continue
   y1,a1=sqrt_mod_p(J),cube_root((J-7)%p)
   if y1 is not None and a1 is not None:xs+=[I,(a1,y1)];break
 return xs
def popax(rprog,rlab,mem,j,X,w):pop60(rprog,rlab,mem,0,2*(X[1]*X[1]%p)%p,0,2*X[1]%p,0,w,ax_bases(j))
def popbt(rprog,rlab,mem,w):
 I=(34<<40)+0x183000
 for o in{((w>>(43*i))&((1<<43)-1))for i in range(6)}:
  J=BT0+(o<<8)
  if J in mem:continue
  for idx,wi in enumerate(SLOTW):tk=o&((1<<41)-1)if idx==5 else o;put(mem,I+64*idx,0);put(mem,I+64*idx+32,BETA*tk*pow(2,43*wi,p)%p)
  a64 .a64_run(rprog,rlab,"tbld",dict(x0=J,x1=I,x2=1,x3=1,x4=7,x5=256,x6=0,x30=SNT),mem);put(mem,J+224,0)
def popbat(rprog,rlab,mem,I,Ys,rho):
 a5=rho;o,J=[],[]
 for yi,Y in enumerate(Ys):
  b,yy=Y;B=yy*yy%p;b3,ae=pow(b,3,p),3*b*b%p
  for q in range(len(I)//4):a0,a2=I[4*q][0],I[4*q+2][0];pop60(rprog,rlab,mem,XT0+q*2560*TIB,2*pow(a0,3,p)%p,0,2*pow(a2,3,p)%p,0,b);A1,A2,A3,A4=I[4*q][1]**2%p,I[4*q+1][1]**2%p,I[4*q+2][1]**2%p,I[4*q+3][1]**2%p;sA,pA,sB,pB=(A1+A2)%p,(A1*A2)%p,(A3+A4)%p,(A3*A4)%p;pop60(rprog,rlab,mem,XT0+(10+q)*2560*TIB,(-sB)%p,(pB-7*sB)%p,(-sA)%p,(pA-7*sA)%p,b3)
  for pr in range(len(I)//2):(a1,y1),(a2,y2)=I[2*pr],I[2*pr+1];c10,c11,c12,ai=xcs(a1);c20,c21,ag,ah=xcs(a2);pop60(rprog,rlab,mem,XT0+(2+pr)*2560*TIB,c11,ai,c21,ah,b3);pop60(rprog,rlab,mem,XT0+(6+pr)*2560*TIB,c12,0,ag,0,ae);ac,ab=y1*y1%p,y2*y2%p;aa,af=(ac+ab)%p,(ac*ab)%p;a8=lkfold(sum((((-aa)%p*((b3>>sh)&m)*pow(2,sh,p)+((af-7*aa)%p if k==0 else 0))%p)for k,(sh,m)in enumerate(zip(W6S,W6M)))+(B*B%p));o.append(a5);J.append(a8);a5=fold(sum(((a8>>64*k)&M64)*shr(a5)[k]for k in range(4)))
 R=pow(a5%p,p-2,p)if a5%p else 0;prod=lambda m_,rr:sum(((m_>>64*k)&M64)*rr[k]for k in range(4))
 for jp in range(len(o)-1,-1,-1):
  yi=(2*jp)//len(I);b,yy=Ys[yi];B=yy*yy%p;pr=jp%(len(I)//2);(a1,y1),(a2,y2)=I[2*pr],I[2*pr+1];A1,A2=y1*y1%p,y2*y2%p;dA=(A2-A1)%p;RR=shr(R);inv_e,R=fold(prod(o[jp],RR)),fold(prod(J[jp],RR));a7=fold(prod(inv_e,shr(fs_op(B,A2))));a9=dA*inv_e;ad=faddw_op(a7,fold(a9%2**256+(a9>>256)*C0))
  for j,X,a4 in((2*pr+1,(a2,y2),ad),(2*pr,(a1,y1),a7)):a6=fold(a4*a4%2**256+(a4*a4>>256)*C0);popax(rprog,rlab,mem,j,X,a6);Mv=lkfold(a4+lookup6(2*(X[1]*X[1]%p)%p,a6));Pv=lookup6(2*X[1]%p,a6);Nv=fold(sum(((Pv>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in range(4)));popbt(rprog,rlab,mem,fs_op(Mv,Nv));popbt(rprog,rlab,mem,faddw_op(Mv,Nv))
def x_words(mem,e,o,l):
 j=l%MX;I=j>>1;X=o[j]
 for k,w in enumerate(corr(X[0])):put(mem,e+256+16*k,w,2)
 put(mem,e+288,X[1]*X[1]%p)
 for k,b_ in enumerate(ax_bases(j)):mem[e+320+8*k]=b_
 if j&1:mem[e+368]=(o[2*I+1][1]**2-o[2*I][1]**2)%p
 if j%4==0:q=j//4;J=(o[4*q+1][1]**2-o[4*q][1]**2)%p;d1=(o[4*q+3][1]**2-o[4*q+2][1]**2)%p;mem[e+376],mem[e+384],mem[e+392]=XT0+q*2560*TIB,2*J,2*d1
 if(j&1)==0:mem[e+400],mem[e+408]=XT0+(2+I)*2560*TIB,XT0+(6+I)*2560*TIB
 if j==0:mem[e+416],mem[e+424]=XT0+10*2560*TIB,XT0+11*2560*TIB
def put_yrec(mem,e,Y):
 b,y=Y;B=y*y%p;mem[e]=0;put(mem,e+YR_B4,pow(b,4,p));put(mem,e+YR_B,B);put(mem,e+YR_B2,B*B%p)
 for k in range(4):put(mem,e+YR_YROW+32*k,y*pow(2,64*k,p)%p)
 put(mem,e+YR_NT,neg_t(b),2);put(mem,e+YR_NT+16,neg_t(BETA*b),2);I,ob3,o=offs6(b),offs6(pow(b,3,p)),offs6(3*b*b%p)
 for k in range(4):mem[e+96+8*k]=ob3[k]
 for k in range(2):mem[e+352+8*k]=ob3[4+k]
 bl=limbs(b,4);mem[e+368],mem[e+376]=bl[0],bl[1]
 for k in range(6):mem[e+384+8*k]=I[k]
 for k in range(6):mem[e+440+8*k]=o[k]
 mem[e+496],mem[e+504]=bl[2],bl[3]
def bwd_model(Xs,Ys,mem,o,a7,R,tm=None,ybase=None):
 if ybase is None:ybase=YR0
 out,prod={},lambda m,rr:sum((m>>64*k&M64)*rr[k]for k in range(4))
 for jp in range(len(o)-1,-1,-1):
  l1,l2=2*jp,2*jp+1;(a1,a9),(a2,a8)=Xs[l1%MX],Xs[l2%MX];b,yy=Ys[l2//MX];B,A1,A2=yy*yy%p,a9*a9%p,a8*a8%p;RR=shr(R);inv_e,R=fold(prod(o[jp],RR)),fold(prod(a7[jp],RR));D2=fs_op(B,A2);dA=(A2-A1)%p;aa=fold(prod(inv_e,shr(D2)));ae=dA*inv_e;af=faddw_op(aa,fold(ae%2**256+(ae>>256)*C0))
  for l,(a,yx),J in((l2,(a2,a8),af),(l1,(a1,a9),aa)):
   a5,ai=Ys[l//MX];AX=yx*yx%p;yr=ybase+512*(l//MX);ad=[get(mem,yr+YR_YROW+32*k)for k in range(4)];I=fold(J*J%2**256+(J*J>>256)*C0);ac=lkfold(J+lookup6(2*AX%p,I));ag=lookup6(2*yx%p,I);ab=fold(prod(ag,ad));w=[fs_op(ac,ab),faddw_op(ac,ab)];bw=[blook(v)for v in w];c0,c1,c2,a4=xcs(a);b4=get(mem,yr+YR_B4);yw=y_words(a5)
   if(l%2)==0:a6=lkfold(sum(((c0*((a5>>sh)&m)*pow(2,sh,p))%p)for sh,m in zip(W6S,W6M))+b4)
   else:
    ah=lkfold(sum(((xcs(a1)[0]*((a5>>sh)&m)*pow(2,sh,p))%p)for sh,m in zip(W6S,W6M))+b4);te=(2*dA)*a5;a6=faddw_op(ah,fold(te%2**256+(te>>256)*C0))
   cs=[a6,lookup6(c1,yw[1],a4),lookup6(c2,yw[2])];fr={}
   for i,(c,yc)in enumerate(zip(cs,(a5,a,None))):
    Rc,cw=shr_chain(c),corr(yc)if yc else(0,0)
    for s in(0,1):fr[s,2*i]=ws(tfp(Rc,w[s],tm),cw[0]);fr[s,2*i+1]=ws(tfp(Rc,bw[s],tm),cw[1])
   out[l]=[fr[s,i]for s in(0,1)for i in range(6)]
 return out,R
def six_leaves(X,Y):
 out=[]
 for t in range(3):o=(pow(BETA,t,p)*X[0]%p,X[1]);out+=[PROD.add(o,Y)[0],PROD.add(o,PROD.neg(Y))[0]]
 return out
def fixup_extra(I):
 out,o={},-1
 for i,(J,_)in enumerate(I):
  if J in("b.lo","b.hs","cbnz"):o=i
  elif J=="stadd":out[i]=i-o
 return out
def flags_dead_after_fixups(J,labels):
 a4=("adc","adcs","sbc","sbcs","csel","cinc","cset","ccmp","ccmn","b.eq","b.ne","b.lo","b.hs")
 def target(t,pc):
  if t[: -1].isdigit()and t[-1:]in("f","b"):c_=labels[t[: -1]];return[q for q in c_ if q<pc][-1]if t[-1]=="b"else[q for q in c_ if q>=pc][0]
  return labels[t][0]
 for i in fixup_extra(J):
  I,seen=[i+1],set()
  while I:
   j=I.pop()
   if j in seen:continue
   seen.add(j);mn,o=J[j]
   if mn in a4:return False
   if mn in("adds","subs","cmp","cmn","tst","ret"):continue
   if mn in("b","tbz","tbnz","cbz","cbnz"):
    I.append(target(o[-1],j+1))
    if mn=="b":continue
   I.append(j+1)
 return True
CBS,CB0,YSP=2304,(34<<40)+0x2000000,(34<<40)+0x4000000
def chain_n(i):rho,q=divmod(i,33);return 33*rho+(0 if q==32 else(16-q//2)*(1-2*(q%2)))
@functools.lru_cache(maxsize=None)
def chain_y(o,i):return PROD.add(o[0],PROD.mul(chain_n(i)%r,o[1]))
def y_generator(a4,I,a5,cert):
 mem,chains={},[]
 for c in range(2):
  C,S=PROD.mul(a5 .rr(1,r),GP),PROD.mul(a5 .rr(1,r),GP);chains.append((C,S));cb=CB0+CBS*c
  for i,P in[(0,C)]+[(i,PROD.mul(i,S))for i in list(range(1,17))+[33]]:put(mem,cb+(1088 if i==33 else 64*i),P[0]);put(mem,cb+(1088 if i==33 else 64*i)+32,P[1])
  mem[cb+2240]=YR0+c*YR_STRIDE
 put(mem,CB0-128,1)
 for a in range(YSP-4096,YSP+1024,8):mem[a]=0
 put(mem,YSP+352,BETA);mem[YSP+384],mem[YSP+424]=C0,CB0;counts=[]
 for a7 in(2,1):mem[YSP+416]=CB0+CBS*a7;counts.append(a64 .a64_run(a4,I,"yround",dict(sp=YSP,x30=SNT),mem,limit=10**7))
 bad=J=0
 for c,rd in((0,0),(1,0),(0,1)):
  for q in range(33):i=33*rd+q;(b,y),e=chain_y(chains[c],i),YR0+c*YR_STRIDE+512*i;B=y*y%p;ok=(get(mem,e+YR_B4)%p==pow(b,4,p)and get(mem,e+YR_B)==B and get(mem,e+YR_B2)==B*B%p);ok&=all(get(mem,e+YR_YROW+32*k)==y*pow(2,64*k,p)%p for k in range(4));ok&=get(mem,e+YR_NT,2)==neg_t(b)and get(mem,e+YR_NT+16,2)==neg_t(BETA*b);ob,ob3,a6=offs6(b),offs6(pow(b,3,p)),offs6(3*b*b%p);ok&=([mem[e+96+8*k]for k in range(4)]+[mem[e+352+8*k]for k in range(2)])==ob3;ok&=[mem[e+384+8*k]for k in range(6)]==ob;o=limbs(b,4);ok&=[mem[e+440+8*k]for k in range(6)]==a6 and(mem[e+368],mem[e+376],mem[e+496],mem[e+504])==(o[0],o[1],o[2],o[3]);bad+=not ok;J+=1
 a9=[(get(mem,CB0+CBS*c),get(mem,CB0+CBS*c+32))for c in range(2)];ok=bad==0 and a9==[chain_y(chains[0],33*2+32),chain_y(chains[1],33+32)];a8=counts==[cert.YR_FIX+2*cert.YR_BLK,cert.YR_FIX+cert.YR_BLK];return dict(records=J,mismatches=bad,centres_ok=ok,counts=counts,counts_ok=a8),mem,chains[0],chains[1]
def leaf_batch(a9,labels,Xs,Ys,q2,ah,cert,rprog=None,rlab=None,W2=0x123456789,KN=5<<26,override=None,chain=None,rho=None,ybase=YR0):
 assert q2%32==0 and ybase%YR_STRIDE==0
 if rprog is None:rprog,rlab=a64 .a64_load(open(os.path.join(HERE,"kernels","route.s")).read())
 fx=fixup_extra(a9);ac={k:dict.fromkeys(fx,0)for k in("fwd","bwd")};nY,SP=q2//MX,A_SP;sb=EB_END-ENTRY*q2;LS=([(AXYX0,AXY_END),(XT0,1<<64),(ybase,ybase+512*nY),(sb,EB_END)],[(sb,EB_END)]);mem=_TraceLC1Dict(rprog,rlab);rho=rho or random.Random(q2).rr(1,p);popbat(rprog,rlab,mem,Xs,Ys[:nY],rho)
 if ah is not None:
  for a in range(ybase-64,ybase+512*nY,8):mem[a]=ah.get(a,0)
 else:
  for a in range(ybase-64,ybase,8):mem[a]=0
  for yi in range(nY):put_yrec(mem,ybase+512*yi,Ys[yi])
 for l in range(q2):x_words(mem,sb+ENTRY*l,Xs,l)
 regs=dict(x0=q2,x1=ybase+64,x3=sb,sp=SP,x30=SNT);regs.update(ars(rho,("x24","x25","x26","x27")));out={};nf=a64 .a64_run(a9,labels,"fwd",regs,mem,out_regs=out,visits=ac["fwd"],ls64=LS);words,aq=[out["x0"]],rgi(out,("x24","x25","x26","x27"));au=[get(mem,sb+ENTRY*(2*j)+ACC)for j in range(q2//2)];es=[get(mem,sb+ENTRY*(2*j))for j in range(q2//2)];R0=pow(aq,p-2,p);am=ybase+512*(nY-1);regs=dict(x3=EB_END,x1=ybase-32,sp=am,x30=SNT);regs.update(ars(R0,("x0","x2","x14","x15")));out={};nb=a64 .a64_run(a9,labels,"bwd",regs,mem,out_regs=out,visits=ac["bwd"],ls64=LS);words.append(out["x5"]);ae=rho
 for l in range(q2):ae=ae*(Ys[l//MX][1]**2-Xs[l%MX][1]**2)%p
 a8=rgi(out,("x0","x2","x14","x15"));sp_kept=out["sp"]==am;tm={};model,mR=bwd_model(Xs,Ys,mem,au,es,R0,tm,ybase);J,a4,a6,a5,leaves=int(mR!=a8),0,0,[[99,-99],[99,-99]],[]
 for l in range(q2):
  X,Y=Xs[l%MX],Ys[l//MX];e,cs=sb+ENTRY*l,components(X,Y);J+=sum(get(mem,e+o,2)!=w for o,w in zip(FRAC_OFF[0]+FRAC_OFF[1],model[l]))
  for i,(dh,dl)in enumerate(frac_errs(mem,e,X,Y)):g=i%6//4;a5[0]=[min(a5[0][0],dh),max(a5[0][1],dh)];a5[1]=[min(a5[1][0],dl),max(a5[1][1],dl)];a4+=not(g-2<=dh<=9 and g-2<=dl<=9)
  ao=six_leaves(X,Y);a6+=sum(dft_leaf(cs[3*s:3*s+3],t,X[0],Y[0])!=ao[2*t+s]for s in(0,1)for t in range(3));leaves.append(ao)
 for(l,s),(x,t,rg)in(override or{}).items():
  f6,c3=represent(x,t,Xs[l%MX][0],Ys[l//MX][0],rg)
  for o,v in zip(FRAC_OFF[s],f6):put(mem,sb+ENTRY*l+o,v,2)
  for t_ in range(3):leaves[l][2*t_+s]=(c3[0]+BT[t_]*c3[1]+BT[2*t_%3]*c3[2])%p
 ai=machine.Port1(mem);regs=dict(x2=W2,x3=sb+KEYS0,x26=KN,sp=SP,x30=SNT);out={};nk=a64 .a64_run(a9,labels,"keys",regs,mem,out_regs=out,dev=ai.dev,ls64=([(sb,EB_END)],[]));sp_kept&=out["sp"]==SP;a7=machine.Sender1(cert);meter,ad=machine.keys_reference(sys.modules[__name__],mem,sb,q2,W2,KN,a7,cert);aa=[leaves[l][2*t+s]for l in range(q2)for s in(0,1)for t in range(3)];ab=[in_W(x)for x in aa];af=[in_D(x)for x in aa];ak=sum(not w and f[0]!=exk(x)[1]for w,f,x in zip(ab,ad,aa));al=sum(f[1:]!=class_label(x)for w,d,f,x in zip(ab,af,ad,aa)if not w and not d);ap=0
 for k,(d,ag)in enumerate(a7 .seq if chain else()):l=(ag[3]>>2)&((1<<24)-1);s=k//3%2;Y,X=chain_y(chain,l//MX),Xs[l%MX];ar=[class_label(PROD.add((pow(BETA,t,p)*X[0]%p,X[1]),PROD.neg(Y)if s else Y)[0])for t in range(3)];ap+=l!=k//6 or Y!=Ys[l//MX]or(not ab[k]and not af[k]and(d,ag[0]&M42,ag[0],ag[1])not in ar)
 I={k:[sum(v for i,v in ac[k].items()if fx[i]==e)for e in(3,4,6)]for k in ac};aj={k:3*I[k][0]+4*I[k][1]+6*I[k][2]for k in ac};return dict(fwd=nf,bwd=nb,keys=nk,acc_ok=aq%p==ae,frac_mis=J,frac_out=a4,frac_err=a5,tfp_terms=tm,id_mis=a6,fixups_3_4=I,fixup_extra=aj,record_decode_mismatches=ap,fixup_words_ok=words==[C0*sum(I[k])for k in("fwd","bwd")]and all(sum(I[k])<=w>>32<=sum(I[k])*C0>>32 for w,k in zip(words,("fwd","bwd"))),fb_form=nf==cert.fwd_count(q2)+aj["fwd"]and nb==cert.bwd_count(q2)+aj["bwd"],records_ok=ai.blocks==a7 .msgs and machine.state1_ok(mem,a7),records=len(a7 .seq),key_paths=meter,in_W=sum(ab),in_D=sum(d and not w for w,d in zip(ab,af)),wrong_index=ak,keys_noncanon=al,keys_formula=nk==cert.KEY_FIX+cert.KEY_INC*q2+meter+a7 .extra,final_R_times_rho_is_one=a8*rho%p==1,sp_kept=sp_kept,reg_trace_ok=all(any(lo<=a<hi for rl,_,lo,hi,_ in cert.REGISTRY if rl=='G')for a in mem.touched))
def macro_counts(a5,entry,a4):
 o,J,I={},None,[ln.split("//")[0].strip()for ln in a5 .split("\n")]
 for s in I:
  if s.startswith(".macro"):J=s.split()[1].rstrip(",");o[J]=[]
  elif s.startswith(".endm"):J=None
  elif J is not None:o[J].append(s)
 out=dict.fromkeys(a4,0)
 def walk(seq):
  for s in seq:
   if not s:continue
   w=s.split()[0].rstrip(",")
   if w in out:out[w]+=1
   if w in o:walk(o[w])
 start=I.index(entry+":");walk(I[start+1:I.index("ret",start)]);return out
def sqrt_mod_p(v):s_=pow(v,(p+1)//4,p);return s_ if s_*s_%p==v%p else None
def degenerate_batch(I,labels,Xs,Ys):
 rprog,rlab=a64 .a64_load(open(os.path.join(HERE,"kernels","route.s")).read());q2,J=32,4;o=EB_END-ENTRY*q2;a5=Ys*J;LS=([(AXYX0,AXY_END),(XT0,1<<64),(YR0,YR0+512*J),(o,EB_END)],[(o,EB_END)]);mem=_TraceLC1Dict(rprog,rlab);popbat(rprog,rlab,mem,Xs,a5,7)
 for a in range(YR0-64,YR0,8):mem[a]=0
 for yi in range(J):put_yrec(mem,YR0+512*yi,a5[yi])
 for l in range(q2):x_words(mem,o+ENTRY*l,Xs,l)
 regs=dict(x0=q2,x1=YR0+64,x3=o,x24=7,x25=0,x26=0,x27=0,sp=A_SP,x30=SNT);out={};a64 .a64_run(I,labels,"fwd",regs,mem,out_regs=out,ls64=LS);a4=rgi(out,("x24","x25","x26","x27"));regs=dict(x3=EB_END,x1=YR0-32,sp=YR0+512*(J-1),x30=SNT);regs.update(ars(pow(a4,p-2,p),("x0","x2","x14","x15")));a64 .a64_run(I,labels,"bwd",regs,mem,ls64=LS);good=0
 for l in range(q2):
  X=Xs[l%MX]
  if(Ys[0][1]**2-X[1]**2)%p:good+=sum(-1<=dh<=8 and-1<=dl<=9 for dh,dl in frac_errs(mem,o+ENTRY*l,X,Ys[0]))
 return a4,good
def edge_batch(a4,labels,Xs,Ys,I,cert,J,o,a5):
 over={}
 for l in range(32):
  for s in(0,1):
   while True:
    x=edge_leaf(J,o[(2*l+s)%len(o)])
    if len(over)>=4 or(in_W(x)if a5=="W"else in_D(x)and not in_W(x)):break
   over[(l,s)]=(x,J.rr(3),J)
 return leaf_batch(a4,labels,Xs,Ys,32,I,cert,override=over)
def tbld_production_tables(ac,I):
 ak=open(os.path.join(HERE,"kernels","route.s")).read();rprog,rlab=a64 .a64_load(ak);ar,au=a64 .a64_load(ak.replace("tst x6,#1\ncsel x18,x3,x2,ne\nlsr x6,x6,#1","cmp x8,#2\ncsel x18,x3,x2,eq"));ap,az=a64 .a64_load(ac+sites.SITE_KERNELS+sites.T_KERNELS);ah,ab,sb=0xDEADBEEFCAFEF00D,0x300000,EB_END-ENTRY*MX;b0=([(AXYX0,AXY_END),(XT0,1<<64),(YR0,YR0+512),(sb,EB_END)],[(sb,EB_END)]);a6=lambda ts:sum(t<<sh for t,sh in zip(ts,W6S))
 def psn(m,b,nw):
  for i in range(nw):m[b+8*i]=ah
 def bp60(m,tb,c1,a1,c2,a2,b1=rprog,b2=rlab):
  aq=0
  for k,sh in enumerate(W6S):psn(m,tb+D6O[k],64);dk=ab+128*k;put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,a1 if k==0 else 0);put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,a2 if k==0 else 0);aq+=a64 .a64_run(b1,b2,"tbld",dict(x0=tb+D6O[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k in(2,5)),x30=SNT),m)
  return aq
 def run_h(b5,regs,m):o={};a64 .a64_run(ap,az,b5,dict(regs,sp=YR0,x30=SNT),m,out_regs=o,ls64=b0);return o
 A1,Bv,cd=I.rr(1,p//2),I.rr(1,p),I.rr(1,1<<63);A2=A1+cd;inv_e,ep=I.rr(1,p),sb+576;mc={};put(mc,YR0+288,Bv);put(mc,ep+224,A2);mc[ep+304]=cd;oc=run_h("t_c1b",dict(x3=ep,**ars(inv_e,("x25","x26","x27","x28"))),mc);w1=inv_e*(Bv-A2)%p;w2=(w1+cd*inv_e)%p;bday_c1b_ok=(get(mc,ep-64)%p,rgi(oc,("x8","x9","x10","x11"))%p)==(w1,w2);ma,b7,b6,am={},I.rr(1,p),I.rr(1,p),ax_bases(0);ae={};na=0;q=am[5]+544;ma[q]=ae[q]=ah
 for k,sh in enumerate(W6S):
  dk=ab+128*k;put(ma,dk,b7*pow(2,sh,p)%p);put(ma,dk+32,0);put(ma,dk+64,b6*pow(2,sh,p)%p);put(ma,dk+96,0)
  for a_ in range(dk,dk+128,8):ae[a_]=ma[a_]
  na+=a64 .a64_run(rprog,rlab,"tbld",dict(x0=am[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5)),x30=SNT),ma);a64 .a64_run(ar,au,"tbld",dict(x0=am[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5)),x30=SNT),ae)
 axyx_ok=(na==6*26+20*(4*128+2*16)and ma[q]==ah and ae[q]!=ah);mp,af,ag,aw,ax,ai,ad={},I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),a6([5,3,1,6,7,1]);ay,ob3=bp60(mp,XT0,af,ag,aw,ax),offs6(ad);put(mp,YR0+64,ai);put(mp,YR0+96,ai);mp[sb+416]=XT0
 for k,v in enumerate(ob3[:4]):mp[YR0+96+8*k]=v
 for k,v in enumerate(ob3[4:]):mp[YR0+352+8*k]=v
 ol=run_h("t_lkp",dict(x1=YR0+64,x3=sb),mp);av=dict(mp);bp60(av,XT0,af,ag,aw,ax,ar,au);b3=run_h("t_lkp",dict(x1=YR0+64,x3=sb),av);a8=(ay==6*26+20*(4*16+2*4)and rgi(ol,("x16","x17","x18","x19"))%p==(af*ad+ag+ai)%p and rgi(b3,("x16","x17","x18","x19"))%p!=(af*ad+ag+ai)%p);J,wb,aj,ao={},a6([2,6,1,5,3,1]),a6([3,4,1,2,6,1]),I.rr(1,p);put(J,YR0+256,ao);put(J,YR0+288,0);a7=limbs(wb,4);J[YR0+368],J[YR0+376],J[YR0+496],J[YR0+504]=a7[0],a7[1],a7[2],a7[3]
 for k,v in enumerate(ob3[:4]):J[YR0+96+8*k]=v
 for k,v in enumerate(ob3[4:]):J[YR0+352+8*k]=v
 for k,v in enumerate(offs6(wb)):J[YR0+384+8*k]=v
 for k,v in enumerate(offs6(aj)):J[YR0+440+8*k]=v
 a4=[(I.rr(1,p),I.rr(1,p),I.rr(1,1<<63),I.rr(1,1<<63))for _ in range(2)];a5=[[(I.rr(1,p),I.rr(1,p)if mi==0 else 0,I.rr(1,p),I.rr(1,p)if mi==0 else 0)for mi in range(2)]for _ in range(4)]
 for q in range(2):tb=XT0+q*2560*TIB;bp60(J,tb,a4[q][0],0,a4[q][1],0);al=sb+ENTRY*(4*q);J[al+376],J[al+384],J[al+392]=tb,2*a4[q][2],2*a4[q][3]
 for pr in range(4):
  for mi in range(2):tb=XT0+(2+4*mi+pr)*2560*TIB;bp60(J,tb,*a5[pr][mi]);J[sb+ENTRY*(2*pr)+400+8*mi]=tb
 run_h("t_efg",dict(x3=EB_END-448,x0=1,x2=2,x14=3,x15=4),J);aa=all(get(J,sb+ENTRY*(4*q+2*s_)+64)%p==(a4[q][s_]*wb+ao)%p and get(J,sb+ENTRY*(4*q+2*s_)+96)%p==(a4[q][s_]*wb+ao+2*a4[q][2+s_]*wb)%p for q in range(2)for s_ in(0,1))and all(get(J,sb+ENTRY*(2*pr)+128+64*mi)%p==(a5[pr][mi][0]*wi+a5[pr][mi][1])%p and get(J,sb+ENTRY*(2*pr)+160+64*mi)%p==(a5[pr][mi][2]*wi+a5[pr][mi][3])%p for pr in range(4)for mi,wi in enumerate((ad,aj)));a9={};psn(a9,BT0,256)
 for idx,wi in enumerate(SLOTW):put(a9,ab+64*idx,BETA*pow(2,43*wi,p)%p);put(a9,ab+64*idx+32,0)
 a64 .a64_run(rprog,rlab,"tbld",dict(x0=BT0,x1=ab,x2=8,x3=2,x4=7,x5=256,x6=1<<4,x30=SNT),a9);u2=sum(t<<(43*i)for i,t in enumerate([6,5,7,5,4,1]));b4=run_h("t_blk",dict(**ars(u2,("x28","x29","x11","x1")),**ars(1,("x0","x2","x5","x6"))),a9);bt_ok=(rgi(b4,("x24","x25","x26","x27"))%p!=BETA*u2%p);ok=bday_c1b_ok and axyx_ok and a8 and aa and bt_ok;return dict(ok=ok,bday_c1b_ok=bday_c1b_ok,axyx_ok=axyx_ok,pair49_ok=a8 and aa,bt_ok=bt_ok)
def part_a_kernels(a8,cert):
 I={};aj=open(os.path.join(HERE,"kernels","leaf.s")).read();a4,labels=a64 .a64_load(aj);a6,a7=a64 .a64_load(open(os.path.join(HERE,"kernels","secp_fe.s")).read());rprog,rlab=a64 .a64_load(open(os.path.join(HERE,"kernels","route.s")).read());yg,ad,ap,av=y_generator(a6,a7,a8,cert);I["y_generator"]=yg;Xs=make_bday_xs(a8);ac,ar=defaultdict(int),{};af,aw=34*(1<<40)+0x20000,34*(1<<40)+0x40000
 for q in range(4):put(ac,af+128*q,Xs[2*q][1]*Xs[2*q][1]%p);put(ac,af+128*q+32,Xs[2*q+1][1]*Xs[2*q+1][1]%p);put(ac,af+128*q+64,2*q+1,2);put(ac,af+128*q+80,2*q+2,2)
 a64 .a64_run(rprog,rlab,"psel",dict(x0=af,x1=4,x2=aw,sp=34*(1<<40)+0x18000,x30=SNT),ac,out_regs=ar);I["bday_xs_psel_ok"]=(ar["x0"]==1 and((Xs[0][1]*Xs[0][1]%p)>>63)!=((Xs[1][1]*Xs[1][1]%p)>>63)and[get(ac,aw+128*q+32*s)for q in range(4)for s in(0,1)]==[X[1]*X[1]%p for X in Xs]);Ys=[chain_y(ap,i)for i in range(12)];runs={}
 for q2 in(32,64):runs[q2]=leaf_batch(a4,labels,Xs,Ys,q2,ad,cert,rprog=rprog,rlab=rlab,chain=ap)
 am={(0,0):(planted_x(a8,(1<<41)-4,0)[0],0,a8),(0,1):(planted_x(a8,(1<<41)-3,1)[0],1,a8)};J=leaf_batch(a4,labels,Xs,[chain_y(av,i)for i in range(4)],32,ad,cert,rprog=rprog,rlab=rlab,override=am,chain=av,ybase=YR0+YR_STRIDE);I["yslot_lc1"]=(J["records_ok"]and J["keys_formula"]and J["fb_form"]and J["sp_kept"]and J["reg_trace_ok"]);I["registry_trace_ok"]=all(d["reg_trace_ok"]for d in runs.values())and J["reg_trace_ok"];I["runs"]={q:dict(fwd=d["fwd"],bwd=d["bwd"],keys=d["keys"],fixup_extra=d["fixup_extra"],records_ok=d["records_ok"])for q,d in runs.items()};au={q:{k:d[k]-d["fixup_extra"][k]for k in("fwd","bwd")}for q,d in runs.items()};I["fb_hold"]=all(d["fb_form"]and d["fixup_words_ok"]for d in runs.values());I["counts_ok"]=all(au[q]["fwd"]==cert.fwd_count(q)and au[q]["bwd"]==cert.bwd_count(q)for q in runs);I["fixups"]=dict(sites=len(fixup_extra(a4)),flags_dead_at_every_join=flags_dead_after_fixups(a4,labels));al=[o for mn,o in a4 if mn in("tbz","tbnz")];I["loop_b34"]=[o for o in al if o[0]=="x3"]==[["x3","#34","1b"]]and all(o[1]=="#11"for o in al if o[0]=="x5");I["keys_hold"]=all(d["keys_formula"]for d in runs.values());I["tfp_terms"]=dict(largest={k:max(d["tfp_terms"][k]for d in runs.values())for k in TERMS},bounds=TERMS,sums_exact=all(d["tfp_terms"]["sum_ok"]for d in runs.values()));I["frac_err"]=[[min(d["frac_err"][g][0]for d in runs.values()),max(d["frac_err"][g][1]for d in runs.values())]for g in(0,1)];I["all_correct"]=(yg["mismatches"]==0 and yg["centres_ok"]and yg["counts_ok"]and I["bday_xs_psel_ok"]and I["yslot_lc1"]and all(d["acc_ok"]and d["frac_mis"]==0 and d["frac_out"]==0 and d["id_mis"]==0 and d["record_decode_mismatches"]==0 and d["records_ok"]and d["keys_noncanon"]==0 and d["wrong_index"]==0 and d["final_R_times_rho_is_one"]and d["sp_kept"]and d["reg_trace_ok"]for d in runs.values()));az=random.Random(SEED+49979687)
 for ay,aq,ao in(("window_batch",("Xm","diff"),"W"),("boundary_batch",("D",),"D")):w=edge_batch(a4,labels,Xs,Ys,ad,cert,az,aq,ao);I[ay]=dict(in_W=w["in_W"],in_D=w["in_D"],ok=(w[{"W":"in_W","D":"in_D"}[ao]]>=4 and w["wrong_index"]==0 and w["keys_noncanon"]==0 and w["records"]==6*32 and w["records_ok"]and w["keys_formula"]and w["fb_form"]and w["fixup_words_ok"]))
 gp,gl=a64 .a64_load(aj.replace("b.lo 9f\n","cbnz xzr,9f\n").replace("b.hs 9f\n","cbnz xzr,9f\n"));g=leaf_batch(gp,gl,Xs,Ys,32,ad,cert,rprog=rprog,rlab=rlab,chain=ap);ag=sum(cert.FIX_SITES_QUAD.values());ah=sum(cert.FIX_EXTRA[m]*cert.FIX_SITES_QUAD[m]for m in cert.FIX_SITES_QUAD);I["forced_fixups"]=dict(fixup_extra=g["fixup_extra"],ok=(g["frac_mis"]==0 and g["frac_out"]==0 and g["id_mis"]==0 and g["fb_form"]and g["fixup_words_ok"]and g["records_ok"]and sum(map(sum,g["fixups_3_4"].values()))==8*ag and sum(g["fixup_extra"].values())==8*ah));ai,ok=set(),True
 for _ in range(3):a=a8 .rr(1,p);mem={};put(mem,0x1000,a);ai.add(a64 .a64_run(a6,a7,"fe_inv",dict(x0=0x2000,x1=0x1000,sp=0x100000,x30=SNT),mem,limit=10**6));ok&=get(mem,0x2000)%p==pow(a,p-2,p)
 I["fe_inv"]=dict(counts=sorted(ai),correct=ok and ai=={cert.INV});ak=("SHX","FDX","F5","FS","FADDW","LKFOLD4");counts={k:macro_counts(aj,k,ak+("LKP6_PAIR","BWDPAIR"))for k in("fwd","bwd","keys")};a9,aa=2*counts["fwd"]["LKP6_PAIR"],counts["bwd"]["BWDPAIR"];a5={}
 for n in ak:a5[n]=(counts["fwd"][n]//2+counts["bwd"][n])//4
 a5["BLK"]=a5 .pop("LKFOLD4");I["site_counts"]=dict(per_x_quad=a5,pair_copies=dict(fwd=a9,bwd=aa));I["sites_cert"]=(a5==cert.FIX_SITES_QUAD and a9==16 and aa==8 and sum(counts["keys"].values())==0);ax=random.Random(SEED+15485863);I["site_checks"]=sites.site_checks(sys.modules[__name__],ax,a6,a7);I["tbld_production"]=tbld_production_tables(aj,random.Random(SEED+27182818));b1=(BETA*Ys[0][0]%p,Ys[0][1]);b0,good=degenerate_batch(a4,labels,[Xs[0],b1,Xs[1],PROD.neg(Ys[0])]+Xs[2:6],[Ys[0]]);I["d_zero_batch"]=dict(acc_zero=b0%p==0,fractions_within_bounds=good,all_other_fractions_wrong=good==0);ae,_=a64 .a64_load(open(os.path.join(HERE,"kernels","route.s")).read());ab=len(a4)+len(ae)+len(a6);I["prog_insns"]=dict(leaf=len(a4),route=len(ae),secp_fe=len(a6),total=ab,fits_2_18_bytes=4*ab<2**18);return I
def lattice_point(o,I):c1=(2*(FA*o+(FA-FB)*I)+p)//(2*p);J=(2*((FA-FB)*o-FB*I)+p)//(2*p);return c1*FB+J*(FA-FB),c1*(FA-FB)-J*FA
def planted_x(o,d1,d2,v=None,j=None):
 while True:
  J=v
  if J is None:
   ha=(o.getrandbits(23)<<41)|(o.rr(16)<<20)|o.rr(200,(1<<20)-200);la=(d1<<23)|o.rr(200,(1<<23)-200);hb=(o.getrandbits(14)<<50)|(o.rr(2)<<20)|o.rr(200,(1<<20)-200);lb=(d2<<23)|o.rr(200,(1<<23)-200);I,uu=lattice_point(((ha<<64|la)*p)>>128,((hb<<64|lb)*p)>>128)
   if not(0<I<p and uu==BETA*I%p and exk(I)[0]==I):continue
   rf=label(I)
   if not(rf[0]==d1 and rf[1]==d2 and(rf[3]&((1<<50)-1))<16):continue
   J=I
  m=o.rr(3)if j is None else j;x=J*pow(BETA,(3-m)%3,p)%p
  if exk(x)[1]==m and not in_W(x)and not in_D(x):return x,J
def load_certificate():
 I=io.StringIO()
 with contextlib.redirect_stdout(I):import certificate as cert
 o=[" ".join(ln.split())for ln in I.getvalue().splitlines()];return cert,[ln for ln in o if" ok"in ln or"FAIL"in ln or ln.startswith(("q1 ","E_total","declared:"))]
def part_c_dependent_pairs():
 pc=43;a9=[(x,y)for x in range(pc)for y in range(pc)if(y*y-x**3-7)%pc==0];rc=len(a9)+1;J=Curve(pc,rc);bc=next(b for b in range(2,pc)if pow(b,3,pc)==1);aa=lambda P,t:None if P is None else(P[0]*pow(bc,t,pc)%pc,P[1]);G0=a9[0];a4=[J.mul(k,G0)for k in range(rc)];o=lambda P,Q:(P is None and Q is None)or(P is not None and Q is not None and pow(P[0],3,pc)==pow(Q[0],3,pc));X1,X2=a4[3],a4[7];a7=sum(1 for af in a4 for a8 in a4 if o(J.add(X1,af),J.add(X2,a8))for ag in a4 if o(J.add(ag,af),J.add(X2,a8)));ad=(a7==36*(rc-1)+1 and a7*((rc-1)**2)<36*(rc**3));I=0
 for t3 in range(3):
  for s3 in(1,-1):
   for t4 in range(3):
    for s4 in(1,-1):
     hits=0
     for A in a4:
      L1=J.add(X1,A)
      for a8 in a4:
       if not o(L1,J.add(X2,a8)):continue
       L4=J.add(aa(X2,t4),a8 if s4>0 else J.neg(a8))
       for B in a4:
        if B==A:continue
        L3=J.add(aa(X1,t3),B if s3>0 else J.neg(B))
        if o(L3,L4):hits+=1
     if hits>I:I=hits
 ae=next(l for l in range(2,rc)if pow(l,3,rc)==1 and aa(G0,1)==J.mul(l,G0));ab=all(aa(a4[k],t)==a4[k*pow(ae,t,rc)%rc]for k in range(rc)for t in range(3));cls=[min(j for j in range(rc)if o(a4[j],a4[k]))for k in range(rc)]
 def hist(t,s_):
  h=Counter()
  for x in range(rc):
   for y in range(rc):h[(cls[(x+y)%rc],cls[(pow(ae,t,rc)*x+s_*y)%rc])]+=1
  return h
 a5=0
 for t3 in range(3):
  for s3 in(1,-1):
   for t4 in range(3):
    for s4 in(1,-1):
     if(t3,s3,t4,s4)==(0,1,0,1):continue
     h1,h2=hist(t3,s3),hist(t4,s4);a5=max(a5,sum(v*h2 .get(k_,0)for k_,v in h1 .items()))
 ac=(ab and I<=36*rc and I*((rc-1)**2)<36*(rc**2*(rc-1))and a5*((rc-1)**2)<36*rc**4);a6=max(max(Counter((D[0]-pow(bc,m,pc)*A[0])%pc for A in a9 if(D:=J.add(A,T))is not None).values())for T in a9 for m in range(3));ok=ad and ac and a6<=4;return dict(ok=ok,curve_p=pc,curve_r=rc,case1_hits=a7,case2a_max=a5,case2b_max=I,lemma_l_deg4_max_fiber=a6)
def part_c_vandermonde_cs(cert):
 def g16(a,b):
  r_=0
  for i in range(4):
   if b>>i&1:r_^=a<<i
  for i in range(6,3,-1):
   if r_>>i&1:r_^=0x13<<(i-4)
  return r_
 mg,rp,J=16,4,64;H=[[((c&3)^g16(c>>2&15,m))&3 for m in range(mg)]for c in range(J)];ac=Counter(rp*m+H[c][m]for c in range(J)for m in range(mg));sm=sum(1 for a in range(J)for b in range(a+1,J)for m in range(mg)if H[a][m]!=H[b][m]);sb=(J*J//(rp*rp)-J)//2;a6=sm+sb*(mg*(mg-1)//2)+(J*(J-1)//2-sb)*mg*(mg-1);a5=sum((mg-s)+mg*(mg-1)-s*(s-1)//2 for a in range(J)for b in range(a+1,J)for s in[sum(H[a][m]==H[b][m]for m in range(mg))]);o=61;a4=[(x,y)for x in range(o)for y in range(o)if(y*y-x**3-7)%o==0];rc=len(a4)+1;I=Curve(o,rc);bc=next(b for b in range(2,o)if pow(b,3,o)==1);a8=lambda P,t:(P[0]*pow(bc,t,o)%o,P[1]);G0=a4[0];lc=next(k for k in range(1,rc)if I.mul(k,G0)==a8(G0,1));un=[1,rc-1,lc,rc-lc,lc*lc%rc,rc-lc*lc%rc];et={(a-b)%rc for a in un for b in un if(a-b)%rc};br={z*a*pow(b,-1,rc)%rc for z in un for a in et for b in et}|set(un);X0,X1=G0,I.mul(next(k for k in range(2,rc)if k not in br),G0);ds=[(t1,s1,t2,s2)for t1 in range(3)for s1 in(1,-1)for t2 in range(3)for s2 in(1,-1)]
 def bvec(A,B,t1,s1,t2,s2):
  v=i=0
  for Y1 in a4:
   L1=I.add(a8(A,t1),Y1 if s1>0 else I.neg(Y1))
   for Y2 in a4:
    L2=I.add(a8(B,t2),Y2 if s2>0 else I.neg(Y2))
    if L1 and L2 and pow(L1[0],3,o)==pow(L2[0],3,o):v|=1<<i
    i+=1
  return v
 vs=[bvec(X0,X1,*d)for d in ds];aa=max((vs[i]&vs[j]).bit_count()for i in range(len(ds))for j in range(i+1,len(ds)));a9=(vs[0]&bvec(X1,X0,0,1,0,1)).bit_count();ab=max(max(Counter((pow(D[0],3,o)-pow(Q[0],3,o))%o for Q in a4 if(D:=I.add(T,I.neg(Q)))).values())for T in a4[:12]);p0=[cert.pool_of_gen(m)for m in range(1024)];a7=sum(cert.pool_of_gen(m)==cert.pool_of_gen((37<<14)+(cert.gf1024_mul(37,5)&15)*1024+m)for m in range(1024));ok=set(ac.values())=={J//rp}and a5==a6 and a9>36 and aa<=30 and ab<=12 and len(set(p0))==1024 and a7==64;return dict(ok=ok,vandermonde_exact_pairs=a5,prod_cluster0_pools=len(set(p0)),prod_dv37_shared=a7,old_swapped_hits=a9,swap_free_gamma_max_hits=aa,deg12_morphism_max_fiber=ab)
def part_d(cert):return dict(passed=not cert.FAILED,failed=cert.FAILED,checks_passed=len(CERT_LINES),tail_summary=CERT_LINES[-3:],claim="%d.%015d"%divmod(int(cert.CLAIM_E*10**15),10**15),q1=cert.Q1GD,per_leaf_declared=format(float(cert.DECLARED["per_leaf"]),".4f"),mu=float(cert.TFULL["mu"]),opened_success_lower_bound=float(cert.LO_RUN),gate_prob=float(cert.P_GATE),success_lower_bound=float(cert.GATED_LO),per_leaf=float(cert.PER_LEAF_OPEN),groups=cert.G_GROUPS)
def main():
 global CERT_LINES;a7=random.Random(SEED);cert,CERT_LINES=load_certificate();OUT["input"]=cfg;OUT["budget"]=dict(seconds=BUDGET,reserve=RESERVE,minimums=dict(class_rule_random=MIN_RANDOM,class_rule_near_tie=MIN_NEAR,class_labels=MIN_CLASSES));I=OUT["A_identities"]=part_a_identities(a7);a4=(I["beta3"]and I["lam3"]and I["lamG_betax"]and I["no_x0"]and I["win_nz"]and I["label_constants"]and I["class_rule"]["mismatches_outside_window_set"]==0 and I["class_labels"]["ok"]and I["c132_ok"])
 if DO_KERNELS:o=OUT["A_kernels"]=part_a_kernels(a7,cert);J=OUT["B_machine"]=machine.part_b(sys.modules[__name__],cert);OUT["counts_match_certificate"]=cm=(o["fb_hold"]and o["counts_ok"]and o["fixups"]["flags_dead_at_every_join"]and o["loop_b34"]and o["keys_hold"]and o["tfp_terms"]["sums_exact"]and o["all_correct"]and o["window_batch"]["ok"]and o["boundary_batch"]["ok"]and o["forced_fixups"]["ok"]and o["fe_inv"]["correct"]and o["sites_cert"]and o["site_checks"]["all_true"]and o["tbld_production"]["ok"]and o["d_zero_batch"]["acc_zero"]and o["d_zero_batch"]["all_other_fractions_wrong"]and o["prog_insns"]["fits_2_18_bytes"]and J["phase0"]["ok"]and J["formulas_ok"]and J["block_paths"]["ok"]and J["grouped_paths"]["ok"]and J["ls64"]["ok"]and J["va55"]["ok"]and J["encoding"]["ok"]and J["phase4_round"]["ok"]and J["y_ring"]["ok"]and J["gate"]["ok"]and J["wait_dispatch"]["ok"]and J["registry_trace"]["ok"])
 else:
  cm=True
 a6=OUT["C_dependent_pairs"]=part_c_dependent_pairs();a5=OUT["C_vandermonde_cs"]=part_c_vandermonde_cs(cert);a8=OUT["D_certificate"]=part_d(cert);OUT["all_checks_pass"]=bool(a4 and cm and(not DO_KERNELS or o["registry_trace_ok"])and a6["ok"]and a5["ok"]and a8["passed"]);OUT["seconds"]=round(time.time()-T0,2);json.dump(OUT,sys.stdout,separators=(",",": "));print()
 if not OUT["all_checks_pass"]:sys.exit(1)
if __name__=="__main__":main()
