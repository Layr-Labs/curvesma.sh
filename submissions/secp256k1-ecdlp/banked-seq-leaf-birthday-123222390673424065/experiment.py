import contextlib,functools,io,json,os,random,sys,time;Rd=random.Random;Rd.rr,Rd.gb=Rd.randrange,Rd.getrandbits;from collections import defaultdict as dfd;from fractions import Fraction;T0=time.time();cfg=json.load(sys.stdin);SEED=cfg.get("seed",1);BUDGET=cfg.get("seconds",12);RESERVE=min(14.5,BUDGET-0.3);T_END=T0+BUDGET-RESERVE;MRA,MNE,MCL=300,100,300;rn,en,dc=range,enumerate,dict
def pkr():
 here=os.path.dirname(os.path.abspath(__file__))
 for base in(here,os.path.join(here,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(base,"kernels","a64.py")):return base
 raise FileNotFoundError("kernels/a64.py not found")
HERE=pkr();sys.path.insert(0,os.path.join(HERE,"kernels"));import a64;a64.a64_load=functools.lru_cache(maxsize=None)(a64.a64_load);sys.path.insert(0,HERE);import sites,machine;ME=sys.modules[__name__];M64=2**64-1;M41,M42=2**40-1,2**41-1;PART=105;M128=2**128;p=2**256-2**32-977;r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16);GX=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",16);GY=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8",16);BETA=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16);LAM=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16);C0=2**32+977;FB=int("e4437ed6010e88286f547fa90abfe4c3",16);FA=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16);FRAC_OFF=((64,80,128,144,32,48),(96,112,160,176,192,208));ENTRY,ACC,KEYS0,MX,Q2M=512,32,64,6,2**24-1;YEB,YR_B,YR_B2,YRW,YR_NT=368,288,64,128,320;YR_K=((256,264,272,280,448,456),(408,416,424,432,440),(464,472,480,488,496,504));W6S=(0,43,86,128,171,214);W6M=tuple((1<<b)-1 for b in(43,43,42,43,43,42));W6E=tuple(en(zip(W6S,W6M)));W6V=tuple((k,s,m,m)for k,(s,m)in W6E);W5V=tuple((k,43*k,2**43-1,2**(84 if k>3 else 43)-1)for k in rn(5));BT=(1,BETA,BETA*BETA%p);O,KO=(2,2,0),(1,0,2);MH_A,MH_B,ML=20,21,23;EH,EL=((-241,237),(-136,133),(-361,360)),((-399,401),(-234,228),(-592,588));OUT,CLS,CE={},[],None
def mm(a,v):a[:]=min(a[0],v),max(a[1],v)
def ksrc(n):return open(os.path.join(HERE,"kernels",n+".s")).read()
def kload(n):return a64.a64_load(ksrc(n))
PL={}
def pln(m):return getattr(m,"pl",0)+getattr(getattr(m,"near",m),"nv",0)
def arun(a,b,e,r,m,*x,**k):
 # one kernel call on memory m; PL counts, by entry, the loads that returned poison (a64.Bank.pl) or never-written
 # near-memory garbage (machine.PZ.nv) during the call
 c=pln(m);n=a64.a64_run(a,b,e,dc(r,x30=SNT),m,*x,**k);c=pln(m)-c
 if c:PL[e]=PL.get(e,0)+c
 return n
class Curve:
 def __init__(U,q,order):U.q,U.r=q,order
 def add(U,P,Q):
  q=U.q
  if P is None:return Q
  if Q is None:return P
  if P[0]==Q[0]:
   if(P[1]+Q[1])%q==0:return None
   o=3*P[0]*P[0]*pow(2*P[1],-1,q)%q
  else:o=(Q[1]-P[1])*pow(Q[0]-P[0],-1,q)%q
  x=(o*o-P[0]-Q[0])%q;return(x,(o*(P[0]-x)-P[1])%q)
 def neg(U,P):return None if P is None else(P[0],(-P[1])%U.q)
 def mul(U,k,P):
  k%=U.r
  if not k or P is None:return None
  q,X,Y,Z=U.q,P[0],P[1],1
  def dbl(X,Y,Z):S,M=4*X*Y*Y%q,3*X*X%q;o=(M*M-2*S)%q;return o,(M*(S-o)-8*pow(Y,4,q))%q,2*Y*Z%q
  for b in bin(k)[3:]:
   if X is not None:
    if Y==0:X=None
    else:X,Y,Z=dbl(X,Y,Z)
   if b=="1":
    if X is None:X,Y,Z=P[0],P[1],1
    else:
     Z2=Z*Z%q;U2,S2=P[0]*Z2%q,P[1]*Z2*Z%q
     if X==U2:
      if(Y+S2)%q==0 or Y==0:X=None
      else:X,Y,Z=dbl(X,Y,Z)
     else:
      H,R=(U2-X)%q,(S2-Y)%q;H2=H*H%q;H3=H*H2%q;V=X*H2%q;X=(R*R-H3-2*V)%q;Y=(R*(V-X)-Y*H3)%q;Z=Z*H%q
  if X is None or Z%q==0:return None
  J=pow(Z,-1,q);I=J*J%q;return(X*I%q,Y*I*J%q)
PROD=Curve(p,r);GP=(GX,GY)
def limbs(v,n=4):return[(v>>(64*k))&M64 for k in rn(n)]
def put(N,a,v,n=4):
 for k,w in en(limbs(v,n)):N[a+8*k]=w
def get(N,a,n=4):return sum(N[a+8*k]<<(64*k)for k in rn(n))
def rgi(regs,o):return sum(regs[n]<<(64*k)for k,n in en(o))
def ars(v,o):return{n:w for n,w in zip(o,limbs(v,len(o)))}
def exk(x):v=[x,BETA*x%p,BETA*BETA%p*x%p];j=min(rn(3),key=lambda t:v[t]);return v[j],j
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
 for _ in rn(3):R+=[(R[-1]%2**192*2**64+(R[-1]>>192)*C0)%2**256]
 return R
def rowsf(w):return(w+((w>>192)*C0>>64))%2**256
def frac_e(r):return r+((r*C0+2**255)>>256)
def blkf(w):return sum(frac_e(BETA*((w>>43*i)&(1<<43)-1)*pow(2,43*i,p)%p)for i in rn(6))%2**256
def ktfp(c1,c2,w,y1):
 y0=rowsf(w);h=lambda a,b:a*b>>64;m,n,d,a,f,g=limbs(c1%2**192,3),limbs(c2,2),limbs((c1-c2)%2**192,3),limbs(y0),limbs(y1),limbs((y0+y1)%2**256)
 def pr(x,y):t=len(x)>2;return(h(x[0],y[2])+x[0]*y[3]+h(x[1],y[1])+x[1]*y[2]+(x[2]*y[1]if t else 0))&M64,(h(x[0],y[1])+x[0]*y[2]+h(x[1],y[0])+x[1]*y[1]+(x[2]*y[0]if t else 0))&M64
 AH,AL=pr(m,a);BH,BL=pr(d,g);P=pr(n,f);return(AH+P[0]&M64)<<64|(AL+P[1]&M64),((BH-AH)&M64)<<64|(BL-AL)&M64
def ksum(z):return sum(x[0]for x in z),sum(x[1]for x in z)
def kents(cl,ce,v,cc):return[CE.eis_red(CE.EIS[CE.KCL[cl][k]],(ce*((v>>sh)&m)*pow(2,sh,p)+(cc if k==0 else 0))%p)for k,sh,_,m in FW("EFG".index(cl))]
def kcomp(a,b):return[ksum(kents("E",2*a**3,b,0)+[CE.eis_red(CE.EIS["EB"],pow(b,4,p))]),ksum(kents("F",-2*a,pow(b,3,p),-a**4)),ksum(kents("G",BT[2]*a*a,3*b*b%p,0))]
def fracs(N,e,s):return[get(N,e+o,2)for o in FRAC_OFF[s]]
def components(X,Y):
 a,b=X[0],Y[0];E,F,G=(b**4+2*a**3*b)%p,(a**4+2*a*b**3)%p,3*a*a*b*b%p;W=[]
 for w in(pow(Y[1]+X[1],-2,p),pow(Y[1]-X[1],-2,p)):W+=[E*w%p,F*w%p,G*w%p]
 return W
def werr(t,v):return[((t>>64*k)-(v>>64*k)+2**63)%2**64-2**63 for k in(1,0)]
def fre(N,e,X,Y):
 cs,W=components(X,Y),[]
 for s in(0,1):o=((cs[3*s]-Y[0])%p,-cs[3*s+1]%p,BT[2]*cs[3*s+2]%p);W+=[werr(t128(BT[i%2]*o[i//2]%p),v)for i,v in en(fracs(N,e,s))]
 return W
def leaf_w(o,t):b0,b1,a0,a1,g0,g1=o;n=wneg;return((ws(a0,a1,g0,n(b0),n(b1)),ws(b0,g1,n(a0))),(ws(n(a0),n(b0),n(b1),n(g0),n(g1)),ws(b0,g0,n(a1))),(ws(a0,a1,b0,n(g0),n(g1)),ws(b1,g0,n(a0))))[t]
def record_fields(w):I,o=w[0]>>64,w[1]>>64;s=(I+o)&M64;j=(2 if o+s>M64 else 1)if I>=o else(2 if I+s>M64 else 0);n=ws(w[0],w[1])^(M128-1);return(j,)+rfields(*((w[0],w[1]),(w[1],n),(n,w[0]))[j])
def in_W(x,o=0):x=BT[o]*x%p;v=[x<<64,(BETA*x%p)<<64,(BETA*BETA%p*x%p)<<64];o=p<<64;lo,hi=min(EH[0][0],EH[1][0]),max(EH[0][1],EH[1][1]);return(any(w<hi*p or w>=o+lo*p for w in v[:2])or v[2]<(1-2*lo)*p or v[2]>o-(3+2*hi)*p or abs(v[0]-v[1])<(hi+1-lo)*p or any(3*lo*p<=w-v[2]<(3*hi+3)*p for w in v[:2]))
def in_D(x,o=0):x=BT[o]*x%p;v=(x,BETA*x%p,BETA*BETA%p*x%p);return any(u%2**mb<hi or u%2**mb>=2**mb+lo for w,I,o in zip(v,EH,EL)for u,mb,(lo,hi)in((t128(w)>>64,MH_A,I),(t128(w)&M64,ML,o)))
def dft_leaf(o,t,a,b):return(o[0]-b+pow(BETA,t,p)*(o[1]-a)+pow(BETA,2*t,p)*o[2])%p
def inverse_dft(L):o,I=pow(3,-1,p),BETA*BETA%p;return((L[0]+L[1]+L[2])*o%p,(L[0]+I*L[1]+BETA*L[2])*o%p,(L[0]+BETA*L[1]+I*L[2])*o%p)
def weak(o):return o.rr(p)if o.random()<0.85 else o.rr(C0)+p
def rpr(x,t,a,b,o,raw=False):
 w,c,k=weak(o)or 1,[o.rr(p)for _ in rn(3)],o.rr(3);c[k]=0;v=(x+b+BT[t]*a)*pow(w,-1,p)-sum(BT[t*i%3]*c[i]for i in rn(3));c[k]=v*pow(BT[t*k%3],-1,p)%p;y1=blkf(w);W=[]
 for i,(ci,y)in en(zip((c[0],-c[1],BT[2]*c[2]),(b,None if raw else a,None))):
  cl=CE.KCL["EFG"[i]]+(("EB",)if i==0 else());ks=[o.rr(p)for _ in rn(len(cl)-1)];A,B=ktfp(*ksum([CE.eis_red(CE.EIS[q],u)for q,u in zip(cl,ks)]+[CE.eis_red(CE.EIS[cl[-1]],(ci-sum(ks))%p)]),w,y1);cw=(corr(y)if i<1 else(t128(y),t128(BETA*y%p)))if y else(0,0);W+=[ws(A,cw[0]),ws(B,cw[1])]
 return W,((c[0]*w-b)%p,(c[1]*w-a)%p,c[2]*w%p)
def edl(o,kind):
 m=o.rr(3)
 if kind=="D":b=o.rr(2);mb=(MH_A,ML)[b];u=o.gb(64-mb)<<mb|o.choice((o.rr(70),2**mb-1-o.rr(70),o.rr(2**mb)));v=-(-(u<<64|o.gb(64)if b==0 else o.gb(64)<<64|u)*p//M128);return v*pow(BETA,-m,p)%p or 1
 off=o.rr(160<<192);w=off if o.random()<0.5 else p-1-off
 if kind=="Xm":return w*pow(BETA,-m,p)%p or 1
 return w*pow((1-BETA,1-BETA*BETA,BETA-BETA*BETA)[m]%p,-1,p)%p or 1
def deadline(I):o=time.time();return o+I*max(0.0,T_END-o)
def ntx(o):I=o.choice([1,BETA,BETA*BETA,1-BETA,1-BETA*BETA,BETA-BETA*BETA]);dd=o.rr(-(1<<o.rr(150,250)),1<<o.rr(150,250));return dd*pow(I%p,-1,p)%p or 1
def ick(x,o):t,a,b=o.rr(3),o.rr(1,p),o.rr(1,p);j=record_fields(leaf_w(rpr(x,t,a,b,o)[0],t))[0];return in_W(x,O[t]),(j+O[t])%3!=exk(x)[1]
def class_labels(o):
 W=dc(classes=0,edge_classes=0,members=[0,0,0],in_W=0,in_D=0,wrong_index=0,mismatches=0,errors_h=[[999,-999]for _ in rn(3)],errors_l=[[999,-999]for _ in rn(3)]);dl=deadline(0.05)
 while W["classes"]<MCL or(W["classes"]<4*MCL and time.time()<dl):
  I=W["classes"]%3==2;x0=edl(o,("D","Xm","diff")[W["classes"]//3%3])if I else o.rr(1,p);v=exk(x0)[0];W["classes"]+=1;W["edge_classes"]+=I;a4=label(v)
  for m in rn(3):
   x=v*pow(BETA,(3-m)%3,p)%p;t,a,b=o.rr(3),o.rr(1,p),o.rr(1,p);w=leaf_w(rpr(x,t,a,b,o)[0],t)
   for k,c in en((w[0],w[1],ws(w[0],w[1])^(M128-1))):eh,el=werr(t128(BT[(O[t]+k)%3]*x%p),c);mm(W["errors_h"][k],eh);mm(W["errors_l"][k],el)
   f=record_fields(w)
   if in_W(x,O[t]):W["in_W"]+=1;continue
   if(f[0]+O[t])%3!=m:W["wrong_index"]+=1;continue
   W["members"][m]+=1
   if in_D(x,O[t]):W["in_D"]+=1
   elif f[1:]!=a4:W["mismatches"]+=1
 J=(W["mismatches"]==0 and W["wrong_index"]==0 and W["in_D"]>0 and W["in_W"]>0 and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(W["errors_h"],EH))and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(W["errors_l"],EL)));return dc(ok=J,classes=W["classes"],in_W=W["in_W"],in_D=W["in_D"],errors_h=W["errors_h"],errors_l=W["errors_l"])
def pai(ac):
 I={};I["beta3"]=pow(BETA,3,p)==1 and BETA!=1;I["lam3"]=pow(LAM,3,r)==1 and LAM!=1;I["lamG_betax"]=PROD.mul(LAM,GP)==(BETA*GX%p,GY);I["no_x0"]=pow(7,(p-1)//2,p)==p-1;I["win_nz"]=all(c%p for c in(1-BETA,1-BETA*BETA,BETA-BETA*BETA,3));I["label_constants"]=(FA*FA-FA*FB+FB*FB==p and FB*BETA%p==FA-FB and(FA-FB)*BETA%p==p-FA);bad=drops=aa=a8=0;a6={"random":0,"near_tie":0};dl=deadline(0.10)
 while a8<26000 and(a8<MRA or time.time()<dl):x,kind=(ac.rr(1,p),"random")if a8<20000 else(ntx(ac),"near_tie");a6[kind]+=1;a8+=1;a7,a9=ick(x,ac);drops+=a7;aa+=a7 and a9;bad+=a9 and not a7
 extra=max(0,MNE-a6["near_tie"]);ab=Rd(SEED+7919)
 for _ in rn(extra):a7,a9=ick(ntx(ab),ab);a6["near_tie"]+=1;drops+=a7;aa+=a7 and a9;bad+=a9 and not a7
 I["class_rule"]=dc(checked=a6,in_W=drops,mismatches_outside_window_set=bad);I["class_labels"]=class_labels(Rd(SEED+104729));J=Rd(SEED+918273);o,c132_ok=1<<64,True
 def _in_r(i_,j_,d_):ai=max(0,i_-d_);return 0<=i_<o and((ai<=j_<=min(o-1,o+d_-2*i_))or(max(ai,o-d_-i_)<=j_<=min(o-1,2*o+d_-2*i_)))
 for ak in rn(600):ae,ad=J.choice([0,1,o-1,o//2,J.gb(64)]),J.choice([0,1,o-1,J.gb(64)]);aj=record_fields((ae<<64,ad<<64))[0];ag=(ae,ad,(~((ae+ad)%o))&M64);c132_ok&=_in_r(ag[aj],ag[(aj+1)%3],0);af=ntx(J)if ak%3==0 else edl(J,("Xm","diff","D")[ak%3]);t_=J.rr(3);al=record_fields(leaf_w(rpr(af,t_,J.rr(1,p),J.rr(1,p),J)[0],t_))[0]+O[t_];ah=pow(BETA,al,p)*af%p;a5,a4=t64(ah),t64(BETA*ah%p);c132_ok&=(_in_r(a5,a4,90)or a5<64 or a5>=o-64 or a4<64 or a4>=o-64 or min(abs(a5+a4-o),abs(a5+a4-2*o),abs(2*a5+a4-o),abs(2*a5+a4-2*o),abs(a5-a4))<=64)
 I["c132_ok"]=c132_ok;return I
SNT=0x5E5E5E5E5E5E5E5E;TIB=1<<40;A_SP=machine.SP;EB_BASE=42<<40;EBE=EB_BASE+(1<<37);YR0,YRS,XT0,BT0=40<<40,1<<31,2**64-32768*TIB,0xfff8<<48;YRP=(32<<40)+YRS;V63=1<<63;XT_END=2**64-2048*TIB;AXYX0=4096*TIB;AXY_END=24576*TIB;AXYX_HI=16384*TIB;D6O=tuple(c<<48 for c in(0,2,8,4,6,9));D7O=tuple(c<<50 for c in(0,2,8,4,6,9));D5O,D5F=(tuple(c<<s for c in(0,2,4,6,8))for s in(48,50));XTK=tuple(XT0+k*10240*TIB for k in rn(3));XTP=(19456*TIB,22016*TIB)
def ax_bases(j):return(AXYX0+((4*j+0)<<49),AXYX0+((4*j+1)<<49),AXYX_HI+((2*j+0)<<48),AXYX0+((4*j+2)<<49),AXYX0+((4*j+3)<<49),AXYX_HI+((2*j+1)<<48))
SLOTW=(0,1,3,2,4,5,3,3)
def FW(i):return(W6V,W5V,W6V)[i]
def DW(i):return(D7O,D5F,D7O)[i]
def offs6(v):return[D5O[k]+(((v>>o)&m)<<6)for k,o,m,_ in W5V]
def offs7(v,i):return[XTK[i]+DW(i)[k]+(((v>>o)&m)<<8)for k,o,m,_ in FW(i)]
def lookup6(o,w,I=0):return lkfold(sum(((o*((w>>J)&m)*pow(2,J,p)+(I if k==0 else 0))%p)for k,J,_,m in W5V))
def pop60(rprog,rlab,N,tb,c1,a5,c2,a4,w,bs=None):
 o=(34<<40)+0x180000
 for k,J,mi,m in(W6V if bs else W5V):
  I=(w>>J)&m;ad=(bs[k]if bs else tb+D5O[k])+((I&mi)<<6)
  if ad in N:continue
  put(N,o,0);put(N,o+32,(c1*I*pow(2,J,p)+(a5 if k==0 else 0))%p);put(N,o+64,0);put(N,o+96,(c2*I*pow(2,J,p)+(a4 if k==0 else 0))%p);arun(rprog,rlab,"tbld",dc(x0=ad,x1=o,x2=1,x3=1,x4=2,x5=64,x6=0),N)
def cube_root(v):c_=pow(v,(p+2)//9,p);return c_ if c_*c_*c_%p==v%p else None
def mbx(o,ks=None):
 xs=[]
 for q in rn(4):
  while True:
   k_=o.rr(1,r);I=PROD.mul(k_,GP);A0=I[1]*I[1]%p;lo=A0&((1<<63)-1)
   if q==0 and lo<(1<<62):continue
   d=o.rr((1<<63)-lo,1<<63)if q==0 else o.rr(1,1<<63);J=A0+d
   if J>=p:continue
   y1,a1=sqrt_mod_p(J),cube_root((J-7)%p)
   if y1 is not None and a1 is not None:
    xs+=[I,(a1,y1)]
    if ks is not None:ks+=[k_]
    break
 return xs
def popax(rprog,rlab,N,j,X,w):pop60(rprog,rlab,N,0,2*(X[1]*X[1]%p)%p,0,2*X[1]%p,0,w,ax_bases(j))
def popbt(rprog,rlab,N,w):
 for o in{((w>>(43*i))&((1<<43)-1))for i in rn(6)}:
  J=BT0+(o<<8)
  if J in N:continue
  for idx,wi in en(SLOTW):put(N,J+32*idx,frac_e(BETA*o*pow(2,43*wi,p)%p))
def eiw(c):return limbs(c["W"][1])+limbs(p-c["U"][1])+limbs(c["kS"])+limbs(c["kT"])+limbs(pow(p,-1,2**192),3)+limbs(c["U"][0]%2**192,3)+limbs(c["W"][0]%2**192,3)+limbs(c["U"][1],2)+limbs(c["W"][1],2)+limbs(c["o"][0]%2**192,3)+limbs(c["o"][1],2)+limbs(c["C"],3)+limbs(c["X2"],2)+[len(c["c"])]+sum((limbs(l1%2**192,3)+limbs(l2%2**128,2)for l1,l2 in c["c"]),[])
KOFF=(0,40,80,128,168,208)
def popk(N,Xs,b):
 for i,(v,ce,cc)in en(((b,lambda a:2*a**3,lambda a:0),(pow(b,3,p),lambda a:-2*a,lambda a:-a**4),(3*b*b%p,lambda a:BT[2]*a*a,lambda a:0))):
  for k,sh,mi,m in FW(i):
   t=(v>>sh)&m;ad=XTK[i]+DW(i)[k]+((t&mi)<<8)
   if ad in N:continue
   for j in rn(6):a=Xs[j][0];c1,c2=CE.eis_red(CE.EIS[CE.KCL["EFG"[i]][k]],(ce(a)*t*pow(2,sh,p)+(cc(a)if k==0 else 0))%p);put(N,ad+KOFF[j],c1%2**192,3);put(N,ad+KOFF[j]+24,c2,2)
   N[ad+120]=N[ad+248]=0
def popbat(rprog,rlab,N,I,Ys,rho):
 a5=rho;o,J=[],[];A=[X[1]*X[1]%p for X in I[:MX]];sA=[(A[2*q]+A[2*q+1])%p for q in rn(3)];pA=[A[2*q]*A[2*q+1]%p for q in rn(3)]
 for yi,Y in en(Ys):
  b,yy=Y;B=yy*yy%p;b3=pow(b,3,p);popk(N,I,b);pop60(rprog,rlab,N,XTP[0],(-sA[1])%p,(pA[1]-7*sA[1])%p,(-sA[0])%p,(pA[0]-7*sA[0])%p,b3);pop60(rprog,rlab,N,XTP[1],0,0,(-sA[2])%p,(pA[2]-7*sA[2])%p,b3)
  for pr in rn(3):aa,af=sA[pr],pA[pr];a8=lkfold(sum((((-aa)%p*((b3>>sh)&m)*pow(2,sh,p)+((af-7*aa)%p if k==0 else 0))%p)for k,sh,_,m in W5V)+(B*B%p));o+=[a5];J+=[a8];a5=fold(sum(((a8>>64*k)&M64)*shr(a5)[k]for k in rn(4)))
 R=pow(a5%p,p-2,p)if a5%p else 0;prod=lambda m_,rr:sum(((m_>>64*k)&M64)*rr[k]for k in rn(4))
 for jp in rn(len(o)-1,-1,-1):
  yi=(2*jp)//MX;b,yy=Ys[yi];B=yy*yy%p;pr=jp%3;(a1,y1),(a2,y2)=I[2*pr],I[2*pr+1];A1,A2=y1*y1%p,y2*y2%p;dA=(A2-A1)%p;RR=shr(R);inv_e,R=fold(prod(o[jp],RR)),fold(prod(J[jp],RR));a7=fold(prod(inv_e,shr(fs_op(B,A2))));ad=lkfold(a7+dA*inv_e)
  for j,X,a4 in((2*pr+1,(a2,y2),ad),(2*pr,(a1,y1),a7)):a6=fold(a4*a4%2**256+(a4*a4>>256)*C0);popax(rprog,rlab,N,j,X,a6);Mv=lkfold(a4+lookup6(2*(X[1]*X[1]%p)%p,a6));Pv=lookup6(2*X[1]%p,a6);Nv=fold(sum(((Pv>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in rn(4)));popbt(rprog,rlab,N,fs_op(Mv,Nv));popbt(rprog,rlab,N,faddw_op(Mv,Nv))
def x_words(N,e,o,l):
 # ENTRY_X's 16 words; the 0x5A pair (kept for even j) stands for the ST64B's two spare registers
 j=l%MX;X=o[j]
 for k,w in en((t128(X[0]),t128(BETA*X[0]%p))):put(N,e+256+16*k,w,2)
 put(N,e+288,X[1]*X[1]%p)
 for k,b_ in en(ax_bases(j)+(0x5A5A5A5A5A5A5A5A,)*2):N[e+320+8*k]=b_
 if j&1:N[e+368]=(o[j][1]**2-o[j-1][1]**2)%p%2**64;N[e+376]=e+320
def pyr(N,e,Y):
 b,y=Y;B=y*y%p;N[e+56]=pow(b,3,p)>>215;c1,c2=CE.eis_red(CE.EIS["EB"],pow(b,4,p));put(N,e+YEB,c1%2**192,3);put(N,e+YEB+24,c2,2);put(N,e+YR_B,B);put(N,e+YR_B2,B*B%p)
 for k in rn(4):put(N,e+YRW+32*k,y*pow(2,64*k,p)%p)
 put(N,e+YR_NT,neg_t(b),2);put(N,e+YR_NT+16,neg_t(BETA*b),2);ob3=offs6(pow(b,3,p))
 for k in rn(4):N[e+96+8*k]=ob3[k]
 N[e+352],N[e+360]=ob3[4],0
 for i,v in en((b,pow(b,3,p),3*b*b%p)):
  for a,x in zip(YR_K[i],offs7(v,i)):N[e+a]=x
def bwd_model(Xs,Ys,N,o,a7,R,ybase=None):
 if ybase is None:ybase=YR0
 W,prod={},lambda m,rr:sum((m>>64*k&M64)*rr[k]for k in rn(4))
 for jp in rn(len(o)-1,-1,-1):
  l1,l2=2*jp,2*jp+1;(a1,a9),(a2,a8)=Xs[l1%MX],Xs[l2%MX];b,yy=Ys[l2//MX];B,A1,A2=yy*yy%p,a9*a9%p,a8*a8%p;RR=shr(R);inv_e,R=fold(prod(o[jp],RR)),fold(prod(a7[jp],RR));D2=fs_op(B,A2);dA=(A2-A1)%p;aa=fold(prod(inv_e,shr(D2)));af=lkfold(aa+dA*inv_e)
  for l,(a,yx),J in((l2,(a2,a8),af),(l1,(a1,a9),aa)):
   a5=Ys[l//MX][0];AX=yx*yx%p;yr=ybase+512*(l//MX);ad=[get(N,yr+YRW+32*k)for k in rn(4)];I=fold(J*J%2**256+(J*J>>256)*C0);ac=lkfold(J+lookup6(2*AX%p,I));ag=lookup6(2*yx%p,I);ab=fold(prod(ag,ad));cs=kcomp(a,a5);fr=[]
   for w in(fs_op(ac,ab),faddw_op(ac,ab)):
    y1=blkf(w)
    for i,c in en(cs):A_,B_=ktfp(*c,w,y1);cw=corr(a5)if i==0 else(0,0);fr+=[ws(A_,cw[0]),ws(B_,cw[1])]
   W[l]=fr
 return W,R
def six_leaves(X,Y):
 W=[]
 for t in rn(3):o=(pow(BETA,t,p)*X[0]%p,X[1]);W+=[PROD.add(o,Y)[0],PROD.add(o,PROD.neg(Y))[0]]
 return W
def fixup_extra(I):
 # a metered fixup ends with ldr t,[sp]; add t,t,x4; str t,[sp]; its extra count is the distance from its branch
 W,o={},-1
 for i,(J,a)in en(I):
  if J in("b.lo","b.hs","cbnz"):o=i
  elif J=="str"and a[1]=="[sp]"and I[i-1][0]=="add"and I[i-1][1][2:]==["x4"]:W[i]=i-o
 return W
def fdaf(J,Lb):
 a4=("adc","adcs","sbc","sbcs","csel","cinc","cset","ccmp","ccmn","b.eq","b.ne","b.lo","b.hs")
 def target(t,pc):
  if t[:-1].isdigit()and t[-1:]in("f","b"):c_=Lb[t[:-1]];return[q for q in c_ if q<pc][-1]if t[-1]=="b"else[q for q in c_ if q>=pc][0]
  return Lb[t][0]
 for i in fixup_extra(J):
  I,seen=[i+1],set()
  while I:
   j=I.pop()
   if j in seen:continue
   seen.add(j);mn,o=J[j]
   if mn in a4:return False
   if mn in("adds","subs","cmp","cmn","tst","ret"):continue
   if mn in("b","tbz","tbnz","cbz","cbnz"):
    I+=[target(o[-1],j+1)]
    if mn=="b":continue
   I+=[j+1]
 return True
CBS,CB0,YSP=2304,(36<<40)+0x2000000,(34<<40)+0x4000000
def chain_n(i):rho,q=divmod(i,33);return 33*rho+(0 if q==32 else(16-q//2)*(1-2*(q%2)))
@functools.lru_cache(maxsize=None)
def cy(o,i):return PROD.add(o[0],PROD.mul(chain_n(i)%r,o[1]))
def ygen(a4,I,a5,K,rt,gp):
 N,chains=machine.mach(ME,K),[];machine.gsel(ME,rt,K,N,gp)
 for c in rn(2):
  C,S=PROD.mul(a5.rr(1,r),GP),PROD.mul(a5.rr(1,r),GP);chains+=[(C,S)];cb=CB0+CBS*c
  for i,P in[(0,C)]+[(i,PROD.mul(i,S))for i in list(rn(1,17))+[33]]:put(N,cb+(1088 if i==33 else 64*i),P[0]);put(N,cb+(1088 if i==33 else 64*i)+32,P[1])
  N[cb+2240]=YRP+c*YRS
 put(N,CB0-128,1)
 put(N,YSP+352,BETA);N[YSP+384],N[YSP+424]=C0,CB0;counts=[];N[YSP+1016]=1
 for k,v in en(eiw(CE.EIS["EB"])):N[YSP+1024+8*k]=v
 for a7 in(2,1):N[YSP+416]=CB0+CBS*a7;counts+=[arun(a4,I,"yround",dc(sp=YSP),N,limit=10**7)]
 bad=J=0
 for c,rd in((0,0),(1,0),(0,1)):
  for q in rn(33):i=33*rd+q;(b,y),e=cy(chains[c],i),YRP+c*YRS+512*i;B=y*y%p;c1,c2=CE.eis_red(CE.EIS["EB"],pow(b,4,p));ok=(get(N,e+YEB,3)==c1%2**192 and get(N,e+YEB+24,2)==c2 and get(N,e+YR_B)==B and get(N,e+YR_B2)==B*B%p);ok&=all(get(N,e+YRW+32*k)==y*pow(2,64*k,p)%p for k in rn(4));ok&=get(N,e+YR_NT,2)==neg_t(b)and get(N,e+YR_NT+16,2)==neg_t(BETA*b);ob3=offs6(pow(b,3,p));ok&=([N[e+96+8*k]for k in rn(4)]+[N[e+352]])==ob3 and N[e+56]==pow(b,3,p)>>215;ok&=all([N[e+a]for a in YR_K[i_]]==offs7(v,i_)for i_,v in en((b,pow(b,3,p),3*b*b%p)));M_=machine.PZ();pyr(M_,0,(b,y));ok&=all(N[e+a]==v for a,v in M_.items());bad+=not ok;J+=1
 a9=[(get(N,CB0+CBS*c),get(N,CB0+CBS*c+32))for c in rn(2)];ok=bad==0 and a9==[cy(chains[0],33*2+32),cy(chains[1],33+32)];N[YSP+1016]=0;counts+=[arun(a4,I,"yround",dc(sp=YSP),N,limit=10**7)];a8=counts==[K.YR_FIX+2*K.YR_BLK,K.YR_FIX+K.YR_BLK,K.YR_FIX+K.YR_BLK_MINE];return dc(records=J,mismatches=bad,centres_ok=ok,counts=counts,counts_ok=a8),N,chains[0],chains[1]
def ysl(o,s):
 # a curve point (b, y) whose slice b^3 >> 215 is s, so that its Y record belongs to sub-rings of slice s
 while True:
  v=s<<215|o.gb(215);b,y=cube_root(v),sqrt_mod_p((v+7)%p)
  if v<p and b is not None and y is not None:return b,y
def ctx(K):return dc(m=machine.mach(ME,K),ht=machine.HT(K,K.SLOT_BITS),done=set(),cq=[],ysig=set())
def ysto(c,K,rt,gp,cls,KN,n,src=YRP):
 # production: ystore copies n producer records from W_S(gp) into W_Y(sigma) at each bank's cursor; a bank's six
 # header words, cursor through the second bwd save word, are written first (YQ_INIT). Returns the count check.
 N=c["m"];machine.gsel(ME,rt,K,N,gp);ss=[N[src+512*i+56]for i in rn(n)]
 for s_ in ss:
  if(s_,cls)not in c["ysig"]:c["ysig"].add((s_,cls));machine.sel(ME,rt,K,N,"Y",s_<<15|cls);arun(*rt,"tinit",dc(x0=6,x1=YR0-64,x2=YR0+64),N)
 machine.gsel(ME,rt,K,N,gp);o={};nn=arun(*rt,"ystore",dc(x15=src+64,x11=n,x26=cls,x24=KN,x25=1<<42,x14=K.SEL["Y"],x13=YR0-64,x12=K.YQ_LIM),N,out_regs=o,ls64=([(src,src+512*n)],[(YR0,K.YQ_LIM+448)]));return nn==K.YSTORE*n+K.YSTORE_FIX and o["x24"]==KN+(n<<42),ss
def lfb(a9,Lb,Xs,Ys,q2,K,rt,c=None,W2=0x123456789,KN=V63|5<<26,override=None,rho=None,pool=7,mi=3,cls=5,prod=True,stop=False,j0=0,lcap=None,u0=0,xw=1):
 # One chunk of virtual generator g=1024 sigma+mi on machine c. Production: the producer's records (model frames
 # pyr in W_S(producer); already stored when prod is false) go through ystore into W_Y(sigma). The member's visit
 # selects W_T(pool, slice), W_S(g) and W_Y(sigma) by bsel, then fwd, bwd and keys run on the frames where ystore
 # put them (frame j0 on); the two words bwd overwrites below the chunk base are saved and restored by cp2.
 assert q2%6==0;c=c or ctx(K);N=c["m"];fx=fixup_extra(a9);ac={k:dc.fromkeys(fx,0)for k in("fwd","bwd")};nY,SP=q2//MX,A_SP;sb=EBE-ENTRY*q2;sl=pow(Ys[0][0],3,p)>>215;assert all(pow(Y[0],3,p)>>215==sl for Y in Ys[:nY])and j0%K.CH==0;sig=sl<<15|cls;g=1024*sig+mi;ybase=YR0+512*j0;rho=rho or Rd(q2).rr(1,p);yok=True
 if prod:
  machine.gsel(ME,rt,K,N,g^1)
  for yi in rn(nY):pyr(N,YRP+512*yi,Ys[yi])
  yok,_=ysto(c,K,rt,g^1,cls,KN,nY)
 sc=[machine.sel(ME,rt,K,N,"T",pool,sl),machine.gsel(ME,rt,K,N,g),machine.sel(ME,rt,K,N,"Y",sig)];LS=([(AXYX0,AXY_END),(XT0,1<<64),(ybase,ybase+512*nY),(sb,EBE)],[(sb,EBE)]);popbat(rt[0],rt[1],N,Xs,Ys[:nY],rho)
 for l in rn(q2*xw):x_words(N,sb+ENTRY*l,Xs,l)
 regs=dc(x0=q2,x1=ybase+64,x3=sb,sp=SP);regs.update(ars(rho,("x24","x25","x26","x27")));W={};nf=arun(a9,Lb,"fwd",regs,N,out_regs=W,visits=ac["fwd"],ls64=LS);words,aq=[W["x0"]],rgi(W,("x24","x25","x26","x27"));au=[get(N,sb+ENTRY*(2*j)+ACC)for j in rn(q2//2)];es=[get(N,sb+ENTRY*(2*j))for j in rn(q2//2)];R0=pow(aq,p-2,p);am=ybase+512*(nY-1);sv=SP+64;sc+=[arun(*rt,"cp2",dc(x1=ybase-32,x2=sv),N)];regs=dc(x3=EBE,x1=ybase-32,sp=am);regs.update(ars(R0,("x0","x2","x14","x15")));W={};nb=arun(a9,Lb,"bwd",regs,N,out_regs=W,visits=ac["bwd"],ls64=LS);words+=[W["x5"]];sc+=[arun(*rt,"cp2",dc(x1=sv,x2=ybase-32),N)];ae=rho
 if stop:return aq,N,sb
 for l in rn(q2):ae=ae*(Ys[l//MX][1]**2-Xs[l%MX][1]**2)%p
 a8=rgi(W,("x0","x2","x14","x15"));sp_kept=W["sp"]==am;model,mR=bwd_model(Xs,Ys,N,au,es,R0,ybase);J,a4,a6,a5,leaves=int(mR!=a8),0,0,[[999,-999],[999,-999]],[]
 for l in rn(q2):
  X,Y=Xs[l%MX],Ys[l//MX];e,cs=sb+ENTRY*l,components(X,Y);J+=sum(get(N,e+o,2)!=w for o,w in zip(FRAC_OFF[0]+FRAC_OFF[1],model[l]))
  for i,(dh,dl)in en(fre(N,e,X,Y)):cb=K.KB2[i%6//2][i%2];mm(a5[0],dh);mm(a5[1],dl);a4+=not(cb[0][0]<=dh<=cb[0][1]and cb[1][0]<=dl<=cb[1][1])
  ao=six_leaves(X,Y);a6+=sum(dft_leaf(cs[3*s:3*s+3],t,X[0],Y[0])!=ao[2*t+s]for s in(0,1)for t in rn(3));leaves+=[ao]
 for(l,s),(x,t,rg)in(override or{}).items():
  f6,c3=rpr(x,t,Xs[l%MX][0],Ys[l//MX][0],rg,True)
  for o,v in zip(FRAC_OFF[s],f6):put(N,sb+ENTRY*l+o,v,2)
  for t_ in rn(3):leaves[l][2*t_+s]=(c3[0]+BT[t_]*c3[1]+BT[2*t_%3]*c3[2])%p
 ht=c["ht"];i0,nc=len(c["cq"]),len(ht.cand);ad,kok=machine.krun(ME,K,rt,(a9,Lb),N,ht,c["done"],g,sb,q2,W2,c["cq"],lcap=lcap);kq=c["cq"][i0:];nk=kq[-1][1];R=machine.kref(ME,N,sb,q2,W2,K)[2];meter=kq[-1][2]-K.KEY_FIX-K.KEY_INC*q2;ft=[(t,leaves[l][2*t+s])for l in rn(q2)for s in(0,1)for t in KO];aa=[x for _,x in ft];ab=[in_W(x,O[t])for t,x in ft];af=[in_D(x,O[t])for t,x in ft];ak=sum(not w and(f[0]+O[t])%3!=exk(x)[1]for w,f,(t,x)in zip(ab,ad,ft));al=sum(f[1:]!=class_label(x)for w,d,f,x in zip(ab,af,ad,aa)if not w and not d);ap=0
 for k,(d,k2,k3,w2,w3)in en(R):l=6*(((w3>>42)&(2**21-1))-u0)+w2-W2;s=k//3%2;Y,X=Ys[l//MX],Xs[l%MX];ar=[class_label(PROD.add((pow(BETA,t,p)*X[0]%p,X[1]),PROD.neg(Y)if s else Y)[0])for t in rn(3)];ap+=l!=k//6 or w3>>63!=1 or(not ab[k]and not af[k]and(d,k2&M42,k2,k3)not in ar)
 I={k:[sum(v for i,v in ac[k].items()if fx[i]==e)for e in(5,6,8)]for k in ac};aj={k:5*I[k][0]+6*I[k][1]+8*I[k][2]for k in ac};return dc(fwd=nf,bwd=nb,keys=nk,acc_ok=aq%p==ae,frac_mis=J,frac_out=a4,frac_err=a5,id_mis=a6,fixups_5_6_8=I,fixup_extra=aj,record_decode_mismatches=ap,fixup_words_ok=words==[C0*sum(I[k])for k in("fwd","bwd")]and all(sum(I[k])<=w>>32<=sum(I[k])*C0>>32 for w,k in zip(words,("fwd","bwd"))),fb_form=nf==K.fwd_count(q2)+aj["fwd"]and nb==K.bwd_count(q2)+aj["bwd"],records_ok=kok and machine.htok(N,K,ht)and yok and all(x==K.BSEL for x in sc[:3])and sc[3:]==[3,3],records=len(R),cand=ht.cand[nc:],key_paths=meter,in_W=sum(ab),in_D=sum(d and not w for w,d in zip(ab,af)),wrong_index=ak,keys_noncanon=al,keys_formula=all(x==y for _,x,y in kq),final_R_times_rho_is_one=a8*rho%p==1,sp_kept=sp_kept,g=g,sig=sig)
def mcc(a5,entry,a4):
 o,J,I={},None,[ln.split("//")[0].strip()for ln in a5.split("\n")]
 for s in I:
  if s.startswith(".macro"):J=s.split()[1].rstrip(",");o[J]=[]
  elif s.startswith(".endm"):J=None
  elif J is not None:o[J].append(s)
 W=dc.fromkeys(a4,0)
 def walk(seq):
  for s in seq:
   if not s:continue
   w=s.split()[0].rstrip(",")
   if w in W:W[w]+=1
   if w in o:walk(o[w])
 start=I.index(entry+":");walk(I[start+1:I.index("ret",start)]);return W
def sqrt_mod_p(v):s_=pow(v,(p+1)//4,p);return s_ if s_*s_%p==v%p else None
def dgb(I,Lb,Xs,Ys,K,rt):
 a4,N,o=lfb(I,Lb,Xs,Ys*4,24,K,rt,rho=7,stop=True);good=0
 for l in rn(24):
  X=Xs[l%MX]
  if(Ys[0][1]**2-X[1]**2)%p:good+=sum(abs(dh)<4096 and abs(dl)<4096 for dh,dl in fre(N,o+ENTRY*l,X,Ys[0]))
 return a4,good
def edge_batch(a4,Lb,Xs,Ys,K,rt,J,o,a5):
 over={}
 for l in rn(24):
  for s in(0,1):
   t=J.rr(3)
   while True:
    x=edl(J,o[(2*l+s)%len(o)])
    if len(over)>=4 or(in_W(x,O[t])if a5=="W"else in_D(x,O[t])and not in_W(x,O[t])):break
   over[(l,s)]=(x,t,J)
 return lfb(a4,Lb,Xs,Ys,24,K,rt,override=over)
def tbpt(ac,I):
 ak=ksrc("bank");rprog,rlab=a64.a64_load(ak);ar,au=a64.a64_load(ak.replace("tst x6,#1\ncsel x18,x3,x2,ne\nlsr x6,x6,#1","cmp x8,#2\ncsel x18,x3,x2,eq"));ap,az=a64.a64_load(ac+sites.SITE_KERNELS+sites.T_KERNELS);fp,fl=kload("secp_fe");ah,ab,sb=0xDEADBEEFCAFEF00D,0x300000,EBE-ENTRY*MX;b0=([(AXYX0,AXY_END),(XT0,1<<64),(YR0,YR0+512),(sb,EBE)],[(sb,EBE)]);a6=lambda ts:sum(t<<43*k for k,t in en(ts))
 def psn(m,b,nw):
  for i in rn(nw):m[b+8*i]=ah
 def bp60(m,tb,c1,a1,c2,a2,b1=rprog,b2=rlab):
  aq=0
  for k,sh,_,_ in W5V:psn(m,tb+D5O[k],64);dk=ab+128*k;put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,a1 if k==0 else 0);put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,a2 if k==0 else 0);aq+=arun(b1,b2,"tbld",dc(x0=tb+D5O[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=0),m)
  return aq
 def rh(b5,regs,m):o={};arun(ap,az,b5,dc(regs,sp=YR0),m,out_regs=o,ls64=b0);return o
 A1,Bv,cd=I.rr(1,p//2),I.rr(1,p),I.rr(1,1<<63);A2=A1+cd;inv_e,ep=I.rr(1,p),sb+576;mc={};put(mc,YR0+288,Bv);put(mc,ep+224,A2);mc[ep+304]=cd;mc[ep+312]=ep+256;oc=rh("t_c1b",dc(x3=ep,**ars(inv_e,("x25","x26","x27","x28"))),mc);w1=inv_e*(Bv-A2)%p;w2=(w1+cd*inv_e)%p;bday_c1b_ok=(get(mc,ep-64)%p,rgi(oc,("x8","x9","x10","x11"))%p,oc["x29"])==(w1,w2,ep+256);ma,b7,b6,am={},I.rr(1,p),I.rr(1,p),ax_bases(0);ae={};na=0;q=am[5]+544;ma[q]=ae[q]=ah
 for k,sh in en(W6S):
  dk=ab+128*k;put(ma,dk,b7*pow(2,sh,p)%p);put(ma,dk+32,0);put(ma,dk+64,b6*pow(2,sh,p)%p);put(ma,dk+96,0)
  for a_ in rn(dk,dk+128,8):ae[a_]=ma[a_]
  na+=arun(rprog,rlab,"tbld",dc(x0=am[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5))),ma);arun(ar,au,"tbld",dc(x0=am[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5))),ae)
 axyx_ok=(na==6*26+20*(4*128+2*16)and ma[q]==ah and ae[q]!=ah);mp,af,ag,aw,ax,ai,ad={},I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),a6([5,3,1,6,7]);ay,ob3=bp60(mp,XTP[0],af,ag,aw,ax),offs6(ad);put(mp,YR0+64,ai);put(mp,YR0+96,ai)
 for k,v in en(ob3[:4]):mp[YR0+96+8*k]=v
 mp[YR0+352],mp[YR0+360]=ob3[4],0
 ol=rh("t_lkp",dc(x1=YR0+64,x3=sb),mp);av=dc(mp);bp60(av,XTP[0],af,ag,aw,ax,ar,au);b3=rh("t_lkp",dc(x1=YR0+64,x3=sb),av);a8=(ay==5*26+20*5*16 and rgi(ol,("x16","x17","x18","x19"))%p==(af*ad+ag+ai)%p and rgi(b3,("x16","x17","x18","x19"))%p!=(af*ad+ag+ai)%p);PB,SK,kok=0x600000,0x500000,True
 def tk(m,cl,dst,starts,steps,nrow):
  for i,v in en(eiw(CE.EIS[cl])):m[PB+8*i]=v
  for j in rn(6):put(m,ab+64*j,steps[j]);put(m,ab+64*j+32,starts[j])
  return arun(fp,fl,"tbldk",dc(x0=dst,x1=ab,x2=nrow,x3=PB,sp=SK),m)
 for cl in CE.EIS:
  me,us,ss={},[I.rr(p)for _ in rn(6)],[I.rr(p)for _ in rn(6)];n=tk(me,cl,sb,ss,us,3);rk=[CE.eis_red(CE.EIS[cl],(ss[j]+t*us[j])%p)for t in rn(3)for j in rn(6)]
  kok&=n==CE.tbldk_count(3,cl)and[(get(me,sb+256*t+KOFF[j],3),get(me,sb+256*t+KOFF[j]+24,2))for t in rn(3)for j in rn(6)]==[(c1%2**192,c2)for c1,c2 in rk]and all(me[sb+256*t+o]==0 for t in rn(3)for o in(120,248))
 Xa=[I.rr(1,p)for _ in rn(6)];wv=[I.gb(256)for _ in rn(3)];bs=CE.eis_red(CE.EIS["EB"],I.rr(p));mk={};ex={}
 for i,(cl,ce,cc)in en((("E",lambda a:2*a**3,lambda a:0),("F",lambda a:-2*a,lambda a:-a**4),("G",lambda a:BT[2]*a*a,lambda a:0))):
  for k,sh,mi,m in FW(i):t=(wv[i]>>sh)&m;tk(mk,CE.KCL[cl][k],XTK[i]+DW(i)[k]+((t&mi)<<8),[(ce(Xa[j])*t*pow(2,sh,p)+(cc(Xa[j])if k==0 else 0))%p for j in rn(6)],[0]*6,1)
  for j in rn(6):ex[i,j]=ksum(kents(cl,ce(Xa[j]),wv[i],cc(Xa[j]))+([bs]if i==0 else[]))
  for a,o in zip(YR_K[i],offs7(wv[i],i)):mk[YR0+a]=o
 put(mk,YR0+YEB,bs[0]%2**192,3);put(mk,YR0+YEB+24,bs[1],2);m2=dc(mk);ea=XTK[0]+D7O[2]+(((wv[0]>>W6S[2])&W6M[2])<<8)+24;m2[ea]=(m2[ea]+1)%2**64
 def lines(m):return{(i,j):(get(m,sb+ENTRY*j+(384,448,192)[i]+8*(j%3==1),3),get(m,sb+ENTRY*j+(384,448,192)[i]+8*(j%3==1)+24,2))for i in rn(3)for j in rn(6)}
 oe=rh("t_efgk",dc(x3=sb+5*ENTRY+64,x0=1,x2=2,x14=3,x15=4),mk);rh("t_efgk",dc(x3=sb+5*ENTRY+64,x0=1,x2=2,x14=3,x15=4),m2);want={k:(v[0]%2**192,v[1])for k,v in ex.items()};got,bad=lines(mk),lines(m2)
 efgk_ok=got==want and bad[0,0]!=want[0,0]and all(bad[k]==want[k]for k in want if k!=(0,0))and[oe[r]for r in("x0","x2","x14","x15","x3","x4")]==[1,2,3,4,sb+5*ENTRY+64,C0];ok=bday_c1b_ok and axyx_ok and a8 and kok and efgk_ok;return dc(ok=ok,bday_c1b_ok=bday_c1b_ok,axyx_ok=axyx_ok,lookup_pair_ok=a8,tbldk_ok=kok,efgk_ok=efgk_ok)
def good(d):return(d["acc_ok"]and d["frac_mis"]==0 and d["frac_out"]==0 and d["id_mis"]==0 and d["record_decode_mismatches"]==0 and d["records_ok"]and d["keys_noncanon"]==0 and d["wrong_index"]==0 and d["final_R_times_rho_is_one"]and d["sp_kept"]and d["keys_formula"]and d["fb_form"]and d["fixup_words_ok"])
def yrec(N,K,sig,j):return[N.bank("Y",(sig,0))[0].get(YR0+512*j+64+8*w)for w in rn(56)]
def pak(a8,K):
 I={};aj=ksrc("leaf");a4,Lb=a64.a64_load(aj);a6,a7=kload("secp_fe");rt=kload("bank");gp=a8.rr(1<<66);yg,YN,ap,av=ygen(a6,a7,a8,K,rt,gp);I["y_generator"]=yg;ka,kb=[],[];Xs,Xb=mbx(a8,ka),mbx(a8,kb);ac,ar=machine.mach(ME,K),{};af,aw=34*2**40+0x20000,34*2**40+0x40000
 for q in rn(4):put(ac,af+128*q,Xs[2*q][1]*Xs[2*q][1]%p);put(ac,af+128*q+32,Xs[2*q+1][1]*Xs[2*q+1][1]%p);put(ac,af+128*q+64,2*q+1,2);put(ac,af+128*q+80,2*q+2,2);put(ac,af+128*q+96,q,4)
 arun(*rt,"psel",dc(x0=af,x1=4,x2=aw,sp=34*2**40+0x18000),ac,out_regs=ar);I["bday_xs_psel_ok"]=(ar["x0"]==1 and((Xs[0][1]*Xs[0][1]%p)>>63)!=((Xs[1][1]*Xs[1][1]%p)>>63)and[get(ac,aw+128*q+32*s)for q in rn(4)for s in(0,1)]==[X[1]*X[1]%p for X in Xs]);s0=a8.rr(1<<41);Ys=[ysl(a8,s0)for _ in rn(12)];runs,cx={},{}
 for q2 in(30,36,42):cx[q2]=ctx(K);runs[q2]=lfb(a4,Lb,Xs,Ys,q2,K,rt,c=cx[q2])
 # two members of one sub-ring with different table banks read the same stored frames in place
 c=ctx(K);r1=lfb(a4,Lb,Xs,Ys,42,K,rt,c=c);sig=r1["sig"];fr=[yrec(c["m"],K,sig,j)for j in rn(7)];M_=machine.PZ();[pyr(M_,512*j,Ys[j])for j in rn(7)];r2=lfb(a4,Lb,Xb,Ys,42,K,rt,c=c,pool=8,mi=4,prod=False,W2=0x223456789);tb=[k for(w,k)in c["m"].banks if c["m"].names[w]=="T"];I["y_bank_members"]=dc(ok=r1==runs[42]and good(r2)and fr==[yrec(c["m"],K,sig,j)for j in rn(7)]==[[(V63|5<<26)+(j<<42)if w==37 else M_[512*j+64+8*w]for w in rn(56)]for j in rn(7)]and sorted(tb)==[(7,s0),(8,s0)]and r1["g"]!=r2["g"],records=r1["records"]+r2["records"],table_banks=len(tb))
 c["m"][YR0-64]=YR0+64;r3=lfb(a4,Lb,Xs,Ys[7:],24,K,rt,c=c,W2=9);yb=I["y_bank_members"];yb["ok"]&=good(r3)and fr[4:]==[yrec(c["m"],K,sig,j)for j in rn(4,7)];yb["records"]+=r3["records"]
 # yround's own records through ystore, each to the sub-ring of its slice; one of them consumed in place
 c=dc(ctx(K),m=YN);KY=V63|a8.gb(41)<<1;ok,ss=ysto(c,K,rt,gp,9,KY,3);P=YN.bank("S",(gp%2**64,gp>>64))[0];ok&=len(set(ss))==3 and all(yrec(YN,K,ss[i]<<15|9,0)==[KY+(i<<42)if w==37 else P[YRP+512*i+64+8*w]for w in rn(56)]for i in rn(3));ry=lfb(a4,Lb,Xs,[cy(ap,0)],6,K,rt,c=c,cls=9,prod=False);I["yround_ystore"]=dc(ok=ok and good(ry)and ry["sig"]==ss[0]<<15|9,records=3,consumed=ry["records"])
 # a second chunk: the member's fwd and bwd run chunk 0 first, and chunk 1 runs on the entries they left, with no
 # set-up write of the entries in between; the cursor stands at frame CH-1 (a set-up write for CH-1 earlier
 # records); the member's bwd on chunk 1 overwrites two words of frame CH-1, which the driver saves and restores
 c=ctx(K);N=c["m"];lfb(a4,Lb,Xs,Ys,24,K,rt,c=c,stop=True);Yc=[ysl(a8,s0)for _ in rn(5)];sg=s0<<15|5;machine.sel(ME,rt,K,N,"Y",sg);N[YR0-64]=YR0+64+512*(K.CH-1);c["ysig"].add((s0,5));g1=1024*sg+2;machine.gsel(ME,rt,K,N,g1);[pyr(N,YRP+512*j,Yc[j])for j in rn(5)];ok,_=ysto(c,K,rt,g1,5,V63,5);f0=yrec(N,K,sg,K.CH-1);rc=lfb(a4,Lb,Xs,Yc[1:],24,K,rt,c=c,prod=False,j0=K.CH,u0=1,xw=0);I["chunk1"]=dc(ok=ok and good(rc)and f0==yrec(N,K,sg,K.CH-1)and None not in f0 and N.bank("Y",(sg,0))[0][YR0-64]==YR0+64+512*(K.CH+4),records=rc["records"])
 # runs[36]'s six Ys again as a second chunk on its machine; each insert must list its leaf against the first's ids
 rb=lfb(a4,Lb,Xs,Ys[5::-1],36,K,rt,c=cx[36],W2=9,KN=V63,rho=1<<250|M64,mi=4,cls=6);I["same_leaves"]=dc(ok=good(runs[36])and good(rb)and rb["fixup_extra"]["fwd"]>0 and rb["cand"]==[(0x123456789+l%6,V63|5<<26|5-l//6<<42,9+l%6,V63|l//6<<42)for l in rn(36)for _ in rn(6)],pairs=len(rb["cand"]))
 I["planted_collision"]=pcol(a8,K,a4,Lb,rt,Xs,Xb,ka[0],kb[0])
 I["runs"]={q:dc(fwd=d["fwd"],bwd=d["bwd"],keys=d["keys"],fixup_extra=d["fixup_extra"],records_ok=d["records_ok"])for q,d in runs.items()};au={q:{k:d[k]-d["fixup_extra"][k]for k in("fwd","bwd")}for q,d in runs.items()};I["fb_hold"]=all(d["fb_form"]and d["fixup_words_ok"]for d in runs.values());I["counts_ok"]=all(au[q]["fwd"]==K.fwd_count(q)and au[q]["bwd"]==K.bwd_count(q)for q in runs);I["fixups"]=dc(sites=len(fixup_extra(a4)),flags_dead_at_every_join=fdaf(a4,Lb));al=[o for mn,o in a4 if mn in("tbz","tbnz")];I["loop_b34"]=al==[["x5","#9","3f"],["x3","#37","1b"],["x0","#%d"%K.SLOT_BITS,"9f"],["x0","#%d"%K.SLOT_M,"9f"]];I["keys_hold"]=all(d["keys_formula"]for d in runs.values());I["frac_err"]=[[min(d["frac_err"][g][0]for d in runs.values()),max(d["frac_err"][g][1]for d in runs.values())]for g in(0,1)];I["all_correct"]=(yg["mismatches"]==0 and yg["centres_ok"]and yg["counts_ok"]and I["bday_xs_psel_ok"]and all(good(d)for d in runs.values()));az=Rd(SEED+49979687)
 for ay,aq,ao in(("window_batch",("Xm","diff"),"W"),("boundary_batch",("D",),"D")):w=edge_batch(a4,Lb,Xs,Ys,K,rt,az,aq,ao);I[ay]=dc(in_W=w["in_W"],in_D=w["in_D"],ok=(w[{"W":"in_W","D":"in_D"}[ao]]>=4 and w["wrong_index"]==0 and w["keys_noncanon"]==0 and w["records"]==6*24 and w["records_ok"]and w["keys_formula"]and w["fb_form"]and w["fixup_words_ok"]))
 gp_,gl=a64.a64_load(aj.replace("b.lo 9f\n","cbnz xzr,9f\n").replace("b.hs 9f\n","cbnz xzr,9f\n"));g=lfb(gp_,gl,Xs,Ys,24,K,rt);ag=sum(K.FIX_SITES_QUAD.values());ah=sum(K.FIX_EXTRA[m]*K.FIX_SITES_QUAD[m]for m in K.FIX_SITES_QUAD);I["forced_fixups"]=dc(fixup_extra=g["fixup_extra"],ok=(g["frac_mis"]==0 and g["frac_out"]==0 and g["id_mis"]==0 and g["fb_form"]and g["fixup_words_ok"]and g["records_ok"]and sum(map(sum,g["fixups_5_6_8"].values()))==6*ag and sum(g["fixup_extra"].values())==6*ah));ai,ok=set(),True
 N=machine.mach(ME,K)
 for _ in rn(3):a=a8.rr(1,p);put(N,A_SP+256,a);ai.add(arun(a6,a7,"fe_inv",dc(x0=A_SP+512,x1=A_SP+256,sp=A_SP),N,limit=10**6));ok&=get(N,A_SP+512)%p==pow(a,p-2,p)
 I["fe_inv"]=dc(counts=sorted(ai),correct=ok and ai=={K.INV});ak=("SHX","FDX","F5","FS","FADDW","LKFOLD4","DFX","ROWSG");counts={k:mcc(aj,k,ak+("FWD6_ONE","BWDPK01","BWDPK00","BWDPK10"))for k in("fwd","bwd","keys")};a9,aa=counts["fwd"]["FWD6_ONE"],sum(counts["bwd"][k]for k in("BWDPK01","BWDPK00","BWDPK10"));a5={}
 for n in ak:a5[n]=(counts["fwd"][n]//2+counts["bwd"][n])//3
 a5["BLK"],a5["DF"],a5["RW"]=a5.pop("LKFOLD4"),a5.pop("DFX"),a5.pop("ROWSG");I["site_counts"]=dc(per_x_quad=a5,pair_copies=dc(fwd=a9,bwd=aa));I["sites_cert"]=(a5==K.FIX_SITES_QUAD and a9==4 and aa==6 and sum(counts["keys"].values())==0);ax=Rd(SEED+15485863);I["site_checks"]=sites.site_checks(ME,ax,a6,a7);I["tbld_production"]=tbpt(aj,Rd(SEED+27182818));b1=(BETA*Ys[0][0]%p,Ys[0][1]);b0,good_=dgb(a4,Lb,[Xs[0],b1,Xs[1],PROD.neg(Ys[0])]+Xs[2:4],[Ys[0]],K,rt);I["d_zero_batch"]=dc(acc_zero=b0%p==0,fractions_within_bounds=good_,all_other_fractions_wrong=good_==0);ae=rt[0];ab=len(a4)+len(ae)+len(a6);I["prog_insns"]=dc(leaf=len(a4),bank=len(ae),secp_fe=len(a6),total=ab,fits_2_18_bytes=4*ab<2**18,stated=f"{len(a4):,} / {len(ae):,} / {len(a6):,} ({ab:,} insns={4*ab:,}"in K.d8);I["insert_listing"]=ilst(a4,Lb,K);return I
def ilst(a4,Lb,K):
 # the insert listing printed in proof.md section 2 is the assembled text of keys' first emit site (11 instructions)
 i=next(i for i,(mn,o)in en(a4)if mn=="lsr"and o[2]=="#23"and i>Lb["keys"][0]);t=["%s %s"%(mn,",".join(o))for mn,o in a4[i:i+K.INS]];return dc(ok=all(x in K.d8 for x in t)and t[3]==t[5]=="dsb sy"and t[4].startswith("stp x0,")and t[8].startswith("cbnz x0,")and len(t)==K.INS,text=t)
def pcol(o,K,a4,Lb,rt,Xa,Xb,k1,k2):
 # Planted collision, end to end: x is secret, Q=[x]P, X=[b]Q+[c]P for the two generators' first X, Y1=[y1]P and
 # Y2=X1+Y1-X2=[y2]P, so the leaves x(X1+Y1) and x(X2+Y2) coincide. Both go through ystore, fwd, bwd and keys on one
 # machine; the second insert lists the pair. The resolve step (Phase 4: both leaves' scalars from their ids, then
 # x) is a Python model here; RESOLVE and REQ_SERVE are analytic allowances, not kernels.
 x=o.rr(1,r);Q=PROD.mul(x,GP);b1,b2=o.rr(1,r),o.rr(1,r);c1,c2=(k1-b1*x)%r,(k2-b2*x)%r
 while True:
  y1=o.rr(1,r);y2=(k1+y1-k2)%r;Y1,Y2=PROD.mul(y1,GP),PROD.mul(y2,GP);L=PROD.add(Xa[0],Y1)
  if Y2 and L and not any(in_W(L[0],t)or in_D(L[0],t)for t in(0,2)):break
 c=ctx(K);ra=lfb(a4,Lb,Xa,[Y1],6,K,rt,c=c,W2=0x1000,KN=V63|11<<1);rb=lfb(a4,Lb,Xb,[Y2],6,K,rt,c=c,pool=8,mi=4,W2=0x2000,KN=V63|13<<1);cand=rb["cand"];ids=cand==[(0x1000,V63|11<<1,0x2000,V63|13<<1)];got=None
 if ids:
  # ids -> (generator, j, u): here j=0, u=0 on both sides; scalars of the two leaves are b x + c + y
  for e in(1,r-1):
   for m in rn(3):
    f=e*pow(LAM,m,r)%r;d=(b1-f*b2)%r
    if d and PROD.mul((f*(c2+y2)-c1-y1)*pow(d,-1,r)%r,GP)==Q:got=(f*(c2+y2)-c1-y1)*pow(d,-1,r)%r
 return dc(ok=ids and got==x and good(ra)and good(rb)and PROD.add(Xb[0],Y2)==L,listed=len(cand),resolved_x_is_planted=got==x)
def lpt(o,I):c1=(2*(FA*o+(FA-FB)*I)+p)//(2*p);J=(2*((FA-FB)*o-FB*I)+p)//(2*p);return c1*FB+J*(FA-FB),c1*(FA-FB)-J*FA
def planted_x(o,d1,d2,v=None,j=None,fx=None):
 # x whose class label has the given d1, d2; fx=(d4lo, klo, hbh) also fixes the bank's three d4 bits, the slot
 # (K3 mod 2^47 = klo) and, unless hbh is None, the K3 bits above bit 49
 while True:
  J=v
  if J is None:
   d4l,kl,hb_=fx or(o.rr(8),o.rr(16),None);ha=(o.gb(20)<<44)|(d4l<<41)|(kl<<20)|o.rr(200,(1<<20)-200);la=(d1<<23)|o.rr(200,(1<<23)-200);hb=((o.gb(14)if hb_ is None else hb_)<<50)|(o.rr(2)<<20)|o.rr(200,(1<<20)-200);lb=(d2<<23)|o.rr(200,(1<<23)-200);I,uu=lpt(((ha<<64|la)*p)>>128,((hb<<64|lb)*p)>>128)
   if not(0<I<p and uu==BETA*I%p and exk(I)[0]==I):continue
   rf=label(I)
   if not(rf[0]==d1 and rf[1]==d2 and rf[3]&((1<<50)-1)==kl and rf[2]>>41&7==d4l and(hb_ is None or rf[3]>>50==hb_)):continue
   J=I
  m=o.rr(3)if j is None else j;x=J*pow(BETA,(3-m)%3,p)%p
  if exk(x)[1]==m and not any(in_W(x,o)or in_D(x,o)for o in(0,2)):return x,J
def lcert():
 I=io.StringIO()
 with contextlib.redirect_stdout(I):import certificate as K
 o=[" ".join(ln.split())for ln in I.getvalue().splitlines()];return K,[ln for ln in o if" ok"in ln or"FAIL"in ln or ln.startswith(("q1 ","E_total","declared:"))]
def part_d(K):return dc(passed=not K.FAILED,failed=K.FAILED,checks_passed=len(CLS),tail_summary=CLS[-3:],claim="%d.%015d"%divmod(int(K.CLAIM_E*10**15),10**15),q1=K.Q1GD,per_leaf_declared=format(float(K.DECLARED["per_leaf"]),".4f"),mu=float(K.TFULL["mu"]),opened_success_lower_bound=float(K.LO_RUN),gate_prob=float(K.P_GATE),success_lower_bound=float(K.GATED_LO),per_leaf=float(K.PER_LEAF_OPEN),groups=K.G_GROUPS)
def main():
 global CLS,CE;a7=Rd(SEED);K,CLS=lcert();CE=K;OUT["input"]=cfg;I=OUT["A_identities"]=pai(a7);a4=(all(I[k]for k in("beta3","lam3","lamG_betax","no_x0","win_nz","label_constants"))and I["class_rule"]["mismatches_outside_window_set"]==0 and I["class_labels"]["ok"]and I["c132_ok"])
 o=OUT["A_kernels"]=pak(a7,K);J=OUT["B_machine"]=machine.part_b(ME,K);OUT["counts_match_certificate"]=cm=(all(o[k]for k in("fb_hold","counts_ok","loop_b34","keys_hold","all_correct","sites_cert"))and o["fixups"]["flags_dead_at_every_join"]and(EH,EL)==tuple(map(tuple,(K.ERRH,K.ERRL)))and all(o[k]["ok"]for k in("window_batch","boundary_batch","forced_fixups","tbld_production","y_bank_members","yround_ystore","chunk1","same_leaves","planted_collision","insert_listing"))and o["fe_inv"]["correct"]and o["site_checks"]["all_true"]and o["d_zero_batch"]["acc_zero"]and o["d_zero_batch"]["all_other_fractions_wrong"]and o["prog_insns"]["fits_2_18_bytes"]and o["prog_insns"]["stated"]and J["formulas_ok"]and all(J[k]["ok"]for k in("occupied","select","zero","y_bank","phase0","gate","ls64","va55","encoding","keys_edges","meter","registry_trace","poison_loads")))
 a8=OUT["D_certificate"]=part_d(K);OUT["all_checks_pass"]=bool(a4 and cm and a8["passed"]);OUT["seconds"]=round(time.time()-T0,2);json.dump(OUT,sys.stdout,separators=(",",": "));print()
 if not OUT["all_checks_pass"]:sys.exit(1)
if __name__=="__main__":main()
