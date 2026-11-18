EB="EB";D="D";fQ="fwd";fP="bwd";X0="x0"
import contextlib,functools,io,json,os,random,sys,time;Rd=random.Random;Rd.rr,Rd.gb=Rd.randrange,Rd.getrandbits;T0=time.time();cfg=json.load(sys.stdin);SEED=cfg.get("seed",1);BUDGET=cfg.get("seconds",12);b4=min(14.5,BUDGET-0.3);cG=T0+BUDGET-b4;h9,ie,hH=300,100,300;rn,en,dc=range,enumerate,dict
def h5():
 here=os.path.dirname(os.path.abspath(__file__))
 for base in(here,os.path.join(here,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(base,"kernels","a64.py")):return base
 raise FileNotFoundError("kernels/a64.py not found")
HERE=h5();sys.path.insert(0,os.path.join(HERE,"kernels"));import a64;a64.a64_load=functools.lru_cache(maxsize=None)(a64.a64_load);sys.path.insert(0,HERE);import sites,machine;M64=2**64-1;fY=2**42-1;by=2**128;p=2**256-2**32-977;r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16);GX=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",16);GY=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8",16);an=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16);hu=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16);C0=2**32+977;FB=int("e4437ed6010e88286f547fa90abfe4c3",16);FA=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16);aH=((64,80,128,144,288,304),(96,112,160,176,192,208));aA,h7,aR,MX=512,32,64,6;gC,cC,bV,hn,aZ=368,288,64,128,320;bQ=((256,264,272,280,448,456),(408,416,424,432,440),(464,472,480,488,496,504));gU=(0,43,86,128,171,214);hJ=tuple((1<<b)-1 for b in(43,43,42,43,43,42));iq=tuple(en(zip(gU,hJ)));gH=tuple((k,s,m,m)for k,(s,m)in iq);gw=tuple((k,43*k,2**43-1,2**(84 if k>3 else 43)-1)for k in rn(5));BT=(1,an,an*an%p);O,KO=(2,2,0),(1,0,2);dp,ML=20,23;EH,EL=((-241,237),(-136,133),(-361,360)),((-399,401),(-234,228),(-592,588));OUT,gJ,CE={},[],None
def mm(a,v):a[:]=min(a[0],v),max(a[1],v)
def aT(n):return open(os.path.join(HERE,"kernels",n+".s")).read()
def aI(n):return a64.a64_load(aT(n))
def at(a,b,e,r,*x,**k):return a64.a64_run(a,b,e,dc(r,x30=iu),*x,**k)
class cH:
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
 def neg(U,P):return None if P is None else(P[0],-P[1]%U.q)
 def mul(U,k,P):
  k%=U.r
  if not k or P is None:return None
  q,X,Y,Z=U.q,P[0],P[1],1
  def hw(X,Y,Z):S,M=4*X*Y*Y%q,3*X*X%q;o=(M*M-2*S)%q;return o,(M*(S-o)-8*pow(Y,4,q))%q,2*Y*Z%q
  for b in bin(k)[3:]:
   if X is not None:
    if Y==0:X=None
    else:X,Y,Z=hw(X,Y,Z)
   if b=="1":
    if X is None:X,Y,Z=P[0],P[1],1
    else:
     Z2=Z*Z%q;U2,S2=P[0]*Z2%q,P[1]*Z2*Z%q
     if X==U2:
      if(Y+S2)%q==0 or Y==0:X=None
      else:X,Y,Z=hw(X,Y,Z)
     else:H,R=(U2-X)%q,(S2-Y)%q;H2=H*H%q;H3=H*H2%q;V=X*H2%q;X=(R*R-H3-2*V)%q;Y=(R*(V-X)-Y*H3)%q;Z=Z*H%q
  if X is None or Z%q==0:return None
  J=pow(Z,-1,q);I=J*J%q;return(X*I%q,Y*I*J%q)
aE=cH(p,r);GP=GX,GY
def limbs(v,n=4):return[v>>64*k&M64 for k in rn(n)]
def put(N,a,v,n=4):
 for k,w in en(limbs(v,n)):N[a+8*k]=w
def get(N,a,n=4):return sum(N[a+8*k]<<64*k for k in rn(n))
def fL(regs,o):return sum(regs[n]<<64*k for k,n in en(o))
def fO(v,o):return{n:w for n,w in zip(o,limbs(v,len(o)))}
def go(x):v=[x,an*x%p,an*an%p*x%p];j=min(rn(3),key=lambda t:v[t]);return v[j],j
def aF(v):return(v<<128)//p
def hE(v):return(v<<64)//p
def ws(*v):return sum(x>>64 for x in v)%2**64<<64|sum(x&M64 for x in v)%2**64
def bP(x):return ws(x^by-1,1<<64|1)
def bk(a,b):o,I,J,lb=a>>64,a&M64,b>>64,b&M64;d1,d2=I>>24,lb>>23;return d1,(o>>41&1)<<41|d2,o>>41<<41|d2,J>>21<<21|o>>20&(1<<21)-1
def label(v):return bk(aF(v),aF(an*v%p))
def aM(x):return label(go(x)[0])
def aS(v):return bP(aF(v%p))
def cK(v):return aS(v),aS(an*v)
def fold(z):return(z%2**256+(z>>256)*C0)%2**256
def lkfold(z):z=z%2**256+(z>>256)*C0;return z-2**256+C0 if z>>256 else z
def bT(a,b):return(a-b if a>=b else a-b+p)%2**256
def uf(z,l):m=2**(64*l);return z%2**256-z%m+(z+(z>>256)*C0)%m
def fg(v,k):return v%2**(64*k)+(v>>256)*C0>=2**(64*k)
def s6(o,w):return sum(o*(w>>s&m)*pow(2,s,p)%p for _,s,m,_ in gH)
def p1m(J,T,f,g=None):
 if g is None:x=f;f,g=lambda b:s6(2*x%p,b),lambda b:s6(2*x*x%p,b)
 q=J*J;m=[(q>>256+64*k&M64)*C0 for k in rn(4)];w=[(m[i]&M64)+(m[i-1]>>64)for i in(1,2,3)];V=q%2**256+(m[0]&M64)+sum((x&M64)<<64*k for k,x in zip((1,2,3),w))+(m[3]>>64<<256);b=uf(V,2);S,Q=f(b),J+g(b);c,a=uf(S,1),uf(Q,1);Z=sum((c>>64*k&M64)*T[k]for k in rn(4));z=uf(Z,2);d,e=a-z,a+z
 return b,d%2**256-d%2**64+(d-C0*(d<0))%2**64,e%2**256-e%2**64+(e+C0*(e>>256))%2**64,[x>>64 for x in w]+[fg(V,2),fg(S,1),fg(Q,1),fg(Z,2),d<0 and d%2**64<C0,e>>256 and e%2**64+C0>=2**64]
def gF(c):
 R=[c]
 for _ in rn(3):R+=[(R[-1]%2**192*2**64+(R[-1]>>192)*C0)%2**256]
 return R
def bX(w):return(w+((w>>192)*C0>>64))%2**256
def bg(r):return r+(r*C0+2**255>>256)
def bq(w):return sum(bg(an*(w>>43*i&(1<<43)-1)*pow(2,43*i,p)%p)for i in rn(6))%2**256
def bS(c1,c2,w,y1):
 y0=bX(w);h=lambda a,b:a*b>>64;m,n,d,a,f,g=limbs(c1%2**192,3),limbs(c2,2),limbs((c1-c2)%2**192,3),limbs(y0),limbs(y1),limbs((y0+y1)%2**256)
 def pr(x,y):t=len(x)>2;return(h(x[0],y[2])+x[0]*y[3]+h(x[1],y[1])+x[1]*y[2]+(x[2]*y[1]if t else 0))&M64,h(x[0],y[1])+x[0]*y[2]+h(x[1],y[0])+x[1]*y[1]+(x[2]*y[0]if t else 0)&M64
 AH,AL=pr(m,a);BH,BL=pr(d,g);P=pr(n,f);return(AH+P[0]&M64)<<64|AL+P[1]&M64,(BH-AH&M64)<<64|BL-AL&M64
def bc(z):return sum(x[0]for x in z),sum(x[1]for x in z)
def aY(cl,ce,v,cc):return[CE.eis_red(CE.EIS[CE.KCL[cl][k]],(ce*(v>>sh&m)*pow(2,sh,p)+(cc if k==0 else 0))%p)for k,sh,_,m in FW("EFG".index(cl))]
def cJ(a,b):return[bc(aY("E",2*a**3,b,0)+[CE.eis_red(CE.EIS[EB],pow(b,4,p))]),bc(aY("F",-2*a,pow(b,3,p),-a**4)),bc(aY("G",BT[2]*a*a,3*b*b%p,0))]
def bY(N,e,s):return[get(N,e+o,2)for o in aH[s]]
def components(X,Y):
 a,b=X[0],Y[0];E,F,G=(b**4+2*a**3*b)%p,(a**4+2*a*b**3)%p,3*a*a*b*b%p;W=[]
 for w in(pow(Y[1]+X[1],-2,p),pow(Y[1]-X[1],-2,p)):W+=[E*w%p,F*w%p,G*w%p]
 return W
def bu(t,v):return[((t>>64*k)-(v>>64*k)+2**63)%2**64-2**63 for k in(1,0)]
def hI(N,e,X,Y,k):
 cs,W=components(X,Y),[]
 for s in(0,1):o=((cs[3*s]-Y[0])%p,-cs[3*s+1]%p,BT[2]*cs[3*s+2]%p);W+=[bu(aF(BT[i%2]*o[i//2]%p),ws(v,get(N,k+16*i,2)*(i<2)))for i,v in en(bY(N,e,s))]
 return W
def aW(o,t):b0,b1,a0,a1,g0,g1=o;n=bP;return((ws(a0,a1,g0,n(b0),n(b1)),ws(b0,g1,n(a0))),(ws(n(a0),n(b0),n(b1),n(g0),n(g1)),ws(b0,g0,n(a1))),(ws(a0,a1,b0,n(g0),n(g1)),ws(b1,g0,n(a0))))[t]
def aB(w):I,o=w[0]>>64,w[1]>>64;s=I+o&M64;j=(2 if o+s>M64 else 1)if I>=o else(2 if I+s>M64 else 0);n=ws(w[0],w[1])^by-1;return(j,)+bk(*((w[0],w[1]),(w[1],n),(n,w[0]))[j])
def in_W(x,o=0):x=BT[o]*x%p;v=[x<<64,an*x%p<<64,an*an%p*x%p<<64];o=p<<64;lo,hi=min(EH[0][0],EH[1][0]),max(EH[0][1],EH[1][1]);return(any(w<hi*p or w>=o+lo*p for w in v[:2])or v[2]<(1-2*lo)*p or v[2]>o-(3+2*hi)*p or abs(v[0]-v[1])<(hi+1-lo)*p or any(3*lo*p<=w-v[2]<(3*hi+3)*p for w in v[:2]))
def in_D(x,o=0):x=BT[o]*x%p;v=(x,an*x%p,an*an%p*x%p);return any(u%2**mb<hi or u%2**mb>=2**mb+lo for w,I,o in zip(v,EH,EL)for u,mb,(lo,hi)in((aF(w)>>64,dp,I),(aF(w)&M64,ML,o)))
def bR(o,t,a,b):return(o[0]-b+pow(an,t,p)*(o[1]-a)+pow(an,2*t,p)*o[2])%p
def bh(L):o,I=pow(3,-1,p),an*an%p;return((L[0]+L[1]+L[2])*o%p,(L[0]+I*L[1]+an*L[2])*o%p,(L[0]+an*L[1]+I*L[2])*o%p)
def cM(o):return o.rr(p)if o.random()<0.85 else o.rr(C0)+p
def gR(x,t,a,b,o,raw=False):
 w,c,k=cM(o)or 1,[o.rr(p)for _ in rn(3)],o.rr(3);c[k]=0;v=(x+b+BT[t]*a)*pow(w,-1,p)-sum(BT[t*i%3]*c[i]for i in rn(3));c[k]=v*pow(BT[t*k%3],-1,p)%p;y1=bq(w);W=[]
 for i,(ci,y)in en(zip((c[0],-c[1],BT[2]*c[2]),(b,None if raw else a,None))):cl=CE.KCL["EFG"[i]]+((EB,)if i==0 else());ks=[o.rr(p)for _ in rn(len(cl)-1)];A,B=bS(*bc([CE.eis_red(CE.EIS[q],u)for q,u in zip(cl,ks)]+[CE.eis_red(CE.EIS[cl[-1]],(ci-sum(ks))%p)]),w,y1);cw=(cK(y)if i<1 else(aF(y),aF(an*y%p)))if y and(i or not raw)else(0,0);W+=[ws(A,cw[0]),ws(B,cw[1])]
 return W,((c[0]*w-b)%p,(c[1]*w-a)%p,c[2]*w%p)
def g9(o,kind):
 m=o.rr(3)
 if kind==D:b=o.rr(2);mb=(dp,ML)[b];u=o.gb(64-mb)<<mb|o.choice((o.rr(70),2**mb-1-o.rr(70),o.rr(2**mb)));v=-(-(u<<64|o.gb(64)if b==0 else o.gb(64)<<64|u)*p//by);return v*pow(an,-m,p)%p or 1
 off=o.rr(160<<192);w=off if o.random()<0.5 else p-1-off
 if kind=="Xm":return w*pow(an,-m,p)%p or 1
 return w*pow((1-an,1-an*an,an-an*an)[m]%p,-1,p)%p or 1
def bd(I):o=time.time();return o+I*max(0.0,cG-o)
def hg(o):I=o.choice([1,an,an*an,1-an,1-an*an,an-an*an]);dd=o.rr(-(1<<o.rr(150,250)),1<<o.rr(150,250));return dd*pow(I%p,-1,p)%p or 1
def hG(x,o):t,a,b=o.rr(3),o.rr(1,p),o.rr(1,p);j=aB(aW(gR(x,t,a,b,o)[0],t))[0];return in_W(x,O[t]),(j+O[t])%3!=go(x)[1]
def class_labels(o):
 W=dc(classes=0,edge_classes=0,members=[0,0,0],in_W=0,in_D=0,wrong_index=0,mismatches=0,errors_h=[[999,-999]for _ in rn(3)],errors_l=[[999,-999]for _ in rn(3)]);dl=bd(0.05)
 while W["classes"]<hH or(W["classes"]<4*hH and time.time()<dl):
  I=W["classes"]%3==2;x0=g9(o,(D,"Xm","diff")[W["classes"]//3%3])if I else o.rr(1,p);v=go(x0)[0];W["classes"]+=1;W["edge_classes"]+=I;a4=label(v)
  for m in rn(3):
   x=v*pow(an,(3-m)%3,p)%p;t,a,b=o.rr(3),o.rr(1,p),o.rr(1,p);w=aW(gR(x,t,a,b,o)[0],t)
   for k,c in en((w[0],w[1],ws(w[0],w[1])^by-1)):eh,el=bu(aF(BT[(O[t]+k)%3]*x%p),c);mm(W["errors_h"][k],eh);mm(W["errors_l"][k],el)
   f=aB(w)
   if in_W(x,O[t]):W["in_W"]+=1;continue
   if(f[0]+O[t])%3!=m:W["wrong_index"]+=1;continue
   W["members"][m]+=1
   if in_D(x,O[t]):W["in_D"]+=1
   elif f[1:]!=a4:W["mismatches"]+=1
 J=W["mismatches"]==0 and W["wrong_index"]==0 and W["in_D"]>0 and W["in_W"]>0 and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(W["errors_h"],EH))and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(W["errors_l"],EL));return dc(ok=J,classes=W["classes"],in_W=W["in_W"],in_D=W["in_D"],errors_h=W["errors_h"],errors_l=W["errors_l"])
def hU(ac):
 I={};I["beta3"]=pow(an,3,p)==1 and an!=1;I["lam3"]=pow(hu,3,r)==1 and hu!=1;I["lamG_betax"]=aE.mul(hu,GP)==(an*GX%p,GY);I["no_x0"]=pow(7,(p-1)//2,p)==p-1;I["win_nz"]=all(c%p for c in(1-an,1-an*an,an-an*an,3));I["label_constants"]=FA*FA-FA*FB+FB*FB==p and FB*an%p==FA-FB and(FA-FB)*an%p==p-FA;bad=drops=a8=0;a6={"random":0,"near_tie":0};dl=bd(0.10)
 while a8<26000 and(a8<h9 or time.time()<dl):x,kind=(ac.rr(1,p),"random")if a8<20000 else(hg(ac),"near_tie");a6[kind]+=1;a8+=1;a7,a9=hG(x,ac);drops+=a7;bad+=a9 and not a7
 extra=max(0,ie-a6["near_tie"]);ab=Rd(SEED+7919)
 for _ in rn(extra):a7,a9=hG(hg(ab),ab);a6["near_tie"]+=1;drops+=a7;bad+=a9 and not a7
 I["class_rule"]=dc(checked=a6,in_W=drops,mismatches_outside_window_set=bad);I["class_labels"]=class_labels(Rd(SEED+104729));J=Rd(SEED+918273);o,c132_ok=1<<64,True
 def bU(i_,j_,d_):ai=max(0,i_-d_);return 0<=i_<o and(ai<=j_<=min(o-1,o+d_-2*i_)or max(ai,o-d_-i_)<=j_<=min(o-1,2*o+d_-2*i_))
 for ak in rn(600):ae,ad=J.choice([0,1,o-1,o//2,J.gb(64)]),J.choice([0,1,o-1,J.gb(64)]);aj=aB((ae<<64,ad<<64))[0];ag=(ae,ad,~((ae+ad)%o)&M64);c132_ok&=bU(ag[aj],ag[(aj+1)%3],0);af=hg(J)if ak%3==0 else g9(J,("Xm","diff",D)[ak%3]);t_=J.rr(3);al=aB(aW(gR(af,t_,J.rr(1,p),J.rr(1,p),J)[0],t_))[0]+O[t_];ah=pow(an,al,p)*af%p;a5,a4=hE(ah),hE(an*ah%p);c132_ok&=bU(a5,a4,90)or a5<64 or a5>=o-64 or a4<64 or a4>=o-64 or min(abs(a5+a4-o),abs(a5+a4-2*o),abs(2*a5+a4-o),abs(2*a5+a4-2*o),abs(a5-a4))<=64
 I["c132_ok"]=c132_ok;return I
iu=0x5E5E5E5E5E5E5E5E;gq=1<<40;cA=(34<<40)+0x7000000;b8=42<<40;fS=b8+(1<<37);fK,gI,gG,hx=40<<40,1<<31,2**64-32768*gq,0xfff8<<48;aK=4096*gq;cY=24576*gq;di=16384*gq;g_=tuple(c<<50 for c in(0,2,8,4,6,9));gn,hS=(tuple(c<<s for c in(0,2,4,6,8))for s in(49,50));gM=tuple(gG+k*10240*gq for k in rn(3));fT=(19456*gq,19968*gq)
def aG(j):return(aK+(4*j+0<<49),aK+(4*j+1<<49),di+(2*j+0<<48),aK+(4*j+2<<49),aK+(4*j+3<<49),di+(2*j+1<<48))
bD=(0,1,3,2,4,5,3)
def FW(i):return(gH,gw,gH)[i]
def DW(i):return(g_,hS,g_)[i]
def a3(v):return[gn[k]+((v>>o&m)<<6)for k,o,m,_ in gw]
def bw(v,i):return[gM[i]+DW(i)[k]+((v>>o&m)<<8)for k,o,m,_ in FW(i)]
class hM(dc):
 def __init__(U,rprog,rlab):super().__init__();U.ih,U.h8,U.touched=rprog,rlab,set()
 def __getitem__(U,a):U.touched.add(a);return super().__getitem__(a)
 def __setitem__(U,a,v):U.touched.add(a);return super().__setitem__(a,v)
 def get(U,a,d=None):
  if a in U:U.touched.add(a)
  return super().get(a,d)
 def __missing__(U,a):
  if 0<=a<16*machine.H1 and a%8==0:d=a//8;b=machine.blk1(d);at(U.ih,U.h8,"init1",dc(x0=1,x1=a,x2=b+32,x3=machine.B1),U);U.touched.add(a);return super().__getitem__(a)
  for lo,hi,first,o in machine.bB[2:]:
   if lo<=a<hi:U.touched.add(a);super().__setitem__(a,0);return 0
  raise KeyError(hex(a))
def bC(rprog,rlab,N,tb,c1,a5,c2,a4,w,bs=None):
 o=(34<<40)+0x180000
 for k,J,mi,m in(gH if bs else gw):
  I=w>>J&m;ad=(bs[k]if bs else tb+gn[k])+((I&mi)<<6)
  if ad in N:continue
  put(N,o,0);put(N,o+32,(c1*I*pow(2,J,p)+(a5 if k==0 else 0))%p);put(N,o+64,0);put(N,o+96,(c2*I*pow(2,J,p)+(a4 if k==0 else 0))%p);at(rprog,rlab,"tbld",dc(x0=ad,x1=o,x2=1,x3=1,x4=2,x5=64,x6=0),N)
def bM(v):c_=pow(v,(p+2)//9,p);return c_ if c_*c_*c_%p==v%p else None
def ir(o):
 xs=[]
 for q in rn(4):
  while True:
   I=aE.mul(o.rr(1,r),GP);A0=I[1]*I[1]%p;lo=A0&(1<<63)-1
   if q==0 and lo<1<<62:continue
   d=o.rr((1<<63)-lo,1<<63)if q==0 else o.rr(1,1<<63);J=A0+d
   if J>=p:continue
   y1,a1=aV(J),bM((J-7)%p)
   if y1 is not None and a1 is not None:xs+=[I,(a1,y1)];break
 return xs
def bW(rprog,rlab,N,j,X,w):bC(rprog,rlab,N,0,2*(X[1]*X[1]%p)%p,0,2*X[1]%p,0,w,aG(j))
def bo(rprog,rlab,N,w):
 for o in{w>>43*i&((1<<43)-1)for i in rn(6)}:
  J=hx+(o<<8)
  if J in N:continue
  for idx,wi in en(bD):put(N,J+32*idx,bg(an*o*pow(2,43*wi,p)%p))
  put(N,J+224,0)
def hL(c):return limbs(c["W"][1])+limbs(p-c["U"][1])+limbs(c["kS"])+limbs(c["kT"])+limbs(pow(p,-1,2**192),3)+limbs(c["U"][0]%2**192,3)+limbs(c["W"][0]%2**192,3)+limbs(c["U"][1],2)+limbs(c["W"][1],2)+limbs(c["o"][0]%2**192,3)+limbs(c["o"][1],2)+limbs(c["C"],3)+limbs(c["X2"],2)+[len(c["c"])]+sum((limbs(l1%2**192,3)+limbs(l2%2**128,2)for l1,l2 in c["c"]),[])
bA=(0,40,80,128,168,208)
def cN(N,Xs,b):
 for i,(v,ce,cc)in en(((b,lambda a:2*a**3,lambda a:0),(pow(b,3,p),lambda a:-2*a,lambda a:-a**4),(3*b*b%p,lambda a:BT[2]*a*a,lambda a:0))):
  for k,sh,mi,m in FW(i):
   t=v>>sh&m;ad=gM[i]+DW(i)[k]+((t&mi)<<8)
   if ad in N:continue
   for j in rn(6):a=Xs[j][0];c1,c2=CE.eis_red(CE.EIS[CE.KCL["EFG"[i]][k]],(ce(a)*t*pow(2,sh,p)+(cc(a)if k==0 else 0))%p);put(N,ad+bA[j],c1%2**192,3);put(N,ad+bA[j]+24,c2,2)
   N[ad+120]=N[ad+248]=0
def bJ(rprog,rlab,N,I,Ys,rho):
 a5=rho;o,J=[],[];A=[X[1]*X[1]%p for X in I[:MX]];sA=[(A[2*q]+A[2*q+1])%p for q in rn(3)];pA=[A[2*q]*A[2*q+1]%p for q in rn(3)]
 for yi,Y in en(Ys):
  b,yy=Y;B=yy*yy%p;b3=pow(b,3,p);cN(N,I,b);bC(rprog,rlab,N,fT[0],-sA[1]%p,(pA[1]-7*sA[1])%p,-sA[0]%p,(pA[0]-7*sA[0])%p,b3);bC(rprog,rlab,N,fT[1],0,0,-sA[2]%p,(pA[2]-7*sA[2])%p,b3)
  for pr in rn(3):aa,af=sA[pr],pA[pr];a8=lkfold(sum(((-aa%p*(b3>>sh&m)*pow(2,sh,p)+((af-7*aa)%p if k==0 else 0))%p)for k,sh,_,m in gw)+B*B%p);o+=[a5];J+=[a8];a5=fold(sum((a8>>64*k&M64)*gF(a5)[k]for k in rn(4)))
 R=pow(a5%p,p-2,p)if a5%p else 0;prod=lambda m_,rr:sum((m_>>64*k&M64)*rr[k]for k in rn(4))
 for jp in rn(len(o)-1,-1,-1):
  yi=2*jp//MX;b,yy=Ys[yi];B=yy*yy%p;pr=jp%3;(a1,y1),(a2,y2)=I[2*pr],I[2*pr+1];A1,A2=y1*y1%p,y2*y2%p;dA=(A2-A1)%p;RR=gF(R);inv_e,R=fold(prod(o[jp],RR)),fold(prod(J[jp],RR));a7=fold(prod(inv_e,gF(bT(B,A2))));ad=lkfold(a7+dA*inv_e)
  for j,X,a4 in((2*pr+1,(a2,y2),ad),(2*pr,(a1,y1),a7)):b,u,v,_=p1m(a4,[yy*pow(2,64*k,p)%p for k in rn(4)],X[1]);bW(rprog,rlab,N,j,X,b);bo(rprog,rlab,N,u);bo(rprog,rlab,N,v)
def bj(N,e,o,l):
 j=l%MX;X=o[j]
 for k,w in en((aF(X[0]),aF(an*X[0]%p))):put(N,e+256+16*k,w,2)
 put(N,e+32,X[1]*X[1]%p)
 for k,b_ in en(aG(j)+(0x5A5A5A5A5A5A5A5A,)*2):N[e+320+8*k]=b_
 if j&1:N[e+368]=(o[j][1]**2-o[j-1][1]**2)%p;N[e+376]=e+320
def hv(N,e,Y):
 b,y=Y;B=y*y%p;N[e]=0;c1,c2=CE.eis_red(CE.EIS[EB],pow(b,4,p));put(N,e+gC,c1%2**192,3);put(N,e+gC+24,c2,2);put(N,e+cC,B);put(N,e+bV,B*B%p)
 for k in rn(4):put(N,e+hn+32*k,y*pow(2,64*k,p)%p)
 put(N,e+aZ,aS(b),2);put(N,e+aZ+16,aS(an*b),2);gt=a3(pow(b,3,p))
 for k in rn(4):N[e+96+8*k]=fT[0]+gt[k]
 N[e+352],N[e+360]=fT[0]+gt[4],0
 for i,v in en((b,pow(b,3,p),3*b*b%p)):
  for a,x in zip(bQ[i],bw(v,i)):N[e+a]=x
def bN(Xs,Ys,N,o,a7,R,ybase=None):
 if ybase is None:ybase=fK
 W,prod={},lambda m,rr:sum((m>>64*k&M64)*rr[k]for k in rn(4))
 for jp in rn(len(o)-1,-1,-1):
  l1,l2=2*jp,2*jp+1;(a1,a9),(a2,a8)=Xs[l1%MX],Xs[l2%MX];b,yy=Ys[l2//MX];B,A1,A2=yy*yy%p,a9*a9%p,a8*a8%p;RR=gF(R);inv_e,R=fold(prod(o[jp],RR)),fold(prod(a7[jp],RR));D2=bT(B,A2);dA=(A2-A1)%p;aa=fold(prod(inv_e,gF(D2)));af=lkfold(aa+dA*inv_e)
  for l,(a,yx),J in((l2,(a2,a8),af),(l1,(a1,a9),aa)):
   a5=Ys[l//MX][0];yr=ybase+512*(l//MX);ad=[get(N,yr+hn+32*k)for k in rn(4)];cs=cJ(a,a5);fr=[]
   for w in p1m(J,ad,yx)[1:3]:
    y1=bq(w)
    for i,c in en(cs):A_,B_=bS(*c,w,y1);fr+=[A_,B_]
   W[l]=fr
 return W,R
def bz(X,Y):
 W=[]
 for t in rn(3):o=(pow(an,t,p)*X[0]%p,X[1]);W+=[aE.add(o,Y)[0],aE.add(o,aE.neg(Y))[0]]
 return W
def fixup_extra(I):
 W,o={},-1
 for i,(J,_)in en(I):
  if J in("b.lo","b.hs","cbnz","cbhs"):o=i
  elif J=="stadd":W[i]=i-o
 return W
def cP(J,Lb):
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
   if mn in("b","tbz","tbnz","cbz","cbnz","cbhs"):
    I+=[target(o[-1],j+1)]
    if mn=="b":continue
   I+=[j+1]
 return True
gQ,gE,f3=2304,(34<<40)+0x2000000,(34<<40)+0x4000000
def bZ(i):rho,q=divmod(i,33);return 33*rho+(0 if q==32 else(16-q//2)*(1-2*(q%2)))
@functools.lru_cache(maxsize=None)
def cy(o,i):return aE.add(o[0],aE.mul(bZ(i)%r,o[1]))
def cL(a4,I,a5,K):
 N,chains={},[]
 for c in rn(2):
  C,S=aE.mul(a5.rr(1,r),GP),aE.mul(a5.rr(1,r),GP);chains+=[(C,S)];cb=gE+gQ*c
  for i,P in[(0,C)]+[(i,aE.mul(i,S))for i in list(rn(1,17))+[33]]:put(N,cb+(1088 if i==33 else 64*i),P[0]);put(N,cb+(1088 if i==33 else 64*i)+32,P[1])
  N[cb+2240]=fK+c*gI
 put(N,gE-128,1);put(N,f3+352,an);N[f3+384],N[f3+424]=C0,gE;counts=[];N[f3+1016]=1
 for k,v in en(hL(CE.EIS[EB])):N[f3+1024+8*k]=v
 for a7 in 2,1:N[f3+416]=gE+gQ*a7;counts+=[at(a4,I,"yround",dc(sp=f3),N,limit=10**7)]
 bad=J=0
 for c,rd in((0,0),(1,0),(0,1)):
  for q in rn(33):i=33*rd+q;(b,y),e=cy(chains[c],i),fK+c*gI+512*i;B=y*y%p;c1,c2=CE.eis_red(CE.EIS[EB],pow(b,4,p));ok=get(N,e+gC,3)==c1%2**192 and get(N,e+gC+24,2)==c2 and get(N,e+cC)==B and get(N,e+bV)==B*B%p;ok&=all(get(N,e+hn+32*k)==y*pow(2,64*k,p)%p for k in rn(4));ok&=get(N,e+aZ,2)==aS(b)and get(N,e+aZ+16,2)==aS(an*b);gt=a3(pow(b,3,p));ok&=[N[e+96+8*k]for k in rn(4)]+[N[e+352]]==[fT[0]+v for v in gt]and N[e+56]==pow(b,3,p)>>215;ok&=all([N[e+a]for a in bQ[i_]]==bw(v,i_)for i_,v in en((b,pow(b,3,p),3*b*b%p)));bad+=not ok;J+=1
 a9=[(get(N,gE+gQ*c),get(N,gE+gQ*c+32))for c in rn(2)];ok=bad==0 and a9==[cy(chains[0],33*2+32),cy(chains[1],33+32)];N[f3+1016]=0;counts+=[at(a4,I,"yround",dc(sp=f3),N,limit=10**7)];a8=counts==[K.YR_FIX+2*K.YR_BLK,K.YR_FIX+K.YR_BLK,K.YR_FIX+K.YR_BLK_MINE];return dc(records=J,mismatches=bad,centres_ok=ok,counts=counts,counts_ok=a8),N,chains[0],chains[1]
def gP(a9,Lb,Xs,Ys,q2,ah,K,rprog=None,rlab=None,W2=0x123456789,KN=5<<26,override=None,chain=None,rho=None,ybase=fK):
 assert q2%6==0 and ybase%gI==0
 if rprog is None:rprog,rlab=aI("route")
 fx=fixup_extra(a9);ac={k:dc.fromkeys(fx,0)for k in(fQ,fP)};nY,SP=q2//MX,cA;sb=fS-aA*q2;LS=([(aK,cY),(gG,1<<64),(ybase,ybase+512*nY),(sb,fS)],[(sb,fS)]);N=hM(rprog,rlab);rho=rho or Rd(q2).rr(1,p);bJ(rprog,rlab,N,Xs,Ys[:nY],rho)
 if ah is not None:
  for a in rn(ybase-64,ybase+512*nY,8):N[a]=ah.get(a,0)if(a-ybase)%512>63 else 7
 else:
  for a in rn(ybase-64,ybase,8):N[a]=0
  for yi in rn(nY):hv(N,ybase+512*yi,Ys[yi])
 for yi in rn(nY):N[ybase+512*yi+360]=KN+(yi<<42)
 for l in rn(q2):bj(N,sb+aA*l,Xs,l)
 regs=dc(x0=q2,x1=ybase+64,x3=sb,sp=SP);regs.update(fO(rho,("x24","x25","x26","x27")));W={};nf=at(a9,Lb,fQ,regs,N,out_regs=W,visits=ac[fQ],ls64=LS);words,aq=[W[X0]],fL(W,("x24","x25","x26","x27"));au=[get(N,sb+aA*(2*j)+h7)for j in rn(q2//2)];es=[get(N,sb+aA*(2*j))for j in rn(q2//2)];R0=pow(aq,p-2,p);am=ybase+512*(nY-1);regs=dc(x3=fS,x1=ybase-32,sp=am);regs.update(fO(R0,(X0,"x2","x14","x15")));W={};nb=at(a9,Lb,fP,regs,N,out_regs=W,visits=ac[fP],ls64=LS);words+=[W["x5"]];ae=rho
 for l in rn(q2):ae=ae*(Ys[l//MX][1]**2-Xs[l%MX][1]**2)%p
 a8=fL(W,(X0,"x2","x14","x15"));sp_kept=W["sp"]==am;model,mR=bN(Xs,Ys,N,au,es,R0,ybase);J,a4,a6,a5,leaves=int(mR!=a8),0,0,[[999,-999],[999,-999]],[]
 for l in rn(q2):
  X,Y=Xs[l%MX],Ys[l//MX];e,cs=sb+aA*l,components(X,Y);J+=sum(get(N,e+o,2)!=w for o,w in zip(aH[0]+aH[1],model[l]))
  for i,(dh,dl)in en(hI(N,e,X,Y,sb+aA*(l-l%6+1)+224)):cb=K.ii[i%6//2][i%2];mm(a5[0],dh);mm(a5[1],dl);a4+=not(cb[0][0]<=dh<=cb[0][1]and cb[1][0]<=dl<=cb[1][1])
  ao=bz(X,Y);a6+=sum(bR(cs[3*s:3*s+3],t,X[0],Y[0])!=ao[2*t+s]for s in(0,1)for t in rn(3));leaves+=[ao]
 for(l,s),(x,t,rg)in(override or{}).items():
  f6,c3=gR(x,t,Xs[l%MX][0],Ys[l//MX][0],rg,True)
  for o,v in zip(aH[s],f6):put(N,sb+aA*l+o,v,2)
  for t_ in rn(3):leaves[l][2*t_+s]=(c3[0]+BT[t_]*c3[1]+BT[2*t_%3]*c3[2])%p
 ai=machine.aD(N);regs=dc(x2=W2+6,x3=sb+aR,sp=SP);W={};nk=at(a9,Lb,"keys",regs,N,out_regs=W,dev=ai.dev,ls64=([(sb,fS)],[]));sp_kept&=W["sp"]==SP;a7=machine.aN(K);meter,ad=machine.bF(sys.modules[__name__],N,sb,q2,W2,KN,a7,K);ft=[(t,leaves[l][2*t+s])for l in rn(q2)for s in(0,1)for t in KO];aa=[x for _,x in ft];ab=[in_W(x,O[t])for t,x in ft];af=[in_D(x,O[t])for t,x in ft];ak=sum(not w and(f[0]+O[t])%3!=go(x)[1]for w,f,(t,x)in zip(ab,ad,ft));al=sum(f[1:]!=aM(x)for w,d,f,x in zip(ab,af,ad,aa)if not w and not d);ap=0
 for k,(d,ag)in en(a7.seq if chain else()):l=6*(ag[3]>>42&2**21-1)+ag[2]-W2;s=k//3%2;Y,X=cy(chain,l//MX),Xs[l%MX];ar=[aM(aE.add((pow(an,t,p)*X[0]%p,X[1]),aE.neg(Y)if s else Y)[0])for t in rn(3)];ap+=l!=k//6 or Y!=Ys[l//MX]or(not ab[k]and not af[k]and(d,ag[0]&fY,ag[0],ag[1])not in ar)
 I={k:[sum(v for i,v in ac[k].items()if fx[i]==e)for e in(3,4,6)]for k in ac};aj={k:3*I[k][0]+4*I[k][1]+6*I[k][2]for k in ac};return dc(fwd=nf,bwd=nb,keys=nk,acc_ok=aq%p==ae,frac_mis=J,frac_out=a4,frac_err=a5,id_mis=a6,fixups_3_4=I,fixup_extra=aj,record_decode_mismatches=ap,fixup_words_ok=words==[C0*sum(I[k])for k in(fQ,fP)]and all(sum(I[k])<=w>>32<=sum(I[k])*C0>>32 for w,k in zip(words,(fQ,fP))),fb_form=nf==K.fwd_count(q2)+aj[fQ]and nb==K.bwd_count(q2)+aj[fP],records_ok=ai.blocks==a7.msgs and machine.aJ(N,a7),records=len(a7.seq),key_paths=meter,in_W=sum(ab),in_D=sum(d and not w for w,d in zip(ab,af)),wrong_index=ak,keys_noncanon=al,keys_formula=nk==K.KEY_FIX+K.KEY_INC*q2+meter+a7.extra,final_R_times_rho_is_one=a8*rho%p==1,sp_kept=sp_kept,reg_trace_ok=all(any(lo<=a<hi for rl,_,lo,hi,_ in K.REGISTRY if rl=='G')for a in N.touched))
def hQ(a5,entry,a4):
 o,J,I={},None,[ln.split("//")[0].strip()for ln in a5.split("\n")]
 for s in I:
  if s.startswith(".macro"):J=s.split()[1].rstrip(",");o[J]=[]
  elif s.startswith(".endm"):J=None
  elif J is not None:o[J].append(s)
 W=dc.fromkeys(a4,0)
 def cB(seq):
  for s in seq:
   if not s:continue
   w=s.split()[0].rstrip(",")
   if w in W:W[w]+=1
   if w in o:cB(o[w])
 start=I.index(entry+":");cB(I[start+1:I.index("ret",start)]);return W
def aV(v):s_=pow(v,(p+1)//4,p);return s_ if s_*s_%p==v%p else None
def hT(I,Lb,Xs,Ys):
 rprog,rlab=aI("route");q2,J=24,4;o=fS-aA*q2;a5=Ys*J;LS=([(aK,cY),(gG,1<<64),(fK,fK+512*J),(o,fS)],[(o,fS)]);N=hM(rprog,rlab);bJ(rprog,rlab,N,Xs,a5,7)
 for a in rn(fK-64,fK,8):N[a]=0
 for yi in rn(J):hv(N,fK+512*yi,a5[yi])
 for l in rn(q2):bj(N,o+aA*l,Xs,l)
 regs=dc(x0=q2,x1=fK+64,x3=o,x24=7,x25=0,x26=0,x27=0,sp=cA);W={};at(I,Lb,fQ,regs,N,out_regs=W,ls64=LS);a4=fL(W,("x24","x25","x26","x27"));regs=dc(x3=fS,x1=fK-32,sp=fK+512*(J-1));regs.update(fO(pow(a4,p-2,p),(X0,"x2","x14","x15")));at(I,Lb,fP,regs,N,ls64=LS);good=0
 for l in rn(q2):
  X=Xs[l%MX]
  if(Ys[0][1]**2-X[1]**2)%p:good+=sum(abs(dh)<4096 and abs(dl)<4096 for dh,dl in hI(N,o+aA*l,X,Ys[0],o+aA*(l-l%6+1)+224))
 return a4,good
def bG(a4,Lb,Xs,Ys,I,K,J,o,a5):
 over={}
 for l in rn(24):
  for s in(0,1):
   t=J.rr(3)
   while True:
    x=g9(J,o[(2*l+s)%len(o)])
    if len(over)>=4 or(in_W(x,O[t])if a5=="W"else in_D(x,O[t])and not in_W(x,O[t])):break
   over[l,s]=x,t,J
 return gP(a4,Lb,Xs,Ys,24,I,K,override=over)
def cQ(ac,I):
 ak=aT("route");rprog,rlab=a64.a64_load(ak);ar,au=a64.a64_load(ak.replace("tst x6,#1\ncsel x18,x3,x2,ne\nlsr x6,x6,#1","cmp x8,#2\ncsel x18,x3,x2,eq"));ap,az=a64.a64_load(ac+sites.aP+sites.aX);fp,fl=aI("secp_fe");ah,ab,sb=0xDEADBEEFCAFEF00D,0x300000,fS-aA*MX;b0=([(aK,cY),(gG,1<<64),(fK,fK+512),(sb,fS)],[(sb,fS)]);a6=lambda ts:sum(t<<43*k for k,t in en(ts))
 def h_(m,b,nw):
  for i in rn(nw):m[b+8*i]=ah
 def b9(m,tb,c1,a1,c2,a2,b1=rprog,b2=rlab):
  aq=0
  for k,sh,_,_ in gw:h_(m,tb+gn[k],64);dk=ab+128*k;put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,a1 if k==0 else 0);put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,a2 if k==0 else 0);aq+=at(b1,b2,"tbld",dc(x0=tb+gn[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=0),m)
  return aq
 def rh(b5,regs,m):o={};at(ap,az,b5,dc(regs,sp=fK),m,out_regs=o,ls64=b0);return o
 A1,Bv,cd=I.rr(1,p//2),I.rr(1,p),I.rr(1,1<<63);A2=A1+cd;inv_e,ep=I.rr(1,p),sb+576;mc={};put(mc,fK+288,Bv);put(mc,ep-32,A2);mc[ep+304]=cd;mc[ep+312]=ep+256;oc=rh("t_c1b",dc(x3=ep+128,**fO(inv_e,("x25","x26","x27","x28"))),mc);w1=inv_e*(Bv-A2)%p;w2=(w1+cd*inv_e)%p;bday_c1b_ok=(get(mc,ep-64)%p,fL(oc,("x8","x9","x10","x11"))%p,oc["x29"])==(w1,w2,ep+256);ma,b7,b6,am={},I.rr(1,p),I.rr(1,p),aG(0);ae={};na=0;q=am[5]+544;ma[q]=ae[q]=ah
 for k,sh in en(gU):
  dk=ab+128*k;put(ma,dk,b7*pow(2,sh,p)%p);put(ma,dk+32,0);put(ma,dk+64,b6*pow(2,sh,p)%p);put(ma,dk+96,0)
  for a_ in rn(dk,dk+128,8):ae[a_]=ma[a_]
  na+=at(rprog,rlab,"tbld",dc(x0=am[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5))),ma);at(ar,au,"tbld",dc(x0=am[k],x1=dk,x2=64,x3=8,x4=2,x5=64,x6=3*(k in(2,5))),ae)
 axyx_ok=na==6*26+20*(4*128+2*16)and ma[q]==ah and ae[q]!=ah;mp,af,ag,aw,ax,ai,ad={},I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),I.rr(1,p),a6([5,3,1,6,7]);ay,gt=b9(mp,fT[0],af,ag,aw,ax),a3(ad);put(mp,fK+64,ai);put(mp,fK+96,ai)
 for k,v in en(gt[:4]):mp[fK+96+8*k]=fT[0]+v
 mp[fK+352],mp[fK+360]=fT[0]+gt[4],0;ol=rh("t_lkp",dc(x1=fK+64,x3=sb),mp);av=dc(mp);b9(av,fT[0],af,ag,aw,ax,ar,au);b3=rh("t_lkp",dc(x1=fK+64,x3=sb),av);a8=ay==5*26+20*5*16 and fL(ol,("x16","x17","x18","x19"))%p==(af*ad+ag+ai)%p and fL(b3,("x16","x17","x18","x19"))%p!=(af*ad+ag+ai)%p;PB,SK,he=0x600000,0x500000,True
 def tk(m,cl,hY,starts,steps,cO):
  for i,v in en(hL(CE.EIS[cl])):m[PB+8*i]=v
  for a in rn(SK-4096,SK,8):m.setdefault(a,0)
  for j in rn(6):put(m,ab+64*j,steps[j]);put(m,ab+64*j+32,starts[j])
  return at(fp,fl,"tbldk",dc(x0=hY,x1=ab,x2=cO,x3=PB,sp=SK),m)
 for cl in CE.EIS:me,us,ss={},[I.rr(p)for _ in rn(6)],[I.rr(p)for _ in rn(6)];n=tk(me,cl,sb,ss,us,3);rk=[CE.eis_red(CE.EIS[cl],(ss[j]+t*us[j])%p)for t in rn(3)for j in rn(6)];he&=n==CE.tbldk_count(3,cl)and[(get(me,sb+256*t+bA[j],3),get(me,sb+256*t+bA[j]+24,2))for t in rn(3)for j in rn(6)]==[(c1%2**192,c2)for c1,c2 in rk]and all(me[sb+256*t+o]==0 for t in rn(3)for o in(120,248))
 Xa=[I.rr(1,p)for _ in rn(6)];wv=[I.gb(256)for _ in rn(3)];bs=CE.eis_red(CE.EIS[EB],I.rr(p));mk={};ex={}
 for i,(cl,ce,cc)in en((("E",lambda a:2*a**3,lambda a:0),("F",lambda a:-2*a,lambda a:-a**4),("G",lambda a:BT[2]*a*a,lambda a:0))):
  for k,sh,mi,m in FW(i):t=wv[i]>>sh&m;tk(mk,CE.KCL[cl][k],gM[i]+DW(i)[k]+((t&mi)<<8),[(ce(Xa[j])*t*pow(2,sh,p)+(cc(Xa[j])if k==0 else 0))%p for j in rn(6)],[0]*6,1)
  for j in rn(6):ex[i,j]=bc(aY(cl,ce(Xa[j]),wv[i],cc(Xa[j]))+([bs]if i==0 else[]))
  for a,o in zip(bQ[i],bw(wv[i],i)):mk[fK+a]=o
 put(mk,fK+gC,bs[0]%2**192,3);put(mk,fK+gC+24,bs[1],2);m2=dc(mk);ea=gM[0]+g_[2]+((wv[0]>>gU[2]&hJ[2])<<8)+24;m2[ea]=(m2[ea]+1)%2**64
 def lines(m):return{(i,j):(get(m,sb+aA*j+(384,448,192)[i]+8*(j%3==1),3),get(m,sb+aA*j+(384,448,192)[i]+8*(j%3==1)+24,2))for i in rn(3)for j in rn(6)}
 oe=rh("t_efgk",dc(x3=sb+2752,x0=1,x2=2,x14=3,x15=4),mk);rh("t_efgk",dc(x3=sb+2752,x0=1,x2=2,x14=3,x15=4),m2);bH={k:(v[0]%2**192,v[1])for k,v in ex.items()};hP,bad=lines(mk),lines(m2);efgk_ok=hP==bH and bad[0,0]!=bH[0,0]and all(bad[k]==bH[k]for k in bH if k!=(0,0))and[oe[r]for r in(X0,"x2","x14","x15","x3","x4")]==[1,2,3,4,sb+2752,C0];ok=bday_c1b_ok and axyx_ok and a8 and he and efgk_ok;return dc(ok=ok,bday_c1b_ok=bday_c1b_ok,axyx_ok=axyx_ok,lookup_pair_ok=a8,tbldk_ok=he,efgk_ok=efgk_ok)
def hR(a8,K):
 I={};aj=aT("leaf");a4,Lb=a64.a64_load(aj);a6,a7=aI("secp_fe");rprog,rlab=aI("route");yg,ad,ap,av=cL(a6,a7,a8,K);I["y_generator"]=yg;Xs=ir(a8);ac,ar=machine.dd(int),{};af,aw=34*2**40+0x20000,34*2**40+0x40000
 for q in rn(4):put(ac,af+128*q,Xs[2*q][1]*Xs[2*q][1]%p);put(ac,af+128*q+32,Xs[2*q+1][1]*Xs[2*q+1][1]%p);put(ac,af+128*q+64,2*q+1,2);put(ac,af+128*q+80,2*q+2,6)
 at(rprog,rlab,"psel",dc(x0=af,x1=4,x2=aw,sp=34*2**40+0x18000),ac,out_regs=ar);I["bday_xs_psel_ok"]=ar[X0]==1 and(Xs[0][1]*Xs[0][1]%p)>>63!=(Xs[1][1]*Xs[1][1]%p>>63)and[get(ac,aw+128*q+32*s)for q in rn(4)for s in(0,1)]==[X[1]*X[1]%p for X in Xs];Ys=[cy(ap,i)for i in rn(12)];runs={}
 for q2 in(30,36,42):runs[q2]=gP(a4,Lb,Xs,Ys,q2,ad,K,rprog=rprog,rlab=rlab,chain=ap)
 am={(0,0):(aC(a8,(1<<40)-4,0)[0],0,a8),(0,1):(aC(a8,(1<<40)-3,1)[0],1,a8)};J=gP(a4,Lb,Xs,[cy(av,i)for i in rn(4)],24,ad,K,rprog=rprog,rlab=rlab,override=am,chain=av,ybase=fK+gI);I["yslot_lc1"]=J["records_ok"]and J["keys_formula"]and J["fb_form"]and J["sp_kept"]and J["reg_trace_ok"];I["registry_trace_ok"]=all(d["reg_trace_ok"]for d in runs.values())and J["reg_trace_ok"];I["runs"]={q:dc(fwd=d[fQ],bwd=d[fP],keys=d["keys"],fixup_extra=d["fixup_extra"],records_ok=d["records_ok"])for q,d in runs.items()};au={q:{k:d[k]-d["fixup_extra"][k]for k in(fQ,fP)}for q,d in runs.items()};I["fb_hold"]=all(d["fb_form"]and d["fixup_words_ok"]for d in runs.values());I["counts_ok"]=all(au[q][fQ]==K.fwd_count(q)and au[q][fP]==K.bwd_count(q)for q in runs);I["fixups"]=dc(sites=len(fixup_extra(a4)),flags_dead_at_every_join=cP(a4,Lb));al=[o for mn,o in a4 if mn in("tbz","tbnz")];I["loop_b34"]=[o for o in al if o[1]!="#12"]==[["x3","#37","1b"],["x5","#9","4f"]]+[["x3","#37","1b"]]*7;I["keys_hold"]=all(d["keys_formula"]for d in runs.values());I["frac_err"]=[[min(d["frac_err"][g][0]for d in runs.values()),max(d["frac_err"][g][1]for d in runs.values())]for g in(0,1)];I["all_correct"]=yg["mismatches"]==0 and yg["centres_ok"]and yg["counts_ok"]and I["bday_xs_psel_ok"]and I["yslot_lc1"]and all(d["acc_ok"]and d["frac_mis"]==0 and d["frac_out"]==0 and d["id_mis"]==0 and d["record_decode_mismatches"]==0 and d["records_ok"]and d["keys_noncanon"]==0 and d["wrong_index"]==0 and d["final_R_times_rho_is_one"]and d["sp_kept"]and d["reg_trace_ok"]for d in runs.values());az=Rd(SEED+49979687)
 for ay,aq,ao in(("window_batch",("Xm","diff"),"W"),("boundary_batch",(D,),D)):w=bG(a4,Lb,Xs,Ys,ad,K,az,aq,ao);I[ay]=dc(in_W=w["in_W"],in_D=w["in_D"],ok=w[{"W":"in_W",D:"in_D"}[ao]]>=4 and w["wrong_index"]==0 and w["keys_noncanon"]==0 and w["records"]==6*24 and w["records_ok"]and w["keys_formula"]and w["fb_form"]and w["fixup_words_ok"])
 gp,gl=a64.a64_load(aj.replace("b.lo 9f\n","cbnz xzr,9f\n").replace("b.hs 9f\n","cbnz xzr,9f\n").replace("cbhs \\a,\\s,","cbnz xzr,").replace("subs \\q,\\q,#-1","cmn \\q,xzr"));g=gP(gp,gl,Xs,Ys,24,ad,K,rprog=rprog,rlab=rlab,chain=ap);ag=sum(K.cT.values());ah=sum(K.FIX_EXTRA[m]*K.cT[m]for m in K.cT);I["forced_fixups"]=dc(fixup_extra=g["fixup_extra"],ok=g["frac_mis"]==0 and g["frac_out"]==0 and g["id_mis"]==0 and g["fb_form"]and g["fixup_words_ok"]and g["records_ok"]and sum(map(sum,g["fixups_3_4"].values()))==6*ag and sum(g["fixup_extra"].values())==6*ah);ai,ok=set(),True
 for _ in rn(3):a=a8.rr(1,p);N={};put(N,0x1000,a);ai.add(at(a6,a7,"fe_inv",dc(x0=0x2000,x1=0x1000,sp=0x100000),N,limit=10**6));ok&=get(N,0x2000)%p==pow(a,p-2,p)
 I["fe_inv"]=dc(counts=sorted(ai),correct=ok and ai=={K.INV});ak=("SHX","FDX","FS","LKFOLD4","DFX","ROWSG");counts={k:hQ(aj,k,ak+("FWD6_ONE","BWDPK01","BWDPK00","BWDPK10"))for k in(fQ,fP,"keys")};a9,aa=counts[fQ]["FWD6_ONE"],sum(counts[fP][k]for k in("BWDPK01","BWDPK00","BWDPK10"));a5={}
 for n in ak:a5[n]=(counts[fQ][n]//2+counts[fP][n])//3
 a5["BLK"],a5["DF"],a5["RW"]=a5.pop("LKFOLD4"),a5.pop("DFX"),a5.pop("ROWSG");I["site_counts"]=dc(per_x_quad=a5,pair_copies=dc(fwd=a9,bwd=aa));I["sites_cert"]=a5==K.cT and a9==4 and aa==6 and sum(counts["keys"].values())==0;ax=Rd(SEED+15485863);I["site_checks"]=sites.site_checks(sys.modules[__name__],ax,a6,a7);I["tbld_production"]=cQ(aj,Rd(SEED+27182818));b1=(an*Ys[0][0]%p,Ys[0][1]);b0,good=hT(a4,Lb,[Xs[0],b1,Xs[1],aE.neg(Ys[0])]+Xs[2:4],[Ys[0]]);I["d_zero_batch"]=dc(acc_zero=b0%p==0,fractions_within_bounds=good,all_other_fractions_wrong=good==0);ae,_=aI("route");ab=len(a4)+len(ae)+len(a6);I["prog_insns"]=dc(leaf=len(a4),route=len(ae),secp_fe=len(a6),total=ab,fits_2_18_bytes=4*ab<2**18,stated=f"{len(a4):,} / {len(ae):,} / {len(a6):,} ({ab:,} insns = {4*ab:,}"in K.d8);return I
def ik(o,I):c1=(2*(FA*o+(FA-FB)*I)+p)//(2*p);J=(2*((FA-FB)*o-FB*I)+p)//(2*p);return c1*FB+J*(FA-FB),c1*(FA-FB)-J*FA
def aC(o,d1,d2,v=None,j=None):
 while True:
  J=v
  if J is None:
   ha=o.gb(22)<<42|d2>>41<<41|o.rr(16)<<20|o.rr(200,(1<<20)-200);la=d1<<24|o.rr(200,(1<<24)-200);hb=o.gb(14)<<50|o.rr(2)<<20|o.rr(200,(1<<20)-200);lb=(d2&(1<<41)-1)<<23|o.rr(200,(1<<23)-200);I,uu=ik((ha<<64|la)*p>>128,(hb<<64|lb)*p>>128)
   if not(0<I<p and uu==an*I%p and go(I)[0]==I):continue
   rf=label(I)
   if not(rf[0]==d1 and rf[1]==d2 and(rf[3]&(1<<50)-1)<16):continue
   J=I
  m=o.rr(3)if j is None else j;x=J*pow(an,(3-m)%3,p)%p
  if go(x)[1]==m and not any(in_W(x,o)or in_D(x,o)for o in(0,2)):return x,J
def cI():
 I=io.StringIO()
 with contextlib.redirect_stdout(I):import certificate as K
 o=[" ".join(ln.split())for ln in I.getvalue().splitlines()];return K,[ln for ln in o if" ok"in ln or"FAIL"in ln or ln.startswith(("q1 ","E_total","declared:"))]
def cD(K):return dc(passed=not K.dg,failed=K.dg,checks_passed=len(gJ),tail_summary=gJ[-3:],claim="%d.%015d"%divmod(int(K.dD*10**15),10**15),q1=K.cV,per_leaf_declared=format(float(K.dH["per_leaf"]),".4f"),mu=float(K.cS["mu"]),opened_success_lower_bound=float(K.cZ),gate_prob=float(K.dy),success_lower_bound=float(K.df),per_leaf=float(K.dq),groups=K.cR)
def main():
 global gJ,CE;a7=Rd(SEED);K,gJ=cI();CE=K;OUT["input"]=cfg;I=OUT["A_identities"]=hU(a7);a4=(all(I[k]for k in("beta3","lam3","lamG_betax","no_x0","win_nz","label_constants"))and I["class_rule"]["mismatches_outside_window_set"]==0 and I["class_labels"]["ok"]and I["c132_ok"]);o=OUT["A_kernels"]=hR(a7,K);J=OUT["B_machine"]=machine.ca(sys.modules[__name__],K);OUT["counts_match_certificate"]=cm=(all(o[k]for k in("fb_hold","counts_ok","loop_b34","keys_hold","all_correct","sites_cert"))and o["fixups"]["flags_dead_at_every_join"]and(EH,EL)==tuple(map(tuple,(K.ERRH,K.ERRL)))and all(o[k]["ok"]for k in("window_batch","boundary_batch","forced_fixups","tbld_production"))and o["fe_inv"]["correct"]and o["site_checks"]["all_true"]and o["d_zero_batch"]["acc_zero"]and o["d_zero_batch"]["all_other_fractions_wrong"]and o["prog_insns"]["fits_2_18_bytes"]and o["prog_insns"]["stated"]and J["formulas_ok"]and all(J[k]["ok"]for k in("phase0","block_paths","grouped_paths","ls64","va55","encoding","phase4_round","y_ring","gate","wait_dispatch","rendezvous","keys_edges","registry_trace")));a8=OUT["D_certificate"]=cD(K);OUT["all_checks_pass"]=bool(a4 and cm and o["registry_trace_ok"]and a8["passed"]);OUT["seconds"]=round(time.time()-T0,2);json.dump(OUT,sys.stdout,separators=(",",": "));print()
 if not OUT["all_checks_pass"]:sys.exit(1)
if __name__=="__main__":main()
