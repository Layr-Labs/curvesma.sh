import contextlib,functools,io,json,os,random,sys,time;random.Random.rr,random.Random.gb=random.Random.randrange,random.Random.getrandbits;from collections import defaultdict;T0=time.time();cfg=json.load(sys.stdin);SEED=cfg.get("seed",1);BUDGET=cfg.get("seconds",12);RESERVE=min(14.5,BUDGET-0.3);T_END=T0+BUDGET-RESERVE;MIN_RANDOM,MIN_NEAR,MIN_CLASSES=300,100,300
def package_root():
 here=os.path.dirname(os.path.abspath(__file__))
 for base in(here,os.path.join(here,"files"),"/snapshot/files"):
  if os.path.isfile(os.path.join(base,"kernels","a64.py")):return base
 raise FileNotFoundError("kernels/a64.py not found")
HERE=package_root();sys.path.insert(0,os.path.join(HERE,"kernels"));import a64;a64.load=functools.lru_cache(maxsize=None)(a64.load);sys.path.insert(0,HERE);import sites,machine;M64=2**64-1;M41,M42=2**40-1,2**41-1;PART=105;M128=2**128;p=2**256-2**32-977;r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",16);GX=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",16);GY=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8",16);BETA=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",16);LAM=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72",16);C0=2**32+977;FB=int("e4437ed6010e88286f547fa90abfe4c3",16);FA=int("114ca50f7a8e2f3f657c1108d9d44cfd9",16);FRAC_OFF=((64,80,128,144,32,48),(96,112,160,176,192,208));ENTRY,ACC,KEYS0,MX=512,32,64,6;YR_B,YR_B2,YR_YROW,YR_NT=288,64,128,320;W6S=(0,52,103,154,205);W6M=tuple((1<<b)-1 for b in(52,51,51,51,51));W6=tuple(enumerate(zip(W6S,W6M)));BT=(1,BETA,BETA*BETA%p);O,KO=(2,2,0),(1,0,2);MH_A,MH_B,ML=20,21,23;OUT,CERT_LINES,CE={},[],None
class Curve:
 def __init__(ac,q,order):ac.q,ac.r=q,order
 def add(ac,P,Q):
  q=ac.q
  if P is None:return Q
  if Q is None:return P
  if P[0]==Q[0]:
   if(P[1]+Q[1])%q==0:return None
   bi=3*P[0]*P[0]*pow(2*P[1],-1,q)%q
  else:bi=(Q[1]-P[1])*pow(Q[0]-P[0],-1,q)%q
  x=(bi*bi-P[0]-Q[0])%q;return(x,(bi*(P[0]-x)-P[1])%q)
 def neg(ac,P):return None if P is None else(P[0],(-P[1])%ac.q)
 def mul(ac,k,P):
  k%=ac.r
  if not k or P is None:return None
  q,X,Y,Z=ac.q,P[0],P[1],1
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
     else:H,R=(U2-X)%q,(S2-Y)%q;H2=H*H%q;H3=H*H2%q;V=X*H2%q;X=(R*R-H3-2*V)%q;Y=(R*(V-X)-Y*H3)%q;Z=Z*H%q
  if X is None or Z%q==0:return None
  iz=pow(Z,-1,q);fa=iz*iz%q;return(X*fa%q,Y*fa*iz%q)
PROD=Curve(p,r);GP=(GX,GY)
def limbs(v,n=4):return[(v>>(64*k))&M64 for k in range(n)]
def put(mem,a,v,n=4):
 for k,w in enumerate(limbs(v,n)):mem[a+8*k]=w
def get(mem,a,n=4):return sum(mem[a+8*k]<<(64*k)for k in range(n))
def R4(a):return tuple("x%d"%i for i in range(a,a+4))
def rgi(regs,ax):return sum(regs[n]<<(64*k)for k,n in enumerate(ax))
def ars(v,ax):return{n:w for n,w in zip(ax,limbs(v,len(ax)))}
def exk(x):
 v=[x,BETA*x%p,BETA*BETA%p*x%p];j=min(range(3),key=lambda t:v[t]);return v[j],j
def t128(v):return(v<<128)//p
def t64(v):return(v<<64)//p
def ws(*v):return sum(x>>64 for x in v)%2**64<<64|sum(x&M64 for x in v)%2**64
def wneg(x):return ws(x^(M128-1),1<<64|1)
def rfields(a,b):ha,la,hb,lb=a>>64,a&M64,b>>64,b&M64;d1,d2=la>>23,lb>>23;return d1,d2,(ha>>41<<41)|d2,(hb>>21<<21)|((ha>>20)&((1<<21)-1))
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
def sqv(a):z=a*a;return z%2**256+(z>>256)*C0
def sqw(a):return fold(sqv(a))
def msum(fe,w,fd=0):return sum((fe*((w>>sh)&m)*pow(2,sh,p)+(fd if k==0 else 0))%p for k,(sh,m)in W6)
def rowsf(w):return(w+((w>>192)*C0>>64))%2**256
def frac_e(r):return r+((r*C0+2**255)>>256)
def blkf(w):return sum(frac_e(BETA*((w>>sh)&m)*pow(2,sh,p)%p)for sh,m in zip(W6S,W6M))%2**256
def ktfp(c1,c2,w,y1):
 y0=rowsf(w);h=lambda a,b:a*b>>64;m,n,d,a,f,g=limbs(c1%2**192,3),limbs(c2,2),limbs((c1-c2)%2**192,3),limbs(y0),limbs(y1),limbs((y0+y1)%2**256)
 def pr(x,y):t=len(x)>2;return(h(x[0],y[2])+x[0]*y[3]+h(x[1],y[1])+x[1]*y[2]+(x[2]*y[1]if t else 0))&M64,(h(x[0],y[1])+x[0]*y[2]+h(x[1],y[0])+x[1]*y[1]+(x[2]*y[0]if t else 0))&M64
 AH,AL=pr(m,a);BH,BL=pr(d,g);P=pr(n,f);return(AH+P[0]&M64)<<64|(AL+P[1]&M64),((BH-AH)&M64)<<64|(BL-AL)&M64
def ksum(z):return sum(x[0]for x in z),sum(x[1]for x in z)
def ered(c,v):return CE.eis_red(CE.EIS[c],v)
def kents(cl,ce,v,cc):return[ered(cl,(ce*((v>>sh)&m)*pow(2,sh,p)+(cc if k==0 else 0))%p)for k,(sh,m)in W6]
def kcomp(a,b):return[ksum(kents("E",2*a**3,b,0)+[ered("EB",pow(b,4,p))]),ksum(kents("FG",-2*a,pow(b,3,p),-a**4)),ksum(kents("FG",BT[2]*a*a,3*b*b%p,0))]
def fracs(mem,e,s):return[get(mem,e+o,2)for o in FRAC_OFF[s]]
def components(X,Y):
 a,b=X[0],Y[0];E,F,G=(b**4+2*a**3*b)%p,(a**4+2*a*b**3)%p,3*a*a*b*b%p;out=[]
 for w in(pow(Y[1]+X[1],-2,p),pow(Y[1]-X[1],-2,p)):out+=[E*w%p,F*w%p,G*w%p]
 return out
def werr(t,v):return[((t>>64*k)-(v>>64*k)+2**63)%2**64-2**63 for k in(1,0)]
def frac_errs(mem,e,X,Y):
 cs,out=components(X,Y),[]
 for s in(0,1):ex=((cs[3*s]-Y[0])%p,-cs[3*s+1]%p,BT[2]*cs[3*s+2]%p);out+=[werr(t128(BT[i%2]*ex[i//2]%p),v)for i,v in enumerate(fracs(mem,e,s))]
 return out
def leaf_w(f6,t):b0,b1,a0,a1,g0,g1=f6;n=wneg;return((ws(a0,a1,g0,n(b0),n(b1)),ws(b0,g1,n(a0))),(ws(n(a0),n(b0),n(b1),n(g0),n(g1)),ws(b0,g0,n(a1))),(ws(a0,a1,b0,n(g0),n(g1)),ws(b1,g0,n(a0))))[t]
def record_fields(w):f0,f1=w[0]>>64,w[1]>>64;s=(f0+f1)&M64;j=(2 if f1+s>M64 else 1)if f0>=f1 else(2 if f0+s>M64 else 0);n=ws(w[0],w[1])^(M128-1);return(j,)+rfields(*((w[0],w[1]),(w[1],n),(n,w[0]))[j])
def in_W(x,o=0):
 x=BT[o]*x%p;v=[x<<64,(BETA*x%p)<<64,(BETA*BETA%p*x%p)<<64];top=p<<64;lo,hi=min(ERR_H[0][0],ERR_H[1][0]),max(ERR_H[0][1],ERR_H[1][1]);return(any(w<hi*p or w>=top+lo*p for w in v[:2])or v[2]<(1-2*lo)*p or v[2]>top-(3+2*hi)*p or abs(v[0]-v[1])<(hi+1-lo)*p or any(3*lo*p<=w-v[2]<(3*hi+3)*p for w in v[:2]))
def in_D(x):v=(x,BETA*x%p,BETA*BETA%p*x%p);return any(u%2**mb<hi or u%2**mb>=2**mb+lo for w,eh,el in zip(v,ERR_H,ERR_L)for u,mb,(lo,hi)in((t128(w)>>64,MH_A,eh),(t128(w)&M64,ML,el)))
def dft_leaf(c3,t,a,b):return(c3[0]-b+pow(BETA,t,p)*(c3[1]-a)+pow(BETA,2*t,p)*c3[2])%p
def inverse_dft(L):i3,b2=pow(3,-1,p),BETA*BETA%p;return((L[0]+L[1]+L[2])*i3%p,(L[0]+b2*L[1]+BETA*L[2])*i3%p,(L[0]+BETA*L[1]+b2*L[2])*i3%p)
def weak(rg):return rg.rr(p)if rg.random()<0.85 else rg.rr(C0)+p
def represent(x,t,a,b,rg,raw=False):
 w,c,k=weak(rg)or 1,[rg.rr(p)for _ in range(3)],rg.rr(3);c[k]=0;v=(x+b+BT[t]*a)*pow(w,-1,p)-sum(BT[t*i%3]*c[i]for i in range(3));c[k]=v*pow(BT[t*k%3],-1,p)%p;y1=blkf(w);out=[]
 for i,(ci,y)in enumerate(zip((c[0],-c[1],BT[2]*c[2]),(b,None if raw else a,None))):cl="E"if i==0 else"FG";ks=[rg.rr(p)for _ in range(5-(i>0))];A,B=ktfp(*ksum([ered(cl,u)for u in ks]+[ered("EB"if i==0 else cl,(ci-sum(ks))%p)]),w,y1);cw=(corr(y)if i<1 else(t128(y),t128(BETA*y%p)))if y else(0,0);out+=[ws(A,cw[0]),ws(B,cw[1])]
 return out,((c[0]*w-b)%p,(c[1]*w-a)%p,c[2]*w%p)
def edge_leaf(rg,kind):
 m=rg.rr(3)
 if kind=="D":b=rg.rr(2);mb=(MH_A,ML)[b];u=rg.gb(64-mb)<<mb|rg.choice((rg.rr(70),2**mb-1-rg.rr(70),rg.rr(2**mb)));v=-(-(u<<64|rg.gb(64)if b==0 else rg.gb(64)<<64|u)*p//M128);return v*pow(BETA,-m,p)%p or 1
 off=rg.rr(160<<192);w=off if rg.random()<0.5 else p-1-off
 if kind=="Xm":return w*pow(BETA,-m,p)%p or 1
 return w*pow((1-BETA,1-BETA*BETA,BETA-BETA*BETA)[m]%p,-1,p)%p or 1
def deadline(fc):dp=time.time();return dp+fc*max(0.0,T_END-dp)
def near_tie_x(rg):mlt=rg.choice([1,BETA,BETA*BETA,1-BETA,1-BETA*BETA,BETA-BETA*BETA]);dd=rg.rr(-(1<<rg.rr(150,250)),1<<rg.rr(150,250));return dd*pow(mlt%p,-1,p)%p or 1
def index_check(x,rg):t,a,b=rg.rr(3),rg.rr(1,p),rg.rr(1,p);j=record_fields(leaf_w(represent(x,t,a,b,rg)[0],t))[0];return in_W(x,O[t]),(j+O[t])%3!=exk(x)[1]
def class_labels(rg):
 out=dict(classes=0,edge_classes=0,members=[0,0,0],in_W=0,in_D=0,wrong_index=0,mismatches=0,errors_h=[[999,-999]for _ in range(3)],errors_l=[[999,-999]for _ in range(3)]);dl=deadline(0.05)
 while out["classes"]<MIN_CLASSES or(out["classes"]<4*MIN_CLASSES and time.time()<dl):
  i0=out["classes"]%3==2;x0=edge_leaf(rg,("D","Xm","diff")[out["classes"]//3%3])if i0 else rg.rr(1,p);v=exk(x0)[0];out["classes"]+=1;out["edge_classes"]+=i0;lab=label(v)
  for m in range(3):
   x=v*pow(BETA,(3-m)%3,p)%p;t,a,b=rg.rr(3),rg.rr(1,p),rg.rr(1,p);w=leaf_w(represent(x,t,a,b,rg)[0],t)
   for k,c in enumerate((w[0],w[1],ws(w[0],w[1])^(M128-1))):eh,el=werr(t128(BT[(O[t]+k)%3]*x%p),c);out["errors_h"][k]=[min(out["errors_h"][k][0],eh),max(out["errors_h"][k][1],eh)];out["errors_l"][k]=[min(out["errors_l"][k][0],el),max(out["errors_l"][k][1],el)]
   f=record_fields(w)
   if in_W(x,O[t]):out["in_W"]+=1;continue
   if(f[0]+O[t])%3!=m:out["wrong_index"]+=1;continue
   out["members"][m]+=1
   if in_D(x):out["in_D"]+=1
   elif f[1:]!=lab:out["mismatches"]+=1
 e6=(out["mismatches"]==0 and out["wrong_index"]==0 and out["in_D"]>0 and out["in_W"]>0 and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors_h"],ERR_H))and all(lo<=e[0]and e[1]<=hi for e,(lo,hi)in zip(out["errors_l"],ERR_L)));return dict(ok=e6,classes=out["classes"],in_W=out["in_W"],in_D=out["in_D"],errors_h=out["errors_h"],errors_l=out["errors_l"])
def part_a_identities(J):
 aa={};aa["beta_cube"]=pow(BETA,3,p)==1 and BETA!=1;aa["lambda_cube"]=pow(LAM,3,r)==1 and LAM!=1;aa["lambda_G_equals_beta_x"]=PROD.mul(LAM,GP)==(BETA*GX%p,GY);aa["no_x0_point"]=pow(7,(p-1)//2,p)==p-1;aa["win_nonzero"]=all(c%p for c in(1-BETA,1-BETA*BETA,BETA-BETA*BETA,3));aa["label_constants"]=(FA*FA-FA*FB+FB*FB==p and FB*BETA%p==FA-FB and(FA-FB)*BETA%p==p-FA);bad=drops=wrong=g9=0;ay={"random":0,"near_tie":0};dl=deadline(0.10)
 while g9<26000 and(g9<MIN_RANDOM or time.time()<dl):x,kind=(J.rr(1,p),"random")if g9<20000 else(near_tie_x(J),"near_tie");ay[kind]+=1;g9+=1;ah,ai=index_check(x,J);drops+=ah;wrong+=ah and ai;bad+=ai and not ah
 extra=max(0,MIN_NEAR-ay["near_tie"]);i1=random.Random(SEED+7919)
 for _ in range(extra):ah,ai=index_check(near_tie_x(i1),i1);ay["near_tie"]+=1;drops+=ah;wrong+=ah and ai;bad+=ai and not ah
 aa["class_rule"]=dict(checked=ay,in_W=drops,mismatches_outside_window_set=bad);aa["class_labels"]=class_labels(random.Random(SEED+104729));aw=random.Random(SEED+918273);ag,c132_ok=1<<64,True
 def _in_r(i_,j_,d_):bt=max(0,i_-d_);return 0<=i_<ag and((bt<=j_<=min(ag-1,ag+d_-2*i_))or(max(bt,ag-d_-i_)<=j_<=min(ag-1,2*ag+d_-2*i_)))
 for dq in range(600):cv,cx=aw.choice([0,1,ag-1,ag//2,aw.gb(64)]),aw.choice([0,1,ag-1,aw.gb(64)]);dr=record_fields((cv<<64,cx<<64))[0];dt=(cv,cx,(~((cv+cx)%ag))&M64);c132_ok&=_in_r(dt[dr],dt[(dr+1)%3],0);du=near_tie_x(aw)if dq%3==0 else edge_leaf(aw,("Xm","diff","D")[dq%3]);t_=aw.rr(3);jx_=record_fields(leaf_w(represent(du,t_,aw.rr(1,p),aw.rr(1,p),aw)[0],t_))[0]+O[t_];dv=pow(BETA,jx_,p)*du%p;bu,bv=t64(dv),t64(BETA*dv%p);w_=CE.WRAPT_C;c132_ok&=(_in_r(bu,bv,CE.PAD_C)or min(bu,bv,ag-1-bu,ag-1-bv)<w_ or min(abs(bu+bv-ag),abs(bu+bv-2*ag))<w_-2)
 aa["c132_ok"]=c132_ok;return aa
SNT=0x5E5E5E5E5E5E5E5E;TIB=1<<40;A_SP=(34<<40)+0x7000000;EB_BASE=(1<<45)+(1<<40);EB_END=EB_BASE+(1<<34);YR0,YR_STRIDE=(1<<45)+(1<<31),1<<31;BT0,BT0_END,BT1=1<<58,2*(1<<58),0xfc00000000000000;AXYX0,AXYX_HI=3*(1<<58),9*(1<<58);XT_END=63*(1<<58);D6O=tuple(c<<57 for c in(0,2,3,4,5));D7O=tuple(c<<57 for c in(0,4,6,8,10));XTK=tuple(c<<57 for c in(42,54,66,78,90,102));XTP=(114<<57,120<<57)
def kf(n):return open(os.path.join(HERE,"kernels",n)).read()
def ax_bases(j):return(AXYX0+(j<<58),AXYX_HI+((4*j+0)<<57),AXYX_HI+((4*j+1)<<57),AXYX_HI+((4*j+2)<<57),AXYX_HI+((4*j+3)<<57))
def offs6(v):return[D6O[k]+(((v>>sh)&m)<<6)for k,(sh,m)in W6]
def offs7(v):return[D7O[k]+(((v>>sh)&m)<<7)for k,(sh,m)in W6]
class _TraceLC1Dict(dict):
 def __init__(ac,rprog,rlab):super().__init__();ac._rp,ac._rl,ac.touched=rprog,rlab,set()
 def __getitem__(ac,a):ac.touched.add(a);return super().__getitem__(a)
 def __setitem__(ac,a,v):ac.touched.add(a);return super().__setitem__(a,v)
 def get(ac,a,d=None):
  if a in ac:ac.touched.add(a)
  return super().get(a,d)
 def __missing__(ac,a):
  if 0<=a<16*machine.H1 and a%8==0:d=a//8;b=machine.blk1(d);a64.run(ac._rp,ac._rl,"init3",dict(x0=1,x1=a,x2=b+32,x3=machine.B1,x30=SNT),ac);ac.touched.add(a);return super().__getitem__(a)
  for lo,hi,first,stride in machine.LC1_TABLES[2:]:
   if lo<=a<hi:ac.touched.add(a);super().__setitem__(a,0);return 0
  raise KeyError(hex(a))
def pop60(rprog,rlab,mem,tb,c1,bj,c2,cj,w,bs=None):
 ab=(34<<40)+0x180000
 for k,(sh,m)in W6:
  tk=(w>>sh)&m;ad=(bs[k]if bs else tb+D6O[k])+(tk<<6)
  if ad in mem:continue
  put(mem,ab,0);put(mem,ab+32,(c1*tk*pow(2,sh,p)+(bj if k==0 else 0))%p);put(mem,ab+64,0);put(mem,ab+96,(c2*tk*pow(2,sh,p)+(cj if k==0 else 0))%p);a64.run(rprog,rlab,"tbld",dict(x0=ad,x1=ab,x2=1,x3=1,x4=2,x5=64,x6=0,x30=SNT),mem)
def cube_root(v):c_=pow(v,(p+2)//9,p);return c_ if c_*c_*c_%p==v%p else None
def make_bday_xs(rg):
 xs=[]
 for q in range(4):
  while True:
   X0=PROD.mul(rg.rr(1,r),GP);A0=X0[1]*X0[1]%p;lo=A0&((1<<63)-1)
   if q==0 and lo<(1<<62):continue
   d=rg.rr((1<<63)-lo,1<<63)if q==0 else rg.rr(1,1<<63);A1=A0+d
   if A1>=p:continue
   y1,a1=sqrt_mod_p(A1),cube_root((A1-7)%p)
   if y1 is not None and a1 is not None:xs+=[X0,(a1,y1)];break
 return xs
def popax(rprog,rlab,mem,j,X,w):pop60(rprog,rlab,mem,0,2*(X[1]*X[1]%p)%p,0,2*X[1]%p,0,w,ax_bases(j))
def popbt(rprog,rlab,mem,w):
 a=BT0+((w&W6M[0])<<6)
 if a not in mem:put(mem,a,frac_e(BETA*(w&W6M[0])%p));put(mem,a+32,0)
 for i in range(1,5):
  t=(w>>W6S[i])&W6M[i];a=BT1+(t<<7)
  if a not in mem:
   for k in range(4):put(mem,a+32*k,frac_e(BETA*t*pow(2,W6S[k+1],p)%p))
def eis_words(c):return limbs(c["W"][1])+limbs(p-c["U"][1])+limbs(c["kS"])+limbs(c["kT"])+limbs(pow(p,-1,2**192),3)+limbs(c["U"][0]%2**192,3)+limbs(c["W"][0]%2**192,3)+limbs(c["U"][1],2)+limbs(c["W"][1],2)+limbs(c["o"][0]%2**192,3)+limbs(c["o"][1],2)+limbs(c["C"],3)+limbs(c["X2"],2)+[len(c["c"])]+sum((limbs(l1%2**192,3)+limbs(l2%2**128,2)for l1,l2 in c["c"]),[])
def popk(mem,Xs,b):
 for i,(v,cl,ce,cc)in enumerate(((b,"E",lambda a:2*a**3,lambda a:0),(pow(b,3,p),"FG",lambda a:-2*a,lambda a:-a**4),(3*b*b%p,"FG",lambda a:BT[2]*a*a,lambda a:0))):
  for g in(0,1):
   for k,(sh,m)in W6:
    t=(v>>sh)&m;ad=XTK[2*i+g]+D7O[k]+(t<<7)
    if ad in mem:continue
    for j in range(3):a=Xs[3*g+j][0];c1,c2=ered(cl,(ce(a)*t*pow(2,sh,p)+(cc(a)if k==0 else 0))%p);put(mem,ad+40*j,c1%2**192,3);put(mem,ad+40*j+24,c2,2)
    mem[ad+120]=0
def popbat(rprog,rlab,mem,Xs,Ys,rho):
 ap=rho;bk,ck=[],[];A=[X[1]*X[1]%p for X in Xs[:MX]];sA=[(A[2*q]+A[2*q+1])%p for q in range(3)];pA=[A[2*q]*A[2*q+1]%p for q in range(3)]
 for yi,Y in enumerate(Ys):
  b,yy=Y;B=yy*yy%p;b3=pow(b,3,p);popk(mem,Xs,b);pop60(rprog,rlab,mem,XTP[0],(-sA[1])%p,(pA[1]-7*sA[1])%p,(-sA[0])%p,(pA[0]-7*sA[0])%p,b3);pop60(rprog,rlab,mem,XTP[1],0,0,(-sA[2])%p,(pA[2]-7*sA[2])%p,b3)
  for pr in range(3):i2=lkfold(sum((((-sA[pr])%p*((b3>>sh)&m)*pow(2,sh,p)+((pA[pr]-7*sA[pr])%p if k==0 else 0))%p)for k,(sh,m)in W6)+(B*B%p));bk.append(ap);ck.append(i2);ap=fold(sum(((i2>>64*k)&M64)*shr(ap)[k]for k in range(4)))
 R=pow(ap%p,p-2,p)if ap%p else 0;prod=lambda m_,rr:sum(((m_>>64*k)&M64)*rr[k]for k in range(4))
 for jp in range(len(bk)-1,-1,-1):
  yi=(2*jp)//MX;b,yy=Ys[yi];B=yy*yy%p;pr=jp%3;(a1,y1),(a2,y2)=Xs[2*pr],Xs[2*pr+1];A1,A2=y1*y1%p,y2*y2%p;dA=(A2-A1)%p;RR=shr(R);inv_e,R=fold(prod(bk[jp],RR)),fold(prod(ck[jp],RR));g7=fold(prod(inv_e,shr(fs_op(B,A2))));h5=lkfold(g7+dA*inv_e)
  for j,X,aq in((2*pr+1,(a2,y2),h5),(2*pr,(a1,y1),g7)):cy=sqw(aq);popax(rprog,rlab,mem,j,X,cy);Mv=lkfold(aq+msum(2*(X[1]*X[1]%p)%p,cy));Pv=lkfold(msum(2*X[1]%p,cy));Nv=fold(sum(((Pv>>64*k)&M64)*(yy*pow(2,64*k,p)%p)for k in range(4)));popbt(rprog,rlab,mem,fs_op(Mv,Nv));popbt(rprog,rlab,mem,faddw_op(Mv,Nv))
def x_words(mem,e,Xs,l):
 j=l%MX;X=Xs[j]
 for k,w in enumerate((t128(X[0]),t128(BETA*X[0]%p))):put(mem,e+256+16*k,w,2)
 put(mem,e+288,X[1]*X[1]%p)
 for k,b_ in enumerate(ax_bases(j)+(0x5A5A5A5A5A5A5A5A,)*3):mem[e+320+8*k]=b_
 if j&1:mem[e+368]=(Xs[j][1]**2-Xs[j-1][1]**2)%p
def put_yrec(mem,e,Y):
 b,y=Y;B=y*y%p;c1,c2=ered("EB",pow(b,4,p));put(mem,e,c1%2**192,3);put(mem,e+24,c2,2);put(mem,e+YR_B,B);put(mem,e+YR_B2,B*B%p)
 for k in range(4):put(mem,e+YR_YROW+32*k,y*pow(2,64*k,p)%p)
 put(mem,e+YR_NT,neg_t(b),2);put(mem,e+YR_NT+16,neg_t(BETA*b),2);az=offs6(pow(b,3,p))
 for k in range(4):mem[e+96+8*k]=az[k]
 mem[e+352]=az[4]
 for v,ads in((b,(384,392,400,408,416)),(pow(b,3,p),(480,488,496,504,424)),(3*b*b%p,(440,448,456,464,472))):
  for a,o in zip(ads,offs7(v)):mem[e+a]=o
def bwd_model(Xs,Ys,mem,bk,ck,R,ybase=None):
 if ybase is None:ybase=YR0
 out,prod={},lambda m,rr:sum((m>>64*k&M64)*rr[k]for k in range(4))
 for jp in range(len(bk)-1,-1,-1):
  l1,l2=2*jp,2*jp+1;(a1,cz),(a2,da)=Xs[l1%MX],Xs[l2%MX];b,yy=Ys[l2//MX];B,A1,A2=yy*yy%p,cz*cz%p,da*da%p;RR=shr(R);inv_e,R=fold(prod(bk[jp],RR)),fold(prod(ck[jp],RR));D2=fs_op(B,A2);dA=(A2-A1)%p;g7=fold(prod(inv_e,shr(D2)));h5=lkfold(g7+dA*inv_e)
  for l,(a,yx),aq in((l2,(a2,da),h5),(l1,(a1,cz),g7)):
   bz=Ys[l//MX][0];AX=yx*yx%p;yr=ybase+512*(l//MX);e7=[get(mem,yr+YR_YROW+32*k)for k in range(4)];b7=sqw(aq);i4=lkfold(aq+msum(2*AX%p,b7));j2=lkfold(msum(2*yx%p,b7));i5=fold(prod(j2,e7));cs=kcomp(a,bz);fr=[]
   for w in(fs_op(i4,i5),faddw_op(i4,i5)):
    y1=blkf(w)
    for i,c in enumerate(cs):A,Bv=ktfp(*c,w,y1);cw=corr(bz)if i==0 else(0,0);fr+=[ws(A,cw[0]),ws(Bv,cw[1])]
   out[l]=fr
 return out,R
def six_leaves(X,Y):
 out=[]
 for t in range(3):PX=(pow(BETA,t,p)*X[0]%p,X[1]);out+=[PROD.add(PX,Y)[0],PROD.add(PX,PROD.neg(Y))[0]]
 return out
def fixup_extra(ae):
 out,i6={},-1
 for i,(mn,_)in enumerate(ae):
  if mn in("b.lo","b.hs","cbnz"):i6=i
  elif mn=="stadd":out[i]=i-i6
 return out
def flags_dead_after_fixups(ae,labels):
 e8=("adc","adcs","sbc","sbcs","csel","cinc","cset","b.eq","b.ne","b.lo","b.hs","b.lt","b.ge")
 def target(t,pc):
  if t[:-1].isdigit()and t[-1:]in("f","b"):
   c_=labels[t[:-1]];return[q for q in c_ if q<pc][-1]if t[-1]=="b"else[q for q in c_ if q>=pc][0]
  return labels[t][0]
 for i in fixup_extra(ae):
  h0,seen=[i+1],set()
  while h0:
   j=h0.pop()
   if j in seen:continue
   seen.add(j);mn,o=ae[j]
   if mn in e8:return False
   if mn in("adds","subs","cmp","cmn","tst","ret"):continue
   if mn in("b","tbz","tbnz","cbz","cbnz"):
    h0.append(target(o[-1],j+1))
    if mn=="b":continue
   h0.append(j+1)
 return True
CBS,CB0,YSP=2304,(34<<40)+0x2000000,(34<<40)+0x4000000
def chain_n(i):rho,q=divmod(i,33);return 33*rho+(0 if q==32 else(16-q//2)*(1-2*(q%2)))
@functools.lru_cache(maxsize=None)
def chain_y(ch,i):return PROD.add(ch[0],PROD.mul(chain_n(i)%r,ch[1]))
def y_generator(bh,ca,J,CT):
 mem,chains={},[]
 for c in range(2):
  C,S=PROD.mul(J.rr(1,r),GP),PROD.mul(J.rr(1,r),GP);chains.append((C,S));cb=CB0+CBS*c
  for i,P in[(0,C)]+[(i,PROD.mul(i,S))for i in list(range(1,17))+[33]]:put(mem,cb+(1088 if i==33 else 64*i),P[0]);put(mem,cb+(1088 if i==33 else 64*i)+32,P[1])
  mem[cb+2240]=YR0+c*YR_STRIDE
 put(mem,CB0-128,1)
 for a in range(YSP-4096,YSP+1024,8):mem[a]=0
 put(mem,YSP+352,BETA);mem[YSP+384],mem[YSP+424]=C0,CB0;counts=[]
 mem[YSP+1016]=1
 for k,v in enumerate(eis_words(CE.EIS["EB"])):mem[YSP+1024+8*k]=v
 for nch in(2,1):mem[YSP+416]=CB0+CBS*nch;counts.append(a64.run(bh,ca,"yround",dict(sp=YSP,x30=SNT),mem,limit=10**7))
 bad=i7=0
 for c,rd in((0,0),(1,0),(0,1)):
  for q in range(33):i=33*rd+q;(b,y),e=chain_y(chains[c],i),YR0+c*YR_STRIDE+512*i;B=y*y%p;c1,c2=ered("EB",pow(b,4,p));ok=(get(mem,e,3)==c1%2**192 and get(mem,e+24,2)==c2 and get(mem,e+YR_B)==B and get(mem,e+YR_B2)==B*B%p);ok&=all(get(mem,e+YR_YROW+32*k)==y*pow(2,64*k,p)%p for k in range(4));ok&=get(mem,e+YR_NT,2)==neg_t(b)and get(mem,e+YR_NT+16,2)==neg_t(BETA*b);ok&=([mem[e+96+8*k]for k in range(4)]+[mem[e+352]])==offs6(pow(b,3,p))and[mem[e+a]for a in(384,392,400,408,416)]==offs7(b)and[mem[e+a]for a in(480,488,496,504,424)]==offs7(pow(b,3,p))and[mem[e+a]for a in(440,448,456,464,472)]==offs7(3*b*b%p);bad+=not ok;i7+=1
 cen=[(get(mem,CB0+CBS*c),get(mem,CB0+CBS*c+32))for c in range(2)];ok=bad==0 and cen==[chain_y(chains[0],33*2+32),chain_y(chains[1],33+32)];mem[YSP+1016]=0;counts.append(a64.run(bh,ca,"yround",dict(sp=YSP,x30=SNT),mem,limit=10**7));cok=counts==[CT.YR_FIX+2*CT.YR_BLK,CT.YR_FIX+CT.YR_BLK,CT.YR_FIX+CT.YR_BLK_MINE];return dict(records=i7,mismatches=bad,centres_ok=ok,counts=counts,counts_ok=cok),mem,chains[0],chains[1]
def leaf_batch(ae,labels,Xs,Ys,q2,bd,CT,rprog=None,rlab=None,W2=0x123456789,KN=5<<26,override=None,chain=None,rho=None,ybase=YR0):
 assert q2%24==0 and ybase%YR_STRIDE==0
 if rprog is None:rprog,rlab=a64.load(kf("route.s"))
 fx=fixup_extra(ae);ce={k:dict.fromkeys(fx,0)for k in("fwd","bwd")};nY,SP=q2//MX,A_SP;sb=EB_END-ENTRY*q2;LS=([(AXYX0,XT_END),(BT0,BT0_END),(BT1,1<<64),(ybase,ybase+512*nY),(sb,EB_END)],[(sb,EB_END)]);mem=_TraceLC1Dict(rprog,rlab);rho=rho or random.Random(q2).rr(1,p);popbat(rprog,rlab,mem,Xs,Ys[:nY],rho)
 if bd is not None:
  for a in range(ybase-64,ybase+512*nY,8):mem[a]=bd.get(a,0)
 else:
  for a in range(ybase-64,ybase,8):mem[a]=0
  for yi in range(nY):put_yrec(mem,ybase+512*yi,Ys[yi])
 for l in range(q2):x_words(mem,sb+ENTRY*l,Xs,l)
 regs=dict(x0=q2,x1=ybase+64,x3=sb,sp=SP,x30=SNT);regs.update(ars(rho,R4(24)));out={};nf=a64.run(ae,labels,"fwd",regs,mem,out_regs=out,visits=ce["fwd"],ls64=LS);words,ap=[out["x0"]],rgi(out,R4(24));j3=[get(mem,sb+ENTRY*(2*j)+ACC)for j in range(q2//2)];es=[get(mem,sb+ENTRY*(2*j))for j in range(q2//2)];R0=pow(ap,p-2,p);top=ybase+512*(nY-1);regs=dict(x3=EB_END,x1=(34<<40)+0x10008,sp=top,x30=SNT);regs.update(ars(R0,("x0","x2","x14","x15")));out={};nb=a64.run(ae,labels,"bwd",regs,mem,out_regs=out,visits=ce["bwd"],ls64=LS);words.append(out["x5"]);h7=rho
 for l in range(q2):h7=h7*(Ys[l//MX][1]**2-Xs[l%MX][1]**2)%p
 b9=rgi(out,("x0","x2","x14","x15"));sp_kept=out["sp"]==top;model,mR=bwd_model(Xs,Ys,mem,j3,es,R0,ybase);d3,d5,d6,ai,leaves=int(mR!=b9),0,0,[[999,-999],[999,-999]],[]
 for l in range(q2):
  X,Y=Xs[l%MX],Ys[l//MX];e,cs=sb+ENTRY*l,components(X,Y);d3+=sum(get(mem,e+o,2)!=w for o,w in zip(FRAC_OFF[0]+FRAC_OFF[1],model[l]))
  for i,(dh,dl)in enumerate(frac_errs(mem,e,X,Y)):cb=CT.KB2[i%6//2][i%2];ai[0]=[min(ai[0][0],dh),max(ai[0][1],dh)];ai[1]=[min(ai[1][0],dl),max(ai[1][1],dl)];d5+=not(cb[0][0]<=dh<=cb[0][1]and cb[1][0]<=dl<=cb[1][1])
  ed=six_leaves(X,Y);d6+=sum(dft_leaf(cs[3*s:3*s+3],t,X[0],Y[0])!=ed[2*t+s]for s in(0,1)for t in range(3));leaves.append(ed)
 for(l,s),(x,t,rg)in(override or{}).items():
  f6,c3=represent(x,t,Xs[l%MX][0],Ys[l//MX][0],rg,True)
  for o,v in zip(FRAC_OFF[s],f6):put(mem,sb+ENTRY*l+o,v,2)
  for t_ in range(3):leaves[l][2*t_+s]=(c3[0]+BT[t_]*c3[1]+BT[2*t_%3]*c3[2])%p
 eb=machine.Port1(mem);regs=dict(x2=W2,x3=sb+KEYS0,x26=KN,sp=SP,x30=SNT);out={};nk=a64.run(ae,labels,"keys",regs,mem,out_regs=out,dev=eb.dev,ls64=([(sb,EB_END)],[]));sp_kept&=out["sp"]==SP;bm=machine.Sender1(CT);meter,d7=machine.keys_reference(sys.modules[__name__],mem,sb,q2,W2,KN,bm,CT);ft=[(t,leaves[l][2*t+s])for l in range(q2)for s in(0,1)for t in KO];flat=[x for _,x in ft];ah=[in_W(x,O[t])for t,x in ft];h8=[in_D(x)for x in flat];wrong=sum(not w and(f[0]+O[t])%3!=exk(x)[1]for w,f,(t,x)in zip(ah,d7,ft));e9=sum(f[1:]!=class_label(x)for w,d,f,x in zip(ah,h8,d7,flat)if not w and not d);ec=0
 for k,(d,cn)in enumerate(bm.seq if chain else()):l=(cn[3]>>2)&((1<<24)-1);s=k//3%2;Y,X=chain_y(chain,l//MX),Xs[l%MX];j4=[class_label(PROD.add((pow(BETA,t,p)*X[0]%p,X[1]),PROD.neg(Y)if s else Y)[0])for t in range(3)];ec+=l!=k//6 or Y!=Ys[l//MX]or(not ah[k]and not h8[k]and(d,cn[0]&M42,cn[0],cn[1])not in j4)
 g5={k:[sum(v for i,v in ce[k].items()if fx[i]==e)for e in(5,6,8)]for k in ce};db={k:5*g5[k][0]+6*g5[k][1]+8*g5[k][2]for k in ce};return dict(fwd=nf,bwd=nb,keys=nk,acc_ok=ap%p==h7,fraction_mismatches=d3,fractions_outside_bounds=d5,fraction_errors=ai,identity_mismatches=d6,fixups_3_4=g5,fixup_extra=db,record_decode_mismatches=ec,fixup_words_ok=words==[C0*sum(g5[k])for k in("fwd","bwd")]and all(sum(g5[k])<=w>>32<=sum(g5[k])*C0>>32 for w,k in zip(words,("fwd","bwd"))),fwd_bwd_formulas=nf==CT.fwd_count(q2)+db["fwd"]and nb==CT.bwd_count(q2)+db["bwd"],records_ok=eb.blocks==bm.msgs and machine.state1_ok(mem,bm),records=len(bm.seq),key_paths=meter,in_W=sum(ah),in_D=sum(d and not w for w,d in zip(ah,h8)),wrong_index=wrong,keys_not_canonical=e9,keys_formula=nk==CT.KEY_FIX+CT.KEY_INC*q2+meter+bm.extra,final_R_times_rho_is_one=b9*rho%p==1,sp_kept=sp_kept,reg_trace_ok=all(any(lo<=a<hi for rl,_,lo,hi,_ in CT.REGISTRY if rl=='G')for a in mem.touched))
def macro_counts(bo,entry,ax):
 a9,cf,bn={},None,[ln.split("//")[0].strip()for ln in bo.split("\n")]
 for s in bn:
  if s.startswith(".macro"):cf=s.split()[1].rstrip(",");a9[cf]=[]
  elif s.startswith(".endm"):cf=None
  elif cf is not None:a9[cf].append(s)
 out=dict.fromkeys(ax,0)
 def walk(seq):
  for s in seq:
   if not s:continue
   w=s.split()[0].rstrip(",")
   if w in out:out[w]+=1
   if w in a9:walk(a9[w])
 start=bn.index(entry+":");walk(bn[start+1:bn.index("ret",start)]);return out
def sqrt_mod_p(v):s_=pow(v,(p+1)//4,p);return s_ if s_*s_%p==v%p else None
def degenerate_batch(ae,labels,Xs,Ys):
 rprog,rlab=a64.load(kf("route.s"));q2,nY=24,4;sb=EB_END-ENTRY*q2;ee=Ys*nY;LS=([(AXYX0,XT_END),(BT0,1<<64),(YR0,YR0+512*nY),(sb,EB_END)],[(sb,EB_END)]);mem=_TraceLC1Dict(rprog,rlab);popbat(rprog,rlab,mem,Xs,ee,7)
 for a in range(YR0-64,YR0,8):mem[a]=0
 for yi in range(nY):put_yrec(mem,YR0+512*yi,ee[yi])
 for l in range(q2):x_words(mem,sb+ENTRY*l,Xs,l)
 regs=dict(x0=q2,x1=YR0+64,x3=sb,x24=7,x25=0,x26=0,x27=0,sp=A_SP,x30=SNT);out={};a64.run(ae,labels,"fwd",regs,mem,out_regs=out,ls64=LS);ap=rgi(out,R4(24));regs=dict(x3=EB_END,x1=YR0-32,sp=YR0+512*(nY-1),x30=SNT);regs.update(ars(pow(ap,p-2,p),("x0","x2","x14","x15")));a64.run(ae,labels,"bwd",regs,mem,ls64=LS);good=0
 for l in range(q2):
  X=Xs[l%MX]
  if(Ys[0][1]**2-X[1]**2)%p:good+=sum(abs(dh)<4096 and abs(dl)<4096 for dh,dl in frac_errs(mem,sb+ENTRY*l,X,Ys[0]))
 return ap,good
def edge_batch(ae,labels,Xs,Ys,bd,CT,rg,ay,co):
 over={}
 for l in range(24):
  for s in(0,1):
   t=rg.rr(3)
   while True:
    x=edge_leaf(rg,ay[(2*l+s)%len(ay)])
    if len(over)>=4 or(in_W(x,O[t])if co=="W"else in_D(x)and not in_W(x,O[t])):break
   over[(l,s)]=(x,t,rg)
 return leaf_batch(ae,labels,Xs,Ys,24,bd,CT,override=over)
def tbld_production_tables(ff,J):
 i8=kf("route.s");rprog,rlab=a64.load(i8);eg,ei=a64.load(i8.replace("tst x6,#1\ncsel x18,x3,x2,ne\nlsr x6,x6,#1","cmp x8,#2\ncsel x18,x3,x2,eq"));fp,fl=a64.load(kf("secp_fe.s"));f2,j6=a64.load(ff+sites.SITE_KERNELS+sites.T_KERNELS);cq,ab,sb=0xDEADBEEFCAFEF00D,0x300000,EB_END-ENTRY*MX;j5=([(AXYX0,XT_END),(BT0,BT0_END),(BT1,1<<64),(YR0,YR0+512),(sb,EB_END)],[(sb,EB_END)]);b1=lambda ts:sum(t<<sh for t,sh in zip(ts,W6S))
 def psn(m,b,nw):
  for i in range(nw):m[b+8*i]=cq
 def run_h(ent,regs,m):o={};a64.run(f2,j6,ent,dict(regs,sp=YR0,x30=SNT),m,out_regs=o,ls64=j5);return o
 A1,Bv,cd=J.rr(1,p//2),J.rr(1,p),J.rr(1,1<<63);A2=A1+cd;inv_e,ep=J.rr(1,p),sb+576;mc={};put(mc,YR0+288,Bv);put(mc,ep+224,A2);mc[ep+304]=cd;oc=run_h("t_c1b",dict(x3=ep,**ars(inv_e,R4(25))),mc);w1=inv_e*(Bv-A2)%p;w2=(w1+cd*inv_e)%p;bday_c1b_ok=(get(mc,ep-64)%p,rgi(oc,R4(8))%p)==(w1,w2);ma,cax,cyx,de={},J.rr(1,p),J.rr(1,p),ax_bases(0);cp={};na=0;q=de[4]+544;ma[q]=cp[q]=cq
 for k,sh in enumerate(W6S):
  dk=ab+128*k;put(ma,dk,cax*pow(2,sh,p)%p);put(ma,dk+32,0);put(ma,dk+64,cyx*pow(2,sh,p)%p);put(ma,dk+96,0)
  for a_ in range(dk,dk+128,8):cp[a_]=ma[a_]
  na+=a64.run(rprog,rlab,"tbld",dict(x0=de[k],x1=dk,x2=16,x3=8,x4=2,x5=64,x6=3*(k>0),x30=SNT),ma);a64.run(eg,ei,"tbld",dict(x0=de[k],x1=dk,x2=16,x3=8,x4=2,x5=64,x6=3*(k>0),x30=SNT),cp)
 axyx_ok=(na==5*26+20*(2*16+4*2*8)and ma[q]==cq and cp[q]!=cq);mp,cr,bj,ej,cj,df,cu={},J.rr(1,p),J.rr(1,p),J.rr(1,p),J.rr(1,p),J.rr(1,p),b1([5,1,0,1,1])
 def bp50(m,tb,c1,a1,c2,a2,dm=rprog,dn=rlab):
  au=0
  for k,sh in enumerate(W6S):psn(m,tb+D6O[k],64);dk=ab+128*k;put(m,dk,c1*pow(2,sh,p)%p);put(m,dk+32,a1 if k==0 else 0);put(m,dk+64,c2*pow(2,sh,p)%p);put(m,dk+96,a2 if k==0 else 0);au+=a64.run(dm,dn,"tbld",dict(x0=tb+D6O[k],x1=dk,x2=8,x3=2,x4=2,x5=64,x6=3*(k>0),x30=SNT),m)
  return au
 j7,az=bp50(mp,XTP[0],cr,bj,ej,cj),offs6(cu);put(mp,YR0+64,df);put(mp,YR0+96,0);mp[sb+416]=XTP[0]
 for k,v in enumerate(az[:4]):mp[YR0+96+8*k]=v
 mp[YR0+352]=az[4];ol=run_h("t_lkp",dict(x1=YR0+64,x3=sb),mp);ef=dict(mp);bp50(ef,XTP[0],cr,bj,ej,cj,eg,ei);olm=run_h("t_lkp",dict(x1=YR0+64,x3=sb),ef);d8=(j7==5*26+20*(2*8+4*2*2)and rgi(ol,R4(16))%p==(cr*cu+bj+df)%p and rgi(olm,R4(16))%p!=(cr*cu+bj+df)%p); PB,SK,kok=0x600000,0x500000,True
 def tk(m,cl,dst,starts,steps,nrow):
  for i,v in enumerate(eis_words(CE.EIS[cl])):m[PB+8*i]=v
  for a in range(SK-4096,SK,8):m.setdefault(a,0)
  for j in range(3):put(m,ab+64*j,steps[j]);put(m,ab+64*j+32,starts[j])
  return a64.run(fp,fl,"tbldk",dict(x0=dst,x1=ab,x2=nrow,x3=PB,sp=SK,x30=SNT),m)
 for cl in("E","FG"):
  me,us,ss={},[J.rr(p)for _ in range(3)],[J.rr(p)for _ in range(3)];n=tk(me,cl,sb,ss,us,3);rk=[ered(cl,(ss[j]+t*us[j])%p)for t in range(3)for j in range(3)]
  kok&=n==CE.tbldk_count(3,cl)and[(get(me,sb+128*t+40*j,3),get(me,sb+128*t+40*j+24,2))for t in range(3)for j in range(3)]==[(c1%2**192,c2)for c1,c2 in rk]and all(me[sb+128*t+120]==0 for t in range(3))
 Xa=[J.rr(1,p)for _ in range(6)];wv=[J.gb(256)for _ in range(3)];bs=ered("EB",J.rr(p));mk={};ex={}
 for i,(cl,ce,cc)in enumerate((("E",lambda a:2*a**3,lambda a:0),("FG",lambda a:-2*a,lambda a:-a**4),("FG",lambda a:BT[2]*a*a,lambda a:0))):
  for g in(0,1):
   for k,(sh,m)in W6:t=(wv[i]>>sh)&m;tk(mk,cl,XTK[2*i+g]+D7O[k]+(t<<7),[(ce(Xa[3*g+j])*t*pow(2,sh,p)+(cc(Xa[3*g+j])if k==0 else 0))%p for j in range(3)],[0,0,0],1)
  for j in range(6):ex[i,j]=ksum(kents(cl,ce(Xa[j]),wv[i],cc(Xa[j]))+([bs]if i==0 else[]))
 for v,ads in zip(wv,((384,392,400,408,416),(480,488,496,504,424),(440,448,456,464,472))):
  for a,o in zip(ads,offs7(v)):mk[YR0+a]=o
 put(mk,YR0,bs[0]%2**192,3);put(mk,YR0+24,bs[1],2);m2=dict(mk);ea=XTK[0]+D7O[2]+(((wv[0]>>W6S[2])&W6M[2])<<7)+24;m2[ea]=(m2[ea]+1)%2**64
 def lines(m):return{(i,j):(get(m,sb+ENTRY*j+(384,448,192)[i]+8*(j%3==1),3),get(m,sb+ENTRY*j+(384,448,192)[i]+8*(j%3==1)+24,2))for i in range(3)for j in range(6)}
 oe=run_h("t_efgk",dict(x3=sb+5*ENTRY+64,x0=1,x2=2,x14=3,x15=4),mk);run_h("t_efgk",dict(x3=sb+5*ENTRY+64,x0=1,x2=2,x14=3,x15=4),m2);want={k:(v[0]%2**192,v[1])for k,v in ex.items()};got,bad=lines(mk),lines(m2)
 efgk_ok=got==want and bad!=want and bad[0,0]!=want[0,0]and all(bad[k]==want[k]for k in want if k!=(0,0))and[oe[r]for r in("x0","x2","x14","x15","x3","x4")]==[1,2,3,4,sb+5*ENTRY+64,C0];ok=bday_c1b_ok and axyx_ok and d8 and kok and efgk_ok;return dict(ok=ok,bday_c1b_ok=bday_c1b_ok,axyx_ok=axyx_ok,lookup_pair_ok=d8,tbldk_ok=kok,efgk_ok=efgk_ok)
def part_a_kernels(J,CT):
 aa={};bo=kf("leaf.s");ae,labels=a64.load(bo);bh,ca=a64.load(kf("secp_fe.s"));rprog,rlab=a64.load(kf("route.s"));yg,bd,dj,em=y_generator(bh,ca,J,CT);aa["y_generator"]=yg;Xs=make_bday_xs(J);bp,en=defaultdict(int),{};cg,eo=34*(1<<40)+0x20000,34*(1<<40)+0x40000
 for q in range(4):put(bp,cg+128*q,Xs[2*q][1]*Xs[2*q][1]%p);put(bp,cg+128*q+32,Xs[2*q+1][1]*Xs[2*q+1][1]%p);put(bp,cg+128*q+64,2*q+1,2);put(bp,cg+128*q+80,2*q+2,2)
 a64.run(rprog,rlab,"psel",dict(x0=cg,x1=4,x2=eo,sp=34*(1<<40)+0x18000,x30=SNT),bp,out_regs=en);aa["bday_xs_psel_ok"]=(en["x0"]==1 and((Xs[0][1]*Xs[0][1]%p)>>63)!=((Xs[1][1]*Xs[1][1]%p)>>63)and[get(bp,eo+128*q+32*s)for q in range(4)for s in(0,1)]==[X[1]*X[1]%p for X in Xs]);Ys=[chain_y(dj,i)for i in range(12)];runs={}
 for q2 in(24,48):runs[q2]=leaf_batch(ae,labels,Xs,Ys,q2,bd,CT,rprog=rprog,rlab=rlab,chain=dj)
 f3={(0,0):(planted_x(J,(1<<41)-4,0)[0],0,J),(0,1):(planted_x(J,(1<<41)-3,1)[0],1,J)};a5=leaf_batch(ae,labels,Xs,[chain_y(em,i)for i in range(4)],24,bd,CT,rprog=rprog,rlab=rlab,override=f3,chain=em,ybase=YR0+YR_STRIDE);aa["yslot_lc1_ok"]=(a5["records_ok"]and a5["keys_formula"]and a5["fwd_bwd_formulas"]and a5["sp_kept"]and a5["reg_trace_ok"]);aa["registry_trace_ok"]=all(d["reg_trace_ok"]for d in runs.values())and a5["reg_trace_ok"];aa["runs"]={q:dict(fwd=d["fwd"],bwd=d["bwd"],keys=d["keys"],fixup_extra=d["fixup_extra"],records_ok=d["records_ok"])for q,d in runs.items()};ek={q:{k:d[k]-d["fixup_extra"][k]for k in("fwd","bwd")}for q,d in runs.items()};aa["fb_hold"]=all(d["fwd_bwd_formulas"]and d["fixup_words_ok"]for d in runs.values());aa["counts_ok"]=all(ek[q]["fwd"]==CT.fwd_count(q)and ek[q]["bwd"]==CT.bwd_count(q)for q in runs);aa["fixups"]=dict(sites=len(fixup_extra(ae)),flags_dead_at_every_join=flags_dead_after_fixups(ae,labels));j0=[o for mn,o in ae if mn in("tbz","tbnz")];aa["loop_bit34"]=[o for o in j0 if o[0]=="x3"]==[["x3","#34","1b"]]and all(o[1]=="#14"for o in j0 if o[0]!="x3");aa["keys_hold"]=all(d["keys_formula"]for d in runs.values());aa["fraction_errors"]=[[min(d["fraction_errors"][g][0]for d in runs.values()),max(d["fraction_errors"][g][1]for d in runs.values())]for g in(0,1)];aa["all_correct"]=(yg["mismatches"]==0 and yg["centres_ok"]and yg["counts_ok"]and aa["bday_xs_psel_ok"]and aa["yslot_lc1_ok"]and all(d["acc_ok"]and d["fraction_mismatches"]==0 and d["fractions_outside_bounds"]==0 and d["identity_mismatches"]==0 and d["record_decode_mismatches"]==0 and d["records_ok"]and d["keys_not_canonical"]==0 and d["wrong_index"]==0 and d["final_R_times_rho_is_one"]and d["sp_kept"]and d["reg_trace_ok"]for d in runs.values()));j8=random.Random(SEED+49979687)
 for name,ay,co in(("window_batch",("Xm","diff"),"W"),("boundary_batch",("D",),"D")):
  w=edge_batch(ae,labels,Xs,Ys,bd,CT,j8,ay,co);aa[name]=dict(in_W=w["in_W"],in_D=w["in_D"],ok=(w[{"W":"in_W","D":"in_D"}[co]]>=4 and w["wrong_index"]==0 and w["keys_not_canonical"]==0 and w["records"]==6*24 and w["records_ok"]and w["keys_formula"]and w["fwd_bwd_formulas"]and w["fixup_words_ok"]))
 gp,gl=a64.load(bo.replace("b.lo 9f\n","cbnz xzr,9f\n").replace("b.hs 9f\n","cbnz xzr,9f\n"));g=leaf_batch(gp,gl,Xs,Ys,24,bd,CT,rprog=rprog,rlab=rlab,chain=dj);f4=sum(CT.FIX_SITES_QUAD.values());f5=sum(CT.FIX_EXTRA[m]*CT.FIX_SITES_QUAD[m]for m in CT.FIX_SITES_QUAD);aa["forced_fixups"]=dict(fixup_extra=g["fixup_extra"],ok=(g["fraction_mismatches"]==0 and g["fractions_outside_bounds"]==0 and g["identity_mismatches"]==0 and g["fwd_bwd_formulas"]and g["fixup_words_ok"]and g["records_ok"]and sum(map(sum,g["fixups_3_4"].values()))==6*f4 and sum(g["fixup_extra"].values())==6*f5));h9,ok=set(),True
 for _ in range(3):a=J.rr(1,p);mem={};put(mem,0x1000,a);h9.add(a64.run(bh,ca,"fe_inv",dict(x0=0x2000,x1=0x1000,sp=0x100000,x30=SNT),mem,limit=10**6));ok&=get(mem,0x2000)%p==pow(a,p-2,p)
 aa["fe_inv"]=dict(counts=sorted(h9),correct=ok and h9=={CT.INV});ax=("SHX","FDX","F5","FS","FADDW","LKFOLD4","DFX","ROWSG");counts={k:macro_counts(bo,k,ax+("FWD6_ONE","BWDPK"))for k in("fwd","bwd","keys")};e1,e2=counts["fwd"]["FWD6_ONE"],counts["bwd"]["BWDPK"];a6={}
 for n in ax:a6[n]=(counts["fwd"][n]//2+counts["bwd"][n])//3
 a6["BLK"],a6["DF"],a6["RW"]=a6.pop("LKFOLD4"),a6.pop("DFX"),a6.pop("ROWSG");aa["site_counts"]=dict(per_x_quad=a6,pair_copies=dict(fwd=e1,bwd=e2));aa["sites_match"]=(a6==CT.FIX_SITES_QUAD and e1==4 and e2==6 and sum(counts["keys"].values())==0);j9=random.Random(SEED+15485863);aa["site_checks"]=sites.site_checks(sys.modules[__name__],j9,bh,ca);aa["tbld_production"]=tbld_production_tables(bo,random.Random(SEED+27182818));phY=(BETA*Ys[0][0]%p,Ys[0][1]);k0,good=degenerate_batch(ae,labels,[Xs[0],phY,Xs[1],PROD.neg(Ys[0])]+Xs[2:4],[Ys[0]]);aa["d_zero_batch"]=dict(acc_zero=k0%p==0,fractions_within_bounds=good,all_other_fractions_wrong=good==0);e3,_=a64.load(kf("route.s"));e0=len(ae)+len(e3)+len(bh);aa["prog_insns"]=dict(leaf=len(ae),route=len(e3),secp_fe=len(bh),total=e0,fits_2_18_bytes=4*e0<2**18,stated=f"{len(ae):,} / {len(e3):,} / {len(bh):,} ({e0:,} insns = {4*e0:,}"in CT._PF);return aa
def lattice_point(tv,tu):c1=(2*(FA*tv+(FA-FB)*tu)+p)//(2*p);c2=(2*((FA-FB)*tv-FB*tu)+p)//(2*p);return c1*FB+c2*(FA-FB),c1*(FA-FB)-c2*FA
def planted_x(rg,d1,d2,v=None,j=None):
 while True:
  vc=v
  if vc is None:
   ha=(rg.gb(23)<<41)|(rg.rr(16)<<20)|rg.rr(200,(1<<20)-200);la=(d1<<23)|rg.rr(200,(1<<23)-200);hb=(rg.gb(14)<<50)|(rg.rr(2)<<20)|rg.rr(200,(1<<20)-200);lb=(d2<<23)|rg.rr(200,(1<<23)-200);vv,uu=lattice_point(((ha<<64|la)*p)>>128,((hb<<64|lb)*p)>>128)
   if not(0<vv<p and uu==BETA*vv%p and exk(vv)[0]==vv):continue
   rf=label(vv)
   if not(rf[0]==d1 and rf[1]==d2 and(rf[3]&((1<<50)-1))<16):continue
   vc=vv
  m=rg.rr(3)if j is None else j;x=vc*pow(BETA,(3-m)%3,p)%p
  if exk(x)[1]==m and not in_W(x)and not in_D(x):return x,vc
def load_certificate():
 er=io.StringIO()
 with contextlib.redirect_stdout(er):import certificate as CT
 bn=[" ".join(ln.split())for ln in er.getvalue().splitlines()];return CT,[ln for ln in bn if" ok"in ln or"FAIL"in ln or ln.startswith(("q1 ","E_total","declared:"))]
def part_d(CT):
 return dict(passed=not CT.FAILED,failed=CT.FAILED,checks_passed=len(CERT_LINES),checks=CERT_LINES,claim="%d.%015d"%divmod(int(CT.CLAIM_E*10**15),10**15),q1=CT.Q1GD,per_leaf_declared=format(float(CT.DECLARED["per_leaf"]),".4f"),mu=float(CT.TFULL["mu"]),opened_success_lower_bound=float(CT.LO_RUN),gate_prob=float(CT.P_GATE),success_lower_bound=float(CT.GATED_LO),per_leaf=float(CT.PER_LEAF_OPEN),groups=CT.G_GROUPS)
def main():
 global CERT_LINES,CE,ERR_H,ERR_L
 J=random.Random(SEED);CT,CERT_LINES=load_certificate();CE=CT;ERR_H,ERR_L=CT.ERRH,CT.ERRL;OUT["input"]=cfg;g4=OUT["A_identities"]=part_a_identities(J);g1=(g4["beta_cube"]and g4["lambda_cube"]and g4["lambda_G_equals_beta_x"]and g4["no_x0_point"]and g4["win_nonzero"]and g4["label_constants"]and g4["class_rule"]["mismatches_outside_window_set"]==0 and g4["class_labels"]["ok"]and g4["c132_ok"])
 af=OUT["A_kernels"]=part_a_kernels(J,CT);ao=OUT["B_machine"]=machine.part_b(sys.modules[__name__],CT);OUT["counts_match_certificate"]=cm=(af["fb_hold"]and af["counts_ok"]and af["fixups"]["flags_dead_at_every_join"]and af["loop_bit34"]and af["keys_hold"]and af["all_correct"]and all(af[k]["ok"]for k in("window_batch","boundary_batch","forced_fixups","tbld_production"))and af["fe_inv"]["correct"]and af["sites_match"]and af["site_checks"]["all_true"]and af["d_zero_batch"]["acc_zero"]and af["d_zero_batch"]["all_other_fractions_wrong"]and af["prog_insns"]["fits_2_18_bytes"]and af["prog_insns"]["stated"]and ao["formulas_ok"]and all(ao[k]["ok"]for k in("phase0","block_paths","grouped_paths","ls64","flat64","encoding","phase4_round","y_ring","gate","wait_dispatch","rendezvous","keys_edges","registry_trace")))
 D_c=OUT["D_certificate"]=part_d(CT);OUT["all_checks_pass"]=bool(g1 and cm and af["registry_trace_ok"]and D_c["passed"]);OUT["seconds"]=round(time.time()-T0,2);json.dump(OUT,sys.stdout,separators=(",",": "));print()
 if not OUT["all_checks_pass"]:sys.exit(1)
if __name__=="__main__":main()
