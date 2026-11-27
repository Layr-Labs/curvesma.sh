#!/usr/bin/env python3
"""Standalone Docker-Safe Verifier for stnfs-q15-hybrid-split-eta1015-v8 (stdlib only, 0 external files, 0 offline MC constants)."""
import bisect,cmath,json,random,sys,time
from fractions import Fraction
from math import ceil,cos,exp,floor,gamma,gcd,lgamma,log,log2,pi,sin,sqrt
RN=range
X_BLS=-0xD201000000010000
ABS_X,LOG2_ABS_X=abs(X_BLS),log2(abs(X_BLS))
P_BLS=0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_BLS=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P0_X=0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb
P0_Y=0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1
Q0_X=(0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e)
Q0_Y=(0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be)
C_COST,V_A_SH,SLAB_BASE=20000.0,0.9351*(2.0**56),1<<56
REF={
"pkg":"stnfs-q15-hybrid-split-eta1015-v8","tag":"CAND8_HYBRID_RQMC_DISC_TOWER",
"Bf":71.45,"Bg":71.81,"eta":1.015,"cush":-0.0125,"disc_pad":0.0015,
"qa":0.65,"rbr":0.82,"Ra":1719.4681182270274,"R0":1355.3079808905318,"sl":1.75,
"lam":340599.43367449316,"lf":408.0,"lg":688.0,"oe":10.0,"w2":16,"w3":32,
"b":{"pglo":28.0,"pw1":55.0,"pw2":58.5,"pghi0":62.5,"pghi1":60.7,"pflo":32.0,"pfhi0":62.1,"pfhi1":60.3,"qg":3.4,"qf":3.8},
"s1b":{"p1":(127<<54)+1,"g1":3,"p2":(114<<54)+1,"g2":7,"p3":(108<<54)+1,"g3":5,"bfly":11.125,"pw":12.0,"sc":12.0,"c1":58.0,"c2":90.0,"saf":1.10,"dm":1.67},
"rep":{"ent":21.50,"bkt":5.5},
"la":{"m":512,"n":32,"p":4,"plg":512,"nl":0.00084,"comp":0.047,"pm1_f":0.92,"x4":11.90,"x1":14.90,"rv":960.0,"ov":1.06},
"cr":{"relation_collection":107.9510,"replay_a_tables":104.4575,"linear_algebra":103.3200,"both_individual_logs_exp":96.6500,"ecm_split_expected":96.0500,"meter_safety_buffer":93.4001,"band_split_and_trees":90.0001,"canonical_a_domain":84.0001,"schirokauer_maps":80.0001,"final_log_backsub_and_gcd":72.0001,"matrix_generator_sge":63.1490,"stage0_setup_and_iso":45.0001,"final_verification_Q":40.0001},
"Kg":7194374076228249696,"num":5339,"den":50,"ram_ub":118.0,
}
class ExactDickman:
 def __init__(s,um=30.0,N=400):
  s.um,s.N=um,N;h,t=1.0/N,int(round(um*N))+1;v=[1.0]*(N+1)+[0.0]*(t-N-1)
  for i in RN(N+1,min(2*N+1,t)):v[i]=1.0-log(i*h)
  if t>2*N+1:
   w=sum(v[N+2:2*N+1])
   for i in RN(2*N+1,t):
    if i%N==0:w=sum(v[i-N+1:i])
    v[i]=(h*(0.5*v[i-N]+w))/(i*h-0.5*h);w+=v[i]-v[i-N+1]
  s.v=v;s.lv=[log(x)if x>1e-300 else-700.0 for x in v]
 def rho(s,u):
  if u<=1.0:return 1.0
  if u>=s.um:return 1e-300
  p=u*s.N;i=int(p)
  if i+1>=len(s.v):return s.v[-1]
  f=p-i;return exp((1.0-f)*s.lv[i]+f*s.lv[i+1])
DK=ExactDickman(30.0,400)
def _primes(n):
 ip=bytearray(b"\x01")*(n+1);ip[0]=ip[1]=0
 for i in RN(2,int(n**0.5)+1):
  if ip[i]:ip[i*i:n+1:i]=bytearray(len(RN(i*i,n+1,i)))
 return[i for i in RN(2,n+1)if ip[i]]
def pi_dus(B):L=B*log(2.0);return((2.0**B)/L)*(1.0+1.0/L+2.334/(L*L))
def U_ub(Bf,Bg):return(pi_dus(Bf)+pi_dus(Bg))/12.0+36.0*pi_dus(0.5*Bf)+12.0*pi_dus(0.5*Bg)
def ppb_dus(lo,hi):
 th=(2.0**hi)*(1.0+0.01/max(1.0,hi*log(2.0)));tl=0.0 if lo<=1.0 else(2.0**lo)*(1.0-0.01/max(1.0,lo*log(2.0)))
 return max(1.0,(th-tl)/log(2.0))
def _g1a(p1,p2):
 if p1 is None:return p2
 if p2 is None:return p1
 x1,y1,x2,y2=p1[0],p1[1],p2[0],p2[1]
 if x1==x2:
  if(y1+y2)%P_BLS==0:return None
  l=(3*x1*x1*pow(2*y1,-1,P_BLS))%P_BLS
 else:l=((y2-y1)*pow(x2-x1,-1,P_BLS))%P_BLS
 x3=(l*l-x1-x2)%P_BLS;return(x3,(l*(x1-x3)-y1)%P_BLS)
def check_stage0():
 assert(X_BLS-1)**2*(X_BLS**4-X_BLS**2+1)//3+X_BLS==P_BLS and X_BLS**4-X_BLS**2+1==R_BLS
 l=X_BLS%R_BLS;assert P_BLS%R_BLS==l and(pow(l,4,R_BLS)-pow(l,2,R_BLS)+1)%R_BLS==0 and pow(l,6,R_BLS)==R_BLS-1 and pow(l,12,R_BLS)==1 and pow(l,4,R_BLS)!=1
 rp,ap,k=None,(P0_X,P0_Y),R_BLS
 while k:
  if k&1:rp=_g1a(rp,ap)
  ap=_g1a(ap,ap);k>>=1
 m2=lambda a,b:((a[0]*b[0]-a[1]*b[1])%P_BLS,(a[0]*b[1]+a[1]*b[0])%P_BLS);q3=m2(m2(Q0_X,Q0_X),Q0_X)
 assert(P0_Y*P0_Y-(P0_X**3+4))%P_BLS==0 and rp is None and m2(Q0_Y,Q0_Y)==((q3[0]+4)%P_BLS,(q3[1]+4)%P_BLS)and P_BLS%15==7 and[j for j in RN(1,5)if pow(P_BLS,j,5)==1][0]==4 and pow(X_BLS%P_BLS,(P_BLS-1)//3,P_BLS)!=1
 assert any((X_BLS**3-X_BLS**2+(c*c%P_BLS)*X_BLS-c)%P_BLS==0 for g in RN(2,20)for cd in(pow(g,(P_BLS-1)//3,P_BLS),)if cd!=1 for c in(cd,(cd*cd)%P_BLS))and X_BLS**6-2*X_BLS**5+2*X_BLS**3+X_BLS+1==3*P_BLS and(R_BLS-1)%12==0 and ABS_X==(1<<16)*906349*254760293
 return{"ok":True,"p_mod_15":7,"ord_5_p":4,"n_sm_chi1":4,"abs_x":"2^16*906349*254760293"}
def _tr(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def _pm(a,b,p):
 a,b=_tr(list(a)),_tr(list(b));iv,db=pow(b[-1],p-2,p),len(b)-1
 while len(a)-1>=db and not(len(a)==1 and a[0]==0):
  c,sh=(a[-1]*iv)%p,len(a)-1-db
  for k in RN(db+1):a[sh+k]=(a[sh+k]-c*b[k])%p
  a.pop()
  if not a:return[0]
  _tr(a)
 return _tr(a)
def _pg(u,v,p):
 u,v=_tr(list(u)),_tr(list(v))
 while not(len(v)==1 and v[0]==0):u,v=v,_pm(u,v,p)
 iv=pow(u[-1],p-2,p);return[(x*iv)%p for x in u]
def _sr(F,p):
 d=len(F)-1
 if d==0:return[]
 if d==1:return[((-F[0])*pow(F[1],p-2,p))%p]
 def ml(x,y):
  r=[0]*(len(x)+len(y)-1)
  for i,vx in enumerate(x):
   if vx:
    for j,vy in enumerate(y):r[i+j]=(r[i+j]+vx*vy)%p
  return _pm(r,F,p)
 a=1
 while True:
  rs,bs,e=[1],[a%p,1],(p-1)//2
  while e:
   if e&1:rs=ml(rs,bs)
   bs=ml(bs,bs);e>>=1
  rs[0]=(rs[0]-1)%p;g=_pg(F,rs,p)
  if 1<=len(g)-1<d:
   F2,d1=list(F),len(g)-1;q=[0]*(d-d1+1)
   for k in RN(d-d1,-1,-1):
    c=F2[k+d1];q[k]=c
    for j in RN(d1+1):F2[k+j]=(F2[k+j]-c*g[j])%p
   return _sr(g,p)+_sr(_tr(q),p)
  a+=1
def _fcr(w,p,e3):
 w2,r0,r1,r2=(w*w)%p,0,1,0
 for bt in bin(p)[3:]:
  c0,c1,c2,c3,c4=r0*r0,2*r0*r1,2*r0*r2+r1*r1,2*r1*r2,r2*r2
  c3+=c4;c2-=c4*w2;c1+=c4*w;r0=(c0+c3*w)%p;r1=(c1-c3*w2)%p;r2=(c2+c3)%p
  if bt=="1":r0,r1,r2=(r2*w)%p,(r0-r2*w2)%p,(r1+r2)%p
 r1=(r1-1)%p;Q=[(-w)%p,w2,p-1,1];ny=0
 if r0==0 and r1==0 and r2==0:
  for rt in _sr(Q,p):
   if pow(rt,e3,p)==1:ny+=3
 else:
  g=_pg(Q,[r0,r1,r2],p)
  if len(g)==2 and pow(((-g[0])*pow(g[1],p-2,p))%p,e3,p)==1:ny+=3
 return ny
def check_chebotarev(p_lim=18000):
 zm=lambda x,y:(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1]);zn=lambda x:x[0]*x[0]-x[0]*x[1]+x[1]*x[1];zd=lambda d,n:all(c%zn(d)==0 for c in zm(n,(d[0]-d[1],-d[1])))
 dw,dw2,p19,p67=(41,24),(17,-24),(-5,-2),(-9,-2);pw,ww=991,113;wr=[t for t in RN(pw)if(t**3-t**2+(ww*ww%pw)*t-ww)%pw==0];kc=pow((wr[0]*pow(wr[1],pw-2,pw))%pw,(pw-1)//3,pw)
 assert zm(p19,p67)==dw and zn(dw)==1273==19*67 and gcd(1273,15)==1 and not zd(p19,dw2)and not zd(p67,dw2)and(118*118+118+1)%151==0 and all((t**3-t**2+32*t-118)%151!=0 for t in RN(151))and pw%45==1 and(ww*ww+ww+1)%pw==0 and wr==[373,802,808]and kc==877 and pow(kc,3,pw)==1
 s3=((0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0));G=[((a,b,c),s)for a in RN(3)for b in RN(3)for c in RN(3)for s in s3];gid={g:i for i,g in enumerate(G)}
 gm=lambda a,b:(((a[0][b[1][0]]+b[0][0])%3,(a[0][b[1][1]]+b[0][1])%3,(a[0][b[1][2]]+b[0][2])%3),(a[1][b[1][0]],a[1][b[1][1]],a[1][b[1][2]]))
 mt=[[gid[gm(a,b)]for b in G]for a in G];it=[next(j for j in RN(162)if mt[i][j]==0)for i in RN(162)];sn,cl=set(),[]
 for i in RN(162):
  if i not in sn:cc=frozenset(mt[mt[h][i]][it[h]]for h in RN(162));sn|=cc;cl.append(cc)
 hz=frozenset(i for i,(v,_)in enumerate(G)if sum(v)%3==0);ns,q={frozenset([0])},[frozenset([0])]
 while q:
  N=q.pop()
  for cc in cl:
   if cc.issubset(N):continue
   M=set(N)|set(cc);fr=set(cc)
   while fr:nx={mt[a][b]for a in M for b in fr}-M;M|=nx;fr=nx
   Mf=frozenset(M)
   if Mf not in ns:ns.add(Mf);q.append(Mf)
 sh=[N for N in ns if N<hz];mx=[N for N in sh if not any(N<M for M in sh)]
 assert len(G)==162 and len(cl)==22 and sorted(len(N)for N in ns)==[1,3,9,27,27,54,81,162]and len(mx)==1 and len(mx[0])==27
 rc,c0,c1,r0,r1={0:0,3:0,6:0,9:0},0,0,[],[]
 for v,s in G:
  fx=sum(1 for i in RN(3)for c in RN(3)if s[i]==i and(c+v[i])%3==c);rc[fx]+=1
  if sum(v)%3==0:r0.append(fx);c0+=int(fx>0)
  else:r1.append(fx);c1+=int(fx>0)
 p0,p3,p6,p9=(Fraction(rc[k],162)for k in(0,3,6,9));dok,d0,dnz=1-p0,Fraction(c0,54),Fraction(c1,108);df=Fraction(1,8)*(Fraction(1,3)*(1-(1-d0)**2)+Fraction(2,3)*(1-(1-dnz)**2))
 m0,m1=sum(r0)/len(r0),sum(r1)/len(r1);v0,v1=sum((x-m0)**2 for x in r0)/len(r0),sum((x-m1)**2 for x in r1)/len(r1)
 assert(p0,p3,p6,p9)==(Fraction(58,81),Fraction(13,54),Fraction(1,27),Fraction(1,162))and dok==Fraction(23,81)and df==Fraction(2131,34992)and abs(((1.0/3.0)*16.0*(2.0*v0+(2.0*m0-2.0)**2)+(2.0/3.0)*16.0*(2.0*v1+(2.0*m1-2.0)**2))/8.0-12.0)<1e-9
 ps=_primes(p_lim);p15=[p for p in ps if p%15==1];assert len(p15)>200;bm=fp=ca=0
 for p in p15:
  e3=(p-1)//3;w=next(pow(g,e3,p)for g in RN(2,p)if pow(g,e3,p)!=1);ph=False
  for wv in(w,(w*w)%p):
   ny=_fcr(wv,p,e3)
   if ny%3:bm+=1
   if ny:ca+=1;ph=True
  if ph:fp+=1
 od,odn=ca/(2.0*len(p15)),fp/float(len(ps));assert bm==0 and abs(od-float(dok))<0.04 and abs(odn-float(df))<0.01
 return{"ok":True,"delta_Q":"41+24*w","kummer_p":pw,"galois_order":162,"d_ok_f":str(dok),"dens_f":str(df),"var_f_root":12.0,"var_g_root":16.0,"split_p":len(p15),"live_d_ok":round(od,5),"live_dens":round(odn,5)}
def _prc(pl,qe,z,o,ad,sb,ml,iv):
 def tf(a):
  while len(a)>1 and a[-1]==z:a.pop()
  return a
 def pm(a,b):
  r=[z]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   if x!=z:
    for j,y in enumerate(b):
     if y!=z:r[i+j]=ad(r[i+j],ml(x,y))
  return tf(r)
 def md(a,m):
  a,dm,il=list(a),len(m)-1,iv(m[-1])
  while len(a)-1>=dm:
   c,sh=ml(a[-1],il),len(a)-1-dm
   if c!=z:
    for i in RN(dm+1):a[sh+i]=sb(a[sh+i],ml(c,m[i]))
   a.pop()
  return tf(a)
 rs,b,e=[o],md([z,o],pl),qe
 while e:
  if e&1:rs=md(pm(rs,b),pl)
  b=md(pm(b,b),pl);e>>=1
 if len(rs)<2:rs+=[z]*(2-len(rs))
 rs[1]=sb(rs[1],o);tf(rs)
 if rs==[z]:return len(pl)-1
 a,c=list(pl),rs
 while any(x!=z for x in c):
  a,c=c,md(a,c)
  if c==[z]:break
 return len(a)-1
def _fast_murphy(Pmax=500000):
 ps=_primes(Pmax);af=ag=af5=ag5=0.0;ck={};l2=log2
 for p in ps:
  er,lp=1.0/(p-1.0),l2(p)
  if p<=300:
   if p==2:ef,eg=0.0,8.0*(1.0/17.0)*(sum(3.0*k*(15.0/(16.0**k))for k in RN(1,6))+16.0/(16.0**5))
   elif p==3:ef,eg=4.0*(3.0/82.0),4.0*(1.75/82.0)
   elif p==5:ef,eg=0.0,2.0*3.0*25.0/(25.0**2-1.0)
   else:
    cu,fd=p%15,1
    while cu!=1:cu=(cu*p)%15;fd+=1
    gc,qs=8//fd,p**fd;pf,xm=float(qs),X_BLS%p
    rg=1 if xm==0 else((3 if pow(xm,(qs-1)//3,p)==1 else 0)if(qs-1)%3==0 else 1)
    eg=fd*gc*rg*pf/(pf*pf-1.0);op=(0,1,lambda a,b:(a+b)%p,lambda a,b:(a-b)%p,lambda a,b:(a*b)%p,lambda a:pow(a,p-2,p))
    if fd==1:
     w1=next(pow(g,(p-1)//3,p)for g in RN(2,p)if pow(g,(p-1)//3,p)!=1)
     ef=sum(4*_prc([(-w)%p,0,0,(w*w)%p,0,0,(-1)%p,0,0,1],p,*op)for w in(w1,(w1*w1)%p))*pf/(pf*pf-1.0)
    elif fd==2 or(fd==4 and p<=13):
     if p%3==1:
      w1=next(pow(g,(p-1)//3,p)for g in RN(2,p)if pow(g,(p-1)//3,p)!=1)
      rf=sum((4//fd)*_prc([(-w)%p,0,0,(w*w)%p,0,0,(-1)%p,0,0,1],qs,*op)for w in(w1,(w1*w1)%p))
     else:
      fi=lambda a:(((a[0]-a[1])*pow((a[0]*a[0]-a[0]*a[1]+a[1]*a[1])%p,p-2,p))%p,((-a[1])*pow((a[0]*a[0]-a[0]*a[1]+a[1]*a[1])%p,p-2,p))%p)
      pf2=[(0,(-1)%p),(0,0),(0,0),((-1)%p,(-1)%p),(0,0),(0,0),((-1)%p,0),(0,0),(0,0),(1,0)]
      rf=gc*_prc(pf2,qs,(0,0),(1,0),lambda a,b:((a[0]+b[0])%p,(a[1]+b[1])%p),lambda a,b:((a[0]-b[0])%p,(a[1]-b[1])%p),lambda a,b:((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0]-a[1]*b[1])%p),fi)
     ef=fd*rf*pf/(pf*pf-1.0)
    else:ef=8.0/pf
  else:
   r15=p%15
   if r15==1:
    pf,xm,e3=float(p),X_BLS%p,(p-1)//3;rg=1 if xm==0 else(3 if pow(xm,e3,p)==1 else 0);eg=8.0*rg*pf/(pf*pf-1.0);w1,g=1,2
    while True:
     c=pow(g,e3,p)
     if c!=1:w1=c;break
     g+=1
    ef=4.0*(_fcr(w1,p,e3)+_fcr((w1*w1)%p,p,e3))*pf/(pf*pf-1.0)
   elif r15 in(11,14):qs=float(p*p);eg,ef=24.0/qs,8.0/qs
   elif r15==4:qs=float(p*p);eg,ef=(24.0 if pow(X_BLS%p,(p-1)//3,p)==1 else 0.0)/qs,8.0/qs
   elif r15 in(2,8):qs=float(p*p);qs*=qs;eg,ef=24.0/qs,8.0/qs
   else:qs=float(p*p);qs*=qs;eg,ef=(24.0 if pow(X_BLS%p,(p-1)//3,p)==1 else 0.0)/qs,8.0/qs
  af+=(er-ef)*lp;ag+=(er-eg)*lp
  if p<=500000:af5,ag5=af,ag
  if p in(4999,49999,199999,499979,999983,1499977):ck[str(p)]=[round(af,5),round(ag,5)]
 ck["500000"]=[round(af5,5),round(ag5,5)]
 return af,ag,ck
def _cof_ex(Bf,ph,mnf,lf,af):
 M,dx,h=80,(Bf-ph)/80.0,2.0;nb=int(ceil(660.0/h))+1;cu=[0.0]*nb;cu[0]=1.0;tp=ep=0.0;ne=mnf+af
 for k in RN(1,11):
  nx=[0.0]*nb
  for m in RN(M):
   x=ph+(m+0.5)*dx;w,sh=dx/(x*k),int(round(x/h))
   for b in RN(nb-sh):
    if cu[b]>0.0:nx[b+sh]+=cu[b]*w
  cu=nx
  for b in RN(nb):
   sb=b*h
   if sb>ne:break
   cv=cu[b]*DK.rho((ne-sb)/ph);tp+=cv
   if sb>lf:ep+=cv
 tp+=DK.rho(ne/ph);return ep/tp,tp/DK.rho(ne/Bf)
def check_murphy_norms(pmax=500000):
 af,ag,ck=_fast_murphy(pmax);b,Ra,R0,q,Bf,Bg=REF["b"],REF["Ra"],REF["R0"],REF["qa"],REF["Bf"],REF["Bg"]
 afu,agu=ck["500000"]
 sf=abs(log2(DK.rho((542.17+afu+0.5)/Bf)/DK.rho((542.17+afu-0.5)/Bf)))
 sg=abs(log2(DK.rho((670.13+agu+0.5)/Bg)/DK.rho((670.13+agu-0.5)/Bg)))
 It=(log(500000.0)+1.0)/(500000.0*(log(2.0)**2))
 sdr,sdh=sqrt((sf*sf*12.0+sg*sg*16.0)*It),sqrt((sf*sf*19.0+sg*sg*23.0+14.0*sf*sg)*It)
 assert 3.0*sdh<abs(REF["cush"])
 dnf=(72.0/q)*log2((Ra**q+R0**q)/4.0);dng=8.0*LOG2_ABS_X+(24.0/q)*log2((R0**q+(Ra**q)/(ABS_X**(q/3.0)))/4.0)
 assert dnf<656.0 and dng<688.0 and 6*656<=4096 and 5*688<=4096 and 10*408<=4096
 tfm=int(floor(656.0/b["qf"]))-127;assert tfm==45 and(128-tfm)+int(floor(656.0/b["qf"]))<=255
 for ngc in RN(610,689):assert(128-max(1,int(floor((ngc+1.0)/b["qg"]))-127))+int(floor(ngc/b["qg"]))<=255
 w0=1.0/(1.0+REF["rbr"]**8);p0,bn0=_cof_ex(Bf,b["pfhi0"],542.17,REF["lf"],afu);p1,bn1=_cof_ex(Bf,b["pfhi1"],542.17,REF["lf"],afu)
 assert abs(bn0-1.0)<0.01 and abs(bn1-1.0)<0.01;pex=w0*p0+(1.0-w0)*p1;cpad=log2(1.0/(1.0-pex))
 return{"ok":True,"pmax":pmax,"af_500k":afu,"ag_500k":agu,"af_end":round(af,5),"ag_end":round(ag,5),"ck":ck,"tail_sd_status":"H1_random_Frobenius_heuristic_NOT_unconditional","tail_500k_sd_root_h1_heur":round(sdr,6),"tail_500k_sd_cov_h1_heur":round(sdh,6),"tail_500k_3sd_cov_h1_heur":round(3.0*sdh,6),"cush_b":REF["cush"],"det_nf":round(dnf,6),"det_ng":round(dng,6),"p_ex_f":round(pex,6),"cof_pad":round(cpad,6),"bnorm":[round(bn0,5),round(bn1,5)]}
VH=((0.481211825060,0.481211825060,-0.481211825060,-0.481211825060),(0.687654016088,-0.206442191029,0.396282847190,-0.877494672250),(0.206442191029,-0.687654016088,0.877494672250,-0.396282847190),(0.877494672250,-0.396282847190,-0.687654016088,0.206442191029),(0.396282847190,-0.877494672250,-0.206442191029,0.687654016088),(1.083936863278,-1.083936863278,0.189840656161,-0.189840656161),(0.189840656161,-0.189840656161,-1.083936863278,1.083936863278))
VHP=tuple((v,0.5*sum(c*c for c in v))for v in VH)
SOB=((1,0,(1,)),(2,1,(1,1)),(3,1,(1,3,1)),(3,2,(1,1,1)),(4,1,(1,1,3,3)),(4,4,(1,3,5,13)),(5,2,(1,1,5,5,17)),(5,4,(1,1,5,5,5)))
def _sobol9(n):
 B=30;V=[[1<<(B-1-k)for k in RN(B)]]
 for s,a,m in SOB:
  mv=list(m)
  for k in RN(s,B):
   v=mv[k-s]^(mv[k-s]<<s)
   for j in RN(1,s):
    if(a>>(s-1-j))&1:v^=(mv[k-j]<<j)
   mv.append(v)
  V.append([mv[k]<<(B-1-k)for k in RN(B)])
 pts,cu=[],[0]*9
 for i in RN(n):
  pts.append(tuple(cu));c=((i+1)&-(i+1)).bit_length()-1
  for d in RN(9):cu[d]^=V[d][c]
 return pts
def _pdp(Bb,plo,phi,qs,ab,bmin,bst,nb,ns_=0.5,pr_=1e-17,TM=260,km=700.0):
 KM=int(round(km/ns_));mo=[0.0]*(KM+1)
 for ai,bi in((0.561459483567,plo),(phi,Bb)):
  if bi<=ai:continue
  ns=max(20,int((bi-ai)*8.0/ns_));dx=(bi-ai)/ns
  for i in RN(ns):
   x0=ai+i*dx;kb=int(round((x0+0.5*dx)/ns_))
   if 1<=kb<=KM:mo[kb]+=log((x0+dx)/x0)
 nzo=[(j,j*mo[j])for j in RN(1,KM+1)if mo[j]>0.0];Po=[0.0]*(KM+1);Po[0]=exp(-sum(mo[1:]))
 for k in RN(1,KM+1):Po[k]=sum(jm*Po[k-j]for j,jm in nzo if j<=k)/k
 mi,m1={},[0.0]*(TM+1);ns=max(80,int((phi-plo)*16.0/ns_));dx=(phi-plo)/ns
 for i in RN(ns):
  x0=plo+i*dx;xm=x0+0.5*dx;w,kb,tk=log((x0+dx)/x0),int(round(xm/ns_)),int(floor(xm/qs))
  if kb>=1 and 1<=tk<=TM:mi[(kb,tk)]=mi.get((kb,tk),0.0)+w;m1[tk]+=w
 li=sum(m1[1:]);nz1=[(j,j*m1[j])for j in RN(1,TM+1)if m1[j]>0.0];P1=[0.0]*(TM+1);P1[0]=exp(-li)
 for k in RN(1,TM+1):P1[k]=sum(jm*P1[k-j]for j,jm in nz1 if j<=k)/k
 sv,rn=[0.0]*(TM+1),0.0
 for t in RN(TM,-1,-1):rn+=P1[t];sv[t]=rn
 its=sorted([(dk,dt,dk*v)for(dk,dt),v in mi.items()]);P2=[None]*(KM+1);P2[0]=[(0,exp(-li))];TM1=TM+1
 for k in RN(its[0][0],KM+1):
  ik,ds,tc=1.0/k,[0.0]*TM1,False
  for dk,dt,dv in its:
   if dk>k:break
   pr=P2[k-dk]
   if pr is None:continue
   sd,lm=dv*ik,TM-dt
   for tt,p0 in pr:
    if tt>lm:break
    ds[tt+dt]+=sd*p0;tc=True
  if tc:
   rp=[(tt,pv)for tt,pv in enumerate(ds)if pv>pr_]
   if rp:P2[k]=rp
 nzk=[(k,P2[k])for k in RN(KM+1)if P2[k]is not None];rt=[]
 for bi in RN(nb):
  kN=max(1,min(KM,int(round((bmin+(bi+0.5)*bst+ab)/ns_))));sm=[0.0]*TM1
  for k,ki in nzk:
   if k>kN:break
   po=Po[kN-k]
   if po<=1e-300:continue
   for tt,pt in ki:sm[tt]+=pt*po
  tot,tl,rw=sum(sm),0.0,[0.0]*129
  for tt in RN(TM,-1,-1):
   tl+=sm[tt]
   if tt<=128:rw[tt]=(tl/tot)if tot>0.0 else 0.0
  rt.append(rw)
 return sv,rt
def check_live_panjer_rqmc_chords_disc(af,ag,n_shifts=24,m_sobol=1024,nch=1440,full_disc=False):
 b,Bf,Bg,qa,rbr,Ra,R0=REF["b"],REF["Bf"],REF["Bg"],REF["qa"],REF["rbr"],REF["Ra"],REF["R0"]
 w0=1.0/(1.0+rbr**8);w1=1.0-w0;G0,GS,GB,F0,FS,FB,N0,NS,NB=610,2,40,510,4,25,440,2,100
 sgL,rgL,sfL,rfL,sgLc,rgLc,sfLc,rfLc,sgLf,rgLf,sfLf,rfLf=[],[],[],[],[],[],[],[],[],[],[],[]
 tmf=320 if full_disc else 260
 for si,(pgh,pfh)in enumerate(((b["pghi0"],b["pfhi0"]),(b["pghi1"],b["pfhi1"]))):
  sg,rg=_pdp(Bg,b["pglo"],pgh,b["qg"],ag,G0,GS,GB,0.5,1e-17,260,695.0);sf,rf=_pdp(Bf,b["pflo"],pfh,b["qf"],af,N0,NS,NB,0.5,1e-17,260,650.0)
  sgL.append(sg);rgL.append(rg);sfL.append(sf);rfL.append(rf)
  if full_disc or si==0:
   sgc,rgc=_pdp(Bg,b["pglo"],pgh,b["qg"],ag,G0,GS,GB,1.0,1e-13,260,695.0);sfc,rfc=_pdp(Bf,b["pflo"],pfh,b["qf"],af,N0,NS,NB,1.0,1e-13,260,650.0)
   sgf,rgf=_pdp(Bg,b["pglo"],pgh,b["qg"],ag,G0,GS,GB,0.25,1e-20,tmf,695.0);sff,rff=_pdp(Bf,b["pflo"],pfh,b["qf"],af,N0,NS,NB,0.25,1e-20,tmf,650.0)
   sgLc.append(sgc);rgLc.append(rgc);sfLc.append(sfc);rfLc.append(rfc);sgLf.append(sgf);rgLf.append(rgf);sfLf.append(sff);rfLf.append(rff)
  else:
   sgLc.append(sg);rgLc.append(rg);sfLc.append(sf);rfLc.append(rf);sgLf.append(sg);rgLf.append(rg);sfLf.append(sf);rfLf.append(rf)
 dk2,dk8=ExactDickman(30.0,200),ExactDickman(30.0,800);ug=[1.5+0.25*k for k in RN(40)]
 d200=max(abs(dk2.rho(u)/dk8.rho(u)-1.0)for u in ug);d400=max(abs(DK.rho(u)/dk8.rho(u)-1.0)for u in ug);d400l=max(abs(log2(DK.rho(u)/dk8.rho(u)))for u in ug)
 al=2.0/qa;Cd=exp(lgamma(4.0*al)-4.0*lgamma(al))
 atb,ata,k0,l0,k1,l1=3.2,4.2,3.2,5.5,2.0,2.5;ik0,il0,ik1,il1,itb,ita=1.0/k0,1.0/l0,1.0/k1,1.0/l1,1.0/atb,1.0/ata
 iqa,il2,tp=1.0/qa,1.0/log(2.0),2.0*pi;iq3,al4=3.0*iqa,4.0*al
 Ap=[(a**6,a**12,a**18)for a in((Ra/(4.0**iqa))*((0.25 if k==0 else 0.75)**0.125)for k in RN(2))]
 ra0c,ra1c,S32,g0c=(Ra*(0.5**0.125))**3,Ra**3,0.8660254037844386,8.0*LOG2_ABS_X
 lg,cs_,sn_,ex_,XF=log,cos,sin,exp,float(X_BLS);Dr,D2r,D8r=DK.rho,dk2.rho,dk8.rho
 def da(rn):
  while True:
   u0,u1,u2=1.36*rn()-0.68,1.36*rn()-0.68,1.36*rn()-0.68;u3=-(u0+u1+u2)
   if abs(u3)>0.68:continue
   ok=True
   for v,h in VHP:
    if abs(u0*v[0]+u1*v[1]+u2*v[2]+u3*v[3])>h:ok=False;break
   if not ok:continue
   e0,e1,e2,e3=ex_(qa*u0),ex_(qa*u1),ex_(qa*u2),ex_(qa*u3);es=e0+e1+e2+e3
   if rn()<=(4.0/es)**al4:ie=1.0/es;return(e0*ie)**iq3,(e1*ie)**iq3,(e2*ie)**iq3,(e3*ie)**iq3
 def sh_ev(ra,rb,a,b_,cs,A,lbp):
  rb3=rb**3;rb6=rb3*rb3;ng0=g0c+24.0*lg(rb)*il2+lbp;fp=gp=pp=1.0;A6,A12,A18=A
  for m in RN(4):
   bm=b_[m];r6=rb6*bm*bm;r12=r6*r6;pp*=A18+A12*r6+A6*r12+r12*r6;ar=ra*a[m]*cs[2*m];ai=ra*a[m]*cs[2*m+1];bt=rb3*bm;gr=ar-XF*bt;gp*=gr*gr+ai*ai+1e-300
   a2r=ar*ar-ai*ai;a2i=2.0*ar*ai;b2_=bt*bt;b3_=b2_*bt;fr=a2r*ar-a2i*ai-a2r*bt+(-0.5*ar+S32*ai)*b2_+0.5*b3_;fi=a2r*ai+a2i*ar-a2i*bt+(-0.5*ai-S32*ar)*b2_-S32*b3_;fp*=fr*fr+fi*fi+1e-300
  return lg(fp)*il2,lg(gp)*il2,ng0,lg(pp)*il2
 def ev2(a,u,cs):
  tb=u[0]**itb;wtb=(8.0/atb)*(tb**(8.0-atb));ta0=u[1]**ita;wta=(8.0/ata)*(ta0**(8.0-ata))
  t0=(1.0-(1.0-u[2])**il0)**ik0;t1=(1.0-(1.0-u[3])**il1)**ik1;t2=(1.0-(1.0-u[4])**il1)**ik1
  p0=k0*l0*(t0**(k0-1.0))*((1.0-t0**k0)**(l0-1.0));p1=k1*l1*(t1**(k1-1.0))*((1.0-t1**k1)**(l1-1.0));p2=k1*l1*(t2**(k1-1.0))*((1.0-t2**k1)**(l1-1.0))
  ub=(t0*t1,t0*(1.0-t1),(1.0-t0)*t2,(1.0-t0)*(1.0-t2));up=ub[0]*ub[1]*ub[2]*ub[3]
  wb=wtb*Cd*(up**(al-1.0))*(t0*(1.0-t0))/(p0*p1*p2);b_=(ub[0]**iq3,ub[1]**iq3,ub[2]**iq3,ub[3]**iq3);lbp=6.0*iqa*lg(up)*il2
  return wb*(w0*wta+w1),((0,w0*wta*wb,*sh_ev(ra0c*(ta0**3),R0*tb,a,b_,cs,Ap[0],lbp)),(1,w1*wb,*sh_ev(ra1c*((0.5+0.5*u[1])**0.375),R0*rbr*tb,a,b_,cs,Ap[1],lbp)))
 rf_=random.Random(20261127);rnf=rf_.random;v2=[[[0.0]*FB for _ in RN(GB)]for _ in RN(2)];s3=[[[[0.0]*NB for _ in RN(FB)]for _ in RN(GB)]for _ in RN(2)];n_fit,tsf=4000,0.0
 for i in RN(n_fit):
  if i%4==0:a_=da(rnf)
  h=[tp*rnf()for _ in RN(4)];cs=(cs_(h[0]),sn_(h[0]),cs_(h[1]),sn_(h[1]),cs_(h[2]),sn_(h[2]),cs_(h[3]),sn_(h[3]))
  _,sr=ev2(a_,(rnf(),rnf(),rnf(),rnf(),rnf()),cs)
  for sh,wi,nf,ng,ng0,nfp in sr:
   ig=max(0,min(GB-1,int((ng0-G0)//GS)));icf=max(0,min(FB-1,int((nfp-F0)//FS)));iff=max(0,min(NB-1,int((nf-N0)//NS)))
   wt=wi*Dr((nf+af)/Bf)*Dr((ng+ag)/Bg);v2[sh][ig][icf]+=wi;s3[sh][ig][icf][iff]+=wt;tsf+=wt
 ivn,ivt,tfm,lp=1.0/n_fit,1.0/tsf,max(0,int(floor(656.0/b["qf"]))-127),REF["lam"]
 Tg,Tf=[[[128]*FB for _ in RN(GB)]for _ in RN(2)],[[[125]*FB for _ in RN(GB)]for _ in RN(2)]
 def hull(rw,sfs):
  h=[]
  for tf in RN(125,tfm-1,-1):
   x,y=sfs[tf],rw[tf]
   while len(h)>=2 and(h[-1][2]-h[-2][2])*(x-h[-1][1])<=(y-h[-1][2])*(h[-1][1]-h[-2][1]):h.pop()
   h.append((tf,x,y))
  return h
 for sh in RN(2):
  rgs,rfs,sgs,sfs=rgL[sh],rfL[sh],sgL[sh],sfL[sh]
  for ig in RN(GB):
   ngu=G0+(ig+1)*GS+1.0;tgm=max(0,min(128,int(floor(ngu/b["qg"]))-127));vw1=sum(v2[sh][ig])*ivn
   nz1=[(iff,sv)for iff in RN(NB)if(sv:=sum(s3[sh][ig][icf][iff]for icf in RN(FB))*ivt)>0.0]
   h1=hull([sum(sv*rfs[iff][tf]for iff,sv in nz1)if(tf>=tfm and vw1>0.0)else 0.0 for tf in RN(126)],sfs)if vw1>0.0 else[]
   for icf in RN(FB):
    vw=v2[sh][ig][icf]*ivn
    if vw<5e-5:
     if vw1<=0.0:continue
     vwu,hc=vw1,h1
    else:
     vwu=vw;nzf=[(iff,s3[sh][ig][icf][iff]*ivt)for iff in RN(NB)if s3[sh][ig][icf][iff]>0.0]
     hc=hull([sum(sv*rfs[iff][tf]for iff,sv in nzf)if tf>=tfm else 0.0 for tf in RN(126)],sfs)
    bo,btg,btf=-1e300,128,125
    for tg in RN(tgm,129):
     rgv,sgt=rgs[ig][tg],lp*vwu*sgs[tg]
     for tf,sfx,rfy in hc:
      o=rgv*rfy-sgt*sfx
      if o>bo:bo,btg,btf=o,tg,tf
    Tg[sh][ig][icf],Tf[sh][ig][icf]=btg,btf
    assert(128-btg)+int(floor(ngu/b["qg"]))<=255 and(128-btf)+int(floor(656.0/b["qf"]))<=255
 N_B=n_shifts;sob=_sobol9(m_sobol);re_=random.Random(20261128);rne=re_.random;s30=1.0/(1<<30)
 b_wk,b_wr,b_sg,b_s1,b_wc,b_wc_,b_wf,b_wrf=[],[],[],[],[],[],[],[];mfs=mgs=swn=nfm=ngm=0.0;ncm=16 if full_disc else 2;iBf,iBg=1.0/Bf,1.0/Bg
 for si in RN(N_B):
  qs=[re_.getrandbits(30)for _ in RN(9)];js=[rne()*s30 for _ in RN(9)]
  swk=swr=ssg=ss1=swc=swkc=swkf=swrf=0.0;dc=(si<ncm)
  for pt in sob:
   a_=da(rne);u=[(pt[d]^qs[d])*s30+js[d]for d in RN(9)]
   c0,s0,c1,s1,c2,s2,c3,s3_=cs_(tp*u[5]),sn_(tp*u[5]),cs_(tp*u[6]),sn_(tp*u[6]),cs_(tp*u[7]),sn_(tp*u[7]),cs_(tp*u[8]),sn_(tp*u[8])
   for uu,cs in((u,(c0,s0,c1,s1,c2,s2,c3,s3_)),((1.0-u[0],1.0-u[1],1.0-u[2],1.0-u[3],1.0-u[4]),(-c0,-s0,-c1,-s1,-c2,-s2,-c3,-s3_))):
    wcp,sr=ev2(a_,uu,cs);swc+=0.5*wcp
    for sh,ws_,nf,ng,ng0,nfp in sr:
     if nf>nfm:nfm=nf
     if ng>ngm:ngm=ng
     ig=int((ng0-G0)//GS);ig=0 if ig<0 else(GB-1 if ig>=GB else ig)
     icf=int((nfp-F0)//FS);icf=0 if icf<0 else(FB-1 if icf>=FB else icf)
     iff=int((nf-N0)//NS);iff=0 if iff<0 else(NB-1 if iff>=NB else iff)
     iga=int((ng-G0)//GS);iga=0 if iga<0 else(GB-1 if iga>=GB else iga)
     wi=0.5*ws_;uf,ug_=(nf+af)*iBf,(ng+ag)*iBg
     wt=wi*Dr(uf)*Dr(ug_);tg,tf=Tg[sh][ig][icf],Tf[sh][ig][icf];wk=wt*rgL[sh][iga][tg]*rfL[sh][iff][tf];sp=sgL[sh][tg]
     swr+=wt;swk+=wk;ssg+=wi*sp;ss1+=wi*sp*sfL[sh][tf];mfs+=wk*nf;mgs+=wk*ng;swn+=wk
     if dc:
      swkc+=wi*D2r(uf)*D2r(ug_)*rgLc[sh][iga][tg]*rfLc[sh][iff][tf]
      wtf=wi*D8r(uf)*D8r(ug_);swrf+=wtf;swkf+=wtf*rgLf[sh][iga][tg]*rfLf[sh][iff][tf]
  b_wk.append(swk/m_sobol);b_wr.append(swr/m_sobol);b_sg.append(ssg/m_sobol);b_s1.append(ss1/m_sobol);b_wc.append(swc/m_sobol)
  if dc:b_wc_.append(swkc/m_sobol);b_wf.append(swkf/m_sobol);b_wrf.append(swrf/m_sobol)
 mwk,mwr,msg,ms1,ewi=sum(b_wk)/N_B,sum(b_wr)/N_B,sum(b_sg)/N_B,sum(b_s1)/N_B,sum(b_wc)/N_B
 se_Vk=REF["sl"]*(sqrt(sum((x-mwk)**2 for x in b_wk)/(N_B*(N_B-1)))/mwk)*il2
 volu=(pi**4)*(gamma(1.0+2.0/qa)**4)/gamma(1.0+8.0/qa)
 Vkb=log2(V_A_SH*(volu*(R0**8)/(70.3125*3.0))*(1.0+rbr**8));U=U_ub(Bf,Bg)
 Vkl=Vkb+REF["sl"]*(log2(REF["eta"]*U)-(Vkb+log2(mwk*0.97843406*(2.0**REF["cush"]))))
 retL,l2sg,l2s1=mwk/mwr,log2(msg),log2(ms1);l2sf=l2s1-l2sg
 mwkc,mwrc=sum(b_wk[:ncm])/ncm,sum(b_wr[:ncm])/ncm
 dVk_cp=-REF["sl"]*log2(mwkc/(sum(b_wc_)/ncm));dVk_pf=-REF["sl"]*log2((sum(b_wf)/ncm)/mwkc);dret_pf=(sum(b_wf)/sum(b_wrf))-(mwkc/mwrc)
 dpad=REF["disc_pad"];assert abs(ewi-1.0)<0.005 and nfm<656.0 and ngm<688.0 and se_Vk<0.0032 and dpad>=2.5*abs(dVk_pf)
 rc_,hq,agd,NC_B=random.Random(20260505),0.5*qa,2.0/qa,nch//N_B;bL4,bLh,bg4=[],[],[]
 for _ in RN(N_B):
  L4b=Lhb=g4b=0.0
  for R,wsh in((R0,w0),(R0*rbr,w1)):
   Rq,ssp,sgp,w248,w124=R**qa,0.0,0.0,0.0,0.0;ns,tl=max(24,int(round(NC_B*wsh))),int(ceil(R*1.15))+20
   for _ in RN(ns):
    gs=[rc_.gammavariate(agd,1.0)for _ in RN(4)];igs=(Rq*(max(rc_.random(),1e-16)**(qa/8.0)))/sum(gs)
    zr,zi2=[],[]
    for gm_ in gs:rm,th=(gm_*igs)**iqa,tp*rc_.random();zr.append(rm*cs_(th));zi2.append((rm*sn_(th))**2)
    z0r,z1r,z2r,z3r=zr;z0i,z1i,z2i,z3i=zi2;cmn,cmx=10000,-10000
    for tc in RN(-tl,tl+1,16):
     s=((z0r+tc)**2+z0i)**hq+((z1r+tc)**2+z1i)**hq
     if s<=Rq and(s:=s+((z2r+tc)**2+z2i)**hq)<=Rq and s+((z3r+tc)**2+z3i)**hq<=Rq:
      if tc<cmn:cmn=tc
      if tc>cmx:cmx=tc
    tlo,thi=(max(-tl,cmn-32),min(tl,cmx+32))if cmn<=cmx else(-256,256);tmn,tmx,ic=10000,-10000,0
    for t in RN(tlo,thi+1):
     s=((z0r+t)**2+z0i)**hq
     if s<=Rq and(s:=s+((z1r+t)**2+z1i)**hq)<=Rq and(s:=s+((z2r+t)**2+z2i)**hq)<=Rq and s+((z3r+t)**2+z3i)**hq<=Rq:
      if t<tmn:tmn=t
      if t>tmx:tmx=t
      ic+=1
    if ic:sp_=tmx-tmn+1;w=1.0/ic;ssp+=w*sp_;sgp+=w*(sp_-ic);w248+=w*((sp_+247)//248);w124+=w*((sp_+123)//124)
   L4b+=wsh*(ssp/w248);Lhb+=wsh*(ssp/w124);g4b+=wsh*(sgp/ssp)
  bL4.append(L4b);bLh.append(Lhb);bg4.append(g4b)
 L4m,Lhm,g4m=sum(bL4)/N_B,sum(bLh)/N_B,sum(bg4)/N_B
 se4,seh,seg=(sqrt(sum((x-m)**2 for x in arr)/(N_B*(N_B-1)))for arr,m in((bL4,L4m),(bLh,Lhm),(bg4,g4m)))
 L4c,Lhc,g4c=L4m-3.0*se4,Lhm-3.0*seh,g4m+3.0*seg;gpad=log2(1.0/(1.0-g4c));assert 95.0<L4c<120.0 and 68.0<Lhc<88.0 and 0.0<g4c<0.008
 return{"ok":True,"sampler":"unbiased_antithetic_9D_Sobol_RQMC","n_fit":n_fit,"n_shifts":N_B,"n_eval_rqmc":N_B*m_sobol*2,"n_chords":nch,"E_w_is":round(ewi,5),"ret_live":round(retL,6),"l2_sg":round(l2sg,6),"l2_sf":round(l2sf,6),"l2_s1":round(l2s1,6),"mean_nf":round(mfs/swn,3),"mean_ng":round(mgs/swn,3),"obs_max":[round(nfm,2),round(ngm,2)],"Vk_base":round(Vkb,6),"Vk_live":round(Vkl,6),"se_Vk":round(se_Vk,6),"pad_3se_Vk":round(3.0*se_Vk,6),"L248_live":round(L4m,3),"se_L248":round(se4,3),"L248_ch":round(L4c,3),"L124_live":round(Lhm,3),"se_L124":round(seh,3),"L124_ch":round(Lhc,3),"gap4_live":round(g4m,6),"se_gap4":round(seg,6),"gap_ch":round(g4c,6),"gap_pad":round(gpad,6),"discretization_convergence":{"ok":True,"primary_grid":"ExactDickman(N=400)+Panjer(0.5b,prune=1e-17,TM=260)","fine_grid":f"ExactDickman(N=800)+Panjer(0.25b,prune=1e-20,TM={tmf})","fine_shells":"Shell0+Shell1"if full_disc else f"Shell0_0.25b_1e-20_TM{tmf}","crn_eval_pts":ncm*m_sobol*2,"dickman_N200_vs_800_max_rel":float(f"{d200:.3e}"),"dickman_N400_vs_800_max_rel":float(f"{d400:.3e}"),"dickman_N400_vs_800_max_l2":float(f"{d400l:.3e}"),"dVk_coarse1b1e13_to_prim05b1e17":round(dVk_cp,6),"dVk_prim_to_fine025b1e20":round(dVk_pf,6),"dret_prim_to_fine":round(dret_pf,6),"disc_pad_b":dpad,"pad_over_shift_x":round(dpad/max(1e-9,abs(dVk_pf)),2),"scope":"numerical_DP_Dickman_only_H1_exploratory_v1"}}
def _ip64(n):
 d,s=n-1,0
 while d%2==0:d//=2;s+=1
 for a in(2,3,5,7,11,13,17,19,23,29,31,37):
  x=pow(a,d,n)
  if x in(1,n-1):continue
  for _ in RN(s-1):
   x=(x*x)%n
   if x==n-1:break
  else:return False
 return True
_sh=lambda a,w,wb,p:(a*w-((a*wb)>>64)*p)&0xFFFFFFFFFFFFFFFF
def _ntt_r4(vec,p,g):
 N,a,p2,p4,l4,t=len(vec),[x%p for x in vec],p<<1,p<<2,0,len(vec)
 while t>1:t>>=2;l4+=1
 for i in RN(N):
  rv,x=0,i
  for _ in RN(l4):rv=(rv<<2)|(x&3);x>>=2
  if i<rv:a[i],a[rv]=a[rv],a[i]
 J=pow(g,(p-1)//4,p);Jp,m=(J<<64)//p,1
 while m<N:
  wm,w1,W=pow(g,(p-1)//(4*m),p),1,[]
  for _ in RN(m):w2=(w1*w1)%p;w3=(w2*w1)%p;W.append((w1,(w1<<64)//p,w2,(w2<<64)//p,w3,(w3<<64)//p));w1=(w1*wm)%p
  for k in RN(0,N,4*m):
   for j in RN(m):
    w1,w1p,w2,w2p,w3,w3p=W[j];x0=a[k+j]-(p2 if a[k+j]>=p2 else 0)
    x1,x2,x3=_sh(a[k+m+j],w1,w1p,p),_sh(a[k+2*m+j],w2,w2p,p),_sh(a[k+3*m+j],w3,w3p,p)
    t0,t1=(x0+x2)-(p2 if x0+x2>=p2 else 0),(x0+p2-x2)-(p2 if x0+p2-x2>=p2 else 0)
    t2,t3=(x1+x3)-(p2 if x1+x3>=p2 else 0),_sh(x1+p2-x3,J,Jp,p)
    a[k+j],a[k+2*m+j],a[k+m+j],a[k+3*m+j]=t0+t2,t0+p2-t2,t1+t3,t1+p2-t3
    assert max(a[k+j],a[k+m+j],a[k+2*m+j],a[k+3*m+j])<p4
  m<<=2
 return[x%p for x in a]
def _ntt1(vec,p,g,inv=False):
 n,a,p2,j=len(vec),[x%p for x in vec],2*p,0
 for i in RN(1,n):
  b=n>>1
  while j&b:j^=b;b>>=1
  j^=b
  if i<j:a[i],a[j]=a[j],a[i]
 L=2
 while L<=n:
  hf,wl=L>>1,pow(g,(p-1)//L,p)
  if inv:wl=pow(wl,p-2,p)
  ws,wp,cu=[1]*hf,[0]*hf,1
  for k in RN(hf):ws[k],wp[k]=cu,(cu<<64)//p;cu=(cu*wl)%p
  for i in RN(0,n,L):
   for k in RN(hf):u=a[i+k]-(p2 if a[i+k]>=p2 else 0);t=_sh(a[i+k+hf],ws[k],wp[k],p);a[i+k],a[i+k+hf]=u+t,u+p2-t
  L<<=1
 return[(x*pow(n,p-2,p))%p for x in a]if inv else[x%p for x in a]
def _mul2d(A,B,pgs):
 m58=(1<<58)-1
 def t58(x):
  r=[]
  while x:r.append(x&m58);x>>=58
  return r or[0]
 la,lb=t58(A),t58(B);k=max(1,(max(len(la),len(lb))-1).bit_length());k1=(k+1)//2;D1,D2=1<<k1,1<<(k-k1);N1,N2=2*D1,2*D2
 la+=[0]*(D1*D2-len(la));lb+=[0]*(D1*D2-len(lb));cm=[]
 for p,g in pgs:
  gA=[_ntt1(la[r*D1:(r+1)*D1]+[0]*D1,p,g)for r in RN(D2)]+[[0]*N1 for _ in RN(D2)]
  gB=[_ntt1(lb[r*D1:(r+1)*D1]+[0]*D1,p,g)for r in RN(D2)]+[[0]*N1 for _ in RN(D2)]
  gC=[[0]*N1 for _ in RN(N2)]
  for c in RN(N1):
   cA,cB=_ntt1([gA[r][c]for r in RN(N2)],p,g),_ntt1([gB[r][c]for r in RN(N2)],p,g)
   cC=_ntt1([(cA[r]*cB[r])%p for r in RN(N2)],p,g,True)
   for r in RN(N2):gC[r][c]=cC[r]
  cm.append([_ntt1(gC[r],p,g,True)for r in RN(N2)])
 p1,p2,p3=pgs[0][0],pgs[1][0],pgs[2][0];i12,i123,p12=pow(p1%p2,p2-2,p2),pow((p1*p2)%p3,p3-2,p3),p1*p2;acc=[0]*(D1*(N2+1)+4)
 for r in RN(N2):
  bs=r*D1
  for c in RN(N1):a1,a2,a3=cm[0][r][c],cm[1][r][c],cm[2][r][c];x12=a1+p1*(((a2-a1)*i12)%p2);acc[bs+c]+=x12+p12*(((a3-x12)*i123)%p3)
 rs=cy=0
 for i,v in enumerate(acc):s=v+cy;rs|=(s&m58)<<(58*i);cy=s>>58
 return rs|(cy<<(58*len(acc)))
def check_ntt():
 s=REF["s1b"];p1,p2,p3,g1,g2,g3=s["p1"],s["p2"],s["p3"],s["g1"],s["g2"],s["g3"];l2P=log2(p1*p2*p3);h1,h2=l2P-181.0,l2P-(log2(ceil((2.0**67.77)/58.0))+116.0)
 assert all(_ip64(p)and(p-1)%(1<<54)==0 and 8*p<(1<<64)and(1<<58)<p for p in(p1,p2,p3))and all(pow(g1,(p1-1)//q,p1)!=1 for q in(2,127))and all(pow(g2,(p2-1)//q,p2)!=1 for q in(2,3,19))and all(pow(g3,(p3-1)//q,p3)!=1 for q in(2,3))and h1>1.57 and h2>4.66
 rng,pgs=random.Random(0x5400612D),[(p1,g1),(p2,g2),(p3,g3)]
 for p,g in pgs:
  vec=[rng.randrange(p)for _ in RN(64)];got=_ntt_r4(vec,p,g);wN=pow(g,(p-1)//64,p)
  for k in RN(64):
   wk,ac,wj=pow(wN,k,p),0,1
   for j in RN(64):ac=(ac+vec[j]*wj)%p;wj=(wj*wk)%p
   assert got[k]==ac
 sz=[512,1024,2048,4096]
 for nb in sz:A,B=rng.getrandbits(nb)|(1<<(nb-1)),rng.getrandbits(nb)|(1<<(nb-1));assert _mul2d(A,B,pgs)==A*B
 return{"ok":True,"r4_bitexact":True,"log2_P":round(l2P,6),"hr_1d":round(h1,4),"hr_2d":round(h2,4),"mul_2d":sz}
def _h576(A,sm,v):
 r256,r2_256,x2h,m64=R_BLS<<256,(2*R_BLS)<<256,0xAC45A4010001A402,(1<<64)-1
 T=A[0]+(A[2]<<32)+((A[2]*x2h)<<64)+r256-ABS_X*(A[1]+(A[3]<<32)+((A[3]*x2h)<<64))+sum(sm[k]*v[k]for k in RN(4))
 assert T>=0;t8=T>>512;assert t8 in(0,1)
 if t8:T-=r2_256
 assert 0<=T<(1<<512);Th=T>>256;t4,t5,t6,t7=Th&m64,(Th>>64)&m64,(Th>>128)&m64,(Th>>192)&m64
 M=(0x42737A020C0D6393,0x65043EB4BE4BAD71,0x38B5DCB707E08ED3,0x355094EDFEDE377C)
 ph=((t4*M[3])+(t5*M[2]+((t5*M[3])<<64))+(t6*M[1]+((t6*M[2])<<64)+((t6*M[3])<<128))+(t7*M[0]+((t7*M[1])<<64)+((t7*M[2])<<128)+((t7*M[3])<<192)))>>64
 s=(T-((Th<<1)+ph)*R_BLS)&((1<<320)-1);assert 0<=s<6*R_BLS;s2=s-(s>>255)*R_BLS;assert 0<=s2<2*R_BLS
 return s2-R_BLS if s2>=R_BLS else s2
def check_kernels_e2e():
 pl=[((i&31)<<56)|(i*131+5)for i in RN(17)];bo=[0,0,1,4,4,8,17];tx=[SLAB_BASE+j*32*SLAB_BASE+pl[k]for j in RN(6)for k in RN(bo[j],bo[j+1])]
 assert len(tx)==17 and tx[-1]==SLAB_BASE+(5*32+16)*SLAB_BASE+(16*131+5)and(900+200)%1009==91 and(91>>5)==2
 Va,pb=1024,37;sl,sh=bytearray((pb+(i*17)%60)for i in RN(Va)),bytearray((i*29)%50 for i in RN(Va))
 sl[13],sh[13],sl[14],sh[14],sl[255],sh[255],sl[256],sh[256]=80,47,80,48,150,105,60,67
 ex=[i for i in RN(Va)if sl[i]+sh[i]>=128];assert all(sl[i]+sh[i]<=255 for i in RN(Va))
 rx=bytearray(sl);sl[:]=bytes([pb])*Va;gs=[]
 for blk in RN(0,Va,128):
  ac=0
  for w in RN(16):ac|=(int.from_bytes(rx[blk+8*w:blk+8*w+8],"little")+int.from_bytes(sh[blk+8*w:blk+8*w+8],"little"))&((1<<64)-1)
  if ac&0x8080808080808080:
   for w in RN(16):
    m=(int.from_bytes(rx[blk+8*w:blk+8*w+8],"little")+int.from_bytes(sh[blk+8*w:blk+8*w+8],"little"))&0x8080808080808080
    while m:gs.append(blk+8*w+(((m&-m).bit_length()-1)>>3));m&=m-1
  sh[blk:blk+128]=bytes(128)
 assert gs==ex and not any(sh)and all(x==pb for x in sl)
 rng=random.Random(0x576B1750)
 for t in RN(60):
  A,sm,vv=([(1<<304)-1,0,(1<<304)-1,0],[R_BLS-1]*4,[R_BLS-1]*4)if t==0 else(([0,(1<<314)-1,0,(1<<314)-1],[0]*4,[0]*4)if t==1 else(([(1<<300)-1]*4,[R_BLS-1]*4,[R_BLS-1]*4)if t==2 else([rng.getrandbits(314)for _ in RN(4)],[rng.randrange(R_BLS)for _ in RN(4)],[rng.randrange(R_BLS)for _ in RN(4)])))
  assert _h576(A,sm,vv)==(A[0]+X_BLS*A[1]+X_BLS*X_BLS*A[2]+(X_BLS**3)*A[3]+sum(sm[k]*vv[k]for k in RN(4)))%R_BLS
 Vv,r60=[rng.randrange(R_BLS)for _ in RN(256)],R_BLS<<60;z1=X_BLS%R_BLS;z2=(z1*z1)%R_BLS;z3=(z2*z1)%R_BLS;te=c1=cu=0
 for _ in RN(60):
  Ap,Am,rr=[0]*4,[0]*4,0
  for ei in RN(60):
   col=rng.randrange(256)
   tp_=MONO[ei]if ei<12 else(MONO[rng.randrange(12)]if rng.random()<0.70 else tuple((0 if rng.random()<0.35 else(1 if rng.random()<0.5 else-1)*(1 if rng.random()<0.4 else rng.randrange(2,1<<29)))for _ in RN(4)))
   te+=1;rr=(rr+(tp_[0]+tp_[1]*z1+tp_[2]*z2+tp_[3]*z3)*Vv[col])%R_BLS
   for k,ck in enumerate(tp_):
    if ck:
     if abs(ck)==1:c1+=1
     else:cu+=1
     if ck>0:Ap[k]+=ck*Vv[col]
     else:Am[k]+=(-ck)*Vv[col]
  As=[Ap[k]+r60-Am[k]for k in RN(4)];assert all(0<=Ap[k]<(1<<293)and 0<=Am[k]<(1<<293)and 0<=As[k]<(1<<315)for k in RN(4))
  sm,vv=[rng.randrange(R_BLS)for _ in RN(4)],[rng.randrange(R_BLS)for _ in RN(4)]
  assert _h576(As,sm,vv)==(rr+sum(sm[j]*vv[j]for j in RN(4)))%R_BLS
 return{"ok":True,"chord_tx_hits":len(tx),"swar_surv":len(gs),"zero_carry_255":True,"horner576":"60/60","z12_4tuple_spmv":"60/60","comp_per_ent":round((c1+cu)/te,3)}
def check_genuine_tower_dlp_sm():
 p,r,m=257,43,16;r2,ep=r*r,r*r-1;assert p%3==2 and(p+1)%r==0 and m*m+1==p
 zm=lambda x,y:(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1]);zn=lambda x:x[0]*x[0]-x[0]*x[1]+x[1]*x[1]
 wm=lambda x,y,M:((x[0]*y[0]-x[1]*y[1])%M,(x[0]*y[1]+x[1]*y[0]-x[1]*y[1])%M)
 wi=lambda x,M:(lambda iv:(((x[0]-x[1])*iv)%M,((-x[1])*iv)%M))(pow(zn(x)%M,-1,M))
 def wg1(x,y):
  while y!=(0,0):ny=zn(y);nm=zm(x,(y[0]-y[1],-y[1]));q=((2*nm[0]+ny)//(2*ny),(2*nm[1]+ny)//(2*ny));t=zm(q,y);x,y=y,(x[0]-t[0],x[1]-t[1])
  return zn(x)==1
 def fm(x,y,M):
  a,b,c,d=wm(x[0],y[0],M),wm(x[1],y[1],M),wm(x[0],y[1],M),wm(x[1],y[0],M)
  return[((a[0]-b[0])%M,(a[1]-b[1])%M),((c[0]+d[0])%M,(c[1]+d[1])%M)]
 def fp(x,e,M):
  rs,bs=[(1,0),(0,0)],[(z[0]%M,z[1]%M)for z in x]
  while e:
   if e&1:rs=fm(rs,bs,M)
   bs=fm(bs,bs,M);e>>=1
  return rs
 assert fp([(3,2),(-1,5)],ep,r)==[(1,0),(0,0)]
 def smq(pw):assert(pw[0][0]-1)%r==0 and pw[0][1]%r==0 and pw[1][0]%r==0 and pw[1][1]%r==0;return[((pw[0][0]-1)//r)%r,(pw[0][1]//r)%r,(pw[1][0]//r)%r,(pw[1][1]//r)%r]
 smf=lambda a,b:smq(fp([(a[0]%r2,a[1]%r2),((-b[0])%r2,(-b[1])%r2)],ep,r2))
 assert smf((2,0),(-1,-2))[3]!=0
 s1,s2=smf((3,1),(2,-1)),smf((5,-2),(1,3))
 assert smq(fp(fm([(3,1),(-2,1)],[(5,-2),(-1,-3)],r2),ep,r2))==[(s1[i]+s2[i])%r for i in RN(4)]and smq(fp(fp([(3,1),(-2,1)],r,r2),ep,r2))==[0,0,0,0]
 ps=_primes(42);sp3,frs,fri={},{},{}
 for q in ps:
  if q%3==1:w0=next(t for t in RN(1,q)if(t*t+t+1)%q==0);sp3[q]=(min(w0,(q-1-w0)%q),max(w0,(q-1-w0)%q))
  if q>3:
   if q%12==1:iq=next(t for t in RN(1,q)if(t*t+1)%q==0);frs[q]=min(iq,q-iq)
   elif q%3==2:fri[q]=(min(t:=next(t for t in RN(1,q)if(t*t+1)%q==0),q-t),0)if q%4==1 else((x0:=next(t for t in RN(1,q)if(3*t*t-1)%q==0)),(2*x0)%q)
 fa,un=sorted(list(frs)+list(fri)),((1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1))
 def hv(ux,vx,wv,q,eq):
  if(ux+vx*wv)%q:return 0
  wl,qk=wv,q
  for _ in RN(eq-1):qk*=q;wl=(wl-(wl*wl+wl+1)*pow((2*wl+1)%qk,-1,qk))%qk
  ev,pk=0,q
  while ev<eq and(ux+vx*wl)%pk==0:ev+=1;pk*=q
  return ev
 def fg_(vK):
  rm,ux,vx,o=zn(vK),vK[0],vK[1],{}
  if rm==0:return None
  for q in ps:
   if rm%q==0:
    eq=0
    while rm%q==0:rm//=q;eq+=1
    if q==3:continue
    if q in sp3:
     e0,e1=hv(ux,vx,sp3[q][0],q,eq),hv(ux,vx,sp3[q][1],q,eq)
     if e0+e1!=eq:return None
     if e0!=e1:o[("g",q)]=e0-e1
    elif eq%2 or ux%(qe:=q**(eq//2))or vx%qe:return None
  return o if rm==1 else None
 def ff_(a,b,vK):
  rm,ux,vx,o=zn(vK),vK[0],vK[1],{}
  if rm==0:return None
  for q in fa:
   if rm%q==0:
    eq=0
    while rm%q==0:rm//=q;eq+=1
    if q in frs:
     gt=0
     for ix,wv in enumerate(sp3[q]):
      ev=hv(ux,vx,wv,q,eq)
      if ev:
       gt+=ev;bq=(b[0]+b[1]*wv)%q
       if bq==0:return None
       k=("fs",q,((a[0]+a[1]*wv)*pow(bq,-1,q))%q);o[k]=o.get(k,0)+(ev if ix==0 else-ev)
     if gt!=eq:return None
    else:
     if eq%2 or ux%(qe:=q**(eq//2))or vx%qe or(b[0]%q==0 and b[1]%q==0):return None
     if q%4==1:continue
     o[("fi",q)]=o.get(("fi",q),0)+((eq//2)if wm(a,wi(b,q),q)==fri[q]else-(eq//2))
  return o if rm==1 else None
 rels,sn=[],set()
 for br in RN(-4,5):
  for bi in RN(5):
   if br<=0 and bi==0:continue
   b=(br,bi);b2=zm(b,b);mb=(m*br,m*bi)
   for ar in RN(-8,9):
    for ai in RN(-8,9):
     if ar==0 and ai==0:continue
     a=(ar,ai);a2=zm(a,a);nfK=(a2[0]+b2[0],a2[1]+b2[1])
     if zn(nfK)%2==0:continue
     ff=ff_(a,b,nfK)
     if ff is None:continue
     fg=fg_((ar-mb[0],ai-mb[1]))
     if fg is None or not wg1(a,b):continue
     cb=min((zm(u,a),zm(u,b))for u in un)
     if cb not in sn:sn.add(cb);rels.append((a,b,fg,ff))
 while True:
  cn={}
  for _,_,fg,ff in rels:
   for k,v in(*fg.items(),*ff.items()):
    if v%r:cn[k]=cn.get(k,0)+1
  nr=[rl for rl in rels if all(cn.get(k,0)>=2 for k,v in(*rl[2].items(),*rl[3].items())if v%r)]
  if len(nr)==len(rels):break
  rels=nr
 fbi=sorted({k for _,_,fg,ff in rels for k,v in(*fg.items(),*ff.items())if v%r});fbm={k:i for i,k in enumerate(fbi)};nfb=len(fbi);sml=[smf(a,b)[3]for a,b,_,_ in rels]
 def rref(mt,nc):
  pv={}
  for rw in mt:
   cu={c:v%r for c,v in rw.items()if v%r}
   while cu:
    pc=min(cu)
    if pc in pv:
     fc=cu[pc]
     for c,val in pv[pc].items():
      nv=(cu.get(c,0)-fc*val)%r
      if nv:cu[c]=nv
      else:cu.pop(c,None)
    else:iv=pow(cu[pc],-1,r);pv[pc]={c:(v*iv)%r for c,v in cu.items()};break
  fr=[c for c in RN(nc)if c not in pv];bs=[]
  for fc in fr:
   v=[0]*nc;v[fc]=1
   for pc in sorted(pv,reverse=True):v[pc]=(-sum(vv*v[c]for c,vv in pv[pc].items()if c!=pc))%r
   bs.append(v)
  return len(fr),bs
 mno,msm=[],[]
 for(_,_,fg,ff),sv in zip(rels,sml):
  d={}
  for k,v in fg.items():
   if v%r:d[fbm[k]]=(d.get(fbm[k],0)+v)%r
  for k,v in ff.items():
   if v%r:d[fbm[k]]=(d.get(fbm[k],0)-v)%r
  mno.append(d);d2=dict(d)
  if sv%r:d2[nfb]=(-sv)%r
  msm.append(d2)
 nno,_=rref(mno,nfb);nsm,ksm=rref(msm,nfb+1);assert nno==0 and nsm==1 and ksm[0][nfb]!=0;kv=ksm[0]
 def wp(x,e):
  rs,bs=(1,0),(x[0]%p,x[1]%p)
  while e:
   if e&1:rs=wm(rs,bs,p)
   bs=wm(bs,bs,p);e>>=1
  return rs
 cof=(p*p-1)//r;gn=next(wp((t,1),cof)for t in RN(1,p)if wp((t,1),cof)!=(1,0));dl,cu={},(1,0)
 for e in RN(r):dl[cu]=e;cu=wm(cu,gn,p)
 rt=set()
 for(a,b,fg,ff),sv in zip(rels,sml):
  td=dl[wp((a[0]-m*b[0],a[1]-m*b[1]),cof)]
  pg=sum(v*kv[fbm[k]]for k,v in fg.items()if v%r)%r;pf=(sum(v*kv[fbm[k]]for k,v in ff.items()if v%r)+sv*kv[nfb])%r;assert pg==pf
  if td:rt.add((pg*pow(td,-1,r))%r)
  else:assert pg==0
 assert len(rt)==1;gg=next(iter(rt));assert gg!=0;oos=[]
 for ar in RN(9,18):
  for ai in RN(-4,5):
   fg=fg_((ar,ai))
   if fg and all(k in fbm for k,v in fg.items()if v%r)and(td:=dl[wp((ar,ai),cof)]):
    oos.append((td,sum(v*kv[fbm[k]]for k,v in fg.items()if v%r)%r))
    if len(oos)==2:break
  if len(oos)==2:break
 assert len(oos)==2 and(oos[1][1]*pow(oos[0][1],-1,r))%r==(oos[1][0]*pow(oos[0][0],-1,r))%r
 return{"ok":True,"tower":"K=Q(w),p=257,r=43|(p+1),g=Y-16,f=Y^2+1(L_f=Q(z12),rk1)","planted_vector":False,"real_relations":len(rels),"n_fb_orbits":nfb,"n_sm_cols":1,"unit_sm_nonzero":True,"sm_hom_and_rth_zero":True,"nullity_without_sm":nno,"nullity_with_sm":nsm,"sm_coeff_nonzero":True,"all_rels_exact_dlp_gauge":gg,"oos_target_dlog_ratio_exact":True}
MONO=((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,0,1,0),(0,-1,0,1),(-1,0,0,0),(0,-1,0,0),(0,0,-1,0),(0,0,0,-1),(1,0,-1,0),(0,1,0,-1))
MMAP={v:i for i,v in enumerate(MONO)}
def _zm(a,k):
 k%=12;ng=k>=6
 if ng:k-=6
 c0,c1,c2,c3=a
 for _ in RN(k):c0,c1,c2,c3=-c3,c0,c1+c3,c2
 return(-c0,-c1,-c2,-c3)if ng else(c0,c1,c2,c3)
def _sgep(U,eta,sd):
 rng=random.Random(sd);Uf,lBf,lBg=int(round(0.4392*U)),REF["Bf"]*log(2.0),REF["Bg"]*log(2.0);Ug,rf1,rg1=U-Uf,DK.rho(6.65),DK.rho(8.31)
 ws=[min(0.45,(1.173/U)/((j+0.5)/Uf)*(DK.rho(7.65-max(0.02,1.0+log((j+0.5)/Uf)/lBf))/rf1)+(0.35/(2.0+2.5*j)if j<40 else 0.0))for j in RN(Uf)]+[min(0.45,(1.196/U)/((j+0.5)/Ug)*(DK.rho(9.31-max(0.02,1.0+log((j+0.5)/Ug)/lBg))/rg1)+(0.40/(2.0+2.5*j)if j<40 else 0.0))for j in RN(Ug)]
 mw,cdf,ac=sum(ws),[],0.0
 for w in ws:ac+=w/mw;cdf.append(ac)
 rws,Lp=[],exp(-mw)
 for _ in RN(int(round(eta*U))):
  kp,pa=0,1.0
  while pa>Lp:kp+=1;pa*=rng.random()
  d={}
  for _ in RN(max(4,kp-1)):
   c,m=min(U-1,bisect.bisect_left(cdf,rng.random())),MONO[rng.randrange(12)]
   if c in d:
    a=d[c];nv=(a[0]+m[0],a[1]+m[1],a[2]+m[2],a[3]+m[3])
    if nv==(0,0,0,0):del d[c]
    else:d[c]=nv
   else:d[c]=m
  if d:rws.append(d)
 cr={}
 for ri,d in enumerate(rws):
  for c in d:cr.setdefault(c,set()).add(ri)
 def kr(ri):
  d=rws[ri]
  if d is None:return
  for c in d:
   s=cr.get(c)
   if s:
    s.discard(ri)
    if not s:cr.pop(c,None)
  rws[ri]=None
 def psg():
  q=[c for c,s in cr.items()if len(s)==1]
  while q:
   c=q.pop();s=cr.get(c)
   if not s or len(s)!=1:continue
   ri=next(iter(s));d=rws[ri]
   if d is None:continue
   nb=list(d);kr(ri)
   for nc in nb:
    ns=cr.get(nc)
    if ns and len(ns)==1:q.append(nc)
 psg();ar=[i for i,d in enumerate(rws)if d]
 while len(ar)-len(cr)>12:kr(max(ar,key=lambda i:len(rws[i])));psg();ar=[i for i,d in enumerate(rws)if d]
 while True:
  did=False
  for dg in(2,3,4):
   cd=None
   for c,s in list(cr.items()):
    if len(s)!=dg:continue
    rl=list(s);mi=[MMAP.get(rws[r][c],-1)for r in rl]
    if any(k<0 for k in mi):continue
    od=sorted(RN(dg),key=lambda ix:len(rws[rl[ix]]));rp,kp=rl[od[0]],mi[od[0]];lp=len(rws[rp])
    if dg>=3 and any(lp+len(rws[rl[o]])>180 for o in od[1:]):continue
    cd=(rp,kp,[(rl[o],mi[o])for o in od[1:]]);break
   if cd is not None:
    rp,kp,tgs=cd;dp=rws[rp]
    for rt,kt in tgs:
     dt,sh=rws[rt],(kt-kp)%12
     for col,v0 in dp.items():
      sv=_zm(v0,sh)
      if col in dt:
       a=dt[col];nv=(a[0]-sv[0],a[1]-sv[1],a[2]-sv[2],a[3]-sv[3])
       if nv==(0,0,0,0):
        del dt[col];cr[col].discard(rt)
        if not cr[col]:cr.pop(col,None)
       else:dt[col]=nv
      else:dt[col]=(-sv[0],-sv[1],-sv[2],-sv[3]);cr.setdefault(col,set()).add(rt)
     if not dt:kr(rt)
    kr(rp);psg();did=True;break
  if not did:break
 al=[d for d in rws if d];cls=sorted({c for d in al for c in d});mx=ov=nz=c1=cu=0
 for d in al:
  for v in d.values():
   nz+=1
   for ck in v:
    av=abs(ck)
    if av==0:continue
    if av==1:c1+=1
    elif av<(1<<30):cu+=1
    else:ov+=1
    if av>mx:mx=av
 return{"rows":len(al),"cols":len(cls),"c1":c1,"cu":cu,"nz":nz,"mx":mx,"ov":ov}
def check_live_sge(full=False):
 ns,mxc,tov,tc1,tcu,tnz=(12 if full else 6),0,0,0,0,0
 for s in RN(ns):
  r=_sgep(512,REF["eta"],7500+s);assert r["rows"]>=r["cols"]and r["ov"]==0 and r["mx"]<(1<<30)
  mxc=max(mxc,r["mx"]);tov+=r["ov"];tc1+=r["c1"];tcu+=r["cu"];tnz+=r["nz"]
 pm1=tc1/max(1,tc1+tcu);cpe=(tc1+tcu)/max(1,tnz);assert pm1>=0.85
 o={"ok":True,"scope":"finite_scale_Z12_4basis_SGE_algebra_only_NOT_asymptotic_H3_exploratory_v1","seeds_u512":ns,"max_abs":mxc,"ov30":tov,"pm1_frac":round(pm1,4),"comp_per_ent":round(cpe,3)}
 if full:
  r1=_sgep(1024,REF["eta"],9401);assert r1["rows"]>=r1["cols"]and r1["ov"]==0 and r1["mx"]<(1<<30)
  o["u1024_check"]={"rows":r1["rows"],"cols":r1["cols"],"max_abs":r1["mx"],"pm1_frac":round(r1["c1"]/max(1,r1["c1"]+r1["cu"]),4)}
 return o
def _mr4(Bb):
 s,oe=REF["s1b"],REF["oe"];L1,L2=(oe*(2.0**60)*8.0)/4.5,(oe*(2.0**60)*8.0)/9.0
 if Bb<=1024.0:nl=max(1.0,Bb/64.0);return 12.0*nl*nl*s["saf"]
 if Bb<=float(1<<59):
  N=float(1<<int(ceil(log2(max(2.0,2.0*ceil(Bb/64.0))))));assert N<=(1<<54)
  return(3.0*(1.5*N*log2(N)*s["bfly"]+N*(s["pw"]+s["sc"]))+((3.0*N*160.0)if(Bb*2.0>L1)else 0.0)+N*s["c1"])*s["saf"]
 L=int(ceil(Bb/58.0));k=max(1,(L-1).bit_length());k1=(k+1)//2;k2=k-k1;N2=float(1<<(k1+1))*float(1<<(k2+1))
 return(3.0*(N2*(1.0*(k1+1)+1.5*(k2+1))*s["bfly"]+N2*(s["pw"]+s["sc"]))+((9.0*N2*160.0)if(Bb*2.0>L2)else 0.0)+N2*s["c2"])*s["saf"]
def _rr4(nlv,lb):
 if nlv<=1.0:return 40000.0
 dm,kp=REF["s1b"]["dm"],max(1,int(4096.0//lb));ns,sb,nll=max(1.0,nlv/float(kp)),float(kp)*lb,lb/64.0
 isb=sum(12.0*(j*nll)*nll*1.15 for j in RN(1,kp));isd=float(kp)*24.0*(sb/64.0)*nll*1.15
 h,ip,id_,cn,cb=int(ceil(log2(ns))),0.0,0.0,ns,sb
 for _ in RN(h):pr=cn*0.5;mi=_mr4(cb);ip+=pr*mi;id_+=pr*dm*mi;cn,cb=pr,cb*2.0
 return ns*(isb+isd)+ip+3.0*_mr4(ns*sb)+id_+nlv*max(25000.0,40000.0*(lb/636.0))
def _oecm(qb,Nb,_ec={},_ep={}):
 qk,ky=round(qb,2),(round(qb,2),int(ceil(Nb/64.0)))
 if ky in _ec:return _ec[ky]
 if qk not in _ep:
  pl,qe,rh=[],max(1.0,qb-3.5),DK.rho
  for lb2 in RN(16,125):
   l1=lb2*0.5;B1=2.0**l1;p1=rh(qe/l1)
   for m in(10.0,25.0,60.0,120.0,250.0):
    B2=B1*m;dt=(log2(B2)-l1)/40.0
    p2=sum((rh((qe-(l1+(i+0.5)*dt))/l1)*(dt/(l1+(i+0.5)*dt)))if(l1+(i+0.5)*dt)<qe else(dt/(l1+(i+0.5)*dt))for i in RN(40))
    pr=min(1.0,p1+p2)
    if pr>1e-40:pl.append((B1,B2,pr))
  _ep[qk]=pl
 L=max(4,int(ceil(Nb/64.0)));mi=212.0 if L==4 else(13.25*L*L+80.0 if L<=16 else min(13.25*L*L+80.0,2.5*_mr4(L*64.0)));bs=(1e300,None)
 for B1,B2,pr in _ep[qk]:
  c1=1.15*(40.0*mi+(B1/log(2.0))*(11.0*mi+93.0*L)+2.5*max(0.0,B2/log(B2)-B1/log(B1))*(mi+40.0*L))
  if c1/pr<bs[0]:bs=(c1/pr,(B1,B2,pr,c1))
 _ec[ky]=bs;return bs
def _dn(Q,sd,Bft,Bgt,Bfm,Bgm,af,ag):
 ld,lv,CR=2.0*log2(1125.0/16.0)+Q,log2((pi**8)/40320.0),2.0**20
 for _ in RN(8):
  lR=(log2(CR)+ld-lv)/16.0+log2(1.25);ls=lR-1.5
  nf=16.0+72.0*ls-(Q if sd=="f" else 0.0);ng=8.0*LOG2_ABS_X+24.0*ls+2.0-(Q if sd=="g" else 0.0)
  CR=max(1024.0,2.0/max(1e-40,DK.rho((nf+af)/Bft)*DK.rho((ng+ag)/Bgt)))
 cu=(CR*6000.0+(2.0**22)+3.0*(_oecm(Bft,nf)[0]+_oecm(Bgt,ng)[0]))/C_COST
 def ec(Ne,lo,hi):
  if hi<=lo+1e-9:return 0.0
  r0,M=DK.rho(Ne/hi),40;dx=(hi-lo)/M
  return sum((DK.rho((Ne-(lo+(i+0.5)*dx))/hi)/r0)*(dx/(lo+(i+0.5)*dx))for i in RN(M)if(lo+(i+0.5)*dx)<Ne)
 return{"cu":cu,"cf":ec(nf+af,Bfm,Bft),"cg":ec(ng+ag,Bgm,Bgt)}
def _s4d(af,ag):
 Bf,Bg=REF["Bf"],REF["Bg"]
 def br(Nr,Qm,lo,hi):
  r0,dx=DK.rho(Nr/Qm),(hi-lo)/200.0
  return sum((DK.rho((Nr-(lo+(i+0.5)*dx))/Qm)/r0)*(dx/(lo+(i+0.5)*dx))for i in RN(200))
 nsp,ni0,ni1,ni2=br(2572.0,275.0,Bg,275.0),0.5*br(2572.0,275.0,0.5*Bg,137.5),0.5*br(2572.0,275.0,137.5,195.0),0.5*br(2572.0,275.0,195.0,275.0)
 d5=_dn(550.0,"g",390.0,390.0,Bf,Bg,af,ag);ut=ni2*d5["cu"];nfc,ngc=ni2*d5["cf"],ni2*d5["cg"]+ni1;nl=[round(ni2,2),round(nfc+ngc,1)]
 for(Q,Bft,Bgt),ag_ in(((390.0,275.0,275.0),nsp+ni0),((275.0,180.0,180.0),0.0),((180.0,125.0,125.0),0.0),((125.0,94.0,94.0),0.0),((94.0,78.0,78.0),0.0),((78.0,Bf,Bg),0.0)):
  df,dg=_dn(Q,"f",Bft,Bgt,Bf,Bg,af,ag),_dn(Q,"g",Bft,Bgt,Bf,Bg,af,ag)
  ut+=nfc*df["cu"]+ngc*dg["cu"];nfc,ngc=nfc*df["cf"]+ngc*dg["cf"],nfc*df["cg"]+ngc*dg["cg"]+ag_;nl.append(round(nfc+ngc,1))
 assert nfc==0.0 and ngc==0.0
 ub=((1.0/(0.918*DK.rho(2572.0/275.0)))*(2.5*_oecm(275.0,2572.0)[0]))/C_COST
 return{"l0_split":round(nsp,2),"l0_inert_le275":round(ni0,2),"l0_inert_275_390":round(ni1,2),"l0_inert_390_550":round(ni2,2),"nodes":nl,"s4":log2(2.0*(ub+ut)),"ecm":log2(2.0*(ut+0.25*ub)),"desc_only":log2(2.0*ut)}
def _lll(B,dl=0.99):
 n,m,B=len(B),len(B[0]),[r[:]for r in B];bo,mu,B2=[[0.0]*m for _ in RN(n)],[[0.0]*n for _ in RN(n)],[0.0]*n
 def gs(st):
  for i in RN(st,n):
   bi=B[i];v=[float(x)for x in bi]
   for j in RN(i):
    bj=bo[j];s=sum(bi[t]*bj[t]for t in RN(m))/B2[j]if B2[j]>1e-18 else 0.0;mu[i][j]=s
    for t in RN(m):v[t]-=s*bj[t]
   bo[i]=v;B2[i]=sum(x*x for x in v)
 gs(0);k=1
 while k<n:
  for j in RN(k-1,-1,-1):
   if abs(mu[k][j])>0.5:
    q=int(round(mu[k][j]));bk,bj=B[k],B[j]
    for t in RN(m):bk[t]-=q*bj[t]
    for t in RN(j):mu[k][t]-=q*mu[j][t]
    mu[k][j]-=q
  if B2[k]>=(dl-mu[k][k-1]**2)*B2[k-1]:k+=1
  else:B[k],B[k-1]=B[k-1],B[k];gs(k-1);k=max(k-1,1)
 return B
def _h90_q11():
 xt=-56;rt=xt**4-xt**2+1;pt=((xt-1)**2*rt)//3+xt;xp=xt%pt;assert rt==9831361 and pt==10647363907 and pt%15==7 and(pt**6+1)%rt==0
 def f12m(a,b):
  rs=[0]*12
  for j1 in RN(3):
   for j2 in RN(3):
    j=j1+j2;mx=xp if j>=3 else 1;jo=j-3 if j>=3 else j;tp_=[0]*7
    for i1 in RN(4):
     c1=a[4*j1+i1]
     if c1:
      for i2 in RN(4):tp_[i1+i2]+=c1*b[4*j2+i2]
    tp_[0]+=tp_[5];tp_[1]+=tp_[6];c4=tp_[4]
    for i in RN(4):rs[4*jo+i]=(rs[4*jo+i]+(tp_[i]-c4)*mx)%pt
  return rs
 def f12p(a,e):
  rs,bs=[1]+[0]*11,list(a)
  while e>0:
   if e&1:rs=f12m(rs,bs)
   bs=f12m(bs,bs);e>>=1
  return rs
 def fr6(a):
  o=[0]*12
  for j in RN(3):c0,c1,c2,c3=a[4*j:4*j+4];o[4*j:4*j+4]=[(c0-c1)%pt,(-c1)%pt,(c3-c1)%pt,(c2-c1)%pt]
  return o
 g0=f12p([3,5,7,11,13,17,19,23,29,31,37,41],(pt**6-1)*((pt**6+1)//rt));assert f12p(g0,rt)==[1]+[0]*11 and g0!=[1]+[0]*11
 tt=f12p(g0,4192837);assert f12m(tt,fr6(tt))==[1]+[0]*11;A=[[0]*12 for _ in RN(12)]
 for k in RN(12):
  ek=[0]*12;ek[k]=1;im=[(ek[i]-f12m(tt,fr6(ek))[i])%pt for i in RN(12)]
  for rr in RN(12):A[rr][k]=im[rr]
 pv,ri=[],0
 for cl in RN(12):
  p0=next((rr for rr in RN(ri,12)if A[rr][cl]%pt!=0),None)
  if p0 is None:continue
  A[ri],A[p0]=A[p0],A[ri];iv=pow(A[ri][cl],pt-2,pt);A[ri]=[(v*iv)%pt for v in A[ri]]
  for rr in RN(12):
   if rr!=ri and A[rr][cl]:f=A[rr][cl];A[rr]=[(A[rr][cc]-f*A[ri][cc])%pt for cc in RN(12)]
  pv.append(cl);ri+=1
 assert len(pv)==6;fc=[c for c in RN(12)if c not in pv]
 B=[[(pt if j==p_ else 0)for j in RN(12)]for p_ in pv]+[[(1 if j==f_ else((-A[pv.index(j)][f_])%pt if j in pv else 0))for j in RN(12)]for f_ in fc]
 ws=[1]*4+[4]*4+[15]*4;Br=_lll([[r[i]*ws[i]for i in RN(12)]for r in B]);al=[Br[0][i]//ws[i]for i in RN(12)];am=[c%pt for c in al]
 assert am==f12m(tt,fr6(am))and any(am)and f12m(f12p(am,pt**6-1),tt)==[1]+[0]*11
 xc,w3,nv=(abs(xt)**(1.0/3.0))*cmath.exp(1j*pi/3.0),cmath.exp(2j*pi/3.0),1.0+0j
 for m in(1,2,3,4):
  z5=cmath.exp(2j*pi*m/5.0)
  for k3 in RN(3):yv=xc*(w3**k3);nv*=sum(al[4*j+i]*(z5**i)*(yv**j)for j in RN(3)for i in RN(4))
 l2n,l2d=log2(abs(nv)),6.0*log2(pt)+4.0*log2(abs(xt));assert l2n<l2d+2.0
 q=1031;assert q%15==11;z5r=next(t for t in RN(2,q)if(t**4+t**3+t**2+t+1)%q==0);yr=pow(xt%q,pow(3,-1,q-1),q);assert(yr**3-xt)%q==0
 B16,p5=[[(q if j==ix else 0)for j in RN(16)]for ix in(0,4)],[pow(z5r,i,q)for i in RN(4)]
 for i in RN(1,4):r0,r1=[0]*16,[0]*16;r0[i],r0[0]=1,(-p5[i])%q;r1[4+i],r1[4]=1,(-p5[i])%q;B16.extend((r0,r1))
 for i in RN(4):r0,r1=[0]*16,[0]*16;r0[8+i],r0[0]=1,(yr*p5[i])%q;r1[12+i],r1[4]=1,(yr*p5[i])%q;B16.extend((r0,r1))
 v0=_lll(B16)[0];assert sum((v0[i]-yr*v0[8+i])*p5[i]for i in RN(4))%q==0 and sum((v0[4+i]-yr*v0[12+i])*p5[i]for i in RN(4))%q==0
 def okm(u,v):
  def z5m(a,b):
   t=[0]*7
   for i in RN(4):
    for j in RN(4):t[i+j]+=a[i]*b[j]
   t[0]+=t[5];t[1]+=t[6];c4=t[4];return[t[i]-c4 for i in RN(4)]
  p00,p01,p10,p11=z5m(u[:4],v[:4]),z5m(u[:4],v[4:]),z5m(u[4:],v[:4]),z5m(u[4:],v[4:])
  return[p00[i]-p11[i]for i in RN(4)]+[p01[i]+p10[i]-p11[i]for i in RN(4)]
 a3,b3=okm(okm(v0[:8],v0[:8]),v0[:8]),okm(okm(v0[8:],v0[8:]),v0[8:]);gok=[a3[i]-xt*b3[i]for i in RN(8)]
 A8,pv_,sg=[okm(gok,[(1 if j==k else 0)for j in RN(8)])for k in RN(8)],1,1
 for k in RN(7):
  if A8[k][k]==0:sw=next(i for i in RN(k+1,8)if A8[i][k]!=0);A8[k],A8[sw]=A8[sw],A8[k];sg=-sg
  pv2=A8[k][k]
  for i in RN(k+1,8):
   for j in RN(k+1,8):A8[i][j]=(A8[i][j]*pv2-A8[i][k]*A8[k][j])//pv_
  pv_=pv2
 nlg=sg*A8[7][7];assert nlg!=0 and nlg%(q*q)==0
 return{"p":pt,"r":rt,"rank_H90":6,"lift_verified":True,"log2_N_Mg":round(l2n,2),"log2_det_w":round(l2d,2),"q11_inert_descent":{"q":q,"q_mod_15":11,"exact_q2_divides_N_Lg":True,"log2_N_Lg_cof":round(log2(abs(nlg)/(q*q)),2)}}
def check_stage4(af,ag,full=False):
 s4=_s4d(af,ag);assert s4["s4"]<REF["cr"]["both_individual_logs_exp"]and s4["ecm"]<REF["cr"]["ecm_split_expected"]
 h90=_h90_q11();lll=33.0-6.0*log2(12.0)+4.5*log2(5.0)+6.0*log2(3.0)+6.0*log2(P_BLS)+4.0*log2(ABS_X);assert lll<2572.0
 o={"ok":True,"lll_ub":round(lll,4),"term":{"Bf":REF["Bf"],"Bg":REF["Bg"],"out_f":0,"out_g":0},"l0_split":s4["l0_split"],"l0_inert_le275":s4["l0_inert_le275"],"l0_inert_275_390":s4["l0_inert_275_390"],"l0_inert_390_550":s4["l0_inert_390_550"],"nodes":s4["nodes"],"desc_only_l2":round(s4["desc_only"],2),"s4_l2":round(s4["s4"],4),"ecm_l2":round(s4["ecm"],4),"real_bls12_h90_q11":h90}
 if not full:return o
 xd=lambda X,Z,A,N:(((X+Z)**2*(X-Z)**2)%N,((((X+Z)**2-(X-Z)**2)%N)*(((X-Z)**2+A*((X+Z)**2-(X-Z)**2))%N))%N)
 xa=lambda X1,Z1,X2,Z2,X,Z,N:((Z*(((X1-Z1)*(X2+Z2)+(X1+Z1)*(X2-Z2))**2))%N,(X*(((X1-Z1)*(X2+Z2)-(X1+Z1)*(X2-Z2))**2))%N)
 def xm(k,X,Z,A,N):
  r0x,r0z,r1x,r1z=1,0,X%N,Z%N
  for bt in bin(k)[2:]:
   if bt=="0":r1x,r1z=xa(r0x,r0z,r1x,r1z,X,Z,N);r0x,r0z=xd(r0x,r0z,A,N)
   else:r0x,r0z=xa(r0x,r0z,r1x,r1z,X,Z,N);r1x,r1z=xd(r1x,r1z,A,N)
  return r0x,r0z
 t16=[p for p in _primes(200)if p>=101]
 for p in t16:
  x=(5*pow(-7,-1,p))%p;d=((2*x*x-1)*pow(x,-4,p))%p;A=(2*(1+d)*pow(1-d,-1,p))%p
  assert(1+sum(1 if(u**3+A*u*u+u)%p==0 else(2 if pow((u**3+A*u*u+u)%p,(p-1)//2,p)==1 else 0)for u in RN(p)))%16==0
 N1=50023*50033;x1=(5*pow(-7,-1,N1))%N1;A24=pow(1-((2*x1*x1-1)*pow(x1,-4,N1))%N1,-1,N1);X,Z=3,1
 for p in _primes(30):
  pk=p
  while pk*p<=30:pk*=p
  X,Z=xm(pk,X,Z,A24,N1)
 assert gcd(Z,N1)==1;bb=[xm(i,X,Z,A24,N1)for i in RN(1,16)if gcd(i,30)==1];DX,DZ=xm(30,X,Z,A24,N1);gp,gc,w1=(1,0),(DX,DZ),1
 for k in RN(1,13):
  for bx,bz in bb:
   g=gcd((gc[0]*bz-bx*gc[1])%N1,N1)
   if 1<g<N1:w1=g;break
  if w1>1:break
  gp,gc=gc,(xd(DX,DZ,A24,N1)if k==1 else xa(gc[0],gc[1],DX,DZ,gp[0],gp[1],N1))
 X2,Z2=xm(2,4,1,2,35);DX2,DZ2=xm(4,X2,Z2,2,35);w2=gcd((DX2*Z2-X2*DZ2)%35,35)
 assert w1==50023 and gcd(Z2,35)==1 and w2==5;o.update({"t16":f"{len(t16)}/{len(t16)}","w1":w1,"w2":[w2,35//w2]});return o
def _tx64(l):
 if l<=0.0:return 0.0
 if l>45.0:return(69.0/64.0)*l+24.0
 ei,pk=0.0,exp(-l)
 for k in RN(1,160):
  pk*=l/k;q,r=divmod(k,64);ei+=pk*(10.0+q*69.0+(35.0 if r&32 else 0)+(19.0 if r&16 else 0)+(11.0 if r&8 else 0)+(7.0 if r&4 else 0)+(3.0 if r&2 else 0)+(3.0 if r&1 else 0))
  if k>l+25.0 and pk<1e-14:break
 return ei
T16=(5.0,10.0,9.0,13.0,15.0,20.0,19.0,23.0,23.0,28.0,27.0,31.0,33.0,38.0,37.0,41.0)
def _tx16(l):
 if l<=0.0:return 8.0
 if l>45.0:return 11.0+(37.0/16.0)*l+22.0
 p0=exp(-l);e,pk=p0*8.0,p0
 for k in RN(1,160):
  pk*=l/k;q,r=divmod(k,16);e+=pk*(11.0+37.0*q+T16[r])
  if k>l+24.0 and pk<1e-13:break
 return e
def _ihyb(plo,pw1,pw2,phi,L248,L124,dok,rmx,w2=16,w3=32,pp=20.0):
 w1,w2_=1.0/(V_A_SH*L248),1.0/(V_A_SH*L124);ih=ip=nok=nb=ne=0.0
 if pw1>plo:
  dx=(pw1-plo)/100.0
  for i in RN(100):
   x=plo+(i+0.5)*dx;p=2.0**x;d=(p/x)*dx;dk=d*dok;nok+=dk;nb+=dk*(p+2.0);ne+=d*V_A_SH
   ip+=dk*(pp+L248*7.0)*w1;ih+=dk*sum(pr*L248*(_tx64(r*V_A_SH/p)+1.5625*(r*V_A_SH/p))for r,pr in rmx)*w1
 for a,b,W in((pw1,pw2,w2),(pw2,phi,w3)):
  if b<=a:continue
  nb_,dx=(L124-1.0)/float(W)+1.0,(b-a)/60.0
  for i in RN(60):
   x=a+(i+0.5)*dx;p=2.0**x;d=(p/x)*dx;dk=d*dok;nok+=dk;nb+=dk*(p/float(W)+2.0);ne+=d*V_A_SH;eip=eih=0.0
   for r,pr in rmx:lb=float(W)*r*V_A_SH/p;pc=1.0-exp(-nb_*lb);eip+=pr*(pp+pc*6.0);eih+=pr*(nb_*_tx16(lb)-(1.0-pc)*nb_*8.0+nb_*1.5625*lb)
   ip+=dk*eip*w2_;ih+=dk*eih*w2_
 return{"ih":ih,"ip":ip,"nok":nok,"nb":nb,"ne":ne}
def _tp(cof,Pt,Bb,km=8):
 if cof<=Pt:return 1.0,0.0,1.0
 om=lambda v:0.0 if v<1.0 else(1.0/v if v<=2.0 else((1.0+log(v-1.0))/v if v<=3.0 else 0.57))
 dy,md,ma=(cof-Pt)/200.0,DK.rho(cof/Pt),0.0;mt=md
 for i in RN(200):
  y=Pt+(i+0.5)*dy;w=DK.rho((cof-y)/Pt)*(om(y/Pt)/Pt)*dy;mt+=w
  if y<=Bb:md+=w
  elif any(k*Pt<=y<=k*Bb for k in RN(2,km+1)):ma+=w
 return md/max(mt,1e-30),ma/max(mt,1e-30),(md+ma)/max(mt,1e-30)
def check_ledger(lm,lis):
 Bf,Bg,eta,b,rep,la=REF["Bf"],REF["Bg"],REF["eta"],REF["b"],REF["rep"],REF["la"]
 af,ag,cp=lm["af_500k"],lm["ag_500k"],lm["cof_pad"]
 L4,Lh,gp,lsg,lsf,ls1,ret=lis["L248_ch"],lis["L124_ch"],lis["gap_pad"],lis["l2_sg"],lis["l2_sf"],lis["l2_s1"],lis["ret_live"]
 dp=lis["discretization_convergence"]["disc_pad_b"];w0=1.0/(1.0+REF["rbr"]**8);w1=1.0-w0
 fmx=[(3,39.0/46.0),(6,6.0/46.0),(9,1.0/46.0)];dg,df,dng,dnf=1.0/3.0,23.0/81.0,1.0/24.0,2131.0/34992.0
 sg0,sf0=_ihyb(b["pglo"],b["pw1"],b["pw2"],b["pghi0"],L4,Lh,dg,[(3,1.0)],REF["w2"],REF["w3"]),_ihyb(b["pflo"],b["pw1"],b["pw2"],b["pfhi0"],L4,Lh,df,fmx,REF["w2"],REF["w3"])
 sg1,sf1=_ihyb(b["pglo"],b["pw1"],b["pw2"],b["pghi1"],L4,Lh,dg,[(3,1.0)],REF["w2"],REF["w3"]),_ihyb(b["pflo"],b["pw1"],b["pw2"],b["pfhi1"],L4,Lh,df,fmx,REF["w2"],REF["w3"])
 ih,ip=w0*(sg0["ih"]+sf0["ih"])+w1*(sg1["ih"]+sf1["ih"]),w0*(sg0["ip"]+sf0["ip"])+w1*(sg1["ip"]+sf1["ip"])
 def rsh(sg,sf):ne,nb=sg["ne"]+sf["ne"],sg["nb"]+sf["nb"];return(ne*rep["ent"]+nb*rep["bkt"])/C_COST,log2(8.0*(ne+nb)+48.0*(sg["nok"]+sf["nok"]))
 r0,m0=rsh(sg0,sf0);r1,m1=rsh(sg1,sf1);ru=r0+r1
 lVk=lis["Vk_live"]+lis["pad_3se_Vk"]+gp+cp+dp;Vt,LtM,lf,lg=2.0**lVk,0.0,REF["lf"],REF["lg"]
 def ps1b(pfm,pgm):
  nonlocal LtM
  s1=2.0**ls1;c1,c2,mg=pfm+(Bf-pfm)/3.0,pfm+2.0*(Bf-pfm)/3.0,0.5*(pgm+Bg)
  _,af1,_=_tp(lf,pfm,Bf);_,af2,_=_tp(lf,c1,Bf);_,af3,_=_tp(lf,c2,Bf);_,ag1,_=_tp(lg,pgm,Bg);_,ag2,_=_tp(lg,mg,Bg)
  rf=min(0.95,sqrt(ret/0.80));rg=ret/rf
  pfs=min(1.0,(DK.rho((lis["mean_nf"]+af)/Bf)*rf)/(2.0**lsf));pgs=min(1.0,(DK.rho((lis["mean_ng"]+ag)/Bg)*rg)/(2.0**lsg))
  cr=s1*1700.0+s1*pfs*1700.0
  for lo,hi,pin,lb,dn in[(1.0,pfm,s1,656.0,dnf),(pfm,c1,s1*af1,lf,dnf),(c1,c2,s1*af2,lf,dnf),(c2,Bf,s1*af3,lf,dnf),(1.0,pgm,s1*pfs,688.0,dng),(pgm,mg,s1*pfs*ag1,lg,dng),(mg,Bg,s1*pfs*ag2,lg,dng)]:
   Lt=max(1.0,(ppb_dus(lo,min(hi,36.0))if lo<36.0 else 0.0)+(dn*ppb_dus(max(lo,36.0),hi)if hi>36.0 else 0.0))
   LtM=max(LtM,Lt);Mt=Vt*pin;Mb=min(Mt,max(1.0,Lt/lb));nbb,LX=max(1.0,Mt/Mb),Mb*lb;ir=_rr4(Mb,lb)
   if LX<Lt:ir+=max(0.0,Lt/LX-1.0)*2.5*_mr4(LX)
   cr+=(nbb*ir)/Vt
  return cr+s1*pfs*pgs*(110.0*256.0*150.0*4.0+0.05*5.0e7*2.0)
 cs1b=w0*ps1b(b["pfhi0"],b["pghi0"])+w1*ps1b(b["pfhi1"],b["pghi1"])
 sll=52.0/256.0+2.5*(19.0/256.0)+64.0/128.0+200.0*(2.0**lsg);slc=sll+0.004;ipcl,ipcc=ih+ip+sll+cs1b,ih+ip+slc+cs1b
 cl_,cc_=(Vt*ipcl)/C_COST,(Vt*ipcc)/C_COST;U=U_ub(Bf,Bg);Nl,Nc=U*la["nl"],U*la["comp"]
 bsp=la["pm1_f"]*10.0+(1.0-la["pm1_f"])*25.0+(Nl*36.0)/Nc;assert bsp<la["x4"]and 948.0<la["rv"]
 Mb,Nb=la["m"],la["n"];Lsq,dmk=Nl/Mb+Nl/Nb+64.0,Nl/Nb+32.0
 ibw=1.02*(Nb*Lsq*(Nc*la["x4"]+Nl*la["rv"])*la["ov"]+dmk*(Nc*la["x1"]+Nl*la["rv"])*la["ov"])
 imk=1.02*dmk*Nb*Nl*283.0*la["ov"];ilg=20.0*((Mb+Nb)**3)*Lsq*log2(Lsq)*236.0*la["ov"]
 ixf=4.0*Mb*(Mb+Nb)*Lsq*ceil(log2(Lsq))*8.0*la["ov"]
 lau=(ibw+imk+ilg+ixf+0.008*(ibw+imk)+30.0*U*283.0*la["ov"])/C_COST
 sgeu=((12500.0*U+25.0*U*log2(U)+100.0*Nc)*1.15)/C_COST;s4=_s4d(af,ag)
 lr={"relation_collection":log2(cl_),"replay_a_tables":log2(ru),"linear_algebra":log2(lau),"both_individual_logs_exp":s4["s4"],"ecm_split_expected":s4["ecm"],"meter_safety_buffer":93.4,"band_split_and_trees":90.0,"canonical_a_domain":84.0,"schirokauer_maps":80.0,"final_log_backsub_and_gcd":72.0,"matrix_generator_sge":log2(sgeu),"stage0_setup_and_iso":45.0,"final_verification_Q":40.0}
 cr=REF["cr"];assert len(lr)==len(cr)==13 and all(cr[k]>lr[k]for k in lr)and cr["relation_collection"]>=log2(cc_)
 wol,woc=sum(2.0**v for v in lr.values()),sum(2.0**v for v in cr.values());lol,loc=log2(wol),log2(woc)
 pb=ceil((2.0**m0)/(12.0*(2.0**60)));gsB=248.0*V_A_SH*(2.0**lsg)*8.0
 hse=(186.0*(2.0**56))/(2.0**60);s1p=log2((2.0**m0)+248.0*(V_A_SH+262144.0)+2.0*186.0*(2.0**56)+gsB)
 ntB,cb,nn,lv=0.0,4080.0,LtM/4080.0,0
 while nn>=1.0:
  if cb<=(1<<59):ntB+=nn*(float(1<<int(ceil(log2(max(2.0,2.0*ceil(cb/64.0))))))*24.0)
  else:L2=int(ceil(cb/58.0));k=max(1,(L2-1).bit_length());k1=(k+1)//2;ntB+=nn*(float(1<<(k1+1))*float(1<<(k-k1+1))*24.0)
  cb*=2.0;nn*=0.5;lv+=1
 s1bB=ntB+(9.0+lv)*LtM/8.0;s1be=s1bB/(2048.0*(2.0**60));s2r=log2((eta*U)*18.0*12.0)
 lmB=Nc*(la["pm1_f"]*8.0+(1.0-la["pm1_f"])*12.0)+Nl*64.0;lvB=Nl*32.0;lae=(lmB/la["p"]+lvB*(1.0+(Nb+6.0)/la["p"]))/(2.0**60)
 dpi=Lsq*Nb/(Mb+Nb);lgB=(Mb*Nb*Lsq+Mb*(Mb+Nb)*Lsq+4.0*((Mb+Nb)**2)*dpi)*32.0;lge=lgB/(la["plg"]*(2.0**60));s3r=log2(lmB+(Nb+6.0+la["p"])*lvB+lgB)
 assert V_A_SH/(2.0**60)<0.06 and hse<12.0<16.0 and s1be<REF["oe"]<16.0 and lae<6.0<16.0 and lge<8.0<16.0 and max(s1p,log2(s1bB),s2r,s3r)<REF["ram_ub"]<120.0
 pm=woc/(2.0**124.0);pf=pm+2.0**(-48)+2.0**(-100)+2.0**(-132)+2.0**(-200);po=1.0-pf;Kg,t64=REF["Kg"],1<<64
 assert pf<2.0e-5 and po>0.99998 and Fraction(Kg,t64)*Fraction(99998,100000)>=Fraction(39,100)and Fraction(Kg-1,t64)*Fraction(99998,100000)<Fraction(39,100)
 pc=(Kg/float(t64))*po;num,den=REF["num"],REF["den"];lgc,lgl=loc+log2(Kg/float(t64)),lol+log2(Kg/float(t64))
 slk=float(Fraction(num,den))-lgc;assert pc>0.39 and gcd(num,den)==1 and slk>2.0*lis["se_Vk"]
 return{"ok":True,"claim":f"{num}/{den}","eta":eta,"U_l2":round(log2(U),6),"U_ub":float(f"{U:.8e}"),"Vk_ch_live":round(lVk,6),"disc_pad_b":dp,"ih":round(ih,6),"ip":round(ip,6),"c_s1b":round(cs1b,6),"slab_live":round(sll,6),"slab_ch":round(slc,6),"ipc_live":round(ipcl,6),"ipc_ch":round(ipcc,6),"blended_spmv_comp":round(bsp,4),"la_mat_eib":round(lmB/(2.0**60),4),"la_owner_eib":round(lae,4),"P_lingen":la["plg"],"lingen_owner_eib":round(lge,4),"half_slab_eib":round(hse,4),"opened_live":round(lol,6),"opened_ch":round(loc,6),"gated_live":round(lgl,6),"gated_ch":round(lgc,6),"slack_b":round(slk,6),"slack_over_live_SE_x":round(slk/lis["se_Vk"],2),"p_meter_cond_H":float(f"{pm:.7e}"),"p_open_cond_H":round(po,8),"p_cold_cond_H":round(pc,8),"ledger_live":{k:round(v,6)for k,v in lr.items()},"ram_l2_B":{"S1_sh0":round(s1p,4),"S1_sh1":round(m1,4),"P_banks":pb,"S1b":round(log2(s1bB),4),"S1b_eib":round(s1be,4),"S2":round(s2r,4),"S3":round(s3r,4),"S3_lingen":round(log2(lgB),4)}}
def main():
 t0=time.perf_counter();raw=sys.stdin.read().strip();inp=json.loads(raw)if raw else{}
 if not isinstance(inp,dict)or(set(inp)-{"part"}):raise ValueError("Input must be a JSON object with optional key 'part'")
 part=inp.get("part","core")
 if part not in("core","diagnostics","all"):raise ValueError(f"Unsupported part: {part}")
 out={"package_id":REF["pkg"],"config_tag":REF["tag"],"part":part,"reads_external_files":0,"offline_mc_constants":0}
 if part in("core","all"):
  out["stage0"]=check_stage0()
  out["chebotarev"]=check_chebotarev(26000 if part=="all" else 15000)
  lm=check_murphy_norms(1500000 if part=="all" else 500000);out["murphy_norms"]=lm
  lis=check_live_panjer_rqmc_chords_disc(lm["af_500k"],lm["ag_500k"],32 if part=="all" else 24,2048 if part=="all" else 1024,8000 if part=="all" else 1200,full_disc=(part=="all"))
  out["live_panjer_rqmc_chords"]=lis;out["ntt_r4_shoup61"]=check_ntt();out["kernels_e2e"]=check_kernels_e2e()
  out["genuine_tower_dlp_sm"]=check_genuine_tower_dlp_sm();out["live_z12_sge_pipe"]=check_live_sge(full=(part=="all"))
  out["stage4_derived"]=check_stage4(lm["af_500k"],lm["ag_500k"],full=(part=="all"));out["ledger_ram_budget"]=check_ledger(lm,lis)
 elif part=="diagnostics":
  lm=check_murphy_norms(1500000);out["murphy_norms"]=lm
  out["live_panjer_rqmc_chords"]=check_live_panjer_rqmc_chords_disc(lm["af_500k"],lm["ag_500k"],32,2048,8000,full_disc=True)
  out["genuine_tower_dlp_sm"]=check_genuine_tower_dlp_sm();out["live_z12_sge_pipe"]=check_live_sge(full=True)
  out["stage4_derived"]=check_stage4(lm["af_500k"],lm["ag_500k"],full=True)
 out["checks_ok"]=all(v.get("ok",True)for v in out.values()if isinstance(v,dict))
 out["elapsed_seconds"]=round(time.perf_counter()-t0,3)
 print(json.dumps(out,indent=2))
 if not out["checks_ok"]:sys.exit(1)
if __name__=="__main__":main()
