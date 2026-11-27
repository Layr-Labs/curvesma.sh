#!/usr/bin/env python3
"""Standalone Docker-Safe Verifier for stnfs-q15-hybrid-split-eta1015-v7 (stdlib only, 0 external files, 0 offline MC constants)."""
from __future__ import annotations
import bisect,cmath,json,math,random,sys,time
from fractions import Fraction
X_BLS= -0xD201000000010000
ABS_X,LOG2_ABS_X=abs(X_BLS),math.log2(abs(X_BLS))
P_BLS=int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",16)
R_BLS=int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",16)
P0_X=int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",16)
P0_Y=int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1",16)
Q0_X=(int("024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8",16),int("13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e",16))
Q0_Y=(int("0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801",16),int("0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be",16))
C_COST,V_A_SH,SLAB_BASE=20000.0,0.9351*(2.0**56),1<<56
REF:dict={
"pkg":"stnfs-q15-hybrid-split-eta1015-v7","tag":"CAND7_HYBRID_LIVE1M_ETA1015",
"Bf":71.45,"Bg":71.81,"eta":1.015,"cush": -0.0080,
"qa":0.65,"rbr":0.82,"Ra":1719.4681182270274,"R0":1355.3079808905318,"sl":1.75,
"lam":340599.43367449316,"lf":408.0,"lg":688.0,"oe":10.0,"w2":16,"w3":32,
"b":{"pglo":28.0,"pw1":55.0,"pw2":58.5,"pghi0":62.5,"pghi1":60.7,"pflo":32.0,"pfhi0":62.1,"pfhi1":60.3,"qg":3.4,"qf":3.8},
"s1b":{"p1":(127<<54)+1,"g1":3,"p2":(114<<54)+1,"g2":7,"p3":(108<<54)+1,"g3":5,"bfly":11.125,"pw":12.0,"sc":12.0,"c1":58.0,"c2":90.0,"saf":1.10,"dm":1.67},
"rep":{"ent":21.50,"bkt":5.5},
"la":{"m":512,"n":32,"p":4,"plg":512,"nl":0.00080,"comp":0.045,"pm1_f":0.92,"x4":11.90,"x1":14.90,"rv":960.0,"ov":1.06},
"sge":[(16,0.002544,0.000743,0.14680,0.854),(18,0.001166,0.000125,0.04930,0.882),(20,0.000449,0.000044,0.00551,0.970),(22,0.000705,0.000212,0.02547,0.958),(23,0.000529,0.000098,0.00966,0.970),(24,0.000595,0.000210,0.01405,0.967)],
"cr":{
"relation_collection":107.9480,"replay_a_tables":104.4575,"linear_algebra":103.1815,
"both_individual_logs_exp":96.6500,"ecm_split_expected":96.0500,
"meter_safety_buffer":93.4001,"band_split_and_trees":90.0001,"canonical_a_domain":84.0001,
"schirokauer_maps":80.0001,"final_log_backsub_and_gcd":72.0001,"matrix_generator_sge":63.1490,
"stage0_setup_and_iso":45.0001,"final_verification_Q":40.0001,
},
"Kg":7194374076228249696,"num":4271,"den":40,"ram_ub":118.0,
}
class ExactDickman:
 def __init__(self,umax:float=30.0,N:int=200)->None:
  self.umax,self.N=umax,N;h,tot=1.0/N,int(round(umax*N))+1
  v=[1.0]*(N+1)+[0.0]*(tot-N-1)
  for i in range(N+1,min(2*N+1,tot)):v[i]=1.0-math.log(i*h)
  for i in range(2*N+1,tot):v[i]=(h*(0.5*v[i-N]+sum(v[i-N+1:i])))/(i*h-0.5*h)
  self.v=v
 def rho(self,u:float)->float:
  if u<=1.0:return 1.0
  if u>=self.umax:return 1e-300
  pos=u*self.N;idx=int(pos)
  if idx+1>=len(self.v):return self.v[-1]
  f=pos-idx;return math.exp((1.0-f)*math.log(self.v[idx])+f*math.log(self.v[idx+1]))
DK=ExactDickman()
def _primes(n:int)->list[int]:
 ip=bytearray(b"\x01")*(n+1);ip[0]=ip[1]=0
 for i in range(2,int(n**0.5)+1):
  if ip[i]:ip[i*i:n+1:i]=bytearray(len(range(i*i,n+1,i)))
 return[i for i in range(2,n+1)if ip[i]]
def pi_dus(B:float)->float:
 lnB=B*math.log(2.0);return((2.0**B)/lnB)*(1.0+1.0/lnB+2.334/(lnB*lnB))
def U_ub(Bf:float,Bg:float)->float:return(pi_dus(Bf)+pi_dus(Bg))/12.0+36.0*pi_dus(0.5*Bf)+12.0*pi_dus(0.5*Bg)
def ppb_dus(lo:float,hi:float)->float:
 th_hi=(2.0**hi)*(1.0+0.01/max(1.0,hi*math.log(2.0)))
 th_lo=0.0 if lo<=1.0 else(2.0**lo)*(1.0-0.01/max(1.0,lo*math.log(2.0)))
 return max(1.0,(th_hi-th_lo)/math.log(2.0))
def _g1a(p1,p2):
 if p1 is None:return p2
 if p2 is None:return p1
 x1,y1,x2,y2=p1[0],p1[1],p2[0],p2[1]
 if x1==x2:
  if(y1+y2)%P_BLS==0:return None
  lam=(3*x1*x1*pow(2*y1,-1,P_BLS))%P_BLS
 else:lam=((y2-y1)*pow(x2-x1,-1,P_BLS))%P_BLS
 x3=(lam*lam-x1-x2)%P_BLS;return(x3,(lam*(x1-x3)-y1)%P_BLS)
def check_stage0()->dict:
 assert(X_BLS-1)**2*(X_BLS**4-X_BLS**2+1)//3+X_BLS==P_BLS and X_BLS**4-X_BLS**2+1==R_BLS
 lam=X_BLS%R_BLS
 assert P_BLS%R_BLS==lam and(pow(lam,4,R_BLS)-pow(lam,2,R_BLS)+1)%R_BLS==0 and pow(lam,6,R_BLS)==R_BLS-1 and pow(lam,12,R_BLS)==1 and pow(lam,4,R_BLS)!=1
 r_pt,a_pt,k=None,(P0_X,P0_Y),R_BLS
 while k:
  if k&1:r_pt=_g1a(r_pt,a_pt)
  a_pt=_g1a(a_pt,a_pt);k>>=1
 assert(P0_Y*P0_Y-(P0_X**3+4))%P_BLS==0 and r_pt is None
 m2=lambda a,b:((a[0]*b[0]-a[1]*b[1])%P_BLS,(a[0]*b[1]+a[1]*b[0])%P_BLS)
 q3=m2(m2(Q0_X,Q0_X),Q0_X);assert m2(Q0_Y,Q0_Y)==((q3[0]+4)%P_BLS,(q3[1]+4)%P_BLS)
 assert P_BLS%15==7 and[j for j in range(1,5)if pow(P_BLS,j,5)==1][0]==4 and pow(X_BLS%P_BLS,(P_BLS-1)//3,P_BLS)!=1
 assert any((X_BLS**3-X_BLS**2+(c*c%P_BLS)*X_BLS-c)%P_BLS==0 for g in range(2,20)for cand in(pow(g,(P_BLS-1)//3,P_BLS),)if cand!=1 for c in(cand,(cand*cand)%P_BLS))
 assert X_BLS**6-2*X_BLS**5+2*X_BLS**3+X_BLS+1==3*P_BLS and(R_BLS-1)%12==0 and ABS_X==(1<<16)*906349*254760293
 return{"ok":True,"p_mod_15":7,"ord_5_p":4,"n_sm_chi1":4,"abs_x":"2^16*906349*254760293"}
def check_chebotarev(p_lim:int=26000)->dict:
 zwm=lambda x,y:(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1])
 zwn=lambda x:x[0]*x[0]-x[0]*x[1]+x[1]*x[1]
 zwdiv=lambda d,n:all(c%zwn(d)==0 for c in zwm(n,(d[0]-d[1],-d[1])))
 dw,dw2,p19,p67=(41,24),(17,-24),(-5,-2),(-9,-2)
 assert zwm(p19,p67)==dw and zwn(dw)==1273==19*67 and math.gcd(1273,15)==1 and not zwdiv(p19,dw2)and not zwdiv(p67,dw2)
 assert(118*118+118+1)%151==0 and all((t**3-t**2+32*t-118)%151!=0 for t in range(151))
 pw,ww=991,113;assert pw%45==1 and(ww*ww+ww+1)%pw==0
 wr=[t for t in range(pw)if(t**3-t**2+(ww*ww%pw)*t-ww)%pw==0];assert wr==[373,802,808]
 kch=pow((wr[0]*pow(wr[1],pw-2,pw))%pw,(pw-1)//3,pw);assert kch==877 and pow(kch,3,pw)==1
 s3=[(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)]
 G=[((a,b,c),s)for a in range(3)for b in range(3)for c in range(3)for s in s3]
 gmul=lambda a,b:(tuple((a[0][b[1][i]]+b[0][i])%3 for i in range(3)),tuple(a[1][b[1][i]]for i in range(3)))
 def ginv(g):
  si=[0,0,0]
  for i in range(3):si[g[1][i]]=i
  st=tuple(si);return(tuple((-g[0][st[i]])%3 for i in range(3)),st)
 seen,ccl=set(),[]
 for g in G:
  if g not in seen:cc=frozenset(gmul(gmul(h,g),ginv(h))for h in G);seen|=cc;ccl.append(cc)
 assert len(G)==162 and len(ccl)==22
 ide,hz9=((0,0,0),(0,1,2)),frozenset(g for g in G if sum(g[0])%3==0)
 nsubs,q={frozenset([ide])},[frozenset([ide])]
 while q:
  N=q.pop()
  for cc in ccl:
   if cc.issubset(N):continue
   M=set(N)|set(cc)
   while True:
    add={gmul(a,b)for a in M for b in cc}-M
    if not add:break
    M|=add
   Mf=frozenset(M)
   if Mf not in nsubs:nsubs.add(Mf);q.append(Mf)
 assert sorted(len(N)for N in nsubs)==[1,3,9,27,27,54,81,162]
 sh=[N for N in nsubs if N<hz9];mx=[N for N in sh if not any(N<M for M in sh)];assert len(mx)==1 and len(mx[0])==27
 rc,c0a,cnza,r0v,rnzv={0:0,3:0,6:0,9:0},0,0,[],[]
 for v,s in G:
  fx=sum(1 for i in range(3)for c in range(3)if s[i]==i and(c+v[i])%3==c);rc[fx]+=1
  if sum(v)%3==0:r0v.append(fx);c0a+=int(fx>0)
  else:rnzv.append(fx);cnza+=int(fx>0)
 p0,p3,p6,p9=(Fraction(rc[k],162)for k in(0,3,6,9))
 dok,d0,dnz=1-p0,Fraction(c0a,54),Fraction(cnza,108)
 df=Fraction(1,8)*(Fraction(1,3)*(1-(1-d0)**2)+Fraction(2,3)*(1-(1-dnz)**2))
 assert(p0,p3,p6,p9)==(Fraction(58,81),Fraction(13,54),Fraction(1,27),Fraction(1,162))and dok==Fraction(23,81)and df==Fraction(2131,34992)
 m0,m1=sum(r0v)/len(r0v),sum(rnzv)/len(rnzv)
 v0,v1=sum((x-m0)**2 for x in r0v)/len(r0v),sum((x-m1)**2 for x in rnzv)/len(rnzv)
 assert abs(((1.0/3.0)*16.0*(2.0*v0+(2.0*m0-2.0)**2)+(2.0/3.0)*16.0*(2.0*v1+(2.0*m1-2.0)**2))/8.0-12.0)<1e-9
 ps=_primes(p_lim);p15=[p for p in ps if p%15==1];assert len(p15)>300;bm=fp=ca=0
 for p in p15:
  w=next(pow(g,(p-1)//3,p)for g in range(2,p)if pow(g,(p-1)//3,p)!=1);ph,e3=False,(p-1)//3
  for wv in(w,(w*w)%p):
   w2,ny=(wv*wv)%p,0
   for t in range(1,p):
    if(t*(t*(t-1)+w2)-wv)%p==0 and pow(t,e3,p)==1:ny+=3
   if ny%3:bm+=1
   if ny:ca+=1;ph=True
  if ph:fp+=1
 od,odn=ca/(2.0*len(p15)),fp/float(len(ps));assert bm==0 and abs(od-float(dok))<0.04 and abs(odn-float(df))<0.01
 return{"ok":True,"delta_Q":"41+24*w","kummer_p":pw,"galois_order":162,"d_ok_f":str(dok),"dens_f":str(df),"var_f_root":12.0,"var_g_root":16.0,"split_p":len(p15),"live_d_ok":round(od,5),"live_dens":round(odn,5)}
def _prc(poly,qe,z,o,ad,sb,ml,iv):
 def trf(a):
  while len(a)>1 and a[-1]==z:a.pop()
  return a
 def pm(a,b):
  r=[z]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   if x!=z:
    for j,y in enumerate(b):
     if y!=z:r[i+j]=ad(r[i+j],ml(x,y))
  return trf(r)
 def md(a,m):
  a,dm,il=list(a),len(m)-1,iv(m[-1])
  while len(a)-1>=dm:
   c,sh=ml(a[-1],il),len(a)-1-dm
   if c!=z:
    for i in range(dm+1):a[sh+i]=sb(a[sh+i],ml(c,m[i]))
   a.pop()
  return trf(a)
 res,b,e=[o],md([z,o],poly),qe
 while e:
  if e&1:res=md(pm(res,b),poly)
  b=md(pm(b,b),poly);e>>=1
 if len(res)<2:res+=[z]*(2-len(res))
 res[1]=sb(res[1],o);trf(res)
 if res==[z]:return len(poly)-1
 a,c=list(poly),res
 while any(x!=z for x in c):
  a,c=c,md(a,c)
  if c==[z]:break
 return len(a)-1
def _fast_murphy(Pmax:int=1000000):
 ps=_primes(Pmax);af=ag=0.0;ck={}
 def trim(a):
  while len(a)>1 and a[-1]==0:a.pop()
  return a
 def pmod(a,b,p):
  a=trim(list(a));iv=pow(b[-1],p-2,p);db=len(b)-1
  while len(a)-1>=db and not(len(a)==1 and a[0]==0):
   c=(a[-1]*iv)%p;sh=len(a)-1-db
   for k in range(db+1):a[sh+k]=(a[sh+k]-c*b[k])%p
   a.pop()
   if not a:return[0]
   trim(a)
  return trim(a)
 def pgcd(u,v,p):
  u,v=trim(list(u)),trim(list(v))
  while not(len(v)==1 and v[0]==0):u,v=v,pmod(u,v,p)
  iv=pow(u[-1],p-2,p);return[(x*iv)%p for x in u]
 def sroots(F,p):
  deg=len(F)-1
  if deg==0:return[]
  if deg==1:return[((-F[0])*pow(F[1],p-2,p))%p]
  def pm(x,y):
   res=[0]*(len(x)+len(y)-1)
   for i,vx in enumerate(x):
    if vx:
     for j,vy in enumerate(y):res[i+j]=(res[i+j]+vx*vy)%p
   return pmod(res,F,p)
  a_s=1
  while True:
   res,bs,e2=[1],[a_s%p,1],(p-1)//2
   while e2:
    if e2&1:res=pm(res,bs)
    bs=pm(bs,bs);e2>>=1
   res[0]=(res[0]-1)%p;g1=pgcd(F,res,p)
   if 1<=len(g1)-1<deg:
    F2=list(F);d1=len(g1)-1;q_=[0]*(deg-d1+1)
    for k in range(deg-d1,-1,-1):
     c=F2[k+d1];q_[k]=c
     for j in range(d1+1):F2[k+j]=(F2[k+j]-c*g1[j])%p
    return sroots(g1,p)+sroots(trim(q_),p)
   a_s+=1
 for p in ps:
  er,lp=1.0/(p-1.0),math.log2(p)
  if p<=300:
   if p==2:ef,eg=0.0,8.0*(1.0/17.0)*(sum(3.0*k*(15.0/(16.0**k))for k in range(1,6))+16.0/(16.0**5))
   elif p==3:ef,eg=4.0*(3.0/82.0),4.0*(1.75/82.0)
   elif p==5:ef,eg=0.0,2.0*3.0*25.0/(25.0**2-1.0)
   else:
    cur,fd=p%15,1
    while cur!=1:cur=(cur*p)%15;fd+=1
    gc,qs=8//fd,p**fd;pf,xm=float(qs),X_BLS%p
    rg=1 if xm==0 else((3 if pow(xm,(qs-1)//3,p)==1 else 0)if(qs-1)%3==0 else 1)
    eg=fd*gc*rg*pf/(pf*pf-1.0)
    ops=(0,1,lambda a,b:(a+b)%p,lambda a,b:(a-b)%p,lambda a,b:(a*b)%p,lambda a:pow(a,p-2,p))
    if fd==1:
     w1=next(pow(g,(p-1)//3,p)for g in range(2,p)if pow(g,(p-1)//3,p)!=1)
     ef=sum(4*_prc([(-w)%p,0,0,(w*w)%p,0,0,(-1)%p,0,0,1],p,*ops)for w in(w1,(w1*w1)%p))*pf/(pf*pf-1.0)
    elif fd==2 or(fd==4 and p<=13):
     if p%3==1:
      w1=next(pow(g,(p-1)//3,p)for g in range(2,p)if pow(g,(p-1)//3,p)!=1)
      rf=sum((4//fd)*_prc([(-w)%p,0,0,(w*w)%p,0,0,(-1)%p,0,0,1],qs,*ops)for w in(w1,(w1*w1)%p))
     else:
      f2i=lambda a:(((a[0]-a[1])*pow((a[0]*a[0]-a[0]*a[1]+a[1]*a[1])%p,p-2,p))%p,((-a[1])*pow((a[0]*a[0]-a[0]*a[1]+a[1]*a[1])%p,p-2,p))%p)
      pf2=[(0,(-1)%p),(0,0),(0,0),((-1)%p,(-1)%p),(0,0),(0,0),((-1)%p,0),(0,0),(0,0),(1,0)]
      rf=gc*_prc(pf2,qs,(0,0),(1,0),lambda a,b:((a[0]+b[0])%p,(a[1]+b[1])%p),lambda a,b:((a[0]-b[0])%p,(a[1]-b[1])%p),lambda a,b:((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0]-a[1]*b[1])%p),f2i)
     ef=fd*rf*pf/(pf*pf-1.0)
    else:ef=8.0/pf
  else:
   r15=p%15
   if r15==1:
    pf,xm,e3=float(p),X_BLS%p,(p-1)//3
    rg=1 if xm==0 else(3 if pow(xm,e3,p)==1 else 0)
    eg=8.0*rg*pf/(pf*pf-1.0);w1,g=1,2
    while True:
     c=pow(g,e3,p)
     if c!=1:w1=c;break
     g+=1
    rf_tot=0
    for w in(w1,(w1*w1)%p):
     w2=(w*w)%p;r0,r1,r2,b0,b1,b2,ev=1,0,0,0,1,0,p
     while ev:
      if ev&1:
       c0,c1,c2,c3,c4=(r0*b0)%p,(r0*b1+r1*b0)%p,(r0*b2+r1*b1+r2*b0)%p,(r1*b2+r2*b1)%p,(r2*b2)%p
       c3,c2,c1=(c3+c4)%p,(c2-c4*w2)%p,(c1+c4*w)%p
       r0,r1,r2=(c0+c3*w)%p,(c1-c3*w2)%p,(c2+c3)%p
      c0,c1,c2,c3,c4=(b0*b0)%p,(2*b0*b1)%p,(2*b0*b2+b1*b1)%p,(2*b1*b2)%p,(b2*b2)%p
      c3,c2,c1=(c3+c4)%p,(c2-c4*w2)%p,(c1+c4*w)%p
      b0,b1,b2=(c0+c3*w)%p,(c1-c3*w2)%p,(c2+c3)%p;ev>>=1
     r1=(r1-1)%p;Q=[(-w)%p,w2,p-1,1]
     if r0==0 and r1==0 and r2==0:
      for rt in sroots(Q,p):
       if pow(rt,e3,p)==1:rf_tot+=3
     else:
      g_=pgcd(Q,[r0,r1,r2],p)
      if len(g_)==2 and pow(((-g_[0])*pow(g_[1],p-2,p))%p,e3,p)==1:rf_tot+=3
    ef=4.0*rf_tot*pf/(pf*pf-1.0)
   else:
    fd=2 if(r15 in(4,11,14))else 4;qs,gc,xm=float(p**fd),8//fd,X_BLS%p
    rg=3 if(p%3==2 or pow(xm,(p-1)//3,p)==1)else 0
    eg=fd*gc*rg/qs;ef=8.0/qs
  af+=(er-ef)*lp;ag+=(er-eg)*lp
  if p<=1300000:af_13,ag_13=af,ag
  if p in(4999,49999,499979,999983,1999993):ck[str(p)]=[round(af,5),round(ag,5)]
 ck["1300000"]=[round(af_13,5),round(ag_13,5)]
 return af,ag,ck
def _cof_ex(Bf:float,pfhi:float,mnf:float,lf:float,af:float)->tuple[float,float]:
 M,dx,h=80,(Bf-pfhi)/80.0,2.0;nb=int(math.ceil(660.0/h))+1;cur=[0.0]*nb;cur[0]=1.0;tp=ep=0.0;ne=mnf+af
 for k in range(1,11):
  nxt=[0.0]*nb
  for m in range(M):
   x=pfhi+(m+0.5)*dx;w=dx/(x*k);sh=int(round(x/h))
   for b in range(nb-sh):
    if cur[b]>0.0:nxt[b+sh]+=cur[b]*w
  cur=nxt
  for b in range(nb):
   sb=b*h
   if sb>ne:break
   cv=cur[b]*DK.rho((ne-sb)/pfhi);tp+=cv
   if sb>lf:ep+=cv
 tp+=DK.rho(ne/pfhi);return ep/tp,tp/DK.rho(ne/Bf)
def check_murphy_norms(pmax:int=1300000)->dict:
 af,ag,ck=_fast_murphy(pmax);b=REF["b"];Ra,R0,q=REF["Ra"],REF["R0"],REF["qa"];Bf,Bg=REF["Bf"],REF["Bg"]
 af_u,ag_u=ck["1300000"][0],ck["1300000"][1]
 sf=abs(math.log2(DK.rho((542.17+af_u+0.5)/Bf)/DK.rho((542.17+af_u-0.5)/Bf)))
 sg=abs(math.log2(DK.rho((670.13+ag_u+0.5)/Bg)/DK.rho((670.13+ag_u-0.5)/Bg)))
 It=(math.log(1300000.0)+1.0)/(1300000.0*(math.log(2.0)**2))
 sd_r=math.sqrt((sf*sf*12.0+sg*sg*16.0)*It);sd_u=math.sqrt((sf*sf*19.0+sg*sg*23.0+14.0*sf*sg)*It);assert 3.0*sd_u<abs(REF["cush"])
 dnf=(72.0/q)*math.log2((Ra**q+R0**q)/4.0);dng=8.0*LOG2_ABS_X+(24.0/q)*math.log2((R0**q+(Ra**q)/(ABS_X**(q/3.0)))/4.0)
 assert dnf<656.0 and dng<688.0 and 6*656<=4096 and 5*688<=4096 and 10*408<=4096
 tfm=int(math.floor(656.0/b["qf"]))-127;assert tfm==45 and(128-tfm)+int(math.floor(656.0/b["qf"]))<=255
 for ngc in range(610,689):assert(128-max(1,int(math.floor((ngc+1.0)/b["qg"]))-127))+int(math.floor(ngc/b["qg"]))<=255
 w0=1.0/(1.0+REF["rbr"]**8);p0,bn0=_cof_ex(Bf,b["pfhi0"],542.17,REF["lf"],af_u);p1,bn1=_cof_ex(Bf,b["pfhi1"],542.17,REF["lf"],af_u)
 assert abs(bn0-1.0)<0.01 and abs(bn1-1.0)<0.01;pex=w0*p0+(1.0-w0)*p1;cpad=math.log2(1.0/(1.0-pex))
 return{"ok":True,"pmax":pmax,"af_1300k":af_u,"ag_1300k":ag_u,"af_end":round(af,5),"ag_end":round(ag,5),"ck":ck,"tail_1300k_sd_root":round(sd_r,6),"tail_1300k_sd_cov_ub":round(sd_u,6),"tail_1300k_3sd_ub":round(3.0*sd_u,6),"cush_b":REF["cush"],"det_nf":round(dnf,6),"det_ng":round(dng,6),"p_ex_f":round(pex,6),"cof_pad":round(cpad,6),"bnorm":[round(bn0,5),round(bn1,5)]}
VH=((0.481211825060,0.481211825060,-0.481211825060,-0.481211825060),(0.687654016088,-0.206442191029,0.396282847190,-0.877494672250),(0.206442191029,-0.687654016088,0.877494672250,-0.396282847190),(0.877494672250,-0.396282847190,-0.687654016088,0.206442191029),(0.396282847190,-0.877494672250,-0.206442191029,0.687654016088),(1.083936863278,-1.083936863278,0.189840656161,-0.189840656161),(0.189840656161,-0.189840656161,-1.083936863278,1.083936863278))
VHH=tuple(0.5*sum(c*c for c in v)for v in VH)
PERMS=((0,1,2,3),(0,1,3,2),(0,2,1,3),(0,2,3,1),(0,3,1,2),(0,3,2,1),(1,0,2,3),(1,0,3,2),(1,2,0,3),(1,2,3,0),(1,3,0,2),(1,3,2,0),(2,0,1,3),(2,0,3,1),(2,1,0,3),(2,1,3,0),(2,3,0,1),(2,3,1,0),(3,0,1,2),(3,0,2,1),(3,1,0,2),(3,1,2,0),(3,2,0,1),(3,2,1,0))
def _panjer_dp(Bb:float,plo:float,phi:float,qs:float,ab:float,bmin:int,bstep:int,nb:int):
 KM,TM=720,260;mo=[0.0]*(KM+1)
 for a_i,b_i in((0.561459483567,plo),(phi,Bb)):
  if b_i<=a_i:continue
  ns=max(20,int((b_i-a_i)*8.0));dx=(b_i-a_i)/ns
  for i in range(ns):
   x0=a_i+i*dx;kb=int(round(x0+0.5*dx))
   if 1<=kb<=KM:mo[kb]+=math.log((x0+dx)/x0)
 nzo=[(j,j*mo[j])for j in range(1,KM+1)if mo[j]>0.0];Po=[0.0]*(KM+1);Po[0]=math.exp(-sum(mo[1:]))
 for k in range(1,KM+1):Po[k]=sum(jm*Po[k-j]for j,jm in nzo if j<=k)/k
 mi,m1={},[0.0]*(TM+1);ns=max(80,int((phi-plo)*16.0));dx=(phi-plo)/ns
 for i in range(ns):
  x0=plo+i*dx;xm=x0+0.5*dx;w,kb,tk=math.log((x0+dx)/x0),int(round(xm)),int(math.floor(xm/qs))
  if kb>=1 and 1<=tk<=TM:mi[(kb,tk)]=mi.get((kb,tk),0.0)+w;m1[tk]+=w
 li=sum(m1[1:]);nz1=[(j,j*m1[j])for j in range(1,TM+1)if m1[j]>0.0];P1=[0.0]*(TM+1);P1[0]=math.exp(-li)
 for k in range(1,TM+1):P1[k]=sum(jm*P1[k-j]for j,jm in nz1 if j<=k)/k
 sv,run=[0.0]*(TM+1),0.0
 for t in range(TM,-1,-1):run+=P1[t];sv[t]=run
 its=sorted([(dk,dt,dk*v)for(dk,dt),v in mi.items()]);P2=[{}for _ in range(KM+1)];P2[0][0]=math.exp(-li)
 for k in range(its[0][0],KM+1):
  ik,dst=1.0/k,{}
  for dk,dt,dv in its:
   if dk>k:break
   pr=P2[k-dk]
   if not pr:continue
   sd=dv*ik
   for tt,p0 in pr.items():
    if tt+dt<=TM:dst[tt+dt]=dst.get(tt+dt,0.0)+sd*p0
  if dst:P2[k]={tt:pv for tt,pv in dst.items()if pv>1e-13}
 nzk=[(k,list(P2[k].items()))for k in range(KM+1)if P2[k]];rt=[]
 for bi in range(nb):
  kN=max(1,min(KM,int(round(bmin+(bi+0.5)*bstep+ab))));sm=[0.0]*(TM+1)
  for k,ki in nzk:
   if k>kN:break
   po=Po[kN-k]
   if po<=1e-300:continue
   for tt,pt in ki:sm[tt]+=pt*po
  tot,tl,row=sum(sm),0.0,[0.0]*129
  for tt in range(TM,-1,-1):
   tl+=sm[tt]
   if tt<=128:row[tt]=(tl/tot)if tot>0.0 else 0.0
  rt.append(row)
 return sv,rt
def check_live_panjer_is_chords(af:float,ag:float,nfit:int=32000,neval:int=256000,nch:int=8000)->dict:
 b,Bf,Bg,qa,rbr,Ra,R0=REF["b"],REF["Bf"],REF["Bg"],REF["qa"],REF["rbr"],REF["Ra"],REF["R0"]
 w0=1.0/(1.0+rbr**8);w1=1.0-w0;G0,GS,GB,F0,FS,FB,N0,NS,NB=610,2,40,510,4,25,440,2,100
 sgL,rgL,sfL,rfL=[],[],[],[]
 for pgh,pfh in((b["pghi0"],b["pfhi0"]),(b["pghi1"],b["pfhi1"])):
  sg,rg=_panjer_dp(Bg,b["pglo"],pgh,b["qg"],ag,G0,GS,GB);sf,rf=_panjer_dp(Bf,b["pflo"],pfh,b["qf"],af,N0,NS,NB)
  sgL.append(sg);rgL.append(rg);sfL.append(sf);rfL.append(rf)
 alpha,alpha_p,a_tb,a_ta0=2.0/qa,1.55,4.6,5.3
 C_is=math.exp(math.lgamma(4.0*alpha)-4.0*math.lgamma(alpha)-(math.lgamma(4.0*alpha_p)-4.0*math.lgamma(alpha_p)))
 d_alp,iqa,agd,il2,tp=alpha-alpha_p,1.0/qa,2.0/qa,1.0/math.log(2.0),2.0*math.pi
 Ap=[(a**6,a**12,a**18)for a in((Ra/(4.0**iqa))*((0.25 if k==0 else 0.75)**0.125)for k in range(2))]
 DK_rho,mlog,mcos,msin=DK.rho,math.log,math.cos,math.sin
 def run_is(n_a:int,rng:random.Random,Tg=None,Tf=None,v2=None,s3=None,nfit_s:float=1.0):
  gvar,urnd=rng.gammavariate,rng.random;swk=swr=ssg=ss1=swc=mnf=mng=swn=nfm=ngm=0.0;vt=vo=0
  for _ in range(n_a):
   while True:
    vt+=1;g0,g1,g2,g3=gvar(agd,1.0),gvar(agd,1.0),gvar(agd,1.0),gvar(agd,1.0);ta=g0+g1+g2+g3
    u0,u1,u2,u3=iqa*mlog(g0/ta),iqa*mlog(g1/ta),iqa*mlog(g2/ta),iqa*mlog(g3/ta);um=0.25*(u0+u1+u2+u3)
    u0-=um;u1-=um;u2-=um;u3-=um
    if all(abs(u0*VH[k][0]+u1*VH[k][1]+u2*VH[k][2]+u3*VH[k][3])<=VHH[k]for k in range(7)):vo+=1;break
   a3_0,a3_1,a3_2,a3_3=(g0/ta)**(3.0*iqa),(g1/ta)**(3.0*iqa),(g2/ta)**(3.0*iqa),(g3/ta)**(3.0*iqa)
   pm=PERMS[int(urnd()*24)];vk=vr=vg=v1=vc=0.0
   for jb in range(4):
    utb=(jb+urnd())*0.25;tb_rad=utb**(1.0/a_tb);w_tb=(8.0/a_tb)*(tb_rad**(8.0-a_tb))
    va=(pm[jb]+urnd())*0.25;ta0=va**(1.0/a_ta0);w_ta0=(8.0/a_ta0)*(ta0**(8.0-a_ta0))
    ra3_0,ra3_1=(Ra*(0.5**0.125)*ta0)**3,(Ra*((0.5+0.5*va)**0.125))**3
    gb0,gb1,gb2,gb3=gvar(alpha_p,1.0),gvar(alpha_p,1.0),gvar(alpha_p,1.0),gvar(alpha_p,1.0)
    tb=gb0+gb1+gb2+gb3;ub0,ub1,ub2,ub3=gb0/tb,gb1/tb,gb2/tb,gb3/tb
    w_gb=C_is*((ub0*ub1*ub2*ub3)**d_alp);vc+=w_tb*w_gb*(w0*w_ta0+w1)
    b1_0,b1_1,b1_2,b1_3=ub0**iqa,ub1**iqa,ub2**iqa,ub3**iqa;b3_0,b3_1,b3_2,b3_3=b1_0**3,b1_1**3,b1_2**3,b1_3**3
    ang0,ang1,ang2,ang3=tp*urnd(),tp*urnd(),tp*urnd(),tp*urnd()
    c0,s0,c1,s1,c2,s2,c3,s3_=mcos(ang0),msin(ang0),mcos(ang1),msin(ang1),mcos(ang2),msin(ang2),mcos(ang3),msin(ang3)
    for sh,sh_w,ra3,rb1 in((0,w0*w_ta0,ra3_0,R0*tb_rad),(1,w1,ra3_1,R0*rbr*tb_rad)):
     A6,A12,A18=Ap[sh];rb3=rb1**3;lnf=lng=lng0=nfp=0.0
     for a3m,b1m,b3m,cm,sm in((a3_0,b1_0,b3_0,c0,s0),(a3_1,b1_1,b3_1,c1,s1),(a3_2,b1_2,b3_2,c2,s2),(a3_3,b1_3,b3_3,c3,s3_)):
      bm=rb1*b1m;r6=bm**6;r12=r6*r6;nfp+=mlog(A18+A12*r6+A6*r12+r12*r6)*il2
      ar,ai,bm3=ra3*a3m*cm,ra3*a3m*sm,rb3*b3m
      gr=ar-X_BLS*bm3;lng+=mlog(gr*gr+ai*ai+1e-300);g0z= -X_BLS*bm3;lng0+=mlog(g0z*g0z+1e-300)
      ar2,ai2=ar*ar-ai*ai,2.0*ar*ai;ar3,ai3=ar2*ar-ai2*ai,ar2*ai+ai2*ar;b6=bm3*bm3;b9=b6*bm3
      fr=ar3-ar2*bm3+(-0.5*ar+0.8660254037844386*ai)*b6+0.5*b9
      fi=ai3-ai2*bm3+(-0.5*ai-0.8660254037844386*ar)*b6-0.8660254037844386*b9
      lnf+=mlog(fr*fr+fi*fi+1e-300)
     nf,ng,ng0=lnf*il2,lng*il2,lng0*il2
     if nf>nfm:nfm=nf
     if ng>ngm:ngm=ng
     ig=max(0,min(GB-1,int(math.floor((ng0-G0)/GS))));icf=max(0,min(FB-1,int(math.floor((nfp-F0)/FS))));iff=max(0,min(NB-1,int(math.floor((nf-N0)/NS))))
     wt_is=sh_w*w_tb*w_gb;wt=wt_is*DK_rho((nf+af)/Bf)*DK_rho((ng+ag)/Bg)
     if Tg is None:
      v2[sh][ig][icf]+=wt_is/nfit_s;s3[sh][ig][icf][iff]+=wt/nfit_s;swr+=wt/nfit_s
     else:
      ig_a=max(0,min(GB-1,int(math.floor((ng-G0)/GS))));tg,tf=Tg[sh][ig][icf],Tf[sh][ig][icf]
      wk=wt*rgL[sh][ig_a][tg]*rfL[sh][iff][tf];sp=sgL[sh][tg]
      vr+=wt;vk+=wk;vg+=wt_is*sp;v1+=wt_is*sp*sfL[sh][tf];mnf+=wk*nf;mng+=wk*ng;swn+=wk
   if Tg is not None:swk+=vk*0.25;swr+=vr*0.25;ssg+=vg*0.25;ss1+=v1*0.25;swc+=vc*0.25
  return swk/n_a,swr/n_a,ssg/n_a,ss1/n_a,swc/n_a,vt,vo,mnf,mng,swn,nfm,ngm
 v2=[[[0.0]*FB for _ in range(GB)]for _ in range(2)];s3=[[[[0.0]*NB for _ in range(FB)]for _ in range(GB)]for _ in range(2)]
 _,tsf_n,*_=run_is(nfit//4,random.Random(20261127),None,None,v2,s3,float(nfit));tsf=tsf_n*(nfit//4)
 for sh in range(2):
  for ig in range(GB):
   for icf in range(FB):
    for iff in range(NB):s3[sh][ig][icf][iff]/=tsf
 tfm,lp=max(0,int(math.floor(656.0/b["qf"]))-127),REF["lam"]
 Tg,Tf=[[[128]*FB for _ in range(GB)]for _ in range(2)],[[[125]*FB for _ in range(GB)]for _ in range(2)]
 for sh in range(2):
  rgs,rfs,sgs,sfs=rgL[sh],rfL[sh],sgL[sh],sfL[sh]
  for ig in range(GB):
   ngu=G0+(ig+1)*GS+1.0;tgm=max(0,min(128,int(math.floor(ngu/b["qg"]))-127));vw1,rf1=sum(v2[sh][ig]),[0.0]*126
   if vw1>0.0:
    for tf in range(tfm,126):rf1[tf]=sum(s3[sh][ig][icf][iff]*rfs[iff][tf]for icf in range(FB)for iff in range(NB)if s3[sh][ig][icf][iff]>0.0)
   for icf in range(FB):
    vw=v2[sh][ig][icf]
    if vw<5e-5:
     if vw1<=0.0:continue
     vwu,rfw=vw1,rf1
    else:
     vwu,rfw=vw,[0.0]*126;nzf=[(iff,s3[sh][ig][icf][iff])for iff in range(NB)if s3[sh][ig][icf][iff]>0.0]
     for tf in range(tfm,126):rfw[tf]=sum(sv*rfs[iff][tf]for iff,sv in nzf)
    bo,btg,btf= -1e300,128,125
    for tg in range(tgm,129):
     rgv,sgt=rgs[ig][tg],sgs[tg]
     for tf in range(tfm,126):
      o=rgv*rfw[tf]-lp*vwu*sgt*sfs[tf]
      if o>bo:bo,btg,btf=o,tg,tf
    Tg[sh][ig][icf],Tf[sh][ig][icf]=btg,btf
    assert(128-btg)+int(math.floor(ngu/b["qg"]))<=255 and(128-btf)+int(math.floor(656.0/b["qf"]))<=255
 N_B,rng_e=32,random.Random(20261128);N_A=neval//(4*N_B);b_wk,b_wr,b_sg,b_s1,b_wc=[],[],[],[],[]
 vt=vo=0;mnf_s=mng_s=sw_n=nfm=ngm=0.0
 for _ in range(N_B):
  wk_,wr_,sg_,s1_,wc_,vti,voi,mnfi,mngi,swni,nfmi,ngmi=run_is(N_A,rng_e,Tg,Tf)
  b_wk.append(wk_);b_wr.append(wr_);b_sg.append(sg_);b_s1.append(s1_);b_wc.append(wc_)
  vt+=vti;vo+=voi;mnf_s+=mnfi;mng_s+=mngi;sw_n+=swni;nfm=max(nfm,nfmi);ngm=max(ngm,ngmi)
 mwk,mwr,msg,ms1=sum(b_wk)/N_B,sum(b_wr)/N_B,sum(b_sg)/N_B,sum(b_s1)/N_B
 se_wk=math.sqrt(sum((x-mwk)**2 for x in b_wk)/(N_B*(N_B-1)));se_b=(se_wk/mwk)*il2;se_Vk=REF["sl"]*se_b
 vol_unit=(math.pi**4)*(math.gamma(1.0+2.0/qa)**4)/math.gamma(1.0+8.0/qa)
 Vk_base=math.log2(V_A_SH*(vol_unit*(R0**8)/(70.3125*3.0))*(1.0+rbr**8))
 U=U_ub(Bf,Bg);Rreq=math.log2(REF["eta"]*U);Rret=Vk_base+math.log2(mwk*0.97843406*(2.0**REF["cush"]))
 Vk_live=Vk_base+REF["sl"]*(Rreq-Rret);retL,l2sg,l2s1=mwk/mwr,math.log2(msg),math.log2(ms1);l2sf=l2s1-l2sg
 ew_is=sum(b_wc)/N_B;assert abs(ew_is-1.0)<0.015 and nfm<656.0 and ngm<688.0 and se_Vk<0.0045
 rc_,hq,NC_B=random.Random(20260505),0.5*qa,nch//N_B;b_L4,b_Lh,b_g4=[],[],[]
 for _ in range(N_B):
  L4b=Lhb=g4b=0.0
  for R,wsh in((R0,w0),(R0*rbr,w1)):
   Rq,iL,ssp,sgp,w248,w124=R**qa,0.0,0.0,0.0,0.0,0.0;ns,tl=max(30,int(round(NC_B*wsh))),int(math.ceil(R*1.15))+20
   for _ in range(ns):
    gs=[rc_.gammavariate(agd,1.0)for _ in range(4)];gsm=sum(gs);rq=Rq*(max(rc_.random(),1e-16)**(qa/8.0))
    zr,zi=[],[]
    for m in range(4):rm,th=(rq*gs[m]/gsm)**iqa,tp*rc_.random();zr.append(rm*mcos(th));zi.append(rm*msin(th))
    ch=[tc for tc in range(-tl,tl+1,16)if sum(((zr[m]+tc)**2+zi[m]**2)**hq for m in range(4))<=Rq]
    tlo,thi=max(-tl,(ch[0]-32)if ch else-256),min(tl,(ch[-1]+32)if ch else 256);tmn,tmx,ic=10000,-10000,0
    for t in range(tlo,thi+1):
     if sum(((zr[m]+t)**2+zi[m]**2)**hq for m in range(4))<=Rq:
      if t<tmn:tmn=t
      if t>tmx:tmx=t
      ic+=1
    if ic:sp_=tmx-tmn+1;w=1.0/ic;iL+=w;ssp+=w*sp_;sgp+=w*(sp_-ic);w248+=w*((sp_+247)//248);w124+=w*((sp_+123)//124)
   L4b+=wsh*(ssp/w248);Lhb+=wsh*(ssp/w124);g4b+=wsh*(sgp/ssp)
  b_L4.append(L4b);b_Lh.append(Lhb);b_g4.append(g4b)
 L4m=sum(b_L4)/N_B;se_L4=math.sqrt(sum((x-L4m)**2 for x in b_L4)/(N_B*(N_B-1)))
 Lhm=sum(b_Lh)/N_B;se_Lh=math.sqrt(sum((x-Lhm)**2 for x in b_Lh)/(N_B*(N_B-1)))
 g4m=sum(b_g4)/N_B;se_g4=math.sqrt(sum((x-g4m)**2 for x in b_g4)/(N_B*(N_B-1)))
 L4_ch,Lh_ch,gap_ch=L4m-3.0*se_L4,Lhm-3.0*se_Lh,g4m+3.0*se_g4;gap_pad=math.log2(1.0/(1.0-gap_ch));assert 95.0<L4_ch<120.0 and 68.0<Lh_ch<88.0 and 0.0<gap_ch<0.008
 return{"ok":True,"n_fit":nfit,"n_eval_is":neval,"n_chords":nch,"E_w_is":round(ew_is,5),"fvor_live":round(vo/float(vt),6),"ret_live":round(retL,6),"l2_sg":round(l2sg,6),"l2_sf":round(l2sf,6),"l2_s1":round(l2s1,6),"mean_nf":round(mnf_s/sw_n,3),"mean_ng":round(mng_s/sw_n,3),"obs_max":[round(nfm,2),round(ngm,2)],"Vk_base":round(Vk_base,6),"Vk_live":round(Vk_live,6),"se_Vk":round(se_Vk,6),"pad_3se_Vk":round(3.0*se_Vk,6),"L248_live":round(L4m,3),"se_L248":round(se_L4,3),"L248_ch":round(L4_ch,3),"L124_live":round(Lhm,3),"se_L124":round(se_Lh,3),"L124_ch":round(Lh_ch,3),"gap4_live":round(g4m,6),"se_gap4":round(se_g4,6),"gap_ch":round(gap_ch,6),"gap_pad":round(gap_pad,6)}
def _ip64(n:int)->bool:
 d,s=n-1,0
 while d%2==0:d//=2;s+=1
 for a in(2,3,5,7,11,13,17,19,23,29,31,37):
  x=pow(a,d,n)
  if x in(1,n-1):continue
  for _ in range(s-1):
   x=(x*x)%n
   if x==n-1:break
  else:return False
 return True
def _shoup(a:int,w:int,wb:int,p:int)->int:return(a*w-((a*wb)>>64)*p)&0xFFFFFFFFFFFFFFFF
def _ntt_r4(vec:list[int],p:int,g:int)->list[int]:
 N,a,p2,p4,l4,t=len(vec),[x%p for x in vec],p<<1,p<<2,0,len(vec)
 while t>1:t>>=2;l4+=1
 for i in range(N):
  rv,x=0,i
  for _ in range(l4):rv=(rv<<2)|(x&3);x>>=2
  if i<rv:a[i],a[rv]=a[rv],a[i]
 J,Jp,m=pow(g,(p-1)//4,p),(pow(g,(p-1)//4,p)<<64)//p,1
 while m<N:
  wm,w1,W=pow(g,(p-1)//(4*m),p),1,[]
  for _ in range(m):w2=(w1*w1)%p;w3=(w2*w1)%p;W.append((w1,(w1<<64)//p,w2,(w2<<64)//p,w3,(w3<<64)//p));w1=(w1*wm)%p
  for k in range(0,N,4*m):
   for j in range(m):
    w1,w1p,w2,w2p,w3,w3p=W[j];x0=a[k+j]-(p2 if a[k+j]>=p2 else 0)
    x1,x2,x3=_shoup(a[k+m+j],w1,w1p,p),_shoup(a[k+2*m+j],w2,w2p,p),_shoup(a[k+3*m+j],w3,w3p,p)
    t0,t1=(x0+x2)-(p2 if x0+x2>=p2 else 0),(x0+p2-x2)-(p2 if x0+p2-x2>=p2 else 0)
    t2,t3=(x1+x3)-(p2 if x1+x3>=p2 else 0),_shoup(x1+p2-x3,J,Jp,p)
    a[k+j],a[k+2*m+j],a[k+m+j],a[k+3*m+j]=t0+t2,t0+p2-t2,t1+t3,t1+p2-t3
    assert a[k+j]<p4 and a[k+m+j]<p4 and a[k+2*m+j]<p4 and a[k+3*m+j]<p4
  m<<=2
 return[x%p for x in a]
def _ntt1(vec:list[int],p:int,g:int,inv:bool=False)->list[int]:
 n,a,p2,j=len(vec),[x%p for x in vec],2*p,0
 for i in range(1,n):
  b=n>>1
  while j&b:j^=b;b>>=1
  j^=b
  if i<j:a[i],a[j]=a[j],a[i]
 L=2
 while L<=n:
  hf,wl=L>>1,pow(g,(p-1)//L,p)
  if inv:wl=pow(wl,p-2,p)
  ws,wp,cur=[1]*hf,[0]*hf,1
  for k in range(hf):ws[k],wp[k]=cur,(cur<<64)//p;cur=(cur*wl)%p
  for i in range(0,n,L):
   for k in range(hf):u=a[i+k]-(p2 if a[i+k]>=p2 else 0);t=_shoup(a[i+k+hf],ws[k],wp[k],p);a[i+k],a[i+k+hf]=u+t,u+p2-t
  L<<=1
 return[(x*pow(n,p-2,p))%p for x in a]if inv else[x%p for x in a]
def _mul2d(A:int,B:int,pgs:list[tuple[int,int]])->int:
 m58=(1<<58)-1
 def t58(x):
  r=[]
  while x:r.append(x&m58);x>>=58
  return r or[0]
 la,lb=t58(A),t58(B);k=max(1,(max(len(la),len(lb))-1).bit_length());k1=(k+1)//2;D1,D2=1<<k1,1<<(k-k1);N1,N2=2*D1,2*D2
 la+=[0]*(D1*D2-len(la));lb+=[0]*(D1*D2-len(lb));cm=[]
 for p,g in pgs:
  gA=[_ntt1(la[r*D1:(r+1)*D1]+[0]*D1,p,g)for r in range(D2)]+[[0]*N1 for _ in range(D2)]
  gB=[_ntt1(lb[r*D1:(r+1)*D1]+[0]*D1,p,g)for r in range(D2)]+[[0]*N1 for _ in range(D2)]
  gC=[[0]*N1 for _ in range(N2)]
  for c in range(N1):
   cA,cB=_ntt1([gA[r][c]for r in range(N2)],p,g),_ntt1([gB[r][c]for r in range(N2)],p,g)
   cC=_ntt1([(cA[r]*cB[r])%p for r in range(N2)],p,g,True)
   for r in range(N2):gC[r][c]=cC[r]
  cm.append([_ntt1(gC[r],p,g,True)for r in range(N2)])
 p1,p2,p3=pgs[0][0],pgs[1][0],pgs[2][0];i12,i123,p12=pow(p1%p2,p2-2,p2),pow((p1*p2)%p3,p3-2,p3),p1*p2;acc=[0]*(D1*(N2+1)+4)
 for r in range(N2):
  bs=r*D1
  for c in range(N1):a1,a2,a3=cm[0][r][c],cm[1][r][c],cm[2][r][c];x12=a1+p1*(((a2-a1)*i12)%p2);acc[bs+c]+=x12+p12*(((a3-x12)*i123)%p3)
 res=cy=0
 for i,v in enumerate(acc):s=v+cy;res|=(s&m58)<<(58*i);cy=s>>58
 return res|(cy<<(58*len(acc)))
def check_ntt()->dict:
 s=REF["s1b"];p1,p2,p3,g1,g2,g3=s["p1"],s["p2"],s["p3"],s["g1"],s["g2"],s["g3"]
 assert all(_ip64(p)and(p-1)%(1<<54)==0 and 8*p<(1<<64)and(1<<58)<p for p in(p1,p2,p3))
 assert all(pow(g1,(p1-1)//q,p1)!=1 for q in(2,127))and all(pow(g2,(p2-1)//q,p2)!=1 for q in(2,3,19))and all(pow(g3,(p3-1)//q,p3)!=1 for q in(2,3))
 l2P=math.log2(p1*p2*p3);h1,h2=l2P-181.0,l2P-(math.log2(math.ceil((2.0**67.77)/58.0))+116.0);assert h1>1.57 and h2>4.66
 rng,pgs=random.Random(0x5400612D),[(p1,g1),(p2,g2),(p3,g3)]
 for p,g in pgs:
  N=64;vec=[rng.randrange(p)for _ in range(N)];got=_ntt_r4(vec,p,g);wN=pow(g,(p-1)//N,p)
  for k in range(N):
   wk,acc,wj=pow(wN,k,p),0,1
   for j in range(N):acc=(acc+vec[j]*wj)%p;wj=(wj*wk)%p
   assert got[k]==acc
 sz=[512,1024,2048,4096]
 for nb in sz:A,B=rng.getrandbits(nb)|(1<<(nb-1)),rng.getrandbits(nb)|(1<<(nb-1));assert _mul2d(A,B,pgs)==A*B
 return{"ok":True,"r4_bitexact":True,"log2_P":round(l2P,6),"hr_1d":round(h1,4),"hr_2d":round(h2,4),"mul_2d":sz}
def _horner576(A:list[int],sm:list[int],v:list[int])->int:
 r256,r2_256,x2h,m64=R_BLS<<256,(2*R_BLS)<<256,0xAC45A4010001A402,(1<<64)-1
 T=A[0]+(A[2]<<32)+((A[2]*x2h)<<64)+r256-ABS_X*(A[1]+(A[3]<<32)+((A[3]*x2h)<<64))+sum(sm[k]*v[k]for k in range(4))
 assert T>=0;t8=T>>512;assert t8 in(0,1)
 if t8:T-=r2_256
 assert 0<=T<(1<<512);Th=T>>256;t4,t5,t6,t7=Th&m64,(Th>>64)&m64,(Th>>128)&m64,(Th>>192)&m64
 M=[0x42737A020C0D6393,0x65043EB4BE4BAD71,0x38B5DCB707E08ED3,0x355094EDFEDE377C]
 ph=((t4*M[3])+(t5*M[2]+((t5*M[3])<<64))+(t6*M[1]+((t6*M[2])<<64)+((t6*M[3])<<128))+(t7*M[0]+((t7*M[1])<<64)+((t7*M[2])<<128)+((t7*M[3])<<192)))>>64
 s=(T-((Th<<1)+ph)*R_BLS)&((1<<320)-1);assert 0<=s<6*R_BLS;s2=s-(s>>255)*R_BLS;assert 0<=s2<2*R_BLS
 return s2-R_BLS if s2>=R_BLS else s2
def check_kernels_e2e()->dict:
 pool=[((i&31)<<56)|(i*131+5)for i in range(17)];bo=[0,0,1,4,4,8,17];base,stride=SLAB_BASE,1<<56;tx_out=[]
 for j in range(6):
  for k in range(bo[j],bo[j+1]):tx_out.append(base+j*32*stride+pool[k])
 assert len(tx_out)==17 and tx_out[-1]==base+(5*32+16)*stride+(16*131+5)
 rc=(900+200)%1009;j0=rc>>5;assert rc==91 and j0==2
 Va,pb=1024,37;s_lo,s_hi=bytearray(Va),bytearray(Va)
 for i in range(Va):s_lo[i]=pb+(i*17)%60;s_hi[i]=(i*29)%50
 s_lo[13],s_hi[13]=80,47;s_lo[14],s_hi[14]=80,48;s_lo[255],s_hi[255]=150,105;s_lo[256],s_hi[256]=60,67
 exp_sv=[i for i in range(Va)if s_lo[i]+s_hi[i]>=128];assert all(s_lo[i]+s_hi[i]<=255 for i in range(Va))
 rx_lo=bytearray(s_lo);s_lo[:]=bytes([pb])*Va;got_sv=[]
 for blk in range(0,Va,128):
  acc=0
  for w in range(16):
   wl=int.from_bytes(rx_lo[blk+8*w:blk+8*w+8],"little");wh=int.from_bytes(s_hi[blk+8*w:blk+8*w+8],"little")
   acc|=(wl+wh)&((1<<64)-1)
  if acc&0x8080808080808080:
   for w in range(16):
    wl=int.from_bytes(rx_lo[blk+8*w:blk+8*w+8],"little");wh=int.from_bytes(s_hi[blk+8*w:blk+8*w+8],"little")
    m=(wl+wh)&0x8080808080808080
    while m:got_sv.append(blk+8*w+(((m& -m).bit_length()-1)>>3));m&=m-1
  s_hi[blk:blk+128]=bytes(128)
 assert got_sv==exp_sv and not any(s_hi)and all(x==pb for x in s_lo)
 rng=random.Random(0x576B1750)
 for t in range(60):
  A,sm,vv=([(1<<304)-1,0,(1<<304)-1,0],[R_BLS-1]*4,[R_BLS-1]*4)if t==0 else(([0,(1<<314)-1,0,(1<<314)-1],[0]*4,[0]*4)if t==1 else(([(1<<300)-1]*4,[R_BLS-1]*4,[R_BLS-1]*4)if t==2 else([rng.getrandbits(314)for _ in range(4)],[rng.randrange(R_BLS)for _ in range(4)],[rng.randrange(R_BLS)for _ in range(4)])))
  assert _horner576(A,sm,vv)==(A[0]+X_BLS*A[1]+X_BLS*X_BLS*A[2]+(X_BLS**3)*A[3]+sum(sm[k]*vv[k]for k in range(4)))%R_BLS
 Vv,r60=[rng.randrange(R_BLS)for _ in range(256)],R_BLS<<60
 z1=X_BLS%R_BLS;z2=(z1*z1)%R_BLS;z3=(z2*z1)%R_BLS;tot_ent=tot_c1=tot_cu=0
 for _ in range(60):
  Ap,Am,rr=[0,0,0,0],[0,0,0,0],0
  for ei in range(60):
   col=rng.randrange(256)
   if ei<12:tup=MONO[ei]
   elif rng.random()<0.70:tup=MONO[rng.randrange(12)]
   else:tup=tuple((0 if rng.random()<0.35 else(1 if rng.random()<0.5 else-1)*(1 if rng.random()<0.4 else rng.randrange(2,1<<29)))for _ in range(4))
   tot_ent+=1;val_fr=(tup[0]+tup[1]*z1+tup[2]*z2+tup[3]*z3)%R_BLS;rr=(rr+val_fr*Vv[col])%R_BLS
   for k in range(4):
    ck=tup[k]
    if ck==0:continue
    if abs(ck)==1:tot_c1+=1
    else:tot_cu+=1
    if ck>0:Ap[k]+=ck*Vv[col]
    else:Am[k]+=(-ck)*Vv[col]
  As=[Ap[k]+r60-Am[k]for k in range(4)]
  assert all(0<=Ap[k]<(1<<293)and 0<=Am[k]<(1<<293)and 0<=As[k]<(1<<315)for k in range(4))
  sm,vv=[rng.randrange(R_BLS)for _ in range(4)],[rng.randrange(R_BLS)for _ in range(4)]
  assert _horner576(As,sm,vv)==(rr+sum(sm[j]*vv[j]for j in range(4)))%R_BLS
 return{"ok":True,"chord_tx_hits":len(tx_out),"swar_surv":len(got_sv),"zero_carry_255":True,"horner576":"60/60","z12_4tuple_spmv":"60/60","comp_per_ent":round((tot_c1+tot_cu)/tot_ent,3)}
MONO=((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,0,1,0),(0,-1,0,1),(-1,0,0,0),(0,-1,0,0),(0,0,-1,0),(0,0,0,-1),(1,0,-1,0),(0,1,0,-1))
MMAP={v:i for i,v in enumerate(MONO)}
def _zmulm(a,k:int):
 k%=12;ng=k>=6
 if ng:k-=6
 c0,c1,c2,c3=a
 for _ in range(k):c0,c1,c2,c3= -c3,c0,c1+c3,c2
 return(-c0,-c1,-c2,-c3)if ng else(c0,c1,c2,c3)
def _z12_sge_hom1d(U:int,eta:float,sd:int)->dict:
 rng=random.Random(sd);Uf,lBf,lBg=int(round(0.4392*U)),REF["Bf"]*math.log(2.0),REF["Bg"]*math.log(2.0)
 Ug,rf1,rg1=U-Uf,DK.rho(6.65),DK.rho(8.31)
 ws=[min(0.45,(1.173/U)/((j+0.5)/Uf)*(DK.rho(7.65-max(0.02,1.0+math.log((j+0.5)/Uf)/lBf))/rf1)+(0.35/(2.0+2.5*j)if j<40 else 0.0))for j in range(Uf)]+[min(0.45,(1.196/U)/((j+0.5)/Ug)*(DK.rho(9.31-max(0.02,1.0+math.log((j+0.5)/Ug)/lBg))/rg1)+(0.40/(2.0+2.5*j)if j<40 else 0.0))for j in range(Ug)]
 mw,cdf,acc=sum(ws),[],0.0
 for w in ws:acc+=w/mw;cdf.append(acc)
 rows,Lp=[],math.exp(-mw)
 for _ in range(int(round(eta*U))):
  kp,pa=0,1.0
  while pa>Lp:kp+=1;pa*=rng.random()
  d={}
  for _ in range(max(4,kp-1)):
   c,m=min(U-1,bisect.bisect_left(cdf,rng.random())),MONO[rng.randrange(12)]
   if c in d:
    a=d[c];nv=(a[0]+m[0],a[1]+m[1],a[2]+m[2],a[3]+m[3])
    if nv==(0,0,0,0):del d[c]
    else:d[c]=nv
   else:d[c]=m
  if d:rows.append(d)
 cr={}
 for ri,d in enumerate(rows):
  for c in d:cr.setdefault(c,set()).add(ri)
 def kr(ri):
  d=rows[ri]
  if d is None:return
  for c in d:
   s=cr.get(c)
   if s:
    s.discard(ri)
    if not s:cr.pop(c,None)
  rows[ri]=None
 def psg():
  q=[c for c,s in cr.items()if len(s)==1]
  while q:
   c=q.pop();s=cr.get(c)
   if not s or len(s)!=1:continue
   ri=next(iter(s));d=rows[ri]
   if d is None:continue
   nb=list(d.keys());kr(ri)
   for nc in nb:
    ns=cr.get(nc)
    if ns and len(ns)==1:q.append(nc)
 psg();ar=[i for i,d in enumerate(rows)if d]
 while len(ar)-len(cr)>12:kr(max(ar,key=lambda i:len(rows[i])));psg();ar=[i for i,d in enumerate(rows)if d]
 while True:
  did=False
  for deg in(2,3,4):
   cd=None
   for c,s in list(cr.items()):
    if len(s)!=deg:continue
    rl=list(s);mi=[MMAP.get(rows[r][c],-1)for r in rl]
    if any(k<0 for k in mi):continue
    od=sorted(range(deg),key=lambda idx:len(rows[rl[idx]]));rp,kp=rl[od[0]],mi[od[0]];lp=len(rows[rp])
    if deg>=3 and any(lp+len(rows[rl[o]])>180 for o in od[1:]):continue
    cd=(rp,kp,[(rl[o],mi[o])for o in od[1:]]);break
   if cd is not None:
    rp,kp,tgs=cd;dp=rows[rp]
    for rt,kt in tgs:
     dt,sh=rows[rt],(kt-kp)%12
     for col,v0 in dp.items():
      sv=_zmulm(v0,sh)
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
 al=[d for d in rows if d];cols=sorted({c for d in al for c in d});mx=ov=nz=c1=cu=0
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
 z1=X_BLS%R_BLS;z2=(z1*z1)%R_BLS;z3=(z2*z1)%R_BLS;n_fb=len(cols);ncols=n_fb+4;cm={c:i for i,c in enumerate(cols)}
 v_true=[rng.randrange(1,R_BLS)for _ in range(ncols)];inv_last=pow(v_true[-1],-1,R_BLS);piv,rk={},0
 for d in al:
  cur={cm[c]:val for c,(c0,c1v,c2v,c3v)in d.items()if(val:=(c0+c1v*z1+c2v*z2+c3v*z3)%R_BLS)}
  s_dot=sum(val*v_true[idx]for idx,val in cur.items())%R_BLS
  for j in range(3):sv=rng.randrange(1,R_BLS);cur[n_fb+j]=sv;s_dot=(s_dot+sv*v_true[n_fb+j])%R_BLS
  s3=(-s_dot*inv_last)%R_BLS
  if s3:cur[n_fb+3]=s3
  while cur:
   pc=min(cur)
   if pc in piv:
    fac=cur[pc]
    for c,pv in piv[pc].items():
     nv=(cur.get(c,0)-fac*pv)%R_BLS
     if nv:cur[c]=nv
     else:cur.pop(c,None)
   else:iv=pow(cur[pc],-1,R_BLS);piv[pc]={c:(v*iv)%R_BLS for c,v in cur.items()};rk+=1;break
 nullity=ncols-rk;rec_ok=False
 if nullity==1:
  fcol=next(c for c in range(ncols)if c not in piv);v_rec=[0]*ncols;v_rec[fcol]=1
  for pc in sorted(piv.keys(),reverse=True):v_rec[pc]=(-sum(v*v_rec[c]for c,v in piv[pc].items()if c!=pc))%R_BLS
  sc=(v_true[fcol]*pow(v_rec[fcol],-1,R_BLS))%R_BLS;rec_ok=all((v_rec[i]*sc)%R_BLS==v_true[i]for i in range(ncols))
 return{"rows":len(al),"cols":ncols,"c1":c1,"cu":cu,"nz":nz,"mx":mx,"ov":ov,"nullity":nullity,"kernel_recovered":rec_ok}
def _genuine_sm_toy_demo()->dict:
 xt= -56;rt=xt**4-xt**2+1;r2=rt*rt
 def r3m(a,b,mod):
  c=[0]*5
  for i in range(3):
   for j in range(3):c[i+j]=(c[i+j]+a[i]*b[j])%mod
  return[(c[0]+c[3]*xt)%mod,(c[1]+c[4]*xt)%mod,c[2]%mod]
 def r3p(a,e,mod):
  res,bs=[1,0,0],[x%mod for x in a]
  while e:
   if e&1:res=r3m(res,bs,mod)
   bs=r3m(bs,bs,mod);e>>=1
  return res
 eps=(rt-1)if pow(xt%rt,(rt-1)//3,rt)==1 else(rt**3-1)
 def sm(a):
  p=r3p(a,eps,r2);assert(p[0]-1)%rt==0 and p[1]%rt==0 and p[2]%rt==0
  return[((p[0]-1)//rt)%rt,(p[1]//rt)%rt,(p[2]//rt)%rt]
 a1,a2=[17,5,-3],[23,-11,7];s1,s2,s12,sr=sm(a1),sm(a2),sm(r3m(a1,a2,r2)),sm(r3p(a1,rt,r2))
 assert s12==[(s1[i]+s2[i])%rt for i in range(3)]and sr==[0,0,0]and any(s1)
 return{"ok":True,"toy_r":rt,"homomorphism_verified":True,"rth_power_zero":True,"nontrivial":True}
def check_live_sge(full:bool=False)->dict:
 ns,fl,mxc,tov,tc1,tcu,tnz=(16 if full else 8),0,0,0,0,0,0
 for s in range(ns):
  r=_z12_sge_hom1d(512,REF["eta"],7500+s)
  if r["rows"]<r["cols"]or r["nullity"]!=1 or not r["kernel_recovered"]:fl+=1
  mxc=max(mxc,r["mx"]);tov+=r["ov"];tc1+=r["c1"];tcu+=r["cu"];tnz+=r["nz"]
 pm1_f=tc1/max(1,tc1+tcu);cpe=(tc1+tcu)/max(1,tnz)
 assert fl==0 and tov==0 and mxc<(1<<30)and pm1_f>=0.85
 u23,u24,la=REF["sge"][-2],REF["sge"][-1],REF["la"]
 assert la["nl"]>u23[1]+2.0*u23[2]and la["nl"]>u24[1]and la["comp"]>3.0*u24[3]
 out={"ok":True,"scope":"finite-scale 1D homogeneous kernel (nullity==1) & 4-basis expansion check (eta=1.015) + real toy Schirokauer map","seeds":ns,"fails":fl,"nullity_all_1":True,"kernel_recovered_all":True,"max_abs":mxc,"ov30":tov,"pm1_frac":round(pm1_f,4),"comp_per_ent":round(cpe,3),"la_hr_u23_x":round((la["nl"]*la["comp"])/(u23[1]*u23[3]),2),"la_hr_u24_x":round((la["nl"]*la["comp"])/(u24[1]*u24[3]),2),"toy_schirokauer_map":_genuine_sm_toy_demo()}
 if full:
  r1=_z12_sge_hom1d(1024,REF["eta"],9401);assert r1["nullity"]==1 and r1["kernel_recovered"]and r1["ov"]==0
  out["u1024_check"]={"rows":r1["rows"],"cols":r1["cols"],"nullity":r1["nullity"],"kernel_recovered":r1["kernel_recovered"],"max_abs":r1["mx"]}
 return out
def _mul_r4(Bb:float)->float:
 s,oe=REF["s1b"],REF["oe"];L1,L2=(oe*(2.0**60)*8.0)/4.5,(oe*(2.0**60)*8.0)/9.0
 if Bb<=1024.0:nl=max(1.0,Bb/64.0);return 12.0*nl*nl*s["saf"]
 if Bb<=float(1<<59):
  N=float(1<<int(math.ceil(math.log2(max(2.0,2.0*math.ceil(Bb/64.0))))));assert N<=(1<<54)
  return(3.0*(1.5*N*math.log2(N)*s["bfly"]+N*(s["pw"]+s["sc"]))+((3.0*N*160.0)if(Bb*2.0>L1)else 0.0)+N*s["c1"])*s["saf"]
 L=int(math.ceil(Bb/58.0));k=max(1,(L-1).bit_length());k1=(k+1)//2;k2=k-k1;N2=float(1<<(k1+1))*float(1<<(k2+1))
 return(3.0*(N2*(1.0*(k1+1)+1.5*(k2+1))*s["bfly"]+N2*(s["pw"]+s["sc"]))+((9.0*N2*160.0)if(Bb*2.0>L2)else 0.0)+N2*s["c2"])*s["saf"]
def _rem_r4(nlv:float,lb:float)->float:
 if nlv<=1.0:return 40000.0
 dm,kp=REF["s1b"]["dm"],max(1,int(4096.0//lb));ns,sb,nll=max(1.0,nlv/float(kp)),float(kp)*lb,lb/64.0
 isb=sum(12.0*(j*nll)*nll*1.15 for j in range(1,kp));isd=float(kp)*24.0*(sb/64.0)*nll*1.15
 h,ip,id_,cn,cb=int(math.ceil(math.log2(ns))),0.0,0.0,ns,sb
 for _ in range(h):pr=cn*0.5;mir=_mul_r4(cb);ip+=pr*mir;id_+=pr*dm*mir;cn,cb=pr,cb*2.0
 return ns*(isb+isd)+ip+3.0*_mul_r4(ns*sb)+id_+nlv*max(25000.0,40000.0*(lb/636.0))
def _opt_ecm(qb:float,Nb:float,_ec:dict={},_ep:dict={}):
 qk,key=round(qb,2),(round(qb,2),int(math.ceil(Nb/64.0)))
 if key in _ec:return _ec[key]
 if qk not in _ep:
  pl,qe=[],max(1.0,qb-3.5)
  for lb4 in range(32,249):
   B1,l1=2.0**(lb4*0.25),lb4*0.25;p1=DK.rho(qe/l1)
   for m in(10.0,25.0,60.0,120.0,250.0):
    B2=B1*m;dt=(math.log2(B2)-l1)/120.0
    p2=sum((DK.rho((qe-(l1+(i+0.5)*dt))/l1)*(dt/(l1+(i+0.5)*dt)))if(l1+(i+0.5)*dt)<qe else(dt/(l1+(i+0.5)*dt))for i in range(120))
    pr=min(1.0,p1+p2)
    if pr>1e-40:pl.append((B1,B2,pr))
  _ep[qk]=pl
 L=max(4,int(math.ceil(Nb/64.0)));mir=212.0 if L==4 else(13.25*L*L+80.0 if L<=16 else min(13.25*L*L+80.0,2.5*_mul_r4(L*64.0)));best=(1e300,None)
 for B1,B2,pr in _ep[qk]:
  c1=1.15*(40.0*mir+(B1/math.log(2.0))*(11.0*mir+93.0*L)+2.5*max(0.0,B2/math.log(B2)-B1/math.log(B1))*(mir+40.0*L))
  if c1/pr<best[0]:best=(c1/pr,(B1,B2,pr,c1))
 _ec[key]=best;return best
def _dnode_derived(Q:float,side:str,Bft:float,Bgt:float,Bfm:float,Bgm:float,af:float,ag:float)->dict:
 l2_det,l2_v16,C_R=2.0*math.log2(1125.0/16.0)+Q,math.log2((math.pi**8)/40320.0),2.0**20
 for _ in range(8):
  l2_Rq=(math.log2(C_R)+l2_det-l2_v16)/16.0+math.log2(1.25);l2_rs=l2_Rq-1.5
  nf=16.0+72.0*l2_rs-(Q if side=="f" else 0.0);ng=8.0*LOG2_ABS_X+24.0*l2_rs+2.0-(Q if side=="g" else 0.0)
  pj=max(1e-40,DK.rho((nf+af)/Bft)*DK.rho((ng+ag)/Bgt));C_R=max(1024.0,2.0/pj)
 cu=(C_R*6000.0+(2.0**22)+3.0*(_opt_ecm(Bft,nf)[0]+_opt_ecm(Bgt,ng)[0]))/C_COST
 def ec(Ne,lo,hi):
  if hi<=lo+1e-9:return 0.0
  r0,M=DK.rho(Ne/hi),40;dx=(hi-lo)/M
  return sum((DK.rho((Ne-(lo+(i+0.5)*dx))/hi)/r0)*(dx/(lo+(i+0.5)*dx))for i in range(M)if(lo+(i+0.5)*dx)<Ne)
 return{"nf":nf,"ng":ng,"cu":cu,"cf":ec(nf+af,Bfm,Bft),"cg":ec(ng+ag,Bgm,Bgt)}
def _s4_derived(af:float,ag:float)->dict:
 Bf,Bg=REF["Bf"],REF["Bg"]
 def br0(Nr,Qmax,Qlo,Qhi,M=200):
  r0,dx=DK.rho(Nr/Qmax),(Qhi-Qlo)/M
  return sum((DK.rho((Nr-(Qlo+(i+0.5)*dx))/Qmax)/r0)*(dx/(Qlo+(i+0.5)*dx))for i in range(M))
 n_sp=2.0*0.5*br0(2572.0,275.0,Bg,275.0);n_in0=1.0*0.5*br0(2572.0,275.0,0.5*Bg,137.5)
 n_in1=1.0*0.5*br0(2572.0,275.0,137.5,195.0);n_in2=1.0*0.5*br0(2572.0,275.0,195.0,275.0)
 d550=_dnode_derived(550.0,"g",390.0,390.0,Bf,Bg,af,ag);ut1=n_in2*d550["cu"]
 nf_cnt,ng_cnt=n_in2*d550["cf"],n_in2*d550["cg"]+n_in1;nd_log=[round(n_in2,2),round(nf_cnt+ng_cnt,1)]
 layers=[((390.0,275.0,275.0),n_sp+n_in0),((275.0,180.0,180.0),0.0),((180.0,125.0,125.0),0.0),((125.0,94.0,94.0),0.0),((94.0,78.0,78.0),0.0),((78.0,Bf,Bg),0.0)]
 for(Q,Bft,Bgt),add_g in layers:
  df=_dnode_derived(Q,"f",Bft,Bgt,Bf,Bg,af,ag);dg=_dnode_derived(Q,"g",Bft,Bgt,Bf,Bg,af,ag)
  ut1+=nf_cnt*df["cu"]+ng_cnt*dg["cu"]
  nf_cnt,ng_cnt=nf_cnt*df["cf"]+ng_cnt*dg["cf"],nf_cnt*df["cg"]+ng_cnt*dg["cg"]+add_g
  nd_log.append(round(nf_cnt+ng_cnt,1))
 assert nf_cnt==0.0 and ng_cnt==0.0
 ub1=((1.0/(0.918*DK.rho(2572.0/275.0)))*(2.5*_opt_ecm(275.0,2572.0)[0]))/C_COST
 return{"l0_split":round(n_sp,2),"l0_inert_le275":round(n_in0,2),"l0_inert_275_390":round(n_in1,2),"l0_inert_390_550":round(n_in2,2),"nodes":nd_log,"s4":math.log2(2.0*(ub1+ut1)),"ecm":math.log2(2.0*(ut1+0.25*ub1)),"desc_only":math.log2(2.0*ut1)}
def _lll(B:list[list[int]],delta:float=0.99)->list[list[int]]:
 n,B=len(B),[r[:]for r in B]
 def gso():
  bo,mu,B2=[[0.0]*len(B[0])for _ in range(n)],[[0.0]*n for _ in range(n)],[0.0]*n
  for i in range(n):
   v=[float(x)for x in B[i]]
   for j in range(i):
    mu[i][j]=sum(B[i][t]*bo[j][t]for t in range(len(v)))/B2[j]if B2[j]>1e-18 else 0.0
    for t in range(len(v)):v[t]-=mu[i][j]*bo[j][t]
   bo[i]=v;B2[i]=sum(x*x for x in v)
  return mu,B2
 mu,B2=gso();k=1
 while k<n:
  for j in range(k-1,-1,-1):
   if abs(mu[k][j])>0.5:
    q=int(round(mu[k][j]))
    for t in range(len(B[0])):B[k][t]-=q*B[j][t]
    mu,B2=gso()
  if B2[k]>=(delta-mu[k][k-1]**2)*B2[k-1]:k+=1
  else:B[k],B[k-1]=B[k-1],B[k];mu,B2=gso();k=max(k-1,1)
 return B
def _real_bls12_h90_and_q11_descent()->dict:
 xt= -56;rt=xt**4-xt**2+1;pt=((xt-1)**2*rt)//3+xt;xp=xt%pt
 assert rt==9831361 and pt==10647363907 and pt%15==7 and(pt**6+1)%rt==0
 def f12m(a,b):
  res=[0]*12
  for j1 in range(3):
   for j2 in range(3):
    j=j1+j2;mx=xp if j>=3 else 1;jo=j-3 if j>=3 else j;tmp=[0]*7
    for i1 in range(4):
     c1=a[4*j1+i1]
     if c1:
      for i2 in range(4):tmp[i1+i2]+=c1*b[4*j2+i2]
    tmp[0]+=tmp[5];tmp[1]+=tmp[6];c4=tmp[4]
    for i in range(4):res[4*jo+i]=(res[4*jo+i]+(tmp[i]-c4)*mx)%pt
  return res
 def f12p(a,e):
  res,bs=[1]+[0]*11,list(a)
  while e>0:
   if e&1:res=f12m(res,bs)
   bs=f12m(bs,bs);e>>=1
  return res
 def fr6(a):
  o=[0]*12
  for j in range(3):c0,c1,c2,c3=a[4*j:4*j+4];o[4*j:4*j+4]=[(c0-c1)%pt,(-c1)%pt,(c3-c1)%pt,(c2-c1)%pt]
  return o
 g0=f12p([3,5,7,11,13,17,19,23,29,31,37,41],(pt**6-1)*((pt**6+1)//rt))
 assert f12p(g0,rt)==[1]+[0]*11 and g0!=[1]+[0]*11
 tt=f12p(g0,4192837);assert f12m(tt,fr6(tt))==[1]+[0]*11
 A=[[0]*12 for _ in range(12)]
 for k in range(12):
  ek=[0]*12;ek[k]=1;img=[(ek[i]-f12m(tt,fr6(ek))[i])%pt for i in range(12)]
  for rr in range(12):A[rr][k]=img[rr]
 pv,ri=[],0
 for col in range(12):
  p0=next((rr for rr in range(ri,12)if A[rr][col]%pt!=0),None)
  if p0 is None:continue
  A[ri],A[p0]=A[p0],A[ri];iv=pow(A[ri][col],pt-2,pt);A[ri]=[(v*iv)%pt for v in A[ri]]
  for rr in range(12):
   if rr!=ri and A[rr][col]:f=A[rr][col];A[rr]=[(A[rr][cc]-f*A[ri][cc])%pt for cc in range(12)]
  pv.append(col);ri+=1
 assert len(pv)==6;fc=[c for c in range(12)if c not in pv]
 B=[[(pt if j==p_ else 0)for j in range(12)]for p_ in pv]+[[(1 if j==f_ else((-A[pv.index(j)][f_])%pt if j in pv else 0))for j in range(12)]for f_ in fc]
 ws=[1,1,1,1,4,4,4,4,15,15,15,15];Br=_lll([[r[i]*ws[i]for i in range(12)]for r in B],0.99)
 al=[Br[0][i]//ws[i]for i in range(12)];am=[c%pt for c in al]
 assert am==f12m(tt,fr6(am))and any(am)and f12m(f12p(am,pt**6-1),tt)==[1]+[0]*11
 xcr,w3,nv=(abs(xt)**(1.0/3.0))*cmath.exp(1j*math.pi/3.0),cmath.exp(2j*math.pi/3.0),1.0+0j
 for m in(1,2,3,4):
  z5=cmath.exp(2j*math.pi*m/5.0)
  for k3 in range(3):yv=xcr*(w3**k3);nv*=sum(al[4*j+i]*(z5**i)*(yv**j)for j in range(3)for i in range(4))
 l2n,l2d=math.log2(abs(nv)),6.0*math.log2(pt)+4.0*math.log2(abs(xt));assert l2n<l2d+2.0
 q=1031;assert q%15==11;z5r=next(t for t in range(2,q)if(t**4+t**3+t**2+t+1)%q==0)
 yr=pow(xt%q,pow(3,-1,q-1),q);assert(yr**3-xt)%q==0
 B16,p5=[[(q if j==idx else 0)for j in range(16)]for idx in(0,4)],[pow(z5r,i,q)for i in range(4)]
 for i in range(1,4):
  r0,r1=[0]*16,[0]*16;r0[i],r0[0]=1,(-p5[i])%q;r1[4+i],r1[4]=1,(-p5[i])%q;B16.extend((r0,r1))
 for i in range(4):
  r0,r1=[0]*16,[0]*16;r0[8+i],r0[0]=1,(yr*p5[i])%q;r1[12+i],r1[4]=1,(yr*p5[i])%q;B16.extend((r0,r1))
 v0=_lll(B16,0.99)[0]
 assert sum((v0[i]-yr*v0[8+i])*p5[i]for i in range(4))%q==0 and sum((v0[4+i]-yr*v0[12+i])*p5[i]for i in range(4))%q==0
 def okm(u,v):
  def z5m(a,b):
   t=[0]*7
   for i in range(4):
    for j in range(4):t[i+j]+=a[i]*b[j]
   t[0]+=t[5];t[1]+=t[6];c4=t[4];return[t[i]-c4 for i in range(4)]
  p00,p01,p10,p11=z5m(u[:4],v[:4]),z5m(u[:4],v[4:]),z5m(u[4:],v[:4]),z5m(u[4:],v[4:])
  return[p00[i]-p11[i]for i in range(4)]+[p01[i]+p10[i]-p11[i]for i in range(4)]
 a3,b3=okm(okm(v0[:8],v0[:8]),v0[:8]),okm(okm(v0[8:],v0[8:]),v0[8:]);gok=[a3[i]-xt*b3[i]for i in range(8)]
 A8,prev,sgn=[okm(gok,[(1 if j==k else 0)for j in range(8)])for k in range(8)],1,1
 for k in range(7):
  if A8[k][k]==0:
   sw=next(i for i in range(k+1,8)if A8[i][k]!=0);A8[k],A8[sw]=A8[sw],A8[k];sgn= -sgn
  piv=A8[k][k]
  for i in range(k+1,8):
   for j in range(k+1,8):A8[i][j]=(A8[i][j]*piv-A8[i][k]*A8[k][j])//prev
  prev=piv
 nlg=sgn*A8[7][7];assert nlg!=0 and nlg%(q*q)==0
 return{"p":pt,"r":rt,"rank_H90":6,"lift_verified":True,"log2_N_Mg":round(l2n,2),"log2_det_w":round(l2d,2),"q11_inert_descent":{"q":q,"q_mod_15":11,"exact_q2_divides_N_Lg":True,"log2_N_Lg_cof":round(math.log2(abs(nlg)/(q*q)),2)}}
def check_stage4(af:float,ag:float,full:bool=False)->dict:
 s4=_s4_derived(af,ag);assert s4["s4"]<REF["cr"]["both_individual_logs_exp"]and s4["ecm"]<REF["cr"]["ecm_split_expected"]
 h90=_real_bls12_h90_and_q11_descent()
 lll=33.0-6.0*math.log2(12.0)+4.5*math.log2(5.0)+6.0*math.log2(3.0)+6.0*math.log2(P_BLS)+4.0*math.log2(ABS_X);assert lll<2572.0
 out={"ok":True,"lll_ub":round(lll,4),"term":{"Bf":REF["Bf"],"Bg":REF["Bg"],"out_f":0,"out_g":0},"l0_split":s4["l0_split"],"l0_inert_le275":s4["l0_inert_le275"],"l0_inert_275_390":s4["l0_inert_275_390"],"l0_inert_390_550":s4["l0_inert_390_550"],"nodes":s4["nodes"],"desc_only_l2":round(s4["desc_only"],2),"s4_l2":round(s4["s4"],4),"ecm_l2":round(s4["ecm"],4),"real_bls12_h90_q11":h90}
 if not full:return out
 xd=lambda X,Z,A24,N:(((X+Z)**2*(X-Z)**2)%N,((((X+Z)**2-(X-Z)**2)%N)*(((X-Z)**2+A24*((X+Z)**2-(X-Z)**2))%N))%N)
 xa=lambda X1,Z1,X2,Z2,Xd,Zd,N:((Zd*(((X1-Z1)*(X2+Z2)+(X1+Z1)*(X2-Z2))**2))%N,(Xd*(((X1-Z1)*(X2+Z2)-(X1+Z1)*(X2-Z2))**2))%N)
 def xm(k,X,Z,A24,N):
  r0x,r0z,r1x,r1z=1,0,X%N,Z%N
  for bt in bin(k)[2:]:
   if bt=="0":r1x,r1z=xa(r0x,r0z,r1x,r1z,X,Z,N);r0x,r0z=xd(r0x,r0z,A24,N)
   else:r0x,r0z=xa(r0x,r0z,r1x,r1z,X,Z,N);r1x,r1z=xd(r1x,r1z,A24,N)
  return r0x,r0z
 nx,dx=5,-7;t16=[p for p in _primes(200)if p>=101]
 for p in t16:
  x=(nx*pow(dx,-1,p))%p;d=((2*x*x-1)*pow(x,-4,p))%p;A=(2*(1+d)*pow(1-d,-1,p))%p
  assert(1+sum(1 if(u**3+A*u*u+u)%p==0 else(2 if pow((u**3+A*u*u+u)%p,(p-1)//2,p)==1 else 0)for u in range(p)))%16==0
 N1=50023*50033;x1=(nx*pow(dx,-1,N1))%N1;A24=pow(1-((2*x1*x1-1)*pow(x1,-4,N1))%N1,-1,N1);X,Z=3,1
 for p in _primes(30):
  pk=p
  while pk*p<=30:pk*=p
  X,Z=xm(pk,X,Z,A24,N1)
 assert math.gcd(Z,N1)==1;bb=[xm(i,X,Z,A24,N1)for i in range(1,16)if math.gcd(i,30)==1];DX,DZ=xm(30,X,Z,A24,N1);gp,gc,w1=(1,0),(DX,DZ),1
 for k in range(1,13):
  for bx,bz in bb:
   g=math.gcd((gc[0]*bz-bx*gc[1])%N1,N1)
   if 1<g<N1:w1=g;break
  if w1>1:break
  gp,gc=gc,(xd(DX,DZ,A24,N1)if k==1 else xa(gc[0],gc[1],DX,DZ,gp[0],gp[1],N1))
 X2,Z2=xm(2,4,1,2,35);DX2,DZ2=xm(4,X2,Z2,2,35);w2=math.gcd((DX2*Z2-X2*DZ2)%35,35)
 assert w1==50023 and math.gcd(Z2,35)==1 and w2==5;out.update({"t16":f"{len(t16)}/{len(t16)}","w1":w1,"w2":[w2,35//w2]});return out
def _tx64(lam:float)->float:
 if lam<=0.0:return 0.0
 if lam>45.0:return(69.0/64.0)*lam+24.0
 ei,pk=0.0,math.exp(-lam)
 for k in range(1,160):
  pk*=lam/k;q,r=divmod(k,64)
  ei+=pk*(10.0+q*69.0+(35.0 if r&32 else 0)+(19.0 if r&16 else 0)+(11.0 if r&8 else 0)+(7.0 if r&4 else 0)+(3.0 if r&2 else 0)+(3.0 if r&1 else 0))
  if k>lam+25.0 and pk<1e-14:break
 return ei
T16=(5.0,10.0,9.0,13.0,15.0,20.0,19.0,23.0,23.0,28.0,27.0,31.0,33.0,38.0,37.0,41.0)
def _tx_chord_x16(lam:float)->float:
 if lam<=0.0:return 8.0
 if lam>45.0:return 11.0+(37.0/16.0)*lam+22.0
 p0=math.exp(-lam);e=p0*8.0;pk=p0
 for k in range(1,160):
  pk*=lam/k;q16,r=divmod(k,16);e+=pk*(11.0+37.0*q16+T16[r])
  if k>lam+24.0 and pk<1e-13:break
 return e
def _integ_hyb(plo,pw1,pw2,phi,L248,L124,dok,rmx,w2=16,w3=32,probe_p=20.0):
 wpc1,wpc2=1.0/(V_A_SH*L248),1.0/(V_A_SH*L124);ih=ip=nok=nb=ne=0.0
 if pw1>plo:
  dx=(pw1-plo)/100.0
  for i in range(100):
   x=plo+(i+0.5)*dx;p=2.0**x;d=(p/x)*dx;dk=d*dok;nok+=dk;nb+=dk*(p+2.0);ne+=d*V_A_SH
   ip+=dk*(probe_p+L248*7.0)*wpc1;ih+=dk*sum(pr*L248*(_tx64(r*V_A_SH/p)+1.5625*(r*V_A_SH/p))for r,pr in rmx)*wpc1
 for(a,b,W)in((pw1,pw2,w2),(pw2,phi,w3)):
  if b<=a:continue
  n_b,dx=(L124-1.0)/float(W)+1.0,(b-a)/60.0
  for i in range(60):
   x=a+(i+0.5)*dx;p=2.0**x;d=(p/x)*dx;dk=d*dok;nok+=dk;nb+=dk*(p/float(W)+2.0);ne+=d*V_A_SH
   eip=eih=0.0
   for r,pr in rmx:
    lb=float(W)*r*V_A_SH/p;pc=1.0-math.exp(-n_b*lb)
    eip+=pr*(probe_p+pc*6.0);eih+=pr*(n_b*_tx_chord_x16(lb)-(1.0-pc)*n_b*8.0+n_b*1.5625*lb)
   ip+=dk*eip*wpc2;ih+=dk*eih*wpc2
 return{"ih":ih,"ip":ip,"nok":nok,"nb":nb,"ne":ne}
def _tprob(cof,Pt,Bb,km=8):
 if cof<=Pt:return 1.0,0.0,1.0
 om=lambda v:0.0 if v<1.0 else(1.0/v if v<=2.0 else((1.0+math.log(v-1.0))/v if v<=3.0 else 0.57))
 dy,md,ma=(cof-Pt)/200.0,DK.rho(cof/Pt),0.0;mt=md
 for i in range(200):
  y=Pt+(i+0.5)*dy;w=DK.rho((cof-y)/Pt)*(om(y/Pt)/Pt)*dy;mt+=w
  if y<=Bb:md+=w
  elif any(k*Pt<=y<=k*Bb for k in range(2,km+1)):ma+=w
 return md/max(mt,1e-30),ma/max(mt,1e-30),(md+ma)/max(mt,1e-30)
def check_ledger(lm:dict,lis:dict)->dict:
 Bf,Bg,eta=REF["Bf"],REF["Bg"],REF["eta"];b,rep,la=REF["b"],REF["rep"],REF["la"]
 af,ag,cpad=lm["af_1300k"],lm["ag_1300k"],lm["cof_pad"]
 L248,L124,gap_pad,l2sg,l2sf,l2s1,ret=lis["L248_ch"],lis["L124_ch"],lis["gap_pad"],lis["l2_sg"],lis["l2_sf"],lis["l2_s1"],lis["ret_live"]
 dokg,dokf,dng,dnf=1.0/3.0,23.0/81.0,1.0/24.0,2131.0/34992.0;w0=1.0/(1.0+REF["rbr"]**8);w1=1.0-w0
 fmx=[(3,39.0/46.0),(6,6.0/46.0),(9,1.0/46.0)]
 sg0=_integ_hyb(b["pglo"],b["pw1"],b["pw2"],b["pghi0"],L248,L124,dokg,[(3,1.0)],REF["w2"],REF["w3"])
 sf0=_integ_hyb(b["pflo"],b["pw1"],b["pw2"],b["pfhi0"],L248,L124,dokf,fmx,REF["w2"],REF["w3"])
 sg1=_integ_hyb(b["pglo"],b["pw1"],b["pw2"],b["pghi1"],L248,L124,dokg,[(3,1.0)],REF["w2"],REF["w3"])
 sf1=_integ_hyb(b["pflo"],b["pw1"],b["pw2"],b["pfhi1"],L248,L124,dokf,fmx,REF["w2"],REF["w3"])
 ih=w0*(sg0["ih"]+sf0["ih"])+w1*(sg1["ih"]+sf1["ih"]);ip=w0*(sg0["ip"]+sf0["ip"])+w1*(sg1["ip"]+sf1["ip"])
 def rsh(sg,sf):
  ne,nb=sg["ne"]+sf["ne"],sg["nb"]+sf["nb"]
  return(ne*rep["ent"]+nb*rep["bkt"])/C_COST,math.log2(8.0*(ne+nb)+48.0*(sg["nok"]+sf["nok"]))
 rep0,mem0=rsh(sg0,sf0);rep1,mem1=rsh(sg1,sf1);repu=rep0+rep1
 l2Vk=lis["Vk_live"]+lis["pad_3se_Vk"]+gap_pad+cpad;Vt,LtM,lf,lg=2.0**l2Vk,0.0,REF["lf"],REF["lg"]
 def ps1b(pfm,pgm):
  nonlocal LtM
  s1=2.0**l2s1;c1,c2,mg=pfm+(Bf-pfm)/3.0,pfm+2.0*(Bf-pfm)/3.0,0.5*(pgm+Bg)
  _,af1,_=_tprob(lf,pfm,Bf);_,af2,_=_tprob(lf,c1,Bf);_,af3,_=_tprob(lf,c2,Bf);_,ag1,_=_tprob(lg,pgm,Bg);_,ag2,_=_tprob(lg,mg,Bg)
  rf=min(0.95,math.sqrt(ret/0.80));rg=ret/rf
  pfs=min(1.0,(DK.rho((lis["mean_nf"]+af)/Bf)*rf)/(2.0**l2sf));pgs=min(1.0,(DK.rho((lis["mean_ng"]+ag)/Bg)*rg)/(2.0**l2sg))
  cr=s1*1700.0+s1*pfs*1700.0
  for lo,hi,pin,lb,dn in[(1.0,pfm,s1,656.0,dnf),(pfm,c1,s1*af1,lf,dnf),(c1,c2,s1*af2,lf,dnf),(c2,Bf,s1*af3,lf,dnf),(1.0,pgm,s1*pfs,688.0,dng),(pgm,mg,s1*pfs*ag1,lg,dng),(mg,Bg,s1*pfs*ag2,lg,dng)]:
   Lt=max(1.0,(ppb_dus(lo,min(hi,36.0))if lo<36.0 else 0.0)+(dn*ppb_dus(max(lo,36.0),hi)if hi>36.0 else 0.0))
   LtM=max(LtM,Lt);Mt=Vt*pin;Mb=min(Mt,max(1.0,Lt/lb));nbb,LX=max(1.0,Mt/Mb),Mb*lb;ir=_rem_r4(Mb,lb)
   if LX<Lt:ir+=max(0.0,Lt/LX-1.0)*2.5*_mul_r4(LX)
   cr+=(nbb*ir)/Vt
  return cr+s1*pfs*pgs*(110.0*256.0*150.0*4.0+0.05*5.0e7*2.0)
 cs1b=w0*ps1b(b["pfhi0"],b["pghi0"])+w1*ps1b(b["pfhi1"],b["pghi1"])
 sl_base=52.0/256.0+2.5*(19.0/256.0)+64.0/128.0
 sll=sl_base+200.0*(2.0**l2sg);slc=sll+0.004;ipcl,ipcc=ih+ip+sll+cs1b,ih+ip+slc+cs1b
 c_l,c_c=(Vt*ipcl)/C_COST,(Vt*ipcc)/C_COST;U=U_ub(Bf,Bg);Nl,Nc=U*la["nl"],U*la["comp"]
 u23,u24=REF["sge"][-2],REF["sge"][-1]
 assert la["nl"]>u23[1]+2.0*u23[2]and la["nl"]>u24[1]and la["comp"]>3.0*u24[3]
 bsp=la["pm1_f"]*10.0+(1.0-la["pm1_f"])*25.0+(Nl*36.0)/Nc;assert bsp<la["x4"]and(76.0+872.0)<la["rv"]
 Mb,Nb=la["m"],la["n"];Lsq,dmk=Nl/Mb+Nl/Nb+64.0,Nl/Nb+32.0
 ibw=1.02*(Nb*Lsq*(Nc*la["x4"]+Nl*la["rv"])*la["ov"]+dmk*(Nc*la["x1"]+Nl*la["rv"])*la["ov"])
 imk=1.02*dmk*Nb*Nl*283.0*la["ov"];ilg=20.0*((Mb+Nb)**3)*Lsq*math.log2(Lsq)*236.0*la["ov"]
 ilg_xf=4.0*Mb*(Mb+Nb)*Lsq*math.ceil(math.log2(Lsq))*8.0*la["ov"]
 lau=(ibw+imk+ilg+ilg_xf+0.008*(ibw+imk)+30.0*U*283.0*la["ov"])/C_COST
 sgeu=((12500.0*U+25.0*U*math.log2(U)+100.0*Nc)*1.15)/C_COST;s4=_s4_derived(af,ag)
 lr={"relation_collection":math.log2(c_l),"replay_a_tables":math.log2(repu),"linear_algebra":math.log2(lau),"both_individual_logs_exp":s4["s4"],"ecm_split_expected":s4["ecm"],"meter_safety_buffer":93.4,"band_split_and_trees":90.0,"canonical_a_domain":84.0,"schirokauer_maps":80.0,"final_log_backsub_and_gcd":72.0,"matrix_generator_sge":math.log2(sgeu),"stage0_setup_and_iso":45.0,"final_verification_Q":40.0}
 cr=REF["cr"];assert len(lr)==len(cr)==13 and all(cr[k]>lr[k]for k in lr)and cr["relation_collection"]>=math.log2(c_c)
 wol,woc=sum(2.0**v for v in lr.values()),sum(2.0**v for v in cr.values());lol,loc=math.log2(wol),math.log2(woc)
 pb=math.ceil((2.0**mem0)/(12.0*(2.0**60)));gsB=248.0*V_A_SH*(2.0**l2sg)*8.0
 half_slab_eib=(186.0*(2.0**56))/(2.0**60)
 s1p=math.log2((2.0**mem0)+248.0*(V_A_SH+262144.0)+2.0*186.0*(2.0**56)+gsB)
 ntB,cb,nn,lv=0.0,4080.0,LtM/4080.0,0
 while nn>=1.0:
  if cb<=(1<<59):ntB+=nn*(float(1<<int(math.ceil(math.log2(max(2.0,2.0*math.ceil(cb/64.0))))))*24.0)
  else:L2=int(math.ceil(cb/58.0));k=max(1,(L2-1).bit_length());k1=(k+1)//2;ntB+=nn*(float(1<<(k1+1))*float(1<<(k-k1+1))*24.0)
  cb*=2.0;nn*=0.5;lv+=1
 s1bB=ntB+(9.0+lv)*LtM/8.0;s1b_oe=s1bB/(2048.0*(2.0**60));s2r=math.log2((eta*U)*18.0*12.0)
 lmB=Nc*(la["pm1_f"]*8.0+(1.0-la["pm1_f"])*12.0)+Nl*64.0;lvB=Nl*32.0
 la_oe=(lmB/la["p"]+lvB*(1.0+(Nb+6.0)/la["p"]))/(2.0**60)
 d_pi=Lsq*Nb/(Mb+Nb)
 lg_totB=(Mb*Nb*Lsq+Mb*(Mb+Nb)*Lsq+4.0*((Mb+Nb)**2)*d_pi)*32.0
 lg_oe=lg_totB/(la["plg"]*(2.0**60));s3r=math.log2(lmB+(Nb+6.0+la["p"])*lvB+lg_totB)
 assert V_A_SH/(2.0**60)<0.06 and half_slab_eib<12.0<16.0 and s1b_oe<REF["oe"]<16.0 and la_oe<6.0<16.0 and lg_oe<8.0<16.0 and max(s1p,math.log2(s1bB),s2r,s3r)<REF["ram_ub"]<120.0
 pm=woc/(2.0**124.0);pf=pm+2.0**(-48)+2.0**(-100)+2.0**(-132)+2.0**(-200);po=1.0-pf;Kg,t64=REF["Kg"],1<<64
 assert pf<2.0e-5 and po>0.99998 and Fraction(Kg,t64)*Fraction(99998,100000)>=Fraction(39,100)and Fraction(Kg-1,t64)*Fraction(99998,100000)<Fraction(39,100)
 pc=(Kg/float(t64))*po;num,den=REF["num"],REF["den"];lgc,lgl=loc+math.log2(Kg/float(t64)),lol+math.log2(Kg/float(t64))
 slack=float(Fraction(num,den))-lgc
 assert pc>0.39 and math.gcd(num,den)==1 and slack>2.0*lis["se_Vk"]
 return{
 "ok":True,"claim":f"{num}/{den}","eta":eta,"U_l2":round(math.log2(U),6),"U_ub":float(f"{U:.8e}"),
 "Vk_ch_live":round(l2Vk,6),"ih":round(ih,6),"ip":round(ip,6),"c_s1b":round(cs1b,6),
 "slab_live":round(sll,6),"slab_ch":round(slc,6),"ipc_live":round(ipcl,6),"ipc_ch":round(ipcc,6),
 "blended_spmv_comp":round(bsp,4),"la_mat_eib":round(lmB/(2.0**60),4),"la_owner_eib":round(la_oe,4),
 "P_lingen":la["plg"],"lingen_owner_eib":round(lg_oe,4),"half_slab_eib":round(half_slab_eib,4),
 "opened_live":round(lol,6),"opened_ch":round(loc,6),"gated_live":round(lgl,6),"gated_ch":round(lgc,6),
 "slack_b":round(slack,6),"slack_over_live_SE_x":round(slack/lis["se_Vk"],2),
 "p_meter_cond_H":float(f"{pm:.7e}"),"p_open_cond_H":round(po,8),"p_cold_cond_H":round(pc,8),
 "ledger_live":{k:round(v,6)for k,v in lr.items()},
 "ram_l2_B":{"S1_sh0":round(s1p,4),"S1_sh1":round(mem1,4),"P_banks":pb,"S1b":round(math.log2(s1bB),4),"S1b_eib":round(s1b_oe,4),"S2":round(s2r,4),"S3":round(s3r,4),"S3_lingen":round(math.log2(lg_totB),4)},
 }
def main()->None:
 t0=time.perf_counter();raw=sys.stdin.read().strip();inp=json.loads(raw)if raw else{}
 if not isinstance(inp,dict)or(set(inp.keys())-{"part"}):raise ValueError("Input must be a JSON object with optional key 'part'")
 part=inp.get("part","core")
 if part not in("core","diagnostics","all"):raise ValueError(f"Unsupported part: {part}")
 out:dict={"package_id":REF["pkg"],"config_tag":REF["tag"],"part":part,"reads_external_files":0,"offline_mc_constants":0}
 if part in("core","all"):
  out["stage0"]=check_stage0()
  out["chebotarev"]=check_chebotarev(26000)
  lm=check_murphy_norms(2000000 if part=="all" else 1300000);out["murphy_norms"]=lm
  lis=check_live_panjer_is_chords(lm["af_1300k"],lm["ag_1300k"],32000,512000,8000)
  out["live_panjer_is_chords"]=lis
  out["ntt_r4_shoup61"]=check_ntt()
  out["kernels_e2e"]=check_kernels_e2e()
  out["live_strict_z12_sge"]=check_live_sge(full=(part=="all"))
  out["stage4_derived"]=check_stage4(lm["af_1300k"],lm["ag_1300k"],full=(part=="all"))
  out["ledger_ram_budget"]=check_ledger(lm,lis)
 elif part=="diagnostics":
  lm=check_murphy_norms(2000000);out["murphy_norms"]=lm
  out["live_panjer_is_chords"]=check_live_panjer_is_chords(lm["af_1300k"],lm["ag_1300k"],32000,512000,8000)
  out["live_strict_z12_sge"]=check_live_sge(full=True)
  out["stage4_derived"]=check_stage4(lm["af_1300k"],lm["ag_1300k"],full=True)
 out["checks_ok"]=all(v.get("ok",True)for v in out.values()if isinstance(v,dict))
 out["elapsed_seconds"]=round(time.perf_counter()-t0,3)
 print(json.dumps(out,indent=2))
 if not out["checks_ok"]:sys.exit(1)
if __name__=="__main__":
 main()
