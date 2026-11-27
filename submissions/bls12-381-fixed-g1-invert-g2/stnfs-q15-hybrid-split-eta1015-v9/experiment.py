#!/usr/bin/env python3
"""v9 verifier."""
import bisect,cmath,json,random,sys,time
from fractions import Fraction
from math import ceil,cos,exp,floor,gamma,gcd,lgamma,log,log2,pi,sin,sqrt
RN,BA,EN,IS,FR,RR,FL,RD,SR,L2,LG,EX,SQ,CE,F_,LN,SM,MX,MN,AB,PW,LS,ST,I,T_,F0,N_=range,bytearray,enumerate,isinstance,Fraction,random.Random,float,round,sorted,log2,log,exp,sqrt,ceil,floor,len,sum,max,min,abs,pow,list,set,int,True,False,None
XB=-0xD201000000010000
AX,L2X=AB(XB),L2(AB(XB))
PB=0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
RB=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P0_X=0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb
P0_Y=0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1
Q0_X=(0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e)
Q0_Y=(0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be)
CC0,VA0,SB0=20000.,.995*(2.**56),1<<56
R_={
"pkg":"stnfs-q15-hybrid-split-eta1015-v9","tag":"CAND9_HYBRID_RQMC_RSA768PARETO_PLA",
"Bf":71.25,"Bg":71.61,"eta":1.015,"cush":-.0075,"disc_pad":.0020,
"qa":.65,"rbr":.82,"Ra":1732.8673461468502,"R0":1355.3079808905318,"sl":1.75,
"lam":340599.43367449316,"lf":408.,"lg":688.,"oe":10.,"w2":16,"w3":32,
"b":{"pglo":28.,"pw1":55.,"pw2":58.5,"pghi0":62.5,"pghi1":60.7,"pflo":32.,"pfhi0":62.1,"pfhi1":60.3,"qg":3.4,"qf":3.8},
"s1b":{"p1":(127<<54)+1,"g1":3,"p2":(114<<54)+1,"g2":7,"p3":(108<<54)+1,"g3":5,"bfly":11.125,"pw":12.,"sc":12.,"c1":58.,"c2":90.,"saf":1.10,"dm":1.67},
"rep":{"ent":21.50,"bkt":5.5},
"la":{"m":1024,"n":64,"p":4,"plg":1024,"e_core":64,"x4":11.45,"x1":14.45,"rv":960.,"ov":1.06},
"cr":{"relation_collection":107.9425,"replay_a_tables":104.5398,"linear_algebra":103.5665,"both_individual_logs_exp":96.6500,"ecm_split_expected":96.0500,"meter_safety_buffer":93.4001,"band_split_and_trees":90.0001,"canonical_a_domain":84.0001,"schirokauer_maps":80.0001,"final_log_backsub_and_gcd":72.0001,"matrix_generator_sge":62.9600,"stage0_setup_and_iso":45.0001,"final_verification_Q":40.0001},
"Kg":7194374076228249696,"num":21357,"den":200,"ram_ub":118.,
}
class EDk:
 def __init__(s,um=30.,N=400):
  s.um,s.N,h,t=um,N,1./N,I(RD(um*N))+1;v=[1.]*(N+1)+[1.-LG(i*h)for i in RN(N+1,MN(2*N+1,t))]+[0.]*MX(0,t-2*N-1)
  if t>2*N+1:
   w=SM(v[N+2:2*N+1])
   for i in RN(2*N+1,t):
    if i%N==0:w=SM(v[i-N+1:i])
    v[i]=(h*(.5*v[i-N]+w))/(i*h-.5*h);w+=v[i]-v[i-N+1]
  s.v,s.lv=v,[LG(x)if x>1e-300 else-700. for x in v]
 def rho(s,u):
  if u<=1.:return 1.
  if u>=s.um:return 1e-300
  p=u*s.N;i=I(p);return s.v[-1]if i+1>=LN(s.v)else EX((1.-(p-i))*s.lv[i]+(p-i)*s.lv[i+1])
DK=EDk(30.,400);D0=DK.rho
def _primes(n):
 ip=BA(b"\x01")*(n+1);ip[0]=ip[1]=0
 for i in RN(2,I(n**.5)+1):
  if ip[i]:ip[i*i:n+1:i]=BA((n-i*i)//i+1)
 return[i for i in RN(2,n+1)if ip[i]]
def pi_dus(B):L=B*LG(2.);return(2.**B/L)*(1.+1./L+2.334/(L*L))
def U_ub(Bf,Bg):return(pi_dus(Bf)+pi_dus(Bg))/12.+36.*pi_dus(.5*Bf)+12.*pi_dus(.5*Bg)
def ppb_dus(lo,hi):th=2.**hi*(1.+.01/MX(1.,hi*LG(2.)));tl=0. if lo<=1. else 2.**lo*(1.-.01/MX(1.,lo*LG(2.)));return MX(1.,(th-tl)/LG(2.))
def _g1a(p1,p2):
 if p1 is N_:return p2
 if p2 is N_:return p1
 x1,y1,x2,y2=p1[0],p1[1],p2[0],p2[1]
 if x1==x2:
  if(y1+y2)%PB==0:return N_
  l=(3*x1*x1*PW(2*y1,-1,PB))%PB
 else:l=((y2-y1)*PW(x2-x1,-1,PB))%PB
 x3=(l*l-x1-x2)%PB;return(x3,(l*(x1-x3)-y1)%PB)
def chk_s0():
 assert(XB-1)**2*(XB**4-XB**2+1)//3+XB==PB and XB**4-XB**2+1==RB;l=XB%RB;assert PB%RB==l and(PW(l,4,RB)-PW(l,2,RB)+1)%RB==0 and PW(l,6,RB)==RB-1 and PW(l,12,RB)==1!=PW(l,4,RB);rp,ap,k=N_,(P0_X,P0_Y),RB
 while k:
  if k&1:rp=_g1a(rp,ap)
  ap=_g1a(ap,ap);k>>=1
 m2=lambda a,b:((a[0]*b[0]-a[1]*b[1])%PB,(a[0]*b[1]+a[1]*b[0])%PB);q3=m2(m2(Q0_X,Q0_X),Q0_X);assert(P0_Y*P0_Y-(P0_X**3+4))%PB==0 and rp is N_ and m2(Q0_Y,Q0_Y)==((q3[0]+4)%PB,(q3[1]+4)%PB)and PB%15==7 and[j for j in RN(1,5)if PW(PB,j,5)==1][0]==4 and PW(XB%PB,(PB-1)//3,PB)!=1 and any((XB**3-XB**2+(c*c%PB)*XB-c)%PB==0 for g in RN(2,20)if(cd:=PW(g,(PB-1)//3,PB))!=1 for c in(cd,(cd*cd)%PB))and XB**6-2*XB**5+2*XB**3+XB+1==3*PB and(RB-1)%12==0 and AX==(1<<16)*906349*254760293
 return{"ok":T_,"p_mod_15":7,"ord_5_p":4,"n_sm_chi1":4,"abs_x":"2^16*906349*254760293"}
def _tr(a):
 while LN(a)>1 and a[-1]==0:a.pop()
 return a
def _pm(a,b,p):
 a,b=_tr(LS(a)),_tr(LS(b));iv,db=PW(b[-1],-1,p),LN(b)-1
 while LN(a)-1>=db and a!=[0]:
  c,sh=(a[-1]*iv)%p,LN(a)-1-db
  for k in RN(db+1):a[sh+k]=(a[sh+k]-c*b[k])%p
  a.pop();a=a or[0];_tr(a)
 return _tr(a)
def _pg(u,v,p):
 u,v=_tr(LS(u)),_tr(LS(v))
 while v!=[0]:u,v=v,_pm(u,v,p)
 iv=PW(u[-1],-1,p);return[(x*iv)%p for x in u]
def _sr(F,p):
 d=LN(F)-1
 if d==0:return[]
 if d==1:return[((-F[0])*PW(F[1],p-2,p))%p]
 def ml(x,y):
  r=[0]*(LN(x)+LN(y)-1)
  for i,vx in EN(x):
   if vx:
    for j,vy in EN(y):r[i+j]=(r[i+j]+vx*vy)%p
  return _pm(r,F,p)
 a=1
 while T_:
  rs,bs,e=[1],[a%p,1],(p-1)//2
  while e:
   if e&1:rs=ml(rs,bs)
   bs=ml(bs,bs);e>>=1
  rs[0]=(rs[0]-1)%p;g=_pg(F,rs,p)
  if 1<=LN(g)-1<d:
   F2,d1=LS(F),LN(g)-1;q=[0]*(d-d1+1)
   for k in RN(d-d1,-1,-1):
    c=F2[k+d1];q[k]=c
    for j in RN(d1+1):F2[k+j]=(F2[k+j]-c*g[j])%p
   return _sr(g,p)+_sr(_tr(q),p)
  a+=1
def _fcr(w,p,e3):
 w2,r0,r1,r2,mk=(w*w)%p,0,1,0,1<<(p.bit_length()-2)
 while mk:
  c4=r2*r2;c3=(2*r1*r2+c4)%p;c2=(2*r0*r2+r1*r1-c4*w2)%p;c1=(2*r0*r1+c4*w)%p;c0=(r0*r0)%p;r0=(c0+c3*w)%p;r1=(c1-c3*w2)%p;r2=(c2+c3)%p
  if p&mk:r0,r1,r2=(r2*w)%p,(r0-r2*w2)%p,(r1+r2)%p
  mk>>=1
 r1=(r1-1)%p
 if r0==0 and r1==0 and r2==0:
  return SM(3 for rt in _sr([(-w)%p,w2,p-1,1],p)if PW(rt,e3,p)==1)
 if r2:
  iv=PW(r2,-1,p);a2=(p-1-iv*r1)%p;a1=(w2-iv*r0)%p;c=(a2*iv)%p;b1=(a1-c*r1)%p;b0=((-w)-c*r0)%p
  if b1==0:return 0
  rt=((-b0)*PW(b1,-1,p))%p;return 3 if((r0+rt*(r1+rt*r2))%p==0 and PW(rt,e3,p)==1)else 0
 if r1:
  rt=((-r0)*PW(r1,-1,p))%p;return 3 if((((rt-1)*rt+w2)*rt-w)%p==0 and PW(rt,e3,p)==1)else 0
 return 0
def chk_chb(pl=18000):
 zm=lambda x,y:(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1]);zn=lambda x:x[0]*x[0]-x[0]*x[1]+x[1]*x[1];zd=lambda d,n:all(c%zn(d)==0 for c in zm(n,(d[0]-d[1],-d[1])));dw,dw2,p19,p67=(41,24),(17,-24),(-5,-2),(-9,-2);pw,ww=991,113;wr=[t for t in RN(pw)if(t**3-t**2+(ww*ww%pw)*t-ww)%pw==0];kc=PW((wr[0]*PW(wr[1],-1,pw))%pw,(pw-1)//3,pw)
 assert zm(p19,p67)==dw and zn(dw)==1273==19*67 and gcd(1273,15)==1 and not zd(p19,dw2)and not zd(p67,dw2)and(118*118+118+1)%151==0 and all((t**3-t**2+32*t-118)%151 for t in RN(151))and pw%45==1 and(ww*ww+ww+1)%pw==0 and wr==[373,802,808]and kc==877 and PW(kc,3,pw)==1;s3=((0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0));G=[((a,b,c),s)for a in RN(3)for b in RN(3)for c in RN(3)for s in s3];gid={g:i for i,g in EN(G)}
 gm=lambda a,b:(((a[0][b[1][0]]+b[0][0])%3,(a[0][b[1][1]]+b[0][1])%3,(a[0][b[1][2]]+b[0][2])%3),(a[1][b[1][0]],a[1][b[1][1]],a[1][b[1][2]]));mt=[[gid[gm(a,b)]for b in G]for a in G];it=[next(j for j in RN(162)if mt[i][j]==0)for i in RN(162)];sn,cl=ST(),[]
 for i in RN(162):
  if i not in sn:cc=frozenset(mt[mt[h][i]][it[h]]for h in RN(162));sn|=cc;cl+=[cc]
 hz=frozenset(i for i,(v,_)in EN(G)if SM(v)%3==0);ns,q={frozenset([0])},[frozenset([0])]
 while q:
  N=q.pop()
  for cc in cl:
   if cc<=N:continue
   M=ST(N)|ST(cc);fr=ST(cc)
   while fr:nx={mt[a][b]for a in M for b in fr}-M;M|=nx;fr=nx
   if(Mf:=frozenset(M))not in ns:ns.add(Mf);q+=[Mf]
 sh=[N for N in ns if N<hz];mx=[N for N in sh if not any(N<M for M in sh)];assert LN(G)==162 and LN(cl)==22 and SR(LN(N)for N in ns)==[1,3,9,27,27,54,81,162]and LN(mx)==1 and LN(mx[0])==27;rc,c0,c1,r0,r1={0:0,3:0,6:0,9:0},0,0,[],[]
 for v,s in G:
  fx=SM(3 for i in RN(3)if s[i]==i and v[i]==0);rc[fx]+=1
  if SM(v)%3==0:r0+=[fx];c0+=I(fx>0)
  else:r1+=[fx];c1+=I(fx>0)
 p0,p3,p6,p9=(FR(rc[k],162)for k in(0,3,6,9));dok,d0,dnz=1-p0,FR(c0,54),FR(c1,108);df=FR(1,8)*(FR(1,3)*(1-(1-d0)**2)+FR(2,3)*(1-(1-dnz)**2));m0,m1=SM(r0)/54.,SM(r1)/108.;v0,v1=SM((x-m0)**2 for x in r0)/54.,SM((x-m1)**2 for x in r1)/108.;assert(p0,p3,p6,p9)==(FR(58,81),FR(13,54),FR(1,27),FR(1,162))and dok==FR(23,81)and df==FR(2131,34992)and AB(((1./3.)*16.*(2.*v0+(2.*m0-2.)**2)+(2./3.)*16.*(2.*v1+(2.*m1-2.)**2))/8.-12.)<1e-9
 ps=_primes(pl);p15=[p for p in ps if p%15==1];assert LN(p15)>200;bm=fp=ca=0
 for p in p15:
  e3=(p-1)//3;w=next(c for g in RN(2,p)if(c:=PW(g,e3,p))!=1);ph=F0
  for wv in(w,(w*w)%p):
   ny=_fcr(wv,p,e3);bm+=I(ny%3!=0)
   if ny:ca+=1;ph=T_
  fp+=I(ph)
 od,odn=ca/(2.*LN(p15)),fp/FL(LN(ps));assert bm==0 and AB(od-FL(dok))<.04 and AB(odn-FL(df))<.01;return{"ok":T_,"delta_Q":"41+24*w","kummer_p":pw,"galois_order":162,"d_ok_f":str(dok),"dens_f":str(df),"var_f_root":12.,"var_g_root":16.,"split_p":LN(p15),"live_d_ok":RD(od,5),"live_dens":RD(odn,5)}
def _prc(pl,qe,z,o,ad,sb,ml,iv):
 def tf(a):
  while LN(a)>1 and a[-1]==z:a.pop()
  return a
 def pm(a,b):
  r=[z]*(LN(a)+LN(b)-1)
  for i,x in EN(a):
   if x!=z:
    for j,y in EN(b):
     if y!=z:r[i+j]=ad(r[i+j],ml(x,y))
  return tf(r)
 def md(a,m):
  a,dm,il=LS(a),LN(m)-1,iv(m[-1])
  while LN(a)-1>=dm:
   c,sh=ml(a[-1],il),LN(a)-1-dm
   if c!=z:
    for i in RN(dm+1):a[sh+i]=sb(a[sh+i],ml(c,m[i]))
   a.pop()
  return tf(a)
 rs,b,e=[o],md([z,o],pl),qe
 while e:
  if e&1:rs=md(pm(rs,b),pl)
  b=md(pm(b,b),pl);e>>=1
 if LN(rs)<2:rs+=[z]*(2-LN(rs))
 rs[1]=sb(rs[1],o);tf(rs)
 if rs==[z]:return LN(pl)-1
 a,c=LS(pl),rs
 while any(x!=z for x in c):
  a,c=c,md(a,c)
  if c==[z]:break
 return LN(a)-1
def _fast_murphy(Pm=1500000):
 ps=_primes(Pm);af=ag=af15=ag15=0.;ck={};l2=L2
 for p in ps:
  er,lp=1./(p-1.),l2(p)
  if p<=300:
   if p==2:ef,eg=0.,8./17.*(SM(3.*k*15./16.**k for k in RN(1,6))+16./16.**5)
   elif p==3:ef,eg=12./82.,7./82.
   elif p==5:ef,eg=0.,150./624.
   else:
    cu,fd=p%15,1
    while cu!=1:cu=(cu*p)%15;fd+=1
    gc,qs=8//fd,p**fd;pf,xm=FL(qs),XB%p;rg=1 if xm==0 else((3 if PW(xm,(qs-1)//3,p)==1 else 0)if(qs-1)%3==0 else 1);eg=fd*gc*rg*pf/(pf*pf-1.);op=(0,1,lambda a,b:(a+b)%p,lambda a,b:(a-b)%p,lambda a,b:(a*b)%p,lambda a:PW(a,-1,p))
    if fd==1:w1=next(c for g in RN(2,p)if(c:=PW(g,(p-1)//3,p))!=1);ef=SM(4*_prc([-w%p,0,0,w*w%p,0,0,p-1,0,0,1],p,*op)for w in(w1,(w1*w1)%p))*pf/(pf*pf-1.)
    elif fd==2 or(fd==4 and p<=13):
     if p%3==1:w1=next(c for g in RN(2,p)if(c:=PW(g,(p-1)//3,p))!=1);rf=SM((4//fd)*_prc([-w%p,0,0,w*w%p,0,0,p-1,0,0,1],qs,*op)for w in(w1,(w1*w1)%p))
     else:fi=lambda a:(lambda iv:(((a[0]-a[1])*iv)%p,((-a[1])*iv)%p))(PW((a[0]*a[0]-a[0]*a[1]+a[1]*a[1])%p,-1,p));rf=gc*_prc([(0,p-1),(0,0),(0,0),(p-1,p-1),(0,0),(0,0),(p-1,0),(0,0),(0,0),(1,0)],qs,(0,0),(1,0),lambda a,b:((a[0]+b[0])%p,(a[1]+b[1])%p),lambda a,b:((a[0]-b[0])%p,(a[1]-b[1])%p),lambda a,b:((a[0]*b[0]-a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0]-a[1]*b[1])%p),fi)
     ef=fd*rf*pf/(pf*pf-1.)
    else:ef=8./pf
   af+=(er-ef)*lp;ag+=(er-eg)*lp
  else:
   r15=p%15
   if r15==1:pf,xm,e3=FL(p),XB%p,(p-1)//3;rg=1 if xm==0 else(3 if PW(xm,e3,p)==1 else 0);eg=8.*rg*pf/(pf*pf-1.);w1=next(c for g in RN(2,p)if(c:=PW(g,e3,p))!=1);ef=4.*(_fcr(w1,p,e3)+_fcr((w1*w1)%p,p,e3))*pf/(pf*pf-1.);af+=(er-ef)*lp;ag+=(er-eg)*lp
   elif r15==4:qs=FL(p*p);af+=(er-8./qs)*lp;ag+=(er-(24. if PW(XB%p,(p-1)//3,p)==1 else 0.)/qs)*lp
   elif r15 in(11,14):qs=FL(p*p);af+=(er-8./qs)*lp;ag+=(er-24./qs)*lp
   elif p<=2000:qs=FL(p*p)**2;af+=(er-8./qs)*lp;ag+=(er-(24. if(r15 in(2,8)or PW(XB%p,(p-1)//3,p)==1)else 0.)/qs)*lp
   else:d=er*lp;af+=d;ag+=d
  if p<=1500000:af15,ag15=af,ag
  if p in(4999,49999,199999,499979,999983,1499977,1999993):ck[str(p)]=[RD(af,5),RD(ag,5)]
 ck["1500000"]=[RD(af15,5),RD(ag15,5)];return af,ag,ck
def _cof_ex(Bf,ph,mnf,lf,af):
 dx=(Bf-ph)/80.;nb=331;cu=[1.]+[0.]*330;tp=ep=0.;ne=mnf+af
 for k in RN(1,11):
  nx=[0.]*nb
  for m in RN(80):
   x=ph+(m+.5)*dx;w,sh=dx/(x*k),I(RD(x/2.))
   for b in RN(nb-sh):
    if cu[b]:nx[b+sh]+=cu[b]*w
  cu=nx
  for b in RN(nb):
   if(sb:=b*2.)>ne:break
   cv=cu[b]*D0((ne-sb)/ph);tp+=cv;ep+=cv if sb>lf else 0.
 tp+=D0(ne/ph);return ep/tp,tp/D0(ne/Bf)
def chk_mrp(pmax=1500000):
 af,ag,ck=_fast_murphy(pmax);b,Ra,R0,q,Bf,Bg=R_["b"],R_["Ra"],R_["R0"],R_["qa"],R_["Bf"],R_["Bg"];afu,agu=ck["1500000"];sf=AB(L2(D0((542.17+afu+.5)/Bf)/D0((542.17+afu-.5)/Bf)));sg=AB(L2(D0((670.13+agu+.5)/Bg)/D0((670.13+agu-.5)/Bg)));It=(LG(1500000.)+1.)/(1500000.*(LG(2.)**2));sdr,sdh=SQ((sf*sf*12.+sg*sg*16.)*It),SQ((sf*sf*19.+sg*sg*23.+14.*sf*sg)*It);assert 3.*sdh<AB(R_["cush"]);dnf=(72./q)*L2((Ra**q+R0**q)/4.);dng=8.*L2X+(24./q)*L2((R0**q+(Ra**q)/(AX**(q/3.)))/4.)
 assert dnf<656. and dng<688. and 6*656<=4096 and 5*688<=4096 and 10*408<=4096 and I(F_(656./b["qf"]))-127==45 and 83+I(F_(656./b["qf"]))<=255 and all((128-MX(1,I(F_((n+1.)/b["qg"]))-127))+I(F_(n/b["qg"]))<=255 for n in RN(610,689))
 w0=1./(1.+R_["rbr"]**8);p0,bn0=_cof_ex(Bf,b["pfhi0"],542.17,R_["lf"],afu);p1,bn1=_cof_ex(Bf,b["pfhi1"],542.17,R_["lf"],afu);assert AB(bn0-1.)<.01 and AB(bn1-1.)<.01;pex=w0*p0+(1.-w0)*p1;cpad=L2(1./(1.-pex))
 return{"ok":T_,"pmax":pmax,"af_1500k":afu,"ag_1500k":agu,"af_end":RD(af,5),"ag_end":RD(ag,5),"ck":ck,"tail_sd_status":"H1_random_Frobenius_NOT_uncond","tail_1500k_sd_root_h1_heur":RD(sdr,6),"tail_1500k_sd_cov_h1_heur":RD(sdh,6),"tail_1500k_3sd_cov_h1_heur":RD(3.*sdh,6),"cush_b":R_["cush"],"det_nf":RD(dnf,6),"det_ng":RD(dng,6),"p_ex_f":RD(pex,6),"cof_pad":RD(cpad,6),"bnorm":[RD(bn0,5),RD(bn1,5)]}
_a,_b,_c,_d,_e,_f,_g=.48121182506,.687654016088,-.206442191029,.39628284719,-.87749467225,1.083936863278,.189840656161
VH=((_a,_a,-_a,-_a),(_b,_c,_d,_e),(-_c,-_b,-_e,-_d),(-_e,-_d,-_b,-_c),(_d,_e,_c,_b),(_f,-_f,_g,-_g),(_g,-_g,-_f,_f))
VHP=tuple((v,.5*SM(c*c for c in v))for v in VH)
SOB=((1,0,(1,)),(2,1,(1,1)),(3,1,(1,3,1)),(3,2,(1,1,1)),(4,1,(1,1,3,3)),(4,4,(1,3,5,13)),(5,2,(1,1,5,5,17)),(5,4,(1,1,5,5,5)))
def _sobol9(n):
 B=30;V=[[1<<(B-1-k)for k in RN(B)]]
 for s,a,m in SOB:
  mv=LS(m)
  for k in RN(s,B):
   v=mv[k-s]^(mv[k-s]<<s)
   for j in RN(1,s):
    if(a>>(s-1-j))&1:v^=(mv[k-j]<<j)
   mv+=[v]
  V+=[[mv[k]<<(B-1-k)for k in RN(B)]]
 pts,cu=[],[0]*9
 for i in RN(n):
  pts+=[tuple(cu)];c=((i+1)&-(i+1)).bit_length()-1
  for d in RN(9):cu[d]^=V[d][c]
 return pts
def _pdp(Bb,plo,phi,qs,ab,bmin,bst,nb,ns_=.5,pr_=1e-17,TM=260,km=700.):
 KM=I(RD(km/ns_));mo=[0.]*(KM+1)
 for ai,bi in((.561459483567,plo),(phi,Bb)):
  if bi<=ai:continue
  ns=MX(20,I((bi-ai)*8./ns_));dx=(bi-ai)/ns
  for i in RN(ns):
   x0=ai+i*dx
   if 1<=(kb:=I(RD((x0+.5*dx)/ns_)))<=KM:mo[kb]+=LG((x0+dx)/x0)
 nzo=[(j,j*mo[j])for j in RN(1,KM+1)if mo[j]>0.];Po=[0.]*(KM+1);Po[0]=EX(-SM(mo[1:]))
 for k in RN(1,KM+1):Po[k]=SM(jm*Po[k-j]for j,jm in nzo if j<=k)/k
 mi,m1={},[0.]*(TM+1);ns=MX(80,I((phi-plo)*16./ns_));dx=(phi-plo)/ns
 for i in RN(ns):
  x0=plo+i*dx;xm=x0+.5*dx;w,kb,tk=LG((x0+dx)/x0),I(RD(xm/ns_)),I(F_(xm/qs))
  if kb>=1 and 1<=tk<=TM:mi[(kb,tk)]=mi.get((kb,tk),0.)+w;m1[tk]+=w
 li=SM(m1[1:]);nz1=[(j,j*m1[j])for j in RN(1,TM+1)if m1[j]>0.];P1=[0.]*(TM+1);P1[0]=EX(-li)
 for k in RN(1,TM+1):P1[k]=SM(jm*P1[k-j]for j,jm in nz1 if j<=k)/k
 sv,rn=[0.]*(TM+1),0.
 for t in RN(TM,-1,-1):rn+=P1[t];sv[t]=rn
 its=SR([(dk,dt,dk*v)for(dk,dt),v in mi.items()]);P2=[N_]*(KM+1);P2[0]=[(0,EX(-li))];TM1=TM+1
 for k in RN(its[0][0],KM+1):
  ik,ds,tc=1./k,[0.]*TM1,F0
  for dk,dt,dv in its:
   if dk>k:break
   if(pr:=P2[k-dk])is N_:continue
   sd,lm=dv*ik,TM-dt
   for tt,p0 in pr:
    if tt>lm:break
    ds[tt+dt]+=sd*p0;tc=T_
  if tc and(rp:=[(tt,pv)for tt,pv in EN(ds)if pv>pr_]):P2[k]=rp
 nzk=[(k,P2[k])for k in RN(KM+1)if P2[k]is not N_];rt=[]
 for bi in RN(nb):
  kN=MX(1,MN(KM,I(RD((bmin+(bi+.5)*bst+ab)/ns_))));sm=[0.]*TM1
  for k,ki in nzk:
   if k>kN:break
   if(po:=Po[kN-k])<=1e-300:continue
   for tt,pt in ki:sm[tt]+=pt*po
  tot,tl,rw=SM(sm),0.,[0.]*129
  for tt in RN(TM,-1,-1):
   tl+=sm[tt]
   if tt<=128:rw[tt]=(tl/tot)if tot>0. else 0.
  rt+=[rw]
 return sv,rt
def chk_rqmc(af,ag,nS=24,mS=1024,nch=1440,fD=F0):
 b,Bf,Bg,qa,rbr,Ra,R0=R_["b"],R_["Bf"],R_["Bg"],R_["qa"],R_["rbr"],R_["Ra"],R_["R0"];pglo,pflo,qg,qf=b["pglo"],b["pflo"],b["qg"],b["qf"];w0=1./(1.+rbr**8);w1=1.-w0;G0,GS,GB,F0_,FS,FB,N0,NS,NB=610,2,40,510,4,25,440,2,100;gL,rL,fL,kL,gc,rc,fc,kc,gf,rF,ff,kf=([]for _ in RN(12));tmf=320 if fD else 260
 for si,(pgh,pfh)in EN(((b["pghi0"],b["pfhi0"]),(b["pghi1"],b["pfhi1"]))):
  sg,rg=_pdp(Bg,pglo,pgh,qg,ag,G0,GS,GB,.5,1e-17,260,695.);sf,rf=_pdp(Bf,pflo,pfh,qf,af,N0,NS,NB,.5,1e-17,260,650.);gL+=[sg];rL+=[rg];fL+=[sf];kL+=[rf]
  if fD or si==0:sgc,rgc=_pdp(Bg,pglo,pgh,qg,ag,G0,GS,GB,1.,1e-13,260,695.);sfc,rfc=_pdp(Bf,pflo,pfh,qf,af,N0,NS,NB,1.,1e-13,260,650.);sgf,rgf=_pdp(Bg,pglo,pgh,qg,ag,G0,GS,GB,.25,1e-20,tmf,695.);sff,rff=_pdp(Bf,pflo,pfh,qf,af,N0,NS,NB,.25,1e-20,tmf,650.);gc+=[sgc];rc+=[rgc];fc+=[sfc];kc+=[rfc];gf+=[sgf];rF+=[rgf];ff+=[sff];kf+=[rff]
  else:gc+=[sg];rc+=[rg];fc+=[sf];kc+=[rf];gf+=[sg];rF+=[rg];ff+=[sf];kf+=[rf]
 dk2,dk8=EDk(30.,200),EDk(30.,800);ug=[1.5+.25*k for k in RN(40)];d200=MX(AB(dk2.rho(u)/dk8.rho(u)-1.)for u in ug);d400=MX(AB(D0(u)/dk8.rho(u)-1.)for u in ug);d400l=MX(AB(L2(D0(u)/dk8.rho(u)))for u in ug);al=2./qa;Cd=EX(lgamma(4.*al)-4.*lgamma(al));atb,ata,k0,l0,k1,l1=3.2,4.2,3.2,5.5,2.,2.5;ik0,il0,ik1,il1,itb,ita=1./k0,1./l0,1./k1,1./l1,1./atb,1./ata;iqa,il2,tp=1./qa,1./LG(2.),2.*pi;iq3,al4=3.*iqa,4.*al;Ap=[(a**6,a**12,a**18)for a in((Ra/(4.**iqa))*((.25 if k==0 else .75)**.125)for k in RN(2))]
 r0c,r1c,S32,g0c=(Ra*(.5**.125))**3,Ra**3,.8660254037844386,8.*L2X;lg,cs_,sn_,ex_,XF=LG,cos,sin,EX,FL(XB);Dr,D2r,D8r=DK.rho,dk2.rho,dk8.rho
 def da(rn):
  while T_:
   u0,u1,u2=1.36*rn()-.68,1.36*rn()-.68,1.36*rn()-.68;u3=-(u0+u1+u2)
   if AB(u3)>.68 or any(AB(u0*v[0]+u1*v[1]+u2*v[2]+u3*v[3])>h for v,h in VHP):continue
   e0,e1,e2,e3=ex_(qa*u0),ex_(qa*u1),ex_(qa*u2),ex_(qa*u3);es=e0+e1+e2+e3
   if rn()<=(4./es)**al4:ie=1./es;return(e0*ie)**iq3,(e1*ie)**iq3,(e2*ie)**iq3,(e3*ie)**iq3
 def se_(ra,rb,a,b_,cs,A,lbp):
  rb3=rb**3;rb6=rb3*rb3;ng0=g0c+24.*lg(rb)*il2+lbp;fp=gp=pp=1.;A6,A12,A18=A
  for m in RN(4):bm=b_[m];r6=rb6*bm*bm;r12=r6*r6;pp*=A18+A12*r6+A6*r12+r12*r6;ar=ra*a[m]*cs[2*m];ai=ra*a[m]*cs[2*m+1];bt=rb3*bm;gr=ar-XF*bt;gp*=gr*gr+ai*ai+1e-300;a2r=ar*ar-ai*ai;a2i=2.*ar*ai;b2_=bt*bt;b3_=b2_*bt;fr=a2r*ar-a2i*ai-a2r*bt+(-.5*ar+S32*ai)*b2_+.5*b3_;fi=a2r*ai+a2i*ar-a2i*bt+(-.5*ai-S32*ar)*b2_-S32*b3_;fp*=fr*fr+fi*fi+1e-300
  return lg(fp)*il2,lg(gp)*il2,ng0,lg(pp)*il2
 def ev2(a,u,cs):
  tb=u[0]**itb;wtb=(8./atb)*(tb**(8.-atb));ta0=u[1]**ita;wta=(8./ata)*(ta0**(8.-ata));t0=(1.-(1.-u[2])**il0)**ik0;t1=(1.-(1.-u[3])**il1)**ik1;t2=(1.-(1.-u[4])**il1)**ik1;p0=k0*l0*(t0**(k0-1.))*((1.-t0**k0)**(l0-1.));p1=k1*l1*(t1**(k1-1.))*((1.-t1**k1)**(l1-1.));p2=k1*l1*(t2**(k1-1.))*((1.-t2**k1)**(l1-1.));ub=(t0*t1,t0*(1.-t1),(1.-t0)*t2,(1.-t0)*(1.-t2));up=ub[0]*ub[1]*ub[2]*ub[3];wb=wtb*Cd*(up**(al-1.))*(t0*(1.-t0))/(p0*p1*p2);b_=(ub[0]**iq3,ub[1]**iq3,ub[2]**iq3,ub[3]**iq3);lbp=6.*iqa*lg(up)*il2
  return wb*(w0*wta+w1),((0,w0*wta*wb,*se_(r0c*(ta0**3),R0*tb,a,b_,cs,Ap[0],lbp)),(1,w1*wb,*se_(r1c*((.5+.5*u[1])**.375),R0*rbr*tb,a,b_,cs,Ap[1],lbp)))
 rf_=RR(20261127);rnf=rf_.random;v2=[[[0.]*FB for _ in RN(GB)]for _ in RN(2)];s3=[[[[0.]*NB for _ in RN(FB)]for _ in RN(GB)]for _ in RN(2)];nF,tsf=4000,0.
 for i in RN(nF):
  if i%4==0:a_=da(rnf)
  h=[tp*rnf()for _ in RN(4)];cs=(cs_(h[0]),sn_(h[0]),cs_(h[1]),sn_(h[1]),cs_(h[2]),sn_(h[2]),cs_(h[3]),sn_(h[3]));_,sr=ev2(a_,(rnf(),rnf(),rnf(),rnf(),rnf()),cs)
  for sh,wi,nf,ng,ng0,nfp in sr:ig=MX(0,MN(GB-1,I((ng0-G0)//GS)));icf=MX(0,MN(FB-1,I((nfp-F0_)//FS)));iff=MX(0,MN(NB-1,I((nf-N0)//NS)));wt=wi*Dr((nf+af)/Bf)*Dr((ng+ag)/Bg);v2[sh][ig][icf]+=wi;s3[sh][ig][icf][iff]+=wt;tsf+=wt
 ivn,ivt,tfm,lp=1./nF,1./tsf,MX(0,I(F_(656./qf))-127),R_["lam"];Tg,Tf=[[[128]*FB for _ in RN(GB)]for _ in RN(2)],[[[125]*FB for _ in RN(GB)]for _ in RN(2)]
 def hull(rw,sfs):
  h=[]
  for tf in RN(125,tfm-1,-1):
   x,y=sfs[tf],rw[tf]
   while LN(h)>=2 and(h[-1][2]-h[-2][2])*(x-h[-1][1])<=(y-h[-1][2])*(h[-1][1]-h[-2][1]):h.pop()
   h+=[(tf,x,y)]
  return h
 for sh in RN(2):
  rgs,rfs,sgs,sfs=rL[sh],kL[sh],gL[sh],fL[sh]
  for ig in RN(GB):
   ngu=G0+(ig+1)*GS+1.;tgm=MX(0,MN(128,I(F_(ngu/qg))-127));vw1=SM(v2[sh][ig])*ivn;nz1=[(iff,sv)for iff in RN(NB)if(sv:=SM(s3[sh][ig][icf][iff]for icf in RN(FB))*ivt)>0.];h1=hull([SM(sv*rfs[iff][tf]for iff,sv in nz1)if(tf>=tfm and vw1>0.)else 0. for tf in RN(126)],sfs)if vw1>0. else[]
   for icf in RN(FB):
    vw=v2[sh][ig][icf]*ivn
    if vw<5e-5:
     if vw1<=0.:continue
     vwu,hc=vw1,h1
    else:vwu=vw;nzf=[(iff,s3[sh][ig][icf][iff]*ivt)for iff in RN(NB)if s3[sh][ig][icf][iff]>0.];hc=hull([SM(sv*rfs[iff][tf]for iff,sv in nzf)if tf>=tfm else 0. for tf in RN(126)],sfs)
    bo,btg,btf=-1e300,128,125
    for tg in RN(tgm,129):
     rgv,sgt=rgs[ig][tg],lp*vwu*sgs[tg]
     for tf,sfx,rfy in hc:
      o=rgv*rfy-sgt*sfx
      if o>bo:bo,btg,btf=o,tg,tf
    Tg[sh][ig][icf],Tf[sh][ig][icf]=btg,btf;assert(128-btg)+I(F_(ngu/qg))<=255 and(128-btf)+I(F_(656./qf))<=255
 NB_=nS;sob=_sobol9(mS);re_=RR(20261128);rne=re_.random;s30=1./(1<<30);bk,br,bg,b1,bc,bC,bf,bR=([]for _ in RN(8));mfs=mgs=swn=nfm=ngm=0.;ncm=16 if fD else 2;iBf,iBg=1./Bf,1./Bg
 for si in RN(NB_):
  qs=[re_.getrandbits(30)for _ in RN(9)];js=[rne()*s30 for _ in RN(9)];swk=swr=ssg=ss1=swc=skc=skf=srf=0.;dc=(si<ncm)
  for pt in sob:
   a_=da(rne);u=[(pt[d]^qs[d])*s30+js[d]for d in RN(9)];c0,s0,c1,s1,c2,s2,c3,s3_=cs_(tp*u[5]),sn_(tp*u[5]),cs_(tp*u[6]),sn_(tp*u[6]),cs_(tp*u[7]),sn_(tp*u[7]),cs_(tp*u[8]),sn_(tp*u[8])
   for uu,cs in((u,(c0,s0,c1,s1,c2,s2,c3,s3_)),((1.-u[0],1.-u[1],1.-u[2],1.-u[3],1.-u[4]),(-c0,-s0,-c1,-s1,-c2,-s2,-c3,-s3_))):
    wcp,sr=ev2(a_,uu,cs);swc+=.5*wcp
    for sh,ws_,nf,ng,ng0,nfp in sr:
     if nf>nfm:nfm=nf
     if ng>ngm:ngm=ng
     ig=MX(0,MN(GB-1,I((ng0-G0)//GS)));icf=MX(0,MN(FB-1,I((nfp-F0_)//FS)));iff=MX(0,MN(NB-1,I((nf-N0)//NS)));iga=MX(0,MN(GB-1,I((ng-G0)//GS)));wi=.5*ws_;uf,ug_=(nf+af)*iBf,(ng+ag)*iBg;wt=wi*Dr(uf)*Dr(ug_);tg,tf=Tg[sh][ig][icf],Tf[sh][ig][icf];wk=wt*rL[sh][iga][tg]*kL[sh][iff][tf];sp=gL[sh][tg];swr+=wt;swk+=wk;ssg+=wi*sp;ss1+=wi*sp*fL[sh][tf];mfs+=wk*nf;mgs+=wk*ng;swn+=wk
     if dc:skc+=wi*D2r(uf)*D2r(ug_)*rc[sh][iga][tg]*kc[sh][iff][tf];wtf=wi*D8r(uf)*D8r(ug_);srf+=wtf;skf+=wtf*rF[sh][iga][tg]*kf[sh][iff][tf]
  bk+=[swk/mS];br+=[swr/mS];bg+=[ssg/mS];b1+=[ss1/mS];bc+=[swc/mS]
  if dc:bC+=[skc/mS];bf+=[skf/mS];bR+=[srf/mS]
 mwk,mwr,msg,ms1,ewi=SM(bk)/NB_,SM(br)/NB_,SM(bg)/NB_,SM(b1)/NB_,SM(bc)/NB_;seV=R_["sl"]*(SQ(SM((x-mwk)**2 for x in bk)/(NB_*(NB_-1)))/mwk)*il2;volu=(pi**4)*(gamma(1.+2./qa)**4)/gamma(1.+8./qa);Vkb=L2(VA0*(volu*(R0**8)/(70.3125*3.))*(1.+rbr**8));U=U_ub(Bf,Bg);Vkl=Vkb+R_["sl"]*(L2(R_["eta"]*U)-(Vkb+L2(mwk*.97843406*(2.**R_["cush"]))));retL,l2sg,l2s1=mwk/mwr,L2(msg),L2(ms1);l2sf=l2s1-l2sg;mc,mr=SM(bk[:ncm])/ncm,SM(br[:ncm])/ncm
 dcp=-R_["sl"]*L2(mc/(SM(bC)/ncm));dpf=-R_["sl"]*L2((SM(bf)/ncm)/mc);drp=(SM(bf)/SM(bR))-(mc/mr);dpad=R_["disc_pad"];assert AB(ewi-1.)<.005 and nfm<656. and ngm<688. and seV<.0032 and dpad>=2.5*AB(dpf);rc_,hq,agd,NC_B=RR(20260505),.5*qa,2./qa,nch//NB_;l4,lh,g4=([]for _ in RN(3))
 for _ in RN(NB_):
  a4=ah=ag_=0.
  for R,wsh in((R0,w0),(R0*rbr,w1)):
   Rq,sp,gp_,w4,wh=R**qa,0.,0.,0.,0.;ns,tl=MX(24,I(RD(NC_B*wsh))),I(CE(R*1.15))+20
   for _ in RN(ns):
    gs=[rc_.gammavariate(agd,1.)for _ in RN(4)];igs=(Rq*(MX(rc_.random(),1e-16)**(qa/8.)))/SM(gs);zr,zi2=[],[]
    for gm_ in gs:rm,th=(gm_*igs)**iqa,tp*rc_.random();zr+=[rm*cs_(th)];zi2+=[(rm*sn_(th))**2]
    z0r,z1r,z2r,z3r=zr;z0i,z1i,z2i,z3i=zi2;cmn,cmx=10000,-10000
    for tc in RN(-tl,tl+1,16):
     s=((z0r+tc)**2+z0i)**hq+((z1r+tc)**2+z1i)**hq
     if s<=Rq and(s:=s+((z2r+tc)**2+z2i)**hq)<=Rq and s+((z3r+tc)**2+z3i)**hq<=Rq:cmn,cmx=MN(cmn,tc),MX(cmx,tc)
    tlo,thi=(MX(-tl,cmn-32),MN(tl,cmx+32))if cmn<=cmx else(-256,256);tmn,tmx,ic=10000,-10000,0
    for t in RN(tlo,thi+1):
     s=((z0r+t)**2+z0i)**hq
     if s<=Rq and(s:=s+((z1r+t)**2+z1i)**hq)<=Rq and(s:=s+((z2r+t)**2+z2i)**hq)<=Rq and s+((z3r+t)**2+z3i)**hq<=Rq:tmn,tmx,ic=MN(tmn,t),MX(tmx,t),ic+1
    if ic:sp_=tmx-tmn+1;w=1./ic;sp+=w*sp_;gp_+=w*(sp_-ic);w4+=w*((sp_+247)//248);wh+=w*((sp_+123)//124)
   a4+=wsh*(sp/w4);ah+=wsh*(sp/wh);ag_+=wsh*(gp_/sp)
  l4+=[a4];lh+=[ah];g4+=[ag_]
 m4,mh,mg_=SM(l4)/NB_,SM(lh)/NB_,SM(g4)/NB_;se4,seh,seg=(SQ(SM((x-m)**2 for x in arr)/(NB_*(NB_-1)))for arr,m in((l4,m4),(lh,mh),(g4,mg_)));c4,ch,cg=m4-3.*se4,mh-3.*seh,mg_+3.*seg;gpad=L2(1./(1.-cg));assert 95.<c4<120. and 68.<ch<88. and 0.<cg<.008
 return{"ok":T_,"sampler":"antithetic_9D_Sobol_RQMC","n_fit":nF,"n_shifts":NB_,"n_eval_rqmc":NB_*mS*2,"n_chords":nch,"E_w_is":RD(ewi,5),"ret_live":RD(retL,6),"l2_sg":RD(l2sg,6),"l2_sf":RD(l2sf,6),"l2_s1":RD(l2s1,6),"mean_nf":RD(mfs/swn,3),"mean_ng":RD(mgs/swn,3),"obs_max":[RD(nfm,2),RD(ngm,2)],"Vk_base":RD(Vkb,6),"Vk_live":RD(Vkl,6),"se_Vk":RD(seV,6),"pad_3se_Vk":RD(3.*seV,6),"L248_live":RD(m4,3),"se_L248":RD(se4,3),"L248_ch":RD(c4,3),"L124_live":RD(mh,3),"se_L124":RD(seh,3),"L124_ch":RD(ch,3),"gap4_live":RD(mg_,6),"se_gap4":RD(seg,6),"gap_ch":RD(cg,6),"gap_pad":RD(gpad,6),"discretization_convergence":{"ok":T_,"primary_grid":"EDk(400)+Panjer(0.5b,1e-17,TM260)","fine_grid":f"EDk(800)+Panjer(0.25b,1e-20,TM{tmf})","fine_shells":"Shell0+Shell1"if fD else f"Shell0_0.25b_1e-20_TM{tmf}","crn_eval_pts":ncm*mS*2,"dickman_N200_vs_800_max_rel":FL(f"{d200:.3e}"),"dickman_N400_vs_800_max_rel":FL(f"{d400:.3e}"),"dickman_N400_vs_800_max_l2":FL(f"{d400l:.3e}"),"dVk_coarse1b1e13_to_prim05b1e17":RD(dcp,6),"dVk_prim_to_fine025b1e20":RD(dpf,6),"dret_prim_to_fine":RD(drp,6),"disc_pad_b":dpad,"pad_over_shift_x":RD(dpad/MX(1e-9,AB(dpf)),2),"scope":"DP_Dickman_only_H1_exploratory_v1"}}
def _ip64(n):
 d,s=n-1,0
 while d%2==0:d//=2;s+=1
 for a in(2,3,5,7,11,13,17,19,23,29,31,37):
  if(x:=PW(a,d,n))in(1,n-1):continue
  for _ in RN(s-1):
   if(x:=(x*x)%n)==n-1:break
  else:return F0
 return T_
_sh=lambda a,w,wb,p:(a*w-((a*wb)>>64)*p)&(1<<64)-1
def _ntt_r4(vc,p,g):
 N,a,p2,p4,l4,t=LN(vc),[x%p for x in vc],p<<1,p<<2,0,LN(vc)
 while t>1:t>>=2;l4+=1
 for i in RN(N):
  rv,x=0,i
  for _ in RN(l4):rv=(rv<<2)|(x&3);x>>=2
  if i<rv:a[i],a[rv]=a[rv],a[i]
 J=PW(g,(p-1)//4,p);Jp,m=(J<<64)//p,1
 while m<N:
  wm,w1,W=PW(g,(p-1)//(4*m),p),1,[]
  for _ in RN(m):w2=(w1*w1)%p;w3=(w2*w1)%p;W+=[(w1,(w1<<64)//p,w2,(w2<<64)//p,w3,(w3<<64)//p)];w1=(w1*wm)%p
  for k in RN(0,N,4*m):
   for j in RN(m):w1,p1,w2,q2,w3,q3=W[j];x0=a[k+j]-(p2 if a[k+j]>=p2 else 0);x1,x2,x3=_sh(a[k+m+j],w1,p1,p),_sh(a[k+2*m+j],w2,q2,p),_sh(a[k+3*m+j],w3,q3,p);t0,t1=(x0+x2)-(p2 if x0+x2>=p2 else 0),(x0+p2-x2)-(p2 if x0+p2-x2>=p2 else 0);t2,t3=(x1+x3)-(p2 if x1+x3>=p2 else 0),_sh(x1+p2-x3,J,Jp,p);a[k+j],a[k+2*m+j],a[k+m+j],a[k+3*m+j]=t0+t2,t0+p2-t2,t1+t3,t1+p2-t3;assert MX(a[k+j],a[k+m+j],a[k+2*m+j],a[k+3*m+j])<p4
  m<<=2
 return[x%p for x in a]
def _ntt1(vc,p,g,iv=F0):
 n,a,p2,j=LN(vc),[x%p for x in vc],2*p,0
 for i in RN(1,n):
  b=n>>1
  while j&b:j^=b;b>>=1
  j^=b
  if i<j:a[i],a[j]=a[j],a[i]
 L=2
 while L<=n:
  hf,wl=L>>1,PW(g,(p-1)//L,p)
  if iv:wl=PW(wl,p-2,p)
  ws,wp,cu=[1]*hf,[0]*hf,1
  for k in RN(hf):ws[k],wp[k]=cu,(cu<<64)//p;cu=(cu*wl)%p
  for i in RN(0,n,L):
   for k in RN(hf):u=a[i+k]-(p2 if a[i+k]>=p2 else 0);t=_sh(a[i+k+hf],ws[k],wp[k],p);a[i+k],a[i+k+hf]=u+t,u+p2-t
  L<<=1
 return[(x*PW(n,p-2,p))%p for x in a]if iv else[x%p for x in a]
def _mul2d(A,B,pgs):
 m5=(1<<58)-1
 def t5(x):
  r=[]
  while x:r+=[x&m5];x>>=58
  return r or[0]
 la,lb=t5(A),t5(B);k=MX(1,(MX(LN(la),LN(lb))-1).bit_length());k1=(k+1)//2;D1,D2=1<<k1,1<<(k-k1);N1,N2=2*D1,2*D2;la+=[0]*(D1*D2-LN(la));lb+=[0]*(D1*D2-LN(lb));cm=[]
 for p,g in pgs:
  gA=[_ntt1(la[r*D1:(r+1)*D1]+[0]*D1,p,g)for r in RN(D2)]+[[0]*N1 for _ in RN(D2)];gB=[_ntt1(lb[r*D1:(r+1)*D1]+[0]*D1,p,g)for r in RN(D2)]+[[0]*N1 for _ in RN(D2)];gC=[[0]*N1 for _ in RN(N2)]
  for c in RN(N1):
   cA,cB=_ntt1([gA[r][c]for r in RN(N2)],p,g),_ntt1([gB[r][c]for r in RN(N2)],p,g);cC=_ntt1([(cA[r]*cB[r])%p for r in RN(N2)],p,g,T_)
   for r in RN(N2):gC[r][c]=cC[r]
  cm+=[[_ntt1(gC[r],p,g,T_)for r in RN(N2)]]
 p1,p2,p3=pgs[0][0],pgs[1][0],pgs[2][0];i12,i123,p12=PW(p1%p2,p2-2,p2),PW((p1*p2)%p3,p3-2,p3),p1*p2;acc=[0]*(D1*(N2+1)+4)
 for r in RN(N2):
  bs=r*D1
  for c in RN(N1):a1,a2,a3=cm[0][r][c],cm[1][r][c],cm[2][r][c];x12=a1+p1*(((a2-a1)*i12)%p2);acc[bs+c]+=x12+p12*(((a3-x12)*i123)%p3)
 rs=cy=0
 for i,v in EN(acc):s=v+cy;rs|=(s&m5)<<(58*i);cy=s>>58
 return rs|(cy<<(58*LN(acc)))
def chk_ntt():
 s=R_["s1b"];p1,p2,p3,g1,g2,g3=s["p1"],s["p2"],s["p3"],s["g1"],s["g2"],s["g3"];l2P=L2(p1*p2*p3);h1,h2=l2P-181.,l2P-(L2(CE((2.**67.77)/58.))+116.);assert all(_ip64(p)and(p-1)%(1<<54)==0 and 8*p<(1<<64)and(1<<58)<p for p in(p1,p2,p3))and all(PW(g1,(p1-1)//q,p1)!=1 for q in(2,127))and all(PW(g2,(p2-1)//q,p2)!=1 for q in(2,3,19))and all(PW(g3,(p3-1)//q,p3)!=1 for q in(2,3))and h1>1.57 and h2>4.66;rng,pgs=RR(0x5400612D),[(p1,g1),(p2,g2),(p3,g3)]
 for p,g in pgs:
  vc=[rng.randrange(p)for _ in RN(64)];gt=_ntt_r4(vc,p,g);wN=PW(g,(p-1)//64,p)
  for k in RN(64):assert gt[k]==SM(vc[j]*PW(wN,j*k,p)for j in RN(64))%p
 sz=[512,1024,2048,4096]
 for nb in sz:A,B=rng.getrandbits(nb)|(1<<(nb-1)),rng.getrandbits(nb)|(1<<(nb-1));assert _mul2d(A,B,pgs)==A*B
 return{"ok":T_,"r4_bitexact":T_,"log2_P":RD(l2P,6),"hr_1d":RD(h1,4),"hr_2d":RD(h2,4),"mul_2d":sz}
def _h576(A,sm,v):
 r256,x2h,m64=RB<<256,0xAC45A4010001A402,(1<<64)-1;T=A[0]+(A[2]<<32)+((A[2]*x2h)<<64)+r256-AX*(A[1]+(A[3]<<32)+((A[3]*x2h)<<64))+SM(sm[k]*v[k]for k in RN(4));assert T>=0 and(t8:=T>>512)in(0,1);T-=t8*(r256<<1);assert 0<=T<(1<<512);Th=T>>256;t4,t5,t6,t7=Th&m64,(Th>>64)&m64,(Th>>128)&m64,(Th>>192)&m64;M=(0x42737A020C0D6393,0x65043EB4BE4BAD71,0x38B5DCB707E08ED3,0x355094EDFEDE377C);ph=((t4*M[3])+(t5*M[2]+((t5*M[3])<<64))+(t6*M[1]+((t6*M[2])<<64)+((t6*M[3])<<128))+(t7*M[0]+((t7*M[1])<<64)+((t7*M[2])<<128)+((t7*M[3])<<192)))>>64;s=(T-((Th<<1)+ph)*RB)&((1<<320)-1);assert 0<=s<6*RB;s2=s-(s>>255)*RB;assert 0<=s2<2*RB;return s2-RB if s2>=RB else s2
def chk_krn():
 p_=[((i&31)<<56)|(i*131+5)for i in RN(17)];b_=[0,0,1,4,4,8,17];t_=[SB0+j*32*SB0+p_[k]for j in RN(6)for k in RN(b_[j],b_[j+1])];assert LN(t_)==17 and t_[-1]==SB0+(5*32+16)*SB0+(16*131+5)and(900+200)%1009==91 and(91>>5)==2;Va,pb=1024,37;sl,sh=BA((pb+(i*17)%60)for i in RN(Va)),BA((i*29)%50 for i in RN(Va));sl[13],sh[13],sl[14],sh[14],sl[255],sh[255],sl[256],sh[256]=80,47,80,48,150,105,60,67;e_=[i for i in RN(Va)if sl[i]+sh[i]>=128];assert all(sl[i]+sh[i]<=255 for i in RN(Va));r_=BA(sl);sl[:]=bytes([pb])*Va;g_=[];u6=lambda b,o:I.from_bytes(b[o:o+8],"little");M8=0x8080808080808080
 for bk in RN(0,Va,128):
  ac=0
  for w in RN(16):ac|=(u6(r_,bk+8*w)+u6(sh,bk+8*w))&((1<<64)-1)
  if ac&M8:
   for w in RN(16):
    m=(u6(r_,bk+8*w)+u6(sh,bk+8*w))&M8
    while m:g_+=[bk+8*w+(((m&-m).bit_length()-1)>>3)];m&=m-1
  sh[bk:bk+128]=bytes(128)
 assert g_==e_ and not any(sh)and all(x==pb for x in sl);rng=RR(0x576B1750)
 for t in RN(60):A,sm,vv=([(1<<304)-1,0,(1<<304)-1,0],[RB-1]*4,[RB-1]*4)if t==0 else(([0,(1<<314)-1,0,(1<<314)-1],[0]*4,[0]*4)if t==1 else(([(1<<300)-1]*4,[RB-1]*4,[RB-1]*4)if t==2 else([rng.getrandbits(314)for _ in RN(4)],[rng.randrange(RB)for _ in RN(4)],[rng.randrange(RB)for _ in RN(4)])));assert _h576(A,sm,vv)==(A[0]+XB*A[1]+XB*XB*A[2]+(XB**3)*A[3]+SM(sm[k]*vv[k]for k in RN(4)))%RB
 Vv,r60=[rng.randrange(RB)for _ in RN(256)],RB<<60;z1=XB%RB;z2,z3=(z1*z1)%RB,(z1*z1*z1)%RB;te=c1=cu=0
 for _ in RN(60):
  Ap,Am,rr=[0]*4,[0]*4,0
  for ei in RN(60):
   col=rng.randrange(256);tp_=MONO[ei]if ei<12 else(MONO[rng.randrange(12)]if rng.random()<.70 else tuple((0 if rng.random()<.35 else(1 if rng.random()<.5 else-1)*(1 if rng.random()<.4 else rng.randrange(2,1<<29)))for _ in RN(4)));te+=1;rr=(rr+(tp_[0]+tp_[1]*z1+tp_[2]*z2+tp_[3]*z3)*Vv[col])%RB
   for k,ck in EN(tp_):
    if ck:c1+=I(AB(ck)==1);cu+=I(AB(ck)!=1);(Ap if ck>0 else Am)[k]+=AB(ck)*Vv[col]
  As=[Ap[k]+r60-Am[k]for k in RN(4)];assert all(0<=Ap[k]<(1<<293)and 0<=Am[k]<(1<<293)and 0<=As[k]<(1<<315)for k in RN(4));sm,vv=[rng.randrange(RB)for _ in RN(4)],[rng.randrange(RB)for _ in RN(4)];assert _h576(As,sm,vv)==(rr+SM(sm[j]*vv[j]for j in RN(4)))%RB
 return{"ok":T_,"chord_tx_hits":LN(t_),"swar_surv":LN(g_),"zero_carry_255":T_,"horner576":"60/60","z12_4tuple_spmv":"60/60","comp_per_ent":RD((c1+cu)/te,3)}
def _rref(mt,nc,md):
 pv={}
 for rw in mt:
  cu={c:t for c,v in rw.items()if(t:=v%md)}
  while cu:
   if(pc:=MN(cu))in pv:
    fc=cu[pc]
    for c,vv in pv[pc].items():
     if(nv:=(cu.get(c,0)-fc*vv)%md):cu[c]=nv
     else:cu.pop(c,N_)
   else:iv=PW(cu[pc],-1,md);pv[pc]={c:(v*iv)%md for c,v in cu.items()};break
 return pv
def chk_twr():
 p,r,m=257,43,16;r2,ep=r*r,r*r-1;assert p%3==2 and(p+1)%r==0 and m*m+1==p;zm=lambda x,y:(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]-x[1]*y[1]);zn=lambda x:x[0]*x[0]-x[0]*x[1]+x[1]*x[1];wm=lambda x,y,M:((x[0]*y[0]-x[1]*y[1])%M,(x[0]*y[1]+x[1]*y[0]-x[1]*y[1])%M);wi=lambda x,M:(lambda iv:(((x[0]-x[1])*iv)%M,((-x[1])*iv)%M))(PW(zn(x)%M,-1,M))
 def wg1(x,y):
  while y!=(0,0):ny=zn(y);nm=zm(x,(y[0]-y[1],-y[1]));q=((2*nm[0]+ny)//(2*ny),(2*nm[1]+ny)//(2*ny));t=zm(q,y);x,y=y,(x[0]-t[0],x[1]-t[1])
  return zn(x)==1
 def fm(x,y,M):a,b,c,d=wm(x[0],y[0],M),wm(x[1],y[1],M),wm(x[0],y[1],M),wm(x[1],y[0],M);return[((a[0]-b[0])%M,(a[1]-b[1])%M),((c[0]+d[0])%M,(c[1]+d[1])%M)]
 def fp(x,e,M):
  rs,bs=[(1,0),(0,0)],[(z[0]%M,z[1]%M)for z in x]
  while e:
   if e&1:rs=fm(rs,bs,M)
   bs=fm(bs,bs,M);e>>=1
  return rs
 assert fp([(3,2),(-1,5)],ep,r)==[(1,0),(0,0)]
 def smq(pw):assert(pw[0][0]-1)%r==0 and pw[0][1]%r==0 and pw[1][0]%r==0 and pw[1][1]%r==0;return[((pw[0][0]-1)//r)%r,(pw[0][1]//r)%r,(pw[1][0]//r)%r,(pw[1][1]//r)%r]
 smf=lambda a,b:smq(fp([(a[0]%r2,a[1]%r2),((-b[0])%r2,(-b[1])%r2)],ep,r2));assert smf((2,0),(-1,-2))[3]!=0;s1,s2=smf((3,1),(2,-1)),smf((5,-2),(1,3));assert smq(fp(fm([(3,1),(-2,1)],[(5,-2),(-1,-3)],r2),ep,r2))==[(s1[i]+s2[i])%r for i in RN(4)]and smq(fp(fp([(3,1),(-2,1)],r,r2),ep,r2))==[0,0,0,0];ps=_primes(42);sp3,frs,fri=({}for _ in RN(3))
 for q in ps:
  if q%3==1:w0=next(t for t in RN(1,q)if(t*t+t+1)%q==0);sp3[q]=(MN(w0,(q-1-w0)%q),MX(w0,(q-1-w0)%q))
  if q>3:
   if q%12==1:iq=next(t for t in RN(1,q)if(t*t+1)%q==0);frs[q]=MN(iq,q-iq)
   elif q%3==2:fri[q]=(MN(t:=next(t for t in RN(1,q)if(t*t+1)%q==0),q-t),0)if q%4==1 else((x0:=next(t for t in RN(1,q)if(3*t*t-1)%q==0)),(2*x0)%q)
 fa,un=SR(LS(frs)+LS(fri)),((1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1))
 def hv(ux,vx,wv,q,eq):
  if(ux+vx*wv)%q:return 0
  wl,qk=wv,q
  for _ in RN(eq-1):qk*=q;wl=(wl-(wl*wl+wl+1)*PW((2*wl+1)%qk,-1,qk))%qk
  ev,pk=0,q
  while ev<eq and(ux+vx*wl)%pk==0:ev+=1;pk*=q
  return ev
 def fg_(vK):
  rm,ux,vx,o=zn(vK),vK[0],vK[1],{}
  if rm==0:return N_
  for q in ps:
   if rm%q==0:
    eq=0
    while rm%q==0:rm//=q;eq+=1
    if q==3:continue
    if q in sp3:
     e0,e1=hv(ux,vx,sp3[q][0],q,eq),hv(ux,vx,sp3[q][1],q,eq)
     if e0+e1!=eq:return N_
     if e0!=e1:o[("g",q)]=e0-e1
    elif eq%2 or ux%(qe:=q**(eq//2))or vx%qe:return N_
  return o if rm==1 else N_
 def ff_(a,b,vK):
  rm,ux,vx,o=zn(vK),vK[0],vK[1],{}
  if rm==0:return N_
  for q in fa:
   if rm%q==0:
    eq=0
    while rm%q==0:rm//=q;eq+=1
    if q in frs:
     gt=0
     for ix,wv in EN(sp3[q]):
      if(ev:=hv(ux,vx,wv,q,eq)):
       gt+=ev;bq=(b[0]+b[1]*wv)%q
       if bq==0:return N_
       k=("fs",q,((a[0]+a[1]*wv)*PW(bq,-1,q))%q);o[k]=o.get(k,0)+(ev if ix==0 else-ev)
     if gt!=eq:return N_
    else:
     if eq%2 or ux%(qe:=q**(eq//2))or vx%qe or(b[0]%q==0 and b[1]%q==0):return N_
     if q%4==1:continue
     o[("fi",q)]=o.get(("fi",q),0)+((eq//2)if wm(a,wi(b,q),q)==fri[q]else-(eq//2))
  return o if rm==1 else N_
 rl,sn=[],ST()
 for br in RN(-4,5):
  for bi in RN(5):
   if br<=0 and bi==0:continue
   b=(br,bi);b2=zm(b,b);mb=(m*br,m*bi)
   for ar in RN(-8,9):
    for ai in RN(-8,9):
     if ar==0 and ai==0:continue
     a=(ar,ai);a2=zm(a,a);nfK=(a2[0]+b2[0],a2[1]+b2[1])
     if zn(nfK)%2==0 or(ff:=ff_(a,b,nfK))is N_ or(fg:=fg_((ar-mb[0],ai-mb[1])))is N_ or not wg1(a,b):continue
     if(cb:=MN((zm(u,a),zm(u,b))for u in un))not in sn:sn.add(cb);rl+=[(a,b,fg,ff)]
 while T_:
  cn={}
  for _,_,fg,ff in rl:
   for k,v in(*fg.items(),*ff.items()):
    if v%r:cn[k]=cn.get(k,0)+1
  nr=[x for x in rl if all(cn.get(k,0)>=2 for k,v in(*x[2].items(),*x[3].items())if v%r)]
  if LN(nr)==LN(rl):break
  rl=nr
 fi=SR({k for _,_,fg,ff in rl for k,v in(*fg.items(),*ff.items())if v%r});fm_={k:i for i,k in EN(fi)};nf=LN(fi);sm=[smf(a,b)[3]for a,b,_,_ in rl]
 def rref(mt,nc):
  pv=_rref(mt,nc,r);fr=[c for c in RN(nc)if c not in pv];bs=[]
  for fc in fr:
   v=[0]*nc;v[fc]=1
   for pc in SR(pv,reverse=T_):v[pc]=(-SM(vv*v[c]for c,vv in pv[pc].items()if c!=pc))%r
   bs+=[v]
  return LN(fr),bs
 m0,m1=[],[]
 for(_,_,fg,ff),sv in zip(rl,sm):
  d={};[d.__setitem__(fm_[k],(d.get(fm_[k],0)+sg*v)%r)for sg,mp in((1,fg),(-1,ff))for k,v in mp.items()if v%r];m0+=[d];d2=dict(d)
  if sv%r:d2[nf]=(-sv)%r
  m1+=[d2]
 n0,_=rref(m0,nf);n1,ks=rref(m1,nf+1);assert n0==0 and n1==1 and ks[0][nf];kv=ks[0]
 def wp(x,e):
  rs,bs=(1,0),(x[0]%p,x[1]%p)
  while e:
   if e&1:rs=wm(rs,bs,p)
   bs=wm(bs,bs,p);e>>=1
  return rs
 cof=(p*p-1)//r;gn=next(wp((t,1),cof)for t in RN(1,p)if wp((t,1),cof)!=(1,0));dl,cu={},(1,0)
 for e in RN(r):dl[cu]=e;cu=wm(cu,gn,p)
 rt=ST()
 for(a,b,fg,ff),sv in zip(rl,sm):
  td=dl[wp((a[0]-m*b[0],a[1]-m*b[1]),cof)];pg=SM(v*kv[fm_[k]]for k,v in fg.items()if v%r)%r;pf=(SM(v*kv[fm_[k]]for k,v in ff.items()if v%r)+sv*kv[nf])%r;assert pg==pf
  if td:rt.add((pg*PW(td,-1,r))%r)
  else:assert pg==0
 assert LN(rt)==1;gg=next(iter(rt));assert gg;os=[]
 for ar in RN(9,18):
  for ai in RN(-4,5):
   if(fg:=fg_((ar,ai)))and all(k in fm_ for k,v in fg.items()if v%r)and(td:=dl[wp((ar,ai),cof)]):
    os+=[(td,SM(v*kv[fm_[k]]for k,v in fg.items()if v%r)%r)]
    if LN(os)==2:break
  if LN(os)==2:break
 assert LN(os)==2 and(os[1][1]*PW(os[0][1],-1,r))%r==(os[1][0]*PW(os[0][0],-1,r))%r;return{"ok":T_,"tower":"K=Q(w),p=257,r=43|(p+1),g=Y-16,f=Y^2+1(L_f=Q(z12),rk1)","planted_vector":F0,"real_relations":LN(rl),"n_fb_orbits":nf,"n_sm_cols":1,"unit_sm_nonzero":T_,"sm_hom_and_rth_zero":T_,"nullity_without_sm":n0,"nullity_with_sm":n1,"sm_coeff_nonzero":T_,"all_rels_exact_dlp_gauge":gg,"oos_target_dlog_ratio_exact":T_}
Z4=(0,0,0,0);MONO=((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,0,1,0),(0,-1,0,1),(-1,0,0,0),(0,-1,0,0),(0,0,-1,0),(0,0,0,-1),(1,0,-1,0),(0,1,0,-1))
MMAP={v:i for i,v in EN(MONO)}
def _zm(a,k):
 k%=12;ng=k>=6
 if ng:k-=6
 c0,c1,c2,c3=a
 for _ in RN(k):c0,c1,c2,c3=-c3,c0,c1+c3,c2
 return(-c0,-c1,-c2,-c3)if ng else(c0,c1,c2,c3)
def _lcsge(U,eta,ltf,ltg,sd,rk=F0):
 rng=RR(sd);Uf=I(RD(.43918*U));Ug=U-Uf
 def sln(Us):
  x=LG(12.*Us*LG(12.*Us+3.))
  for _ in RN(10):x=LG(12.*Us*x)
  return x
 def sul(lB,rt,tl):
  lo,hi=1.6,12.
  for _ in RN(36):md=.5*(lo+hi);lo,hi=(md,hi)if(eta*rt/lB)*(D0(md-1.)/D0(md))<tl else(lo,md)
  return .5*(lo+hi)
 lnf,lng=sln(Uf),sln(Ug);uf,ug=sul(lnf,U/Uf,ltf),sul(lng,U/Ug,ltg);rf,rg=D0(uf),D0(ug);ws=[]
 for Uc,lB,us,rs in((Uf,lnf,uf,rf),(Ug,lng,ug,rg)):
  for j in RN(Uc):
   tg=12.*(j+.5);pj=MX(31.,tg*LG(tg+31.));pj=MX(31.,tg*LG(pj));pj=MX(31.,tg*LG(pj));ws+=[MN(.45,(12./pj)*(D0(us-MN(1.,LG(pj)/lB))/rs))]
 tw=SM(ws);sc=U/tw;pr=[w*sc for w in ws];al=[0]*U;sm=[i for i,v in EN(pr)if v<1.];lg=[i for i,v in EN(pr)if v>=1.]
 while sm and lg:s_,l_=sm.pop(),lg.pop();al[s_]=l_;pr[l_]+=pr[s_]-1.;(sm if pr[l_]<1. else lg).append(l_)
 for i in lg+sm:pr[i]=1.
 R=I(RD(eta*U));rw=[N_]*R;cr={};rn=rng.random;rg_=rng.randrange;eL=EX(-tw)
 for ri in RN(R):
  k,p_=0,1.
  while p_>eL:k+=1;p_*=rn()
  d={}
  for _ in RN(MX(4,k-1)):
   ix=rg_(U);c=ix if rn()<pr[ix]else al[ix];m=MONO[rg_(12)]
   if c in d:
    a=d[c];nv=(a[0]+m[0],a[1]+m[1],a[2]+m[2],a[3]+m[3])
    if nv==Z4:del d[c]
    else:d[c]=nv
   else:d[c]=m
  rw[ri]=d
  for c in d:cr.setdefault(c,ST()).add(ri)
 def kr(ri):
  if(d:=rw[ri])is N_:return
  for c in d:
   if(s:=cr.get(c)):
    s.discard(ri)
    if not s:cr.pop(c,N_)
  rw[ri]=N_
 def pg():
  q=[c for c,s in cr.items()if LN(s)==1]
  while q:
   if not(s:=cr.get(q.pop()))or LN(s)!=1 or(d:=rw[ri:=next(iter(s))])is N_:continue
   nb=LS(d);kr(ri);q+=[nc for nc in nb if(ns:=cr.get(nc))and LN(ns)==1]
 pg();et=16 if U>=1024 else 12;ar=[i for i,d in EN(rw)if d]
 while LN(ar)-LN(cr)>et:
  bt=MX(1,MN((LN(ar)-LN(cr)-et)//4,32));ar.sort(key=lambda i:LN(rw[i]),reverse=T_)
  for ri in ar[:bt]:
   if rw[ri]is not N_:kr(ri)
  pg();ar=[i for i,d in EN(rw)if d]
 bk={g:ST()for g in RN(2,7)}
 for c,s in cr.items():
  if 2<=LN(s)<=6:bk[LN(s)].add(c)
 def ub(c):
  ln=LN(s)if(s:=cr.get(c))else 0
  for g in RN(2,7):bk[g].discard(c)
  if 2<=ln<=6:bk[ln].add(c)
 while T_:
  did=F0
  for dg in(2,3,4,5,6):
   cd=N_
   for c in LS(bk[dg]):
    if not(s:=cr.get(c))or LN(s)!=dg:bk[dg].discard(c);continue
    rl=LS(s);mi=[MMAP.get(rw[r][c],-1)for r in rl]
    if any(k<0 for k in mi):continue
    od=SR(RN(dg),key=lambda ix:LN(rw[rl[ix]]));rp,kp=rl[od[0]],mi[od[0]];lp=LN(rw[rp])
    if dg>=3 and any(lp+LN(rw[rl[o]])>140 for o in od[1:]):continue
    cd=(rp,kp,[(rl[o],mi[o])for o in od[1:]]);break
   if cd is not N_:
    rp,kp,tg_=cd;dp=rw[rp];tc=ST(dp)
    for rt,kt in tg_:
     dt,sh=rw[rt],(kt-kp)%12
     for col,v0 in dp.items():
      sv=_zm(v0,sh)
      if col in dt:
       a=dt[col];nv=(a[0]-sv[0],a[1]-sv[1],a[2]-sv[2],a[3]-sv[3])
       if nv==Z4:
        del dt[col];cr[col].discard(rt)
        if not cr[col]:cr.pop(col,N_)
       else:dt[col]=nv
      else:dt[col]=(-sv[0],-sv[1],-sv[2],-sv[3]);cr.setdefault(col,ST()).add(rt)
     if not dt:kr(rt)
    kr(rp);pg()
    for col in tc:ub(col)
    did=T_;break
  if not did:
   ad=F0
   for c,s in cr.items():
    if 2<=LN(s)<=6 and c not in bk[LN(s)]:bk[LN(s)].add(c);ad=T_
   if not ad:break
 al=[d for d in rw if d];cls=SR({c for d in al for c in d});mx=ov=nz=c1=cu=0
 for d in al:
  for v in d.values():
   nz+=1
   for ck in v:
    if(av:=AB(ck)):c1+=I(av==1);cu+=I(1<av<(1<<30));ov+=I(av>=(1<<30));mx=MX(mx,av)
 rk1=dmin=N_
 if rk and LN(al)>LN(cls)>=4:
  q=1000033;z1=next(PW(g,(q-1)//12,q)for g in RN(2,100)if(PW(g,(q-1)//3,q)-PW(g,(q-1)//6,q)+1)%q==0);z2,z3=(z1*z1)%q,(z1*z1*z1)%q;cm={c:i for i,c in EN(cls)};fm=[{cm[c]:t for c,v in d.items()if(t:=(v[0]+v[1]*z1+v[2]*z2+v[3]*z3)%q)}for d in al]
  while T_:
   cc={}
   for ri,d in EN(fm):
    for c in d:cc.setdefault(c,[]).append(ri)
   c2=next((c for c,rs in cc.items()if LN(rs)<=2),N_)
   if c2 is N_:break
   rs=cc[c2]
   if LN(rs)==1:fm[rs[0]]={}
   else:
    r0,r1=rs;fc=(fm[r1][c2]*PW(fm[r0][c2],-1,q))%q
    for c,v in fm[r0].items():
     nv=(fm[r1].get(c,0)-fc*v)%q
     if nv:fm[r1][c]=nv
     else:fm[r1].pop(c,N_)
    fm[r0]={}
   fm=[d for d in fm if d]
  fc_=SR({c for d in fm for c in d});fmp={c:i for i,c in EN(fc_)};nc=LN(fc_);dmin=MN(SM(1 for d in fm if c in d)for c in fc_)if fc_ else 0
  if LN(fm)>nc>=4:
   vt=[rng.randrange(1,q)for _ in RN(nc+4)];iv3=PW(vt[nc+3],-1,q);mt=[]
   for d in fm:r0={fmp[c]:v for c,v in d.items()};s0,s1,s2=rng.randrange(q),rng.randrange(q),rng.randrange(q);r0[nc],r0[nc+1],r0[nc+2],r0[nc+3]=s0,s1,s2,((-(SM(v*vt[c]for c,v in r0.items())+s0*vt[nc]+s1*vt[nc+1]+s2*vt[nc+2]))*iv3)%q;mt+=[r0]
   rk1=(dmin>=3 and LN(_rref(mt,nc+4,q))==nc+3)
 return{"nl":LN(cls),"rows":LN(al),"c1":c1,"cu":cu,"nz":nz,"mx":mx,"ov":ov,"rk1":rk1,"dmin":dmin}
def chk_h3(af,ag,mnf,mng,full=F0):
 Bf,Bg,eta,b=R_["Bf"],R_["Bg"],R_["eta"],R_["b"];lf,lg=Bf*LG(2.),Bg*LG(2.);ur_=(2.**Bf/lf)/((2.**Bf/lf)+(2.**Bg/lg));ug_=1.-ur_;uf,ug=(mnf+af)/Bf,(mng+ag)/Bg;rf,rg=D0(uf),D0(ug);fs=((2.**b["pfhi0"]/b["pfhi0"])+(2.**b["pghi0"]/b["pghi0"]))/((2.**Bf/Bf)+(2.**Bg/Bg));M,cs=500,[]
 for B,ur,us,rs in((Bf,ur_,uf,rf),(Bg,ug_,ug,rg)):dy=(B-2.)/M;ws=[((2.**(y-B))*(B/y)*LG(2.)*dy,eta*(1./ur)*(2.**(B-y))/(B*LG(2.))*(D0(us-y/B)/rs))for y in(2.+(i+.5)*dy for i in RN(M))];sw=SM(x[0]for x in ws);cs.extend((ur*wr/sw,lm)for wr,lm in ws)
 ltf,ltg=eta*(1./ur_)/lf*(D0(uf-1.)/rf),eta*(1./ug_)/lg*(D0(ug-1.)/rg);w0=SM(w*lm for w,lm in cs)/eta;phi=.8
 for _ in RN(80):phi=EX(-SM(w*lm*EX(-phi*lm)for w,lm in cs)/eta)
 R1,C1=eta*phi,SM(w*(1.-(1.+phi*lm)*EX(-phi*lm))for w,lm in cs);E1=R1-C1;assert E1>.15;sch=((10,4),(11,4),(12,3),(13,3),(14,2))if full else((10,3),(11,3),(12,2),(13,2));lad,rkp,rkt,mxc,tov=[],0,0,0,0
 for p2,ns in sch:
  U=1<<p2;rs=[_lcsge(U,eta,ltf,ltg,7500+100*p2+s,rk=(p2<=11))for s in RN(ns)]
  for r in rs:assert r["rows"]>=r["nl"]and r["ov"]==0 and r["mx"]<(1<<30);mxc=MX(mxc,r["mx"]);tov+=r["ov"]
  if p2<=11:rkp+=SM(1 for r in rs if r["rk1"]);rkt+=ns
  tc=SM(r["c1"]+r["cu"]for r in rs);lad+=[{"log2U":p2,"seeds":ns,"nl_U":RD(SM(r["nl"]for r in rs)/(ns*U),5),"comp_U":RD(tc/(ns*U),5),"pm1":RD(SM(r["c1"]for r in rs)/MX(1,tc),4),"max_abs":MX(r["mx"]for r in rs)}]
 assert rkp==rkt and tov==0
 def ols(xs,ys):
  n=LN(xs);mx,my=SM(xs)/n,SM(ys)/n;bl=SM((x-mx)*(y-my)for x,y in zip(xs,ys))/SM((x-mx)**2 for x in xs);al=my-bl*mx;return EX(al),-bl,1.-SM((y-(al+bl*x))**2 for x,y in zip(xs,ys))/SM((y-my)**2 for y in ys)
 xs=[LG(LG(1<<r["log2U"]))for r in lad];Al,bl,r2l=ols(xs,[LG(r["nl_U"])for r in lad]);Ac,bc,r2c=ols(xs,[LG(r["comp_U"])for r in lad]);U_p=U_ub(Bf,Bg);lnU=LG(U_p);nlo,nco=Al*(lnU**(-bl)),Ac*(lnU**(-bc));nrs=.004037*((24.589/lnU)**3);wrc=RD(144.18*(4./3.),2);nlc=RD(MX(1.10*nrs,10.*nlo),6);ncc=RD(MX(5.*nco,wrc*nlc),6);assert r2l>.85 and nlc<.5*fs;yc=Bg+L2(nlc)
 for _ in RN(12):yc=Bg+L2(nlc*(yc/Bg))
 dst={Z4:1.};ec1,ecu={},{}
 for k in RN(1,9):
  nd={}
  for v,pv in dst.items():
   p12=pv/12.
   for m in MONO:nv=(v[0]+m[0],v[1]+m[1],v[2]+m[2],v[3]+m[3]);nd[nv]=nd.get(nv,0.)+p12
  dst=nd;ec1[k]=SM(pv*SM(1 for x in v if AB(x)==1)for v,pv in dst.items());ecu[k]=SM(pv*SM(1 for x in v if AB(x)>=2)for v,pv in dst.items())
 ps=[p for p in _primes(20000)if p>=17];pm=1.;y0=L2(20000.)
 for Lm in(2.,3.,5.,10.,20.,40.):
  t1=tu=0.
  for B,us,rs in((Bf,uf,rf),(Bg,ug,rg)):
   for qt in(1/16.,1/16.,1/25.,1/25.,1/81.,1/81.,1/121.,1/121.):
    lm=Lm*MN(.95,qt*(D0(us-4./B)/rs));pk=EX(-lm);sc_=pk
    for k in RN(1,9):pk*=lm/k;sc_+=pk;t1+=pk*ec1[k];tu+=pk*ecu[k]
    tu+=4.*MX(0.,1.-sc_)
   for p in ps:lm=Lm*MN(.95,(12./p)*(D0(us-L2(p)/B)/rs));em=EX(-lm);p1=lm*em;p2=.5*lm*p1;p3=lm*p2/3.;p4=MX(0.,1.-em-p1-p2-p3);t1+=(p1*ec1[1]+p2*ec1[2]+p3*ec1[3]+p4*1.5)/12.;tu+=(p1*ecu[1]+p2*ecu[2]+p3*ecu[3]+p4*2.5)/12.
   dy=(yc-y0)/100.
   for i in RN(100):y=y0+(i+.5)*dy;p=2.**y;do=(p/(12.*y))*dy;lm=Lm*(12./p)*(D0(us-y/B)/rs);em=EX(-lm);p1=lm*em;p2=.5*lm*p1;p3=MX(0.,1.-em-p1-p2);t1+=do*(p1*ec1[1]+p2*ec1[2]+p3*1.8);tu+=do*(p1*ecu[1]+p2*ecu[2]+p3*2.)
  pm=MN(pm,t1/(t1+tu))
 p1c=MN(.92,RD(pm-.05,4));assert pm>.97 and p1c==.92;Nl=U_p*nlc;ec=R_["la"]["e_core"];l2r,l2n=L2(FL(RB)),L2(Nl);pbw=2.*Nl/FL(RB);psm=1./FL(RB);ld=-(ec-3)*l2r;lpj=-2.*l2r;ls43=L2(9.)-5.*l2n;ls54=L2(2764.8)-6.*l2n;ls33=L2(6.)-3.*l2n-l2r;ls23=L2(3.)-l2n-2.*l2r;pla=pbw+psm+2.**ls43+2.**ls54+2.**ls33+2.**ls23+2.**lpj;assert L2(pla)<-200.
 return{"ok":T_,"scope":"H3_RSA768_full_Pareto_derived_pLA","prod_2core":{"ltf":RD(ltf,4),"ltg":RD(ltg,4),"w0":RD(w0,2),"fb_frac":RD(fs,6),"phi":RD(phi,5),"R1_U":RD(R1,5),"C1_U":RD(C1,5),"E1_U":RD(E1,5)},"ladder":lad,"ols":{"nl_exp":RD(bl,3),"nl_R2":RD(r2l,4),"nl_pred":FL(f"{nlo:.3e}"),"nc_exp":RD(bc,3),"nc_R2":RD(r2c,4),"nc_pred":FL(f"{nco:.3e}"),"rsa768_scaled":RD(nrs,6),"rsa768_w_comp":wrc},"derived_caps":{"nl_cap":nlc,"comp_cap":ncc,"y_cut_derived":RD(yc,3),"pm1_min_exact":RD(pm,4),"pm1_cap":p1c,"z12_k2_pm1_exact":RD(ec1[2]/(ec1[2]+ecu[2]),4),"max_abs":mxc,"ov30":tov},"la_prob_4part":{"e_core":ec,"fq_sparse_core_dmin3_rk1":f"{rkp}/{rkt}","l2_p_bw_uncond":RD(L2(pbw),2),"l2_p_dense_rank":RD(ld,1),"l2_p_sq_proj":RD(lpj,1),"l2_p_sparse_k4s3":RD(ls43,1),"l2_p_sparse_k5s4":RD(ls54,1),"l2_p_sparse_k3s3":RD(ls33,1),"l2_p_sparse_k2s3":RD(ls23,1),"l2_p_sm_heur":RD(L2(psm),2),"p_la_total":FL(f"{pla:.4e}"),"l2_p_la":RD(L2(pla),2)}}
def _mr4(Bb):
 s,oe=R_["s1b"],R_["oe"];M1,M2=(oe*(2.**60)*8.)/4.5,(oe*(2.**60)*8.)/9.
 if Bb<=1024.:return 12.*MX(1.,Bb/64.)**2*s["saf"]
 if Bb<=FL(1<<59):N=FL(1<<I(CE(L2(MX(2.,2.*CE(Bb/64.))))));assert N<=(1<<54);return(3.*(1.5*N*L2(N)*s["bfly"]+N*(s["pw"]+s["sc"]))+(3.*N*160. if Bb*2.>M1 else 0.)+N*s["c1"])*s["saf"]
 L=I(CE(Bb/58.));k=MX(1,(L-1).bit_length());k1=(k+1)//2;k2=k-k1;N2=FL(1<<(k1+1))*FL(1<<(k2+1));return(3.*(N2*(k1+1+1.5*(k2+1))*s["bfly"]+N2*(s["pw"]+s["sc"]))+(9.*N2*160. if Bb*2.>M2 else 0.)+N2*s["c2"])*s["saf"]
def _rr4(nlv,lb):
 if nlv<=1.:return 40000.
 dm,kp=R_["s1b"]["dm"],MX(1,I(4096.//lb));ns,sb,nll=MX(1.,nlv/kp),kp*lb,lb/64.;isb=SM(12.*j*nll*nll*1.15 for j in RN(1,kp));isd=kp*24.*(sb/64.)*nll*1.15;h,ip,id_,cn,cb=I(CE(L2(ns))),0.,0.,ns,sb
 for _ in RN(h):pr=cn*.5;mi=_mr4(cb);ip+=pr*mi;id_+=pr*dm*mi;cn,cb=pr,cb*2.
 return ns*(isb+isd)+ip+3.*_mr4(ns*sb)+id_+nlv*MX(25000.,40000.*(lb/636.))
def _oecm(qb,Nb,_ec={},_ep={}):
 qk,ky=RD(qb,2),(RD(qb,2),I(CE(Nb/64.)))
 if ky in _ec:return _ec[ky]
 if qk not in _ep:
  pl,qe,rh=[],MX(1.,qb-3.5),DK.rho
  for lb2 in RN(16,125):
   l1=lb2*.5;B1=2.**l1;p1=rh(qe/l1)
   for m in(10.,25.,60.,120.,250.):
    B2=B1*m;dt=(L2(B2)-l1)/40.;pr=MN(1.,p1+SM((rh((qe-(l1+(i+.5)*dt))/l1)*(dt/(l1+(i+.5)*dt)))if(l1+(i+.5)*dt)<qe else(dt/(l1+(i+.5)*dt))for i in RN(40)))
    if pr>1e-40:pl+=[(B1,B2,pr)]
  _ep[qk]=pl
 L=MX(4,I(CE(Nb/64.)));mi=212. if L==4 else(13.25*L*L+80. if L<=16 else MN(13.25*L*L+80.,2.5*_mr4(L*64.)));bs=(1e300,N_)
 for B1,B2,pr in _ep[qk]:
  c1=1.15*(40.*mi+(B1/LG(2.))*(11.*mi+93.*L)+2.5*MX(0.,B2/LG(B2)-B1/LG(B1))*(mi+40.*L))
  if c1/pr<bs[0]:bs=(c1/pr,(B1,B2,pr,c1))
 _ec[ky]=bs;return bs
def _dn(Q,sd,Bft,Bgt,Bfm,Bgm,af,ag):
 ld,lv,CR=2.*L2(1125./16.)+Q,L2((pi**8)/40320.),2.**20
 for _ in RN(8):ls=(L2(CR)+ld-lv)/16.+L2(1.25)-1.5;nf=16.+72.*ls-(Q if sd=="f" else 0.);ng=8.*L2X+24.*ls+2.-(Q if sd=="g" else 0.);CR=MX(1024.,2./MX(1e-40,D0((nf+af)/Bft)*D0((ng+ag)/Bgt)))
 cu=(CR*6000.+(2.**22)+3.*(_oecm(Bft,nf)[0]+_oecm(Bgt,ng)[0]))/CC0;ec=lambda Ne,lo,hi:0. if hi<=lo+1e-9 else SM((D0((Ne-(lo+(i+.5)*(hi-lo)/40.))/hi)/D0(Ne/hi))*(((hi-lo)/40.)/(lo+(i+.5)*(hi-lo)/40.))for i in RN(40)if(lo+(i+.5)*(hi-lo)/40.)<Ne);return{"cu":cu,"cf":ec(nf+af,Bfm,Bft),"cg":ec(ng+ag,Bgm,Bgt)}
def _s4d(af,ag):
 Bf,Bg=R_["Bf"],R_["Bg"];br=lambda Nr,Qm,lo,hi:SM((D0((Nr-(lo+(i+.5)*(hi-lo)/200.))/Qm)/D0(Nr/Qm))*(((hi-lo)/200.)/(lo+(i+.5)*(hi-lo)/200.))for i in RN(200));nsp,ni0,ni1,ni2=br(2572.,275.,Bg,275.),.5*br(2572.,275.,.5*Bg,137.5),.5*br(2572.,275.,137.5,195.),.5*br(2572.,275.,195.,275.);d5=_dn(550.,"g",390.,390.,Bf,Bg,af,ag);ut=ni2*d5["cu"];nfc,ngc=ni2*d5["cf"],ni2*d5["cg"]+ni1;nl=[RD(ni2,2),RD(nfc+ngc,1)]
 for(Q,Bft,Bgt),ag_ in(((390.,275.,275.),nsp+ni0),((275.,180.,180.),0.),((180.,125.,125.),0.),((125.,94.,94.),0.),((94.,78.,78.),0.),((78.,Bf,Bg),0.)):df,dg=_dn(Q,"f",Bft,Bgt,Bf,Bg,af,ag),_dn(Q,"g",Bft,Bgt,Bf,Bg,af,ag);ut+=nfc*df["cu"]+ngc*dg["cu"];nfc,ngc=nfc*df["cf"]+ngc*dg["cf"],nfc*df["cg"]+ngc*dg["cg"]+ag_;nl+=[RD(nfc+ngc,1)]
 assert nfc==0.==ngc;ub=((1./(.918*D0(2572./275.)))*(2.5*_oecm(275.,2572.)[0]))/CC0;return{"l0_split":RD(nsp,2),"l0_inert_le275":RD(ni0,2),"l0_inert_275_390":RD(ni1,2),"l0_inert_390_550":RD(ni2,2),"nodes":nl,"s4":L2(2.*(ub+ut)),"ecm":L2(2.*(ut+.25*ub)),"desc_only":L2(2.*ut)}
def _lll(B,dl=.99):
 n,m,B=LN(B),LN(B[0]),[r[:]for r in B];bo,mu,B2=[[0.]*m for _ in RN(n)],[[0.]*n for _ in RN(n)],[0.]*n
 def gs(st):
  for i in RN(st,n):
   bi=B[i];v=[FL(x)for x in bi]
   for j in RN(i):
    bj=bo[j];mu[i][j]=(s:=SM(bi[t]*bj[t]for t in RN(m))/B2[j]if B2[j]>1e-18 else 0.)
    for t in RN(m):v[t]-=s*bj[t]
   bo[i]=v;B2[i]=SM(x*x for x in v)
 gs(0);k=1
 while k<n:
  for j in RN(k-1,-1,-1):
   if AB(mu[k][j])>.5:
    q=I(RD(mu[k][j]));bk,bj=B[k],B[j]
    for t in RN(m):bk[t]-=q*bj[t]
    for t in RN(j):mu[k][t]-=q*mu[j][t]
    mu[k][j]-=q
  if B2[k]>=(dl-mu[k][k-1]**2)*B2[k-1]:k+=1
  else:B[k],B[k-1]=B[k-1],B[k];gs(k-1);k=MX(k-1,1)
 return B
def _h90_q11():
 xt=-56;rt=xt**4-xt**2+1;pt=((xt-1)**2*rt)//3+xt;xp=xt%pt;assert rt==9831361 and pt==10647363907 and pt%15==7 and(pt**6+1)%rt==0
 def fm(a,b):
  rs=[0]*12
  for j1 in RN(3):
   for j2 in RN(3):
    j=j1+j2;mx,jo=(xp,j-3)if j>=3 else(1,j);tp_=[0]*7
    for i1 in RN(4):
     if(c1:=a[4*j1+i1]):
      for i2 in RN(4):tp_[i1+i2]+=c1*b[4*j2+i2]
    tp_[0]+=tp_[5];tp_[1]+=tp_[6];c4=tp_[4]
    for i in RN(4):rs[4*jo+i]=(rs[4*jo+i]+(tp_[i]-c4)*mx)%pt
  return rs
 def fp(a,e):
  rs,bs=[1]+[0]*11,LS(a)
  while e:
   if e&1:rs=fm(rs,bs)
   bs=fm(bs,bs);e>>=1
  return rs
 f6=lambda a:[v for j in RN(3)for c0,c1,c2,c3 in(a[4*j:4*j+4],)for v in((c0-c1)%pt,(-c1)%pt,(c3-c1)%pt,(c2-c1)%pt)]
 g0=fp([3,5,7,11,13,17,19,23,29,31,37,41],(pt**6-1)*((pt**6+1)//rt));assert fp(g0,rt)==[1]+[0]*11!=g0;tt=fp(g0,4192837);assert fm(tt,f6(tt))==[1]+[0]*11;A=[[0]*12 for _ in RN(12)]
 for k in RN(12):
  ek=[1 if i==k else 0 for i in RN(12)];im=[(ek[i]-fm(tt,f6(ek))[i])%pt for i in RN(12)]
  for rr in RN(12):A[rr][k]=im[rr]
 pv,ri=[],0
 for cl in RN(12):
  if(p0:=next((rr for rr in RN(ri,12)if A[rr][cl]%pt),N_))is N_:continue
  A[ri],A[p0]=A[p0],A[ri];iv=PW(A[ri][cl],-1,pt);A[ri]=[(v*iv)%pt for v in A[ri]]
  for rr in RN(12):
   if rr!=ri and(f:=A[rr][cl]):A[rr]=[(A[rr][cc]-f*A[ri][cc])%pt for cc in RN(12)]
  pv+=[cl];ri+=1
 assert LN(pv)==6;fc=[c for c in RN(12)if c not in pv];B=[[(pt if j==p_ else 0)for j in RN(12)]for p_ in pv]+[[(1 if j==f_ else((-A[pv.index(j)][f_])%pt if j in pv else 0))for j in RN(12)]for f_ in fc];ws=[1]*4+[4]*4+[15]*4;Br=_lll([[r[i]*ws[i]for i in RN(12)]for r in B]);al=[Br[0][i]//ws[i]for i in RN(12)];am=[c%pt for c in al];assert am==fm(tt,f6(am))and any(am)and fm(fp(am,pt**6-1),tt)==[1]+[0]*11;xc,w3,nv=(AB(xt)**(1./3.))*cmath.exp(1j*pi/3.),cmath.exp(2j*pi/3.),1.+0j
 for m in(1,2,3,4):
  z5=cmath.exp(2j*pi*m/5.)
  for k3 in RN(3):yv=xc*(w3**k3);nv*=SM(al[4*j+i]*(z5**i)*(yv**j)for j in RN(3)for i in RN(4))
 l2n,l2d=L2(AB(nv)),6.*L2(pt)+4.*L2(AB(xt));assert l2n<l2d+2.;q=1031;assert q%15==11;z5r=next(t for t in RN(2,q)if(t**4+t**3+t**2+t+1)%q==0);yr=PW(xt%q,PW(3,-1,q-1),q);assert(yr**3-xt)%q==0;b6,p5=[[(q if j==ix else 0)for j in RN(16)]for ix in(0,4)],[PW(z5r,i,q)for i in RN(4)]
 for i in RN(1,4):r0,r1=[0]*16,[0]*16;r0[i],r0[0]=1,(-p5[i])%q;r1[4+i],r1[4]=1,(-p5[i])%q;b6+=[r0,r1]
 for i in RN(4):r0,r1=[0]*16,[0]*16;r0[8+i],r0[0]=1,(yr*p5[i])%q;r1[12+i],r1[4]=1,(yr*p5[i])%q;b6+=[r0,r1]
 v0=_lll(b6)[0];assert SM((v0[i]-yr*v0[8+i])*p5[i]for i in RN(4))%q==0 and SM((v0[4+i]-yr*v0[12+i])*p5[i]for i in RN(4))%q==0
 def om(u,v):
  def zm(a,b):
   t=[0]*7
   for i in RN(4):
    for j in RN(4):t[i+j]+=a[i]*b[j]
   t[0]+=t[5];t[1]+=t[6];c4=t[4];return[t[i]-c4 for i in RN(4)]
  p00,p01,p10,p11=zm(u[:4],v[:4]),zm(u[:4],v[4:]),zm(u[4:],v[:4]),zm(u[4:],v[4:]);return[p00[i]-p11[i]for i in RN(4)]+[p01[i]+p10[i]-p11[i]for i in RN(4)]
 a3,b3=om(om(v0[:8],v0[:8]),v0[:8]),om(om(v0[8:],v0[8:]),v0[8:]);gok=[a3[i]-xt*b3[i]for i in RN(8)];A8,pv_,sg=[om(gok,[1 if j==k else 0 for j in RN(8)])for k in RN(8)],1,1
 for k in RN(7):
  if A8[k][k]==0:sw=next(i for i in RN(k+1,8)if A8[i][k]);A8[k],A8[sw]=A8[sw],A8[k];sg=-sg
  pv2=A8[k][k]
  for i in RN(k+1,8):
   for j in RN(k+1,8):A8[i][j]=(A8[i][j]*pv2-A8[i][k]*A8[k][j])//pv_
  pv_=pv2
 nlg=sg*A8[7][7];assert nlg and nlg%(q*q)==0;return{"p":pt,"r":rt,"rank_H90":6,"lift_verified":T_,"log2_N_Mg":RD(l2n,2),"log2_det_w":RD(l2d,2),"q11_inert_descent":{"q":q,"q_mod_15":11,"exact_q2_divides_N_Lg":T_,"log2_N_Lg_cof":RD(L2(AB(nlg)/(q*q)),2)}}
def chk_s4(af,ag,full=F0):
 s4=_s4d(af,ag);assert s4["s4"]<R_["cr"]["both_individual_logs_exp"]and s4["ecm"]<R_["cr"]["ecm_split_expected"];h90=_h90_q11();lll=33.-6.*L2(12.)+4.5*L2(5.)+6.*L2(3.)+6.*L2(PB)+4.*L2(AX);assert lll<2572.
 o={"ok":T_,"lll_ub":RD(lll,4),"term":{"Bf":R_["Bf"],"Bg":R_["Bg"],"out_f":0,"out_g":0},"l0_split":s4["l0_split"],"l0_inert_le275":s4["l0_inert_le275"],"l0_inert_275_390":s4["l0_inert_275_390"],"l0_inert_390_550":s4["l0_inert_390_550"],"nodes":s4["nodes"],"desc_only_l2":RD(s4["desc_only"],2),"s4_l2":RD(s4["s4"],4),"ecm_l2":RD(s4["ecm"],4),"real_bls12_h90_q11":h90}
 if not full:return o
 xd=lambda X,Z,A,N:(((X+Z)**2*(X-Z)**2)%N,(((X+Z)**2-(X-Z)**2)*((X-Z)**2+A*((X+Z)**2-(X-Z)**2)))%N);xa=lambda a,b,c,d,X,Z,N:((Z*((a-b)*(c+d)+(a+b)*(c-d))**2)%N,(X*((a-b)*(c+d)-(a+b)*(c-d))**2)%N)
 def xm(k,X,Z,A,N):
  a,b,c,d=1,0,X%N,Z%N
  for bt in bin(k)[2:]:(c,d),(a,b)=(xa(a,b,c,d,X,Z,N),xd(a,b,A,N))if bt=="0" else(xd(c,d,A,N),xa(a,b,c,d,X,Z,N))
  return a,b
 t16=[p for p in _primes(200)if p>=101]
 for p in t16:x=(5*PW(-7,-1,p))%p;d=((2*x*x-1)*PW(x,-4,p))%p;A=(2*(1+d)*PW(1-d,-1,p))%p;assert(1+SM(1 if(v:=(u**3+A*u*u+u)%p)==0 else(2 if PW(v,(p-1)//2,p)==1 else 0)for u in RN(p)))%16==0
 N1=50023*50033;x1=(5*PW(-7,-1,N1))%N1;A24=PW(1-((2*x1*x1-1)*PW(x1,-4,N1))%N1,-1,N1);X,Z=3,1
 for p in _primes(30):
  pk=p
  while pk*p<=30:pk*=p
  X,Z=xm(pk,X,Z,A24,N1)
 assert gcd(Z,N1)==1;bb=[xm(i,X,Z,A24,N1)for i in RN(1,16)if gcd(i,30)==1];DX,DZ=xm(30,X,Z,A24,N1);gp,gc,w1=(1,0),(DX,DZ),1
 for k in RN(1,13):
  if(w1:=next((g for bx,bz in bb if 1<(g:=gcd((gc[0]*bz-bx*gc[1])%N1,N1))<N1),1))>1:break
  gp,gc=gc,(xd(DX,DZ,A24,N1)if k==1 else xa(gc[0],gc[1],DX,DZ,gp[0],gp[1],N1))
 X2,Z2=xm(2,4,1,2,35);DX2,DZ2=xm(4,X2,Z2,2,35);w2=gcd((DX2*Z2-X2*DZ2)%35,35);assert w1==50023 and gcd(Z2,35)==1 and w2==5;o.update({"t16":f"{LN(t16)}/{LN(t16)}","w1":w1,"w2":[w2,35//w2]});return o
def _tx64(l):
 if l<=0.:return 0.
 if l>45.:return(69./64.)*l+24.
 ei,pk=0.,EX(-l)
 for k in RN(1,160):
  pk*=l/k;q,r=divmod(k,64);ei+=pk*(10.+q*69.+(35 if r&32 else 0)+(19 if r&16 else 0)+(11 if r&8 else 0)+(7 if r&4 else 0)+(3 if r&2 else 0)+(3 if r&1 else 0))
  if k>l+25. and pk<1e-14:break
 return ei
T16=(5,10,9,13,15,20,19,23,23,28,27,31,33,38,37,41)
def _tx16(l):
 if l<=0.:return 8.
 if l>45.:return 33.+(37./16.)*l
 p0=EX(-l);e,pk=p0*8.,p0
 for k in RN(1,160):
  pk*=l/k;q,r=divmod(k,16);e+=pk*(11.+37.*q+T16[r])
  if k>l+24. and pk<1e-13:break
 return e
def _ihyb(plo,pw1,pw2,phi,L4,Lh,dk_,rmx,w2=16,w3=32,pp=20.):
 w1,w2_=1./(VA0*L4),1./(VA0*Lh);ih=ip=nok=nb=ne=0.
 if pw1>plo:
  dx=(pw1-plo)/100.
  for i in RN(100):x=plo+(i+.5)*dx;p=2.**x;d=(p/x)*dx;dk=d*dk_;nok+=dk;nb+=dk*(p+2.);ne+=d*VA0;ip+=dk*(pp+L4*7.)*w1;ih+=dk*SM(pr*L4*(_tx64(r*VA0/p)+1.5625*r*VA0/p)for r,pr in rmx)*w1
 for a,b,W in((pw1,pw2,w2),(pw2,phi,w3)):
  if b<=a:continue
  nb_,dx=(Lh-1.)/W+1.,(b-a)/60.
  for i in RN(60):
   x=a+(i+.5)*dx;p=2.**x;d=(p/x)*dx;dk=d*dk_;nok+=dk;nb+=dk*(p/W+2.);ne+=d*VA0;eip=eih=0.
   for r,pr in rmx:lb=W*r*VA0/p;pc=1.-EX(-nb_*lb);eip+=pr*(pp+pc*6.);eih+=pr*(nb_*_tx16(lb)-(1.-pc)*nb_*8.+nb_*1.5625*lb)
   ip+=dk*eip*w2_;ih+=dk*eih*w2_
 return{"ih":ih,"ip":ip,"nok":nok,"nb":nb,"ne":ne}
def _tp(cof,Pt,Bb,km=8):
 if cof<=Pt:return 1.,0.,1.
 om=lambda v:0. if v<1. else(1./v if v<=2. else((1.+LG(v-1.))/v if v<=3. else .57));dy,md,ma=(cof-Pt)/200.,D0(cof/Pt),0.;mt=md
 for i in RN(200):
  y=Pt+(i+.5)*dy;w=D0((cof-y)/Pt)*(om(y/Pt)/Pt)*dy;mt+=w;md+=w if y<=Bb else 0.;ma+=w if(y>Bb and any(k*Pt<=y<=k*Bb for k in RN(2,km+1)))else 0.
 return md/MX(mt,1e-30),ma/MX(mt,1e-30),(md+ma)/MX(mt,1e-30)
def chk_ldg(lm,lis,lh3):
 Bf,Bg,eta,b,rep,la=R_["Bf"],R_["Bg"],R_["eta"],R_["b"],R_["rep"],R_["la"];af,ag,cp=lm["af_1500k"],lm["ag_1500k"],lm["cof_pad"];L4,Lh,gp,lsg,lsf,ls1,ret=lis["L248_ch"],lis["L124_ch"],lis["gap_pad"],lis["l2_sg"],lis["l2_sf"],lis["l2_s1"],lis["ret_live"];dp=lis["discretization_convergence"]["disc_pad_b"];w0=1./(1.+R_["rbr"]**8);w1=1.-w0;fmx=[(3,39./46.),(6,6./46.),(9,1./46.)];dg,df,dng,dnf=1./3.,23./81.,1./24.,2131./34992.;w2,w3=R_["w2"],R_["w3"]
 sg0,sf0,sg1,sf1=(_ihyb(b[lo],b["pw1"],b["pw2"],b[hi],L4,Lh,dk,rm,w2,w3)for lo,hi,dk,rm in(("pglo","pghi0",dg,[(3,1.)]),("pflo","pfhi0",df,fmx),("pglo","pghi1",dg,[(3,1.)]),("pflo","pfhi1",df,fmx)));ih,ip=w0*(sg0["ih"]+sf0["ih"])+w1*(sg1["ih"]+sf1["ih"]),w0*(sg0["ip"]+sf0["ip"])+w1*(sg1["ip"]+sf1["ip"])
 rsh=lambda sg,sf:(((sg["ne"]+sf["ne"])*rep["ent"]+(sg["nb"]+sf["nb"])*rep["bkt"])/CC0,L2(8.*(sg["ne"]+sf["ne"]+sg["nb"]+sf["nb"])+48.*(sg["nok"]+sf["nok"])));r0,m0=rsh(sg0,sf0);r1,m1=rsh(sg1,sf1);ru=r0+r1;lV=lis["Vk_live"]+lis["pad_3se_Vk"]+gp+cp+dp;Vt,LtM,lf,lg=2.**lV,0.,R_["lf"],R_["lg"]
 def p1b(pfm,pgm):
  nonlocal LtM;s1=2.**ls1;c1,c2,mg=pfm+(Bf-pfm)/3.,pfm+2.*(Bf-pfm)/3.,.5*(pgm+Bg);af1,af2,af3,ag1,ag2=(_tp(l_,p_,B_)[1]for l_,p_,B_ in((lf,pfm,Bf),(lf,c1,Bf),(lf,c2,Bf),(lg,pgm,Bg),(lg,mg,Bg)));rf=MN(.95,SQ(ret/.80));rg=ret/rf;pfs,pgs=MN(1.,(D0((lis["mean_nf"]+af)/Bf)*rf)/(2.**lsf)),MN(1.,(D0((lis["mean_ng"]+ag)/Bg)*rg)/(2.**lsg));cr=s1*1700.*(1.+pfs)
  for lo,hi,pin,lb,dn in((1.,pfm,s1,656.,dnf),(pfm,c1,s1*af1,lf,dnf),(c1,c2,s1*af2,lf,dnf),(c2,Bf,s1*af3,lf,dnf),(1.,pgm,s1*pfs,688.,dng),(pgm,mg,s1*pfs*ag1,lg,dng),(mg,Bg,s1*pfs*ag2,lg,dng)):
   Lt=MX(1.,(ppb_dus(lo,MN(hi,36.))if lo<36. else 0.)+(dn*ppb_dus(MX(lo,36.),hi)if hi>36. else 0.));LtM=MX(LtM,Lt);Mt=Vt*pin;Mb=MN(Mt,MX(1.,Lt/lb));nbb,LX=MX(1.,Mt/Mb),Mb*lb;ir=_rr4(Mb,lb)+(MX(0.,Lt/LX-1.)*2.5*_mr4(LX)if LX<Lt else 0.);cr+=(nbb*ir)/Vt
  return cr+s1*pfs*pgs*(110.*256.*150.*4.+.05*5.0e7*2.)
 c1=w0*p1b(b["pfhi0"],b["pghi0"])+w1*p1b(b["pfhi1"],b["pghi1"]);sl=52./256.+2.5*(19./256.)+64./128.+200.*(2.**lsg);sc=sl+.004;il,ic=ih+ip+sl+c1,ih+ip+sc+c1;cl_,cc_=(Vt*il)/CC0,(Vt*ic)/CC0;U=U_ub(Bf,Bg);dc_=lh3["derived_caps"];n_,c_,p1=dc_["nl_cap"],dc_["comp_cap"],dc_["pm1_cap"];Nl,Nc=U*n_,U*c_;bsp=p1*10.+(1.-p1)*25.+(Nl*36.)/Nc;assert bsp<la["x4"]and 948.<la["rv"];Mb,Nb=la["m"],la["n"];Lq,dm=Nl/Mb+Nl/Nb+64.,Nl/Nb+32.
 iw=1.02*(Nb*Lq*(Nc*la["x4"]+Nl*la["rv"])+dm*(Nc*la["x1"]+Nl*la["rv"]))*la["ov"];ik=1.02*dm*Nb*Nl*283.*la["ov"];ig=20.*((Mb+Nb)**3)*Lq*L2(Lq)*236.*la["ov"];ix=4.*Mb*(Mb+Nb)*Lq*CE(L2(Lq))*8.*la["ov"];lu=(iw+ik+ig+ix+.008*(iw+ik)+30.*U*283.*la["ov"])/CC0;su=((12500.*U+25.*U*L2(U)+100.*Nc)*1.15)/CC0;s4=_s4d(af,ag)
 lr={"relation_collection":L2(cl_),"replay_a_tables":L2(ru),"linear_algebra":L2(lu),"both_individual_logs_exp":s4["s4"],"ecm_split_expected":s4["ecm"],"meter_safety_buffer":93.4,"band_split_and_trees":90.,"canonical_a_domain":84.,"schirokauer_maps":80.,"final_log_backsub_and_gcd":72.,"matrix_generator_sge":L2(su),"stage0_setup_and_iso":45.,"final_verification_Q":40.};cr=R_["cr"];assert LN(lr)==LN(cr)==13 and all(cr[k]>lr[k]for k in lr)and cr["relation_collection"]>=L2(cc_)
 wl,wc=SM(2.**v for v in lr.values()),SM(2.**v for v in cr.values());ll,lc=L2(wl),L2(wc);pb=CE((2.**m0)/(12.*(2.**60)));g0=248.*VA0*(2.**lsg)*8.;he=(186.*(2.**56))/(2.**60);sp=L2((2.**m0)+248.*(VA0+262144.)+2.*186.*(2.**56)+g0);ntB,cb,nn,lv=0.,4080.,LtM/4080.,0
 while nn>=1.:
  if cb<=(1<<59):ntB+=nn*(FL(1<<I(CE(L2(MX(2.,2.*CE(cb/64.))))))*24.)
  else:L_=I(CE(cb/58.));k=MX(1,(L_-1).bit_length());k1=(k+1)//2;ntB+=nn*(FL(1<<(k1+1))*FL(1<<(k-k1+1))*24.)
  cb*=2.;nn*=.5;lv+=1
 sB=ntB+(9.+lv)*LtM/8.;se=sB/(2048.*(2.**60));s2=L2((eta*U)*18.*12.);mB=Nc*(p1*8.+(1.-p1)*12.)+Nl*64.;vB=Nl*32.;le=(mB/la["p"]+vB*(1.+(Nb+6.)/la["p"]))/(2.**60);dpi=Lq*Nb/(Mb+Nb);gB=(Mb*Nb*Lq+Mb*(Mb+Nb)*Lq+4.*((Mb+Nb)**2)*dpi)*32.;ge=gB/(la["plg"]*(2.**60));s3=L2(mB+(Nb+6.+la["p"])*vB+gB);assert VA0/(2.**60)<.07 and he<12.<16. and se<R_["oe"]<16. and le<10.<16. and ge<8.<16. and MX(sp,L2(sB),s2,s3)<R_["ram_ub"]<120.
 pm=wc/(2.**124.);p_la=lh3["la_prob_4part"]["p_la_total"];pf=pm+2.**(-48)+2.**(-100)+2.**(-132)+p_la;po=1.-pf;Kg,t64=R_["Kg"],1<<64;assert pf<2.0e-5 and po>.99998 and FR(Kg,t64)*FR(99998,100000)>=FR(39,100)and FR(Kg-1,t64)*FR(99998,100000)<FR(39,100);pc=(Kg/FL(t64))*po;num,den=R_["num"],R_["den"];gc,gl=lc+L2(Kg/FL(t64)),ll+L2(Kg/FL(t64));slk=FL(FR(num,den))-gc;assert pc>.39 and gcd(num,den)==1 and slk>2.*lis["se_Vk"]
 return{"ok":T_,"claim":f"{num}/{den}","eta":eta,"U_l2":RD(L2(U),6),"U_ub":FL(f"{U:.8e}"),"Vk_ch_live":RD(lV,6),"disc_pad_b":dp,"ih":RD(ih,6),"ip":RD(ip,6),"c_s1b":RD(c1,6),"slab_live":RD(sl,6),"slab_ch":RD(sc,6),"ipc_live":RD(il,6),"ipc_ch":RD(ic,6),"blended_spmv_comp":RD(bsp,4),"la_mat_eib":RD(mB/(2.**60),4),"la_owner_eib":RD(le,4),"P_lingen":la["plg"],"lingen_owner_eib":RD(ge,4),"half_slab_eib":RD(he,4),"opened_live":RD(ll,6),"opened_ch":RD(lc,6),"gated_live":RD(gl,6),"gated_ch":RD(gc,6),"slack_b":RD(slk,6),"slack_over_live_SE_x":RD(slk/lis["se_Vk"],2),"p_meter_cond_H":FL(f"{pm:.7e}"),"p_open_cond_H":RD(po,8),"p_cold_cond_H":RD(pc,8),"ledger_live":{k:RD(v,6)for k,v in lr.items()},"la_sens_gated_live":{f"{m}x":RD(L2(wl+(m-1.)*lu)+L2(Kg/FL(t64)),6)for m in(1.,1.15,1.25,1.5)},"ram_l2_B":{"S1_sh0":RD(sp,4),"S1_sh1":RD(m1,4),"P_banks":pb,"S1b":RD(L2(sB),4),"S1b_eib":RD(se,4),"S2":RD(s2,4),"S3":RD(s3,4),"S3_lingen":RD(L2(gB),4)}}
def main():
 t0=time.perf_counter();raw=sys.stdin.read().strip();inp=json.loads(raw)if raw else{}
 if not IS(inp,dict)or(ST(inp)-{"part"}):raise ValueError("Bad input")
 if(pt:=inp.get("part","core"))not in("core","diagnostics","all"):raise ValueError(f"Bad part: {pt}")
 fl=(pt!="core");out={"package_id":R_["pkg"],"config_tag":R_["tag"],"part":pt,"reads_external_files":0,"offline_mc_constants":0}
 if pt!="diagnostics":out["stage0"]=chk_s0();out["chebotarev"]=chk_chb(26000 if fl else 15000)
 lm=chk_mrp(2000000 if fl else 1500000);out["murphy_norms"]=lm;lis=chk_rqmc(lm["af_1500k"],lm["ag_1500k"],32 if fl else 24,2048 if fl else 1024,8000 if fl else 1200,fD=fl);out["live_panjer_rqmc_chords"]=lis
 if pt!="diagnostics":out["ntt_r4_shoup61"]=chk_ntt();out["kernels_e2e"]=chk_krn()
 out["genuine_tower_dlp_sm"]=chk_twr();lh3=chk_h3(lm["af_1500k"],lm["ag_1500k"],lis["mean_nf"],lis["mean_ng"],full=fl);out["h3_extrap_and_la_prob"]=lh3;out["stage4_derived"]=chk_s4(lm["af_1500k"],lm["ag_1500k"],full=fl)
 if pt!="diagnostics":out["ledger_ram_budget"]=chk_ldg(lm,lis,lh3)
 out["checks_ok"]=all(v.get("ok",T_)for v in out.values()if IS(v,dict));out["elapsed_seconds"]=RD(time.perf_counter()-t0,3);print(json.dumps(out,indent=2))
 if not out["checks_ok"]:sys.exit(1)
if __name__=="__main__":main()
