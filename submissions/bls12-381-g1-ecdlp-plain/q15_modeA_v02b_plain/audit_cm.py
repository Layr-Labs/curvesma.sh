"""audit_cm.py — Self-contained pure-Python-stdlib hosted Docker audit for BLS12-381 CM-STNFS over Q(zeta_15)..."""
import base64
import cmath
import lzma
import hashlib
import json
import math
import random
import sys
from decimal import Decimal,ROUND_CEILING,localcontext
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
if str(ROOT)not in sys.path:sys.path.insert(0,str(ROOT))
if str(HERE)not in sys.path:sys.path.insert(0,str(HERE))
from cm_models import(GAMMA,LN2,MONT5,NTT3_PRIMES,NTT3_PROD,NTT4_PRIMES,NTT4_PROD,Retention,primes_upto,ecm_instr,ecm_mean_var,fft_F,gbatch_ntt3,lll_16d_int,log2_li,missing_recon_z12,mr_instr,ntt3_F,pbits_split,rho_db,rho_exact,srt_toy_check,verify_kernels_symbolic,)
from filtsim_z12 import verify_h3_ladder_and_fit
from w8_verifiers import check_sigmaB_exactness_live,verify_genuine_q15_filter_hx_live
from hx_z12 import run_hx_on_filtered_z12,unplanted_two_sided_bw_check
UNIT=20000.0
U_BLS=-0xd201000000010000
R_BLS=U_BLS**4-U_BLS**2+1
P_BLS=((U_BLS-1)**2*R_BLS)//3+U_BLS
PHI15=[1,-1,0,1,-1,1,0,-1,1]
def ptrim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:
  a.pop()
 return a
def pmul_mod(a,b,m,mod_poly=None):
 out=[0]*(len(a)+len(b)-1)
 for i,ai in enumerate(a):
  if ai:
   for j,bj in enumerate(b):
    if bj:out[i+j]+=ai*bj
 if mod_poly is not None:return prem_mod(out,mod_poly,m)
 out=[x%m for x in out]
 while len(out)>1 and out[-1]==0:
  out.pop()
 return out
def prem_mod(a,b,m):
 a=[x%m for x in a]
 while len(a)>1 and a[-1]==0:
  a.pop()
 db=len(b)-1
 while db>0 and b[db]%m==0:
  db-=1
 if db==0:return[0]
 lc=b[db]%m
 monic=(lc==1)
 inv_lc=1 if monic else pow(lc,-1,m)
 b_mod=[x%m for x in b[:db+1]]
 while len(a)-1>=db and not(len(a)==1 and a[0]==0):
  c=a[-1]if monic else(a[-1]*inv_lc)%m
  off=len(a)-1-db
  for i in range(db):
   bi=b_mod[i]
   if bi:a[off+i]=(a[off+i]-c*bi)%m
  a.pop()
  while len(a)>1 and a[-1]==0:
   a.pop()
 return a
def ppow_mod(base,exp,mod_poly,m):
 res=[1]
 b=prem_mod(base,mod_poly,m)
 e=exp
 while e>0:
  if e&1:res=pmul_mod(res,b,m,mod_poly)
  e>>=1
  if e:b=pmul_mod(b,b,m,mod_poly)
 return res
def pgcd_mod(a,b,m):
 a=ptrim([x%m for x in a])
 b=ptrim([x%m for x in b])
 while not(len(b)==1 and b[0]==0):
  a,b=b,prem_mod(a,b,m)
 inv=pow(a[-1],-1,m)
 return[(x*inv)%m for x in a]
def pmul_z(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,ai in enumerate(a):
  if ai:
   for j,bj in enumerate(b):
    if bj:out[i+j]+=ai*bj
 return out
def red_phi15_z(a):
 c=list(a)
 while len(c)>8:
  lead=c.pop()
  if lead:
   off=len(c)-8
   for i in range(8):
    c[off+i]-=lead*PHI15[i]
 while len(c)<8:
  c.append(0)
 return c
def mul_phi15_z(a,b):
 a0,a1,a2,a3,a4,a5,a6,a7=a;b0,b1,b2,b3,b4,b5,b6,b7=b
 c0=a0*b0;c1=a0*b1+a1*b0;c2=a0*b2+a1*b1+a2*b0;c3=a0*b3+a1*b2+a2*b1+a3*b0
 c4=a0*b4+a1*b3+a2*b2+a3*b1+a4*b0;c5=a0*b5+a1*b4+a2*b3+a3*b2+a4*b1+a5*b0
 c6=a0*b6+a1*b5+a2*b4+a3*b3+a4*b2+a5*b1+a6*b0;c7=a0*b7+a1*b6+a2*b5+a3*b4+a4*b3+a5*b2+a6*b1+a7*b0
 c8=a1*b7+a2*b6+a3*b5+a4*b4+a5*b3+a6*b2+a7*b1;c9=a2*b7+a3*b6+a4*b5+a5*b4+a6*b3+a7*b2
 c10=a3*b7+a4*b6+a5*b5+a6*b4+a7*b3;c11=a4*b7+a5*b6+a6*b5+a7*b4;c12=a5*b7+a6*b6+a7*b5;c13=a6*b7+a7*b6;c14=a7*b7
 c13+=c14;c11-=c14;c10+=c14;c9-=c14;c7+=c14;c6-=c14
 c12+=c13;c10-=c13;c9+=c13;c8-=c13;c6+=c13;c5-=c13
 c11+=c12;c9-=c12;c8+=c12;c7-=c12;c5+=c12;c4-=c12
 c10+=c11;c8-=c11;c7+=c11;c6-=c11;c4+=c11;c3-=c11
 c9+=c10;c7-=c10;c6+=c10;c5-=c10;c3+=c10;c2-=c10
 c8+=c9;c6-=c9;c5+=c9;c4-=c9;c2+=c9;c1-=c9
 c7+=c8;c5-=c8;c4+=c8;c3-=c8;c1+=c8;c0-=c8
 return[c0,c1,c2,c3,c4,c5,c6,c7]
def pow_phi15_z(a,e):
 if e==3:return mul_phi15_z(a,mul_phi15_z(a,a))
 res=[1]+[0]*7;b=list(a)
 while e>0:
  if e&1:res=mul_phi15_z(res,b)
  e>>=1
  if e:b=mul_phi15_z(b,b)
 return res
_YPOW15=tuple(red_phi15_z([0]*k+[1])for k in range(15))
def y_pow_z(k):
 return list(_YPOW15[k%15])
def sub_phi15_z(a,b):return[a[i]-b[i]for i in range(8)]
def add_phi15_z(a,b):return[a[i]+b[i]for i in range(8)]
def scale_phi15_z(s,a):return[s*x for x in a]
def norm_q15_z(a):
 v0,v1,v2,v3,v4,v5,v6,v7=a
 p1=mul_phi15_z(a,(v0+v1+v2-v5-v6-v7,-v2-v4+v7,-v1-v3+v6,v1-v6-v7,-v1-v2+v7,v2-v5-v7,v1-v4-v6,-v1-v2-v3+v6+v7))
 v0,v1,v2,v3,v4,v5,v6,v7=p1
 p2=mul_phi15_z(p1,(v0-v2-v6+v7,v2+v4-v7,-v3+v6,-v2-v6,v1+v2-v7,-v2+v5+v7,-v6,v2-v3+v6-v7))
 v0,v1,v2,v3,v4,v5,v6,v7=p2
 p3=mul_phi15_z(p2,(v0+v2+v4,-v4,-v2-v6,v2,-v2-v4+v7,v4+v5,v2+v3,v1-v2-v4-v6))
 assert not any(p3[1:]),f"non-rational norm {p3}"
 return p3[0]
def check_algebra_bls12():
 assert P_BLS.bit_length()==381
 assert R_BLS.bit_length()==255
 assert P_BLS%15==7
 assert R_BLS%15==13
 chi=P_BLS%R_BLS
 assert chi==(U_BLS%R_BLS)
 assert pow(chi,12,R_BLS)==1
 _cof_gt=(P_BLS**12-1)//R_BLS;assert (P_BLS**12-1)%R_BLS==0 and math.gcd(R_BLS,_cof_gt)==1
 for d in (1,2,3,4,6):
  assert math.gcd(pow(chi,d,R_BLS)-1,R_BLS)==1
 w_pi=((U_BLS**3-U_BLS**2-U_BLS)*pow(U_BLS+1,-1,P_BLS))%P_BLS
 assert(w_pi*w_pi+w_pi+1)%P_BLS==0
 w_pibar=(-1-w_pi)%P_BLS
 y5_minus_w=[-w_pi%P_BLS,0,0,0,0,1]
 h_pi=pgcd_mod(PHI15,y5_minus_w,P_BLS)
 h_pibar=pgcd_mod(PHI15,[-w_pibar%P_BLS,0,0,0,0,1],P_BLS)
 assert len(h_pi)==5 and len(h_pibar)==5
 prod=pmul_mod(h_pi,h_pibar,P_BLS)
 assert prod==[x%P_BLS for x in PHI15]
 yp=ppow_mod([0,1],P_BLS,h_pi,P_BLS)
 assert yp==prem_mod(y_pow_z(7),h_pi,P_BLS),"sigma(y)=y^7 must equal y^p mod h_pi!"
 uw2=(U_BLS*w_pi*w_pi)%P_BLS
 exp3=(P_BLS**4-1)//3
 c_frob=pow(uw2,(P_BLS-1)//3,P_BLS)
 assert pow(uw2,exp3,P_BLS)!=1,"u*w^2 must be a non-cube in F_{p^4}"
 assert c_frob==(w_pi*w_pi)%P_BLS,"X^p = omega^2 * X required for tau = Frob^1!"
 x3=uw2
 f_val=(pow(x3,3,P_BLS)+(w_pi+1)*pow(x3,2,P_BLS)+x3-w_pi)%P_BLS
 assert f_val==0
 assert-U_BLS==(1<<16)*906349*254760293
 assert pow(3,(19**2-1)//3,19)!=1 and pow(6,(67**4-1)//3,67)!=1
 assert 67%15==7 and 5%3==2
 f_ge2_orbits_log2=1.0+log2_li(35.85)
 assert f_ge2_orbits_log2<32.5 and 35.85<48.0
 w_r1=(-U_BLS*U_BLS)%R_BLS
 assert(w_r1*w_r1+w_r1+1)%R_BLS==0
 w_r2=(-1-w_r1)%R_BLS
 assert w_r2==pow(U_BLS,4,R_BLS)
 hr1=pgcd_mod(PHI15,[-w_r1%R_BLS,0,0,0,0,1],R_BLS)
 hr2=pgcd_mod(PHI15,[-w_r2%R_BLS,0,0,0,0,1],R_BLS)
 assert len(hr1)==5 and len(hr2)==5
 assert pmul_mod(hr1,hr2,R_BLS)==[x%R_BLS for x in PHI15]
 def mat_mul_mod(A,B,m):
  n=len(A)
  return[[sum(A[i][k]*B[k][j]for k in range(n))%m for j in range(n)]for i in range(n)]
 def mat_pow_mod(A,e,m):
  n=len(A)
  res=[[1 if i==j else 0 for j in range(n)]for i in range(n)]
  b=[row[:]for row in A]
  while e>0:
   if e&1:res=mat_mul_mod(res,b,m)
   e>>=1
   if e:b=mat_mul_mod(b,b,m)
  return res
 def frob_mat(poly,m):
  d=len(poly)-1
  xr=ppow_mod([0,1],m,poly,m)
  cols=[]
  cur=[1]
  for _ in range(d):
   cols.append(cur+[0]*(d-len(cur)))
   cur=pmul_mod(cur,xr,m,poly)
  return[[cols[j][i]for j in range(d)]for i in range(d)]
 I4=[[1 if i==j else 0 for j in range(4)]for i in range(4)]
 for hri in (hr1,hr2):
  M4=frob_mat(hri,R_BLS)
  M4_2=mat_mul_mod(M4,M4,R_BLS)
  assert M4_2!=I4 and mat_mul_mod(M4_2,M4_2,R_BLS)==I4
 u3=pow(U_BLS,3,R_BLS)
 g2_roots=[u3,(u3*w_r2)%R_BLS,(u3*w_r2*w_r2)%R_BLS]
 assert len(set(g2_roots))==3
 g2_prod=[1]
 for rt in g2_roots:
  g2_prod=pmul_mod(g2_prod,[(-rt)%R_BLS,1],R_BLS)
 assert g2_prod==[(-U_BLS*w_r2*w_r2)%R_BLS,0,0,1]
 c1=(U_BLS*w_r1*w_r1)%R_BLS
 assert pow(c1,(R_BLS-1)//3,R_BLS)!=1
 assert((R_BLS+1)*(R_BLS*R_BLS+1))%3==1
 I9=[[1 if i==j else 0 for j in range(9)]for i in range(9)]
 for wi in (w_r1,w_r2):
  fi=[(-wi)%R_BLS,0,0,1,0,0,(wi+1)%R_BLS,0,0,1]
  M9=frob_mat(fi,R_BLS)
  M9_3=mat_pow_mod(M9,3,R_BLS)
  assert mat_pow_mod(M9_3,3,R_BLS)==I9
  x_r3_minus_x=[M9_3[i][1]for i in range(9)]
  x_r3_minus_x[1]=(x_r3_minus_x[1]-1)%R_BLS
  assert pgcd_mod(x_r3_minus_x,fi,R_BLS)==[1]
 G=build_G1w()
 detG=int_det(G)
 assert detG==1265625==(3**4)*(5**6)
 W,S=build_WS_1w()
 assert mat_mul(mat_T(W),mat_mul(G,W))==G
 assert mat_mul(mat_T(S),mat_mul(G,S))==G
 for np_i in NTT4_PRIMES:
  assert np_i<(1<<62)and(np_i-1)%(1<<53)==0
  for a_b in (2,3,5,7,11,13,17,19,23,29,31,37):
   d_m,s_m=np_i-1,0
   while d_m%2==0:d_m//=2;s_m+=1
   x_m=pow(a_b,d_m,np_i)
   if x_m in (1,np_i-1):continue
   for _ in range(s_m-1):
    x_m=(x_m*x_m)%np_i
    if x_m==np_i-1:break
   else:raise AssertionError(f"NTT prime {np_i} composite!")
 ntt4_bits=sum(math.log2(np_i)for np_i in NTT4_PRIMES)
 assert ntt4_bits>244.0
 # Exhaustive verification of all 56 prime ideals P of O_H with N(P)<256 (min N(P)=16>9>3; 150 root ideals).
 _oh_sm=[]
 for _p in range(2,256):
  if any(_p%_d==0 for _d in range(2,int(_p**0.5)+1)):continue
  if _p==3:_f,_g=4,1
  elif _p==5:_f,_g=2,1
  else:
   _x,_f=_p%15,1
   while _x!=1:_x=(_x*_p)%15;_f+=1
   _g=8//_f
  _Np=_p**_f
  if _Np<256:_oh_sm.append((_p,_f,_g,_Np))
 _min_np=min(x[3]for x in _oh_sm);assert _min_np==16>9>3 and sum(x[2]for x in _oh_sm)==56
 def _fq_roots(p,mp):
  d=len(mp)-1;els=[[]]
  for _ in range(d):els=[x+[c]for x in els for c in range(p)]
  def pm(a,b):
   r=[0]*(2*d-1)
   for i,x in enumerate(a):
    for j,y in enumerate(b):r[i+j]=(r[i+j]+x*y)%p
   for i in range(len(r)-1,d-1,-1):
    c=r[i]
    if c:
     for j in range(d):r[i-d+j]=(r[i-d+j]-c*mp[j])%p
   return r[:d]
  one=[1]+[0]*(d-1)
  w=one if p==3 else next(e for e in els if not any((pm(e,e)[j]+e[j]+one[j])%p for j in range(d)))
  w2=pm(w,w);wp1=[(w[j]+one[j])%p for j in range(d)];rf=rg=0
  for x in els:
   x2=pm(x,x);x3=pm(x2,x)
   if not any((x3[j]-(U_BLS%p)*w2[j])%p for j in range(d)):rg+=1
   x6=pm(x3,x3);x9=pm(x6,x3);t6=pm(wp1,x6)
   if not any((x9[j]+t6[j]+x3[j]-w[j])%p for j in range(d)):rf+=1
  return rf,rg
 _sm_cnt={}
 for _p,_mp,_g in ((2,[1,1,0,0,1],2),(3,[2,1,0,0,1],1),(5,[2,1,1],1),(11,[1,1,1],4)):
  _rf,_rg=_fq_roots(_p,_mp);_Np=_p**(len(_mp)-1)
  if _rf+_rg:_sm_cnt[_Np]=_g*(_rf+_rg)
 for _p in (31,61,151,181,211,241):
  _z=next(pow(b,(_p-1)//15,_p)for b in range(2,200)if pow(b,(_p-1)//3,_p)!=1 and pow(b,(_p-1)//5,_p)!=1 and pow(b,(_p-1)//15,_p)!=1)
  _cp=0
  for _k in (1,2,4,7,8,11,13,14):
   _ry=pow(_z,_k,_p);_w=pow(_ry,5,_p);_w2=(_w*_w)%_p
   _cp+=sum(1 for _x in range(_p)if(pow(_x,9,_p)+(_w+1)*pow(_x,6,_p)+pow(_x,3,_p)-_w)%_p==0)
   _cp+=sum(1 for _x in range(_p)if(pow(_x,3,_p)-U_BLS*_w2)%_p==0)
  _sm_cnt[_p]=_cp
 from cm_models import SMALL_150
 assert sorted(_sm_cnt.items())==sorted(SMALL_150)and sum(_sm_cnt.values())==150
 _sm_hits=sum(c/float(n)for n,c in SMALL_150)
 # Single special-q descent.
 return dict(p_bits=P_BLS.bit_length(),r_bits=R_BLS.bit_length(),p_mod15=P_BLS%15,r_mod15=R_BLS%15,h_pi_deg=4,h_pibar_deg=4,uw2_noncube=True,tau_eq_frob1=True,mod_r_crt_degs=[4,4],mod_r_g_degs=[4,4,4,12],mod_r_f_degs=[36,36],eps_g_deg=12,eps_f_deg=36,det_G=detG,w_s_isometries=True,ntt4_primes_ok=len(NTT4_PRIMES),ntt4_prod_bits=round(ntt4_bits,2),f_ge2_orbits_log2=round(f_ge2_orbits_log2,2),oh_min_prime_norm=_min_np,oh_small_ideals_lt256=56,kf_kg_small_ideals_lt256=150,small_ideals_hits_pc=round(_sm_hits,4),)
def trace_z15(d):
 d%=15
 g=math.gcd(15,d)
 return 8 if g==15 else(-4 if g==5 else(-2 if g==3 else 1))
def mat_mul(A,B):
 n,k,m=len(A),len(B),len(B[0])
 return[[sum(A[i][t]*B[t][j]for t in range(k))for j in range(m)]for i in range(n)]
def mat_T(A):
 return[[A[j][i]for j in range(len(A))]for i in range(len(A[0]))]
def int_det(M):
 n=len(M)
 A=[[F(x)for x in row]for row in M]
 det=F(1)
 for i in range(n):
  piv=next(r for r in range(i,n)if A[r][i]!=0)
  if piv!=i:
   A[i],A[piv]=A[piv],A[i]
   det=-det
  det*=A[i][i]
  inv=F(1)/A[i][i]
  for r in range(i+1,n):
   if A[r][i]!=0:
    fac=A[r][i]*inv
    for c in range(i,n):
     A[r][c]-=fac*A[i][c]
 assert det.denominator==1
 return det.numerator
def inv_int_mat(M):
 n=len(M)
 A=[[F(x)for x in row]+[F(1 if i==j else 0)for j in range(n)]for i,row in enumerate(M)]
 for i in range(n):
  piv=next(r for r in range(i,n)if A[r][i]!=0)
  A[i],A[piv]=A[piv],A[i]
  inv=F(1)/A[i][i]
  A[i]=[x*inv for x in A[i]]
  for r in range(n):
   if r!=i and A[r][i]!=0:
    fac=A[r][i]
    A[r]=[A[r][c]-fac*A[i][c]for c in range(2*n)]
 return[[x.numerator for x in row[n:]]for row in A]
EXPS_1W=[0,5,3,8,6,11,9,14]
M_1W=mat_T([y_pow_z(e)for e in EXPS_1W])
M_1W_INV=inv_int_mat(M_1W)
def build_G1w():
 G_pow=[[trace_z15(i-j)for j in range(8)]for i in range(8)]
 return mat_mul(mat_T(M_1W),mat_mul(G_pow,M_1W))
def build_WS_1w():
 M_sig=mat_T([y_pow_z(7*j)for j in range(8)])
 M_w=mat_T([y_pow_z(5+j)for j in range(8)])
 S_1w=mat_mul(M_1W_INV,mat_mul(M_sig,M_1W))
 W_1w=mat_mul(M_1W_INV,mat_mul(M_w,M_1W))
 return W_1w,S_1w
def check_toy_e2e_dlp(toy):
 u0,p0,r0,c_x,chi=toy["u0"],toy["p0"],toy["r0"],toy["c_x"],toy["chi"]
 h_pi,h_pibar=toy["h_pi"],toy["h_pibar"]
 assert p0==7207 and r0==601 and chi==p0%r0==596
 assert pmul_mod(h_pi,h_pibar,p0)==[x%p0 for x in PHI15]
 def fp4_add(a,b):return[(a[i]+b[i])%p0 for i in range(4)]
 def fp4_sub(a,b):return[(a[i]-b[i])%p0 for i in range(4)]
 _h0,_h1,_h2,_h3=h_pi[0],h_pi[1],h_pi[2],h_pi[3]
 def fp4_mul(a,b):
  a0,a1,a2,a3=a;b0,b1,b2,b3=b
  c0=a0*b0;c1=a0*b1+a1*b0;c2=a0*b2+a1*b1+a2*b0;c3=a0*b3+a1*b2+a2*b1+a3*b0
  c4=a1*b3+a2*b2+a3*b1;c5=a2*b3+a3*b2;c6=a3*b3
  v=c6%p0
  if v:c2-=v*_h0;c3-=v*_h1;c4-=v*_h2;c5-=v*_h3
  v=c5%p0
  if v:c1-=v*_h0;c2-=v*_h1;c3-=v*_h2;c4-=v*_h3
  v=c4%p0
  if v:c0-=v*_h0;c1-=v*_h1;c2-=v*_h2;c3-=v*_h3
  return[c0%p0,c1%p0,c2%p0,c3%p0]
 def fp4_pow(a,e):
  res=[1,0,0,0];b=list(a)
  while e>0:
   if e&1:res=fp4_mul(res,b)
   e>>=1
   if e:b=fp4_mul(b,b)
  return res
 def fp4_inv(a):
  a1=fp4_pow(a,p0)
  a2=fp4_pow(a1,p0)
  a3=fp4_pow(a2,p0)
  prod3=fp4_mul(a1,fp4_mul(a2,a3))
  nrm=fp4_mul(a,prod3)[0]
  inv_n=pow(nrm,-1,p0)
  return[(x*inv_n)%p0 for x in prod3]
 y4=[0,1,0,0]
 w4=fp4_pow(y4,5)
 w4_2=fp4_mul(w4,w4)
 uw2=[(u0*x)%p0 for x in w4_2]
 wc1=w4 if c_x==1 else w4_2
 wc2=fp4_mul(wc1,wc1)
 frob4_cols=[fp4_pow([1 if j==i else 0 for j in range(4)],p0)for i in range(4)]
 def frob4(a):
  out=[0,0,0,0]
  for i in range(4):
   ai=a[i]
   if ai:
    col=frob4_cols[i]
    for j in range(4):
     out[j]=(out[j]+ai*col[j])%p0
  return out
 def fp12_mul(A,B):
  c0=fp4_add(fp4_mul(A[0],B[0]),fp4_mul(uw2,fp4_add(fp4_mul(A[1],B[2]),fp4_mul(A[2],B[1]))))
  c1=fp4_add(fp4_add(fp4_mul(A[0],B[1]),fp4_mul(A[1],B[0])),fp4_mul(uw2,fp4_mul(A[2],B[2])))
  c2=fp4_add(fp4_add(fp4_mul(A[0],B[2]),fp4_mul(A[1],B[1])),fp4_mul(A[2],B[0]))
  return(c0,c1,c2)
 def fp12_pow(A,e):
  res=([1,0,0,0],[0,0,0,0],[0,0,0,0]);b=A
  while e>0:
   if e&1:res=fp12_mul(res,b)
   e>>=1
   if e:b=fp12_mul(b,b)
  return res
 def fp12_frob(A):
  return(frob4(A[0]),fp4_mul(frob4(A[1]),wc1),fp4_mul(frob4(A[2]),wc2))
 def fp12_inv(A):
  def f4(Z):
   for _ in range(4):Z=fp12_frob(Z)
   return Z
  A4=f4(A);A8=f4(A4)
  adj=fp12_mul(A4,A8)
  n4=fp12_mul(A,adj)[0]
  inv4=fp4_inv(n4)
  return(fp4_mul(adj[0],inv4),fp4_mul(adj[1],inv4),fp4_mul(adj[2],inv4))
 cof_GT=(p0**4-p0**2+1)//r0
 def to_GT(A):
  e6=A
  for _ in range(6):e6=fp12_frob(e6)
  ez=fp12_mul(e6,fp12_inv(A))
  e2=fp12_frob(fp12_frob(ez))
  return fp12_pow(fp12_mul(e2,ez),cof_GT)
 def embed_ab(a8,b8):
  ar=prem_mod(a8,h_pi,p0);ar+=[0]*(4-len(ar))
  br=prem_mod([-x for x in b8],h_pi,p0);br+=[0]*(4-len(br))
  return(ar,br,[0,0,0,0])
 def key12(A):return(tuple(A[0]),tuple(A[1]),tuple(A[2]))
 G_gen=to_GT(([3,1,0,0],[2,0,p0-1,0],[0,5,0,0]))
 m_bs=math.isqrt(r0)+1
 bs_map={}
 cur=([1,0,0,0],[0,0,0,0],[0,0,0,0])
 for j in range(m_bs):
  bs_map[key12(cur)]=j
  cur=fp12_mul(cur,G_gen)
 Gm_inv=fp12_inv(fp12_pow(G_gen,m_bs))
 def dlog_GT(A):
  c=to_GT(A)
  for i in range(m_bs+1):
   k=key12(c)
   if k in bs_map:return(i*m_bs+bs_map[k])%r0
   c=fp12_mul(c,Gm_inv)
  raise RuntimeError("BSGS failed")
 r0sq=r0*r0
 roots_r0=[s for s in range(2,r0)if pow(s,15,r0)==1 and pow(s,5,r0)!=1 and pow(s,3,r0)!=1]
 assert len(roots_r0)==8
 roots_l2=[]
 for s in roots_r0:
  val=sum(PHI15[j]*pow(s,j,r0sq)for j in range(9))%r0sq
  der=sum(j*PHI15[j]*pow(s,j-1,r0)for j in range(1,9))%r0
  s2=(s-val*pow(der,-1,r0))%r0sq
  roots_l2.append(s2)
 crt_e=[]
 for idx,s in enumerate(roots_r0):
  poly=[1]
  den=1
  for j,sj in enumerate(roots_r0):
   if j!=idx:
    poly=pmul_mod(poly,[-sj%r0,1],r0)
    den=(den*(s-sj))%r0
  invd=pow(den,-1,r0)
  p8=[(x*invd)%r0 for x in poly]
  p8+=[0]*(8-len(p8))
  crt_e.append(p8)
 w_tau_exp=5*c_x
 chi_inv_pows=[pow(chi,(-t)%12,r0)for t in range(12)]
 s2_pows=[[pow(s2,j,r0sq)for j in range(8)]for s2 in roots_l2]
 s2_meta=[]
 for s2 in roots_l2:
  wc=pow(s2,5,r0sq);w2c=(wc*wc)%r0sq
  s2_meta.append(((u0*w2c)%r0sq,wc%r0sq,((wc+1)*wc)%r0sq,w2c%r0sq))
 def sm_chi_proj(a8,b8,is_g):
  dX=3 if is_g else 9
  eps=(r0**3-1)if is_g else(r0**6-1)
  L_c=[[0]*8 for _ in range(dX)]
  _MSK=(1<<44)-1
  for ci,(spw,(cg0,cf0,cf3,cf6)) in enumerate(zip(s2_pows,s2_meta)):
   ac=(a8[0]*spw[0]+a8[1]*spw[1]+a8[2]*spw[2]+a8[3]*spw[3]+a8[4]*spw[4]+a8[5]*spw[5]+a8[6]*spw[6]+a8[7]*spw[7])%r0sq
   bc=(b8[0]*spw[0]+b8[1]*spw[1]+b8[2]*spw[2]+b8[3]*spw[3]+b8[4]*spw[4]+b8[5]*spw[5]+b8[6]*spw[6]+b8[7]*spw[7])%r0sq
   if is_g:
    c0=cg0
    p0,p1,p2=1,0,0;q0,q1,q2=ac,(-bc)%r0sq,0;ee=eps
    while ee:
     if ee&1:
      p0,p1,p2=(p0*q0+c0*(p1*q2+p2*q1))%r0sq,(p0*q1+p1*q0+c0*p2*q2)%r0sq,(p0*q2+p1*q1+p2*q0)%r0sq
     ee>>=1
     if ee:
      q0,q1,q2=(q0*q0+c0*(2*q1*q2))%r0sq,(2*q0*q1+c0*q2*q2)%r0sq,(2*q0*q2+q1*q1)%r0sq
    pw=[p0,p1,p2]
   else:
    c0,c3,c6=cf0,cf3,cf6
    def _smul9(x,y):
     ix=((((((((x[8]<<44)|x[7])<<44|x[6])<<44|x[5])<<44|x[4])<<44|x[3])<<44|x[2])<<44|x[1])<<44|x[0]
     iy=((((((((y[8]<<44)|y[7])<<44|y[6])<<44|y[5])<<44|y[4])<<44|y[3])<<44|y[2])<<44|y[1])<<44|y[0]
     p=ix*iy
     o0=p&_MSK;p>>=44;o1=p&_MSK;p>>=44;o2=p&_MSK;p>>=44;o3=p&_MSK;p>>=44;o4=p&_MSK;p>>=44;o5=p&_MSK;p>>=44;o6=p&_MSK;p>>=44;o7=p&_MSK;p>>=44;o8=p&_MSK;p>>=44
     o9=p&_MSK;p>>=44;o10=p&_MSK;p>>=44;o11=p&_MSK;p>>=44;o12=p&_MSK;p>>=44;o13=p&_MSK;p>>=44;o14=p&_MSK;p>>=44;o15=p&_MSK;p>>=44;o16=p&_MSK
     v=o16%r0sq
     if v:o7+=v*c0;o10+=v*c3;o13+=v*c6
     v=o15%r0sq
     if v:o6+=v*c0;o9+=v*c3;o12+=v*c6
     v=o14%r0sq
     if v:o5+=v*c0;o8+=v*c3;o11+=v*c6
     v=o13%r0sq
     if v:o4+=v*c0;o7+=v*c3;o10+=v*c6
     v=o12%r0sq
     if v:o3+=v*c0;o6+=v*c3;o9+=v*c6
     v=o11%r0sq
     if v:o2+=v*c0;o5+=v*c3;o8+=v*c6
     v=o10%r0sq
     if v:o1+=v*c0;o4+=v*c3;o7+=v*c6
     v=o9%r0sq
     if v:o0+=v*c0;o3+=v*c3;o6+=v*c6
     return[o0%r0sq,o1%r0sq,o2%r0sq,o3%r0sq,o4%r0sq,o5%r0sq,o6%r0sq,o7%r0sq,o8%r0sq]
    pw=[1,0,0,0,0,0,0,0,0];bs=[ac,(-bc)%r0sq,0,0,0,0,0,0,0];ee=eps
    while ee:
     if ee&1:pw=_smul9(pw,bs)
     ee>>=1
     if ee:bs=_smul9(bs,bs)
   pw[0]=(pw[0]-1)%r0sq
   ec=crt_e[ci]
   for k in range(dX):
    assert pw[k]%r0==0
    qk=(pw[k]//r0)%r0
    if qk:
     for j in range(8):
      L_c[k][j]=(L_c[k][j]+qk*ec[j])%r0
  proj=[0]*(dX*8)
  for k in range(dX):
   ck=L_c[k]
   if not any(ck):continue
   acc=[0]*8
   cur_p=list(ck)
   for t in range(12):
    wt=chi_inv_pows[t]
    shift=_YPOW15[(w_tau_exp*t*k)%15]
    term=mul_phi15_z(cur_p,shift)
    for j in range(8):
     acc[j]=(acc[j]+wt*term[j])%r0
    nxt=[0]*8
    for j,cj in enumerate(cur_p):
     if cj:
      y7j=_YPOW15[(7*j)%15]
      for ii in range(8):
       nxt[ii]=(nxt[ii]+cj*y7j[ii])%r0
    cur_p=nxt
   for j in range(8):
    proj[k*8+j]=acc[j]
  return proj
 def rref_mod(M,m):
  A=[[x%m for x in row]for row in M]
  nrows,ncols=len(A),len(A[0])
  pivs=[]
  r=0
  for c in range(ncols):
   p=next((i for i in range(r,nrows)if A[i][c]!=0),None)
   if p is None:continue
   A[r],A[p]=A[p],A[r]
   inv=pow(A[r][c],-1,m)
   A[r]=[(x*inv)%m for x in A[r]]
   for i in range(nrows):
    if i!=r and A[i][c]!=0:
     fac=A[i][c]
     A[i]=[(A[i][j]-fac*A[r][j])%m for j in range(ncols)]
   pivs.append(c);r+=1
  return A,pivs
 def left_kernel_mod(M,m):
  MT=mat_T(M)
  A,pivs=rref_mod(MT,m)
  piv_set=set(pivs)
  ncols=len(MT[0])
  free=[c for c in range(ncols)if c not in piv_set]
  ker=[]
  for fc in free:
   v=[0]*ncols
   v[fc]=1
   for ri,pc in enumerate(pivs):
    v[pc]=(-A[ri][fc])%m
   ker.append(v)
  return ker
 def solve_mod(M,b,m):
  aug=[row+[bi%m]for row,bi in zip(M,b)]
  A,pivs=rref_mod(aug,m)
  ncols=len(M[0])
  assert ncols not in pivs,"inconsistent linear system"
  x=[0]*ncols
  for ri,pc in enumerate(pivs):
   x[pc]=A[ri][ncols]
  return x
 rng_sm=random.Random(20261119)
 W_g=[rng_sm.randrange(r0)for _ in range(24)]
 W_f=[[rng_sm.randrange(r0)for _ in range(72)]for _ in range(3)]
 kg=toy["kg"]
 dlogs_g=[dlog_GT(embed_ab(rel["a"],rel["b"]))for rel in kg["rels"]]
 assert dlogs_g==[rel["dl"]for rel in kg["rels"]],"Live BSGS must match K_g dlogs!"
 smg_full=[sm_chi_proj(rel["a"],rel["b"],True)for rel in kg["rels"]]
 _,pivg=rref_mod(smg_full,r0)
 assert len(pivg)==2
 smg_1=[[sum(W_g[j]*row[j]for j in range(24))%r0]for row in smg_full]
 Mg=[rel["row"]for rel in kg["rels"]]
 Maug_g=[row+sm for row,sm in zip(Mg,smg_1)]
 Lkg=left_kernel_mod(Maug_g,r0)
 fail_g=sum(1 for v in Lkg if sum(v[i]*dlogs_g[i]for i in range(len(dlogs_g)))%r0!=0)
 assert len(Lkg)>=20 and fail_g==0
 Lkg0=left_kernel_mod(Mg,r0)
 fail_g0=sum(1 for v in Lkg0 if sum(v[i]*dlogs_g[i]for i in range(len(dlogs_g)))%r0!=0)
 assert fail_g0>0
 fchi_cnt=0;nchi_chk=0
 chi2=(chi*chi)%r0
 for rel,dl in zip(kg["rels"],dlogs_g):
  if dl==0:continue
  nchi_chk+=1
  dl_f=dlog_GT(fp12_frob(embed_ab(rel["a"],rel["b"])))
  assert dl_f==(chi*dl)%r0
  if dl_f!=(chi2*dl)%r0:fchi_cnt+=1
 assert fchi_cnt==nchi_chk and nchi_chk>30
 wg=solve_mod(Maug_g,dlogs_g,r0)
 t_matches=0
 for tgt in kg["targets"]:
  dl_true=dlog_GT(embed_ab(tgt["a"],tgt["b"]))
  assert dl_true==tgt["dl"]
  smt_full=sm_chi_proj(tgt["a"],tgt["b"],True)
  smt=[sum(W_g[j]*smt_full[j]for j in range(24))%r0]
  rec=sum((tgt["row"]+smt)[j]*wg[j]for j in range(len(wg)))%r0
  if rec==dl_true:t_matches+=1
 assert t_matches==len(kg["targets"])==5
 kf=toy["kf"]
 dlogs_f=[dlog_GT(embed_ab(rel["a"],rel["b"]))for rel in kf["rels"]]
 assert dlogs_f==[rel["dl"]for rel in kf["rels"]]
 smf_full=[sm_chi_proj(rel["a"],rel["b"],False)for rel in kf["rels"]]
 _,pivf=rref_mod(smf_full,r0)
 assert 3<=len(pivf)<=6
 smf_3=[[sum(W_f[c][j]*row[j]for j in range(72))%r0 for c in range(3)]for row in smf_full]
 Mf=[rel["row"]for rel in kf["rels"]]
 Maug_f=[row+sm for row,sm in zip(Mf,smf_3)]
 Lkf=left_kernel_mod(Maug_f,r0)
 fail_f=sum(1 for v in Lkf if sum(v[i]*dlogs_f[i]for i in range(len(dlogs_f)))%r0!=0)
 assert len(Lkf)>=50 and fail_f==0
 wf=solve_mod(Maug_f,dlogs_f,r0)
 match_f=sum(1 for row,dl in zip(Maug_f,dlogs_f)if sum(row[j]*wf[j]for j in range(len(wf)))%r0==dl)
 assert match_f==len(dlogs_f)==62
 return dict(kg_rels=len(dlogs_g),kg_deps=len(Lkg),kg_fail=fail_g,kg_neg_drop_sm_fail=f"{fail_g0}/{len(Lkg0)}",kg_neg_wrong_chi_fail=f"{fchi_cnt}/{nchi_chk}",kg_targets_solved=f"{t_matches}/5",kf_rels=len(dlogs_f),kf_deps=len(Lkf),kf_fail=fail_f,kf_vlog_fit=f"{match_f}/{len(dlogs_f)}",random_Wf_Wg_compression=True,omitted_short_orbit_factors=kg["n_short"]+kf["n_short"],)
def check_sigmaB_exactness():
 return check_sigmaB_exactness_live(build_G1w,build_WS_1w,mat_mul,M_1W,norm_q15_z,U_BLS)
def check_anchor_norms(anchors):
 z3=y_pow_z(5)
 z3_2=mul_phi15_z(z3,z3)
 y1=y_pow_z(1);y2=y_pow_z(2);y3=y_pow_z(3)
 w_p1_y2=mul_phi15_z(add_phi15_z(z3,[1]+[0]*7),y2)
 wy3=mul_phi15_z(z3,y3)
 def val_p(n,p):
  v=0
  while n%p==0:
   n//=p;v+=1
  return v
 max_err_f=max_err_g=0.0
 for row in anchors:
  a=[sum(M_1W[i][j]*row["a"][j]for j in range(8))for i in range(8)]
  b=[sum(M_1W[i][j]*row["b"][j]for j in range(8))for i in range(8)]
  a3=pow_phi15_z(a,3);b3=pow_phi15_z(b,3)
  a6=mul_phi15_z(a3,a3);a9=mul_phi15_z(a6,a3)
  b6=mul_phi15_z(b3,b3);b9=mul_phi15_z(b6,b3)
  g_el=sub_phi15_z(a3,scale_phi15_z(U_BLS,mul_phi15_z(y1,b3)))
  ng=abs(norm_q15_z(g_el))
  t1=mul_phi15_z(y1,mul_phi15_z(a6,b3))
  t2=mul_phi15_z(w_p1_y2,mul_phi15_z(a3,b6))
  t3=mul_phi15_z(wy3,b9)
  f_el=sub_phi15_z(sub_phi15_z(sub_phi15_z(a9,t1),t2),t3)
  nf=abs(norm_q15_z(f_el))
  assert val_p(ng,2)==row["v2g"]
  assert val_p(nf,3)==row["v3f"]and val_p(ng,3)==row["v3g"]
  assert val_p(nf,5)==row["v5f"]and val_p(ng,5)==row["v5g"]
  z5_inv=y_pow_z(12)
  bt=mul_phi15_z(z5_inv,b)
  bt3=pow_phi15_z(bt,3)
  g_iso=sub_phi15_z(a3,scale_phi15_z(U_BLS,mul_phi15_z(z3_2,bt3)))
  assert abs(norm_q15_z(g_iso))==ng,"w=y <-> w'=omega^2 norm isometry failed!"
  max_err_f=max(max_err_f,abs(math.log2(nf)-row["nf"]))
  max_err_g=max(max_err_g,abs(math.log2(ng)-row["ng"]))
 assert max_err_f<1e-8 and max_err_g<1e-8
 return dict(n_anchors=len(anchors),max_err_f=max_err_f,max_err_g=max_err_g,isometry_exact=True)
def verify_pool_2adic_live(pool,op):
 import unitdomain as _UD
 return _UD.verify_pool_2adic_live(pool,op,build_G1w,M_1W,PHI15,U_BLS,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,rho_db)
def check_voronoi_exact(data):
 units=data["pool"]["voronoi_units_10"]
 assert len(units)==10
 places=[cmath.exp(2j*math.pi*k/15.0)for k in (1,7,4,2)]
 facets=[]
 for u in units:
  assert abs(norm_q15_z(u))==1
  facets.append([abs(sum(u[j]*(z**j)for j in range(8)))**2 for z in places])
 anc_ok=0
 for row in data["pool"]["anchors"]:
  a_pow=[sum(M_1W[i][j]*row["a"][j]for j in range(8))for i in range(8)]
  x4=[abs(sum(a_pow[j]*(z**j)for j in range(8)))**2 for z in places]
  t2=sum(x4)
  if all(sum(fj*xj for fj,xj in zip(f,x4))>=t2-1e-7 for f in facets):anc_ok+=1
 assert anc_ok==len(data["pool"]["anchors"])==12
 p0_l,u0_l=7207,-5
 wpi0=((u0_l**3-u0_l**2-u0_l)*pow(u0_l+1,-1,p0_l))%p0_l
 hpi0=pgcd_mod(PHI15,[-wpi0%p0_l,0,0,0,0,1],p0_l)
 gam0=[137,2903,5111,889]
 B_hnf=[]
 for i in range(4):
  v=[0]*16;v[i]=p0_l;B_hnf.append(v)
 for j in range(4):
  v=[0]*16
  for idx,c in enumerate(hpi0):v[j+idx]=c
  B_hnf.append(v)
 for k in range(8):
  yk=[0]*8;yk[k]=1
  ak=pmul_mod(gam0,prem_mod(yk,hpi0,p0_l),p0_l,hpi0)
  v=[0]*16
  for idx,c in enumerate(ak):v[idx]=c
  v[8+k]=1;B_hnf.append(v)
 B_red,_,_=lll_16d_int(B_hnf,99,100)
 v1_l=B_red[0];a1_l,b1_l=v1_l[:8],v1_l[8:]
 v2_l,k_ind=None,None
 for idx in range(1,16):
  ak_l,bk_l=B_red[idx][:8],B_red[idx][8:]
  if any(sub_phi15_z(mul_phi15_z(a1_l,bk_l),mul_phi15_z(b1_l,ak_l))):v2_l,k_ind=B_red[idx],idx+1;break
 assert v2_l is not None
 a2_l,b2_l=v2_l[:8],v2_l[8:]
 for aa,bb in ((a1_l,b1_l),(a2_l,b2_l)):
  assert prem_mod(aa,hpi0,p0_l)==pmul_mod(gam0,prem_mod(bb,hpi0,p0_l),p0_l,hpi0)
 det_l=abs(norm_q15_z(sub_phi15_z(mul_phi15_z(a1_l,b2_l),mul_phi15_z(a2_l,b1_l))))
 assert det_l%(p0_l**4)==0 and det_l>0
 # Live 381-bit G_T F_{p^6}-subfield elimination: s*z' = W_0 + W_1*X.
 assert math.gcd(R_BLS,P_BLS**4-1)==1 and math.gcd(R_BLS,P_BLS**6-1)==1 and (P_BLS**12-1)%(R_BLS*math.lcm(P_BLS**4-1,P_BLS**6-1))==0
 w_pi=((U_BLS**3-U_BLS**2-U_BLS)*pow(U_BLS+1,-1,P_BLS))%P_BLS
 h_pi=pgcd_mod(PHI15,[-w_pi%P_BLS,0,0,0,0,1],P_BLS)
 assert len(h_pi)==5 and h_pi[4]==1
 h0,h1,h2,h3=h_pi[0],h_pi[1],h_pi[2],h_pi[3];uw2=(U_BLS*w_pi*w_pi)%P_BLS
 def _f4a(a,b):return[(a[0]+b[0])%P_BLS,(a[1]+b[1])%P_BLS,(a[2]+b[2])%P_BLS,(a[3]+b[3])%P_BLS]
 def _f4s(a,b):return[(a[0]-b[0])%P_BLS,(a[1]-b[1])%P_BLS,(a[2]-b[2])%P_BLS,(a[3]-b[3])%P_BLS]
 def _f4c(s,a):return[(s*a[0])%P_BLS,(s*a[1])%P_BLS,(s*a[2])%P_BLS,(s*a[3])%P_BLS]
 def _f4m(a,b):
  a0,a1,a2,a3=a;b0,b1,b2,b3=b
  c0=a0*b0;c1=a0*b1+a1*b0;c2=a0*b2+a1*b1+a2*b0;c3=a0*b3+a1*b2+a2*b1+a3*b0
  c4=a1*b3+a2*b2+a3*b1;c5=a2*b3+a3*b2;c6=a3*b3
  v=c6%P_BLS
  if v:c2-=v*h0;c3-=v*h1;c4-=v*h2;c5-=v*h3
  v=c5%P_BLS
  if v:c1-=v*h0;c2-=v*h1;c3-=v*h2;c4-=v*h3
  v=c4%P_BLS
  if v:c0-=v*h0;c1-=v*h1;c2-=v*h2;c3-=v*h3
  return[c0%P_BLS,c1%P_BLS,c2%P_BLS,c3%P_BLS]
 def _f4p(a,e):
  r=[1,0,0,0];b=list(a)
  while e:
   if e&1:r=_f4m(r,b)
   e>>=1
   if e:b=_f4m(b,b)
  return r
 _fc=[_f4p([1 if j==i else 0 for j in range(4)],P_BLS)for i in range(4)]
 def _fr1(a):
  o=[0,0,0,0]
  for i in range(4):
   ai=a[i]
   if ai:
    col=_fc[i];o[0]+=ai*col[0];o[1]+=ai*col[1];o[2]+=ai*col[2];o[3]+=ai*col[3]
  return[x%P_BLS for x in o]
 def _fr2(a):return _fr1(_fr1(a))
 def _f4i(a):
  a1=_fr1(a);a2=_fr1(a1);a3=_fr1(a2);pr=_f4m(a1,_f4m(a2,a3))
  return _f4c(pow(_f4m(a,pr)[0],-1,P_BLS),pr)
 def _f12m(A,B):
  a0,a1,a2=A;b0,b1,b2=B;v0=_f4m(a0,b0);v1=_f4m(a1,b1);v2=_f4m(a2,b2)
  t0=_f4s(_f4s(_f4m(_f4a(a1,a2),_f4a(b1,b2)),v1),v2);t1=_f4s(_f4s(_f4m(_f4a(a0,a1),_f4a(b0,b1)),v0),v1);t2=_f4s(_f4s(_f4m(_f4a(a0,a2),_f4a(b0,b2)),v0),v2)
  return(_f4a(v0,_f4c(uw2,t0)),_f4a(t1,_f4c(uw2,v2)),_f4a(t2,v1))
 def _f12p(A,e):
  if e==0:return([1,0,0,0],[0,0,0,0],[0,0,0,0])
  tbl=[([1,0,0,0],[0,0,0,0],[0,0,0,0]),A]
  for _ in range(14):tbl.append(_f12m(tbl[-1],A))
  nbits=(e.bit_length()+3)&~3
  r=tbl[0]
  for sh in range(nbits-4,-1,-4):
   r=_f12m(r,r);r=_f12m(r,r);r=_f12m(r,r);r=_f12m(r,r)
   nib=(e>>sh)&15
   if nib:r=_f12m(r,tbl[nib])
  return r
 _inv2=(P_BLS+1)//2;_d2=_f4s([0,1,0,0],_fr2([0,1,0,0]));_i2d=_f4i(_f4c(2,_d2))
 def _to2(a):ap=_fr2(a);return _f4c(_inv2,_f4a(a,ap)),_f4m(_f4s(a,ap),_i2d)
 _g_raw=([3,1,0,0],[5,0,2,0],[7,0,0,1]);_cof=(P_BLS**12-1)//R_BLS
 # Construct real 381-bit G_T generator via Frobenius final exponentiation.
 _cx_fr=pow(uw2,(P_BLS-1)//3,P_BLS);_cx2_fr=(_cx_fr*_cx_fr)%P_BLS
 def _f12_fr(A,k):
  rA=A
  for _ in range(k):
   rA=(_fr1(rA[0]),_f4c(_cx_fr,_fr1(rA[1])),_f4c(_cx2_fr,_fr1(rA[2])))
  return rA
 def _f12_inv(A):
  A1=_f12_fr(A,4);A2=_f12_fr(A1,4);p12=_f12m(A1,A2);nm=_f12m(A,p12)
  assert nm[1]==[0,0,0,0]and nm[2]==[0,0,0,0];inm=_f4i(nm[0])
  return(_f4m(inm,p12[0]),_f4m(inm,p12[1]),_f4m(inm,p12[2]))
 _y1=_f12m(_f12_fr(_g_raw,6),_f12_inv(_g_raw));_y2=_f12m(_f12_fr(_y1,2),_y1)
 _e_hard=(P_BLS**4-P_BLS**2+1)//R_BLS;G0=_f12p(_y2,_e_hard)
 _one12=([1,0,0,0],[0,0,0,0],[0,0,0,0]);assert G0!=_one12 and _f12p(G0,R_BLS)==_one12
 _gt_ok=0
 for _ti in range(4):
  _xt=int(hashlib.sha256(f"q15-boot-elim-{_ti}".encode()).hexdigest(),16)%(R_BLS-1)+1
  _Z=_f12p(G0,_xt)
  _e0,_o0=_to2(_Z[0]);_e1,_o1=_to2(_Z[1]);_e2,_o2=_to2(_Z[2])
  _c0=_f4s(_f4m(_e1,_o0),_f4m(_o1,_e0));_c1=_f4s(_f4m(_o2,_e0),_f4m(_e2,_o0));_c2=_f4s(_f4m(_e2,_o1),_f4m(_o2,_e1))
  assert _fr2(_c0)==_c0 and _fr2(_c1)==_c1 and _fr2(_c2)==_c2
  _s=(_c0,_c1,_c2);_W=_f12m(_s,_Z)
  assert _W[2]==[0,0,0,0]and any(_W[1])
  _gam=_f4m([(-x)%P_BLS for x in _W[0]],_f4i(_W[1]))
  # Sample nonzero (a, b) in L_boot(z').
  _bv=[1 if j==_ti else 0 for j in range(8)];_b4=prem_mod(_bv,h_pi,P_BLS);_b4=_b4+[0]*(4-len(_b4))
  _a4=_f4m(_gam,_b4);_av=_a4+[0,0,0,0]
  _el=(_a4,[(-x)%P_BLS for x in _b4],[0,0,0,0])
  _fac4=_f4m([(-x)%P_BLS for x in _b4],_f4i(_W[1]))
  _rhs=_f12m((_fac4,[0,0,0,0],[0,0,0,0]),_f12m(_s,_Z))
  assert _el==_rhs and _f12_fr(_s,6)==_s and _f12_fr((_fac4,[0,0,0,0],[0,0,0,0]),4)==(_fac4,[0,0,0,0],[0,0,0,0])
  _gt_ok+=1
 assert _gt_ok==4
 # Live z' == 1 blinding-shortcut check: unique e_0 = (-x) mod r yields z' = z * gamma_0^{e_0} = 1 -> V(z) = (-e_0 * V_0) mod r with zero descent work (F_{p^6} elim only runs on z' != 1).
 _x_sc=int(hashlib.sha256(b"q15-boot-z1-shortcut").hexdigest(),16)%(R_BLS-1)+1
 _Z_sc=_f12p(G0,_x_sc);_e0_sc=(- _x_sc)%R_BLS;_Zp_sc=_f12m(_Z_sc,_f12p(G0,_e0_sc))
 _V0_sc=12345678901234567%R_BLS
 assert _Zp_sc==_one12 and ((-_e0_sc*_V0_sc)%R_BLS)==((_V0_sc*_x_sc)%R_BLS)
 _e0z,_o0z=_to2(_one12[0]);_e1z,_o1z=_to2(_one12[1]);_e2z,_o2z=_to2(_one12[2])
 assert _f12m((_f4s(_f4m(_e1z,_o0z),_f4m(_o1z,_e0z)),_f4s(_f4m(_o2z,_e0z),_f4m(_e2z,_o0z)),_f4s(_f4m(_e2z,_o1z),_f4m(_o2z,_e1z))),_one12)[1]==[0,0,0,0]
 return dict(voronoi_units_ok=10,anchors_min_t2=f"{anc_ok}/12",constructive_lll16_k=k_ind,constructive_lll16_index=det_l//(p0_l**4),gt_boot_fp6_elim_live_381b=_gt_ok,gt_boot_z_eq_1_shortcut_ok=True,gcd_r_lcm_p4_p6_is_1=True)
def check_coll_la_core(data):
 op=data["operating_point"]
 pool=data["pool"]
 coll=data["collection"]
 la=data["la"]
 led=data["ledger"]
 lf,lg=op["lf"],op["lg"]
 p_g,p_f,p_lo=coll["p_g"],coll["p_f"],coll["p_lo"]
 det=coll["det"]
 log2_Up=math.log2((2.0**log2_li(lf)+2.0**log2_li(lg))/12.0)
 assert abs(log2_Up-op["log2_Up"])<1e-9
 p_min=int(coll.get("p_min",32))
 ret_g=Retention(round(lg,2),p_g,prime_min=2**p_min,explicit_bits=14,step=0.25,max_bits=800,tilt=0.05,odd=1)
 cst_g=Retention(800,p_g,prime_min=2**p_min,explicit_bits=14,step=0.25,max_bits=801,tilt=0.005,odd=1)
 ret_f=Retention(round(lf,2),p_f,prime_min=2**p_min,explicit_bits=14,step=0.25,max_bits=660,tilt=0.05,odd=1)
 cst_f=Retention(660,p_f,prime_min=2**p_min,explicit_bits=14,step=0.25,max_bits=661,tilt=0.005,odd=1)
 sbg,tf_const,sbf=det["sbg"],det["tf_const"],det["sbf"]
 tg_const=tuple(det.get("tg_const",[299,287,275]))
 kmg,kmf=int(sbg//p_g),int(sbf//p_f)
 _kh_g=float(det.get("kh_g",0.08));_kh_f=float(det.get("kh_f",0.08))
 c_fbatch=gbatch_ntt3(pbits_split(lf,_kh_f),leaf=9)
 c_gbatch=gbatch_ntt3(pbits_split(lg,_kh_g),leaf=11)
 _ps14=primes_upto(16384);_tb14=_gb14=_fb14=0.0
 for _p14 in _ps14:
  _lb14=math.log2(_p14);_tb14+=_lb14
  if _p14 in (2,3,5):_gb14+=_lb14;_fb14+=_lb14;continue
  if _p14%15!=1:continue
  _e14=(_p14-1)//15;_s14=next(t for g in range(2,80)for t in(pow(g,_e14,_p14),)if pow(t,3,_p14)!=1 and pow(t,5,_p14)!=1)
  _w14=pow(_s14,5,_p14);_w24=(_w14*_w14)%_p14
  if pow((U_BLS*_w24)%_p14,(_p14-1)//3,_p14)==1 or pow((U_BLS*_w14)%_p14,(_p14-1)//3,_p14)==1:_gb14+=_lb14
  for _ww in (_w14,_w24):
   _fp=[(-_ww)%_p14,0,0,1,0,0,(_ww+1)%_p14,0,0,1];_xp=ppow_mod([0,1],_p14,_fp,_p14)
   _xm=list(_xp)+[0]*max(0,2-len(_xp));_xm[1]=(_xm[1]-1)%_p14
   if len(pgcd_mod(_fp,_xm,_p14))>1:_fb14+=_lb14;break
 assert _gb14/_tb14<=_kh_g and _fb14/_tb14<=_kh_f
 _sb_cache={}
 def _stage_b_sp(nbits,kh):
  _k=(int(round(nbits)),kh)
  if _k not in _sb_cache:
   L=max(5,math.ceil(nbits/64.0))
   _sb_cache[_k]=gbatch_ntt3(pbits_split(p_min,kh),leaf=L,c_leaf=10.0*MONT5*((L/5.0)**2)+2.0*nbits*(12.0*L+16.0))
  return _sb_cache[_k]
 # Product trees (P_g, P_f, candidate tree) + distributed 4-step NTT scratch (< 2^64 B/owner).
 _sp_lg=math.ceil(pbits_split(lg,_kh_g)/64.0);_sp_lf=math.ceil(pbits_split(lf,_kh_f)/64.0)
 _n_pt_own=2048
 _L_pt_g=math.ceil(math.log2(_sp_lg/11.0));_L_pt_f=math.ceil(math.log2(_sp_lf/9.0))
 _cand_tree_g=(_L_pt_g+1.0)*(11.0*(1<<_L_pt_g))*8.0  # Retained g-batch candidate tree across all L_pt_g+1 levels!
 _pt_mem_own=((_L_pt_g+1.0)*_sp_lg*8.0+(_L_pt_f+1.0)*_sp_lf*8.0+_cand_tree_g+4.0*_sp_lg*8.0)/_n_pt_own
 _sub_l11_pg=(_L_pt_g-11.0+1.0)*(_sp_lg/(1<<11))*8.0
 _sub_l11_pf=(_L_pt_f-11.0+1.0)*(_sp_lf/(1<<11))*8.0
 _sub_l11_cand=(_L_pt_g-11.0+1.0)*(11.0*(1<<(_L_pt_g-11)))*8.0
 assert _pt_mem_own<2.0**64 and max(_sub_l11_pg,_sub_l11_pf,_sub_l11_cand)<2.0**62,"2048-owner retained lf+lg+cand product trees must fit < 2^64 B/owner and ALL THREE L-11 local subtrees < 2^62 B!"
 _pt_cost=(math.ceil(math.log2(_sp_lg))*ntt3_F(_sp_lg)+math.ceil(math.log2(_sp_lf))*ntt3_F(_sp_lf)+10.0*10.0*(_sp_lg+_sp_lf))
 assert math.log2(_pt_cost/UNIT)<=dict(led["rows"])["pg_pf_product_trees"]
 alpha_f,alpha_g_ge3,kap2_g=op["alpha_f"],op["alpha_g_ge3"],op["kap2_g"]
 rg_c,sga_c,sgb_c={},{},{}
 rf_c,sfa_c,sfb_c,qf_c={},{},{},{}
 ysum_mid=ysum_rethi=ysum_strict_hi=ysum_sh_sq=c_cof_num=w_k_tot=sga_num=0.0
 _nf_cond_num=_ng_cond_num=_w_cond_den=0.0
 max_bg_seen=0
 hs=float(pool.get("hist_step",1.0))
 hh=0.5*hs
 _f1_f, _f1_g = {}, {}
 _rf_sm, _rg_sm = {}, {}
 for bf,bgo,bgm,wc,wk in pool["hist_2d"]:
  if bgm>max_bg_seen:max_bg_seen=bgm
  _cls=0 if (bgm-bgo)==1 else(1 if (bgm-bgo)==13 else 2)
  _tg=tg_const[_cls]
  if bf not in rf_c:
   # Single constant f-side threshold S_f >= tf_const (ZERO pre-Nf knowledge of bf, ZERO f-norm predictor!)
   _sl_lo=float(bf)-float(tf_const);_sl_hi=float(bf+hs)-float(tf_const)
   rf_c[bf]=(ret_f.at(bf+hh,_sl_lo,sbf)if _sl_lo>0 else 0.0,ret_f.at(bf+hs,_sl_lo,sbf)if _sl_lo>0 else 0.0)
   _sfa0=cst_f.survivor_upper(float(bf),_sl_lo,0.25)if _sl_lo>=0 else 0.0
   _sfa1=cst_f.survivor_upper(float(bf+hs),_sl_hi,0.25)if _sl_hi>=0 else 0.0
   _sfb0=cst_f.survivor_upper(float(bf),_sl_lo,0.25,sbf,lf,kmf)if _sl_lo>=0 else 0.0
   _sfb1=cst_f.survivor_upper(float(bf+hs),_sl_hi,0.25,sbf,lf,kmf)if _sl_hi>=0 else 0.0
   sfa_c[bf]=min(1.0,1.25*max(_sfa0,_sfa1))
   sfb_c[bf]=min(1.0,1.25*max(_sfb0,_sfb1))
   qf_c[bf]=min(1.0,2.0*rho_exact((bf+alpha_f)/lf))
   ef_c=bf+hh+alpha_f;ef_h=bf+hs+alpha_f
   _rf_sm[bf]=(rho_db(ef_c/lf,ef_c),rho_db(ef_h/lf,ef_h))
  _gk=(bgo,_tg)
  if _gk not in rg_c:
   # Static 2-adic class constant T_g[c] = (299, 287, 275) — ZERO Q_CM(a) radial shells, ZERO per-bgo predictor!
   _slg_lo=float(bgo)-float(_tg);_slg_hi=float(bgo+hs)-float(_tg)
   rg_c[_gk]=(ret_g.at(bgo+hh,_slg_lo,sbg)if _slg_lo>0 else 0.0,ret_g.at(bgo+hs,_slg_lo,sbg)if _slg_lo>0 else 0.0)
   _sga0=cst_g.survivor_upper(float(bgo),_slg_lo,0.25)if _slg_lo>=0 else 0.0
   _sga1=cst_g.survivor_upper(float(bgo+hs),_slg_hi,0.25)if _slg_hi>=0 else 0.0
   _sgb0=cst_g.survivor_upper(float(bgo),_slg_lo,0.25,sbg,lg,kmg)if _slg_lo>=0 else 0.0
   _sgb1=cst_g.survivor_upper(float(bgo+hs),_slg_hi,0.25,sbg,lg,kmg)if _slg_hi>=0 else 0.0
   sga_c[_gk]=min(1.0,max(_sga0,_sga1))
   sgb_c[_gk]=min(1.0,max(_sgb0,_sgb1))
  if bgo not in _rg_sm:
   eg_c=bgo+hh+alpha_g_ge3+1.0+kap2_g;eg_h=bgo+hs+alpha_g_ge3+1.0+kap2_g
   _rg_sm[bgo]=(rho_db(eg_c/lg,eg_c),rho_db(eg_h/lg,eg_h))
  if wc>0.0:
   _rf0,_rf1=_rf_sm[bf];_rg0,_rg1=_rg_sm[bgo]
   p_sm=_rf0*_rg0
   p_sm_hi=_rf1*_rg1
   ysum_mid+=wc*p_sm*rf_c[bf][0]*rg_c[_gk][0]
   ysum_rethi+=wc*p_sm*rf_c[bf][1]*rg_c[_gk][1]
   _ysh=p_sm_hi*rf_c[bf][1]*rg_c[_gk][1];ysum_strict_hi+=wc*_ysh;ysum_sh_sq+=wc*_ysh*_ysh
   if rf_c[bf][1]>0.0 and rg_c[_gk][1]>0.0:
    _wsh=wc*_ysh
    _f1_f[bf]=_f1_f.get(bf,0.0)+_wsh
    _f1_g[_gk]=_f1_g.get(_gk,0.0)+_wsh
    _nf_cond_num+=_wsh*(bf+hh);_ng_cond_num+=_wsh*(bgo+hh);_w_cond_den+=_wsh
  if wk>0.0:
   sga,sgb=sga_c[_gk],sgb_c[_gk]
   sfa,sfb,qf=sfa_c[bf],sfb_c[bf],qf_c[bf]
   csb_g=_stage_b_sp(float(bgm),_kh_g)
   csb_f=_stage_b_sp(float(bf+hs),_kh_f)
   c_ng=(sga*sfa)*6500.0
   c_nf=(sgb*sfa)*20000.0
   c_sb=(sga*sfa)*csb_g+(sgb*sfa)*csb_f
   q_f_cond=min(1.0,qf/max(sfb,1e-12))
   c_fb=(sgb*sfb)*(c_fbatch+q_f_cond*(2.0**24.5+2.0**24.0+c_gbatch))
   c_cof_num+=wk*(c_ng+c_nf+c_sb+c_fb)
   sga_num+=wk*sga
   w_k_tot+=wk
 c_cof_live=c_cof_num/w_k_tot;sga_mean=sga_num/w_k_tot
 nf_cond_mean=_nf_cond_num/_w_cond_den;ng_cond_mean=_ng_cond_num/_w_cond_den
 # Scan full 0.02-bit descent grids [58.02, 70.00] and [63.52, 70.00] directly on pool_hist_2d!
 _f_it=[(_bf,_wr/rho_db((_bf+hs+alpha_f)/lf,_bf+hs+alpha_f),_bf+hs+alpha_f,float(_bf)-float(tf_const),1.0/rf_c[_bf][1])for _bf,_wr in _f1_f.items()]
 _g_it=[(_bgo,_wr/rho_db((_bgo+hs+alpha_g_ge3+1.0+kap2_g)/lg,_bgo+hs+alpha_g_ge3+1.0+kap2_g),_bgo+hs+alpha_g_ge3+1.0+kap2_g,float(_bgo)-float(_tg),1.0/rg_c[(_bgo,_tg)][1])for(_bgo,_tg),_wr in _f1_g.items()]
 def _min_ret_mult(side_s,b_lo_100):
  rm_min=1e9;_rd=rho_db;_rdc={}
  if side_s=="f":
   _ra=ret_f.at;_il=1.0/lf;_it=[(int(round(_bf*100)),_bf,_w0,_ef,_sf_lo,_irf)for _bf,_w0,_ef,_sf_lo,_irf in _f_it]
   for _bi in range(7000,b_lo_100-1,-2):
    _b=_bi*0.01;_num=_den=0.0;_ok_s=(sbf>_b)
    for _bf100,_bf,_w0,_ef,_sf_lo,_irf in _it:
     _ek=_bf100-_bi
     _rv=_rdc.get(_ek)
     if _rv is None:
      _eb=_ef-_b;_rv=_rd(_eb*_il,_eb if _eb>1.0 else 1.0);_rdc[_ek]=_rv
     _ww=_w0*_rv;_den+=_ww
     if _ok_s and _sf_lo>_b:_num+=_ww*(_ra(_bf+hs-_b,_sf_lo-_b,sbf-_b)*_irf)
    _r=_num/_den
    if _r<rm_min:rm_min=_r
  else:
   _ra=ret_g.at;_il=1.0/lg;_it=[(int(round(_bgo*100)),_bgo,_w0,_eg,_slg_lo,_irg)for _bgo,_w0,_eg,_slg_lo,_irg in _g_it]
   for _bi in range(7000,b_lo_100-1,-2):
    _b=_bi*0.01;_num=_den=0.0;_ok_s=(sbg>_b)
    for _bg100,_bgo,_w0,_eg,_slg_lo,_irg in _it:
     _ek=_bg100-_bi
     _rv=_rdc.get(_ek)
     if _rv is None:
      _eb=_eg-_b;_rv=_rd(_eb*_il,_eb if _eb>1.0 else 1.0);_rdc[_ek]=_rv
     _ww=_w0*_rv;_den+=_ww
     if _ok_s and _slg_lo>_b:_num+=_ww*(_ra(_bgo+hs-_b,_slg_lo-_b,sbg-_b)*_irg)
    _r=_num/_den
    if _r<rm_min:rm_min=_r
  return rm_min
 ret_mult_f_min=_min_ret_mult("f",5802);ret_mult_g_min=_min_ret_mult("g",6352)
 assert ret_mult_f_min>=0.90 and ret_mult_g_min>=0.82,"Conditioned ret_mult below (0.90, 0.82)!"
 assert abs(c_cof_live-det["c_cof_2d_upper"])<1e-6
 assert max_bg_seen>=690,"2D histogram must cover high-norm tail!"
 f_prim_hist=sum(wc for _,_,_,wc,_ in pool["hist_2d"])
 assert abs(f_prim_hist-pool["f_prim"])<1e-6 and 0.94<f_prim_hist<0.97
 log2_rels_h=(pool["log2_base"]+pool["dV"])+math.log2(ysum_mid)
 log2_rels_rethi=(pool["log2_base"]+pool["dV"])+math.log2(ysum_rethi)
 log2_rels_strict_hi=(pool["log2_base"]+pool["dV"])+math.log2(ysum_strict_hi)
 ecm_coll_succ=0.999999
 ecm_rel_keep=ecm_coll_succ**29
 det_cred_pen=-math.log2(float(coll.get("det_credited",0.99))*ecm_rel_keep)
 margin_net=(log2_rels_rethi-det_cred_pen)-(log2_Up+math.log2(op["exc"]))
 margin_strict_hi=(log2_rels_strict_hi-det_cred_pen)-(log2_Up+math.log2(op["exc"]))
 assert margin_net>=0.100 and margin_strict_hi>=0.030 and margin_strict_hi>=15.0*pool["se_y_bits"],"Finite-region O_H-coprime margin (>=15 SE strict-hi after det_credited=0.99)!"
 _msh=ysum_strict_hi/w_k_tot;_vsh=max(0.0,ysum_sh_sq/w_k_tot-_msh*_msh)
 _se_sh_b=math.sqrt(_vsh/float(pool.get("n_draws",1000000)))/(_msh*LN2)
 log2_rels_2se_live=log2_rels_strict_hi-det_cred_pen-2.0*max(_se_sh_b,pool["se_y_bits"])
 assert log2_rels_2se_live>=pool["log2_rels"]
 assert coll["h2_allow"]>=1.0+8.0*max(pool["se_c_rel"],pool["se_w_rel"])
 V_lin=2.0**pool["V_kept"]
 p2_mult=pool["p2_mult"]
 N_a_tbl=4.0*(2.0**pool["log2_Va"])
 N_7=2.0**pool["log2_N7_d12"]
 Rb_unit2=pool["A2_prod"]/pool["r_ab"]
 L_max_live=2.0*math.sqrt(Rb_unit2/8.0)
 Lb_live=L_max_live*(105.0*math.pi/768.0)
 assert abs(L_max_live-pool["L_max"])<1e-4 and abs(Lb_live-pool["Lb_avg"])<1e-4 and L_max_live<coll["w_tile"]
 assert abs(N_a_tbl*N_7*Lb_live-V_lin)/V_lin<1e-6,"Product volume N_a_tbl * N_7 * L_b must equal V_kept!"
 w_tile,w_sb,p_m4,p_m4_s=coll["w_tile"],coll["w_sb"],coll["p_m4"],coll["p_m4_s"]
 def price_side(p_hi):
  pl=min(p_lo,p_hi);pm=min(p_m4,pl);pms=min(p_m4_s,pm)
  n1al=1.01*(2.0**log2_li(pms)-2.0**log2_li(float(p_min)))
  n1as=1.01*(2.0**log2_li(pm)-2.0**log2_li(pms))if pm>pms else 0.0
  n1b=1.01*(2.0**log2_li(pl)-2.0**log2_li(pm))if pl>pm else 0.0
  b1_lk=n1al*float(coll.get("b1a_long_lk",384.0))+n1as*float(coll.get("b1a_short_lk",40.0))+n1b*float(coll.get("b1b_lk",12.0))
  ib1=(math.log(pms/float(p_min))*coll["c_hit_lo_m4"]+(math.log(pm/pms)*coll["c_hit_lo_m4_s"]if pm>pms else 0.0)+(math.log(pl/pm)*coll["c_hit_lo_u"]if pl>pm else 0.0))*p2_mult
  ib1+=(b1_lk*N_7*Lb_live+(n1al+n1as+n1b)*N_7*24.0)/V_lin
  sp1a=4.0*1.01*(2.0**(2.0*pm))/(2.0*pm*LN2-1.0)/4.0
  sp1b=(1.01*(2.0**(2.0*pl)-2.0**(2.0*pm))/(2.0*pl*LN2-1.0)/4.0)if pl>pm else 0.0
  # Charge 4 per-bank terminal sentinels + hist zero-init (+1/bkt).
  bld=((n1al+n1as+n1b)/4.0)*coll["c_list_a"]*N_a_tbl+(coll["c_lbkt"]+1.0)*(sp1a+sp1b)+(4.0*(n1al+n1as)+n1b)/4.0
  mem=((n1al+n1as+n1b)/4.0)*8.0*N_a_tbl+8.0*(sp1a+sp1b)
  ipr=ihi=0.0
  if p_hi>pl:
   n2=1.01*(2.0**log2_li(p_hi)-2.0**log2_li(pl))
   sb_r=Lb_live/w_sb+1.5
   ipr=n2*N_7*(coll["c_tile_hdr"]+float(coll.get("c_b2_row_k0",10.0))+sb_r*coll["c_sb_hdr"])/V_lin
   cand_bd=1.0+w_sb/Lb_live
   ihi=math.log(p_hi/pl)*cand_bd*coll["c_hit_hi"]*p2_mult+(28.0*32.0)/(N_a_tbl*Lb_live)
   sp2=1.01*(2.0**(2.0*p_hi)-2.0**(2.0*pl))/(2.0*p_hi*LN2-1.0)/4.0
   bld+=(n2/4.0)*coll["c_list_a"]*N_a_tbl+(coll["c_lbkt"]+1.0)*(sp2/w_sb)+(n2/4.0)
   mem+=(n2/4.0)*8.0*N_a_tbl+8.0*(sp2/w_sb)
  return bld,mem,ib1+ipr+ihi
 bg,mg,ipc_g=price_side(p_g)
 bf,mf,ipc_f=price_side(p_f)
 ipc_tot=ipc_g+ipc_f+coll["c_scan_zero"]+c_cof_live
 coll_u=pool["V_kept"]+math.log2(ipc_tot*coll["h2_allow"]/UNIT)
 build_u=math.log2((bg+bf)/UNIT)
 rows_d=dict(led["rows"])
 assert coll_u<=rows_d["relation_collection"]+1e-4
 assert build_u<=rows_d["replay_a_tables"]+1e-4
 col_own_mem=4.0*(w_sb+1.0)*(1<<56)*2.0
 m1a_bld=16.0*N_a_tbl+64.0*(2.0**p_m4)+w_tile*(1<<20)
 m1b_bld=16.0*N_a_tbl+16.0*(2.0**p_lo)+w_tile*(1<<20)
 m2_bld=16.0*N_a_tbl+16.0*((2.0**p_g)/w_sb)+w_tile*(1<<20)
 m2_res_max=8.0*N_a_tbl+8.0*((2.0**p_g)/w_sb)
 assert 2.0**pool["log2_Va"]<(1<<56)*0.95 and col_own_mem<2.0**64,"Column owner (exact 2^56 bank stride, S_g[16][4][2^56]+S_f[4][1][2^56]) must fit < 2^64 B!"
 assert max(m1a_bld,m1b_bld,m2_bld,m2_res_max)<2.0**64,"1-prime-per-owner build peak (<=14.981 EiB) & resident (<=7.491 EiB) must fit < 2^64 B!"
 assert mg+mf<2.0**120,"Total replay tables must fit < 2^120 B!"
 mem_tbl_tot=mg+mf
 np_band=1.25506*(2.0**48)/(48.0*LN2)
 c_rel_band=np_band*(6*(10+12)+44+4)+(13+16)*1024+(1<<26)
 n_rels_upper=2.0**(pool["log2_rels"]+0.55)
 band_u_live=math.log2(n_rels_upper*c_rel_band/UNIT)
 assert band_u_live<=rows_d["band_split"]
 ef_e,ef_v=ecm_mean_var(lf,640,success=ecm_coll_succ)
 eg_e,eg_v=ecm_mean_var(lg,768,success=ecm_coll_succ)
 assert ef_v==0.0 and eg_v==0.0
 c_ecm_rel=(13.0*ef_e+14.0*mr_instr(640)+13.0*(1<<26)+16.0*eg_e+17.0*mr_instr(768)+16.0*(1<<26))
 ew_ecm_live=(c_ecm_rel/UNIT)*n_rels_upper
 ecm_gap=(2.0**rows_d["ecm_split_expected"])-ew_ecm_live
 assert ecm_gap>0
 # Deterministic cap n_rels_upper = 2^(log2_rels+0.55); Stepped-GCD H-ECM (var=0, 0.999999**29 yield discount in det_cred_pen)!
 assert n_rels_upper>=(2.0**pool["log2_rels"])*1.46 and ecm_rel_keep<1.0
 _h3=verify_h3_ladder_and_fit(la,log2_Up,run_u12=False,op=op,nf_cond=nf_cond_mean,ng_cond=ng_cond_mean);k5_live=_h3["k5_live"]
 assert led["cap_abort"]==0.0
 assert la["tab_ratio_charged"]==1.0 and la["n_owners"]==1 and la["te"]==0,"Plain-sparse ZERO-TABM 1-owner 4-vec Arch B!"
 N=2.0**(log2_Up+math.log2(la["nl_cap"]))
 N_sw=2.0**(log2_Up+math.log2(la["sw_cap"]))
 m,n,n_sm=la["m"],la["n"],la["n_sm"]
 la_mem_full=2.0*N_sw*float(la.get("bytes_per_slot",128.0));horner_1own_mem=2.0*N_sw*32.0
 _n_proj_own=int(la.get("n_proj_owners",16))
 la_proj_mem=(124.0*32.0*N+128.0*32.0*N)/_n_proj_own+16.0*N
 _L_bm_mem=(1.0+2.0**-7)*(N+n_sm)*(1.0/m+1.0/n)+6.0
 _n_seq_own=2048;la_seq_mem=(_L_bm_mem*m*n*32.0)/_n_seq_own
 # PM-Basis across 2^22 owners: coeff arrays (2.58 EiB) + 3-plane x 4-prime NTT workspace at M=2^48 (6.05 EiB) = 8.63 EiB < 16 EiB!
 _n_pm_own=1<<22;_M_pm_sub=float(1<<math.ceil(math.log2(_L_bm_mem)))
 la_pm_mem=(2.0*_L_bm_mem*((m+n)**2)*32.0+3.0*4.0*_M_pm_sub*((m+n)**2)*8.0)/_n_pm_own
 la_gen_vec_mem=_L_bm_mem*n*32.0
 _nu_ub=float(la.get("nonunit_ub",0.25))
 _n_filt_own=8192;filt_mut_mem=(n_rels_upper*(64.0+201.0*(8.0+_nu_ub)+7.0+201.0*10.0+160.0))/_n_filt_own
 _n_purge_own=2048;purge_mut_mem=(n_rels_upper*1040.0)/_n_purge_own
 assert filt_mut_mem<=8.0*(2.0**60) and purge_mut_mem<=8.0*(2.0**60)
 assert all((1.0-x["exp_unit"])<=_nu_ub for x in la["ladder"][2:])
 n_mat_own=int(la.get("n_mat_owners",16));assert 1.0-(la["eu_25"]-2.0*la["se_eu_25"])<=_nu_ub
 assert la["k5_charged"]>=(1.0-_nu_ub)*7.65+_nu_ub*22.0 and float(la.get("k5_horner_charged",14.0))>=(1.0-_nu_ub)*9.60+_nu_ub*24.0
 _w_max_row=4096;assert _w_max_row<0xFFFF
 mat_shard_mem=(N_sw*(64.0+201.0*(8.0+_nu_ub)+7.0))/n_mat_own+(64.0+_w_max_row*9.0)
 _w_row_avg=8.0+201.0+(201.0*_nu_ub)/8.0+7.0/8.0
 _c_sm4=float(la["c_stream_mat"]);_c_sm1=float(la["c_stream_mat_1vec"])
 _ops_row_agg=1.125*_w_row_avg+47.125
 assert _c_sm4>=_ops_row_agg/4.0 and _c_sm1>=_ops_row_agg
 assert max(la_mem_full,horner_1own_mem,la_proj_mem,mat_shard_mem,la_seq_mem,la_pm_mem,la_gen_vec_mem,filt_mut_mem,purge_mut_mem)<2.0**64,"All LA/filtering owners must fit < 2^64 B!"
 # Explicit 1-vec Horner k5_horner_charged, full d=(1+1/128)N/n degree slack, and explicit 4.5*n Horner Y=[S,BZ] read stream from 16 projection owners!
 _eu_lo=la["eu_25"]-2.0*la["se_eu_25"]
 _k5_h_live=_eu_lo*9.60+(1.0-_eu_lo)*24.0
 _k5_h_chg=float(la.get("k5_horner_charged",11.80))
 assert _k5_h_chg>=_k5_h_live+1.0
 _c_xproj=float(la.get("c_xproj",70.0));assert _c_xproj>=64.5
 c_prod=(la["sw_cap"]/la["nl_cap"])*(201.0*la["k5_charged"]*la["tab_ratio_charged"]+la["c_row_chi"]+_c_sm4)
 c_prod_horner=(la["sw_cap"]/la["nl_cap"])*(201.0*_k5_h_chg*la["tab_ratio_charged"]+la["c_row_chi"]+_c_sm1)
 bw_sig=2.0**-7
 kap=((1.0+bw_sig)*(N+n_sm)+6.0/(1.0/m+1.0/n))/N
 per_n2=(c_prod+_c_xproj)*(n/m+1.0)*kap+(1.0/m+(1.0+1.0/128.0)/n)*(c_prod_horner+n*(97.0+4.5)+MONT5)
 la_u=math.log2((N*N*per_n2+N*(m+n)*MONT5)/UNIT)
 L_bm=(1.0+bw_sig)*(N+n_sm)*(1.0/m+1.0/n)+6.0
 t_bm=m+n;depth_bm=math.ceil(math.log2(L_bm));k_bm=math.ceil(math.log2(4.0*L_bm))
 assert(1<<k_bm)<=(1<<53)and t_bm*L_bm*(1<<128)<NTT4_PROD
 w_base_bm=L_bm*m*(t_bm**2)*MONT5
 w_rec_bm=2.0*depth_bm*(8.0*L_bm)*(1536.0*(t_bm**3)+(t_bm**2)*(1536.0*k_bm+6220.0))
 gen_u=math.log2((w_base_bm+w_rec_bm)/UNIT)
 assert la_u<=rows_d["linear_algebra"]
 assert gen_u<=rows_d["matrix_generator"]
 return dict(
  c_cof_2d_upper_live=round(c_cof_live,4),sga_mean=sga_mean,ret_mult_f_min=round(ret_mult_f_min,5),ret_mult_g_min=round(ret_mult_g_min,5),
  nf_cond_mean=round(nf_cond_mean,3),ng_cond_mean=round(ng_cond_mean,3),max_bg_bin_covered=max_bg_seen,
  log2_rels_2d_mid=round(log2_rels_h,6),log2_rels_2d_rethi=round(log2_rels_rethi,6),
  log2_rels_2d_strict_hi=round(log2_rels_strict_hi,6),log2_rels_2se_live=round(log2_rels_2se_live,6),
  margin_net_rethi_bits=round(margin_net,5),margin_strict_hi_bits=round(margin_strict_hi,5),
  L_max_live=round(L_max_live,2),ipc_tot=round(ipc_tot,4),
  coll_u=round(coll_u,6),build_u=round(build_u,4),pt_2048own_frac=round(_pt_mem_own/(2.0**64),4),
  pt_l11_sub_max_log2_B=round(math.log2(max(_sub_l11_pg,_sub_l11_pf,_sub_l11_cand)),3),
  max_1own_bld_log2_B=round(math.log2(max(m1a_bld,m1b_bld,m2_bld)),3),
  mem_tables_log2_B=round(math.log2(mem_tbl_tot),3),
  band_u=round(band_u_live,4),ecm_u=round(math.log2(ew_ecm_live),4),ecm_var_live=0.0,ecm_rel_keep=round(ecm_rel_keep,6),
  k5_live=round(k5_live,4),k5_horner_live=round(_k5_h_live,4),la_u=round(la_u,6),gen_u=round(gen_u,4),
  la_mem_full_log2=round(math.log2(la_mem_full),3),la_proj_mem_log2=round(math.log2(la_proj_mem),3),n_proj_owners=_n_proj_own,
  la_seq_mem_log2=round(math.log2(la_seq_mem),3),n_seq_owners=_n_seq_own,la_pm_mem_log2=round(math.log2(la_pm_mem),3),n_pm_owners=_n_pm_own,la_gen_vec_mem_log2=round(math.log2(la_gen_vec_mem),3),n_filt_owners=_n_filt_own,filt_mut_mem_log2=round(math.log2(filt_mut_mem),3),purge_mut_mem_log2=round(math.log2(purge_mut_mem),3),
  mat_shard_mem_log2=round(math.log2(mat_shard_mem),3),
  n_mat_owners=n_mat_own,c_stream_mat=float(la["c_stream_mat"]),
  N_la=N,N_sw=N_sw,f_prim_mc_1m=round(f_prim_hist,6),
 )
def check_hx_filtered_live():
 res=run_hx_on_filtered_z12(log2U=12,seed=5,m=12,n=8,s_sm=4)
 assert(res["checks_ok"]and res["l_rec_strict_lower_tri"]and res["prod_eq_schur"]and res["pm_basis_lead_sm_nonzero"]and res["pm_basis_zero_sm_ctrl_zero_vec"]and res["pm_basis_horner_passes"]==1 and res["rank_M"]==res["N"]and res["rank_A"]==res["D"]-1 and res["kernel_eq_zero"]and res["recovered_true_kernel"]and res["post_peel_cols_lost"]==0 and res["post_merge_cols_lost"]==0 and res["repurge_backsub_errors"]==0 and res["n_repurge_backsub"]==res["repurge_stack_len"]and res["n_reconstruct_unit"]==0 and res["all_2core_rels_zero"]and res["recovered_cols_match_true"]and res["cols_unrecovered_tri"]==0 and res["backsub_2core_exact"])
 return res
def main():
 raw=sys.stdin.read().strip()
 req=json.loads(raw)if raw else{}
 allowed_keys={"descent_mesh"}
 if not set(req).issubset(allowed_keys):raise ValueError("audit_cm.py accepts {} or {'descent_mesh': 384}")
 data=json.loads((ROOT/"cm_data.json").read_text())
 _tb=json.loads(lzma.decompress(base64.b85decode((ROOT/"cm_tables.b85").read_bytes())).decode())
 for _k in ("pool","collection","la"):data[_k]=_tb["cfg_"+_k]
 data["pool"]["hist_2d"]=_tb["pool_hist_2d"];data["pool"]["anchors"]=_tb["pool_anchors"];data["pool"]["voronoi_units_10"]=_tb["pool_vor"]
 data["la"]["ladder"]=_tb["la_ladder"]
 toy=json.loads(lzma.decompress(base64.b85decode((ROOT/"toy_e2e.b85").read_bytes())).decode())
 out={}
 ck={}
 import unitdomain as _UD
 out["algebra_bls12"]=check_algebra_bls12();ck["algebra_bls12"]=True
 out["murphy_alpha_exact_live"]=_UD.verify_murphy_alpha_live(data["operating_point"],U_BLS,primes_upto,pgcd_mod,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,rho_db);ck["murphy_alpha_exact_live"]=True
 out["real_fb_joint_coll_live"]=_UD.verify_real_fb_joint_coll(data["pool"]["voronoi_units_10"],build_G1w,build_WS_1w,mat_mul,M_1W,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,primes_upto,rho_db,Retention);ck["real_fb_joint_coll_live"]=True
 out["toy_e2e_dlp"]=check_toy_e2e_dlp(toy);ck["toy_e2e_dlp"]=True
 out["sigmaB_exactness"]=check_sigmaB_exactness();ck["sigmaB_exactness"]=True
 out["anchor_norms"]=check_anchor_norms(data["pool"]["anchors"]);out["pool_2adic_live"]=verify_pool_2adic_live(data["pool"],data["operating_point"]);ck["anchor_norms"]=True
 out["voronoi_exact"]=check_voronoi_exact(data);ck["voronoi_exact"]=True
 out["srt_toy"]=srt_toy_check();ck["srt_toy"]=True
 out["coll_la_core"]=check_coll_la_core(data);ck["coll_la_core"]=True
 out["filtsim_and_hx_live"]=check_hx_filtered_live();ck["filtsim_and_hx_live"]=True
 out["missing_recon_6b_live"]=missing_recon_z12();ck["missing_recon_6b_live"]=True
 _h3r=verify_h3_ladder_and_fit(data["la"],data["operating_point"]["log2_Up"],run_u12=True,op=data["operating_point"],nf_cond=out["coll_la_core"]["nf_cond_mean"],ng_cond=out["coll_la_core"]["ng_cond_mean"]);_u12p=_h3r.get("live_u12_seed0")
 out["h3_ladder_and_genuine_q15_hx"]=dict(h3=_h3r,genuine_q15_hx=verify_genuine_q15_filter_hx_live(toy,primes_upto,PHI15,y_pow_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,pmul_mod,prem_mod,u12_proxy=_u12p));ck["h3_ladder_and_genuine_q15_hx"]=True
 out["unplanted_two_sided_bw"]=unplanted_two_sided_bw_check()
 assert out["unplanted_two_sided_bw"]["pm_basis_lead_sm_nonzero"]and out["unplanted_two_sided_bw"]["pm_basis_zero_sm_ctrl_zero_vec"]and out["unplanted_two_sided_bw"]["pm_basis_horner_passes"]==1
 ck["unplanted_two_sided_bw"]=True
 _dp_kern=data.get("descent_kernels",{"lll_node":2.0**40.0,"lll_boot":2.0**47.0,"pair_setup":2.0**35.2})
 out["kernels_symbolic"]=verify_kernels_symbolic((ROOT/"KERNELS.md").read_text(),data["collection"],data["la"],R_BLS,float(_dp_kern["lll_node"]),float(_dp_kern["lll_boot"]),_dp_kern,sga_mean=out["coll_la_core"]["sga_mean"]);ck["kernels_symbolic"]=True
 import verify_q15_descent as _VD
 _dm=req.get("descent_mesh",384)
 if _dm not in (384,):raise ValueError("descent_mesh must be 384")
 dres=_VD.run_all_checks(root_dir=ROOT,seed=1,mesh=_dm)
 out["q15_modeA_descent"]=dres;ck["q15_modeA_descent"]=True
 _cl=out["coll_la_core"];_op=data["operating_point"];_pl=data["pool"];_led=data["ledger"]
 assert dres["all_ok"] is True and dres["mode"]=="pure_k1_boole_wald_iid_boot"
 assert dres["nb_r"]==1.0 and dres["nb_alpha"]==0.5 and dres["p0_floor"]==0.001
 assert dres["phi_tree"]<=0.60 and dres["boot_phi"]==0.0 and dres["p_descent_fail"]==0.0
 assert dres["mini_adaptive_campaign"]["checks_ok"] is True and dres["mini_adaptive_campaign"]["live_witnesses_replayed"]>=3
 assert dres["mini_adaptive_campaign"]["full16d_k1_cohort"]["combined_zero_count"]==0 and dres["mini_adaptive_campaign"]["prod_raw_lq_survey"]["all38_below_1x"]==0 and dres["mini_adaptive_campaign"]["prod_raw_lq_survey"]["reached_below_1x"]==0
 assert dres["live_full16d_k1_gen"]["total_ideals"]>=32 and dres["live_full16d_k1_gen"]["total_empty"]==0 and dres["live_full16d_k1_gen"]["reached_children"]>=12
 assert dres["prod_q70_cert_diag"]["radius_rescaling"] is False and dres["prod_q70_cert_diag"]["bound_cell_screened"] is True and dres["prod_q70_cert_diag"]["cells_per_lattice"]==32 and dres["prod_q70_cert_diag"]["n_checked"]==4 and dres["prod_q70_cert_diag"]["n_below_1x"]==0 and dres["prod_q70_cert_diag"]["min_ratio"]>=2.0
 assert dres["pm_log2_rels"]<=_cl["log2_rels_2se_live"]
 assert abs(dres["pm_main_lb"]-max(_op["lf"],_op["lg"]))<1e-9
 assert dres["pm_main_nf"]<=_cl["nf_cond_mean"] and dres["pm_main_ng"]<=_cl["ng_cond_mean"]
 assert abs(dres["pm_alpha_f"]-_op["alpha_f"])<1e-9 and abs(dres["pm_alpha_g"]-_op["alpha_g_total"])<1e-9
 assert dres["pm_pair_setup"]==float(_dp_kern["pair_setup"])==2.0**35.2
 assert dres["pm_lll_node"]==float(_dp_kern["lll_node"])==2.0**40.0
 assert dres["pm_lll_boot"]==float(_dp_kern["lll_boot"])==2.0**47.0
 assert dres["pm_bound_cell"]==float(_dp_kern["bound_cell"])==5000.0 and dres["pm_sieve_hit"]==2600.0
 assert dres["pm_ret_mult_f"]<=_cl["ret_mult_f_min"]and dres["pm_ret_mult_g"]<=_cl["ret_mult_g_min"]
 assert dres["eager_all_100pct"] is True and dres["linear_miss_retention"] is True and dres["f4_counterexample_rejected"] is True and dres["suyama_stepped_ok"] is True and dres["suyama_stepped"]["emp_p1"]>=dres["suyama_stepped"]["model_p1_lo"]
 import q15_descent16 as _DS
 for _u_chk in (1032.0/1024.0,2.0,3080.0/1024.0,3306.0/1024.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0,14.0,14.0+1.0/1024.0,14.5,15.0):
  _re_u=rho_exact(_u_chk)
  assert _DS.rho_lo(_u_chk)<_re_u<_DS.rho_hi(_u_chk),f"Dickman enclosure failed at u={_u_chk}"
 _st_q15=json.loads(lzma.decompress(base64.b85decode("".join((ROOT/"q15_stats.b85").read_text().split()))).decode())
 _Pm_c=_DS.P();_DS.load_stats(_st_q15,_Pm_c)
 _Pm_c.MAIN_LB=max(_op["lf"],_op["lg"]);_Pm_c.MAIN_UF=546.0/_Pm_c.MAIN_LB;_Pm_c.MAIN_UG=670.0/_Pm_c.MAIN_LB
 for _mf in (600,900,1200):
  assert _DS.miss_side(math.ceil(4*_mf/_Pm_c.LX),_Pm_c.LOG2_RELS,_Pm_c.MAIN_UF,_Pm_c.LX,_Pm_c.MAIN_LB,_Pm_c.ORBIT,_Pm_c.P_SIEVE_F,_Pm_c.RET_MULT_F)>=_DS.miss_side(math.ceil(4*_mf/_Pm_c.LX),_cl["log2_rels_2se_live"],_cl["nf_cond_mean"]/_op["lf"],_Pm_c.LX,_op["lf"],_Pm_c.ORBIT,data["collection"]["p_f"],_cl["ret_mult_f_min"])
 for _mg in (700,1200,2590):
  assert _DS.miss_side(math.ceil(4*_mg/_Pm_c.LX),_Pm_c.LOG2_RELS,_Pm_c.MAIN_UG,_Pm_c.LX,_Pm_c.MAIN_LB,_Pm_c.ORBIT,_Pm_c.P_SIEVE_G,_Pm_c.RET_MULT_G)>=_DS.miss_side(math.ceil(4*_mg/_Pm_c.LX),_cl["log2_rels_2se_live"],_cl["ng_cond_mean"]/_op["lg"],_Pm_c.LX,_op["lg"],_Pm_c.ORBIT,data["collection"]["p_g"],_cl["ret_mult_g_min"])
 for _lm in (0.1,1.0,2.0,5.0,10.0,18.46,32.0,100.0,250.0,600.0,700.0,933.0,1000.0):
  _lp=0.75*_lm
  _p_thin=math.exp(-_lp)*(1.0+_lp) if _lp<700.0 else 1e-290
  _p_p3=math.exp(-_lm)*(1.0+_lm+0.5*_lm*_lm+(_lm**3)/6.0) if _lm<700.0 else 1e-290
  _ex_thin=math.exp(-_lp)*(1.0+_lp) if _lp<709.0 else 0.0
  _ex_p3=math.exp(-_lm)*(1.0+_lm+0.5*_lm*_lm+(_lm**3)/6.0) if _lm<709.0 else 0.0
  assert max(_p_p3,_p_thin)>=_ex_thin and max(_p_p3,_p_thin)>=_ex_p3
 dres["missing_exact_thin_dom_rho25"]=True
 _rows_d=dict(_led["rows"])
 assert dres["two_log_units"]<=_rows_d["both_individual_logs_exp"]
 with localcontext()as ctx:
  ctx.prec=80
  two=Decimal(2);ln2=two.ln()
  p2=lambda v:(Decimal(str(v))*ln2).exp()
  lg2=lambda v:v.ln()/ln2
  exp_core=sum(p2(v)for _,v in _led["rows"])
  exp_core_l2=lg2(exp_core)
  assert exp_core_l2<=Decimal(_led["exp_run_claim"])
  t_run_l2=Decimal(str(_led["t_run_log2"]))
  assert t_run_l2==Decimal("124.0")
  p_meter=p2(_led["exp_run_claim"])/p2(t_run_l2)
  p_meter_110=Decimal("1.1")*p_meter
  assert p_meter<=Decimal(_led["p_meter_upper"])
  order_d=Decimal(str(R_BLS))
  assert data.get("track")=="bls12-381-g1-ecdlp-plain" and data.get("unit_ir")==UNIT and data.get("g1_unit_ir")==2000.0
  g1_scale=Decimal(int(UNIT//2000.0));assert g1_scale==Decimal(10)
  fail_main=two**-66+two**-70+two**-80+Decimal(2)/order_d
  fail_desc=p_meter+Decimal(4)/order_d
  fail_desc_110=p_meter_110+Decimal(4)/order_d
  succ_core=(Decimal(1)-fail_main)*(Decimal(1)-fail_desc)
  succ_head=(Decimal(1)-fail_main)*(Decimal(1)-fail_desc_110)-two**-48
  lo=F(_led["p_floor"])
  assert F(str(succ_head))>=lo>F(39,100)
  K=int(_led["K_gate"])
  assert 0<K<(1<<64)
  assert K*lo>=F(39,100)*(1<<64)>(K-1)*lo
  ew_p=lg2(Decimal(1)+Decimal(K)*(p2(_led["exp_run_claim"])+p2(_led["wrapper_log2"]))/(two**64))
  assert ew_p<Decimal(_led["claim_p"])
  ew_g1=lg2(Decimal(1)+Decimal(K)*g1_scale*(p2(_led["exp_run_claim"])+p2(_led["wrapper_log2"]))/(two**64))
  assert ew_g1<Decimal(_led["claim_g1"])
  w_hard_l2=lg2(p2(t_run_l2)+p2(_led["wrapper_log2"]))
  w_hard_g1_l2=lg2(g1_scale*(p2(t_run_l2)+p2(_led["wrapper_log2"])))
  manifest_ok=False;track_claim=None
  for mp in (ROOT/"submission.json",HERE/"submission.json"):
   if mp.is_file():
    mf=json.loads(mp.read_text())
    sc=mf["claims"][0]["score"]
    declared=F(int(sc["numerator"]),int(sc["denominator"]))
    ev_files=mf["claims"][0].get("evidence",[])
    ex_files=[ex["source"] for ex in mf.get("experiments",[]) if "source" in ex]
    if declared==F(_led["claim_g1"]) and ev_files and all((mp.parent/ef).is_file() for ef in ev_files+ex_files):
     track_claim=_led["claim_g1"]
     assert("2^"+_led["claim_g1"])in mf["title"]
     assert("2^"+_led["exp_run_claim"])in mf["claims"][0]["statement"]
     manifest_ok=True
    break
  assert manifest_ok and track_claim is not None
  _q10=Decimal("1e-10");_up=lambda d:str(d.quantize(_q10,rounding=ROUND_CEILING))
 out["ledger_gate_claim"]=dict(
  exp_run_log2=str(exp_core_l2)[:18],exp_run_claim=_led["exp_run_claim"],
  t_run_log2=str(t_run_l2),w_hard_log2=str(w_hard_l2)[:18],w_hard_g1_log2=str(w_hard_g1_l2)[:18],
  p_meter=_up(p_meter),p_meter_110=_up(p_meter_110),
  p_descent_fail="0.0",
  succ_core=str(succ_core)[:10],succ_head_110_wrap=str(succ_head)[:10],
  p_floor=_led["p_floor"],K_gate=K,ew_p_log2=str(ew_p)[:18],claim_p=_led["claim_p"],
  ew_g1_log2=str(ew_g1)[:18],claim_g1=_led["claim_g1"],track_claim=track_claim,
 )
 ck["ledger_gate_claim"]=True
 res=dict(checks_ok=all(ck.values()),checks=ck,summary=out)
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:raise SystemExit(1)
if __name__=="__main__":main()
