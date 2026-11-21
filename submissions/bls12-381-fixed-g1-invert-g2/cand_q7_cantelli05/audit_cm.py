"""audit_cm.py — Self-contained pure-Python-stdlib hosted Docker audit for BLS12-381 CM-STNFS over Q(zeta_15)..."""
import base64
import cmath
import lzma
import hashlib
import json
import math
import random
import sys
from decimal import Decimal,localcontext
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
if str(ROOT)not in sys.path:sys.path.insert(0,str(ROOT))
if str(HERE)not in sys.path:sys.path.insert(0,str(HERE))
from cm_models import(GAMMA,LN2,MONT5,NTT4_PRIMES,NTT4_PROD,Retention,primes_upto,ecm_instr,ecm_mean_var,eval_descent_live,fft_F,gbatch_cost,lll_16d_int,log2_li,mr_instr,rho_db,rho_exact,srt_toy_check,stage_b_analytic,verify_gt_boot_20,verify_kernels_symbolic,verify_prod_sq24,)
from filtsim_z12 import verify_h3_ladder_and_fit
from w7_verifiers import verify_genuine_q15_filter_hx_live,verify_q15_k1_strict_h5_live
from hx_z12 import run_hx_on_filtered_z12
from pm_adaptive import unplanted_two_sided_bw_check,verify_b16_dsieve_live,verify_prod_cond_q15,verify_q15_93_b16_live,verify_q15_enc_bias_live
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
def y_pow_z(k):
 k%=15;c=[0]*(k+1);c[k]=1;return red_phi15_z(c)
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
 # Exhaustive verification of ALL 56 prime ideals P of O_H=Z[zeta_15] with N(P)<256:
 # Proves unconditionally (1) min N(P)=16 > 9 = deg(f_9) > 3 = deg(g_3), hence ZERO fixed prime-ideal divisors P!=Q on any saturated M=L_q (|N(det)|==q);
 # (2) EXACT finite-field root count of f_9 and g_3 across all 56 O_H ideals <256 = 150 ideals (sum 1/N(Q)=1.7059), ALL sieved explicitly in dsieve16_wc_pc (ZERO skipped primes from 2 to 2^lx!).
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
  c=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[1]*b[1]+a[2]*b[0],a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0],a[1]*b[3]+a[2]*b[2]+a[3]*b[1],a[2]*b[3]+a[3]*b[2],a[3]*b[3]]
  for k in (6,5,4):
   v=c[k]%p0
   if v:o=k-4;c[o]-=v*_h0;c[o+1]-=v*_h1;c[o+2]-=v*_h2;c[o+3]-=v*_h3
  return[c[0]%p0,c[1]%p0,c[2]%p0,c[3]%p0]
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
 def sm_chi_proj(a8,b8,is_g):
  dX=3 if is_g else 9
  eps=(r0**3-1)if is_g else(r0**6-1)
  L_c=[[0]*8 for _ in range(dX)]
  for ci,s2 in enumerate(roots_l2):
   w_c=pow(s2,5,r0sq)
   w2_c=(w_c*w_c)%r0sq
   ac=sum(a8[j]*pow(s2,j,r0sq)for j in range(8))%r0sq
   bc=sum(b8[j]*pow(s2,j,r0sq)for j in range(8))%r0sq
   c0=(u0*w2_c)%r0sq if is_g else w_c%r0sq
   c3=(((w_c+1)*w_c)%r0sq)if not is_g else 0;c6=w2_c%r0sq if not is_g else 0
   _MSK=(1<<44)-1
   def _smul(x,y):
    ix=0
    for k in range(dX-1,-1,-1):ix=(ix<<44)|x[k]
    iy=0
    for k in range(dX-1,-1,-1):iy=(iy<<44)|y[k]
    p=ix*iy
    o=[(p>>(44*k))&_MSK for k in range(2*dX-1)]
    if is_g:
     for k in (4,3):v=o[k]%r0sq;o[k-3]+=v*c0
    else:
     for k in range(16,8,-1):
      v=o[k]%r0sq
      if v:o[k-9]+=v*c0;o[k-6]+=v*c3;o[k-3]+=v*c6
    return[o[k]%r0sq for k in range(dX)]
   pw=[1]+[0]*(dX-1);bs=[ac,(-bc)%r0sq]+[0]*(dX-2);ee=eps
   while ee:
    if ee&1:pw=_smul(pw,bs)
    ee>>=1
    if ee:bs=_smul(bs,bs)
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
    shift=y_pow_z((w_tau_exp*t*k)%15)
    term=mul_phi15_z(cur_p,shift)
    for j in range(8):
     acc[j]=(acc[j]+wt*term[j])%r0
    nxt=[0]*8
    for j,cj in enumerate(cur_p):
     if cj:
      y7j=y_pow_z(7*j)
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
 G=build_G1w()
 W,S=build_WS_1w()
 def mv(M,v):return[sum(M[i][j]*v[j]for j in range(8))for i in range(8)]
 def q_1w(v):return sum(v[i]*G[i][j]*v[j]for i in range(8)for j in range(8))
 I8=[[1 if i==j else 0 for j in range(8)]for i in range(8)]
 G12_6=[]
 cur_s=I8
 for _ in range(4):
  cur_w=cur_s
  for _ in range(3):
   assert all(cur_w[i][j]==0 for i in range(2,8)for j in range(2))
   G12_6.append([row[2:]for row in cur_w[2:]])
   cur_w=mat_mul(W,cur_w)
  cur_s=mat_mul(S,cur_s)
 assert len(G12_6)==12
 def in_D12(b):
  U=b[2]-b[4]-b[6]
  V=b[3]-b[5]-b[7]
  if U==0 and V==0:
   t=tuple(b[2:])
   if not any(t):return False
   return t==min(tuple(sum(M6[i][j]*t[j]for j in range(6))for i in range(6))for M6 in G12_6)
  if not(U>V and V>=0):return False
  if b[2]!=0 or b[3]!=0:return(b[2]>0)or(b[2]==0 and b[3]>0)
  d1x,d1y=b[4]-b[6],b[5]-b[7]
  if d1x!=0 or d1y!=0:return(d1x>0)or(d1x==0 and d1y>0)
  return False
 R=[[0.0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(i+1):
   s=sum(R[k][i]*R[k][j]for k in range(j))
   R[j][i]=math.sqrt(G[i][i]-s)if i==j else(G[i][j]-s)/R[j][j]
 pts=[]
 cur=[0]*8
 def fp_rec(k,rem):
  s=sum(R[k][j]*cur[j]for j in range(k+1,8))
  d=R[k][k];rad=math.sqrt(max(0.0,rem))/d;ctr=-s/d
  lo=math.ceil(ctr-rad-1e-11);hi=math.floor(ctr+rad+1e-11)
  if k==0:
   for v in range(lo,hi+1):
    cur[0]=v
    q=q_1w(cur)
    if 0<q<=30:pts.append(tuple(cur))
  else:
   for v in range(lo,hi+1):
    cur[k]=v;diff=d*v+s;nrem=rem-diff*diff
    if nrem>=-1e-9:fp_rec(k-1,nrem)
 fp_rec(7,30.0)
 seen=set();n_orb12=exact1=dup=missed=0
 for b in pts:
  if b in seen:continue
  orb=set();cs=list(b)
  for _ in range(4):
   cw=list(cs)
   for _ in range(3):
    orb.add(tuple(cw));cw=mv(W,cw)
   cs=mv(S,cs)
  seen|=orb
  if len(orb)==12:
   n_orb12+=1
   h=sum(1 for x in orb if in_D12(x))
   if h==1:exact1+=1
   elif h>1:dup+=1
   else:missed+=1
 assert dup==0 and missed==0 and len(pts)==2940 and n_orb12==239 and exact1==239
 def _pows_1w(s,p):
  s3=pow(s,3,p);s5=pow(s,5,p);s6=(s3*s3)%p;s9=(s6*s3)%p
  return(1,s5,s3,(s3*s5)%p,s6,(s6*s5)%p,s9,(s9*s5)%p)
 def eval_1w_pows(v,pw,p):
  return(v[0]*pw[0]+v[1]*pw[1]+v[2]*pw[2]+v[3]*pw[3]+v[4]*pw[4]+v[5]*pw[5]+v[6]*pw[6]+v[7]*pw[7])%p
 rng=random.Random(20261119)
 R2_a,R2_b=80.0,70.0;SCALES=(1.00,1.50,2.00)
 a_orbits=[];seen_a=set()
 while len(a_orbits)<40:
  v0=[rng.randint(-2,2)for _ in range(8)]
  if not any(v0):continue
  v1=mv(S,v0);v2=mv(S,v1);v3=mv(S,v2)
  orb=[tuple(x)for x in (v0,v1,v2,v3)]
  if len(set(orb))!=4 or any(x in seen_a for x in orb):continue
  qa=q_1w(v0)/R2_a
  a0_p=[sum(M_1W[i][j]*v0[j]for j in range(8))for i in range(8)]
  _na0=abs(norm_q15_z(a0_p));_v2=(_na0&-_na0).bit_length()-1 if _na0 else 8
  cls=0 if _v2==0 else(1 if _v2==4 else 2)
  for _vk in (v1,v2,v3):
   _akp=[sum(M_1W[i][j]*_vk[j]for j in range(8))for i in range(8)]
   _nak=abs(norm_q15_z(_akp));_v2k=(_nak&-_nak).bit_length()-1 if _nak else 8
   assert(0 if _v2k==0 else(1 if _v2k==4 else 2))==cls
  if qa<=SCALES[cls]:
   seen_a.update(orb)
   a_orbits.append(([v0,v1,v2,v3],cls))
 b_rows=[];seen_t=set()
 while len(b_rows)<35:
  tail=[0]+[rng.randint(-1,1)for _ in range(7)]
  if not in_D12(tail):continue
  tt=tuple(tail)
  if tt in seen_t:continue
  pts_r=[]
  for b0 in range(-6,7):
   bv=list(tail);bv[0]=b0
   assert in_D12(bv)
   if q_1w(bv)<=R2_b:pts_r.append((b0,bv))
  if len(pts_r)>=4:seen_t.add(tt);b_rows.append((tail,pts_r))
 def is_prime(n):return n>=2 and all(n%d!=0 for d in range(2,math.isqrt(n)+1))
 orbits4=[]
 for p in range(31,350):
  if p%15!=1 or not is_prime(p):continue
  s_roots=[s for s in range(2,p)if pow(s,15,p)==1 and pow(s,5,p)!=1 and pow(s,3,p)!=1]
  seen_s=set()
  for s0 in s_roots:
   if s0 in seen_s:continue
   sk=[s0,pow(s0,7,p),pow(s0,4,p),pow(s0,13,p)]
   seen_s.update(sk)
   w=pow(s0,5,p);w2=(w*w)%p
   rhs_g=(U_BLS*w2)%p
   for r in range(1,p):
    if pow(r,3,p)==rhs_g:orbits4.append(("g",p,sk,r));break
   c6=(-w2)%p;c3=((-w-1)*w)%p;c0=(-w)%p
   for r in range(1,p):
    r3=pow(r,3,p)
    if(pow(r3,3,p)+c6*r3*r3+c3*r3+c0)%p==0:orbits4.append(("f",p,sk,r));break
  if len(orbits4)>=12:break
 import struct
 W_SB=32;tot_hits=mis=0;MASK=(1<<64)-1;MASK59=(1<<59)-1;OID_M=(1<<57)-1
 for side,p,sk,r in orbits4:
  s0=sk[0];r_inv=pow(r,-1,p);n_sb=(p+W_SB-1)//W_SB
  pw0=_pows_1w(s0,p);hsb=[0]*n_sb;kbuf=[]
  for m in range(4):
   for vs,_ in a_orbits:
    ka=(eval_1w_pows(vs[m],pw0,p)*r_inv)%p
    hsb[ka>>5]+=1;kbuf.append(ka)
  off=[0]
  for c in hsb:off.append(off[-1]+c*8)
  fba=bytearray(off[-1]);pos=list(off[:-1]);idx_kb=0
  for m in range(4):
   x_oidm=m<<57
   for _ in range(len(a_orbits)):
    ka=kbuf[idx_kb];idx_kb+=1;sbi=ka>>5
    struct.pack_into("<Q",fba,pos[sbi],((ka&31)<<59)|x_oidm);pos[sbi]+=8;x_oidm+=1
  flat=bytes(fba)
  for k in range(4):
   s_k=sk[k];pwk=_pows_1w(s_k,p);xk57=k<<57
   a_evals=[(eval_1w_pows(vs[0],pwk,p),eval_1w_pows(vs[1],pwk,p),eval_1w_pows(vs[2],pwk,p),eval_1w_pows(vs[3],pwk,p))for(vs,cls) in a_orbits]
   for tail,pts_r in b_rows:
    L_b=len(pts_r);k0=eval_1w_pows(pts_r[0][1],pwk,p)
    rep_hits=set()
    ranges=[(k0,k0+L_b-1)]if k0+L_b<=p else[(k0,p-1),(0,(k0+L_b-1)%p)]
    for lo,hi in ranges:
     for sbi in range(lo>>5,(hi>>5)+1):
      for ofs in range(off[sbi],off[sbi+1],8):
       rec,=struct.unpack_from("<Q",flat,ofs)
       ka=(sbi<<5)|(rec>>59);t=(ka-k0)%p
       if t<L_b:
        rot=((rec-xk57)&MASK)&MASK59
        rep_hits.add((rot&OID_M,rot>>57,t))
    gt_hits=set()
    for t,(b0,bv) in enumerate(pts_r):
     rbs=(r*eval_1w_pows(bv,pwk,p))%p
     for oid,ev4 in enumerate(a_evals):
      if ev4[0]==rbs:gt_hits.add((oid,0,t))
      if ev4[1]==rbs:gt_hits.add((oid,1,t))
      if ev4[2]==rbs:gt_hits.add((oid,2,t))
      if ev4[3]==rbs:gt_hits.add((oid,3,t))
    if rep_hits!=gt_hits:mis+=len(rep_hits^gt_hits)
    tot_hits+=len(gt_hits)
 assert mis==0 and tot_hits>1000
 return dict(d12_pts=len(pts),d12_orb12=n_orb12,d12_exact1=exact1,d12_dups=dup,d12_missed=missed,tables=len(orbits4),ideals=4*len(orbits4),gt_hits=tot_hits,arch="Arch2-Prod-8B-flat-W32",mismatches=mis)
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
 import random
 A2,r_ab,s_vals=pool["A2_prod"],pool["r_ab"],pool["s_vals"]
 R2_a,R2_b=A2*r_ab,A2/r_ab;G=build_G1w()
 L=[[0.0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(i+1):
   s=G[i][j]-sum(L[i][k]*L[j][k]for k in range(j))
   L[i][j]=math.sqrt(s)if i==j else s/L[j][j]
 Um=[[L[j][i]for j in range(8)]for i in range(8)]
 def _su(rhs):
  x=[0.0]*8
  for i in range(7,-1,-1):x[i]=(rhs[i]-sum(Um[i][j]*x[j]for j in range(i+1,8)))/Um[i][i]
  return x
 places=[cmath.exp(2j*math.pi*k/15.0)for k in (1,7,4,2)]
 facets=[[abs(sum(u[j]*(z**j)for j in range(8)))**2 for z in places]for u in pool["voronoi_units_10"]]
 z3=y_pow_z(5);y1=y_pow_z(1);y2=y_pow_z(2);y3=y_pow_z(3)
 wp1_y2=mul_phi15_z(add_phi15_z(z3,[1]+[0]*7),y2);wy3=mul_phi15_z(z3,y3)
 def _pr2(a,mp):
  r=[x&1 for x in a]
  for i in range(7,3,-1):
   if r[i]:
    for j in range(5):r[i-4+j]^=mp[j]
  return r[:4]
 P1,P2=[1,1,0,0,1],[1,0,0,1,1]
 def _cls(ap):
  d1=not any(_pr2(ap,P1));d2=not any(_pr2(ap,P2))
  if not d1 and not d2:return 0
  if d1 and d2:return 2
  na=abs(norm_q15_z([x&3 for x in ap]))
  return 1 if(na and((na&-na).bit_length()-1)==4)else 2
 def _q1w(v):return sum(v[i]*math.sumprod(G[i],v)for i in range(8))
 def _tp(v):return[math.sumprod(M_1W[i],v)for i in range(8)]
 lf_op,lg_op,af_op,ag_op,k2_op=float(op["lf"]),float(op["lg"]),float(op["alpha_f"]),float(op["alpha_g_ge3"]),float(op["kap2_g"])
 smax=max(s_vals);rng=random.Random(pool["seed"]);nfs,ngs,ys,ysp,pr_fl=[],[],[],[],[];cc=[0,0,0];nvt=nvo=n_pr=0;ND=400
 while len(nfs)<ND:
  ua=[rng.gauss(0.0,1.0)for _ in range(8)];ra=math.sqrt(math.sumprod(ua,ua))
  da=_su([x/ra for x in ua]);dap=_tp(da)
  x4=[abs(sum(dap[j]*(z**j)for j in range(8)))**2 for z in places];t2=sum(x4);nvt+=1
  if not all(math.sumprod(f,x4)>=t2-1e-9 for f in facets):continue
  nvo+=1
  qat=smax*(rng.random()**0.25);sa=math.sqrt(qat*R2_a)
  a1w=[round(x*sa)for x in da]
  if not any(a1w):continue
  ap=_tp(a1w);c=_cls(ap);sc=s_vals[c]
  if _q1w(a1w)>sc*R2_a:continue
  ub=[rng.gauss(0.0,1.0)for _ in range(8)];rb=math.sqrt(math.sumprod(ub,ub));sb=math.sqrt(R2_b*(rng.random()**0.25))
  b1w=[round(x*sb)for x in _su([x/rb for x in ub])]
  if not any(b1w)or _q1w(b1w)>R2_b:continue
  bp=_tp(b1w);a3,b3=pow_phi15_z(ap,3),pow_phi15_z(bp,3);a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3)
  a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
  ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_BLS,mul_phi15_z(y1,b3)))))
  nf=abs(norm_q15_z(sub_phi15_z(sub_phi15_z(sub_phi15_z(a9,mul_phi15_z(y1,mul_phi15_z(a6,b3))),mul_phi15_z(wp1_y2,mul_phi15_z(a3,b6))),mul_phi15_z(wy3,b9))))
  if ng==0 or nf==0:continue
  v2g=(ng&-ng).bit_length()-1
  assert(c==0 and v2g==0)or(c==1 and v2g==12)or(c==2 and v2g>=24)
  l2f=math.log2(nf);l2g=math.log2(ng);ef_l=l2f+af_op;eg_l=l2g-v2g+ag_op+1.0+k2_op
  nfs.append(l2f);ngs.append(l2g);cc[c]+=1
  _yv=rho_db(ef_l/lf_op,ef_l)*rho_db(eg_l/lg_op,eg_l);ys.append(_yv)
  _g0=math.gcd(abs(norm_q15_z(ap)),abs(norm_q15_z(bp)));_cp=True
  if _g0>1:
   _t=_g0;_d=2;_ps=[]
   while _d*_d<=_t and _d<=10000:
    if _t%_d==0:
     _ps.append(_d)
     while _t%_d==0:_t//=_d
    _d+=1 if _d==2 else 2
   if _t>1:_ps.append(_t)
   for _pp in _ps:
    _g1=pgcd_mod(ap,PHI15,_pp)
    if len(_g1)>1 and len(pgcd_mod(_g1,bp,_pp))>1:_cp=False;break
  pr_fl.append(_cp)
  if _cp:n_pr+=1;ysp.append(_yv)
 nf_m=sum(nfs)/ND;ng_m=sum(ngs)/ND;fv=nvo/nvt;fv_ex=3.0/(5.0*math.sqrt(5.0))
 assert abs(fv-fv_ex)<=3.0*math.sqrt(fv_ex*(1.0-fv_ex)/nvt)
 nf_se=math.sqrt(sum((x-nf_m)**2 for x in nfs)/(ND*(ND-1)))
 ng_se=math.sqrt(sum((x-ng_m)**2 for x in ngs)/(ND*(ND-1)))
 z_nf=abs(nf_m-pool["nf_mean"])/nf_se;z_ng=abs(ng_m-pool["ng_mean"])/ng_se
 assert z_nf<=2.5 and z_ng<=2.5
 w_cls=(225.0/256.0,225.0/2048.0,23.0/2048.0);tws=sum(wc*(sc**4)for wc,sc in zip(w_cls,s_vals))
 for c in range(3):
  ec=w_cls[c]*(s_vals[c]**4)/tws;assert abs(cc[c]/ND-ec)/math.sqrt(ec*(1.0-ec)/ND)<=3.0
 vol8=lambda R2:(math.pi**4/24.0)*(R2**4)/1125.0
 va_l2=math.log2(tws*vol8(R2_a)*fv_ex/120.0)
 vb_l2=math.log2(vol8(R2_b)/3.0)
 v16_l2=math.log2(tws*vol8(R2_a)*vol8(R2_b)*fv_ex/360.0)
 y_p_all=[_yv if _f else 0.0 for _yv,_f in zip(ys,pr_fl)]
 ym=sum(y_p_all)/ND;yse=math.sqrt(sum((x-ym)**2 for x in y_p_all)/(ND*(ND-1)))
 l2_rels_unret=v16_l2+math.log2(ym);se_rels_b=yse/(ym*math.log(2.0))
 f_prim_l=n_pr/ND;z_prim=abs(f_prim_l-pool["f_prim"])/math.sqrt(pool["f_prim"]*(1.0-pool["f_prim"])/ND)
 assert z_prim<=2.5
 hh=0.5*float(pool.get("hist_step",2.0))
 hy_m=sum(wc*rho_db((bf+hh+af_op)/lf_op,bf+hh+af_op)*rho_db((bgo+hh+ag_op+1.0+k2_op)/lg_op,bgo+hh+ag_op+1.0+k2_op)for bf,bgo,_,wc,_ in pool["hist_2d"]if wc>0)
 l2_rels_hist_unret=pool["V_kept"]+math.log2(hy_m)
 z_rels=abs(l2_rels_unret-l2_rels_hist_unret)/se_rels_b
 assert abs(va_l2-pool["log2_Va"])<1e-4 and abs(vb_l2-pool["log2_Vb_all"])<1e-4 and abs(v16_l2-pool["V_kept"])<1e-4
 assert z_rels<=2.5
 return dict(n_draws=ND,v2g_exact=True,nf_live=round(nf_m,3),ng_live=round(ng_m,3),z_nf=round(z_nf,2),z_ng=round(z_ng,2),cls_shares=[round(x/ND,4)for x in cc],f_vor_exact=round(fv_ex,6),f_vor_mc=round(fv,5),log2_Va_live=round(va_l2,4),log2_Vb_live=round(vb_l2,4),V_kept_live=round(v16_l2,4),log2_rels_unret_live=round(l2_rels_unret,4),log2_rels_unret_hist=round(l2_rels_hist_unret,4),z_rels=round(z_rels,2),prim_live_400=f"{n_pr}/{ND}",z_prim=round(z_prim,2),prim_yield_ratio=round((sum(ysp)/len(ysp))/(sum(ys)/ND),4))
def check_voronoi_and_boot_witnesses(data,boot_w):
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
 w_pi=((U_BLS**3-U_BLS**2-U_BLS)*pow(U_BLS+1,-1,P_BLS))%P_BLS
 h_pi=pgcd_mod(PHI15,[-w_pi%P_BLS,0,0,0,0,1],P_BLS)
 bt_res=verify_gt_boot_20(boot_w,data,P_BLS,R_BLS,U_BLS,h_pi,w_pi,pmul_mod,prem_mod,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,)
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
 return dict(voronoi_units_ok=10,anchors_min_t2=f"{anc_ok}/12",**bt_res,constructive_lll16_k=k_ind,constructive_lll16_index=det_l//(p0_l**4),)
def check_pool_coll_la_descent_ledger(data,descent_mesh=None,dlive_out=None):
 op=data["operating_point"]
 pool=data["pool"]
 coll=data["collection"]
 la=data["la"]
 desc=data["descent"]
 led=data["ledger"]
 lf,lg=op["lf"],op["lg"]
 p_g,p_f,p_lo=coll["p_g"],coll["p_f"],coll["p_lo"]
 det=coll["det"]
 log2_Up=math.log2((2.0**log2_li(lf)+2.0**log2_li(lg))/12.0)
 assert abs(log2_Up-op["log2_Up"])<1e-9
 ret_g=Retention(round(lg,2),int(round(p_g)),prime_min=2**24,explicit_bits=14,step=0.5,max_bits=800,tilt=0.05,odd=1)
 cst_g=Retention(800,int(round(p_g)),prime_min=2**24,explicit_bits=14,step=0.5,max_bits=801,tilt=0.005,odd=1)
 ret_f=Retention(round(lf,2),int(round(p_f)),prime_min=2**24,explicit_bits=14,step=0.5,max_bits=660,tilt=0.05,odd=1)
 cst_f=Retention(660,int(round(p_f)),prime_min=2**24,explicit_bits=14,step=0.5,max_bits=661,tilt=0.005,odd=1)
 sag,sbg,saf,sbf=det["sag"],det["sbg"],det["saf"],det["sbf"]
 kmg,kmf=int(sbg//p_g),int(sbf//p_f)
 c_fbatch=gbatch_cost(1.000028*(2.0**lf)/LN2,leaf=9)
 c_batch2=gbatch_cost(1.000028*(2.0**p_g)/LN2,leaf=11)
 c_ecm_g=ecm_instr(lg,sbg+8.0,0.999)
 c_mr_g=mr_instr(sbg+8.0)
 alpha_f,alpha_g_ge3,kap2_g=op["alpha_f"],op["alpha_g_ge3"],op["kap2_g"]
 rg_c,sga_c,sgb_c={},{},{}
 rf_c,sfa_c,sfb_c,qf_c={},{},{},{}
 ysum_mid=ysum_rethi=ysum_strict_hi=c_cof_num=w_k_tot=0.0
 max_bg_seen=0
 hs=float(pool.get("hist_step",2.0))
 hh=0.5*hs
 for bf,bgo,bgm,wc,wk in pool["hist_2d"]:
  if bgm>max_bg_seen:max_bg_seen=bgm
  if bf not in rf_c:
   rf_c[bf]=(ret_f.at(bf+hh,saf,sbf),ret_f.at(bf+hs,saf,sbf))
   sfa_c[bf]=min(1.0,1.25*cst_f.survivor_upper(float(bf),saf,0.25))
   sfb_c[bf]=min(1.0,1.25*cst_f.survivor_upper(float(bf),saf,0.25,sbf,lf,kmf))
   qf_c[bf]=min(1.0,2.0*rho_exact((bf+alpha_f)/lf))
  if bgo not in rg_c:
   rg_c[bgo]=(ret_g.at(bgo+hh,sag,sbg),ret_g.at(bgo+hs,sag,sbg))
   sga_c[bgo]=cst_g.survivor_upper(float(bgo),sag,0.25)
   sgb_c[bgo]=cst_g.survivor_upper(float(bgo),sag,0.25,sbg,lg,kmg)
  if wc>0.0:
   ef_c=bf+hh+alpha_f;eg_c=bgo+hh+alpha_g_ge3+1.0+kap2_g
   ef_h=bf+hs+alpha_f;eg_h=bgo+hs+alpha_g_ge3+1.0+kap2_g
   p_sm=rho_db(ef_c/lf,ef_c)*rho_db(eg_c/lg,eg_c)
   p_sm_hi=rho_db(ef_h/lf,ef_h)*rho_db(eg_h/lg,eg_h)
   ysum_mid+=wc*p_sm*rf_c[bf][0]*rg_c[bgo][0]
   ysum_rethi+=wc*p_sm*rf_c[bf][1]*rg_c[bgo][1]
   ysum_strict_hi+=wc*p_sm_hi*rf_c[bf][1]*rg_c[bgo][1]
  if wk>0.0:
   sga,sgb=sga_c[bgo],sgb_c[bgo]
   sfa,sfb,qf=sfa_c[bf],sfb_c[bf],qf_c[bf]
   csb_g=stage_b_analytic(float(bgm))
   csb_f=stage_b_analytic(float(bf+hs))
   c_ng=sga*6500.0
   c_nf=sga*20000.0
   c_sb=(sga*sfa)*csb_g+(sgb*sfa)*csb_f
   q_f_cond=min(1.0,qf/max(sfb,1e-12))
   c_fb=(sgb*sfb)*(c_fbatch+q_f_cond*(2.0**24.5+2.0**24.0+c_batch2+c_mr_g+c_ecm_g))
   c_cof_num+=wk*(c_ng+c_nf+c_sb+c_fb)
   w_k_tot+=wk
 c_cof_live=c_cof_num/w_k_tot
 assert abs(c_cof_live-det["c_cof_2d_upper"])<1e-6
 assert max_bg_seen>=690,"2D histogram must cover high-norm tail!"
 f_prim_hist=sum(wc for _,_,_,wc,_ in pool["hist_2d"])
 assert abs(f_prim_hist-pool["f_prim"])<1e-6 and 0.94<f_prim_hist<0.97
 log2_rels_h=(pool["log2_base"]+pool["dV"])+math.log2(ysum_mid)
 log2_rels_rethi=(pool["log2_base"]+pool["dV"])+math.log2(ysum_rethi)
 log2_rels_strict_hi=(pool["log2_base"]+pool["dV"])+math.log2(ysum_strict_hi)
 margin_net=(log2_rels_rethi-0.009)-(log2_Up+math.log2(op["exc"]))
 margin_strict_hi=log2_rels_strict_hi-(log2_Up+math.log2(op["exc"]))
 assert margin_net>=0.350 and margin_strict_hi>=0.250 and margin_strict_hi>=40.0*pool["se_y_bits"],"Generous finite-region O_H-coprime margin (>=+0.24b strict-hi, >=+0.37b net)!"
 assert log2_rels_rethi-0.009-2.0*pool["se_y_bits"]>=desc["params"]["log2_rels"]
 assert coll["h2_allow"]>=1.0+2.0*max(pool["se_c_rel"],pool["se_w_rel"])
 V_lin=2.0**pool["V_kept"]
 p2_mult=pool["p2_mult"]
 N_a_tbl=4.0*(2.0**pool["log2_Va"])
 N_7=2.0**pool["log2_N7_d12"]
 Rb_unit2=pool["A2_prod"]/pool["r_ab"]
 L_max_live=2.0*math.sqrt(Rb_unit2/8.0)
 Lb_live=L_max_live*(105.0*math.pi/768.0)
 assert abs(L_max_live-pool["L_max"])<1e-4 and abs(Lb_live-pool["Lb_avg"])<1e-4 and L_max_live<coll["w_tile"]
 assert abs(N_a_tbl*N_7*Lb_live-V_lin)/V_lin<1e-6,"Product volume N_a_tbl * N_7 * L_b must equal V_kept!"
 w_tile,w_sb,p_m4=coll["w_tile"],coll["w_sb"],coll["p_m4"]
 def price_side(p_hi):
  pl=min(p_lo,p_hi);pm=min(p_m4,pl)
  n1a=1.01*(2.0**log2_li(pm)-2.0**log2_li(24.0))
  n1b=1.01*(2.0**log2_li(pl)-2.0**log2_li(pm))if pl>pm else 0.0
  b1_lk=n1a*44.0+n1b*34.0
  ib1=(math.log(pm/24.0)*coll["c_hit_lo_m4"]+(math.log(pl/pm)*coll["c_hit_lo_u"]if pl>pm else 0.0))*p2_mult
  ib1+=(b1_lk*N_7*Lb_live+(n1a+n1b)*N_7*24.0)/V_lin
  sp1a=4.0*1.01*(2.0**(2.0*pm))/(2.0*pm*LN2-1.0)/4.0
  sp1b=(1.01*(2.0**(2.0*pl)-2.0**(2.0*pm))/(2.0*pl*LN2-1.0)/4.0)if pl>pm else 0.0
  bld=((n1a+n1b)/4.0)*coll["c_list_a"]*N_a_tbl+coll["c_lbkt"]*(sp1a+sp1b)
  mem=((n1a+n1b)/4.0)*8.0*N_a_tbl+8.0*(sp1a+sp1b)
  ipr=ihi=0.0
  if p_hi>pl:
   n2=1.01*(2.0**log2_li(p_hi)-2.0**log2_li(pl))
   sb_r=Lb_live/w_sb+1.5
   ipr=n2*N_7*(coll["c_tile_hdr"]+sb_r*coll["c_sb_hdr"])/V_lin
   cand_bd=1.0+w_sb/Lb_live
   ihi=math.log(p_hi/pl)*(cand_bd*coll["c_bd_chk"]+coll["c_hit_hi"])*p2_mult
   sp2=1.01*(2.0**(2.0*p_hi)-2.0**(2.0*pl))/(2.0*p_hi*LN2-1.0)/4.0
   bld+=(n2/4.0)*coll["c_list_a"]*N_a_tbl+coll["c_lbkt"]*(sp2/w_sb)
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
 col_own_mem=4.0*(1<<57)*4.0
 m1a_bld=16.0*N_a_tbl+64.0*(2.0**p_m4)+w_tile*(1<<20)
 m1b_bld=16.0*N_a_tbl+16.0*(2.0**p_lo)+w_tile*(1<<20)
 m2_bld=16.0*N_a_tbl+16.0*((2.0**p_g)/w_sb)+w_tile*(1<<20)
 assert 2.0**pool["log2_Va"]<(1<<57)*0.90 and col_own_mem<2.0**64,"Column owner (exact 2^57 bank stride) must fit < 2^64 B!"
 assert max(m1a_bld,m1b_bld,m2_bld)<2.0**64,"1-owner simultaneous build (kbuf+out+off+pos+staging) must fit < 2^64 B!"
 assert mg+mf<2.0**120,"Total replay tables must fit < 2^120 B!"
 mem_tbl_tot=mg+mf
 np_band=1.25506*(2.0**48)/(48.0*LN2)
 c_rel_band=np_band*(6*(10+12)+44+4)+(13+16)*1024+(1<<26)
 n_rels_upper=2.0**(pool["log2_rels"]+0.55)
 band_u_live=math.log2(n_rels_upper*c_rel_band/UNIT)
 assert band_u_live<=rows_d["band_split"]
 ef_e,ef_v=ecm_mean_var(lf,640,0.999)
 eg_e,eg_v=ecm_mean_var(lg,768,0.999)
 c_ecm_rel=(13.0*ef_e+14.0*mr_instr(640)+13.0*(1<<26)+16.0*eg_e+17.0*mr_instr(768)+16.0*(1<<26))
 ew_ecm_live=(c_ecm_rel/UNIT)*n_rels_upper
 sd_rel=(13.0*math.sqrt(ef_v)+16.0*math.sqrt(eg_v))/UNIT
 ecm_var_live=n_rels_upper*(sd_rel*sd_rel)
 ecm_gap=(2.0**rows_d["ecm_split_expected"])-ew_ecm_live
 assert ecm_gap>0
 ecm_cheb_p=ecm_var_live/(ecm_gap*ecm_gap)
 assert ecm_cheb_p<2.0**-70
 _h3=verify_h3_ladder_and_fit(la,log2_Up,run_u12=False);k5_live=_h3["k5_live"]
 assert led["cap_abort"]==0.0
 assert desc["mode"]=="pure_k1_boole_wald_iid_boot" and desc["params"]["p0_mode"]=="cantelli05" and desc["params"]["boot_rels"]==1
 assert la["tab_ratio_charged"]==1.0 and la["n_owners"]==1 and la["te"]==0,"Plain-sparse ZERO-TABM 1-owner 4-vec Arch B!"
 N=2.0**(log2_Up+math.log2(la["nl_cap"]))
 N_sw=2.0**(log2_Up+math.log2(la["sw_cap"]))
 m,n,n_sm=la["m"],la["n"],la["n_sm"]
 la_mem_full=2.0*N_sw*128.0;horner_1own_mem=2.0*N_sw*32.0
 assert la_mem_full<2.0**64 and horner_1own_mem<2.0**64,"Single 64-bit LA owner must fit < 2^64 B!"
 c_prod=(la["sw_cap"]/la["nl_cap"])*(201.0*la["k5_charged"]*la["tab_ratio_charged"]+la["c_row_chi"])
 bw_sig=2.0**-7
 kap=((1.0+bw_sig)*(N+n_sm)+6.0/(1.0/m+1.0/n))/N
 per_n2=(c_prod+60.0)*(n/m+1.0)*kap+(1.0/m+1.0/n)*(c_prod+n*97.0+270.0)
 la_u=math.log2((N*N*per_n2+N*(m+n)*270.0)/UNIT)
 L_bm=(1.0+bw_sig)*(N+n_sm)*(1.0/m+1.0/n)+6.0
 t_bm=m+n;depth_bm=math.ceil(math.log2(L_bm));k_bm=math.ceil(math.log2(4.0*L_bm))
 assert(1<<k_bm)<=(1<<53)and t_bm*L_bm*(1<<128)<NTT4_PROD
 w_base_bm=L_bm*m*(t_bm**2)*MONT5
 w_rec_bm=2.0*depth_bm*(8.0*L_bm)*(1536.0*(t_bm**3)+(t_bm**2)*(1536.0*k_bm+6220.0))
 gen_u=math.log2((w_base_bm+w_rec_bm)/UNIT)
 assert la_u<=rows_d["linear_algebra"]
 assert gen_u<=rows_d["matrix_generator"]
 assert len(desc["plan_struct"])==2*desc["n_bins"]==78
 m_node=384 if descent_mesh==384 else(96 if descent_mesh==96 else(64 if descent_mesh==64 else 32))
 m_boot=370 if descent_mesh==384 else 185
 dlive=eval_descent_live(desc,mesh_node=m_node,mesh_boot=m_boot)
 if dlive_out is not None:dlive_out["dlive"]=dlive
 _udup_live,_yr_min,_fp_h_n,_ds_nph=verify_b16_dsieve_live(data,desc,dlive,m_node,U_BLS,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int,)
 assert dlive["peel_h"]<0.50 and dlive["exc_crit"]<0.95
 two_log_exp_live=dlive["two_log_exp_units"]
 pr_b_log2_live=dlive["boot_pr_log2"]
 p_desc_loss_live=dlive["p_desc_total_loss"]
 assert two_log_exp_live<desc["two_log_upper"]
 assert float(desc["params"]["p_unrec_leaf"])==0.0 and p_desc_loss_live==0.0 and float(desc["p_desc_loss_upper"])==0.0
 assert dlive["min_mu"]>=30.0 and dlive["max_p0"]<0.01 and 0.0<dlive["Phi_tree"]<=0.60 and dlive["N_tree_1log_det"]<=4096
 assert dlive["p70_tree_build_u"]<=rows_d["pg_pf_product_trees"]-1.0
 assert desc["boot"]["cells_b"]==1 and desc["params"]["ball_pair"]==0 and desc["audit_11"]["ball_pair_charged"]==0
 assert desc["audit_11"]["sv_max"]<desc["boot"]["m_sv"]
 with localcontext()as ctx:
  ctx.prec=80
  two=Decimal(2);ln2=two.ln()
  p2=lambda v:(Decimal(str(v))*ln2).exp()
  lg2=lambda v:v.ln()/ln2
  exp_core=sum(p2(v)for _,v in led["rows"])
  exp_core_l2=lg2(exp_core)
  assert exp_core_l2<=Decimal(led["exp_run_claim"])
  t_run_l2=Decimal(str(led["t_run_log2"]))
  assert t_run_l2==Decimal("123.0")
  p_meter=p2(led["exp_run_claim"])/p2(t_run_l2)
  assert p_meter<=Decimal(led["p_meter_upper"])
  order_d=Decimal(str(R_BLS))
  p_rej_1=(Decimal(1)-order_d/(two**255))**64
  n_rej_ent=Decimal(20+(n-4)*(math.ceil(N)+4))
  p_rej_tot=n_rej_ent*p_rej_1
  assert p_rej_tot<two**-150
  fail_core=(two**-66+Decimal(str(led["cap_abort"]))+two**-70+two**-80+Decimal(str(desc["p_desc_loss_upper"]))+p_meter+Decimal(2)/order_d+Decimal(4)/order_d)
  fail_hr=fail_core+two**-48
  succ_core=Decimal(1)-fail_core
  succ_hr=Decimal(1)-fail_hr
  lo=F(led["p_floor"])
  assert F(str(succ_hr))>=lo>F(39,100)
  K=int(led["K_gate"])
  assert 0<K<(1<<64)
  assert K*lo>=F(39,100)*(1<<64)>(K-1)*lo
  ew_p=lg2(Decimal(1)+Decimal(K)*(p2(led["exp_run_claim"])+p2(led["wrapper_log2"]))/(two**64))
  ew_g1=lg2(Decimal(1)+Decimal(K)*Decimal(10)*(p2(led["exp_run_claim"])+p2(led["wrapper_log2"]))/(two**64))
  assert ew_p<Decimal(led["claim_p"])
  assert ew_g1<Decimal(led["claim_g1"])
  manifest_ok=True;track_claim=None
  for mp in (ROOT/"submission.json",HERE/"submission.json"):
   if mp.is_file():
    mf=json.loads(mp.read_text())
    sc=mf["claims"][0]["score"]
    declared=F(int(sc["numerator"]),int(sc["denominator"]))
    if declared==F(led["claim_p"]):
     track_claim=led["claim_p"]
     assert("2^"+led["claim_p"]) in mf["title"]
     assert("2^"+led["exp_run_claim"]) in mf["claims"][0]["statement"]
    elif declared==F(led["claim_g1"]):
     track_claim=led["claim_g1"]
     assert("2^"+led["claim_g1"]) in mf["title"]
     assert("2^"+led["exp_run_claim"]) in mf["claims"][0]["statement"]
    else:manifest_ok=False
    break
  assert manifest_ok
 return dict(c_cof_2d_upper_live=round(c_cof_live,4),max_bg_bin_covered=max_bg_seen,log2_rels_2d_mid=round(log2_rels_h,6),log2_rels_2d_rethi=round(log2_rels_rethi,6),log2_rels_2d_strict_hi=round(log2_rels_strict_hi,6),margin_net_rethi_bits=round(margin_net,5),margin_strict_hi_bits=round(margin_strict_hi,5),L_max_live=round(L_max_live,2),ipc_tot=round(ipc_tot,4),coll_u=round(coll_u,6),build_u=round(build_u,4),max_1own_bld_log2_B=round(math.log2(max(m1a_bld,m1b_bld,m2_bld)),3),mem_tables_log2_B=round(math.log2(mem_tbl_tot),3),ecm_cheb_log2=round(math.log2(ecm_cheb_p),2),k5_live=round(k5_live,4),la_u=round(la_u,6),gen_u=round(gen_u,4),la_mem_full_log2=round(math.log2(la_mem_full),3),peel_h=round(dlive["peel_h"],6),exc_crit=round(dlive["exc_crit"],4),desc_nodes_live=len(desc["plan_struct"]),desc_mesh=f"{m_node}/{m_boot}",boot_pr_log2_live=round(pr_b_log2_live,6),cells_b16_exact_live=True,b16_udup_live=_udup_live,b16_yr_min_live=_yr_min,dsieve16_fp_bf_hits=_fp_h_n,dsieve16_nodes_per_hit_R64=_ds_nph,f_prim_mc_200k=round(f_prim_hist,6),desc_two_log_exp_live=round(two_log_exp_live,6),desc_p0_mode=dlive["p0_mode"],desc_c_mu=dlive["c_mu"],desc_N_tree_1log_det=int(dlive["N_tree_1log_det"]),desc_N_tree_2log_det=int(dlive["N_tree_2log_det"]),desc_max_p0=float(f"{dlive['max_p0']:.6f}"),desc_Phi_tree=round(dlive["Phi_tree"],6),desc_wald_mult=round(dlive["wald_mult"],4),desc_p_total_loss=0.0,desc_p71_tree_build_u=round(dlive["p71_tree_build_u"],3),desc_min_mu=round(dlive["min_mu"],2),unified_meter_log2=str(t_run_l2),p_meter=str(p_meter)[:12],h3_fixed_instance_cap_abort=led["cap_abort"],exp_run_log2=str(exp_core_l2)[:18],succ_with_wrap=str(succ_hr)[:10],p_floor=led["p_floor"],K_gate=K,ew_p_log2=str(ew_p)[:18],ew_g1_log2=str(ew_g1)[:18],track_claim=track_claim,)
def check_hx_filtered_live():
 res=run_hx_on_filtered_z12(log2U=12,seed=5,m=12,n=8,s_sm=4)
 assert(res["checks_ok"]and res["l_rec_strict_lower_tri"]and res["prod_eq_schur"]and res["rank_M"]==res["N"]and res["rank_A"]==res["D"]-1 and res["kernel_eq_zero"]and res["recovered_true_kernel"]and res["post_peel_cols_lost"]==0 and res["post_merge_cols_lost"]==0 and res["repurge_backsub_errors"]==0 and res["n_repurge_backsub"]==res["repurge_stack_len"]and res["n_reconstruct_unit"]==0 and res["all_2core_rels_zero"]and res["recovered_cols_match_true"]and res["cols_unrecovered_tri"]==0 and res["backsub_2core_exact"])
 return res
def main():
 raw=sys.stdin.read().strip()
 req=json.loads(raw)if raw else{}
 allowed_keys={"descent_mesh"}
 if not set(req).issubset(allowed_keys):raise ValueError("audit_cm.py accepts {} or {'descent_mesh': 384}")
 descent_mesh=req.get("descent_mesh")
 if descent_mesh is not None and descent_mesh not in (32,64,96,384):raise ValueError("descent_mesh must be 32, 64, 96, or 384")
 data=json.loads((ROOT/"cm_data.json").read_text())
 _tb=json.loads(lzma.decompress(base64.b85decode((ROOT/"cm_tables.b85").read_bytes())).decode())
 for _k in ("pool","collection","la","descent"):data[_k]=_tb["cfg_"+_k]
 data["pool"]["hist_2d"]=_tb["pool_hist_2d"];data["pool"]["anchors"]=_tb["pool_anchors"];data["pool"]["voronoi_units_10"]=_tb["pool_vor"]
 data["prod_sq24"]=_tb["prod_sq24"];data["prod_enc12"]=_tb["prod_enc12"];data["q15_93_b16"]=_tb["q15_93_b16"];data["q15_cond7e"]=_tb["q15_cond7e"]
 data["descent"]["plan_struct"]=_tb["descent_plan_struct"];data["descent"]["qbins"]=_tb["descent_qbins"]
 data["descent"]["qbins_lo"]=_tb["descent_qbins_lo"];data["descent"]["qbins_sub"]=_tb["descent_qbins_sub"];data["descent"]["cells_rank2"]=_tb["descent_cells_rank2"];data["la"]["ladder"]=_tb["la_ladder"];data["descent"]["audit_11"]=_tb["desc_audit11"]
 toy=json.loads(lzma.decompress(base64.b85decode((ROOT/"toy_e2e.b85").read_bytes())).decode())
 bw_p=ROOT/"boot_witness.json" if(ROOT/"boot_witness.json").is_file()else(ROOT/"boot_witness_20.json")
 boot_w=json.loads(bw_p.read_text())
 out={}
 ck={}
 import unitdomain as _UD
 out["algebra_bls12"]=check_algebra_bls12();ck["algebra_bls12"]=True
 out["murphy_alpha_exact_live"]=_UD.verify_murphy_alpha_live(data["operating_point"],U_BLS,primes_upto,pgcd_mod,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,rho_db);ck["murphy_alpha_exact_live"]=True
 out["real_fb_joint_coll_live"]=_UD.verify_real_fb_joint_coll(data["pool"]["voronoi_units_10"],build_G1w,build_WS_1w,mat_mul,M_1W,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,primes_upto,rho_db,Retention);ck["real_fb_joint_coll_live"]=True
 out["toy_e2e_dlp"]=check_toy_e2e_dlp(toy);ck["toy_e2e_dlp"]=True
 out["sigmaB_exactness"]=check_sigmaB_exactness();ck["sigmaB_exactness"]=True
 out["anchor_norms"]=check_anchor_norms(data["pool"]["anchors"]);out["pool_2adic_live"]=verify_pool_2adic_live(data["pool"],data["operating_point"]);ck["anchor_norms"]=True
 out["voronoi_and_boot_live"]=check_voronoi_and_boot_witnesses(data,boot_w);ck["voronoi_and_boot_live"]=True
 out["srt_toy"]=srt_toy_check();ck["srt_toy"]=True
 dlive_box={}
 out["ledger_gate_claim"]=check_pool_coll_la_descent_ledger(data,descent_mesh=descent_mesh,dlive_out=dlive_box);ck["ledger_gate_claim"]=True
 out["filtsim_and_hx_live"]=check_hx_filtered_live();ck["filtsim_and_hx_live"]=True
 _h3r=verify_h3_ladder_and_fit(data["la"],data["operating_point"]["log2_Up"],run_u12=True);_u12p=_h3r.get("live_u12_seed0");out["h3_ladder_and_genuine_q15_hx"]=dict(h3=_h3r,genuine_q15_hx=verify_genuine_q15_filter_hx_live(toy,primes_upto,PHI15,y_pow_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,pmul_mod,prem_mod,u12_proxy=_u12p));ck["h3_ladder_and_genuine_q15_hx"]=True
 out["unplanted_two_sided_bw"]=unplanted_two_sided_bw_check();ck["unplanted_two_sided_bw"]=True
 out["kernels_symbolic"]=verify_kernels_symbolic((ROOT/"KERNELS.md").read_text(),data["collection"],data["la"],R_BLS,float(data["descent"]["params"]["lll_node"]),float(data["descent"]["params"]["lll_boot"]),data["descent"]["params"]);ck["kernels_symbolic"]=True
 out["prod_cond_q15_live"]=verify_prod_cond_q15(data["prod_enc12"],data["descent"],dlive_box["dlive"],U_BLS,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,lll_16d_int,vor_units=data["pool"]["voronoi_units_10"],pgcd_mod=pgcd_mod);ck["prod_cond_q15_live"]=True
 out["prod_sq24_live"]=verify_prod_sq24(data["prod_sq24"],data["descent"],dlive_box["dlive"],U_BLS,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,vor_units=data["pool"]["voronoi_units_10"]);ck["prod_sq24_live"]=True
 out["q15_93_b16_live"]=verify_q15_93_b16_live(data["q15_93_b16"],primes_upto,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int);ck["q15_93_b16_live"]=True
 out["q15_enc_bias_live"]=verify_q15_enc_bias_live(data["q15_cond7e"],data["descent"],U_BLS,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,lll_16d_int);ck["q15_enc_bias_live"]=True
 out["q15_k1_strict_h5_live"]=verify_q15_k1_strict_h5_live(None,primes_upto,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,lll_16d_int,desc_params=data["descent"]["params"]);ck["q15_k1_strict_h5_live"]=True
 res=dict(checks_ok=all(ck.values()),checks=ck,summary=out)
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:raise SystemExit(1)
if __name__=="__main__":main()
