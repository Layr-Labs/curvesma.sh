"""v20_lat_b2.py — BN254 C6 v26 pure-stdlib live verifier for Lattice Band 2 (LSE STADDH, 5-op LUT4 radius lookup, finite-difference k=1, 4-plane cone, direct live score nodes, 32x8=256 Haldane-HT) and the P_min=24 SRT."""
import math,random,struct
FP_C0=21.0;FP_C1=82.0;FP_C2=99.0
FP_CL=[21,77,94,79,77,75,73,71,69,67,65,94]
PHI_F3=((1,0,2,2,0,1),(1,2,0,2,1,2),(1,2,2,1,2,2))
WHEEL_OPS=dict(chord_enter=0,emit9_step1=9,emit9_step2=9,emit9_step3=9,emit9_last=7,wheel_m0_check=2,wheel_m0_disp=4,wheel_chord_rot=3,wheel_seg_disp0=5,wheel_skip1_stub=2,wheel_skip2_stub=2)
AMAX_MR=(416,387,362,276,276,276)
BMAX_MR=(553,514,481,367,367,367)
MA_MR=tuple(2*m+1 for m in AMAX_MR)
MB_MR=tuple(2*m+1 for m in BMAX_MR)
WA_MR=[1]
for _m in MA_MR[:-1]:WA_MR.append(WA_MR[-1]*_m)
WB_MR=[1]
for _m in MB_MR[:-1]:WB_MR.append(WB_MR[-1]*_m)
PROD_MA_MR=WA_MR[-1]*MA_MR[-1]
PROD_MB_MR=WB_MR[-1]*MB_MR[-1]
A_OFF_MR=sum(AMAX_MR[i]*WA_MR[i]for i in range(6))
B_OFF_MR=sum(BMAX_MR[i]*WB_MR[i]for i in range(6))
def _mr_da(a):return sum(a[i]*WA_MR[i]for i in range(6))&((1<<64)-1)
def _mr_db5(b):return(sum(b[i]*WB_MR[i]for i in range(6))<<5)&((1<<64)-1)
def _mr_pack_a(a):return(A_OFF_MR+sum(a[i]*WA_MR[i]for i in range(6)))&((1<<64)-1)
def _mr_pack_b(b):return B_OFF_MR+sum(b[i]*WB_MR[i]for i in range(6))
HT_SEED_BASE=1
def _st_f3(v):
 a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5=v
 t0=(a0+2*a2+2*a3+a5)%3;t1=(a0+2*a1+2*a3+a4+2*a5)%3;t2=(a0+2*a1+2*a2+a3+2*a4+2*a5)%3
 t3=(b0+2*b2+2*b3+b5)%3;t4=(b0+2*b1+2*b3+b4+2*b5)%3;t5=(b0+2*b1+2*b2+b3+2*b4+2*b5)%3
 return t0+3*(t1+3*(t2+3*(t3+3*(t4+3*t5))))
def _add_st_f3(s1,s2):
 out=0;p3=1
 for _ in range(6):
  out+=((s1%3+s2%3)%3)*p3;s1//=3;s2//=3;p3*=3
 return out
def _bad_st_f3(st):
 t0=st%3;st//=3;t1=st%3;st//=3;t2=st%3;st//=3;t3=st%3;st//=3;t4=st%3;t5=st//3
 return(t0==0 and t3==0)or(t1==0 and t4==0)or(t2==0 and t5==0)
def _bad_f3(v):
 a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5=v
 return(((a0+2*a2+2*a3+a5)%3==0 and(b0+2*b2+2*b3+b5)%3==0)or((a0+2*a1+2*a3+a4+2*a5)%3==0 and(b0+2*b1+2*b3+b4+2*b5)%3==0)or((a0+2*a1+2*a2+a3+2*a4+2*a5)%3==0 and(b0+2*b1+2*b2+b3+2*b4+2*b5)%3==0))
def wheel_prices_from_kernels(path):
 s=count_kernel_segments(path)
 return{k:s.get(k,-1)for k in WHEEL_OPS}
def _fpw(v,o=0):
 return sum(FP_CL[i+o]*x for i,x in enumerate(v))
from fractions import Fraction as Frac
from math import sqrt,ceil,floor,log2
import model_bn254_v6 as M
import descent_bn254_v6 as DESC
M64=(1<<64)-1
I64_MIN= -(1<<63)
I64_MAX=(1<<63)-1
Q32=1<<32
Q38=1<<38
def _chk_u64(x):
 if x<0 or x>M64:
  raise OverflowError(f"u64 overflow: {x}")
 return x
def _chk_i64(x):
 if x<I64_MIN or x>I64_MAX:
  raise OverflowError(f"i64 overflow: {x}")
 return x
def count_kernel_segments(path):
 """Instruction lines inside each '//@seg NAME' ... '//@end' block of kernels_v26.s (labels, directives, comments, blanks excluded)."""
 import re
 ins=re.compile(r"^\s+[a-z][a-z0-9.]*\b");segs={};cur=None
 for raw in open(path,encoding="ascii"):
  st=raw.strip()
  if st.startswith("//@seg"):cur=st.split()[1];segs.setdefault(cur,0);continue
  if st.startswith("//@end"):cur=None;continue
  line=raw.split("//",1)[0].rstrip()
  if cur is None or not line.strip()or line.endswith(":")or line.strip().startswith("."):continue
  if ins.match(line):segs[cur]+=1
 return segs
def fp_prices_from_kernels(path):
 s=count_kernel_segments(path)
 out=[s["k0_body"]+s["k0_back"],s["top_2"]+s["body_2"]+s["k1_exit"]]
 for m in range(2,11):out.append(s["top_%d"%(m+1)]+s["body_%d"%(m+1)]+s["top_%d"%m])
 return out,s
def _umulh64(a,b):
 _chk_u64(a);_chk_u64(b)
 return(a*b)>>64
_SQRT_LUT=[0]*512+[ceil(sqrt(idx+1)*(2.0**52.5))+(1<<31)for idx in range(512,2048)]
def _isqrt_u40_lut7_q32(r2):
 """v25's 7-op 16 KiB LUT (kept only as the reference that LUT4 must not exceed)."""
 if r2<1 or r2>=(1<<40):
  raise OverflowError(f"u40 range error: {r2}")
 c=64-r2.bit_length()
 s=c&~1
 idx=(r2<<s)>>53
 return _chk_u64(_SQRT_LUT[idx]>>(s>>1))
def _lut4_buckets():
 """Every reachable LUT4 index i=(c<<11)|top11 of r2'=umulh(rem,inv) in [0,2^40): (i, lo, hi) of its r2' bucket."""
 out=[]
 for c in range(24,65):
  if c==64:out.append(((64<<11),0,0));continue
  for top in range(1024,2048):
   if c<=53:lo_v,hi_v=top<<(53-c),((top+1)<<(53-c))-1
   else:
    if top&((1<<(c-53))-1):continue
    lo_v=hi_v=top>>(c-53)
   out.append(((c<<11)|top,lo_v,hi_v))
 return out
_PB=dict(n=0,fail=0,maxb=-999.0,n2=0,fail2=0,maxb2=-999.0)
_LUT4=[0]*((64<<11)+1)
for _i,_lo,_hi in _lut4_buckets():_LUT4[_i]=math.isqrt((_hi+1)<<52)+1
def _isqrt_lut4(r2):
 """v26 5-op kernel (umulh; clz; lsl; extr #53; ldr): r2=r2'+1 with r2'=umulh(rem,inv); returns LUT4[(clz(r2')<<11)|((r2'<<clz)>>53)] > sqrt(r2)*2^26."""
 if r2<1 or r2>=(1<<40)+1:
  raise OverflowError(f"u40 range error: {r2}")
 r2p=r2-1;c=64-r2p.bit_length();t=(r2p<<(c&63))&M64
 return _LUT4[(c<<11)|(t>>53)]
def _staddh_u16(mem,off,w_lp):
 """Exact ARMv8.1-A FEAT_LSE STADDH Ws, [Xn] semantics: *(uint16_t*)Xn = (*(uint16_t*)Xn + (Ws & 0xFFFF)) & 0xFFFF."""
 mem[off]=(mem[off]+(w_lp&0xFFFF))&0xFFFF
_SP=math.sumprod
def _mv6(Mat,v):
 return(_SP(Mat[0],v),_SP(Mat[1],v),_SP(Mat[2],v),_SP(Mat[3],v),_SP(Mat[4],v),_SP(Mat[5],v))
def _qf6(G,v):
 return _SP(v,(_SP(G[0],v),_SP(G[1],v),_SP(G[2],v),_SP(G[3],v),_SP(G[4],v),_SP(G[5],v)))
def _eval_poly(cs,s,p):
 r=0
 for c in reversed(cs):
  r=(r*s+c)%p
 return r
def _eval_zk(v,s,p):
 inv3=pow(3,-1,p)
 s2=(s*s)%p
 s3=(s2*s)%p
 s4=(s3*s)%p
 s5=(s4*s)%p
 tot=0
 zk=M.ZK_NUM
 for k in range(6):
  vk=v[k]
  if vk:
   z=zk[k]
   ek=(z[0]+z[1]*s+z[2]*s2+z[3]*s3+z[4]*s4+z[5]*s5)%p
   tot=(tot+vk*ek)%p
 return(tot*inv3)%p
def _sqrt_mod(n,p):
 n%=p
 if n==0:return 0
 if pow(n,(p-1)//2,p)!=1:return-1
 if p%4==3:
  r=pow(n,(p+1)//4,p)
  return min(r,p-r)
 q,s=p-1,0
 while(q&1)==0:
  q>>=1;s+=1
 z=2
 while pow(z,(p-1)//2,p)!=p-1:
  z+=1
 m,c,t,r=s,pow(z,q,p),pow(n,q,p),pow(n,(q+1)//2,p)
 while t!=1:
  tt,i=t,0
  while tt!=1:
   tt=(tt*tt)%p;i+=1
  b=pow(c,1<<(m-i-1),p)
  m=i;c=(b*b)%p;t=(t*c)%p;r=(r*b)%p
 return min(r,p-r)
def _find_all_roots_h(p):
 def pmul(a,b):
  la,lb=len(a),len(b)
  r=[0]*(la+lb-1)
  for i in range(la):
   xv=a[i]
   if xv:
    for j in range(lb):r[i+j]+=xv*b[j]
  return[x%p for x in r]
 def pmod(a,m):
  a=a[:];lm=len(m)
  inv_m=1 if m[-1]==1 else pow(m[-1],-1,p)
  while len(a)>=lm:
   c=a[-1]if inv_m==1 else(a[-1]*inv_m)%p
   off=len(a)-lm
   for i in range(lm):a[off+i]=(a[off+i]-c*m[i])%p
   while a and a[-1]==0:a.pop()
  return a
 def pdiv(a,b):
  a=a[:];lb=len(b);inv_b=pow(b[-1],-1,p)
  q=[0]*(len(a)-lb+1)
  while len(a)>=lb:
   c=(a[-1]*inv_b)%p;off=len(a)-lb;q[off]=c
   for i in range(lb):a[off+i]=(a[off+i]-c*b[i])%p
   while a and a[-1]==0:a.pop()
  while q and q[-1]==0:q.pop()
  return q
 def ppow(b,e,m):
  if len(m)==7 and m==h:
   res=M._ppow_h(b,e,p)
   while res and res[-1]==0:res.pop()
   return res
  r=[1]
  while e:
   if e&1:r=pmod(pmul(r,b),m)
   b=pmod(pmul(b,b),m);e>>=1
  return r
 def pgcd(a,b):
  while b:a,b=b,pmod(a,b)
  if a:
   inv=pow(a[-1],-1,p);a=[(xv*inv)%p for xv in a]
  return a
 h=[c%p for c in M.H_COEFFS]
 xp=ppow([0,1],p,h)
 if len(xp)<2:xp+=[0]*(2-len(xp))
 xp[1]=(xp[1]-1)%p
 while xp and xp[-1]==0:xp.pop()
 g=pgcd(h,xp)
 if len(g)<4:return[]
 rng_cz=random.Random(p)
 stk=[g];roots=[]
 while stk:
  f=stk.pop()
  if len(f)==2:
   roots.append((-f[0]*pow(f[1],-1,p))%p)
  elif len(f)==3:
   inv2=pow(2*f[2],-1,p)
   sq=_sqrt_mod((f[1]*f[1]-4*f[2]*f[0])%p,p)
   roots.append(((-f[1]+sq)*inv2)%p)
   roots.append(((-f[1]-sq)*inv2)%p)
  else:
   while True:
    a=rng_cz.randrange(p)
    t=ppow([a,1],(p-1)//2,f)
    if t:t[0]=(t[0]-1)%p
    while t and t[-1]==0:t.pop()
    d=pgcd(f,t)
    if 1<len(d)<len(f):
     stk.append(d);stk.append(pdiv(f,d));break
 roots.sort()
 return roots
def _find_root_h(p):
 rs=_find_all_roots_h(p)
 return rs[0]if rs else None
SH_A=(20,10,0,50,40,30)
BI_A_STEP=(512,512,512,512,512,512)
BI_A_CAN=(512,512,512,0,512,512)
MK_A=(0x3FF,)*6
SH_B=(51,30,20,10,0,41)
BI_B_STEP=(1024,1024,512,512,512,512)
BI_B_CAN=(0,1024,512,512,512,512)
MK_B=(0x7FF,0x7FF,0x3FF,0x3FF,0x3FF,0x3FF)
A_ORD=(3,4,5,0,1,2)
B_ORD=(0,5,1,2,3,4)
BIAS_A_STEP_X2=sum((2*BI_A_STEP[i])<<SH_A[i]for i in range(6))
A3_BIAS_SHIFT=512<<50
BIAS_B_STEP_X2=sum((2*BI_B_STEP[i])<<SH_B[i]for i in range(6))
B0_BIAS_SHIFT=1024<<51
def _pack_v(v,sh,bi):
 p=0
 for i in range(6):
  p|=(v[i]+bi[i])<<sh[i]
 return p
def _pack_d(d,sh):
 acc=0
 for i in range(6):
  acc=(acc+(d[i]<<sh[i]))&M64
 return acc
def _unpack_v(p,sh,bi,mk):
 return tuple(((p>>sh[i])&mk[i])-bi[i]for i in range(6))
def _cone_of(v):
 x0,x1,x2=v[3],-v[4],v[5]
 sx=x0+x1+x2
 if sx==0:return-1
 if sx>0:sa_pos=1
 else:sa_pos=0;x0,x1,x2= -x0,-x1,-x2
 if x0>x1 and x0>x2:k,bk0=0,v[6]
 elif x1>x0 and x1>x2:k,bk0=1,v[6]+v[8]+v[10]
 elif x2>x0 and x2>x1:k,bk0=2,v[6]+v[7]+v[10]
 else:return-1
 if bk0==0:return-1
 return(k<<2)|(sa_pos<<1)|int(bk0>0)
def _canon_cone_ref(a,b,S,S2):
 sx=a[3]-a[4]+a[5]
 sas=(1,)if sx>0 else((-1,)if sx<0 else(1,-1))
 best_key=None;best_ab=None
 for ta,tb in((a,b),(_mv6(S,a),_mv6(S,b)),(_mv6(S2,a),_mv6(S2,b))):
  sb=1
  for i in B_ORD:
   if tb[i]!=0:
    sb=1 if tb[i]>0 else-1;break
  cb_t=tuple(sb*x for x in tb)
  bkey=(cb_t[0],cb_t[1],cb_t[2],cb_t[3],cb_t[4],cb_t[5])
  for sa in sas:
   ca_t=ta if sa==1 else(-ta[0],-ta[1],-ta[2],-ta[3],-ta[4],-ta[5])
   key=(sa*sx,ca_t[3],ca_t[3],ca_t[4],ca_t[5],ca_t[0],ca_t[1],ca_t[2])+bkey
   if best_key is None or key>best_key:
    best_key=key
    best_ab=(ca_t,cb_t)
 return best_ab
def _floor_div(num,den):
 if den<0:num,den= -num,-den
 return num//den
def _ceil_div(num,den):
 if den<0:num,den= -num,-den
 return-((-num)//den)
def _certify_min_gsym(G_sym,G_SC):
 """Exact 6D Fraction LDL^T certificate proving min_{0!=v in Z^6} v^T G_sym v == G_sym[0][0] == 753729."""
 D_f=[Frac(0)]*6
 L_f=[[Frac(1 if i==j else 0)for j in range(6)]for i in range(6)]
 for c in range(6):
  d=Frac(G_sym[c][c])-sum(L_f[c][k]*L_f[c][k]*D_f[k]for k in range(c))
  assert d>0;D_f[c]=d
  for r in range(c+1,6):
   L_f[r][c]=(Frac(G_sym[r][c])-sum(L_f[r][k]*L_f[c][k]*D_f[k]for k in range(c)))/d
 T_lim=5*G_SC
 v=[0]*6;best=[10**18];cnt=[0]
 def enum_q(idx,rem):
  sh=sum(L_f[j][idx]*v[j]for j in range(idx+1,6))
  rad2=rem/D_f[idx]
  r_hi=Frac(math.isqrt(rad2.numerator*rad2.denominator)+1,rad2.denominator)
  di=D_f[idx]
  for xv in range((-sh-r_hi).__floor__(),(-sh+r_hi).__ceil__()+1):
   d=xv+sh;nrem=rem-di*d*d
   if nrem>=0:
    v[idx]=xv
    if idx==0:
     if any(v):
      q=_qf6(G_sym,v);cnt[0]+=1
      if q<best[0]:best[0]=q
    else:
     enum_q(idx-1,nrem)
 enum_q(5,Frac(T_lim))
 return best[0],cnt[0]
def verify_v20_lat_b2_and_srt(seed=20261121):
 T4=M.load_data()["t4_geom"]
 G_sym,D_mats,UMATS,S,S2,G_SC=M.build_exact_int_geom(T4)
 smat_cyclic_ok=(S[3]==[0,0,0,0,-1,0]and S[4]==[0,0,0,0,0,-1]and S[5]==[0,0,0,1,0,0])
 m6_exact,m6_pts=_certify_min_gsym(G_sym,G_SC)
 min_b0_sq=1000*m6_exact
 b_star_min_f=Frac(min_b0_sq)*(Frac(13,20)**11)
 T12_prod_max=int(round(1764*(595.0**2)*G_SC))
 uncond_inv_max=int((Frac(Q38*T12_prod_max)/b_star_min_f).__ceil__())
 uncond_r2_max=((Q38*uncond_inv_max)>>64)+1
 uncond_rk_max=math.isqrt(int((Frac(T12_prod_max)/b_star_min_f).__ceil__()))+1
 uncond_xmax0=int((Frac(3,2)**11*uncond_rk_max).__ceil__())
 uncond_ck_max=11*((1<<31)+1)*uncond_xmax0
 uncond_rad_q32_max=_isqrt_lut4(uncond_r2_max)
 uncond_d19_max=(uncond_rad_q32_max+Q32)>>13
 uncond_d2_max=uncond_d19_max*uncond_d19_max
 _b4=sorted(_lut4_buckets(),key=lambda e:e[1])
 lut_mono_ok=all(_LUT4[_b4[j][0]]<=_LUT4[_b4[j+1][0]]for j in range(len(_b4)-1))and all(_b4[j][2]+1==_b4[j+1][1]for j in range(len(_b4)-1))and _b4[0][1]==0 and _b4[-1][2]==(1<<40)-1
 lut4_bucket_viol=sum(1 for _i,_lo,_hi in _b4 if _LUT4[_i]*_LUT4[_i]<=(_hi+1)<<52)
 lut4_gt_lut7=sum(1 for _i,_lo,_hi in _b4 if _LUT4[_i]>_isqrt_u40_lut7_q32(_lo+1))
 lut4_max_over=max(_LUT4[_i]/(sqrt(_lo+1.0)*67108864.0)-1.0 for _i,_lo,_hi in _b4)
 uncond_64bit_zero_exc_ok=(
  m6_exact==G_sym[0][0]==753729
  and lut_mono_ok
  and b_star_min_f>32*G_SC
  and uncond_inv_max<=M64
  and uncond_r2_max<(1<<40)
  and uncond_ck_max<(1<<63)
  and uncond_d19_max<(1<<32)
  and uncond_d2_max<=M64
 )
 rng_sq=random.Random(seed^0x51A7)
 isqrt_errs=0;isqrt_max_over=0.0
 for n_t in list(range(1,2049))+[(1<<k)-1 for k in range(8,41)]+[(1<<k)for k in range(8,40)]+[1+M.rand_below(rng_sq.getrandbits,(1<<40)-1)for _ in range(19000)]:
  rad_t=_isqrt_lut4(n_t);ex_t=sqrt(n_t)*(1<<26)
  if rad_t<=ex_t:isqrt_errs+=1
  ov_t=(rad_t-ex_t)/ex_t
  if ov_t>isqrt_max_over:isqrt_max_over=ov_t
 staddh_errs=0
 mem_lse=[0]*256;mem_ref=[0]*256
 for i_t in range(4096):
  off_t=rng_sq.randrange(256)
  val_t=0xFFFF if i_t<16 else(0x8000 if i_t<32 else rng_sq.randrange(65536))
  _staddh_u16(mem_lse,off_t,val_t)
  w_r=mem_ref[off_t];w_r=(w_r+(val_t&0xFFFF))&0xFFFF;mem_ref[off_t]=w_r
  if mem_lse[off_t]!=mem_ref[off_t]:staddh_errs+=1
 orbits=[]
 p_cand=101
 while len(orbits)<6:
  if M._isp(p_cand)and p_cand not in(2,3,19):
   rts=[s for s in range(p_cand)if _eval_poly(M.H_COEFFS,s,p_cand)==0]
   if len(rts)>=3:
    used=set();inv3=pow(3,-1,p_cand)
    for s0 in rts:
     if s0 in used:continue
     s1=(_eval_poly(M.SIGMA_NUM,s0,p_cand)*inv3)%p_cand
     s2=(_eval_poly(M.SIGMA_NUM,s1,p_cand)*inv3)%p_cand
     if len({s0,s1,s2})!=3:continue
     used.update((s0,s1,s2))
     ws=_eval_poly(M.W_COEFFS,s0,p_cand)
     assert _eval_poly(M.W_COEFFS,s1,p_cand)==ws and _eval_poly(M.W_COEFFS,s2,p_cand)==ws
     uw=(-(((M.U_BN254+1)%p_cand)*ws+2))%p_cand
     if uw==0:continue
     rho=_sqrt_mod(pow(uw,-1,p_cand),p_cand)
     if rho>0:
      orbits.append((p_cand,(s0,s1,s2),rho))
      if len(orbits)==6:break
  p_cand+=2
 def in_jvor(a,b,Wa,Wb):
  for fi in range(6):
   q=Wa*_qf6(D_mats[fi],a)+Wb*_qf6(D_mats[fi],b)
   if q<0:return False
   if q==0:
    ua,ub=_mv6(UMATS[fi],a),_mv6(UMATS[fi],b)
    if _canon_cone_ref(ua,ub,S,S2)<_canon_cone_ref(a,b,S,S2):
     return False
  return True
 def stab_of(a,b):
  c0=_canon_cone_ref(a,b,S,S2);st=0
  for ta,tb in((a,b),(_mv6(S,a),_mv6(S,b)),(_mv6(S2,a),_mv6(S2,b))):
   if _canon_cone_ref(ta,tb,S,S2)==c0:
    if ta in(a,tuple(-x for x in a))and tb in(b,tuple(-x for x in b)):st+=1
  return st
 Wa_t,Wb_t=1,1;T12_t=18*G_SC;pts6=[]
 for a5 in range(-2,3):
  for a4 in range(-2,3):
   for a3 in range(-2,3):
    for a2 in range(-2,3):
     for a1 in range(-2,3):
      for a0 in range(-3,4):
       v=(a0,a1,a2,a3,a4,a5)
       if not any(v):continue
       qv=_qf6(G_sym,v)
       if qv<=T12_t:pts6.append((v,qv))
 raw_reg={};jvor_reg={}
 for va,qa in pts6:
  for vb,qb in pts6:
   if qa+qb>T12_t:continue
   ca,cb=_canon_cone_ref(va,vb,S,S2)
   if(ca,cb)not in raw_reg:
    st=stab_of(ca,cb)
    raw_reg[(ca,cb)]=st
    if in_jvor(ca,cb,Wa_t,Wb_t):jvor_reg[(ca,cb)]=st
 lat_mismatches=0;tot_raw_hits=0;tot_jvor_hits=0
 for p,(s0,s1,s2),rho in orbits:
  bf_raw={};bf_jvor={}
  for(ca,cb),st in raw_reg.items():
   mult=0
   for sk in(s0,s1,s2):
    av=_eval_zk(ca,sk,p);rb=(rho*_eval_zk(cb,sk,p))%p
    if av==rb:mult+=1
    if av==(p-rb)%p:mult+=1
   if st==3:mult//=3
   if mult:
    bf_raw[(ca,cb)]=mult
    if(ca,cb)in jvor_reg:bf_jvor[(ca,cb)]=mult
  inv3=pow(3,-1,p);sp=[1]
  for _ in range(5):sp.append((sp[-1]*s0)%p)
  es=[(sum(M.ZK_NUM[k][j]*sp[j]for j in range(6))*inv3)%p for k in range(6)]
  lf=es+[((-rho)*e)%p for e in es]
  inv0=pow(lf[0],-1,p)
  B=[[p if c==0 else 0 for c in range(12)]]
  for i in range(1,12):
   c0=(-lf[i]*inv0)%p
   if c0>p//2:c0-=p
   B.append([c0 if c==0 else(1 if c==i else 0)for c in range(12)])
  DESC.bareiss_lll_12d(B,G_sym,Wa_t,Wb_t,False,9,10)
  GB=[[sum(B[i][u]*G_sym[u][v]*B[j][v]for u in range(6)for v in range(6))
    +sum(B[i][6+u]*G_sym[u][v]*B[j][6+v]for u in range(6)for v in range(6))
    for j in range(12)]for i in range(12)]
  D_fp=[1]+[0]*12;lam_fp=[[0]*12 for _ in range(12)]
  for i in range(12):
   for j in range(i+1):
    u=GB[i][j]
    for k in range(j):u=(D_fp[k+1]*u-lam_fp[i][k]*lam_fp[j][k])//D_fp[k]
    if j<i:lam_fp[i][j]=u
    else:D_fp[i+1]=u
  mu_lo=[[_chk_i64((lam_fp[j][k]<<32)//D_fp[k+1])if j>k else 0 for k in range(12)]for j in range(12)]
  b2_lo=[_chk_u64(min(M64,(D_fp[k+1]<<64)//(D_fp[k]*T12_t)))for k in range(12)]
  inv_b2=[_chk_u64(((D_fp[k]*T12_t<<38)+D_fp[k+1]-1)//D_fp[k+1])for k in range(12)]
  Xm=[0]*12
  for k in range(11,-1,-1):
   rk=math.isqrt(((D_fp[k]*T12_t)+D_fp[k+1]-1)//D_fp[k+1])+1
   Xm[k]=rk+sum((Xm[j]+1)//2 for j in range(k+1,12))
  ec_arr=[sum(Xm[j]for j in range(k+1,12))+1 for k in range(12)]
  fp_raw={};fp_jvor={};x=[0]*12
  def enum64(k,rem,az,cur):
   ck=0
   for j in range(k+1,12):ck=_chk_i64(ck-_chk_i64(mu_lo[j][k]*x[j]))
   r2=_chk_u64(_umulh64(rem,inv_b2[k])+1)
   rad=_isqrt_lut4(r2)
   ec=ec_arr[k]
   lo=(1 if k==0 else 0)if az else(_chk_i64(ck-ec-rad+Q32-1)>>32)
   hi=_chk_i64(ck+ec+rad)>>32
   Bk=B[k];b2k=b2_lo[k]
   for xv in range(lo,hi+1):
    if k==0:
     nxt=[cur[c]+xv*Bk[c]for c in range(12)]
     av,bv=tuple(nxt[:6]),tuple(nxt[6:])
     if any(av)and any(bv)and _qf6(G_sym,av)+_qf6(G_sym,bv)<=T12_t:
      ca,cb=_canon_cone_ref(av,bv,S,S2)
      fp_raw[(ca,cb)]=fp_raw.get((ca,cb),0)+1
      if in_jvor(ca,cb,Wa_t,Wb_t):fp_jvor[(ca,cb)]=fp_jvor.get((ca,cb),0)+1
    else:
     ad=_chk_u64(abs(_chk_i64((xv<<32)-ck)))
     d_lo=ad-ec if ad>ec else 0
     d19=d_lo>>13
     sub=_umulh64(_chk_u64(d19*d19),b2k)
     if sub<=rem:
      x[k]=xv
      enum64(k-1,rem-sub,az and(xv==0),[cur[c]+xv*Bk[c]for c in range(12)])
  enum64(11,Q38,True,[0]*12)
  if fp_raw!=bf_raw or fp_jvor!=bf_jvor:lat_mismatches+=1
  tot_raw_hits+=sum(fp_raw.values());tot_jvor_hits+=sum(fp_jvor.values())
 def _canon_sign_b(tb):
  sb=1
  for i in B_ORD:
   if tb[i]!=0:sb=1 if tb[i]>0 else-1;break
  return tuple(sb*x for x in tb)
 def _zerosplit_4plane(a,b):
  """Exact 4-plane B_TAB[0..3] lookup: a_tab[pack_a] = ((2*ai)<<5)|k (<=54 bits for all ai<2^48; k in {0,1,2} if S*a!=a else k=3),
  B_TAB[k][pack_b] = canon_sign(S^k b) for k in {0,1,2}, and canon_S_sign(b) for k=3."""
  a1,a2=_mv6(S,a),_mv6(S2,a)
  sx=a[3]-a[4]+a[5]
  sas=(1,)if sx>0 else((-1,)if sx<0 else(1,-1))
  if a1==a:
   # S-fixed a: k=3, ca = best +-a, cb = max over S^j of canon_sign(S^j b)
   ca_best=max(((sa*sx,sa*a[3],sa*a[3],sa*a[4],sa*a[5],sa*a[0],sa*a[1],sa*a[2]),tuple(sa*x for x in a))for sa in sas)[1]
   cb_best=max(_canon_sign_b(tb)for tb in(b,_mv6(S,b),_mv6(S2,b)))
   return 3,ca_best,cb_best
  best_key=None;best_k=0;best_ca=None
  for k_i,ta in((0,a),(1,a1),(2,a2)):
   for sa in sas:
    ca_t=ta if sa==1 else(-ta[0],-ta[1],-ta[2],-ta[3],-ta[4],-ta[5])
    key=(sa*sx,ca_t[3],ca_t[3],ca_t[4],ca_t[5],ca_t[0],ca_t[1],ca_t[2])
    if best_key is None or key>best_key:best_key=key;best_k=k_i;best_ca=ca_t
  tb=b if best_k==0 else(_mv6(S,b)if best_k==1 else _mv6(S2,b))
  return best_k,best_ca,_canon_sign_b(tb)
 def execute_chord_split(base,B0,l0,h0,T12_chk=None,verify_ref=True,p_div=None,s0_div=None,rho_div=None):
  """Zero-Split 4-plane B_TAB[0..3] + 2-word Z/2^64Z mixed-radix emitter: 1 segment, 0 clips, 0 splits, 0 boundary leaves."""
  if l0>h0:
   return dict(lvs=0,fast=0,bnd=0,segs=0,bspl=0,fspl=0,cross=0,clips=0,mm=0,ae=0,div_err=0,max_sa=0)
  length=h0-l0+1
  vl=[base[c]+l0*B0[c]for c in range(12)]
  pa=(_mr_pack_a(vl[:6]))&M64;pb5=(_mr_pack_b(vl[6:])<<5)&M64
  da=_mr_da(B0[:6]);db5=_mr_db5(B0[6:])
  mm=ae=div_err=0
  _zke=None
  if p_div is not None:
   _i3=pow(3,-1,p_div);_sp_d=[1]
   for _ in range(5):_sp_d.append((_sp_d[-1]*s0_div)%p_div)
   _zke=[(sum(M.ZK_NUM[k][j]*_sp_d[j]for j in range(6))*_i3)%p_div for k in range(6)]
   _zke_r=[(-rho_div*e)%p_div for e in _zke];_zke12=_zke+_zke_r;_sp_div=math.sumprod
  for u in range(length):
   if verify_ref or p_div is not None:
    v=[vl[c]+u*B0[c]for c in range(12)]
    if _zke is not None and _sp_div(v,_zke12)%p_div!=0:div_err+=1
    av,bv=tuple(v[:6]),tuple(v[6:])
    if any(abs(av[i])>AMAX_MR[i]for i in range(6))or any(abs(bv[i])>BMAX_MR[i]for i in range(6)):
     ae+=1
    else:
     pa_exp=_mr_pack_a(av);pb_exp=_mr_pack_b(bv)
     if pa!=pa_exp or(pb5&31)!=0 or(pb5>>5)!=pb_exp or pa>=PROD_MA_MR or pb_exp>=PROD_MB_MR or pb_exp>=(1<<59):
      ae+=1
     if verify_ref:
      k_p,ca_g,cb_g=_zerosplit_4plane(av,bv)
      # Verify exact ARM64 `extr x6, x0, x5, #5` plane-index formation!
      x6_sim=(((k_p<<64)|pb5)>>5)&M64
      if x6_sim!=((k_p<<59)|pb_exp):ae+=1
      ca_r,cb_r=_canon_cone_ref(av,bv,S,S2)
      if ca_g!=ca_r or cb_g!=cb_r:mm+=1
   pa=(pa+da)&M64;pb5=(pb5+db5)&M64
  return dict(lvs=length,fast=length,bnd=0,segs=1,bspl=0,fspl=0,cross=0,clips=0,mm=mm,ae=ae,div_err=div_err,max_sa=1)
 tie_cases=[
  ((1,0,0,0,0,0),(0,1,0,0,0,0)),
  ((5,0,0,0,0,0),(1,2,-3,4,-1,2)),
  ((-7,0,0,0,0,0),(0,0,1,-2,3,-4)),
  ((1,1,1,0,0,0),(0,1,-1,2,0,-1)),
 ]
 for v_seed,b_seed in[((0,1,0,0,0,0),(0,1,0,0,0,0)),((0,0,0,1,0,0),(1,-2,3,0,1,-1)),((1,0,1,0,1,0),(0,0,1,0,0,-2))]:
  v1=_mv6(S,v_seed);v2=_mv6(S2,v_seed)
  a_fix=tuple(v_seed[i]+v1[i]+v2[i]for i in range(6))
  assert _mv6(S,a_fix)==a_fix
  tie_cases.append((a_fix,b_seed))
 tie_errs=0
 for ta,tb in tie_cases:
  outs=set();ca_i,cb_i=ta,tb
  for _ in range(3):
   for sa_i in(1,-1):
    for sb_i in(1,-1):
     outs.add(_canon_cone_ref(tuple(sa_i*x for x in ca_i),tuple(sb_i*x for x in cb_i),S,S2))
   ca_i,cb_i=_mv6(S,ca_i),_mv6(S,cb_i)
  if len(outs)!=1:tie_errs+=1
  if ta[3]==0 and ta[4]==0 and ta[5]==0:
   rt=execute_chord_split(list(ta)+list(tb),[1,0,0,0,0,0,0,1,0,0,0,0],0,4,verify_ref=True)
   if rt["mm"]or rt["ae"]:tie_errs+=1
 rng=random.Random(seed)
 T12_prod=int(round(1764*(589.43**2)*G_SC))
 omega=[math.pi**(d/2.0)/math.gamma(d/2.0+1.0)for d in range(13)]
 det_G12_sqrt=sqrt((1764.0**6)*(1000.0**6)*(M._det_mat(T4["GA"])**2))*(G_SC**6)
 pbits_list=(49,51,53,55,57,59,61)
 comp_tot_lvs=0;comp_tot_nodes=0;comp_mm=0;comp_ae=0;comp_div_err=0;comp_max_sa=0;wheel_tree_state_errs=0;wheel_tree_leaves_checked=0;swar21_cur_errs=0
 comp_node_ratios_obs=[];comp_node_ratios_pred=[]
 comp_fp_ir_obs=[];comp_fp_ir_pred=[]
 prod_exec_lvs=0;prod_exec_mm=0;prod_exec_ae=0;prod_exec_div_err=0;prod_exec_clips=0
 prod_anal_ir_vals=[];prod_min_b2_bits=999.0
 prod_max_inv_bits=0;prod_max_r2_bits=0;prod_max_d2_bits=0;prod_max_ck_bits=0
 inc_c0_errs=0
 for pbits_t in pbits_list:
  pp=(1<<pbits_t)+1
  while True:
   if M._isp(pp):
    s0p=_find_root_h(pp)
    if s0p is not None:
     wsp=_eval_poly(M.W_COEFFS,s0p,pp)
     uwp=(-(((M.U_BN254+1)%pp)*wsp+2))%pp
     if uwp!=0:
      rhop=_sqrt_mod(pow(uwp,-1,pp),pp)
      if rhop>0:break
   pp+=2
  inv3p=pow(3,-1,pp);spp=[1]
  for _ in range(5):spp.append((spp[-1]*s0p)%pp)
  esp=[(sum(M.ZK_NUM[k][j]*spp[j]for j in range(6))*inv3p)%pp for k in range(6)]
  lfp=esp+[((-rhop)*e)%pp for e in esp]
  inv0p=pow(lfp[0],-1,pp)
  Bp=[[pp if c==0 else 0 for c in range(12)]]
  for i in range(1,12):
   c0=(-lfp[i]*inv0p)%pp
   if c0>pp//2:c0-=pp
   Bp.append([c0 if c==0 else(1 if c==i else 0)for c in range(12)])
  DESC.bareiss_lll_12d(Bp,G_sym,1764,1000,False,9,10,(99,100))
  Bpa=[(_mv6(G_sym,r[:6]),_mv6(G_sym,r[6:]))for r in Bp]
  GBp=[[1764*math.sumprod(Bpa[i][0],Bp[j][:6])+1000*math.sumprod(Bpa[i][1],Bp[j][6:])for j in range(12)]for i in range(12)]
  Dp=[1]+[0]*12;lamp=[[0]*12 for _ in range(12)]
  for i in range(12):
   for j in range(i+1):
    u=GBp[i][j]
    for k in range(j):u=(Dp[k+1]*u-lamp[i][k]*lamp[j][k])//Dp[k]
    if j<i:lamp[i][j]=u
    else:Dp[i+1]=u
  min_b2=min(Dp[k+1]/Dp[k]for k in range(12))
  prod_min_b2_bits=min(prod_min_b2_bits,log2(min_b2))
  assert min_b2>=float(b_star_min_f)*0.999999
  T_sc=int((4200.0*pp*det_G12_sqrt/omega[12])**(1.0/6.0))
  mu_s=[[_chk_i64((lamp[j][k]<<32)//Dp[k+1])if j>k else 0 for k in range(12)]for j in range(12)]
  b2_s=[_chk_u64(min(M64,(Dp[k+1]<<64)//(Dp[k]*T_sc)))for k in range(12)]
  inv_s=[_chk_u64(((Dp[k]*T_sc<<38)+Dp[k+1]-1)//Dp[k+1])for k in range(12)]
  Xm_s=[0]*12
  for k in range(11,-1,-1):
   rk=math.isqrt(((Dp[k]*T_sc)+Dp[k+1]-1)//Dp[k+1])+1
   Xm_s[k]=rk+sum((Xm_s[j]+1)//2 for j in range(k+1,12))
  ec_s=[sum(Xm_s[j]for j in range(k+1,12))+1 for k in range(12)]
  nk_s=[0]*12;fl_s=[0]*12;xs=[0]*12;orb_lvs=[0]
  mu_s_col=[[mu_s[j][k]for j in range(k+1,12)]for k in range(12)]
  A1_fd_s=(Dp[2]<<38)//(Dp[1]*T_sc);d2r1_fd_s=2*A1_fd_s
  _sp_f=math.sumprod
  _B21=[(_mr_da(r[:6]),_mr_db5(r[6:]))for r in Bp];_BI21=(A_OFF_MR&M64,(B_OFF_MR<<5)&M64)
  _stb=[_st_f3(Bp[j])for j in range(12)]
  _step_k=[None]+[[[_add_st_f3(st,_add_st_f3(_stb[j],_stb[j])if r_m==2 else(_stb[j]if r_m==1 else 0))for st in range(729)]for r_m in range(3)]for j in range(1,12)]
  _step1_tab=[(( (1 if _bad_st_f3(st)else 0)|(2 if _bad_st_f3(_add_st_f3(st,_stb[0]))else 0)|(4 if _bad_st_f3(_add_st_f3(st,_add_st_f3(_stb[0],_stb[0])))else 0) )<<10)|_step_k[1][1][st]for st in range(729)]
  def full_tree(k,rem,az,cur,c0_passed=None,st_cur=0,bmask_passed=None,c21=None):
   nonlocal comp_mm,comp_ae,comp_div_err,comp_max_sa,inc_c0_errs,wheel_tree_state_errs,wheel_tree_leaves_checked,swar21_cur_errs
   nk_s[k]+=1
   ck= -_sp_f(mu_s_col[k],xs[k+1:])if k<11 else 0
   if ck<I64_MIN or ck>I64_MAX:raise OverflowError(f"i64 ck overflow: {ck}")
   if k==0 and c0_passed is not None and c0_passed!=ck:inc_c0_errs+=1
   r2=((rem*inv_s[k])>>64)+1
   if r2<1 or r2>=(1<<40):raise OverflowError(f"u40 r2 overflow: {r2}")
   rad=_isqrt_lut4(r2)
   ec=ec_s[k]
   lo=(1 if k==0 else 0)if az else((ck-ec-rad+Q32-1)>>32)
   hi=(ck+ec+rad)>>32
   if lo>hi:return
   Bk=Bp[k];B21k=_B21[k]
   if k==0:
    if c21 is not None:
     pa_l=(c21[0]+lo*B21k[0])&M64;pb5_l=(c21[1]+lo*B21k[1])&M64
     vl_chk=[cur[c]+lo*Bk[c]for c in range(12)]
     if pa_l!=_mr_pack_a(vl_chk[:6])or(pb5_l&31)!=0 or(pb5_l>>5)!=_mr_pack_b(vl_chk[6:]):swar21_cur_errs+=1
    if bmask_passed is not None:
     for u_w in range(lo,hi+1):
      wheel_tree_leaves_checked+=1
      if bool((bmask_passed>>(u_w%3))&1)!=_bad_f3([cur[c]+u_w*Bk[c]for c in range(12)]):
       wheel_tree_state_errs+=1
    rc=execute_chord_split(cur,Bk,lo,hi,T12_chk=T_sc,verify_ref=(orb_lvs[0]<400),p_div=pp,s0_div=s0p,rho_div=rhop)
    orb_lvs[0]+=rc["lvs"];comp_mm+=rc["mm"];comp_ae+=rc["ae"];comp_div_err+=rc["div_err"]
    if rc["max_sa"]>comp_max_sa:comp_max_sa=rc["max_sa"]
   elif k==1:
    R1=hi-lo+1
    pad1=((2*R1*ec*A1_fd_s)>>32)+2+R1*R1
    d0_q32=(lo<<32)-ck
    d19_0=abs(d0_q32)>>13
    f0=((d19_0*d19_0)*b2_s[1])>>64
    rem0=rem+pad1-f0
    dr1=(((2*lo+1)*Q32-2*ck)*A1_fd_s)>>32
    c0_cur= -_sp_f(mu_s_col[0][1:],xs[2:])-mu_s[1][0]*lo
    x23_sim=((hi-lo)<<12)|_step_k[1][lo%3][st_cur]
    for xv in range(lo,hi+1):
     ent14=_step1_tab[x23_sim&0x3FF]
     x23_sim=(x23_sim&~0x3FF)|(ent14&0x3FF)
     fl_s[0]+=(rem0<0)
     if rem0>=0:
      xs[1]=xv
      full_tree(0,rem0,az and(xv==0),[cur[c]+xv*Bk[c]for c in range(12)],c0_cur,0,ent14>>10,((c21[0]+xv*B21k[0])&M64,(c21[1]+xv*B21k[1])&M64))
     x23_sim=(x23_sim-(1<<12))&M64
     rem0-=dr1;dr1+=d2r1_fd_s;c0_cur-=mu_s[1][0]
   else:
    b2k=b2_s[k]
    st_child=_step_k[k][lo%3][st_cur]
    for xv in range(lo,hi+1):
     df_x=(xv<<32)-ck
     ad= -df_x if df_x<0 else df_x
     d_lo=ad-ec if ad>ec else 0
     d19=d_lo>>13
     sub=((d19*d19)*b2k)>>64
     fl_s[k-1]+=(sub>rem)
     if sub<=rem:
      xs[k]=xv
      full_tree(k-1,rem-sub,az and(xv==0),[cur[c]+xv*Bk[c]for c in range(12)],None,st_child,None,((c21[0]+xv*B21k[0])&M64,(c21[1]+xv*B21k[1])&M64))
     st_child=_step_k[k][1][st_child]
  full_tree(11,Q38,True,[0]*12,None,0,None,_BI21)
  comp_tot_lvs+=orb_lvs[0];comp_tot_nodes+=sum(nk_s)
  comp_node_ratios_obs.append(round(sum(nk_s)/max(1,orb_lvs[0]),3))
  pred_lvls=[2.0**(log2(omega[12-k]/omega[12])+0.5*(log2(Dp[k])-k*log2(T_sc)))for k in range(1,12)]
  comp_node_ratios_pred.append(round(sum(pred_lvls),3))
  comp_fp_ir_obs.append(round((_fpw([a+b for a,b in zip(nk_s,fl_s)]))/max(1,orb_lvs[0]),2))
  comp_fp_ir_pred.append(round(_fpw(pred_lvls),2))
  mu_p=[[_chk_i64((lamp[j][k]<<32)//Dp[k+1])if j>k else 0 for k in range(12)]for j in range(12)]
  b2_p=[_chk_u64(min(M64,(Dp[k+1]<<64)//(Dp[k]*T12_prod)))for k in range(12)]
  _mb=max(log2(Dp[k+1])-log2(Dp[k])-log2(T12_prod)for k in range(12));_PB['n2']+=1;_PB['fail2']+=(_mb>=0);_PB['maxb2']=max(_PB['maxb2'],_mb)
  inv_p=[_chk_u64(((Dp[k]*T12_prod<<38)+Dp[k+1]-1)//Dp[k+1])for k in range(12)]
  for ib in inv_p:prod_max_inv_bits=max(prod_max_inv_bits,ib.bit_length())
  Xm_p=[0]*12
  for k in range(11,-1,-1):
   rk=math.isqrt(((Dp[k]*T12_prod)+Dp[k+1]-1)//Dp[k+1])+1
   Xm_p[k]=rk+sum((Xm_p[j]+1)//2 for j in range(k+1,12))
  ec_p=[sum(Xm_p[j]for j in range(k+1,12))+1 for k in range(12)]
  Bnorm_f=[Dp[k+1]/Dp[k]for k in range(12)];R_f=sqrt(float(T12_prod))
  p_lvs=0;p_chords=0;p_fast=0;p_bnd=0;p_segs=0;p_bspl=0;p_fspl=0;p_cross=0;p_clips=0
  while p_lvs<300:
   g=[rng.gauss(0.0,1.0)for _ in range(11)]
   n2=sum(z*z for z in g);rad_g=R_f*(rng.random()**(1.0/11.0))/sqrt(n2)
   xp=[0]*12;rem_p=Q38;ok_s=True
   for k in range(11,0,-1):
    ck=0
    for j in range(k+1,12):ck=_chk_i64(ck-_chk_i64(mu_p[j][k]*xp[j]))
    prod_max_ck_bits=max(prod_max_ck_bits,abs(ck).bit_length())
    xp[k]=int(round(ck/float(Q32)+(g[k-1]*rad_g)/sqrt(Bnorm_f[k])))
    ad=_chk_u64(abs(_chk_i64((xp[k]<<32)-ck)))
    d_lo=ad-ec_p[k]if ad>ec_p[k]else 0
    d19=d_lo>>13;d2=_chk_u64(d19*d19)
    prod_max_d2_bits=max(prod_max_d2_bits,d2.bit_length())
    sub=_umulh64(d2,b2_p[k])
    if sub>rem_p:ok_s=False;break
    rem_p-=sub
   if not ok_s:continue
   ck0=0
   for j in range(1,12):ck0=_chk_i64(ck0-_chk_i64(mu_p[j][0]*xp[j]))
   r20=_chk_u64(_umulh64(rem_p,inv_p[0])+1)
   prod_max_r2_bits=max(prod_max_r2_bits,r20.bit_length())
   rad0=_isqrt_lut4(r20)
   l0=_chk_i64(ck0-ec_p[0]-rad0+Q32-1)>>32
   h0=_chk_i64(ck0+ec_p[0]+rad0)>>32
   if l0>h0:continue
   base_v=[sum(xp[j]*Bp[j][c]for j in range(1,12))for c in range(12)]
   rc=execute_chord_split(base_v,Bp[0],l0,h0,T12_chk=T12_prod,verify_ref=True,p_div=pp,s0_div=s0p,rho_div=rhop)
   if rc["lvs"]==0:continue
   p_chords+=1;p_lvs+=rc["lvs"];p_fast+=rc["fast"];p_bnd+=rc["bnd"]
   p_segs+=rc["segs"];p_bspl+=rc["bspl"];p_fspl+=rc["fspl"];p_cross+=rc["cross"];p_clips+=rc["clips"]
   prod_exec_mm+=rc["mm"];prod_exec_ae+=rc["ae"];prod_exec_div_err+=rc["div_err"]
  prod_exec_lvs+=p_lvs;prod_exec_clips+=p_clips
  ratios_12=[2.0**(log2(omega[12-k]/omega[12])+0.5*(log2(Dp[k])-k*log2(T12_prod)))for k in range(1,12)]
  n0_an=ratios_12[0];n1p_an=sum(ratios_12[1:])
  fp_ir=_fpw(ratios_12)
  sp_per_ch=float(WHEEL_OPS["chord_enter"]+WHEEL_OPS["wheel_m0_check"]+WHEEL_OPS["wheel_chord_rot"]+WHEEL_OPS["wheel_seg_disp0"])
  sp_ir=n0_an*sp_per_ch
  lf_ir=float(WHEEL_OPS["emit9_step1"])+0.015
  tot_ir=fp_ir+sp_ir+lf_ir+M.C_LAT_APPLY
  prod_anal_ir_vals.append(round(tot_ir,2))
 wsum=sum(1.0/pb for pb in pbits_list)
 prod_inv_bit_mean_ir=round(sum((1.0/pb)*ir for pb,ir in zip(pbits_list,prod_anal_ir_vals))/wsum,2)
 prod_max_ir=max(prod_anal_ir_vals)
 comp_max_rel_err_pct=round(100.0*max(abs(o-p)/min(o,p)for o,p in zip(comp_node_ratios_obs,comp_node_ratios_pred)),2)
 comp_fp_ir_max_rel_err_pct=round(100.0*max(abs(o-p)/min(o,p)for o,p in zip(comp_fp_ir_obs,comp_fp_ir_pred)),2)
 N_STRATA=32;K_DRAW=8;N_ORB48=N_STRATA*K_DRAW
 orb48_ir=[];orb48_lvs=0;strata_ht=[];tot_rej_tries=0
 _hb=HT_SEED_BASE
 p_lat_f=M.OP["p_lat"];lx2_f=M.OP["lx2"];_sp=math.sumprod
 for s_idx in range(N_STRATA):
  lo_b=p_lat_f+(lx2_f-p_lat_f)*s_idx/N_STRATA
  hi_b=p_lat_f+(lx2_f-p_lat_f)*(s_idx+1)/N_STRATA
  L_i=int(2.0**lo_b)|1;H_i=int(2.0**hi_b)
  if(H_i&1)==0:H_i-=1
  M_s=(H_i-L_i)//2+1;S_T=0;s_draws=[]
  for _kd in range(K_DRAW):
   rng_orb=random.Random((_hb<<20)|(s_idx<<8)|_kd)
   tries=0
   while True:
    tries+=1
    pp4=L_i+2*M.rand_below(rng_orb.getrandbits,M_s)
    if(pp4&7)not in(1,3)or not M._isp(pp4):continue
    xp4=M._ppow_h([0,1,0,0,0,0],pp4,pp4);xp4[1]=(xp4[1]-1)%pp4
    gp4=M.H_COEFFS if not any(xp4)else M._pgcd(xp4,M.H_COEFFS,pp4)
    if len(gp4)<=1:continue
    rs4=M._roots_mod(gp4,pp4)
    if not rs4:continue
    vp4=[];u_mod4=M.U_BN254%pp4
    for sr4 in rs4:
     ws4=_eval_poly(M.W_COEFFS,sr4,pp4)
     gl4=(-(((M.U_BN254+1)%pp4)*ws4+2))%pp4
     if gl4!=0:
      rh_pos=_sqrt_mod(pow(gl4,-1,pp4),pp4)
      if rh_pos>0:vp4.append((sr4,rh_pos))
    if vp4:
     m_p4=len(vp4)//3
     s04,rh_pos=rng_orb.choice(vp4)
     rh4=rh_pos if rng_orb.randrange(2)==0 else pp4-rh_pos
     break
   S_T+=tries;tot_rej_tries+=tries
   i34=pow(3,-1,pp4);sp4=[1]
   for _ in range(5):sp4.append((sp4[-1]*s04)%pp4)
   es4=[(sum(M.ZK_NUM[k][j]*sp4[j]for j in range(6))*i34)%pp4 for k in range(6)]
   lf4=es4+[((-rh4)*e)%pp4 for e in es4];i04=pow(lf4[0],-1,pp4)
   B4=[[pp4 if c==0 else 0 for c in range(12)]]
   for r1 in range(1,12):
    c0=(-lf4[r1]*i04)%pp4
    if c0>pp4//2:c0-=pp4
    B4.append([c0 if c==0 else(1 if c==r1 else 0)for c in range(12)])
   _,D4,lm4=DESC.bareiss_lll_12d(B4,G_sym,1764,1000,False,9,10,(99,100),True)
   r12_4=[2.0**(log2(omega[12-k]/omega[12])+0.5*(log2(D4[k])-k*log2(T12_prod)))for k in range(1,12)]
   n0_4,n1p_4=r12_4[0],sum(r12_4[1:])
   mu4=[[(lm4[j][k]<<32)//D4[k+1]if j>k else 0 for k in range(12)]for j in range(12)]
   b24=[min(M64,(D4[k+1]<<64)//(D4[k]*T12_prod))for k in range(12)]
   _mb=max(log2(D4[k+1])-log2(D4[k])-log2(T12_prod)for k in range(12));_PB['n']+=1;_PB['fail']+=(_mb>=0);_PB['maxb']=max(_PB['maxb'],_mb)
   in4=[((D4[k]*T12_prod<<38)+D4[k+1]-1)//D4[k+1]for k in range(12)]
   Xm4=[0]*12
   for k in range(11,-1,-1):
    rk4=math.isqrt(((D4[k]*T12_prod)+D4[k+1]-1)//D4[k+1])+1
    Xm4[k]=rk4+sum((Xm4[j]+1)//2 for j in range(k+1,12))
   ec4=[sum(Xm4[j]for j in range(k+1,12))+1 for k in range(12)]
   Bn4=[D4[k+1]/D4[k]for k in range(12)];R_f4=sqrt(float(T12_prod))
   inv_sB4=[1.0/sqrt(Bn4[k])for k in range(12)]
   mu4_c=[[mu4[j][k]for j in range(k+1,12)]for k in range(12)]
   B4_c=[[B4[j][c]for j in range(1,12)]for c in range(12)]
   c4_ch=c4_lv=0;c4_sp_w=c4_lf_w=c4_ap_w=0.0
   _gs=rng_orb.gauss;_rd=rng_orb.random;_iq32=1.0/float(Q32);_b0=B4[0];_in0=in4[0];_ec0=ec4[0];_mc0=mu4_c[0]
   _WO=WHEEL_OPS;_FL=0.015;_CH_ENT=float(_WO["chord_enter"]);_EM9=float(_WO["emit9_step1"]);_ELAST=float(_WO["emit9_last"])
   _M0C=float(_WO["wheel_m0_check"]);_M0D=float(_WO["wheel_m0_disp"]);_CROT=float(_WO["wheel_chord_rot"])
   _SDISP0=float(_WO["wheel_seg_disp0"]);_SK1=float(_WO["wheel_skip1_stub"]);_SK2=float(_WO["wheel_skip2_stub"])
   while c4_lv<160:
    gv=[_gs(0.0,1.0)for _ in range(11)]
    rg=R_f4*(_rd()**(1.0/11.0))/sqrt(_sp(gv,gv))
    xp4=[0]*12;rm4=Q38;ok4=True
    for k in range(11,0,-1):
     ck=-_sp(mu4_c[k],xp4[k+1:])if k<11 else 0
     xk=int(round(ck*_iq32+gv[k-1]*rg*inv_sB4[k]));xp4[k]=xk
     ad=abs((xk<<32)-ck);ec_k=ec4[k]
     if ad>ec_k:
      d19=(ad-ec_k)>>13;sub=((d19*d19)*b24[k])>>64
      if sub>rm4:ok4=False;break
      rm4-=sub
    if not ok4:continue
    ck0=-_sp(_mc0,xp4[1:])
    r20=((rm4*_in0)>>64)+1
    rad0=math.isqrt(r20<<52)+1
    l0=(ck0-_ec0-rad0+Q32-1)>>32;h0=(ck0+_ec0+rad0)>>32
    if l0>h0:continue
    x1p=xp4[1:]
    bv4=[_sp(x1p,B4_c[c])for c in range(12)]
    bmask=(1 if _bad_f3(bv4)else 0)|(2 if _bad_f3([bv4[c]+_b0[c]for c in range(12)])else 0)|(4 if _bad_f3([bv4[c]+2*_b0[c]for c in range(12)])else 0)
    lv=h0-l0+1
    c4_ch+=1;c4_lv+=lv
    if bmask==7:
     c4_sp_w+=_M0C+_CH_ENT+_CROT+_SDISP0
    elif bmask==0:
     chunks=(lv+127)//128
     c4_sp_w+=_M0C+_CH_ENT+_M0D+(chunks-1)*(_M0D+20.0)
     c4_lf_w+=(lv-1)*(_EM9+_FL)+(_ELAST+_FL)
     c4_ap_w+=lv*M.C_LAT_APPLY
    else:
     c4_sp_w+=_M0C+_CH_ENT+_CROT
     cmask=((bmask|(bmask<<3))>>(l0%3))&7
     is_p2=(cmask in(3,5,6))
     full3,rem3=divmod(lv,3)
     kept_i=full3*(1 if is_p2 else 2)+sum(1 for t_i in range(rem3)if not((cmask>>t_i)&1))
     chunks=(lv+127)//128
     skip_c=0.0
     if kept_i>=1:
      if(cmask&1)==0:skip_c=0.0
      elif(cmask&2)==0:skip_c=_SK1
      else:skip_c=_SK2
     c4_lf_w+=_SDISP0+skip_c+(chunks-1)*(_SDISP0+skip_c+20.0)
     if kept_i>=1:
      c4_lf_w+=(kept_i-1)*(_EM9+_FL)+(_ELAST+_FL)
      c4_ap_w+=kept_i*M.C_LAT_APPLY
   orb48_lvs+=c4_lv
   ir4=_fpw(r12_4)+n0_4*(c4_sp_w/c4_ch)+(c4_lf_w+c4_ap_w)/c4_lv
   orb48_ir.append(round(ir4,2));s_draws.append((pp4,m_p4,ir4))
  strata_ht.append((M_s,S_T,s_draws))
 Y_ht=X_ht=0.0
 for M_s,S_T,s_draws in strata_ht:
  N_hat=M_s*(K_DRAW-1.0)/(S_T-1.0)
  Y_ht+=N_hat*sum(mp*ir/pv for pv,mp,ir in s_draws)/K_DRAW
  X_ht+=N_hat*sum(mp/pv for pv,mp,_ in s_draws)/K_DRAW
 orb48_wmean=Y_ht/X_ht;var_ht=0.0
 for M_s,S_T,s_draws in strata_ht:
  N_hat=M_s*(K_DRAW-1.0)/(S_T-1.0)
  us_ht=[(mp/pv)*(ir-orb48_wmean)for pv,mp,ir in s_draws]
  ub_ht=sum(us_ht)/K_DRAW
  s2_ht=sum((u-ub_ht)**2 for u in us_ht)/(K_DRAW-1)
  var_N=(N_hat*N_hat)/float(K_DRAW-2)
  var_ht+=(N_hat*N_hat+var_N)*(s2_ht/K_DRAW)+var_N*(ub_ht*ub_ht)
 orb48_se=sqrt(var_ht)/X_ht
 orb48_wmean_p3se=round(orb48_wmean+3.0*orb48_se,3)
 extrap_factor=1.0+max(0.03,comp_fp_ir_max_rel_err_pct/100.0)
 orb48_extrap_ub=round(extrap_factor*orb48_wmean_p3se,3)
 GA=T4["GA"]
 aug=[GA[i][:]+[1.0 if i==j else 0.0 for j in range(6)]for i in range(6)]
 for col in range(6):
  piv=max(range(col,6),key=lambda r:abs(aug[r][col]))
  aug[col],aug[piv]=aug[piv],aug[col]
  pv=aug[col][col];aug[col]=[xv/pv for xv in aug[col]]
  for r in range(6):
   if r!=col and aug[r][col]:
    fac=aug[r][col];aug[r]=[aug[r][c]-fac*aug[col][c]for c in range(12)]
 max_cu=[sqrt(aug[i][6+i])for i in range(6)]
 max_a_prod=[round(589.43*sqrt(1.0+M.SC_GUARD_DELTA)*c,1)for c in max_cu]
 max_b_prod=[round(589.43*sqrt(1.764*(1.0+M.SC_GUARD_DELTA))*c,1)for c in max_cu]
 swar_bounds_ok=(all(max_a_prod[i]<AMAX_MR[i]+1.0 for i in range(6))and all(max_b_prod[i]<BMAX_MR[i]+1.0 for i in range(6))and PROD_MB_MR<(1<<59))
 Na_t=1024;B_st=16;n_b0=16;n_ch_slot=6
 a_vecs_raw=[]
 for a3i in range(0,4):
  for a4i in range(-4,5):
   for a5i in range(-4,5):
    for a0i in range(-2,3):
     av_t=(a0i,0,0,a3i,a4i,a5i)
     if any(av_t):
      qa_t=sum(av_t[u]*G_sym[u][v]*av_t[v]for u in range(6)for v in range(6))
      a_vecs_raw.append((qa_t,av_t))
     if len(a_vecs_raw)>=Na_t:break
    if len(a_vecs_raw)>=Na_t:break
   if len(a_vecs_raw)>=Na_t:break
  if len(a_vecs_raw)>=Na_t:break
 a_vecs_raw.sort(key=lambda x:x[0])
 a_sorted=[av for _,av in a_vecs_raw]
 a_keys=[_pack_v(av,SH_A,BI_A_CAN)for av in a_sorted]
 def _atab_sim_idx(pa_v):
  a3=(pa_v>>50)&0x1FF
  a4o=((pa_v>>40)&0x3FF)-508
  a5o=((pa_v>>30)&0x3FF)-508
  a0o=((pa_v>>20)&0x3FF)-510
  if((pa_v>>10)&0x3FF)!=512 or(pa_v&0x3FF)!=512 or not(0<=a3<4 and 0<=a4o<9 and 0<=a5o<9 and 0<=a0o<5):
   return-1
  return(a3<<11)|(a4o<<7)|(a5o<<3)|a0o
 a_tab_arr=[-1]*8192
 for _ai_i,_pa_i in enumerate(a_keys,1):
  _idx_i=_atab_sim_idx(_pa_i)
  assert _idx_i>=0 and a_tab_arr[_idx_i]==-1
  a_tab_arr[_idx_i]=_ai_i
 def _pack_a16(av):
  w0=((av[0]&0xFFFF)|((av[1]&0xFFFF)<<16)|((av[2]&0xFFFF)<<32)|((av[3]&0xFFFF)<<48))
  w1=((av[4]&0xFFFF)|((av[5]&0xFFFF)<<16))
  return(_chk_u64(w0),_chk_u64(w1))
 def _sxth16(v):
  v&=0xFFFF
  return v-65536 if v>=32768 else v
 a_packed16=[None]+[_pack_a16(av)for av in a_sorted]
 def _barrett_kta_18op(w0,w1,c_k,p_val,bias_val,inv_p):
  a0=_sxth16(w0);a1=_sxth16(w0>>16);a2=_sxth16(w0>>32);a3=_sxth16(w0>>48)
  a4=_sxth16(w1);a5=_sxth16(w1>>16)
  s_acc=_chk_i64(a0*c_k[0])
  s_acc=_chk_i64(s_acc+_chk_i64(a1*c_k[1]))
  s_acc=_chk_i64(s_acc+_chk_i64(a2*c_k[2]))
  s_acc=_chk_i64(s_acc+_chk_i64(a3*c_k[3]))
  s_acc=_chk_i64(s_acc+_chk_i64(a4*c_k[4]))
  s_acc=_chk_i64(s_acc+_chk_i64(a5*c_k[5]))
  sp=_chk_u64(s_acc+bias_val)
  q_est=_umulh64(sp,inv_p)
  rem=_chk_u64(sp-_chk_u64(q_est*p_val))
  assert rem<2*p_val
  return rem-p_val if rem>=p_val else rem
 barrett_kta_errs=0;barrett_kta_tested=0;barrett_max_sp_bits=0
 for pb_b2 in(42,46,48):
  p_b2=(1<<pb_b2)+1
  while True:
   if M._isp(p_b2):
    s_b2=_find_root_h(p_b2)
    if s_b2 is not None:
     ws_b2=_eval_poly(M.W_COEFFS,s_b2,p_b2)
     uw_b2=(-(((M.U_BN254+1)%p_b2)*ws_b2+2))%p_b2
     if uw_b2!=0:
      rho_b2=_sqrt_mod(pow(uw_b2,-1,p_b2),p_b2)
      if rho_b2>0:break
   p_b2+=2
  inv3_b=pow(3,-1,p_b2);sp_b=[1]
  for _ in range(5):sp_b.append((sp_b[-1]*s_b2)%p_b2)
  e_b=[(sum(M.ZK_NUM[k][j]*sp_b[j]for j in range(6))*inv3_b)%p_b2 for k in range(6)]
  neg_inv_rho=(-pow(rho_b2,-1,p_b2))%p_b2
  c_b=[(ek*neg_inv_rho)%p_b2 for ek in e_b]
  bias_b=_chk_u64(3072*p_b2)
  inv_p_b=_chk_u64((1<<64)//p_b2)
  barrett_max_sp_bits=max(barrett_max_sp_bits,(2*bias_b).bit_length())
  for ai in range(Na_t,0,-1):
   w0,w1=a_packed16[ai]
   got_k=_barrett_kta_18op(w0,w1,c_b,p_b2,bias_b,inv_p_b)
   av=a_sorted[ai-1]
   ref_k=(sum(av[k]*e_b[k]for k in range(6))*neg_inv_rho)%p_b2
   if got_k!=ref_k:barrett_kta_errs+=1
   barrett_kta_tested+=1
 p_l=127
 while True:
  if M._isp(p_l):
   s_l=_find_root_h(p_l)
   if s_l is not None:
    ws_l=_eval_poly(M.W_COEFFS,s_l,p_l)
    uw_l=(-(((M.U_BN254+1)%p_l)*ws_l+2))%p_l
    if uw_l!=0:
     rho_l=_sqrt_mod(pow(uw_l,-1,p_l),p_l)
     if rho_l>0:break
  p_l+=2
 inv3_l=pow(3,-1,p_l);sp_l=[1]
 for _ in range(5):sp_l.append((sp_l[-1]*s_l)%p_l)
 e_l=[(sum(M.ZK_NUM[k][j]*sp_l[j]for j in range(6))*inv3_l)%p_l for k in range(6)]
 neg_ir_l=(-pow(rho_l,-1,p_l))%p_l
 c_l=[(ek*neg_ir_l)%p_l for ek in e_l]
 bias_l=_chk_u64(3072*p_l);inv_p_l=_chk_u64((1<<64)//p_l)
 W_l=16;nb_l=(p_l+W_l-1)//W_l
 kta_l=[0]*(Na_t+1)
 head_l=[0]*nb_l;nxt_l=[0]*(Na_t+1);M47=(1<<47)-1
 for ai in range(Na_t,0,-1):
  w0,w1=a_packed16[ai]
  k=_barrett_kta_18op(w0,w1,c_l,p_l,bias_l,inv_p_l)
  kta_l[ai]=k;sb=k//W_l
  nxt_l[ai]=head_l[sb];head_l[sb]=((k&(W_l-1))<<47)|ai
 pref_errs=0
 for sb in range(nb_l):
  namax=rng.randint(40,Na_t)
  got_s=set();rec=head_l[sb]
  while rec:
   ai=rec&M47
   if ai>namax:break
   got_s.add((rec>>47,ai));rec=nxt_l[ai]
  ref_s={(kta_l[ai]&(W_l-1),ai)for ai in range(1,namax+1)if kta_l[ai]//W_l==sb}
  if got_s!=ref_s:pref_errs+=1
 off_x4=[0]*(nb_l+1)
 for ai in range(1,Na_t+1):off_x4[kta_l[ai]//W_l+1]+=1
 for sb in range(nb_l):off_x4[sb+1]+=off_x4[sb]
 flat_x4=[0]*Na_t;end_x4=off_x4[-1]
 for ai in range(Na_t,0,-1):
  k=kta_l[ai];sb=k//W_l
  off_x4[sb+1]-=1;flat_x4[off_x4[sb+1]]=(ai<<17)|(k&(W_l-1))
 off_x4=off_x4[1:]+[end_x4]
 x4_errs=0
 ROW_STRIDE_SIM=1<<11 # power-of-2 row stride (2048 halfwords = 4096 B >= 2*(Na_t+1))
 for sb in range(nb_l):
  namax=rng.randint(40,Na_t)
  namax_sh17=(namax<<17)|15
  for sgn_v in(1,-1):
   got_sc=[0]*(W_l*ROW_STRIDE_SIM);ref_sc=[0]*(W_l*ROW_STRIDE_SIM)
   st_i,en_i=off_x4[sb],off_x4[sb+1];p_i=st_i
   base_row=0 if sgn_v>0 else (W_l-1)*ROW_STRIDE_SIM
   while p_i+4<=en_i:
    r0,r1,r2,r3=flat_x4[p_i],flat_x4[p_i+1],flat_x4[p_i+2],flat_x4[p_i+3]
    if r3>namax_sh17:break
    for rv in(r0,r1,r2,r3):
     jr=rv&0xF
     r_off=(base_row+(jr<<11))if sgn_v>0 else(base_row-(jr<<11))
     _staddh_u16(got_sc,r_off+((rv>>16)>>1),720)
    p_i+=4
   while p_i<en_i:
    rv=flat_x4[p_i];p_i+=1
    if rv>namax_sh17:break
    jr=rv&0xF
    r_off=(base_row+(jr<<11))if sgn_v>0 else(base_row-(jr<<11))
    _staddh_u16(got_sc,r_off+((rv>>16)>>1),720)
   for ai in range(1,namax+1):
    if kta_l[ai]//W_l==sb:
     jr=kta_l[ai]&(W_l-1)
     if sgn_v<0:jr=(W_l-1)-jr
     ref_sc[jr*ROW_STRIDE_SIM+ai]=(ref_sc[jr*ROW_STRIDE_SIM+ai]+720)&0xFFFF
   if got_sc!=ref_sc:x4_errs+=1
 N_SC_SIM=B_st;N_LIST_NODES=4
 fabric_1hop_addr_errs=0
 for sid_t in(0,1,(1<<34)+7,(1<<35)-1):
  if _chk_u64(sid_t<<15)>=(1<<50):fabric_1hop_addr_errs+=1
  if _chk_u64(sid_t<<14)>=(1<<49):fabric_1hop_addr_errs+=1
 for rb_t in(0,1,32768,(1<<28)-1):
  bt_e=_chk_u64((((1<<35)-1)<<29)|rb_t)
  if(bt_e>>29)!=((1<<35)-1)or(((bt_e<<35)&M64)!=_chk_u64(rb_t<<35)):fabric_1hop_addr_errs+=1
  if _chk_u64((((bt_e<<35)&M64)+((1<<49)-2))&M64)!=(rb_t<<35)+((1<<49)-2):fabric_1hop_addr_errs+=1
  if _chk_u64((rb_t<<35)+((1<<49)-2))>=(1<<64):fabric_1hop_addr_errs+=1
 b_tab_packing_errs=0
 N_ROWS_SIM=N_SC_SIM*n_b0*8
 b_tab_sim=[0]*(N_ROWS_SIM+1) # flat array indexed by b_id; index N_ROWS_SIM = out-of-E_b over-radius dummy sink -> C[0] row 0!
 c_row_cnt=[1]+[0]*(N_SC_SIM-1) # C[0] reserves row 0 for out-of-E_b dummy sink; real C[0] rows start at row 1!
 for b_id in range(N_ROWS_SIM):
  if(b_id&7)==0:
   best_s=min(range(N_SC_SIM),key=lambda s:c_row_cnt[s])
   cur_s=best_s
  r_idx=c_row_cnt[cur_s];c_row_cnt[cur_s]=r_idx+1
  b_tab_sim[b_id]=_chk_u64((cur_s<<29)|r_idx)
 b_tab_sim[N_ROWS_SIM]=0 # out-of-E_b dummy sink: s_id=0, row_idx=0 (strictly disjoint from all real C[0] rows >= 1!)
 if any(b_tab_sim[b_id]==0 for b_id in range(N_ROWS_SIM)):b_tab_packing_errs+=1
 c_live_scores=[[[0]*(Na_t+1)for _ in range(c_row_cnt[s])]for s in range(N_SC_SIM)]
 score_ref=[[[0]*(Na_t+1)for _ in range(c_row_cnt[s])]for s in range(N_SC_SIM)]
 PAGE_CAP=2046;PTR_AI_MASK=(1<<35)-1
 stage_off=[16]*N_SC_SIM
 link_E_to_C=[bytearray(16384)for _ in range(N_SC_SIM)]
 flushed_full_pages=0;flushed_tail_pages=0;streamed_pages=0
 _pk_u64=struct.Struct("<Q").pack_into;_uk_u64=struct.Struct("<Q").unpack_from
 _pk_rec=struct.Struct("<QQ").pack_into;_uk_rec=struct.Struct("<QQ").unpack_from
 def _flush_1hop_to_C(sid_v,w_lp_v,cnt_v):
  nonlocal streamed_pages
  lbuf=link_E_to_C[sid_v]
  _pk_rec(lbuf,0,cnt_v|(w_lp_v<<16),0)
  hdr0,_=_uk_rec(lbuf,0)
  c_cnt=hdr0&0xFFFF;c_wlp=(hdr0>>16)&0xFFFF
  assert 1<=c_cnt<=PAGE_CAP
  streamed_pages+=1;off_b=16;c_sc=c_live_scores[sid_v]
  for _ in range(c_cnt):
   ptr_rec=_uk_u64(lbuf,off_b)[0];off_b+=8
   _staddh_u16(c_sc[ptr_rec>>35],(ptr_rec&PTR_AI_MASK)>>1,c_wlp)
 _rb=M.rand_below;_gb=rng.getrandbits
 for b_idx in range(224):
  w_lp_ticks=752+b_idx
  n_recs_b=37000 if b_idx in(0,223)else 18
  for rec_i in range(n_recs_b):
   ai=1+_rb(_gb,Na_t)
   b_id=N_ROWS_SIM if(rec_i%500==0)else _rb(_gb,N_ROWS_SIM)
   pa=a_keys[ai-1]
   br_e=b_tab_sim[b_id]
   sid_e=br_e>>29;r_idx_e=br_e&((1<<29)-1)
   if b_id<N_ROWS_SIM:
    score_ref[sid_e][r_idx_e][ai]=(score_ref[sid_e][r_idx_e][ai]+w_lp_ticks)&0xFFFF
   ai_e=a_tab_arr[_atab_sim_idx(pa)]
   ptr_rec=_chk_u64(((br_e<<35)&M64)+(ai_e<<1))
   if ptr_rec>=(1<<64):fabric_1hop_addr_errs+=1
   _pk_u64(link_E_to_C[sid_e],stage_off[sid_e],ptr_rec)
   stage_off[sid_e]+=8
   if stage_off[sid_e]&16384:
    if stage_off[sid_e]!=16384:fabric_1hop_addr_errs+=1
    _flush_1hop_to_C(sid_e,w_lp_ticks,PAGE_CAP)
    stage_off[sid_e]=16;flushed_full_pages+=1
  for sid_e in range(N_SC_SIM):
   cnt=(stage_off[sid_e]-16)>>3
   if cnt:
    _flush_1hop_to_C(sid_e,w_lp_ticks,cnt)
    stage_off[sid_e]=16;flushed_tail_pages+=1
 link_E_to_C=None;a_tab_arr=None;b_tab_sim=None
 supertile_errs=0
 for s_i in range(N_SC_SIM):
  r_start=1 if s_i==0 else 0
  for r_i in range(r_start,c_row_cnt[s_i]):
   if c_live_scores[s_i][r_i]!=score_ref[s_i][r_i]:supertile_errs+=1
 zero_frag_exact_ok=(supertile_errs==0 and fabric_1hop_addr_errs==0 and streamed_pages==flushed_full_pages+flushed_tail_pages)
 n_packed_nodes_sim=N_SC_SIM
 list_nodes=[bytearray(16384)for _ in range(N_LIST_NODES)]
 link_L_to_C=[bytearray(16384)for _ in range(N_LIST_NODES)]
 list_node_cnts=[0]*N_LIST_NODES
 for ai_l in range(1,Na_t+1):
  ln_k=(kta_l[ai_l]//W_l)%N_LIST_NODES;cnt_k=list_node_cnts[ln_k]
  _pk_u64(list_nodes[ln_k],cnt_k*8,(ai_l<<17)|(kta_l[ai_l]&(W_l-1)))
  list_node_cnts[ln_k]=cnt_k+1
 list_1hop_errs=0
 got_l_sc=[0]*(W_l*ROW_STRIDE_SIM);ref_l_sc=[0]*(W_l*ROW_STRIDE_SIM)
 for ai_l in range(1,Na_t+1):
  idx_l=(kta_l[ai_l]&(W_l-1))*ROW_STRIDE_SIM+ai_l
  ref_l_sc[idx_l]=(ref_l_sc[idx_l]+720)&0xFFFF
 for ln_k in range(N_LIST_NODES):
  link_L_to_C[ln_k][:]=list_nodes[ln_k]
  for j_l in range(list_node_cnts[ln_k]):
   rv_l=_uk_u64(link_L_to_C[ln_k],j_l*8)[0]
   _staddh_u16(got_l_sc,((rv_l&0xF)<<11)+((rv_l>>16)>>1),720)
 if got_l_sc!=ref_l_sc:list_1hop_errs+=1
 FRAG_TILES=256;frag_stress_nodes=N_SC_SIM;frag_stress_errs=0
 _g_loads=[0]*N_SC_SIM;_g_items=[rng.randint(1,64)for _ in range(FRAG_TILES)]
 for _it in _g_items:
  _bi=min(range(N_SC_SIM),key=lambda s:_g_loads[s]);_g_loads[_bi]+=_it
 if max(_g_loads)>sum(_g_items)/float(N_SC_SIM)+max(_g_items)+1e-9:frag_stress_errs+=1
 GB_t=T4["GB"];g00_t=GB_t[0][0]
 S5_t=[[GB_t[r][c]-GB_t[r][0]*GB_t[0][c]/g00_t for c in range(1,6)]for r in range(1,6)]
 delta5_t=sqrt(max(sum((0.5 if(m>>i)&1 else-0.5)*S5_t[i][j]*(0.5 if(m>>j)&1 else-0.5)for i in range(5)for j in range(5))for m in range(32)))
 V5_u=(8.0*math.pi*math.pi)/15.0/sqrt(M._det_mat(S5_t))
 aug5=[S5_t[i][:]+[1.0 if i==j else 0.0 for j in range(5)]for i in range(5)]
 for c5 in range(5):
  pv5=max(range(c5,5),key=lambda r:abs(aug5[r][c5]))
  aug5[c5],aug5[pv5]=aug5[pv5],aug5[c5];dv5=aug5[c5][c5]
  aug5[c5]=[x/dv5 for x in aug5[c5]]
  for r5 in range(5):
   if r5!=c5 and aug5[r5][c5]:
    fc5=aug5[r5][c5];aug5[r5]=[aug5[r5][j]-fc5*aug5[c5][j]for j in range(10)]
 inv5_d=[sqrt(aug5[i][5+i])for i in range(5)]
 discrete_chord_bound_test_errs=0
 for R_t in(4.0,8.0,12.0):
  b5=[int(floor(R_t*inv5_d[i]))for i in range(5)];R2_t=R_t*R_t;cnt5=[0]
  def _e5(ix,v5):
   if ix==5:
    if sum(v5[i]*S5_t[i][j]*v5[j]for i in range(5)for j in range(5))<=R2_t+1e-12:cnt5[0]+=1
    return
   for x5 in range(-b5[ix],b5[ix]+1):v5[ix]=x5;_e5(ix+1,v5)
  _e5(0,[0]*5)
  if cnt5[0]>V5_u*((R_t+delta5_t)**5):discrete_chord_bound_test_errs+=1
 def _pmul_h_z(a,b):
  t=[0]*11
  for i in range(6):
   for j in range(6):t[i+j]+=a[i]*b[j]
  for k in range(10,5,-1):
   c=t[k]
   if c:
    for j in range(6):t[k-6+j]-=c*M.H_COEFFS[j]
  return t[:6]
 def _det6(M6):
  A=[r[:]for r in M6];sg=1;pv=1
  for k in range(5):
   if A[k][k]==0:
    p0=next((r for r in range(k+1,6)if A[r][k]!=0),-1)
    if p0<0:return 0
    A[k],A[p0]=A[p0],A[k];sg= -sg
   p1=A[k][k]
   for i in range(k+1,6):
    for j in range(k+1,6):A[i][j]=(A[i][j]*p1-A[i][k]*A[k][j])//pv
   pv=p1
  return sg*A[5][5]
 def _exact_ng(a,b,u_val=M.U_BN254):
  a3=[sum(a[k]*M.ZK_NUM[k][j]for k in range(6))for j in range(6)]
  b3=[sum(b[k]*M.ZK_NUM[k][j]for k in range(6))for j in range(6)]
  a32=_pmul_h_z(a3,a3)
  g9=[(u_val+1)*x+2*z+y for x,z,y in zip(_pmul_h_z(M.W_COEFFS,a32),a32,_pmul_h_z(b3,b3))]
  cols=[];cur=g9[:]
  for j in range(6):
   cols.append(cur)
   if j<5:
    nxt=[0]+cur[:5];c=cur[5]
    if c:
     for m in range(6):nxt[m]-=c*M.H_COEFFS[m]
    cur=nxt
  dv=abs(_det6([[cols[c][r]for c in range(6)]for r in range(6)]))
  assert dv%(9**6)==0
  return dv//(9**6)
 def _bernstein_strip(P0,xs):
  pt=[list(xs)]
  while len(pt[-1])>1:
   pv=pt[-1]
   pt.append([pv[i]*pv[i+1]if i+1<len(pv)else pv[i]for i in range(0,len(pv),2)])
  rem=[P0%pt[-1][0]]
  for lvl in range(len(pt)-2,-1,-1):
   rem=[rem[i>>1]%v for i,v in enumerate(pt[lvl])]
  out=[]
  for x_v,r_v in zip(xs,rem):
   cur=x_v;g=math.gcd(cur,r_v)
   while g>1:
    cur//=g;g=math.gcd(cur,g)
   out.append((cur,x_v//cur))
  return out
 gen_norms=[]
 while len(gen_norms)<64:
  av=[rng.randint(-220,220)for _ in range(6)]
  bv=[rng.randint(-300,300)for _ in range(6)]
  if any(av)and any(bv)and M.is_ideal_primitive_oh(av,bv):
   gen_norms.append(_exact_ng(av,bv))
 ps14=M._primes_upto(1<<14);P14=1
 for p_v in ps14:P14*=p_v
 srt_out=[]
 for off in range(0,64,16):srt_out.extend(_bernstein_strip(P14,gen_norms[off:off+16]))
 srt_errs=0
 for idx_n,(cof_s,sm_s)in enumerate(srt_out):
  cur=gen_norms[idx_n]
  for p_v in ps14:
   while cur%p_v==0:cur//=p_v
  if cur!=cof_s or(gen_norms[idx_n]//cur)!=sm_s:srt_errs+=1
 def _is_prob_prime(n):
  if n<2:return False
  for p in(2,3,5,7,11,13,17,19,23,29,31,37):
   if n==p:return True
   if n%p==0:return False
  d=n-1;s=0
  while(d&1)==0:d>>=1;s+=1
  for a in(2,3,5,7,11,13,17,19,23,29,31,37):
   if a>=n:continue
   x=pow(a,d,n)
   if x in(1,n-1):continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 def _pollard_rho(n):
  if(n&1)==0:return 2
  for c in range(1,50):
   x=2;y=2;d=1
   while d==1:
    x=(x*x+c)%n;y=(y*y+c)%n;y=(y*y+c)%n
    d=math.gcd(abs(x-y),n)
   if d!=n:return d
  return n
 def _factor_full(n,small_ps):
  facs={}
  for p in small_ps:
   if p*p>n:break
   while n%p==0:facs[p]=facs.get(p,0)+1;n//=p
  stk=[n]if n>1 else[]
  while stk:
   m=stk.pop()
   if m==1:continue
   if _is_prob_prime(m):facs[m]=facs.get(m,0)+1
   else:d=_pollard_rho(m);stk.append(d);stk.append(m//d)
  return facs
 B_STRIP_E2E=1<<8;B_SIEVE_E2E=1<<20
 ps_e2e_small=M._primes_upto(1<<10)
 ps_e2e_strip=M._primes_upto(B_STRIP_E2E);P_e2e_strip=1
 for p_v in ps_e2e_strip:P_e2e_strip*=p_v
 e2e_norms=[]
 while len(e2e_norms)<320:
  av=[rng.randint(-4,4)for _ in range(6)]
  bv=[rng.randint(-6,6)for _ in range(6)]
  if any(av)and any(bv)and M.is_ideal_primitive_oh(av,bv):
   ng_e=_exact_ng(av,bv,u_val=5)
   if ng_e>=(1<<25):e2e_norms.append(ng_e)
 e2e_srt=[]
 for off in range(0,len(e2e_norms),32):e2e_srt.extend(_bernstein_strip(P_e2e_strip,e2e_norms[off:off+32]))
 e2e_smooth_sqf=0;e2e_smooth_p2=0;e2e_1lp=0;e2e_multi=0
 fn_0lp=0;fp_0lp=0;min_m_sm=1e9;min_g_1lp=1e9
 for ng_e,(c_post,_)in zip(e2e_norms,e2e_srt):
  nb=log2(ng_e);facs=_factor_full(c_post,ps_e2e_small)
  score_ticks=sum(floor(16.0*log2(p_v))for p_v in facs if B_STRIP_E2E<=p_v<=B_SIEVE_E2E)
  has_p2=any(e_v>=2 for p_v,e_v in facs.items()if B_STRIP_E2E<=p_v<=B_SIEVE_E2E)
  cap_nb=floor(ceil(nb)/8.0);rb=(cap_nb+1)/16.0
  res=(log2(c_post)if c_post>1 else 0.0)-score_ticks/16.0
  lps=[p_v for p_v in facs if p_v>B_SIEVE_E2E]
  passed=(res<=rb+0.25)
  if not lps:
   if not has_p2:
    e2e_smooth_sqf+=1
    mrg=(rb+0.25)-res
    if mrg<min_m_sm:min_m_sm=mrg
    if not passed:fn_0lp+=1
   else:
    e2e_smooth_p2+=1
  else:
   if len(lps)==1 and facs[lps[0]]==1:e2e_1lp+=1
   else:e2e_multi+=1
   gap=res-(rb+0.25)
   if gap<min_g_1lp:min_g_1lp=gap
   if passed:fp_0lp+=1
 ok=(smat_cyclic_ok and uncond_64bit_zero_exc_ok and isqrt_errs==0 and isqrt_max_over<=0.0012 and lut4_bucket_viol==0 and lut4_gt_lut7==0
  and staddh_errs==0 and x4_errs==0 and lat_mismatches==0
  and tie_errs==0 and barrett_kta_errs==0 and barrett_kta_tested>=3000 and barrett_max_sp_bits<64
  and pref_errs==0 and supertile_errs==0 and fabric_1hop_addr_errs==0 and b_tab_packing_errs==0 and discrete_chord_bound_test_errs==0 and list_1hop_errs==0 and zero_frag_exact_ok and frag_stress_errs==0 and flushed_full_pages>=20 
  and streamed_pages==flushed_full_pages+flushed_tail_pages and inc_c0_errs==0 and swar_bounds_ok
  and comp_tot_lvs>=10000 and comp_mm==0 and comp_ae==0 and comp_div_err==0 and comp_max_sa<=4 and wheel_tree_state_errs==0 and wheel_tree_leaves_checked>=10000 and swar21_cur_errs==0
  and comp_max_rel_err_pct<=5.0 and comp_fp_ir_max_rel_err_pct<=5.0
  and prod_exec_lvs>=2000 and prod_exec_mm==0 and prod_exec_ae==0 and prod_exec_div_err==0
  and prod_max_inv_bits<64 and prod_max_r2_bits<40 and prod_max_d2_bits<64 and prod_max_ck_bits<63
  and srt_errs==0
  and e2e_smooth_sqf>=10 and e2e_1lp>=100 and fn_0lp==0 and fp_0lp==0)
 return dict(
  ok=ok,
  smat_cyclic_a345_ok=smat_cyclic_ok,
  zero_float_fp_enum=True,
  strict_u64_i64_umulh_checked=True,
  zero_exception_64bit_q38_proved=uncond_64bit_zero_exc_ok,
  cert_min_gsym_exact=m6_exact,
  cert_min_gsym_pts_le_5gsc=m6_pts,
  uncond_min_bk_star_bits=round(log2(float(b_star_min_f)),3),
  uncond_inv_max_bits=uncond_inv_max.bit_length(),
  uncond_r2_q12_max_bits=uncond_r2_max.bit_length(),
  uncond_ck_max_bits=uncond_ck_max.bit_length(),
  uncond_d2_max_bits=uncond_d2_max.bit_length(),
  isqrt_u40_lut7_errs=isqrt_errs,lut4_buckets=len(_b4),lut4_bucket_violations=lut4_bucket_viol,lut4_buckets_above_lut7=lut4_gt_lut7,lut4_max_overshoot_pct=round(lut4_max_over*100.0,4),
  isqrt_u40_lut7_max_over_pct=round(isqrt_max_over*100.0,4),
  staddh_16b_semantic_errs=staddh_errs,
  b2_x4_prefix_staddh_errs=x4_errs,
  orbits_6ideal_checked=len(orbits),
  lattice_b2_mismatches=lat_mismatches,
  lattice_b2_raw_hits=tot_raw_hits,
  lattice_b2_jvor_hits=tot_jvor_hits,
  complete_root_to_leaf_orbits=len(pbits_list),
  complete_root_to_leaf_leaves=comp_tot_lvs,
  complete_root_to_leaf_nodes=comp_tot_nodes,
  complete_root_to_leaf_mismatches=comp_mm,
  complete_root_to_leaf_apply_errs=comp_ae,
  complete_root_to_leaf_div_errs=comp_div_err,
  complete_root_to_leaf_max_sa_cones=comp_max_sa,
  wheel_tree_leaves_checked=wheel_tree_leaves_checked,
  wheel_tree_state_errs=wheel_tree_state_errs,
  swar21_cur_errs=swar21_cur_errs,
  complete_scaled_nodes_per_leaf_obs=comp_node_ratios_obs,
  complete_scaled_nodes_per_leaf_pred=comp_node_ratios_pred,
  complete_scaled_max_rel_err_pct=comp_max_rel_err_pct,
  complete_scaled_fp_ir_obs=comp_fp_ir_obs,
  complete_scaled_fp_ir_pred=comp_fp_ir_pred,
  complete_scaled_fp_ir_max_rel_err_pct=comp_fp_ir_max_rel_err_pct,
  prod_exec_leaves=prod_exec_lvs,
  prod_exec_mismatches=prod_exec_mm,
  prod_exec_apply_errs=prod_exec_ae,
  prod_exec_div_errs=prod_exec_div_err,
  prod_exec_endpoint_clips=prod_exec_clips,
  prod_min_bk_star_bits=round(prod_min_b2_bits,2),
  prod_max_inv_b2_bits=prod_max_inv_bits,
  prod_max_r2_q12_bits=prod_max_r2_bits,
  prod_max_d2_bits=prod_max_d2_bits,
  prod_max_ck_bits=prod_max_ck_bits,
  prod_12lvl_analytic_ir_per_leaf=prod_anal_ir_vals,
  prod_inv_bit_weighted_mean_ir=prod_inv_bit_mean_ir,
  prod_max_ir=prod_max_ir,
  prod_256orb_strata=N_STRATA,
  prod_256orb_draws_per_stratum=K_DRAW,
  prod_256orb_count=N_ORB48,
  prod_256orb_leaves_sampled=orb48_lvs,
  prod_256orb_min_ir=min(orb48_ir),
  prod_256orb_max_ir=max(orb48_ir),
  prod_256orb_wmean_ht_ir=round(orb48_wmean,3),
  prod_256orb_se_strat_ir=round(orb48_se,3),
  prod_256orb_wmean_plus_3se_ir=orb48_wmean_p3se,
  prod_256orb_extrap_factor=round(extrap_factor,4),
  prod_256orb_extrap_upper_bound_ir=orb48_extrap_ub,
  prod_256orb_exact_rejection_ht=True,
  prod_256orb_uniform_random_root_choice=True,prod_256orb_run=_hb,prod_256orb_strata_ht=strata_ht,
  beta_ge1_census=dict(orbits_256=_PB['n'],ge1_256=_PB['fail'],max_log2_beta_256=round(_PB['maxb'],3),orbits_prod=_PB['n2'],ge1_prod=_PB['fail2'],max_log2_beta_prod=round(_PB['maxb2'],3)),prod_256orb_avg_rej_tries=round(tot_rej_tries/float(N_ORB48),1),
  prod_budget_c_lat_leaf=M.C_LAT_LEAF,
  swar_tie_cases_checked=len(tie_cases),
  swar_tie_canon_errs=tie_errs,
  barrett_kta_18op_tested=barrett_kta_tested,
  barrett_kta_18op_errs=barrett_kta_errs,
  barrett_max_sp_bits=barrett_max_sp_bits,
  desc_ai_prefix_chain_prime=p_l,
  desc_ai_prefix_chain_errs=pref_errs,
  score_node_replay_errs=supertile_errs,
  fabric_1hop_addr_errs=fabric_1hop_addr_errs,
  b_tab_packing_errs=b_tab_packing_errs,
  discrete_chord_bound_test_errs=discrete_chord_bound_test_errs,
  fabric_1hop_list_node_errs=list_1hop_errs,
  fabric_1hop_per_node_bytearray_zero_dict=True,
  fabric_1hop_zero_frag_exact_ok=zero_frag_exact_ok,
  score_nodes_sim=n_packed_nodes_sim,
  greedy_bound_subbands=FRAG_TILES,
  greedy_bound_nodes=frag_stress_nodes,
  greedy_node_bound_errs=frag_stress_errs,
  fabric_1hop_list_nodes_sim=N_LIST_NODES,
  logp_bands_1_16b=224,
  lat_8b_direct_ptr_page_cap=PAGE_CAP,
  emit_full_pages_flushed=flushed_full_pages,
  emit_tail_pages_flushed=flushed_tail_pages,
  emit_pages_streamed=streamed_pages,
  zero_hash_table_apply=True,
  inc_c0_exact_errs=inc_c0_errs,
  swar_bounds_ok=swar_bounds_ok,
  max_a_prod_coords=max_a_prod,
  max_b_prod_coords=max_b_prod,
  srt_p24_genuine_norms_checked=len(gen_norms),
  srt_batch_mismatches=srt_errs,
  e2e_factored_c6_norms_checked=len(e2e_norms),
  e2e_factored_smooth_sqfree=e2e_smooth_sqf,
  e2e_factored_smooth_p2=e2e_smooth_p2,
  e2e_factored_1lp=e2e_1lp,
  e2e_factored_multi_lp=e2e_multi,
  gate_0lp_fn=fn_0lp,
  gate_0lp_fp=fp_0lp,
  gate_0lp_min_smooth_margin_b=round(min_m_sm,4),
  gate_0lp_min_1lp_gap_b=round(min_g_1lp,2),
 )
