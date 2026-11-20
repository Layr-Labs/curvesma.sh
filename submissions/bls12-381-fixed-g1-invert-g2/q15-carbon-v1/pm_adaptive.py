"""Self-contained classic PM-Basis live verifier & Q(zeta_15) real-2-core-FB adaptive H5 verifier."""
import base64,cmath,json,lzma,math,random
from decimal import Context,Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
def _padd(a,b,p):
 n=max(len(a),len(b));o=[0]*n
 for i,v in enumerate(a):o[i]=v
 for i,v in enumerate(b):o[i]=(o[i]+v)%p
 while o and o[-1]==0:o.pop()
 return o
def _pmul(a,b,p,cap=None):
 if not a or not b:return []
 la,lb=len(a),len(b);n=la+lb-1 if cap is None else min(la+lb-1,cap);o=[0]*n
 for i,ai in enumerate(a):
  if ai:
   for j in range(min(lb,n-i)):o[i+j]=(o[i+j]+ai*b[j])%p
 while o and o[-1]==0:o.pop()
 return o
def _pmat_mul(A,B,p,cap=None):
 r,k,c=len(A),len(A[0]),len(B[0])
 C=[[[] for _ in range(c)] for _ in range(r)]
 for i in range(r):
  for t in range(k):
   ait=A[i][t]
   if not ait:continue
   bt=B[t]
   for j in range(c):
    if bt[j]:C[i][j]=_padd(C[i][j],_pmul(ait,bt[j],p,cap),p)
 return C
def _mbasis1(F0,s,p):
 m,t=len(F0),len(F0[0]);order=sorted(range(t),key=lambda j:(s[j],j))
 pivs={};Pcols=[[0]*t for _ in range(t)];dinc=[0]*t
 for cj in order:
  vec=[F0[r][cj]%p for r in range(m)];pc=[0]*t;pc[cj]=1
  for r in range(m):
   v=vec[r]
   if not v:continue
   if r in pivs:
    pv,pp=pivs[r]
    for rr in range(r,m):
     if pv[rr]:vec[rr]=(vec[rr]-v*pv[rr])%p
    for cc in range(t):
     if pp[cc]:pc[cc]=(pc[cc]-v*pp[cc])%p
   else:
    inv=pow(v,-1,p);pivs[r]=([(x*inv)%p for x in vec],[(x*inv)%p for x in pc]);dinc[cj]=1;break
  Pcols[cj]=pc
 P=[[[] for _ in range(t)] for _ in range(t)]
 for j in range(t):
  for i in range(t):
   v=Pcols[j][i]
   if v:P[i][j]=[0,v] if dinc[j] else [v]
 return P,[s[j]+dinc[j] for j in range(t)]
def pm_basis(F,sigma,s,p):
 """Giorgi-Jeannerod-Villard (2003) / Thome (2002) divide-and-conquer PM-Basis over F_p[z]."""
 m,t=len(F),len(F[0])
 if sigma==1:return _mbasis1([[F[i][j][0] if F[i][j] else 0 for j in range(t)] for i in range(m)],s,p)
 s1=sigma//2;s2=sigma-s1
 P1,sm=pm_basis([[F[i][j][:s1] for j in range(t)] for i in range(m)],s1,s,p)
 FP1=_pmat_mul(F,P1,p,cap=sigma)
 G=[[FP1[i][j][s1:sigma] for j in range(t)] for i in range(m)]
 P2,so=pm_basis(G,s2,sm,p)
 return _pmat_mul(P1,P2,p),so
def pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_ref):
 """Run live divide-and-conquer PM-Basis on S_k = X A^k Y (k>=1), recover exact ker(A) via Horner."""
 m=len(seq[0]);n=len(seq[0][0]);D=len(A_sq);L=hr+hc-1
 S=seq[1:1+L]
 M=[[[] for _ in range(n+m)] for _ in range(m)]
 for i in range(m):
  for j in range(n):
   poly=[S[k][i][j] for k in range(L)]
   while poly and poly[-1]==0:poly.pop()
   M[i][j]=poly
  M[i][n+i]=[r-1]
 P_b,s_out=pm_basis(M,L,[0]*n+[1]*m,r)
 def mvec(P,v):return [sum(P[i][j]*v[j] for j in range(len(v)))%r for i in range(len(P))]
 v_pm=None;deg_used=None
 for col in sorted(range(n+m),key=lambda j:s_out[j]):
  if s_out[col]>hc-1:continue
  q=[P_b[i][col] for i in range(n)]
  dg=max((len(poly)-1 for poly in q),default=-1)
  if dg<0:continue
  acc=[0]*D
  for t in range(dg+1):
   qt=[q[bi][t] if t<len(q[bi]) else 0 for bi in range(n)]
   yc=[sum(Y[ri][bi]*qt[bi] for bi in range(n))%r for ri in range(D)]
   acc=[(a+y)%r for a,y in zip(mvec(A_sq,acc),yc)]
  if any(acc) and all(x==0 for x in mvec(A_sq,acc)):
   v_pm=acc;deg_used=dg;break
 assert v_pm is not None,"PM-Basis failed to recover ker(A)!"
 idx0=next(i for i,x in enumerate(v_ref) if x)
 ratio=(v_pm[idx0]*pow(v_ref[idx0],-1,r))%r
 assert ratio!=0 and all((v_ref[i]*ratio)%r==v_pm[i] for i in range(D))
 return dict(pm_basis_order=L,pm_basis_cols=n+m,pm_basis_deg=deg_used,pm_basis_matches_true_ker=True)
def _is_prime_det(n):
 if n<2:return False
 for p in (2,3,5,7,11,13,17,19,23,29,31,37):
  if n==p:return True
  if n%p==0:return False
 d=n-1;s=0
 while d%2==0:d//=2;s+=1
 for a in (2,3,5,7,11,13,17,19,23,29,31,37):
  x=pow(a,d,n)
  if x==1 or x==n-1:continue
  for _ in range(s-1):
   x=(x*x)%n
   if x==n-1:break
  else:return False
 return True
def _fp_exact_int(Br,dg,lm,R2,collect=False):
 z=[0]*16;nodes=0;hits=0;out=[]
 def _rec(k,rem_sc,nz):
  nonlocal nodes,hits
  num_k=-sum(lm[j][k]*z[j] for j in range(k+1,16));dk1=dg[k+1];dk=dg[k]
  rhs=rem_sc*dk;rad=math.isqrt(rhs)
  lo=(num_k-rad+dk1-1)//dk1;hi=(num_k+rad)//dk1
  if not nz:lo=max(0,lo)
  for zk in range(lo,hi+1):
   z[k]=zk;nz2=nz or (zk>0);diff=zk*dk1-num_k;rnum=rhs-diff*diff
   if rnum<0:continue
   if k==0:
    if not nz2:continue
    hits+=1
    if collect:
     v=[sum(z[i]*Br[i][t] for i in range(16)) for t in range(16)]
     for x in v:
      if x!=0:
       if x<0:v=[-t for t in v]
       break
     out.append(tuple(v))
   else:
    assert rnum%dk1==0;nodes+=1;_rec(k-1,rnum//dk1,nz2)
 _rec(15,R2*dg[16],False)
 return hits,nodes,out
def _fp_q32_int(Br,dg,lm,R2,collect=False):
 S=1<<32;MU=[[0]*16 for _ in range(16)]
 for j in range(16):
  for k in range(j):
   dk1=dg[k+1];num=lm[j][k]*S
   MU[j][k]=(num+dk1//2)//dk1 if num>=0 else -((-num+dk1//2)//dk1)
 L2_LO=[(dg[k+1]*S)//dg[k] for k in range(16)]
 INV_HI=[math.isqrt((dg[k]*S*S+dg[k+1]-1)//dg[k+1])+1 for k in range(16)]
 z=[0]*16;nodes=0;hits=0;out=[]
 def _rec(k,rem_q,nz,asum):
  nonlocal nodes,hits
  cq=-sum(MU[j][k]*z[j] for j in range(k+1,16));ec=(asum+1)>>1
  sr=math.isqrt(rem_q)+1;rq=((sr*INV_HI[k])>>16)+1+ec
  lo=(cq-rq+S-1)>>32;hi=(cq+rq)>>32
  if not nz:lo=max(0,lo)
  l2k=L2_LO[k]
  for zk in range(lo,hi+1):
   z[k]=zk;nz2=nz or (zk>0);da=abs((zk<<32)-cq);dl=da-ec if da>ec else 0
   sq=(((dl*dl)>>32)*l2k)>>32
   if sq>rem_q:continue
   if k==0:
    if not nz2:continue
    v=[sum(z[i]*Br[i][t] for i in range(16)) for t in range(16)]
    if sum(x*x for x in v)<=R2:
     hits+=1
     if collect:
      for x in v:
       if x!=0:
        if x<0:v=[-t for t in v]
        break
      out.append(tuple(v))
   else:nodes+=1;_rec(k-1,rem_q-sq,nz2,asum+abs(zk))
 _rec(15,R2<<32,False,0)
 return hits,nodes,out
def verify_b16_dsieve_live(data,desc,dlive,m_node,U_bls,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 _Rm=160;_M=_Rm*_Rm
 def _ksq_b(poly,bpc):
  _x=int.from_bytes(b"".join(v.to_bytes(bpc,"little") for v in poly),"little");_p=_x*_x
  _nb=max((_M+1)*bpc,(_p.bit_length()+7)//8);_bs=_p.to_bytes(_nb,"little")
  return [int.from_bytes(_bs[i*bpc:(i+1)*bpc],"little") for i in range(_M+1)]
 _t1=[0]*(_M+1);_t1[0]=1
 for _k in range(1,_Rm+1):_t1[_k*_k]=2
 _t16=_ksq_b(_ksq_b(_ksq_b(_ksq_b(_t1,4),8),8),16)
 for _rk,_cv in desc["cells_rank2"].items():
  _R=int(_rk);assert int(_cv)==(sum(_t16[:_R*_R+1])-1)//2,f"B_16(R)/+-1 mismatch R={_R}"
 assert 4*max(int(desc["cells_rank2"][str(r["R"])]) for r in dlive["rows"])<(1<<111)
 _yp=[y_pow_z(k) for k in range(15)]
 _ug=[[1,0,0,1,0,0,0,0],[-1,0,1,0,0,0,0,0],[-1,1,0,0,0,0,0,0]]
 _ui=[[1,0,-1,0,0,0,1,-1],[0,-1,0,0,-1,1,-1,0],[0,-1,-1,0,-1,0,0,-1]]
 for _a,_b in zip(_ug,_ui):assert mul_phi15_z(_a,_b)==[1]+[0]*7
 _ym=[tuple(tuple(mul_phi15_z(y,[1 if i==j else 0 for i in range(8)])[r] for j in range(8)) for r in range(8)) for y in _yp]
 _pc=[{0:[1]+[0]*7} for _ in range(3)]
 def _upow(_idx,_n):
  _dd=_pc[_idx]
  if _n in _dd:return _dd[_n]
  _rv=mul_phi15_z(_upow(_idx,_n-1),_ug[_idx]) if _n>0 else mul_phi15_z(_upow(_idx,_n+1),_ui[_idx])
  _dd[_n]=_rv;return _rv
 _CL=Context(prec=40,rounding=ROUND_FLOOR);_CH=Context(prec=40,rounding=ROUND_CEILING)
 def _iv_add(a,b):return (_CL.add(a[0],b[0]),_CH.add(a[1],b[1]))
 def _iv_mul(a,b):
  a0,a1=a;b0,b1=b
  return (min(_CL.multiply(a0,b0),_CL.multiply(a0,b1),_CL.multiply(a1,b0),_CL.multiply(a1,b1)),
   max(_CH.multiply(a0,b0),_CH.multiply(a0,b1),_CH.multiply(a1,b0),_CH.multiply(a1,b1)))
 def _iv_ln(a):assert a[0]>0;return (_CL.ln(a[0]),_CH.ln(a[1]))
 def _iv_cs(m_red):
  with localcontext() as c:
   c.prec=45
   pi_m=Decimal("3.1415926535897932384626433832795028841971")
   th=Decimal(2*m_red)*pi_m/Decimal(15);cv=Decimal(1);sv=th;tc=Decimal(1);ts=th;t2=-th*th
   for k in range(1,28):tc=tc*t2/((2*k-1)*(2*k));cv+=tc;ts=ts*t2/((2*k)*(2*k+1));sv+=ts
   eps=Decimal("1e-35")
  clo,slo=_CL.subtract(cv,eps),_CL.subtract(sv,eps);chi,shi=_CH.add(cv,eps),_CH.add(sv,eps)
  return (clo,chi),(slo,shi)
 _zp_iv=[];_Z_TAB=[]
 _S40=1<<40
 for _kk in (1,7,4,2):
  _row_iv=[];_row_t=[]
  for _j in range(8):
   _mr=(_kk*_j)%15
   if _mr>7:_mr-=15
   _ci,_si=_iv_cs(_mr);_row_iv.append((_ci,_si))
   with localcontext() as c:
    c.prec=40
    _cn=int((_ci[0]*_S40+Decimal("0.5")).to_integral_value(rounding=ROUND_FLOOR))
    _sn=int((_si[0]*_S40+Decimal("0.5")).to_integral_value(rounding=ROUND_FLOOR))
    assert abs(_ci[0]*_S40-_cn)<=Decimal("0.51") and abs(_ci[1]*_S40-_cn)<=Decimal("0.51")
    assert abs(_si[0]*_S40-_sn)<=Decimal("0.51") and abs(_si[1]*_S40-_sn)<=Decimal("0.51")
   _row_t.append((_cn,_sn))
  _zp_iv.append(_row_iv);_Z_TAB.append(_row_t)
 def _iv_norm2(_vec,_k):
  _re=(Decimal(0),Decimal(0));_im=(Decimal(0),Decimal(0))
  for _j in range(8):
   _v=_vec[_j]
   if _v:
    _dv=(Decimal(_v),Decimal(_v))
    _re=_iv_add(_re,_iv_mul(_dv,_zp_iv[_k][_j][0]))
    _im=_iv_add(_im,_iv_mul(_dv,_zp_iv[_k][_j][1]))
  def _sq(iv):
   lo=Decimal(0) if iv[0]<=0<=iv[1] else min(_CL.multiply(iv[0],iv[0]),_CL.multiply(iv[1],iv[1]))
   return (lo,max(_CH.multiply(iv[0],iv[0]),_CH.multiply(iv[1],iv[1])))
  return _iv_add(_sq(_re),_sq(_im))
 _Lg_iv=[[_iv_mul((Decimal("0.5"),Decimal("0.5")),_iv_ln(_iv_norm2(u,k))) for k in range(4)] for u in _ug]
 _Mm_iv=[[_Lg_iv[i][k] for i in range(3)] for k in range(3)]
 with localcontext() as _c45:
  _c45.prec=45
  _Mm0=[[(_Mm_iv[r][c][0]+_Mm_iv[r][c][1])/Decimal(2) for c in range(3)] for r in range(3)]
  _a,_b,_c=_Mm0[0];_d,_e,_f=_Mm0[1];_g,_h,_ii=_Mm0[2]
  _dt=_a*(_e*_ii-_f*_h)-_b*(_d*_ii-_f*_g)+_c*(_d*_h-_e*_g)
  assert abs(Decimal(8)*abs(_dt)-Decimal("4.6618207772966835"))<Decimal("1e-12")
  _Mi0=[[(_e*_ii-_f*_h)/_dt,(_c*_h-_b*_ii)/_dt,(_b*_f-_c*_e)/_dt],
   [(_f*_g-_d*_ii)/_dt,(_a*_ii-_c*_g)/_dt,(_c*_d-_a*_f)/_dt],
   [(_d*_h-_e*_g)/_dt,(_b*_g-_a*_h)/_dt,(_a*_e-_b*_d)/_dt]]
  _mi_inf=max(sum(abs(_Mi0[r][c]) for c in range(3)) for r in range(3))
 _eps_E=Decimal(0)
 for _r in range(3):
  _rs=Decimal(0)
  for _cl in range(3):
   _s=(Decimal(0),Decimal(0))
   for _j in range(3):_s=_iv_add(_s,_iv_mul(_Mm_iv[_r][_j],(_Mi0[_j][_cl],_Mi0[_j][_cl])))
   _tgt=Decimal(1) if _r==_cl else Decimal(0)
   _rs=_CH.add(_rs,max(abs(_CH.subtract(_s[0],_tgt)),abs(_CH.subtract(_s[1],_tgt))))
  if _rs>_eps_E:_eps_E=_rs
 assert _eps_E<Decimal("1e-30")
 _mi_rad=_CH.divide(_CH.multiply(_mi_inf,_eps_E),_CL.subtract(Decimal(1),_eps_E))
 _Mi_iv=[]
 for _r in range(3):
  _rw=[]
  for _cl in range(3):
   _lo=_CL.subtract(_Mi0[_r][_cl],_mi_rad);_hi=_CH.add(_Mi0[_r][_cl],_mi_rad)
   _rw.append((_lo,_hi))
  _Mi_iv.append(_rw)
 _ep_lo=[{} for _ in range(3)];_SC=Decimal(10**9);_SC4=10**36
 def _epow_lo(_idx,_n):
  if _n not in _ep_lo[_idx]:
   _out=[]
   for _k in range(4):
    _iv=_Lg_iv[_idx][_k]
    _ex=_CL.multiply(Decimal(2*_n),_iv[0] if _n>=0 else _iv[1])
    _out.append(int(_CL.multiply(_CL.exp(_ex),_SC).to_integral_value(rounding=ROUND_FLOOR)))
   _ep_lo[_idx][_n]=_out
  return _ep_lo[_idx][_n]
 _d15001=Decimal("15.001")
 def _cnt_orb(_va,_vb,_R2):
  _eks=[_iv_mul((Decimal(2),Decimal(2)),_iv_add(_iv_norm2(_va,k),_iv_norm2(_vb,k))) for k in range(4)]
  _bd_hi=_CH.multiply(_d15001,Decimal(_R2));_bd_i=int(_CH.multiply(_bd_hi,Decimal(_SC4)).to_integral_value(rounding=ROUND_CEILING))
  _eks_i=[int(_CL.multiply(_ek[0],_SC).to_integral_value(rounding=ROUND_FLOOR)) for _ek in _eks]
  _Bk=[_CH.multiply(_CH.ln(_CH.divide(_bd_hi,_ek[0])),Decimal("0.5")) for _ek in _eks]
  _SB=Decimal(0)
  for _b in _Bk:_SB=_CH.add(_SB,_b)
  _lo=[_CL.minus(_CH.subtract(_SB,_Bk[k])) for k in range(3)]
  _nlo=[];_nhi=[]
  for _i in range(3):
   _iv_s=(Decimal(0),Decimal(0))
   for _k in range(3):_iv_s=_iv_add(_iv_s,_iv_mul(_Mi_iv[_i][_k],(_lo[_k],_Bk[_k])))
   _nlo.append(int(_iv_s[0].to_integral_value(rounding=ROUND_FLOOR)))
   _nhi.append(int(_iv_s[1].to_integral_value(rounding=ROUND_CEILING)))
  _mo=0
  for _n0 in range(_nlo[0],_nhi[0]+1):
   _e0=_epow_lo(0,_n0)
   for _n1 in range(_nlo[1],_nhi[1]+1):
    _e1=_epow_lo(1,_n1);_e01=[_eks_i[k]*_e0[k]*_e1[k] for k in range(4)]
    for _n2 in range(_nlo[2],_nhi[2]+1):
     _e2=_epow_lo(2,_n2)
     if _e01[0]*_e2[0]+_e01[1]*_e2[1]+_e01[2]*_e2[2]+_e01[3]*_e2[3]>_bd_i:continue
     _u=mul_phi15_z(mul_phi15_z(_upow(0,_n0),_upow(1,_n1)),_upow(2,_n2))
     a0,a1,a2,a3,a4,a5,a6,a7=mul_phi15_z(_u,_va);b0,b1,b2,b3,b4,b5,b6,b7=mul_phi15_z(_u,_vb)
     for _mt in _ym:
      _s=0
      for r0,r1,r2,r3,r4,r5,r6,r7 in _mt:
       _xa=r0*a0+r1*a1+r2*a2+r3*a3+r4*a4+r5*a5+r6*a6+r7*a7;_s+=_xa*_xa
       if _s>_R2:break
       _xb=r0*b0+r1*b1+r2*b2+r3*b3+r4*b4+r5*b5+r6*b6+r7*b7;_s+=_xb*_xb
       if _s>_R2:break
      else:_mo+=1
  assert _mo>=1;return _mo
 def _is_cp(_va,_vb):
  _g0=math.gcd(abs(norm_q15_z(_va)),abs(norm_q15_z(_vb)))
  if _g0<=1:return True
  _t=_g0;_d=2;_ps=[]
  while _d*_d<=_t and _d<=10000:
   if _t%_d==0:
    _ps.append(_d)
    while _t%_d==0:_t//=_d
   _d+=1 if _d==2 else 2
  if _t>1:_ps.append(_t)
  for _pp in _ps:
   _g1=pgcd_mod(_va,PHI15,_pp)
   if len(_g1)>1 and len(pgcd_mod(_g1,_vb,_pp))>1:return False
  return True
 _rng=random.Random(20261120)
 _p_ds=31;assert _is_prime_det(_p_ds);_c_ds=[1]+[_rng.randrange(_p_ds) for _ in range(15)]
 _B_ds=[[0]*16 for _ in range(16)];_B_ds[0][0]=_p_ds
 for _j in range(1,16):_B_ds[_j][0]=(-_c_ds[_j])%_p_ds;_B_ds[_j][_j]=1
 _Br_ds,_dg_ds,_lm_ds=lll_16d_int(_B_ds,9,10);assert _dg_ds[16]==_p_ds*_p_ds
 _fp_cnt,_fp_nd,_fp_h=_fp_exact_int(_Br_ds,_dg_ds,_lm_ds,4,collect=True)
 _q32_cnt,_q32_nd,_q32_h=_fp_q32_int(_Br_ds,_dg_ds,_lm_ds,4,collect=True)
 _bf_h=[];_vb_ds=[0]*16
 def _bf_rec(_idx,_rm,_nz):
  if _idx==16:
   if _nz and sum(_c_ds[j]*_vb_ds[j] for j in range(16))%_p_ds==0:_bf_h.append(tuple(_vb_ds))
   return
  _r=math.isqrt(_rm)
  for _x in range(-_r,_r+1):
   if not _nz and _x<0:continue
   _vb_ds[_idx]=_x;_bf_rec(_idx+1,_rm-_x*_x,_nz or (_x>0))
 _bf_rec(0,4,False)
 assert _fp_cnt==_q32_cnt==len(_fp_h)==len(_q32_h)==len(set(_q32_h))==len(_bf_h)>400 and set(_fp_h)==set(_q32_h)==set(_bf_h)
 _bc=[[0]*5 for _ in range(17)];_bc[0]=[1]*5
 for _d in range(1,17):
  for _s in range(5):_bc[_d][_s]=sum(_bc[_d-1][_s-_x*_x] for _x in range(-math.isqrt(_s),math.isqrt(_s)+1))
 _s16_2=(_bc[16][4]-1)//2;_off2=(_bc[16][4]+1)//2;_rks=set()
 for _hv in _q32_h:
  _rk=0;_rm=4
  for _d in range(16):
   for _x in range(-math.isqrt(_rm),_hv[_d]):_rk+=_bc[15-_d][_rm-_x*_x]
   _rm-=_hv[_d]*_hv[_d]
  _idx=_rk-_off2;assert 0<=_idx<_s16_2;_rks.add(_idx)
 assert len(_rks)==len(_q32_h)
 def _find_prime_Q(_bp,_sd_seed):
  _rp=random.Random(_sd_seed);_p=((1<<(_bp-1))|(_rp.getrandbits(_bp-2)<<1)|1);_p+=(1-_p)%15
  while True:
   if _is_prime_det(_p) and (_p%9 in (4,7)):
    for _g in range(2,50):
     _s=pow(_g,(_p-1)//15,_p)
     if pow(_s,3,_p)!=1 and pow(_s,5,_p)!=1:
      _rhs=(U_bls*pow(_s,10,_p))%_p
      if pow(_rhs,(_p-1)//3,_p)==1:
       _t=pow(_rhs,(2*_p+1)//9 if _p%9==4 else (_p+2)//9,_p)
       if pow(_t,3,_p)==_rhs:return _p,_s,_t
      break
   _p+=15
 def _struct_L(_a1,_b1,_a2,_b2,_p,_s,_t):
  _sp=[pow(_s,j,_p) for j in range(8)]
  _ev=lambda poly:sum(poly[j]*_sp[j] for j in range(8))%_p
  _u1=(_ev(_a1)-_t*_ev(_b1))%_p;_u2=(_ev(_a2)-_t*_ev(_b2))%_p
  _c=[(_u1*_sp[j])%_p for j in range(8)]+[(_u2*_sp[j])%_p for j in range(8)]
  _k0=next(i for i,v in enumerate(_c) if v!=0);_iv=pow(_c[_k0],-1,_p);_c=[(v*_iv)%_p for v in _c]
  _Bb=[[0]*16 for _ in range(16)];_Bb[0][_k0]=_p;_rr=1
  for _j in range(16):
   if _j==_k0:continue
   _Bb[_rr][_k0]=(-_c[_j])%_p;_Bb[_rr][_j]=1;_rr+=1
  return _Bb
 _alpha=100.0/74.0;_vb64=[1.0]+[(math.pi**(d/2.0)/math.gamma(d/2.0+1.0))*(64.05**d) for d in range(1,17)]
 def _b_int_64(_ls):
  _inv=[1.0/x for x in _ls];_dp=[1.0];_tot=0.0
  for _k in range(15,0,-1):
   _u=_inv[_k];_ndp=[0.0]*(len(_dp)+1)
   for _s,_v in enumerate(_dp):_ndp[_s]+=_v;_ndp[_s+1]+=_v*_u
   _dp=_ndp;_tot+=0.5*(sum(_vb64[_s]*_dp[_s] for _s in range(len(_dp)))+1.0)
  return _tot
 _ds_disc={};_gate_ratios=[]
 for _yi,_lbl in ((0,"g"),(12,"f")):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16)
  _Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_qq;_rp=[pow(_ry,i,_qq) for i in range(8)]
  for _i in range(1,8):_Bq[_i][0]=(-_rp[_i])%_qq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=(_rt*_rp[_i])%_qq;_Bq[8+_i][8+_i]=1
  _Rr,_,_=lll_16d_int(_Bq,99,100);_a1,_b1=_Rr[0][:8],_Rr[0][8:]
  _a2,_b2=next((v[:8],v[8:]) for v in _Rr[1:] if any(sub_phi15_z(mul_phi15_z(_a1,v[8:]),mul_phi15_z(_b1,v[:8]))))
  _p83,_s83,_t83=_find_prime_Q(84,1084+_yi*10)
  _Br83,_dg83,_lm83=lll_16d_int(_struct_L(_a1,_b1,_a2,_b2,_p83,_s83,_t83),9,10);assert _dg83[16]==_p83*_p83
  _hq83,_ndq83,_=_fp_q32_int(_Br83,_dg83,_lm83,64*64,collect=False)
  _ls83=[math.sqrt(_dg83[i+1]/_dg83[i]) for i in range(16)];_bi83=_b_int_64(_ls83)
  assert _hq83>400 and _ndq83<=_bi83
  if _lbl=="g":
   _h83,_nd83,_=_fp_exact_int(_Br83,_dg83,_lm83,64*64,collect=False);assert _h83==_hq83 and _nd83==_ndq83
  _ds_disc[f"{_lbl}_p84_R64"]=dict(hits=_hq83,nodes_q32=_ndq83,B_int_ub=round(_bi83,1),actual_over_ub=round(_ndq83/_bi83,4))
  for _bp in (16,48,60,70):
   _pb,_sb,_tb=_find_prime_Q(_bp,2000+_yi*100+_bp)
   _,_dgb,_=lll_16d_int(_struct_L(_a1,_b1,_a2,_b2,_pb,_sb,_tb),9,10);assert _dgb[16]==_pb*_pb
   _ls=[math.sqrt(_dgb[i+1]/_dgb[i]) for i in range(16)]
   _cwc=(_pb**(1.0/16.0))*(_alpha**(60.0/16.0));_lswc=[_cwc*(_alpha**(-j/2.0)) for j in range(16)]
   _gr=_b_int_64(_ls)/_b_int_64(_lswc);assert _gr<0.50;_gate_ratios.append(round(_gr,3))
 def _bnorm(re,im,rad,sh):
  mx=max(abs(re),abs(im),rad)
  if mx>=(1<<48):
   k=mx.bit_length()-48;mk=(1<<k)-1
   return (re>>k if re>=0 else -((-re)>>k)),(im>>k if im>=0 else -((-im)>>k)),((rad+mk)>>k)+2,sh+k
  return re,im,rad,sh
 def _badd(x,y):
  r1,i1,e1,s1=x;r2,i2,e2,s2=y
  if s1<s2:
   d=s2-s1
   if d>60:return _bnorm(r2,i2,e2+((abs(r1)+abs(i1)+e1+((1<<d)-1))>>d),s2)
   mk=(1<<d)-1;r1=r1>>d if r1>=0 else -((-r1)>>d);i1=i1>>d if i1>=0 else -((-i1)>>d);e1=((e1+mk)>>d)+2;s1=s2
  elif s2<s1:return _badd(y,x)
  return _bnorm(r1+r2,i1+i2,e1+e2,s1)
 def _bmul(x,y):
  r1,i1,e1,s1=x;r2,i2,e2,s2=y
  return _bnorm(r1*r2-i1*i2,r1*i2+i1*r2,(abs(r1)+abs(i1))*e2+(abs(r2)+abs(i2))*e1+e1*e2,s1+s2)
 with localcontext() as _c30:
  _c30.prec=30;_ln2=Decimal(2).ln()
  _LUT=[int(((Decimal(64+t+1)/Decimal(64)).ln()/_ln2*Decimal(1<<16)).to_integral_value(rounding=ROUND_CEILING)) for t in range(64)]
 def _bl2(x):
  re,im,rad,sh=x;M=math.isqrt(re*re+im*im)+1+rad;bl=M.bit_length()-1
  t=((M>>(bl-6))&63) if bl>=6 else ((M<<(6-bl))&63)
  return ((bl+sh)<<16)+_LUT[t]
 _w_el=y_pow_z(5);_wp1=add_phi15_z(_w_el,[1]+[0]*7);_w2=y_pow_z(10)
 _cells_bc=[([_rng.randint(-50,50) for _ in range(8)],[_rng.randint(-50,50) for _ in range(8)]) for _ in range(20)]
 for _k in range(5):
  _bb=[_rng.randint(-5,5) for _ in range(8)];_bb[0]=1
  _aa=[int(round(abs(U_bls)**(1.0/3.0)))*_bb[j] for j in range(8)];_cells_bc.append((_aa,_bb))
 _mgf=_mgg=0.0
 for _aa,_bb in _cells_bc:
  if not any(_aa) or not any(_bb):continue
  _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
  _ng_ex=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
  _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
  _nf_ex=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
  _bg_q=_bf_q=0
  for _k in range(4):
   _zt=_Z_TAB[_k]
   _sa=_bnorm(sum(_aa[j]*_zt[j][0] for j in range(8)),sum(_aa[j]*_zt[j][1] for j in range(8)),sum(abs(c) for c in _aa),-40)
   _sb=_bnorm(sum(_bb[j]*_zt[j][0] for j in range(8)),sum(_bb[j]*_zt[j][1] for j in range(8)),sum(abs(c) for c in _bb),-40)
   _a3b=_bmul(_bmul(_sa,_sa),_sa);_b3b=_bmul(_bmul(_sb,_sb),_sb)
   _wb=_bnorm(_zt[5][0],_zt[5][1],1,-40);_w2b=_bmul(_wb,_wb);_ub=_bnorm(U_bls,0,0,0)
   _tgb=_bmul(_bmul(_ub,_w2b),_b3b);_gv=_badd(_a3b,(-_tgb[0],-_tgb[1],_tgb[2],_tgb[3]))
   _a6b=_bmul(_a3b,_a3b);_b6b=_bmul(_b3b,_b3b);_a9b=_bmul(_a6b,_a3b);_b9b=_bmul(_b6b,_b3b)
   _wp1b=_badd(_wb,(1<<40,0,0,-40));_tfb=_bmul(_wb,_b9b)
   _fv=_badd(_badd(_badd(_a9b,_bmul(_wp1b,_bmul(_a6b,_b3b))),_bmul(_a3b,_b6b)),(-_tfb[0],-_tfb[1],_tfb[2],_tfb[3]))
   _bg_q+=2*_bl2(_gv);_bf_q+=2*_bl2(_fv)
  _bf=_bf_q/65536.0;_bg=_bg_q/65536.0
  if _ng_ex>0:assert _bg>=math.log2(_ng_ex);_mgg=max(_mgg,_bg-math.log2(_ng_ex))
  if _nf_ex>0:assert _bf>=math.log2(_nf_ex);_mgf=max(_mgf,_bf-math.log2(_nf_ex))
 assert _mgf<0.25 and _mgg<0.25
 _yr_min=1e99;_eff_u=[];_ss_yr=dlive["_ss_eval"]
 for _yi in (0,5,12,17):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_sn="g" if _sd==0 else "f"
  _qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16);_nb=_qq.bit_length()-1
  _prw=next(r for r in dlive["rows"] if r["side"]==_sn and r["lo"]<=_nb<r["hi"])
  _Rp=_prw["R"];_Rp2=_Rp*_Rp;_Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_qq;_rp=[pow(_ry,i,_qq) for i in range(8)]
  for _i in range(1,8):_Bq[_i][0]=(-_rp[_i])%_qq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=(_rt*_rp[_i])%_qq;_Bq[8+_i][8+_i]=1
  _Rr,_,_=lll_16d_int(_Bq,99,100);_a1,_b1=_Rr[0][:8],_Rr[0][8:]
  _a2,_b2=next((v[:8],v[8:]) for v in _Rr[1:] if any(sub_phi15_z(mul_phi15_z(_a1,v[8:]),mul_phi15_z(_b1,v[:8]))))
  _nfs,_ngs,_wts=[],[],[]
  for _ in range(12):
   while True:
    _g=[_rng.gauss(0,1) for _ in range(16)];_rn=math.sqrt(sum(x*x for x in _g))
    _rd=_Rp*(_rng.random()**(1.0/16.0))/_rn;_vv=[int(round(x*_rd)) for x in _g]
    if 0<sum(x*x for x in _vv)<=_Rp2:break
   _al,_bt=_vv[:8],_vv[8:]
   if not _is_cp(_al,_bt):_wts.append(0.0);_nfs.append(9999.0);_ngs.append(9999.0);continue
   _wts.append(1.0/_cnt_orb(_al,_bt,_Rp2))
   _aa=add_phi15_z(mul_phi15_z(_al,_a1),mul_phi15_z(_bt,_a2));_bb=add_phi15_z(mul_phi15_z(_al,_b1),mul_phi15_z(_bt,_b2))
   _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
   _ng=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
   _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
   _nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
   _nfs.append(math.log2(_nf)-(math.log2(_qq) if _sn=="f" else 0.0));_ngs.append(math.log2(_ng)-(math.log2(_qq) if _sn=="g" else 0.0))
  _pf=_ss_yr(_nfs,float(_prw["n"]),_prw["iA"],m_node,float(desc["params"]["alpha_f"]),_prw["kA"],_prw["kB"])
  _pg=_ss_yr(_ngs,float(_prw["n"]),_prw["iA"],m_node,float(desc["params"]["alpha_g"]),_prw["kA"],_prw["kB"])
  _es=sum(a*b for a,b in zip(_pf,_pg));_esm=sum(w*a*b for w,a,b in zip(_wts,_pf,_pg))
  _yr=(_esm/len(_wts))/_prw["pr"]
  if _es>0:_eff_u.append(round(_esm/_es,4))
  if _yr<_yr_min:_yr_min=_yr
 assert _yr_min>=10.0*float(desc["params"]["udup_node"])
 return round(min(_eff_u),4),round(_yr_min,2),len(_q32_h),dict(discrete_R64=_ds_disc,gate_8primes_max=max(_gate_ratios),gate_8primes=_gate_ratios,bound_int_max_gap_bits=dict(f=round(_mgf,4),g=round(_mgg,4)))
def verify_q15_93_b16_live(camp93,primes_upto,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 U_TOY=-5;W_EL=y_pow_z(5);WP1=add_phi15_z(W_EL,[1]+[0]*7);W2=y_pow_z(10)
 out={}
 for grp in ("g","f","enc"):
  rows=camp93[grp];assert len(rows)==40
  cnts=[r[4] if grp=="enc" else r[3] for r in rows]
  mus=[r[5] if grp=="enc" else r[4] for r in rows]
  mc=sum(cnts)/40.0;vc=sum((x-mc)**2 for x in cnts)/40.0;fano=vc/mc
  emp0=cnts.count(0)/40.0;mod0=sum(1.0/(1.0+0.5*m) for m in mus)/40.0
  wpgf=max((sum(t**c for c in cnts)/40.0)/(sum(1.0/(1.0+0.5*m*(1.0-t)) for m in mus)/40.0) for t in [i*0.05 for i in range(17)])
  assert fano<1.5 and emp0<mod0 and wpgf<0.90
  out[grp]=dict(n=40,mean_M=round(mc,2),mean_mu_raw=round(sum(mus)/40.0,2),fano=round(fano,2),zeros=cnts.count(0),emp0=round(emp0,4),h5_empty=round(mod0,4),worst_pgf=round(wpgf,4))
 P_SM=math.prod(primes_upto(256));P_MD=math.prod(p for p in primes_upto(8192) if p>256);NA=1<<24;NA3=NA**3
 pts=[((i,1),) for i in range(16)]
 for i in range(16):
  for j in range(i+1,16):
   pts.append(((i,1),(j,-1)));pts.append(((i,1),(j,1)))
   for k in range(j+1,16):
    pts.append(((i,1),(j,-1),(k,-1)));pts.append(((i,1),(j,-1),(k,1)))
    pts.append(((i,1),(j,1),(k,-1)));pts.append(((i,1),(j,1),(k,1)))
 _rng=random.Random(7)
 def _rho(n):
  if n%2==0:return 2
  while True:
   c=_rng.randrange(1,n);y=_rng.randrange(2,n);r=1;q=1;g=1
   while g==1:
    x=y
    for _ in range(r):y=(y*y+c)%n
    k=0
    while k<r and g==1:
     ys=y;m=min(64,r-k)
     for _ in range(m):y=(y*y+c)%n;q=(q*abs(x-y))%n
     g=math.gcd(q,n);k+=m
    r<<=1
   if g==n:
    while True:
     ys=(ys*ys+c)%n;g=math.gcd(abs(x-ys),n)
     if g>1:break
   if g!=n:return g
 def _ok_sw(n):
  if n<=0:return False
  for P in (P_SM,P_MD):
   while True:
    g=math.gcd(n,P)
    if g==1:break
    n//=g
   if n==1:return True
  if n>NA3:return False
  if _is_prime_det(n):return n<=NA
  st=[n];fc=0
  while st:
   m=st.pop()
   if m==1:continue
   if _is_prime_det(m):
    if m>NA:return False
    fc+=1
    if fc>3:return False
   else:d=_rho(m);st.append(d);st.append(m//d)
  return True
 wit_ok=0
 for sd in ("f",):
  p,ry,rt,exp_cnt,_=camp93[sd][0]
  B=[[0]*16 for _ in range(16)];B[0][0]=p;rp=[pow(ry,i,p) for i in range(8)]
  for i in range(1,8):B[i][0]=(-rp[i])%p;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*rp[i])%p;B[8+i][8+i]=1
  Rr,_,_=lll_16d_int(B,99,100);a1,b1=Rr[0][:8],Rr[0][8:]
  a2,b2=next((v[:8],v[8:]) for v in Rr[1:] if any(sub_phi15_z(mul_phi15_z(a1,v[8:]),mul_phi15_z(b1,v[:8]))))
  seen=set()
  ba=[mul_phi15_z(y_pow_z(i),a1) for i in range(8)]+[mul_phi15_z(y_pow_z(i),a2) for i in range(8)]
  bb_t=[mul_phi15_z(y_pow_z(i),b1) for i in range(8)]+[mul_phi15_z(y_pow_z(i),b2) for i in range(8)]
  for nz in pts:
   if len(nz)==1:i0,s0=nz[0];aa=ba[i0];bb=bb_t[i0]
   elif len(nz)==2:
    (i0,s0),(i1,s1)=nz;a0,a1=ba[i0],ba[i1];b0,b1=bb_t[i0],bb_t[i1]
    aa=[a0[t]+a1[t] if s1>0 else a0[t]-a1[t] for t in range(8)];bb=[b0[t]+b1[t] if s1>0 else b0[t]-b1[t] for t in range(8)]
   else:
    (i0,s0),(i1,s1),(i2,s2)=nz;a0,a1,a2=ba[i0],ba[i1],ba[i2];b0,b1,b2=bb_t[i0],bb_t[i1],bb_t[i2]
    aa=[a0[t]+s1*a1[t]+s2*a2[t] for t in range(8)];bb=[b0[t]+s1*b1[t]+s2*b2[t] for t in range(8)]
   a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3)
   ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_TOY,mul_phi15_z(W2,b3)))))
   if ng==0:continue
   if sd=="g":assert ng%p==0;ng//=p
   if not _ok_sw(ng):continue
   a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
   nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_EL,b9))))
   if nf==0:continue
   if sd=="f":assert nf%p==0;nf//=p
   if not _ok_sw(nf):continue
   g0=math.gcd(abs(norm_q15_z(aa)),abs(norm_q15_z(bb)));cp=True
   if g0>1:
    t_=g0;d_=2;ps_=[]
    while d_*d_<=t_ and d_<=10000:
     if t_%d_==0:
      ps_.append(d_)
      while t_%d_==0:t_//=d_
     d_+=1 if d_==2 else 2
    if t_>1:ps_.append(t_)
    for pp_ in ps_:
     g1_=pgcd_mod(aa,PHI15,pp_)
     if len(g1_)>1 and len(pgcd_mod(g1_,bb,pp_))>1:cp=False;break
   if cp:seen.add((nf,ng))
  assert len(seen)==exp_cnt;wit_ok+=1
 out["live_regenerated_93_b16_balls"]=wit_ok;out["checks_ok"]=True
 return out
def verify_prod_cond_q15(enc12,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,lll_16d_int):
 assert len(enc12)==12
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 log2_u8=float(desc["params"]["log2_u8"]);S8=math.log2(desc["cells_rank2"]["8"])
 qb_f=[x[1] for x in desc["qbins"]];qb_g=[x[2] for x in desc["qbins"]]
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4) for _ in range(8)];bt=[rng.randint(-4,4) for _ in range(8)]
  if any(al) and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
 cfs,cgs,pr_rats=[],[],[]
 _ss=dlive["_ss_eval"];_mn=dlive["mesh_node"];_af=float(desc["params"]["alpha_f"]);_ag=float(desc["params"]["alpha_g"])
 for sid,q_hex,pa,pb,ry_hex,rt_hex in enc12:
  q,ry,rt=int(q_hex,16),int(ry_hex,16),int(rt_hex,16)
  pa3,pb3=pow_fn(pa,3),pow_fn(pb,3)
  if sid==0:pn=abs(norm_fn(sub_fn(pa3,scale_fn(U_bls,mul_fn(w2_el,pb3)))))
  else:
   pa6,pb6=mul_fn(pa3,pa3),mul_fn(pb3,pb3);pa9,pb9=mul_fn(pa6,pa3),mul_fn(pb6,pb3)
   pn=abs(norm_fn(sub_fn(add_fn(add_fn(pa9,mul_fn(wp1_el,mul_fn(pa6,pb3))),mul_fn(pa3,pb6)),mul_fn(w_el,pb9))))
  assert pn>0 and pn%q==0
  ry_p=[pow(ry,i,q) for i in range(8)]
  assert (sum(pa[i]*ry_p[i] for i in range(8))-rt*sum(pb[i]*ry_p[i] for i in range(8)))%q==0
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  for i in range(1,8):B[i][0]=(-ry_p[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*ry_p[i])%q;B[8+i][8+i]=1
  R,dsq,lsq=lll_16d_int(B,99,100);assert dsq[16]==q*q
  for kk in range(1,16):assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  a1,b1=R[0][:8],R[0][8:]
  a2,b2=next((v[:8],v[8:]) for v in R[1:] if any(sub_fn(mul_fn(a1,v[8:]),mul_fn(b1,v[:8]))))
  sname="g" if sid==0 else "f";nb=q.bit_length()-1;lq=math.log2(q)
  prow=next(r for r in dlive["rows"] if r["side"]==sname and r["lo"]<=nb<r["hi"])
  S_bin=math.log2(desc["cells_rank2"][str(prow["R"])]);c_act=(lq+S_bin)/16.0
  c_val=(lq+S8)/16.0;nfas,ngas=[],[]
  for al,bt in ab_pts:
   a=add_fn(mul_fn(al,a1),mul_fn(bt,a2));b=add_fn(mul_fn(al,b1),mul_fn(bt,b2))
   a3,b3=pow_fn(a,3),pow_fn(b,3)
   ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
   a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
   nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
   assert (ng if sid==0 else nf)%q==0
   cf=math.log2(nf)-72.0*c_val;cg=math.log2(ng)-24.0*c_val-log2_u8
   cfs.append(cf);cgs.append(cg)
   nfas.append(72.0*c_act+cf-(lq if sname=="f" else 0.0));ngas.append(24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0))
  _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,_af,prow["kA"],prow["kB"])
  _pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,_ag,prow["kA"],prow["kB"])
  pr_rats.append((sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts))/prow["pr"])
 cfs.sort();cgs.sort();n=len(cfs)
 sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
 sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
 assert all(sf6[i]<qb_f[i] for i in range(6)) and all(sg6[i]<qb_g[i] for i in range(6)) and min(pr_rats)>=10.0
 _prs=sorted(pr_rats)
 return dict(n_enc_live=12,parent_divisibility_live="12/12",live_bareiss_99_100_all12=True,
  presupplied_vectors=0,enc_sext_f=sf6,enc_sext_g=sg6,
  min_headroom_f=round(min(qb_f[i]-sf6[i] for i in range(6)),2),
  min_headroom_g=round(min(qb_g[i]-sg6[i] for i in range(6)),2),
  joint_pr_over_charged_min=round(_prs[0],2),joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),
  all_sextiles_below=True,checks_ok=True)
def unplanted_two_sided_bw_check():
 import random
 p,r,m=467,233,381
 f_c,g_c=[1,1,-5,0,1],[-10,-1,1]
 assert sum(ci*pow(m,i,p) for i,ci in enumerate(f_c))%p==0
 assert sum(ci*pow(m,i,p) for i,ci in enumerate(g_c))%p==0
 eps_f,eps_g,r2=r-1,r*r-1,r*r
 def peval(c,x,mod):
  s=0
  for ci in reversed(c):
   s=(s*x+ci)%mod
  return s
 def pmul(a,b,mod_p,q):
  d=len(mod_p)-1
  out=[0]*(len(a)+len(b)-1)
  for i,ai in enumerate(a):
   if ai:
    for j,bj in enumerate(b):
     out[i+j]=(out[i+j]+ai*bj)%q
  for k in range(len(out)-1,d-1,-1):
   ck=out[k]
   if ck:
    for j in range(d):
     out[k-d+j]=(out[k-d+j]-ck*mod_p[j])%q
  return out[:d]
 def ppow(base,exp,mod_p,q):
  d=len(mod_p)-1
  res=[1]+[0]*(d-1)
  b=[x%q for x in base]+[0]*(d-len(base))
  e=exp
  while e:
   if e & 1:res=pmul(res,b,mod_p,q)
   e >>=1
   if e:b=pmul(b,b,mod_p,q)
  return res
 primes=[q for q in range(2,115) if all(q%d for d in range(2,int(q**0.5)+1))]
 def unram(c,q):
  for rt in range(q):
   if peval(c,rt,q)==0:
    if sum(i*c[i]*pow(rt,i-1,q) for i in range(1,len(c)))%q==0:return False
  return True
 fb_f,fb_g={},{}
 for q in primes:
  if unram(f_c,q):
   for rt in range(q):
    if peval(f_c,rt,q)==0:fb_f[(q,rt)]=len(fb_f)
  if unram(g_c,q):
   for rt in range(q):
    if peval(g_c,rt,q)==0:fb_g[(q,rt)]=len(fb_g)
 pf_act=sorted({q for (q,_) in fb_f})
 pg_act=sorted({q for (q,_) in fb_g})
 def sfac(val,plist):
  n=abs(val)
  if n==0:return None
  fac={}
  for q in plist:
   if n%q==0:
    e=0
    while n%q==0:
     n//=q
     e+=1
    fac[q]=e
    if n==1:return fac
  return fac if n==1 else None
 rels=[]
 for b in range(1,55):
  for a in range(-450,451):
   if a==0 or math.gcd(abs(a),b)!=1 or (a-b*m)%p==0:continue
   fg=sfac(a*a-a*b-10*b*b,pg_act)
   if fg is None:continue
   ff=sfac(a**4-5*a*a*b*b+a*b**3+b**4,pf_act)
   if ff is None:continue
   vf={fb_f[(q,(a*pow(b,-1,q))%q)]:e for q,e in ff.items()}
   vg={fb_g[(q,(a*pow(b,-1,q))%q)]:e for q,e in fg.items()}
   rels.append((a,b,vf,vg))
 nf_cnt=len(fb_f)
 cur=list(range(len(rels)))
 while True:
  ccnt={}
  for ri in cur:
   for k in rels[ri][2]:
    ccnt[k]=ccnt.get(k,0)+1
   for k in rels[ri][3]:
    ccnt[nf_cnt+k]=ccnt.get(nf_cnt+k,0)+1
  nxt=[ri for ri in cur if all(ccnt[k]>=2 for k in rels[ri][2]) and all(ccnt[nf_cnt+k]>=2 for k in rels[ri][3])]
  if len(nxt)==len(cur):break
  cur=nxt
 used=sorted({k for ri in cur for k in rels[ri][2]} | {nf_cnt+k for ri in cur for k in rels[ri][3]})
 cmap={c:i for i,c in enumerate(used)}
 N_2core=len(used)
 def smap(a,b,mod_p,eps):
  el=ppow([a%r2,(-b)%r2],eps,mod_p,r2)
  assert el[0]%r==1 and all(x%r==0 for x in el[1:])
  return [((el[0]-1)//r)%r]+[(x//r)%r for x in el[1:]]
 rows_M,rows_S=[],[]
 for ri in cur:
  a,b,vf,vg=rels[ri]
  row=[0]*N_2core
  for k,e in vf.items():
   row[cmap[k]]=e%r
  for k,e in vg.items():
   row[cmap[nf_cnt+k]]=(-e)%r
  rows_M.append(row)
  rows_S.append(smap(a,b,f_c,eps_f)[1:4]+smap(a,b,g_c,eps_g)[1:2])
 m_rows=[row[:] for row in rows_M];s_rows=[row[:] for row in rows_S]
 alive_r=[True]*len(cur);alive_c=[True]*N_2core
 merge_stack=[]
 while True:
  piv_c=next((c for c in range(N_2core) if alive_c[c] and sum(1 for i in range(len(cur)) if alive_r[i] and m_rows[i][c]!=0)==2),-1)
  if piv_c<0:break
  r_pair=[i for i in range(len(cur)) if alive_r[i] and m_rows[i][piv_c]!=0]
  rp,ro=r_pair[0],r_pair[1]
  merge_stack.append((piv_c,m_rows[rp][:],s_rows[rp][:]))
  fac=(m_rows[ro][piv_c]*pow(m_rows[rp][piv_c],-1,r))%r
  m_rows[ro]=[(m_rows[ro][j]-fac*m_rows[rp][j])%r for j in range(N_2core)]
  s_rows[ro]=[(s_rows[ro][j]-fac*s_rows[rp][j])%r for j in range(4)]
  alive_r[rp]=False;alive_c[piv_c]=False
 core_cols=[c for c in range(N_2core) if alive_c[c]]
 N,D=len(core_cols),len(core_cols)+4
 active=[i for i in range(len(cur)) if alive_r[i]]
 cc={c:sum(1 for i in active if m_rows[i][c]!=0) for c in core_cols}
 while len(active)>D:
  cands=[i for i in active if all(cc[c]>=3 for c in core_cols if m_rows[i][c]!=0)] or [i for i in active if all(cc[c]>=2 for c in core_cols if m_rows[i][c]!=0)]
  drop=max(cands,key=lambda i:(sum(1 for c in core_cols if m_rows[i][c]!=0),i))
  active.remove(drop)
  for c in core_cols:
   if m_rows[drop][c]!=0:cc[c]-=1
 def rref_null(mat):
  mr=[row[:] for row in mat]
  R,C=len(mr),len(mr[0])
  pivs,rk=[],0
  for c in range(C):
   piv=next((i for i in range(rk,R) if mr[i][c]%r!=0),-1)
   if piv<0:continue
   mr[rk],mr[piv]=mr[piv],mr[rk]
   inv=pow(mr[rk][c],-1,r)
   rr=mr[rk]
   for j in range(c,C):
    rr[j]=(rr[j]*inv)%r
   for i in range(R):
    if i!=rk and mr[i][c]:
     fv=mr[i][c]
     ri_r=mr[i]
     for j in range(c,C):
      ri_r[j]=(ri_r[j]-fv*rr[j])%r
   pivs.append(c)
   rk+=1
  pset=set(pivs)
  vecs=[]
  for fc in [c for c in range(C) if c not in pset]:
   v=[0]*C
   v[fc]=1
   for ri_i,pc in enumerate(pivs):
    v[pc]=(-mr[ri_i][fc])%r
   vecs.append(v)
  return rk,vecs
 M_sq=[[m_rows[i][c] for c in core_cols] for i in active]
 S_sq=[s_rows[i][:] for i in active]
 A_sq=[M_sq[i]+S_sq[i] for i in range(D)]
 rk_M,_=rref_null(M_sq)
 rk_A,ker_A=rref_null(A_sq)
 rk_3sm,ker_3sm=rref_null([M_sq[i]+S_sq[i][:3] for i in range(D)])
 assert rk_M==N and rk_A==D-1 and len(ker_A)==1 and len(ker_3sm)==0
 def mmul(P,Q):
  R,K,C=len(P),len(P[0]),len(Q[0])
  out=[[0]*C for _ in range(R)]
  for i in range(R):
   Pi,oi=P[i],out[i]
   for k in range(K):
    v=Pi[k]
    if v:
     Qk=Q[k]
     for j in range(C):
      oi[j]=(oi[j]+v*Qk[j])%r
  return out
 def mvec(P,v):
  return [sum(P[i][j]*v[j] for j in range(len(v)))%r for i in range(len(P))]
 rng=random.Random(20261119)
 m_bw=n_bw=4
 bnds=[(k*D)//4 for k in range(5)]
 X_good=[[0]*D for _ in range(m_bw)]
 for k in range(m_bw):
  for j in range(bnds[k],bnds[k+1]):
   X_good[k][j]=rng.randrange(1,r)
 Y=[[rng.randrange(1,r) for _ in range(n_bw)] for _ in range(D)]
 hr,hc=(D//m_bw)+2,(D//n_bw)+2
 AY,seq=[row[:] for row in Y],[]
 for _ in range(hr+hc+2):
  seq.append(mmul(X_good,AY))
  AY=mmul(A_sq,AY)
 H_rows=[[v for j in range(hc) for v in seq[i+j+1][ri_b]] for i in range(hr) for ri_b in range(m_bw)]
 rk_H,ker_H=rref_null(H_rows)
 v_rec=None
 for c_vec in ker_H:
  acc=[0]*D
  for j in range(hc-1,-1,-1):
   cj=c_vec[j*n_bw:(j+1)*n_bw]
   yc=[sum(Y[ri_i][bi]*cj[bi] for bi in range(n_bw))%r for ri_i in range(D)]
   acc=[(a_i+y_i)%r for a_i,y_i in zip(mvec(A_sq,acc),yc)]
  if any(acc) and all(x==0 for x in mvec(A_sq,acc)):
   v_rec=acc
   break
 assert v_rec is not None
 O_rows,cX=[],[row[:] for row in X_good]
 for _ in range(hr):
  cX=mmul(cX,A_sq)
  O_rows.extend(cX)
 rk_OA,_=rref_null(O_rows)
 C_cols,cY=[],[row[:] for row in Y]
 for _ in range(hc):
  for bi in range(n_bw):
   C_cols.append([cY[ri_i][bi] for ri_i in range(D)])
  cY=mmul(A_sq,cY)
 rk_Cd,_=rref_null(C_cols)
 X_bad=[[((k+1)*X_good[0][j])%r for j in range(D)] for k in range(m_bw)]
 O_bad,cXb=[],[row[:] for row in X_bad]
 for _ in range(hr):
  cXb=mmul(cXb,A_sq)
  O_bad.extend(cXb)
 rk_OA_bad,_=rref_null(O_bad)
 vlog_full={c:v_rec[j] for j,c in enumerate(core_cols)}
 sm_rec=v_rec[N:]
 for piv_c,pm,ps in reversed(merge_stack):
  acc=(sum(pm[j]*vlog_full[j] for j in range(N_2core) if j!=piv_c and pm[j]!=0)+sum(ps[s]*sm_rec[s] for s in range(4)))%r
  vlog_full[piv_c]=(-acc*pow(pm[piv_c],-1,r))%r
 while True:
  prog=False
  for idx in range(len(cur)):
   unk=[j for j in range(N_2core) if rows_M[idx][j]!=0 and j not in vlog_full]
   if len(unk)==1:
    cn=unk[0]
    acc=(sum(rows_M[idx][j]*vlog_full[j] for j in range(N_2core) if j!=cn and rows_M[idx][j]!=0)+sum(rows_S[idx][s]*sm_rec[s] for s in range(4)))%r
    vlog_full[cn]=(-acc*pow(rows_M[idx][cn],-1,r))%r
    prog=True
  if not prog:break
 assert len(vlog_full)==N_2core
 gen=next(g for g in range(2,p) if pow(g,(p-1)//r,p)!=1 and pow(g,2,p)!=1)
 g_r,dlog,ag=pow(gen,(p-1)//r,p),{},1
 for e in range(r):
  dlog[ag]=e
  ag=(ag*g_r)%p
 ratios=set()
 for idx,ri in enumerate(cur):
  a,b,vf,vg=rels[ri]
  lhs=(sum(e*vlog_full[cmap[k]] for k,e in vf.items())+sum(rows_S[idx][j]*sm_rec[j] for j in range(3)))%r
  rhs=(sum(e*vlog_full[cmap[nf_cnt+k]] for k,e in vg.items())-rows_S[idx][3]*sm_rec[3])%r
  assert lhs==rhs
  dl=dlog[pow((a-b*m)%p,(p-1)//r,p)]
  if dl!=0:ratios.add((lhs*pow(dl,-1,r))%r)
 pm_res=pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_rec)
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==D and rk_OA_bad<D-1 and len(ratios)==1
 return dict(p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(merge_stack),N=N,D=D,rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=rk_OA_bad,**pm_res,unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(ratios),checks_ok=True)
