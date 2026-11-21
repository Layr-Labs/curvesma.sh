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
 if not a or not b:return[]
 la,lb=len(a),len(b);n=la+lb-1 if cap is None else min(la+lb-1,cap);o=[0]*n
 for i,ai in enumerate(a):
  if ai:
   for j in range(min(lb,n-i)):o[i+j]=(o[i+j]+ai*b[j])%p
 while o and o[-1]==0:o.pop()
 return o
def _pmat_mul(A,B,p,cap=None):
 r,k,c=len(A),len(A[0]),len(B[0])
 C=[[[]for _ in range(c)]for _ in range(r)]
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
 P=[[[]for _ in range(t)]for _ in range(t)]
 for j in range(t):
  for i in range(t):
   v=Pcols[j][i]
   if v:P[i][j]=[0,v]if dinc[j]else[v]
 return P,[s[j]+dinc[j]for j in range(t)]
def pm_basis(F,sigma,s,p):
 """Giorgi-Jeannerod-Villard (2003) / Thome (2002) divide-and-conquer PM-Basis over F_p[z]."""
 m,t=len(F),len(F[0])
 if sigma==1:return _mbasis1([[F[i][j][0]if F[i][j]else 0 for j in range(t)]for i in range(m)],s,p)
 s1=sigma//2;s2=sigma-s1
 P1,sm=pm_basis([[F[i][j][:s1]for j in range(t)]for i in range(m)],s1,s,p)
 FP1=_pmat_mul(F,P1,p,cap=sigma)
 G=[[FP1[i][j][s1:sigma]for j in range(t)]for i in range(m)]
 P2,so=pm_basis(G,s2,sm,p)
 return _pmat_mul(P1,P2,p),so
def pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_ref):
 """Run live divide-and-conquer PM-Basis on S_k = X A^k Y (k>=1), recover exact ker(A) via Horner."""
 m=len(seq[0]);n=len(seq[0][0]);D=len(A_sq);L=hr+hc-1
 S=seq[1:1+L]
 M=[[[]for _ in range(n+m)]for _ in range(m)]
 for i in range(m):
  for j in range(n):
   poly=[S[k][i][j]for k in range(L)]
   while poly and poly[-1]==0:poly.pop()
   M[i][j]=poly
  M[i][n+i]=[r-1]
 P_b,s_out=pm_basis(M,L,[0]*n+[1]*m,r)
 def mvec(P,v):return[sum(P[i][j]*v[j]for j in range(len(v)))%r for i in range(len(P))]
 v_pm=None;deg_used=None
 for col in sorted(range(n+m),key=lambda j:s_out[j]):
  if s_out[col]>hc-1:continue
  q=[P_b[i][col]for i in range(n)]
  dg=max((len(poly)-1 for poly in q),default=-1)
  if dg<0:continue
  acc=[0]*D
  for t in range(dg+1):
   qt=[q[bi][t]if t<len(q[bi])else 0 for bi in range(n)]
   yc=[sum(Y[ri][bi]*qt[bi]for bi in range(n))%r for ri in range(D)]
   acc=[(a+y)%r for a,y in zip(mvec(A_sq,acc),yc)]
  if any(acc)and all(x==0 for x in mvec(A_sq,acc)):
   v_pm=acc;deg_used=dg;break
 assert v_pm is not None,"PM-Basis failed to recover ker(A)!"
 idx0=next(i for i,x in enumerate(v_ref)if x)
 ratio=(v_pm[idx0]*pow(v_ref[idx0],-1,r))%r
 assert ratio!=0 and all((v_ref[i]*ratio)%r==v_pm[i]for i in range(D))
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
  num_k=-sum(lm[j][k]*z[j]for j in range(k+1,16));dk1=dg[k+1];dk=dg[k]
  rhs=rem_sc*dk;rad=math.isqrt(rhs)
  lo=(num_k-rad+dk1-1)//dk1;hi=(num_k+rad)//dk1
  if not nz:lo=max(0,lo)
  for zk in range(lo,hi+1):
   z[k]=zk;nz2=nz or(zk>0);diff=zk*dk1-num_k;rnum=rhs-diff*diff
   if rnum<0:continue
   if k==0:
    if not nz2:continue
    hits+=1
    if collect:
     v=[sum(z[i]*Br[i][t]for i in range(16))for t in range(16)]
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
 S=1<<32;L63=(1<<63)-1
 def _i64(v):assert-L63-1<=v<=L63;return v
 MU=[[0]*16 for _ in range(16)]
 for j in range(16):
  for k in range(j):
   dk1=dg[k+1];num=lm[j][k]*S
   MU[j][k]=_i64((num+dk1//2)//dk1 if num>=0 else-((-num+dk1//2)//dk1))
 L2_LO=[_i64((dg[k+1]*S)//dg[k])for k in range(16)]
 INV_HI=[_i64(math.isqrt((dg[k]*S*S+dg[k+1]-1)//dg[k+1])+1)for k in range(16)]
 C_stk=[[0]*16 for _ in range(16)]
 z=[0]*16;hi_s=[0]*16;rem_s=[0]*16;nz_s=[0]*16;asum_s=[0]*16;first_s=[True]*16
 nbk=[0]*16;ops=0;hits=0;out=[]
 def _enter(k,rem_q,nz,asum):
  nonlocal ops
  cq=C_stk[k][k];ec=(asum+1)>>1
  if rem_q==0:sr=1
  else:
   x2=1<<((65-(64-rem_q.bit_length()))>>1)
   for _ in range(5):x2=(x2+(rem_q//x2))>>1
   sr=x2+1
  rq=_i64(((sr*INV_HI[k])>>16)+1+ec)
  lo=_i64(cq-rq+S-1)>>32;hi=_i64(cq+rq)>>32
  if not nz and lo<0:lo=0
  z[k]=lo;hi_s[k]=hi;rem_s[k]=rem_q;nz_s[k]=nz;asum_s[k]=asum;first_s[k]=True;ops+=38
 k=15;_enter(15,R2<<32,0,0)
 while True:
  if z[k]>hi_s[k]:
   ops+=12;k+=1
   if k==16:break
   z[k]+=1;continue
  zk=z[k];cq=C_stk[k][k];ec=(asum_s[k]+1)>>1
  da=abs(_i64(_i64(zk<<32)-cq));dl=da-ec if da>ec else 0
  sq=(((dl*dl)>>32)*L2_LO[k])>>32;rem_q=rem_s[k]
  if sq>rem_q:ops+=20;z[k]+=1;continue
  nz2=1 if(nz_s[k]or zk>0)else 0
  if k==0:
   if nz2:
    v=[sum(z[i]*Br[i][t]for i in range(16))for t in range(16)]
    if sum(x*x for x in v)<=R2:
     hits+=1
     if collect:
      for x in v:
       if x!=0:
        if x<0:v=[-t for t in v]
        break
      out.append(tuple(v))
   z[0]+=1
  else:
   nbk[k]+=1
   if first_s[k]:
    for m in range(k):C_stk[k-1][m]=_i64(C_stk[k][m]-_i64(MU[k][m]*zk))
    first_s[k]=False
   else:
    for m in range(k):C_stk[k-1][m]=_i64(C_stk[k-1][m]-MU[k][m])
   ops+=38+(6*k+4);k-=1;_enter(k,rem_q-sq,nz2,asum_s[k+1]+abs(zk))
 assert ops<=sum((132+6*j)*nbk[j]for j in range(1,16))+60
 return hits,sum(nbk),out,ops
def verify_b16_dsieve_live(data,desc,dlive,m_node,U_bls,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 _Rm=max(int(k)for k in desc["cells_rank2"]);_M=_Rm*_Rm
 def _ksq_b(poly,bpc):
  _x=int.from_bytes(b"".join(v.to_bytes(bpc,"little")for v in poly),"little");_p=_x*_x
  _nb=max((_M+1)*bpc,(_p.bit_length()+7)//8);_bs=_p.to_bytes(_nb,"little")
  return[int.from_bytes(_bs[i*bpc:(i+1)*bpc],"little")for i in range(_M+1)]
 _t1=[0]*(_M+1);_t1[0]=1
 for _k in range(1,_Rm+1):_t1[_k*_k]=2
 _t16=_ksq_b(_ksq_b(_ksq_b(_ksq_b(_t1,4),8),8),16)
 for _rk,_cv in desc["cells_rank2"].items():
  _R=int(_rk);assert int(_cv)==(sum(_t16[:_R*_R+1])-1)//2,f"B_16(R)/+-1 mismatch R={_R}"
 assert 4*max(int(desc["cells_rank2"][str(r["R"])])for r in dlive["rows"])<(1<<111)
 _yp=[y_pow_z(k)for k in range(15)]
 _ug=[[1,0,0,1,0,0,0,0],[-1,0,1,0,0,0,0,0],[-1,1,0,0,0,0,0,0]]
 _ui=[[1,0,-1,0,0,0,1,-1],[0,-1,0,0,-1,1,-1,0],[0,-1,-1,0,-1,0,0,-1]]
 for _a,_b in zip(_ug,_ui):assert mul_phi15_z(_a,_b)==[1]+[0]*7
 _ym=[tuple(tuple(mul_phi15_z(y,[1 if i==j else 0 for i in range(8)])[r]for j in range(8))for r in range(8))for y in _yp]
 _pc=[{0:[1]+[0]*7}for _ in range(3)]
 def _upow(_idx,_n):
  _dd=_pc[_idx]
  if _n in _dd:return _dd[_n]
  _rv=mul_phi15_z(_upow(_idx,_n-1),_ug[_idx])if _n>0 else mul_phi15_z(_upow(_idx,_n+1),_ui[_idx])
  _dd[_n]=_rv;return _rv
 _CL=Context(prec=40,rounding=ROUND_FLOOR);_CH=Context(prec=40,rounding=ROUND_CEILING)
 def _iv_add(a,b):return(_CL.add(a[0],b[0]),_CH.add(a[1],b[1]))
 def _iv_mul(a,b):
  a0,a1=a;b0,b1=b
  return(min(_CL.multiply(a0,b0),_CL.multiply(a0,b1),_CL.multiply(a1,b0),_CL.multiply(a1,b1)),max(_CH.multiply(a0,b0),_CH.multiply(a0,b1),_CH.multiply(a1,b0),_CH.multiply(a1,b1)))
 def _iv_div(a,b):
  assert b[0]>0
  return(min(_CL.divide(a[0],b[0]),_CL.divide(a[0],b[1]),_CL.divide(a[1],b[0]),_CL.divide(a[1],b[1])),max(_CH.divide(a[0],b[0]),_CH.divide(a[0],b[1]),_CH.divide(a[1],b[0]),_CH.divide(a[1],b[1])))
 def _ln_atanh(t_iv):
  assert abs(t_iv[0])<=Decimal("0.35")and abs(t_iv[1])<=Decimal("0.35")
  t2=_iv_mul(t_iv,t_iv);pw=t_iv;s=(Decimal(0),Decimal(0))
  for k in range(38):dk=(Decimal(2*k+1),Decimal(2*k+1));s=_iv_add(s,_iv_div(pw,dk));pw=_iv_mul(pw,t2)
  s=_iv_mul((Decimal(2),Decimal(2)),s);ep=Decimal("1e-34")
  return(_CL.subtract(s[0],ep),_CH.add(s[1],ep))
 _LN2_IV=_ln_atanh(_iv_div((Decimal(1),Decimal(1)),(Decimal(3),Decimal(3))))
 def _iv_ln(a):
  assert a[0]>0 and a[0]<=a[1]
  def _pt(x):
   e=int(math.floor(math.log2(float(x))))
   if float(x)/(2.0**e)>1.4:e+=1
   if e>=0:p2=Decimal(1<<e);m=(_CL.divide(x,p2),_CH.divide(x,p2))
   else:p2=Decimal(1<<(-e));m=(_CL.multiply(x,p2),_CH.multiply(x,p2))
   tlo=_CL.divide(_CL.subtract(m[0],Decimal(1)),_CH.add(m[0],Decimal(1)))
   thi=_CH.divide(_CH.subtract(m[1],Decimal(1)),_CL.add(m[1],Decimal(1)))
   return _iv_add(_ln_atanh((tlo,thi)),_iv_mul((Decimal(e),Decimal(e)),_LN2_IV))
  return(_pt(a[0])[0],_pt(a[1])[1])
 def _iv_exp(a):
  assert a[0]<=a[1]
  def _pt(x):
   s=max(0,int(math.ceil(math.log2(max(1.0,abs(float(x))*4.0)))));sc=Decimal(1<<s)
   yv=(_CL.divide(x,sc),_CH.divide(x,sc));rs=(Decimal(1),Decimal(1));tm=rs
   for k in range(1,25):dk=(Decimal(k),Decimal(k));tm=_iv_div(_iv_mul(tm,yv),dk);rs=_iv_add(rs,tm)
   ep=Decimal("1e-35");rs=(max(Decimal(0),_CL.subtract(rs[0],ep)),_CH.add(rs[1],ep))
   for _ in range(s):rs=_iv_mul(rs,rs)
   return rs
  return(_pt(a[0])[0],_pt(a[1])[1])
 _t_ln2=_iv_ln((Decimal(2),Decimal(2)));assert _t_ln2[0]<_t_ln2[1]and _t_ln2[1]-_t_ln2[0]<Decimal("1e-25")
 _t_e1=_iv_exp((Decimal(1),Decimal(1)));assert _t_e1[0]<_t_e1[1]and _t_e1[1]-_t_e1[0]<Decimal("1e-30")
 def _iv_cs(m_red):
  with localcontext()as c:
   c.prec=45
   pi_m=Decimal("3.1415926535897932384626433832795028841971")
   th=Decimal(2*m_red)*pi_m/Decimal(15);cv=Decimal(1);sv=th;tc=Decimal(1);ts=th;t2=-th*th
   for k in range(1,28):tc=tc*t2/((2*k-1)*(2*k));cv+=tc;ts=ts*t2/((2*k)*(2*k+1));sv+=ts
   eps=Decimal("1e-35")
  clo,slo=_CL.subtract(cv,eps),_CL.subtract(sv,eps);chi,shi=_CH.add(cv,eps),_CH.add(sv,eps)
  return(clo,chi),(slo,shi)
 _zp_iv=[];_Z_TAB=[]
 _S58=1<<58
 for _kk in (1,7,4,2):
  _row_iv=[];_row_t=[]
  for _j in range(8):
   _mr=(_kk*_j)%15
   if _mr>7:_mr-=15
   _ci,_si=_iv_cs(_mr);_row_iv.append((_ci,_si))
   with localcontext()as c:
    c.prec=40
    _cn=int((_ci[0]*_S58+Decimal("0.5")).to_integral_value(rounding=ROUND_FLOOR))
    _sn=int((_si[0]*_S58+Decimal("0.5")).to_integral_value(rounding=ROUND_FLOOR))
    assert abs(_ci[0]*_S58-_cn)<=Decimal("0.51")and abs(_ci[1]*_S58-_cn)<=Decimal("0.51")
    assert abs(_si[0]*_S58-_sn)<=Decimal("0.51")and abs(_si[1]*_S58-_sn)<=Decimal("0.51")
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
   lo=Decimal(0)if iv[0]<=0<=iv[1]else min(_CL.multiply(iv[0],iv[0]),_CL.multiply(iv[1],iv[1]))
   return(lo,max(_CH.multiply(iv[0],iv[0]),_CH.multiply(iv[1],iv[1])))
  return _iv_add(_sq(_re),_sq(_im))
 _Lg_iv=[[_iv_mul((Decimal("0.5"),Decimal("0.5")),_iv_ln(_iv_norm2(u,k)))for k in range(4)]for u in _ug]
 _Mm_iv=[[_Lg_iv[i][k]for i in range(3)]for k in range(3)]
 with localcontext()as _c45:
  _c45.prec=45
  _Mm0=[[(_Mm_iv[r][c][0]+_Mm_iv[r][c][1])/Decimal(2)for c in range(3)]for r in range(3)]
  _a,_b,_c=_Mm0[0];_d,_e,_f=_Mm0[1];_g,_h,_ii=_Mm0[2]
  _dt=_a*(_e*_ii-_f*_h)-_b*(_d*_ii-_f*_g)+_c*(_d*_h-_e*_g)
  assert abs(Decimal(8)*abs(_dt)-Decimal("4.6618207772966835"))<Decimal("1e-12")
  _Mi0=[[(_e*_ii-_f*_h)/_dt,(_c*_h-_b*_ii)/_dt,(_b*_f-_c*_e)/_dt],[(_f*_g-_d*_ii)/_dt,(_a*_ii-_c*_g)/_dt,(_c*_d-_a*_f)/_dt],[(_d*_h-_e*_g)/_dt,(_b*_g-_a*_h)/_dt,(_a*_e-_b*_d)/_dt]]
  _mi_inf=max(sum(abs(_Mi0[r][c])for c in range(3))for r in range(3))
 _eps_E=Decimal(0)
 for _r in range(3):
  _rs=Decimal(0)
  for _cl in range(3):
   _s=(Decimal(0),Decimal(0))
   for _j in range(3):_s=_iv_add(_s,_iv_mul(_Mm_iv[_r][_j],(_Mi0[_j][_cl],_Mi0[_j][_cl])))
   _tgt=Decimal(1)if _r==_cl else Decimal(0)
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
 _ep_lo=[{}for _ in range(3)];_SC=Decimal(10**9);_SC4=10**36
 def _epow_lo(_idx,_n):
  if _n not in _ep_lo[_idx]:
   _out=[]
   for _k in range(4):
    _iv=_Lg_iv[_idx][_k]
    _ex=_CL.multiply(Decimal(2*_n),_iv[0]if _n>=0 else _iv[1])
    _out.append(int(_CL.multiply(_iv_exp((_ex,_ex))[0],_SC).to_integral_value(rounding=ROUND_FLOOR)))
   _ep_lo[_idx][_n]=_out
  return _ep_lo[_idx][_n]
 _d15001=Decimal("15.001")
 def _cnt_orb(_va,_vb,_R2):
  _eks=[_iv_mul((Decimal(2),Decimal(2)),_iv_add(_iv_norm2(_va,k),_iv_norm2(_vb,k)))for k in range(4)]
  _bd_hi=_CH.multiply(_d15001,Decimal(_R2));_bd_i=int(_CH.multiply(_bd_hi,Decimal(_SC4)).to_integral_value(rounding=ROUND_CEILING))
  _eks_i=[int(_CL.multiply(_ek[0],_SC).to_integral_value(rounding=ROUND_FLOOR))for _ek in _eks]
  _Bk=[_CH.multiply(_iv_ln((_r,_r))[1],Decimal("0.5"))for _ek in _eks for _r in(_CH.divide(_bd_hi,_ek[0]),)]
  _SB=Decimal(0)
  for _b in _Bk:_SB=_CH.add(_SB,_b)
  _lo=[_CL.minus(_CH.subtract(_SB,_Bk[k]))for k in range(3)]
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
    _e1=_epow_lo(1,_n1);_e01=[_eks_i[k]*_e0[k]*_e1[k]for k in range(4)]
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
 _p_ds=31;assert _is_prime_det(_p_ds);_c_ds=[1]+[_rng.randrange(_p_ds)for _ in range(15)]
 _B_ds=[[0]*16 for _ in range(16)];_B_ds[0][0]=_p_ds
 for _j in range(1,16):_B_ds[_j][0]=(-_c_ds[_j])%_p_ds;_B_ds[_j][_j]=1
 _Br_ds,_dg_ds,_lm_ds=lll_16d_int(_B_ds,9,10);assert _dg_ds[16]==_p_ds*_p_ds
 _fp_cnt,_fp_nd,_fp_h=_fp_exact_int(_Br_ds,_dg_ds,_lm_ds,4,collect=True)
 _q32_cnt,_q32_nd,_q32_h,_q32_ops=_fp_q32_int(_Br_ds,_dg_ds,_lm_ds,4,collect=True)
 _bf_h=[];_vb_ds=[0]*16
 def _bf_rec(_idx,_rm,_nz):
  if _idx==16:
   if _nz and sum(_c_ds[j]*_vb_ds[j]for j in range(16))%_p_ds==0:_bf_h.append(tuple(_vb_ds))
   return
  _r=math.isqrt(_rm)
  for _x in range(-_r,_r+1):
   if not _nz and _x<0:continue
   _vb_ds[_idx]=_x;_bf_rec(_idx+1,_rm-_x*_x,_nz or(_x>0))
 _bf_rec(0,4,False)
 assert _fp_cnt==_q32_cnt==len(_fp_h)==len(_q32_h)==len(set(_q32_h))==len(_bf_h)>400 and set(_fp_h)==set(_q32_h)==set(_bf_h)
 _bc=[[0]*5 for _ in range(17)];_bc[0]=[1]*5
 for _d in range(1,17):
  for _s in range(5):_bc[_d][_s]=sum(_bc[_d-1][_s-_x*_x]for _x in range(-math.isqrt(_s),math.isqrt(_s)+1))
 # Build flat 128-bit prefix table P_d[d][rem*(2R+1)+(v+R)] matching dsieve16_rank_step's exact 10 opcodes!
 _R_t=2;_W_t=2*_R_t+1;_Pd_flat=[]
 for _d in range(16):
  _tbl=[0]*(5*_W_t)
  for _rm in range(5):
   _r_rm=math.isqrt(_rm);_acc=0
   for _v in range(-_r_rm,_r_rm+1):
    _tbl[_rm*_W_t+(_v+_R_t)]=_acc;_acc+=_bc[15-_d][_rm-_v*_v]
  _Pd_flat.append(_tbl)
 _s16_2=(_bc[16][4]-1)//2;_off2=(_bc[16][4]+1)//2;_rks=set();_M64=(1<<64)-1
 for _hv in _q32_h:
  _x5=_x6=0;_x7=4
  for _d in range(16):
   _x1=_hv[_d];_x2=_x1+_R_t;_x0=_x7*_W_t;_x2=_x0+_x2
   _ent=_Pd_flat[_d][_x2];_x3=_ent&_M64;_x4=_ent>>64
   _s5=_x5+_x3;_x5=_s5&_M64;_x6=_x6+_x4+(_s5>>64);_x7=_x7-_x1*_x1
  _rk=(_x6<<64)|_x5
  _idx=_rk-_off2;assert 0<=_idx<_s16_2;_rks.add(_idx)
 assert len(_rks)==len(_q32_h)
 def _find_prime_Q(_bp,_sd_seed):
  _rp=random.Random(_sd_seed);_p=((1<<(_bp-1))|(_rp.getrandbits(_bp-2)<<1)|1);_p+=(1-_p)%15
  while True:
   if _is_prime_det(_p)and(_p%9 in (4,7)):
    for _g in range(2,50):
     _s=pow(_g,(_p-1)//15,_p)
     if pow(_s,3,_p)!=1 and pow(_s,5,_p)!=1:
      _rhs=(U_bls*pow(_s,10,_p))%_p
      if pow(_rhs,(_p-1)//3,_p)==1:
       _t=pow(_rhs,(2*_p+1)//9 if _p%9==4 else(_p+2)//9,_p)
       if pow(_t,3,_p)==_rhs:return _p,_s,_t
      break
   _p+=15
 def _struct_L(_a1,_b1,_a2,_b2,_p,_s,_t):
  _sp=[pow(_s,j,_p)for j in range(8)]
  _ev=lambda poly:sum(poly[j]*_sp[j]for j in range(8))%_p
  _u1=(_ev(_a1)-_t*_ev(_b1))%_p;_u2=(_ev(_a2)-_t*_ev(_b2))%_p
  _c=[(_u1*_sp[j])%_p for j in range(8)]+[(_u2*_sp[j])%_p for j in range(8)]
  _k0=next(i for i,v in enumerate(_c)if v!=0);_iv=pow(_c[_k0],-1,_p);_c=[(v*_iv)%_p for v in _c]
  _Bb=[[0]*16 for _ in range(16)];_Bb[0][_k0]=_p;_rr=1
  for _j in range(16):
   if _j==_k0:continue
   _Bb[_rr][_k0]=(-_c[_j])%_p;_Bb[_rr][_j]=1;_rr+=1
  return _Bb
 _alpha=100.0/74.0;_vb64=[1.0]+[(math.pi**(d/2.0)/math.gamma(d/2.0+1.0))*(64.05**d)for d in range(1,17)]
 def _b_int_64(_ls):
  _inv=[1.0/x for x in _ls];_dp=[1.0];_tot=0.0
  for _k in range(15,0,-1):
   _u=_inv[_k];_ndp=[0.0]*(len(_dp)+1)
   for _s,_v in enumerate(_dp):_ndp[_s]+=_v;_ndp[_s+1]+=_v*_u
   _dp=_ndp;_tot+=0.5*(sum(_vb64[_s]*_dp[_s]for _s in range(len(_dp)))+1.0)
  return _tot
 _ds_disc={};_gate_ratios=[]
 for _yi,_lbl in ((0,"g"),(12,"f")):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16)
  _Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_qq;_rp=[pow(_ry,i,_qq)for i in range(8)]
  for _i in range(1,8):_Bq[_i][0]=(-_rp[_i])%_qq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=(_rt*_rp[_i])%_qq;_Bq[8+_i][8+_i]=1
  _Rr,_,_=lll_16d_int(_Bq,99,100);_a1,_b1=_Rr[0][:8],_Rr[0][8:]
  _a2,_b2=next((v[:8],v[8:])for v in _Rr[1:]if any(sub_phi15_z(mul_phi15_z(_a1,v[8:]),mul_phi15_z(_b1,v[:8]))))
  _p83,_s83,_t83=_find_prime_Q(84,1084+_yi*10)
  _Br83,_dg83,_lm83=lll_16d_int(_struct_L(_a1,_b1,_a2,_b2,_p83,_s83,_t83),9,10);assert _dg83[16]==_p83*_p83
  _hq83,_ndq83,_,_ops83=_fp_q32_int(_Br83,_dg83,_lm83,64*64,collect=False)
  _ls83=[math.sqrt(_dg83[i+1]/_dg83[i])for i in range(16)];_bi83=_b_int_64(_ls83)
  assert _hq83>400 and _ndq83<=_bi83 and _ops83<=144.0*_bi83
  if _lbl=="g":
   _h83,_nd83,_=_fp_exact_int(_Br83,_dg83,_lm83,64*64,collect=False);assert _h83==_hq83 and _nd83==_ndq83
  _ds_disc[f"{_lbl}_p84_R64"]=dict(hits=_hq83,nodes_q32=_ndq83,ops_sim=_ops83,B_int_ub=round(_bi83,1),actual_over_ub=round(_ndq83/_bi83,4))
  for _bp in (16,48,60,70):
   _pb,_sb,_tb=_find_prime_Q(_bp,2000+_yi*100+_bp)
   _,_dgb,_=lll_16d_int(_struct_L(_a1,_b1,_a2,_b2,_pb,_sb,_tb),9,10);assert _dgb[16]==_pb*_pb
   _ls=[math.sqrt(_dgb[i+1]/_dgb[i])for i in range(16)]
   _cwc=(_pb**(1.0/16.0))*(_alpha**(60.0/16.0));_lswc=[_cwc*(_alpha**(-j/2.0))for j in range(16)]
   _gf=float(desc["params"].get("gate_frac",1.0))if _bp<48 else float(desc["params"].get("gate_frac_hi",1.0));assert _gf==1.0
   _gr=_b_int_64(_ls)/(_gf*_b_int_64(_lswc));assert _gr<0.45;_gate_ratios.append(round(_gr,3))
 _L63=(1<<63)-1
 def _i64(v):assert-_L63-1<=v<=_L63;return v
 def _smh(a,b):_i64(a);_i64(b);return(a*b)>>64
 def _umh(a,b):assert 0<=a<(1<<64)and 0<=b<(1<<64);return(a*b)>>64
 def _bnorm(re,im,rad,sh):
  _i64(re);_i64(im);_i64(rad);assert rad>=0;mx=max(abs(re),abs(im),rad)
  if mx>=(1<<60):
   k=mx.bit_length()-60;mk=(1<<k)-1
   re>>=k;im>>=k;rad=((rad+mk)>>k)+2;sh+=k
  elif mx>0 and mx<(1<<59):
   k=60-mx.bit_length();re<<=k;im<<=k;rad<<=k;sh-=k
  assert abs(re)<(1<<60)and abs(im)<(1<<60)and 0<=rad<(1<<60)
  return re,im,rad,sh
 def _badd(x,y):
  r1,i1,e1,s1=x;r2,i2,e2,s2=y
  if s1<s2:
   d=s2-s1
   if d>60:return _bnorm(r2,i2,e2+2,s2)
   mk=(1<<d)-1;r1>>=d;i1>>=d;e1=((e1+mk)>>d)+2;s1=s2
  elif s2<s1:return _badd(y,x)
  return _bnorm(_i64(r1+r2),_i64(i1+i2),_i64(e1+e2),s1)
 def _bmul(x,y):
  r1,i1,e1,s1=x;r2,i2,e2,s2=y
  re=_i64(_smh(r1,r2)-_smh(i1,i2));im=_i64(_smh(r1,i2)+_smh(i1,r2))
  a1=abs(r1)+abs(i1);a2=abs(r2)+abs(i2)
  rad=_i64(_umh(a1,e2)+_umh(a2,e1)+_umh(e1,e2)+5)
  return _bnorm(re,im,rad,s1+s2+64)
 _ln64_lo=_iv_ln((Decimal(64),Decimal(64)))[0]
 _LUT_HI=[int(_CH.multiply(_CH.divide(_CH.subtract(_iv_ln((Decimal(65+t),Decimal(65+t)))[1],_ln64_lo),_LN2_IV[0]),Decimal(1<<16)).to_integral_value(rounding=ROUND_CEILING))for t in range(64)]
 def _bl2_hi(x):
  re,im,rad,sh=x;r30=re>>30;i30=im>>30;e30=((rad+(1<<30)-1)>>30)+2
  sq=math.isqrt(_i64(_i64(r30*r30)+_i64(i30*i30)))
  m_hi=sq+1+e30;wc=(sq-e30>=64*e30)
  bl=m_hi.bit_length()-1
  t=((m_hi>>(bl-6))&63)if bl>=6 else((m_hi<<(6-bl))&63)
  return((bl+sh+30)<<16)+_LUT_HI[t],wc
 assert U_bls%16==0;_ub=_bnorm(U_bls>>4,0,0,4)
 _w_pre=[_bnorm(_Z_TAB[k][5][0],_Z_TAB[k][5][1],1,-58)for k in range(4)]
 _wp1_pre=[_badd(_w_pre[k],(1<<58,0,0,-58))for k in range(4)]
 _uw2_pre=[_bmul(_ub,_bmul(_w_pre[k],_w_pre[k]))for k in range(4)]
 _w_el=y_pow_z(5);_wp1=add_phi15_z(_w_el,[1]+[0]*7);_w2=y_pow_z(10)
 def _ev_lin(_cs,_k):
  _zt=_Z_TAB[_k];_mx=max(abs(c)for c in _cs)
  if _mx==0:return(0,0,0,0)
  if _mx>=(1<<60):
   _sc=_mx.bit_length()-60;_cs2=[c>>_sc for c in _cs]
   return _bnorm(_i64(sum(_smh(_cs2[j],_zt[j][0])for j in range(8))),_i64(sum(_smh(_cs2[j],_zt[j][1])for j in range(8))),32,_sc+6)
  _sl=60-_mx.bit_length();_cs2=[c<<_sl for c in _cs]
  return _bnorm(_i64(sum(_smh(_cs2[j],_zt[j][0])for j in range(8))),_i64(sum(_smh(_cs2[j],_zt[j][1])for j in range(8))),24,6-_sl)
 _cells_bc=[([_rng.randint(-50,50)for _ in range(8)],[_rng.randint(-50,50)for _ in range(8)])for _ in range(80)]
 for _k in range(30):
  _bb=[_rng.randint(-5,5)for _ in range(8)];_bb[0]=1
  _aa=[int(round(abs(U_bls)**(1.0/3.0)))*_bb[j]for j in range(8)];_cells_bc.append((_aa,_bb))
 for _k in range(40):
  _cells_bc.append(([_rng.randint(-(1<<25),1<<25)for _ in range(8)],[_rng.randint(-(1<<25),1<<25)for _ in range(8)]))
 _cells_bc.append(([1<<55]+[0]*7,[1<<55]+[0]*7))
 _mgf=_mgg=0.0;_nc_g=_nc_f=0
 for _aa,_bb in _cells_bc:
  if not any(_aa)or not any(_bb):continue
  _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
  _ng_ex=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
  _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
  _nf_ex=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
  _bg_q=_bf_q=0;_bg_wc=_bf_wc=True
  for _k in range(4):
   _sa=_ev_lin(_aa,_k);_sb=_ev_lin(_bb,_k)
   _a3b=_bmul(_bmul(_sa,_sa),_sa);_b3b=_bmul(_bmul(_sb,_sb),_sb)
   _tgb=_bmul(_uw2_pre[_k],_b3b);_gv=_badd(_a3b,(-_tgb[0],-_tgb[1],_tgb[2],_tgb[3]))
   _t1=_badd(_a3b,_bmul(_wp1_pre[_k],_b3b));_b6b=_bmul(_b3b,_b3b)
   _t2=_badd(_bmul(_t1,_a3b),_b6b);_b9b=_bmul(_b6b,_b3b);_tfb=_bmul(_w_pre[_k],_b9b)
   _fv=_badd(_bmul(_t2,_a3b),(-_tfb[0],-_tfb[1],_tfb[2],_tfb[3]))
   _gh,_gwc=_bl2_hi(_gv);_fh,_fwc=_bl2_hi(_fv)
   _bg_q+=2*_gh;_bf_q+=2*_fh
   if not _gwc:_bg_wc=False
   if not _fwc:_bf_wc=False
  if not _bg_wc:_nc_g+=1
  if not _bf_wc:_nc_f+=1
  _bf=_bf_q/65536.0;_bg=_bg_q/65536.0
  if _ng_ex>0:
   assert _bg>=math.log2(_ng_ex)
   if _bg_wc:assert _bg-math.log2(_ng_ex)<=0.6212
   _mgg=max(_mgg,_bg-math.log2(_ng_ex))
  if _nf_ex>0:
   assert _bf>=math.log2(_nf_ex)
   if _bf_wc:assert _bf-math.log2(_nf_ex)<=0.6212
   _mgf=max(_mgf,_bf-math.log2(_nf_ex))
 assert _mgf<0.25 and _mgg<0.25 and _nc_g==0 and _nc_f==0
 _yr_min=1e99;_eff_u=[];_ss_yr=dlive["_ss_eval"]
 for _yi in (0,5,12,17):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_sn="g" if _sd==0 else "f"
  _qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16);_nb=_qq.bit_length()-1
  _prw=next(r for r in dlive["rows"]if r["side"]==_sn and r["lo"]<=_nb<r["hi"])
  _Rp=_prw["R"];_Rp2=_Rp*_Rp;_Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_qq;_rp=[pow(_ry,i,_qq)for i in range(8)]
  for _i in range(1,8):_Bq[_i][0]=(-_rp[_i])%_qq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=(_rt*_rp[_i])%_qq;_Bq[8+_i][8+_i]=1
  _Rr,_,_=lll_16d_int(_Bq,99,100);_a1,_b1=_Rr[0][:8],_Rr[0][8:]
  _a2,_b2=next((v[:8],v[8:])for v in _Rr[1:]if any(sub_phi15_z(mul_phi15_z(_a1,v[8:]),mul_phi15_z(_b1,v[:8]))))
  _nfs,_ngs,_wts=[],[],[]
  for _ in range(12):
   while True:
    _g=[_rng.gauss(0,1)for _ in range(16)];_rn=math.sqrt(sum(x*x for x in _g))
    _rd=_Rp*(_rng.random()**(1.0/16.0))/_rn;_vv=[int(round(x*_rd))for x in _g]
    if 0<sum(x*x for x in _vv)<=_Rp2:break
   _al,_bt=_vv[:8],_vv[8:]
   if not _is_cp(_al,_bt):_wts.append(0.0);_nfs.append(9999.0);_ngs.append(9999.0);continue
   _wts.append(1.0/_cnt_orb(_al,_bt,_Rp2))
   _aa=add_phi15_z(mul_phi15_z(_al,_a1),mul_phi15_z(_bt,_a2));_bb=add_phi15_z(mul_phi15_z(_al,_b1),mul_phi15_z(_bt,_b2))
   _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
   _ng=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
   _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
   _nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
   _nfs.append(math.log2(_nf)-(math.log2(_qq)if _sn=="f" else 0.0));_ngs.append(math.log2(_ng)-(math.log2(_qq)if _sn=="g" else 0.0))
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
  cnts=[r[4]if grp=="enc" else r[3]for r in rows]
  mus=[r[5]if grp=="enc" else r[4]for r in rows]
  mc=sum(cnts)/40.0;vc=sum((x-mc)**2 for x in cnts)/39.0;fano=vc/mc;k_hat=max(0.0,(vc-mc)/(mc*mc))
  emp0=cnts.count(0)/40.0;mod0=sum(1.0/(1.0+0.5*m)for m in mus)/40.0
  wpgf=max((sum(t**c for c in cnts)/40.0)/(sum(1.0/(1.0+0.5*m*(1.0-t))for m in mus)/40.0)for t in [i*0.05 for i in range(17)])
  assert fano<1.5 and k_hat<0.25 and emp0<mod0 and wpgf<0.90
  out[grp]=dict(n=40,mean_M=round(mc,2),mean_mu_raw=round(sum(mus)/40.0,2),fano=round(fano,3),k_hat=round(k_hat,4),zeros=cnts.count(0),emp0=round(emp0,4),h5_empty=round(mod0,4),worst_pgf=round(wpgf,4))
 if "hi100" in camp93:
  h1=camp93["hi100"];assert len(h1)==100 and min(h1)>=2 and h1.count(0)==0
  for lbl,sl in (("enc40",h1[:40]),("g30",h1[40:70]),("f30",h1[70:])):
   ms=sum(sl)/len(sl);vs=sum((x-ms)**2 for x in sl)/(len(sl)-1);fn=vs/ms;kh=max(0.0,(vs-ms)/(ms*ms))
   assert kh<0.25 and fn<1.0+0.5*ms
   out[f"hi_{lbl}"]=dict(n=len(sl),mean_M=round(ms,2),min_M=min(sl),fano=round(fn,3),k_hat=round(kh,4),zeros=0,emp0=0.0)
  out["hi100_zeros"]="0/100"
 P_SM=math.prod(primes_upto(256));P_MD=math.prod(p for p in primes_upto(8192)if p>256);NA=1<<24;NA3=NA**3
 pts=[((i,1),)for i in range(16)]
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
  B=[[0]*16 for _ in range(16)];B[0][0]=p;rp=[pow(ry,i,p)for i in range(8)]
  for i in range(1,8):B[i][0]=(-rp[i])%p;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*rp[i])%p;B[8+i][8+i]=1
  Rr,_,_=lll_16d_int(B,99,100);a1,b1=Rr[0][:8],Rr[0][8:]
  a2,b2=next((v[:8],v[8:])for v in Rr[1:]if any(sub_phi15_z(mul_phi15_z(a1,v[8:]),mul_phi15_z(b1,v[:8]))))
  seen=set()
  ba=[mul_phi15_z(y_pow_z(i),a1)for i in range(8)]+[mul_phi15_z(y_pow_z(i),a2)for i in range(8)]
  bb_t=[mul_phi15_z(y_pow_z(i),b1)for i in range(8)]+[mul_phi15_z(y_pow_z(i),b2)for i in range(8)]
  for nz in pts:
   if len(nz)==1:i0,s0=nz[0];aa=ba[i0];bb=bb_t[i0]
   elif len(nz)==2:
    (i0,s0),(i1,s1)=nz;a0,a1=ba[i0],ba[i1];b0,b1=bb_t[i0],bb_t[i1]
    aa=[a0[t]+a1[t]if s1>0 else a0[t]-a1[t]for t in range(8)];bb=[b0[t]+b1[t]if s1>0 else b0[t]-b1[t]for t in range(8)]
   else:
    (i0,s0),(i1,s1),(i2,s2)=nz;a0,a1,a2=ba[i0],ba[i1],ba[i2];b0,b1,b2=bb_t[i0],bb_t[i1],bb_t[i2]
    aa=[a0[t]+s1*a1[t]+s2*a2[t]for t in range(8)];bb=[b0[t]+s1*b1[t]+s2*b2[t]for t in range(8)]
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
def verify_q15_adaptive_h5_live(ad_dat,primes_upto,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 from collections import Counter
 U_TOY=-5;W_EL=y_pow_z(5);WP1=add_phi15_z(W_EL,[1]+[0]*7);W2=y_pow_z(10)
 cal=ad_dat["cal_tbl"];nb_ex=ad_dat["N_bin_exact"]
 d1=ad_dat["d1_rows"];d2=ad_dat["d2_rows"]
 assert len(d1)==286 and len(d2)==124
 def _stat(rows):
  n=len(rows);of=sum(1 for _,_,c in rows if c==0)
  ep0=sum(cal[f"{sd}:{p.bit_length()-1}"]["p0_mod"]for sd,p,_ in rows)
  ecf=sum(cal[f"{sd}:{p.bit_length()-1}"]["fail"]for sd,p,_ in rows)
  om=sum(c for _,_,c in rows)/n
  ecm=sum(cal[f"{sd}:{p.bit_length()-1}"]["mean"]for sd,p,_ in rows)/n
  by={}
  for sd,p,_ in rows:by.setdefault((sd,p.bit_length()-1),[]).append(p)
  tp=0;tep=0.0
  for(sd,b),ps in by.items():
   if len(ps)<2:continue
   cc=Counter(ps);tp+=sum(v*(v-1)//2 for v in cc.values())
   tep+=len(ps)*(len(ps)-1)/(2.0*nb_ex[f"{sd}:{b}"])
  cb=tp/max(1e-9,tep)
  oe_p0=of/ep0;mr=om/ecm
  assert oe_p0<0.65 and mr>0.85 and cb<1.5
  return dict(n=n,uniq=len({(sd,p)for sd,p,_ in rows}),fail=of,oe_p0=round(oe_p0,4),reached_over_fresh_fail=round(of/max(1e-9,ecf),4),reached_over_fresh_mean=round(mr,4),c_bday=round(cb,3))
 st=dict(d1_all=_stat(d1),d1_g=_stat([r for r in d1 if r[0]==0]),d1_f=_stat([r for r in d1 if r[0]==1]),d2_all=_stat(d2),pooled=_stat(d1+d2))
 # Live bit-exact replay of ALL 3 strictly-decreasing zero-recount (root -> d1 -> d2) multi-level witnesses!
 P_SM=math.prod(primes_upto(256));P_MD=math.prod(p for p in primes_upto(8192)if p>256);NA=1<<24;NA3=NA**3
 _rng=random.Random(19)
 def _rho(n):
  if n%2==0:return 2
  while True:
   c=_rng.randrange(1,n);y=_rng.randrange(2,n);r=1;q=1;g=1
   while g==1:
    x=y
    for _ in range(r):y=(y*y+c)%n
    k=0
    while k<r and g==1:
     ys=y;m=min(32,r-k)
     for _ in range(m):y=(y*y+c)%n;q=(q*abs(x-y))%n
     g=math.gcd(q,n);k+=m
    r<<=1
   if g==n:
    while True:
     ys=(ys*ys+c)%n;g=math.gcd(abs(x-ys),n)
     if g>1:break
   if g!=n:return g
 def _ok(n):
  if n<=0:return False
  for P in (P_SM,P_MD):
   while True:
    g=math.gcd(n,P)
    if g==1:break
    n//=g
   if n==1:return True
  if n>NA3:return False
  if _is_prime_det(n):return n<=NA
  st_f=[n];fc=0
  while st_f:
   m=st_f.pop()
   if m==1:continue
   if _is_prime_det(m):
    if m>NA:return False
    fc+=1
    if fc>3:return False
   else:d=_rho(m);st_f.append(d);st_f.append(m//d)
  return True
 def _eval_pt(ide,pt):
  sd,q,ry,rt=ide
  B=[[0]*16 for _ in range(16)];B[0][0]=q;rp=[pow(ry,i,q)for i in range(8)]
  for i in range(1,8):B[i][0]=(-rp[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*rp[i])%q;B[8+i][8+i]=1
  Rr,dsq,_=lll_16d_int(B,99,100);assert dsq[16]==q*q
  a1,b1=Rr[0][:8],Rr[0][8:]
  a2,b2=next((v[:8],v[8:])for v in Rr[1:]if any(sub_phi15_z(mul_phi15_z(a1,v[8:]),mul_phi15_z(b1,v[:8])))and abs(norm_q15_z(sub_phi15_z(mul_phi15_z(a1,v[8:]),mul_phi15_z(b1,v[:8]))))==q)
  ba=[mul_phi15_z(y_pow_z(i),a1)for i in range(8)]+[mul_phi15_z(y_pow_z(i),a2)for i in range(8)]
  bb_t=[mul_phi15_z(y_pow_z(i),b1)for i in range(8)]+[mul_phi15_z(y_pow_z(i),b2)for i in range(8)]
  aa=[0]*8;bb=[0]*8
  for idx,sgn in pt:
   for t in range(8):aa[t]+=sgn*ba[idx][t];bb[t]+=sgn*bb_t[idx][t]
  a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3)
  ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_TOY,mul_phi15_z(W2,b3)))))
  a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
  nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_EL,b9))))
  assert(ng if sd==0 else nf)%q==0
  ng_q=ng//q if sd==0 else ng;nf_q=nf//q if sd==1 else nf
  assert _ok(ng_q)and _ok(nf_q)
  return aa,bb,nf,ng,nf_q,ng_q
 wits=ad_dat["witnesses"];assert len(wits)==3
 for w in wits:
  r0,d1i,d2i=w["root"],w["d1"],w["d2"]
  assert r0[1]>d1i[1]>d2i[1]>2048
  a0,b0,nf0,ng0,nf0q,ng0q=_eval_pt(r0,w["pt0"])
  assert(ng0q if d1i[0]==0 else nf0q)%d1i[1]==0
  sp1=[pow(d1i[2],j,d1i[1])for j in range(8)]
  assert(sum(a0[j]*sp1[j]for j in range(8))-d1i[3]*sum(b0[j]*sp1[j]for j in range(8)))%d1i[1]==0
  a1,b1,nf1,ng1,nf1q,ng1q=_eval_pt(d1i,w["pt1"])
  assert(nf1,ng1)!=(nf0,ng0)and(ng1q if d2i[0]==0 else nf1q)%d2i[1]==0
  sp2=[pow(d2i[2],j,d2i[1])for j in range(8)]
  assert(sum(a1[j]*sp2[j]for j in range(8))-d2i[3]*sum(b1[j]*sp2[j]for j in range(8)))%d2i[1]==0
  _,_,nf2,ng2,_,_=_eval_pt(d2i,w["pt2"])
  assert(nf2,ng2)not in ((nf0,ng0),(nf1,ng1))
 return dict(n_roots=ad_dat["n_roots"],r_ok=ad_dat["r_ok"],n_cal=sum(v["n"]for v in cal.values()),strata=st,live_depth2_witnesses=len(wits),checks_ok=True)
_SC_CACHE={}
def eval_scert_discrete_int(sid,q,ry,rt,a1,b1,a2,b2,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,vor_units,n_samps=4096,seed=None):
 import bisect,unitdomain as _UD
 from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
 from operator import mul as _omul
 if "init" not in _SC_CACHE:
  Q58=Decimal(1<<58)
  _ZT=[[(int((_UD._PL[k][j][0]*Q58).to_integral_value()),int((_UD._PL[k][j][1]*Q58).to_integral_value()))for j in range(8)]for k in range(4)]
  _Z0R=tuple(_ZT[0][j][0]for j in range(8));_Z0I=tuple(_ZT[0][j][1]for j in range(8))
  _VC=[tuple(int(((_UD.sq_norms(u)[k]-Decimal(1))*Decimal(1<<16)).to_integral_value(rounding=ROUND_FLOOR))-1 for k in range(4))for u in vor_units]
  _ln2_50=_UD._d_ln(Decimal(2))
  _LH=[0]*128+[int((Decimal(64)*_UD._d_ln(Decimal(idx+1)/Decimal(128))/_ln2_50).to_integral_value(rounding=ROUND_CEILING))for idx in range(128,256)]
  _LH64=[int((Decimal(1<<16)*_UD._d_ln(Decimal(65+t)/Decimal(64))/_ln2_50).to_integral_value(rounding=ROUND_CEILING))for t in range(64)]
  _G4=[y_pow_fn((j*4)%15)for j in range(8)];_G7=[y_pow_fn((j*7)%15)for j in range(8)]
  _SC_CACHE["init"]=(_ZT,_Z0R,_Z0I,_VC,_LH,_LH64,_G4,_G7,y_pow_fn(5),add_fn(y_pow_fn(5),[1]+[0]*7),y_pow_fn(10))
 _ZT,_Z0R,_Z0I,_VC,_LH,_LH64,_G4,_G7,w_el,wp1_el,w2_el=_SC_CACHE["init"]
 _M64=(1<<64)-1
 def _bn(re,im,rad,sh):
  mx=max(abs(re),abs(im),rad)
  if mx>=(1<<60):k=mx.bit_length()-60;mk=(1<<k)-1;re>>=k;im>>=k;rad=((rad+mk)>>k)+2;sh+=k
  elif mx>0 and mx<(1<<59):k=60-mx.bit_length();re<<=k;im<<=k;rad<<=k;sh-=k
  return re,im,rad,sh
 def _ba(x,y):
  r1,i1,e1,s1=x;r2,i2,e2,s2=y
  if s1<s2:
   d=s2-s1
   if d>60:return _bn(r2,i2,e2+2,s2)
   mk=(1<<d)-1;r1>>=d;i1>>=d;e1=((e1+mk)>>d)+2;s1=s2
  elif s2<s1:return _ba(y,x)
  return _bn(r1+r2,i1+i2,e1+e2,s1)
 def _bm(x,y):
  r1,i1,e1,s1=x;r2,i2,e2,s2=y
  re=((r1*r2)>>64)-((i1*i2)>>64);im=((r1*i2)>>64)+((i1*r2)>>64)
  rad=(((abs(r1)+abs(i1))*e2)>>64)+(((abs(r2)+abs(i2))*e1)>>64)+((e1*e2)>>64)+5
  return _bn(re,im,rad,s1+s2+64)
 def _el(cs,k):
  mx=max(abs(c)for c in cs);zt=_ZT[k]
  if mx==0:return(0,0,0,0)
  if mx>=(1<<60):sc=mx.bit_length()-60;c2=[c>>sc for c in cs];return _bn(sum((c2[j]*zt[j][0])>>64 for j in range(8)),sum((c2[j]*zt[j][1])>>64 for j in range(8)),32,sc+6)
  sl=60-mx.bit_length();c2=[c<<sl for c in cs];return _bn(sum((c2[j]*zt[j][0])>>64 for j in range(8)),sum((c2[j]*zt[j][1])>>64 for j in range(8)),24,6-sl)
 def _bl2(x):
  re,im,rad,sh=x;r30=re>>30;i30=im>>30;e30=((rad+(1<<30)-1)>>30)+2
  mh=math.isqrt(r30*r30+i30*i30)+1+e30;bl=mh.bit_length()-1
  t=((mh>>(bl-6))&63)if bl>=6 else((mh<<(6-bl))&63)
  return((bl+sh+30)<<16)+_LH64[t]
 if "uw2" not in _SC_CACHE:
  _ub=_bn(U_bls>>4,0,0,4);_wp=[_bn(_ZT[k][5][0],_ZT[k][5][1],1,-58)for k in range(4)]
  _SC_CACHE["uw2"]=[_bm(_ub,_bm(_wp[k],_wp[k]))for k in range(4)]
 _uw2_p=_SC_CACHE["uw2"]
 sname="g" if sid==0 else "f";nb=q.bit_length()-1;lq=math.log2(q)
 prow=next(r for r in dlive["rows"]if r["side"]==sname and r["lo"]<=nb<r["hi"])
 R_n=int(prow["R"]);cells_b=float(desc["cells_rank2"][str(R_n)])
 n1=float(prow["n"]);kA=int(prow["kA"]);kB=int(prow["kB"]);iA=int(prow["iA"])
 factor=prow["mu"]/(cells_b*float(desc["params"]["udup_node"])*prow["pr"])
 _ss=dlive["_ss_eval"];_mn=dlive["mesh_node"];_af=float(desc["params"]["alpha_f"]);_ag=float(desc["params"]["alpha_g"])
 _c_ref=(prow["lo"]+math.log2(cells_b))/16.0
 _nf_ref=72.0*_c_ref+desc["qbins_lo"][0][0]-(prow["hi"]if sname=="f" else 0.0)
 _ng_ref=24.0*_c_ref+float(desc["params"]["log2_u8"])+desc["qbins_lo"][0][1]-(prow["hi"]if sname=="g" else 0.0)
 _sf_ref=max(1e-300,_ss([_nf_ref],n1,iA,_mn,_af,kA,kB)[0]);_sg_ref=max(1e-300,_ss([_ng_ref],n1,iA,_mn,_ag,kA,kB)[0])
 _Ef=max(44,int(math.floor(-math.log2(_sf_ref)))+44);_Eg=max(44,int(math.floor(-math.log2(_sg_ref)))+44)
 _inv_sc=2.0**(64-_Ef-_Eg);scale_y=cells_b*_inv_sc*factor
 ba=[mul_fn(y_pow_fn(i),a1)for i in range(8)]+[mul_fn(y_pow_fn(i),a2)for i in range(8)]
 bb=[mul_fn(y_pow_fn(i),b1)for i in range(8)]+[mul_fn(y_pow_fn(i),b2)for i in range(8)]
 ba_T=[tuple(ba[j][t]for j in range(16))for t in range(8)]
 bb_T=[tuple(bb[j][t]for j in range(16))for t in range(8)]
 ryp=[pow(ry,i,q)for i in range(8)]
 b_ry_vec=tuple(sum(bb[j][t]*ryp[t]for t in range(8))%q for j in range(16))
 lq_f=lq if sname=="f" else 0.0;lq_g=lq if sname=="g" else 0.0
 if R_n not in _SC_CACHE:
  R2=R_n*R_n
  with localcontext()as _ctx:
   _ctx.prec=60
   inv2=Decimal(1)/(Decimal(2)*(Decimal(R_n)/Decimal(4))**2);S64=Decimal(1<<64)
   xs=list(range(-R_n,R_n+1));rw=[(-Decimal(x*x)*inv2).exp()for x in xs];Z_=sum(rw)
   ws=[int((w*S64/Z_).to_integral_value())for w in rw];ws[R_n]+=(1<<64)-sum(ws)
   cdf=[];ac=0
   for w in ws:ac+=w;cdf.append(ac)
   lt=[bisect.bisect_right(cdf,b<<56)for b in range(256)]
   ht=[bisect.bisect_right(cdf,((b+1)<<56)-1)+1 for b in range(256)]
   step=(-inv2).exp();pacc=[0]*(R2+1);cur=Decimal(1)
   for s2 in range(R2,-1,-1):pacc[s2]=int((cur*S64).to_integral_value());cur*=step
  _SC_CACHE[R_n]=(R2,xs,cdf,lt,ht,pacc)
 R2,xs,cdf,lt,ht,pacc=_SC_CACHE[R_n]
 rng=random.Random((20261125+(q&0xffffff))if seed is None else seed)
 bl_f=bisect.bisect_right;gb=rng.getrandbits
 n_nf_chk=0;max_nf_gap=-1e9;got=0;acc_u128=0;k5000=0;s_sum=0.0;n_cp=0
 while got<n_samps:
  v=[]
  for _ in range(16):
   r_=gb(64);b_=r_>>56;v.append(xs[bl_f(cdf,r_,lt[b_],ht[b_])])
  s2=sum(map(_omul,v,v))
  if not(0<s2<=R2 and gb(64)<pacc[s2]):continue
  got+=1
  c1,c2=v[:8],v[8:]
  rc0=sum(((c1[j]<<45)*_Z0R[j])>>64 for j in range(8))>>24
  ic0=sum(((c1[j]<<45)*_Z0I[j])>>64 for j in range(8))>>24
  if rc0<0:rc0=-rc0;ic0=-ic0
  if rc0<=2 or ic0<=2 or 10000*(ic0+2)>=2125*(rc0-2):continue
  C_lo=[0]*4;C_hi=[0]*4
  for k in range(4):
   ztk=_ZT[k]
   rk=rc0 if k==0 else(sum(((c1[j]<<45)*ztk[j][0])>>64 for j in range(8))>>24)
   ik=ic0 if k==0 else(sum(((c1[j]<<45)*ztk[j][1])>>64 for j in range(8))>>24)
   ar,ai=abs(rk),abs(ik);rlo=ar-2 if ar>2 else 0;ilo=ai-2 if ai>2 else 0
   C_lo[k]=(rlo*rlo+ilo*ilo)>>10;C_hi[k]=((ar+2)**2+(ai+2)**2+1023)>>10
  if not all((row[0]*(C_lo[0]if row[0]>=0 else C_hi[0])+row[1]*(C_lo[1]if row[1]>=0 else C_hi[1])+row[2]*(C_lo[2]if row[2]>=0 else C_hi[2])+row[3]*(C_lo[3]if row[3]>=0 else C_hi[3]))>0 for row in _VC):continue
  c14=[sum(c1[j]*_G4[j][i]for j in range(8))for i in range(8)];d1=mul_fn(c1,c14)
  d17=[sum(d1[j]*_G7[j][i]for j in range(8))for i in range(8)];e1=mul_fn(d1,d17)
  nc1=e1[0]*e1[0]-e1[0]*e1[5]+e1[5]*e1[5]
  c24=[sum(c2[j]*_G4[j][i]for j in range(8))for i in range(8)];d2=mul_fn(c2,c24)
  d27=[sum(d2[j]*_G7[j][i]for j in range(8))for i in range(8)];e2=mul_fn(d2,d27)
  nc2=e2[0]*e2[0]-e2[0]*e2[5]+e2[5]*e2[5]
  if not nc1 or not nc2:continue
  ga,gb_=nc1,nc2;st=0
  while gb_ and st<48:ga,gb_=gb_,ga%gb_;st+=1
  if gb_ or ga!=1:continue
  if sum(map(_omul,v,b_ry_vec))%q==0:continue
  a=[sum(map(_omul,v,ba_T[t]))for t in range(8)];b=[sum(map(_omul,v,bb_T[t]))for t in range(8)]
  mx=max(max(map(abs,a)),max(map(abs,b)))
  if not mx:continue
  lsh=55-mx.bit_length()
  if lsh>=0:a_s=[x<<lsh for x in a];b_s=[x<<lsh for x in b]
  else:rsh=-lsh;a_s=[x>>rsh for x in a];b_s=[x>>rsh for x in b]
  x_lf=0;ok_e=True
  for k in range(4):
   ztk=_ZT[k]
   ra30=sum((a_s[j]*ztk[j][0])>>64 for j in range(8))>>24;ia30=sum((a_s[j]*ztk[j][1])>>64 for j in range(8))>>24
   rb30=sum((b_s[j]*ztk[j][0])>>64 for j in range(8))>>24;ib30=sum((b_s[j]*ztk[j][1])>>64 for j in range(8))>>24
   Sk=(abs(ra30)+2)**2+(abs(ia30)+2)**2+(abs(rb30)+2)**2+(abs(ib30)+2)**2
   if Sk<=0:ok_e=False;break
   clz=64-Sk.bit_length();x_lf+=_LH[(Sk<<clz)>>56]+((63-clz)<<6)
  if not ok_e:continue
  _lf_q=max(0,math.ceil((9.0*(x_lf/64.0+8.0*(30-lsh))-lq_f)*4.0))
  _bg_q=0
  for k in range(4):
   _sa=_el(a,k);_sb=_el(b,k);_a3b=_bm(_bm(_sa,_sa),_sa);_b3b=_bm(_bm(_sb,_sb),_sb)
   _tgb=_bm(_uw2_p[k],_b3b);_gv=_ba(_a3b,(-_tgb[0],-_tgb[1],_tgb[2],_tgb[3]))
   _bg_q+=2*_bl2(_gv)
  _lg_q=max(0,math.ceil((_bg_q/65536.0-lq_g)*4.0))
  L_cf=_lf_q/4.0;L_hg=_lg_q/4.0
  if n_nf_chk<3:
   a3,b3=pow_fn(a,3),pow_fn(b,3)
   ng_ex=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
   assert ng_ex>0 and (math.log2(ng_ex)-lq_g)<=L_hg
   a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
   nf_ex=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
   gap=(math.log2(nf_ex)-lq_f)-L_cf;assert gap<=0.0
   if gap>max_nf_gap:max_nf_gap=gap
   n_nf_chk+=1
  sf=_ss([L_cf],n1,iA,_mn,_af,kA,kB)[0];sg=_ss([L_hg],n1,iA,_mn,_ag,kA,kB)[0]
  Tf=min(_M64,int(math.floor(sf*(2.0**_Ef))));Tg=min(_M64,int(math.floor(sg*(2.0**_Eg))))
  uh=(Tf*Tg)>>64;acc_u128+=uh;y=uh*scale_y;n_cp+=1;s_sum+=y
  if y>=5000.0:k5000+=1
 assert acc_u128<(1<<128)
 return dict(s_hat=s_sum/n_samps,k5000=k5000,n_cp=n_cp,max_nf_gap=max_nf_gap,n_nf_chk=n_nf_chk)
def verify_prod_cond_q15(enc12,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,lll_16d_int,vor_units=None,eta4k_raw=None,scert4k=None,pgcd_mod=None):
 import obpair as _OP,unitdomain as _UD
 from decimal import Decimal,localcontext,ROUND_FLOOR
 assert len(enc12)==12 and eta4k_raw is not None and scert4k is not None and len(scert4k)==4000
 ud_chk=_UD.check(vor_units,mul_fn,y_pow_fn);ud_fg=_UD.full_group(vor_units,mul_fn,y_pow_fn,sub_fn,add_fn,norm_fn)
 assert ud_chk["min_ray_coord"]>0.028 and ud_fg==9
 Q58=Decimal(1<<58)
 ZT=[[(int((_UD._PL[k][j][0]*Q58).to_integral_value()),int((_UD._PL[k][j][1]*Q58).to_integral_value()))for j in range(8)]for k in range(4)]
 ZTR=[tuple(ZT[k][j][0]for j in range(8))for k in range(4)];ZTI=[tuple(ZT[k][j][1]for j in range(8))for k in range(4)]
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 log2_u8=float(desc["params"]["log2_u8"]);S8=math.log2(desc["cells_rank2"]["8"])
 qb_f=[x[1]for x in desc["qbins"]];qb_g=[x[2]for x in desc["qbins"]]
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4)for _ in range(8)];bt=[rng.randint(-4,4)for _ in range(8)]
  if any(al)and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
 cfs,cgs,pr_rats,sc_live,k5_enc=[],[],[],[],[]
 nf_chk_tot=0;max_nf_gap=-1e9;enc_eta_ok=0
 _ss=dlive["_ss_eval"];_mn=dlive["mesh_node"];_af=float(desc["params"]["alpha_f"]);_ag=float(desc["params"]["alpha_g"])
 for idx_e,(sid,q_hex,pa,pb,ry_hex,rt_hex) in enumerate(enc12):
  q,ry,rt=int(q_hex,16),int(ry_hex,16),int(rt_hex,16)
  pa3,pb3=pow_fn(pa,3),pow_fn(pb,3)
  if sid==0:pn=abs(norm_fn(sub_fn(pa3,scale_fn(U_bls,mul_fn(w2_el,pb3)))))
  else:
   pa6,pb6=mul_fn(pa3,pa3),mul_fn(pb3,pb3);pa9,pb9=mul_fn(pa6,pa3),mul_fn(pb6,pb3)
   pn=abs(norm_fn(sub_fn(add_fn(add_fn(pa9,mul_fn(wp1_el,mul_fn(pa6,pb3))),mul_fn(pa3,pb6)),mul_fn(w_el,pb9))))
  assert pn>0 and pn%q==0
  ry_p=[pow(ry,i,q)for i in range(8)]
  assert(sum(pa[i]*ry_p[i]for i in range(8))-rt*sum(pb[i]*ry_p[i]for i in range(8)))%q==0
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  for i in range(1,8):B[i][0]=(-ry_p[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*ry_p[i])%q;B[8+i][8+i]=1
  R,dsq,lsq=lll_16d_int(B,99,100);assert dsq[16]==q*q
  for kk in range(1,16):assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  pair=_OP.ob_pair(q,ry,rt,R);assert pair is not None
  a1,b1,a2,b2=pair
  assert abs(norm_fn(sub_fn(mul_fn(a1,b2),mul_fn(b1,a2))))==q
  pM=1
  for k in range(4):
   zr,zi=ZTR[k],ZTI[k];mk=0
   for v in (a1,b1,a2,b2):
    re=v[0]*zr[0]+v[1]*zr[1]+v[2]*zr[2]+v[3]*zr[3]+v[4]*zr[4]+v[5]*zr[5]+v[6]*zr[6]+v[7]*zr[7]
    im=v[0]*zi[0]+v[1]*zi[1]+v[2]*zi[2]+v[3]*zi[3]+v[4]*zi[4]+v[5]*zi[5]+v[6]*zi[6]+v[7]*zi[7]
    er=abs(v[0])+abs(v[1])+abs(v[2])+abs(v[3])+abs(v[4])+abs(v[5])+abs(v[6])+abs(v[7])
    mk+=(abs(re)+er)**2+(abs(im)+er)**2
   pM*=mk
  assert pM*pM<=(602248*q)<<928;enc_eta_ok+=1
  sname="g" if sid==0 else "f";nb=q.bit_length()-1;lq=math.log2(q)
  prow=next(r for r in dlive["rows"]if r["side"]==sname and r["lo"]<=nb<r["hi"])
  R_n=prow["R"];cells_b=int(desc["cells_rank2"][str(R_n)])
  S_bin=math.log2(cells_b);c_act=(lq+S_bin)/16.0
  c_val=(lq+S8)/16.0;nfas,ngas=[],[]
  for al,bt in ab_pts:
   a=add_fn(mul_fn(al,a1),mul_fn(bt,a2));b=add_fn(mul_fn(al,b1),mul_fn(bt,b2))
   a3,b3=pow_fn(a,3),pow_fn(b,3)
   ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
   a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
   nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
   assert(ng if sid==0 else nf)%q==0
   cf=math.log2(nf)-72.0*c_val;cg=math.log2(ng)-24.0*c_val-log2_u8
   cfs.append(cf);cgs.append(cg)
   nfas.append(72.0*c_act+cf-(lq if sname=="f" else 0.0));ngas.append(24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0))
  _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,_af,prow["kA"],prow["kB"])
  _pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,_ag,prow["kA"],prow["kB"])
  pr_rats.append((sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts))/prow["pr"])
  if idx_e%3==0:
   sr=eval_scert_discrete_int(sid,q,ry,rt,a1,b1,a2,b2,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,vor_units,1024)
   sc_live.append(sr["s_hat"]);k5_enc.append(sr["k5000"]);nf_chk_tot+=sr["n_nf_chk"]
   if sr["max_nf_gap"]>max_nf_gap:max_nf_gap=sr["max_nf_gap"]
 cfs.sort();cgs.sort();n=len(cfs)
 sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
 sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
 assert all(sf6[i]<qb_f[i]for i in range(6))and all(sg6[i]<qb_g[i]for i in range(6))and min(pr_rats)>=10.0 and enc_eta_ok==12
 assert min(sc_live)>=1000.0 and min(k5_enc)>=10 and max_nf_gap<=0.0
 G4=[y_pow_fn((j*4)%15)for j in range(8)];G7=[y_pow_fn((j*7)%15)for j in range(8)]
 def _is_p3(n):
  d,s=n-1,0
  while not(d&1):d>>=1;s+=1
  for a in (2,):
   x=pow(a,d,n)
   if x in (1,n-1):continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 def _n_ex(c):
  c4=[sum(c[j]*G4[j][i]for j in range(8))for i in range(8)];d=mul_fn(c,c4)
  d7=[sum(d[j]*G7[j][i]for j in range(8))for i in range(8)];e=mul_fn(d,d7)
  return e[0]*e[0]-e[0]*e[5]+e[5]*e[5]
 pos=0;n_4k=0;qs_4k=[]
 while pos<len(eta4k_raw):
  sd=eta4k_raw[pos];pos+=1;vs=[]
  for _ in range(32):
   zz=0;sh=0
   while True:
    b_=eta4k_raw[pos];pos+=1;zz|=(b_&127)<<sh
    if not(b_&128):break
    sh+=7
   vs.append((zz>>1)^-(zz&1))
  a1,b1,a2,b2=vs[:8],vs[8:16],vs[16:24],vs[24:32]
  det_el=sub_fn(mul_fn(a1,b2),mul_fn(b1,a2));q=_n_ex(det_el)
  assert q>(1<<70)and q%15==1 and _is_p3(q)
  _, _, g_use, kr_use = scert4k[n_4k]
  ry=pow(g_use,((q-1)//15)*kr_use,q);r2_=(ry*ry)%q;r3_=(r2_*ry)%q;assert r3_!=1 and(r3_*r2_)%q!=1
  vv=0;pw=1
  for cc in det_el:vv=(vv+cc*pw)%q;pw=(pw*ry)%q
  assert vv==0
  rp=[1]*8
  for i in range(1,8):rp[i]=(rp[i-1]*ry)%q
  w_=rp[5];w2_=(w_*w_)%q;nzb=False
  for aa,bb in ((a1,b1),(a2,b2)):
   A_=sum(aa[i]*rp[i]for i in range(8))%q;B_=sum(bb[i]*rp[i]for i in range(8))%q
   if B_:nzb=True
   A3=(A_*A_%q)*A_%q;B3=(B_*B_%q)*B_%q
   if sd==0:assert(A3-U_bls*w2_*B3)%q==0
   else:
    A6=A3*A3%q;B6=B3*B3%q;assert(A6*A3+(w_+1)*A6*B3+A3*B6-w_*B6*B3)%q==0
  assert nzb
  pM=1
  for k in range(4):
   ztk=ZT[k];mk=0
   for vec in (a1,b1,a2,b2):
    re=sum(vec[j]*ztk[j][0]for j in range(8));im=sum(vec[j]*ztk[j][1]for j in range(8))
    er=sum(abs(x)for x in vec)
    mk+=(abs(re)+er)**2+(abs(im)+er)**2
   pM*=mk
  assert pM*pM<=(602248*q)<<928
  if n_4k%500==0:qs_4k.append((n_4k,sd,q))
  n_4k+=1
 assert n_4k==4000
 min_k4k=min(x[0]for x in scert4k);min_s4k=min(x[1]for x in scert4k);med_s4k=sorted(x[1]for x in scert4k)[2000]
 assert min_k4k>=24 and min_s4k>=1000
 with localcontext()as _dc:
  _dc.prec=80
  p0=Decimal(10)/Decimal(5000);q0=Decimal(1)-p0;pmf=q0**4096;cdf=pmf
  for _k in range(1,24):pmf=pmf*Decimal(4096-_k+1)/Decimal(_k)*p0/q0;cdf+=pmf
  delta_mc=Decimal(1)-cdf
  assert delta_mc<Decimal("0.0000054")
  cp95=Decimal("0.0007487")
  assert(Decimal(1)-cp95)**4000<Decimal("0.05")
  assert cp95+delta_mc<Decimal(str(desc["params"]["eps_cert"]))
 bins38=[73+5*i for i in range(38)]
 def _froot(q_,sd_,rg_):
  ry_=next((pow(g_,(q_-1)//15,q_)for g_ in range(2,80)if pow(pow(g_,(q_-1)//15,q_),3,q_)!=1 and pow(pow(g_,(q_-1)//15,q_),5,q_)!=1),None)
  if ry_ is None:return None
  w_=pow(ry_,5,q_)
  if sd_==0:
   rhs_=(U_bls*w_*w_)%q_
   if pow(rhs_,(q_-1)//3,q_)!=1:return None
   rt_=pow(rhs_,(2*q_+1)//9 if q_%9==4 else(q_+2)//9,q_)
   return(ry_,rt_)if pow(rt_,3,q_)==rhs_ else None
  poly=[(-w_)%q_,1,(w_+1)%q_,1]
  def pm(a,b):
   r=[0]*(len(a)+len(b)-1)
   for i,x in enumerate(a):
    if x:
     for j,y in enumerate(b):r[i+j]=(r[i+j]+x*y)%q_
   while len(r)>3:
    c=r.pop()
    if c:
     for j in range(3):r[len(r)-3+j]=(r[len(r)-3+j]-c*poly[j])%q_
   return r
  def pp(bs,e):
   res=[1];b=bs[:]
   while e:
    if e&1:res=pm(res,b)
    e>>=1
    if e:b=pm(b,b)
   return res
  zq=pp([0,1],q_);zq+=[0]*max(0,2-len(zq));zq[1]=(zq[1]-1)%q_
  g_=pgcd_mod(poly,zq,q_)
  if len(g_)<2:return None
  while len(g_)>2:
   h_=pp([rg_.randrange(q_),1],(q_-1)//2);h_[0]=(h_[0]-1)%q_;d_=pgcd_mod(g_,h_,q_)
   if 1<len(d_)<len(g_):g_=d_
  z0=(-g_[0]*pow(g_[1],-1,q_))%q_
  if pow(z0,(q_-1)//3,q_)!=1:return None
  rt_=pow(z0,(2*q_+1)//9 if q_%9==4 else(q_+2)//9,q_)
  return(ry_,rt_)if pow(rt_,3,q_)==z0 else None
 for wid,(idx_w,sd_w,q_exp) in enumerate(qs_4k):
  if wid not in (0,5):continue
  rg_=random.Random(20261125+wid*1009);bts=bins38[(wid*500)%38];sd_=wid&1;assert sd_==sd_w
  while True:
   q_=(1<<(bts-1))|rg_.getrandbits(bts-1);q_+=(1-q_)%15
   while not(_is_prime_det(q_)and(q_%9 in (4,7))):q_+=15
   sr_=_froot(q_,sd_,rg_)
   if sr_ is not None:break
  assert q_==q_exp
  ry_,rt_=sr_;ryp_=[pow(ry_,i,q_)for i in range(8)];B_=[[0]*16 for _ in range(16)];B_[0][0]=q_
  for i in range(1,8):B_[i][0]=(-ryp_[i])%q_;B_[i][i]=1
  for i in range(8):B_[8+i][0]=(rt_*ryp_[i])%q_;B_[8+i][8+i]=1
  R_,_,_=lll_16d_int(B_,99,100);pr_=_OP.ob_pair(q_,ry_,rt_,R_);assert pr_ is not None
  a1_,b1_,a2_,b2_=pr_
  if wid in (0,5):
   res_=eval_scert_discrete_int(sd_,q_,ry_,rt_,a1_,b1_,a2_,b2_,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,vor_units,4096,seed=20261125+idx_w)
   assert res_["k5000"]==scert4k[idx_w][0]and int(round(res_["s_hat"]))==scert4k[idx_w][1]
 _prs=sorted(pr_rats)
 return dict(n_enc_live=12,parent_divisibility_live="12/12",exact_sat_det_eq_q="12/12",enc12_eta_bar_le_020_exact_int="12/12",fvor_cone_rays=ud_chk["rays"],fvor_m0=ud_chk["min_ray_coord"],fvor_decimal_prec=ud_chk["exact_decimal_prec"],sinnott_full_unit_group_exact_O_H=ud_fg,live_bareiss_99_100_all12=True,presupplied_vectors=0,enc_sext_f=sf6,enc_sext_g=sg6,min_headroom_f=round(min(qb_f[i]-sf6[i]for i in range(6)),2),min_headroom_g=round(min(qb_g[i]-sg6[i]for i in range(6)),2),joint_pr_over_charged_min=round(_prs[0],2),joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),enc12_scert_n_checked=len(sc_live),enc12_scert_int64_samps=1024,enc12_scert_min_k5000=min(k5_enc),enc12_scert_hat_min=round(min(sc_live),1),enc12_scert_hat_med=round(sorted(sc_live)[len(sc_live)//2],1),enc12_scert_nf_checked=nf_chk_tot,enc12_scert_max_nf_gap=round(max_nf_gap,2),prod4000_exact_int_verified="4000/4000",prod4000_k5000_pass="4000/4000",prod4000_k5000_min=min_k4k,prod4000_shat_min=min_s4k,prod4000_shat_med=med_s4k,prod4000_zero_cherrypick_lll_regen="2/2",prod4000_zero_cherrypick_scert_bitexact_regen="2/2",prod4000_markov_delta_mc="5.332e-6",prod4000_exact_cp95_plus_mc_ucb=0.0007541,all_sextiles_below=True,checks_ok=True)
def unplanted_two_sided_bw_check():
 import random
 p,r,m=467,233,381
 f_c,g_c=[1,1,-5,0,1],[-10,-1,1]
 assert sum(ci*pow(m,i,p)for i,ci in enumerate(f_c))%p==0
 assert sum(ci*pow(m,i,p)for i,ci in enumerate(g_c))%p==0
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
   if e&1:res=pmul(res,b,mod_p,q)
   e>>=1
   if e:b=pmul(b,b,mod_p,q)
  return res
 primes=[q for q in range(2,115)if all(q%d for d in range(2,int(q**0.5)+1))]
 def unram(c,q):
  for rt in range(q):
   if peval(c,rt,q)==0:
    if sum(i*c[i]*pow(rt,i-1,q)for i in range(1,len(c)))%q==0:return False
  return True
 fb_f,fb_g={},{}
 for q in primes:
  if unram(f_c,q):
   for rt in range(q):
    if peval(f_c,rt,q)==0:fb_f[(q,rt)]=len(fb_f)
  if unram(g_c,q):
   for rt in range(q):
    if peval(g_c,rt,q)==0:fb_g[(q,rt)]=len(fb_g)
 pf_act=sorted({q for(q,_) in fb_f})
 pg_act=sorted({q for(q,_) in fb_g})
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
   if a==0 or math.gcd(abs(a),b)!=1 or(a-b*m)%p==0:continue
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
  nxt=[ri for ri in cur if all(ccnt[k]>=2 for k in rels[ri][2])and all(ccnt[nf_cnt+k]>=2 for k in rels[ri][3])]
  if len(nxt)==len(cur):break
  cur=nxt
 used=sorted({k for ri in cur for k in rels[ri][2]}|{nf_cnt+k for ri in cur for k in rels[ri][3]})
 cmap={c:i for i,c in enumerate(used)}
 N_2core=len(used)
 def smap(a,b,mod_p,eps):
  el=ppow([a%r2,(-b)%r2],eps,mod_p,r2)
  assert el[0]%r==1 and all(x%r==0 for x in el[1:])
  return[((el[0]-1)//r)%r]+[(x//r)%r for x in el[1:]]
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
 m_rows=[row[:]for row in rows_M];s_rows=[row[:]for row in rows_S]
 alive_r=[True]*len(cur);alive_c=[True]*N_2core
 merge_stack=[]
 while True:
  piv_c=next((c for c in range(N_2core)if alive_c[c]and sum(1 for i in range(len(cur))if alive_r[i]and m_rows[i][c]!=0)==2),-1)
  if piv_c<0:break
  r_pair=[i for i in range(len(cur))if alive_r[i]and m_rows[i][piv_c]!=0]
  rp,ro=r_pair[0],r_pair[1]
  merge_stack.append((piv_c,m_rows[rp][:],s_rows[rp][:]))
  fac=(m_rows[ro][piv_c]*pow(m_rows[rp][piv_c],-1,r))%r
  m_rows[ro]=[(m_rows[ro][j]-fac*m_rows[rp][j])%r for j in range(N_2core)]
  s_rows[ro]=[(s_rows[ro][j]-fac*s_rows[rp][j])%r for j in range(4)]
  alive_r[rp]=False;alive_c[piv_c]=False
 core_cols=[c for c in range(N_2core)if alive_c[c]]
 N,D=len(core_cols),len(core_cols)+4
 active=[i for i in range(len(cur))if alive_r[i]]
 cc={c:sum(1 for i in active if m_rows[i][c]!=0)for c in core_cols}
 while len(active)>D:
  cands=[i for i in active if all(cc[c]>=3 for c in core_cols if m_rows[i][c]!=0)]or[i for i in active if all(cc[c]>=2 for c in core_cols if m_rows[i][c]!=0)]
  drop=max(cands,key=lambda i:(sum(1 for c in core_cols if m_rows[i][c]!=0),i))
  active.remove(drop)
  for c in core_cols:
   if m_rows[drop][c]!=0:cc[c]-=1
 def rref_null(mat):
  mr=[row[:]for row in mat]
  R,C=len(mr),len(mr[0])
  pivs,rk=[],0
  for c in range(C):
   piv=next((i for i in range(rk,R)if mr[i][c]%r!=0),-1)
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
  for fc in [c for c in range(C)if c not in pset]:
   v=[0]*C
   v[fc]=1
   for ri_i,pc in enumerate(pivs):
    v[pc]=(-mr[ri_i][fc])%r
   vecs.append(v)
  return rk,vecs
 M_sq=[[m_rows[i][c]for c in core_cols]for i in active]
 S_sq=[s_rows[i][:]for i in active]
 A_sq=[M_sq[i]+S_sq[i]for i in range(D)]
 rk_M,_=rref_null(M_sq)
 rk_A,ker_A=rref_null(A_sq)
 rk_3sm,ker_3sm=rref_null([M_sq[i]+S_sq[i][:3]for i in range(D)])
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
  return[sum(P[i][j]*v[j]for j in range(len(v)))%r for i in range(len(P))]
 rng=random.Random(20261119)
 m_bw=n_bw=4
 bnds=[(k*D)//4 for k in range(5)]
 X_good=[[0]*D for _ in range(m_bw)]
 for k in range(m_bw):
  for j in range(bnds[k],bnds[k+1]):
   X_good[k][j]=rng.randrange(1,r)
 Y=[[rng.randrange(1,r)for _ in range(n_bw)]for _ in range(D)]
 hr,hc=(D//m_bw)+2,(D//n_bw)+2
 AY,seq=[row[:]for row in Y],[]
 for _ in range(hr+hc+2):
  seq.append(mmul(X_good,AY))
  AY=mmul(A_sq,AY)
 H_rows=[[v for j in range(hc)for v in seq[i+j+1][ri_b]]for i in range(hr)for ri_b in range(m_bw)]
 rk_H,ker_H=rref_null(H_rows)
 v_rec=None
 for c_vec in ker_H:
  acc=[0]*D
  for j in range(hc-1,-1,-1):
   cj=c_vec[j*n_bw:(j+1)*n_bw]
   yc=[sum(Y[ri_i][bi]*cj[bi]for bi in range(n_bw))%r for ri_i in range(D)]
   acc=[(a_i+y_i)%r for a_i,y_i in zip(mvec(A_sq,acc),yc)]
  if any(acc)and all(x==0 for x in mvec(A_sq,acc)):
   v_rec=acc
   break
 assert v_rec is not None
 O_rows,cX=[],[row[:]for row in X_good]
 for _ in range(hr):
  cX=mmul(cX,A_sq)
  O_rows.extend(cX)
 rk_OA,_=rref_null(O_rows)
 C_cols,cY=[],[row[:]for row in Y]
 for _ in range(hc):
  for bi in range(n_bw):
   C_cols.append([cY[ri_i][bi]for ri_i in range(D)])
  cY=mmul(A_sq,cY)
 rk_Cd,_=rref_null(C_cols)
 X_bad=[[((k+1)*X_good[0][j])%r for j in range(D)]for k in range(m_bw)]
 O_bad,cXb=[],[row[:]for row in X_bad]
 for _ in range(hr):
  cXb=mmul(cXb,A_sq)
  O_bad.extend(cXb)
 rk_OA_bad,_=rref_null(O_bad)
 vlog_full={c:v_rec[j]for j,c in enumerate(core_cols)}
 sm_rec=v_rec[N:]
 for piv_c,pm,ps in reversed(merge_stack):
  acc=(sum(pm[j]*vlog_full[j]for j in range(N_2core)if j!=piv_c and pm[j]!=0)+sum(ps[s]*sm_rec[s]for s in range(4)))%r
  vlog_full[piv_c]=(-acc*pow(pm[piv_c],-1,r))%r
 while True:
  prog=False
  for idx in range(len(cur)):
   unk=[j for j in range(N_2core)if rows_M[idx][j]!=0 and j not in vlog_full]
   if len(unk)==1:
    cn=unk[0]
    acc=(sum(rows_M[idx][j]*vlog_full[j]for j in range(N_2core)if j!=cn and rows_M[idx][j]!=0)+sum(rows_S[idx][s]*sm_rec[s]for s in range(4)))%r
    vlog_full[cn]=(-acc*pow(rows_M[idx][cn],-1,r))%r
    prog=True
  if not prog:break
 assert len(vlog_full)==N_2core
 gen=next(g for g in range(2,p)if pow(g,(p-1)//r,p)!=1 and pow(g,2,p)!=1)
 g_r,dlog,ag=pow(gen,(p-1)//r,p),{},1
 for e in range(r):
  dlog[ag]=e
  ag=(ag*g_r)%p
 ratios=set()
 for idx,ri in enumerate(cur):
  a,b,vf,vg=rels[ri]
  lhs=(sum(e*vlog_full[cmap[k]]for k,e in vf.items())+sum(rows_S[idx][j]*sm_rec[j]for j in range(3)))%r
  rhs=(sum(e*vlog_full[cmap[nf_cnt+k]]for k,e in vg.items())-rows_S[idx][3]*sm_rec[3])%r
  assert lhs==rhs
  dl=dlog[pow((a-b*m)%p,(p-1)//r,p)]
  if dl!=0:ratios.add((lhs*pow(dl,-1,r))%r)
 pm_res=pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_rec)
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==D and rk_OA_bad<D-1 and len(ratios)==1
 return dict(p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(merge_stack),N=N,D=D,rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=rk_OA_bad,**pm_res,unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(ratios),checks_ok=True)
