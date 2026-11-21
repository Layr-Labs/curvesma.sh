"""Self-contained classic PM-Basis live verifier & Q(zeta_15) single-special-q Pure-k=1 Wald-i.i.d.-boot Boole H5 verifier."""
import base64,cmath,json,lzma,math,random
from decimal import Context,Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from fractions import Fraction as F
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
 z=[0]*16;nodes=0;hits=0;out=[];_sp=math.sumprod;_isq=math.isqrt
 lm_cols=[[lm[j][k]for j in range(k+1,16)]for k in range(16)]
 BrT=[[Br[i][t]for i in range(16)]for t in range(16)]
 def _rec(k,rem_sc,nz):
  nonlocal nodes,hits
  num_k=-_sp(lm_cols[k],z[k+1:])if k<15 else 0;dk1=dg[k+1];dk=dg[k]
  rhs=rem_sc*dk;rad=_isq(rhs)
  lo=(num_k-rad+dk1-1)//dk1;hi=(num_k+rad)//dk1
  if not nz and lo<0:lo=0
  for zk in range(lo,hi+1):
   z[k]=zk;nz2=nz or(zk>0);diff=zk*dk1-num_k;rnum=rhs-diff*diff
   if rnum<0:continue
   if k==0:
    if not nz2:continue
    hits+=1
    if collect:
     v=[_sp(z,BrT[t])for t in range(16)]
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
 S=1<<32;L63=(1<<63)-1;N63=-(1<<63);_sp=math.sumprod
 m_act=16
 while m_act>0 and dg[m_act]>R2*dg[m_act-1]:m_act-=1
 if m_act==0:return 0,0,[],8
 for k in range(m_act):
  assert dg[k]<=642*dg[k+1]and dg[k+1]<=642*R2*dg[k]and _sp(Br[k],Br[k])<=3100*R2
 MU=[[0]*m_act for _ in range(m_act)]
 for j in range(m_act):
  for k in range(j):
   dk1=dg[k+1];num=lm[j][k]*S
   mv=(num+dk1//2)//dk1 if num>=0 else-((-num+dk1//2)//dk1)
   assert N63<=mv<=L63;MU[j][k]=mv
 L2_LO=[(dg[k+1]*S)//dg[k]for k in range(m_act)]
 INV_HI=[math.isqrt((dg[k]*S*S+dg[k+1]-1)//dg[k+1])+1 for k in range(m_act)]
 assert all(N63<=v<=L63 for v in L2_LO)and all(N63<=v<=L63 for v in INV_HI)
 C_stk=[[0]*m_act for _ in range(m_act)];V_stk=[[0]*16 for _ in range(m_act+1)];v_cur=[0]*16
 z=[0]*m_act;hi_s=[0]*m_act;rem_s=[0]*m_act;nz_s=[0]*m_act;asum_s=[0]*m_act;first_s=[True]*m_act
 nbk=[0]*m_act;ops=8;hits=0;out=[];b0=Br[0]
 def _enter(k,rem_q,nz,asum):
  nonlocal ops,v_cur
  cq=C_stk[k][k];ec=(asum+1)>>1
  if rem_q==0:sr=1
  else:
   x2=1<<((65-(64-rem_q.bit_length()))>>1)
   for _ in range(5):x2=(x2+(rem_q//x2))>>1
   sr=x2+1
  rq=((sr*INV_HI[k])>>16)+1+ec;c1=cq-rq+S-1;c2=cq+rq
  assert N63<=rq<=L63 and N63<=c1<=L63 and N63<=c2<=L63
  lo=c1>>32;hi=c2>>32
  if not nz and lo<0:lo=0
  z[k]=lo;hi_s[k]=hi;rem_s[k]=rem_q;nz_s[k]=nz;asum_s[k]=asum;first_s[k]=True;ops+=55
  if k==0 and lo<=hi:
   lm1=lo-1;base0=V_stk[0]
   v_cur=[base0[t]+lm1*b0[t]for t in range(16)];ops+=40
   assert N63<=min(v_cur)and max(v_cur)<=L63
 k=m_act-1;_enter(k,R2<<32,0,0)
 while True:
  if z[k]>hi_s[k]:
   ops+=12;k+=1
   if k==m_act:break
   z[k]+=1;continue
  zk=z[k];nz2=1 if(nz_s[k]or zk>0)else 0
  if k==0:
   v_cur=[v_cur[t]+b0[t]for t in range(16)]
   assert N63<=min(v_cur)and max(v_cur)<=L63
   if nz2 and _sp(v_cur,v_cur)<=R2:
    hits+=1
    if collect:
     vc=list(v_cur)
     for x in vc:
      if x!=0:
       if x<0:vc=[-t for t in vc]
       break
     out.append(tuple(vc))
   else:ops+=72
   z[0]+=1;continue
  cq=C_stk[k][k];ec=(asum_s[k]+1)>>1
  z32=zk<<32;df=z32-cq
  assert N63<=z32<=L63 and N63<=df<=L63
  da=-df if df<0 else df;dl=da-ec if da>ec else 0
  sq=(((dl*dl)>>32)*L2_LO[k])>>32;rem_q=rem_s[k]
  if sq>rem_q:ops+=20;z[k]+=1;continue
  nbk[k]+=1;vk1=V_stk[k];rk=Br[k]
  vn=[vk1[t]+zk*rk[t]for t in range(16)]
  assert N63<=min(vn)and max(vn)<=L63;V_stk[k-1]=vn
  ckm1=C_stk[k-1];muk=MU[k]
  if first_s[k]:
   ck=C_stk[k]
   for m in range(k):
    p=muk[m]*zk;v=ck[m]-p
    assert N63<=p<=L63 and N63<=v<=L63;ckm1[m]=v
   first_s[k]=False
  else:
   for m in range(k):
    v=ckm1[m]-muk[m]
    assert N63<=v<=L63;ckm1[m]=v
  ops+=31+40+(6*k+4);k-=1;_enter(k,rem_q-sq,nz2,asum_s[k+1]+abs(zk))
 assert ops<=sum((172+6*j+(260 if j==1 else 0))*nbk[j]for j in range(1,m_act))+140+(260 if m_act==1 else 0)
 return hits,sum(nbk),out,ops
def verify_b16_dsieve_live(data,desc,dlive,m_node,U_bls,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 _Rm=max(int(k)for k in desc["cells_rank2"]);_M=_Rm*_Rm
 def _ksq_b(poly,bpc):
  _x=int.from_bytes(b"".join(v.to_bytes(bpc,"little")for v in poly),"little");_p=_x*_x
  _nb=max((_M+1)*bpc,(_p.bit_length()+7)//8);_bs=_p.to_bytes(_nb,"little")
  return[int.from_bytes(_bs[i*bpc:(i+1)*bpc],"little")for i in range(_M+1)]
 _t1=[0]*(_M+1);_t1[0]=1
 for _k in range(1,_Rm+1):_t1[_k*_k]=2
 _t16=_ksq_b(_ksq_b(_ksq_b(_ksq_b(_t1,2),4),8),14)
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
 def _cnt_orb_R(_va,_vb,_QinvT,_qq,_Rop2,_R2):
  _eks=[_iv_mul((Decimal(2),Decimal(2)),_iv_add(_iv_norm2(_va,k),_iv_norm2(_vb,k)))for k in range(4)]
  _bd_hi=_CH.multiply(_d15001,Decimal(_Rop2*_R2));_bd_i=int(_CH.multiply(_bd_hi,Decimal(_SC4)).to_integral_value(rounding=ROUND_CEILING))
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
  _mo=0;_qq2=_qq*_qq;_R2q2=_R2*_qq2;_sp=math.sumprod
  _e2_rng=[(_n2,_epow_lo(2,_n2))for _n2 in range(_nlo[2],_nhi[2]+1)]
  _e1_rng=[(_n1,_epow_lo(1,_n1))for _n1 in range(_nlo[1],_nhi[1]+1)]
  for _n0 in range(_nlo[0],_nhi[0]+1):
   _e0=_epow_lo(0,_n0)
   _ek00,_ek01,_ek02,_ek03=_eks_i[0]*_e0[0],_eks_i[1]*_e0[1],_eks_i[2]*_e0[2],_eks_i[3]*_e0[3]
   for _n1,_e1 in _e1_rng:
    _a0,_a1,_a2,_a3=_ek00*_e1[0],_ek01*_e1[1],_ek02*_e1[2],_ek03*_e1[3]
    for _n2,_e2 in _e2_rng:
     if _a0*_e2[0]+_a1*_e2[1]+_a2*_e2[2]+_a3*_e2[3]>_bd_i:continue
     _u=mul_phi15_z(mul_phi15_z(_upow(0,_n0),_upow(1,_n1)),_upow(2,_n2))
     a0,a1,a2,a3,a4,a5,a6,a7=mul_phi15_z(_u,_va);b0,b1,b2,b3,b4,b5,b6,b7=mul_phi15_z(_u,_vb)
     for _ in range(15):
      _w=(a0,a1,a2,a3,a4,a5,a6,a7,b0,b1,b2,b3,b4,b5,b6,b7)
      _s=0
      for _qt in _QinvT:
       _qj=_sp(_w,_qt);_s+=_qj*_qj
       if _s>_R2q2:break
      else:_mo+=1
      # Exact multiplication by y mod Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1 (6 adds vs 64 muls!)
      a0,a1,a2,a3,a4,a5,a6,a7=-a7,a0+a7,a1,a2-a7,a3+a7,a4-a7,a5,a6+a7
      b0,b1,b2,b3,b4,b5,b6,b7=-b7,b0+b7,b1,b2-b7,b3+b7,b4-b7,b5,b6+b7
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
 # Unconditional m-truncated (9/10)-Lovasz Q32 i64 proof test on diag(1,...,1,p70) (m=15), diag(1,p70,...,p70) (m=1), & diag(p70,...,p70) (m=0):
 _p70d=(1<<70)+15;_B_d15=[[(_p70d if i==15 else 1)if i==j else 0 for j in range(16)]for i in range(16)]
 _Br_d15,_dg_d15,_lm_d15=lll_16d_int(_B_d15,9,10);assert _dg_d15[16]==_p70d*_p70d
 _hd15_2,_,_hd15_h,_=_fp_q32_int(_Br_d15,_dg_d15,_lm_d15,4,collect=True)
 _fe15_2,_,_fe15_h=_fp_exact_int(_Br_d15,_dg_d15,_lm_d15,4,collect=True)
 assert _hd15_2==_fe15_2>300 and set(_hd15_h)==set(_fe15_h)and all(v[15]==0 for v in _hd15_h)
 _B_d15_64=[[(_p70d if i==15 else 64)if i==j else 0 for j in range(16)]for i in range(16)]
 _Br_d64,_dg_d64,_lm_d64=lll_16d_int(_B_d15_64,9,10)
 _hd15_128,_,_,_=_fp_q32_int(_Br_d64,_dg_d64,_lm_d64,128*128,collect=False)
 assert _hd15_128==_hd15_2
 _B_d1=[[(1 if i==0 else _p70d)if i==j else 0 for j in range(16)]for i in range(16)]
 _Br_d1,_dg_d1,_lm_d1=lll_16d_int(_B_d1,9,10)
 _hd1_2,_,_hd1_h,_=_fp_q32_int(_Br_d1,_dg_d1,_lm_d1,4,collect=True);assert _hd1_2==2 and set(_hd1_h)=={tuple([1]+[0]*15),tuple([2]+[0]*15)}
 _B_d0=[[_p70d if i==j else 0 for j in range(16)]for i in range(16)]
 _Br_d0,_dg_d0,_lm_d0=lll_16d_int(_B_d0,9,10)
 _hd0_128,_,_,_=_fp_q32_int(_Br_d0,_dg_d0,_lm_d0,128*128,collect=True);assert _hd0_128==0
 # Exact Fraction verification of Newton-Maclaurin without-replacement inequality e_t(f)/C(N,t) <= (e_1(f)/N)^t:
 for _fv in ([1]*8+[0]*8,[F(1,3),F(4,5),F(1,2),0,1,F(2,7),F(9,10),0,F(3,4),F(1,5)]):
  _Nf=len(_fv);_el=[F(1)]+[F(0)]*_Nf
  for _x in _fv:
   for _t in range(_Nf,0,-1):_el[_t]+=_el[_t-1]*_x
  _m1=_el[1]/_Nf
  for _t in range(1,_Nf+1):assert _el[_t]/math.comb(_Nf,_t)<=_m1**_t
 assert F(math.comb(96,24),math.comb(192,24))<F(1,1<<24)
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
 def _struct_L_R(_R,_p,_s,_t):
  _sp=[pow(_s,j,_p)for j in range(8)]
  _c=[(sum(_R[j][k]*_sp[k]for k in range(8))-_t*sum(_R[j][8+k]*_sp[k]for k in range(8)))%_p for j in range(16)]
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
 _ds_disc={};_gate_ratios=[];_R_cache={}
 for _yi,_lbl in ((0,"g"),(12,"f")):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16)
  _Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_qq;_rp=[pow(_ry,i,_qq)for i in range(8)]
  for _i in range(1,8):_Bq[_i][0]=(-_rp[_i])%_qq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=(_rt*_rp[_i])%_qq;_Bq[8+_i][8+_i]=1
  _Rr,_dsq_q,_=lll_16d_int(_Bq,99,100);assert _dsq_q[16]==_qq*_qq;_R_cache[_yi]=_Rr
  if _yi==0:_f16_ba=[_Rr[j][:8]for j in range(16)]
  _p83,_s83,_t83=_find_prime_Q(84,1084+_yi*10)
  _Br83,_dg83,_lm83=lll_16d_int(_struct_L_R(_Rr,_p83,_s83,_t83),9,10);assert _dg83[16]==_p83*_p83
  _hq83,_ndq83,_,_ops83=_fp_q32_int(_Br83,_dg83,_lm83,64*64,collect=False)
  _ls83=[math.sqrt(_dg83[i+1]/_dg83[i])for i in range(16)];_bi83=_b_int_64(_ls83)
  _ir83=140.0;_dp83=[1.0]
  for _k83 in range(15,0,-1):
   _u83=1.0/_ls83[_k83];_nd83_l=[0.0]*(len(_dp83)+1)
   for _s83,_v83 in enumerate(_dp83):_nd83_l[_s83]+=_v83;_nd83_l[_s83+1]+=_v83*_u83
   _dp83=_nd83_l;_nk83=0.5*(sum(_vb64[_s83]*_dp83[_s83]for _s83 in range(len(_dp83)))+1.0)
   _ir83+=(172.0+6.0*_k83+(260.0 if _k83==1 else 0.0))*_nk83
  assert _hq83>400 and _ndq83<=_bi83 and _ops83<=_ir83
  if _lbl=="g":
   _h83,_nd83,_=_fp_exact_int(_Br83,_dg83,_lm83,64*64,collect=False);assert _h83==_hq83 and _nd83==_ndq83
  _ds_disc[f"{_lbl}_p84_R64"]=dict(hits=_hq83,nodes_q32=_ndq83,ops_sim=_ops83,B_int_ub=round(_bi83,1),actual_over_ub=round(_ndq83/_bi83,4))
  for _bp in (16,48,60,70):
   _pb,_sb,_tb=_find_prime_Q(_bp,2000+_yi*100+_bp)
   _,_dgb,_=lll_16d_int(_struct_L_R(_Rr,_pb,_sb,_tb),9,10);assert _dgb[16]==_pb*_pb
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
  ep1=e30+1;m_hi=sq+ep1;wc=(sq-ep1>=64*ep1)
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
   if _bg_wc:assert _bg-math.log2(_ng_ex)<=0.5343
   _mgg=max(_mgg,_bg-math.log2(_ng_ex))
  if _nf_ex>0:
   assert _bf>=math.log2(_nf_ex)
   if _bf_wc:assert _bf-math.log2(_nf_ex)<=0.5343
   _mgf=max(_mgf,_bf-math.log2(_nf_ex))
 assert _mgf<0.25 and _mgg<0.25 and _nc_g==0 and _nc_f==0
 # Live literal -57.3b near-cancellation witness: verifies nc==1 (sq - e30 < 64*e30) fires & routes to exact NORMS (0 false negatives)!
 _anc=[-15760645,-18238986,-10524099,103233,763938,0,0,0];_bnc=[-17,0,0,0,0,0,0,0]
 _a3n,_b3n=pow_phi15_z(_anc,3),pow_phi15_z(_bnc,3)
 _ng_nc=abs(norm_q15_z(sub_phi15_z(_a3n,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3n)))))
 _bg_nc_q=0;_gwc_nc=True
 for _k in range(4):
  _sa=_ev_lin(_anc,_k);_sb=_ev_lin(_bnc,_k)
  _a3b=_bmul(_bmul(_sa,_sa),_sa);_b3b=_bmul(_bmul(_sb,_sb),_sb)
  _tgb=_bmul(_uw2_pre[_k],_b3b);_gv=_badd(_a3b,(-_tgb[0],-_tgb[1],_tgb[2],_tgb[3]))
  _gh,_gwc=_bl2_hi(_gv);_bg_nc_q+=2*_gh
  if not _gwc:_gwc_nc=False
 _bg_nc=_bg_nc_q/65536.0;_nc_gap=round(_bg_nc-math.log2(_ng_nc),2)
 assert(not _gwc_nc)and _bg_nc>=math.log2(_ng_nc)and _nc_gap>10.0
 _yr_min=1e99;_eff_u=[];_ss_yr=dlive["_ss_eval"]
 _bl4=[]
 for _qp in range(2049,8192,2):
  if _qp%15!=1 or not _is_prime_det(_qp):continue
  _e15=(_qp-1)//15;_g0=next(pow(g,_e15,_qp)for g in range(2,80)if pow(pow(g,_e15,_qp),3,_qp)!=1 and pow(pow(g,_e15,_qp),5,_qp)!=1)
  _wp=pow(_g0,5,_qp)
  _rts=[_rx for _rx in range(1,_qp)if(pow((_rx*_rx*_rx)%_qp,3,_qp)+(_wp+1)*pow((_rx*_rx*_rx)%_qp,2,_qp)+(_rx*_rx*_rx)%_qp-_wp)%_qp==0]
  if len(_rts)>=2:_bl4.append((_qp,_g0,_rts[0],_rts[1]))
  if len(_bl4)==4:break
 assert len(_bl4)==4
 for _bi,_yi in enumerate((0,5,12,17)):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_sn="g" if _sd==0 else "f"
  _qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16);_nb=_qq.bit_length()-1
  _qp,_ryp_b,_rtp_b,_rtp_o=_bl4[_bi];_Mq=_qq*_qp;_iq=pow(_qq,-1,_qp)
  _prw=next(r for r in dlive["rows"]if r["side"]==_sn and r["lo"]<=_nb<r["hi"])
  _Rp=_prw["R"];_Rp2=_Rp*_Rp
  _rp=[pow(_ry,i,_qq)for i in range(8)];_rpb=[pow(_ryp_b,i,_qp)for i in range(8)]
  _Yq=[u+_qq*(((v-u)*_iq)%_qp)for u,v in zip(_rp,_rpb)]
  _Tq=[u+_qq*(((v-u)*_iq)%_qp)for u,v in zip([(_rt*u)%_qq for u in _rp],[(_rtp_b*v)%_qp for v in _rpb])]
  _Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_Mq
  for _i in range(1,8):_Bq[_i][0]=(-_Yq[_i])%_Mq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=_Tq[_i]%_Mq;_Bq[8+_i][8+_i]=1
  _Rr,_dsq_y,_=lll_16d_int(_Bq,99,100);assert _dsq_y[16]==_Mq*_Mq
  # Verify live that any vector in L_{Q, Q'} also vanishing at sibling root _rtp_o != _rtp_b above same O_H-prime is non-coprime in O_H:
  _c_oth=[(sum(_Rr[j][k]*_rpb[k]for k in range(8))-_rtp_o*sum(_Rr[j][8+k]*_rpb[k]for k in range(8)))%_qp for j in range(16)]
  _k0o=next(i for i,v in enumerate(_c_oth)if v!=0);_ivo=pow(_c_oth[_k0o],-1,_qp)
  _j1o=next(i for i in range(16)if i!=_k0o);_vo=[0]*16;_vo[_j1o]=1;_vo[_k0o]=(-_c_oth[_j1o]*_ivo)%_qp
  _a_o=[sum(_vo[j]*_Rr[j][t]for j in range(16))for t in range(8)];_b_o=[sum(_vo[j]*_Rr[j][8+t]for j in range(16))for t in range(8)]
  assert(sum(_a_o[k]*_rpb[k]for k in range(8))-_rtp_o*sum(_b_o[k]*_rpb[k]for k in range(8)))%_qp==0 and not _is_cp(_a_o,_b_o)
  _nfs,_ngs,_wts=[],[],[]
  _aug=[[F(_Rr[i][j])for j in range(16)]+[F(1 if i==j else 0)for j in range(16)]for i in range(16)]
  for _c in range(16):
   _pv=next(r for r in range(_c,16)if _aug[r][_c]!=0);_aug[_c],_aug[_pv]=_aug[_pv],_aug[_c]
   _iv=F(1)/_aug[_c][_c];_aug[_c]=[x*_iv for x in _aug[_c]]
   for _r in range(16):
    if _r!=_c and _aug[_r][_c]!=0:
     _fc=_aug[_r][_c];_aug[_r]=[_aug[_r][j]-_fc*_aug[_c][j]for j in range(32)]
  _QinvT=[[int(_aug[i][16+j]*_Mq)for i in range(16)]for j in range(16)]
  _Rop2=max(sum(abs(math.sumprod(_Rr[i],_Rr[j]))for j in range(16))for i in range(16));_rp_y=[pow(_ry,i,_qq)for i in range(8)]
  for _ in range(12):
   while True:
    _g=[_rng.gauss(0,1)for _ in range(16)];_rn=math.sqrt(sum(x*x for x in _g))
    _rd=_Rp*(_rng.random()**(1.0/16.0))/_rn;_vv=[int(round(x*_rd))for x in _g]
    if 0<sum(x*x for x in _vv)<=_Rp2:break
   _aa=[sum(_vv[j]*_Rr[j][t]for j in range(16))for t in range(8)]
   _bb=[sum(_vv[j]*_Rr[j][8+t]for j in range(16))for t in range(8)]
   if(sum(_aa[j]*_rp_y[j]for j in range(8))%_qq==0 and sum(_bb[j]*_rp_y[j]for j in range(8))%_qq==0)or not _is_cp(_aa,_bb):
    _wts.append(0.0);_nfs.append(9999.0);_ngs.append(9999.0);continue
   _wts.append(1.0/_cnt_orb_R(_aa,_bb,_QinvT,_Mq,_Rop2,_Rp2))
   _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
   _ng=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
   _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
   _nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
   assert _nf%_qp==0 and(_ng if _sn=="g" else(_nf//_qp))%_qq==0
   _nfs.append(math.log2(_nf)-((math.log2(_qq)if _sn=="f" else 0.0)+math.log2(_qp)));_ngs.append(math.log2(_ng)-(math.log2(_qq)if _sn=="g" else 0.0))
  _pf=_ss_yr(_nfs,float(_prw["n"]),_prw["iA"],m_node,float(desc["params"]["alpha_f"]),_prw["kA"],_prw["kB"])
  _pg=_ss_yr(_ngs,float(_prw["n"]),_prw["iA"],m_node,float(desc["params"]["alpha_g"]),_prw["kA"],_prw["kB"])
  _es=sum(a*b for a,b in zip(_pf,_pg));_esm=sum(w*a*b for w,a,b in zip(_wts,_pf,_pg))
  _yr=(_esm/len(_wts))/_prw["pr"]
  if _es>0:_eff_u.append(round(_esm/_es,4))
  if _yr<_yr_min:_yr_min=_yr
 assert _yr_min>=10.0*float(desc["params"]["udup_node"])
 # Live residue-field kernel lattice check for BOTH F_16 = (Z/2Z)^4 ideals P1=(2,y^4+y+1) & P2=(2,y^4+y^3+1) (N(P)=16, f=4):
 from itertools import combinations as _combs,product as _prods
 def _f16_red(poly,mp):
  v=0
  for i,c in enumerate(poly):
   if c&1:
    r=1
    for _ in range(i):
     r<<=1
     if r&16:r^=mp
    v^=r
  return v
 _f16_res={}
 _Re68=68.05;_vb68=[1.0]+[(math.pi**(d/2.0)/math.gamma(d/2.0+1.0))*(_Re68**d)for d in range(1,17)]
 _al=100.0/74.0;_cw=(16.0**(1.0/16.0))*(_al**(60.0/16.0));_iw=[1.0/(_cw*(_al**(-j/2.0)))for j in range(16)]
 _wc16=0.0;_dpw=[1.0]
 for _k in range(15,0,-1):
  _u=_iw[_k];_nd=[0.0]*(len(_dpw)+1)
  for _s,_v in enumerate(_dpw):_nd[_s]+=_v;_nd[_s+1]+=_v*_u
  _dpw=_nd;_wc16+=0.5*(sum(_vb68[_s]*_dpw[_s]for _s in range(len(_dpw)))+1.0)
 for _nm,_mp in (("P1_f16",19),("P2_f16",25)):
  _cf=[_f16_red(_f16_ba[j],_mp)for j in range(16)]
  _M=[[( _cf[j]>>r)&1 for j in range(16)]for r in range(4)]
  _piv=[];_rr=0
  for _c in range(16):
   _pr=next((i for i in range(_rr,4)if _M[i][_c]),None)
   if _pr is None:continue
   _M[_rr],_M[_pr]=_M[_pr],_M[_rr]
   for i in range(4):
    if i!=_rr and _M[i][_c]:
     for j in range(16):_M[i][j]^=_M[_rr][j]
   _piv.append(_c);_rr+=1
  assert len(_piv)==4;_ps=set(_piv);_BL=[]
  for _pc in _piv:
   _v=[0]*16;_v[_pc]=2;_BL.append(_v)
  for _c in range(16):
   if _c in _ps:continue
   _v=[0]*16;_v[_c]=1
   for _ri,_pc in enumerate(_piv):
    if _M[_ri][_c]:_v[_pc]=1
   _BL.append(_v)
  _Rl,_dsq,_lsq=lll_16d_int(_BL,9,10);assert _dsq[16]==256
  _qcnt,_qnod,_qhits,_=_fp_q32_int(_Rl,_dsq,_lsq,4,collect=True)
  _bf=set()
  for i in range(16):_v=[0]*16;_v[i]=2;_bf.add(tuple(_v))
  for _wk in (1,2,3,4):
   for _ix in _combs(range(16),_wk):
    _xv=0
    for _id in _ix:_xv^=_cf[_id]
    if _xv:continue
    for _sg in _prods((-1,1),repeat=_wk):
     if _sg[0]<0:continue
     _v=[0]*16
     for _id,_s in zip(_ix,_sg):_v[_id]=_s
     _bf.add(tuple(_v))
  _qs={tuple(v if next(x for x in v if x)>0 else[-x for x in v])for v in _qhits}
  assert _qs==_bf
  _ls=[math.sqrt(_dsq[k+1]/_dsq[k])for k in range(16)];_bi=0.0;_dp=[1.0]
  for _k in range(15,0,-1):
   _u=1.0/_ls[_k];_nd=[0.0]*(len(_dp)+1)
   for _s,_v in enumerate(_dp):_nd[_s]+=_v;_nd[_s+1]+=_v*_u
   _dp=_nd;_bi+=0.5*(sum(_vb68[_s]*_dp[_s]for _s in range(len(_dp)))+1.0)
  _rat=round(_bi/_wc16,4);assert _rat<0.60
  _f16_res[_nm]=dict(det=16,r2_hits=len(_bf),mismatches=0,b_int_over_wc_R68=_rat)
 # Live 2-adic g-side prime-power valuation recovery check: P_i^k | a (k=1..6, weights 12,12,12,12,12,4 -> exact 64b per P_i, 128b total -> 0 deficit!)
 _w2_wt=(12,12,12,12,12,4);_b1=[1]+[0]*7;_b1_3=pow_phi15_z(_b1,3);_p1g=[1,1,0,0,1,0,0,0]
 for _k in range(1,8):
  _ng2=abs(norm_q15_z(sub_phi15_z(pow_phi15_z([1<<_k]+[0]*7,3),scale_phi15_z(U_bls,mul_phi15_z(_w2,_b1_3)))))
  assert((_ng2&-_ng2).bit_length()-1)==2*sum(_w2_wt[:min(_k,6)])
  _ngp=abs(norm_q15_z(sub_phi15_z(pow_phi15_z(pow_phi15_z(_p1g,_k),3),scale_phi15_z(U_bls,mul_phi15_z(_w2,_b1_3)))))
  assert((_ngp&-_ngp).bit_length()-1)==sum(_w2_wt[:min(_k,6)])
 return round(min(_eff_u),4),round(_yr_min,2),len(_q32_h),dict(discrete_R64=_ds_disc,gate_8primes_max=max(_gate_ratios),gate_8primes=_gate_ratios,f16_small_ideals=_f16_res,p2_adic_128b_zero_deficit=True,bound_int_max_gap_bits=dict(f=round(_mgf,4),g=round(_mgg,4)),near_cancel_witness=dict(nc_triggered=True,gap_bits=_nc_gap,routed_to_exact_norms=True))
def verify_q15_93_b16_live(camp93,primes_upto,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 U_TOY=-5;W_EL=y_pow_z(5);WP1=add_phi15_z(W_EL,[1]+[0]*7);W2=y_pow_z(10);out={}
 for grp in ("g","f","enc"):
  rows=camp93[grp];assert len(rows)==40
  cnts=[r[4]if grp=="enc" else r[3]for r in rows];mus=[r[5]if grp=="enc" else r[4]for r in rows]
  mc=sum(cnts)/40.0;vc=sum((x-mc)**2 for x in cnts)/39.0;fano=vc/mc;k_hat=max(0.0,(vc-mc)/(mc*mc))
  emp0=cnts.count(0)/40.0;mod0_c=sum(1.0/(1.0+0.5*m)for m in mus)/40.0
  assert fano<1.5 and k_hat<0.25 and emp0<mod0_c
  out[grp]=dict(n=40,p0_mode="cantelli05",c_mu=1.0,mean_M=round(mc,2),mean_mu_raw=round(sum(mus)/40.0,2),fano=round(fano,3),k_hat=round(k_hat,4),zeros=cnts.count(0),emp0=round(emp0,4),cantelli_empty=round(mod0_c,4))
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
 p0_f,ry0_f,rt0_f,_,_=camp93["f"][0]
 B0=[[0]*16 for _ in range(16)];B0[0][0]=p0_f;rp0=[pow(ry0_f,i,p0_f)for i in range(8)]
 for i in range(1,8):B0[i][0]=(-rp0[i])%p0_f;B0[i][i]=1
 for i in range(8):B0[8+i][0]=(rt0_f*rp0[i])%p0_f;B0[8+i][8+i]=1
 Rr0,dsq0,_=lll_16d_int(B0,99,100);assert dsq0[16]==p0_f*p0_f
 seen0=set()
 for nz in pts[:512]:
  aa=[0]*8;bb=[0]*8
  for idx,sgn in nz:
   row=Rr0[idx]
   for t in range(8):aa[t]+=sgn*row[t];bb[t]+=sgn*row[8+t]
  a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3)
  ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_TOY,mul_phi15_z(W2,b3)))))
  if ng==0 or not _ok_sw(ng):continue
  a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
  nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_EL,b9))))
  if nf==0:continue
  assert nf%p0_f==0;nf//=p0_f
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
  if cp:seen0.add((nf,ng))
 assert len(seen0)>=5
 out["direct_16vec_Z_scan_only"]=True;out["live_regenerated_directR_cnt"]=len(seen0);out["checks_ok"]=True
 return out
def verify_q15_enc_bias_live(c7e,desc,U_bls,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,lll_16d_int):
 pA=c7e["partA_ubls_enc_vs_rand"];pC=c7e["partC_smoothness_weighted"]
 assert pA["n_enc"]==120 and pA["discards"]==0 and pA["ks_cf_p"]>=0.10 and pA["ks_cg_p"]>=0.05
 assert pC["n_enc"]==120 and pC["discards"]==0 and abs(pC["diff_hf"])<=2.0*pC["se_diff_hf"]
 W_EL=y_pow_z(5);WP1=add_phi15_z(W_EL,[1]+[0]*7);W2=y_pow_z(10)
 log2_u8=float(desc["params"]["log2_u8"]);S8=math.log2(desc["cells_rank2"]["8"])
 qb_f=[x[1]for x in desc["qbins"]];qb_g=[x[2]for x in desc["qbins"]]
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<12:
  al=[rng.randint(-4,4)for _ in range(8)];bt=[rng.randint(-4,4)for _ in range(8)]
  if any(al)and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
 live_A_ok=0
 for wlist in (pA["live_witnesses_enc"],pA["live_witnesses_rand"]):
  for sd,p,ry,rt in wlist:
   rp=[pow(ry,i,p)for i in range(8)];B=[[0]*16 for _ in range(16)];B[0][0]=p
   for i in range(1,8):B[i][0]=(-rp[i])%p;B[i][i]=1
   for i in range(8):B[8+i][0]=(rt*rp[i])%p;B[8+i][8+i]=1
   R,dsq,_=lll_16d_int(B,99,100);assert dsq[16]==p*p
   c_v=(math.log2(p)+S8)/16.0
   for al,bt in ab_pts:
    vv=al+bt;a=[sum(vv[j]*R[j][t]for j in range(16))for t in range(8)];b=[sum(vv[j]*R[j][8+t]for j in range(16))for t in range(8)]
    a3,b3=pow_phi15_z(a,3),pow_phi15_z(b,3)
    ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_bls,mul_phi15_z(W2,b3)))))
    a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
    nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_EL,b9))))
    assert(ng if sd==0 else nf)%p==0
    assert math.log2(nf)-72.0*c_v<qb_f[-1]and math.log2(ng)-24.0*c_v-log2_u8<qb_g[-1]
   live_A_ok+=1
 assert live_A_ok==8
 return dict(partA_ubls_enc_vs_rand={k:v for k,v in pA.items()if not k.startswith("live_")},partA_live_reverified="8/8",partC_smoothness_weighted=pC,checks_ok=True)

def verify_prod_cond_q15(enc12,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,lll_16d_int,vor_units=None,pgcd_mod=None):
 import unitdomain as _UD
 assert len(enc12)==12
 ud_chk=_UD.check(vor_units,mul_fn,y_pow_fn);ud_fg=_UD.full_group(vor_units,mul_fn,y_pow_fn,sub_fn,add_fn,norm_fn)
 assert ud_chk["min_ray_coord"]>0.028 and ud_fg==9
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 log2_u8=float(desc["params"]["log2_u8"]);S8=math.log2(desc["cells_rank2"]["8"])
 qb_f=[x[1]for x in desc["qbins"]];qb_g=[x[2]for x in desc["qbins"]]
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4)for _ in range(8)];bt=[rng.randint(-4,4)for _ in range(8)]
  if any(al)and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
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
  ry_p=[pow(ry,i,q)for i in range(8)]
  assert(sum(pa[i]*ry_p[i]for i in range(8))-rt*sum(pb[i]*ry_p[i]for i in range(8)))%q==0
  Tq=[(rt*u)%q for u in ry_p]
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  for i in range(1,8):B[i][0]=(-ry_p[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=Tq[i]%q;B[8+i][8+i]=1
  R,dsq,lsq=lll_16d_int(B,99,100);assert dsq[16]==q*q
  for kk in range(1,16):assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  sname="g" if sid==0 else "f";nb=q.bit_length()-1;lq=math.log2(q)
  prow=next(r for r in dlive["rows"]if r["side"]==sname and r["lo"]<=nb<r["hi"])
  R_n=prow["R"];cells_b=int(desc["cells_rank2"][str(R_n)]);S_bin=math.log2(cells_b)
  c_act=(lq+S_bin)/16.0;c_val=(lq+S8)/16.0;nfas,ngas=[],[]
  for al,bt in ab_pts:
   vv=al+bt
   a=[sum(vv[j]*R[j][t]for j in range(16))for t in range(8)]
   b=[sum(vv[j]*R[j][8+t]for j in range(16))for t in range(8)]
   a3,b3=pow_fn(a,3),pow_fn(b,3)
   ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
   a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
   nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,pb3 if False else b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
   assert(ng if sid==0 else nf)%q==0
   cf=math.log2(nf)-72.0*c_val;cg=math.log2(ng)-24.0*c_val-log2_u8
   cfs.append(cf);cgs.append(cg)
   nfas.append(72.0*c_act+cf-(lq if sname=="f" else 0.0));ngas.append(24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0))
  _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,_af,prow["kA"],prow["kB"])
  _pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,_ag,prow["kA"],prow["kB"])
  pr_rats.append((sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts))/prow["pr"])
 cfs.sort();cgs.sort();n=len(cfs)
 sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
 sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2)for i in range(6)]
 assert all(sf6[i]<qb_f[i]for i in range(6))and all(sg6[i]<qb_g[i]for i in range(6))and min(pr_rats)>=10.0
 _prs=sorted(pr_rats)
 return dict(n_enc_live=12,parent_divisibility_live="12/12",Z_basis_det_eq_q="12/12",fvor_cone_rays=ud_chk["rays"],fvor_m0=ud_chk["min_ray_coord"],fvor_decimal_prec=ud_chk["exact_decimal_prec"],sinnott_full_unit_group_exact_O_H=ud_fg,live_bareiss_99_100_all12=True,presupplied_vectors=0,enc_sext_f=sf6,enc_sext_g=sg6,min_headroom_f=round(min(qb_f[i]-sf6[i]for i in range(6)),2),min_headroom_g=round(min(qb_g[i]-sg6[i]for i in range(6)),2),modeled_joint_pr_over_charged_min=round(_prs[0],2),modeled_joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),all_sextiles_below=True,checks_ok=True)

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
