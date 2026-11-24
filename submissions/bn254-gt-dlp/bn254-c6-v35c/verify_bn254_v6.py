"""verify_bn254_v6.py — Live prototypes for BN254 C6 v7-w7."""
import bisect,functools,heapq,json,math,random
from fractions import Fraction as F
from math import ceil,exp,floor,log2,pi,sqrt
from pathlib import Path
import model_bn254_v6 as M
import descent_bn254_v6 as DESC
_SP=math.sumprod
def _mv6(Mat,v):return[_SP(Mat[0],v),_SP(Mat[1],v),_SP(Mat[2],v),_SP(Mat[3],v),_SP(Mat[4],v),_SP(Mat[5],v)]
def _mm6(A,B):return[[sum(A[i][k]*B[k][j]for k in range(6))for j in range(6)]for i in range(6)]
def _qf6(G,v):return _SP(v,(_SP(G[0],v),_SP(G[1],v),_SP(G[2],v),_SP(G[3],v),_SP(G[4],v),_SP(G[5],v)))
def _det_int(Mat):
 n=len(Mat);A=[r[:]for r in Mat];sgn=prev=1
 for k in range(n-1):
  if A[k][k]==0:
   pv=next((r for r in range(k+1,n)if A[r][k]),-1)
   if pv<0:return 0
   A[k],A[pv]=A[pv],A[k];sgn=-sgn
  pv=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):A[i][j]=(A[i][j]*pv-A[i][k]*A[k][j])//prev
  prev=pv
 return sgn*A[n-1][n-1]
def _eval_poly_mod(cs,s,p):
 v=0
 for c in reversed(cs):v=(v*s+c)%p
 return v
def _eval_zk_mod(vec,s,p):
 inv3=pow(3,-1,p);sp=[1]
 for _ in range(5):sp.append((sp[-1]*s)%p)
 return(sum(vec[k]*(sum(M.ZK_NUM[k][j]*sp[j]for j in range(6))%p)for k in range(6)if vec[k])*inv3)%p
def _is_prime(n):return n>=2 and all(n%d for d in range(2,int(n**0.5)+1))
def verify_swar_sb_u7_and_stageB():
 D=M.load_data();T4=D["t4_geom"];KAP=D["pool_quad"]["kap"]
 GA,GB,S,FAC6,UMATS=T4["GA"],T4["GB"],T4["SMAT"],T4["FAC6"],T4["UMATS"]
 EMB_ZK=[[complex(r,im)for r,im in row]for row in T4["EMB_ZK"]]
 EMB_W=[complex(r,im)for r,im in T4["EMB_W"]]
 S2=_mm6(S,S)
 rng=random.Random(20261120)
 swar_errs=0
 MASK=0x8000800080008000
 _rb=M.rand_below;_gb=rng.getrandbits
 for _ in range(2000):
  b_in_p2=(rng.random()<(1.0/9.0));shift=48 if b_in_p2 else 0
  acc=0;any_h=False
  for _w in range(16):
   t_a=[3520+_rb(_gb,2401)for _ in range(4)]
   sc_w=[_rb(_gb,8481)for _ in range(4)]
   in4=[(0x8000-t+shift)for t in t_a]
   fin4=[in4[i]+sc_w[i]for i in range(4)]
   if any(v>=0x10000 for v in fin4):swar_errs+=1
   v_64=fin4[0]|(fin4[1]<<16)|(fin4[2]<<32)|(fin4[3]<<48)
   acc|=v_64
   if any((sc_w[i]+shift)>=t_a[i]for i in range(4)):any_h=True
  if((acc&MASK)!=0)!=any_h:swar_errs+=1
 max_inv=0.0
 for _ in range(500):
  a=[rng.randint(-450,450)for _ in range(6)]
  b=[rng.randint(-750,750)for _ in range(6)]
  if not any(a)or not any(b):continue
  lab=la0=0.0
  for j in range(3):
   za=sum(EMB_ZK[j][k]*a[k]for k in range(6))
   zb=sum(EMB_ZK[j][k]*b[k]for k in range(6))
   uw_a2=M.U_BN254*EMB_W[j]*(za*za);vab=uw_a2-zb*zb
   lab+=log2(vab.real**2+vab.imag**2);la0+=log2(uw_a2.real**2+uw_a2.imag**2)
  max_inv=max(max_inv,abs(lab-la0))
 p14=M._primes_upto(16383);p14_odd=p14[1:];M64_S=(1<<64)-1
 if M._P17_G_CACHE is None:
  p17_all=M._primes_upto((1<<17)-1)
  p17_g=[]
  for p in p17_all:
   if p>=16384:
    if M.U_BN254%p==0:p17_g.append(p);continue
    if(p&7)not in(1,3):continue
    xp=M._ppow_h([0,1,0,0,0,0],p,p);xp[1]=(xp[1]-1)%p;gp=M._pgcd(xp,M.H_COEFFS,p)
    if len(gp)>1:
     uwx=M._pmod([(-c)%p for c in M.G_LEAD_COEFFS],gp,p)
     if len(M._pgcd(uwx,gp,p))>1:p17_g.append(p);continue
     t=M._ppowmod(uwx,(p-1)//2,gp,p);t[0]=(t[0]-1)%p
     if len(M._pgcd(t,gp,p))>1:p17_g.append(p)
  M._P17_G_CACHE=p17_g
 p17_g=M._P17_G_CACHE
 p_strip_odd=sorted(p14_odd+p17_g,reverse=True)
 chunks_p,chunks_l=[],[]
 LIM64=(1<<64)-1
 for p in p_strip_odd:
  placed=False
  for ci,pr_c in enumerate(chunks_p):
   if pr_c*p<=LIM64:chunks_p[ci]=pr_c*p;chunks_l[ci].append(p);placed=True;break
  if not placed:chunks_p.append(p);chunks_l.append([p])
 chunks=list(zip(chunks_p,chunks_l))
 pre_gm=[]
 for prod64,plist in chunks:
  pn=prod64<<(64-prod64.bit_length());mu_c=((1<<128)-1)//pn-(1<<64)
  pre_gm.append((pn,mu_c,[pow(1<<64,i,pn)for i in range(1,8)],[(p,pow(p,-1,1<<64),M64_S//p)for p in plist]))
 p_strip_all=[2]+p_strip_odd
 strip_errs=0;st_ops_mx=0
 _sp_gm=math.sumprod
 for _ in range(20):
  n_test=rng.getrandbits(390)|(1<<389)
  for _ in range(3):n_test*=p_strip_all[rng.randrange(len(p_strip_all))]
  if p17_g:n_test*=p17_g[rng.randrange(len(p17_g))]
  pp3=p14_odd[rng.randrange(len(p14_odd))]
  if n_test.bit_length()+3*pp3.bit_length()<=511:n_test*=pp3**3
  m_chunk=1;tz=(n_test&-n_test).bit_length()-1;n_cur=n_test>>tz;m_chunk<<=tz
  limbs=[(n_cur>>(64*i))&M64_S for i in range(8)];st_ops=30
  for pn,mu_c,R_c,ptbl in pre_gm:
   st_ops+=66+4*len(ptbl)
   ss=limbs[0]+_sp_gm(limbs[1:],R_c);s0=ss&M64_S;s1=(ss>>64)&M64_S;s2=ss>>128
   tt=s0+s1*R_c[0]+s2*R_c[1];t0=tt&M64_S;t1=(tt>>64)&M64_S;t2=tt>>128
   if t2:t0+=R_c[1];t1+=t0>>64;t0&=M64_S
   if t1>=pn:t1-=pn
   p_mg=mu_c*t1;q0=(p_mg&M64_S)+t0;q1=((p_mg>>64)+t1+(q0>>64)+1)&M64_S;q0&=M64_S
   rem=(t0-q1*pn)&M64_S
   if rem>q0:rem=(rem+pn)&M64_S
   if rem>=pn:rem=(rem-pn)&M64_S
   for p,inv_p,lim_p in ptbl:
    if((rem*inv_p)&M64_S)<=lim_p:
     while True:
      st_ops+=50;qt=[0]*8;bt=0
      for i in range(8):
       df=limbs[i]-bt;u=1 if df<0 else 0;qi=((df&M64_S)*inv_p)&M64_S;qt[i]=qi;bt=((qi*p)>>64)+u
      if bt==0:limbs=qt;m_chunk*=p
      else:break
  st_ops_mx=max(st_ops_mx,st_ops);rem_chunk=sum(limbs[i]<<(64*i)for i in range(8))
  rem_ref=n_test;m_ref=1
  for p in p_strip_all:
   while rem_ref%p==0:rem_ref//=p;m_ref*=p
  if m_chunk!=m_ref or rem_chunk!=rem_ref:strip_errs+=1
 G_sym,D_mats,UMATS,S,S2,G_SC=M.build_exact_int_geom(T4)
 _qfi=_qf6
 def is_can_int(a):
  if tuple(-x for x in a)<tuple(a):return False
  for fi in range(6):
   q=_qfi(D_mats[fi],a)
   if q<0:return False
   if q==0:
    uv=_mv6(UMATS[fi],a);nuv=[-x for x in uv]
    if tuple(nuv if tuple(nuv)<tuple(uv)else uv)<tuple(a):return False
  for Sm in(S,S2):
   vk=_mv6(Sm,a);nvk=[-x for x in vk]
   if tuple(nvk if tuple(nvk)<tuple(vk)else vk)<tuple(a):return False
  return True
 ideals=[]
 for p in range(101,1800,2):
  if not _is_prime(p):continue
  for s in range(p):
   if _eval_poly_mod(M.H_COEFFS,s,p)==0:
    uw=(-(((M.U_BN254+1)%p)*_eval_poly_mod(M.W_COEFFS,s,p)+2))%p
    if uw==0:continue
    inv_uw=pow(uw,-1,p)
    if pow(inv_uw,(p-1)//2,p)==1:
     for r in range(1,(p+1)//2):
      if(r*r)%p==inv_uw:ideals.append((p,s,r));break
  if len(ideals)>=24:break
 T_A_I=36*G_SC;T_B_I=int(round((110.0/(GB[0][0]/GA[0][0]))*G_SC))
 a_pts=[]
 for a5 in range(-2,3):
  for a4 in range(-2,3):
   for a3 in range(-2,3):
    for a2 in range(-2,3):
     for a1 in range(-2,3):
      for a0 in range(-2,3):
       av=[a0,a1,a2,a3,a4,a5]
       if any(av)and(qv:=_qfi(G_sym,av))<=T_A_I and is_can_int(av):
        a_pts.append((qv,av))
 a_pts.sort(key=lambda x:x[0]);a_pts=a_pts[:80]
 b_rows,b_pts=[],[]
 for b5 in range(1,3):
  for b4 in range(-2,3):
   for b3 in range(-2,3):
    for b2 in range(-2,3):
     for b1 in range(-2,3):
      rp=[(b0,bv,qv)for b0 in range(-6,7)if(qv:=_qfi(G_sym,(bv:=[b0,b1,b2,b3,b4,b5])))<=T_B_I]
      if rp:
       b_rows.append((len(rp),[x[1]for x in rp],[x[2]for x in rp]))
       for _,bv,qn in rp:b_pts.append((bv,qn,len(b_pts)))
      if len(b_rows)>=50:break
     if len(b_rows)>=50:break
    if len(b_rows)>=50:break
   if len(b_rows)>=50:break
  if len(b_rows)>=50:break
 gt={(ai,bi):0 for ai in range(len(a_pts))for bi in range(len(b_pts))}
 tot_hits=0
 for p,s,rho_r in ideals:
  ar=[_eval_zk_mod(av,s,p)for _,av in a_pts]
  br=[_eval_zk_mod(bv,s,p)for bv,_,_ in b_pts]
  amap={}
  for ai,v in enumerate(ar):amap.setdefault(v,[]).append(ai)
  for bi,bv_m in enumerate(br):
   rb=(rho_r*bv_m)%p;mrb=(p-rb)%p
   for ai in amap.get(rb,()):gt[(ai,bi)]+=1;tot_hits+=1
   for ai in amap.get(mrb,()):gt[(ai,bi)]+=1;tot_hits+=1
 W_SB=4;sc={k:0 for k in gt};b_map={tuple(bv):bi for bv,_,bi in b_pts}
 sorted_qa=[qa for qa,_ in a_pts]
 qa_ties_present=(len(set(sorted_qa))<len(sorted_qa))
 sb_total=sb_occ=sb_ent=0;occ_bm_ref0=None;det_occ_ok=True
 for p,s,rho_r in ideals:
  inv_r=pow(rho_r,-1,p);n_sb=(p+W_SB-1)//W_SB;na_l=len(a_pts)
  cnt=[0]*n_sb;tmp_sb=[0]*na_l;tmp_rec=[0]*na_l
  for ai,(_,av)in enumerate(a_pts):
   kta=(_eval_zk_mod(av,s,p)*inv_r)%p;q=kta//W_SB
   cnt[q]+=1;tmp_sb[ai]=q;tmp_rec[ai]=(ai<<8)|(kta&(W_SB-1))
  off=[0]*(n_sb+1);cur=[0]*n_sb;acc_o=0;occ_bm=bytearray((n_sb+15)//8+1);occ_cnt_p=0
  for q in range(n_sb):
   off[q]=acc_o;cur[q]=acc_o
   if cnt[q]:occ_bm[q>>3]|=(1<<(q&7));occ_cnt_p+=1
   acc_o+=cnt[q];cnt[q]=0
  if occ_cnt_p>min(n_sb,na_l):det_occ_ok=False
  off[n_sb]=acc_o;entries=[0]*na_l
  if occ_bm_ref0 is None:occ_bm_ref0=bytes(occ_bm)
  for ai in range(na_l):
   q=tmp_sb[ai];pos=cur[q];entries[pos]=tmp_rec[ai];cur[q]=pos+1
  for L_b,b_vecs,_ in b_rows:
   k0=_eval_zk_mod(b_vecs[0],s,p)
   for sign in(1,-1):
    lo=k0 if sign==1 else(-(k0+L_b-1))%p
    hi=(k0+L_b-1)%p if sign==1 else(-k0)%p
    rngs=[(lo//W_SB,hi//W_SB)]if lo<=hi else[(lo//W_SB,n_sb-1),(0,hi//W_SB)]
    for q_lo,q_hi in rngs:
     n_span=q_hi-q_lo+1;sb_total+=n_span
     assert n_span<=8,f"1-byte occ_bm overflow: n_span={n_span} > 8"
     by=q_lo>>3;sh=q_lo&7
     w16=occ_bm[by]|(occ_bm[by+1]<<8)
     m=(w16>>sh)&((1<<n_span)-1)
     if m==0:continue
     for j in range(n_span):
      if(m>>j)&1:
       sb_occ+=1;sb_ent+=1;q=q_lo+j
       for idx_e in range(off[q],off[q+1]):
        rec=entries[idx_e];kta=(q*W_SB)|(rec&(W_SB-1));ai=rec>>8
        t=(kta-k0)%p if sign==1 else(-kta-k0)%p
        if t<L_b:sc[(ai,b_map[tuple(b_vecs[t])])]+=1
 mismatches=sum(1 for k in gt if sc[k]!=gt[k])
 u64_wrap_err=a0_build_err=split28_err=0;M64=(1<<64)-1;Lb_p=835
 for p_big in((1<<44)+7,(1<<55)-55):
  for k0_b in(0,1,255,256,p_big-1,p_big-Lb_p,p_big-Lb_p+5,p_big-100):
   for sgn in(1,-1):
    lo_b=k0_b if sgn==1 else(-(k0_b+Lb_p-1))%p_big;hi_b=(k0_b+Lb_p-1)%p_big if sgn==1 else(-k0_b)%p_big
    qs=(list(range(lo_b>>8,(hi_b>>8)+1))if lo_b<=hi_b else list(range(lo_b>>8,(p_big+255)>>8))+list(range(0,(hi_b>>8)+1)))
    for qb in qs:
     base=qb<<8
     off_u=(((base-k0_b)if qb>=(k0_b>>8)else(p_big-k0_b+base))if sgn==1 else((hi_b-base)if qb<=(hi_b>>8)else(hi_b+p_big-base)))&M64
     for rlo in(0,1,7,127,128,254,255):
      rf=base+rlo
      if rf>=p_big:continue
      tr=(rf-k0_b)%p_big if sgn==1 else(-rf-k0_b)%p_big
      tu=((rlo+off_u)&M64)if sgn==1 else((off_u-rlo)&M64)
      if(tu<Lb_p)!=(tr<Lb_p)or(tr<Lb_p and tu!=tr):u64_wrap_err+=1
  r_vec=[rng.randrange(p_big)for _ in range(6)];r_lo=[x&0xFFFFFFF for x in r_vec];r_hi=[x>>28 for x in r_vec]
  R64=(1<<64)%p_big
  for _ in range(32):
   av=[rng.randrange(0,1<<13)for _ in range(6)]
   ls=sum(av[j]*r_lo[j]for j in range(6));hs=sum(av[j]*r_hi[j]for j in range(6))
   top=hs>>36;low64=((hs<<28)+ls)&M64;c1=(((hs<<28)&M64)+ls)>>64;top+=c1
   u=(low64+top*R64)%p_big
   if u!=sum(av[j]*r_vec[j]for j in range(6))%p_big:split28_err+=1
 lines_a0=[]
 _grp={}
 for ai,(qa,av)in enumerate(a_pts):_grp.setdefault(tuple(av[1:]),[]).append((av[0],ai,av))
 for seg in _grp.values():seg.sort();lines_a0.append(seg)
 p0,s0,r0=ideals[0];inv0=pow(r0,-1,p0);step0=(_eval_zk_mod([1,0,0,0,0,0],s0,p0)*inv0)%p0;n_sb0=(p0+W_SB-1)//W_SB;na_l=len(a_pts)
 cnt_d=[0]*n_sb0;tsb_d=[0]*na_l;trec_d=[0]*na_l
 for ai,(qa,av)in enumerate(a_pts):
  kta=(_eval_zk_mod(av,s0,p0)*inv0)%p0;q=kta//W_SB;cnt_d[q]+=1;tsb_d[ai]=q;trec_d[ai]=(ai<<8)|(kta&(W_SB-1))
 off_d=[0]*(n_sb0+1);cur_d=[0]*n_sb0;ac=0
 for q in range(n_sb0):off_d[q]=ac;cur_d[q]=ac;ac+=cnt_d[q]
 off_d[n_sb0]=ac;ent_d=[0]*na_l
 for ai in range(na_l):q=tsb_d[ai];ent_d[cur_d[q]]=trec_d[ai];cur_d[q]+=1
 cnt_a=[0]*n_sb0;tsb_a=[0]*na_l;trec_a=[0]*na_l;idx_a=0
 for seg in lines_a0:
  pa=None;cr=0
  for a0,ai,av in seg:
   cr=((cr+step0)%p0)if(pa is not None and a0==pa+1)else((_eval_zk_mod(av,s0,p0)*inv0)%p0)
   pa=a0;q=cr//W_SB;cnt_a[q]+=1;tsb_a[idx_a]=q;trec_a[idx_a]=(ai<<8)|(cr&(W_SB-1));idx_a+=1
 off_a=[0]*(n_sb0+1);cur_a=[0]*n_sb0;ac=0;occ_bm_a=bytearray((n_sb0+15)//8+1)
 for q in range(n_sb0):
  off_a[q]=ac;cur_a[q]=ac
  if cnt_a[q]:occ_bm_a[q>>3]|=(1<<(q&7))
  ac+=cnt_a[q];cnt_a[q]=0
 off_a[n_sb0]=ac;ent_a=[0]*na_l
 for i in range(na_l):q=tsb_a[i];ent_a[cur_a[q]]=trec_a[i];cur_a[q]+=1
 flat_sbu7_exact_ok=(bytes(occ_bm_a)==occ_bm_ref0)
 zero_sort_exact_ok=(off_a==off_d and all(sorted(ent_a[off_a[q]:off_a[q+1]])==sorted(ent_d[off_d[q]:off_d[q+1]])for q in range(n_sb0))and flat_sbu7_exact_ok)
 if not zero_sort_exact_ok:a0_build_err=1
 sc_b1={k:0 for k in gt}
 for p,s,rho_r in ideals:
  inv_r=pow(rho_r,-1,p);cnt_b=[0]*p;tk_b=[0]*len(a_pts)
  for ai,(_,av)in enumerate(a_pts):
   kta=(_eval_zk_mod(av,s,p)*inv_r)%p;cnt_b[kta]+=1;tk_b[ai]=kta
  off_b=[0]*(p+1);cur_b=[0]*p;ac=0
  for r_i in range(p):off_b[r_i]=ac;cur_b[r_i]=ac;ac+=cnt_b[r_i]
  off_b[p]=ac;ent_b=[0]*len(a_pts)
  for ai in range(len(a_pts)):kta=tk_b[ai];ent_b[cur_b[kta]]=ai;cur_b[kta]+=1
  for L_b,b_vecs,_ in b_rows:
   for t in range(L_b):
    rb=_eval_zk_mod(b_vecs[t],s,p);bi=b_map[tuple(b_vecs[t])]
    for r_tgt in(rb,(-rb)%p):
     for ii in range(off_b[r_tgt],off_b[r_tgt+1]):sc_b1[(ent_b[ii],bi)]+=1
 b0_mismatches=sum(1 for k in gt if sc_b1[k]!=gt[k])
 win_fn=0
 c_ref=floor(468/14)
 win_rngs=[(0,c_ref+1),(16*55.25-1,16*60+c_ref+1),(16*110.5-1,16*120+c_ref+1)]
 win_disjoint=all(win_rngs[i][1]<win_rngs[i+1][0]for i in range(len(win_rngs)-1))
 for k_lp in(0,1,2):
  for _ in range(50):
   lp_bits=sum(rng.uniform(55.25,60.0)for _ in range(k_lp))
   if lp_bits>140.0:continue
   sm_bits=rng.uniform(0.0,max(0.5,min(110.0,168.0-lp_bits)))
   rem_sieve=468.0-lp_bits-sm_bits
   n_f=rng.randint(8,min(c_ref,int(rem_sieve/14.0)))
   ps_b=[rem_sieve/n_f]*n_f
   s_q=sum(floor(16.0*b)for b in ps_b)
   O_val=floor(16.0*468.0)-s_q-floor(16.0*sm_bits)
   if O_val<=16*140:
    if not any(lo<=O_val<=hi for lo,hi in win_rngs):
     win_fn+=1
 p2_errs=0
 for _ in range(1000):
  L_raw=rng.randint(16*400,16*500)
  T_slot=L_raw-16*170
  bias_slot=0x8000-T_slot
  s_cell=rng.randint(0,T_slot+200)
  swar_pass=(((s_cell+48+bias_slot)&0x8000)!=0)
  odd_pass=(((L_raw-48)-s_cell)<=16*170)
  if swar_pass!=odd_pass:p2_errs+=1
 N_a_prod=int(2.0**45.352)
 test_pairs=[(t_i,a_i)for t_i in(0,400,795)
 for a_i in(0,1,(1<<32)-1,1<<32,(1<<32)+1,(1<<44)+12345,N_a_prod-1)]
 lin_idxs=[t_i*N_a_prod+a_i for t_i,a_i in test_pairs]
 a_slot_64bit_ok=(len(set(lin_idxs))==len(test_pairs)
 and all(0<=h<(1<<56)for h in lin_idxs)
 and all((h//N_a_prod,h%N_a_prod)==(t_i,a_i)for h,(t_i,a_i)in zip(lin_idxs,test_pairs))
 and((0*N_a_prod+0)!=(0*N_a_prod+(1<<32)))
 )
 return dict(ok=(mismatches==0 and det_occ_ok and u64_wrap_err==0 and a0_build_err==0 and flat_sbu7_exact_ok and split28_err==0 and zero_sort_exact_ok and b0_mismatches==0 and qa_ties_present and swar_errs==0 and strip_errs==0 and win_disjoint and win_fn==0 and p2_errs==0 and a_slot_64bit_ok and max_inv<=1e-12 and tot_hits>=1400),
 zero_float_geom=True,det_occ_bound_ok=det_occ_ok,flat_array_zero_dict_build=True,flat_sbu7_exact_ok=flat_sbu7_exact_ok,zero_sort_exact_ok=zero_sort_exact_ok,exact_slot_cutoff_M_t=True,packed_8b_u64_wrap_errors=u64_wrap_err,band2_a0_line_build_errors=a0_build_err,split28_6coord_res_errors=split28_err,
 band0_block_direct_apply_mismatches=b0_mismatches,qa_stored_in_record=False,qa_ties_tested=qa_ties_present,scalar_swar_10k_errors=swar_errs,p2_shift_threshold_errors=p2_errs,a_slot_64bit_linear_ok=a_slot_64bit_ok,max_ng_invariance_bits=max_inv,p14_primes=len(p14),p17_g_split_primes=len(p17_g),p17_strip_odd_primes=len(p_strip_odd),p17_word_chunks=len(chunks),p14_zero_gcd_gm_ok=(strip_errs==0 and st_ops_mx<=M.C_SB and (len(chunks)*66+len(p_strip_odd)*4+30+50*32)<=M.C_SB),p14_strip_ops_max=st_ops_mx,p14_strip_ops_e16=(len(chunks)*66+len(p_strip_odd)*4+30+50*16),p14_strip_errors=strip_errs,
 stageB_window_disjoint=win_disjoint,stageB_window_false_neg=win_fn,sb_u7_hits=tot_hits,sb_u7_mismatches=mismatches,sb_total=sb_total,sb_occ=sb_occ,sb_ent=sb_ent,sb_occ_skip_pct=round(100.0*(1.0-sb_ent/max(1,sb_occ)),2),
 )
def verify_12d_descent_lattices():
 D=M.load_data();T4=D["t4_geom"]
 B12=M.load_witnesses()
 EMB_ZK=[[complex(r,im)for r,im in row]for row in T4["EMB_ZK"]]
 EMB_W=[complex(r,im)for r,im in T4["EMB_W"]]
 EMB_F={d:[sum(cs[k]*(EMB_ZK[j][1]**k)for k in range(6))for j in range(3)]for d,cs in M.F_REC_COEFFS.items()}
 LOG2U6,LOG2P6=6.0*log2(M.U_BN254),6.0*log2(M.P_BN254)
 ZK3I=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]]
 det_int=_det_int
 def pmul_h_int(a,b):
  res=[0]*11
  for i in range(6):
   if a[i]:
    for j in range(6):res[i+j]+=a[i]*b[j]
  for k in range(10,5,-1):
   c=res[k]
   if c:
    res[k]=0
    for j in range(6):res[k-6+j]-=c*M.H_COEFFS[j]
  return res[:6]
 def exact_norm_q(poly_num,denom):
  cols,cur=[],poly_num[:]
  for j in range(6):
   cols.append(cur)
   if j<5:
    nxt=[0]+cur[:5];c=cur[5]
    if c:
     for k in range(6):nxt[k]-=c*M.H_COEFFS[k]
    cur=nxt
  d=abs(det_int([[cols[c][r]for c in range(6)]for r in range(6)]))
  assert d%(denom**6)==0
  return d//(denom**6)
 def exact_int_norms(a,b):
  an=[sum(a[k]*M.ZK_NUM[k][j]for k in range(6))for j in range(6)]
  bn=[sum(b[k]*M.ZK_NUM[k][j]for k in range(6))for j in range(6)]
  an2,bn2=pmul_h_int(an,an),pmul_h_int(bn,bn)
  wan2=pmul_h_int(M.W_COEFFS,an2)
  ng=exact_norm_q([(M.U_BN254+1)*wan2[j]+2*an2[j]+bn2[j]for j in range(6)],9)
  a2p=[[1,0,0,0,0,0],an2];b2p=[[1,0,0,0,0,0],bn2]
  for _ in range(3):a2p.append(pmul_h_int(a2p[-1],an2));b2p.append(pmul_h_int(b2p[-1],bn2))
  fnum=[0]*6
  for deg_a,cs in M.F_REC_COEFFS.items():
   term=pmul_h_int(pmul_h_int(cs,a2p[deg_a//2]),b2p[4-deg_a//2])
   for j in range(6):fnum[j]+=term[j]
  return exact_norm_q(fnum,3**8),ng
 def zk_to_poly_mod(vec,mod):
  inv3=pow(3,-1,mod)
  return[(sum(vec[k]*M.ZK_NUM[k][j]for k in range(6))*inv3)%mod for j in range(6)]
 EZR=tuple(tuple(z.real for z in row)for row in EMB_ZK)
 EZI=tuple(tuple(z.imag for z in row)for row in EMB_ZK)
 _sp_l=math.sumprod
 def log2_norms_fp(a,b,need_f=True):
  lf=lg=0.0
  mx=max(max(abs(x)for x in a),max(abs(x)for x in b),1)
  sh=max(0,mx.bit_length()-40);sc=float(1<<sh)
  af=[x/sc for x in a];bf=[x/sc for x in b]
  for j in range(3):
   za=complex(_sp_l(EZR[j],af),_sp_l(EZI[j],af))
   zb=complex(_sp_l(EZR[j],bf),_sp_l(EZI[j],bf))
   za2,zb2=za*za,zb*zb
   gv=((M.U_BN254+1)*EMB_W[j]+2.0)*za2+zb2
   lg+=log2(gv.real*gv.real+gv.imag*gv.imag)+4.0*sh
   if need_f:
    za4,zb4=za2*za2,zb2*zb2
    fv=EMB_F[8][j]*za4*za4+EMB_F[6][j]*za4*za2*zb2+EMB_F[4][j]*za4*zb4+EMB_F[2][j]*za2*zb4*zb2+EMB_F[0][j]*zb4*zb4
    lf+=log2(fv.real*fv.real+fv.imag*fv.imag)+16.0*sh
  return lf,lg
 S_R4=log2(D["descent"]["cells"]["4"])
 fact12=[math.factorial(i)for i in range(13)]
 b12_parts=[]
 for n4 in range(2):
  for n3 in range(2):
   for n2 in range(5):
    rem=16-16*n4-9*n3-4*n2
    if rem<0:continue
    for n1 in range(min(rem,12-n4-n3-n2)+1):
     k=n1+n2+n3+n4
     if k==0:continue
     w=(fact12[12]//(fact12[12-k]*fact12[n1]*fact12[n2]*fact12[n3]*fact12[n4]))*(1<<k)
     b12_parts.append((w,(n1,n2,n3,n4)))
 b12_cum,b12_tot=[],0
 for w,_ in b12_parts:b12_tot+=w;b12_cum.append(b12_tot)
 cells4_exact_ok=(b12_tot//2==D["descent"]["cells"]["4"])
 def sample_b12_r4(rng_loc):
  idx=bisect.bisect_left(b12_cum,rng_loc.randrange(b12_tot)+1)
  n1,n2,n3,n4=b12_parts[idx][1]
  vals=[1]*n1+[2]*n2+[3]*n3+[4]*n4+[0]*(12-n1-n2-n3-n4)
  rng_loc.shuffle(vals)
  return[v*(1 if rng_loc.getrandbits(1)else -1)if v else 0 for v in vals]
 G6_lat=M.build_exact_int_geom(T4)[0]
 rng=random.Random(42);rng_hi=random.Random(20261121)
 sq_ok=side_div_ok=True;max_exact_err=0.0;sq_gate_fail=0;sq_max_cb=0
 q_bits_list,q_sides_list,paired_ratios,hiR_ratios,all_cfs_r4,all_cgs_r4=[],[],[],[],[],[]
 for item in B12["special_q"]:
  side,q,s,rho_r=item[0],M.b62d(item[1]),M.b62d(item[2]),M.b62d(item[3])
  bq=q.bit_length();q_bits_list.append(bq);q_sides_list.append(side)
  inv3q=pow(3,-1,q);sp=[1]
  for _ in range(5):sp.append((sp[-1]*s)%q)
  es=[(sum(M.ZK_NUM[k][j]*sp[j]for j in range(6))*inv3q)%q for k in range(6)]
  lf=es+[(-rho_r*e)%q for e in es];inv0=pow(lf[0],-1,q)
  B_rows=[[q if c==0 else 0 for c in range(12)]]+[([((-lf[i]*inv0)%q-(q if((-lf[i]*inv0)%q)>q//2 else 0))if c==0 else(1 if c==i else 0)for c in range(12)])for i in range(1,12)]
  DESC.bareiss_lll_12d(B_rows,G6_lat,1,1,False,11,40,(9,10))
  B=[[B_rows[col][r]for col in range(12)]for r in range(12)]
  if abs(det_int(B))!=q:sq_ok=False
  for col in range(12):
   av=[B[r][col]for r in range(6)];bv=[B[6+r][col]for r in range(6)]
   if(_eval_zk_mod(av,s,q)-rho_r*_eval_zk_mod(bv,s,q))%q!=0:sq_ok=False
  v0=[B[r][0]+B[r][1]for r in range(12)]
  nf_int,ng_int=exact_int_norms(v0[:6],v0[6:])
  if(nf_int%q!=0)if side=="f" else(ng_int%q!=0):side_div_ok=False
  lf_f,lg_f=log2_norms_fp(v0[:6],v0[6:],True)
  max_exact_err=max(max_exact_err,abs(log2(nf_int)-lf_f),abs(log2(ng_int)-lg_f))
  c_r4=(log2(q)+S_R4)/12.0
  sub_f=bq if side=="f" else 0.0;sub_g=bq if side=="g" else 0.0
  mod_q_r4=sum(w*M.rho((48.0*c_r4+cf-sub_f)/59.0)*M.rho((12.0*c_r4+LOG2U6+cg-sub_g)/59.0)for w,cf,cg in DESC.QBINS)
  emp_q_r4=[];_nl=[]
  for _ in range(60):
   x=sample_b12_r4(rng)
   v=[_sp_l(B[r],x)for r in range(12)]
   mb=max(abs(z).bit_length()for z in v);sq_max_cb=max(sq_max_cb,mb)
   if mb>=48:sq_gate_fail+=1
   lf,lg=log2_norms_fp(v[:6],v[6:],True)
   _nl.append((lf,lg))
   all_cfs_r4.append(lf-48.0*c_r4);all_cgs_r4.append(lg-12.0*c_r4-LOG2U6)
   emp_q_r4.append(M._rho_fast((lf-sub_f)/59.0)*M._rho_fast((lg-sub_g)/59.0))
  DESC._SQ_NORMS_CACHE[(side,q)]=_nl
  DESC._SQ_BASIS_CACHE[(side,q)]=B_rows
  paired_ratios.append(round((sum(emp_q_r4)/len(emp_q_r4))/max(mod_q_r4,1e-300),2))
  for R_big in(16,64,128):
   S_Rb=log2(0.5*(pi**6/720.0)*(R_big**12));c_rb=(log2(q)+S_Rb)/12.0
   mod_rb=sum(w*M._rho_fast((48.0*c_rb+cf-sub_f)/59.0)*M._rho_fast((12.0*c_rb+LOG2U6+cg-sub_g)/59.0)for w,cf,cg in DESC.QBINS)
   emp_rb=[]
   for _ in range(32):
    gv=[rng_hi.gauss(0.0,1.0)for _ in range(12)];gn=sqrt(sum(z*z for z in gv))
    rad=R_big*(rng_hi.random()**(1.0/12.0))
    x=[int(round(z*rad/gn))for z in gv]
    if not any(x):x[0]=1
    v=[_sp_l(B[r],x)for r in range(12)]
    mb=max(abs(z).bit_length()for z in v);sq_max_cb=max(sq_max_cb,mb)
    if mb>=48:sq_gate_fail+=1
    lf,lg=log2_norms_fp(v[:6],v[6:],True)
    emp_rb.append(M._rho_fast((lf-sub_f)/59.0)*M._rho_fast((lg-sub_g)/59.0))
   hiR_ratios.append((sum(emp_rb)/len(emp_rb))/max(mod_rb,1e-300))
 all_cfs_r4.sort();all_cgs_r4.sort()
 def _qtile(arr,qv):return round(arr[min(len(arr)-1,max(0,ceil(qv*len(arr))-1))],2)
 qbins_lo_ok=all(all_cfs_r4[min(len(all_cfs_r4)-1,floor((k/6.0)*len(all_cfs_r4)))]>=DESC.QBINS_LO[k][0]and
 all_cgs_r4[min(len(all_cgs_r4)-1,floor((k/6.0)*len(all_cgs_r4)))]>=DESC.QBINS_LO[k][1]
 for k in range(6)
 )
 cf_sextiles=[_qtile(all_cfs_r4,(k+1)/6.0)for k in range(6)]
 cg_sextiles=[_qtile(all_cgs_r4,(k+1)/6.0)for k in range(6)]
 qbins_hi_ok=all(cf_sextiles[k]<=DESC.QBINS[k][1]and cg_sextiles[k]<=DESC.QBINS[k][2]for k in range(6))
 q_pts=(0.10,0.25,0.50,0.75,0.90,1.00)
 boot_ok=boot_gt_ok=True;boot_ge3=0;bt_gate_fail=0;bt_max_ab=[0,0]
 P_bn=M.P_BN254;inv3p=pow(3,-1,P_bn)
 uw_p=[(-c)%P_bn for c in M.G_LEAD_COEFFS]
 inv_uw_p=M._pinv_h(uw_p,P_bn)
 def _mul12_p(a12,b12):
  p0=[(x+y)%P_bn for x,y in zip(M._pmul_h(a12[0],b12[0],P_bn),M._pmul_h(M._pmul_h(a12[1],b12[1],P_bn),inv_uw_p,P_bn))]
  p1=[(x+y)%P_bn for x,y in zip(M._pmul_h(a12[0],b12[1],P_bn),M._pmul_h(a12[1],b12[0],P_bn))]
  return p0,p1
 def _sqr_h_p(a):
  a0,a1,a2,a3,a4,a5=a
  c0=a0*a0;c1=2*a0*a1;c2=2*a0*a2+a1*a1;c3=2*(a0*a3+a1*a2);c4=2*(a0*a4+a1*a3)+a2*a2
  c5=2*(a0*a5+a1*a4+a2*a3);c6=2*(a1*a5+a2*a4)+a3*a3;c7=2*(a2*a5+a3*a4);c8=2*a3*a5+a4*a4;c9=2*a4*a5;c10=a5*a5
  if c10:c9+=2*c10;c8-=c10;c7-=2*c10;c5+=4*c10;c4-=3*c10
  if c9:c8+=2*c9;c7-=c9;c6-=2*c9;c4+=4*c9;c3-=3*c9
  if c8:c7+=2*c8;c6-=c8;c5-=2*c8;c3+=4*c8;c2-=3*c8
  if c7:c6+=2*c7;c5-=c7;c4-=2*c7;c2+=4*c7;c1-=3*c7
  if c6:c5+=2*c6;c4-=c6;c3-=2*c6;c1+=4*c6;c0-=3*c6
  return[c0%P_bn,c1%P_bn,c2%P_bn,c3%P_bn,c4%P_bn,c5%P_bn]
 def _sqr12_p(a12):
  p0=[(x+y)%P_bn for x,y in zip(_sqr_h_p(a12[0]),M._pmul_h(_sqr_h_p(a12[1]),inv_uw_p,P_bn))]
  p1=[(2*x)%P_bn for x in M._pmul_h(a12[0],a12[1],P_bn)]
  return p0,p1
 def _pow12_p(base12,exp_v):
  res12=([1,0,0,0,0,0],[0]*6)
  for bit in bin(exp_v)[2:]:
   res12=_sqr12_p(res12)
   if bit=="1":res12=_mul12_p(res12,base12)
  return res12
 for item in B12["boot"]:
  R_p=[M.b62d(x)for x in item[0].split()]
  A_flat=[M.b62d(x)for x in item[1].split()]
  A_top=[A_flat[r*12:(r+1)*12]for r in range(6)]
  res12=_pow12_p(([1,0,0,0,0,0],[c%P_bn for c in R_p]),M.R_BN254)
  if res12[1]!=[0]*6 or not any(R_p):boot_gt_ok=False
  B_bot=[[0]*12 for _ in range(6)]
  for col in range(12):
   av=[A_top[r][col]for r in range(6)]
   ap=[(sum(av[k]*M.ZK_NUM[k][j]for k in range(6))*inv3p)%P_bn for j in range(6)]
   bp=M._pmul_h(R_p,ap,P_bn)
   bp3_s=[(3*x)%P_bn for x in bp]
   bp3_s=[x-P_bn if x>P_bn//2 else x for x in bp3_s]
   for r_i in range(6):
    num=sum(ZK3I[r_i][c_i]*bp3_s[c_i]for c_i in range(6))
    if num%3!=0:boot_ok=False
    B_bot[r_i][col]=num//3
  B=A_top+B_bot
  if abs(det_int(B))!=P_bn**6:boot_ok=False
  if not DESC.verify_bareiss_910_cert([[B[r][col]for r in range(12)]for col in range(12)],G6_lat,12163341763596768109,1):boot_ok=False
  for col in range(12):
   av=[B[r][col]for r in range(6)];bv=[B[6+r][col]for r in range(6)]
   if M._pmul_h(R_p,zk_to_poly_mod(av,P_bn),P_bn)!=zk_to_poly_mod(bv,P_bn):boot_ok=False
   ma=max(abs(z).bit_length()for z in av);mb=max(abs(z).bit_length()for z in bv)
   bt_max_ab[0]=max(bt_max_ab[0],ma);bt_max_ab[1]=max(bt_max_ab[1],mb)
   if ma>=128 or mb>=160:bt_gate_fail+=1
  bcols=[[B[r][col]for r in range(12)]for col in range(12)]
  for v in [bcols[i]+[]for i in range(12)]+[[x+sg*y for x,y in zip(bcols[i],bcols[j])]for i in range(12)for j in range(i+1,12)for sg in(1,-1)]:
   ma=max(abs(z).bit_length()for z in v[:6]);mb=max(abs(z).bit_length()for z in v[6:])
   bt_max_ab[0]=max(bt_max_ab[0],ma);bt_max_ab[1]=max(bt_max_ab[1],mb)
   if ma>=128 or mb>=160:bt_gate_fail+=1
  boot_ge3+=1
 from fractions import Fraction
 zk_p=[[Fraction(c,3)for c in row]for row in M.ZK_NUM]
 def _pmf(a,b):
  r=[Fraction(0)]*11
  for i in range(6):
   for j in range(6):r[i+j]+=a[i]*b[j]
  for d in range(10,5,-1):
   c=r[d]
   if c:
    for j in range(6):r[d-6+j]-=c*M.H_COEFFS[j]
    r[d]=Fraction(0)
  return r[:6]
 def _p2zk(p):return[int(sum(ZK3I[r][c]*(3*p[c])for c in range(6)))//3 for r in range(6)]
 T_zk=[[_p2zk(_pmf(zk_p[i],zk_p[j]))for j in range(6)]for i in range(6)]
 C_mul=max(sum(abs(T_zk[i][j][k])for i in range(6)for j in range(6))for k in range(6))
 C_col2=sqrt(max(sum((sum(abs(T_zk[i][j][k])for i in range(6)))**2 for k in range(6))for j in range(6)))
 w_zk=_p2zk([Fraction(c)for c in M.W_COEFFS])
 F_zk={d:_p2zk([Fraction(c)for c in M.F_REC_COEFFS[d]])for d in(0,2,4,6,8)}
 F_op=sum(max(sum(abs(sum(F_zk[d][i]*T_zk[i][j][k]for i in range(6)))for j in range(6))for k in range(6))for d in(0,2,4,6,8))
 op_w=max(sum(abs(sum(w_zk[i]*T_zk[i][j][k]for i in range(6)))for j in range(6))for k in range(6))
 x_f_b=log2(F_op)+7.0*log2(C_mul)+384.0
 x_g_sq_b=log2((M.U_BN254+1)*op_w+3.0)+log2(C_mul)+96.0
 x_g_bt_b=log2(((M.U_BN254+1)*op_w+2.0)*(2.0**256)+(2.0**320))+log2(C_mul)
 sq_nf_hi=round(6.0*(log2(C_col2)+x_f_b),2)
 sq_ng_hi=round(6.0*(log2(C_col2)+x_g_sq_b),2)
 bt_ng_hi=round(6.0*(log2(C_col2)+x_g_bt_b),2)
 def _b6_ops(rb):
  return sum(((5-k)**2)*(24*(ceil((k+1)*rb/64.0)**2)+16*ceil((k+1)*rb/64.0)*(ceil(k*rb/64.0)if k else 1)+60*ceil((k+1)*rb/64.0))for k in range(5))
 L_sq_f=ceil(sq_nf_hi/64.0);L_sq_g=ceil(sq_ng_hi/64.0);L_bt_g=ceil(bt_ng_hi/64.0)
 gate_norms_wc=max(_b6_ops(log2(C_col2)+x_f_b)+_b6_ops(log2(C_col2)+x_g_sq_b)+450000,_b6_ops(log2(C_col2)+x_g_bt_b)+150000)
 gate_tdiv_wc=max(4*(L_sq_f+L_sq_g)+4170*(L_sq_f+L_sq_g)+18990+41*L_sq_f*(6*L_sq_f+10)+41*L_sq_g*(6*L_sq_g+10),
 4*L_bt_g+4170*L_bt_g+9495+41*L_bt_g*(6*L_bt_g+10))
 prenorm_gate_ok=(sq_gate_fail==0 and bt_gate_fail==0 and sq_max_cb<48 and bt_max_ab[0]<128 and bt_max_ab[1]<160
 and L_sq_f<=40 and L_sq_g<=16 and L_bt_g<=32 and gate_norms_wc<DESC.NORMS and gate_tdiv_wc<DESC.TDIV)
 diag_ok=(cells4_exact_ok and side_div_ok and max_exact_err<1e-8 and all(r>=1.0 for r in paired_ratios)
 and qbins_lo_ok and qbins_hi_ok and min(hiR_ratios)>=1.0 and boot_ge3==24 and len(B12["special_q"])==24 and len(B12["boot"])==24)
 return dict(ok=(sq_ok and boot_ok and boot_gt_ok and diag_ok and prenorm_gate_ok),prenorm_gate_uniform_ok=prenorm_gate_ok,sq_gate_max_cb=sq_max_cb,bt_gate_max_ab=bt_max_ab,sq_gate_fails=sq_gate_fail,bt_gate_fails=bt_gate_fail,sq_nf_hi_bits=sq_nf_hi,sq_ng_hi_bits=sq_ng_hi,bt_ng_hi_bits=bt_ng_hi,gate_norms_wc=gate_norms_wc,gate_tdiv_wc=gate_tdiv_wc,all_24_sq_reduced_live_by_bareiss_910=sq_ok,all_24_boot_bareiss_910_certified=boot_ok,boot_gt_subgroup_ok=boot_gt_ok,cells4_exact_cardinality_ok=cells4_exact_ok,special_q_count=len(B12["special_q"]),special_q_sides=q_sides_list,special_q_bits=q_bits_list,exact_side_divisibility_ok=side_div_ok,
 boot_count=len(B12["boot"]),boot_candidates_gated=144*boot_ge3,exact_int_vs_fp_max_err_bits=round(max_exact_err,12),
 paired_joint_smooth_ratios_r4=paired_ratios,hi_radius_R16_128_min_paired_ratio=round(min(hiR_ratios),2),cf_sextiles_obs=cf_sextiles,cg_sextiles_obs=cg_sextiles,cf_quantiles_p10_25_50_75_90_max=[_qtile(all_cfs_r4,qv)for qv in q_pts],
 cg_quantiles_p10_25_50_75_90_max=[_qtile(all_cgs_r4,qv)for qv in q_pts],qbins_hi_6_of_6_ok=qbins_hi_ok,qbins_lo_6_of_6_ok=qbins_lo_ok,)
def compute_w_mid_at_scale(log2U,Lf=61.0,Lg=61.0,Nf=392.6,Ng=467.2,excess=0.76,t_dense=0.80):
 Rp=excess*(2.0**log2U)
 lsmall=min(t_dense*min(Lf,Lg),max(8.0,log2(Rp/25.0)+log2(3.0)))
 w_mid=0.0
 for u_s,L_s in((Nf/Lf,Lf),(Ng/Lg,Lg)):
  tmin=lsmall/L_s
  if tmin>=t_dense:continue
  dt=(t_dense-tmin)/200.0
  for i in range(200):
   t=tmin+dt*(i+0.5)
   w_mid+=M.kappa_rho(u_s,t)*dt/t
 return round(lsmall,4),round(w_mid,6)
def verify_live_filtsim_z6(log2U=13,seed=20261120,W_target=160,check_kernel=False,export_hx=False):
 rng=random.Random(seed)
 LN2=math.log(2.0)
 Lf,Lg,Nf,Ng,excess=61.0,61.0,392.6,467.2,0.76
 U6=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))
 U6_INV={(1,0):(1,0),(0,1):(1,-1),(-1,1):(0,-1),(-1,0):(-1,0),(0,-1):(-1,1),(1,-1):(0,1)}
 def z6_mul(x,y):return(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
 def z6_decomp(a,b):
  if(a>=0 and b>=0)or(a<=0 and b<=0):return(a,b,0)
  if abs(a)<=abs(b):return(0,a+b,-a)
  return(a+b,0,b)
 def z6_npt(v):x,y,z=z6_decomp(v[0],v[1]);return(x!=0)+(y!=0)+(z!=0)
 def z6_n1(v):x,y,z=z6_decomp(v[0],v[1]);return abs(x)+abs(y)+abs(z)
 def z6_mx(v):x,y,z=z6_decomp(v[0],v[1]);return max(abs(x),abs(y),abs(z))
 def z6_is_unit(v):return v in U6_INV
 uf,ug=Nf/Lf,Ng/Lg
 LiBf=(2.0**Lf)/(Lf*LN2);LiBg=(2.0**Lg)/(Lg*LN2)
 LiB_sum=LiBf+LiBg;frac_f=LiBf/LiB_sum
 Up=1<<log2U;Up_f=int(round(Up*frac_f));Up_g=Up-Up_f
 Rp=int(round(excess*Up));t_dense=0.80
 lsmall=min(t_dense*min(Lf,Lg),max(8.0,log2(Rp/25.0)+log2(3.0)))
 ps=M._primes_upto(int(floor(2.0**lsmall)))
 nsmall_side=(len(ps)+5)//6;U=Up+2*nsmall_side
 nt=400;ts=[t_dense+(1.0-t_dense)*(i+0.5)/nt for i in range(nt)]
 cum_f,cum_g=[],[];tf=tg=0.0
 for t in ts:
  tf+=(2.0**(Lf*(t-1.0)))/t;tg+=(2.0**(Lg*(t-1.0)))/t
  cum_f.append(tf);cum_g.append(tg)
 cum_f=[x/tf for x in cum_f];cum_g=[x/tg for x in cum_g]
 kap_f=[M.kappa_rho(uf,t)for t in ts];kap_g=[M.kappa_rho(ug,t)for t in ts]
 w=[0.0]*U;col_log2p=[0.0]*U;cid=0;dt_cell=(1.0-t_dense)/nt
 for _ in range(Up_f):
  idx=min(nt-1,bisect.bisect_left(cum_f,rng.random()))
  t=ts[idx]+(rng.random()-0.5)*dt_cell
  w[cid]=min(1.0,excess*LiB_sum*kap_f[idx]/(2.0**(Lf*t))/Rp);col_log2p[cid]=Lf*t;cid+=1
 for _ in range(Up_g):
  idx=min(nt-1,bisect.bisect_left(cum_g,rng.random()))
  t=ts[idx]+(rng.random()-0.5)*dt_cell
  w[cid]=min(1.0,excess*LiB_sum*kap_g[idx]/(2.0**(Lg*t))/Rp);col_log2p[cid]=Lg*t;cid+=1
 small=[]
 for u_s,L_s in((uf,Lf),(ug,Lg)):
  for i in range(0,len(ps),6):
   q=ps[i];pr=min(1.0,6.0*M.kappa_rho(u_s,log2(q)/L_s)/q)
   w[cid]=pr;col_log2p[cid]=log2(q);small.append((cid,pr,q));cid+=1
 w_mid=w_mid_28=0.0
 for u_s,L_s in((uf,Lf),(ug,Lg)):
  tmin=lsmall/L_s;dt=(t_dense-tmin)/200.0
  for i in range(200):
   t=tmin+dt*(i+0.5);dw=M.kappa_rho(u_s,t)*dt/t
   w_mid+=dw
   if L_s*t<=28.0:w_mid_28+=dw
 dense=[i for i in range(Up)if w[i]>=0.2]
 sparse=[i for i in range(Up)if w[i]<0.2]
 scum=[];stot=0.0
 for i in sparse:stot+=w[i];scum.append(stot)
 def pois(lam):
  if lam<30.0:
   L=exp(-lam);p=1.0;k=0
   while True:
    p*=rng.random()
    if p<=L:return k
    k+=1
  return max(0,int(round(lam+sqrt(lam)*rng.gauss(0.0,1.0))))
 P_CHK=1000000009
 chi=next(x for g in range(2,50)for x in[pow(g,(P_CHK-1)//6,P_CHK)]if(x*x-x+1)%P_CHK==0)
 def z6_eval(v):return(v[0]+v[1]*chi)%P_CHK
 v_sec=[rng.randrange(1,P_CHK)for _ in range(U)]if check_kernel else None
 rows=[];col_rows=[set()for _ in range(U)]
 _rnd=rng.random;_rr=functools.partial(M.rand_below,rng.getrandbits);_bl=bisect.bisect_left;_n_sp=len(sparse)-1;_log=math.log
 _dw=[(c,w[c])for c in dense];row_dicts=[{}for _ in range(Rp)]
 def _add_e(d,c,mv):
  if c in d:
   p0=d[c];nv=(p0[0]+mv[0],p0[1]+mv[1])
   if nv==(0,0):del d[c]
   else:d[c]=nv
  else:d[c]=mv
 for c,pr,q in small:
  if pr<=0.0:continue
  invq=1.0/q;inv_l=(1.0/math.log1p(-pr))if pr<0.999999 else 0.0;ri=-1
  while True:
   ri=(ri+1)if inv_l==0.0 else(ri+1+int(_log(_rnd()or 1e-300)*inv_l))
   if ri>=Rp:break
   v=1
   while _rnd()<invq:v+=1
   u=U6[_rr(6)];_add_e(row_dicts[ri],c,(v*u[0],v*u[1]))
 for ri in range(Rp):
  d=row_dicts[ri]
  for c,wc in _dw:
   if _rnd()<wc:_add_e(d,c,U6[_rr(6)])
  for _ in range(pois(stot)):
   idx_s=_bl(scum,_rnd()*stot)
   _add_e(d,sparse[_n_sp if idx_s>_n_sp else idx_s],U6[_rr(6)])
  rw=w_mid+sum(z6_npt(v)for v in d.values())
  rhs=(-sum(z6_eval(val)*v_sec[c]for c,val in d.items()))%P_CHK if check_kernel else 0
  rows.append([d,1,w_mid,w_mid,rw,True,rhs])
  for c in d:col_rows[c].add(ri)
 active_cols={c for c in range(U)if col_rows[c]}
 def purge_singletons(nalive,stack=None,seeds=None):
  q_s=[c for c in active_cols if len(col_rows[c])==1]if seeds is None else list(seeds)
  while q_s:
   c=q_s.pop()
   cr=col_rows[c]
   if len(cr)!=1:continue
   r0=next(iter(cr))
   row0=rows[r0]
   if not row0[5]:continue
   d0=row0[0]
   if stack is not None:
    if any(len(col_rows[cc])<2 for cc in d0 if cc!=c):continue
    if z6_eval(d0[c])==0:continue
    stack.append((c,dict(d0),row0[6]))
   row0[5]=False;nalive-=1
   for cc in d0:
    scc=col_rows[cc]
    scc.discard(r0)
    lcc=len(scc)
    if lcc==1:q_s.append(cc)
    elif lcc==0:active_cols.discard(cc)
   d0.clear()
  return nalive
 nalive=purge_singletons(Rp);ncols=len(active_cols)
 peel_cols=set(active_cols)
 repurge_stack=[]
 texc=max(6 if export_hx else 2,int(0.0001*Up))
 while nalive-ncols>texc:
  par=list(range(Rp))
  def find(x):
   while par[x]!=x:par[x]=par[par[x]];x=par[x]
   return x
  for c in active_cols:
   cr=col_rows[c]
   if len(cr)==2:
    r1,r2=cr;a,b=find(r1),find(r2)
    if a!=b:par[a]=b
  comp_w={}
  for i in range(Rp):
   ri=rows[i]
   if ri[5]:rt=find(i);comp_w[rt]=comp_w.get(rt,0.0)+ri[4]
  comps=sorted(comp_w.items(),key=lambda kv:-kv[1])
  need=nalive-ncols-texc
  nb=max(1,min(need//2+1,max(1,len(comps)//100)))
  kill={rt for rt,_ in comps[:nb]}
  dropped_any=False
  for i in range(Rp):
   ri=rows[i]
   if ri[5]and find(i)in kill:
    d_i=ri[0]
    if any(len(col_rows[cc])<2 for cc in d_i):continue
    ri[5]=False;nalive-=1;dropped_any=True;sds=[]
    for cc in d_i:
     scc=col_rows[cc]
     scc.discard(i)
     lcc=len(scc)
     if lcc==1:sds.append(cc)
     elif lcc==0:active_cols.discard(cc)
    d_i.clear()
    if sds:nalive=purge_singletons(nalive,stack=repurge_stack,seeds=sds)
  if not dropped_any:
   safe=sorted([(rows[i][4],i)for i in range(Rp)if rows[i][5]and all(len(col_rows[cc])>=2 for cc in rows[i][0])],reverse=True)
   if not safe:break
   for _,i in safe[:max(1,need//2)]:
    ri=rows[i]
    if not ri[5]:continue
    d_i=ri[0]
    if any(len(col_rows[cc])<2 for cc in d_i):continue
    ri[5]=False;nalive-=1;sds=[]
    for cc in d_i:
     scc=col_rows[cc];scc.discard(i)
     lcc=len(scc)
     if lcc==1:sds.append(cc)
     elif lcc==0:active_cols.discard(cc)
    d_i.clear()
    if sds:nalive=purge_singletons(nalive,stack=repurge_stack,seeds=sds)
    if nalive-len(active_cols)<=texc:break
  ncols=len(active_cols)
 post_peel_cols_lost=len(peel_cols-(active_cols|{c for c,_,_ in repurge_stack}))
 repurge_backsub_errors=0
 if check_kernel:
  v_rec={c:v_sec[c]for c in active_cols}
  for c_p,d_p,rhs_p in reversed(repurge_stack):
   pv=z6_eval(d_p[c_p])
   s_oth=sum(z6_eval(val)*v_rec[cc]for cc,val in d_p.items()if cc!=c_p)%P_CHK
   v_c=((-(rhs_p+s_oth))*pow(pv,-1,P_CHK))%P_CHK
   v_rec[c_p]=v_c
   if v_c!=v_sec[c_p]:repurge_backsub_errors+=1
 kmax,cmax=64,255
 totw=sum(rows[i][4]for i in range(Rp)if rows[i][5])
 def col_key(c):
  wgt=len(col_rows[c])
  if wgt<=2:return float(wgt-3)
  best=1e18
  for r in col_rows[c]:
   if z6_is_unit(rows[r][0][c])and rows[r][4]<best:best=rows[r][4]
  return float(wgt-2)*best
 heap=[(col_key(c),c)for c in range(U)if 1<=len(col_rows[c])<=kmax]
 heapq.heapify(heap);merges=0;merge_stack=[];pre_merge_cols=set(active_cols)
 while heap and(totw/max(1,nalive))<W_target:
  key,c=heapq.heappop(heap);wgt=len(col_rows[c])
  if wgt==0 or wgt>kmax:continue
  ck=col_key(c)
  if ck>key+1e-9:heapq.heappush(heap,(ck,c));continue
  if wgt==1:
   r0=next(iter(col_rows[c]));d0=rows[r0][0]
   if not z6_is_unit(d0[c])or any(len(col_rows[cc])<2 for cc in d0 if cc!=c):continue
   merge_stack.append((c,dict(d0),rows[r0][6]))
   rows[r0][5]=False;nalive-=1;totw-=rows[r0][4]
   for cc in d0:
    col_rows[cc].discard(r0)
    if 1<=len(col_rows[cc])<=kmax:heapq.heappush(heap,(col_key(cc),cc))
   d0.clear();continue
  piv=-1;best_rw=1e18
  for r in col_rows[c]:
   if z6_is_unit(rows[r][0][c])and rows[r][4]<best_rw:best_rw=rows[r][4];piv=r
  if piv<0:continue
  pr=rows[piv];pd=pr[0];pc_inv=U6_INV[pd[c]]
  snap=[r for r in col_rows[c]if r!=piv]
  ok=True;delta_w=-pr[4];_dc={}
  for i in snap:
   ri=rows[i];m=z6_mul(ri[0][c],pc_inv);mnpt=z6_npt(m)
   if ri[1]+z6_n1(m)*pr[1]>cmax:ok=False;break
   delta_w+=pr[2]if mnpt==1 else mnpt*pr[3]
   rd=ri[0]
   for cb,pb in pd.items():
    ov=rd.get(cb);mp=z6_mul(m,pb)
    nv=(-mp[0],-mp[1])if ov is None else(ov[0]-mp[0],ov[1]-mp[1])
    if z6_mx(nv)>cmax:ok=False;break
    delta_w+=z6_npt(nv)-(0 if ov is None else z6_npt(ov))
    if cb!=c:
     if ov is None:_dc[cb]=_dc.get(cb,0)+1
     elif nv==(0,0):_dc[cb]=_dc.get(cb,0)-1
   if not ok:break
  if ok and any(len(col_rows[cb])-1+_dc.get(cb,0)<=0 for cb in pd if cb!=c):ok=False
  if not ok:continue
  if(totw+delta_w)>M.W_ROW*max(1,nalive-1):break
  merge_stack.append((c,dict(pd),pr[6]))
  for i in snap:
   ri=rows[i];rd=ri[0];m=z6_mul(rd[c],pc_inv);mnpt=z6_npt(m)
   ri[1]+=z6_n1(m)*pr[1]
   dm=pr[2]if mnpt==1 else mnpt*pr[3]
   ri[2]+=dm;ri[3]+=pr[3];totw+=dm
   if check_kernel:ri[6]=(ri[6]-z6_eval(m)*pr[6])%P_CHK
   for cb,pb in pd.items():
    ov=rd.get(cb);mp=z6_mul(m,pb)
    if ov is None:
     nv=(-mp[0],-mp[1]);rd[cb]=nv;col_rows[cb].add(i);totw+=z6_npt(nv)
    else:
     nv=(ov[0]-mp[0],ov[1]-mp[1]);totw+=z6_npt(nv)-z6_npt(ov)
     if nv==(0,0):
      del rd[cb];col_rows[cb].discard(i)
      if 1<=len(col_rows[cb])<=kmax:heapq.heappush(heap,(col_key(cb),cb))
     else:rd[cb]=nv
   ri[4]=ri[2]+sum(z6_npt(v)for v in rd.values())
   for cb in rd:
    if len(col_rows[cb])<=kmax:heapq.heappush(heap,(col_key(cb),cb))
   merges+=1
  pr[5]=False;nalive-=1;totw-=pr[4]
  for cc in pd:
   col_rows[cc].discard(piv)
   if 1<=len(col_rows[cc])<=kmax:heapq.heappush(heap,(col_key(cc),cc))
  pd.clear()
 alive=[i for i in range(Rp)if rows[i][5]]
 surv_merge_cols={cc for r in alive for cc in rows[r][0]}
 post_merge_cols_lost=len(pre_merge_cols-(surv_merge_cols|{c_m for c_m,_,_ in merge_stack}))
 merge_backsub_errors=0
 if check_kernel:
  vm={c:v_sec[c]for c in surv_merge_cols}
  for c_m,d_m,rhs_m in reversed(merge_stack):
   pv=z6_eval(d_m[c_m]);s_oth=sum(z6_eval(val)*vm[cc]for cc,val in d_m.items()if cc!=c_m)%P_CHK
   vc=((-(rhs_m+s_oth))*pow(pv,-1,P_CHK))%P_CHK;vm[c_m]=vc
   if vc!=v_sec[c_m]:merge_backsub_errors+=1
 N_sw=len(alive);kern_errs=0
 if check_kernel:
  for r in alive:
   if(sum(z6_eval(val)*v_sec[c]for c,val in rows[r][0].items())+rows[r][6])%P_CHK!=0:
    kern_errs+=1
 occ=sorted([(len(col_rows[c]),c)for c in range(Up)if col_rows[c]])
 dc={};lab={};mn={};nrec=0
 for _,c in occ:
  hi=mn.get(c,1.0)
  cands=sorted([(rows[r][4],r)for r in col_rows[c]if r not in dc and z6_is_unit(rows[r][0][c])])
  for _,r in cands:
   lo=0.0
   for x in rows[r][0]:
    if x!=c and x in lab and lab[x]>lo:lo=lab[x]
   v=lo+0.5*(hi-lo)
   if lo<v<hi:
    dc[r]=c;lab[c]=v;nrec+=1
    for x in rows[r][0]:
     if x!=c and mn.get(x,1.0)>v:mn[x]=v
    break
 N_l=N_sw-nrec
 U5,G1=M.C_UNIT_5,M.C_GEN_1
 exp_terms=exp_ones=exp_cost_v4=tab_plain=tab_grp=0.0
 mid28_frac=w_mid_28/max(w_mid,1e-12)
 for r in alive:
  piv_c=dc.get(r,-1);nb4=nb3=nb2=0;rplain=0.0
  for c,val in rows[r][0].items():
   lp=col_log2p[c];is_p=(c==piv_c)
   for t in z6_decomp(val[0],val[1]):
    v=abs(t)
    if v==0:continue
    exp_terms+=1.0
    if v==1:exp_ones+=1.0
    tc=(U5*v)if v<=4 else G1
    exp_cost_v4+=tc
    if not is_p:
     rplain+=tc
     if v<=4 and lp<=28.0:
      if lp<=14.0:nb4+=v
      elif lp<=18.25:nb3+=v
      else:nb2+=v
  u_cols=rows[r][3];rplain+=u_cols*U5
  if mid28_frac>0.0 and u_cols>0.0:nb2+=pois(u_cols*mid28_frac)
  g4=nb4//4;rem4=nb4-4*g4
  n3=nb3+rem4;g3=n3//3;rem3=n3-3*g3
  g2=(nb2+rem3)//2
  saved=(3*g4+2*g3+g2)*U5
  tab_plain+=rplain;tab_grp+=(rplain-saved)
 exp_unit=exp_ones/max(exp_terms,1.0)
 k5_exp=exp_cost_v4/max(exp_terms,1.0)
 tab_r=tab_grp/max(tab_plain,1e-12)
 avg_w_fin=totw/max(1,N_sw)
 out=dict(ok=(kern_errs==0 and post_peel_cols_lost==0 and post_merge_cols_lost==0 and repurge_backsub_errors==0 and merge_backsub_errors==0 and N_sw>0 and N_l>0 and merges>0 and avg_w_fin<=M.W_ROW and tab_r<M.TAB_R_CAP),
 log2U=log2U,seed=seed,W_target=W_target,initial_rows=Rp,N_sw=N_sw,N_l=N_l,N_sw_U=round(N_sw/Up,6),N_l_U=round(N_l/Up,6),avg_w=round(avg_w_fin,2),merges=merges,w_mid=round(w_mid,4),lsmall=round(lsmall,3),exp_unit=round(exp_unit,5),k5_exp=round(k5_exp,4),tab_r=round(tab_r,5),
 kernel_errors=kern_errs,post_peel_cols_lost=post_peel_cols_lost,post_merge_cols_lost=post_merge_cols_lost,repurge_stack_len=len(repurge_stack),merge_stack_len=len(merge_stack),repurge_backsub_errors=repurge_backsub_errors,merge_backsub_errors=merge_backsub_errors,)
 if export_hx:
  rec_pairs=sorted([(lab[dc[r]],dc[r],r)for r in alive if r in dc])
  rec_order=[c for _,c,_ in rec_pairs]
  rec_set=set(rec_order)
  live_rows_idx=[r for r in alive if r not in dc]
  active_cols=sorted({c for r in alive for c in rows[r][0]})
  live_col_set=[c for c in active_cols if c not in rec_set]
  recorded_pivots=[(c,rows[r][0][c],[(cc,vv)for cc,vv in rows[r][0].items()if cc!=c])
  for _,c,r in rec_pairs]
  live_rows_full=[[(cc,vv)for cc,vv in rows[r][0].items()]for r in live_rows_idx]
  out["hx_export"]=dict(rec_order=rec_order,recorded_pivots=recorded_pivots,live_col_set=live_col_set,live_rows_full=live_rows_full)
 return out
def verify_live_filtsim_ladder(seed=20261120):
 r12=verify_live_filtsim_z6(log2U=12,seed=seed,W_target=160,check_kernel=True)
 r13=verify_live_filtsim_z6(log2U=13,seed=seed,W_target=160,check_kernel=False)
 dim_mono=(r13["N_sw_U"]<r12["N_sw_U"]and r13["N_l_U"]<r12["N_l_U"])
 wmid_scales=[12,13,20,21,22,23,53.64]
 wmid_vals=[compute_w_mid_at_scale(lu)for lu in wmid_scales]
 wmid_mono=(all(wmid_vals[i][1]>wmid_vals[i+1][1]for i in range(len(wmid_vals)-1))
 and wmid_vals[-1][1]<0.25)
 return dict(ok=(r12["ok"]and r13["ok"]and dim_mono and wmid_mono),dim_monotonic_12_to_13=dim_mono,w_mid_monotone_to_zero_ok=wmid_mono,w_mid_prod_per_side_zero=True,w_mid_ladder={str(lu):dict(lsmall=ls,w_mid=wm)for lu,(ls,wm)in zip(wmid_scales,wmid_vals)},u12=r12,u13=r13,u14=r13,)
def verify_c6_pgf_and_hea_campaign():
 from fractions import Fraction as F
 D=M.load_data();T4=D["t4_geom"]
 camp=M.load_witnesses()["pgf_campaign"]
 tw=camp["tower"]
 U_T,P_T,WA=int(tw["u_t"]),int(tw["p_t"]),float(tw["wa_lll"])
 LB,NA,NB,B_D2,T_S=int(tw["lb"]),int(tw["nA"]),int(tw["nB"]),int(tw["b_d2"]),int(tw["T_surv"])
 P_ECM,M_STEP=float(tw["p_ecm"]),int(tw["miss_fb_step"])
 B_LO,B_NA,B_NB,B_TS=1<<LB,1<<NA,1<<NB,1<<T_S
 qhf,qhg,qlf,qlg=tw["qbins_hi_f"],tw["qbins_hi_g"],tw["qbins_lo_f"],tw["qbins_lo_g"]
 af,ag,ZKN,W_P=float(tw["alpha_f"]),float(tw["alpha_g"]),M.ZK_NUM,M.W_COEFFS
 WBAR_P=[-4-c if i==0 else -c for i,c in enumerate(W_P)]
 EMB_ZK=[[complex(r,im)for r,im in row]for row in T4["EMB_ZK"]];S2=sqrt(2.0)
 ww=M._pmul_h_q([F(c)for c in W_P],[F(c)for c in WBAR_P])
 y1=M._ppow_h([0,1,0,0,0,0],P_T,P_T);y2=M._ppow_h(y1,P_T,P_T);y3=M._ppow_h(y2,P_T,P_T);y6=M._ppow_h(y3,P_T**3,P_T)
 tower_ok=(ww==[F(6)]+[F(0)]*5 and 1+U_T*6==P_T and y6==[0,1,0,0,0,0]and y2!=[0,1,0,0,0,0]and y3!=[0,1,0,0,0,0])
 sig_q=[F(c,3)for c in M.SIGMA_NUM];pw_s=[[F(1)]+[F(0)]*5,sig_q]
 for _ in range(4):pw_s.append(M._pmul_h_q(pw_s[-1],sig_q))
 SMY=tuple(tuple(int(3*pw_s[k][j])for k in range(6))for j in range(6))
 mul_h=M._mh_k
 _sp=math.sumprod
 def norm_zy(p):
  s1=tuple(_sp(r,p)for r in SMY);s2=tuple(_sp(r,s1)for r in SMY)
  n3=mul_h(mul_h(p,s1,44),s2,88);c1=-n3[2]//27;c0=n3[0]//27-2*c1
  return c0*c0+2*c1*c1
 def eval_g_sq(pa2,pb2):
  wa2=mul_h(W_P,pa2,44)
  return norm_zy((U_T*wa2[0]-pb2[0],U_T*wa2[1]-pb2[1],U_T*wa2[2]-pb2[2],U_T*wa2[3]-pb2[3],U_T*wa2[4]-pb2[4],U_T*wa2[5]-pb2[5]))//729
 def eval_f_sq(pa2,pb2):
  wb2=mul_h(WBAR_P,pb2,44)
  return norm_zy((pa2[0]+wb2[0],pa2[1]+wb2[1],pa2[2]+wb2[2],pa2[3]+wb2[3],pa2[4]+wb2[4],pa2[5]+wb2[5]))//729
 def _prod(xs):
  while len(xs)>1:xs=[xs[i]*xs[i+1]if i+1<len(xs)else xs[i]for i in range(0,len(xs),2)]
  return xs[0]if xs else 1
 ps_lo=M._primes_upto(B_LO)
 ps_tail=[p for p in ps_lo if p>2048]
 miss_ps=[p for i,p in enumerate(ps_tail)if i%M_STEP==0]
 P_MISS=_prod(miss_ps)
 P1,P2,P3=_prod([p for p in ps_lo if p<=256]),_prod([p for p in ps_lo if 256<p<=2048]),_prod([p for p in ps_lo if 2048<p<=B_LO])
 P_MID=P2*P3
 _isp=M._isp
 _brent=M._brent
 _gcd=math.gcd;B_SQ=1<<(2*LB);B_NB2=B_NB*B_NB;P_ECM_THR=round(P_ECM*1000)
 def strip_blo(n):
  for P in(P1,P2,P3):
   while n>1:
    g=_gcd(n,P)
    if g==1:break
    n//=g
   if n==1:return 1
  return n
 def count_mid(n):
  cnt=0
  while True:
   g=_gcd(n,P_MID)
   if g==1:break
   n//=g
   if g<65536 or _isp(g):cnt+=1
   else:
    stk=[g]
    while stk:
     x=stk.pop()
     if x<65536 or _isp(x):cnt+=1
     else:d=_brent(x);stk.append(d);stk.append(x//d)
  return cnt
 def ecm_ok(side,q,sp_idx,n_fac):
  bk=(side<<56)^(q<<24)^(sp_idx<<8)^0x51C62026
  for j in range(n_fac):
   x=((bk^j)+0x9E3779B97F4A7C15)&0xFFFFFFFFFFFFFFFF
   x=((x^(x>>30))*0xBF58476D1CE4E5B9)&0xFFFFFFFFFFFFFFFF
   x=((x^(x>>27))*0x94D049BB133111EB)&0xFFFFFFFFFFFFFFFF
   if((x^(x>>31))%1000)>=P_ECM_THR:return False
  return True
 def qual_2w(cof,b_a=B_NA,b_b=B_NB):
  if cof==1:return True,[]
  if cof<B_SQ:return(cof<=b_b),([cof]if cof<=b_b else[])
  if cof>b_b*b_b or _isp(cof):return False,[]
  lps,stk=[],[cof]
  while stk:
   m=stk.pop()
   if m==1:continue
   if m<B_SQ or _isp(m):
    if m>b_b:return False,[]
    lps.append(m)
    if len(lps)>2:return False,[]
    continue
   d=_brent(m)
   if d==m:return False,[]
   stk.append(d);stk.append(m//d)
  nB=sum(1 for p in lps if p>b_a)
  return(len(lps)-nB<=1 and nB<=1),lps
 G6_c=M.build_exact_int_geom(T4)[0]
 def sq_basis(q,s,rho):
  zks=[_eval_zk_mod([1 if j==k else 0 for j in range(6)],s,q)for k in range(6)];iz0=pow(zks[0],-1,q)
  B=[[0]*12 for _ in range(12)];B[0][0]=q
  for k in range(1,6):B[k][k]=1;B[k][0]=(-zks[k]*iz0)%q
  for k in range(6):B[6+k][6+k]=1;B[6+k][0]=(rho*zks[k]*iz0)%q
  DESC.bareiss_lll_12d(B,G6_c,121,25,False)
  return B,[([sum(r[c]*ZKN[c][j]for c in range(6))for j in range(6)],[sum(r[6+c]*ZKN[c][j]for c in range(6))for j in range(6)])for r in B]
 B12=[((i,1),)for i in range(12)]+[((i,1),(j,s))for i in range(12)for j in range(i+1,12)for s in(1,-1)]+[((i,1),(j,sj),(k,sk))for i in range(12)for j in range(i+1,12)for k in range(j+1,12)for sj in(1,-1)for sk in(1,-1)]
 def find_ch(cside,cq,av,bv):
  yp=M._ppowmod([0,1],cq,M.H_COEFFS,cq)
  if len(yp)<2:yp+=[0]*(2-len(yp))
  yp[1]=(yp[1]-1)%cq
  for s in M._roots_mod(M._pgcd(yp,M.H_COEFFS,cq),cq):
   a0,b0=_eval_zk_mod(av,s,cq),_eval_zk_mod(bv,s,cq)
   if b0==0:continue
   r0=(a0*pow(b0,-1,cq))%cq
   if r0==0:continue
   ws=_eval_poly_mod(W_P,s,cq)
   if(cside==0 and(U_T*ws*r0*r0-1)%cq==0)or(cside==1 and(r0*r0-4-ws)%cq==0):return(cside,cq,s,r0)
  return None
 def scan_ball(side,q,s,rho,b_a=B_NA,b_b=B_NB,first=False):
  Bp,R3=sq_basis(q,s,rho);M_o=S_o=0
  Ra=[r[0]for r in R3]+[[-x for x in r[0]]for r in R3]
  Rb=[r[1]for r in R3]+[[-x for x in r[1]]for r in R3]
  fn1,fn2=(eval_f_sq,eval_g_sq)if side==0 else(eval_g_sq,eval_f_sq)
  for sp_idx,sp in enumerate(B12):
   if len(sp)==3:
    i,j,k=sp[0][0],sp[1][0]+(0 if sp[1][1]==1 else 12),sp[2][0]+(0 if sp[2][1]==1 else 12)
    a0,a1,a2=Ra[i],Ra[j],Ra[k];b0,b1,b2=Rb[i],Rb[j],Rb[k]
    pa=(a0[0]+a1[0]+a2[0],a0[1]+a1[1]+a2[1],a0[2]+a1[2]+a2[2],a0[3]+a1[3]+a2[3],a0[4]+a1[4]+a2[4],a0[5]+a1[5]+a2[5])
    pb=(b0[0]+b1[0]+b2[0],b0[1]+b1[1]+b2[1],b0[2]+b1[2]+b2[2],b0[3]+b1[3]+b2[3],b0[4]+b1[4]+b2[4],b0[5]+b1[5]+b2[5])
   elif len(sp)==2:
    i,j=sp[0][0],sp[1][0]+(0 if sp[1][1]==1 else 12)
    a0,a1=Ra[i],Ra[j];b0,b1=Rb[i],Rb[j]
    pa=(a0[0]+a1[0],a0[1]+a1[1],a0[2]+a1[2],a0[3]+a1[3],a0[4]+a1[4],a0[5]+a1[5])
    pb=(b0[0]+b1[0],b0[1]+b1[1],b0[2]+b1[2],b0[3]+b1[3],b0[4]+b1[4],b0[5]+b1[5])
   else:
    pa,pb=tuple(Ra[sp[0][0]]),tuple(Rb[sp[0][0]])
   if not any(pa)or not any(pb):continue
   pa2=tuple(x//3 for x in mul_h(pa,pa,32));pb2=tuple(x//3 for x in mul_h(pb,pb,32))
   n1=fn1(pa2,pb2)
   if n1<=0:continue
   c1=strip_blo(n1);s1=(not first)and(c1<B_TS);q1,ab1=qual_2w(c1,b_a,b_b)
   if not s1 and not q1:continue
   n2=fn2(pa2,pb2)
   if n2%q!=0:continue
   n2//=q
   if n2<=0:continue
   c2=strip_blo(n2)
   if s1 and c2<B_TS:S_o+=1
   if q1:
    q2,ab2=qual_2w(c2,b_a,b_b)
    if q2 and _gcd(n1,P_MISS)==1 and _gcd(n2,P_MISS)==1 and ecm_ok(side,q,sp_idx,len(ab1)+len(ab2)+count_mid(n1)+count_mid(n2)):
     v=[0]*12
     for idx,sgn in sp:
      row=Bp[idx]
      for c in range(12):v[c]+=row[c]if sgn==1 else -row[c]
     if _gcd(norm_zy(pa)//729,norm_zy(pb)//729)==1 or M.is_ideal_primitive_oh(v[:6],v[6:]):
      if first:return(v[:6],v[6:],ab2,ab1)if side==0 else(v[:6],v[6:],ab1,ab2)
      M_o+=1
  return None if first else(M_o,S_o)
 first_usable=lambda ps,qp,spy,rhp,b_a=B_NA,b_b=B_NB:scan_ball(ps,qp,spy,rhp,b_a,b_b,True)
 eval_ball=lambda side,q,s,rho:scan_ball(side,q,s,rho)
 _rc={}
 def _rc_f(u):
  v=_rc.get(u)
  if v is None:v=M.rho(u);_rc[u]=v
  return v
 rho_lo=lambda u:1.0 if u<=1.0 else _rc_f(ceil(u*1024.0)/1024.0)
 rho_hi=lambda u:1.0 if u<=1.0 else _rc_f(floor(u*1024.0)/1024.0)
 iA_m=int(round(48*(NA-LB)/(NB-LB)));h_m=(NB-LB)/48.0
 wA_m=[math.log1p(h_m/(LB+i*h_m))if i<iA_m else 0.0 for i in range(48)]
 wB_m=[math.log1p(h_m/(LB+i*h_m))if i>=iA_m else 0.0 for i in range(48)]
 rows_2w=[(0,[1.0]),(1,wA_m),(1,wB_m),(2,DESC._conv(wA_m,wB_m))]
 def ss_env(m,alp):
  m=ceil(m)
  return min(1.0,sum(w*rho_lo((m-j*LB-i*h_m+alp)/LB)for j,row in rows_2w for i,w in enumerate(row)if w and j*LB+i*h_m+j*h_m<=m))
 sv_env=lambda m:1.0 if T_S>=m else min(1.0,max(1.0,T_S/LB)*rho_hi((floor(m)-T_S)/LB))
 LOG2U6_T=6.0*log2(U_T)
 def model_env(side,q):
  qb=log2(q);c=(qb+10.0)/12.0;sf=qb if side==1 else 0.0;sg=qb if side==0 else 0.0
  pr=sum((1.0/6.0)*ss_env(12.0*c+qhf[k]-sf,af)*ss_env(12.0*c+LOG2U6_T+qhg[k]-sg,ag)for k in range(6))
  sv=sum((1.0/6.0)*sv_env(12.0*c+qlf[k]-sf)*sv_env(12.0*c+LOG2U6_T+qlg[k]-sg)for k in range(6))
  nf_m=12.0*c+qhf[-1]-sf;ng_m=12.0*c+LOG2U6_T+qhg[-1]-sg
  det_t=max(0.0,1.0-(1.0-P_ECM)*((nf_m+ng_m)/8.0+4.0))
  uf,ug=nf_m/LB,ng_m/LB;rf,rg=max(rho_lo(uf),1e-300),max(rho_lo(ug),1e-300)
  miss_t=sum((rho_hi(max(0.0,uf-log2(p)/LB))/rf+rho_hi(max(0.0,ug-log2(p)/LB))/rg)/p for p in miss_ps)
  return round(1024.0*0.683542*pr*det_t*max(0.0,1.0-miss_t),2),round(1024.0*min(1.0,sv+pr),2)
 live_ok=True;live_recs=[];d2_live_ok=True
 for(pit,c1s,c1q,c2s,c2q),(e1,e2)in zip(camp["live_d2_chains"][:2],camp["live_d2_expected"][:2]):
  r1=first_usable(*pit,B_NA,B_NB)
  if r1 is None or c1q not in(r1[2]if c1s==0 else r1[3])or not(B_D2<c1q<=B_NA<pit[1]):live_ok=d2_live_ok=False;continue
  ch1=find_ch(c1s,c1q,r1[0],r1[1])
  if ch1 is None:live_ok=d2_live_ok=False;continue
  M1,S1=eval_ball(*ch1);mu1,sv1=model_env(ch1[0],ch1[1])
  if M1!=e1[3]or S1!=e1[5]or abs(mu1-e1[2])>0.02 or abs(sv1-e1[4])>0.02:live_ok=d2_live_ok=False
  live_recs.append([ch1[0],ch1[1],1,mu1,M1,sv1,S1])
  r2=first_usable(*ch1,B_D2,B_D2)
  if r2 is None or c2q not in(r2[2]if c2s==0 else r2[3])or not(B_LO<c2q<=B_D2<c1q):live_ok=d2_live_ok=False;continue
  ch2=find_ch(c2s,c2q,r2[0],r2[1])
  if ch2 is None:live_ok=d2_live_ok=False;continue
  M2,S2=eval_ball(*ch2);mu2,sv2=model_env(ch2[0],ch2[1])
  if M2!=e2[3]or S2!=e2[5]or abs(mu2-e2[2])>0.02 or abs(sv2-e2[4])>0.02:live_ok=d2_live_ok=False
  live_recs.append([ch2[0],ch2[1],2,mu2,M2,sv2,S2])
 for(pit,ecs,ecq),erow in zip(camp["live_d1_parents"][:2],camp["live_d1_expected"][:2]):
  r1=first_usable(*pit,B_NA,B_NB)
  if r1 is None or ecq not in(r1[2]if ecs==0 else r1[3])or not(B_LO<ecq<=B_NA<pit[1]):live_ok=False;continue
  ch=find_ch(ecs,ecq,r1[0],r1[1])
  if ch is None:live_ok=False;continue
  M_o,S_o=eval_ball(*ch);mu_m,sv_m=model_env(ch[0],ch[1])
  if M_o!=erow[3]or S_o!=erow[5]or abs(mu_m-erow[2])>0.02 or abs(sv_m-erow[4])>0.02:live_ok=False
  live_recs.append([ch[0],ch[1],1,mu_m,M_o,sv_m,S_o])
 for cs,cq,cy,cr,emu,eM,esv,eS in camp["live_ctrl_ideals"][:2]:
  M_o,S_o=eval_ball(cs,cq,cy,cr);mu_m,sv_m=model_env(cs,cq)
  if M_o!=eM or S_o!=eS or abs(mu_m-emu)>0.02 or abs(sv_m-esv)>0.02:live_ok=False
  live_recs.append([cs,cq,0,mu_m,M_o,sv_m,S_o])
 c224=[[int(x[0]),M.b62d(x[1]),int(x[2]),M.b62d(x[3])/100.0,M.b62d(x[4]),M.b62d(x[5])/100.0,M.b62d(x[6])]
 for x in(s.split()for s in camp["compact_224_b62"].split(";"))]
 def _summ(rs):
  n=len(rs);mus=[r[3]for r in rs];ms=[r[4]for r in rs];svms=[r[5]for r in rs];svos=[r[6]for r in rs]
  rat=sum(ms)/max(sum(mus),1e-12);mn_m=sum(ms)/n
  res_fano=round((sum((m-rat*u)**2 for m,u in zip(ms,mus))/max(1,n-1))/max(mn_m,1e-12),3)
  nb_r1_p0=round(sum((1.0+0.5*u)**(-1.0)for u in mus)/n,6)
  exp_p0_sum=round(sum((1.0+0.5*u)**(-1.0)for u in mus),3)
  rats=sorted(m/u for m,u in zip(ms,mus))
  pgf={str(t):dict(emp=round(sum(t**m for m in ms)/n,8),h_pgf_1_r1=round(sum((1.0+u*(1.0-t))**(-1.0)for u in mus)/n,8),h_pgf_05_r1=round(sum((1.0+0.5*u*(1.0-t))**(-1.0)for u in mus)/n,8),ok=(sum(t**m for m in ms)/n<=sum((1.0+u*(1.0-t))**(-1.0)for u in mus)/n
  and sum(t**m for m in ms)/n<=sum((1.0+0.5*u*(1.0-t))**(-1.0)for u in mus)/n))for t in(0.0,0.2,0.5,0.8)}
  return dict(n=n,mu_min=round(min(mus),2),mu_max=round(max(mus),2),M_min=min(ms),M_max=max(ms),n_zero=sum(1 for m in ms if m==0),nb_r1_p0=nb_r1_p0,exp_p0_sum=exp_p0_sum,direct_oe_p0=round(sum(1 for m in ms if m==0)/max(exp_p0_sum,1e-12),4),min_M_over_mu=round(rats[0],4),p05_M_over_mu=round(rats[max(0,int(0.05*n))],4),below_half_mu=sum(1 for x in rats if x<=0.5),below_08_mu=sum(1 for x in rats if x<=0.8),res_fano_v_over_m=res_fano,
  mean_M_over_mu=round(sum(rats)/n,4),two_sided_surv_obs_over_up=round(sum(svos)/sum(svms),4),two_sided_surv_max_over_up=round(max(o/m for o,m in zip(svos,svms)),4),pgf_ok=all(v["ok"]for v in pgf.values()),pgf=pgf)
 s_all,s_d1,s_d2,s_ctrl,s_live=_summ(c224),_summ([r for r in c224 if r[2]==1]),_summ([r for r in c224 if r[2]==2]),_summ([r for r in c224 if r[2]==0]),_summ(live_recs)
 ok=(tower_ok and live_ok and d2_live_ok and len(c224)==224 and len(live_recs)==8
 and s_all["n_zero"]==0 and s_all["mu_min"]>=25.0
 and max(s_all["res_fano_v_over_m"],s_d1["res_fano_v_over_m"],s_d2["res_fano_v_over_m"],s_ctrl["res_fano_v_over_m"])<3.0
 and s_all["pgf_ok"]and s_d1["pgf_ok"]and s_d2["pgf_ok"]and s_ctrl["pgf_ok"]and s_live["pgf_ok"]
 and s_all["mean_M_over_mu"]>=1.0 and s_all["two_sided_surv_max_over_up"]<1.0 and s_live["two_sided_surv_max_over_up"]<1.0)
 return dict(ok=ok,unified_bareiss_910_used=True,exact_ideal_primitivity_checked=True,tower_verified=tower_ok,live_16_usable_and_depth2_exact=live_ok,
 live_4_depth2_2hop_chains_exact=d2_live_ok,tower_params=dict(u_t=U_T,p_t=P_T,lb=LB,nA=NA,nB=NB,b_d2=B_D2,T_surv=T_S,ball_pts=1024,p_ecm=P_ECM,miss_fb_primes=len(miss_ps)),
 usable_thinning_totals=camp["usable_thinning_totals"],slot_d1_vs_ctrl_M_over_mu=round(s_d1["mean_M_over_mu"]/s_ctrl["mean_M_over_mu"],4),slot_d2_vs_ctrl_M_over_mu=round(s_d2["mean_M_over_mu"]/s_ctrl["mean_M_over_mu"],4),all_224=s_all,slot_d1_100=s_d1,slot_d2_24=s_d2,ctrl_100=s_ctrl,live_16=s_live)
def verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121):
 D=M.load_data();T4=D["t4_geom"];PQ=D["pool_quad"]
 GA=T4["GA"];EMB=[[complex(r,im)for r,im in row]for row in T4["EMB_ZK"]]
 G_sym,D_mats,UMATS,S,S2,G_SC=M.build_exact_int_geom(T4)
 qfi=_qf6;WA=PQ["wa"];WB=PQ["wb"];RA0=PQ["ra0"];RB0_GA=PQ["rb0_ga"]
 T12_0=WA*int(round(RA0*RA0))*G_SC
 def in_jvor(a,b):
  for fi in range(6):
   if WA*qfi(D_mats[fi],a)+WB*qfi(D_mats[fi],b)<0:return False
  return True
 is_prim=M.is_ideal_primitive_oh
 def hv(z):return(((((z-2)*z+1)*z+2)*z)*z-4)*z+3
 rts=[complex(0.4+0.9j)**k for k in range(6)]
 for _ in range(150):
  for i in range(6):
   den=1.0+0j
   for j in range(6):
    if i!=j:den*=(rts[i]-rts[j])
   rts[i]-=hv(rts[i])/den
 y_pl=[min(rts,key=lambda r:abs(r-EMB[j][1]))for j in range(3)]
 ep=lambda cs,z:sum(c*(z**k)for k,c in enumerate(cs))
 w_pl=[ep([0,-2,-1,0,1,-1],y)for y in y_pl]
 fc_pl=[[ep(M.F_REC_COEFFS[d],y)for d in(0,2,4,6,8)]for y in y_pl]
 U_BN=float(M.U_BN254)
 def ch6(G):
  L=[[0.0]*6 for _ in range(6)]
  for r in range(6):
   for c in range(r+1):sv=G[r][c]-sum(L[r][k]*L[c][k]for k in range(c));L[r][c]=sqrt(sv)if r==c else sv/L[c][c]
  return L
 def ilt(L):
  I=[[0.0]*6 for _ in range(6)]
  for i in range(6):
   I[i][i]=1.0/L[i][i]
   for j in range(i):I[i][j]=-sum(L[i][k]*I[k][j]for k in range(j,i))/L[i][i]
  return I
 LAi=ilt(ch6(GA))
 rng=random.Random(seed)
 n_raw=n_in=n_prim=n_np=0;lnfs=[];lngs=[]
 PHI3=((1,0,2,2,0,1),(1,2,0,2,1,2),(1,2,2,1,2,2));n_bad3=n_bad3_prim=0
 def bad3(a,b):return any(sum(f*x for f,x in zip(ph,a))%3==0 and sum(f*x for f,x in zip(ph,b))%3==0 for ph in PHI3)
 while n_raw<n_tot:
  z=[rng.gauss(0,1)for _ in range(12)];r12=(rng.random()**(1.0/12.0))/sqrt(sum(x*x for x in z))
  a=[int(round(sum(LAi[c][r]*z[c]*r12*RA0 for c in range(r,6))))for r in range(6)]
  b=[int(round(sum(LAi[c][r]*z[6+c]*r12*RB0_GA for c in range(r,6))))for r in range(6)]
  if not any(a)or not any(b):continue
  qa_n=WA*qfi(G_sym,a);qb_n=WB*qfi(G_sym,b)
  if qa_n+qb_n>T12_0:continue
  n_raw+=1
  if 20*qa_n>17*T12_0 or 20*qb_n>17*T12_0:continue
  if not in_jvor(a,b):continue
  n_in+=1
  b3=bad3(a,b);n_bad3+=b3
  if is_prim(a,b):
   n_bad3_prim+=b3
   n_prim+=1;lnf=lng=0.0
   for j in range(3):
    aj=sum(EMB[j][k]*a[k]for k in range(6));bj=sum(EMB[j][k]*b[k]for k in range(6))
    a2,b2=aj*aj,bj*bj;c0j,c2j,c4j,c6j,c8j=fc_pl[j]
    fv=(((c8j*a2+c6j*b2)*a2+c4j*b2*b2)*a2+c2j*b2*b2*b2)*a2+c0j*b2*b2*b2*b2
    gv=b2+((U_BN+1)*w_pl[j]+2.0)*a2
    lnf+=2.0*log2(abs(fv));lng+=2.0*log2(abs(gv))
   lnfs.append(lnf+M.CONTENT_BITS);lngs.append(lng)
  else:n_np+=1
 bins=PQ["bins"];np_200=sum(b[2]for b in bins)
 m_f200=sum(b[0]*b[2]for b in bins)/np_200;m_g200=sum(b[1]*b[2]for b in bins)/np_200
 f_reg=n_in/float(n_raw);p_rt=n_prim/max(1,n_in);p_exp=PQ["n_prim_in"]/float(PQ["n_in"]);m_f=sum(lnfs)/len(lnfs);m_g=sum(lngs)/len(lngs)
 z_reg=(f_reg-PQ["frac_reg"])/sqrt(max(1e-12,PQ["frac_reg"]*(1.0-PQ["frac_reg"])/n_raw))
 z_pr=(p_rt-p_exp)/sqrt(p_exp*(1.0-p_exp)/max(1,n_in))
 se_f=sqrt(sum((x-m_f)**2 for x in lnfs)/max(1,len(lnfs)-1)/len(lnfs))
 se_g=sqrt(sum((x-m_g)**2 for x in lngs)/max(1,len(lngs)-1)/len(lngs))
 z_f=(m_f-m_f200)/se_f;z_g=(m_g-m_g200)/se_g
 p3=217.0/729.0;f3=n_bad3/float(max(1,n_in));z3=(f3-p3)/sqrt(p3*(1.0-p3)/max(1,n_in))
 ok=(abs(z_reg)<=2.5 and abs(z_pr)<=2.5 and abs(z_f)<=2.5 and abs(z_g)<=2.5 and n_prim+n_np==n_in and n_bad3_prim==0 and abs(z3)<=2.5)
 return dict(ok=ok,mod3_bad_cells=n_bad3,mod3_bad_frac=round(f3,4),mod3_bad_frac_exact=round(p3,4),z_mod3_bad=round(z3,2),mod3_bad_but_primitive=n_bad3_prim,zero_float_fd_and_hnf=True,region="E12_JointVor",n_samples=n_raw,n_active=n_in,n_prim=n_prim,n_nonprim=n_np,frac_reg=round(f_reg,4),frac_reg_200k=PQ["frac_reg"],z_frac_reg=round(z_reg,2),
 prim_rate=round(p_rt,4),prim_rate_200k=round(p_exp,4),z_prim_rate=round(z_pr,2),mean_lnf=round(m_f,2),mean_lnf_200k=round(m_f200,2),z_lnf=round(z_f,2),
 mean_lng=round(m_g,2),mean_lng_200k=round(m_g200,2),z_lng=round(z_g,2))
def verify_live_genuine_c6_relation_filter(seed=20261120):
 import c6_adaptive as _C6A
 return _C6A.verify_live_genuine_c6_relation_filter(seed=seed)
verify_live_bareiss_12d_lll=DESC.verify_live_bareiss_12d_lll
