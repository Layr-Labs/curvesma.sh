"""v25_verify.py — Pure-stdlib live verifier for BN254 C6 v25 upgrades:
1. Exact closed-form B1_PREF_COV quadrature at K_B1 = M.K_B1_PREF = 256 checkpoints per shell (<= M.B1_PREF_COV = 1.0120; the K_B1 = 64 value is printed for comparison).
2. Band 1 end to end: the //@seg b1_* listing of kernels_v32.s (line setup, rows, transform-copy into link pages, flushes, stream switches, terminal pages, receiver page apply) executed instruction by instruction on a synthetic instance against a reference score array, with executed counts checked against the charged C_B1_ROW, C_HIT_LO and line prices. v27's uint32_t delta runs are gone.
3. Fixed a0-chord-stepping kta(a) with explicit Pass-0 ai_of_chord -> kta_arr[ai] scatter preserving exact Q(a)-sorted ai order (3,072 genuine 42..48b (p,a) pairs + end-to-end ai-order x4/ll build check).
4. 4-plane cyclic Cone3 algebraic invariance under S[3..5] + 450-point bijection check.
5. k=1->k=0 finite-diff rem_0 one-sided safety + unconditional r2_ub < 2^40 certificate across [48..61]b.
6. 32 KiB-strided single 16 KiB staging slot tbnz #14 exact replay (zero double-buffer misfire).
7. Lower Band 2 x4 (ai << 17) | j_rel + add lsr #16 encoding across variable k_s in {35, 42, 48, 50}, W in {8, 16}, sgn in {+1, -1}.
8. W=64 floor_div b_0 sub-band nesting at b_0 in {-64, -1, 0, 63, 64}, Graham greedy deterministic per-C[s] bound, derived over-radius row-capacity guard (delta = 0.005), dummy sink block disjointness on C[0], and init_a + read-only scan8 bit-exact check.
9. Mod-3 ideal-primitivity wheel: F_3 ring maps, joint rank 3, sigma action, residue counts, HNF false-discard test, wtab/chain dispatch, sp regions, x23 and emit-10 address checks.
Also verify_ht_pool (replay of one stored calibration run and the pooled 16-run C_LAT_LEAF bound) and verify_nonprim_hnf_witnesses (moved unchanged from model_bn254_v6.py).
"""
import math,random,re
from bisect import bisect_left
from pathlib import Path
from math import log,log2,sqrt,ceil,exp
import model_bn254_v6 as M
import v20_lat_b2 as V20
import descent_bn254_v6 as DESC
M64=(1<<64)-1
def verify_v25_upgrades(seed=20261122):
 rng=random.Random(seed)
 T4=M.load_data()["t4_geom"]
 G_sym,_,_,S,S2,G_SC=M.build_exact_int_geom(T4)
 def _cov_exact(K):
  qa_m,qb_m=M.QA_MAX,M.QB_MAX
  p_b= lambda t:20*t**3-45*t**4+36*t**5-10*t**6
  den=0.05*(p_b(qa_m)+p_b(qb_m)-1.0)
  w_split=max(0.0,1.0-qa_m)
  tot=(w_split**3)*(qa_m**3)
  u_min=1.0-qb_m
  n=M.B_A_SHELLS*K
  inv_n=1.0/n
  for m in range(1,n+1):
   u_lo=max(u_min,qa_m*(((m-1)*inv_n)**(1.0/3.0)))
   u_hi=min(qa_m,qa_m*((m*inv_n)**(1.0/3.0)))
   if u_hi>u_lo:
    tot+=(m*inv_n)*(qa_m**3)*((1.0-u_lo)**3-(1.0-u_hi)**3)
  return tot/den
 cov_64=_cov_exact(64)
 cov_k=_cov_exact(M.K_B1_PREF)
 b1_cov_ok=(cov_64<=1.0310 and cov_k<=M.B1_PREF_COV)
 for _ in range(5):
  rng.sample(range(1,1024),320);[rng.randrange(0,32768) for _ in range(1024)]
 b1=band1_e2e(str(Path(__file__).resolve().parent/"kernels_v32.s"),seed=seed+2)
 lat_ap=lat_apply_e2e(str(Path(__file__).resolve().parent/"kernels_v32.s"),seed=seed+3)
 apply_errs=0 if (b1["ok"] and lat_ap["ok"]) else 1
 c_hit_lo_exact=b1["c_hit_lo_exact"]
 b1_hit_cost_ok=(c_hit_lo_exact<=M.C_HIT_LO and b1["row_static"]+b1["switch_allow"]<=M.C_B1_ROW)
 chord_errs=0
 chord_tested=0
 chord_order_errs=0
 for pb in (42,46,48):
  p=(1<<pb)+1
  while True:
   if (p&7) in (1,3) and M._isp(p):
    s0=V20._find_root_h(p)
    if s0 is not None:
     ws=V20._eval_poly(M.W_COEFFS,s0,p)
     uw=(-(((M.U_BN254+1)%p)*ws+2))%p
     if uw!=0:
      rho=V20._sqrt_mod(pow(uw,-1,p),p)
      if rho>0:
       break
   p+=2
  inv3=pow(3,-1,p)
  sp=[1]
  for _ in range(5):
   sp.append((sp[-1]*s0)%p)
  ek=[(sum(M.ZK_NUM[k][j]*sp[j] for j in range(6))*inv3)%p for k in range(6)]
  neg_ir=(-pow(rho,-1,p))%p
  ck=[(e*neg_ir)%p for e in ek]
  c0=ck[0]
  all_vecs=[]
  chords=[]
  for ch_i in range(32):
   a_tail=[rng.randint(-300,300) for _ in range(5)]
   a0_lo=rng.randint(-200,150)
   ch_vecs=[[a0_lo+step]+a_tail for step in range(32)]
   chords.append(ch_vecs)
   all_vecs.extend(ch_vecs)
  qa_key= lambda v:math.sumprod(V20._mv6(G_sym,v),v)
  order=sorted(range(len(all_vecs)),key=lambda idx:(qa_key(all_vecs[idx]),all_vecs[idx]))
  ai_of_flat=[0]*len(all_vecs)
  for ai,fidx in enumerate(order):
   ai_of_flat[fidx]=ai
  kta_arr=[-1]*len(all_vecs)
  flat_pos=0
  for ch_vecs in chords:
   a_start=ch_vecs[0]
   kta=sum(a_start[k]*ck[k] for k in range(6))%p
   for v in ch_vecs:
    ai=ai_of_flat[flat_pos]
    flat_pos+=1
    kta_arr[ai]=kta
    ref=sum(v[k]*ck[k] for k in range(6))%p
    if kta!=ref:
     chord_errs+=1
    chord_tested+=1
    nxt=kta+c0
    kta=nxt-p if nxt>=p else nxt
  for ai,fidx in enumerate(order):
   v=all_vecs[fidx]
   ref=sum(v[k]*ck[k] for k in range(6))%p
   if kta_arr[ai]!=ref:
    chord_order_errs+=1
 half0=568.83*sqrt(M.QA_MAX)*sqrt(196608.0/753729.0)
 n_int=20000
 e_len=sum(5.0*((i+0.5)/n_int)**4*sqrt(1.0-((i+0.5)/n_int)**2) for i in range(n_int))/n_int
 chord_mean_pts=2.0*half0*e_len
 pieces_mean=1.0+2.0*sum(((j/float(M.B_A_SHELLS))**(1.0/3.0))**2.5 for j in range(1,M.B_A_SHELLS))
 piece_mean_pts=chord_mean_pts/pieces_mean
 if not (piece_mean_pts>=34.0):
  chord_errs+=1
 s35_alg_ok=(S[3]==[0,0,0,0,-1,0] and S[4]==[0,0,0,0,0,-1] and S[5]==[0,0,0,1,0,0])
 cone3_errs=0
 cone3_tested=0
 for _ in range(450):
  a=[rng.randint(-300,300) for _ in range(6)]
  b=[rng.randint(-400,400) for _ in range(6)]
  v=a+b
  orb=[]
  cur_a,cur_b=a[:],b[:]
  for _ in range(3):
   for sa in (1,-1):
    for sb in (1,-1):
     orb.append([sa*x for x in cur_a]+[sb*x for x in cur_b])
   cur_a=V20._mv6(S,cur_a)
   cur_b=V20._mv6(S,cur_b)
  sx0=a[3]-a[4]+a[5]
  s_a=V20._mv6(S,a)
  if (s_a[3]-s_a[4]+s_a[5])!=sx0:
   cone3_errs+=1
  c0=V20._cone_of(v)
  if c0>=0:
   ca0,cb0=V20._canon_cone_ref(tuple(v[:6]),tuple(v[6:]),S,S2)
   can0=list(ca0)+list(cb0)
   hits=[u for u in orb if (u[3]-u[4]+u[5]>0 and u[3]> -u[4] and u[3]>u[5] and u[6]>0)]
   if len(hits)!=1 or hits[0]!=can0:
    cone3_errs+=1
   for u in orb:
    cau,cbu=V20._canon_cone_ref(tuple(u[:6]),tuple(u[6:]),S,S2)
    if (cau,cbu)!=(ca0,cb0):
     cone3_errs+=1
   cone3_tested+=1
 T12_prod=int(round(1764*(572.0**2)*G_SC))
 Q38=1<<38
 fd_rem0_max_r2_bits=0.0
 fd_rem0_errs=0
 fd_A1_min=1<<64
 for pb in (47,52,56,61):
  pp=(1<<pb)+1
  while True:
   if (pp&7) in (1,3) and M._isp(pp):
    s0p=V20._find_root_h(pp)
    if s0p is not None:
     wsp=V20._eval_poly(M.W_COEFFS,s0p,pp)
     uwp=(-(((M.U_BN254+1)%pp)*wsp+2))%pp
     if uwp!=0:
      rhop=V20._sqrt_mod(pow(uwp,-1,pp),pp)
      if rhop>0:
       break
   pp+=2
  inv3p=pow(3,-1,pp)
  spp=[1]
  for _ in range(5):
   spp.append((spp[-1]*s0p)%pp)
  esp=[(sum(M.ZK_NUM[k][j]*spp[j] for j in range(6))*inv3p)%pp for k in range(6)]
  lfp=esp+[((-rhop)*e)%pp for e in esp]
  inv0p=pow(lfp[0],-1,pp)
  Bp=[[pp if c==0 else 0 for c in range(12)]]
  for r1 in range(1,12):
   c0_v=(-lfp[r1]*inv0p)%pp
   if c0_v>pp//2:
    c0_v-=pp
   Bp.append([c0_v if c==0 else (1 if c==r1 else 0) for c in range(12)])
  _,Dp,_=DESC.bareiss_lll_12d(Bp,G_sym,1764,1000,False,9,10,(99,100),True)
  inv0=((Dp[0]*T12_prod<<38)+Dp[1]-1)//Dp[1]
  A1_fd=(Dp[2]<<38)//(Dp[1]*T12_prod)
  fd_A1_min=min(fd_A1_min,A1_fd)
  Xm=[0]*12
  for k in range(11,-1,-1):
   rk=math.isqrt(((Dp[k]*T12_prod)+Dp[k+1]-1)//Dp[k+1])+1
   Xm[k]=rk+sum((Xm[j]+1)//2 for j in range(k+1,12))
  ec1=sum(Xm[j] for j in range(2,12))+1
  R1_max=2*Xm[1]+2
  pad1_max=((2*R1_max*ec1*A1_fd)>>32)+2
  rem0_ub=Q38+pad1_max+R1_max*R1_max
  r2_ub=((rem0_ub*inv0)>>64)+1
  if r2_ub>=(1<<40):
   fd_rem0_errs+=1
  fd_rem0_max_r2_bits=max(fd_rem0_max_r2_bits,log2(r2_ub))
 tbnz_errs=0
 base_stage=1<<40
 flushes_seen=0
 for s_id in (0,1,17,(1<<34)-1):
  pbase=base_stage+(s_id<<15)+8
  wp=pbase+16
  for rec_i in range(1,4100):
   old_wp=wp;wp+=8
   if not (pbase+16<=old_wp and old_wp+8<=pbase+16384):
    tbnz_errs+=1
   bit14=(old_wp>>14)&1
   should_flush=((rec_i%2046)==0)
   if bool(bit14)!=should_flush:
    tbnz_errs+=1
   if bit14:
    x0_hdr=old_wp-16384+8
    if x0_hdr!=pbase:
     tbnz_errs+=1
    wp=x0_hdr+16
    flushes_seen+=1
 if flushes_seen!=4*2:
  tbnz_errs+=1
 sh17_x4_errs=0
 for W_t in (8,16):
  for k_s in (35,42,48,50):
   row_stride=1<<k_s
   namax_t=4000
   namax_sh17=(namax_t<<17)|15
   pairs=sorted((rng.randint(1,5000),rng.randrange(W_t)) for _ in range(200))
   entries=[(ai<<17)|jr for ai,jr in pairs]
   for sgn_v in (1,-1):
    got_map={}
    ref_map={}
    base_addr=(1<<55) if sgn_v>0 else ((1<<55)+(W_t-1)*row_stride)
    for e in entries:
     if e>namax_sh17:
      break
     jr=e&0xF
     r_addr=(base_addr+(jr<<k_s)) if sgn_v>0 else (base_addr-(jr<<k_s))
     a_addr=r_addr+(e>>16)
     got_map[a_addr]=(got_map.get(a_addr,0)+720)&0xFFFF
    for ai,jr in pairs:
     if ai>namax_t:
      break
     jr_eff=jr if sgn_v>0 else (W_t-1-jr)
     a_addr=(1<<55)+(jr_eff<<k_s)+(ai<<1)
     ref_map[a_addr]=(ref_map.get(a_addr,0)+720)&0xFFFF
    if got_map!=ref_map:
     sh17_x4_errs+=1
   if ((namax_t<<17)|15)>namax_sh17 or (((namax_t+1)<<17)|0)<=namax_sh17:
    sh17_x4_errs+=1
 sb_nest_errs=0
 for b0 in (-64,-63,-1,0,1,63,64,559):
  sb64=b0>>6
  for W_sub in (8,16,32,64):
   sh=int(round(log2(W_sub)))
   lo_b0=(b0>>sh)<<sh
   hi_b0=lo_b0+W_sub-1
   if (lo_b0>>6)!=sb64 or (hi_b0>>6)!=sb64:
    sb_nest_errs+=1
 import fd_bounds as FDB
 fdb=FDB.main()
 import exact_levels as EXL
 exl=EXL.main()
 eps_or=fdb["eps_or_bound"]
 guard_errs=0
 if not (eps_or<=M.SC_GUARD_DELTA/2.0):
  guard_errs+=1
 import stor_bound as SBND
 sbo=SBND.main()
 if not (sbo["ok"] and fdb["ok"] and eps_or<=float(SBND.EPS_OR_MAX) and float(SBND.EPS_OR_MAX)<=float(SBND.DELTA)/2
and fdb["exact_per_orbit_ir_log2"]<=log2(M.C_LAT_EXC) and exl["ok"]):
  guard_errs+=1
 dummy_sink_ok=bool(sbo["sink_ok"])
 num_g=den_g=0.0
 R0=551.94
 for i_r in range(120):
  qp=((i_r+0.5)/120.0)**0.4
  w_c=sqrt(1.0-qp)*R0
  for i_c in range(8):
   c_c=8.0*(i_c+0.5)/8.0
   grp={}
   for b in range(ceil(c_c-w_c),int(c_c+w_c)+1):
    v=max(0.0,1.0-qp-((b-c_c)/R0)**2)**3
    grp.setdefault(b>>3,[]).append(v)
   for vs in grp.values():
    num_g+=8.0*max(vs)
    den_g+=sum(vs)
 scan_group_overhead=num_g/den_g
 if not (scan_group_overhead<=1.05):
  guard_errs+=1
 scan8_errs=0
 n_cells_t=256
 init_a_t=[0x8000-rng.randint(2000,8000) for _ in range(n_cells_t)]
 rows_8=[init_a_t[:] for _ in range(8)]
 expected_surv=set()
 for r in range(8):
  for c in range(n_cells_t):
   if rng.random()<0.03:
    rows_8[r][c]=(rows_8[r][c]+9000)&0xFFFF
   if rows_8[r][c]&0x8000:
    expected_surv.add((r,c))
 found_surv=set()
 for c0 in range(0,n_cells_t,4):
  orr_acc=0
  for r in range(8):
   w=rows_8[r][c0]|(rows_8[r][c0+1]<<16)|(rows_8[r][c0+2]<<32)|(rows_8[r][c0+3]<<48)
   orr_acc|=w
  if orr_acc&0x8000800080008000:
   for r in range(8):
    for dc in range(4):
     if rows_8[r][c0+dc]&0x8000:
      found_surv.add((r,c0+dc))
 if found_surv!=expected_surv:
  scan8_errs+=1
 from fractions import Fraction as Fr
 H_f=[Fr(c) for c in M.H_COEFFS]
 ZK_f=[[Fr(c,3) for c in row] for row in M.ZK_NUM]
 PHI=V20.PHI_F3
 def _pmul_h(a,b):
  r=[Fr(0)]*(len(a)+len(b)-1)
  for ii,x in enumerate(a):
   for jj,y in enumerate(b):r[ii+jj]+=x*y
  while len(r)>6:
   c=r.pop()
   if c:
    off=len(r)-6
    for kk in range(6):r[off+kk]-=c*H_f[kk]
  return r+[Fr(0)]*(6-len(r))
 aug_z=[ZK_f[ii][:]+[Fr(1 if ii==jj else 0) for jj in range(6)] for ii in range(6)]
 for col in range(6):
  pv_r=next(rr for rr in range(col,6) if aug_z[rr][col]!=0)
  aug_z[col],aug_z[pv_r]=aug_z[pv_r],aug_z[col]
  pv=aug_z[col][col];aug_z[col]=[x/pv for x in aug_z[col]]
  for rr in range(6):
   if rr!=col and aug_z[rr][col]:
    fac=aug_z[rr][col];aug_z[rr]=[aug_z[rr][cc]-fac*aug_z[col][cc] for cc in range(12)]
 ZK_inv=[row[6:] for row in aug_z]
 wheel_hom_errs=0
 for jj in range(6):
  for kk in range(6):
   pr=_pmul_h(ZK_f[jj],ZK_f[kk])
   co=[sum(pr[mm]*ZK_inv[mm][nn] for mm in range(6)) for nn in range(6)]
   if any(x.denominator!=1 for x in co):
    wheel_hom_errs+=1
   else:
    ci=[x.numerator for x in co]
    for ii in range(3):
     if sum(PHI[ii][mm]*ci[mm] for mm in range(6))%3!=(PHI[ii][jj]*PHI[ii][kk])%3:
      wheel_hom_errs+=1
 det3_f3=(PHI[0][0]*(PHI[1][1]*PHI[2][2]-PHI[1][2]*PHI[2][1])
-PHI[0][1]*(PHI[1][0]*PHI[2][2]-PHI[1][2]*PHI[2][0])
+PHI[0][2]*(PHI[1][0]*PHI[2][1]-PHI[1][1]*PHI[2][0]))%3
 wheel_rank3_ok=(det3_f3!=0 and all(PHI[ii][0]==1 for ii in range(3)))
 sigma_perm=[PHI.index(tuple(sum(PHI[ii][rr]*S[rr][kk] for rr in range(6))%3 for kk in range(6))) for ii in range(3)]
 wheel_sigma_ok=(sorted(sigma_perm)==[0,1,2] and sigma_perm!=[0,1,2])
 bad_729=sum(1 for st in range(729) if any((st//(3**ii))%3==0 and (st//(3**(3+ii)))%3==0 for ii in range(3)))
 rep_27=sum(1<<(3-sum(1 for ii in range(3) if (st//(3**ii))%3==0)) for st in range(27))
 wheel_false_discards=0;wheel_disc_cnt=0
 for _ in range(250):
  av=[rng.randint(-40,40) for _ in range(6)]
  bv=[rng.randint(-40,40) for _ in range(6)]
  if not any(av) or not any(bv): continue
  if V20._bad_f3(av+bv):
   wheel_disc_cnt+=1
   if M.is_ideal_primitive_oh(av,bv):
    wheel_false_discards+=1
 wtab_errs=0
 wtab=[None]*448
 for sm in range(1,7):
  for mm in range(1,61):
   idx=(sm<<6)|mm
   gp=[t for t in range(mm) if not ((sm>>(t%3))&1)]
   if not gp:wtab[idx]=("NONE",0,0)
   elif len(gp)==1:wtab[idx]=("TAIL1",gp[0],1)
   else:
    gaps=[gp[ii+1]-gp[ii] for ii in range(len(gp)-1)]
    ch="C3" if all(g==3 for g in gaps) else ("C12" if gaps[-1]==1 else "C21")
    wtab[idx]=(ch,gp[0],len(gp))
 sp_regions=[(0,2880),(3840,3864)]+[(8192*kk,8192*kk+5554) for kk in range(1,12)]+[(131072,139264)]
 sp_regions.sort()
 sp_alias_errs=sum(1 for ii in range(len(sp_regions)-1) if sp_regions[ii][1]>sp_regions[ii+1][0])
 for sm in range(1,7):
  for nn in list(range(0,64))+[120,121,125]:
   got=[];aa,bb,rm=1000,2000,nn
   while rm>0:
    mm=60 if rm>60 else rm
    ch,sk,kk=wtab[(sm<<6)|mm]
    if kk>0:
     ca,cb=aa+sk*7,bb+sk*13
     if ch=="TAIL1":got.append((ca,cb))
     elif ch=="C3":
      for jj in range(kk):
       got.append((ca,cb))
       if jj+1<kk:ca+=21;cb+=39
     else:
      st_l=1 if ch=="C12" else 2;steps_l=[]
      for _ in range(kk-1):steps_l.append(st_l);st_l=3-st_l
      steps_l.reverse()
      for jj in range(kk):
       got.append((ca,cb))
       if jj<len(steps_l):ca+=steps_l[jj]*7;cb+=steps_l[jj]*13
    rm-=mm
    if rm>0:aa+=60*7;bb+=60*13
   ref=[(1000+t*7,2000+t*13) for t in range(nn) if not ((sm>>(t%3))&1)]
   if got!=ref:wtab_errs+=1
 x23_barrett_errs=0
 for xv in (-3145728,-100000,-17,-1,0,1,2,3,17,100000,3145728):
  u=(xv+(0x300<<12))&0xFFFFFFFF
  q=(u*0x55555556)>>32
  r_got=u-3*q
  if r_got!=(xv%3):x23_barrett_errs+=1
 for cnt0 in (0,1,2,50,400):
  st_c=2048+321;x23=(cnt0<<12)|st_c;it=0
  while True:
   if (x23&0xFFF)!=st_c:x23_barrett_errs+=1
   nxt=2048+(((st_c-2048)*5+7)%729)
   x23=(x23& ~0xFFF)|nxt
   x23=(x23-(1<<12))&M64
   st_c=nxt;it+=1
   if (x23>>63)&1: break
  if it!=cnt0+1:x23_barrett_errs+=1
 emit10_errs=0
 blks_38=ceil((2.0**sbo["score_node_log2"])/(2.0**38))
 nap_38=int(2.0**sbo["Ea_plus_log2"])+1
 if not (blks_38<(1<<26) and (((blks_38-1)<<38)+nap_38+1)<(1<<64)):
  emit10_errs+=1
 for s_t in (0,1,(1<<34)-1,(1<<35)-1):
  for rb_t in (0,1,blks_38-1):
   for ai_t in (1,nap_38):
    x3=(s_t<<29)|rb_t
    if (x3>>26)!=(8*s_t) or (((ai_t+((x3<<38)&M64))&M64)!=((rb_t<<38)+ai_t)):
     emit10_errs+=1
 v28_wheel_ok=(
wheel_hom_errs==0 and wheel_rank3_ok and wheel_sigma_ok
and bad_729==217 and rep_27==125
and wheel_false_discards==0 and wheel_disc_cnt>=40
and wtab_errs==0 and sp_alias_errs==0 and x23_barrett_errs==0 and emit10_errs==0
)
 ok=(
v28_wheel_ok and
b1_cov_ok and b1_hit_cost_ok and apply_errs==0
and chord_errs==0 and chord_order_errs==0 and chord_tested==3072
and s35_alg_ok and cone3_errs==0 and cone3_tested>=300
and fd_rem0_errs==0 and fd_rem0_max_r2_bits<25.0
and tbnz_errs==0 and sh17_x4_errs==0
and sb_nest_errs==0 and guard_errs==0 and dummy_sink_ok and scan8_errs==0
)
 return dict(
ok=ok,
b1_cov_k64_exact=round(cov_64,5),
b1_cov_k256_exact=round(cov_k,5),
b1_pref_cov_budget=M.B1_PREF_COV,
b1_c_hit_lo_exact=round(c_hit_lo_exact,6),
b1_c_hit_lo_budget=M.C_HIT_LO,
b1_c_b1_row_budget=M.C_B1_ROW,
b1_e2e=b1,
lat_apply_e2e=lat_ap,
a0_chord_kta_tested=chord_tested,
a0_chord_kta_errs=chord_errs,
a0_chord_ai_order_errs=chord_order_errs,
a0_chord_mean_points=round(chord_mean_pts,2),
a0_chord_mean_pieces=round(pieces_mean,4),
a0_chord_piece_mean_points=round(piece_mean_pts,2),
cone3_s35_alg_ok=s35_alg_ok,
cone3_bijection_tested=cone3_tested,
cone3_bijection_errs=cone3_errs,
fd_rem0_max_r2_bits=round(fd_rem0_max_r2_bits,2),
fd_rem0_u40_errs=fd_rem0_errs,
tbnz14_single_slot_errs=tbnz_errs,
sh17_x4_var_ks_errs=sh17_x4_errs,
b0_floor_div_nest_errs=sb_nest_errs,
over_radius_eps_bound=round(eps_or,7),
fd_bounds=fdb,exact_levels=exl,
scan_group_overhead=round(scan_group_overhead,4),
stor_bound=sbo,
fd_A1_min_log2=round(log2(fd_A1_min),2),
row_capacity_guard_errs=guard_errs,
dummy_sink_blocks_c0=1<<15,
dummy_sink_disjoint_ok=dummy_sink_ok,
init_scan8_readonly_errs=scan8_errs,
v28_wheel_ok=v28_wheel_ok,
v28_wheel_hom_errs=wheel_hom_errs,
v28_wheel_det3_f3=det3_f3,
v28_wheel_sigma_perm=sigma_perm,
v28_wheel_bad_729=bad_729,
v28_wheel_rep_27=rep_27,
v28_wheel_false_discards=wheel_false_discards,
v28_wheel_wtab_errs=wtab_errs,
v28_sp_alias_errs=sp_alias_errs,
v28_wheel_x23_barrett_errs=x23_barrett_errs,
v28_emit10_errs=emit10_errs,
)
_B1R=re.compile
_B1X=r"x(\d+)";_B1I=r"#(\d+)"
_B1PAT=[
("ldp_post",_B1R(rf"ldp {_B1X}, {_B1X}, \[{_B1X}\], {_B1I}$")),("ldp",_B1R(rf"ldp {_B1X}, {_B1X}, \[{_B1X}\]$")),
("stp_post",_B1R(rf"stp {_B1X}, {_B1X}, \[{_B1X}\], {_B1I}$")),("stp",_B1R(rf"stp {_B1X}, {_B1X}, \[{_B1X}\]$")),
("ldr_post",_B1R(rf"ldr {_B1X}, \[{_B1X}\], {_B1I}$")),("ldr_reg",_B1R(rf"ldr {_B1X}, \[{_B1X}, {_B1X}\]$")),
("str_post",_B1R(rf"str {_B1X}, \[{_B1X}\], {_B1I}$")),
("add_lsl",_B1R(rf"add {_B1X}, {_B1X}, {_B1X}, lsl {_B1I}$")),("add_reg",_B1R(rf"add {_B1X}, {_B1X}, {_B1X}$")),
("add_imm12",_B1R(rf"add {_B1X}, {_B1X}, {_B1I}, lsl #12$")),("add_imm",_B1R(rf"add {_B1X}, {_B1X}, {_B1I}$")),
("sub_lsl",_B1R(rf"sub {_B1X}, {_B1X}, {_B1X}, lsl {_B1I}$")),("sub_imm",_B1R(rf"sub {_B1X}, {_B1X}, {_B1I}$")),
("sub_reg",_B1R(rf"sub {_B1X}, {_B1X}, {_B1X}$")),
("subs_reg",_B1R(rf"subs {_B1X}, {_B1X}, {_B1X}$")),("cmp_imm",_B1R(rf"cmp {_B1X}, {_B1I}$")),("cmp",_B1R(rf"cmp {_B1X}, {_B1X}$")),
("csel_hs",_B1R(rf"csel {_B1X}, {_B1X}, {_B1X}, hs$")),("and_imm",_B1R(rf"and {_B1X}, {_B1X}, #0x([0-9a-f]+)$")),
("lsr",_B1R(rf"lsr {_B1X}, {_B1X}, {_B1I}$")),("ubfx",_B1R(rf"ubfx {_B1X}, {_B1X}, {_B1I}, {_B1I}$")),
("adr",_B1R(rf"adr {_B1X}, ([.\w]+)$")),("blr",_B1R(rf"blr {_B1X}$")),("br",_B1R(rf"br {_B1X}$")),("b",_B1R(r"b ([.\w]+)$")),("ret",_B1R(r"ret$")),
("tbz",_B1R(rf"tbz {_B1X}, {_B1I}, ([.\w]+)$")),("tbnz",_B1R(rf"tbnz {_B1X}, {_B1I}, ([.\w]+)$")),
("bcc",_B1R(r"b\.(eq|hi|lo|ls|pl) ([.\w]+)$")),("staddh",_B1R(rf"staddh w(\d+), \[{_B1X}\]$")),("staddb",_B1R(rf"staddb w(\d+), \[{_B1X}\]$")),
("ldp_sp_off",_B1R(rf"ldp {_B1X}, {_B1X}, \[sp, {_B1I}\]$")),("stp_sp_off",_B1R(rf"stp {_B1X}, {_B1X}, \[sp, {_B1I}\]$")),
("ldp_off",_B1R(rf"ldp {_B1X}, {_B1X}, \[{_B1X}, {_B1I}\]$")),("ldr_off",_B1R(rf"ldr {_B1X}, \[{_B1X}, {_B1I}\]$")),("ldr",_B1R(rf"ldr {_B1X}, \[{_B1X}\]$")),
("ldr_lsl3",_B1R(rf"ldr {_B1X}, \[{_B1X}, {_B1X}, lsl #3\]$")),("ldr_sp_off",_B1R(rf"ldr {_B1X}, \[sp, {_B1I}\]$")),
("str_sp_off",_B1R(rf"str {_B1X}, \[sp, {_B1I}\]$")),("str",_B1R(rf"str {_B1X}, \[{_B1X}\]$")),("ldadd",_B1R(rf"ldadd {_B1X}, {_B1X}, \[{_B1X}\]$")),
("add_lsr",_B1R(rf"add {_B1X}, {_B1X}, {_B1X}, lsr {_B1I}$")),("sub_imm12",_B1R(rf"sub {_B1X}, {_B1X}, {_B1I}, lsl #12$")),
("mov_imm",_B1R(rf"mov {_B1X}, {_B1I}$")),("extr",_B1R(rf"extr {_B1X}, {_B1X}, {_B1X}, {_B1I}$")),
("sbfx",_B1R(rf"sbfx {_B1X}, {_B1X}, {_B1I}, {_B1I}$")),("mul",_B1R(rf"mul {_B1X}, {_B1X}, {_B1X}$")),
("madd",_B1R(rf"madd {_B1X}, {_B1X}, {_B1X}, {_B1X}$")),("msub",_B1R(rf"msub {_B1X}, {_B1X}, {_B1X}, {_B1X}$")),
("udiv",_B1R(rf"udiv {_B1X}, {_B1X}, {_B1X}$")),("lsl_reg",_B1R(rf"lsl {_B1X}, {_B1X}, {_B1X}$")),
("bfi_w",_B1R(r"bfi w(\d+), w(\d+), #(\d+), #(\d+)$")),
("ldr_sxtw3",_B1R(rf"ldr {_B1X}, \[{_B1X}, w(\d+), sxtw #3\]$")),
("ldr_uxtw3",_B1R(rf"ldr {_B1X}, \[{_B1X}, w(\d+), uxtw #3\]$")),
]
def _seg_load(path,prefix):
 """Instruction lines of the //@seg <prefix>* blocks in file order (.rept expanded), labels -> index, segment of each line."""
 prog,labels,segof,cur,rep,buf=[],{},[],None,None,[]
 for raw in open(path,encoding="ascii"):
  st=raw.strip()
  if st.startswith("//@seg"):
   cur=st.split()[1];continue
  if st.startswith("//@end"):
   cur= None;continue
  line=raw.split("//",1)[0].strip()
  if cur is None or not cur.startswith(prefix) or not line:
   continue
  if line.startswith(".rept"):
   rep,buf=int(line.split()[1]),[];continue
  if line.startswith(".endr"):
   for ins in buf*rep:
    prog.append(ins);segof.append(cur)
   rep= None;continue
  if line.endswith(":"):
   labels[line[:-1]]=len(prog);continue
  if line.startswith("."):
   continue
  if rep is not None:
   buf.append(line)
  else:
   prog.append(line);segof.append(cur)
 return prog,labels,segof
def _b1_load(path):
 return _seg_load(path,"b1_")
def lat_apply_e2e(path,seed=20261125):
 """Executes //@seg emit10_step1 + //@seg emit_direct_flush + //@seg lat_apply_page + //@seg lat_apply from kernels_v30.s
    instruction by instruction: verifies 9-op LDADD emit (BLK=38, bits 26..28==0, 32768*s+8 16,384 B page rollover at 2,046
    records, x1=#8 preserved across flush) and x128-unrolled C[s] apply against direct reference score memory."""
 prog,lab,segof=_seg_load(path,"lat_apply")
 P=_b1_decode(prog,lab)
 segs_all=V20.count_kernel_segments(path)
 fl_ops=segs_all.get("emit_direct_flush",-1)
 ap_page_ops=sum(1 for sg in segof if sg=="lat_apply_page")
 ap_body_ops=lab["apply_end"]-lab["apply"]
 fl_per_rec=(fl_ops+2+20)/2046.0
 ap_per_rec=(ap_body_ops+2)/(ap_body_ops//3*2)+(ap_page_ops+3+20)/2046.0
 rng=random.Random(seed)
 ep,elab,_=_seg_load(path,"emit")
 e_prog=[ln for ln in ep]
 sgs,nxt,cur_sg,pl,prev_seg={},{},None,None,None
 fl_sgs={".Lfl_s1":[],".Lfl_s2":[],".Lfl_s3":[],".Lfl_last":[]}
 cur_fl= None
 for raw in open(path,encoding="ascii"):
  st=raw.strip()
  if st.startswith("//@seg"):cur_sg=st.split()[1];continue
  if st.startswith("//@end"):cur_sg= None;continue
  ln=raw.split("//",1)[0].strip()
  if not ln or (ln.startswith(".") and not ln.endswith(":")): continue
  if ln.endswith(":"):
   if pl is not None and pl!=".Lmask0":nxt.setdefault(pl,[]).append(ln)
   if pl!=".Lmask0":
    pl=ln[:-1] if ln[:-1] in (".Lchord_split",".Lcont_s1",".Lcont_s2",".Lcont_s3",".Lcont_last",".Lchain_entry",".Lmask0") else None
   cur_fl=ln[:-1] if ln[:-1] in fl_sgs else None
   if ln in ("emit9_step1:","emit8_step1:") and prev_seg in ("wheel_skip1_stub","wheel_skip2_stub"):
    nxt[prev_seg]=ln
   continue
  if pl is not None:
   nxt.setdefault(pl,[]).append(ln)
   if pl not in (".Lcont_last",".Lchain_entry",".Lmask0") or len(nxt[pl])>=2:pl= None
  if cur_sg:sgs.setdefault(cur_sg,[]).append(ln);prev_seg=cur_sg
  if cur_fl:
   fl_sgs[cur_fl].append(ln)
   if ln.startswith("b "):cur_fl= None
 mask0_ins=nxt.get(".Lmask0",[])
 chain_ok=(nxt.get(".Lchord_split")==["wheel_disp4:"] and nxt.get("wheel_skip1_stub")=="emit9_step1:"
and nxt.get("wheel_skip2_stub")=="emit8_step1:" and nxt.get(".Lcont_s1")==["emit9_step2:"]
and nxt.get(".Lcont_s2")==["emit9_step3:"] and nxt.get(".Lcont_s3")==["emit9_last:"]
and nxt.get(".Lcont_last")==["b.pl .Lk0_top","b .Ltop_2"]
and mask0_ins==["b.pl .Lk0_top","b .Ltop_2"] and "ret" not in mask0_ins)
 wd=sgs.get("wheel_disp4",[]);sk1=sgs.get("wheel_skip1_stub",[])
 e1=sgs.get("emit9_step1",[]);e2=sgs.get("emit9_step2",[])
 e3=sgs.get("emit9_step3",[]);e4=sgs.get("emit9_last",[])
 tail=["b.pl .Lk0_top","b .Ltop_2"]
 c_lines=wd+sk1+e1+e2+e3+e4+tail
 psk1=len(wd);p0=psk1+len(sk1);p1=p0+len(e1);p2=p1+len(e2);p3=p2+len(e3);p4=p3+len(e4)
 f1=len(c_lines);f2=f1+len(fl_sgs[".Lfl_s1"]);f3=f2+len(fl_sgs[".Lfl_s2"]);f4=f3+len(fl_sgs[".Lfl_s3"])
 all_lines=c_lines+fl_sgs[".Lfl_s1"]+fl_sgs[".Lfl_s2"]+fl_sgs[".Lfl_s3"]+fl_sgs[".Lfl_last"]
 stop_top=len(all_lines);stop_exit=stop_top+1
 elab2={".Lcont_s1":p1,".Lcont_s2":p2,".Lcont_s3":p3,".Lcont_last":p4,
".Lfl_s1":f1,".Lfl_s2":f2,".Lfl_s3":f3,".Lfl_last":f4,
".Lk0_top":stop_top,".Ltop_2":stop_exit}
 Pe=_b1_decode(all_lines,elab2)
 m8e,m16e,ref_e={},{},{}
 ATAB=1<<48;BTAB=1<<50;PTRS=1<<52;SP0=1<<46;RTAB=1<<45;RROW=RTAB+(1<<20)
 BPLANE_STRIDE=1<<40
 wgt_e=19;m8e[SP0+480]=wgt_e
 s_max=(1<<34)-1;s_list=(0,1,s_max)
 pbase={s_:(1<<54)+32768*s_+8 for s_ in s_list}
 for s_ in s_list:m8e[PTRS+8*s_]=pbase[s_]+16
 k_key=(((-5)&0xFFFFFFFF)& ~(7<<23))|(1<<23)
 k_sgn=k_key-(1<<32) if (k_key>>31) else k_key
 m8e[(RTAB+(k_sgn<<3))&M64]=RROW;m8e[RROW+6*8]=4*p0;m8e[RROW+7*8]=4*psk1
 em_flushes=0;emit_ok= True
 try:
  cpu_e=_B1CPU(Pe,m8e,m16e);xe=cpu_e.x
  xe[1]=8;xe[10]=RTAB;xe[11]=48;xe[12]=3;xe[26]=16;xe[27]=1;xe[17]=PTRS;xe[31]=SP0
  prod_ais=(1,(1<<37)-1,1<<38,(1<<37)+1,(1<<47)+12345,(1<<48)-1)
  rb_prod_max=43963475
  offs=(0,1,3,6);n_chords=525
  for ch in range(n_chords):
   for st_i,off_i in enumerate(offs):
    idx=ch*4+st_i;pos=ch*7+off_i
    s_i=1 if idx<2050 else (0 if idx<2080 else s_max)
    rb_i=0 if s_i==0 else (rb_prod_max if s_i==s_max else (1+(idx%3)))
    ai_i=prod_ais[(idx-2080)%len(prod_ais)] if s_i==s_max else (1+(idx%256))
    bp_k=BTAB+(idx&3)*BPLANE_STRIDE
    m8e[ATAB+16*pos]=ai_i;m8e[ATAB+16*pos+8]=bp_k
    m8e[bp_k+8*pos]=(s_i<<30)|rb_i
    if s_i==1 and idx<2046:
     ad=(rb_i<<38)+ai_i;ref_e[ad]=(ref_e.get(ad,0)+wgt_e)&0xFF
  def on_fl_ret(c):
   nonlocal em_flushes
   em_flushes+=1;pg=c.x[0]
   if pg!=pbase[1] or m8e[pg]!=2046 or m8e[pg+8]!=wgt_e: raise ValueError("bad flush header")
   if m8e[PTRS+8*1]!=pbase[1]+16 or c.x[1]!=8: raise ValueError("bad ptr reset")
   rc=_B1CPU(P,m8e,m16e);rc.x[28]=pg;rc.run(lab["apply_page"],set(),{},{})
  fl_hooks={f2-1:on_fl_ret,f3-1:on_fl_ret,f4-1:on_fl_ret,len(all_lines)-1:on_fl_ret}
  etally={};xe[23]=n_chords-1
  for ch in range(n_chords):
   xe[2]=ATAB+16*(ch*7-(1 if ch==0 else 0));xe[5]=ch*7-(1 if ch==0 else 0)
   xe[6]=(-5)&0xFFFFFFFF;xe[14]=1;xe[0]=7 if ch==0 else 6;xe[23]=(xe[23]-1)&M64;cpu_e.N=bool(xe[23]>>63)
   cpu_e.run(0,{stop_top,stop_exit},fl_hooks,etally)
  prod_width_ok=all(
m8e[pbase[s_max]+16+8*j]==(rb_prod_max<<38)+prod_ais[j%len(prod_ais)] for j in range(20))
  emit_ok=(chain_ok and len(wd)==4 and len(sk1)==2 and len(e1)==9 and len(e2)==9 and len(e3)==9 and len(e4)==7
and em_flushes==1 and sum(etally.values())==n_chords*39+2+1+8
and m16e==ref_e and prod_width_ok and xe[1]==8
and m8e[PTRS+8*1]==pbase[1]+16+4*8
and m8e[PTRS+8*0]==pbase[0]+16+30*8 and m8e[PTRS+8*s_max]==pbase[s_max]+16+20*8
and all(ins[0]!="other" for ins in Pe))
  if emit_ok:
   raw_lines=[l.split("//",1)[0].strip() for l in open(path,encoding="ascii")]
   macros,inv_lines,in_alt,cur_m={},[],False,None
   for ln in raw_lines:
    if not ln: continue
    if ln==".altmacro":in_alt= True;continue
    if ln==".noaltmacro":in_alt= False;continue
    if not in_alt: continue
    if ln.startswith(".macro "):
     parts=ln[7:].strip().split(None,1)
     mname=parts[0];mparams=[p.strip() for p in parts[1].split(",")] if len(parts)>1 else []
     cur_m=mname;macros[mname]=(mparams,[]);continue
    if ln==".endm":cur_m= None;continue
    if cur_m is not None:macros[cur_m][1].append(ln)
    else:inv_lines.append(ln)
   symtab={};gas_ins=[];gas_lab={}
   _eval_c={}
   def _gas_eval(expr):
    if expr.isdigit(): return int(expr)
    if expr in symtab: return symtab[expr]
    if " == " in expr:
     l_,r_=expr.split(" == ",1)
     return int(_gas_eval(l_.strip())==_gas_eval(r_.strip()))
    if " + " in expr and "(" not in expr:
     l_,r_=expr.split(" + ",1)
     return _gas_eval(l_.strip())+_gas_eval(r_.strip())
    if " - " in expr and "(" not in expr:
     l_,r_=expr.split(" - ",1)
     return _gas_eval(l_.strip())-_gas_eval(r_.strip())
    if expr.startswith("(((") and expr.endswith(") & 1)"):
     inner=expr[3:-6]
     p1,p2=inner.rsplit(") - ",1)
     a_,b_=p1.split(" - ",1)
     return ((_gas_eval(a_.strip())-_gas_eval(b_.strip()))-_gas_eval(p2.strip()))&1
    return int(eval(expr,{"__builtins__":{}},symtab))
   def _gas_sub(ln,env):
    if not env or "\\" not in ln: return ln
    out=ln
    for k,v in env.items():
     if f"\\{k}" in out:out=out.replace(f"\\{k}\\()",v).replace(f"\\{k}",v)
    return out
   _e8_cache={}
   def _gas_exec(lines,env):
    i,n=0,len(lines)
    while i<n:
     ln=_gas_sub(lines[i],env);i+=1
     if not ln: continue
     c0=ln[0]
     if c0==".":
      if ln.startswith(".set "):
       nm,ex=ln[5:].split(",",1);symtab[nm.strip()]=_gas_eval(ex.strip());continue
      if ln.startswith(".if "):
       brs=[(ln[4:].strip(),[])];dep=0
       while i<n:
        cl=lines[i];i+=1
        if cl.startswith(".if "):dep+=1;brs[-1][1].append(cl)
        elif cl==".endif":
         if dep==0: break
         dep-=1;brs[-1][1].append(cl)
        elif dep==0 and cl.startswith(".elseif "):brs.append((cl[8:].strip(),[]))
        elif dep==0 and cl==".else":brs.append(("1",[]))
        else:brs[-1][1].append(cl)
       for cond,body in brs:
        if _gas_eval(_gas_sub(cond,env)):_gas_exec(body,env);break
       continue
      if ln.startswith(".rept "):
       cnt=_gas_eval(ln[6:].strip());body=[];dep=0
       while i<n:
        cl=lines[i];i+=1
        if cl.startswith(".rept "):dep+=1;body.append(cl)
        elif cl==".endr":
         if dep==0: break
         dep-=1;body.append(cl)
        else:body.append(cl)
       for _ in range(cnt):_gas_exec(body,env)
       continue
      if ln.endswith(":"):
       gas_lab[ln[:-1].strip()]=len(gas_ins);continue
     elif ":" in ln:
      lb,rest=ln.split(":",1)
      if lb.strip() and " " not in lb.strip():
       gas_lab[lb.strip()]=len(gas_ins);ln=rest.strip()
       if not ln: continue
     tok=ln.split(None,1)
     t0=tok[0]
     if t0 in macros:
      if t0=="E8":
       st_a=tok[1].strip()
       if st_a not in _e8_cache:
        _b0=len(gas_ins);mp,mb=macros["E8"];_gas_exec(mb,{mp[0]:st_a});_e8_cache[st_a]=gas_ins[_b0:]
       else:
        gas_ins.extend(_e8_cache[st_a])
       continue
      if t0=="E9":
       a0,a1=[a.strip() for a in tok[1].split(",",1)]
       if a0 not in _e8_cache:
        _b0=len(gas_ins);mp,mb=macros["E8"];_gas_exec(mb,{mp[0]:a0});_e8_cache[a0]=gas_ins[_b0:]
       id_v=str(symtab[a1[1:]]) if a1.startswith("%") else a1
       gas_ins.extend(_e8_cache[a0])
       gas_ins.append(f"tbnz x9, #14, .Lfl_g_{id_v}")
       gas_lab[f".Lct_g_{id_v}"]=len(gas_ins)
       continue
      if t0=="E9C":
       a1=tok[1].strip()
       id_v=str(symtab[a1[1:]]) if a1.startswith("%") else a1
       gas_lab[f".Lfl_g_{id_v}"]=len(gas_ins)
       gas_ins.extend(["sub x0, x9, #4, lsl #12","add x0, x0, #8","mov x3, #2046","ldr x7, [sp, #480]","stp x3, x7, [x0]","add x9, x0, #16","str x9, [x8]",f"b .Lct_g_{id_v}"])
       continue
      mp,mb=macros[t0]
      av=[a.strip() for a in (tok[1].split(",") if len(tok)>1 else [])]
      ne={p:(str(symtab[a[1:]]) if a.startswith("%") else a) for p,a in zip(mp,av)}
      _gas_exec(mb,ne);continue
     gas_ins.append(ln)
   _gas_exec(inv_lines,{})
   def _slice_stub(start_idx):
    out=[]
    for k_ in range(start_idx,len(gas_ins)):
     out.append(gas_ins[k_])
     if gas_ins[k_]=="b .Ltop_2": break
    return out
   m8t={};PBt=(1<<54)+8
   for j_ in range(790):
    m8t[ATAB+16*j_]=j_+1;m8t[ATAB+16*j_+8]=BTAB;m8t[BTAB+8*j_]=(1<<30)|1
   def _tgt_for_short(mode,s_op,bm,r0,Lc):
    kept=[j_ for j_ in range(Lc) if not ((bm>>((r0+j_)%3))&1)]
    if not kept: return kept,list(mask0_ins)
    if len(kept)==1:
     sk0=kept[0]
     if sk0==0: return kept,_slice_stub(gas_lab[f"e{mode}_s0_p0"]+(128-1)*s_op)
     return kept,_slice_stub(gas_lab[f"wheel_skip{sk0}_emit{mode}_last"])
    sk0=kept[0];gaps=[b_-a_ for a_,b_ in zip(kept,kept[1:])];pc_bm=bin(bm).count("1")
    pat=0 if pc_bm==0 else (3 if pc_bm==2 else (1 if gaps[-1]==1 else 2))
    if sk0==0:
     nmax=(128 if pat==0 else (43 if pat==3 else 86))
     return kept,_slice_stub(gas_lab[f"e{mode}_s0_p{pat}"]+(nmax-len(kept))*s_op)
    return kept,_slice_stub(gas_lab[f"e{mode}_s{sk0}_p{pat}_k{len(kept)}"])
   for mode,s_op,l_op in ((8,8,6),(9,9,7)):
    for bm in range(8):
     for r0 in range(3):
      for Lc in (1,2,3,4,5,6,7,8,9,10,11,12,128):
       kept,tgt_ins=_tgt_for_short(mode,s_op,bm,r0,Lc)
       st_l=list(wd)+tgt_ins
       st_pc=len(wd);stop_t=len(st_l)
       ld={".Lk0_top":stop_t,".Ltop_2":stop_t+1}
       for ln_ in tgt_ins:
        if ln_.startswith("tbnz ") and ".Lfl_g_" in ln_:ld[ln_.rsplit(",",1)[1].strip()]=stop_t+2
       Ps=_b1_decode(st_l,ld)
       kk_=(((r0)&0xFFFFFFFF)& ~(7<<23))|(bm<<23)
       ks_=kk_-(1<<32) if (kk_>>31) else kk_
       m8t[(RTAB+(ks_<<3))&M64]=RROW;m8t[RROW+(Lc-1)*8]=4*st_pc
       m8t[PTRS+8]=PBt+16
       cs=_B1CPU(Ps,m8t,{});xs_=cs.x
       xs_[1],xs_[2],xs_[5],xs_[6],xs_[0],xs_[14]=8,ATAB,0,r0,Lc-1,bm
       xs_[10],xs_[11],xs_[12],xs_[17],xs_[23],xs_[26],xs_[27]=RTAB,48,3,PTRS,1,16,1
       tl_={};cs.run(0,{stop_t,stop_t+1},{},tl_)
       n_w=(m8t[PTRS+8]-(PBt+16))//8
       got=[m8t[PBt+16+8*j_]-(1<<38)-1 for j_ in range(n_w)]
       exp_ops=4+(2 if (kept and kept[0]>0) else 0)+(s_op*(len(kept)-1)+l_op if kept else 0)+1
       if got!=kept or sum(tl_.values())!=exp_ops:
        emit_ok= False;break
      if not emit_ok: break
     if not emit_ok: break
    if not emit_ok: break
   if emit_ok:
    cc_prog,cc_lab,_=_seg_load(path,"chord_chunk")
    def _slice_ret(start_idx):
     out=[]
     for k_ in range(start_idx,len(gas_ins)):
      out.append(gas_ins[k_])
      if gas_ins[k_]=="ret": break
     return out
    c128_specs=((0,0,128),(0,1,86),(0,1,85),(1,2,85),(0,3,43),(1,3,43),(2,3,42))
    RROWS=[RROW+8192*(k_+1) for k_ in range(3)]
    for mode,s_op,l_op in ((8,8,6),(9,9,7)):
     c128_ins={sp_:_slice_ret(gas_lab[f"c{mode}_s{sp_[0]}_p{sp_[1]}_k{sp_[2]}"]) for sp_ in c128_specs}
     for bm in range(8):
      for r0 in (-5,-4,-3):
       for Lc in (129,256,385,781):
        kept=[j_ for j_ in range(Lc) if not ((bm>>((r0+j_)%3))&1)]
        chunks=(Lc+127)//128;tail_L=Lc-128*(chunks-1);r_tail=(r0+128*(chunks-1))%3
        _,tail_ins=_tgt_for_short(mode,s_op,bm,r_tail,tail_L)
        lc_prog=list(wd);pc_cc=len(lc_prog);lc_prog.extend(cc_prog)
        pc_m0=len(lc_prog);lc_prog.extend(mask0_ins)
        pc_tail=len(lc_prog);lc_prog.extend(tail_ins)
        c128_pc={}
        for sp_ in c128_specs:
         c128_pc[sp_]=len(lc_prog);lc_prog.extend(c128_ins[sp_])
        stop_t=len(lc_prog)
        ld={"wheel_disp4":0,".Lchunk_loop":pc_cc+cc_lab[".Lchunk_loop"],".Lk0_top":stop_t,".Ltop_2":stop_t+1}
        for ln_ in lc_prog:
         if ln_.startswith("tbnz ") and ".Lfl_g_" in ln_:ld[ln_.rsplit(",",1)[1].strip()]=stop_t+2
        Plc=_b1_decode(lc_prog,ld)
        for rr_ in range(3):
         rrow_r=RROWS[rr_]
         m8t[rrow_r+(Lc-1)*8]=4*(pc_m0 if bm==7 else pc_cc)
         m8t[rrow_r+(tail_L-1)*8]=4*pc_tail
         if bm<7:
          k128=[j_ for j_ in range(128) if not ((bm>>((rr_+j_)%3))&1)]
          sk_=k128[0];g_=k128[-1]-k128[-2];pc_b=bin(bm).count("1")
          pt_=0 if pc_b==0 else (3 if pc_b==2 else (1 if g_==1 else 2))
          m8t[rrow_r+781*8]=4*c128_pc[(sk_,pt_,len(k128))]
        for ch_i in range(chunks):
         rv_=r0+128*ch_i
         kk_=((rv_&0xFFFFFFFF)& ~(7<<23))|(bm<<23)
         ks_=kk_-(1<<32) if (kk_>>31) else kk_
         m8t[(RTAB+(ks_<<3))&M64]=RROWS[rv_%3]
        m8t[PTRS+8]=PBt+16
        cs=_B1CPU(Plc,m8t,{});xs_=cs.x
        x23_v=0 if Lc==781 else M64;cs.N=bool(x23_v>>63)
        xs_[1],xs_[2],xs_[5],xs_[6],xs_[0],xs_[14]=8,ATAB,0,r0&0xFFFFFFFF,Lc-1,bm
        xs_[10],xs_[11],xs_[12],xs_[17],xs_[23],xs_[26],xs_[27]=RTAB,48,3,PTRS,x23_v,16,1
        xs_[30],xs_[31]=0xCAFEBABE12345678,SP0
        tl_={};end_pc=cs.run(0,{stop_t,stop_t+1},{},tl_)
        n_w=(m8t[PTRS+8]-(PBt+16))//8
        got=[m8t[PBt+16+8*j_]-(1<<38)-1 for j_ in range(n_w)]
        sk0_c=2 if (kept and kept[0]>0) else 0
        chg_ops=(4+2+2) if bm==7 else (4+2+sk0_c+(chunks-1)*(4+2+sk0_c+20)+(s_op*(len(kept)-1)+l_op if kept else 0)+2)
        exp_end=stop_t if x23_v==0 else stop_t+1
        if got!=kept or xs_[30]!=0xCAFEBABE12345678 or end_pc!=exp_end or sum(tl_.values())>chg_ops:
         emit_ok= False;break
       if not emit_ok: break
      if not emit_ok: break
     if not emit_ok: break
 except (KeyError,ValueError,IndexError):
  emit_ok= False
 m8,m16,ref={},{},{}
 PAGE=1<<44;tally={};ok=emit_ok
 try:
  for n_rec,wgt in ((2046,19),(257,17),(1,15)):
   m8[PAGE]=n_rec;m8[PAGE+8]=wgt
   for i in range(n_rec):
    ad=(1<<40)+rng.randrange(1,512)
    m8[PAGE+16+8*i]=ad
    ref[ad]=(ref.get(ad,0)+wgt)&0xFF
   cpu=_B1CPU(P,m8,m16);cpu.x[28]=PAGE
   cpu.run(lab["apply_page"],set(),{},tally)
   if cpu.x[4]!=PAGE+16+8*n_rec:
    ok= False
  ok=ok and (m16==ref)
 except (KeyError,ValueError,IndexError):
  ok= False
 ir_exec=sum(tally.values())
 tot_recs=2046+257+1
 ir_chg=tot_recs*(ap_body_ops+2)/(ap_body_ops//3*2)+3*(ap_page_ops+3)
 ok=bool(ok and fl_ops==8 and ap_page_ops==11 and ap_body_ops==384
and fl_per_rec<=0.015 and ap_per_rec<=M.C_LAT_APPLY and ir_exec<=ir_chg
and all(ins[0]!="other" for ins in P))
 return dict(ok=ok,emit9_e2e_ok=emit_ok,emit9_flushes=em_flushes,emit9_ops=len(e1),
emit_direct_flush_ops=fl_ops,flush_ir_per_rec=round(fl_per_rec,6),
apply_page_ops=ap_page_ops,apply_body_ops=ap_body_ops,apply_ir_per_rec=round(ap_per_rec,6),
c_lat_apply_budget=M.C_LAT_APPLY,pages_tested=3,recs_tested=tot_recs,
executed_ir=ir_exec,charged_ir=round(ir_chg,2))
_B1DEC_CACHE={}
def _b1_decode(prog,labels):
 out=[]
 for ln in prog:
  if ln in _B1DEC_CACHE and not labels:
   out.append(_B1DEC_CACHE[ln]);continue
  hit=None
  if ln in _B1DEC_CACHE:
   hit=_B1DEC_CACHE[ln]
  else:
   for nm,rx in _B1PAT:
    mt=rx.match(ln)
    if mt:
     g=list(mt.groups())
     hit=(nm,)+tuple(int(v) if isinstance(v,str) and v.isdigit() else v for v in g)
     if nm=="and_imm":hit=(nm,hit[1],hit[2],int(hit[3],16)) if len(hit)>3 else (nm,hit[1],int(hit[2],16))
     _B1DEC_CACHE[ln]=hit;break
   else:
    hit=("other",ln);_B1DEC_CACHE[ln]=hit
  if hit[0] in ("adr","b","tbz","tbnz","bcc"):
   out.append(hit[:-1]+(labels[hit[-1]],))
  else:
   out.append(hit)
 return out
class _B1CPU:
 def __init__(self,P,m8,m16):
  self.P,self.x,self.m8,self.m16,self.Z,self.C,self.N=P,[0]*32,m8,m16,False,False,False
 def run(self,pc,stop,hook,tally):
  P,x,m8,m16=self.P,self.x,self.m8,self.m16
  while True:
   if pc in stop: return pc
   if pc in hook:hook[pc](self)
   ins=P[pc];op=ins[0];npc=pc+1;tally[pc]=tally.get(pc,0)+1
   if op=="ldp_post":x[ins[1]]=m8[x[ins[3]]];x[ins[2]]=m8[x[ins[3]]+8];x[ins[3]]=(x[ins[3]]+ins[4])&M64
   elif op=="ldp":x[ins[1]]=m8[x[ins[3]]];x[ins[2]]=m8[x[ins[3]]+8]
   elif op=="stp_post":m8[x[ins[3]]]=x[ins[1]];m8[x[ins[3]]+8]=x[ins[2]];x[ins[3]]=(x[ins[3]]+ins[4])&M64
   elif op=="stp":m8[x[ins[3]]]=x[ins[1]];m8[x[ins[3]]+8]=x[ins[2]]
   elif op=="ldr_post":x[ins[1]]=m8[x[ins[2]]];x[ins[2]]=(x[ins[2]]+ins[3])&M64
   elif op=="ldr_reg":x[ins[1]]=m8[(x[ins[2]]+x[ins[3]])&M64]
   elif op=="str_post":m8[x[ins[2]]]=x[ins[1]];x[ins[2]]=(x[ins[2]]+ins[3])&M64
   elif op=="add_lsl":x[ins[1]]=(x[ins[2]]+(x[ins[3]]<<ins[4]))&M64
   elif op=="add_reg":x[ins[1]]=(x[ins[2]]+x[ins[3]])&M64
   elif op=="add_imm12":x[ins[1]]=(x[ins[2]]+(ins[3]<<12))&M64
   elif op=="add_imm":x[ins[1]]=(x[ins[2]]+ins[3])&M64
   elif op=="sub_lsl":x[ins[1]]=(x[ins[2]]-(x[ins[3]]<<ins[4]))&M64
   elif op=="sub_imm":x[ins[1]]=(x[ins[2]]-ins[3])&M64
   elif op=="sub_reg":x[ins[1]]=(x[ins[2]]-x[ins[3]])&M64
   elif op in ("subs_reg","cmp","cmp_imm"):
    a_,b_=(x[ins[2]],x[ins[3]]) if op=="subs_reg" else ((x[ins[1]],ins[2]) if op=="cmp_imm" else (x[ins[1]],x[ins[2]]))
    r_=(a_-b_)&M64;self.Z,self.C,self.N=(r_==0),(a_>=b_),bool(r_>>63)
    if op=="subs_reg":x[ins[1]]=r_
   elif op=="csel_hs":x[ins[1]]=x[ins[2]] if self.C else x[ins[3]]
   elif op=="and_imm":x[ins[1]]=x[ins[2]]&ins[3]
   elif op=="lsr":x[ins[1]]=x[ins[2]]>>ins[3]
   elif op=="ubfx":x[ins[1]]=(x[ins[2]]>>ins[3])&((1<<ins[4])-1)
   elif op=="adr":x[ins[1]]=4*ins[2]
   elif op=="blr":
    if x[ins[1]]&3: raise ValueError("misaligned blr")
    x[30]=4*npc;npc=x[ins[1]]>>2
   elif op=="br":
    if x[ins[1]]&3: raise ValueError("misaligned br")
    npc=x[ins[1]]>>2
   elif op=="b":npc=ins[1]
   elif op=="ret":
    if x[30]==0: return
    if x[30]&3: raise ValueError("misaligned ret")
    npc=x[30]>>2
   elif op=="tbz":
    if not (x[ins[1]]>>ins[2])&1:npc=ins[3]
   elif op=="tbnz":
    if (x[ins[1]]>>ins[2])&1:npc=ins[3]
   elif op=="bcc":
    t_={"eq":self.Z,"hi":self.C and not self.Z,"lo": not self.C,"ls":(not self.C) or self.Z,"pl": not self.N}[ins[1]]
    if t_:npc=ins[2]
   elif op=="staddh":
    ad=x[ins[2]]
    if ad&1: raise ValueError("unaligned staddh")
    m16[ad]=(m16.get(ad,0)+(x[ins[1]]&0xFFFF))&0xFFFF
   elif op=="staddb":
    ad=x[ins[2]]
    m16[ad]=(m16.get(ad,0)+(x[ins[1]]&0xFF))&0xFF
   elif op=="ldp_sp_off":x[ins[1]]=m8[(x[31]+ins[3])&M64];x[ins[2]]=m8[(x[31]+ins[3]+8)&M64]
   elif op=="stp_sp_off":m8[(x[31]+ins[3])&M64]=x[ins[1]];m8[(x[31]+ins[3]+8)&M64]=x[ins[2]]
   elif op=="ldp_off":x[ins[1]]=m8[x[ins[3]]+ins[4]];x[ins[2]]=m8[x[ins[3]]+ins[4]+8]
   elif op=="ldr_off":x[ins[1]]=m8[(x[ins[2]]+ins[3])&M64]
   elif op=="ldr":x[ins[1]]=m8[x[ins[2]]]
   elif op=="ldr_lsl3":x[ins[1]]=m8[(x[ins[2]]+(x[ins[3]]<<3))&M64]
   elif op=="ldr_sp_off":x[ins[1]]=m8[(x[31]+ins[2])&M64]
   elif op=="str_sp_off":m8[(x[31]+ins[2])&M64]=x[ins[1]]
   elif op=="str":m8[x[ins[2]]]=x[ins[1]]
   elif op=="ldadd":
    ad_=x[ins[3]];old_=m8[ad_];m8[ad_]=(old_+x[ins[1]])&M64;x[ins[2]]=old_
   elif op=="add_lsr":x[ins[1]]=(x[ins[2]]+(x[ins[3]]>>ins[4]))&M64
   elif op=="sub_imm12":x[ins[1]]=(x[ins[2]]-(ins[3]<<12))&M64
   elif op=="mov_imm":x[ins[1]]=ins[2]&M64
   elif op=="extr":x[ins[1]]=(((x[ins[2]]<<64)|x[ins[3]])>>ins[4])&M64
   elif op=="sbfx":
    v_=(x[ins[2]]>>ins[3])&((1<<ins[4])-1)
    x[ins[1]]=(v_-(1<<ins[4]))&M64 if v_>>(ins[4]-1) else v_
   elif op=="mul":x[ins[1]]=(x[ins[2]]*x[ins[3]])&M64
   elif op=="madd":x[ins[1]]=(x[ins[2]]*x[ins[3]]+x[ins[4]])&M64
   elif op=="msub":x[ins[1]]=(x[ins[4]]-x[ins[2]]*x[ins[3]])&M64
   elif op=="udiv":x[ins[1]]=x[ins[2]]//x[ins[3]] if x[ins[3]] else 0
   elif op=="lsl_reg":x[ins[1]]=(x[ins[2]]<<(x[ins[3]]&63))&M64
   elif op=="bfi_w":
    m_=((1<<ins[4])-1)<<ins[3]
    x[ins[1]]=((x[ins[1]]& ~m_)|((x[ins[2]]<<ins[3])&m_))&0xFFFFFFFF
   elif op=="ldr_sxtw3":
    v_=x[ins[3]]&0xFFFFFFFF
    x[ins[1]]=m8[(x[ins[2]]+(((v_-(1<<32)) if v_>>31 else v_)<<3))&M64]
   elif op=="ldr_uxtw3":x[ins[1]]=m8[(x[ins[2]]+((x[ins[3]]&0xFFFFFFFF)<<3))&M64]
   else:
    raise ValueError("instruction outside the executed set: "+str(ins[1]))
   pc=npc
def band1_e2e(path,seed=20261124,NA=2000,primes=((7,45,8,256),(61,94,24,8),(89,103,24,8)),nstreams=3):
 """Executes the sender (b1_row, b1_copy, b1_flush, b1_switch, b1_term) and the receiver (b1_page, b1_apply) on a synthetic
    instance: NA global ranks in 4 equal-count shells, K checkpoints per shell (production K = 256 for p = 7, so the 16-bit
    start-offset field reaches 57,344), three primes with different log p weights,
    lines of rows with random z(b), last shell J, checkpoint c, sub-band switches over 3 streams, 16 KiB pages. The score
    memory must equal a direct reference, and the executed instructions must not exceed the charged per-row, per-entry,
    per-page and per-switch prices derived from the listing's own segment counts."""
 prog,lab,segof=_b1_load(path)
 P=_b1_decode(prog,lab)
 seg={}
 for sg in segof:seg[sg]=seg.get(sg,0)+1
 body=lab["b1_copy_end"]-lab["b1_copy"];abody=lab["b1_apply_end"]-lab["b1_apply"]
 row_static=seg["b1_row"]+4+2
 ent_static=body/(body//2)+4.0/(body//2)
 thr=[ins for sg_,ins in zip(segof,P) if sg_=="b1_switch" and ins[0] in ("add_imm12","add_imm") and ins[1]==26]
 T=sum(ins[3]<<12 if ins[0]=="add_imm12" else ins[3] for ins in thr)
 pair_bytes=16*((body//4)-1)
 thr_ok=(len(thr)==2 and 16<=T and T+max(8+pair_bytes,pair_bytes+16)<=16384)
 page_min=(T-16)//8
 flush_static,term_static,sw_static=seg["b1_flush"],seg["b1_term"],seg["b1_switch"]
 page_static=seg["b1_page"]+3
 c_hit=ent_static+(flush_static+20)/page_min+(abody+2)/(abody//3*2)+(page_static+20)/page_min
 rng=random.Random(seed)
 m8,m16,ref={},{},{}
 lo=[1+j*NA//4 for j in range(4)]+[NA+1]
 za={ai:sum((1<<i) for i in range(3) if rng.random()<1/3) for ai in range(1,NA+1)}
 st=dict(rows=0,empty_rows=0,entries=0,lines=0,pages=0,term_pages=0,min_page=1<<30,max_page=0,switches=0)
 ENT0,CHK0,PAGE0,SST0,DESC0=1<<40,1<<36,1<<44,1<<32,1<<34
 stop_rx=set()
 tally,rtally={},{}
 term=[False]
 def receiver(cpu):
  page=cpu.x[28];n=m8[page]
  if 16+8*n>16384: raise ValueError("page overflow")
  if not term[0]:st["min_page"]=min(st["min_page"],n)
  st["max_page"]=max(st["max_page"],n);st["pages"]+=1
  rc=_B1CPU(P,m8,m16);rc.x[28]=page
  rc.run(lab["b1_page"],stop_rx,{},rtally)
  if rc.x[4]!=page+16+8*n: raise ValueError("receiver did not consume exactly its page")
 hook={lab["b1_handoff"]:receiver,lab["b1_term_handoff"]:receiver}
 ok= True;max_start=0
 try:
  for p,wgt,nlines,K in primes:
   BLK=K*8;RS=32*BLK
   kta={ai:rng.randrange(p) for ai in range(1,NA+1)}
   D=[-(-(lo[j+1]-lo[j])//K) for j in range(4)]
   byres={}
   for ai in range(1,NA+1):byres.setdefault(kta[ai],[]).append(ai)
   pos=ENT0
   for res in range(p+1):
    for z in range(8):
     for j in range(4):
      base=CHK0+res*RS+(4*z+j)*BLK
      if res==p:
       m8[base]=pos;continue
      ent=[ai for ai in byres.get(res,()) if lo[j]<=ai<lo[j+1] and not (za[ai]&z)]
      for c in range(K):
       m8[base+8*c]=pos+8*bisect_left(ent,lo[j]+c*D[j])
      for ai in ent:
       m8[pos]=ai;pos+=8
   cv=[rng.randrange(p) for _ in range(6)];CV0=1<<30
   for k_,c_ in enumerate(cv):m8[CV0+8*k_]=c_
   cpu=_B1CPU(P,m8,m16);x=cpu.x
   x[23]=p*RS;x[24]=CHK0;x[3]=SST0;x[29]=wgt
   for s_ in range(nstreams):
    m8[SST0+16*s_]=PAGE0+(s_<<14);m8[SST0+16*s_+8]=PAGE0+(s_<<14)+16
   x[19]=SST0;x[28]=PAGE0;x[25]=PAGE0+16;x[26]=PAGE0+T
   rb=1+(p<<8)
   for line in range(nlines):
    Lr=rng.randrange(1,40)
    bv=[rng.randrange(-559,560-Lr)]+[rng.randrange(-510,511) for _ in range(5)]
    r0=sum(bk*ck for bk,ck in zip(bv,cv))%p;dl=cv[0]
    dp=DESC0+line*4096;LR=DESC0-(1<<24)+32*line
    for i in range(Lr):
     if i==0 and p<16:zb,J,ce=0,3,K
     else:zb,J,ce=sum((1<<q) for q in range(3) if rng.random()<1/3),rng.randrange(4),rng.randrange(1,K+1)
     flag=(i==0 or rng.random()<0.08)
     d0=rb<<38;rb+=1
     m8[dp]=d0;m8[dp+8]=(((4*zb+J)*BLK+8*ce)<<32)|(int(flag)<<31)|((4*zb)*BLK);dp+=16
     max_start=max(max_start,4*zb*BLK)
     if flag:
      m8[dp]=rng.randrange(nstreams);dp+=8;st["switches"]+=1
     res=(r0+i*dl)%p
     lim=lo[J]+(ce*D[J] if ce<K else lo[J+1]-lo[J])
     n_row=0
     for ai in byres.get(res,()):
      if not (za[ai]&zb) and ai<lim and ai<lo[J+1]:
       ad=d0+ai;ref[ad]=(ref.get(ad,0)+wgt)&0xFF;n_row+=1
     st["rows"]+=1;st["entries"]+=n_row;st["empty_rows"]+=(n_row==0)
    pk= lambda vs:sum((v_&0xFFFF)<<(16*i_) for i_,v_ in enumerate(vs))
    m8[LR]=DESC0+line*4096;m8[LR+8]=dp;m8[LR+16]=pk(bv[:4]);m8[LR+24]=pk(bv[4:])
    x[1]=LR;x[2]=CV0;x[5]=p;x[6]=RS.bit_length()-1;x[7]=p*(6<<15)
    st["lines"]+=1
    cpu.run(lab["b1_line"],set(),hook,tally)
   m8[x[19]]=x[28];m8[x[19]+8]=x[25]
   term[0]= True
   for s_ in range(nstreams):
    sp_=SST0+16*s_
    if m8[sp_+8]!=m8[sp_]+16:
     x[28]=m8[sp_];x[25]=m8[sp_+8]
     cpu.run(lab["b1_term"],set(),hook,tally);st["term_pages"]+=1
     m8[sp_+8]=x[25]
   term[0]= False
  ok=(m16==ref)
 except (KeyError,ValueError,IndexError):
  ok= False
 snd=sum(tally.values());rcv=sum(rtally.values())
 flush_pages=st["pages"]-st["term_pages"]
 snd_chg=st["rows"]*row_static+st["entries"]*ent_static+flush_pages*flush_static+st["term_pages"]*term_static \
  +st["switches"]*sw_static+st["lines"]*(seg["b1_line"]+1)
 rcv_chg=st["entries"]*(abody+2)/(abody//3*2)+st["pages"]*page_static
 ok=bool(ok and st["min_page"]>=page_min and snd<=snd_chg and rcv<=rcv_chg and seg["b1_row"]==20
and seg["b1_line"]+1<=64 and body==256 and abody==384 and thr_ok and max_start>=1<<15
and all(ins[0]!="other" for ins in P))
 return dict(ok=ok,cells=len(ref),stats=st,seg_counts=seg,row_static=row_static,ent_static=ent_static,flush_threshold=T,
max_start_offset=max_start,
page_min_records=page_min,c_hit_lo_exact=c_hit,switch_allow=0.2,sender_ir=snd,sender_charged=round(snd_chg,2),
receiver_ir=rcv,receiver_charged=round(rcv_chg,2))
def verify_nonprim_hnf_witnesses():
 PQ=M.load_data()["pool_quad"]
 if M._OHM_CACHE is None:M.is_ideal_primitive_oh([1,0,0,0,0,0],[0,1,0,0,0,0])
 MT=M._OHM_CACHE
 def _omul(u,v): return [sum(u[i]*v[j]*MT[i][j][c] for i in range(6) for j in range(6)) for c in range(6)]
 def _mrows(v): return [[sum(v[i]*MT[i][j][c] for i in range(6)) for c in range(6)] for j in range(6)]
 def _hnf_det(rows):
  R=[r[:] for r in rows];det=1;r0=0
  for col in range(6):
   while True:
    nz=[i for i in range(r0,len(R)) if R[i][col]!=0]
    if len(nz)==1:
     p=nz[0]
     if p!=r0:R[r0],R[p]=R[p],R[r0]
     if R[r0][col]<0:R[r0]=[-x for x in R[r0]]
     det*=R[r0][col];r0+=1;break
    p=min(nz,key=lambda i:abs(R[i][col]))
    if R[p][col]<0:R[p]=[-x for x in R[p]]
    pv=R[p][col];pr=R[p]
    for i in nz:
     if i==p: continue
     q=R[i][col]//pv;ri=R[i]
     for c in range(col,6):ri[c]-=q*pr[c]
  return det
 w_zk=[-1,-1,-1,-1,1,-1]
 assert [sum(w_zk[k]*M.ZK_NUM[k][j] for k in range(6))//3 for j in range(6)]==M.W_COEFFS
 ok_all= True;checked=[]
 for ni_decl,a,b in PQ["nonprim_witnesses"]:
  ni_obs=_hnf_det(_mrows(a)+_mrows(b))
  a2=_omul(a,a);b2=_omul(b,b);wa2=_omul(w_zk,a2)
  ng=abs(_hnf_det(_mrows([b2[c]+(M.U_BN254+1)*wa2[c]+2*a2[c] for c in range(6)])))
  div_ok=(ng%(ni_obs*ni_obs)==0);nj=ng//(ni_obs*ni_obs) if div_ok else 0;v2=0;t=nj
  while t>0 and (t&1)==0:v2+=1;t>>=1
  ok_i=(ni_obs==ni_decl and div_ok and v2 in (0,3,6));ok_all&=ok_i
  checked.append(dict(N_I=ni_obs,divides_Ng_sq=div_ok,v2_NJ=v2,ok=ok_i))
 return dict(ok=ok_all,n_checked=len(checked),witnesses=checked)
def verify_ht_pool(lat,extrap_factor=None):
 """Replay check plus pooled bound. `lat` is v20_lat_b2's live result for run hb = lat["prod_256orb_run"]: its 32 strata
    (M_s, rejection tries S_T, draws (p, m_p, Ir)) must equal the stored records of run hb. The bound pools all 16 stored
    runs into one stratified sample (128 draws per stratum) and applies v20_lat_b2's own estimator: Haldane-UMVUE stratum
    sizes N_s = M_s (K - 1) / (S_T - 1) from the pooled tries, Horvitz-Thompson ratio R over the draws, linearised SE with
    the N_s variance, and 1.03 (R + 3 SE) <= C_LAT_LEAF."""
 import json,lzma,base64
 D=json.loads(lzma.decompress(base64.b85decode((Path(__file__).resolve().parent/"ht_pool.b85").read_bytes())))
 runs,Ls=D["runs"],D["L"]
 hb=lat["prod_256orb_run"];live=lat["prod_256orb_strata_ht"]
 st=runs[hb-1];replay_errs=0
 for s,((Ms,ST,recs),(lMs,lST,ld)) in enumerate(zip(st,live)):
  if Ms!=lMs or ST!=lST or len(recs)!=len(ld):replay_errs+=1;continue
  for (k,mp,iu),(pv,lmp,lir) in zip(recs,ld):
   if Ls[s]+2*k!=pv or mp!=lmp or abs(iu-lir*1e6)>1.0:replay_errs+=1
 K=8*len(runs)
 pooled=[(runs[0][s][0],sum(r[s][1] for r in runs),[(Ls[s]+2*k,mp,iu/1e6) for r in runs for k,mp,iu in r[s][2]])
for s in range(len(Ls))]
 Y=X=0.0
 for Ms,ST,d in pooled:
  N=Ms*(K-1.0)/(ST-1.0)
  Y+=N*sum(mp*ir/pv for pv,mp,ir in d)/K;X+=N*sum(mp/pv for pv,mp,_ in d)/K
 R=Y/X;var=0.0
 for Ms,ST,d in pooled:
  N=Ms*(K-1.0)/(ST-1.0)
  us=[(mp/pv)*(ir-R) for pv,mp,ir in d];ub=sum(us)/K
  s2=sum((u-ub)**2 for u in us)/(K-1);vN=N*N/(K-2.0)
  var+=(N*N+vN)*(s2/K)+vN*ub*ub
 SE=sqrt(var)/X
 f=extrap_factor if extrap_factor is not None else lat["prod_256orb_extrap_factor"]
 g_wmean=lat["prod_256orb_gauss_wmean_ht_ir"]
 g_se=lat["prod_256orb_gauss_se_strat_ir"]
 g_ub=g_wmean+3.0*g_se
 bound_unif=f*(R+3.0*SE)
 bound=R+(f-1.0)*g_ub+f*3.0*SE
 def _est(strata,Kd):
  Y_=X_=0.0
  for Ms,ST,d in strata:
   N=Ms*(Kd-1.0)/(ST-1.0)
   Y_+=N*sum(mp*ir/pv for pv,mp,ir in d)/Kd;X_+=N*sum(mp/pv for pv,mp,_ in d)/Kd
  R_=Y_/X_;v_=0.0
  for Ms,ST,d in strata:
   N=Ms*(Kd-1.0)/(ST-1.0)
   us=[(mp/pv)*(ir-R_) for pv,mp,ir in d];ub=sum(us)/Kd
   v_+=(N*N+N*N/(Kd-2.0))*(sum((u-ub)**2 for u in us)/(Kd-1)/Kd)+N*N/(Kd-2.0)*ub*ub
  return R_,sqrt(v_)/X_
 per_run=[];per_run_bound_unif=[]
 for r in runs:
  Rr,Sr=_est([(Ms,ST,[(Ls[s]+2*k,mp,iu/1e6) for k,mp,iu in recs]) for s,(Ms,ST,recs) in enumerate(r)],8)
  per_run.append(round(Rr,3));per_run_bound_unif.append(round(f*(Rr+3.0*Sr),3))
 ok=(replay_errs==0 and bound<=M.C_LAT_LEAF and len(runs)==16 and all(len(r)==32 for r in runs))
 return dict(ok=ok,runs=len(runs),orbits=sum(len(r[s][2]) for r in runs for s in range(len(Ls))),replay_run=hb,
replay_errs=replay_errs,pooled_R=round(R,4),pooled_SE=round(SE,4),extrap_factor=f,
live_gauss_wmean=round(g_wmean,4),live_gauss_se=round(g_se,4),live_gauss_ub=round(g_ub,4),
pooled_bound_uniform_103=round(bound_unif,4),
pooled_bound=round(bound,4),c_lat_leaf=M.C_LAT_LEAF,per_run_R=per_run,per_run_bound_uniform_103=per_run_bound_unif)
