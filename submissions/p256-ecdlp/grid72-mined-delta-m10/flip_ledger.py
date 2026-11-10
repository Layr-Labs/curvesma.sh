"""Exact arithmetic only. Does not independently prove the supplied symbolic inequalities.
Flipped-grid grouped early-abort E-ledger (G groups of YB Y-chains, MX permanent X's, 175-bit fraction keys): proof.md sections 3-7, 15."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551# protected subgroup order
P_FIELD=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
assert r<P_FIELD and r.bit_length()==256 and r%2==1
NODE_BITS=72;SIDE_BITS=36
NODES=1<<NODE_BITS;SIDE=1<<SIDE_BITS
CAND_CAP=1<<27;REQ_CAP=1<<29;RELAY_CAP=1<<34
ID_RECOVERY=1<<22
MX=1024;MP=MX//2;YB=1024;G=1600# permanent X's per node, X-pairs; Y-chains per group (one batch = one step of YB chains)
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KBH=F(513,512)*F(1<<184,r-1)
KEYLO=17;KEY=192-KEYLO;BAD=33# key = bits [17,192) of the label (175 bits); 33 residues mod 2^17 within the label error range [-21,12]
CLOW=28;C29=F(2)+F(2*BAD,1<<CLOW)# conservative for every routing/hash-cell pattern, lowest fixed fraction bit28
B45=(F(2)+F(2*BAD,1<<KEYLO))*F(1<<211,r-1)# stored K3' subset:45 bits, lowest fixed bit17
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512)if f==36 else C29 if f==72 else F(2)+F(2*BAD,1<<KEYLO)if f==KEY else 2)
KBC=lambda f,h:F(1<<(256-f),r-1)*C29
HB=58;CAP_PG=2300032;SLOTS0=(1<<51)+(1<<56)# hop1 slots follow the 2^56-byte pair tables at [2^51,2^51+2^56)
VISIT=4;HIT=7;VER=18;RECV=6+F(27,128);H2S=7+F(1,32)+F(10,128);HOP1=F(291,128)# hop2 blocks: proof-tails.md Section12, RECEIVER_PROTOCOL section3
GF=1277*SIDE+28;HOLES=G*SIDE# per-group fixed work; one 32-byte hole per odd final block at a home
FWDP=F(152)
# 1925/2: 1929/2 minus the two mov x17,#32 that used to sit in the pair.
# +1/4: the pair-end pointer is reloaded once per four unrolled pairs, on every iteration.
# The four-instruction borrow stays inside this count on every pair. No stopping-weighted credit.
BWDP_LONG=F(1861,2)+F(1,4)
BWDP=BWDP_LONG
FWDY,BWDY,FWDC,BWDC=36,13,17,24+99+4# per Y and per call: backward Y pays two more (str end, mov x19,#32)
YCF,YCB,YPP,YTS=123,561,1055,12# per Y lane, tail-free (measured): ychain_f, ychain_b, yprep; their tail sites
CY=YCF+YCB+YPP+4*YTS# per-Y chain step and Y entry with every tail charged
INV=24725;DRV=68;YFIX=11+7+10# fixed inverse allowance (two per batch); flip_group's own work per batch (executed); fixed per call of the chain/prep kernels
TAILW=12;TAILE=F(7,2);FIRST=8# three sites per incidence (worst 12, H_T mean at most 6); two positions per batch charged in full
FOLDW=4;FOLDMEAN=F(108*((1<<256)-P_FIELD),r-1);assert 9*((1<<256)-P_FIELD)<P_FIELD
ARITHM=(FWDP+BWDP)/2;ARITH=ARITHM+TAILW+FOLDW
TBROW,TBFIX=140,87# tbld per row (8 slots, slot 7 zero) and per table (measured); 2^37 rows
GENMUL=23441408# COLD_ACCOUNTING generator item
def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def ln_interval(x,terms=64):
 assert 1<=x<=2
 t=(x-1)/(x+1);t2=t*t;power=t;lo=F(0)
 for j in range(terms):
  lo+=2*power/(2*j+1);power*=t2
 tail=2*power/((2*terms+1)*(1-t2))
 return lo,lo+tail
def score_bounds(total):
 e=124;x=F(total,1000*(1<<e))
 while x>2:e+=1;x=x/2
 while x<1:e-=1;x=x*2
 a,b=ln_interval(x);c,d=ln_interval(F(2))
 return F(e)+a/d,F(e)+b/c
import math
EXP_TERMS=80;FACT80=math.factorial(EXP_TERMS)
# Phase 0 and the public X's (proof.md sections 2, 5, 6; proof-appendices.md PHASE0)
M_C=7<<101# candidates [s]P,7*2^29 per node
LQ,LN=128,4;NCLS=LQ*(LQ//LN);NEED=LQ*LQ;SHARE=2*LN*LN;assert MX==2*LN*LQ# line design: classes of different slopes share SHARE X's
QP=F(1<<66,r-1)# PHASE0, Lemma M
C1=M_C*M_C*QP
F2X=(2*C1+C1*C1+M_C*C1*F(1<<68,r-1))/(MX*(MX-1));F2=-(-F2X.numerator//F2X.denominator)# two leaves through one Y point
PATS=9**4;F4X=32*PATS*sum(C1**k for k in range(1,5));F4=-(-F4X.numerator//F4X.denominator);P_REL=F(1,1<<150);assert F4X/(r-1)<P_REL
SQP=math.isqrt(P_FIELD)+1;S_A=F(r-1,2)
P_OR_LO=((1<<64)-1)*(F(P_FIELD,4)-3*SQP-(1<<64))/(S_A*S_A);P_OR_HI=((1<<64)-1)*(F(P_FIELD,4)+3*SQP+4)/(S_A*S_A)
EN_LO,EN_HI=M_C*(M_C-1)*P_OR_LO,M_C*(M_C-1)*P_OR_HI;VAR_N=EN_HI+M_C**3*QP*QP;TS=1<<10
E_STR=M_C*M_C*F(2*((1<<154)+1)<<64,r-1)*QP# close pairs split by a2^102 cell boundary
E_Z=M_C*M_C*F(1,r-1)+M_C**3*QP*QP# equal x, or two close pairs sharing a candidate
MCAP,HREC=1<<16,1<<35# hop1 mining slot cap; mining records per home
OVF1,OVF2=3*(M_C//NODES)*F(2<<220,r-1)/MCAP,3*M_C*F(2<<184,r-1)/HREC;assert OVF1<F(1,1<<15)and OVF2<F(2,3)and HREC*F(58,100)>1000;P_OVF=F(1,1<<800)# binomial tails (3nq/cap)^cap
MHB,HCAP,LCAP,PCAP=34,48,64,1<<20
QD,QH=F(2<<(184-MHB),r-1),F(2<<184,r-1)
P_CH=math.comb(M_C,HCAP+2)*QD**(HCAP+1);P_LS=(M_C*M_C*QP)**LCAP*QH**(LCAP-1)/math.factorial(LCAP);P_PC=VAR_N/(PCAP-EN_HI)**2
P_AB0X=VAR_N/(EN_LO-NEED-TS)**2+E_STR/TS+E_Z+P_OVF+P_CH+P_LS+P_PC;P_AB0=F(-(-(P_AB0X.numerator<<80)//P_AB0X.denominator),1<<80)
assert EN_LO>NEED+TS and PCAP>EN_HI and SQP*SQP>P_FIELD and M_C//NODES<1<<33
MINE_VISIT=32;MINE_REC=2048;MINE_REC_ITEMS=120+12+24+4+HCAP*MINE_VISIT+2
MINE_NODE=1<<48;MINE_COLL=PCAP*(20*24+8+8*40+80)+(1<<15)*8
MINE_NODE_ITEMS=2*GF+(1<<(MHB-14))*1026+1+8*LCAP**2+LCAP*4*40+2*40+MINE_COLL+(1<<15)*88+(1<<16)+8*1024
assert MINE_REC_ITEMS<=MINE_REC and MINE_NODE_ITEMS<MINE_NODE
MINE_W=M_C*(GENMUL+MINE_REC)+NODES*MINE_NODE
assert M_C//NODES*32+PCAP*64+(1<<15)*32+LCAP*16<=1<<56 and HREC*64<1<<57# Phase 0 scratch fits the unbuilt table region
FIXED=4*MX# the leaves through the witness's two Y points
def calc(L,g=G,caph=None):
 """Success lower bound after the first g groups (chains of L steps); caps fixed by the full run when given."""
 NY=g*YB*L;n=2*MX*NY;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
 npc=NODES//NCLS;keep=lambda sh:1-F(sh*(sh+1),2*MX*MX)
 pr_same=NCLS*F(npc*(npc-1),2);pr_slope=LQ*F((LQ//LN)*(LQ//LN-1),2)*npc*npc;pr_diff=F(LQ*(LQ-1),2)*(LQ//LN)**2*npc*npc
 assert pr_same+pr_slope+pr_diff==F(NODES*(NODES-1),2)
 gam_p=n*n*(pr_same*keep(MX)+pr_slope+pr_diff*keep(SHARE))# Gamma': pairs with both X's shared and j>=j' dropped (section5)
 NYT=G*YB*L;WC=4*MX*L# leaves on the witness's two chains
 cs=gam*12*(NYT+MX*L)*S*(P1*P1+F(36,(r-1)**2))
 db=YB*MX*F(2,r-1)+L*YB*F(3,r-1)+F(MX,r-1)
 TWP,TWH=1<<20,1<<40# Markov allowances for the witness chain's own records
 mp=2*MX*YB*L*KB(36)+1;vp=mp+1+YB*L*(2*MX)**2*F2*KB(36)**2
 assert CAP_PG-FIXED-TWP>mp
 dp=vp/(CAP_PG-FIXED-TWP-mp)**2+RHO*2*MX*L*KB(36)/TWP
 mh=NODES*2*MX*NYT*KBH+1
 if caph is None:caph=-(-mh.numerator//mh.denominator)+(1<<41)+HOLES
 assert caph-HOLES-129-FIXED-TWH>mh
 dh=(mh+1+NODES*NYT*(2*MX)**2*F2*KBH**2)/(caph-HOLES-129-FIXED-TWH-mh)**2+RHO*WC*KBH/TWH
 intr=RHO*(S+FIXED*F2)*KB(KEY)
 Mk=2*KB(KEY);TI=1<<20;TIV=1<<25
 e_o=RHO*Mk*((FIXED+WC)*(S+FIXED+WC)+FIXED*FIXED*F4)# first moment (section6)
 h=RHO*S*S*KB(72)*Mk;gg=RHO*n*S*Mk;v=RHO*F(S,SIDE)*S*KB(36)*Mk# pair counts of free leaves (section6)
 T2=lambda m:2*m*RHO*Mk*S+2*m*2*MX+2*m*2*MX*S*RHO*Mk+m*(2*MX)**2# shared leaf (third free; third on its partner's Y point), one same-Y couple, four leaves on two Y points (no product bound: charged in full)
 th,tg,tv=T2(h),T2(gg),T2(v)
 dc=(h*h+h+th)/(CAND_CAP-4-TI)**2+e_o/TI;dr=(gg*gg+gg+tg)/(F(REQ_CAP-8,2)-TI)**2+e_o/TI;dv=(v*v+v+tv)/(F(RELAY_CAP-108,18)-TIV)**2+e_o/TIV
 per=db+dp+dh+dc+dr+dv
 loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
 glob=NODES*G*L*db+F(NODES*(2*G*YB+1)+M_C,1<<248)+P_AB0+P_REL
 dlab=gam*F(4*BAD<<(256-KEYLO),(r-1)**2)# bad-label collision mass (section6)
 mu_g=gam_p*P1-dlab# good collisions of Gamma'
 M=(mu_g.numerator<<200)//mu_g.denominator;mu_lo=F(M,1<<200)
 acc=0;t=1;fk=1# e_hi = sum_{k<=80} (-mu_lo)^k/k!, summed exactly over the common denominator 2^(200*80)*80!
 for k in range(EXP_TERMS+1):
  if k:t*=-M;fk*=k
  acc+=t*(1<<(200*(EXP_TERMS-k)))*(FACT80//fk)
 ED=(1<<(200*EXP_TERMS))*FACT80;EH=-(-acc*(1<<160)//ED)
 assert acc*(1<<160)<=EH*ED<acc*(1<<160)+ED
 e_hi=F(EH,1<<160)
 assert 0<mu_lo<16
 id_loss=6*mu*(1+F(1,r))*RHO*F2*B45
 success=1-e_hi-cs-loss-glob-id_loss
 return dict(n=n,S=S,mu=mu,mu_g=mu_g,th=th,tg=tg,tv=tv,success=success,caph=caph,mp=mp,cs=cs,intr=intr,loss=loss,global_loss=glob,
                bad_label_mass=dlab,id_recovery_loss=id_loss,h=h,g=gg,v=v,d_cl=dc,d_rq=dr,d_rl=dv,d_batch=db,dp=dp,dh=dh,e_o=e_o)
with localcontext()as ctx:
 ctx.prec=80
 THRESH=F(39,100)
 def fits(m):# home records (32 bytes) end below the heads at7*2^61
  mh=NODES*2*MX*(G*YB*m)*KBH+1;ch=-(-mh.numerator//mh.denominator)+(1<<41)+HOLES
  return SLOTS0+32+32*CAP_PG*SIDE+48*ch<=7<<61
 lo,hi=1<<20,1<<30
 while hi-lo>1:
  m=(lo+hi)//2
  if fits(m):lo=m
  else:hi=m
 L=67500000;assert fits(L);G0=G# fixed run geometry; every global round fits the47-bit Y index
 assert G*L<1<<37# every group stores the global round rho=g*L+u; no home-index epoch bit is needed
 YMAX=YB*(G*L-1)+(YB-1);assert YMAX<1<<47 and((YMAX>>28)<<45)<1<<64
 assert SIDE_BITS+28==64 and SIDE_BITS+9+19==64 and MX==2*(1<<9)
 assert ((SIDE-1)|(((YMAX&((1<<28)-1))<<36)))<1<<64
 assert ((SIDE-1)|((MX//2-1)<<36)|((YMAX>>28)<<45))<1<<64# w2 stores source row36+Ylo28; w3 stores source col36+Xpair9+Yhi19
 a=calc(L);caph=a['caph'];assert fits(L)and caph==a['caph']
 K0=-((-39*(1<<64)*a['success'].denominator)//(100*a['success'].numerator));assert F(K0,1<<64)*a['success']>=THRESH>F(K0-1,1<<64)*a['success']
 GATE_K=K0;RHO_G=F(GATE_K,1<<64);GATE_COST=64# section17 gate: run iff u<GATE_K; least64-bit integer for this changed run
 assert 0<GATE_K<1<<64 and RHO_G*a['success']>=THRESH
 up=lambda x:F(-((-x.numerator<<120)//x.denominator),1<<120)
 PHI=[F(1)]+[min(F(1),up(max(F(0),1-calc(L,g,caph)['success'])))for g in range(1,G)]# PHI[g]=1-prefix success (section15)
 for g in range(1,G):PHI[g]=min(PHI[g],PHI[g-1])# Pr(T>g-1) is nonincreasing in g, so an earlier prefix's bound also bounds it
 assert all(PHI[g]>=PHI[g+1]for g in range(G-1))and all(0<=p<=1 for p in PHI)
 log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
 CELL128=KBC(HB+NODE_BITS,0)
 CELL136=KBC(HB+NODE_BITS+6,0)# a K2 hit: all 64 bits of K2 and the home
 PAIRM=CELL128+F(1,r-1)+F(HIT,VISIT)*(CELL136+F(1,r-1))# main pairs: one fixed cell, +1/(r-1) for an identity leaf; the K2-hit cost is folded in per visit
 PAIRM_SAME=2*CELL128+F(1,r-1)+F(HIT,VISIT)*(2*CELL136+F(1,r-1))# same-batch pairs: two-cell union
 S=a['S'];s=NODES*2*MX*YB*L# leaves per group over all nodes
 INC=YB*MX# incidences per batch
 sx=NODES*L*(2*INC)*(2*INC-1)//2# same-batch leaf pairs inside one group (section15)
 leaf=NODES*L*(ARITHM*INC+(FWDY+BWDY+CY)*YB+FWDC+BWDC+YFIX+2*INV+DRV)
 route=s*(HOP1+H2S+RECV)
 WG=[];walk_main=F(0);gp=0
 for g in range(G):
  evis=(s*gp+s*(s-1)//2-sx)*PAIRM
  WG.append(leaf+route+VISIT*evis);walk_main+=PHI[g]*VISIT*evis;gp+=s
 assert gp==S
 PAIRS=S*(S-1)//2
 SAMEY=NODES*G*L*YB*(2*MX)*(2*MX-1)//2# unordered leaf pairs through one Y point: Lemma M's factor F2 (section15)
 SXV=VISIT*PAIRM_SAME*(sx*G+(F2-1)*SAMEY)
 EK2=VER*(PAIRS+(F2-1)*SAMEY)*2*KB(KEY)# K3' (full 175-bit key) matches
 DEF=2*(VISIT+VER)*a['d_batch']*caph*S# defective walkers, and defective priors (each visited by at most CAP_H later records)
 rows={
      'leaf_route_walks_E':sum(PHI[g]*WG[g]for g in range(G)),
      'same_batch_walks':SXV,'k2_verifications':EK2,'defective_walks':DEF,
      'group_fixed':sum(PHI)*NODES*GF,# per executed group (section15)
      'hash_head_zeroing':NODES*((1<<(HB-14))*1026+1),
      'node_init':NODES*(1<<20),
      'gen_points_headers':NODES*(MX*GENMUL+G*YB*2*GENMUL+GENMUL),# X's, chain starts and steps, R_g
      'scalar_sampling':(NODES*(2*G*YB+1)+M_C)*256,
      'phase0':MINE_W,
      'resolution':NODES*CAND_CAP*42410496,
      'source_service':NODES*2*CAND_CAP*((1<<10)+GENMUL+(1<<16)+ID_RECOVERY),# lifetime: the cumulative candidate cap gives at most 2*N*CAND_CAP requests over the run; each served request regenerates its Y
      'phase4_stages':sum(PHI)*NODES*(4*(1<<41)+(1<<37)),
      'barriers':sum(PHI)*NODES*4*(1<<10),
      'tables':NODES*(3*MP*(TBROW*(1<<37)+TBFIX+(1<<12))+MX*(1<<17)),# tbld, table constants, X areas
      'phase_control':NODES*max(1<<20,1024*G),
      'validation_output':(1<<28)+NODES*64+G*2*(SIDE+NODES)*40,
      'correction_tails':sum(PHI)*NODES*L*INC*TAILE+NODES*G*L*INC*TAILW*a['d_batch'],# section16: non-defective PHI-weighted; defective batches unweighted at the worst
      'fold_tails':sum(PHI)*NODES*L*INC*FOLDMEAN+NODES*G*L*INC*FOLDW*a['d_batch'],
      'first_prefix_tail':sum(PHI)*NODES*L*FIRST,
    }
 E_RUN=sum(rows.values());E_exact=GATE_COST+RHO_G*E_RUN;E=-(-E_exact.numerator//E_exact.denominator)# gated (section17)
 TW=NODES*G*(L*(INC*(TAILW+FOLDW)+FIRST)+GF+4*(1<<41)+(1<<37)+4*(1<<10));NT=('leaf_route_walks_E','correction_tails','first_prefix_tail','fold_tails','group_fixed','phase4_stages','barriers')
 WORST_exact=sum(WG)+TW+sum(v for k,v in rows.items()if k not in NT);WORST=-(-WORST_exact.numerator//WORST_exact.denominator)
 RUN_WORST=sum(WG[g]-VISIT*(s*s*g+s*(s-1)//2-sx)*PAIRM for g in range(G))+S*caph*(VISIT+VER)+TW+sum(v for k,v in rows.items()if k not in NT+('same_batch_walks','k2_verifications','defective_walks'))
 RUN_WORST=-(-RUN_WORST.numerator//RUN_WORST.denominator)
 assert RUN_WORST>WORST
 def log2_upper(x):
  e=x.numerator.bit_length()-x.denominator.bit_length();y=x/F(2)**e
  while y>=2:y/=2;e+=1
  while y<1:y*=2;e-=1
  return e+ln_interval(y)[1]/ln_interval(F(2))[0]
 RW_UP=log2_upper(F(RUN_WORST,1000));RUN_WORST_LOG=-((-RW_UP.numerator*10**4)//RW_UP.denominator)
 s_lo,s_hi=score_bounds(E)
 CLAIM=F(-((-s_hi.numerator*10**15)//s_hi.denominator),10**15)# least 15-decimal ceiling of the certified upper bound
 # disclosed sensitivity (proof.md section4): the same rows with the leaf passes' LS64 accesses as LDP/STP (+27/2 forward, +48 backward per pair)
 E_NL=E_exact+RHO_G*sum(PHI[g]for g in range(G))*NODES*L*INC*F(27+96,4);E_NL=-(-E_NL.numerator//E_NL.denominator)
 nl_hi=score_bounds(E_NL)[1];CLAIM_NL=F(-((-nl_hi.numerator*10**15)//nl_hi.denominator),10**15)
 assert s_hi<=CLAIM and s_lo>CLAIM-F(1,10**15) and s_hi<score_bounds(WORST)[0]
 assert PHI[-1]<F(8,1000)and sum(PHI)<F(35,100)*G and log2_upper(F(WORST,1000))<F(128,1)# printed bounds (section15)
 assert RW_UP<=F(RUN_WORST_LOG,10**4)and a['success']>F(992,1000)
 setup=sum(rows[k]for k in['group_fixed','hash_head_zeroing','node_init','gen_points_headers','scalar_sampling','phase_control','barriers','phase4_stages','tables','phase0'])
 assert setup<=1000*(1<<120)and log2_upper(F(setup,1000))<F(120,1)# printed preprocessing bound (section7)
 slotbytes=32*CAP_PG*SIDE
 assert SLOTS0+32+slotbytes+48*caph<=7<<61 and(7<<61)+(8<<HB)==1<<64
 assert CAP_PG%128==0
 print(json.dumps({'research_only':True,'L':L,'mx':MX,'yb':YB,'node_bits':NODE_BITS,'side_bits':SIDE_BITS,'candidate_cap':CAND_CAP,'request_cap':REQ_CAP,'relay_cap':RELAY_CAP,'id_recovery_instructions':ID_RECOVERY,'groups':G,'cap_pg':CAP_PG,'g0':G0,'hash_bits':HB,'key_bits':KEY,'bad_residues':BAD,
       'chains_per_node':G*YB,'y_per_node':G*YB*L,'bad_label_mass':str(dec(a['bad_label_mass']))[:12],'id_recovery_loss':str(dec(a['id_recovery_loss'])),'partial_key_bits':45,'lowest_cell_bit':CLOW,'exp_terms':EXP_TERMS,'d_cl':str(dec(a['d_cl']))[:12],'mu':str(dec(a['mu']))[:12],'mu_g':str(dec(a['mu_g']))[:12],'loads':{k:'%.6e'%float(a[k])for k in('h','g','v')},'drops':{k:'%.6e'%float(a[k])for k in('d_cl','d_rq','d_rl','dp','dh')},'second_moment_rest':{k:'%.3e'%float(a[k])for k in('th','tg','tv')},
       'success_formula':str(dec(a['success']))[:24],'gate_k':GATE_K,'gated_success':str(dec(RHO_G*a['success']))[:24],'E_run_log2_units':str(log(F(E_RUN,1000)))[:20],'hop1_group_drop':'%.6e'%float(a['dp']),'home_drop':'%.6e'%float(a['dh']),'intruder_rate':'%.6e'%float(a['intr']),'bad_label_mass_e':'%.6e'%float(a['bad_label_mass']),
       'phi':[str(dec(PHI[g]))[:10]for g in(0,1,G//2,G-1)],'E_factor_groups':str(dec(sum(PHI)/G))[:12],
       'walk_instructions_per_leaf_E':str(dec(walk_main/S))[:10],'E_total_instructions':str(E),'E_log2_units':str(log(F(E,1000)))[:40],
       'phi1_log2_units':str(log(F(WORST,1000)))[:20],'run_worst_log2_units_ceil4':str(dec(F(RUN_WORST_LOG,10**4))),'claim':'%d/10^15'%(CLAIM*10**15),
       'score_enclosure':[str(dec(s_lo))[:40],str(dec(s_hi))[:40]],'claim_no_ls64':'%d/10^15'%(CLAIM_NL*10**15),'gap_no_ls64':'%d/10^15'%((CLAIM_NL-CLAIM)*10**15),'fold_tail_mean':str(FOLDMEAN),'fold_tail_worst':FOLDW,'arith_per_incidence':str(ARITHM),'per_y':CY+FWDY+BWDY,
       'tails':{'worst_per_incidence':TAILW,'charged_mean':float(TAILE),'first_per_batch':FIRST},'label_interval':[-21,12],'rows':{k:str(v)for k,v in rows.items()},'cap_h':caph,'memory_end':1<<64,
       'setup_log2_units':str(log(F(setup,1000)))[:12],'scope':'Grouped E-ledger; coefficients in proof.md.'}))
