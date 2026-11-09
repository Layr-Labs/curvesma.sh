from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
P_FIELD=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
assert r<P_FIELD and r.bit_length()==256 and r%2==1
NODE_BITS=72;SIDE_BITS=36
NODES=1<<NODE_BITS;SIDE=1<<SIDE_BITS
CAND_CAP=1<<28;REQ_CAP=1<<29;RELAY_CAP=1<<34
ID_RECOVERY=1<<22
MX=1024;MP=MX//2;YB=1024;G=1600# permanent X's per node, X-pairs; Y-chains per group (one batch = one step of YB chains)
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KBH=F(513,512)*F(1<<184,r-1)
KEYLO=17;KEY=192-KEYLO;BAD=33
CLOW=28;C29=F(2)+F(2*BAD,1<<CLOW)
B45=(F(2)+F(2*BAD,1<<KEYLO))*F(1<<211,r-1)# stored K3' subset:45 bits, lowest fixed bit17
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512)if f==36 else C29 if f==72 else F(2)+F(2*BAD,1<<KEYLO)if f==KEY else 2)
KBC=lambda f,h:F(1<<(256-f),r-1)*C29
HB=58;CAP_PG=90112;SLOTS0=(1<<51)+(1<<56)
VISIT=4;VER=7;RECV=12;H2S=F(129,32);HOP1=F(291,128)
FWDP=F(152)
BWDP_LONG=F(1925,2)
BWDP=BWDP_LONG
FWDY,BWDY,FWDC,BWDC=36,7,17,27-MP//2
YCF,YCB,YPP,YTS=123,561,945,11
CY=YCF+YCB+YPP+4*YTS
INV=24725;DRV=68;YFIX=11+7+10# fixed inverse allowance (two per batch); flip_group's own work per batch (executed); fixed per call of the chain/prep kernels
TAILW=12;TAILE=F(7,2);FIRST=8
FOLDW=4;FOLDMEAN=F(108*((1<<256)-P_FIELD),r-1);assert 9*((1<<256)-P_FIELD)<P_FIELD
ARITHM=(FWDP+BWDP)/2;ARITH=ARITHM+TAILW+FOLDW
TBROW,TBFIX=140,87
GENMUL=23441408
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
def calc(L,g=G,caph=None):
 NY=g*YB*L;n=2*MX*NY;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
 NYT=G*YB*L
 cs=gam*12*(NYT+MX*L)*S*(P1*P1+F(36,(r-1)**2))
 db=YB*MX*F(2,r-1)+L*YB*F(3,r-1)+F(MX,r-1)
 mp=2*MX*YB*L*KB(36)+1;vp=mp+1+60*(YB*L+MX)*KB(36)+4*(L*L+MX*L)*KB(36)
 assert CAP_PG-12>mp
 dp=vp/(CAP_PG-12-mp)**2
 mh=NODES*2*MX*NYT*KBH+1
 if caph is None:caph=-(-mh.numerator//mh.denominator)+(1<<41)
 assert caph-12>mh
 dh=(mh+1+60*(NYT+MX)*KBH+4*(L*L+MX*L)*KBH)/(caph-12-mh)**2
 intr=RHO*(S+18)*KB(KEY)
 h=RHO*(S*S*KB(72)*(2*KB(KEY))+(4*S+20)*(2*KB(KEY)))
 gg=2*RHO*n*S*(2*KB(KEY))
 v=18*RHO*(F(S,SIDE)*S*KB(36)*(2*KB(KEY))+(4*S+20)*(2*KB(KEY)))
 dc=h/(CAND_CAP-4);dr=gg/(REQ_CAP-8);dv=v/(RELAY_CAP-108)
 per=db+dp+dh+dc+dr+dv
 loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
 glob=NODES*G*L*db+NODES*F(MX,r)+F(NODES*(MX+2*G*YB+1),1<<248)
 M=(mu.numerator<<200)//mu.denominator;mu_lo=F(M,1<<200)
 acc=0;t=1;fk=1
 for k in range(EXP_TERMS+1):
  if k:t*=-M;fk*=k
  acc+=t*(1<<(200*(EXP_TERMS-k)))*(FACT80//fk)
 ED=(1<<(200*EXP_TERMS))*FACT80;EH=-(-acc*(1<<160)//ED)
 assert acc*(1<<160)<=EH*ED<acc*(1<<160)+ED
 e_hi=F(EH,1<<160)
 assert 0<mu_lo<8
 dlab=gam*F(4*BAD<<(256-KEYLO),(r-1)**2)
 id_loss=6*mu*(1+F(1,r))*RHO*B45
 success=1-e_hi-cs-loss-glob-dlab-id_loss
 return dict(n=n,S=S,mu=mu,success=success,caph=caph,mp=mp,cs=cs,intr=intr,loss=loss,global_loss=glob,
                label_loss=dlab,id_recovery_loss=id_loss,h=h,g=gg,v=v,d_cl=dc,d_rq=dr,d_rl=dv,d_batch=db,dp=dp,dh=dh)
with localcontext()as ctx:
 ctx.prec=80
 THRESH=F(39,100)
 def fits(m):
  mh=NODES*2*MX*(G*YB*m)*KBH+1;ch=-(-mh.numerator//mh.denominator)+(1<<41)
  return SLOTS0+32+32*CAP_PG*SIDE+32*ch<=7<<61
 lo,hi=1<<20,1<<30
 while hi-lo>1:
  m=(lo+hi)//2
  if fits(m):lo=m
  else:hi=m
 L=55000000;assert fits(L);G0=G
 assert G*L<1<<37
 YMAX=YB*(G*L-1)+(YB-1);assert YMAX<1<<47 and((YMAX>>28)<<45)<1<<64
 assert SIDE_BITS+28==64 and SIDE_BITS+9+19==64 and MX==2*(1<<9)
 assert ((SIDE-1)|(((YMAX&((1<<28)-1))<<36)))<1<<64
 assert ((SIDE-1)|((MX//2-1)<<36)|((YMAX>>28)<<45))<1<<64
 a=calc(L);caph=a['caph'];assert fits(L)and caph==a['caph']
 K0=-((-39*(1<<64)*a['success'].denominator)//(100*a['success'].numerator));assert F(K0,1<<64)*a['success']>=THRESH>F(K0-1,1<<64)*a['success']
 GATE_K=K0;RHO_G=F(GATE_K,1<<64);GATE_COST=64
 assert 0<GATE_K<1<<64 and RHO_G*a['success']>=THRESH
 up=lambda x:F(-((-x.numerator<<120)//x.denominator),1<<120)
 PHI=[F(1)]+[min(F(1),up(max(F(0),1-calc(L,g,caph)['success'])))for g in range(1,G)]# PHI[g]=1-prefix success (section15)
 assert all(PHI[g]>=PHI[g+1]for g in range(G-1))and all(0<=p<=1 for p in PHI)
 log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
 CELL128=KBC(HB+NODE_BITS,0)
 PAIRM=CELL128+F(1,r-1)
 PAIRM_SAME=2*CELL128+F(1,r-1)
 S=a['S'];s=NODES*2*MX*YB*L# leaves per group over all nodes
 INC=YB*MX
 sx=NODES*L*(2*INC)*(2*INC-1)//2
 leaf=NODES*L*(ARITHM*INC+(FWDY+BWDY+CY)*YB+FWDC+BWDC+YFIX+2*INV+DRV)
 route=s*(HOP1+H2S+RECV)
 WG=[];walk_main=F(0);gp=0
 for g in range(G):
  evis=(s*gp+s*(s-1)//2-sx)*PAIRM
  WG.append(leaf+route+VISIT*evis);walk_main+=PHI[g]*VISIT*evis;gp+=s
 assert gp==S
 PAIRS=S*(S-1)//2
 SXV=VISIT*PAIRM_SAME*sx*G
 EK2=VER*PAIRS*2*KB(KEY)# K3' (full 175-bit key) matches
 DEF=2*(VISIT+VER)*a['d_batch']*caph*S# defective walkers, and defective priors (each visited by at most CAP_H later records)
 rows={
      'leaf_route_walks_E':sum(PHI[g]*WG[g]for g in range(G)),
      'same_batch_walks':SXV,'k2_verifications':EK2,'defective_walks':DEF,
      'group_fixed':sum(PHI)*NODES*(1236*SIDE+27),# per executed group (section15)
      'hash_head_zeroing':NODES*((1<<(HB-14))*1026+1),
      'node_init':NODES*(1<<20),
      'gen_points_headers':NODES*(MX*GENMUL+G*YB*2*GENMUL+GENMUL),# X's, chain starts and steps, R_g
      'scalar_sampling':NODES*(MX+2*G*YB+1)*256,
      'resolution':NODES*CAND_CAP*42410496,
      'source_service':NODES*REQ_CAP*((1<<10)+GENMUL+(1<<16)+ID_RECOVERY),# each capped request regenerates its Y
      'phase4_stages':sum(PHI)*NODES*(4*(1<<41)+(1<<37)),
      'barriers':sum(PHI)*NODES*4*(1<<10),
      'tables':NODES*(4*MP*(TBROW*(1<<37)+TBFIX+(1<<12))+MX*(1<<17)),# tbld, table constants, X areas
      'phase_control':NODES*max(1<<20,1024*G),
      'validation_output':(1<<28)+NODES*64+G*2*(SIDE+NODES)*40,
      'correction_tails':sum(PHI)*NODES*L*INC*TAILE+NODES*G*L*INC*TAILW*a['d_batch'],# section16: non-defective PHI-weighted; defective batches unweighted at the worst
      'fold_tails':sum(PHI)*NODES*L*INC*FOLDMEAN+NODES*G*L*INC*FOLDW*a['d_batch'],
      'first_prefix_tail':sum(PHI)*NODES*L*FIRST,
    }
 E_RUN=sum(rows.values());E_exact=GATE_COST+RHO_G*E_RUN;E=-(-E_exact.numerator//E_exact.denominator)
 TW=NODES*G*(L*(INC*(TAILW+FOLDW)+FIRST)+1236*SIDE+27+4*(1<<41)+(1<<37)+4*(1<<10));NT=('leaf_route_walks_E','correction_tails','first_prefix_tail','fold_tails','group_fixed','phase4_stages','barriers')
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
 CLAIM=F(-((-s_hi.numerator*10**15)//s_hi.denominator),10**15)
 E_NL=E_exact+RHO_G*sum(PHI[g]for g in range(G))*NODES*L*INC*F(27+96,4);E_NL=-(-E_NL.numerator//E_NL.denominator)
 nl_hi=score_bounds(E_NL)[1];CLAIM_NL=F(-((-nl_hi.numerator*10**15)//nl_hi.denominator),10**15)
 assert s_hi<=CLAIM and s_lo>CLAIM-F(1,10**15) and s_hi<score_bounds(WORST)[0]
 assert PHI[-1]<F(8,1000)and sum(PHI)<F(35,100)*G and log2_upper(F(WORST,1000))<F(128,1)
 assert RW_UP<=F(RUN_WORST_LOG,10**4)and a['success']>F(992,1000)
 setup=sum(rows[k]for k in['group_fixed','hash_head_zeroing','node_init','gen_points_headers','scalar_sampling','phase_control','barriers','phase4_stages','tables'])
 assert setup<=1000*(1<<119)and log2_upper(F(setup,1000))<F(119,1)
 slotbytes=32*CAP_PG*SIDE
 assert SLOTS0+32+slotbytes+32*caph<=7<<61 and(7<<61)+(8<<HB)==1<<64
 assert CAP_PG%128==0
 print(json.dumps({'research_only':True,'L':L,'mx':MX,'yb':YB,'node_bits':NODE_BITS,'side_bits':SIDE_BITS,'candidate_cap':CAND_CAP,'request_cap':REQ_CAP,'relay_cap':RELAY_CAP,'id_recovery_instructions':ID_RECOVERY,'groups':G,'cap_pg':CAP_PG,'g0':G0,'hash_bits':HB,'key_bits':KEY,'bad_residues':BAD,
       'chains_per_node':G*YB,'y_per_node':G*YB*L,'label_loss':str(dec(a['label_loss']))[:12],'id_recovery_loss':str(dec(a['id_recovery_loss'])),'partial_key_bits':45,'lowest_cell_bit':CLOW,'exp_terms':EXP_TERMS,'d_cl':str(dec(a['d_cl']))[:12],
       'success_formula':str(dec(a['success']))[:24],'gate_k':GATE_K,'gated_success':str(dec(RHO_G*a['success']))[:24],'E_run_log2_units':str(log(F(E_RUN,1000)))[:20],'hop1_group_drop':'%.6e'%float(a['dp']),'home_drop':'%.6e'%float(a['dh']),'intruder_rate':'%.6e'%float(a['intr']),'label_loss_e':'%.6e'%float(a['label_loss']),
       'phi':[str(dec(PHI[g]))[:10]for g in(0,1,G//2,G-1)],'E_factor_groups':str(dec(sum(PHI)/G))[:12],
       'walk_instructions_per_leaf_E':str(dec(walk_main/S))[:10],'E_total_instructions':str(E),'E_log2_units':str(log(F(E,1000)))[:40],
       'phi1_log2_units':str(log(F(WORST,1000)))[:20],'run_worst_log2_units_ceil4':str(dec(F(RUN_WORST_LOG,10**4))),'claim':'%d/10^15'%(CLAIM*10**15),
       'score_enclosure':[str(dec(s_lo))[:40],str(dec(s_hi))[:40]],'claim_no_ls64':'%d/10^15'%(CLAIM_NL*10**15),'gap_no_ls64':'%d/10^15'%((CLAIM_NL-CLAIM)*10**15),'fold_tail_mean':str(FOLDMEAN),'fold_tail_worst':FOLDW,'arith_per_incidence':str(ARITHM),'per_y':CY+FWDY+BWDY,
       'tails':{'worst_per_incidence':TAILW,'charged_mean':float(TAILE),'first_per_batch':FIRST},'label_interval':[-21,12],'rows':{k:str(v)for k,v in rows.items()},'cap_h':caph,'memory_end':1<<64,
       'setup_log2_units':str(log(F(setup,1000)))[:12],'scope':'Grouped E-ledger; coefficients in proof.md.'}))
