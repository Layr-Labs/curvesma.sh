"""Exact arithmetic only. Does not independently prove the supplied symbolic inequalities.
Flipped-grid grouped early-abort E-ledger (G groups of YB Y-chains, MX permanent X's, 176-bit fraction keys): proof.md sections 3-7, 15."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551# protected subgroup order
P_FIELD=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
assert r<P_FIELD and r.bit_length()==256 and r%2==1
NODES=1<<70;SIDE=1<<35
MX=1024;MP=MX//2;YB=1024;G=526# permanent X's per node, X-pairs; Y-chains per group (one batch = one step of YB chains)
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KBH=F(513,512)*F(1<<186,r-1)
KEY=176;BAD=33# key = bits [16,192) of the label; 33 residues mod 2^16 within the label error range [-21,12]
C29=F(2)+F(2*BAD,1<<29)# counting bound 2(1+33*2^-a) at lowest fixed bit a=29 (FRACTION_LABEL.md section3)
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512)if f==35 else C29 if f==70 else F(2)+F(2*BAD,1<<16)if f==KEY else 2)
KBC=lambda f,h:F(1<<(256-f),r-1)*C29
HB=60;CAP_PG=5*(1<<17)
VISIT=4;VER=5;RECV=11;H2S=F(129,32);HOP1=F(291,128)
FWDP=F(359,2);BWDP=F(1981,2)# per X-pair, tail-free (measured): flip_fwd, flip_bwd
FWDY,BWDY,FWDC,BWDC=12,11,18,24# per Y and per call (measured)
YCF,YCB,YPP,YTS=123,561,945,11# per Y lane, tail-free (measured): ychain_f, ychain_b, yprep; their tail sites
CY=YCF+YCB+YPP+4*YTS# per-Y chain step and Y entry with every tail charged
INV=39479;DRV=68;YFIX=11+7+10# fixed inverse allowance (two per batch); flip_group's own work per batch (executed); fixed per call of the chain/prep kernels
TAILW=12;TAILE=6;FIRST=8# three sites per incidence (worst 12, H_T mean at most 6); two positions per batch charged in full
ARITHM=(FWDP+BWDP)/2;ARITH=ARITHM+TAILW
TBROW,TBFIX=140,126# tbld per row (8 slots) and per table (measured)
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
 e=125;x=F(total,1000*(1<<e))
 if x>2:e=126;x=x/2
 a,b=ln_interval(x);c,d=ln_interval(F(2))
 return F(e)+a/d,F(e)+b/c
import math
FACT40=math.factorial(40)
def calc(L,g=G,caph=None):
 """Success lower bound after the first g groups (chains of L steps); caps fixed by the full run when given."""
 NY=g*YB*L;n=2*MX*NY;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
 NYT=G*YB*L
 cs=gam*12*(NYT+MX*L)*S*(P1*P1+F(36,(r-1)**2))
 db=YB*MX*F(2,r-1)+L*YB*F(3,r-1)+F(MX,r-1)# per batch: a zero x-difference, a degenerate chain point or step of its group up to its round, or an identity X
 mp=2*MX*YB*L*KB(35)+1;vp=mp+1+60*(YB*L+MX)*KB(35)+4*(L*L+MX*L)*KB(35)# witness incidences; witness-chain pairs (section6)
 assert CAP_PG-12>mp
 dp=vp/(CAP_PG-12-mp)**2
 mh=NODES*2*MX*NYT*KBH+1
 if caph is None:caph=-(-mh.numerator//mh.denominator)+(1<<50)
 assert caph-12>mh
 dh=(mh+1+60*(NYT+MX)*KBH+4*(L*L+MX*L)*KBH)/(caph-12-mh)**2
 intr=RHO*(S+12)*KB(KEY)
 h=RHO*(S*S*KB(70)*(2*KB(KEY))+(2*S+20)*(2*KB(KEY)))
 gg=2*RHO*n*S*(2*KB(KEY))
 v=18*RHO*(F(S,SIDE)*S*KB(35)*(2*KB(KEY))+(2*S+20)*(2*KB(KEY)))
 dc=h/((1<<24)-4);dr=gg/((1<<29)-8);dv=v/((1<<30)-16)
 per=db+dp+dh+dc+dr+dv
 loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
 glob=NODES*G*L*db+NODES*F(MX,r)+F(NODES*(MX+2*G*YB+1),1<<248)
 M=(mu.numerator<<200)//mu.denominator;mu_lo=F(M,1<<200)
 acc=0;t=1;fk=1# e_hi = sum_{k<=40} (-mu_lo)^k/k!, summed exactly over the common denominator 2^8000*40!
 for k in range(41):
  if k:t*=-M;fk*=k
  acc+=t*(1<<(200*(40-k)))*(FACT40//fk)
 e_hi=F(acc,(1<<8000)*FACT40)
 assert 0<mu_lo<2
 dlab=gam*F(4*BAD<<240,(r-1)**2)
 success=1-e_hi-cs-loss-glob-dlab
 return dict(n=n,S=S,mu=mu,success=success,caph=caph,mp=mp,cs=cs,intr=intr,loss=loss,global_loss=glob,
                label_loss=dlab,h=h,g=gg,v=v,d_cl=dc,d_rq=dr,d_rl=dv,d_batch=db,dp=dp,dh=dh)
with localcontext()as ctx:
 ctx.prec=80
 THRESH=F(39,100)
 lo,hi=1<<26,1<<28
 while hi-lo>1:
  m=(hi+lo)//2
  if calc(m)['success']>=THRESH:hi=m
  else:lo=m
 L=hi;a=calc(L);caph=a['caph']
 assert calc(L)['success']>=THRESH and calc(L-1)['success']<THRESH
 up=lambda x:F(-((-x.numerator<<120)//x.denominator),1<<120)
 PHI=[F(1)]+[min(F(1),up(max(F(0),1-calc(L,g,caph)['success'])))for g in range(1,G)]# PHI[g]=1-prefix success (section15)
 assert all(PHI[g]>=PHI[g+1]for g in range(G-1))and all(0<=p<=1 for p in PHI)
 log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
 CELL130=KBC(HB+70,0)
 PAIRM=CELL130+F(1,r-1)# main pairs: one fixed cell, +1/(r-1) for an identity leaf
 PAIRM_SAME=2*CELL130+F(1,r-1)# same-batch pairs: two-cell union
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
 SXV=VISIT*PAIRM_SAME*sx*G
 EK2=VER*PAIRS*2*KBC(134,0)
 DEF=2*(VISIT+VER)*a['d_batch']*caph*S# defective walkers, and defective priors (each visited by at most CAP_H later records)
 rows={
      'leaf_route_walks_E':sum(PHI[g]*WG[g]for g in range(G)),
      'same_batch_walks':SXV,'k2_verifications':EK2,'defective_walks':DEF,
      'group_fixed':NODES*G*(1236*SIDE+27),
      'hash_head_zeroing':NODES*((1<<(HB-14))*1026+1),
      'node_init':NODES*(1<<20),
      'gen_points_headers':NODES*(MX*GENMUL+G*YB*2*GENMUL+GENMUL),# X's, chain starts and steps, R_g
      'scalar_sampling':NODES*(MX+2*G*YB+1)*256,
      'resolution':NODES*(1<<24)*42410496,
      'source_service':NODES*(1<<29)*((1<<10)+GENMUL+(1<<16)),# each capped request regenerates its Y
      'phase4_stages':NODES*G*4*(1<<41),
      'barriers':NODES*(4*G)*(1<<10),
      'tables':NODES*(4*MP*(TBROW*(1<<32)+TBFIX+(1<<12))+MX*(1<<17)),# tbld, table constants, X areas
      'phase_control':NODES*(1<<20),
      'validation_output':(1<<28)+NODES*64+G*2*(SIDE+NODES)*40,
      'correction_tails':NODES*G*L*INC*(TAILE+a['d_batch']*(TAILW-TAILE)),
      'first_prefix_tail':NODES*G*L*FIRST,# section16 positions, no PHI
    }
 E_exact=sum(rows.values());E=-(-E_exact.numerator//E_exact.denominator)
 TW=NODES*G*L*INC*TAILW;NT=('leaf_route_walks_E','correction_tails')
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
 # disclosed sensitivity (proof.md section4): the same rows with the leaf passes' LS64 accesses as LDP/STP (+11 forward, +51 backward per pair)
 E_NL=E_exact+sum(PHI[g]for g in range(G))*NODES*L*INC*F(11+51,2);E_NL=-(-E_NL.numerator//E_NL.denominator)
 nl_hi=score_bounds(E_NL)[1];CLAIM_NL=F(-((-nl_hi.numerator*10**15)//nl_hi.denominator),10**15)
 assert s_hi<=CLAIM and s_lo>CLAIM-F(1,10**15) and s_hi<score_bounds(WORST)[0]
 assert PHI[-1]<F(61115,10**5)and sum(PHI)<F(85733,10**5)*G and log2_upper(F(WORST,1000))<F(1258410,10**4)# printed bounds (section15)
 assert RUN_WORST_LOG==1782022 and RW_UP<=F(RUN_WORST_LOG,10**4)and a['success']>F(39000000016,10**11)
 setup=sum(rows[k]for k in['group_fixed','hash_head_zeroing','node_init','gen_points_headers','scalar_sampling','phase_control','barriers','phase4_stages','tables'])
 assert setup<=1000*(1<<119)and log2_upper(F(setup,1000))<F(11651305,10**5)# printed preprocessing bound (section7)
 slotbytes=32*CAP_PG*SIDE
 assert(1<<52)+32+slotbytes+40*caph<=1<<63 and(1<<63)+(8<<HB)==1<<64
 assert CAP_PG%128==0
 print(json.dumps({'research_only':True,'L':L,'mx':MX,'yb':YB,'groups':G,'cap_pg':CAP_PG,'hash_bits':HB,'key_bits':KEY,'bad_residues':BAD,
       'chains_per_node':G*YB,'y_per_node':G*YB*L,'label_loss':str(dec(a['label_loss']))[:12],'d_cl':str(dec(a['d_cl']))[:12],
       'success_formula':str(dec(a['success']))[:24],'hop1_group_drop':'%.6e'%float(a['dp']),'home_drop':'%.6e'%float(a['dh']),'intruder_rate':'%.6e'%float(a['intr']),'label_loss_e':'%.6e'%float(a['label_loss']),
       'phi':[str(dec(PHI[g]))[:10]for g in(0,1,G//2,G-1)],'E_factor_groups':str(dec(sum(PHI)/G))[:12],
       'walk_instructions_per_leaf_E':str(dec(walk_main/S))[:10],'E_total_instructions':str(E),'E_log2_units':str(log(F(E,1000)))[:40],
       'phi1_log2_units':str(log(F(WORST,1000)))[:20],'run_worst_log2_units_ceil4':str(dec(F(RUN_WORST_LOG,10**4))),'claim':'%d/10^15'%(CLAIM*10**15),
       'score_enclosure':[str(dec(s_lo))[:40],str(dec(s_hi))[:40]],'claim_no_ls64':'%d/10^15'%(CLAIM_NL*10**15),'gap_no_ls64':'%d/10^15'%((CLAIM_NL-CLAIM)*10**15),'arith_per_incidence':str(ARITHM),'per_y':CY+FWDY+BWDY,
       'tails':{'worst_per_incidence':TAILW,'charged_mean':TAILE,'first_per_batch':FIRST},'label_interval':[-21,12],'rows':{k:str(v)for k,v in rows.items()},'cap_h':caph,'memory_end':1<<64,
       'setup_log2_units':str(log(F(setup,1000)))[:12],'scope':'Grouped E-ledger; coefficients in proof.md.'}))
