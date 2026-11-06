"""Exact arithmetic only. Does not independently prove the supplied symbolic inequalities.
Grouped early-abort E-ledger (G groups, hash-chain detection, 176-bit fraction keys): proof.md sections 3, 4, 6, 7, 15."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551# protected subgroup order
P_FIELD=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff# field modulus, distinct
assert r<P_FIELD and r.bit_length()==256 and r%2==1
NODES=1<<70;SIDE=1<<35;Q2=1<<28
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KBH=F(513,512)*F(1<<186,r-1)# home70 capacity bound
KEY=176;BAD=29# key = bits [16,192) of the 192-bit fraction label; 29 = residues mod 2^16 within the error range [-19,10]
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512)if f==35 else F(2)+F(1,1<<23)if f==70 else F(2)+F(2*BAD,1<<16)if f==KEY else 2)
KBC=lambda f,h:F(1<<(256-f),r-1)*(F(2)+F(1,1<<23))# one f-bit cell whose lowest fixed bit is >= 29
G=526;HB=57;KS=1<<29;CAP_PG=5*(1<<17);CAP_P=25*(1<<18)
VISIT=4;VER=5;RECV=11;H2S=F(129,32);HOP1=F(291,128);ARITH=F(2614,4)# section4, section12
TAILW=20;TAILE=10;ARITHM=ARITH-TAILW# section16: minimal path + correction tails (5 sites, worst 20, H_T mean<=10)
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
 while x>2:e+=1;x=x/2
 a,b=ln_interval(x);c,d=ln_interval(F(2))
 return F(e)+a/d,F(e)+b/c
def split(n,k):
 base,rem=divmod(n,k);return[base+1 if g<rem else base for g in range(k)]
def group_sizes(q):# section17: groups 0..G1-1 hold rows [0,KS), the rest rows [KS,q)
 G1=(G*KS+q//2)//q;qs=split(KS,G1)+split(q-KS,G-G1)
 ms=[0]
 for x in qs:ms.append(ms[-1]+x)
 return qs,ms
def calc(q,caps=None):
 """Success lower bound for a run (or prefix) of q X-rows per node; caps fixed by the full run when given."""
 n=2*q*Q2;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
 cs=gam*12*(q+Q2)*S*(P1*P1+F(36,(r-1)**2))
 db=Q2*F(2,r-1)+F(1,r)
 qg=caps['qg']if caps else max(group_sizes(q)[0])
 mp=2*qg*Q2*KB(35)+1;vp=mp+1+60*(qg+Q2)*KB(35)
 assert CAP_PG-12>mp
 dp=vp/(CAP_PG-12-mp)**2
 mh=S*KBH+1;caph=caps['caph']if caps else-(-mh.numerator//mh.denominator)+(1<<50)
 assert caph-12>mh
 dh=(mh+1+60*(q+Q2)*KBH)/(caph-12-mh)**2
 intr=RHO*(S+12)*KB(KEY)
 h=RHO*(S*S*KB(70)*(2*KB(KEY))+(2*S+20)*(2*KB(KEY)))
 g=2*RHO*n*S*(2*KB(KEY))
 v=18*RHO*(F(S,SIDE)*S*KB(35)*(2*KB(KEY))+(2*S+20)*(2*KB(KEY)))
 dc=h/((1<<24)-4);dr=g/((1<<29)-8);dv=v/((1<<30)-16)
 per=db+dp+dh+dc+dr+dv
 loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
 assert F((1<<256)-r,1<<256)<F(1,1<<31)
 glob=NODES*q*db+F(NODES*(q+Q2+1),1<<248)
 mu_lo=F((mu.numerator<<200)//mu.denominator,1<<200)
 term=F(1);e_hi=term
 for k in range(1,41):
  term*=(-mu_lo)/k;e_hi+=term
 assert 0<mu_lo<4
 dlab=gam*F(4*BAD<<240,(r-1)**2)# a pair is a collision whose class fraction T(x) lies in the BAD residues mod 2^16
 success=1-e_hi-cs-loss-glob-dlab
 return dict(n=n,S=S,mu=mu,success=success,caph=caph,mp=mp,cs=cs,intr=intr,loss=loss,global_loss=glob,
                label_loss=dlab,h=h,g=g,v=v,d_cl=dc,d_rq=dr,d_rl=dv,d_batch=db,dp=dp,dh=dh)
GATE_COST=64# section17: one fresh 64-bit gate word u; run iff u<GATE_K
with localcontext()as ctx:
 ctx.prec=80
 THRESH=F(39,100)
 def pre(m):n=2*m*Q2;S=NODES*n;mh=S*KBH+1;return F(S*S-NODES*n*n,2)*P1,-(-mh.numerator//mh.denominator)+(1<<50)
 fits=lambda m:(lambda mu,c:mu<4 and(1<<52)+32+32*CAP_PG*SIDE+40*c<=15<<60)(*pre(m))# records below the heads at15*2^60
 lo,hi=KS,KS*2
 while hi-lo>1:
  m=(lo+hi)//2
  if fits(m):lo=m
  else:hi=m
 hi=lo;assert fits(hi)and not fits(hi+1)and KS<hi<2*KS# largest q1 that fits; k mod2^29 plus the group-boundary bit identifies every row
 a=calc(hi);q1=hi;caps=dict(caph=a['caph'],qg=max(group_sizes(hi)[0]))
 GATE_K=-((-39*(1<<64)*a['success'].denominator)//(100*a['success'].numerator));RHO_G=F(GATE_K,1<<64);assert 0<GATE_K<(1<<64)
 assert calc(q1,caps)['success']==a['success']and RHO_G*a['success']>=THRESH and F(GATE_K-1,1<<64)*a['success']<THRESH# least K
 qs,ms=group_sizes(q1);assert sum(qs)==q1 and max(qs)==caps['qg']and min(qs)>=1 and KS in ms and q1-KS<=KS
 # PHI[g]=1-prefix success (section15)
 up=lambda x:F(-((-x.numerator<<120)//x.denominator),1<<120)
 PHI=[F(1)]+[min(F(1),up(max(F(0),1-calc(ms[g],caps)['success'])))for g in range(1,G)]
 assert all(PHI[g]>=PHI[g+1]for g in range(G-1))and all(0<=p<=1 for p in PHI)
 log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
 CELL130=KBC(HB+70,0)# 129-bit cell: home 70 + K2 index 59
 PAIRM=CELL130+F(1,r)# main pairs: one fixed cell, +1/r for X=O
 PAIRM_SAME=2*CELL130+F(1,r)# same-row pairs: two-cell union
 hot=ARITH/2+HOP1+H2S+RECV
 assert hot==F(44039,128)
 assert ARITHM==F(2534,4)and TAILE*2==TAILW# each of 5 tail sites charged at probability 1/2
 S=a['S'];sg=[NODES*2*q*Q2 for q in qs]
 WG=[];samex=[];gp=0;walk_main=F(0)
 for g,q in enumerate(qs):
  s=sg[g]
  leaf=NODES*q*(ARITHM*Q2+39479+1024+(1<<17))
  route=s*(HOP1+H2S+RECV)
  sx=NODES*q*(2*Q2)*(2*Q2-1)//2
  evis=(s*gp+s*(s-1)//2-sx)*PAIRM
  WG.append(leaf+route+VISIT*evis);samex.append(sx);walk_main+=PHI[g]*VISIT*evis;gp+=s
 assert gp==S and sum(sg)==S
 PAIRS=S*(S-1)//2
 SXV=VISIT*PAIRM_SAME*sum(samex)
 EK2=VER*PAIRS*2*KBC(134,0)# all pairs, two-cell union
 DEF=(VISIT+VER)*a['d_batch']*a['caph']*S
 rows={
      'leaf_route_walks_E':sum(PHI[g]*WG[g]for g in range(G)),
      'same_X_walks':SXV,'k2_verifications':EK2,'defective_walks':DEF,
      'group_fixed':NODES*G*(1236*SIDE+27),
      'hash_head_zeroing':NODES*((1<<(HB-14))*1026+1),
      'node_init':NODES*(1<<20),
      'gen_points_headers':NODES*(q1*23441408+(Q2)*23396352+23441408),
      'scalar_sampling':NODES*(q1+Q2+1)*256,
      'resolution':NODES*(1<<24)*42410496,
      'source_service':NODES*(1<<29)*(1<<10),
      'phase4_stages':NODES*G*4*(1<<41),
      'barriers':NODES*(4*G)*(1<<10),
      'fraction_tables':NODES*(q1*(1<<17)+Q2*(1<<16)),# per-X area and per-Y fraction entry (COLD_ACCOUNTING.md)
      'range_tree_driver':NODES*8192*Q2,
      'norm_sort_search':NODES*(1280*q1+3968*Q2+4352),
      'phase_control':NODES*(1<<20),
      'validation_output':(1<<28)+NODES*64+G*2*(SIDE+NODES)*40,
      'correction_tails':NODES*q1*Q2*(TAILE+a['d_batch']*(TAILW-TAILE)),
      'first_prefix_tail':NODES*q1*4,# section16 exception, no PHI
    }
 E_RUN=sum(rows.values());E_exact=GATE_COST+RHO_G*E_RUN;E=-(-E_exact.numerator//E_exact.denominator)# integer ceiling; gated (section17)
 TW=NODES*q1*Q2*TAILW;NT=('leaf_route_walks_E','correction_tails')
 WORST_exact=sum(WG)+TW+sum(v for k,v in rows.items()if k not in NT);WORST=-(-WORST_exact.numerator//WORST_exact.denominator)
 # all groups run, whole-chain walks
 RUN_WORST=sum(WG[g]-VISIT*(sg[g]*(sum(sg[:g]))+sg[g]*(sg[g]-1)//2-samex[g])*PAIRM for g in range(G))+S*a['caph']*(VISIT+VER)+TW+sum(v for k,v in rows.items()if k not in NT+('same_X_walks','k2_verifications','defective_walks'))
 assert RUN_WORST.denominator==1 and RUN_WORST>WORST
 def log2_upper(x):# certified log2 upper bound
  e=x.numerator.bit_length()-x.denominator.bit_length();y=x/F(2)**e
  while y>=2:y/=2;e+=1
  while y<1:y*=2;e-=1
  return e+ln_interval(y)[1]/ln_interval(F(2))[0]
 assert PHI[-1]<F(13682,10**5)and sum(PHI[g]*qs[g]for g in range(G))<F(59526,10**5)*q1 and log2_upper(F(WORST,1000))<F(1270259,10**4)
 RW_UP=log2_upper(RUN_WORST/1000);RUN_WORST_LOG=-((-RW_UP.numerator*10**4)//RW_UP.denominator)# outward 4-decimal ceiling
 assert RW_UP<=F(RUN_WORST_LOG,10**4)and RUN_WORST_LOG==1802438
 s_lo,s_hi=score_bounds(E)
 CLAIM=F(125111070365338857,10**15)
 assert s_hi<=CLAIM and s_lo>CLAIM-F(1,10**15)
 assert s_hi<score_bounds(WORST)[0]
 assert F(1,(1<<34)-108)<=F(1,(1<<30)-16)
 setup=sum(rows[k]for k in['group_fixed','hash_head_zeroing','node_init','gen_points_headers','scalar_sampling','range_tree_driver','norm_sort_search','phase_control','barriers','phase4_stages','fraction_tables'])
 assert setup<=1000*(1<<119)
 slotbytes=32*CAP_PG*SIDE;caph=a['caph']
 fixed=(1<<51)+2*(1<<39)+2*(1<<41)+hi*96+Q2*96+Q2*256+512+(1<<34)+(1<<44)+(8192*SIDE+32)+(1<<38)+8*SIDE+64*SIDE+4096+4096# + Y fraction entries, X area
 assert fixed<(1<<52)
 assert(1<<52)+32+slotbytes+40*caph<=15<<60# slots, then 40-byte records
 memory_end=(15<<60)+(8<<HB);assert memory_end==1<<64# heads at 15*2^60
 assert(1<<52)+32+slotbytes+32*(2*max(qs)*Q2)<(1<<64)
 assert CAP_PG%128==0
 assert fixed+(1<<30)+2048+192*Q2+32*hi+(1<<40)<(1<<52)
 print(json.dumps({'research_only':True,'q1':hi,'q2':Q2,'groups':G,'group_sizes':[qs[0],qs[-1]],'cap_pg':CAP_PG,'hash_bits':HB,'key_bits':KEY,'label_loss':str(dec(a['label_loss']))[:12],'d_cl':str(dec(a['d_cl']))[:12],
       'success_formula':str(dec(a['success']))[:24],'gate_k':GATE_K,'gated_success':str(dec(RHO_G*a['success']))[:24],'E_run_instructions':str(-(-E_RUN.numerator//E_RUN.denominator)),'hop1_group_drop':str(dec(a['dp']))[:12],'home_drop':str(dec(a['dh']))[:12],'intruder_rate':str(dec(a['intr']))[:12],
       'hot_per_leaf_fixed':str(hot),'phi':[str(dec(PHI[g]))[:10]for g in(0,1,G//2,G-1)],'E_factor_groups':str(dec(sum(PHI[g]*qs[g]for g in range(G))/q1))[:12],
       'walk_instructions_per_leaf_E':str(dec(walk_main/S))[:10],'E_total_instructions':str(E),'E_log2_units':str(log(F(E,1000)))[:40],
       'phi1_log2_units':str(log(F(WORST,1000)))[:20],'run_worst_log2_units_ceil4':str(dec(F(RUN_WORST_LOG,10**4))),'claim':str(CLAIM),'score_enclosure':[str(dec(s_lo))[:40],str(dec(s_hi))[:40]],
       'tails':{'worst_per_incidence':TAILW,'charged_mean':TAILE,'arith_min':str(ARITHM)},'rows':{k:str(v)for k,v in rows.items()},'layout_checks':True,'memory_end':memory_end,'cap_h':a['caph'],
       'scope':'Grouped E-ledger; coefficients in proof.md.'}))
