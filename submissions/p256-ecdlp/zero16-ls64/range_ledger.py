"""Exact arithmetic only. Does not independently prove the supplied symbolic inequalities.
Grouped early-abort E-ledger (G groups, hash-chain detection): proof.md sections 4, 6, 7, 15."""
from fractions import Fraction as F
from decimal import Decimal,localcontext
import json
r=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551  # protected subgroup order
P_FIELD=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff  # field modulus, distinct
assert r<P_FIELD and r.bit_length()==256 and r%2==1
NODES=1<<70;SIDE=1<<35;Q2=1<<28
P1=F(2*(r-1)+1,r*r);RHO=F(r,r-1)**4
KBH=F(513,512)*F(1<<186,r-1)  # home70 capacity bound; pair loads keep crude KB(70)
KB=lambda f:F(1<<(256-f),r-1)*(F(513,512) if f==35 else F(2)+F(1,1<<30) if f==70 else F(2)+F(1,1<<14) if f in (144,182) else 2)
KBC=lambda f,h:F((2<<(256-f))+(4<<(224-f+h)),r-1)  # one f-bit cell, h fixed positions in [224,256)
G=512;HB=60;CAP_PG=5*(1<<17);CAP_P=25*(1<<18)
VISIT=6;VER=4;RECV=13;H2S=F(145,32);HOP1=F(291,128);ARITH=F(2930,4)  # LS64 transport (section4)
TAILW=30;TAILE=15;ARITHM=ARITH-TAILW  # section16: minimal path + correction tails (worst 30, H_T mean<=15)
def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def ln_interval(x, terms=64):
 assert 1 <= x <= 2
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
def group_sizes(q):
 base,rem=divmod(q,G);qs=[base+1 if g<rem else base for g in range(G)]
 ms=[0]
 for x in qs:ms.append(ms[-1]+x)
 return qs,ms
def calc(q,caps=None):
 """Success lower bound for a run (or prefix) of q X-rows per node; caps fixed by the full run when given."""
 n=2*q*Q2;S=NODES*n;gam=F(S*S-NODES*n*n,2);mu=gam*P1
 cs=gam*12*(q+Q2)*S*(P1*P1+F(36,(r-1)**2))
 db=Q2*F(2,r-1)+F(1,r)
 qg=caps['qg'] if caps else -(-q//G)
 mp=2*qg*Q2*KB(35)+1;vp=mp+1+60*(qg+Q2)*KB(35)
 assert CAP_PG-12>mp
 dp=vp/(CAP_PG-12-mp)**2
 mh=S*KBH+1;caph=caps['caph'] if caps else -(-mh.numerator//mh.denominator)+(1<<50)
 assert caph-12>mh
 dh=(mh+1+60*(q+Q2)*KBH)/(caph-12-mh)**2
 intr=RHO*(S+12)*KB(182)
 h=RHO*(S*S*KB(70)*(2*KB(182))+(2*S+20)*(2*KB(182)))
 g=2*RHO*n*S*(2*KB(182))
 v=18*RHO*(F(S,SIDE)*S*KB(35)*(2*KB(182))+(2*S+20)*(2*KB(182)))
 dc=h/((1<<24)-4);dr=g/((1<<29)-8);dv=v/((1<<30)-16)
 inc=q*Q2;mean=inc*(F(1,2)+F(1,1<<116))
 cap=caps['carry_cap'] if caps else -(-mean.numerator//mean.denominator)+(1<<50)
 var=F(inc*(q+Q2-1),4);anchored=4*(q+Q2)
 assert cap-mean>anchored
 dm=var/(cap-mean-anchored)**2
 per=db+dp+dh+dc+dr+dv+dm
 loss=mu*(1+F(1,r))*(F(1,r)+intr+2*per)
 assert F((1<<256)-r,1<<256)<F(1,1<<31)
 glob=NODES*q*db+F(NODES*(q+Q2+1),1<<248)
 mu_lo=F((mu.numerator<<200)//mu.denominator,1<<200)
 term=F(1);e_hi=term
 for k in range(1,41):
  term*=(-mu_lo)/k;e_hi+=term
 assert 0<mu_lo<2
 success=1-e_hi-cs-loss-glob-gam*F(1<<227,r*r)
 return dict(n=n,S=S,mu=mu,success=success,caph=caph,mp=mp,cs=cs,intr=intr,loss=loss,global_loss=glob,
                carry_cap=cap,carry_loss=dm,carry_mean=mean,carry_variance=var,anchored_incidence_bound=anchored,d_batch=db,dp=dp,dh=dh)
with localcontext() as ctx:
 ctx.prec=80
 THRESH=F(390000001,10**9)
 lo,hi=1<<28,1<<29
 while hi-lo>1:
  m=(hi+lo)//2
  if calc(m)['success']>=THRESH:hi=m
  else:lo=m
 assert hi<(1<<29)
 a=calc(hi);q1=hi;caps=dict(caph=a['caph'],carry_cap=a['carry_cap'],qg=-(-hi//G))
 assert calc(q1,caps)['success']>=THRESH and calc(q1-1,caps)['success']<THRESH
 qs,ms=group_sizes(q1);assert sum(qs)==q1 and max(qs)==-(-q1//G) and min(qs)>=1
 # PHI[g] = 1 - prefix success at ms[g] rows (proof.md 15)
 up=lambda x:F(-((-x.numerator<<120)//x.denominator),1<<120)  # dyadic upper rounding
 PHI=[F(1)]+[min(F(1),up(max(F(0),1-calc(ms[g],caps)['success']))) for g in range(1,G)]
 assert all(PHI[g]>=PHI[g+1] for g in range(G-1)) and all(0<=p<=1 for p in PHI)
 log=lambda x:dec(F(x)).ln()/Decimal(2).ln()
 CELL130=KBC(HB+70,0)  # one 130-bit cell: home 70 bits + the 60 hash-index bits of K2
 PAIRM=CELL130+F(1,r)  # main pairs: one fixed cell, +1/r for X=O
 PAIRM_SAME=2*CELL130+F(1,r)  # same-row pairs: two-cell union
 hot=ARITH/2+HOP1+H2S+RECV
 assert hot==F(49415,128)
 assert ARITHM==F(2810,4) and TAILE*2==TAILW  # each of 7 tail sites charged at probability 1/2
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
 EK2=VER*PAIRS*2*KBC(134,0)  # all pairs, two-cell union
 DEF=(VISIT+VER)*a['d_batch']*a['caph']*S
 rows={
      'leaf_route_walks_E':sum(PHI[g]*WG[g] for g in range(G)),
      'same_X_walks':SXV,'k2_verifications':EK2,'defective_walks':DEF,
      'group_fixed':NODES*G*(1280*SIDE+64),
      'hash_head_zeroing':NODES*((1<<(HB-12))*1026+1),
      'node_init':NODES*(1<<20),
      'gen_points_headers':NODES*(q1*23441408+(Q2)*23396352+23441408),
      'scalar_sampling':NODES*(q1+Q2+1)*256,
      'resolution':NODES*(1<<24)*42410496,
      'source_service':NODES*(1<<29)*(1<<10),
      'phase4_stages':NODES*G*4*(1<<41),
      'barriers':NODES*2048*(1<<10),
      'xsum_corrections':NODES*4*a['carry_cap'],
      'range_tree_driver':NODES*8192*Q2,
      'norm_sort_search':NODES*(1280*q1+3968*Q2+4352),
      'phase_control':NODES*(1<<20),
      'validation_output':(1<<28)+NODES*64+G*2*(SIDE+NODES)*40,
      'correction_tails':NODES*q1*Q2*(TAILE+a['d_batch']*(TAILW-TAILE)),  # all rows, no PHI
    }
 E_exact=sum(rows.values());E=-(-E_exact.numerator//E_exact.denominator)  # integer ceiling
 TW=NODES*q1*Q2*TAILW;NT=('leaf_route_walks_E','correction_tails')
 WORST_exact=sum(WG)+TW+sum(v for k,v in rows.items() if k not in NT);WORST=-(-WORST_exact.numerator//WORST_exact.denominator)
 # worst case: all G groups run, every walk visits its whole chain
 RUN_WORST=sum(WG[g]-VISIT*(sg[g]*(sum(sg[:g]))+sg[g]*(sg[g]-1)//2-samex[g])*PAIRM for g in range(G))+S*a['caph']*(VISIT+VER)+TW+sum(v for k,v in rows.items() if k not in NT+('same_X_walks','k2_verifications','defective_walks'))
 assert RUN_WORST.denominator==1 and RUN_WORST>WORST
 def log2_upper(x):  # certified log2 upper bound, rational x>=1
  e=x.numerator.bit_length()-x.denominator.bit_length();y=x/F(2)**e
  while y>=2:y/=2;e+=1
  while y<1:y*=2;e-=1
  return e+ln_interval(y)[1]/ln_interval(F(2))[0]
 RW_UP=log2_upper(RUN_WORST/1000);RUN_WORST_LOG=-((-RW_UP.numerator*10**4)//RW_UP.denominator)  # outward 4-decimal ceiling
 assert RW_UP<=F(RUN_WORST_LOG,10**4) and RUN_WORST_LOG==1783505
 s_lo,s_hi=score_bounds(E)
 CLAIM=F(125888933322107017,10**15)
 assert s_hi<=CLAIM and s_lo>CLAIM-F(1,10**15)
 assert s_hi<score_bounds(WORST)[0]
 assert F(1,(1<<34)-108)<=F(1,(1<<30)-16)
 setup=sum(rows[k] for k in ['group_fixed','hash_head_zeroing','node_init','gen_points_headers','scalar_sampling','range_tree_driver','norm_sort_search','phase_control','barriers','phase4_stages','xsum_corrections'])
 assert setup<=1000*(1<<119)
 slotbytes=32*CAP_PG*SIDE;caph=a['caph']
 fixed=(1<<51)+2*(1<<39)+2*(1<<41)+hi*96+Q2*96+(1<<34)+(1<<44)+(8192*SIDE+32)+(1<<38)+8*SIDE+64*SIDE+4096+4096
 assert fixed<(1<<52)
 memory_end=(1<<52)+64+slotbytes+40*caph+(8<<HB)  # slots start at 2^52+32 (32 mod 64)
 assert memory_end<(1<<64)
 assert (1<<52)+32+slotbytes+32*(2*max(qs)*Q2)<(1<<64)
 assert CAP_PG%128==0
 assert fixed+(1<<30)+2048+192*Q2+32*hi+(1<<40)<(1<<52)
 print(json.dumps({'research_only':True,'q1':hi,'q2':Q2,'groups':G,'group_sizes':[qs[0],qs[-1]],'cap_pg':CAP_PG,'hash_bits':HB,'carry_cap':a['carry_cap'],
       'success_formula':str(dec(a['success']))[:24],'hop1_group_drop':str(dec(a['dp']))[:12],'home_drop':str(dec(a['dh']))[:12],'intruder_rate':str(dec(a['intr']))[:12],
       'hot_per_leaf_fixed':str(hot),'phi':[str(dec(PHI[g]))[:10] for g in (0,1,G//2,G-1)],'E_factor_groups':str(dec(sum(PHI[g]*qs[g] for g in range(G))/q1))[:12],
       'walk_instructions_per_leaf_E':str(dec(walk_main/S))[:10],'E_total_instructions':str(E),'E_log2_units':str(log(F(E,1000)))[:40],
       'phi1_log2_units':str(log(F(WORST,1000)))[:20],'run_worst_log2_units_ceil4':str(dec(F(RUN_WORST_LOG,10**4))),'claim':str(CLAIM),'score_enclosure':[str(dec(s_lo))[:40],str(dec(s_hi))[:40]],
       'tails':{'worst_per_incidence':TAILW,'charged_mean':TAILE,'arith_min':str(ARITHM)},'rows':{k:str(v) for k,v in rows.items()},'layout_checks':True,'memory_end':memory_end,'cap_h':a['caph'],
       'scope':'Grouped E-ledger; coefficients in proof.md.'}))
