"""Exact rational current-cost calculation over a verbatim reference certificate.

Import is silent and does not consume stdin. All claim-bearing calculations use
integers or Fractions; floating approximations only seed searches whose final
bounds are checked with rational intervals. Python caches this module after its
first import. The reference certificate is additionally cached by its own path.

base/ is the preserved predecessor, including its old proof and manifest. Its
checks certify that reference, not the new kernels. The current proof, primitive
checks and memory-map checks must justify the replacement bounds below.
"""
import contextlib
import functools
import importlib.util
import hashlib
import io
import json
import math
import sys
import time
import types
from fractions import Fraction as F
from dataclasses import dataclass, fields
from pathlib import Path

t0=time.time()
BASE_DIR=Path(__file__).resolve().parent/'base'
BASELINE_HASHES={'certificate.py': '9317444f084150d8616bd1e6ea5bcb177ef12c85a7beb8df0efd9106ec5e393a', 'encoding_fuzz.py': '621744e1b692e4b56d240d413d6bc9eb4e2a685440fd9c681e342dcbd4735d1f', 'experiment.py': 'f0565b74740b93ee61866a96aa4705cdcbc5aef5706da45c9a78a58c2276f90c', 'kernels/a64.py': 'fbf39bb8c7a64d2d6f3645f9a596f32d28ea7d64822a22ff80a97be8478422dc', 'kernels/leaf.s': '97a18f3b7322dc2f63aa9b32589c313afc91b8dbd22a3c08da39db763c1cf2c8', 'kernels/route.s': '159cf182bd9117b4fd29bbba5784e3a16b7bb45c63044898313e6e44b7978ebb', 'kernels/secp_fe.s': 'f382d760d310333fe99fbb6ac0710afb3ce22929e49f7c50b23845b880bfe867', 'machine.py': 'aee50573ee6afab76b6e474889777fb08d2e70e4cc9aa89af9b3989eec4f6abc', 'proof.md': '7867bdae79a6fc8a7c472d377646b40694c634e50a043c0df9f116a1334091cd', 'sites.py': 'fcd9efc35261e2d26425b31d422fe5036ba922555bf3a0e4686b808c3c9975f6', 'submission.json': 'e130d1bb752d3a3c757a66d10a5c5f1c10c87d3c57a62da9131572702bddb944', 'success.md': '594a0b3b03ce467ce33191ae81d2246b9787f74ba3b4a1328bbbb6e99a2fb741'}
for _rel,_sha in BASELINE_HASHES.items():
    assert hashlib.sha256((BASE_DIR/_rel).read_bytes()).hexdigest()==_sha, ('reference bytes changed',_rel)
_base_path=BASE_DIR/'certificate.py'
_base_name='_card48_reference_'+hashlib.sha256(str(_base_path).encode()).hexdigest()[:16]
if _base_name not in sys.modules:
    _spec=importlib.util.spec_from_file_location(_base_name,_base_path)
    assert _spec is not None and _spec.loader is not None
    _base=importlib.util.module_from_spec(_spec)
    _transcript=io.StringIO()
    sys.modules[_base_name]=_base
    try:
        with contextlib.redirect_stdout(_transcript):
            _spec.loader.exec_module(_base)
        _base._card48_transcript=_transcript.getvalue()
    except BaseException:
        del sys.modules[_base_name]
        raise
c=BASELINE=sys.modules[_base_name]
BASELINE_TRANSCRIPT=c._card48_transcript
assert not c.FAILED, ('reference certificate failed',tuple(c.FAILED),BASELINE_TRANSCRIPT)
assert c.CLAIM_E==F(123216802973012036,10**15)
assert c.C==1000 and c.THRESH==F(39,100)
baseline_time=time.time()-t0
def cd(a,b):return -(-a//b)
def ceilf(x):
    x=F(x)
    return cd(x.numerator,x.denominator)
def rat(x):
    x=F(x)
    return dict(numerator=str(x.numerator),denominator=str(x.denominator))
def decimal_down(x,n=12):
    scale=10**n;v=F(x)*scale;k=v.numerator//v.denominator
    return ('-' if k<0 else '')+str(abs(k)//scale)+'.'+str(abs(k)%scale).zfill(n)
def decimal_up(x,n=12):
    scale=10**n;k=ceilf(F(x)*scale)
    return ('-' if k<0 else '')+str(abs(k)//scale)+'.'+str(abs(k)%scale).zfill(n)
def grid(x,n=4):return str(c.least_grid_up(F(x),10**n))

TB=2**40
UPPER_BLOCKS=2**55//12288
LOWER_BLOCKS=c.F2-UPPER_BLOCKS
LOWER_H2_BEGIN=32767*TB-LOWER_BLOCKS*12288
RCB_END=LOWER_H2_BEGIN//65536*65536
RELAY_CAP=(RCB_END-256*TB-65536)//64
assert RELAY_CAP==277059750328320
assert RCB_END%65536==0 and RCB_END<=LOWER_H2_BEGIN
@dataclass(frozen=True)
class PrimitiveBounds:
    # Source-count obligations: the current kernel checks must establish them.
    k2b:int=256
    flush2:int=21
    dlvhb:int=512
    dlvhs:int=811
    blockh:int=565
    finalh:int=1153
    new_insert:int=3
    main_visit:int=10
    unweighted_visit:int=12
    relay_save:F=F(129,64)
    home_save:F=F(33,16)
    y_transport:int=129
    offer:int=576
    relay_tables:int=51809096722
    relay_table_init:int=262
    home_table_init:int=518
    program_bytes:int=2**19

INPUTS=PrimitiveBounds()
for _field in fields(INPUTS):
    assert type(getattr(INPUTS,_field.name)) is _field.type, _field.name
    assert getattr(INPUTS,_field.name)>0, _field.name
K2B=INPUTS.k2b
FLUSH2,DLVHB,DLVHS,BLOCKH,FINALH=(INPUTS.flush2,INPUTS.dlvhb,INPUTS.dlvhs,INPUTS.blockh,INPUTS.finalh)
NEW_INSERT,MAIN_VISIT,UNWEIGHTED_VISIT=INPUTS.new_insert,INPUTS.main_visit,INPUTS.unweighted_visit
RELAY_SAVE,HOME_SAVE=INPUTS.relay_save,INPUTS.home_save
Y_TRANSPORT,OFFER=INPUTS.y_transport,INPUTS.offer
TABLES_R=INPUTS.relay_tables
REL_TABLE_INIT,HOME_TABLE_INIT=INPUTS.relay_table_init,INPUTS.home_table_init
PROGRAM_BYTES=INPUTS.program_bytes
assert BLOCKH==FLUSH2+c.EXC_ENTRY+c.DONE_C+DLVHB
assert FINALH==c.FINH_DEST+c.FINH_SENT+c.DONE_C+DLVHS+c.MSG_LINK+(DLVHS-DLVHB)
assert RELAY_SAVE==F(43*6,128) and HOME_SAVE==F(88*6,256)
assert NEW_INSERT<4 and MAIN_VISIT<=UNWEIGHTED_VISIT
assert PROGRAM_BYTES%32==0 and PROGRAM_BYTES>=2**18
SCALE=1<<256

# The only changed success dependency is the relay request capacity. Derive
# its exact change directly, independently of certificate.spx's code.
rl_old=c.RELAY_CAP*9//10
rl_new=RELAY_CAP*9//10
TH_RELAY=RELAY_CAP-18*c.FIXED_LEAVES-rl_new
CN=(F(TH_RELAY)-c.M_VG_MAX)**2
assert TH_RELAY>c.M_VG_MAX
P_MINE_CAPS=(c.P_MINE_HOME_CAP+c.LAM_ORD*(
    c._OTH_H_OUT/F(c.CAND_CAP-2*c.PTS_MINE_CHAIN)
    +2*c._OTH_G_OUT/F(c.REQ_CAP-4*c.PTS_MINE_CHAIN)
    +c._RL_CR_MINE/c.cc
    +c._VAR_R_OUT/(RELAY_CAP-36*c.PTS_MINE_CHAIN-c.cc-c._M_R_OUT)**2))
CK=c.ck+P_MINE_CAPS-c.P_MINE_CAPS
DK=CK-c.ck
def delta_parts(q):
    n=c.LEAVES*q*c.Q2;s=c.NODES*n
    dw=c.RHO*s*c.D_WLT*(2+c.D_M)+c.DWS
    rlcr=18*c.KB_M*(c.FIXED_LEAVES*(s+2*dw)
        +c.LA*c.KB_D1*(c.SPROD_MAX*c.RHO*c.D_WLT*(2+c.D_M)+c.SG_MAX*c.DWS))
    o=c.RHO*c.SPROD_MAX*c.KB_D1*c.KB_M
    vr=648*(o*(c._FDIFF0+2*c.CLUST_BATCH*s*c.KB_M)
             +c.CLUST_BATCH*s*c.RL_SAME_PL)
    drl=rlcr*(F(1,rl_new)-F(1,rl_old))+vr*(1/CN-1/c.cn)
    rows=2*(c.GOOD_NODE_PAIRS*n*n)*c.P1*c.cr*drl
    return drl,rows,DK+rows

# delta(q)=constant+quadratic+cubic, because rlcr andvr are affine inq and
# intensity is quadratic. Use integer Newton coefficients to price all22000
# groups without22000 Fraction gcds. Check the polynomial at disjoint points.
values=[delta_parts(q)[2] for q in range(4)]
df=[]
while values:
    df.append(values[0]);values=[b-a for a,b in zip(values,values[1:])]
DEN=math.lcm(*(v.denominator for v in df))
NUM=[int(v*DEN) for v in df]
def delta_num(q):return sum(a*math.comb(q,k) for k,a in enumerate(NUM))
def delta(q):return F(delta_num(q),DEN)
assert all(delta(q)==delta_parts(q)[2] for q in (4,5,17,10**6,c.Q1GD))
# The same degree3 polynomial in power basis, evaluated with three integer
# multiplies. Division occurs only once per weight; no Fraction gcd per group.
POLY_NUM=[6*NUM[0],6*NUM[1]-3*NUM[2]+2*NUM[3],3*NUM[2]-3*NUM[3],NUM[3]]
POLY_DEN=6*DEN
_poly_gcd=math.gcd(POLY_DEN,*POLY_NUM)
POLY_NUM=[v//_poly_gcd for v in POLY_NUM]
POLY_DEN//=_poly_gcd
def delta_poly_num(q):
    a,b,d,e=POLY_NUM
    return ((e*q+d)*q+b)*q+a
assert all(F(delta_poly_num(q),POLY_DEN)==delta(q)
    for q in (0,1,2,3,4,5,17,10**6,c.Q1GD))

# Crosscheck the independently derived delta against the published full exact
# success formula with just the new capacity dependencies supplied.
env=dict(c.__dict__)
env.update(RELAY_CAP=RELAY_CAP,TH_RELAY=TH_RELAY,cn=CN,P_MINE_CAPS=P_MINE_CAPS,ck=CK)
new_spx=types.FunctionType(c.spx.__code__,env,'new_spx')
check_q=(0,1,c.MQ1[1],c.MQ1[c.G_GROUPS//2],c.MQ1[-2],c.Q1GD)
for q in check_q:
    old=c.spx(q);new=new_spx(q)
    assert old['lo']-new['lo']==delta(q)
    assert new['d_rl']-old['d_rl']==delta_parts(q)[0]
    assert new['glob']-old['glob']==DK
    for key in ('mu','cs','d_p2','d_f','lam_int','d_cl','d_rq','g_def'):
        assert old[key]==new[key]
TFULL=new_spx(c.Q1GD)
LO=TFULL['lo']
assert LO==c.LO_RUN-delta(c.Q1GD)
KG=cd(c.THRESH.numerator*LO.denominator*(1<<128),c.THRESH.denominator*LO.numerator)
PG=F(KG,1<<128)
assert PG*LO>=F(39,100) and F(KG-1,1<<128)*LO<F(39,100)

# For allg, oldPHI>=old exact tail and oldPHI<oldtail+2^-256 when unclipped.
# Thus ceil_256(oldPHI+delta), clipped at1, upper-bounds the new tail and is
# at most one256-bit ulp above a fresh round-up. This extra conservative ulp is
# carried through every work row. It avoids recomputing identical exponentials.
PHIN=[]
OLD_EFFECTIVE_NUM=0
for q,old,sg in zip(c.MQ1[:-1],c.PHI,c.SG):
    assert SCALE%old.denominator==0
    oldn=old.numerator*(SCALE//old.denominator)
    OLD_EFFECTIVE_NUM+=oldn*sg
    newn=min(SCALE,oldn+cd(delta_poly_num(q)*SCALE,POLY_DEN))
    assert oldn<=newn<=SCALE
    PHIN.append(newn)
assert PHIN[0]==SCALE and all(a>=b for a,b in zip(PHIN,PHIN[1:]))
for q in check_q[:-1]:
    if q in c.MQ1[:-1]:
        i=c.MQ1.index(q);tail=1-new_spx(q)['lo'];phi=F(PHIN[i],SCALE)
        assert phi>=min(1,tail)
        if phi<1:assert phi-tail< F(2,SCALE)

def weighted_int(values):
    assert len(values)==len(PHIN) and all(isinstance(v,int) for v in values)
    return F(sum(n*v for n,v in zip(PHIN,values)),SCALE)
SUMPHI=F(sum(PHIN),SCALE)
WS=weighted_int(c.SG)
NEWL=[];NEWR=[];MAINS=[]
positive_core=0
ROUTE_SAVE_PER_LEAF=(RELAY_SAVE+HOME_SAVE)/6
for g,(q,s) in enumerate(zip(c.Q1G,c.SG)):
    assert s%K2B==0
    # Exact ceil difference, rather than the floating forecast's linear form.
    yhot=c.NODES*(cd(q*c.L_CHAIN*Y_TRANSPORT,c.M_G)
                  -cd(q*c.L_CHAIN*c.YRING_TOT,c.M_G))
    yr=q*(OFFER*c.NODES*1023//1024+8*4096*(c.NODES+c.N_REL))
    spill=c.NODES*q*c.CHK*64
    prep=c.NODES*cd(q,c.M_G)*c.L_CHAIN*3
    assert s%36==0 and s%c.K2B==0
    delta_leaf=yhot+yr+spill+prep-(s//36)*19-s//6
    assert c.WGL[g].denominator==1
    nl=c.WGL[g].numerator+delta_leaf
    assert type(nl) is int
    NEWL.append(nl)
    blocks=(s//K2B)*BLOCKH-(s//c.K2B)*c.BLOCKH
    tail=c.N_HOME*(FINALH-c.FINALH)
    split=c.N_REL*64
    null=cd(s,K2B)
    partial=c.N_HOME*2*1024
    assert s%ROUTE_SAVE_PER_LEAF.denominator==0
    routing_credit=(s//ROUTE_SAVE_PER_LEAF.denominator)*ROUTE_SAVE_PER_LEAF.numerator
    assert c.WGR[g].denominator==1
    nr=c.WGR[g].numerator+blocks-routing_credit+tail+split+null+partial
    assert type(nr) is int
    NEWR.append(nr)
    gp=c.NODES*c.LEAVES*c.MQ1[g]*c.Q2
    main=s*gp+s*(s-1)//2-c.SAMEX[g]
    _expected_num=((s*4+c.ALLN*c.STOP_CHECK)*c.PAIRM.denominator
                   +4*main*c.PAIRM.numerator)
    assert c.WGG[g].numerator*c.PAIRM.denominator==_expected_num*c.WGG[g].denominator
    MAINS.append(main)
    # Retain the old deterministic core maximum, ignore all reductions and add
    # every positive new operation. This bounds all exceptional execution too.
    positive_core+=yr+spill+prep+tail+split+null+partial
WMAIN=weighted_int(MAINS)
WLEAF=weighted_int(NEWL)
WROUTE=weighted_int(NEWR)
WHOME=NEW_INSERT*WS+c.ALLN*c.STOP_CHECK*SUMPHI+MAIN_VISIT*c.PAIRM*WMAIN
WK2=c.HCHG_VER*c.PAIRK2*WMAIN
WJ=c.PJ_MEAN*WS+F(c.SJOIN,3)*WS
WEX=SUMPHI*c.EXCH_ROW

# Keep twelve-instruction visits everywhere outside the ordinary weighted main
# row. The proven actual visit cost10<=12; all visit population coefficients
# below are nonnegative. This is intentional conservative overcharging.
visit_extra=UNWEIGHTED_VISIT-4
assert UNWEIGHTED_VISIT>=MAIN_VISIT and visit_extra==8
DEF_MINE=ceilf(c._P_MINE_DEF*(c.W_MINE_DET_MAX+visit_extra*c.S_MINE_TOT*c.CAP_H))
E_MINE=(c.E_MINE-c.S_MINE_TOT+visit_extra*c.E_VIS_MINE+DEF_MINE-c.DEF_MINE
    +cd(c.S_MINE_TOT*BLOCKH,K2B)-cd(c.S_MINE_TOT*c.BLOCKH,c.K2B)
    +c.N_HOME*(FINALH-c.FINALH)+cd(c.S_MINE_TOT,K2B)+c.N_HOME*3*1024)
EXC_WORK=ceilf(c.P_EXC*(c.W_DET_MAX+positive_core+visit_extra*c.SG_FULL*c.CAP_H))
SXV=F(UNWEIGHTED_VISIT,4)*(c.SXV-c.EXC_WORK)+EXC_WORK
DEFV=F(UNWEIGHTED_VISIT,4)*c.DEFV
wlt_visits=c.SG_FULL*(c.RHO*c.SG_FULL*c.D_WLT*(2+c.D_M)*c.KB_C132+c.DWS*c.KB(132))
DEFWLT=c.DEFWLT+visit_extra*wlt_visits
assert all(x>=0 for x in (c.SXV-c.EXC_WORK,c.DEFV,wlt_visits,c.E_VIS_MINE,
    c._P_MINE_DEF*c.S_MINE_TOT*c.CAP_H,c.P_EXC*c.SG_FULL*c.CAP_H))

COLD=dict(
    program_delivery_delta=c.ALLN*((PROGRAM_BYTES-2**18)//32)*2*c.MSG,
    yview_tables=c.NODES*34216000,
    sender_processors=(c.NODES//1024)*4096,
    y_handler_vbar=c.NODES*8,
    clone_arena=c.N_HOME*2**23,
    relay_replacement_page_tables=c.N_REL*1040*TABLES_R,
    relay_source_table=c.N_REL*REL_TABLE_INIT,
    home_source_table=c.N_HOME*HOME_TABLE_INIT,
)
COLDSUM=sum(COLD.values())
ENTRY_DELTA=c.NODES*6*c.CH*9
SEED_COLD=c.NODES*15*2**43
SPLIT_COLD=c.N_REL*128
E_PRE=c.E_PRE+(E_MINE-c.E_MINE)+COLDSUM+ENTRY_DELTA+SEED_COLD+SPLIT_COLD
E_SHARED=c.E_SHARED_REST+E_PRE-c.E_PRE
E_RUN=(E_SHARED+c.E_FIX_WORK+c.E_CAND_WORK+c.E_REQ_WORK
    +WLEAF+WROUTE+WHOME+WK2+WJ+WEX+SXV+c.EK2*c.HCHG_VER+c.ELIST*c.HCHG_LS1
    +DEFV+c.DEFK2+c.DEFL+c.DEFJ+DEFWLT+c.SXK2*c.HCHG_VER)
E_TOTAL=c.GATE_INSNS+PG*E_RUN
CLAIM=c.claim15_of(E_TOTAL/c.C)
assert c.P2(CLAIM)[0]>=E_TOTAL/c.C
assert c.P2(CLAIM-F(1,10**15))[1]<E_TOTAL/c.C
TARGET=F(123186802973,10**9)
assert CLAIM<=TARGET

ROWS=[c.ALLN*c.NODE_INIT+COLDSUM,c.NODES*c.GEN_MUL,(2+2*c.G_GROUPS)*c.ALLN*c.MSG,
    E_MINE,c.NODES*(3*cd(c.Q1GD,c.M_G)+1)*c.COIN_SCALAR,c.NODES*c.MX*c.XSETUP,
    c.NODES*6*c.CH*(c.ENTRY_X+9)+c.NODES*c.HEADER,c.NODES*c.XTBLD_ALL+SEED_COLD,
    c.NODES*c.TBLD,c.N_REL*((c.INIT3_BKT+c.TINIT_ENTRY)*c.F2+2)
        +c.NODES*(c.INIT3_BKT*c.F1+2)+SPLIT_COLD,
    c.ZERO_INIT,c.HR_INIT,WLEAF,c.E_FIX_WORK,WJ,WROUTE,WHOME,SXV,
    WK2+c.SXK2*c.HCHG_VER+c.EK2*c.HCHG_VER,c.ELIST*c.HCHG_LS1,
    DEFV+c.DEFK2+c.DEFL+c.DEFJ+DEFWLT,c.E_CAND_WORK,c.E_REQ_WORK,WEX,
    c.NBAR_ALLOW*c.BARRIER,c.ALLN*c.PHASE_CTRL_G,c.OUTPUT]
assert len(ROWS)==27 and sum(ROWS)==E_RUN

# The predecessor's literal declarations are deliberately reference-only.
# They were checked unchanged above. New declarations bind the new cost below.
REFERENCE_DIFFERS=(NEW_INSERT!=c.HCHG_NEW and MAIN_VISIT!=c.VISIT
    and K2B!=c.K2B and KG!=c.K_GATE and CLAIM!=c.CLAIM_E)
assert REFERENCE_DIFFERS
# A stale gate is a concrete failing control: the old probability no longer
# reaches the threshold after the reduced relay cap is charged.
STALE_GATE_REJECTED=c.P_GATE*LO<c.THRESH
assert STALE_GATE_REJECTED
assert all(type(v) in (int,F) for v in ROWS)
assert all(type(v) is F for v in (LO,PG,E_RUN,E_TOTAL,CLAIM))
assert all(v>=0 for v in ROWS)
LO_RUN,P_GATE,GATED_LO,K_GATE=LO,PG,PG*LO,KG
CLAIM_E=CLAIM
LEDGER_ROWS=ROWS
SUCCESS_CHECKS=dict(
    c19b_mine_caps=P_MINE_CAPS<c.Z2(-19,1),c20_success=LO>=F(9985,10000),
    c21_gate=PG*LO>=F(39,100),
    c30_cap_terms=max(TFULL['d_cl'],TFULL['d_rq'],TFULL['d_rl'])<c.Z2(-152,10),
    c31_global=TFULL['glob']<c.Z2(-18,1),
    c32_weights_monotone=all(a>=b for a,b in zip(PHIN,PHIN[1:])))
assert all(SUCCESS_CHECKS.values())
IMPORT_SECONDS=time.time()-t0
@functools.lru_cache(maxsize=1)
def _report():
    # Optional display grids and large success strings are formed on request.
    eff=F(OLD_EFFECTIVE_NUM,6*SCALE)
    target_run=(c.P2(TARGET)[0]*c.C-c.GATE_INSNS)/PG
    budget=target_run-E_RUN
    rowinfo=[dict(row=i+1,log2_upper_4dp=decimal_up(c.least_grid_up(F(v),10**4),4))
             for i,v in enumerate(ROWS)]
    result=dict(
        scope='Exact current-cost arithmetic over a verbatim predecessor; primitive bounds and the current map are checked separately.',
        baseline_certificate_pass=not c.FAILED,baseline_claim=rat(c.CLAIM_E),
        exact_claim_upper=rat(CLAIM),exact_claim_decimal=decimal_up(CLAIM,15),
        target=rat(TARGET),target_margin_bits_lower=decimal_down(TARGET-CLAIM,15),
        success_lower_12dp=decimal_down(LO),gate_upper_15dp=decimal_up(PG,15),
        exact_gate_numerator=str(KG),gate_denominator=str(1<<128),
        exact_success=rat(LO),exact_gated_success=rat(PG*LO),
        relay_cap=RELAY_CAP,relay_threshold=TH_RELAY,
        full_success_terms={k:dict(lower=decimal_down(TFULL[k],18),upper=decimal_up(TFULL[k],18))
                            for k in ('d_cl','d_rq','d_rl','d_p2','rows','glob')},
        mine_caps_upper_18dp=decimal_up(P_MINE_CAPS,18),
        run_log2_upper_6dp=decimal_up(c.least_grid_up(E_RUN,10**6),6),
        run_units_log2_upper_6dp=decimal_up(c.least_grid_up(E_RUN/c.C,10**6),6),
        total_log2_upper_6dp=decimal_up(c.least_grid_up(E_TOTAL,10**6),6),
        preprocessing_log2_upper_4dp=decimal_up(c.least_grid_up(E_PRE/c.C,10**4),4),
        per_leaf_upper_4dp=decimal_up(E_RUN/c.SG_FULL,4),
        mine_log2_upper_6dp=decimal_up(c.least_grid_up(E_MINE,10**6),6),
        exact_def_mine=str(DEF_MINE),exact_exc_work=str(EXC_WORK),
        program_bytes=PROGRAM_BYTES,
        rows=rowinfo,rows_sum_equals_run=(sum(ROWS)==E_RUN),
        exact_positive_budget=budget>0,
        spare_opened_run_insns_per_effective_inc_lower=decimal_down(budget/eff,15),
        source_table_cold=dict(relay=REL_TABLE_INIT,home=HOME_TABLE_INIT),
        cold_charges={k:str(v) for k,v in COLD.items()},
        visit_costs=dict(weighted_main=MAIN_VISIT,unweighted=UNWEIGHTED_VISIT,
            dominance_proved=(UNWEIGHTED_VISIT>=MAIN_VISIT)),
        reference_costs_superseded=REFERENCE_DIFFERS,
        stale_gate_negative_control_rejected=STALE_GATE_REJECTED,
        unchanged_success_check_status=SUCCESS_CHECKS,
        current_map_obligation='The separate current-map checker must validate the split H2 layout, RCB end, YVIEW, clone arena, source tables and code bound.',
        phi_bound='ceil256(old_upper+exact_delta), clipped1; all22000 weights updated; at most one extra ulp over fresh rounding.',
        baseline_seconds=round(baseline_time,3),total_seconds=round(IMPORT_SECONDS,3),
    )
    return result

def summary(include_exact_success=False):
    """Return JSON-compatible values; omit large success rationals by default."""
    result=dict(_report())
    if not include_exact_success:
        result.pop('exact_success',None)
        result.pop('exact_gated_success',None)
    return result


def __getattr__(name):
    # Backwards-compatible access without computing display values on import.
    if name=='REPORT':return _report()
    raise AttributeError(name)
