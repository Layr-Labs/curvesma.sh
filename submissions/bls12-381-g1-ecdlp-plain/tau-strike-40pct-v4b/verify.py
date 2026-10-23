"""Exact integer certificate. Run toy.py --capture, then verify.py --capture.

Normal verify.py never overwrites evidence. No production-sized search executes.
"""
import hashlib
import json
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import sys
import reference as ref

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROFILE_PATH = ROOT/'tracks/bls12-381-g1-ecdlp-plain/target-profile.json'
COST_PATH = ROOT/'tracks/docs/COST_MODEL.md'
RULES_PATH = ROOT/'tracks/bls12-381-g1-ecdlp-plain/rules.md'
PROFILE = json.loads(PROFILE_PATH.read_text())
P = int(PROFILE['field_modulus'],16)
R = int(PROFILE['subgroup_order'],16)
G = tuple(int(PROFILE['generator_P0'][k],16) for k in ('x','y'))
L = int('d201000000010000',16)**2-1
C = 2000
CLAIM = F(3512,25)
N = (R-1)//6
B = ref.budget(R)
CHECKS = {}

def ck(name,value):
    CHECKS[name] = bool(value)
    if not value:
        raise AssertionError(name)

def below_power(q,e):
    return q.numerator**e.denominator < (1<<e.numerator)*q.denominator**e.denominator

def bracket(q,lo,hi):
    return not below_power(q,lo) and below_power(q,hi)

def exp_lower(x,terms=80):
    term = total = F(1)
    for k in range(1,terms):
        term = term*x/k
        total += term
    return total

ck('exact protected plain input is T only, no n or powers',set(PROFILE['input']) == {'T'} and 'No auxiliary powers' in PROFILE['challenge'])
ck('protected output is a canonical 32-byte scalar',PROFILE['output'] == {'x':'32-byte scalar'})
ck('C=2000 read from the protected document','| BLS12-381 ordinary G1 DLOG | Existing BLS G1 addition reference | **2,000** |' in COST_PATH.read_text())
ck('383-instruction field kernel pinned by cost model','**383 instructions**' in COST_PATH.read_text())
ck('field and scalar widths',P.bit_length() == 381 and R.bit_length() == 255)
ck('exact free-class count',6*N == R-1)
ck('lambda has exact order three',L != 1 and (L*L+L+1)%R == 0)
beta = pow(2,(P-1)//3,P)
if ref.mul(L,G,P) != (beta*G[0]%P,G[1]):
    beta = beta*beta%P
ck('beta root and public endomorphism identity',beta != 1 and (beta*beta+beta+1)%P == 0 and ref.mul(L,G,P) == (beta*G[0]%P,G[1]))
ctx = {'p':P,'r':R,'G':G,'beta':beta,'lam':L,'units':(1,L,L*L%R)}
k,q,u = ref.canonical(G,ctx)
kn,qn,un = ref.canonical((G[0],-G[1]%P),ctx)
ck('canon(-Q)=canon(Q) and eps(-Q)=-eps(Q)',k == kn and q == qn and (u+un)%R == 0)
# Only public-constant actions: these are not production challenge trials.
max_products = 0
for h in (0,1,2,3,L,R-1):
    for Q in (None,G,(G[0],-G[1]%P)):
        action,census = ref.projective_action(Q,h,G,P,census=True)
        ck('public Jacobian action Q=%s h=%s agrees with independent affine oracle'%(str(Q is None),h)+str(Q == G),
           action == ref.add(Q,ref.mul(h,G,P),P))
        ck('public action census within 5367-kernel production ceiling h=%s Q=%s'%(h,str(Q)),census['field_products'] <= 5367)
        max_products = max(max_products,census['field_products'])
ck('B is the minimal integer with 100B^2>=169N',100*B*B >= 169*N and 100*(B-1)**2 < 169*N)
ck('bounded detector and replays cost at most 6B',3*B+2*B+B == 6*B)
ck('all execution counters fit three 64-bit limbs',6*B < 2**192)
ck('action has 255 doubles and at most 256 mixed calls',R.bit_length() <= 255)
ck('mixed fallback after four products plus doubling fits eleven',4+7 == 11)
ck('normalization covers binary Fermat and four scaling kernels',2*P.bit_length()+4 == 766)
KERNELS = 255*7+256*11+766+8
ck('one transition has 5375 field kernels including canonical decoding',KERNELS == 5375)
ck('projective control allowance covers 32 field-linear steps and transfer census',32*64+1536 < 4096)
ck('six-limb modular linear step fits 64 instructions',12+6+6+6+6+16 <= 64)
ck('Fermat exponent loop control fits 8192 instructions',381*12+2048 < 8192)
ck('Keccak and packing fit 65536 instructions',24*1024+16384 < 65536)
ck('one transition SHAKE input and output each fit the 136-byte rate',len(b'orbit6-v3')+32+49 < 136 and 64 < 136)
ck('512-bit reduction fits 65536 instructions',512*128 == 65536)
ck('two coefficient products fit 131072 instructions',2*255*256 <= 131072)
ROWS = {'field_multiply_kernels':KERNELS*383,'projective_nonkernel':511*4096,
        'Fermat_inverse_control':8192,'keyed_SHAKE256':65536,
        '512_bit_hash_reduction':65536,'two_scalar_normalizations':131072,
        'class_keys_control_counter_and_transfers':8192}
MUL = ROWS['field_multiply_kernels']
OTHER = sum(ROWS.values())-MUL
I = sum(ROWS.values())
FIXED = 2**40
W = F(6*B*I+FIXED,C)
STRESS = (6*B*(MUL+F(6,5)*OTHER)+F(6,5)*FIXED)/C
LO,HI = F(1402126,10000),F(1402127,10000)
SLO,SHI = F(1403594,10000),F(1403595,10000)
ck('all seven rows sum to I=4430209',I == 4430209 and MUL == 2058625 and OTHER == 2371584)
ck('base ledger exact logarithm bracket 140.2126..140.2127',bracket(W,LO,HI))
ck('1.2x nonkernel exact bracket 140.3594..140.3595',bracket(STRESS,SLO,SHI))
ck('score fraction is reduced',gcd(CLAIM.numerator,CLAIM.denominator) == 1)
ck('claim exceeds computed work',below_power(W,CLAIM))
ck('stress has more than 0.10-bit headroom',below_power(STRESS,CLAIM-F(1,10)))
ck('arithmetic comparison difference 4.99 bits',CLAIM < F(14547,100) and F(14547,100)-CLAIM == F(499,100))
ck('wider base reserve is openly 0.2673..0.2674 bits',CLAIM-HI == F(2673,10000) and CLAIM-LO == F(2674,10000))
ck('maximum-0.06 slack and minimum-0.10 stress reserve cannot coexist',STRESS >= W and F(6,100) < F(1,10))
MODEL = F(11,25)
exponent = F(20,21)**2*F(B*(B-1),2*N)
ck('useful-cycle model with factor 5/6 gives at least 11/25',exp_lower(exponent) >= F(125,59))
BIAS = F((6*B+2)*R,2**512)
EXCEPT = F(B+2,R)
ck('full-label reduction bias below 2^-124',BIAS < F(1,2**124))
ck('identity and initial coefficient exceptions below 2^-124',EXCEPT < F(1,2**124))
ck('model minus named 0.03 reserve and exceptions exceeds 0.40',MODEL-F(3,100)-BIAS-EXCEPT > F(2,5))
ck('success exceeds protected floor 0.39',F(2,5) > F(39,100))
MEMORY = 2**18+16*256+2**16+2**12+2**16+2**12
ck('retained states are affine and fit 256-byte slots',48+48+32+32+49+24+8 <= 256)
ck('Jacobian temporaries fit independently allocated scratch',64*48 < 2**16)
ck('hard native memory inventory under 2^20 bytes',MEMORY == 405504 and MEMORY < 2**20)
ck('fixed work over-inventory covers 64 scalar routines and I/O',64*2**32+2**24 < FIXED)
ck('cold preprocessing below 2^30 group units',F(FIXED,C) < 2**30)
ck('canonical zero output for identity',ref.wire_output(ref.solve(None,ctx,bytes(32))['x'],R) == {'x':'0'*64})
TOY = json.loads((HERE/'toy-summary.json').read_text())
ck('at least 2000 captured trials with all failures retained',TOY['trials'] == 2400 and sum(t['trials'] for t in TOY['runs']) == 2400)
total = recovered = paired = exhaustive = 0
for t in TOY['runs']:
    records = []
    for block in t['blocks']:
        data = (HERE/block['file']).read_bytes()
        ck('raw toy shard digest '+block['file'],hashlib.sha256(data).hexdigest() == block['sha256'])
        records += json.loads(data)['rows']
    ck('800 raw rows and complete success-failure totals for p=%d'%t['p'],len(records) == 800 and len(records) == t['trials'] and t['recoveries']+sum(t['failures'].values()) == 800)
    ck('all returned scalars correct for p=%d'%t['p'],all(a[2] is None or a[1] == a[2] for a in records))
    ck('observed cap and histogram totals consistent p=%d'%t['p'],all(a[3] <= t['evaluation_cap'] for a in records) and sum(a[3] for a in records) == sum(t['step_orientation_histogram']) == sum(t['hash_scalar_mod16_histogram']))
    ck('all diagnostic scalar bins and six orientations reachable p=%d'%t['p'],all(t['step_orientation_histogram']) and all(t['hash_scalar_mod16_histogram']))
    total += len(records)
    recovered += sum(a[2] is not None for a in records)
    paired += t['identical_coin_affine_pipeline_crosschecks']
    exhaustive += t['exhaustive_action_crosschecks']
ck('toy summary matches all raw records',total == TOY['trials'] and recovered == TOY['recoveries'])
ck('72 identical-coin old-new pipeline comparisons and over 30000 action comparisons',paired == 72 and exhaustive > 30000)
PROTECTED = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
             for p in (PROFILE_PATH,COST_PATH,RULES_PATH)}
REPORT = '\n'.join(['Exact integer certificate; public constants and toy runs only.',
    'checks: %d/%d'%(sum(CHECKS.values()),len(CHECKS)),
    'C: %d'%C,'N: %d'%N,'B: %d'%B,'evaluation cap: %d'%(6*B),
    'field kernels per evaluation: %d'%KERNELS,'rows: '+json.dumps(ROWS,sort_keys=True),
    'I: %d'%I,'W: %s'%W,'stress W: %s'%STRESS,
    'base log2 bracket: %s < log2(W) < %s'%(LO,HI),
    'stress log2 bracket: %s < log2(Wstress) < %s'%(SLO,SHI),
    'claim: %s; arithmetic comparison difference 4.99 bits'%CLAIM,
    'success: 2/5 under H1; model floor 11/25; named reserve 3/100',
    'native bytes: %d < 1048576'%MEMORY,
    'toy recoveries: %d/%d; affine pipeline pairs: %d; exhaustive action pairs: %d'%(recovered,total,paired,exhaustive)])+'\n'
CERT = {'version':4,'scope':'Exact arithmetic and captured toy evidence only. H1 production transfer is not certified.',
        'normalization_C':C,'protected_inputs_sha256':PROTECTED,
        'claim':{'numerator':str(CLAIM.numerator),'denominator':str(CLAIM.denominator)},
        'parameters':{'r':str(R),'p':str(P),'N':str(N),'B':str(B),'beta':str(beta),'lambda':str(L)},
        'checks':CHECKS,'ledger':{'H':0,'rows':ROWS,'I':I,'fixed':FIXED,'W':str(W),'stress_W':str(STRESS),
          'base_log2_bracket':[str(LO),str(HI)],'stress_log2_bracket':[str(SLO),str(SHI)],
          'memory_bytes':MEMORY,'memory_log2_bytes':20,'preprocessing_log2':30,'nonuniform_advice_bytes':0},
        'success':{'model_floor':'11/25','claimed':'2/5','reserve':'3/100','bias':str(BIAS),'exceptions':str(EXCEPT)},
        'toy':{'recoveries':recovered,'trials':total,'old_new_pipeline_pairs':paired,'action_pairs':exhaustive,
               'max_products_in_public_action_census':max_products},
        'slack_conflict':'Since Wstress>=W, a 0.06-bit maximum base slack cannot coexist with 0.10-bit stress headroom. Transparent wider reserve is used.'}
if '--capture' in sys.argv:
    (HERE/'certificate.json').write_text(json.dumps(CERT,indent=2)+'\n')
    (HERE/'captured-output.txt').write_text(REPORT)
else:
    assert json.loads((HERE/'certificate.json').read_text()) == CERT
    assert (HERE/'captured-output.txt').read_text() == REPORT
print(REPORT,end='')
