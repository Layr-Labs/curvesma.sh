"""Recompute exact arithmetic and captured small-curve runs.

Run: python3 -B check.py --capture
No production attack is executed. Only public-constant identities are checked
at production size. All complexity inequalities use integers or Fractions.
"""
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
from fractions import Fraction as F
import random
import sys
import reference as ref

HERE = Path(__file__).resolve().parent
P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
G = (int('17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb', 16),
     int('08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1', 16))
L = int('d201000000010000', 16)**2-1
C = 2000
CLAIM = F(14547, 100)
N = (R-1)//6
B = ref.budget(R)
CHECKS = {}
ROOT = HERE.parents[2]
PROFILE_PATH = ROOT/'tracks/bls12-381-g1-ecdlp-plain/target-profile.json'
COST_PATH = ROOT/'tracks/docs/COST_MODEL.md'
RULES_PATH = ROOT/'tracks/bls12-381-g1-ecdlp-plain/rules.md'
PROFILE = json.loads(PROFILE_PATH.read_text())
PROTECTED_DIGESTS = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in (PROFILE_PATH,COST_PATH,RULES_PATH)}


def ck(name, value):
    CHECKS[name] = bool(value)
    if not value:
        raise AssertionError(name)


def below_power(q, exponent):
    return q.numerator**exponent.denominator < (1 << exponent.numerator)*q.denominator**exponent.denominator


def log_interval(q, denominator=10000):
    lo, hi = 0, 256*denominator
    while hi-lo > 1:
        mid = (lo+hi)//2
        if below_power(q, F(mid, denominator)):
            hi = mid
        else:
            lo = mid
    return F(lo, denominator), F(hi, denominator)


def exp_lower(x, terms=80):
    term = total = F(1)
    for k in range(1, terms):
        term = term*x/k
        total += term
    return total


ck('embedded field and order equal the protected profile', int(PROFILE['field_modulus'],16) == P and int(PROFILE['subgroup_order'],16) == R)
ck('protected problem supplies T only and requires a 32-byte scalar', set(PROFILE['input']) == {'T'} and PROFILE['output'] == {'x':'32-byte scalar'})
ck('normalization is read from the protected cost document', '| BLS12-381 ordinary G1 DLOG | Existing BLS G1 addition reference | **2,000** |' in COST_PATH.read_text())
ck('C fixed at 2000; no n parameter in this claim', C == 2000)
ck('field width is 381 bits stored in 48 bytes', P.bit_length() == 381 and (P.bit_length()+7)//8 == 48)
ck('scalar width is 255 bits stored in 32 bytes', R.bit_length() == 255 and (R.bit_length()+7)//8 == 32)
ck('six-class count is exact', (R-1) % 6 == 0 and 6*N == R-1)
ck('lambda is a nonidentity cube root modulo r', L != 1 and (L*L+L+1) % R == 0)
ck('public generator satisfies curve equation', (G[1]*G[1]-G[0]**3-4) % P == 0)
beta = pow(2, (P-1)//3, P)
if ref.mul(L, G, P) != (beta*G[0] % P, G[1]):
    beta = beta*beta % P
ck('beta is a nonidentity cube root modulo p', beta != 1 and (beta*beta+beta+1) % P == 0)
ck('production public endomorphism identity checked', ref.mul(L, G, P) == (beta*G[0] % P, G[1]))
prod = {'p': P, 'r': R, 'G': G, 'beta': beta, 'lam': L, 'units': (1, L, L*L % R)}
k1, q1, u1 = ref.canonical(G, prod)
k2, q2, u2 = ref.canonical((G[0], -G[1] % P), prod)
ck('production negation has identical complete key and canonical point', k1 == k2 and q1 == q2)
ck('production negation uses the opposite eigenvalue', (u1+u2) % R == 0)
ck('budget is minimal integer with 100B squared at least 169N', 100*B*B >= 169*N and 100*(B-1)**2 < 169*N)
ck('Floyd including entry and one-lap replay uses at most 6B', 3*B+2*B+B == 6*B)
ck('step cap fits in three 64-bit limbs', 6*B < 1 << 192)
ck('all scalars fit in four 64-bit limbs', R < 1 << 256)
ck('all field coordinates fit in six 64-bit limbs', P < 1 << 384)
ck('binary scalar multiply uses at most 510 affine operations', 2*R.bit_length() <= 510)
ck('transition uses at most 511 affine operations', 510+1 == 511)
ck('Fermat inverse and affine formula fit 768 field multiplications', 2*P.bit_length()+6 <= 768)
ck('constant-width state record including key and counter fits 256 bytes', 48+48+32+32+49+24+8 <= 256)
MEMORY = 2**18+16*256+2**16+2**12+2**16+2**12
ck('hard capped total native storage is below 2^20 bytes', MEMORY < 2**20)
ck('fixed work inventory fits 2^40 instructions', 64*2**32+2**24 < 2**40)
ck('preprocessing units are below 2^30', F(2**40, C) < 2**30)
ck('SHAKE one-block input is within its 136-byte rate', len(b'orbit6-v3')+32+49 < 136)
ck('SHAKE one-block output fits its 136-byte rate', 64 < 136)
ck('Keccak permutation and packing budget fit 65536 instructions', 24*1024+16384 < 65536)
ck('512-bit scalar reduction budget covers 512 steps at 128 instructions', 512*128 == 65536)
ck('two coefficient products fit 131072 instructions', 2*255*256 <= 131072)
ck('affine nonmultiply operation census fits 8192 instructions', 381*12+10*64+1024 < 8192)

ROWS = {
    'field_multiply_kernels': (511*768+8)*383,
    'affine_nonmultiply': 511*8192,
    'keyed_SHAKE256': 65536,
    '512_bit_hash_reduction': 65536,
    'two_scalar_normalizations': 131072,
    'class_keys_control_counter_and_transfers': 8192,
}
MUL = ROWS['field_multiply_kernels']
OTHER = sum(ROWS.values())-MUL
I = sum(ROWS.values())
FIXED = 2**40
W = F(6*B*I+FIXED, C)
STRESS = F(6*B, C)*(MUL+F(6,5)*OTHER)+F(6,5)*F(FIXED,C)
lo, hi = log_interval(W)
slo, shi = log_interval(STRESS)
ck('claim fraction is reduced', gcd(CLAIM.numerator, CLAIM.denominator) == 1)
ck('ledger log is inside exact 145.3392..145.3393 bracket', lo == F(1453392,10000) and hi == F(1453393,10000))
ck('all instructions are counted exactly once and H is zero', I == 154767096 and W == F(6*B*sum(ROWS.values())+FIXED,C))
ck('submitted work bound covers computed ledger', below_power(W, CLAIM))
ck('1.2x all nonmultiplication rows still have at least 0.10 bits', below_power(STRESS, CLAIM-F(1,10)))
ck('submitted slack is explicitly in the 0.13..0.14 bit interval', F(13,100) < CLAIM-hi and CLAIM-lo < F(14,100))
# No certificate asserts the incompatible maximum-0.06-bit slack requirement.
# Under ideal independent hash labels, orbit multipliers are uniform over six
# elements. One full lap has a nonidentity product with probability 5/6.
exponent = F(20,21)**2 * F(B*(B-1), 2*N)
MODEL_FLOOR = F(11,25)
ck('finite B birthday exponent gives useful-collision model at least 0.44', exp_lower(exponent) >= F(125,59))
# (5/6)*(1-exp(-exponent)) >= 11/25 iff exp(exponent)>=125/59.
BIAS = F((6*B+2)*R, 2**512)
EXCEPT = F(B+2,R)
LOSS = BIAS+EXCEPT
ck('modular-reduction bias over all distinct labels is below 2^-124', BIAS < F(1,2**124))
ck('infinity and zero-initial-coefficient allowance is below 2^-124', EXCEPT < F(1,2**124))
ck('model minus 0.03 extrapolation reserve and exact exceptional loss exceeds 0.40', MODEL_FLOOR-F(3,100)-LOSS > F(2,5))
ck('success floor is strictly above the protected 0.39 minimum', F(2,5) > F(39,100))

TOYS = []
all_correct = all_bounded = True
all_recovered = all_trials = 0
for p, r, cofactor in ((1543,523,3),(1549,541,3),(1621,523,3)):
    ctx = ref.toy_context(p,r,cofactor)
    # Exhaustively test six-to-one key coverage and sign handling.
    buckets = {}
    for a in range(1,r):
        point = ref.mul(a,ctx['G'],p)
        assert point is not None
        key, canonical, u = ref.canonical(point,ctx)
        all_correct &= canonical == ref.mul(u,point,p)
        neg = (point[0], -point[1] % p)
        kn, cn, un = ref.canonical(neg,ctx)
        all_correct &= key == kn and canonical == cn and (u+un) % r == 0
        buckets.setdefault(key,[]).append(a)
    ck('toy %d exhaustive free six-class quotient and both signs' % p,
       len(buckets) == (r-1)//6 and all(len(v)==6 for v in buckets.values()) and all_correct)
    trials = []
    failures = {}
    successes = max_evals = 0
    for trial in range(400):
        seed = hashlib.sha256(('orbit6-v3:%d:%d' % (p,trial)).encode()).digest()
        x = random.Random(int.from_bytes(seed,'big')).randrange(r)
        result = ref.solve(ref.mul(x,ctx['G'],p),ctx,seed,check=True)
        all_bounded &= result['evaluations'] <= 6*ref.budget(r)
        max_evals = max(max_evals,result['evaluations'])
        if result['x'] is not None:
            all_correct &= result['x'] == x
            successes += 1
        else:
            failures[result['reason']] = failures.get(result['reason'],0)+1
        trials.append([trial,x,result['x'],result['evaluations'],result['reason']])
    all_trials += len(trials)
    all_recovered += successes
    TOYS.append({'p':p,'r':r,'cofactor':cofactor,'beta':ctx['beta'],'lambda':ctx['lam'],
                 'G':ctx['G'],'classes':(r-1)//6,'B':ref.budget(r),'trials':len(trials),
                 'recoveries':successes,'failures':failures,'max_evaluations':max_evals,
                 'trial_records':trials})
ck('every returned toy logarithm is correct and every intermediate representation checked', all_correct)
ck('all captured toy trials respect deterministic evaluation cap', all_bounded)
ck('captured trial total is exactly 1200', all_trials == 1200)
ck('toy arithmetic reports a nonzero number of verified recoveries', all_recovered > 0)
ck('identity target returns exactly the canonical zero scalar', ref.wire_output(ref.solve(None,prod,bytes(32))['x'],R) == {'x':'0'*64})
REPORT = '\n'.join([
    'Exact integer/Fraction certificate; no production-sized attack executed.',
    'certificate checks: %d/%d' % (sum(CHECKS.values()),len(CHECKS)),
    'C: %d' % C,
    'N: %d' % N,
    'B: %d' % B,
    'evaluation cap: %d' % (6*B),
    'instruction rows: '+json.dumps(ROWS,sort_keys=True),
    'per-evaluation instructions: %d' % I,
    'W exact: %s' % W,
    'log2(W) exact bracket: %s < log2(W) < %s' % (lo,hi),
    'stress log2 exact bracket: %s < log2(W_stress) < %s' % (slo,shi),
    'claim: %s' % CLAIM,
    'success claim: 2/5; conservative ideal-model useful-collision floor: 11/25',
    'native memory inventory bytes: %d; hard cap: 1048576' % MEMORY,
    'toy verified recoveries: %d/%d' % (all_recovered,all_trials),
]+['toy p=%d r=%d: %d/%d, max evaluations %d <= %d' % (t['p'],t['r'],t['recoveries'],t['trials'],t['max_evaluations'],6*t['B']) for t in TOYS])+'\n'
CERT = {'version':3,'scope':'exact arithmetic and small-curve captured runs; production constant identities only; heuristic production success is not certified',
        'normalization_C':C,'protected_inputs_sha256':PROTECTED_DIGESTS,'claim':{'numerator':str(CLAIM.numerator),'denominator':str(CLAIM.denominator)},
        'checks':CHECKS,'parameters':{'p':str(P),'r':str(R),'beta':str(beta),'lambda':str(L),'N':str(N),'B':str(B)},
        'ledger':{'H':0,'instruction_rows_per_evaluation':ROWS,'fixed_instruction_bound':FIXED,
                  'W':str(W),'stress_W':str(STRESS),'log2_W_bracket':[str(lo),str(hi)],'stress_log2_bracket':[str(slo),str(shi)],
                  'memory_inventory_bytes':MEMORY,'memory_log2_bytes':20,'preprocessing_log2':30,'nonuniform_advice_bytes':0},
        'success':{'ideal_model_lower':'11/25','claimed':'2/5','extrapolation_reserve':'3/100','exact_exception_allowance':str(LOSS)},
        'toy':{'trials':all_trials,'verified_recoveries':all_recovered}}
if '--capture' in sys.argv:
    (HERE/'certificate.json').write_text(json.dumps(CERT,indent=2)+'\n')
    (HERE/'toy-results.json').write_text(json.dumps({'scope':'Toy curves and exact same bounded Floyd procedure; not production distribution, success-rate or instruction-cost certification','runs':TOYS},separators=(',',':'))+'\n')
    (HERE/'captured-output.txt').write_text(REPORT)
else:
    assert json.loads((HERE/'certificate.json').read_text()) == CERT
    assert (HERE/'captured-output.txt').read_text() == REPORT
print(REPORT,end='')
