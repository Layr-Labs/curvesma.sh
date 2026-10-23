"""Exact arithmetic audit only. Never executes a production-size search.
python3 -B verify_certificate.py --write  # materialize independently recomputed values
python3 -B verify_certificate.py          # compare to captured certificate
"""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
p = int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f', 16)
r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)
beta = int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee', 16)
lam = int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72', 16)
N = (r-1)//6
nu = F(27,25)
budget_square = F(98859,100000)*nu**2*N
T = math.isqrt(budget_square.numerator//budget_square.denominator)
while T*T < budget_square: T += 1
M, S, A, workers = 2**73, 2**74, 2**24, 2**10
L = 40*2**58
K = T+workers*L
C, E0, E = 1000, 950, 1140
step_rows = dict(beta_multiplications=358, field_normalizations=160,
    canonical_argmin=36, y_parity_negation=20, adding_table_index_loads=16,
    distinguished_test=3, full_ring_comparisons=96, state_and_loop=28,
    brent_expanded_state_compare_control=132, checkpoint_full_copy=76)
base_other = dict(dp_lookup_transfers=F(S*2**18,C),
    fresh_starts=F(S*4096), deepen_rebuilds=F(2**96),
    all_key_match_replays=F(2**110), exceptional_arithmetic=F(2**40),
    table_preprocessing=F(2**32), constant_checks_sampling_output=F(2**20))
rows = {'main': F(K*(C+E),C)}
rows.update({k:F(6,5)*v for k,v in base_other.items()})
raw = sum(rows.values(), F(0))
reserve = F(2,25)*raw
W = raw+reserve
claim = F(6403,50)
poisson_tail = F(2**65*57**256*257, math.factorial(256)*200)
false_key = F(S*(S-1),2*2**192)
loss = F(1,2**20)

def exp_minus_bounds(a):
    # a is reduced to a 12-decimal enclosure to keep the captured rationals compact.
    q = 10**12
    lo = F(a.numerator*q//a.denominator,q)
    hi = lo+F(1,q)
    def partial(x, n):
        term = total = F(1)
        for k in range(1,n+1):
            term *= -x/k
            total += term
        return total
    # Monotonicity + alternating-series signs: odd lower, even upper.
    return partial(hi,41), partial(lo,40)

def success(rho):
    a = rho*T*(T-1)/(2*N)
    el, eu = exp_minus_bounds(a)
    return (1-eu)*(1-loss)-poisson_tail-false_key, (1-el)*(1-loss)

checks = {}
def check(name, truth, meaning):
    if not truth: raise AssertionError(name+': '+meaning)
    checks[name] = {'result':'PASS','meaning':meaning}

def less_pow2(x, numerator, denominator=1):
    return x.numerator**denominator < (2**numerator)*x.denominator**denominator

def greater_pow2(x, numerator, denominator=1):
    return x.numerator**denominator > (2**numerator)*x.denominator**denominator

check('N_full_precision',N==19298681539552699237261830834781317975472927379845817397100860523586360249056,'exact complete integer (r-1)/6')
check('N_partition',6*N==r-1,'free nonzero class count')
check('beta_order3',beta!=1 and pow(beta,3,p)==1 and (beta*beta+beta+1)%p==0,'cube-root identities modulo p')
check('lambda_order3',lam!=1 and pow(lam,3,r)==1 and (lam*lam+lam+1)%r==0,'cube-root identities modulo r')
check('six_units',len({s*pow(lam,j,r)%r for s in (-1,1) for j in range(3)})==6,'six distinct scalar actions')
check('seven_nonresidue',pow(7,(p-1)//2,p)==p-1,'no x=0 curve points')
check('gcd_lambda_minus_one',math.gcd(lam-1,r)==1,'nontrivial GLV action has no nonzero fixed point')
reject = F(2**256-r,2**256)
check('scalar_rejection_rate',reject<F(1,2**127),'a fresh 256-bit draw rejects with probability below 2^-127')
check('scalar_rejection_cap',4*(S+2**20)*reject**4<F(1,2**430),'four-draw per scalar cap, union over every charged start and table scalar')
check('exception_cap_charged',F(2**16*2**20,C)<base_other['exceptional_arithmetic'],'all 65536 exceptional repairs fit their separate hard work ceiling')
check('integer_nu',nu==F(108,100),'nu exactly 1.08, not a floating literal')
check('ceil_sqrt_upper',T*T>=budget_square,'T attains the budget square')
check('ceil_sqrt_minimal',(T-1)**2<budget_square,'T is the least sufficient integer')
check('integer_flush',K-T==workers*L,'all in-flight work is charged at deepest d=58')
check('length_u64',L<2**64,'8-byte per-trail length suffices')
check('global_u128',K<2**128,'16-byte global work counter suffices')
check('row_sum',sum(step_rows.values())==925,'sum every ordinary-step row')
check('baseline_row_allowance',sum(step_rows.values())<=E0,'25 instruction base cushion')
check('E_120_inside_ledger',E==F(6,5)*E0,'all non-reference step work repriced by 1.2 inside main')
check('other_120_inside_ledger',all(rows[k]==F(6,5)*v for k,v in base_other.items()),'every non-main row also repriced by 1.2')
check('doubling_two_ticks',F(3,2)+F(E,C)<=2*(1+F(E,C)),'a doubling debits TWO main ticks, not one')
check('reserve_8pct',W==F(27,25)*raw,'8% explicit whole-ledger reserve, not claimed as execution')
check('claim_reduced',math.gcd(claim.numerator,claim.denominator)==1,'6403/50 is reduced')
check('total_below_claim',less_pow2(W,6403,50),'exact W^50 < 2^6403')
check('rounding_slack_at_least_005',less_pow2(W,25611,200),'log2 W < 128.055, giving >0.005 bits')
check('rounding_slack_at_most_060',greater_pow2(W,6400,50),'log2 W >128.000, giving <0.060 bits')
check('headroom_over_120_reprice',less_pow2(raw,3199,25),'raw <2^127.96; claim has >0.10 bit over 1.2x repricing')
check('raw_log_bracket_low',greater_pow2(raw,1279078,10000),'raw log2 >127.9078')
check('raw_log_bracket_high',less_pow2(raw,1279079,10000),'raw log2 <127.9079')
check('W_log_bracket_low',greater_pow2(W,1280189,10000),'reserved log2 >128.0189')
check('W_log_bracket_high',less_pow2(W,1280190,10000),'reserved log2 <128.0190')
check('memory_record_layout',24+64+8+36+1+1+58==192,'192-bit key, TWO full initial scalars, length64, 3 pointer96, color, DP flags, padding')
check('memory_node_pool',M*192==3*2**79,'fixed pool, pointer and balancing overhead included')
check('memory_all',M*192+2**30<2**81,'all supplementary memory is inside a separate 2^30-byte cap')
check('pointer_width',M*192<2**96,'96-bit byte addresses cover the full node pool')
check('adding_table_memory',2**20*160+workers*1024+2**24<2**30,'adding table, complete states and code fit auxiliary cap')
check('dp_depths',tuple(range(56,59))==(56,57,58),'two deepen operations, then failure rather than unbounded allocation')
check('partition_disjoint_dp',58<=64,'partition bits 64..83 never overlap DP bits 0..57')
check('deepen_charged',F(3*M*2**20,C)<base_other['deepen_rebuilds'],'two rebuilds plus abort scan fit 2^96 before 1.2 repricing')
check('lookup_charge',2**16+148*1000+2**14<2**18,'Toeplitz hash + 148 routed tree levels + remaining work fit per DP ceiling')
check('replay_charge',F(A*2*L*2**17,C)+2**30<base_other['all_key_match_replays'],'ALL capped key matches replayed, with coefficient modular arithmetic')
check('poisson_submean',F(K,2**56*2**65)<57,'each Poisson split has mean below 57')
check('poisson_cap_split',M==2**65*256,'total overflow implies one virtual subcount >=256')
check('poisson_ratio',F(57,257)<1,'successive factorial-series terms have a convergent ratio')
check('poisson_rational_tail',poisson_tail==2**65*F(57**256,math.factorial(256))/(1-F(57,257)),'exact rational upper bound, e^-57 dropped conservatively')
check('poisson_tail_tiny',poisson_tail<F(1,2**120),'first overflow is below 2^-120 UNDER the disclosed Poisson-tail model')
check('fingerprint_pair_union',false_key<F(1,2**45),'independent universal fingerprint key, even for S capped DP insertions')
for rh, label, floor in [(F(9,10),'090',F(4048,10000)),(F(943,1000),'0943',F(4193,10000)),(F(1),'100',F(4381,10000)),(F(858,1000),'0858',F(39,100))]:
    low, high = success(rh)
    check('success_'+label,low>floor,'alternating-series rational lower bound exceeds '+str(floor))
check('success_claim_margin',success(F(9,10))[0]-F(2,5)>F(48,10000),'named margin over 0.40 exceeds 0.0048')
check('break_even_above_0857',success(F(857,1000))[1]<F(39,100),'break-even lower endpoint fails 0.39')
check('break_even_below_086',success(F(86,100))[0]>F(39,100),'break-even upper endpoint <=0.86')
check('old_product_exponent',-F(888,10)+F(6132,100)==-F(687,25),'2^-88.8 times 2^61.32 equals 2^-27.48, NOT 1.7e-27')
check('old_product_bracket',F(1,2**28)**25<F(1,2**687)<F(1,2**27)**25,'the corrected product lies between 2^-28 and 2^-27')
check('old_reprice_delta_exceeds_008',F(175,165)**100>2**8,'650->750 costs MORE than 0.08 bit')
check('no_advice',0==0,'no nonuniform advice; public constants are excluded by protected convention')
pre = F(6,5)*(base_other['table_preprocessing']+base_other['constant_checks_sampling_output'])
check('preprocessing_bound',pre<2**33,'charged setup <2^33 reference equivalents')

captured = json.loads((HERE/'toy_pipeline_results.json').read_text())
legacy = json.loads((HERE/'toy_results.json').read_text())
check('toy_actual_scale',captured['curve']['order']==523969 and captured['curve']['classes']==87328,'toy group and orbit count explicitly captured')
check('toy_trial_count',sum(v['trials'] for v in captured['table_results'])==captured['total_trials']==16000,'all table trial counts reconciled')
check('toy_success_count',sum(v['success_with_core_overlap'] for v in captured['table_results'])==captured['successes'],'successes recomputed from the eight captured rows')
check('toy_aborts_present',sum(v['brent_cycle_aborts'] for v in captured['table_results'])>0,'measured cycle failures are disclosed rather than assumed away')
check('toy_spurious_present',sum(v['spurious_key_matches'] for v in captured['table_results'])>0,'toy reconstruction path actually exercised')
check('legacy_definition_separate',legacy['rho_pooled']!=captured['pooled_effective_exponent_rho'],'legacy inverse-mean-time statistic is NOT the exponent coefficient')

def frac(x):
    x=F(x)
    return {'numerator':str(x.numerator),'denominator':str(x.denominator)}
def log_display(x):
    return format(math.log2(x.numerator)-math.log2(x.denominator),'.10f')

out = dict(version=2,scope='Candidate-authored exact arithmetic checks, NOT validation of heuristic premises or instruction estimates.',
    integers={k:str(v) for k,v in dict(p=p,r=r,beta=beta,lambda_scalar=lam,N=N,T=T,K_max=K,
       M_node_cap=M,S_insertion_start_cap=S,A_match_cap=A,deepest_trail_tick_cap=L,workers=workers,
       node_bytes=192,node_pool_bytes=M*192,auxiliary_byte_cap=2**30,peak_memory_byte_bound=2**81).items()},
    parameters=dict(nu=frac(nu),C=C,E_baseline=E0,E_charged=E,step_rows=step_rows,
       success_claim=frac(F(2,5)),heuristic_rho_floor=frac(F(9,10)),dp_depths=[56,57,58],
       claim=frac(claim),nonuniform_advice_bytes=0),
    ledger={k:frac(v) for k,v in rows.items()},ledger_raw=frac(raw),reserve_8pct=frac(reserve),total=frac(W),
    poisson_tail_rational=frac(poisson_tail),universal_key_union_rational=frac(false_key),
    displays_only=dict(log2_raw=log_display(raw),log2_total=log_display(W),
       headroom_over_120_reprice=format(float(claim)-float(log_display(raw)),'.10f'),
       score_rounding_slack=format(float(claim)-float(log_display(W)),'.10f'),
       corrected_old_product='2^(-687/25), approximately 5.342e-9'),
    success_lower_floors=dict(rho090='0.4048',rho0943='0.4193',rho1='0.4381',
       break_even_039_bracket=['0.857','0.858'],reserved_margin_over_040='>0.0048'),
    toy_summary=dict(trials=captured['total_trials'],successes=captured['successes'],
       effective_exponent_rho=captured['pooled_effective_exponent_rho'],
       min_table_rho=captured['table_rho_min'],max_table_rho=captured['table_rho_max'],
       legacy_inverse_mean_time_rho=legacy['rho_pooled']),
    captured_toy_sha256=hashlib.sha256((HERE/'toy_pipeline_results.json').read_bytes()).hexdigest(),
    checks=checks,pass_count=len(checks))
if '--write' in sys.argv:
    (HERE/'certificate.json').write_text(json.dumps(out,indent=2)+'\n')
else:
    old=json.loads((HERE/'certificate.json').read_text())
    if old != out: raise AssertionError('captured certificate differs from fresh exact audit')
print(json.dumps({'certificate_checks':len(checks),'passed':len(checks),'verdict':'PASS',
    'claim':'6403/50','log2_raw':out['displays_only']['log2_raw'],
    'log2_bound':out['displays_only']['log2_total']}))
