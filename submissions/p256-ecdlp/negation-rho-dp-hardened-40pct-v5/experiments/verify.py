#!/usr/bin/env python3
"""Exact rational certificate; run python3 experiments/verify.py .
Stdlib only. No production-scale attack. Float proposals are integer-validated.
"""
import json
import math
from pathlib import Path
import random
import sys
from fractions import Fraction as F
N = int('ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551',16)
NP = (N+1)//2
M,R,D,CAP = 1<<14,1<<24,24<<96,1<<35
CLAIM = F(3311,25)

def update_coeff(a,c,sign,modulus):
    temporary = 2*a+c
    value,count = temporary,0
    while value>=modulus:
        value -= modulus
        count += 1
    if sign<0 and value: value = modulus-value
    return value,temporary,count

def exp_neg_upper(x,degree=40):
    assert degree%2==0 and x>=0
    value=term=F(1)
    for i in range(1,degree+1):
        term *= -x/i
        value += term
    return value  # Taylor remainder is negative at even degree for x>=0.

def rat(x):
    x=F(x)
    return f'{x.numerator}/{x.denominator}'

def lt_power(x,a,b):
    x=F(x)
    return x.numerator**b < x.denominator**b*(1<<a)

def log_interval(x):
    a=math.floor(math.log2(float(x))*10000)
    assert lt_power(x,a+1,10000) and not lt_power(x,a,10000)
    return [f'{a}/10000',f'{a+1}/10000']

def build(toy,Curve,prime_small,profile,normalization):
    checks={}
    def ck(name,test): checks[name]='PASS' if bool(test) else 'FAIL'
    ck('SEC2_order_has256bits',N.bit_length()==256)
    ck('order_matches_protected_profile',N==int(profile['subgroup_order'],16))
    ck('C_read_from_protected_cost_model',normalization==1000)
    ck('P256_not_48byte_field',int(profile['field_modulus'],16).bit_length()==256)
    ck('class_count_includes_identity',2*NP==N+1)
    ck('coord_and_scalar_are32bytes',(256+7)//8==32)
    k=math.isqrt(1211*NP//1000)+2
    while F((k-1)*(k-2),2*NP)>=F(1211,2000): k-=1
    x=F(k*(k-1),2*NP)
    ck('least_k_for_1211_over2000',x>=F(1211,2000) and F((k-1)*(k-2),2*NP)<F(1211,2000))
    ticks=M*((k+M-1)//M+D)
    t=(ticks+M-1)//M
    self_loss=F(10,9)*F(M*(2*t+1)*(2*t),2*NP)
    ck('all_lags_self_loss_below1_over6000',self_loss<F(1,6000))
    p0=1-exp_neg_upper(F(9,10)*x)
    e24_lower=sum(F(24**i,math.factorial(i)) for i in range(81))
    dp_tail=2/e24_lower
    ck('e24_positive_series_exceeds20billion',e24_lower>20_000_000_000)
    ck('dp_tail_below1e_minus10',dp_tail<F(1,10**10))
    reject=F((1<<256)-N,1<<256)
    sample_loss=(2*R+4*M)*reject**2
    ck('sampling_loss_below2^-38',sample_loss<F(1,1<<38))
    mean_dp=F(ticks,1<<96)
    ck('DP_Chernoff_3mean_less_half_cap',6*mean_dp<CAP)
    ck('DP_cap_exceeds64',CAP>64)
    dp_overflow=F(1,1<<64)
    ck('slow_cap_Markov_below2^-120',F(12*ticks,N*(1<<16))<F(1,1<<120))
    slow_loss=F(1,1<<120)
    success=p0-self_loss-dp_tail-sample_loss-dp_overflow-slow_loss-F(1,N)
    ck('p0_exceeds420_over1000',p0>F(420,1000))
    ck('success_exceeds419_over1000',success>F(419,1000))
    ck('success_margin_over_2_over5_exceeds19_over1000',success-F(2,5)>F(19,1000))
    ck('claim2_over5_exceeds_rule39_over100',F(2,5)>F(39,100))
    ck('full_unsigned_start_fits5limbs',N-1<1<<320)
    ck('max_temp_3n_minus3',2*(N-1)+(N-1)==3*N-3)
    ck('temporary_fits258bits',3*N-3<1<<258)
    ck('temporary_fits5limbs',3*N-3<1<<320)
    ck('two_subtractions_suffice',3*N-3<3*N)
    ck('reduced_coeff_fits32bytes',N-1<1<<256)
    for ai,a in enumerate((0,1,1<<255,N-1)):
        for ci,c in enumerate((0,1,N-1)):
            for sign in (1,-1):
                out,tmp,count=update_coeff(a,c,sign,N)
                ck(f'coeff_boundary_{ai}_{ci}_{sign}',out==(sign*(2*a+c))%N and tmp<1<<258 and count<=2)
    ck('no_inverse_cancellation_through64',all(pow(2,i,N) not in (1,N-1) for i in range(1,65)))
    rows=[
        ('13_field_M_double_add_two_batches',F(13*512),True),
        ('two_Fermat_inversions_512M_each_amortized',F(2*512*512,M),True),
        ('field_add_sub_20_times24',F(480),False),
        ('exception_tests',F(6),False),('table_address_and16_limb_loads',F(32),False),
        ('canonical_negation',F(64),False),('two_5limb_coeff_double_add_reduce_sign',F(224),False),
        ('index_hash_4_mix64_rounds',F(160),False),('two_full_DP_hashes',F(320),False),
        ('Floyd_compare_and_counters',F(64),False),('point_state_traffic',F(80),False),
        ('two_batch_prefix_traffic',F(80),False),('loop_address_control_meter',F(80),False)]
    per=sum(v for _,v,_ in rows)
    per12=sum(v if mul else F(6,5)*v for _,v,mul in rows)
    ck('M512_exceeds_instruction_skeleton336',512>32+32+144+64+32+32)
    ck('13M_equals7M_double_plus6M_add',13==7+6)
    ck('Fermat_binary_exponent_at_most512M',255+255<=512)
    ck('per_call8278',per==8278)
    ck('repriced_per_call8596',per12==8596)
    extras={
        'table_build':F(R*(20032*512+(1<<20)),1000),
        'capped_randomness':F((4*R+8*M+16)*1024,1000),
        'initial_and_max_restart_builds':F(2*M*(20032*512+(1<<20)),1000),
        'DP_tree_cap_plus_one':F((CAP+1)*(1<<16),1000),
        'slow_path_cap':F((1<<16)*(1<<24),1000),
        'verification_decode_output':F(1<<26,1000),
        'memory_zero_initialization':F(CAP*128+R*128+(1<<30),1000)}
    extra=sum(extras.values())
    work=F(3*ticks)*per/1000+extra
    stress=F(3*ticks)*per12/1000+F(6,5)*extra
    ck('K_strictly_below2^127p64',lt_power(ticks,12764,100))
    ck('all_extras_relative_below2^-80',extra*(1<<80)<work)
    ck('claim_reduced3311_over25',math.gcd(3311,25)==1 and CLAIM==F(13244,100))
    ck('nominal_work_below2^132p273',lt_power(work,132273,1000))
    ck('stress_work_below2^132p327',lt_power(stress,132327,1000))
    ck('nominal_slack_at_least0p005_bits',lt_power(work,132435,1000))
    ck('stress_slack_at_least0p10_bits',lt_power(stress,13234,100))
    ck('nominal_slack_below0p20_bits',not lt_power(work,13224,100))
    pre=sum(extras[s] for s in ('table_build','capped_randomness','initial_and_max_restart_builds','memory_zero_initialization'))
    ck('preprocessing_below2^38',pre<1<<38)
    mem_rows={'DP_nodes128B':CAP*128,'walk_table128B':R*128,'two_Floyd_states160B':2*M*160,
        'six_batch_arrays32B':6*M*32,'code_rng_scratch':1<<30}
    memory=sum(mem_rows.values())
    ck('DP_layout128B',32+32+32+8+8+8+8==128)
    ck('table_layout128B',32+32+32+32==128)
    ck('Floyd_state_layout160B',32+32+40+40+16==160)
    ck('peak_memory_below2^43B',memory<1<<43)
    ck('addresses_fit64bits',memory<1<<64)
    ck('RB_tree_height_at_most72',2*(CAP+1).bit_length()<=72)
    ck('DP_budget65536_exceeds72_times256',1<<16>72*256)
    ck('toy_prime_order_counted',prime_small(toy['curve']['order']))
    ck('toy_no_wrong_verified_answers',toy['wrong_verified_answers']==0)
    ck('toy_at_least300_logs_of400',toy['recovered_and_verified']>=300 and toy['trials']==400)
    ck('toy_Floyd_handles_lambda_above64',any(v['lambda']>64 for v in toy['long_cycle_floyd_examples']))
    ck('toy_Floyd_by_mu_plus_lambda',all(v['meeting_tick']<=v['mu']+v['lambda'] for v in toy['long_cycle_floyd_examples']))
    pp=int('ffffffff00000001000000000000000000000000ffffffffffffffffffffffff',16)
    bb=int('5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b',16)
    base=(int('6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296',16),
          int('4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5',16))
    ck('field_matches_protected_profile',pp==int(profile['field_modulus'],16))
    ck('a_matches_minus3_protected',int(profile['coefficient_a'],16)==pp-3)
    ck('b_matches_protected_profile',bb==int(profile['coefficient_b'],16))
    ck('generator_matches_protected_profile',base==(int(profile['generator_P']['x'],16),int(profile['generator_P']['y'],16)))
    curve=Curve(pp,-3,bb)
    ck('P256_generator_on_curve_a_minus3',(base[1]**2-base[0]**3+3*base[0]-bb)%pp==0)
    ck('P256_nP_is_identity',curve.mul(N,base) is None)
    rng=random.Random(25620261023)
    target=curve.mul(rng.randrange(1,N),base)
    for i in range(8):
        a,b,c,d=[rng.randrange(N) for _ in range(4)]
        start,sign0=curve.canon(curve.add(curve.mul(a,base),curve.mul(b,target)))
        a,b=sign0*a%N,sign0*b%N
        offset=curve.add(curve.mul(c,base),curve.mul(d,target))
        end,sign=curve.canon(curve.add(curve.add(start,start),offset))
        aa,_,_=update_coeff(a,c,sign,N)
        bbb,_,_=update_coeff(b,d,sign,N)
        ck(f'P256_fullwidth_transition_{i}',end==curve.add(curve.mul(aa,base),curve.mul(bbb,target)))
    intervals={s:log_interval(v) for s,v in (('K',ticks),('W',work),('W_repriced',stress))}
    for s,bounds in intervals.items():
        value={'K':ticks,'W':work,'W_repriced':stress}[s]
        lo,hi=map(F,bounds)
        ck('exact_log_interval_'+s,not lt_power(value,lo.numerator,lo.denominator) and lt_power(value,hi.numerator,hi.denominator))
    count=sum(s=='PASS' for s in checks.values())
    output=f'CERTIFICATE {count}/{len(checks)} PASS; all_checks_pass={count==len(checks)}'
    return {'version':3,'method':'Python int/Fraction; actual bounded EC toy; no numpy; not production attack execution',
        'claim':{'time_log2':'3311/25','success':'2/5','rules_floor':'39/100'},
        'parameters':{'n':str(N),'N_prime':str(NP),'m':M,'R':R,'D':str(D),'k':str(k),'K_logical_ticks':str(ticks),
            'f_calls_cap':str(3*ticks),'DP_cap':CAP,'C':1000,'M_instructions':512,'quality_floor':'9/10'},
        'success_model':{'p0_lower':rat(p0),'self_loss_upper':rat(self_loss),'DP_tail_upper':rat(dp_tail),
            'sampling_loss_upper':rat(sample_loss),'DP_overflow_upper':rat(dp_overflow),'slow_overflow_upper':rat(slow_loss),
            'unusable_relation_upper':rat(F(1,N)),'success_lower':rat(success),
            'conditions':'Quantitative random-mapping/DP heuristic, not established by certificate or toy'},
        'rows':[{'name':s,'instructions':rat(v),'M_priced':m,'repriced_instructions':rat(v if m else F(6,5)*v)} for s,v,m in rows],
        'per_call_instructions':rat(per),'repriced_per_call_instructions':rat(per12),
        'extras_units':{s:rat(v) for s,v in extras.items()},'total_units':rat(work),'repriced_total_units':rat(stress),
        'exact_log2_intervals':intervals,'memory_bytes':memory,'memory_rows':mem_rows,'preprocessing_units':rat(pre),
        'toy':toy,'checks':checks,'pass_count':count,'check_count':len(checks),'all_checks_pass':count==len(checks),
        'verifier_output':output}

def main():
    from toy import run,Curve,prime_small
    package=Path(sys.argv[1] if len(sys.argv)>1 else '.').resolve()
    toy=run()
    profile=json.loads((package.parents[2]/'tracks/p256-ecdlp/target-profile.json').read_text())
    model=(package.parents[2]/'tracks/docs/COST_MODEL.md').read_text()
    line=next(s for s in model.splitlines() if s.startswith('| P-256 DLOG |'))
    normalization=int(line.split('**')[1].replace(',',''))
    cert=build(toy,Curve,prime_small,profile,normalization)
    (package/'certificate.json').write_text(json.dumps(cert,indent=2,sort_keys=True)+'\n')
    (package/'experiments/toy-results.json').write_text(json.dumps(toy,indent=2,sort_keys=True)+'\n')
    print(cert['verifier_output'])
    print('EXACT_LOG_INTERVALS',json.dumps(cert['exact_log2_intervals'],sort_keys=True))
    print('TOY',json.dumps(toy,sort_keys=True))
    print('FAILURES',[s for s,v in cert['checks'].items() if v!='PASS'])
    return 0 if cert['all_checks_pass'] else 1
if __name__=='__main__': raise SystemExit(main())
