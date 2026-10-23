#!/usr/bin/env python3
"""Exact arithmetic and captured-toy consistency, not a production mixing proof."""
import argparse
from fractions import Fraction as F
from math import gcd,isqrt,prod
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import zlib
from ledger import (P,R,N,NP,T0,T,RESERVE,CLAIM,C,LN_LO,LN_HI,ROWS,
    ORDINARY_CEILING,WEXEC,WSTRESS,W,below_pow2,interval,exp_minus_upper,ancillary)
from base_audit import summary,CONFIGS,DENSITIES
ROOT=Path(__file__).resolve().parents[1]


def frac(x):
    x=F(x);return dict(numerator=str(x.numerator),denominator=str(x.denominator))


def prime(n):return n>1 and (n==2 or (n%2 and all(n%d for d in range(3,isqrt(n)+1,2))))


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true')
    p.add_argument('--bootstrap',action='store_true',help='Arithmetic only; documents remain pending')
    p.add_argument('--raw-dir',type=Path);p.add_argument('--rerun',action='store_true')
    a=p.parse_args()
    if a.rerun:
        assert a.raw_dir,'--rerun requires an external --raw-dir'
        subprocess.run([sys.executable,'-B',str(ROOT/'scripts/faithful_audit.py'),'--write',
            '--raw-dir',str(a.raw_dir)],check=True)
        subprocess.run([sys.executable,'-B',str(ROOT/'scripts/algebra_audit.py'),'--write'],check=True)
    checks={}
    def ck(k,b):checks[k]=bool(b)
    ck('C_20000',C==20000)
    ck('N_exact',12*N==R-1)
    ck('lambda_order12',pow(P,12,R)==1 and pow(P,6,R)==R-1 and pow(P,4,R)!=1)
    ck('cyclotomic',(1+pow(P,4,R)+pow(P,8,R))%R==0)
    ck('Np_ceiling',(NP-1)*1023<N*1024<=NP*1023)
    ck('ln_bracket',F(51082562376599,10**14)<LN_LO<LN_HI<F(51082562376600,10**14))
    ck('ln_width',LN_HI-LN_LO<F(1,1<<300))
    ck('T0_ceiling_sqrt',(T0-1)**2<2*NP*LN_HI<=T0*T0)
    ck('T0_true_ln_same_integer',(T0-1)**2<2*NP*LN_LO)
    ck('bounded_uniform_sampling',R>1<<254 and F(1<<12)*F(3,4)**1024<F(1,1<<400))
    ck('stretch_109_over100',RESERVE==F(109,100) and 100*(T-1)<109*T0<=100*T)
    ck('comparison_target_040',1-exp_minus_upper(LN_LO)>F(399999,10**6))
    ck('named_001_margin',F(2,5)-F(39,100)==F(1,100))
    factors={2:32,3:1,11:1,19:1,10177:1,906349:2,2508409:1,125527:1,
             859267:1,2529403:1,52437899:1,254760293:2}
    ck('r_minus1_factors',prod(q**e for q,e in factors.items())==R-1)
    ck('factor_primes',all(prime(q) for q in factors))
    o=R-1
    for q in factors:
        while o%q==0 and pow(2,o//q,R)==1:o//=q
    d=o//gcd(o,12)
    ck('order2_exact',pow(2,o,R)==1 and all(pow(2,o//q,R)!=1 for q in factors if o%q==0))
    ck('square_segment_no_closure',d>2*T)
    ck('unsquared_positive_counts',T<R)
    ck('denominator_loss_small',F(1,R)<F(1,1<<250))
    ck('Hkey_identity_union',6*T*(1<<120)<N)
    ck('Hkey_tie_union',6*T*1024*(1<<240)<P)
    ck('nonreference_sum_8996',sum(ROWS.values())==8996)
    ck('ceiling_9000',sum(ROWS.values())<=ORDINARY_CEILING==9000)
    ck('two_decodes_per_key',ROWS['first_key_two_Montgomery_decodes_2_times_383']==766
        and ROWS['second_key_two_Montgomery_decodes_2_times_383']==766)
    ck('six_dense_worst_case',WEXEC==F(147,20)*T+(1<<100))
    ck('stress_all_rows',WSTRESS==F(381,50)*T+(1<<101))
    ck('envelope_8T',W==8*T+(1<<100) and WEXEC<W)
    ck('score_reduced',gcd(CLAIM.numerator,CLAIM.denominator)==1)
    ck('score_12881_over100',CLAIM==F(12881,100))
    ck('score_below_arithmetic_threshold',CLAIM<F(12882,100))
    ck('nominal_bracket',not below_pow2(W,1287765,10000) and below_pow2(W,1287766,10000))
    ck('claim_slack_at_least_0005',below_pow2(W,128805,1000))
    ck('claim_slack_at_most_006',not below_pow2(W,128750,1000))
    ck('stress_bracket',not below_pow2(WSTRESS,1287063,10000) and below_pow2(WSTRESS,1287064,10000))
    ck('stress_headroom_010',below_pow2(WSTRESS,128710,1000))
    ar,at=ancillary();M=1<<69;G=129
    fp=F(G*M*(M-1),1<<129)
    ck('ancillary_2pow100',at<1<<100)
    ck('ancillary_stress_2pow101',F(6,5)*at<1<<101)
    ck('record_128_bytes',16+32+32+32+16==128)
    ck('memory_2pow78',128*M+(1<<40)<1<<78)
    ck('fingerprint_128_bits',4608*128//8==73728)
    ck('fingerprint_expectation',fp<1<<17)
    ck('fingerprint_Markov_cap',fp/(1<<32)<F(1,1<<15))
    ck('report_budget',4608*16+2*256*256+72*383+(1<<20)<1<<26)
    ck('report_mean_not_cap',9*T<M*(1<<60))
    ck('DP_bits_in_field',70+128<=381)
    ck('setup_table_bound',1024*(2*510*23000+12*(72*383+10000))+(1<<30)<(1<<50)*20000)
    ck('setup_table_memory',1024*12*(576+64)<1<<24)
    toy=json.loads((ROOT/'toy_results.json').read_text())
    ck('seven_settings',len(toy['configurations'])==7)
    ck('final_fixed_seed',toy['seed']==2026102309 and toy['trials_per_config']==25000)
    ck('final_fixed_stretch',toy['reserve']=={'numerator':109,'denominator':100})
    # Pair difference in [-1,1], n=25000, radius .03: exp(-45/4).
    exp_lower=sum(F(45,4)**j/prod(range(1,j+1)) for j in range(61))
    ck('paired_7_tail_union_lt_001',exp_lower>700)
    # Endpoint radius .015 gives the same exponent; 14 endpoints, two tails.
    ck('endpoint_28_tail_union_lt_001',exp_lower>2800)
    fixed_counts=[];cycle_mass=[];tot=0
    for i,cfg in enumerate(toy['configurations']):
        text=''.join((ROOT/f).read_text() for f in cfg['parts']);c=json.loads(text)
        tag=f'toy{i}_';tr=c['trials'];tot+=tr
        ck(tag+'parameters',(c['q'],c['r'],c['branches'])==CONFIGS[i] and tr==25000 and c['seed']==toy['seed']+i)
        ck(tag+'all_branches_live',c['dead_branches']==0 and min(c['exhaustive_branch_histogram'])>0)
        ck(tag+'census_mass',sum(c['exhaustive_branch_histogram'])==c['classes'])
        ck(tag+'literal_metadata',c['literal_tower_crosschecks']==24)
        cr=c['control'];v=c['variants']['orientation_guard'];ss=v['shared'];lab=v['labels']
        diff=F(v['useful_first_repeat_at_stretch']-cr['repeat_at_base'],tr)
        ck(tag+'paired_CDF_positive',diff-F(3,100)>0)
        ck(tag+'paired_exact',v['paired']['useful_repeat_CDF_difference']==
            {'numerator':diff.numerator,'denominator':diff.denominator})
        ck(tag+'endpoint_lower_gt039',F(v['useful_first_repeat_at_stretch'],tr)-F(3,200)>F(39,100))
        ck(tag+'no_structural_cycles',lab['structurally_fruitless']==0)
        ck(tag+'no_fruitless_two_cycles',lab['structurally_fruitless_two_cycles']==0)
        ck(tag+'label_totals',lab['zero_denominator']+lab['verified_recovery_eventual']==tr)
        ck(tag+'all_variants_reported',set(c['variants'])=={'baseline','inverse_guard','orientation_guard'})
        for name,z in c['variants'].items():
            ztag=tag+name+'_'
            ck(ztag+'literal_steps',z['literal_R_adding_twisted_and_guard_checks']==32)
            ck(ztag+'branch_visits',sum(z['visited_branch_histogram'])==z['transitions'])
            ck(ztag+'DP_censor_mass',all(z['shared']['dp'][str(b)]['hit_before_repeat']+
                z['shared']['dp'][str(b)]['no_hit_before_repeat']==tr for b in DENSITIES))
            ck(ztag+'cycle_hist_mass',all(sum(n for _,n in z['shared'][k]['histogram'])==tr
                for k in ('mu','cycle_lambda','first_repeat','floyd_meeting')))
        if a.raw_dir:
            raw=zlib.decompress((a.raw_dir/f'raw-{i}.json.zlib').read_bytes());rawobj=json.loads(raw)
            ck(tag+'raw_SHA256',hashlib.sha256(raw).hexdigest()==c['raw_sha256'])
            ck(tag+'raw_control_summary',summary(rawobj['control'])==cr)
            for name,z in c['variants'].items():
                rr=rawobj['data'][name];ll=rawobj['labels'][name]
                ck(tag+name+'_raw_summaries',summary(rr)==z['shared'] and len(ll)==tr)
                names=('zero_denominator','structurally_fruitless','verified_recovery_eventual',
                    'verified_recovery_within_cap','two_cycles','structurally_fruitless_two_cycles')
                ck(tag+name+'_raw_labels',z['labels']=={n:sum(row[j] for row in ll) for j,n in enumerate(names)})
                ck(tag+name+'_raw_endpoint',sum(int(row[5] and label[2]) for row,label in zip(rr,ll))==
                    z['useful_first_repeat_at_stretch'])
        fixed_counts.append([c['classes'],c['branches'],cr['repeat_at_base'],
            v['useful_first_repeat_at_stretch'],lab['verified_recovery_within_cap'],lab['zero_denominator']])
        cycle_mass.append([c['branches']]+[c['variants'][k]['labels']['structurally_fruitless']
            for k in ('baseline','inverse_guard','orientation_guard')])
    ck('175000_trial_pairs',tot==175000)
    ck('three_variants_525000_trajectories',3*tot==525000)
    alg=json.loads((ROOT/'evidence/algebra.json').read_text())
    ck('two_production_GT_samples',alg['production']['samples']==2)
    ck('24_guard_conjugacy_checks',sum(z['guard_conjugacy_checks'] for z in alg['production']['results'])==24)
    ck('24_Montgomery_key_checks',sum(z['canonical_key_Montgomery_pipeline_checks'] for z in alg['production']['results'])==24)
    ck('316368_alignment_pairs',alg['alignment']['pairs']==316368)
    ck('10947_baseline_cancellations',alg['forced_inverse_cycle']['exhaustive_label_cancellations']==10947)
    ck('no_false_answers',alg['fingerprint']['wrong_answers_accepted']==0 and alg['fingerprint']['false_pairs_rejected']>0)
    for k,v in alg['fingerprint']['miniature_universal_zero_counts'].items():ck('universal_hash_'+k,v==512)
    if not a.bootstrap:
        sub=json.loads((ROOT/'submission.json').read_text());note=(ROOT/'public-note.md').read_text()
        arg=(ROOT/'argument.md').read_text();st=sub['claims'][0]['statement']
        ck('manifest_exact_keys',set(sub)=={'version','title','claims','experiments'} and sub['version']==1)
        ck('claim_exact_keys',set(sub['claims'][0])=={'evidence','score','statement'} and sub['claims'][0]['score']==frac(CLAIM))
        ck('empty_experiments_dependencies',sub['experiments']==[] and json.loads((ROOT/'dependencies.json').read_text())==[])
        ck('self_contained_statement',all(w not in st for w in ('fc8e8','mpjunior','prior','improves','128.82','7da5f')))
        ck('note_10_to14_KiB',10240<=len(note.encode())<=14336)
        ck('full_attribution','GPT 6.1 Sol Max / Hermes' in note)
        ck('documents_core_terms',all(w in note and w in arg for w in ('8996','9000','109/100','Hmix-guard','128.81','0.39')))
        ck('note_lineage','fc8e8f88' in note and '7da5f23a' in note and 'mpjunior92' in note)
    paths=sorted(p for p in ROOT.rglob('*') if p.is_file() and p.suffix!='.pyc')
    ck('file_count_le32',len(paths)<=32)
    ck('files_le16KiB',all(p.stat().st_size<=16384 for p in paths if p.name!='certificate.json'))
    ck('payload_le256KiB',sum(p.stat().st_size for p in paths)<262144)
    ck('no_private_path_leak',all((chr(47)+'Users'+chr(47)) not in p.read_text() for p in paths if p.suffix!='.pyc'))
    bindings={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
        if p.name not in ('certificate.json','certificate-details.json','certificate_run.log','RESTART.md')}
    out=dict(version=1,package='ftwist-escape-rho-40pct-v8',
        scope='Exact arithmetic and captured toy consistency only. Hmix-guard and Hkey remain production hypotheses.',
        documents_pending=a.bootstrap,constants=dict(p=P,r=R,N=N,Np=NP,T0=T0,T=T,
            order2=o,order_mod_aut=d,claim=frac(CLAIM),rows_sum=sum(ROWS.values()),
            ceiling=ORDINARY_CEILING,W_executed=frac(WEXEC),W_scored=frac(W),W_repriced=frac(WSTRESS),
            ancillary_total=frac(at),fingerprint_pair_expectation=frac(fp)),
        log2_intervals=dict(executed=interval(WEXEC),scored=interval(W),repriced=interval(WSTRESS)),
        endpoint_columns=['N','B','control_base','guard_useful_repeat_stretch','guard_Floyd_stretch','guard_den_zero'],
        endpoints=fixed_counts,structural_cycle_columns=['B','baseline','inverse_guard','orientation_guard'],
        structural_cycles=cycle_mass,raw_captures_checked=a.raw_dir is not None,
        checks=checks,all_checks_pass=all(checks.values()),source_sha256=bindings)
    details={k:v for k,v in out.items() if k not in ('checks','all_checks_pass')}
    detail_text=json.dumps(details,separators=(',',':'))+'\n'
    root=dict(version=1,scope=out['scope'],documents_pending=a.bootstrap,
        details='evidence/certificate-details.json',
        details_sha256=hashlib.sha256(detail_text.encode()).hexdigest(),
        checks=checks,all_checks_pass=all(checks.values()))
    text=json.dumps(root,separators=(',',':'))+'\n'
    if a.write:
        (ROOT/'evidence/certificate-details.json').write_text(detail_text)
        (ROOT/'certificate.json').write_text(text)
    print(json.dumps(dict(checks=len(checks),passed=sum(checks.values()),bootstrap=a.bootstrap,
        log2=out['log2_intervals'],failed=[k for k,v in checks.items() if not v]),indent=2))
    if not all(checks.values()):raise SystemExit(2)
if __name__=='__main__':main()
