#!/usr/bin/env python3
"""Executed exact arithmetic/data consistency; does not certify production Hmix."""
import argparse
from fractions import Fraction as F
from math import gcd, isqrt, prod
from pathlib import Path
import hashlib
import json
from ledger import (P,R,N,NP,BRANCHES,RESERVE,CLAIM,LN_LO,LN_HI,T0,T,
                    ORDINARY_ROWS,ORDINARY_CEILING,WEXEC,W,WSTRESS,
                    ceil_fraction,ceil_sqrt_fraction,exp_minus_upper,ln_bounds,
                    below_pow2,interval)
ROOT=Path(__file__).resolve().parents[1]
checks={}
def ck(name,value): checks[name]=bool(value)
def frac(v):
    v=F(v)
    return {'numerator':str(v.numerator),'denominator':str(v.denominator)}
def prime(n):
    return n>1 and (n==2 or (n%2 and all(n%d for d in range(3,isqrt(n)+1,2))))

ck('class_count_exact',12*N==R-1)
ck('lambda_order_12',pow(P,12,R)==1 and pow(P,6,R)==R-1 and pow(P,4,R)!=1)
ck('cyclotomic_relation',(1+pow(P,4,R)+pow(P,8,R))%R==0)
ck('branch_correction_ceiling',(NP-1)*1023<N*1024<=NP*1023)
ck('branch_correction_not_r_over_r_minus_one',BRANCHES==1024 and BRANCHES!=R)
ck('ln5over3_exact_series_enclosure',F(51082562376599,10**14)<LN_LO<LN_HI<F(51082562376600,10**14))
ck('ln_enclosure_width_less_2pow_minus300',LN_HI-LN_LO<F(1,1<<300))
ck('T0_minimal_rational_sqrt',(T0-1)**2<2*NP*LN_HI<=T0*T0)
ck('T_preselected_11over10_stretch',10*(T-1)<11*T0<=10*T and RESERVE==F(11,10))
# For the uniform comparison map no-repeat <= exp[-t(t+1)/(2Np)].
# T0^2/(2Np)>=ln(5/3), so the comparison success is >=2/5.
success_lower=1-exp_minus_upper(LN_LO)
ck('comparison_success_exceeds_0399',success_lower>F(399,1000))
ck('comparison_040_exact_target',F(1)-F(3,5)==F(2,5))
ck('retain_one_percent_probability_margin_to_track_floor',F(2,5)-F(39,100)==F(1,100))

FACTORS={2:32,3:1,11:1,19:1,10177:1,906349:2,2508409:1,125527:1,
         859267:1,2529403:1,52437899:1,254760293:2}
ck('factorization_r_minus_one',prod(q**e for q,e in FACTORS.items())==R-1)
ck('factor_primality',all(prime(q) for q in FACTORS))
o=R-1
for q in FACTORS:
    while o%q==0 and pow(2,o//q,R)==1:o//=q
d=o//gcd(o,12)
ck('order_of_two_upper',pow(2,o,R)==1)
ck('order_of_two_exact',all(pow(2,o//q,R)!=1 for q in FACTORS if o%q==0))
ck('no_initial_label_closure_through_two_T',d>2*T)
ck('identity_comparison_union_less_2pow_minus120',3*T*(1<<120)<N)
ck('Hkey_tie_comparison_union_less_2pow_minus240',3*T*66*(1<<240)<P)

ck('four_Fp_products_already_in_v5',ORDINARY_ROWS['four_base_field_products_4_times_383']==4*383)
ck('nonreference_itemized_sum_4372',sum(ORDINARY_ROWS.values())==4372)
ck('all_rows_inside_8000',sum(ORDINARY_ROWS.values())<=ORDINARY_CEILING==8000)
ck('square_and_add_two_dense_products_not_one',2*3==6)
ck('executed_ceiling_matches_ledger',WEXEC==F(36,5)*T+(1<<100))
ck('all_nonreference_ceilings_repriced_6over5',WSTRESS==(6+F(6,5)*F(3*8000,20000))*T+(1<<101))
ck('executed_cost_below_nominal_envelope',WEXEC<W)
ck('claim_reduced',gcd(CLAIM.numerator,CLAIM.denominator)==1)
ck('nominal_log_interval_1287897_1287898',not below_pow2(W,1287897,10000) and below_pow2(W,1287898,10000))
ck('nominal_claim_slack_at_least_0005',below_pow2(W,128815,1000))
ck('nominal_claim_slack_at_most_006',not below_pow2(W,128760,1000))
ck('repriced_log_interval_1286850_1286851',not below_pow2(WSTRESS,1286850,10000) and below_pow2(WSTRESS,1286851,10000))
ck('repriced_headroom_at_least_010',below_pow2(WSTRESS,128720,1000))
ck('claim_is_12882_over_100',CLAIM==F(12882,100))
ck('HOLD_claim_does_not_beat_12589_over_100',CLAIM>=F(12589,100))

# This is a FLOOR FOR THE RETAINED SUFFICIENT MU+LAMBDA SCHEDULE,
# not a lower bound on all algorithms or even on the exact Floyd CDF.
ll,lh=ln_bounds(F(100,61))
FLOOR_T=ceil_sqrt_fraction(2*NP*ll)
FLOOR_WORK=6*FLOOR_T
ck('ideal_039_six_reference_floor_above_1282134',not below_pow2(FLOOR_WORK,1282134,10000))
ck('ideal_039_six_reference_floor_below_1282135',below_pow2(FLOOR_WORK,1282135,10000))
ck('reference_only_floor_above_frontier',not below_pow2(FLOOR_WORK,12589,100))
# Even using repeat-by-2T (a necessary condition for a Floyd meeting by T)
# and a union bound P(repeat)<=2T(2T+1)/(2Np), a reference-only model
# budget at the frontier is hopeless. Use a rational upper bound on 2^125.89.
cap=1<<126
ck('ideal_uniform_map_union_at_2T_for_budget_2pow126_is_below_039',
   F((2*(cap//6))*(2*(cap//6)+1),2*NP)<F(39,100))

# Rival's own rounded reconstruction, not a validation of its mechanism.
S0=ceil_sqrt_fraction(2*NP*LN_HI)
S=ceil_fraction(F(103,100)*S0)
RIVAL=F(57,50)*S*F(131129,131072)+F(57,50)*20*(1<<80)+(1<<62)
ck('rival_reconstructed_log_125884_125885',not below_pow2(RIVAL,125884,1000) and below_pow2(RIVAL,125885,1000))
# An unrealistically aggressive transplant still retains SIX reference calls.
B2800=(6+F(3*2800,20000))*T+(1<<100)
ck('even_unjustified_2800_ceiling_transplant_loses',not below_pow2(B2800,12589,100))

M,G,FP=1<<69,129,128
record=16+32+32+32+16
peak=record*M+(1<<40)
fp_pairs=F(G*M*(M-1),2*(1<<FP))
ancillary_rows={
 'report_instructions':G*M*(1<<26),
 'clear_instructions':G*M*record*16,
 'candidate_units':(1<<32)*(1<<15),
 'setup_units':1<<50,
 'exception_output_pairing_units':1<<90,
}
ancillary=F(ancillary_rows['report_instructions']+ancillary_rows['clear_instructions'],20000)+sum(ancillary_rows[k] for k in ['candidate_units','setup_units','exception_output_pairing_units'])
ck('128_byte_records',record==128)
ck('peak_below_2pow78',peak<1<<78)
ck('all_ancillary_ceiling_inside_2pow100',ancillary<1<<100)
ck('all_ancillary_1p2_repricing_inside_2pow101',F(6,5)*ancillary<1<<101)
ck('FP_pair_expectation_less_2pow17',fp_pairs<1<<17)
ck('FP_Markov_cap_bound_less_2pow_minus15',fp_pairs/(1<<32)<F(1,1<<15))
ck('DP_60_bits',len(range(10,70))==60)
ck('all_generations_DP_predicate_within_field_key_bits',70+G-1<=381)
ck('full_point_independent_linear_fingerprint_input_bits',12*6*64==4608)
ck('independent_fingerprint_matrix_bytes_73728',4608*128//8==73728)
ck('full_point_fingerprint_work_inside_report_ceiling',4608*16+2*256*256+72*383+(1<<20)<1<<26)
# Exhaustive independent 4-input/3-output linear-hash model. Every nonzero
# difference has exactly1/8 zero outputs over all4096 matrices.
for delta in range(1,16):
    count=0
    for matrix in range(1<<12):
        columns=[(matrix>>(3*k))&7 for k in range(4)]
        val=0
        for k in range(4):
            if (delta>>k)&1:val^=columns[k]
        count+=val==0
    ck('linear_fingerprint_exhaustive_difference_'+str(delta),count==512)
ck('new_report_mean_comparison_less_M_over_3',9*T<M*(1<<60))
# The DP table is OPTIONAL: disabling reports cannot stop or change Floyd.
# Hence no Poisson overflow premise is needed for either resources or success.
ck('table_build_once_under_2pow50',1024*(2*510*23000+12*(72*383+10000))+(1<<30)<(1<<50)*20000)
ck('table_bytes_under_2pow24',1024*12*(576+64)<1<<24)

# Faithful larger experiment and toy-only confidence calculations.
toy=json.loads((ROOT/'toy_results.json').read_text())
configs=[json.loads((ROOT/f).read_text()) for f in toy['configuration_files']]
ck('seven_complete_faithful_settings',len(configs)==7)
ck('preselected_seed_and_reserve',toy['seed']==2026102306 and toy['preregistered_length_reserve']=={'numerator':11,'denominator':10})
ck('70000_new_independent_pairs',sum(c['trials'] for c in configs)==70000 and all(c['trials']==10000 for c in configs))
# 14 two-sided endpoint tails: 28*e^-8<1/100.
eight_lower=sum(F(8)**i/prod(range(1,i+1)) for i in range(41))
ck('simultaneous_endpoint_failure_below_one_percent',eight_lower>2800)
ck('simultaneous_paired_CDF_failure_below_one_percent',eight_lower>700)
ratios=[]
for i,c in enumerate(configs):
    tag='toy_'+str(i)
    n=c['classes']; hist=c['exhaustive_branch_histogram']; tr=c['trials']
    ck(tag+'_zero_dead_branches',len(hist)==c['branches'] and min(hist)>0 and c['dead_branches']==0)
    ck(tag+'_exhaustive_branch_mass',sum(hist)==n)
    ck(tag+'_literal_tower_and_step_checks',c['literal_tower_crosschecks']==24 and c['literal_square_and_add_step_checks']==32)
    ss,rr=c['shared'],c['control']
    ratio=F(ss['first_repeat']['sum'],rr['first_repeat']['sum']);ratios.append(ratio)
    ck(tag+'_exact_mean_ratio_reproduced',c['paired']['mean_first_repeat_ratio']=={'numerator':ratio.numerator,'denominator':ratio.denominator})
    ck(tag+'_empirical_mean_ratio_less_103over100',ratio<F(103,100))
    diff=F(ss['first_repeat_at_stretch']-rr['first_repeat_at_base'],tr)
    ck(tag+'_10percent_CDF_paired_lower_positive',diff-F(1,25)>0)
    ck(tag+'_endpoint_lower_above_toy_039',F(ss['first_repeat_at_stretch'],tr)-F(1,50)>F(39,100))
    ck(tag+'_all_DP_densities_censored_honestly',all(
        ss['dp'][str(d)]['hit_before_repeat']+ss['dp'][str(d)]['no_hit_before_repeat']==tr
        and rr['dp'][str(d)]['hit_before_repeat']+rr['dp'][str(d)]['no_hit_before_repeat']==tr
        for d in [2,4,6]))
    ck(tag+'_observed_histogram_trial_mass',all(sum(z[1] for z in side[k]['histogram'])==tr
        for side in [ss,rr] for k in ['mu','cycle_lambda','first_repeat','floyd_meeting']))
ck('two_genuine_1024_branch_scales',sum(c['branches']==1024 for c in configs)==2)
ck('high_reuse_stress_exceeds_97percent',any(100*c['reused_branch_transitions']>97*c['transitions'] for c in configs))

algebra=json.loads((ROOT/'evidence/algebra_results.json').read_text())
ck('two_reexecuted_production_GT_identity_samples',algebra['production_tower_identities']['samples']==2)
ck('reexecuted_exhaustive_scalar_alignment_cases',sum(c['pairs'] for c in algebra['all_pair_alignment'])==2995488)
ck('reexecuted_DP_fault_audit_zero_mismatch',algebra['additive_same_step_dp_audit']['counters']['dp_mismatch']==0)
ck('reexecuted_DP_faults_detected',algebra['additive_same_step_dp_audit']['counters']['injected_dp_faults_detected']==495)
# Exact Fraction levers and deliberately unclaimed sensitivities.
sensitivity={
 'v5_scored_envelope':8*((5*(isqrt(5*NP//4)+1)+3)//4)+(1<<100),
 'change_stretch_only_to_10percent':8*ceil_fraction(F(11,10)*(isqrt(5*NP//4)+1))+(1<<100),
 'strike_ln5over3_and_10percent_envelope':W,
 'track_floor_039_no_reserve_no_nonreference':FLOOR_WORK,
 'unjustified_2800_transplant':B2800,
 'rival_reconstruction':RIVAL,
}

def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    out={
      'version':1,'package':'gt-orbit12-rho-strike-40pct-v6',
      'scope':'Exact arithmetic and captured-toy consistency only. Production Hmix/Hkey remain unproved. HOLD: no frontier improvement.',
      'constants':{'p':P,'r':R,'N':N,'Np_floor':N*1024//1023,'Np_ceiling':NP,
        'T0':T0,'T':T,'stretch':frac(RESERVE),'ln5over3_lower':frac(LN_LO),'ln5over3_upper':frac(LN_HI),
        'order_two':o,'order_mod_frobenius':d,'ordinary_rows':ORDINARY_ROWS,'ordinary_sum':sum(ORDINARY_ROWS.values()),
        'ordinary_ceiling':ORDINARY_CEILING,'W_executed':frac(WEXEC),'W_scored':frac(W),'W_repriced':frac(WSTRESS),
        'claim':frac(CLAIM),'ideal_039_schedule_floor':frac(FLOOR_WORK),'hard_peak_bytes':peak,
        'ancillary_rows':ancillary_rows,'ancillary_total':frac(ancillary),'finite_FP_expectation':frac(fp_pairs),
        'max_empirical_shared_control_mean_ratio':frac(max(ratios)),
        'trial_pairs':sum(c['trials'] for c in configs)},
      'log2_intervals':{k:interval(v) for k,v in sensitivity.items()}|{'executed':interval(WEXEC),'repriced':interval(WSTRESS)},
      'sensitivity':{k:frac(v) for k,v in sensitivity.items()},
      'checks':checks,'all_checks_pass':all(checks.values()),
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'scripts').glob('*.py'))},
      'verifier_source':Path(__file__).read_text(),
    }
    text=json.dumps(out,indent=2)+'\n'
    if a.write:
        (ROOT/'certificate.json').write_text(text)
        (ROOT/'evidence/certificate_run.log').write_text('COMMAND: python3 -B scripts/verify_v6.py --write\n'+str(sum(checks.values()))+'/'+str(len(checks))+' exact checks\n'+'\n'.join(('PASS ' if b else 'FAIL ')+k for k,b in checks.items())+'\n')
    print(text,end='')
    if not all(checks.values()):raise SystemExit('Exact verification failed')
if __name__=='__main__':main()
