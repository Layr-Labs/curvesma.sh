"""Regenerate exact-integer/Fraction certificate; never runs production search.
python3 -B verify_certificate.py --write  # author-side capture
python3 -B verify_certificate.py          # compare exact captured objects
PASS checks arithmetic and captured evidence, NOT production ERSE premises.
"""
import hashlib
import json
import math
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import ledger as z
HERE=Path(__file__).resolve().parent
checks={}
def ck(name,truth):
    if not truth: raise AssertionError(name)
    checks[name]=True

def lt(x,a,b=1):
    x=F(x); return x.numerator**b<2**a*x.denominator**b

def gt(x,a,b=1):
    x=F(x); return x.numerator**b>2**a*x.denominator**b

ck('N_full_precision',z.N==19298681539552699237261830834781317975472927379845817397100860523586360249056)
ck('class_partition',6*z.N==z.r-1)
ck('beta_order3',z.beta!=1 and (z.beta*z.beta+z.beta+1)%z.p==0)
ck('lambda_order3',z.lam!=1 and (z.lam*z.lam+z.lam+1)%z.r==0)
ck('six_units',len({s*pow(z.lam,j,z.r)%z.r for s in (-1,1) for j in range(3)})==6)
ck('seven_nonresidue',pow(7,(z.p-1)//2,z.p)==z.p-1)
ck('nu_exact_two',z.nu==2)
ck('ceil_sqrt_upper',z.T*z.T>=z.budget_square)
ck('ceil_sqrt_minimal',(z.T-1)**2<z.budget_square)
ck('flush_charged',z.K-z.T==z.workers*z.L)
ck('trail_counter64',z.L<2**64)
ck('global_counter128',z.K<2**128)
ck('step_rows_sum925',sum(z.step_rows.values())==925)
ck('E0_covers_rows',sum(z.step_rows.values())<=z.E0)
ck('step_repriced_inside',z.E==F(6,5)*z.E0)
ck('every_other_row_repriced',all(z.rows[k]==F(6,5)*v for k,v in z.base_other.items()))
ck('doubling_charged_two_ticks',F(3,2)+F(z.E,z.C)<=2*(1+F(z.E,z.C)))
ck('whole_ledger_reserve8pct',z.W==F(27,25)*z.raw)
ck('claim_reduced',z.claim==F(2579,20))
ck('total_below_claim',lt(z.W,2579,20))
ck('roundup_more_than005',lt(z.W,25789,200))
ck('roundup_less_than060',gt(z.W,12889,100))
ck('repriced_headroom_more_than010',lt(z.raw,2577,20))
ck('raw_log_lower',gt(z.raw,1287968,10000))
ck('raw_log_upper',lt(z.raw,1287969,10000))
ck('reserved_log_lower',gt(z.W,1289078,10000))
ck('reserved_log_upper',lt(z.W,1289079,10000))
ck('record192',24+64+8+36+1+1+58==192)
ck('fixed_pool_bytes',z.M*192==3*2**79)
ck('hard_memory',z.M*192+2**30<2**81)
ck('pointer96',z.M*192<2**96)
ck('table_and_states_fit',z.R*160+z.workers*1024+2**24<2**30)
ck('partition_dp_disjoint',58<=64)
ck('lookup_charge',2**16+148*1000+2**14<2**18)
ck('deepen_charge',F(3*z.M*2**20,z.C)<z.base_other['deepen_rebuilds'])
ck('replay_charge',F(z.A*2*z.L*2**17,z.C)+2**30<z.base_other['all_key_match_replays'])
ck('exception_charge',F(2**16*2**20,z.C)<z.base_other['exceptional_arithmetic'])
ck('preprocessing_included',F(6,5)*(z.base_other['table_preprocessing']+z.base_other['constant_checks_sampling_output'])<2**33)
ck('markov_pool_formula',z.count_pool==F(11,10)*F(z.K,2**58*z.M))
ck('markov_events_starts_formula',z.count_events_starts==F(11,10)*F(z.K,2**56*(z.S-z.workers)))
ck('markov_loss_disclosed',F(334,1000)<z.count_pool+z.count_events_starts<F(335,1000))
ck('fingerprint_pair_union',z.false_key<F(1,2**45))
ck('relation_exact_uniform_label',z.relation==F(1,z.r))
ck('sampling_union_below430',z.sampling<F(1,2**430))
ck('exception_union_below118',z.exception<F(1,2**118))
ck('detector_union_below120',z.detector<F(1,2**120))
ck('exp17_gt_2pow24',sum((F(17**k,math.factorial(k)) for k in range(70)),F(0))>2**24)
ck('production_age_window_exponent',20*F(9,10)*z.velocity>17)
ck('production_wait_second_moment',2*(F(10,9)/z.velocity)**2<3)
ck('production_age_second_moment',6*(F(10,9)/z.velocity)**2<8)
ck('tail_union_itemized',z.trail_cap==F(44,2**24))
ck('production_core_law',z.renewal_core>F(74998,100000))
ck('production_eta_sparse',(2**58)**2*2**116<z.N)
ck('production_transport_pair_bound',F(64*2**116,z.N)<F(1,2**128))
# Prove ERSE complete rho >=0.74997 without a floating-point logarithm.
a=F(z.T*(z.T-1),2*z.N)
_,e_model_hi=z.exp_minus_bounds(z.renewal_core*a)
e_floor_lo,_=z.exp_minus_bounds(F(74997,100000)*a)
ck('production_complete_ERSE_floor',e_model_hi+F(64*2**116,z.N)+z.trail_cap<e_floor_lo)
transfer=[]
for d in (20,36,56,58):
    delta=F(64*(1<<d)**2,z.N)+z.trail_cap
    ck(f'production_transfer_d{d}',e_model_hi+delta<e_floor_lo)
    transfer.append(dict(scope='EXTRAPOLATED ERSE model, not measurement',
        table_entries=str(z.R),dp_depth=d,spacing_table_ratio=z.frac(F(1<<d,z.R)),
        endpoint_loss=z.frac(delta),complete_rho_lower=z.frac(F(74997,100000))))
ck('working_quality_below_ERSE',z.rho<F(74997,100000))
ck('success_more_than0414',z.success()>F(414,1000))
ck('success_claim040',z.success()>F(2,5))
ck('all_losses_additive',set(z.losses)=={'pool_markov','event_start_markov','fingerprint_union','relation_failure','sampling_rejections','exception_cap','unforced_detector_abort','collision_transport_trail_cap'})
sub=json.loads((HERE/'submission.json').read_text())
ck('manifest_schema_exact',set(sub)=={'version','title','claims','experiments'})
ck('claim_schema_exact',set(sub['claims'][0])=={'statement','evidence','score'})
ck('submission_score_matches',sub['claims'][0]['score']==z.frac(z.claim))
ck('no_declared_experiments',sub['experiments']==[])
ck('dependencies_empty',json.loads((HERE/'dependencies.json').read_text())==[])
summary=json.loads((HERE/'sweep-summary.json').read_text())
ck('sweep_points24',len(summary['points'])==24)
ck('sweep_trials48000',sum(v['trials'] for v in summary['points'])==48000)
ck('toy_group_scale',summary['group']['classes']==87328 and summary['group']['order']==523969)
ck('toy_fixed_partition',summary['parameters']['partition_shift']>=max(summary['parameters']['depths']))
analysis=[]; evidence_files=['sweep-summary.json']
for R in (8,32,256,1024):
    name=f'branches-R{R}.json'; h=json.loads((HERE/name).read_text()); evidence_files.append(name)
    ck(f'R{R}_all_branches_live',len(h['full_class_histogram'])==R and min(h['full_class_histogram'])>0 and sum(h['full_class_histogram'])==87328)
    for d in (4,6,8,10,12,14):
        low=4+4*((d-4)//4)
        name=f'sweep-R{R}-d{low}-{low+2}.json'
        if name not in evidence_files: evidence_files.append(name)
        pair=json.loads((HERE/name).read_text())
        v=next(v for v in pair if v['dp_bits']==d); label=f'R{R}_d{d}'
        cc=Counter()
        for row in v['table_rows']:
            cc.update(row['counters'])
            assert sum(row['dp_event_histogram'].values())==row['trials']
        ck(label+'_2000_trials',len(v['table_rows'])==8 and sum(r['trials'] for r in v['table_rows'])==v['trials']==2000)
        ck(label+'_counts_reconciled',dict(sorted(cc.items()))==v['counters'] and sum(r['successes'] for r in v['table_rows'])==v['successes'])
        ordinary=cc['ticks']-2*cc['escapes']-cc['rewind_ticks']
        ck(label+'_escape_reward_envelope',F(2*cc['escapes']+cc['rewind_ticks'],ordinary)<=F(8,R))
        AT=F(318*317,2*87328); quality=F(3,4)/(1+F(8,R))**2
        elo,ehi=z.exp_minus_bounds(quality*AT)
        core=F(cc['core_collision_trials'],2000); pipeline=F(v['successes'],2000)
        ck(label+'_core_law_empirical_compatibility',core>=1-elo)
        delta=min(F(1),F(64*(1<<d)**2,87328)+z.trail_cap)
        ck(label+'_transport_law_empirical_compatibility',pipeline>=max(F(0),1-ehi-delta))
        if d in (4,6): ck(label+'_dense_first_moment_110',F(cc['dp_events'],cc['dp_tests'])<=F(11,10*2**d))
        analysis.append(dict(width=R,dp_bits=d,ratio=v['ratio'],trials=2000,successes=v['successes'],
            measured_complete_rho=v['rho'],measured_core_rho=-math.log1p(-float(core))/float(AT),
            nominal_DP_density_ratio=v['measured_dp_event_density']*2**d,
            ERSE_pipeline_probability_lower_display=float(max(F(0),1-ehi-delta))))
values=dict(integers={k:str(getattr(z,k)) for k in ('p','r','beta','lam','N','T','K','M','S','A','L','R','workers')},
    parameters=dict(nu=z.frac(z.nu),rho=z.frac(z.rho),balance=z.frac(z.balance),claim=z.frac(z.claim)),
    ledger={k:z.frac(v) for k,v in z.rows.items()},raw=z.frac(z.raw),reserved=z.frac(z.W),
    losses={k:z.frac(v) for k,v in z.losses.items()},
    success_lower=z.frac(z.success()),displays=dict(log2_raw=z.log_display(z.raw),log2_reserved=z.log_display(z.W),
    success_lower=float(z.success()),loss_sum=float(sum(z.losses.values(),F(0)))))
out=dict(version=1,scope='Exact conditional arithmetic and captured toy reconciliation only; no production-premise certification.',
    evidence_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in sorted(evidence_files)},
    checks=checks,pass_count=len(checks))
for name,obj in (('certificate.json',out),('certificate-values.json',values),('sweep-analysis.json',analysis),('transfer-table.json',transfer)):
    data=json.dumps(obj,separators=(',',':'))+'\n'
    if '--write' in sys.argv: (HERE/name).write_text(data)
    else: assert (HERE/name).read_text()==data, 'Captured artifact mismatch: '+name
print(json.dumps(dict(verdict='PASS',certificate_checks=len(checks),claim='2579/20',T=str(z.T),
    log2_raw=z.log_display(z.raw),log2_reserved=z.log_display(z.W),success_lower=float(z.success()))))
