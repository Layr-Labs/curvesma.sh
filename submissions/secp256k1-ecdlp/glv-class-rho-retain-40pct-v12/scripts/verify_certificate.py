"""Exact author certificate with fail-closed fixed-spacing qualification.
Run --write to capture certificate.json and evidence/certificate_run.log;
run without it to compare regenerated artifacts. Exit 2 means HOLD, not PASS.
Fractions/integers prove bounds. Display statistics in captured JSON are not
used as production probability evidence. This verifier is embedded verbatim.
"""
import hashlib
import json
import math
import re
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from typing import Any
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import ledger as z
from scripts.transports import delta_transport
from scripts.toy_sweep import binding_budget, configs
from scripts.renewal_law import finite_reserve, lower_envelope, path_ratios, tail_rows, PATHS


def rational(v):
    return F(int(v['numerator']),int(v['denominator']))


def below(x,a,b=1):
    x=F(x)
    return x.numerator**b < (1<<a)*x.denominator**b


def above(x,a,b=1):
    x=F(x)
    return x.numerator**b > (1<<a)*x.denominator**b


def retention_audit() -> dict[str, Any]:
    """Certify the SAME fixed-density first moments, conditional on H-rate.

    Let D count all emitted endpoints before the K cap. At most D distinct
    records can be stored, so abort-on-full is contained in {D>=M}. Markov
    gives balance*K/(DP*M), with no independence or replacement premise.
    Terminal/start counts use their separately stated H-rate first moment.
    """
    fill=F(z.K,z.DP*z.M)
    pool=z.balance*fill
    starts=z.balance*F(z.K,z.DP*(z.S-z.workers))
    corrected=z.success()+z.count_pool+z.count_events_starts-pool-starts
    return dict(checks=dict(fixed_spacing_pool_rate_matches=z.count_pool==pool,
        fixed_spacing_start_rate_matches=z.count_events_starts==starts,
        fixed_spacing_fill_matches=z.fill==fill,
        pool_loss_at_most_001=pool<=F(1,100),
        event_start_loss_at_most_0005=starts<=F(1,200),
        fixed_spacing_success_040=corrected>=F(2,5),
        no_evictions_abort_on_full=z.pool_policy=='abort_on_full_no_evictions',
        both_markov_denominators_actual=z.count_pool*z.M==z.count_events_starts*(z.S-z.workers)),
        actual_spacing=str(z.DP),ledger_pool_spacing=str(z.DP),ledger_start_spacing=str(z.DP),
        nominal_fill=z.frac(fill),expected_emitted_upper=z.frac(z.balance*F(z.K,z.DP)),
        fixed_spacing_pool_markov=z.frac(pool),fixed_spacing_start_markov=z.frac(starts),
        corrected_unclipped_union_lower=z.frac(corrected),
        clipped_probability_lower=z.frac(max(F(0),corrected)),
        issue='Conditional capacity-retention arithmetic passes at the actual fixed spacing; no endpoint is evicted or filtered. H-rate/H-window remain disclosed production hypotheses, not theorems from the toy sweep.')


def build_certificate() -> dict[str, Any]:
    checks: dict[str,bool]={}
    def ck(name,value): checks[name]=bool(value)
    ck('production_class_count',6*z.N==z.r-1)
    ck('production_beta_order3',z.beta!=1 and (z.beta*z.beta+z.beta+1)%z.p==0)
    ck('production_lambda_order3',z.lam!=1 and (z.lam*z.lam+z.lam+1)%z.r==0)
    ck('six_distinct_units',len({s*pow(z.lam,j,z.r)%z.r for s in (-1,1) for j in range(3)})==6)
    ck('C1000',z.C==1000)
    ck('nu137over100',z.nu==F(137,100))
    ck('ceil_core_upper',z.T*z.T>=z.budget_square)
    ck('ceil_core_minimal',(z.T-1)**2<z.budget_square)
    ck('drain_charged',z.K-z.T==z.workers*z.L)
    ck('step2432',sum(z.step_rows.values())==2432)
    ck('two_limb_length_24instructions',z.step_rows['two_limb_length_control']==24)
    ck('step_budget8000',sum(z.step_rows.values())<=z.E0==8000)
    ck('step_repriced9600',z.E==F(6,5)*z.E0==9600)
    ck('all_other_rows_repriced',all(z.rows[k]==F(6,5)*v for k,v in z.base_other.items()))
    ck('lookup_repriced_actual_S',z.base_other['dp_lookup_transfers']==F(z.S*2**18,z.C))
    ck('starts_repriced_actual_S',z.base_other['fresh_starts']==z.S*4096)
    ck('retention_maintenance_2pow96',z.base_other.get('retention_maintenance')==F(2**96))
    ck('zero_deepening_rows',not any('deepen' in k for k in z.rows))
    ck('raw_sums',z.raw==sum(z.rows.values(),F(0)))
    ck('reserve27over25',z.W==F(27,25)*z.raw)
    ck('claim_reduced',z.claim==F(13068,100) and math.gcd(z.claim.numerator,z.claim.denominator)==1)
    ck('claim_below_006_window',above(z.W,13062,100))
    ck('claim_below_user133',z.claim<=133)
    ck('claim_below_gate13662',z.claim<F(6831,50))
    ck('total_below_claim',below(z.W,13068,100))
    ck('raw_headroom_010',below(z.raw,13058,100))
    ck('raw_log_lower1305594',above(z.raw,1305594,10000))
    ck('raw_log_upper1305595',below(z.raw,1305595,10000))
    ck('reserved_log_lower1306704',above(z.W,1306704,10000))
    ck('reserved_log_upper1306705',below(z.W,1306705,10000))
    ck('fixed_DP64',z.d==64 and z.DP==2**z.d)
    ck('fixed_branch_disjoint',z.branch_shift>=z.d and z.branch_shift+20<=256)
    ck('trail_cut_actual_spacing',z.L==40*z.DP)
    ck('two_limb_length_fits',z.L<2**128 and z.record_fields['length']==16)
    ck('global_two_limb_ticks_fit',z.K<2**128)
    ck('storage192_bytes',sum(z.record_fields.values())==z.record_bytes==192)
    ck('full_scalar_entropy64bytes',z.record_fields['starting_scalars']==64)
    ck('storage_pool_3times2pow78',z.memory_rows['node_pool']==z.M*192==3*2**78)
    ck('memory_peak_row_sums',z.peak_memory==sum(z.memory_rows.values()))
    ck('storage_total_below2pow80',z.peak_memory<=2**80)
    ck('memory_aux_cap2pow30',z.memory_rows['auxiliary_cap']==2**30)
    ck('M_records2pow72',z.M==2**72)
    ck('S_events2pow85',z.S==2**85)
    ck('tree_height_under146',2*(z.M+1).bit_length()<=146)
    ck('replay_ceiling_covers_longer_L',2*(z.A+2**16)*z.L*2**10+(z.A+2**16)*2**12<=z.base_other['all_key_match_replays'])
    ck('transport_exact',z.transport==F(64*z.DP*z.DP,z.N))
    ck('transport_under1e30',z.transport<F(1,10**30))
    ck('detector_recomputed',z.detector==F(8*z.L*z.L,z.N))
    ck('fingerprint_actual_S',z.false_key==F(z.S*(z.S-1),2*2**192))
    ck('conditional_success_040',z.success()>=F(2,5))
    ck('named_margin_001_over039',F(2,5)-F(39,100)==F(1,100))
    ck('conditional_margin002_over040',z.success()-F(2,5)>F(1,500))
    sub=json.loads((HERE/'submission.json').read_text())
    ck('manifest_exact_keys',set(sub)=={'version','title','claims','experiments'})
    ck('claim_exact_keys',len(sub['claims'])==1 and set(sub['claims'][0])=={'evidence','score','statement'})
    ck('score3267over25',sub['claims'][0]['score']==z.frac(z.claim))
    ck('experiments_empty',sub['experiments']==[])
    ck('dependencies_empty',json.loads((HERE/'dependencies.json').read_text())==[])
    statement=sub['claims'][0]['statement']
    ck('claim_self_contained_no_lineage',not re.search(r'6255b11|136\.62|improv|v9|submission [a-f0-9]',statement,re.I))
    argument=(HERE/'argument.md').read_bytes(); note=(HERE/'public-note.md').read_bytes()
    ck('argument_under8192',len(argument)<8192)
    ck('note_10to14KiB',10240<=len(note)<14336)
    ck('note_no_private_paths',b'/Users/' not in note)
    ck('note_lineage_attribution',b'qwen3.8-flash-next / Hermes' in note)
    ck('note_packaging_attribution',b'GPT 6.1 Sol' in note)
    ck('note_context_lineage',b'6255b11' in note and b'136.62' in note)
    summary=json.loads((HERE/'sweep-summary.json').read_text())
    points=summary['points']
    ck('sweep_bound_to_production_d',summary['parameters']['production_d']==z.d and summary['parameters']['production_spacing']==str(z.DP))
    ck('three_scales',set(summary['groups'])=={'small','medium','large'})
    ck('configs24',len(points)==24)
    ck('trials48000',sum(p['trials'] for p in points)==48000)
    records=[]
    for path in sorted((HERE/'evidence').glob('config-*.json')):
        packed=json.loads(path.read_text())
        for row in packed:
            for t in row['tables']:
                t['counters']={k:v for k,v in zip(row['table_counter_fields'],t.pop('counts')) if v}
        records.extend(packed)
    ck('all_captured_configs',len(records)==24 and {p['id'] for p in records}=={p['id'] for p in points})
    by_id={p['id']:p for p in records}
    for scale,g in summary['groups'].items():
        ck(scale+'_full_enumeration',g['subgroup_points_enumerated']==g['order']-1==6*g['classes'])
        ck(scale+'_GLV_verified',g['subgroup_order_verified'] and g['matching_glv_all_nonidentity'])
        ck(scale+'_addition_fixtures2000',g['independent_group_addition_fixtures']==2000)
        pp=[p for p in points if p['scale']==scale]
        ck(scale+'_configs8',len(pp)==8)
        expected=configs()[scale]
        ck(scale+'_spacing_comoved_configs',[(p['width'],p['dp_bits'],p['variant']) for p in pp]==expected)
        ck(scale+'_transport_ratio1e30',all(rational(p['delta_transport'])/z.transport>10**30 for p in pp))
    for p in points:
        label=p['id']; full=by_id[label]; n=p['classes']; d=p['dp_bits']; T=p['T']
        delta=delta_transport(d,classes=n); budget,envelope=binding_budget(n,d)
        ck(label+'_2000trials',p['trials']==2000==sum(t['trials'] for t in full['tables']))
        ck(label+'_8tables',len(full['tables'])==8)
        ck(label+'_delta_exact',rational(p['delta_transport'])==delta)
        ck(label+'_delta_under030',delta<F(3,10))
        ck(label+'_delta_ratio_exact',rational(p['delta_ratio_to_production'])==delta/z.transport)
        ck(label+'_minimal_binding_budget',T==budget and rational(p['envelope_lower'])==envelope)
        ck(label+'_envelope_positive',envelope>F(1,50))
        ck(label+'_binding_empirical',F(p['recovered'],p['trials'])>=envelope)
        hist=full['branch_histogram']
        ck(label+'_all_branches_live',len(hist)==p['width'] and min(hist)>0 and sum(hist)==n)
        ck(label+'_DP_disjoint',p['partition_shift']>=d)
        ck(label+'_no_wrong_recovery',p['wrong']==0)
        ck(label+'_outcomes_reconcile',p['recovered']+p['no_recovery']==p['trials'])
        cc=Counter()
        for row in full['tables']: cc.update(row['counters'])
        ck(label+'_counters_reconcile',dict(sorted(cc.items()))==full['counters'])
        ck(label+'_success_reconcile',cc['recovered']==p['recovered'])
        ck(label+'_first_repeat_exact_sum',sum(t['repeat_samples'] for t in full['tables'])==p['repeat_samples'] and sum(t['repeat_sum'] for t in full['tables'])==p['repeat_sum'])
        ck(label+'_summary_matches',all(full[k]==v for k,v in p.items()))
    # Every new probability input is recomputed from collision-selected counts.
    rr=[]
    for path in sorted((HERE/'evidence').glob('renewal-*.json')):
        rr.extend(json.loads(path.read_text()))
    ck('renewal24',len(rr)==24 and {r['id'] for r in rr}==set(by_id))
    ck('eta_measured_directed',z.eta==lower_envelope(rr)==rational(summary['renewal_lower_envelope']))
    ck('eta_nonvacuous',F(4,5)<z.eta<=1)
    # exp(10) exceeds this positive finite Taylor sum, so union failure of
    # 120 conditional finite-mixture Hoeffding cells is less than 0.01.
    exp10=sum((F(10**j,math.factorial(j)) for j in range(25)),F(0))
    ck('finite_Hoeffding_union_99pct',F(120,1)/exp10<F(1,100))
    for row in rr:
        name=row['id']; p=by_id[name]; n=row['selected']; eps=finite_reserve(n)
        ck(name+'_select',n==p['counters']['core_repeat_trials'] and 0<n<=p['trials'])
        ck(name+'_complete',0<=row['completed_pairs']<=n)
        ck(name+'_age',sum(v[-1] for v in row['age_hist'])==n)
        ck(name+'_renew_age',sum(v[-1] for v in row['renewal_age_hist'])==n)
        ck(name+'_gap',sum(v[-1] for v in row['inter_escape_gap_hist'])==row['gap_events']==p['counters'].get('escapes',0))
        ck(name+'_finite_reserve',eps*eps*n>=5 and (eps-F(1,10000))**2*n<5)
        ck(name+'_tails',len(row['tails'])==4 and [v[0] for v in row['tails']]==[1,2,4,40] and all(0<=v[1]<=n for v in row['tails']))
        ck(name+'_tail_lower',all(t['lower']>=z.eta*t['null'] for t in tail_rows(row)))
        ck(name+'_completion_lower',F(row['completed_pairs'],n)-eps>=z.eta)
        ck(name+'_paths_exhaustive',set(row['terminal_counts'])==set(PATHS) and sum(row['terminal_counts'].values())==p['counters']['starts'])
        ck(name+'_paths_exact',all(row['terminal_counts'][k]==p['counters'].get('terminal_'+k,0) and row['path_trial_counts'][k]==p['path_trial_counts'].get(k,0) for k in PATHS))
        ck(name+'_selected_paths',sum(row['selected_terminal_counts'].values())==2*n)
        for path,ratio in path_ratios(row,p).items(): ck(name+'_dom_'+path,ratio<=1)
        terminals=sum(row['terminal_counts'].values()); ticks=p['counters']['ticks']
        ck(name+'_terminal_rate',F(terminals*p['spacing'],ticks)<=F(11,10))
        ck(name+'_DP_rate',F(p['counters']['dp_events']*p['spacing'],ticks)<=F(11,10))
    # This is an ordering of the assumed scaled ceilings, NOT an assertion
    # of monotone measured rare-path probabilities or production convergence.
    ordered=[]
    for scale in ('small','medium','large'):
        pp=[p for p in points if p['scale']==scale and p['width']==8]
        ordered.append(max(F(8*(40*p['spacing'])**2,p['classes']) for p in pp))
    ck('scaled_detector_ceiling_order_only',ordered[0]>=ordered[1]>=ordered[2])
    ck('primitive179_pinned',z.instruction_primitives['fp_mul']==179)
    ck('primitive_normalize80',z.instruction_primitives['normalize52_pack']==2*5*3+2*5+2*5+12+9+9)
    ck('primitive_sub96',z.instruction_primitives['modsub52']==5+1+5+5+z.instruction_primitives['normalize52_pack'])
    ck('primitive_cmp5',z.instruction_primitives['cmp256']==4+1)
    ck('primitive_unpack15',z.instruction_primitives['unpack52']==5*3)
    ck('itemized34_unique',len(z.itemized_step_rows)==len({r['name'] for r in z.itemized_step_rows})==34)
    for row in z.itemized_step_rows:
        ck('cost_'+row['name'],isinstance(row['count'],int) and row['count']>=0 and z.step_rows[row['name']]==row['count']*z.instruction_primitives[row['primitive']])
    ck('itemized_sum2432_ceiling8000',sum(r['count']*z.instruction_primitives[r['primitive']] for r in z.itemized_step_rows)==2432<=z.E0)
    ck('probability_not_exponent_discount',z.success()==z.eta*(1-z.exp_minus_bounds(z.rho*z.T*(z.T-1)/(2*z.N))[1])-sum(z.losses.values(),F(0)))
    audit=retention_audit()
    checks.update(audit['checks'])
    evidence_paths=sorted((HERE/'evidence').glob('config-*.json'))+sorted((HERE/'evidence').glob('renewal-*.json'))+[HERE/'sweep-summary.json']
    values=dict(claim=z.frac(z.claim),raw=z.frac(z.raw),reserved=z.frac(z.W),
        success_conditional_lower=z.frac(z.success()),named_success_margin=z.frac(F(1,100)),
        rows={k:z.frac(v) for k,v in z.rows.items()},losses={k:z.frac(v) for k,v in z.losses.items()},
        integers={k:str(getattr(z,k)) for k in ('N','T','K','M','S','DP','R','L','workers')},
        eta=z.frac(z.eta),itemized_instructions=sum(z.step_rows.values()),instruction_primitives=z.instruction_primitives,
        finite_Hoeffding_scope='99% simultaneous for conditional selected mixtures of the fixed captured tables only; not a production confidence bound',
        nu=z.frac(z.nu),record_fields=z.record_fields,memory_rows={k:str(v) for k,v in z.memory_rows.items()},
        peak_memory_bytes=str(z.peak_memory),success_margin_over040=z.frac(z.success()-F(2,5)))
    return dict(version=1,verdict='PASS' if all(checks.values()) else 'HOLD',
        scope='Exact conditional fixed-spacing capacity, work and probability arithmetic plus captured toy reconciliation. H-Rom/H-rate/H-window production transfer remains an explicitly disclosed hypothesis, not an arithmetic theorem.',
        pass_count=sum(checks.values()),check_count=len(checks),checks=checks,
        failed_checks=[k for k,v in checks.items() if not v],
        production_qualification=dict(verdict='HOLD',reason='The observed path ratios pass, but detector ceilings are vacuous at toy sizes and measured rare-path rates do not establish the requested monotone production transfer.'),
        retention_audit={k:v for k,v in audit.items() if k!='checks'},values=values,
        evidence_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in evidence_paths},
        verifier_source=Path(__file__).read_text())


def main():
    out=build_certificate()
    data=json.dumps(out,separators=(',',':'))+'\n'
    log=f"{out['verdict']} {out['pass_count']}/{out['check_count']} conditional arithmetic/capture checks; claim {z.claim}\n"
    log+='PASS selected-pair tails, path ratios and 34-row instruction derivation\n' if not out['failed_checks'] else 'HOLD '+','.join(out['failed_checks'])+'\n'
    log+='HOLD production qualification: detector calibration vacuous; measured monotonicity not observed\n'
    if '--write' in sys.argv:
        (HERE/'certificate.json').write_text(data)
        (HERE/'evidence/certificate_run.log').write_text(log)
    else:
        assert (HERE/'certificate.json').read_text()==data,'stale certificate'
        assert (HERE/'evidence/certificate_run.log').read_text()==log,'stale run log'
    print(json.dumps({k:out[k] for k in ('verdict','pass_count','check_count','failed_checks')}))
    return 0 if out['verdict']=='PASS' else 2
if __name__=='__main__': sys.exit(main())
