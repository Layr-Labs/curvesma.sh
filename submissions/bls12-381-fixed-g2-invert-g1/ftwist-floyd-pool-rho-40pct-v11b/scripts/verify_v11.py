#!/usr/bin/env python3
"""Exact fail-closed v11 verification. Exit 2 means HOLD, not release.
All arithmetic is integer/Fraction; toy captures cannot prove production H-Rom.
"""
import argparse,ast,hashlib,io,json,sys,unittest,zlib
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from ledger import *
from base_audit import CONFIGS,DENSITIES
from tower import Tower
from cap_analysis import evaluate
ROOT=Path(__file__).resolve().parents[1]


def frac(x):
    x=F(x);return dict(numerator=str(x.numerator),denominator=str(x.denominator))


def rf(x):return F(int(x['numerator']),int(x['denominator']))


def under(x,y):
    y=F(y);return below_pow2(x,y.numerator,y.denominator)


def exhaustive():
    checked=blind=0;ok=True
    for n in range(1,5):
        for dest in product(range(n),repeat=n):
            for start in range(n):
                seen={};path=[];u=start
                while u not in seen:
                    seen[u]=len(path);path.append(u);u=dest[u]
                mu=seen[u];cyc=path[mu:];ell=len(cyc);t=h=start;tick=0
                while True:
                    tick+=1;t=dest[t];h=dest[dest[h]]
                    if t==h:break
                    if tick>mu+ell:ok=False;break
                ok &= tick<=mu+ell
                for mask in range(1<<n):
                    checked+=1
                    if not any(mask&(1<<i) for i in cyc):
                        blind+=1;ok &= t==h
    return dict(function_start_DP_partitions=checked,DP_free_partitions=blind,
        all_caught=bool(ok),scope='All functional maps on at most four states; analytic mu+lambda proof covers arbitrary finite periods.')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true')
    ap.add_argument('--raw-dir',type=Path);a=ap.parse_args();ck={}
    ledger_source=(ROOT/'scripts/ledger.py').read_text()
    def check(k,v):ck[k]='PASS' if bool(v) else 'FAIL'
    check('protected_C20000',C==20000)
    check('class_count',12*N==R-1)
    check('aut_order',pow(P,12,R)==1 and pow(P,6,R)==R-1 and pow(P,4,R)!=1)
    check('cyclotomic',(1+pow(P,4,R)+pow(P,8,R))%R==0)
    check('Np_directed_ceiling',(NP-1)*1023<N*1024<=NP*1023)
    check('ln_rational_enclosure',F(51082562376599,10**14)<LN_LO<LN_HI<F(51082562376600,10**14))
    check('ln_width',LN_HI-LN_LO<F(1,1<<300))
    check('T0_unique',(T0-1)**2<2*NP*LN_LO<2*NP*LN_HI<=T0*T0)
    check('T_reserve',100*(T-1)<109*T0<=100*T)
    check('ordinary_sum5602',sum(ROWS.values())==5602<=ORDINARY_CEILING==5700)
    check('Floyd256_added',ROWS['Floyd_key_comparison_control']==256)
    check('expanded_compare132_added',ROWS['brent_expanded_state_compare_control']==132)
    check('checkpoint76_added',ROWS['checkpoint_full_copy']==76)
    check('pool256_retained',ROWS['DP_pool_key_comparison_control']==256)
    check('no_orientation_guard',not any('guard' in k for k in ROWS))
    check('escape_rows',sum(ESCAPE_ROWS.values())==53290<=ESCAPE_INSTRUCTIONS==55000)
    check('single_paid_inplace_add',ESCAPE_REF==1)
    check('no_fresh_start_row',not any('restart' in k or 'fresh_start' in k for k in ESCAPE_ROWS))
    check('opportunity_prediction',ESCAPE_PRED==F(3*T,12*BRANCHES))
    check('four_times_mean_is_conditional_ceiling',ESCAPE_MEAN_CEILING==4*ESCAPE_PRED)
    check('global_allowance_ceiling',8*(ESCAPE_CAP-1)<T<=8*ESCAPE_CAP)
    check('global_allowance_tail',ESCAPE_TAIL==ESCAPE_MEAN_CEILING/ESCAPE_CAP<=F(1,128))
    check('no_per_trail_escape_cap','ESCAPES_PER_TRAIL' not in ledger_source)
    check('charged_escape_tail_is_global',STRICT_ESCAPE_TAIL==ESCAPE_TAIL==ESCAPE_MEAN_CEILING/ESCAPE_CAP)
    check('e4_only_an_ablation',E4_ABLATION_ESCAPE_CAP==4 and 'E4_ABLATION_ESCAPE_CAP' in ledger_source)
    check('execution_envelope',WEXEC==F(771,200)*T+ESCAPE_BUDGET_EXEC+(1<<100))
    check('repriced_all_nonreference',WSTRESS==F(2013,500)*T+ESCAPE_BUDGET_STRESS+(1<<101))
    check('escape_budget_prices_allowance',
        ESCAPE_BUDGET_EXEC==(ESCAPE_REF+F(ESCAPE_INSTRUCTIONS,C))*ESCAPE_CAP and
        ESCAPE_BUDGET_STRESS==(ESCAPE_REF+F(6,5)*F(ESCAPE_INSTRUCTIONS,C))*ESCAPE_CAP)
    check('allowance_covers_mean_ceiling',ESCAPE_CAP>=ESCAPE_MEAN_CEILING)
    check('repriced_below_execution_envelope',WEXEC<WSTRESS<W)
    check('scored_envelope',W==F(29,5)*T+(1<<100))
    check('reduced_claim',CLAIM==F(2567,20))
    check('numeric_ship_threshold',CLAIM<F(12881,100))
    check('claim_window_low',under(W,CLAIM-F(1,200)))
    check('claim_window_high',not under(W,CLAIM-F(3,50)))
    check('repriced_headroom010',under(WSTRESS,CLAIM-F(1,10)))
    check('scored_interval',interval(W)==[1283126,1283127,10000])
    check('repriced_interval',interval(WSTRESS)==[1279666,1279667,10000])
    check('execution_interval',interval(WEXEC)==[1278888,1278889,10000])
    check('fixed_spacing_capacity',DP_SPACING==1<<54 and POOL_CAPACITY==1<<80)
    check('all_report_opportunities',K==3*T+ESCAPE_CAP+1)
    check('same_spacing_pool_Markov',POOL_MARKOV==BALANCE*F(K,DP_SPACING*POOL_CAPACITY))
    check('same_spacing_report_Markov',REPORT_MARKOV==BALANCE*F(K,DP_SPACING*REPORT_CAP))
    check('pool_debit_below0012',POOL_MARKOV<F(12,1000))
    check('transport_exact_arithmetic_only',TRANSPORT==F(64*DP_SPACING**2,N))
    check('full_entropy_layout',24+32+32+24+16==RECORD_BYTES==128)
    check('wide_address_links',96>80 and 192*POOL_CAPACITY+(1<<40)<1<<88)
    check('both_routing_endpoints',2*192*(1<<14)+(1<<20)<1<<26)
    check('independent_fullpoint_fingerprint',FINGERPRINT_BITS==192 and 192*4608//8<1<<17)
    check('false_candidate_cap',FINGERPRINT_TAIL==F(REPORT_CAP*POOL_CAPACITY,(1<<192)*(1<<32))<=F(1,1<<60))
    ar,at=ancillary()
    check('report_work',ar['report_instructions']==REPORT_CAP*(1<<26))
    check('initialization_work',ar['clear_instructions']==POOL_CAPACITY*192*16)
    check('all_ancillary_charged',at<1<<100 and F(6,5)*at<1<<101)
    check('two_jump_tables_setup',2*1024*(2*510*23000+12*(72*383+10000))+(1<<32)<(1<<50)*C)
    check('two_jump_tables_storage',2*1024*12*(576+64)<1<<24)
    check('field_exception_unions',7*T*(1<<119)<N and 7*T*1024*(1<<239)<P)
    check('named_margin',F(2,5)-SUCCESS_FLOOR==NAMED_MODEL_MARGIN==F(1,100))
    pre=CONDITIONAL_MODEL_REFERENCE-ESCAPE_TAIL-POOL_MARKOV-REPORT_MARKOV-TRANSPORT-FINGERPRINT_TAIL
    check('reference_model_precap_margin_only',pre>F(2,5))
    check('correct_postcap_model_arithmetic',MODEL_LOWER==CONDITIONAL_MODEL_REFERENCE-STRICT_ESCAPE_TAIL-POOL_MARKOV-REPORT_MARKOV-TRANSPORT-FINGERPRINT_TAIL)
    check('tail_row_is_global_allowance_not_clipped_mean',
        STRICT_ESCAPE_TAIL==ESCAPE_TAIL and STRICT_ESCAPE_TAIL!=1 and 'min(F(1)' not in ledger_source)
    check('floor_from_global_tail_only',MODEL_LOWER>=F(39,100) and MODEL_LOWER>=F(2,5))
    suite=unittest.defaultTestLoader.discover(str(ROOT/'scripts'),pattern='test_mechanism.py')
    result=unittest.TextTestRunner(stream=io.StringIO()).run(suite)
    check('seven_mechanism_regressions',result.wasSuccessful() and result.testsRun==7)
    source=(ROOT/'scripts/retention.py').read_text();tree=ast.parse(source)
    check('zero_delete_AST',not any(isinstance(n,ast.Delete) for n in ast.walk(tree)))
    check('no_pool_pop_clear',not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ('pop','popitem','clear','update') for n in ast.walk(tree)))
    check('single_immutable_insert',source.count('self.records[key]=record')==1 and 'if key in self.records:return self.records[key]' in source)
    esc=source.split('    def escape_event(self):',1)[1].split('\n\ndef run_walk',1)[0]
    check('escape_cannot_erase_records','self.records' not in esc)
    check('full_pool_aborts','pool_failure=1;return True' in source and 'self.full=True' in source)
    check('Floyd_independent_of_DP','if tort[0]!=hare[0]:continue' in source and 'hare=advance(hare)' in source)
    check('detector_reset_not_point_reset','hare=tort  # detector reset' in source and 'tort=escape(tort)' in source)
    check('per_trail_escapes_unbounded_charged','No per-trail escape cap' in source and
        'if pool.escapes>=escape_allowance' in source and 'escape_R_adds=pool.escapes' in source)
    proof=exhaustive()
    check('all_period4_DP_partitions',proof['all_caught'] and proof['function_start_DP_partitions']==17066)
    check('DP_free_partitions_nonvacuous',proof['DP_free_partitions']>1000)
    toy=json.loads((ROOT/'toy_results.json').read_text());configs=[];parts=[];total=0
    check('seven_configurations',len(toy['configurations'])==7)
    check('2500_trials_each',toy['trials_per_config']==2500 and toy['seed']==2026102310)
    check('runtime_under25minutes',0<toy['elapsed_nanoseconds']<1500*10**9)
    compact=[]
    for i,item in enumerate(toy['configurations']):
        text=''.join((ROOT/p).read_text() for p in item['parts']);parts+=item['parts']
        c=json.loads(text);configs.append(c);tag=f'c{i}_';tr=c['trials'];B=c['branches'];total+=tr
        check(tag+'lossless_concat_SHA',hashlib.sha256(text.encode()).hexdigest()==item['summary_sha256'])
        check(tag+'params',(c['q'],c['r'],B)==CONFIGS[i] and tr==2500 and c['seed']==toy['seed']+i)
        check(tag+'exhaustive_live',min(c['exhaustive_branch_histogram'])>0 and sum(c['exhaustive_branch_histogram'])==c['classes'])
        check(tag+'literal_tower',c['literal_tower_crosschecks']==24 and c['literal_raw_twist_checks']==32)
        check(tag+'fixed_spacings',c['fixed_DP_bits']==[2,4,6])
        check(tag+'three_eval_tick',c['ordinary_cap']==3*c['floyd_cap'])
        check(tag+'production_no_per_trail_cap',c['e4_ablation_escape_cap']==4 and c['escape_allowance']>4)
        check(tag+'max_trail_escapes_within_allowance',0<c['max_trail_escapes']<=c['escape_allowance'])
        check(tag+'all_branches_visited',min(c['live_branch_histogram'])>0)
        check(tag+'all_class_pool',c['pool_capacity']==c['classes']+1)
        check(tag+'label_checks',c['label_invariant_checks']>0)
        check(tag+'shared_fresh',c['variants']['shared']['trials']==c['variants']['fresh']['trials']==tr)
        check(tag+'fixed_correlated',c['fixed_table_control']['starts']==256 and 'Correlated' in c['fixed_table_control']['scope'])
        f=Tower(c['q']);g=tuple(c['generator']);r=c['r'];shift=B.bit_length()-1
        for d in DENSITIES:
            p=c['pipelines'][str(d)];dt=tag+f'd{d}_'
            check(dt+'fraction',rf(p['success_fraction'])==F(p['success'],tr))
            check(dt+'failed_once',p['success']+p['capped_failures']==tr)
            check(dt+'zero_delete',p['deletions']==0)
            check(dt+'zero_escape_erase',p['escape_retention_failures']==0)
            check(dt+'append_only',p['insertions']==p['pool_final'])
            check(dt+'one_jump_charge',p['escape_R_adds']==p['escapes'])
            check(dt+'hard_caps',p['ordinary_steps']<=tr*c['ordinary_cap'] and p['escapes']<=tr*c['escape_allowance'])
            check(dt+'no_toy_overflow',p['pool_cap_failure']==p['report_cap_failure']==0)
            check(dt+'escape_frequency_exact',rf(p['escape_rate'])==F(p['escapes'],max(1,p['ordinary_steps'])) and rf(p['prediction'])==F(1,12*B))
            check(dt+'paired_v10',p['success']-p['v10_success']==p['new_not_v10']-p['v10_not_new'])
            check(dt+'paired_old',p['success']-c['old_floyd']['success']==p['new_not_old_floyd']-p['old_floyd_not_new'])
            check(dt+'all_reachable_blind_caught',p['unrecovered_initial_blind_caught']==p['unrecovered_initial_blind_reachable'])
            check(dt+'blind_events_nonvacuous',p['fruitless_DP_free_caught']>0)
            check(dt+'finite_frequency_below_fourmean',F(p['escapes'],max(1,p['ordinary_steps']))<=F(4,12*B))
            if B==1024:check(dt+'finite_margin040',F(p['success'],tr)-F(3,50)>=F(2,5))
            w=c['no_DP_fruitless_witnesses'][str(d)]
            check(dt+'blind_witness',w['states'][0]==w['states'][1] and all(((key>>shift)&((1<<d)-1))!=0 for key in w['keys']))
            st=tuple(w['states'][0]);literal_ok=label_ok=True
            for e,key in zip(w['edges'],w['keys']):
                u=f.pow(g,e['source_rep']);v=f.mul(u,f.pow(g,e['jump']));_,k,rep=f.metadata(v)
                literal_ok &= f.keys(u)[0]==key and rep==f.pow(g,e['destination_rep']) and pow(c['q'],k,r)==e['alignment']
                label_ok &= st[0]==e['source'] and (st[1]+w['target_x']*st[2])%r==e['source_rep']
                m=e['alignment'];st=(e['destination'],m*(st[1]+e['jump']-w['target_x']*e['beta'])%r,m*(st[2]+e['beta'])%r)
            check(dt+'literal_edges',literal_ok)
            check(dt+'closed_label_lap',label_ok and st==tuple(w['states'][0]))
            due=w['floyd_first_meeting']<=c['floyd_cap'];direct=w['direct_probe']
            check(dt+'direct_DPfree_probe_correct_bound',
                (due and direct['fruitless_meetings']>0 and direct['escapes']==1) or
                (not due and direct['floyd_meetings']==direct['escapes']==0 and direct['ordinary_steps']==c['ordinary_cap']))
        if a.raw_dir:
            raw=zlib.decompress((a.raw_dir/f'paired-{i}.json.zlib').read_bytes());rows=json.loads(raw)
            check(tag+'raw_SHA',hashlib.sha256(raw).hexdigest()==c['paired_raw_sha256'])
            check(tag+'raw_size',len(rows)==tr and all(len(row)==25 for row in rows))
            check(tag+'raw_old',sum(row[0] for row in rows)==c['old_floyd']['success'])
            for j,d in enumerate(DENSITIES):
                p=c['pipelines'][str(d)];b=1+8*j
                cols=('success','escapes','cross_escape_recovery','reports','allowance_failure','v10_success','e4_ablation_success','fruitless_DP_free_caught')
                check(tag+f'd{d}_raw_sums',all(sum(row[b+k] for row in rows)==p[name] for k,name in enumerate(cols)))
        compact.append(dict(N=c['classes'],B=B,old=c['old_floyd']['success'],v10=[c['pipelines'][str(d)]['v10_success'] for d in DENSITIES],v11=[c['pipelines'][str(d)]['success'] for d in DENSITIES]))
    check('17500_independent_sets',total==17500)
    check('split_inventory',set(parts)=={str(p.relative_to(ROOT)) for p in (ROOT/'evidence').glob('config-*.txt')})
    stress=configs[4];p=stress['pipelines']['2']
    check('B8_same_seed_controls',stress['branches']==8 and stress['old_floyd']['success']==613 and p['v10_success']==512)
    check('B8_production_recovery_above40pct',p['success']==1453 and F(p['success'],2500)>=F(2,5))
    check('B8_blind_catch_100pct',p['unrecovered_initial_blind_reachable']==p['unrecovered_initial_blind_caught']==554)
    check('B8_no_DP_events',p['fruitless_DP_free_caught']==3721<=p['fruitless_meetings']==6573)
    check('B8_cross_escape_credit',p['cross_escape_recovery']==778)
    check('B8_escapes_exceed_four_charged',stress['max_trail_escapes']==10 and p['trails_above_four_escapes']==358 and p['escape_R_adds']==p['escapes'])
    check('B8_no_allowance_failure',p['allowance_failure']==0)
    check('B8_e4_ablation_disclosed',p['e4_ablation_success']==1331 and p['e4_ablation_allowance_failures']==358)
    check('B8_cap_removal_gain_disclosed',p['success']-p['e4_ablation_success']==122)
    scale=evaluate()
    check('diagnostic_scope','NOT a measured production distribution' in scale['renewal_diagnostic']['scope'])
    check('diagnostic_low_union',rf(scale['renewal_diagnostic']['success_upper_union'])<SUCCESS_FLOOR)
    check('Floyd_only_not_requalified',not scale['floyd_only_preserving_emitted_records']['qualified'])
    # Both production qualification gates now pass on the honest charged-
    # escape accounting: 1/128 global tail and margin above 2/5.
    check('production_global_tail_qualified',scale['qualification']['global_allowance_tail_nonvacuous'])
    check('production_success_margin040',scale['qualification']['production_success_margin'])
    sub=json.loads((ROOT/'submission.json').read_text());note=(ROOT/'public-note.md').read_text();arg=(ROOT/'argument.md').read_text()
    check('manifest_exact',set(sub)=={'version','title','claims','experiments'} and sub['version']==1)
    check('one_claim_exact',len(sub['claims'])==1 and set(sub['claims'][0])=={'evidence','score','statement'})
    check('reduced_score_strings',sub['claims'][0]['score']==frac(CLAIM))
    check('empty_experiments_deps',sub['experiments']==[] and json.loads((ROOT/'dependencies.json').read_text())==[])
    check('self_contained_statement',all(s not in sub['claims'][0]['statement'] for s in ('19c19d8c','mpjunior92','improves','128.81')))
    check('release_statement_clean',all(s not in sub['claims'][0]['statement'] for s in ('HOLD','unqualified','vacuous','not an achieved bound')))
    check('documents_not_held','HOLD' not in note and 'HOLD' not in arg and 'not submitted' not in note)
    check('note_10to14KiB',10240<=len(note.encode())<14336)
    check('argument_under8000',len(arg.encode())<8000)
    check('actual_design_attribution','GPT 6.1 Sol / Hermes' in note and 'qwen3.8-flash-next / Hermes' in note)
    check('public_prior_lineage',all(s in note for s in ('mpjunior92','ftwist-escape-rho-40pct-v8','ftwist-rare-escape-rho-40pct-v9','ftwist-retain-rho-40pct-v10')))
    excluded={'certificate.json','certificate-details.json','certificate_run.log','intake-audit.json'}
    paths=sorted(p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    source_paths=[p for p in paths if p.name not in excluded]
    check('archive_filecount32',len(source_paths)+3<=32)
    check('evidence_parts_under16KiB',all(p.stat().st_size<=16384 for p in source_paths if 'evidence' in p.parts))
    check('source_files_under64KiB',all(p.stat().st_size<=65536 for p in source_paths))
    check('private_path_clean',all((chr(47)+'Users'+chr(47)) not in p.read_text() for p in source_paths))
    check('secret_shape_clean',all(('y'+'kn'+'_') not in p.read_text() for p in source_paths))
    check('no_float_literals',not any(isinstance(n,ast.Constant) and isinstance(n.value,float) for p in (ROOT/'scripts').glob('*.py') for n in ast.walk(ast.parse(p.read_text()))))
    bindings={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    details=dict(version=1,status='HOLD',scope='Exact arithmetic and captured finite execution only; escape accounting is charged-per-escape with the priced global first-moment allowance.',
        constants=dict(N=N,Np=NP,T0=T0,T=T,per_trail_escape_cap=None,global_priced_allowance=ESCAPE_CAP,escape_instructions=ESCAPE_INSTRUCTIONS,K=K,claim=frac(CLAIM),DP_spacing=DP_SPACING,pool_capacity=POOL_CAPACITY,record_bytes=RECORD_BYTES,storage_bytes=192*POOL_CAPACITY+(1<<40)),
        work=dict(execution_envelope=frac(WEXEC),repriced=frac(WSTRESS),scored=frac(W)),
        log2_intervals=dict(execution=interval(WEXEC),repriced=interval(WSTRESS),scored=interval(W)),
        losses=dict(global_escape_tail=frac(STRICT_ESCAPE_TAIL),pool=frac(POOL_MARKOV),report=frac(REPORT_MARKOV),transport=frac(TRANSPORT),fingerprint=frac(FINGERPRINT_TAIL)),
        finite_results=compact,exhaustive_Floyd=proof,cap_analysis=scale,
        raw_captures_checked=a.raw_dir is not None,source_sha256=bindings)
    dt=json.dumps(details,separators=(',',':'))+'\n'
    check('details_under16KiB',len(dt.encode())<=16384)
    check('archive_byte_budget',sum(p.stat().st_size for p in source_paths)+len(dt.encode())+50000<262144)
    failed=[k for k,v in ck.items() if v!='PASS'];status='HOLD' if failed else 'PASS'
    if details['status']!=status:
        details['status']=status;dt=json.dumps(details,separators=(',',':'))+'\n'
    out=dict(version=1,status=status,checks=dict(ck),all_checks_pass=not failed,
        details='evidence/certificate-details.json',details_sha256=hashlib.sha256(dt.encode()).hexdigest(),
        verifier_source=(ROOT/'scripts/verify_v11.py').read_text())
    text=json.dumps(out,separators=(',',':'))+'\n'
    assert len(text.encode())<=50000
    lines=[json.dumps(dict(status=status,checks=len(ck),passed=len(ck)-len(failed),failed=failed,claim=frac(CLAIM)),separators=(',',':'))]
    if a.write:
        (ROOT/'evidence/certificate-details.json').write_text(dt)
        (ROOT/'certificate.json').write_text(text)
        (ROOT/'evidence/certificate_run.log').write_text('\n'.join(lines)+'\n')
    else:
        assert (ROOT/'evidence/certificate-details.json').read_text()==dt
        assert json.loads((ROOT/'certificate.json').read_text())==out
    print('\n'.join(lines),flush=True)
    if failed:raise SystemExit(2)
if __name__=='__main__':main()
