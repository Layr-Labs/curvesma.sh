"""Full routing plus independent range-scanner tests; no production attack."""
import contextlib,io,json,runpy,sys
from pathlib import Path
base=Path(__file__).resolve().parent
base=base/'files' if (base/'files'/'experiment.py').is_file() else base
sys.path.insert(0,str(base))
def invoke(name,inputs):
    buf=io.StringIO();old=sys.stdin
    try:
        sys.stdin=io.StringIO(json.dumps(inputs))
        with contextlib.redirect_stdout(buf):
            runpy.run_path(str(base/name),run_name='__main__')
    finally:sys.stdin=old
    return json.loads(buf.getvalue())
inputs=json.load(sys.stdin)
seed=int(inputs.get('seed',1))
full=invoke('experiment.py',{'seed':seed,'trials':1500,'control_trials':500,'curve_trials':60,'seconds':12,'kernels':True})
isolated=invoke('research_scan.py',{'seed':seed})
checks={'cost_certificate':full['D_certificate']['passed'],
        'kernel_counts':full['counts_match_certificate'],
        'field_identities':all(full['A_identities'][k] for k in ('beta_cube_root_of_1_mod_p','lambda_cube_root_of_1_mod_r','lambda_G_equals_beta_x','window_forms_nonzero','seven_nonresidue_no_x0_point')),
        'class_rule':full['A_identities']['class_rule']['mismatches']==0,
        'field_kernel_correctness':full['A_kernels']['all_correct'],
        'arithmetic_sites':full['A_kernels']['site_checks']['all_ok'] and full['A_kernels']['sites_match_certificate'],
        'inversion':full['A_kernels']['fe_inv']['correct'],
        'routing_messages':full['B_routing']['messages_equal_emit_reference'],
        'routing_slots_tags':full['B_routing']['slots_equal_reference'] and full['B_routing']['source_tags_ok'],
        'routing_order':full['B_routing']['buckets_pointer_order_and_scan_ok'],
        'routing_candidates':full['B_routing']['candidates_equal_exact_order_reference'],
        'isolated_scanner':isolated['scanner_correct'],
        'isolated_cost_certificate':not isolated['certificate_failed_checks']}
print(json.dumps({'checks_ok':all(checks.values()),'checks':checks,'full_pipeline':full,'isolated':isolated}))
if not all(checks.values()):raise SystemExit(1)
