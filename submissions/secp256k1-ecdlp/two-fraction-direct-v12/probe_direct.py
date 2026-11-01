import json,sys,io,contextlib,random
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root));sys.path.insert(0,str(root/'kernels'))
cfg=json.load(sys.stdin);sys.stdin=io.StringIO(json.dumps(dict(seed=97,trials=400,curve_trials=10,seconds=15)))
import experiment as ex
cert=ex.load_certificate()
a=ex.part_a_kernels(random.Random(97),cert)
assert a["all_correct"] and a["keys_formula_holds"] and a["fwd_bwd_formulas_hold"] and a["sites_match_certificate"]
assert a["fixup_batch"]["ok"] and a["boundary_batch"]["ok"] and a["forced_fixups"]["ok"] and a["site_checks"]["fixups_ok"] and a["fixups"]["flags_dead_at_every_join"]
import routing_check
b=routing_check.check()
b_wide=routing_check.check(cap1=17,cap3=9,q2=8)
assert all(b_wide[k] for k in ["formulas_ok","messages_equal_emit_reference","slots_equal_reference","source_tags_ok","buckets_pointer_order_and_scan_ok","candidates_equal_exact_order_reference"])
assert b['formulas_ok'] and b['messages_equal_emit_reference'] and b['slots_equal_reference'] and b['source_tags_ok'] and b['buckets_pointer_order_and_scan_ok'] and b['candidates_equal_exact_order_reference'],b
with contextlib.redirect_stdout(io.StringIO()):import certificate as m
q=145289002;v=m.success_terms(q);cost=m.ledger(q)
assert v['lo']>=m.THRESH and m.success_terms(q-1)['lo']<m.THRESH
ids=ex.part_a_identities(random.Random(101))
assert ids['class_labels']['ok']
assert all(v for v in ids.values() if isinstance(v,bool))
pr,la=ex.a64.a64_load((root/'kernels/route.s').read_text())
scan=ex.scan_paths(pr,la,random.Random(103),cert);assert scan['ok']
toy=ex.toy.part_c(random.Random(107),ex.deadline,400,10,400,10)
assert not m.FAILED,m.FAILED
assert toy['curve_pipeline_vs_logs']['identical_outcomes']==toy['curve_pipeline_vs_logs']['trials']
s=m.SIDE;u=m.Q2;bcap=v['cap_3'];F=(2**20+5*2048+1024+288*q+160*u+32*(q+u+1)+256*u+2*64+8*s+32*s+16*2**24+8*2**48+16*2**28+64*2**29+16*2**28*64+64*2**34+96*s+2048+2**34)
assert F<2**51+2**44<2**52
end=2**52+32*m.CAP_P*s+40*2**24*bcap
attempt=2**52+32*m.CAP_P*(s-1)+32*6*((q+1)//2)*u
assert end<2**64 and attempt<2**63
assert 12*s+88*8*m.REQ_CAP<=m.EXCH
assert m.pow2_interval(m.Fraction(125142304476874,10**12))[0]>=m.Fraction(cost['total'],1000)
assert cost['preprocessing']<1000*2**111
print(json.dumps(dict(research_only=False,arithmetic=a,routing=b,routing_wide=b_wide,identities=ids,scanner=scan,toy=toy,memory_fixed=F,memory_end=end,attempted_cursor_bound=attempt,conditional_log2=cost['log2_display'],declared_score="125.142304476874",success_lower_bound=float(v["lo"]),q1=q,total=str(cost['total']),scope='Complete finite arithmetic/routing and exact consequences of the supplied analytical proof; not a production execution or independent mathematical certificate.')))
