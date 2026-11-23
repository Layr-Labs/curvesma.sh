"""Current-package checks. Execution belongs in llm-verifier's isolated runner.

The preserved reference certificate is recomputed once. Current sources are
loaded explicitly; the historical manifest and runtime count declarations are
not substituted for current checks. Toy executions establish only their scope.
"""
import importlib.util
import io
import json
import sys
import time
from pathlib import Path

START=time.time()
ROOT=Path('/snapshot/files') if Path('/snapshot/files/base').is_dir() else Path(__file__).parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'kernels'))
import a64
import costs

def load_reference(cfg):
    old_stdin=sys.stdin
    sys.stdin=io.StringIO(json.dumps(cfg))
    name='reference_experiment'
    spec=importlib.util.spec_from_file_location(name,ROOT/'base/experiment.py')
    E=importlib.util.module_from_spec(spec);sys.modules[name]=E
    try:spec.loader.exec_module(E)
    finally:sys.stdin=old_stdin
    E.CE=costs.c
    E.CLS=costs.BASELINE_TRANSCRIPT.splitlines()
    # No second certificate import: all original checks have already run.
    E.lcert=lambda:(costs.c,E.CLS)
    sys.path.insert(0,str(ROOT))
    return E

def full_y(E,K,sources,F):
    program,labels=a64.a64_load(sources['leaf']+'\n'+sources['secp_fe']+'\n'+sources['ytransport'])
    N=F['memory'];initial=dict(N)
    YS=(34<<40)+0x7100000;DONE=(34<<40)+0x7120000;SP=(34<<40)+0x7140000
    def execute(pointer_change=False,link_change=False):
        N.protect=False;N.clear();N.update(initial)
        N.coins=iter(E.limbs(F['rho']));N.view_writes=[];N.view_reads=set()
        if pointer_change:N[F['sb']+176]+=512
        if link_change:
            y=F['ys'][0][0]
            first=E.offs7(y,0)[0]
            N[first+120]=first+128
        N[YS+48]=0x123456789;N[YS+56]=0
        for a,name in zip((64,72,80,88),('fwd','fe_inv','bwd','keys')):N[YS+a]=labels[name][0]
        E.put(N,YS+96,E.r);N[DONE]=0;N.protect=True
        SN=E.machine.gS1(N);port=E.machine.Port1(SN);regs={}
        LS=([(E.AXYX0,E.AXY_END),(E.XT0,1<<64),(F['view'],F['view']+512*F['n']),
             (F['sb'],E.EBE)],[(F['sb'],E.EBE)])
        count=E.arun(program,labels,'yconsume',dict(x0=F['view'],x1=F['n'],x2=DONE,x3=9,sp=SP),
            SN,out_regs=regs,dev=port.dev,ls64=LS)
        errors=sum(N[a]!=v for a,v in F['expected_eb'].items())
        reference_sender=E.machine.Sender1(K)
        E.machine.kref(E,F['reference'],F['sb'],6*F['n'],0x123456789,5<<26,reference_sender,K)
        return dict(count=count,eb_mismatches=errors,expected_eb_words=len(F['expected_eb']),
            pointers_kept=all(N[a]==v for a,v in F['pointers'].items()),
            read_only=not N.view_writes,done=N[DONE],stack_restored=regs['sp']==SP,
            coins_ok=N[YS+56]==0,records_match=port.blocks==reference_sender.msgs,
            sender_state=E.machine.state1_ok(SN,reference_sender))
    good=execute()
    assert good['eb_mismatches']==0 and good['done']==9
    assert all(good[k] for k in ('pointers_kept','read_only','stack_restored','coins_ok','records_match','sender_state'))
    wrong_pointer=execute(pointer_change=True)
    wrong_link=execute(link_change=True)
    assert wrong_pointer['eb_mismatches']>0 and wrong_link['eb_mismatches']>0
    return dict(ok=True,positive=good,wrong_pointer_output_errors=wrong_pointer['eb_mismatches'],
        wrong_link_output_errors=wrong_link['eb_mismatches'])

def main(cfg):
    import leaf_check,map_check,ydriver_check,fault_check
    E=load_reference(cfg);K=costs.c
    sources={name:(ROOT/'kernels'/f'{name}.s').read_text()
             for name in ('leaf','secp_fe','route','ytransport','pointer_init')}
    programs={name:a64.a64_load(src) for name,src in sources.items()}
    counts={name:len(pair[0]) for name,pair in programs.items()}
    faults=fault_check.run(programs['leaf'],programs['route'])
    bad_ops=list(programs['route'][0]);bad_labels=programs['route'][1]
    bad_ops[bad_labels['krq6a'][0]+0x20000//4]=('b',['rsx24'])
    try:fault_check.run(programs['leaf'],(bad_ops,bad_labels))
    except AssertionError as exc:bad_shadow='krq6a_shadow' in str(exc)
    else:bad_shadow=False
    assert bad_shadow
    faults['wrong_shadow_control']=bad_shadow
    out={'scope':'Exact cost recomputation plus finite current-source checks, not a production attack run.',
         'costs':costs.summary(),'program_instructions':counts,'fault_wiring':faults}
    mapping=map_check.audit(counts)
    assert mapping['relay']['cap']==costs.RELAY_CAP
    assert mapping['relay']['tables_charged']==costs.TABLES_R
    assert mapping['program']['allowance']==costs.PROGRAM_BYTES
    peak=sum(mapping['physical'][k]['pa_end']*v for k,v in {'G':K.NODES,'R':K.N_REL,'H':K.N_HOME}.items())
    assert peak<=2**138
    out['memory_peak_reserved_bytes']=str(peak)
    out['map']={k:mapping[k] for k in ('all_checks_pass','checks','negative_controls','relay','clone','yview','program')}
    out['map']['pa_ends']={k:v['pa_end'] for k,v in mapping['physical'].items()}
    out['leaf'],fixture=leaf_check.leaf_check(E,K,sources['leaf'],sources['pointer_init'],n=4)
    out['full_y']=full_y(E,K,sources,fixture)
    out['y_driver']=ydriver_check.run_checks()
    import relay_checks,home_checks
    original=(ROOT/'base/kernels/route.s').read_text()
    out['relay']=relay_checks.run(E,sources['route'],original,mode='compact')
    out['home']=home_checks.run(E,sources['route'],original,mode='compact')
    out['cold']=leaf_check.run_cold(E,K,sources['secp_fe'])
    relay_surplus=K.FLUSH+K.EXC_ENTRY+K.DONE_C1+K.H2S_FULL-43
    assert relay_surplus==12 and relay_surplus>=0
    out['relay_partial_surplus_per128']=relay_surplus
    partial=max(23-(o>0)+K.dsg(0,o)+K.dsg(o,n)+K.dsg(o+n,256-o-n)
                for o in range(257) for n in range(257-o) if (o,n)!=(0,256))
    assert partial==costs.DLVHS==811
    out['inherited_partial_copy_envelope']=partial
    assert out['home']['partial_cost']['charged_boundary']==costs.FINALH+2048
    assert out['home']['partial_cost']['max_nonhash']<=1200
    out['all_checks_pass']=True
    out['seconds']=round(time.time()-START,3)
    print(json.dumps(out,separators=(',',':')))

if __name__=='__main__':main(json.load(sys.stdin))
