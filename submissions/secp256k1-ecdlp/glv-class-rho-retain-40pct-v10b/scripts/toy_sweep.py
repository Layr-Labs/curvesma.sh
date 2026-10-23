"""Deterministic three-scale, binding-envelope toy capture; stdlib only.
python3 -B scripts/toy_sweep.py [--quick]
Actual GLV coordinates choose branch and DP; exponent oracle only adds points.
Bare-map first repeats and full ring-escape DP recovery are separate statistics.
Neither finite toy evidence nor author arithmetic establishes production mixing.
"""
import hashlib
import json
import math
import random
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import ledger as z
from toy_group import scale_groups
from scripts.transports import delta_transport


def binding_budget(N,d):
    delta=delta_transport(d,classes=N)
    assert delta < F(3,10)
    T=16
    while True:
        _,hi=z.exp_minus_bounds(F(7,10)*F(T*(T-1),2*N))
        exact_lower=1-hi-delta
        # Directed decimal rational keeps captures compact; never round up.
        envelope=F(exact_lower.numerator*10**12//exact_lower.denominator,10**12)
        if envelope>F(1,50): return T,envelope
        T+=16


def recover(q,a1,b1,a2,b2):
    den=(a2-a1)%q
    return (b1-b2)*pow(den,-1,q)%q if den else None


def one(g,width,seed,d,trials,T,shift,trans):
    q,N=g.Q,g.N
    jumps,branch,dest,unit,hist=trans
    rng=random.Random(seed*1000003+d*73)
    kl,xs=g.kl,g.xs
    dp=[(x & ((1<<d)-1))==0 for x in xs]
    horizon=40*(1<<d)
    total=Counter(); repeats=[]; cycles=[]
    for _ in range(trials):
        secret=rng.randrange(1,q)
        ja=[rng.randrange(q) for _ in jumps]
        jb=[(k-a*secret)%q for k,a in zip(jumps,ja)]
        # Uninterrupted bare class map: exact first repeat, not a censored mean.
        x=rng.randrange(N); seen={}; t=0
        while x>=0 and x not in seen:
            seen[x]=t; x=dest[x]; t+=1
        if x>=0: repeats.append(t); cycles.append(t-seen[x])
        else: total['bare_identity_trials']+=1
        used=0; serial=0; good=False; any_recovery=False; core_repeat=False
        first={}; pool={}; st=Counter()
        def start():
            nonlocal serial
            serial+=1; st['starts']+=1
            x=rng.randrange(N); a=rng.randrange(q)
            return [x,a,(kl[x]-a*secret)%q,[],None,1,0,0,set(),serial]
        def visit(v):
            nonlocal core_repeat
            if used<T:
                x=v[0]; v[8].add(x)
                if x in first and first[x]!=v[9]: core_repeat=True
                first.setdefault(x,v[9])
        def endpoint(v):
            nonlocal good,any_recovery
            x,a,b=v[:3]; st['dp_events']+=1
            # Exact endpoints; production fingerprint errors are a separate debit.
            for oa,ob,oc in pool.get(x,[]):
                candidate=recover(q,oa,ob,a,b)
                if candidate is None: st['degenerate_relation_pairs']+=1; continue
                if candidate!=secret or g.points[candidate]!=g.points[secret]:
                    st['wrong_recoveries']+=1
                    raise AssertionError('wrong recovered logarithm')
                any_recovery=True
                if v[8].intersection(oc): good=True
            pool.setdefault(x,[]).append((a,b,set(v[8])))
        def tick(v):
            nonlocal used,any_recovery
            if v[7]>=horizon: st['length_aborts']+=1; return False
            visit(v)
            x=v[0]; j=branch[x]; u=unit[x]; xn=dest[x]
            if xn<0: u=1
            v[1]=u*(v[1]+ja[j])%q; v[2]=u*(v[2]+jb[j])%q
            used+=1; v[7]+=1; st['dp_tests']+=1
            if xn<0:
                st['identity_events']+=1
                if v[1]:
                    candidate=(-v[2]*pow(v[1],-1,q))%q
                    assert candidate==secret
                    any_recovery=True
                return False
            v[0]=xn
            if dp[xn]: endpoint(v); return False
            ring=v[3]
            if xn in ring:
                cycle=ring[ring.index(xn):]; minimum=min(cycle)
                while v[0]!=minimum and v[7]<horizon:
                    visit(v); x=v[0]; j=branch[x]; u=unit[x]
                    v[0]=dest[x]
                    assert v[0]>=0
                    v[1]=u*(v[1]+ja[j])%q; v[2]=u*(v[2]+jb[j])%q
                    used+=1; v[7]+=1; st['rewind_ticks']+=1
                if v[7]+2>horizon: st['length_aborts']+=1; return False
                visit(v); x=v[0]; u=g.double_unit[x]
                v[0]=g.double[x]; v[1]=2*u*v[1]%q; v[2]=2*u*v[2]%q
                used+=2; v[7]+=2; st['escapes']+=1; st['dp_tests']+=1
                v[3]=[]; ring=v[3]
                if dp[v[0]]: endpoint(v); return False
            ring.append(v[0])
            if len(ring)>8: ring.pop(0)
            state=(v[0],tuple(ring)); v[6]+=1
            if state==v[4]: st['detector_aborts']+=1; return False
            if v[6]==v[5]: v[4]=state; v[5]*=2; v[6]=0
            return True
        active: list[list | None]=[start() for _ in range(8)]; cursor=0
        while used<T:
            v=active[cursor]
            if v is None: v=start()
            active[cursor]=v if tick(v) else None
            cursor=(cursor+1)%8
        for v in active:
            if v is not None:
                while tick(v): pass
        total.update(st)
        total['ticks']+=used
        total['recovered']+=good
        total['all_recoveries']+=any_recovery
        total['no_recovery']+=not good
        total['core_repeat_trials']+=core_repeat
    return dict(seed=seed,trials=trials,counters=dict(sorted(total.items())),
        repeat_samples=len(repeats),repeat_sum=sum(repeats),repeat_sq_sum=sum(t*t for t in repeats),
        repeat_min=min(repeats) if repeats else None,repeat_max=max(repeats) if repeats else None,
        cycle_sum=sum(cycles),dp_classes=sum(dp),
        trajectory_sha256=hashlib.sha256(json.dumps([repeats,cycles],separators=(',',':')).encode()).hexdigest())


def configs():
    # Co-move the finite sensitivity matrix by one denser exponent step when
    # production spacing widens. Clamp at each scale's exact delta<.30 cap;
    # literal +10 bits is infeasible at these N and is NOT claimed simulated.
    offset=min(1,max(0,z.d-54))
    small=[(8,0,j) for j in range(8)]
    medium=[(w,min(d+offset,1),j) for w in (8,64) for d in (0,1) for j in range(2)]
    large=[(w,min(d+offset,4),0) for w,d in ((8,3),(8,4),(64,2),(64,3),(64,4),(256,2),(256,3),(256,4))]
    return dict(small=small,medium=medium,large=large)


def run():
    groups=scale_groups(); points=[]; per_table=25 if '--quick' in sys.argv else 250
    evidence=HERE/'evidence'; evidence.mkdir(exist_ok=True)
    for scale,settings in configs().items():
        g=groups[scale]; rows=[]
        for ci,(width,d,variant) in enumerate(settings):
            T,envelope=binding_budget(g.N,d)
            tables=[]; hist=None; shift=4
            for tid in range(8):
                seed=91001+100*width+10000*ci+1000000*len(points)+tid
                tr=g.transitions(width,seed,shift)
                hist=tr[4]
                tables.append(one(g,width,seed,d,per_table,T,shift,tr))
            cc=Counter()
            for row in tables: cc.update(row['counters'])
            n=8*per_table; samples=sum(t['repeat_samples'] for t in tables)
            repeat_sum=sum(t['repeat_sum'] for t in tables)
            assert samples>0
            mean_repeat=repeat_sum/samples
            mean_ratio=mean_repeat/math.sqrt(math.pi*g.N/2)
            p=cc['recovered']/n; a=T*(T-1)/(2*g.N)
            point=dict(id=f'{scale}-{ci:02d}',scale=scale,classes=g.N,order=g.Q,width=width,
                variant=variant,dp_bits=d,spacing=1<<d,partition_shift=shift,T=T,
                a=z.frac(F(T*(T-1),2*g.N)),delta_transport=z.frac(delta_transport(d,classes=g.N)),
                delta_ratio_to_production=z.frac(delta_transport(d,classes=g.N)/z.transport),
                envelope_lower=z.frac(envelope),envelope_display=float(envelope),
                trials=n,recovered=cc['recovered'],wrong=cc['wrong_recoveries'],no_recovery=cc['no_recovery'],
                measured_probability=p,observed_above_binding_envelope=F(cc['recovered'],n)>=envelope,
                pipeline_rho=-math.log1p(-p)/a if p<1 else None,
                repeat_samples=samples,repeat_sum=repeat_sum,mean_repeat=mean_repeat,rho_emp=mean_ratio,
                rho_emp_definition='mean bare-map first-repeat length / sqrt(pi*N/2); NOT birthday exponent coefficient',
                branch_histogram=hist,dp_classes=tables[0]['dp_classes'],counters=dict(sorted(cc.items())),tables=tables)
            rows.append(point); points.append(point)
            print(json.dumps({k:point[k] for k in ('id','classes','width','spacing','T','envelope_display','recovered','rho_emp','pipeline_rho')}),flush=True)
        # Pair rows losslessly so the complete archive stays below 32 files.
        for k in range(0,8,2):
            data=json.dumps(rows[k:k+2],separators=(',',':'))+'\n'
            (evidence/f'config-{scale}-{k:02d}-{k+1:02d}.json').write_text(data)
    summary=dict(version=1,scope='Toy-coordinate mechanics and finite binding tests only; no production mixing, retention-rate, or instruction-price certification.',
        groups={name:g.validation for name,g in groups.items()},
        parameters=dict(tables_per_config=8,trials_per_table=per_table,workers=8,ring=8,
            model_rho=z.frac(z.rho),budget_rule='smallest positive multiple of 16 with envelope > 1/50',
            spacing_fixed_per_config=True,full_pool_no_evictions=True,
            production_d=z.d,production_spacing=str(z.DP),
            spacing_comovement='one-bit sensitivity step for d>54, clamped at exact scale feasibility (small0/medium1/large4); not literal production simulation',
            matrix_scope='Small: width8 only, N216; medium: widths8/64, N908; large: widths8/64/256, N87328. Infeasible tiny-scale branch widths excluded, no dead-branch assertion relaxed.'),
        points=[{k:v for k,v in p.items() if k not in ('tables','counters','branch_histogram','rho_emp_definition')} for p in points],
        total_trials=sum(p['trials'] for p in points),
        min_envelope_display=min(p['envelope_display'] for p in points),
        rho_emp_band=[min(p['rho_emp'] for p in points),max(p['rho_emp'] for p in points)],
        pipeline_rho_band=[min(p['pipeline_rho'] for p in points),max(p['pipeline_rho'] for p in points)],
        envelope_compatibility_count=sum(p['observed_above_binding_envelope'] for p in points))
    (HERE/'sweep-summary.json').write_text(json.dumps(summary,separators=(',',':'))+'\n')
    print(json.dumps({k:summary[k] for k in ('total_trials','min_envelope_display','rho_emp_band','pipeline_rho_band','envelope_compatibility_count')}),flush=True)
if __name__=='__main__': run()
