"""Coordinate-faithful, capped end-to-end GLV spacing/table sweep.
Run python3 -B toy_sweep.py [--quick]. Output files are deterministic except time.
Fixed partition bits 14..(14+log2 R-1) never overlap DP bits 0..13.
Eight tables x 250 fresh targets/starts = 2000 trials at EACH ratio point.
8 interleaved walkers; core-overlap credit only; ordinary=1, doubling=2 ticks.
"""
import hashlib
import json
import math
import random
import sys
from collections import Counter,defaultdict
from pathlib import Path
from toy_group import Group,P,B,Q,N
HERE=Path(__file__).resolve().parent
TRIALS=20 if '--quick' in sys.argv else 250
TABLES=8
WIDTHS=(8,32,256,1024)
DEPTHS=(4,6,8,10,12,14)
WORKERS=8
T=318

def one(g,R,seed,d,trials,trans):
    jumps,branch,dest,unit,hist=trans
    rng=random.Random(seed*101+d*10000019)
    xs,kl=g.xs,g.kl
    mask=(1<<d)-1
    horizon=40*(1<<d)
    dp=[not (x&mask) for x in xs]
    out=Counter()
    count_hist=Counter()
    for trial in range(trials):
        secret=rng.randrange(1,Q)
        bs=[rng.randrange(Q) for _ in range(R)]
        aa=[(k-b*secret)%Q for k,b in zip(jumps,bs)]
        used=0; serial=0; good=False; recovered=False; core_collision=False
        pool=defaultdict(list); first={}
        st=Counter()
        def start():
            nonlocal serial
            serial+=1; st['starts']+=1
            x=rng.randrange(N); b=rng.randrange(Q)
            return [x,(kl[x]-b*secret)%Q,b,[],(x,()),1,0,0,set(),serial]
        def visit(v):
            nonlocal core_collision
            if used<T:
                x=v[0]; v[8].add(x)
                if x in first and first[x]!=v[9]: core_collision=True
                first.setdefault(x,v[9])
        def endpoint(v):
            nonlocal good,recovered
            x,a,b=v[:3]; st['dp_events']+=1
            # Separate deterministic toy 12-bit fingerprint, not production Toeplitz.
            key=((xs[x]*0x9e3779b1)^(xs[x]>>19))&4095
            for ox,oa,ob,oc in pool[key]:
                if ox!=x: st['spurious_keys']+=1; continue
                if ob==b: st['degenerate_relation_pairs']+=1; continue
                candidate=(oa-a)*pow((b-ob)%Q,-1,Q)%Q
                assert candidate==secret and g.points[candidate]==g.points[secret]
                recovered=True
                if v[8].intersection(oc): good=True
            pool[key].append((x,a,b,v[8]))
        def advance(v):
            nonlocal used,good,recovered
            if v[7]>=horizon: st['length_aborts']+=1; return False
            visit(v)
            x=v[0]; j=branch[x]; u=unit[x]; xn=dest[x]
            # Identity canonicalisation has unit 1, not the oracle's zero slot.
            if xn<0: u=1
            v[1]=u*(v[1]+aa[j])%Q; v[2]=u*(v[2]+bs[j])%Q
            if g.points[kl[x]][0]==g.points[jumps[j]][0] if jumps[j] else False:
                st['equal_x_events']+=1
            used+=1; v[7]+=1; st['dp_tests']+=1
            if xn<0:
                st['identity_events']+=1
                if v[2]:
                    assert (-v[1]*pow(v[2],-1,Q))%Q==secret
                    # Direct identity recovery is separately reported, not rho credit.
                    recovered=True
                return False
            v[0]=xn
            if dp[xn]: endpoint(v); return False
            ring=v[3]
            if xn in ring:
                cycle=ring[ring.index(xn):]
                minimum=min(cycle)
                while v[0]!=minimum and v[7]<horizon:
                    visit(v); x=v[0]; j=branch[x]; u=unit[x]
                    v[0]=dest[x]; v[1]=u*(v[1]+aa[j])%Q; v[2]=u*(v[2]+bs[j])%Q
                    used+=1; v[7]+=1; st['rewind_ticks']+=1
                if v[7]+2>horizon: st['length_aborts']+=1; return False
                visit(v); x=v[0]; u=g.double_unit[x]
                v[0]=g.double[x]; v[1]=2*u*v[1]%Q; v[2]=2*u*v[2]%Q
                used+=2; v[7]+=2; st['escapes']+=1; st['dp_tests']+=1
                v[3]=[]; ring=v[3]
                if dp[v[0]]: endpoint(v); return False
            ring.append(v[0])
            if len(ring)>8: ring.pop(0)
            state=(v[0],tuple(ring)); v[6]+=1
            if state==v[4]: st['detector_aborts']+=1; return False
            if v[6]==v[5]: v[4]=state; v[5]*=2; v[6]=0
            return True
        active=[start() for _ in range(WORKERS)]
        cursor=0
        while used<T:
            v=active[cursor]
            if v is None: v=start()
            active[cursor]=v if advance(v) else None
            cursor=(cursor+1)%WORKERS
        for v in active:
            if v is not None:
                while advance(v): pass
        out.update(st); out['ticks']+=used
        out['successes']+=good; out['all_recoveries']+=recovered
        out['core_collision_trials']+=core_collision
        out['trials_with_detector_abort']+=bool(st['detector_aborts'])
        out['core_collision_no_recovery']+=bool(core_collision and not good)
        count_hist[st['dp_events']]+=1
    p=out['successes']/trials
    rho=-math.log1p(-p)/(T*(T-1)/(2*N)) if p<1 else None
    return dict(seed=seed,dp_bits=d,nominal_spacing=1<<d,ratio_num=1<<d,ratio_den=R,
        trials=trials,successes=out['successes'],rho=rho,counters=dict(sorted(out.items())),
        dp_class_count=sum(dp),dp_event_histogram={str(k):v for k,v in sorted(count_hist.items())})

def run():
    g=Group()
    summary=[]
    for R in WIDTHS:
        rows=[]; hist=None
        for tid in range(TABLES):
            seed=91001+R*100+tid
            tr=g.transitions(R,seed)
            if hist is None: hist=tr[4]
            assert hist==tr[4]
            for d in DEPTHS: rows.append(one(g,R,seed,d,TRIALS,tr))
            print(json.dumps(dict(progress_width=R,table=tid+1,tables=TABLES)),file=sys.stderr,flush=True)
        # Lossless pairs satisfy both the 32-file intake and small per-file cap.
        records=[]
        for d in DEPTHS:
            rr=[r for r in rows if r['dp_bits']==d]
            n=sum(v['trials'] for v in rr); s=sum(v['successes'] for v in rr)
            p=s/n
            c=Counter()
            for v in rr: c.update(v['counters'])
            record=dict(width=R,dp_bits=d,spacing=1<<d,ratio=(1<<d)/R,
                trials=n,successes=s,probability=p,
                rho=-math.log1p(-p)/(T*(T-1)/(2*N)) if p<1 else None,
                counters=dict(sorted(c.items())),dp_class_count=rr[0]['dp_class_count'],
                measured_dp_event_density=c['dp_events']/c['dp_tests'],
                nominal_density=2.0**-d,
                mean_dp_events=c['dp_events']/n,
                mean_ticks=c['ticks']/n,
                table_rows=rr)
            records.append(record)
            summary.append({k:v for k,v in record.items() if k!='table_rows'})
        for i in range(0,len(records),2):
            pair=records[i:i+2]
            fn=f"sweep-R{R}-d{pair[0]['dp_bits']}-{pair[1]['dp_bits']}.json"
            (HERE/fn).write_text(json.dumps(pair,separators=(',',':'))+'\n')
        (HERE/f'branches-R{R}.json').write_text(json.dumps(dict(width=R,partition_shift=14,
            full_class_histogram=hist,total_classes=sum(hist),minimum=min(hist),maximum=max(hist)))+'\n')
    root=dict(version=1,scope='Toy subgroup only. No production mixing, rare-basin or instruction-price certification.',
        group=g.validation,parameters=dict(widths=WIDTHS,depths=DEPTHS,tables=TABLES,
        trials_per_table=TRIALS,trials_per_ratio=TABLES*TRIALS,workers=WORKERS,
        core_ticks=T,ring=8,partition_shift=14,drain_tick_cap='40*2^d',seed_rule='91001+100*R+table_id'),
        points=[{k:v for k,v in row.items() if k!='counters'} for row in summary])
    (HERE/'sweep-summary.json').write_text(json.dumps(root,separators=(',',':'))+'\n')
    print(json.dumps(dict(points=len(summary),total_trials=sum(x['trials'] for x in summary),
        trend=[[v['width'],v['dp_bits'],v['rho']] for v in summary])))
if __name__=='__main__': run()
