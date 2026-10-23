#!/usr/bin/env python3
"""Same-table, actual-coordinate comparison of three stationary class maps.
The full orientation guard squares a tentative R-add precisely when that
R-add needs a nontrivial canonical Frobenius. No history, retries or free work.
Raw rows are optionally captured outside the upload; summaries are lossless.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
import time
from typing import Any
from base_audit import universe, trace, summary, fraction, CONFIGS, DENSITIES
from ledger import ln_bounds, ceil_sqrt_fraction

ROOT=Path(__file__).resolve().parents[1]
SEED=2026102309
TRIALS=25000
VARIANTS=('baseline','inverse_guard','orientation_guard')
RESERVE=F(109,100)


def transition(i,jumps,reps,keys,cid,align,r,variant):
    B=len(jumps); h=keys[i]%B; raw=(reps[i]+jumps[h])%r
    alpha=align[raw]; j=cid[raw]
    guarded=(variant=='orientation_guard' and alpha!=1) or (
        variant=='inverse_guard' and alpha==r-1 and keys[j]%B==h)
    c=2 if guarded else 1
    raw=c*raw%r
    return cid[raw],c*align[raw]%r,h,int(guarded)


def labels(states,row,jumps,beta,x,D0,reps,keys,cid,align,r,variant):
    mu,ell,_,t=row[:4]; A=(reps[states[0]]-x*D0)%r; D=D0
    prod=1; weights={}; snap=None; squares=0; cycle_squares=0
    for n in range(2*t):
        i=states[n] if n<len(states) else states[mu+(n-mu)%ell]
        j,m,h,sq=transition(i,jumps,reps,keys,cid,align,r,variant)
        weights[h]=(weights.get(h,0)+pow(prod,-1,r))%r
        prod=prod*m%r
        A=m*(A+jumps[h]-x*beta[h])%r; D=m*(D+beta[h])%r
        assert (A+x*D)%r==reps[j]
        squares+=sq
        if t<=n<2*t:cycle_squares+=sq
        if n+1==t:snap=(A,D,prod,dict(weights))
    assert snap is not None
    At,Dt,pt,wt=snap
    den=(Dt-D)%r; num=(A-At)%r
    structural=(pt-prod)%r==0 and all((pt*wt.get(h,0)-prod*weights.get(h,0))%r==0
        for h in weights.keys()|wt.keys())
    assert not structural or den==num==0
    recovered=den!=0 and num*pow(den,-1,r)%r==x
    assert den==0 or recovered
    # all lambda factors are 1 on unguarded orientation-guard edges.
    if variant=='orientation_guard' and cycle_squares==0:
        assert prod==pt
    return [int(den==0),int(structural),int(recovered),int(recovered and row[6]),
            int(ell==2),int(ell==2 and structural),squares,cycle_squares]


def run_config(q,r,B,trials,seed):
    f,g,L,cid,align,keys,reps,meta=universe(q,r,B)
    N=(r-1)//12; np=(N*B+B-2)//(B-1)
    base=ceil_sqrt_fraction(2*np*ln_bounds(F(5,3))[1]); cap=(109*base+99)//100
    rng=random.Random(seed); data={v:[] for v in VARIANTS}; labs={v:[] for v in VARIANTS}
    control=[]; visits={v:Counter() for v in VARIANTS}
    counts={v:Counter() for v in VARIANTS}; literal={v:0 for v in VARIANTS}
    for trial in range(trials):
        jumps=[rng.randrange(r) for _ in range(B)]; beta=[rng.randrange(r) for _ in range(B)]
        x=rng.randrange(r); D0=rng.randrange(r); start=cid[rng.randrange(1,r)]
        for v in VARIANTS:
            used=set()
            def step(i):
                j,m,h,sq=transition(i,jumps,reps,keys,cid,align,r,v)
                visits[v][h]+=1; counts[v]['transitions']+=1
                counts[v]['reused']+=h in used; used.add(h); counts[v]['squares']+=sq
                if literal[v]<32:
                    y=f.pow(g,reps[i]); J=f.pow(g,jumps[h]); raw=f.mul(y,J)
                    if sq:raw=f.mul(raw,raw)
                    assert f.metadata(raw)[2]==f.pow(g,reps[j])
                    # Explicitly noncanonical representative with twisted jump.
                    k=literal[v]%12; yn=f.frob(y,k); got=f.mul(yn,f.pow(J,L[k]))
                    if sq:got=f.mul(got,got)
                    assert f.metadata(got)[2]==f.pow(g,reps[j])
                    literal[v]+=1
                return j
            row,states,_=trace(step,start,keys,B.bit_length()-1,base,cap)
            data[v].append(row)
            labs[v].append(labels(states,row,jumps,beta,x,D0,reps,keys,cid,align,r,v))
        mapping={}
        def ctrl(i):
            if i not in mapping:mapping[i]=rng.randrange(len(keys))
            return mapping[i]
        cr,_,_=trace(ctrl,start,keys,B.bit_length()-1,base,cap); control.append(cr)
    cc=summary(control); variants={}
    names=('zero_denominator','structurally_fruitless','verified_recovery_eventual',
           'verified_recovery_within_cap','two_cycles','structurally_fruitless_two_cycles')
    for v in VARIANTS:
        ss=summary(data[v]); lc={key:sum(z[j] for z in labs[v]) for j,key in enumerate(names)}
        useful_repeat=sum(int(z[2] and row[5]) for row,z in zip(data[v],labs[v]))
        diff=F(useful_repeat-cc['repeat_at_base'],trials)
        variants[v]=dict(shared=ss,labels=lc,useful_first_repeat_at_stretch=useful_repeat,
            paired=dict(mean_repeat_ratio=fraction(F(ss['first_repeat']['sum'],cc['first_repeat']['sum'])),
                useful_repeat_CDF_difference=fraction(diff),toy_only_lower_radius_003=fraction(diff-F(3,100))),
            visited_branch_histogram=[visits[v][h] for h in range(B)],
            **dict(counts[v]),literal_R_adding_twisted_and_guard_checks=literal[v])
    meta.update(seed=seed,trials=trials,base_cap=base,stretched_iteration_cap=cap,
                comparison_Np_ceiling=np,control=cc,variants=variants)
    raw=json.dumps(dict(data=data,labels=labs,control=control),separators=(',',':')).encode()
    meta['raw_sha256']=hashlib.sha256(raw).hexdigest()
    return meta,raw


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true')
    p.add_argument('--trials',type=int,default=TRIALS);p.add_argument('--raw-dir',type=Path)
    a=p.parse_args(); start=time.monotonic()
    if a.raw_dir:a.raw_dir.mkdir(parents=True,exist_ok=True)
    result: dict[str,Any]=dict(version=1,seed=SEED,trials_per_config=a.trials,reserve=fraction(RESERVE),
        scope='Actual tower keys and retained tables. Baseline failures retained. Toy evidence, not a production theorem.',
        variants=list(VARIANTS),configurations=[])
    for i,(q,r,B) in enumerate(CONFIGS):
        row,raw=run_config(q,r,B,a.trials,SEED+i)
        if a.write:
            if a.raw_dir:
                import zlib
                (a.raw_dir/f'raw-{i}.json.zlib').write_bytes(zlib.compress(raw,9))
            text=json.dumps(row,separators=(',',':'))+'\n'; parts=[]
            for j in range(0,len(text),14000):
                name=f'evidence/config-{i}-{j//14000}.txt'
                (ROOT/name).write_text(text[j:j+14000]);parts.append(name)
            result['configurations'].append(dict(parts=parts,q=q,r=r,branches=B))
            (ROOT/'toy_results.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(q=q,r=r,B=B,trials=a.trials,branch_min=row['histogram_min'],
            control_base=row['control']['repeat_at_base'],variants={v:dict(labels=z['labels'],
            useful_repeat=z['useful_first_repeat_at_stretch'],squares=z['squares'],steps=z['transitions'],
            paired=z['paired']) for v,z in row['variants'].items()})),flush=True)
    result['elapsed_seconds']=time.monotonic()-start
    if a.write:(ROOT/'toy_results.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
