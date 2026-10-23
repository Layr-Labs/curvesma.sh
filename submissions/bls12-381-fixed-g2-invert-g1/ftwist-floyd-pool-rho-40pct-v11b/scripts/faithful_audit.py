#!/usr/bin/env python3
"""Faithful capture of the repaired v11 semantics: Floyd stalls, persistent
DP credit, per-escape-charged in-place jumps under the priced global
allowance. The old hard E=4 run is retained as an ablation only.
All labels and branch/DP keys are literal tower-derived, never exponent hashes.
"""
import argparse,hashlib,json,random,time,zlib
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from base_audit import universe,trace,summary,fraction,CONFIGS,DENSITIES
from ledger import ln_bounds,ceil_sqrt_fraction,E4_ABLATION_ESCAPE_CAP
from retention import run_walk
from controls import run_walk as v10_walk,old_floyd
ROOT=Path(__file__).resolve().parents[1]
SEED=2026102310
TRIALS=2500


def run_config(q,r,B,trials,seed):
    f,g,L,cid,align,keys,reps,meta=universe(q,r,B)
    N=(r-1)//12;np=(N*B+B-2)//(B-1)
    t0=ceil_sqrt_fraction(2*np*ln_bounds(F(5,3))[1]);cap=(109*t0+99)//100
    ordinary_cap=3*cap
    # Priced global first-moment allowance for this finite T: ceil(512*mean).
    allowance=(512*ordinary_cap+12*B-1)//(12*B)
    capacity=N+1;shift=B.bit_length()-1;rng=random.Random(seed)
    counts={d:Counter() for d in DENSITIES};visits=Counter();labels_checked=0
    rawrows=[];controlrows=[];oldcounts=Counter();witnesses={};literal=0
    maxesc=[0]
    per_trial=[]
    for trial in range(trials):
        jumps=[rng.randrange(r) for _ in range(B)];beta=[rng.randrange(r) for _ in range(B)]
        offsets=[rng.randrange(r) for _ in range(B)];obeta=[rng.randrange(r) for _ in range(B)]
        x=rng.randrange(r);D0=rng.randrange(r);start=cid[rng.randrange(1,r)]
        s0=(start,(reps[start]-x*D0)%r,D0)
        def basic(i):
            h=keys[i]%B;v=(reps[i]+jumps[h])%r
            return cid[v],align[v],h
        def advance(st):
            nonlocal labels_checked
            i,A,D=st;j,m,h=basic(i)
            out=(j,m*(A+jumps[h]-x*beta[h])%r,m*(D+beta[h])%r)
            assert (out[1]+x*out[2])%r==reps[j]
            labels_checked+=1;visits[h]+=1
            return out
        def escaper():
            erng=random.Random((seed<<32)^trial^0xE5CA9E)
            def esc(st):
                nonlocal labels_checked
                i,A,D=st;h=erng.randrange(B);v=(reps[i]+offsets[h])%r;m=align[v]
                out=(cid[v],m*(A+offsets[h]-x*obeta[h])%r,m*(D+obeta[h])%r)
                assert (out[1]+x*out[2])%r==reps[out[0]]
                labels_checked+=1
                return out
            return esc
        row,states,_=trace(lambda i:basic(i)[0],start,keys,shift,t0,cap)
        rawrows.append(row);cycle=states[row[0]:]
        st=s0
        for _ in range(row[0]):st=advance(st)
        cycle_start=st
        for _ in cycle:st=advance(st)
        structural=st==cycle_start
        mapping={};crng=random.Random((seed<<32)^trial^0xC01720)
        def fresh(i):
            if i not in mapping:mapping[i]=crng.randrange(len(keys))
            return mapping[i]
        crow,_,_=trace(fresh,start,keys,shift,t0,cap);controlrows.append(crow)
        got,es,os=old_floyd(advance,escaper(),s0,cap,r)
        assert got is None or got==x
        oldcounts.update(success=int(got is not None),escapes=es,ordinary_steps=os)
        trial_data=[int(got is not None)]
        for d in DENSITIES:
            isdp=lambda i,d=d:((keys[i]>>shift)&((1<<d)-1))==0
            probe=Counter()
            def observer(tort,hare,steps):
                # Read-only cycle census after a true Floyd event. This is
                # instrumented evidence, never a production escape trigger.
                cur=tort[0];cs=[]
                while cur not in cs:
                    cs.append(cur);cur=basic(cur)[0]
                cyc=cs[cs.index(cur):]
                blind=not any(isdp(i) for i in cyc)
                probe['fruitless_DP_free_caught']+=int(blind)
                probe['fruitless_DP_free_two_caught']+=int(blind and len(cyc)==2)
            z=run_walk(advance,escaper(),s0,isdp,capacity,ordinary_cap,allowance,r,observer=observer)
            old10=v10_walk(advance,escaper(),s0,isdp,capacity,ordinary_cap,allowance,r)
            ablation=run_walk(advance,escaper(),s0,isdp,capacity,ordinary_cap,E4_ABLATION_ESCAPE_CAP,r)
            for out in (z,old10,ablation):
                assert out['recovered_scalar'] is None or out['recovered_scalar']==x
            c=counts[d]
            for k,v in z.items():
                if k!='recovered_scalar':c[k]+=v
            c.update(probe)
            c['capped_failures']+=int(not z['success'])
            c['v10_success']+=old10['success']
            c['e4_ablation_success']+=ablation['success']
            c['e4_ablation_escapes']+=ablation['escapes']
            c['e4_ablation_allowance_failures']+=ablation['allowance_failure']
            c['new_not_v10']+=int(z['success'] and not old10['success'])
            c['v10_not_new']+=int(old10['success'] and not z['success'])
            c['new_not_old_floyd']+=int(z['success'] and got is None)
            c['old_floyd_not_new']+=int(not z['success'] and got is not None)
            c['trails_above_four_escapes']+=int(z['escapes']>4)
            maxesc[0]=max(maxesc[0],z['escapes'])
            stuck=bool(structural and not any(isdp(i) for i in cycle))
            c['initial_fruitless_DP_free']+=int(stuck)
            c['initial_fruitless_DP_free_two']+=int(stuck and len(cycle)==2)
            # Catch-rate denominator contains only reachable first meetings:
            # earlier useful DP recovery is a success, not a missed cycle.
            reachable=stuck and row[3]<=cap and not z['success']
            c['unrecovered_initial_blind_reachable']+=int(reachable)
            c['unrecovered_initial_blind_caught']+=int(reachable and z['floyd_meetings']>0)
            if stuck and str(d) not in witnesses:
                direct=run_walk(advance,escaper(),s0,lambda i:False,capacity,ordinary_cap,1,r)
                assert row[3]>cap or direct['fruitless_meetings']>0
                witnesses[str(d)]=dict(trial=trial,cycle=cycle,keys=[keys[i] for i in cycle],
                    states=[list(cycle_start),list(st)],mu=row[0],cycle_length=len(cycle),
                    floyd_first_meeting=row[3],direct_probe=direct,target_x=x,
                    edges=[dict(source=i,source_rep=reps[i],branch=keys[i]%B,
                        jump=jumps[keys[i]%B],beta=beta[keys[i]%B],destination=basic(i)[0],
                        destination_rep=reps[basic(i)[0]],alignment=basic(i)[1]) for i in cycle])
            trial_data.extend([z[k] for k in ('success','escapes','cross_escape_recovery','reports','allowance_failure')])
            trial_data.extend([old10['success'],ablation['success'],probe['fruitless_DP_free_caught']])
        per_trial.append(trial_data)
        if literal<32:
            i=start;h=keys[i]%B;y=f.pow(g,reps[i]);J=f.pow(g,jumps[h]);j,_,_=basic(i)
            assert f.metadata(f.mul(y,J))[2]==f.pow(g,reps[j])
            k=literal%12;yn=f.frob(y,k);_,kap,_=f.metadata(yn)
            assert f.metadata(f.mul(yn,f.frob(J,-kap)))[2]==f.pow(g,reps[j])
            literal+=1
    fixed_jumps=[rng.randrange(r) for _ in range(B)];fixed=[]
    for _ in range(256):
        def fs(i):return cid[(reps[i]+fixed_jumps[keys[i]%B])%r]
        row,_,_=trace(fs,cid[rng.randrange(1,r)],keys,shift,t0,cap);fixed.append(row)
    result={}
    for d,c in counts.items():
        result[str(d)]=dict(c,success_fraction=fraction(F(c['success'],trials)),
            endpoint_minus006=fraction(F(c['success'],trials)-F(3,50)),
            escape_rate=fraction(F(c['escapes'],max(1,c['ordinary_steps']))),
            prediction=fraction(F(1,12*B)),
            blind_catch_fraction=fraction(F(c['unrecovered_initial_blind_caught'],max(1,c['unrecovered_initial_blind_reachable']))))
    meta.update(seed=seed,trials=trials,base_cap=t0,floyd_cap=cap,
        ordinary_cap=ordinary_cap,escape_allowance=allowance,
        e4_ablation_escape_cap=E4_ABLATION_ESCAPE_CAP,pool_capacity=capacity,
        fixed_DP_bits=list(DENSITIES),old_floyd=dict(oldcounts),pipelines=result,
        variants=dict(shared=summary(rawrows),fresh=summary(controlrows)),
        live_branch_histogram=[visits[h] for h in range(B)],
        max_trail_escapes=maxesc[0],
        label_invariant_checks=labels_checked,literal_raw_twist_checks=literal,
        no_DP_fruitless_witnesses=witnesses,
        fixed_table_control=dict(starts=256,scope='Correlated; descriptive only.',summary=summary(fixed)))
    raw=json.dumps(per_trial,separators=(',',':')).encode()
    meta['paired_raw_sha256']=hashlib.sha256(raw).hexdigest()
    meta['paired_raw_columns']=['old_success']+[f'd{d}_{k}' for d in DENSITIES for k in
        ('success','escapes','cross_escape_recovery','reports','allowance_failure','v10_success','e4_ablation_success','fruitless_DP_free_caught')]
    return meta,raw


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=TRIALS)
    ap.add_argument('--raw-dir',type=Path);a=ap.parse_args();start=time.monotonic_ns()
    if a.raw_dir:a.raw_dir.mkdir(parents=True,exist_ok=True)
    out=dict(version=1,seed=SEED,trials_per_config=a.trials,
        scope='Literal tower coordinate partition, exact labels, shared/fresh controls. Floyd plus zero-delete pool; per-trail escapes are unbounded but each is charged, termination is bounded by the 3T trail cap, and the priced global first-moment allowance fails closed. The hard E=4 run is the historical ablation, not the claimed algorithm. Toy evidence is not production mixing, sparse-DP transport or instruction-price proof.',configurations=[])
    for i,(q,r,B) in enumerate(CONFIGS):
        c,raw=run_config(q,r,B,a.trials,SEED+i)
        text=json.dumps(c,separators=(',',':'))+'\n';parts=[]
        for j in range(0,len(text),14000):
            name=f'evidence/config-{i}-{j//14000}.txt'
            (ROOT/name).parent.mkdir(exist_ok=True)
            (ROOT/name).write_text(text[j:j+14000]);parts.append(name)
        if a.raw_dir:(a.raw_dir/f'paired-{i}.json.zlib').write_bytes(zlib.compress(raw,9))
        out['configurations'].append(dict(index=i,q=q,r=r,branches=B,parts=parts,summary_sha256=hashlib.sha256(text.encode()).hexdigest()))
        out['elapsed_nanoseconds']=time.monotonic_ns()-start
        (ROOT/'toy_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(configurations=len(out['configurations']),trials=a.trials,
        elapsed_nanoseconds=out['elapsed_nanoseconds'],status='CAPTURED')),flush=True)
if __name__=='__main__':main()
