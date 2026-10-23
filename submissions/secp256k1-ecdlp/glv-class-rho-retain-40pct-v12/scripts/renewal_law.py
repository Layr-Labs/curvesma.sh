"""Collision-selected, not uniform-pair, finite renewal diagnostics.
No statistic here is a production theorem. All probability gates use Fractions.
A trial selects its FIRST cross-trail shared-class visit before the core cap.
Both entire selected trails are followed to their actual terminal paths.
"""
from collections import Counter
from fractions import Fraction as F
from math import isqrt

MULTIPLIERS=(1,2,4,40)
PATHS=('dp','identity','detector','exception_cap','trail_limit','pool_full','event_cap')

def frac(x):
    x=F(x)
    return dict(numerator=str(x.numerator),denominator=str(x.denominator))

def rat(v):
    return F(int(v['numerator']),int(v['denominator']))

def bucket(t,s):
    # Dyadic spacing-normalized bins: [0,s),[s,2s),...,[32s,infinity).
    return min(6,(t//s).bit_length())

def finite_reserve(n):
    # eps >= sqrt(5/n); 120*exp(-10)<1/100, proved by rational series.
    u=(5*10**8+n-1)//n
    t=isqrt(u)
    return F(t+(t*t<u),10000)

def joint_null(width,s,h):
    # Null only: DP-before-escape is the actual priority ordering in tick().
    # q = neither DP nor escape; escape hazard is 1/(12*width).
    e=F(1,12*width); dp=F(1,s); q=(1-dp)*(1-e)
    one=1-(1-dp)*e*(1-q**h)/(1-q)
    exact=one*one
    return F((exact.numerator*10**12+exact.denominator-1)//exact.denominator,10**12)

def new_capture():
    return dict(selected=0,completed_pairs=0,ages=Counter(),renewal_ages=Counter(),gaps=Counter(),
        selected_gaps=Counter(),tails=Counter(),paths=Counter(),selected_paths=Counter(),
        selected_nonDP_failure=Counter(),gap_exposure=0,gap_events=0)

def observe(capture,trails,pair,s):
    for rec in trails.values():
        capture['paths'][rec['path']]+=1
        prev=0
        for t in rec['escapes']:
            capture['gaps'][bucket(t-prev,s)]+=1
            capture['gap_events']+=1; prev=t
        capture['gap_exposure']+=rec['length']
    if pair is None: return
    capture['selected']+=1
    r1,r2=[trails[i] for i,_ in pair]
    capture['completed_pairs']+=r1['path']=='dp' and r2['path']=='dp' and r1['endpoint']==r2['endpoint'] and r1['a']!=r2['a']
    records=[]
    for serial,pos in pair:
        rec=trails[serial]
        prev=max([0]+[t for t in rec['escapes'] if t<=pos])
        records.append((rec,pos,prev))
        capture['selected_paths'][rec['path']]+=1
        capture['selected_nonDP_failure'][rec['path']]+=rec['path'] not in ('dp','identity')
        last=0
        for t in rec['escapes']:
            capture['selected_gaps'][bucket(t-last,s)]+=1; last=t
    capture['ages'][tuple(bucket(pos,s) for _,pos,_ in records)]+=1
    capture['renewal_ages'][tuple(bucket(pos-prev,s) for _,pos,prev in records)]+=1
    for m in MULTIPLIERS:
        h=m*s
        good=True
        for rec,pos,_ in records:
            # A DP before h ends exposure successfully, NOT as right censoring.
            # An abort before h is a failure, not a spurious long renewal gap.
            if any(pos<t<=pos+h for t in rec['escapes']): good=False
            if rec['path'] not in ('dp','identity') and rec['length']-pos<=h: good=False
        capture['tails'][m]+=good

def merge(parts):
    out=new_capture()
    for p in parts:
        for k,v in p.items():
            if isinstance(v,Counter): out[k].update(v)
            else: out[k]+=v
    return out

def compact_hist(c):
    return [[*(k if isinstance(k,tuple) else (k,)),v] for k,v in sorted(c.items()) if v]

def finish(c,point):
    n=c['selected']; s=point['spacing']; B=point['width']; trials=point['trials']
    assert n>0
    return dict(id=point['id'],classes=point['classes'],width=B,spacing=s,
        selected=n,completed_pairs=c['completed_pairs'],
        age_hist=compact_hist(c['ages']),renewal_age_hist=compact_hist(c['renewal_ages']),
        inter_escape_gap_hist=compact_hist(c['gaps']),selected_gap_hist=compact_hist(c['selected_gaps']),
        gap_exposure=c['gap_exposure'],gap_events=c['gap_events'],
        tails=[[m,c['tails'][m]] for m in MULTIPLIERS],
        terminal_counts={k:c['paths'][k] for k in PATHS},
        selected_terminal_counts={k:c['selected_paths'][k] for k in PATHS},
        path_trial_counts={k:point['path_trial_counts'].get(k,0) for k in PATHS})

def path_ratios(row,point):
    s=point['spacing']; B=point['width']; K=point['T']+8*40*s
    assumed=dict(detector=min(F(1),F(8*(40*s)**2,point['classes'])),
        exception_cap=min(F(1),F(4*B*K,point['classes']*2**16)),
        trail_limit=F(44,2**24),pool_full=min(F(1),F(11*K,10*s*point['pool_capacity'])),
        event_cap=min(F(1),F(11*K,10*s*(point['event_capacity']-8))))
    return {k:F(row['path_trial_counts'][k],point['trials'])/v for k,v in assumed.items()}

def tail_rows(row):
    n=row['selected']; eps=finite_reserve(n)
    return [dict(window=m,successes=k,samples=n,null=joint_null(row['width'],row['spacing'],m*row['spacing']),
        observed=F(k,n),lower=max(F(0),F(k,n)-eps)) for m,k in row['tails']]

def lower_envelope(rows):
    low=min([t['lower']/t['null'] for row in rows for t in tail_rows(row)]+
        [max(F(0),F(row['completed_pairs'],row['selected'])-finite_reserve(row['selected'])) for row in rows])
    return max(F(0),min(F(1),F(low.numerator*10000//low.denominator,10000)))
