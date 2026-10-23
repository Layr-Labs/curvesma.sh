#!/usr/bin/env python3
"""Unmodified universe, trace and histogram helpers from the v7 audit.
An exhaustive exponent oracle accelerates literal tower transitions only.
"""
from array import array
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import random
from typing import Any
from tower import Tower, exact_order

SEED=2026102307
CONFIGS=[(14009539,12289,64),(2677963,49177,256),
         (159508651,196657,1024),(19066387,786433,1024),
         (19066387,786433,8),(19066387,786433,32),(19066387,786433,128)]
DENSITIES=(2,4,6)
CACHE={}

def fraction(x):
    x=F(x); return {'numerator':x.numerator,'denominator':x.denominator}

def prime(n):
    from math import isqrt
    return n>1 and (n==2 or (n%2 and all(n%d for d in range(3,isqrt(n)+1,2))))

def universe(q,r,B):
    if (q,r) not in CACHE:
        assert prime(q) and prime(r) and exact_order(q%r,r,12)
        f=Tower(q); rng=random.Random(SEED^q)
        assert f.f2pow((1,1),(q*q-1)//2)!=(1,0)
        assert f.f2pow((1,1),(q*q-1)//3)!=(1,0)
        while True:
            g=f.pow(tuple(rng.randrange(q) for _ in range(12)),(q**12-1)//r)
            if g!=f.one:break
        assert f.pow(g,r)==f.one
        L=[pow(q,k,r) for k in range(12)]
        cid=array('i',[-1])*r; cid[0]=0
        align=array('i',[1])*r
        keys=[0]; reps=[0]; literal=0
        y=f.one
        for e in range(1,r):
            y=f.mul(y,g)
            if cid[e]>=0:continue
            K,k,rep=f.metadata(y); ce=e*L[k]%r
            orbit=[e*a%r for a in L]; assert len(set(orbit))==12
            i=len(keys)
            for j,ep in enumerate(orbit):
                assert cid[ep]==-1
                cid[ep]=i; align[ep]=L[(k-j)%12]
                assert ep*align[ep]%r==ce
            if literal<24:
                assert f.pow(g,ce)==rep
                assert tuple((f.frob(y,j)[1]+f.frob(y,j)[7])%q for j in range(12))==f.keys(y)
                literal+=1
            keys.append(K);reps.append(ce)
        assert min(cid)>=0 and len(keys)==1+(r-1)//12
        CACHE[q,r]=(f,g,L,cid,align,keys,reps,literal)
    f,g,L,cid,align,keys,reps,literal=CACHE[q,r]
    hist=Counter(K%B for K in keys[1:]); census=[hist[h] for h in range(B)]
    assert min(census)>0,'DEAD BRANCH'
    meta: dict[str,Any]=dict(q=q,r=r,classes=(r-1)//12,branches=B,
        exhaustive_branch_histogram=census,histogram_min=min(census),histogram_max=max(census),
        dead_branches=0,literal_tower_crosschecks=literal,generator=list(g),
        exhaustive_dp_counts={str(d):sum(((K>>(B.bit_length()-1))&((1<<d)-1))==0 for K in keys[1:]) for d in DENSITIES},
        coordinate_metadata_sha256=hashlib.sha256(json.dumps([keys,reps],separators=(',',':')).encode()).hexdigest())
    return f,g,L,cid,align,keys,reps,meta

def trace(step,start,keys,shift,base,cap):
    seen={}; states=[]; first: list[int|None]=[None]*len(DENSITIES); a=start
    while a not in seen:
        t=len(states); seen[a]=t; states.append(a); K=keys[a]
        for j,d in enumerate(DENSITIES):
            if first[j] is None and ((K>>shift)&((1<<d)-1))==0:first[j]=t
        a=step(a)
    mu=seen[a]; ell=len(states)-mu; meet=ell*max(1,(mu+ell-1)//ell)
    row=[mu,ell,len(states),meet,int(len(states)<=base),int(len(states)<=cap),int(meet<=cap),*first]
    return row,states,a

def distribution(vals):
    s=sorted(vals); width=max(1,(s[-1]+15)//16); c=Counter(v//width for v in vals)
    return dict(count=len(vals),sum=sum(vals),min=s[0],max=s[-1],
        quantiles={str(q):s[(len(s)-1)*q//100] for q in (10,25,50,75,90,95,99)},
        histogram_width=width,histogram=[[k*width,c[k]] for k in sorted(c)])

def summary(rows):
    out: dict[str,Any]={k:distribution([row[i] for row in rows]) for i,k in enumerate(
        ['mu','cycle_lambda','first_repeat','floyd_meeting'])}
    out.update(trials=len(rows),repeat_at_base=sum(r[4] for r in rows),
        repeat_at_stretch=sum(r[5] for r in rows),floyd_at_stretch=sum(r[6] for r in rows),dp={})
    for i,d in enumerate(DENSITIES):
        hits=[r[7+i] for r in rows if r[7+i] is not None]
        out['dp'][str(d)]=dict(hit_before_repeat=len(hits),no_hit_before_repeat=len(rows)-len(hits),
                              first_hit=distribution(hits) if hits else None)
    return out
