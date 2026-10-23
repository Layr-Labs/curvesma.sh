#!/usr/bin/env python3
# toy_hard.py — replicates HARDENED v3 class-walk rules on toy groups.
# Production analogue: GT classes under e->e*p^j mod r, ord_r(p)=12.
#   K = min orbit element (class invariant), kappa = argmin.
#   LADDER (class state ONLY): m = 2 if ctz(K) >= s else 1.
#   STEP: e' = e*p^kappa + m*(a_h + b_h*x) mod r, h = K mod R.
#     (walk tracks A,B with e = A + B*x: A' = A*p^kap + m*a_h, B' = B*p^kap + m*b_h)
#   DP iff K mod 2^d == 0: match in table on fp = K>>d; different (A,B)
#     => x = (A1-A2)/(B2-B1) mod r, verify against ground truth.
#   tau = own steps since own last DP; tau >= taumax -> salted restart
#     (fresh e,A,B; lockstep rounds).
# Audits: (C1) merge-preservation: same class at t => same class at t+1
#   unless one side salted-restarted; DP simultaneity for merged pairs.
#   (C2) cycle-length distribution ladder-off vs exact random-mapping
#   formula P(on cycle exactly k) = (N-1)_{k-1}/N^k, and ladder-on cycles.
import random, math, json
from fractions import Fraction

def is_prime(n):
    if n<2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31):
        if n%q==0: return n==q
    d,s=n-1,0
    while d%2==0: d//=2; s+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31):
        y=pow(a,d,n)
        if y in (1,n-1): continue
        for _ in range(s-1):
            y=y*y%n
            if y==n-1: break
        else: return False
    return True

def order(a,n):
    o,x=1,a%n
    while x!=1:
        x=x*a%n; o+=1
        if o>n: raise ValueError
    return o

def params(rmin):
    r=rmin+((1-rmin)%12)
    while True:
        if r%12==1 and is_prime(r):
            for p in range(2,r):
                if order(p,r)==12: return r,p
        r+=12

def ctz(v,cap=64):
    u=0
    while v%2==0 and u<cap: v//=2; u+=1
    return u

class Map:
    def __init__(s,r,p,R,x,seed):
        s.r,s.p,s.R=r,p,R
        s.P=[pow(p,j,r) for j in range(12)]
        rr=random.Random(seed)
        s.a=[rr.randrange(r) for _ in range(R)]
        s.b=[rr.randrange(r) for _ in range(R)]
        s.x=x
    def canon(s,e):
        v=[e*q%s.r for q in s.P]; K=min(v); return K,v.index(K)
    def step(s,e,s_esc):
        K,k=s.canon(e)
        m=2 if ctz(K,s_esc+8)>=s_esc else 1
        h=K%s.R
        pk=s.P[k]
        return K,(e*pk+m*(s.a[h]+s.b[h]*s.x))%s.r,pk,s.a[h],s.b[h]

def trial(r,p,R,s_esc,d,Mn,steps,taumax,rng):
    x=rng.randrange(1,r)
    c=Map(r,p,R,x,rng.getrandbits(64))
    W=[]
    for _ in range(Mn):
        A0=rng.randrange(1,r); B0=rng.randrange(1,r)
        W.append([(A0+B0*x)%r,A0,B0,0,0])
    tab={}; solved=False
    merges=splits_map=dp_mismatch=false_x=restarts=fp_collision=0
    pK=[None]*Mn; pR=[-1]*Mn; pT=[None]*Mn
    for t in range(steps):
        cK=[None]*Mn
        for i,w in enumerate(W):
            K,e1,pk,ah,bh=c.step(w[0],s_esc)
            m=2 if ctz(K,s_esc+8)>=s_esc else 1
            w[0]=e1
            w[1]=(pk*w[1]+m*ah)%r
            w[2]=(pk*w[2]+m*bh)%r
            w[3]+=1; cK[i]=K
            if dm0(K,d):
                w[3]=0
                Ac=pk*w[1]%r; Bc=pk*w[2]%r   # canonical representative exponents
                for (K0,A0,B0) in tab.get(K>>d,[]):
                    if K0==K:
                        if (A0,B0)!=(Ac,Bc):
                            den=(Bc-B0)%r
                            if den:
                                xh=(A0-Ac)*pow(den,r-2,r)%r
                                if xh==x: solved=True
                                else: false_x+=1
                    else:
                        fp_collision+=1
                tab.setdefault(K>>d,[]).append((K,Ac,Bc))
            if w[3]>=taumax:
                A0=rng.randrange(1,r); B0=rng.randrange(1,r)
                w[0]=(A0+B0*x)%r; w[1]=A0; w[2]=B0
                w[3]=0; w[4]+=1
                cK[i]=None; restarts+=1
        gmap={}
        for i,K in enumerate(pK):
            if K is not None: gmap.setdefault(K,[]).append(i)
        for K,ids in gmap.items():
            for ai in range(len(ids)):
                for bi in range(ai+1,len(ids)):
                    i,j=ids[ai],ids[bi]
                    merges+=1
                    if cK[i]!=cK[j]:
                        ri=W[i][4]>pR[i]; rj=W[j][4]>pR[j]
                        if ri or rj:
                            pass  # salted restart on >=1 side: allowed divergence
                        else:
                            splits_map+=1   # TRUE class-consistency violation
                    elif cK[i]==cK[j] and \
                        ((cK[i] is not None and dm0(cK[i],d))!=
                         (cK[j] is not None and dm0(cK[j],d))):
                        dp_mismatch+=1
        pK=cK; pR=[w[4] for w in W]; pT=[w[3] for w in W]
    return dict(solved=solved,merges=merges,splits_map=splits_map,
                dp_mismatch=dp_mismatch,false_x=false_x,restarts=restarts,
                fp_checks=fp_collision)

def dm0(K,d): return K%(1<<d)==0
def K0dp(K,d): return K%(1<<d)==0

def cycles(c,s_esc,samples=4000):
    r=c.r; out=[]
    for _ in range(samples):
        e=min([rng.randrange(1,r)*q%r for q in c.P])
        seen={}; ks=[]; t=0
        while e not in seen:
            seen[e]=t
            K,k=c.canon(e); ks.append(K)
            m=2 if (s_esc<40 and ctz(K,64)>=s_esc) else 1
            h=K%c.R
            e=(e*c.P[k]+m*(c.a[h]+c.b[h]*c.x))%r
            e=min([e*q%r for q in c.P]); t+=1
        out.append((t-seen[e],ks[seen[e]:]))
    return out

def Ppoint(n,k):
    num=1
    for j in range(1,k): num*=(n-j)
    return Fraction(num,n**k)

rng=random.Random(20261022)
out={"seed":20261022,"run":"toy_hard v3-hardened","configs":[]}
print("== C1 merge/DP audit (class-only ladder, salted tau-restart, lockstep) ==")
TOT=0
for (rmin,R,s_esc,d,Mn,steps,taumax,ntr) in [
        (61,4,2,2,6,400,24,800),
        (109,8,2,3,8,800,32,700),
        (397,8,2,4,10,1500,48,600),
        (1321,8,2,4,12,2000,64,150)]:
    r,p=params(rmin)
    agg=dict(solved=0,merges=0,splits_map=0,dp_mismatch=0,false_x=0,restarts=0,fp_checks=0)
    for _ in range(ntr):
        q=trial(r,p,R,s_esc,d,Mn,steps,taumax,rng); TOT+=1
        for kk in agg: agg[kk]+=q[kk]
    print(f"r={r} p={p} N={(r-1)//12} R={R} s={s_esc} d={d} M={Mn} steps={steps} "
          f"taumax={taumax} trials={ntr} solved={agg['solved']} mergepairs={agg['merges']} "
          f"nonrestart_splits={agg['splits_map']} dp_mismatch={agg['dp_mismatch']} "
          f"false_x={agg['false_x']} restarts={agg['restarts']} fp_collisions={agg['fp_checks']}")
    out["configs"].append(dict(r=r,p=p,N=(r-1)//12,R=R,s=s_esc,d=d,M=Mn,
        steps=steps,taumax=taumax,trials=ntr,**agg))
print("TOTAL_TRIALS",TOT)
out["total_trials"]=TOT
print("== C2 cycle distribution: ladder-off vs exact random-mapping, ladder-on max ==")
for rmin in (109,397):
    r,p=params(rmin); n=(r-1)//12
    c=Map(r,p,8,rng.randrange(1,r),rng.getrandbits(64))
    raw=cycles(c,40);   hard=cycles(c,2)
    rm=sum(L for L,_ in raw)/len(raw); hm=sum(L for L,_ in hard)/len(hard)
    print(f"r={r} N={n} ladderOFF_mean={rm:.2f} theory_sqrt(piN/8)={math.sqrt(math.pi*n/8):.2f}"
          f" ladderON_mean={hm:.2f} ladderON_max={max(L for L,_ in hard)}")
    for L in (1,2,3,4,8,16):
        print(f"  P(cyc={L}) emp={sum(1 for l,_ in raw if l==L)/len(raw):.5f} "
              f"exact={float(Ppoint(n,L)):.5f}")
    lad=[(L,ks) for L,ks in raw if all(ctz(K,64)<2 for K in ks)]
    wasted=sum(L for L,_ in lad)/len(raw)
    pred=sum(float(Ppoint(n,k))*k*Fraction(3,4)**k for k in range(1,140))
    print(f"  ladderless_cycles={len(lad)}/4000 wasted_steps/walk emp={wasted:.4f} "
          f"exact_pred(s=2)={pred:.4f}")
    out.setdefault("cycles",[]).append(dict(r=r,n=n,raw_mean=rm,hard_mean=hm,
        hard_max=max(L for L,_ in hard),
        emp={L:sum(1 for l,_ in raw if l==L)/len(raw) for L in (1,2,3,4,8,16)},
        theory={L:float(Ppoint(n,L)) for L in (1,2,3,4,8,16)},
        ladderless=len(lad),wasted=wasted,pred=pred))
nB=(52435875175126190479447740508185965837690552500527637822603658699938581184513-1)//12
eps=Fraction(0); num=1
for k in range(1,1025):
    if k>1: num*=(nB-(k-1))
    term=Fraction(num*k,nB**k)*Fraction(255,256)**k
    eps+=term
print("PROD expected wasted steps/walk from ladderless cycles (s=8):", f"{float(eps):.3e}")
out["prod_wasted_steps_s8"]=float(eps)
json.dump(out,open("toy_hard_results.json","w"),indent=1)
print("wrote toy_hard_results.json")
