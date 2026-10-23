#!/usr/bin/env python3
"""Exact Cheon checks; --toy runs real EC grids, not a cost benchmark."""
import json,sys,random
from math import isqrt,prod

def cs(n):return isqrt(n)+(isqrt(n)**2<n)
def prime(n):return n>=2 and all(n%j for j in range(2,isqrt(n)+1))
def grid(P,A,b,m,h,r,mul):
    babies={}; Q=P
    for j in range(m):babies[Q]=j; Q=mul(b,Q)
    Q=A; f=pow(b,-m,r)
    for i in range(h):
        if Q in babies:return i*m+babies[Q]
        Q=mul(f,Q)
    return None

def curve():
    p=884600761357; r=193863853
    def add(A,B):
        if A is None:return B
        if B is None:return A
        x,y=A; u,v=B
        if x==u and (y+v)%p==0:return None
        a=(3*x*x*pow(2*y,-1,p) if A==B else (v-y)*pow(u-x,-1,p))%p
        z=(a*a-x-u)%p
        return z,(a*(x-z)-y)%p
    def mul(k,A):
        assert k>=0
        B=None
        while k:
            if k&1:B=add(B,A)
            A=add(A,A); k>>=1
        return B
    return p,r,(241242499517,755159451635),mul

def test():
    p,r,P,mul=curve()
    assert prime(p) and prime(r)
    assert mul(0,P) is None and mul(r,P) is None and mul(r+1,P)==P
    assert P[1]**2%p==(P[0]**3+3)%p
    beta=next(pow(a,(p-1)//3,p) for a in range(2,99) if pow(a,(p-1)//3,p)!=1)
    orbit={(P[0]*pow(beta,i,p)%p,j*P[1]%p) for i in range(3) for j in (-1,1)}
    assert len(orbit)==6 and all(mul(r,Q) is None for Q in orbit)
    for k in range(6):
        assert grid(1,pow(4,k,13),4,2,3,13,lambda a,P:a*P%13)==k
    assert grid(1,pow(4,5,13),4,2,2,13,lambda a,P:a*P%13) is None
    return beta

def toy():
    beta=test(); p,r,P,mul=curve(); d=24367; q=(r-1)//d; g=2
    assert all(pow(g,(r-1)//l,r)!=1 for l in (2,3,7,13,17,59))
    m=(27*cs((2*q+4)//5)+24)//25; t=cs(d); rng=random.Random(20261023); hits=0
    for _ in range(2000):
        a=rng.randrange(d); b=rng.randrange(q); s=pow(g,q*a+b,r)
        T=mul(s,P); A=mul(pow(s,d,r),P)
        k=grid(P,A,pow(g,d,r),m,m,r,mul)
        assert (k is not None)==(b<m*m)
        if k is None:continue
        assert k==b
        j=grid(P,mul(pow(g,-k,r),T),pow(g,q,r),t,t,r,mul)
        assert j is not None and pow(g,q*(j%d)+k,r)==s
        hits+=1
    return dict(p=p,r=r,b=3,P=list(P),beta=beta,d=d,q=q,m=m,m2=t,trials=2000,
                successes=hits,wrong_outputs=0,cycles=0,escapes=0,seed=20261023,
                coverage=[m*m,q],grid_action_cap=2*m+2*t)

def certificate():
    r=int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16)
    fs=[(2,32),(3,1),(11,1),(19,1),(10177,1),(125527,1),(859267,1),(906349,2),
        (2508409,1),(2529403,1),(52437899,1),(254760293,2)]
    d=2**20; q=(r-1)//d; m0=cs((2*q+4)//5); m=(27*m0+24)//25; N=2*m+2048
    E=510*(64*2048+256*64+4096)+1024*2048+2**18
    route=176*129*256; R=N*(E+route+2**16)+2**80+2**20*d
    L=1<<(N-1).bit_length(); peak=288*L+2**64
    checks=dict(factors=prod(l**e for l,e in fs)==r-1 and all(prime(l) for l,e in fs),
      primitive=all(pow(7,(r-1)//l,r)!=1 for l,e in fs),divisor=d*q==r-1,
      coverage=m*m<q and 625*5*m*m>=729*2*q and 200*m*m>=81*q,
      limbs=r<2**255 and (r-1)**2<2**510 and (r-1)**2+(2**256-1)*r<2**512,
      memory=peak<2**129 and N<2**128,
      claim=R<2000*2**134,stress=6*R<5*2000*2**134)
    assert all(checks.values())
    return dict(checks=checks,d=d,q=q,m0=m0,m=m,N=N,E=E,route=route,R=R,C=2000,
       L=L,peak_bytes=peak,hard_cap_bytes=2**129,scalar_limbs=4,product_limbs=8,
       success=[m*m,q],claim=[134,1],stress_work=[6*R,10000],scalar_scratch_limbs=9)

if __name__=='__main__':
    if '--toy' in sys.argv:o=toy()
    elif '--test' in sys.argv:test(); o={'regression':'PASS'}
    else:o=certificate()
    print(json.dumps(o,sort_keys=True,separators=(',',':')))
