#!/usr/bin/env python3
"""Published prime-count inputs plus exact rational consuming arithmetic.
Analytic input only: Dusart, arXiv:1002.0442v1, Thm6.9 eq(6.6).
No numerical test below is a proof of that published theorem.
"""
from fractions import Fraction as F
from pathlib import Path
import json

SOURCE="https://arxiv.org/pdf/1002.0442"
X1=2**82
X2=2**70
NA=11786267987634812636981
M=(NA+1)//2

def small_log_bounds(u,J=64):
    """1<=u<=2. Positive atanh series, with geometric tail majorant."""
    u=F(u)
    assert 1<=u<=2
    if u==1:return F(0),F(0)
    t=(u-1)/(u+1)
    lo=sum((2*t**(2*j+1)/F(2*j+1) for j in range(J)),F(0))
    rem=2*t**(2*J+1)/(F(2*J+1)*(1-t*t))
    return lo,lo+rem

LN2_LO,LN2_HI=small_log_bounds(2)
def log_bounds(x):
    x=F(x); assert x>0
    k=0
    while x>=2:x/=2;k+=1
    while x<1:x*=2;k-=1
    lo,hi=small_log_bounds(x)
    if k>=0:return k*LN2_LO+lo,k*LN2_HI+hi
    return k*LN2_HI+lo,k*LN2_LO+hi

def pi_lower(x):
    x=F(x)
    # Conservative strict versions of Dusart's lower validity threshold.
    assert x>5393
    lo,hi=log_bounds(x)
    assert hi>1
    return x/(hi-1)

def pi_upper(x):
    x=F(x)
    # Conservative strict version of the upper validity threshold.
    assert x>60184
    lo,hi=log_bounds(x)
    assert lo>F(11,10)
    return x/(lo-F(11,10))

def floorq(x):return x.numerator//x.denominator

def finite_consuming_checks():
    base=F(X1,2)
    assert NA%2==1 and 2*M==NA+1
    baby=pi_lower(F(57,100)*X1)-pi_upper(base)
    giant=pi_lower(F(16,25)*X1)-pi_upper(base)
    assert baby>M and giant>2*M
    total_prime_upper=floorq(pi_upper(X2))
    assert total_prime_upper==24896332478097058409
    return {
      "analytic_input":"Dusart Thm6.9 eq(6.6), valid ranges checked; not re-proved",
      "source":SOURCE,"X1":X1,"X2":X2,"NA":NA,"M":M,"two_M":2*M,
      "pi_X2_upper":total_prime_upper,
      "baby_endpoint":"57*X1/100","baby_rank_lower_floor":floorq(baby),
      "baby_rank_required":M,"baby_margin_floor":floorq(baby)-M,
      "giant_endpoint":"16*X1/25","giant_rank_lower_floor":floorq(giant),
      "giant_rank_required":2*M,"giant_margin_floor":floorq(giant)-2*M,
      "status":"PASS"
    }

def small_transcription_spotchecks():
    """Independent small sieve; checks transcription, not large-x validity."""
    N=2**20
    flags=bytearray(b'\x01')*(N+1)
    flags[0]=flags[1]=0
    for p in range(2,1025):
        if flags[p]:
            start=p*p
            flags[start:N+1:p]=b'\x00'*(((N-start)//p)+1)
    rows=[]
    for x in [2**13,2**16,2**18,2**20]:
        exact=sum(flags[:x+1])
        lo,hi=log_bounds(x)
        # Opposite enclosure directions certify the exact published formulas
        # at these checkpoints, rather than only their relaxed consuming bounds.
        assert F(x)/(lo-1)<=exact
        upper_pass=None
        if x>60184:
            assert exact<=F(x)/(hi-F(11,10))
            upper_pass=True
        rows.append({"x":x,"pi_exact_by_sieve":exact,"lower_pass":True,
                     "upper_pass":upper_pass})
    return {"scope":"small-scale constant/transcription checks only; not production prime-count verification",
            "enclosure_direction":"opposite enclosures certify the exact stated formulas at checkpoints",
            "rows":rows,"status":"PASS"}
