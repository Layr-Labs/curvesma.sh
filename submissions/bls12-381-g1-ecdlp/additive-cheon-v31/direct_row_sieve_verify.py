#!/usr/bin/env python3
"""Exact finite direct-row sieve, under declared P1/endpoints; NOT a P2 proof."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
r=52435875175126190479447740508185965837690552500527637822603658699938581184513
d=1<<20;m=(r-1)//d;X=1<<82;X2=1<<70;Z=X*X2
n=67704887904593045852627984777072135;M=5893133993817406318491
B0=F(36*n**4,m);B=1<<128
def sqrt_bounds(v):
 v=F(v);q=isqrt(v.numerator*B*B//v.denominator)
 lo=F(q,B);hi=F(q+1,B)
 assert lo*lo<=v<hi*hi
 return lo,hi
S=isqrt(r);assert S*S<r<(S+1)*(S+1)
q3=isqrt(isqrt(3*B**4))
assert q3**4<=3*B**4<(q3+1)**4
root3lo=F(q3,B)
w=(1+sqrt_bounds(500*(d-1))[1])/d
Qbar=n*M*w
AB=F(57,100)*X;AG=F(64,100)*X
ZB=F(57,100)*Z;ZG=F(64,100)*Z
t0=(r-1)//(2*Z);z0=2*t0//X
tc=isqrt((r-1)//12)
assert 3*X*X<r and 3*X*X*X2<r and z0==1899369
def pi_exact(v):
 sieve=bytearray(b"\x01")*(v+1);sieve[:2]=b"\x00\x00"
 for p in range(2,isqrt(v)+1):
  if sieve[p]:sieve[p*p:v+1:p]=b"\x00"*(((v-p*p)//p)+1)
 return sum(sieve)
piz0=pi_exact(z0);assert piz0==141993
exception=F(24*n*M*(piz0+1),1)/B0
Asame=24*Qbar*(ZB*ZB+ZG*ZG)/r/B0
# Selected original primitive-vector CDF: g=1 event maps + n-weighted cross.
# The 16n/r cross area coefficient and 4n cross constant are conservative.
lead=24*(ZG*2*sqrt_bounds(Qbar*(8*n*AB/r+F(16*n,r)))[1]
        +ZB*2*sqrt_bounds(Qbar*(8*n*AG/r+F(16*n,r)))[1])/B0
content=24*(ZB+ZG)*4*n*(F(1,t0)+F(1,tc))/B0
# Choose an l1-shortest original vector. The diamond overlap argument
# gives ||v||_1<=sqrt(2r); since its norm is integer use floor(sqrt(2r)).
l1max=isqrt(2*r)
assert l1max*l1max<=2*r<(l1max+1)*(l1max+1)
row_direct=24*(ZB+ZG)*Qbar*l1max/r/B0
row_constant=48*Qbar/B0
primes=[2,3,5,7,11,13,17,19,23,29,31,37,41]
rows=[];Q=phi=1;prod=F(1);sum_ap=F();divs=[(1,0,F(1))]
for k,p in enumerate(primes,1):
 Q*=p;phi*=p-1
 ap=F(1) if p==2 else sqrt_bounds(F(2,p))[1]
 prod*=1+ap;sum_ap+=ap
 divs=divs+[(v*p,o+1,s*(1+F(1,p))) for v,o,s in divs]
 if k<2:continue
 assert Q<X//2
 # Check validity at all optimal integration endpoints for every divisor.
 for v,o,sig in divs[1:]:
  for A in (AB,AG):
   C=(8*n*A/r+F(16*n,r))*(1<<o)/v
   T2=Qbar/C
   assert T2>t0*t0 and T2<(F(X*X,8))**2
   assert 12*T2<r*X*X
 rho=F(phi,Q)
 main=Asame*rho*rho
 # Exact 2^b*1_(b>0)=2*sum_(odd subsets)1.
 # Original selected shortest direction is ordinary primitive, so
 # its allocated bad-event F has content g=1 (no sigma/tau inflation).
 # Since ap(2)=1, product(1-ap)=0, giving exact odd sum =prod(1+ap).
 # Sharper chord lower bound I/h-fmax gives HALF the previous loss:
 # 2^(b-1)*1_(b>0)=sum_(odd subsets)1.
 bad=lead*rho*prod/2
 bad+=content*rho*(3**k-(-1)**k)/2
 # Mod6: sum allowed row-parameter classes per period of 6 row indices is4.
 interval_factor=3**(k-2)
 rounding=interval_factor*(F(2,3)*row_direct+4*row_constant)+row_constant
 total=main+bad+rounding+exception
 rows.append(dict(pmax=p,number_primes=k,Q=Q,phi=phi,
    main_exact=str(main),main_decimal=float(main),
    bad_direction_exact=str(bad),bad_direction_decimal=float(bad),
    row_rounding_exact=str(rounding),row_rounding_decimal=float(rounding),
    short_exception_decimal=float(exception),total_exact=str(total),total_decimal=float(total)))
best=min(rows,key=lambda z:F(z["total_exact"]))
assert best["pmax"]==31
total_upper=F(best["total_exact"])
assert total_upper<F(156980,1000000)
out=dict(status="PASS_SHARED_SMALL_BOUND_ONLY_NOT_P2",
 assumptions="Accepted P1(delta_w<=500) BOTH sides, certified block endpoints .57/.64X1; actual prefixes; no H2.",
 scope="Exactly-one-shared-small, two same-side placements. Pruned inner placements and >=2-shared strata accounted separately.",
 proof_inputs="Direct-row chord inequalities and selected-vector divisor CDF proved in finite-sieve-alternative.txt.",
 P1_max_mass_upper_exact=str(w),
 short_exception_exact=str(exception),short_exception_decimal=float(exception),
 bad_direction_leading_coefficient_exact=str(lead),
 bad_direction_content_coefficient_exact=str(content),
 row_direct_coefficient_exact=str(row_direct),
 row_constant_coefficient_exact=str(row_constant),
 scan=rows,best=best,safe_bound="0.156980",
 missing="No bound on coprime excess Delta0; this shared-small estimate alone does not prove route B or original P2.")
Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(dict(status=out["status"],best=best,safe_bound=out["safe_bound"],
 scan=[{k:x[k] for k in ("pmax","main_decimal","bad_direction_decimal","row_rounding_decimal","total_decimal")} for x in rows]),indent=2))

