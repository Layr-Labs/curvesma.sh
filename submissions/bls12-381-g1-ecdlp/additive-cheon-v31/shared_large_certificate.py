#!/usr/bin/env python3
"""Standalone exact certificate for the shared-large theorem.

Analytic input (not proved by this program): Dusart2010, Theorem6.9,
equation6.6: pi(x)>=x/(ln x-1), x>5393; pi(x)<=x/(ln x-1.1),
x>60184. Every use/domain and consuming arithmetic is exposed below.
Class-weight arithmetic premise: delta_B^w,delta_G^w<=500.
All other finite inputs, including primality of r, are checked here.
"""
from fractions import Fraction as F
from math import isqrt, gcd
from pathlib import Path
import json

def prime_by_trial(q):
    if q<2:return False
    if q==2:return True
    if q%2==0:return False
    return all(q%a for a in range(3,isqrt(q)+1,2))

# Elementary Pocklington witness criterion, proved in the companion text.
z=0xd201000000010000
L=z*z-1
r=z**4-z*z+1
assert r==52435875175126190479447740508185965837690552500527637822603658699938581184513
assert z==2**16*906349*254760293
factor_F={2:32,906349:2,254760293:2}
for p in factor_F:assert prime_by_trial(p)
known_F=1
for p,e in factor_F.items():known_F*=p**e
assert known_F==z*z and (r-1)%known_F==0 and known_F*known_F>r
witnesses={2:5,906349:2,254760293:2}
primality_witnesses=[]
for p,a in witnesses.items():
    assert pow(a,r-1,r)==1
    residue=pow(a,(r-1)//p,r)
    assert gcd(residue-1,r)==1
    primality_witnesses.append({'prime_factor':p,'known_exponent':factor_F[p],
      'base':a,'power_residue':residue,'gcd':gcd(residue-1,r)})
# The criterion implies r is prime; no probabilistic primality test is used.

d=2**20
assert (r-1)%d==0
m=(r-1)//d
assert m%6==0
n=isqrt(11*m//120)
if 120*n*n<11*m:n+=1
assert 120*(n-1)**2<11*m<=120*n*n
X1,X2=2**82,2**70
NA=196*n*d//X2+2
M=(NA+1)//2
assert n==67704887904593045852627984777072135
assert M==5893133993817406318491
assert X2<X1//2 and 2*X1*X2*X2<r and 2*X1*X1*X2<r
assert 3*X1*X1<r
Z=X1*X2
baseline=F(36*n**4,m)

# Check the explicit root-lattice bases and centered-cell radius.
assert L*L+L+1==r
Rcell=(L+3)//2
bases=[
 (L,(L,1),(-1,L+1)),
 (-L,(-L,1),(1,L+1)),
 (L+1,(L+1,1),(-1,L)),
 (-L-1,(-L-1,1),(1,L)),
]
for rho,v,w in bases:
    assert (v[0]-rho*v[1])%r==0 and (w[0]-rho*w[1])%r==0
    assert abs(v[0]*w[1]-v[1]*w[0])==r
    assert all(abs(v[j])+abs(w[j])<=2*Rcell for j in [0,1])
    assert any((rho*rho+s*rho+1)%r==0 for s in [1,-1])

# Preserve the exact rounding choices of the original numerical certificate.
ln2lo=F(693147180559945309,10**18)
ln2hi=F(693148,10**6)
s=sum((F(2,(2*j+1)*3**(2*j+1)) for j in range(40)),F(0))
tail=F(2,81*3**81)*F(9,8)
assert ln2lo<s and s+tail<ln2hi
def log_ratio_bounds(u,terms=64):
    u=F(u);assert u>=1
    t=(u-1)/(u+1)
    lo=sum((2*t**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    hi=lo+2*t**(2*terms+1)/(F(2*terms+1)*(1-t*t))
    return lo,hi
def sqrt_up(v,Q=10**100):
    v=F(v)
    u=F(isqrt(v.numerator*Q*Q//v.denominator)+1,Q)
    assert u*u>=v
    return u

# Complete finite checks for the cutoff-free mixed-root k=0 certificate.
def v2(q):
    e=0
    while q%2==0:q//=2;e+=1
    return e
Uints=[1,L,L*L,L+1,L*L+1,L*(L+1)]
assert len(set(Uints))==6 and all(0<u<r for u in Uints)
assert {u%r for u in Uints}=={1,r-1,L%r,(-L)%r,(L+1)%r,(-L-1)%r}
assert v2(L+1)==v2(L*(L+1))==32
assert all(v2(u)==0 for u in [1,L,L*L]) and v2(L*L+1)==1
small_factor=2*5*37*1637
assert all(prime_by_trial(p) for p in [2,5,37,1637])
assert (L*L+1)%small_factor==0 and (L*L+1)//small_factor>1
assert gcd(L*L+1,L)==1

# The only prime-counting analytic input; the program verifies its domains
# and all substitutions/arithmetic, not the theorem's proof.
assert X2>60184
P=F(X2,70*ln2lo-F(11,10))   # outward upper for pi(X2)
Bexp=55
B0=2**Bexp
assert B0>60184
Psmall=F(B0,Bexp*ln2lo-F(11,10))
endpoint_certificates=[]
for alpha,needed in [(F(57,100),M),(F(16,25),2*M)]:
    endpoint=alpha*X1
    assert endpoint>5393 and X1//2>60184
    ulo,uhi=log_ratio_bounds(2*alpha)
    count_lower=endpoint/(81*ln2hi+uhi-1)-F(X1,2)/(81*ln2lo-F(11,10))
    assert count_lower>needed
    endpoint_certificates.append({'coefficient':str(alpha),'endpoint':str(endpoint),
      'required_primes_above_X1_over_2':needed,
      'lower_count_floor':count_lower.numerator//count_lower.denominator})
ZB=F(57,100)*Z
ZG=F(16,25)*Z

mu=sqrt_up(F(501,d))
Qbar=2*n*mu*P
T=isqrt((r-1)*B0*B0//(3*X2*X2))
L0=isqrt(r//(3*X2*X2))
assert L0>0
assert 3*X2*X2*T*T<r*B0*B0
assert 3*X2*X2*(T+1)*(T+1)>=r*B0*B0
assert X2*T<B0**4
assert T<2**56*L0
C1=6*n*mu*X2*sqrt_up(F(2,r))
C0=12*n*mu*Rcell*sqrt_up(F(2,r))
Ibar=Qbar/T+C1*(56*ln2hi)+C0/L0
Qwheel=30
unit_residues=[a for a in range(Qwheel) if gcd(a,Qwheel)==1]
assert len(unit_residues)==8
rho_w=F(16,225)
sigma_w=F(32,15)
tau_w=64

Eretained=24*(rho_w*(Qbar/2)*(ZB*ZB+ZG*ZG)/r+
  sigma_w*(Ibar/2)*(ZB+ZG)+
  sigma_w*Qbar*(ZB+ZG)/isqrt(r)+tau_w*Qbar)
DB=2*ZB*ZB//r+1
DG=2*ZG*ZG//r+1
Elow=48*n*Psmall*(DB+DG)
Ewheel_exception=192*n*M*3
Etotal=Eretained+Elow+Ewheel_exception
ratio=Etotal/baseline
# The former decimal is a rounded display, slightly below the exact rational.
assert ratio>F('0.002332911916573674')
assert ratio<F('0.002332911916573675')
assert ratio<F('0.002333')

components={
 'coordinate_wheel_area':24*rho_w*(Qbar/2)*(ZB*ZB+ZG*ZG)/r/baseline,
 'inverse_minimum_boundary':24*sigma_w*(Ibar/2)*(ZB+ZG)/baseline,
 'direct_minimum_boundary':24*sigma_w*Qbar*(ZB+ZG)/isqrt(r)/baseline,
 'finite_constants':24*tau_w*Qbar/baseline,
 'low_coefficient_split':Elow/baseline,
 'coordinate_wheel_exceptions':Ewheel_exception/baseline,
}
assert sum(components.values(),F(0))==ratio

# Explicit completeness correction for two common large primes and no
# common small prime. This is separate from the exactly-one-large theorem.
Y=X2*X2
Dtwo=2*Y*Y//r+1
Etwo_large=F(96*M*M*Dtwo,1)/baseline
assert Etwo_large<F('0.0000000000000000164')
pruning_mass_over_m=F(6*M*n,m)
Erepeat_ratio=F(48*n*(2*Z*Z//r+1),1)/baseline
assert Erepeat_ratio<F('0.0000000000000000002671')
