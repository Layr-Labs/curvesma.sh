"""Exact rational enclosure of a declared logarithm, no floating point."""
from fractions import Fraction as F
import json

def ln_interval(x, terms=64):
    assert F(1)<=x<=F(2)
    z=(x-1)/(x+1)
    lower=2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
    tail=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lower,lower+tail

def certify(n,d,ceiling):
    x=F(n,d);k=n.bit_length()-d.bit_length()
    while x<2**k:k-=1
    while x>=2**(k+1):k+=1
    l,u=ln_interval(x/F(2**k));a,b=ln_interval(F(2))
    low=F(k)+l/b;high=F(k)+u/a
    assert high<ceiling
    assert low>ceiling-F(1,10**7)
    return {'declared_ceiling':str(ceiling),'certified':True,'width_below_1e_minus_50':high-low<F(1,10**50),'terms':64,'scope':'Certifies only the logarithm of the supplied total, not the algorithm or heuristics.'}
if __name__=='__main__':
 print(json.dumps(certify(201409072310421922867463303485386601261,2000,F(116277611990909175,10**15)),indent=2))
