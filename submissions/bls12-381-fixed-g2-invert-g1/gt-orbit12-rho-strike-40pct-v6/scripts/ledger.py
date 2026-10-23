#!/usr/bin/env python3
"""Exact rational logarithm enclosures and the strike ledger; no float bounds."""
from fractions import Fraction as F
from math import isqrt
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
N=(R-1)//12
BRANCHES=1024
NP=(N*BRANCHES+BRANCHES-2)//(BRANCHES-1)
RESERVE=F(11,10)
CLAIM=F(6441,50)  # 128.82, a HOLD candidate, never a sub-frontier assertion.

def ln_bounds(x, terms=80):
    x=F(x)
    assert 1<=x<=2
    z=(x-1)/(x+1)
    lower=2*sum((z**(2*k+1)/F(2*k+1)) for k in range(terms))
    upper=lower+2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return lower,upper

def ceil_fraction(x):
    return (x.numerator+x.denominator-1)//x.denominator

def ceil_sqrt_fraction(x):
    x=F(x); a=isqrt(x.numerator//x.denominator)
    return a if a*a*x.denominator==x.numerator else a+1

def exp_minus_upper(x, terms=80):
    x=F(x); assert 0<x<1 and terms%2==0
    term=total=F(1)
    for j in range(1,terms+1):
        term*=-x/j
        total+=term
    return total

LN_LO,LN_HI=ln_bounds(F(5,3))
T0=ceil_sqrt_fraction(2*NP*LN_HI)
T=ceil_fraction(RESERVE*T0)
ORDINARY_ROWS={
 'four_base_field_products_4_times_383':1532,
 'nine_six_limb_add_sub_neg_9_times_64':576,
 'minimum_of_twelve_six_limb_keys_and_index':600,
 'two_scalar_double_add_reduce_updates':384,
 'twisted_jump_and_label_load_addressing':256,
 'DP_predicate_and_report_dispatch':128,
 'Floyd_key_comparison_and_control':256,
 'loop_and_meter_updates':128,
 'ordinary_state_load_store_and_scratch':512,
}
ORDINARY_CEILING=8000
WEXEC=(6+F(3*ORDINARY_CEILING,20000))*T+(1<<100)
W=8*T+(1<<100)
WSTRESS=(6+F(6,5)*F(3*ORDINARY_CEILING,20000))*T+(1<<101)

def below_pow2(x,num,den):
    x=F(x)
    return pow(x.numerator,den)<(1<<num)*pow(x.denominator,den)

def interval(x,den=10000):
    x=F(x)
    # Narrow display interval by integer-power binary search only.
    lo=(x.numerator.bit_length()-x.denominator.bit_length()-1)*den
    hi=lo+3*den
    while hi-lo>1:
        mid=(lo+hi)//2
        if below_pow2(x,mid,den):hi=mid
        else:lo=mid
    return [lo,hi,den]
