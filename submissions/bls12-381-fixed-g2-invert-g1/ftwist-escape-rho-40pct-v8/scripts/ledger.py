#!/usr/bin/env python3
"""Exact integer/rational ledger for the orientation-guarded class walk.
Instruction prices are analytical ceilings, not measured full attack assembly.
"""
from fractions import Fraction as F
from math import isqrt
P=0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R=0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
N=(R-1)//12
BRANCHES=1024
NP=(N*BRANCHES+BRANCHES-2)//(BRANCHES-1)
RESERVE=F(109,100)
CLAIM=F(12881,100)
C=20000
FP_MUL=383
FP_ADD=64
SUCCESS_FLOOR=F(39,100)
NAMED_MODEL_MARGIN=F(1,100)


def ln_bounds(x,terms=80):
    x=F(x); assert 1<=x<=2
    z=(x-1)/(x+1)
    lo=2*sum(z**(2*k+1)/F(2*k+1) for k in range(terms))
    return lo,lo+2*z**(2*terms+1)/((2*terms+1)*(1-z*z))


def ceil_fraction(x):
    x=F(x); return (x.numerator+x.denominator-1)//x.denominator


def ceil_sqrt_fraction(x):
    x=F(x); n=isqrt(x.numerator//x.denominator)
    return n if n*n*x.denominator==x.numerator else n+1


def below_pow2(x,num,den):
    x=F(x)
    return pow(x.numerator,den)<(1<<num)*pow(x.denominator,den)


def interval(x,den=10000):
    x=F(x); lo=(x.numerator.bit_length()-x.denominator.bit_length()-1)*den; hi=lo+3*den
    while hi-lo>1:
        mid=(lo+hi)//2
        if below_pow2(x,mid,den):hi=mid
        else:lo=mid
    return [lo,hi,den]


def exp_minus_upper(u,terms=80):
    u=F(u);assert 0<u<1 and terms%2==0
    term=out=F(1)
    for j in range(1,terms+1):term*=-u/j;out+=term
    return out

LN_LO,LN_HI=ln_bounds(F(5,3))
T0=ceil_sqrt_fraction(2*NP*LN_HI)
T=ceil_fraction(RESERVE*T0)
ROWS={
 'first_key_two_Montgomery_decodes_2_times_383':2*FP_MUL,
 'first_key_four_Fp_products_4_times_383':4*FP_MUL,
 'first_key_nine_six_limb_add_sub_neg_9_times_64':9*FP_ADD,
 'first_key_minimum_twelve_six_limb_keys_index':600,
 'two_four_limb_scalar_add_reduce_updates':256,
 'twisted_multiplier_label_load_addressing':256,
 'DP_predicate_report_dispatch':128,
 'Floyd_key_comparison_control':256,
 'loop_hard_work_meters':128,
 'ordinary_state_scratch_traffic':512,
 'orientation_guard_difference_and_branch':128,
 'second_key_two_Montgomery_decodes_2_times_383':2*FP_MUL,
 'second_key_four_Fp_products_4_times_383':4*FP_MUL,
 'second_key_nine_six_limb_add_sub_neg_9_times_64':9*FP_ADD,
 'second_key_minimum_twelve_six_limb_keys_index':600,
 'guard_two_scalar_doublings_reductions':128,
 'guard_extra_state_traffic':256,
}
ORDINARY_CEILING=9000
# One R-add and at most one dense square per call. Three calls per Floyd
# iteration. Both metadata evaluations and all escape work are charged even
# on ordinary calls: no frequency assumption or escape-cap tail is required.
WEXEC=(6+F(3*ORDINARY_CEILING,C))*T+(1<<100)
WSTRESS=(6+F(6,5)*F(3*ORDINARY_CEILING,C))*T+(1<<101)
W=8*T+(1<<100)


def frobenius_rows(tower):
    rows=[]
    for k,constants in enumerate(tower.factors):
        mul=add=0;categories=[]
        for a,b in constants:
            if (a,b)==(1,0):cat='one'
            elif (a,b)==(P-1,0):cat='minus_one';add+=2
            elif b==0:cat='real';mul+=2
            elif a==0 and b in (1,P-1):cat='imaginary_unit';add+=1
            elif a==0:cat='imaginary';mul+=2;add+=1
            elif (a-b)%P==0 or (a+b)%P==0:cat='diagonal';mul+=2;add+=2
            else:cat='general';mul+=3;add+=4
            categories.append(cat)
        if k&1:add+=6
        traffic=0 if k==0 else 512;control=0 if k==0 else 128
        rows.append(dict(k=k,Fp_products=mul,Fp_add_sub_neg=add,traffic=traffic,
            control=control,instructions=mul*FP_MUL+add*FP_ADD+traffic+control,categories=categories))
    return rows


def ancillary():
    M=1<<69;G=129;record=128
    rows={'report_instructions':G*M*(1<<26),'clear_instructions':G*M*record*16,
          'candidate_units':(1<<32)*(1<<15),'setup_units':1<<50,
          'exception_pairing_output_units':1<<90}
    total=F(rows['report_instructions']+rows['clear_instructions'],C)+sum(
        rows[k] for k in ('candidate_units','setup_units','exception_pairing_output_units'))
    return rows,total

if __name__=='__main__':
    import json
    print(json.dumps(dict(T0=T0,T=T,nonreference_sum=sum(ROWS.values()),
        scored_log2=interval(W),executed_log2=interval(WEXEC),repriced_log2=interval(WSTRESS)),indent=2))
