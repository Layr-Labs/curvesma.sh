#!/usr/bin/env python3
"""Exact integer/rational ledger for raw R-adding with paid rare escapes.
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
CLAIM=F(12835,100)
E4_ABLATION_ESCAPE_CAP=4  # historical ablation only; never production semantics
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
    if num>=0:return pow(x.numerator,den)<(1<<num)*pow(x.denominator,den)
    return (pow(x.numerator,den)<<(-num))<pow(x.denominator,den)


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
 'DP_pool_key_comparison_control':256,
 'loop_hard_work_meters':128,
 'ordinary_state_scratch_traffic':512,
 'label_alignment_index_bookkeeping':128,
 'Floyd_key_comparison_control':256,
 'brent_expanded_state_compare_control':132,
 'checkpoint_full_copy':76,
}
ORDINARY_CEILING=5700
# The forced inverse/same-branch opportunity predicts 1/(12B). Complete
# zero-denominator Floyd meetings also include multilap/long-cycle effects.
# Hescape allows FOUR times that mean, not a false exact-rate assertion.
ESCAPE_PRED=F(3*T,12*BRANCHES)
ESCAPE_MEAN_CEILING=4*ESCAPE_PRED
ESCAPE_CAP=ceil_fraction(512*ESCAPE_PRED)
ESCAPE_TAIL=ESCAPE_MEAN_CEILING/ESCAPE_CAP  # global allowance Markov tail
# Escape accounting (repaired): a trail has NO per-trail escape cap. Each
# in-place offset is charged one dense product plus ESCAPE_INSTRUCTIONS
# (ESCAPE_ROWS, one_jump_charge), and termination is bounded by the ordinary
# trail cap 3T. ESCAPE_CAP is the PRICED global allowance, so Markov on the
# conditional first moment bounds the event "a trail exceeds the allowance"
# by ESCAPE_MEAN_CEILING/ESCAPE_CAP = 1/128. That global row is the only
# escape residual subtracted below; the former min(1,mean/5)=1 clipping row
# described an E=4 stopping rule this algorithm never had to use.
STRICT_ESCAPE_TAIL=ESCAPE_TAIL
ESCAPE_REF=1
ESCAPE_ROWS={
 'one_metadata_key':3474,
 'one_label_update_index_in_place':1024,
 'one_full_canonical_Frobenius_confirmation':72*FP_MUL+12000,
 'full_point_compare_fingerprint_and_dispatch':8192,
 'fresh_index_randomness_and_scratch':1024,
}
ESCAPE_INSTRUCTIONS=55000
# Fixed, production-long spacing and capacity; neither is adapted at runtime.
DP_SPACING=1<<54
POOL_CAPACITY=1<<80
RECORD_BYTES=24+32+32+24+16
FINGERPRINT_BITS=192
REPORT_CAP=1<<84
BALANCE=F(11,10)
K=3*T+ESCAPE_CAP+1
# Conditional first moment is for THIS spacing, not a coarser predecessor.
POOL_MARKOV=BALANCE*F(K,DP_SPACING*POOL_CAPACITY)
TRANSPORT=F(64*DP_SPACING*DP_SPACING,N)
FINGERPRINT_TAIL=F(REPORT_CAP*POOL_CAPACITY,(1<<FINGERPRINT_BITS)*(1<<32))
REPORT_MARKOV=BALANCE*F(K,DP_SPACING*REPORT_CAP)
MEASURED_RECOVERY=F(5661,10000)  # historical diagnostic, not v11 premise
CONDITIONAL_MODEL_REFERENCE=F(5661,10000)  # explicitly unqualified H-Rom
MODEL_LOWER=CONDITIONAL_MODEL_REFERENCE-STRICT_ESCAPE_TAIL-POOL_MARKOV-REPORT_MARKOV-TRANSPORT-FINGERPRINT_TAIL
# Pool full => abort. There is NO deletion, replacement or density change.
# P(any full-pool event) <= E[reports]/M under the NAMED report-count law.
# This is not a derivation of that law for a stalled history-dependent walk.
# Both unrepriced executed ceilings and repriced RAW rows are published.
# All Fp products/decodes are nonreference and receive the 6/5 repricing.
ESCAPE_BUDGET_EXEC=(ESCAPE_REF+F(ESCAPE_INSTRUCTIONS,C))*ESCAPE_CAP
ESCAPE_BUDGET_STRESS=(ESCAPE_REF+F(6,5)*F(ESCAPE_INSTRUCTIONS,C))*ESCAPE_CAP
WEXEC=(3+F(3*ORDINARY_CEILING,C))*T+ESCAPE_BUDGET_EXEC+(1<<100)
WSTRESS=(3+F(6,5)*F(3*ORDINARY_CEILING,C))*T+ESCAPE_BUDGET_STRESS+(1<<101)
W=F(29,5)*T+(1<<100)
# There is no separate strict per-trail execution row: the executed ceiling
# already prices ESCAPE_CAP charged escapes, which covers unbounded-per-trail
# escapes up to the priced allowance and fails closed above it.


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
    M=POOL_CAPACITY;G=1;record=RECORD_BYTES
    rows={'report_instructions':REPORT_CAP*(1<<26),'clear_instructions':G*M*(record+64)*16,
          'candidate_units':(1<<32)*(1<<15),'setup_units':1<<50,
          'exception_pairing_output_units':1<<90}
    total=F(rows['report_instructions']+rows['clear_instructions'],C)+sum(
        rows[k] for k in ('candidate_units','setup_units','exception_pairing_output_units'))
    return rows,total

if __name__=='__main__':
    import json
    print(json.dumps(dict(T0=T0,T=T,nonreference_sum=sum(ROWS.values()),
        scored_log2=interval(W),executed_log2=interval(WEXEC),repriced_log2=interval(WSTRESS))))
