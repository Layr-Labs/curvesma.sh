#!/usr/bin/env python3
"""Exact charged-escape residual and a clearly labelled renewal diagnostic."""
from fractions import Fraction as F
from ledger import *

def frac(x):
    x=F(x);return dict(numerator=str(x.numerator),denominator=str(x.denominator))

def evaluate():
    h=F(1,12*BRANCHES)
    # Five geometric forced-basin segments with a deliberately inflated 4x
    # visit count. This is a compatibility countermodel only; it is decoupled
    # from escape stopping, which has no per-trail cap.
    k=5
    second=F(k*k+k*(1-h),1)/(h*h)
    before=min(F(1),8*second/N)
    marked=min(F(1),F(4*k,1)/(h*DP_SPACING))
    return dict(version=1,per_trail_escape_cap=None,
        global_priced_allowance=ESCAPE_CAP,
        expected_forced_opportunities=frac(ESCAPE_PRED),
        conditional_mean_ceiling=frac(ESCAPE_MEAN_CEILING),
        global_allowance_Markov=frac(ESCAPE_TAIL),
        applied_escape_tail=frac(STRICT_ESCAPE_TAIL),
        floor_arithmetic=frac(MODEL_LOWER),
        binding='Each escape is charged and the priced global allowance bounds the escape residual at 1/128; there is no per-trail escape cap to price.',
        renewal_diagnostic=dict(scope='Compatibility countermodel only; NOT a measured production distribution, coordinate-map theorem or universal lower bound.',
            hazard=frac(h),segments=k,visits_multiplier=4,
            useful_repeat_union=frac(before),any_DP_union=frac(marked),
            success_upper_union=frac(min(F(1),before+marked))),
        qualification=dict(global_allowance_tail_nonvacuous=STRICT_ESCAPE_TAIL<F(1,100),
            production_success_margin=MODEL_LOWER>=F(2,5)),
        floyd_only_preserving_emitted_records=dict(
            semantics='Keep the same charged in-place jump and pool; omit the Floyd stall-oracle rows only. This is repricing, not a built second variant.',
            row='Floyd_key_comparison_control=256 per call',
            cost_ceiling=frac(WEXEC-F(3*256*T,C)),
            qualified=False,
            reason='Dropping the stall oracle rows gives up exactly the DP-free trigger this package exists to repair, and loses cross-escape DP credit. Retaining all non-DP visits would require a separately specified and priced table.'))

if __name__=='__main__':
    import json
    print(json.dumps(evaluate(),separators=(',',':')),flush=True)
