"""Exact full-production constants and inherited v8 ledger; only nu changes.
No production search. Fractions are used for every scored operation.
"""
import math
from fractions import Fraction as F
p=int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f',16)
r=int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141',16)
beta=int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee',16)
lam=int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72',16)
N=(r-1)//6
nu=F(2)
c=F(98859,100000)
budget_square=c*nu*nu*N
T=math.isqrt(budget_square.numerator//budget_square.denominator)
if T*T<budget_square: T+=1
M,S,A,workers=2**73,2**74,2**24,1024
L=40*2**58
K=T+workers*L
R=2**20
C,E0,E=1000,950,1140
step_rows=dict(beta_multiplications=358,field_normalizations=160,canonical_argmin=36,
    y_parity_negation=20,adding_table_index_loads=16,distinguished_test=3,
    full_ring_comparisons=96,state_and_loop=28,
    brent_expanded_state_compare_control=132,checkpoint_full_copy=76)
base_other=dict(dp_lookup_transfers=F(S*2**18,C),fresh_starts=F(S*4096),
    deepen_rebuilds=F(2**96),all_key_match_replays=F(2**110),
    exceptional_arithmetic=F(2**40),table_preprocessing=F(2**32),
    constant_checks_sampling_output=F(2**20))
rows={'main':F(K*(C+E),C)}
rows.update({k:F(6,5)*v for k,v in base_other.items()})
raw=sum(rows.values(),F(0))
W=F(27,25)*raw
claim=F(2579,20)
rho=F(7,10)
# Escape-renewal calibration and exact Markov, not a Poisson model.
velocity=1/(1+F(8,R))
renewal_core=F(3,4)*velocity**2
balance=F(11,10)
count_pool=balance*F(K,2**58*M)
count_events_starts=balance*F(K,2**56*(S-workers))
false_key=F(S*(S-1),2*2**192)
relation=F(1,r)
reject=F(2**256-r,2**256)
sampling=4*(S+R)*reject**4
exception=F(4*R*K,N*2**16)
# Four <=L pieces (two prefixes and two suffixes): <=8L^2 unordered pairs.
detector=F(8*L*L,N)
# At production v_R>0.99 and productive DP density >=0.9/s, a 20s
# window has exponent >17. e^-17<2^-24; two ages and two residuals.
trail_cap=F(44,2**24)
losses=dict(pool_markov=count_pool,event_start_markov=count_events_starts,
    fingerprint_union=false_key,relation_failure=relation,sampling_rejections=sampling,
    exception_cap=exception,unforced_detector_abort=detector,
    collision_transport_trail_cap=trail_cap)

def exp_minus_bounds(a):
    """Directed input enclosure and alternating rational series (0<=a<3)."""
    assert 0<=a<3
    q=10**12
    lo=F(a.numerator*q//a.denominator,q); hi=lo+F(1,q)
    def partial(x,n):
        term=total=F(1)
        for k in range(1,n+1): term*=-x/k; total+=term
        return total
    return partial(hi,61),partial(lo,60)

def success(quality=rho):
    _,upper=exp_minus_bounds(quality*T*(T-1)/(2*N))
    return 1-upper-sum(losses.values(),F(0))

def frac(x):
    x=F(x)
    return dict(numerator=str(x.numerator),denominator=str(x.denominator))

def log_display(x):
    x=F(x)
    return f'{math.log2(x.numerator)-math.log2(x.denominator):.10f}'
