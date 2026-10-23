"""Exact v10b fixed-spacing ledger; no adaptive deepening or record eviction.
Joint retune: d=64, M=2^72 records, S=2^85 events, nu=153/125.
Abort on full; both Markov debits use the ACTUAL fixed coordinate-DP density.
C and R/C follow tracks/docs/COST_MODEL.md; memory is reported, not scored.
The finite toy tests do not prove the explicitly conditional production laws.
"""
import math
from fractions import Fraction as F
p=int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f',16)
r=int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141',16)
beta=int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee',16)
lam=int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72',16)
N=(r-1)//6
nu=F(153,125)
c=F(98859,100000)
budget_square=c*nu*nu*N
T=math.isqrt(budget_square.numerator//budget_square.denominator)
if T*T<budget_square: T+=1
M,S,A,workers=2**72,2**85,2**24,1024
R=2**20
d=64
DP=2**d
branch_shift=64
L=40*DP                 # actual fixed-spacing trail cutoff, not an old depth
K=T+workers*L
C,E0,E=1000,950,1140
# Twenty-four newly itemized instructions cover two-limb length increments,
# limits, checkpoint copies and replay bookkeeping. 949 fits the SAME 950 cap.
step_rows=dict(beta_multiplications=358,field_normalizations=160,canonical_argmin=36,
    y_parity_negation=20,adding_table_index_loads=16,distinguished_test=3,
    full_ring_comparisons=96,state_and_loop=28,two_limb_length_control=24,
    brent_expanded_state_compare_control=132,checkpoint_full_copy=76)
base_other=dict(dp_lookup_transfers=F(S*2**18,C),fresh_starts=F(S*4096),
    retention_maintenance=F(2**96),all_key_match_replays=F(2**110),
    exceptional_arithmetic=F(2**40),table_preprocessing=F(2**32),
    constant_checks_sampling_output=F(2**20))
rows={'main':F(K*(C+E),C)}
rows.update({k:F(6,5)*v for k,v in base_other.items()})
raw=sum(rows.values(),F(0))
W=F(27,25)*raw
claim=F(12821,100)
# Each full record retains both independent 256-bit starting labels. Extra
# length limb consumes padding, not a scalar limb; no entropy is discarded.
record_fields=dict(fingerprint=24,starting_scalars=64,length=16,
    pointers=36,flags=2,metadata=50)
record_bytes=sum(record_fields.values())
memory_rows=dict(node_pool=M*record_bytes,auxiliary_cap=2**30)
peak_memory=sum(memory_rows.values())
pool_policy='abort_on_full_no_evictions'
rho=F(7,10)
velocity=1/(1+F(8,R))
renewal_core=F(3,4)*velocity**2
balance=F(11,10)
# H-rate: E[D] <= balance*K/DP for the fixed-density event process.
# At most D records can be live. If D >= M, fail BEFORE deleting a record.
# Markov therefore gives P(pool abort) <= E[D]/M = balance*fill.
fill=F(K,DP*M)
count_pool=balance*fill
# Every start after the initial workers consumes a terminal event; the same
# fixed-spacing first-moment law gives E[D]/(S-workers), never depth 56/58.
count_events_starts=balance*F(K,DP*(S-workers))
false_key=F(S*(S-1),2*2**192)
relation=F(1,r)
reject=F(2**256-r,2**256)
sampling=4*(S+R)*reject**4
exception=F(4*R*K,N*2**16)
# Four <=L pieces: <=8L^2 pairs. Recomputed when d and L move together.
detector=F(8*L*L,N)
# H-window selected-pair transport debit, NOT a record-replacement bound.
transport=F(64*DP*DP,N)
trail_cap=F(44,2**24)
losses=dict(pool_markov=count_pool,event_start_markov=count_events_starts,
    fingerprint_union=false_key,relation_failure=relation,sampling_rejections=sampling,
    exception_cap=exception,unforced_detector_abort=detector,
    endpoint_transport=transport,collision_transport_trail_cap=trail_cap)

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
