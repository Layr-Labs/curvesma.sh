"""Exact v12 fixed-spacing ledger; no adaptive deepening or record eviction.
Joint retune: d=64, M=2^72 records, S=2^85 events, nu=137/100.
Abort on full; both Markov debits use the ACTUAL fixed coordinate-DP density.
C and R/C follow tracks/docs/COST_MODEL.md; memory is reported, not scored.
The finite toy tests do not prove the explicitly conditional production laws.
"""
import math
import json
from pathlib import Path
from fractions import Fraction as F
p=int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f',16)
r=int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141',16)
beta=int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee',16)
lam=int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72',16)
N=(r-1)//6
nu=F(137,100)
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
C,E0,E=1000,8000,9600
# Multiplication internal work is charged only for EXTRA beta calls; the six
# reference products are already in C. Other scalar-64 operations cost one.
instruction_primitives={'insn': 1, 'fp_mul': 179, 'normalize52_pack': 80, 'modsub52': 96, 'unpack52': 15, 'cmp256': 5}
itemized_step_rows=[{'name': 'affine_six_differences', 'routine': 'add/repair', 'primitive': 'modsub52', 'count': 6}, {'name': 'beta_images', 'routine': 'canonicalize', 'primitive': 'fp_mul', 'count': 2}, {'name': 'three_x_normalizations', 'routine': 'canonicalize', 'primitive': 'normalize52_pack', 'count': 3}, {'name': 'y_normalization', 'routine': 'canonicalize', 'primitive': 'normalize52_pack', 'count': 1}, {'name': 'unpack_live_xy', 'routine': 'add', 'primitive': 'unpack52', 'count': 2}, {'name': 'argmin_comparisons', 'routine': 'canonicalize', 'primitive': 'cmp256', 'count': 2}, {'name': 'argmin_word_selects', 'routine': 'canonicalize', 'primitive': 'insn', 'count': 8}, {'name': 'argmin_index_selects', 'routine': 'canonicalize', 'primitive': 'insn', 'count': 2}, {'name': 'y_parity_control', 'routine': 'canonicalize', 'primitive': 'insn', 'count': 3}, {'name': 'y_negation', 'routine': 'canonicalize', 'primitive': 'modsub52', 'count': 1}, {'name': 'y_orientation_select', 'routine': 'canonicalize', 'primitive': 'insn', 'count': 4}, {'name': 'branch_address', 'routine': 'step', 'primitive': 'insn', 'count': 3}, {'name': 'jump_xy_loads', 'routine': 'step', 'primitive': 'insn', 'count': 12}, {'name': 'DP64_load_compare_branch', 'routine': 'step', 'primitive': 'insn', 'count': 3}, {'name': 'eight_ring_equalities', 'routine': 'escape', 'primitive': 'insn', 'count': 72}, {'name': 'cycle_min_cmp_select', 'routine': 'escape', 'primitive': 'insn', 'count': 63}, {'name': 'ring_shift_load_store', 'routine': 'escape', 'primitive': 'insn', 'count': 56}, {'name': 'ring_append_load_store', 'routine': 'escape', 'primitive': 'insn', 'count': 8}, {'name': 'ring_clear_stores', 'routine': 'escape', 'primitive': 'insn', 'count': 32}, {'name': 'expanded_state_compare', 'routine': 'Brent', 'primitive': 'insn', 'count': 164}, {'name': 'checkpoint_copy', 'routine': 'Brent', 'primitive': 'insn', 'count': 82}, {'name': 'Brent_power_length_control', 'routine': 'Brent', 'primitive': 'insn', 'count': 12}, {'name': 'two_limb_length_control', 'routine': 'meters', 'primitive': 'insn', 'count': 24}, {'name': 'point_state_load_store', 'routine': 'step', 'primitive': 'insn', 'count': 16}, {'name': 'loop_and_dispatch', 'routine': 'step', 'primitive': 'insn', 'count': 16}, {'name': 'rewind_minimum_dispatch', 'routine': 'escape', 'primitive': 'insn', 'count': 16}, {'name': 'equal_x_identity_guards', 'routine': 'add/repair', 'primitive': 'insn', 'count': 18}, {'name': 'batch_queue_traffic', 'routine': 'add', 'primitive': 'insn', 'count': 32}, {'name': 'spill_reload_ceiling', 'routine': 'step', 'primitive': 'insn', 'count': 256}, {'name': 'scheduler_control', 'routine': 'step', 'primitive': 'insn', 'count': 64}, {'name': 'replay_tag_write', 'routine': 'step', 'primitive': 'insn', 'count': 4}, {'name': 'global_meter_control', 'routine': 'meters', 'primitive': 'insn', 'count': 8}, {'name': 'ring_flags', 'routine': 'escape', 'primitive': 'insn', 'count': 32}, {'name': 'pointers_addressing', 'routine': 'step', 'primitive': 'insn', 'count': 32}]
step_rows={q["name"]:q["count"]*instruction_primitives[q["primitive"]] for q in itemized_step_rows}
base_other=dict(dp_lookup_transfers=F(S*2**18,C),fresh_starts=F(S*4096),
    retention_maintenance=F(2**96),all_key_match_replays=F(2**110),
    exceptional_arithmetic=F(2**40),table_preprocessing=F(2**32),
    constant_checks_sampling_output=F(2**20))
rows={'main':F(K*(C+E),C)}
rows.update({k:F(6,5)*v for k,v in base_other.items()})
raw=sum(rows.values(),F(0))
W=F(27,25)*raw
claim=F(13068,100)
# Each full record retains both independent 256-bit starting labels. Extra
# length limb consumes padding, not a scalar limb; no entropy is discarded.
record_fields=dict(fingerprint=24,starting_scalars=64,length=16,
    pointers=36,flags=2,metadata=50)
record_bytes=sum(record_fields.values())
memory_rows=dict(node_pool=M*record_bytes,auxiliary_cap=2**30)
peak_memory=sum(memory_rows.values())
pool_policy='abort_on_full_no_evictions'
rho=F(7,10)
# Measured finite lower envelope; production transfer remains H-window.
_capture=json.loads((Path(__file__).parent/"sweep-summary.json").read_text())
_e=_capture["renewal_lower_envelope"]
eta=F(int(_e["numerator"]),int(_e["denominator"]))
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
    return eta*(1-upper)-sum(losses.values(),F(0))

def frac(x):
    x=F(x)
    return dict(numerator=str(x.numerator),denominator=str(x.denominator))

def log_display(x):
    x=F(x)
    return f'{math.log2(x.numerator)-math.log2(x.denominator):.10f}'
