"""Exact integer production ledger; executed only by the isolated runner."""
import math
from decimal import Decimal, getcontext

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
D = 1 << 20
C = 2000
LANES = 1 << 16
SLOTS = 1 << 36
RECORD_CAP = 1 << 34
ENDPOINT_CAP = 1 << 34
CHAIN_CAP = 1 << 96
DP_BITS = 88


def compute():
    q = (R-1)//(6*D)
    assert q*6*D == R-1
    ceil_root = math.isqrt(q)+1
    states = LANES*((4*ceil_root+LANES-1)//LANES)
    per_state = {'seven_reference_affine_additions':7*C,
                 'two_scalar_Montgomery_calls':2*512,
                 'SHA256_one_block':8192,
                 'seven_batch_inversion_shares':7*(1 << 20)//LANES,
                 'all_table_gathers_staging_key_counters_and_control':8192}
    assert sum(per_state.values()) < 32768
    # Cold setup, all exceptional work, endpoint handling, restarts, root recovery
    # and output have explicit component bounds in proof.md totalling below 2^64.
    instructions = states*32768+(1 << 64)
    point_bytes = 2*8*(1 << 32)*96
    record_bytes = SLOTS*128
    scratch_bytes = 1 << 38
    memory = point_bytes+record_bytes+scratch_bytes
    assert memory < 1 << 44
    assert instructions < C*(1 << 123)
    getcontext().prec = 80
    logwork = (Decimal(instructions)/C).ln()/Decimal(2).ln()
    theta = 2.0**-DP_BITS
    k = states-LANES*CHAIN_CAP
    cross_lower = 1-math.exp(-k*k/(16*q))
    duplicate_bound = 64/(1 << 20)
    cycle_bound = 4*(ENDPOINT_CAP+LANES)/(q*theta*theta)
    probe_bound = SLOTS*(math.e/4)**128/(1-math.e/4)
    assert cross_lower > 0.61 and cycle_bound < 2**-20 and probe_bound < 2**-30
    return {'n_minimum':D,'n_scored':D,'d':D,'quotient_size':str(q),
            'maximum_states':str(states),'lanes':LANES,'dp_bits':DP_BITS,
            'chain_cap':str(CHAIN_CAP),'record_cap':RECORD_CAP,'endpoint_cap':ENDPOINT_CAP,
            'hash_slots':SLOTS,'per_state_instruction_components':per_state,
            'per_state_instruction_budget':32768,'other_instruction_budget':str(1 << 64),
            'total_instruction_bound':str(instructions),'C':C,
            'work_log2_diagnostic':str(logwork),'submitted_log2_bound':'123',
            'score':{'numerator':'123','denominator':'1'},
            'memory_bytes':str(memory),'peak_memory_log2_bound':44,
            'preprocessing_log2_bound':54,'advice_bytes':0,
            'success_claim':0.40,'success_is_heuristic':True,
            'conditional_success_diagnostics':{'Hwalk_encounter_lower':cross_lower,
              'duplicate_endpoint_Markov_upper':duplicate_bound,
              'Hwalk_intra_trail_repeat_union_upper':cycle_bound,
              'ideal_index_probe_union_upper':probe_bound,
              'reserved_combined_tail_failure':0.01},
            'scope':'Exact ledger arithmetic; primitive budgets and production heuristic need proof.md review.'}
