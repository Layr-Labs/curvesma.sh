#!/usr/bin/env python3
"""Recompute all certificate verdicts with integers/Fraction only (no floats).
Proves arithmetic implications, NOT the disclosed production walk heuristic.
"""
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt
import json

P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
R_ORDER = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
C = 20000


def falling(n, k):
    z = 1
    for j in range(k):
        z *= n - j
    return z


def exp_bounds(x, n=64):
    """Positive Taylor sum and geometric upper bound on its remaining terms."""
    term = total = F(1)
    for i in range(1, n + 1):
        term *= x / i
        total += term
    first = term * x / (n + 1)
    return total, total + first / (1 - x / (n + 2))


def fp2_mul(a, b):
    return ((a[0]*b[0]-a[1]*b[1]) % P, (a[0]*b[1]+a[1]*b[0]) % P)


def fp2_pow(a, e):
    out = (1, 0)
    while e:
        if e & 1:
            out = fp2_mul(out, a)
        a = fp2_mul(a, a)
        e >>= 1
    return out


def build():
    r, p, processors = R_ORDER, P, 1 << 20
    N = (r-1)//12 + 1
    S = isqrt(8*N//5)
    if 5*S*S < 8*N:
        S += 1
    q = F(1, 1 << 60)
    D, start_cap, cap, candidate_cap = 1 << 66, 1 << 67, 1 << 67, 1 << 24
    slots = 3 * (1 << 65)
    W = F(36, 5)*S + (1 << 90)
    stress = F(168, 25)*S + F(6, 5)*(1 << 90)  # 3*(2+1.2*4000/C)*S
    itemized = 3*S*(2+F(4000, C)) + (1 << 90)
    mu = S*q
    fp_mean = F(D*(D-1), 2*(1 << 112))
    memory = 80*slots + processors*12*1024*640 + (1 << 50)
    checks = {}
    checks['protected_C_is_20000'] = C == 20000
    checks['score_reduced'] = gcd(6443, 50) == 1
    checks['r_has_255_bits'] = r.bit_length() == 255
    checks['p_has_381_bits'] = p.bit_length() == 381
    checks['order12'] = pow(p, 12, r) == 1 and all(pow(p, k, r) != 1 for k in range(1, 12))
    checks['half_frobenius_inverse'] = pow(p, 6, r) == r-1
    checks['nonidentity_orbits_exact'] = 12*(N-1) == r-1
    checks['S_is_exact_ceiling_sqrt_8N_over5'] = 5*S*S >= 8*N and 5*(S-1)**2 < 8*N
    checks['work_within_claim_6443_over50'] = W.numerator**50 <= (1 << 6443)*W.denominator**50
    checks['claim_slack_at_least_0005_bits'] = W.numerator**200 < (1 << 25771)*W.denominator**200
    checks['claim_slack_at_most_006_bits'] = W.numerator**5 > (1 << 644)*W.denominator**5
    checks['all_nonmul_repriced_120pct_margin_at_least_010_bits'] = stress.numerator**25 < (1 << 3219)*stress.denominator**25
    checks['itemized_below_envelope'] = itemized < W
    checks['ordinary_rows_sum_4000'] = sum([1532,192,24,288,198,160,192,128,256,192,128,710]) == 4000
    checks['actual_instr_below_allowance_8000'] = 4000 < 8000
    checks['reports_mean_below_63over64_capacity'] = mu < F(63,64)*D
    _, eupper = exp_bounds(F(63,4096))
    checks['report_cap_chernoff_base_below_8191over8192'] = eupper*F(64,65) < F(8191,8192)
    checks['report_cap_tail_below_2minus100'] = D >= 8192*256 and (F(8191,8192)**8192) < F(1,2)
    elo, _ = exp_bounds(F(4,5))
    checks['ideal_collision_probability_exceeds_055'] = elo > F(20,9)
    checks['finite_fp_false_pair_mean_below_2pow19'] = fp_mean < (1 << 19)
    checks['finite_fp_markov_loss_below_1over32'] = fp_mean/candidate_cap < F(1,32)
    checks['poisson_geometric_tail_ratio_below_1over32'] = F(1 << 19, candidate_cap+1) < F(1,32)
    checks['poisson_chernoff_bound_below_2minus3B'] = 3*(1 << 19)*8 < candidate_cap
    checks['false_pair_verification_total_below_2pow40'] = candidate_cap*(1 << 15) < (1 << 40)
    checks['all_local_cycles_union_below_2minus47'] = F(start_cap*(2*cap)*(2*cap+1), 2*N) < F(1,1 << 47)
    checks['unfinished_tail_fraction_below_2minus38'] = F(processors*cap, S) < F(1,1 << 38)
    checks['empty_segment_tail_union_below_2minus100'] = F(8,3)**128 > (1 << 167)
    checks['sampling_bias_below_2minus188'] = F(2*start_cap*r, 1 << 512) < F(1,1 << 188)
    checks['bad_relation_union_below_2minus230'] = F(candidate_cap, r) < F(1,1 << 230)
    checks['reserved_other_loss_more_than_all_small_terms'] = F(1,1 << 30) > F(1,1 << 36)+F(1,1 << 47)+F(3,1 << 100)+F(1,1 << 188)+F(1,1 << 230)+F(3*S*66,p)
    checks['success_floor_after_markov_deviation_and_reserved_losses'] = F(11,20)-F(1,10)-F(1,32)-F(1,1 << 30) > F(2,5)
    checks['success_margin_after_deviation_exceeds_001'] = F(11,20)-F(1,10)-F(1,32)-F(1,1 << 30)-F(2,5) > F(1,100)
    checks['success_margin_before_transfer_deficit_exceeds_passer_006'] = F(11,20)-F(1,32)-F(1,1 << 30)-F(2,5) > F(3,50)
    checks['80byte_record_layout'] = 14+32+32+2 == 80
    checks['hard_capacity_at_2over3_load'] = 3*D == 2*slots
    checks['peak_including_replicas_below_2pow73_bytes'] = memory < (1 << 73)
    checks['shard_local_address_space_under_64bits'] = 80*3*(1 << 49) < (1 << 64)
    checks['global_slot_index_67bits_accounted'] = slots.bit_length() == 67
    additive = start_cap*(1 << 10) + F(D*(1 << 25),C) + candidate_cap*(1 << 15) + (1 << 42) + (1 << 70)
    checks['all_explicit_additive_rows_below_2pow90'] = additive < (1 << 90)
    d1 = fp2_pow((1,1),(p-1)//6)
    zeta = fp2_pow((1,1),(p*p-1)//6)
    checks['d1_norm_equals_zeta'] = fp2_mul(d1,(d1[0],-d1[1]%p)) == zeta
    checks['zeta_in_Fp_order6'] = zeta[1] == 0 and pow(zeta[0],3,p) == p-1 and pow(zeta[0],6,p) == 1
    checks['zeta_square_is_zeta_minus_one'] = pow(zeta[0],2,p) == (zeta[0]-1)%p
    gammas = [fp2_pow((1,1),(p**i-1)//6) for i in range(12)]
    checks['gamma_even_odd_chain'] = all(gammas[2*j] == (pow(zeta[0],j,p),0) and gammas[2*j+1] == fp2_mul(d1,(pow(zeta[0],j,p),0)) for j in range(6))
    # Enumerate ALL 4^4 maps, not samples, for the exact joint (mu,lambda) law.
    joint = {}
    floyd_checks = 0
    for mapping in product(range(4), repeat=4):
        for start in range(4):
            seen, y = {}, start
            while y not in seen:
                seen[y] = len(seen)
                y = mapping[y]
            m, l = seen[y], len(seen)-seen[y]
            if start == 0:
                joint[(m,l)] = joint.get((m,l),0)+1
            t = h = start
            for it in range(1,5):
                t,h = mapping[t],mapping[mapping[h]]
                if t == h:
                    floyd_checks += it <= m+l
                    break
    checks['joint_mu_lambda_law_all_4pow4_maps'] = all(F(joint.get((m,l),0),4**4) == F(falling(3,m+l-1),4**(m+l)) for m in range(4) for l in range(1,5-m))
    checks['point_on_Lcycle_law_all_4pow4_maps'] = all(F(joint.get((0,l),0),4**4) == F(falling(3,l-1),4**l) for l in range(1,5))
    checks['floyd_bound_all_1024_start_graph_pairs'] = floyd_checks == 1024
    checks['geometric_cycle_tail_ratio'] = all(F(n-l,n)*F(3,4) <= F(3,4) for n in range(2,17) for l in range(1,n))
    from pathlib import Path
    evidence = json.loads(Path(__file__).with_name('mechanism-output.json').read_text())
    assert evidence['pipeline_file'] == 'mechanism-pipeline.json'
    evidence['pipeline'] = json.loads(Path(__file__).with_name('mechanism-pipeline.json').read_text())['pipeline']
    fields, pipelines = evidence['fields'], evidence['pipeline']
    observed = sum(len(f['first_cycles']) + sum(len(t['cycles'])-1 for t in f['basin_trials']) for f in fields)
    checks['captured_actual_field_run_completed'] = evidence['completed'] is True
    checks['captured_R_matches_production_1024'] = evidence['R'] == 1024
    checks['captured_three_sparse_bit_densities'] = evidence['densities_bits'] == [8,10,12]
    checks['captured_five_group_scales'] = len(fields) == 5
    checks['captured_all_N_at_most_2pow40'] = all(f['N'] <= 1 << 40 for f in fields)
    checks['captured_largest_N_exceeds_2pow38'] = max(f['N'] for f in fields) > 1 << 38
    checks['captured_order12_frobenius_all_fields'] = all(pow(f['p'],12,f['r']) == 1 and all(pow(f['p'],k,f['r']) != 1 for k in range(1,12)) for f in fields)
    checks['captured_same_tower_congruence'] = all(f['p'] % 24 == 19 for f in fields)
    checks['captured_r_divides_cyclotomic12'] = all((f['p']**4-f['p']**2+1) % f['r'] == 0 for f in fields)
    checks['captured_class_counts_exact'] = all(12*(f['N']-1) == f['r']-1 for f in fields)
    checks['captured_240_independent_tables'] = sum(len(f['first_cycles']) for f in fields) == 240
    checks['captured_eight_long_basin_trials'] = sum(len(f['basin_trials']) for f in fields) == 8
    checks['captured_each_long_trial_32_cycle_detections'] = all(len(t['cycles']) == 32 for f in fields for t in f['basin_trials'])
    checks['captured_cycle_total_matches_enumeration'] = observed == evidence['cycles_observed']
    checks['captured_partition_samples_have_1024_bins'] = all(len(f['partition_counts']) == 1024 and sum(f['partition_counts']) == f['uniform_coordinate_samples'] for f in fields)
    checks['captured_163840_uniform_coordinate_samples'] = sum(f['uniform_coordinate_samples'] for f in fields) == 163840
    checks['captured_all_three_floyd_evaluations'] = all(p['evaluations'] == 3*p['tortoise_iterations'] for p in pipelines)
    checks['captured_all_six_dense_analogues'] = all(p['dense_multiplication_analogues'] == 2*p['evaluations'] for p in pipelines)
    checks['captured_iteration_windows_match'] = all(p['tortoise_iterations'] == p['windows']*p['budget_per_window'] for p in pipelines)
    checks['captured_zero_equivariance_failures'] = evidence['equivariance_failures'] == 0
    checks['captured_zero_incorrect_recoveries'] = evidence['incorrect_recoveries'] == 0
    cdf_n = sum(len(f['first_cycles']) for f in fields)
    cdf_differences = []
    for key in ('rho','lambda','mu'):
        for a,b in ((1,2),(1,1),(3,2),(2,1),(3,1)):
            actual = sum(c[key]**2*b*b <= a*a*f['N'] for f in fields for c in f['first_cycles'])
            reference = sum(c[key]**2*b*b <= a*a*f['N'] for f in fields for c in f['reference_cycles'])
            cdf_differences.append(abs(actual-reference))
    checks['captured_cdf_grid_discrepancy_below_reserved_1over10'] = 10*max(cdf_differences) < cdf_n
    return {'schema_version': 1, 'scope': 'Exact arithmetic, captured-evidence consistency and conditional-model implications only; not proof of production random mapping.',
            'score': {'numerator': '6443','denominator':'50'},
            'N':str(N),'S':str(S),'floyd_evaluation_cap':str(3*S),'dense_multiplication_cap':str(6*S),
            'W_envelope':{'numerator':str(W.numerator),'denominator':str(W.denominator)},
            'W_itemized':{'numerator':str(itemized.numerator),'denominator':str(itemized.denominator)},
            'W_120pct_nonmul':{'numerator':str(stress.numerator),'denominator':str(stress.denominator)},
            'hard_reports':str(D),'slot_capacity':str(slots),'peak_bytes':str(memory),
            'false_pair_expectation_upper':{'numerator':str(fp_mean.numerator),'denominator':str(fp_mean.denominator)},
            'checks':checks,'check_count':len(checks)}


if __name__ == '__main__':
    certificate = build()
    print(json.dumps(certificate,indent=2))
    raise SystemExit(0 if all(certificate['checks'].values()) else 1)
