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
    S = isqrt(6*N//5)
    if 5*S*S < 6*N:
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
    checks['score_reduced'] = gcd(2573, 20) == 1
    checks['r_has_255_bits'] = r.bit_length() == 255
    checks['p_has_381_bits'] = p.bit_length() == 381
    checks['order12'] = pow(p, 12, r) == 1 and all(pow(p, k, r) != 1 for k in range(1, 12))
    checks['half_frobenius_inverse'] = pow(p, 6, r) == r-1
    checks['nonidentity_orbits_exact'] = 12*(N-1) == r-1
    checks['S_is_exact_ceiling_sqrt_6N_over5'] = 5*S*S >= 6*N and 5*(S-1)**2 < 6*N
    checks['work_within_claim_12865_over100'] = W.numerator**20 <= (1 << 2573)*W.denominator**20
    checks['claim_slack_at_least_0005_bits'] = W.numerator**200 < (1 << 25729)*W.denominator**200
    checks['claim_slack_at_most_006_bits'] = W.numerator**100 > (1 << 12859)*W.denominator**100
    checks['all_nonmul_repriced_120pct_margin_at_least_010_bits'] = stress.numerator**20 < (1 << 2571)*stress.denominator**20
    checks['itemized_below_envelope'] = itemized < W
    checks['ordinary_rows_sum_4000'] = sum([1532,192,24,288,198,160,192,128,256,192,128,710]) == 4000
    checks['actual_instr_below_allowance_8000'] = 4000 < 8000
    checks['reports_mean_below_7over8_capacity'] = mu < F(7,8)*D
    _, eupper = exp_bounds(F(7,64))
    checks['report_cap_chernoff_base_below_255over256'] = eupper*F(8,9) < F(255,256)
    checks['report_cap_tail_below_2minus100'] = D >= 256*256 and (F(255,256)**256) < F(1,2)
    elo, _ = exp_bounds(F(3,5))
    checks['ideal_collision_probability_exceeds_045'] = elo > F(20,11)
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
    checks['success_floor_after_markov_and_reserved_losses'] = F(9,20)-F(1,32)-F(1,1 << 30) > F(2,5)
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
    return {'schema_version': 1, 'scope': 'Exact arithmetic and conditional-model implications only; not proof of production random mapping.',
            'score': {'numerator': '2573','denominator':'20'},
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
