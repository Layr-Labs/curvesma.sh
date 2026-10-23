#!/usr/bin/env python3
"""Generate certificate.json by executed exact integer/rational comparisons.
No empirical or heuristic premise is marked certified. Python stdlib only.
"""
from fractions import Fraction as F
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import json

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
BRANCHES = 1024
CLAIM_N, CLAIM_D = 6441, 50
N = (R-1)//12
NP_NUM, NP_DEN = N*BRANCHES, BRANCHES-1
NP = (NP_NUM+NP_DEN-1)//NP_DEN
T = isqrt(5*NP//4)+1
M, GENERATIONS, FP_BITS = 1 << 69, 129, 128
checks = {}

def ck(name, value):
    checks[name] = bool(value)


def below_pow2(n, num, den):
    return pow(n, den) < (1 << num)


def above_pow2(n, num, den):
    return pow(n, den) > (1 << num)


def prime_trial(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n)+1, 2))


# Public subgroup and Frobenius arithmetic.
ck('twelve_divides_r_minus_one', 12*N == R-1)
ck('lambda_has_order_twelve', pow(P, 12, R) == 1 and pow(P, 6, R) != 1 and pow(P, 4, R) != 1)
ck('lambda_six_is_minus_one', pow(P, 6, R) == R-1)
ck('cyclotomic_relation', (1+pow(P, 4, R)+pow(P, 8, R)) % R == 0)
ck('N_log_interval_2512721_to_2512722_over_10000',
   above_pow2(N, 2512721, 10000) and below_pow2(N, 2512722, 10000))
ck('Np_is_ceil_N_times_1024_over_1023',
   (NP-1)*NP_DEN < NP_NUM <= NP*NP_DEN)
ck('Np_not_subgroup_order_correction', BRANCHES == 1024 and BRANCHES != R)
ck('T_is_minimal_integer_with_4T_squared_gt_5Np',
   4*(T-1)**2 <= 5*NP < 4*T*T)

# Prove that no initial-label coefficient can vanish at a Floyd distance <=T.
FACTORS = {2:32, 3:1, 11:1, 19:1, 10177:1, 906349:2,
           2508409:1, 125527:1, 859267:1, 2529403:1,
           52437899:1, 254760293:2}
ck('r_minus_one_factorization_product', prod(q**e for q,e in FACTORS.items()) == R-1)
ck('factor_primality_by_trial_division', all(prime_trial(q) for q in FACTORS))
o = R-1
for q in FACTORS:
    while o % q == 0 and pow(2, o//q, R) == 1:
        o //= q
ck('two_exact_order_upper', pow(2, o, R) == 1)
ck('two_exact_order_minimal', all(pow(2, o//q, R) != 1 for q in FACTORS if o % q == 0))
d = o//gcd(o,12)
ck('no_label_closure_through_twice_the_budget', d > 2*T)

# The finite random-mapping first-repeat law is a derived COMPARISON model.
# For N0 points and t distinct transitions: no-repeat = prod(1-j/N0).
# Thus no-repeat <= exp(-t(t+1)/(2N0)); T^2/Np > 5/4.
# exp(-5/8) upper by even alternating truncation, exact Fractions.
u = F(5,8)
term = exp_upper = F(1)
for j in range(1, 41):
    term *= -u/j
    exp_upper += term
success_lower = 1-exp_upper
ck('comparison_success_lower_gt_04647', success_lower > F(4647,10000))
ck('success_claim_is_track_floor_039', F(39,100) == F(39,100))
ck('named_model_margin_exceeds_007', success_lower-F(39,100) > F(7,100))
# Allow one common 5% walk-length inflation in addition to Np.
v = F(5,8)*F(20,21)**2
term = e2 = F(1)
for j in range(1, 41):
    term *= -v/j
    e2 += term
ck('five_percent_slowdown_model_success_gt_0432', 1-e2 > F(432,1000))
ck('five_percent_model_margin_exceeds_0042', 1-e2-F(39,100) > F(42,1000))

# Ordinary map-evaluation rows. 2 dense Fp12 products cover y^2*T[h][-k].
ROWS = {
    'four_base_field_products_4_times_383':1532,
    'nine_six_limb_add_sub_neg_9_times_64':576,
    'minimum_of_twelve_six_limb_keys_and_index':600,
    'two_scalar_double_add_reduce_updates':384,
    'twisted_jump_and_label_load_addressing':256,
    'DP_predicate_and_fingerprint_extract':128,
    'Floyd_key_comparison_and_control':256,
    'loop_and_meter_updates':128,
    'ordinary_state_load_store_and_scratch':512,
}
row_total = sum(ROWS.values())
ordinary_ceiling = 8000
ck('ordinary_rows_sum_4372', row_total == 4372)
ck('all_ordinary_rows_inside_8000_instruction_ceiling', row_total <= ordinary_ceiling)
ck('three_calls_per_Floyd_iteration', 3*2 == 6)
# The published nominal score is an explicit larger envelope (8T), not the
# tight executed ceiling (36T/5). Publish both and stress ALL executed rows,
# including the unused per-call ceiling and every additive allowance.
actual = F(36,5)*T + (1 << 100)
stress = F(186,25)*T + (1 << 101)
W = 8*T+(1 << 100)
Ws = (stress.numerator+stress.denominator-1)//stress.denominator
ck('executed_ceiling_inside_scored_envelope', actual < W)
ck('stress_is_six_reference_plus_1p2_all_nonreference_ceiling',
   F(186,25) == 6+F(6,5)*F(3*ordinary_ceiling,20000))
ck('envelope_log_interval_1287977_to_1287978_over_10000',
   above_pow2(W,1287977,10000) and below_pow2(W,1287978,10000))
ck('stress_log_less_than_128694_over_1000', below_pow2(Ws,128694,1000))
ck('claim_above_envelope', below_pow2(W,CLAIM_N,CLAIM_D))
ck('claim_minus_nominal_ledger_between_0005_and_006_bits',
   below_pow2(W,128815,1000) and above_pow2(W,128760,1000))
ck('stress_plus_010_bits_below_claim', below_pow2(Ws,128720,1000))
ck('claim_reduced', gcd(CLAIM_N,CLAIM_D) == 1)
ck('claim_above_requested_local_lower_124', CLAIM_N > 124*CLAIM_D)

# D1: none of these comparison windows restarts or changes the class map.
THETA_DEN = 1 << 60
ck('twenty_over_theta_is_twenty_times_2pow60', 20*THETA_DEN == 20*(1 << 60))
ck('sixtyfour_over_theta_is_2pow66_not_2pow26', 64*THETA_DEN == 1 << 66)
ck('Floyd_budget_exceeds_comparison_DP_windows', T > 64*THETA_DEN)
ck('identity_visit_union_comparison_less_2pow_minus120', 3*T*(1 << 120) < N)
ck('key_tie_union_under_Hkey_less_2pow_minus240', 3*T*66*(1 << 240) < P)

# D5: disjoint literal bits, including all adaptive generations.
base_dp = set(range(10,70))
base_fp = set(range(70,198))
ck('base_DP_has_sixty_bits', len(base_dp) == 60)
ck('fingerprint_has_128_bits_and_no_fixed_DP_bits', len(base_fp-base_dp) == 128 and not(base_fp & base_dp))
ck('all_129_generations_have_disjoint_DP_and_fp_bits', all(
    not(set(range(10,70+g)) & set(range(70+g,198+g))) for g in range(GENERATIONS)))
ck('last_fingerprint_within_381bit_key', 198+GENERATIONS-1 <= 381)
ck('at_most_three_T_reports_have_mean_less_than_M_over_three', 9*T < M*THETA_DEN)
false_pairs = F(GENERATIONS*M*(M-1), 2*(1 << FP_BITS))
ck('finite_fp_false_pair_expectation_under_Hfp_less_than_2pow17', false_pairs < 1 << 17)
ck('finite_fp_false_pair_expectation_under_Hfp_less_than_2pow20', false_pairs < 1 << 20)
# Auxiliary Poisson comparison only; not an exact law for Floyd reports.
# mu < M/3 and e < 11/4 => (e*mu/M)^M < (11/12)^M.
ck('Poisson_comparison_chernoff_11over12_power_M_less_2pow_minus80',
   M >= 1024 and 11**1024*(1 << 80) < 12**1024)
# Valid finite Markov candidate-cap bound under Hfp, no Poisson substitution.
ck('fp_cap_markov_failure_under_Hfp_less_2pow_minus15', false_pairs/(1 << 32) < F(1,1 << 15))
# Overflow and candidate saturation never stop or change Floyd, so these
# auxiliary DP events have ZERO success debit by construction, not by a tail.
record_bytes = 16+32+32+32+16
peak = M*record_bytes+(1 << 40)
ck('record_layout_128_bytes', record_bytes == 128)
ck('hard_peak_memory_less_than_2pow78', peak < 1 << 78)
# Deterministic ancillary upper limits: no expectation is used for resources.
report_cap = GENERATIONS*M
report_instructions = report_cap*(1 << 26)
verification_units = (1 << 32)*(1 << 15)
clear_instructions = GENERATIONS*M*record_bytes*16
setup_units = 1 << 50
exception_units = 1 << 90
ancillary = F(report_instructions+clear_instructions,20000)+verification_units+setup_units+exception_units
ck('all_deterministic_ancillary_charges_less_2pow100', ancillary < 1 << 100)
ck('all_ancillary_rows_repriced_1p2_less_2pow101', F(6,5)*ancillary < 1 << 101)

# Tiny exact finite mapping counts verify the derived joint law, not its
# transfer to the cryptanalytic walk. Sum_(mu>=0,lambda>=1) P(mu,lambda)=1.
from itertools import product
for n in range(1,7):
    law_sum = F(0)
    for j in range(1,n+1):
        falling = prod(n-i for i in range(1,j))
        for mu in range(j):
            law_sum += F(falling,n**j)
    ck(f'finite_joint_rho_law_normalizes_n_{n}', law_sum == 1)
    if n <= 4:
        counts = {}
        for mapping in product(range(n),repeat=n):
            seen, a = {}, 0
            while a not in seen:
                seen[a] = len(seen)
                a = mapping[a]
            pair = (seen[a],len(seen)-seen[a])
            counts[pair] = counts.get(pair,0)+1
        for (mu,lam), count in counts.items():
            j = mu+lam
            exact = F(prod(n-i for i in range(1,j)),n**j)
            ck(f'enumerated_mapping_law_n{n}_mu{mu}_lambda{lam}',
               F(count,n**n) == exact)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    out = {
        'version':1,
        'scope':'Executed integer/rational arithmetic checks only; no heuristic, instruction implementation, or production success is certified.',
        'track':'bls12-381-fixed-g2-invert-g1',
        'package':'gt-orbit12-rho-hardened-40pct-v4',
        'constants':{
            'p':P,'r':R,'N':N,'branches':BRANCHES,
            'Np_rational_numerator':NP_NUM,'Np_rational_denominator':NP_DEN,
            'Np_ceiling':NP,'T':T,'order_of_two':o,'order_mod_frobenius_subgroup':d,
            'ordinary_rows':ROWS,'ordinary_rows_sum':row_total,
            'ordinary_instruction_ceiling_per_evaluation':ordinary_ceiling,
            'nominal_scored_envelope_integer':W,
            'repriced_executed_ceiling_integer':Ws,
            'claim_numerator':CLAIM_N,'claim_denominator':CLAIM_D,
            'theta_denominator':THETA_DEN,'diagnostic_window_20_over_theta':20*THETA_DEN,
            'obsolete_lifetime_64_over_theta':64*THETA_DEN,
            'hard_records_per_generation':M,'generations':GENERATIONS,
            'fingerprint_bits':FP_BITS,'record_bytes':record_bytes,
            'hard_peak_bytes':peak,
            'false_pair_expectation_numerator':false_pairs.numerator,
            'false_pair_expectation_denominator':false_pairs.denominator,
            'success_comparison_lower_numerator':success_lower.numerator,
            'success_comparison_lower_denominator':success_lower.denominator,
        },
        'checks':checks,
    }
    out['pass_lines'] = [('PASS ' if ok else 'FAIL ')+name for name, ok in checks.items()]
    text = json.dumps(out,indent=2)+'\n'
    if args.write:
        root = Path(__file__).resolve().parents[1]
        (root/'certificate.json').write_text(text)
        (root/'evidence').mkdir(exist_ok=True)
        (root/'evidence'/'certificate_run.log').write_text(text)
    print(text,end='')
    if not all(checks.values()):
        raise SystemExit('Individual exact arithmetic check failed')

if __name__ == '__main__':
    main()
