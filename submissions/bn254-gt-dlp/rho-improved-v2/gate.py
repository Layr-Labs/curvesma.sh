"""Gate and expected-work ledger of this submission (Python standard library only).

Reads {"mean_allowance_percent": int, "cap_multiple": int, "seed": int} on stdin (defaults 5, 64) and recomputes,
in exact decimal arithmetic, every figure of proof.md sections 3.4, 4 and 5 that differs from the base
package: the expected number of steps, the step cap and its Markov bound, the failure budget, the gate
constant K, the expected total work and the declared score. The step prices are the counted prices of
machine.md (reproduced by experiment.py part D). It also implements the attack-side counters of machine.md
section 4 (three-word work meter and step counter, two-word stuck and report counters) on 64-bit words with
add-with-carry and word-by-word comparison, and checks them against big integers up to their caps.
Nothing here is a production run.
"""
import json
import sys
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
from fractions import Fraction
from pathlib import Path

getcontext().prec = 80
base = Path("/snapshot/files")
if not base.is_dir():
    base = Path(__file__).parent
try:
    request = json.load(sys.stdin)
except ValueError:
    request = {}
allow = Fraction(100 + int(request.get("mean_allowance_percent", 5)), 100)
cap_multiple = int(request.get("cap_multiple", 64))
profile = json.loads((base / "target-profile.json").read_text())
r = int(profile["subgroup_order"], 16)
p = int(profile["field_modulus"], 16)

LN2 = Decimal(2).ln()


def lg(v):
    return Decimal(v).ln() / LN2


def dec(fr):
    return Decimal(fr.numerator) / Decimal(fr.denominator)


# pi to 70 digits by Machin's formula
def arctan_inv(n):
    x = Decimal(1) / n
    total, term, k, n2 = x, x, 1, n * n
    while abs(term) > Decimal(10) ** -75:
        term = -term / n2
        k += 2
        total += term / k
    return total


PI = 4 * (4 * arctan_inv(5) - arctan_inv(239))

# ---- quantities shared with the base package (proof.md sections 2-5) ----
assert (r - 1) % 12 == 0 and pow(p, 12, r) == 1 and pow(p, 6, r) == r - 1
N = (r - 1) // 12                                   # Frobenius classes
R_TABLE = 1024
S_SQ = Fraction(1, 2)                               # squaring probability
sum_q2 = S_SQ**2 + (1 - S_SQ) ** 2 / R_TABLE        # conservative pair-collision loss
N_eff = Decimal(N) / dec(1 - sum_q2)                # N' = N / (1 - sum q^2)
MUL_STEP, SQ_STEP, C_UNIT = 18191, 6965, 20000      # counted instructions per step (machine.md)
c_avg = dec(Fraction(MUL_STEP + SQ_STEP, 2 * C_UNIT))
c_mul = dec(Fraction(MUL_STEP, C_UNIT))
eps = dec(Fraction(5, 12 * R_TABLE) + Fraction(250, R_TABLE**2))
THETA_BITS, WALKS = 60, 2**20
E_delay = Decimal(WALKS * 20) * Decimal(2) ** THETA_BITS

# ---- this submission ----
sqrtN = N_eff.sqrt()
model_mean = (PI * N_eff / 2).sqrt()                # H1: E[T_c] = sqrt(pi N'/2)
ET = dec(allow) * model_mean                        # charged mean (H-mean)
S_cap = cap_multiple * sqrtN                        # step cap
markov = ET / S_cap                                 # P(T_c > S_cap) <= E[T_c]/S_cap
loss_reports_cap = Decimal(2) ** -8                 # report count above 2^8 times its mean (Markov)
loss_small = Decimal("0.001")                       # detection e^-20, single-walk collisions 2^-20, key ties,
#                                                     false-candidate and restart caps, beta1 = beta2
fail = markov + loss_reports_cap + loss_small
success_claimed = Decimal("0.97")
K = int((Decimal("0.39") / success_claimed * Decimal(2) ** 64).to_integral_value(rounding=ROUND_CEILING))
gate = Decimal(K) / Decimal(2) ** 64
walk = c_avg * ET * (1 + eps)                       # expected walk work up to the first useful collision
delay = c_mul * E_delay                             # detection delay, charged at the multiply price
NONWALK = Decimal(2) ** 66                          # reports, setup, restarts, candidates, output (expected)
reports_expected = ET * Decimal(2) ** -THETA_BITS * (1 + eps) + Decimal(WALKS)
run = walk + delay + NONWALK
total = 1 + gate * run
abort = c_mul * (S_cap * (1 + eps) + E_delay) + NONWALK * 2**8   # every-run meter (not the score)

# ---- attack-side counters (machine.md section 4): exact integer model on 64-bit words ----
MASK = 2**64 - 1


class WordCounter:
    """An n-word little-endian counter updated only with 64-bit adds-with-carry and compared word by word."""

    def __init__(self, words):
        self.w = [0] * words
        self.overflow = False

    def add(self, value):                       # value < 2^64: adds, then adcs/adc with zero
        carry = value
        for i in range(len(self.w)):
            s = self.w[i] + carry
            self.w[i] = s & MASK
            carry = s >> 64
        if carry:
            self.overflow = True

    def value(self):
        return sum(x << (64 * i) for i, x in enumerate(self.w))

    def at_least(self, cap_words):              # cmp / sbcs chain from the low word, then b.hs
        borrow = 0
        for a, b in zip(self.w, cap_words):
            borrow = 1 if a - b - borrow < 0 else 0
        return borrow == 0


def words_of(n, k):
    return [(n >> (64 * i)) & MASK for i in range(k)]


meter_cap_instr = int((abort * C_UNIT).to_integral_value(rounding=ROUND_CEILING))   # cap in instructions
step_cap_int = int(S_cap.to_integral_value(rounding=ROUND_FLOOR))
stuck_threshold = 20 * 2**THETA_BITS
report_cap = int((reports_expected * 2**8).to_integral_value(rounding=ROUND_CEILING))
MAX_SINGLE_CHARGE = 2**48                       # largest amount added to the meter in one addition
import random as _random
_rng = _random.Random(int(request.get("seed", 20261111)))
counter_tests = 0
counters_ok = True
for words, cap in ((3, meter_cap_instr), (3, step_cap_int), (2, stuck_threshold), (2, report_cap)):
    capw = words_of(cap, words)
    counters_ok = counters_ok and cap < 2 ** (64 * words) and cap + MAX_SINGLE_CHARGE < 2 ** (64 * words)
    for start in (0, 1, cap - MAX_SINGLE_CHARGE, cap - 18191, cap - 1, cap, 2**64 - 1, 2**128 - 5,
                  _rng.randrange(cap), _rng.randrange(cap), _rng.randrange(cap)):
        start = max(0, min(start, cap))
        ctr = WordCounter(words)
        ctr.w = words_of(start, words)
        ref = start
        for _ in range(200):
            inc = _rng.choice((1, 6965, 18191, 10000, MAX_SINGLE_CHARGE - 1, _rng.randrange(1, MAX_SINGLE_CHARGE)))
            if ctr.at_least(capw) != (ref >= cap):
                counters_ok = False
            if ref >= cap:                      # the run aborts here; no further addition is made
                break
            ctr.add(inc)
            ref += inc
            counter_tests += 1
            if ctr.value() != ref or ctr.overflow:
                counters_ok = False
BOOKKEEPING = {"T bits and DP test": 12, "table address": 4, "exponent update": 52, "2-cycle test": 10,
               "checkpoint test": 10, "history": 12, "checkpoint period (per-walk counter mod 64)": 6,
               "stuck-walk counter, two words": 7, "global step counter and step cap, three words": 11,
               "work meter, three words": 12, "branch dispatch": 6, "walk loop": 8}
itemised = sum(BOOKKEEPING.values())
out_counters = {
    "meter_cap_instructions_bits": meter_cap_instr.bit_length(), "log2_meter_cap_instructions": float(lg(meter_cap_instr)),
    "step_cap_bits": step_cap_int.bit_length(), "stuck_threshold_bits": stuck_threshold.bit_length(),
    "report_cap_bits": report_cap.bit_length(),
    "meter_fits_192_bits_with_headroom": meter_cap_instr + MAX_SINGLE_CHARGE < 2**192,
    "meter_would_wrap_128_bits": meter_cap_instr >= 2**128,
    "stuck_threshold_needs_two_words": stuck_threshold >= 2**64,
    "word_counter_additions_checked": counter_tests, "word_counters_match_big_integers": counters_ok,
    "bookkeeping_itemised_instructions": itemised, "bookkeeping_charged_instructions": 300,
    "bookkeeping_within_charge": itemised <= 300,
}

try:
    sc = json.loads((base / "submission.json").read_text())["claims"][0]["score"]
    declared = Decimal(sc["numerator"]) / Decimal(sc["denominator"])
except (OSError, KeyError, ValueError):
    declared = Decimal("123.63")
out = {
    "inputs": {"mean_allowance": str(dec(allow)), "cap_multiple": cap_multiple},
    "log2_N_classes": float(lg(N)), "log2_N_effective": float(lg(N_eff)),
    "average_step_units": str(c_avg), "multiply_step_units": str(c_mul), "eps": float(eps),
    "log2_model_mean_steps": float(lg(model_mean)), "log2_charged_mean_steps": float(lg(ET)),
    "log2_step_cap": float(lg(S_cap)), "markov_bound_cap": float(markov),
    "failure_budget_past_gate": float(fail), "success_past_gate_lower": float(1 - fail),
    "success_past_gate_claimed": str(success_claimed), "budget_ok": 1 - fail >= success_claimed,
    "gate_K": str(K), "gate_probability": float(gate),
    "total_success": float(gate * success_claimed), "success_ok": gate * success_claimed >= Decimal("0.39"),
    "K_is_least": Decimal(K - 1) / Decimal(2) ** 64 * success_claimed < Decimal("0.39"),
    "log2_expected_walk_units": float(lg(walk)), "log2_delay_units": float(lg(delay)),
    "log2_expected_reports": float(lg(reports_expected)),
    "nonwalk_bucket_covers_reports": reports_expected / 2 + Decimal(2) ** 60 <= NONWALK,
    "log2_expected_run_units": float(lg(run)), "log2_expected_total_units": float(lg(total)),
    "score_rounded_up_to_hundredth": str((lg(total) * 100).to_integral_value(rounding=ROUND_CEILING) / 100),
    "declared_score": str(declared), "declared_ok": declared >= lg(total),
    "log2_every_run_meter_units": float(lg(abort)),
    "base_claim_for_comparison": {"log2_S_at_model_0.43": float(lg((2 * N_eff * (Decimal(1) / Decimal("0.57")).ln()).sqrt())),
                                  "bits_gained_on_walk": float(lg((2 * N_eff * (Decimal(1) / Decimal("0.57")).ln()).sqrt() * c_avg) - lg(gate * walk))},
    "counters": out_counters,
    "scope": "Arithmetic of the gate and expected-work ledger only. H1, H2 and H-mean are heuristics; no production run.",
}
print(json.dumps(out))
