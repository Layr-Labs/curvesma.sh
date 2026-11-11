"""Run-gate ledger for the tag-tracing rho, plus a toy mean-collision measurement.

Part L  exact ledger: step cap, gate constant K (least integer), success past the gate,
        expected total work and every-run maximum, in exact decimal arithmetic. Prices
        (273 instructions per certified step, 24,576 per window end, 2^110 fixed bucket)
        are those of machine.md, from the supplied listings; nothing is re-priced here.
Part T  toy measurement of the mean total collision time of the serial 8-adding
        distinguished-point walk (selector floor(8y/q), as in production but on a
        one-coordinate field) in prime-order subgroups of small prime fields: mean of
        T_c / sqrt(n) against sqrt(pi*kappa/2) for kappa = 8/7 and the charged 1.25.
Standard library only. Optional JSON on stdin: {"seed": int, "trials": int, "bits": int}.
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from pathlib import Path

try:
    cfg = json.loads(sys.stdin.read() or "{}")
except ValueError:
    cfg = {}
SEED = int(cfg.get("seed", 1))
TRIALS = int(cfg.get("trials", 200))
BITS = int(cfg.get("bits", 20))

base = Path("/snapshot/files")
if not (base / "target-profile.json").is_file():
    base = Path(__file__).resolve().parent
profile = json.loads((base / "target-profile.json").read_text())
r = int(profile["subgroup_order"], 16)

C_UNIT = 20000
STEP_I = 273                 # certified step (ttstep.disasm.txt)
WINDOW_I = 24576             # window end, per 256 certified steps (tower.disasm.txt + driver)
WINDOW = 256
FIXED_LOG2_I = 110           # every other work, metered (machine.md "Fixed work")
GATE_I = 64                  # draw one 64-bit word, compare, branch (generous)
KAPPA = Decimal("1.25")      # charged collision constant (model 8/7)
CAP_MULT = 64                # step cap = 64 sqrt(kappa r)
DP_BITS = 86                 # distinguished density 2^-86
WALK_LEN_LOG2 = 92           # per-walk step cap
WALKS_LOG2 = 50              # walk cap (store sized for 2^50 keys)
EVAL_CAP_LOG2 = 90           # exact-evaluation meter
UNCERT_LOG2 = Decimal(97).ln() / Decimal(2).ln() - 61   # P[uncertified step] under H2


def part_L():
    getcontext().prec = 80
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2
    pi = Decimal("3.14159265358979323846264338327950288419716939937510")
    kr = KAPPA * r
    mean_tc = (pi * kr / 2).sqrt()                     # E[T_c] <= int_0^inf exp(-s^2/(2 kappa r)) ds
    mean_tail = Decimal(2) ** DP_BITS                  # E[steps from the touch to the stored DP], H2
    S_cap = int((CAP_MULT * kr.sqrt()).to_integral_value(rounding="ROUND_CEILING"))
    per_step = Decimal(STEP_I) + Decimal(WINDOW_I) / WINDOW   # 369
    run_max_I = GATE_I + STEP_I * S_cap + WINDOW_I * (S_cap // WINDOW + 1) + 2 ** FIXED_LOG2_I
    p_self = Decimal(2) ** (WALKS_LOG2 + 2 * DP_BITS + 1) / r
    p_no_dp = Decimal(2) ** WALKS_LOG2 * Decimal(-64).exp()
    exp_steps = mean_tc + 1 + mean_tail + (p_self + p_no_dp) * S_cap
    exp_run_I = per_step * exp_steps + WINDOW_I + Decimal(2) ** FIXED_LOG2_I
    # losses past the gate
    exp_walks = exp_steps / Decimal(2) ** DP_BITS + 1
    exp_evals = exp_steps * Decimal(2) ** UNCERT_LOG2 + exp_walks
    losses = {
        "step_cap_markov": exp_steps / S_cap,
        "walk_cap_markov": exp_walks / Decimal(2) ** WALKS_LOG2,
        "exact_eval_meter_markov": exp_evals / Decimal(2) ** EVAL_CAP_LOG2,
        "store_probe_meter_markov_inherited": Decimal(2) ** -16,
        "walk_without_dp_within_2^92_H2": Decimal(2) ** WALKS_LOG2 * Decimal(-64).exp(),
        "self_cycle_before_dp": Decimal(2) ** WALKS_LOG2 * Decimal(2) ** (2 * DP_BITS + 1) / r,
        "degeneracy_and_sampling": Decimal(2) ** -200,
    }
    past = 1 - sum(losses.values())
    sigma = Decimal("0.97")
    assert past >= sigma
    K = int((Decimal("0.39") / sigma * Decimal(2) ** 64).to_integral_value(rounding="ROUND_CEILING"))
    gamma = Decimal(K) / Decimal(2) ** 64
    exp_total_I = GATE_I + gamma * exp_run_I
    claim_log2 = lg(exp_total_I / C_UNIT)
    claim = (claim_log2 * 10 ** 6).to_integral_value(rounding="ROUND_CEILING") / 10 ** 6
    checks = {
        "K_least": gamma * sigma >= Decimal("0.39") and Decimal(K - 1) / Decimal(2) ** 64 * sigma < Decimal("0.39"),
        "past_gate_ge_sigma": past >= sigma,
        "S_cap_bits_le_192": S_cap.bit_length() <= 192,
        "walk_cap_times_len_covers_S_cap": WALKS_LOG2 + WALK_LEN_LOG2 >= S_cap.bit_length(),
        "claim_ge_computed": claim >= claim_log2,
    }
    return {
        "kappa_charged": str(KAPPA), "log2_mean_Tc": float(lg(mean_tc)),
        "mean_Tc_over_sqrt_r": float(mean_tc / Decimal(r).sqrt()),
        "S_cap": str(S_cap), "log2_S_cap": float(lg(S_cap)),
        "per_step_instructions": float(per_step),
        "log2_expected_run_units": float(lg(exp_run_I / C_UNIT)),
        "log2_every_run_max_units": float(lg(Decimal(run_max_I) / C_UNIT)),
        "losses": {k: float(v) for k, v in losses.items()},
        "success_past_gate_lb": float(past), "sigma_claimed": str(sigma),
        "K": str(K), "gate_probability": float(gamma), "success_lb": float(gamma * sigma),
        "log2_expected_total_units": float(claim_log2), "claim_time_log2": str(claim),
        "checks": checks,
    }


# ------------------------------------------------------------------ Part T: toy collision law
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def toy_field(rng, bits):
    while True:
        n = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
        if not is_prime(n):
            continue
        for k in range(2, 4000, 2):
            q = k * n + 1
            if is_prime(q):
                while True:
                    h = pow(rng.randrange(2, q - 1), k, q)
                    if h != 1:
                        return q, n, h


def toy_trial(rng, bits, dp_bits):
    q, n, h = toy_field(rng, bits)
    x = rng.randrange(n)
    z = pow(h, x, q)
    M = []
    for _ in range(8):
        c, d = rng.randrange(n), rng.randrange(n)
        M.append(pow(h, c, q) * pow(z, d, q) % q)
    seen = set()
    total = 0
    dp_lim = q >> dp_bits
    while True:
        a, b = rng.randrange(n), rng.randrange(n)
        y = pow(h, a, q) * pow(z, b, q) % q
        while True:
            if y in seen:
                return total / math.sqrt(n)
            seen.add(y)
            if y < dp_lim:
                break
            y = y * M[(8 * y) // q] % q
            total += 1


def part_T():
    rng = random.Random(SEED)
    vals = [toy_trial(rng, BITS, 4) for _ in range(TRIALS)]
    m = sum(vals) / len(vals)
    sd = math.sqrt(sum((v - m) ** 2 for v in vals) / (len(vals) - 1))
    return {"bits": BITS, "trials": TRIALS, "dp_bits": 4,
            "mean_Tc_over_sqrt_n": m, "standard_error": sd / math.sqrt(len(vals)),
            "model_kappa_8_7": math.sqrt(math.pi * 8 / 7 / 2),
            "charged_kappa_1_25": math.sqrt(math.pi * 1.25 / 2),
            "within_charged": m <= math.sqrt(math.pi * 1.25 / 2)}


L = part_L()
T = part_T()
out = {"seed": SEED, "partL_gate_ledger": L, "partT_toy_mean_collision": T,
       "all_ok": all(L["checks"].values())}
print(json.dumps(out, indent=1, default=str))
sys.exit(0 if out["all_ok"] else 1)
