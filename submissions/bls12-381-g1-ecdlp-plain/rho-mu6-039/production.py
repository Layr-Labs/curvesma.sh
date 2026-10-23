"""Production-constant checks and the exact W bound for proof.md (BLS12-381 G1, plain DLOG).

Checks on the protected constants: n prime, n = 1 mod 3, the GLV relation
phi(P)=(beta*x, y)=[lambda]P on the real generator (and lambda^2), 1+beta+beta^2=0,
then recomputes every ledger term as exact rationals and the rounded-up score.
"""
import json
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

json.load(sys.stdin)
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
n = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
G = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
     0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
C = 2000


def is_prime(m):
    dd, s = m - 1, 0
    while dd % 2 == 0:
        dd //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        x = pow(a, dd, m)
        if x in (1, m - 1):
            continue
        for _ in range(s - 1):
            x = x * x % m
            if x == m - 1:
                break
        else:
            return False
    return True


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0:
            return None
        lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, p) % p
    else:
        lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x3 = (lam * lam - P[0] - Q[0]) % p
    return (x3, (lam * (P[0] - x3) - P[1]) % p)


def smul(k, P):
    out = None
    while k:
        if k & 1:
            out = add(out, P)
        P = add(P, P)
        k >>= 1
    return out


out = {}
assert is_prime(p) and is_prime(n) and n % 3 == 1 and p % 3 == 1
assert (G[1] ** 2 - G[0] ** 3 - 4) % p == 0 and smul(n, G) is None
out["n_prime_and_1_mod_3"] = True
betas = {pow(a, (p - 1) // 3, p) for a in (2, 3, 5)} - {1}
lams = {pow(a, (n - 1) // 3, n) for a in (2, 3, 5, 7)} - {1}
pair = None
for be in betas:
    for la in lams:
        if smul(la, G) == (be * G[0] % p, G[1]):
            pair = (be, la)
beta, lam = pair
assert smul(lam * lam % n, G) == (beta * beta * G[0] % p, G[1])
assert (1 + beta + beta * beta) % p == 0 and (1 + lam + lam * lam) % n == 0
out["beta"] = hex(beta)
out["lambda"] = hex(lam)
out["glv_phi_equals_lambda_on_G"] = True
out["log2_n"] = math.log2(n)

getcontext().prec = 80
SLACK = F(104, 100)
# T = ceil(SLACK * sqrt(2 ln(1/0.61) * n / 6)) additions, computed with 80-digit decimals.
core = (2 * (Decimal(100) / Decimal(61)).ln() * Decimal(n) / 6).sqrt()
T = int((Decimal(SLACK.numerator) / Decimal(SLACK.denominator) * core).to_integral_value(rounding="ROUND_CEILING"))
out["random_mapping_steps_for_0.39"] = str(core)[:40]
out["T_steps"] = T
out["log2_T"] = math.log2(T)
out["predicted_success_at_T_random_mapping"] = 1 - math.exp(-(T * T) / (2 * (n / 6)))

E = F(1) + F(600, C)                     # per step: measured 584 + inversion share 7 + allowance 9
walks, theta_log = 2 ** 32, 40           # parallel walks, distinguished-point rate 2^-40
steps = F(T) * E
tail = F(walks * 2 ** theta_log) * E     # finish trails to their next DP after the budget
dps = -(-T // 2 ** theta_log)
restart = F(dps + walks) * 2000          # fresh start aP0+bT per trail (2 scalar mults, <= 2000 units)
dpstore = F(dps * 600, C)                # DP insert/lookup, 40-byte record
cycle = F(T, 2 ** 10) * F(200, C)        # periodic fruitless-cycle check every 2^10 steps
table = F(2048 * 2000)                   # r=2048 step points R_i = c_i P0 + d_i T
final = F(2 * 20 * 2 ** theta_log) * E + 10 ** 4   # recompute two trails (<= 20/theta each), solve, verify
misc = F(2 ** 20)
W = steps + tail + restart + dpstore + cycle + table + final + misc
out["ledger"] = {k: str(v) for k, v in {
    "per_step_E": E, "steps": steps, "dp_tail": tail, "trail_starts": restart,
    "dp_store": dpstore, "cycle_checks": cycle, "step_table": table, "final_solve": final,
    "decode_constants_output": misc, "W_total": W}.items()}
lw = (Decimal(W.numerator).ln() - Decimal(W.denominator).ln()) / Decimal(2).ln()
num = int((lw * 10**6).to_integral_value(rounding="ROUND_CEILING"))
assert Decimal(num) / Decimal(10**6) - lw > Decimal("1e-30")
out["log2_W"] = str(lw)[:30]
out["submitted_time_log2"] = f"{num}/1000000"
out["log2_memory_bytes"] = math.log2(dps * 48 + 2048 * 64 + walks * 128 + 2 ** 20)
print(json.dumps(out))
