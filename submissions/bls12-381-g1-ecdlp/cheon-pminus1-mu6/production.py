"""Production-constant checks and the exact W(2^20) bound for proof.md.

Verifies on the protected BLS12-381 constants: the factorisation of r-1,
6d | r-1 for d=2^20, -1 and lambda in H = (F_r^*)^d, the GLV relation
phi(P0)=(beta*x, y)=[lambda]P0 on the real generator, m*g >= 0.39*N', and
recomputes every ledger term of proof.md as exact rationals.
"""
import json
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

json.load(sys.stdin)

p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P0 = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
      0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
FACT = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2,
        2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}
LOG_N, D = 20, 1 << 20
C = 2000


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        if n % q == 0:
            return n == q
    dd, s = n - 1, 0
    while dd % 2 == 0:
        dd //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        x = pow(a, dd, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
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
    R = None
    while k:
        if k & 1:
            R = add(R, P)
        P = add(P, P)
        k >>= 1
    return R


out = {}
prod = 1
for q, e in FACT.items():
    assert is_prime(q)
    prod *= q ** e
assert prod == r - 1
out["r_minus_1_factorisation_verified"] = True
out["log2_r"] = math.log2(r)
assert (r - 1) % (6 * D) == 0
N = (r - 1) // D
Np = N // 6
out["six_d_divides_r_minus_1"] = True
out["log2_N_prime"] = math.log2(Np)

zeta = next(g for g in range(2, 100) if all(pow(g, (r - 1) // q, r) != 1 for q in FACT))
out["zeta_generator_of_Fr_star"] = zeta
assert pow(r - 1, N, r) == 1
lams = [pow(zeta, (r - 1) // 3, r), pow(zeta, 2 * (r - 1) // 3, r)]
for l in lams:
    assert pow(l, N, r) == 1
out["minus_one_and_cube_roots_of_unity_in_H"] = True

assert P0[1] ** 2 % p == (P0[0] ** 3 + 4) % p
assert smul(r, P0) is None
betas = [b for b in (pow(a, (p - 1) // 3, p) for a in (2, 3, 5, 7)) if b != 1][:1]
beta = betas[0]
phi = (beta * P0[0] % p, P0[1])
lam = next(l for l in lams if smul(l, P0) == phi)
beta2 = beta * beta % p
assert smul(lam * lam % r, P0) == (beta2 * P0[0] % p, P0[1])
assert (1 + beta + beta2) % p == 0
out["glv_phi_equals_lambda_on_P0"] = True
out["class_x_set_is_x_beta_x_minus_x_minus_beta_x"] = True

need = -(-39 * Np // 100)
m = math.isqrt(need)
if m * m < need:
    m += 1
g = m
assert m * g * 100 >= 39 * Np and m * g <= Np
out["m"] = m
out["g"] = g
out["log2_m"] = math.log2(m)
out["m_times_g_times_100_ge_39_N_prime"] = m * g * 100 >= 39 * Np
out["success_probability_m_g_over_N_prime"] = m * g / Np

# Ledger (target-group-operation equivalents, C = 2000 instructions).
E = F(2) + F(800, C)                        # per BSGS element: measured 734 + inversion share 13 + allowance 53
bsgs = (m + g) * E
tables = 6 * ((1 << 85) - 1) * (F(1) + F(400, C)) + 10**4
slow_path = F(2) ** 38
recover = F(510 * 250, C) + 6 * D * 3
misc = F(2) ** 20
hclear = F((1 << 117) * 10, C)           # clear the 2^117-slot hash table, 10 words per slot
W = bsgs + tables + slow_path + recover + misc + hclear
out["ledger"] = {k: str(v) for k, v in
                 {"per_element_E": E, "bsgs": bsgs, "tables": tables, "slow_path": slow_path,
                  "recovery": recover, "decode_constants_output": misc, "hash_table_clear": hclear, "W_total": W}.items()}
getcontext().prec = 80
lw = (Decimal(W.numerator).ln() - Decimal(W.denominator).ln()) / Decimal(2).ln()
num = int((lw * 10**6).to_integral_value(rounding="ROUND_CEILING"))
assert Decimal(num) / Decimal(10**6) - lw > Decimal("1e-30")
out["log2_W_80_digits"] = str(lw)[:30]
out["submitted_time_log2"] = f"{num}/1000000"
pre = tables + m * E + hclear
lpre = (Decimal(pre.numerator).ln() - Decimal(pre.denominator).ln()) / Decimal(2).ln()
out["log2_preprocessing_tables_plus_baby_phase"] = str(lpre)[:12]
mem = (1 << 117) * 80 + 6 * (1 << 85) * 96 + (1 << 20)
out["log2_memory_bytes"] = math.log2(mem)
out["hash_load_factor"] = m / (1 << 117)
print(json.dumps(out))
