"""BLS12-381 G1 arithmetic on standard integers (affine, y^2 = x^3 + 4), stdlib only.

Used by gen_inputs_v3.py and replica_v3.py. The group law is the same one blst computes, so
the kernel's affine outputs must equal these outputs exactly.
"""
import json


def load_profile(path):
    prof = json.load(open(path))
    p = int(prof["field_modulus"], 16)
    r = int(prof["subgroup_order"], 16)
    g = (int(prof["generator_P0"]["x"], 16), int(prof["generator_P0"]["y"], 16))
    return p, r, g


def on_curve(P, p):
    x, y = P
    return (y * y - x * x * x - 4) % p == 0


def add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def mul(k, P, p):
    R = None
    Q = P
    while k:
        if k & 1:
            R = add(R, Q, p)
        Q = add(Q, Q, p)
        k >>= 1
    return R


def cube_roots_of_unity(m):
    out = []
    for g in range(2, 200):
        c = pow(g, (m - 1) // 3, m)
        if c != 1 and pow(c, 3, m) == 1 and c not in out:
            out.append(c)
            if len(out) == 2:
                break
    return out


def matched_beta_lambda(P0, p, r):
    """Pick beta in Fp and lambda in Fr with (beta*x, y) = [lambda]P0 (the GLV endomorphism)."""
    x0, y0 = P0
    for lam in cube_roots_of_unity(r):
        target = mul(lam, P0, p)
        for beta in cube_roots_of_unity(p):
            if (beta * x0 % p, y0) == target:
                return beta, lam
    raise ValueError("no matching beta, lambda")


def mont_limbs(v, p):
    m = v * (1 << 384) % p
    return [(m >> (64 * i)) & ((1 << 64) - 1) for i in range(6)]


def less_mont(a, b, p):
    la, lb = mont_limbs(a, p), mont_limbs(b, p)
    for i in range(5, -1, -1):
        if la[i] != lb[i]:
            return la[i] < lb[i]
    return False
