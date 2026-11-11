"""Pairing output Q = [x]Q0 on the BN254 G2 twist, with counted Fp operations.

Given the verified logarithm x (g^x = z, proof.md step 4), the algorithm outputs P = P0 and Q = [x]Q0;
for z = 1 it outputs (P0, infinity). e(P0, [x]Q0) = e(P0, Q0)^x = g^x = z by bilinearity.
This script (standard library, own code, under 2 s):
  - checks Q0 and the twist coefficient b' = 3/(9+u) from the protected profile constants, and [r]Q0 = O;
  - computes [x]Q0 by left-to-right double-and-add in Jacobian coordinates over Fp2 (dbl-2009-l with a = 0,
    madd-2007-bl), one Fp2 inversion (norm and a Fermat inversion in Fp) and the affine conversion,
    counting every Fp multiplication, addition and subtraction;
  - compares each result with an independent affine reference, checks it on the twist, and checks
    [r]Q = O on a sample;
  - prices the counts with machine.md's fp_mul 155 + 4 (call) and fp_add/fp_sub 21 + 4, plus 32 per
    routine call, and gives the analytical worst case over all x < r (253 doublings and 253 additions).
Optional stdin JSON: {"seed": int, "cases": int}.
"""
import json
import random
import sys

try:
    cfg = json.loads(sys.stdin.read() or "{}")
except ValueError:
    cfg = {}
SEED, CASES = int(cfg.get("seed", 1)), int(cfg.get("cases", 24))
rng = random.Random(SEED)

p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
r = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
Q0 = ((0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed,
       0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2),
      (0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa,
       0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b))
MUL_I, ADD_I, CALL_I = 155 + 4, 21 + 4, 32


class F2:
    """Fp2 = Fp[u]/(u^2 + 1) with counted Fp operations (Karatsuba product: 3 mul, 5 add/sub)."""
    def __init__(self):
        self.m = self.a = self.calls = 0

    def mul(self, x, y):
        self.m += 3; self.a += 5; self.calls += 1
        t0, t1 = x[0] * y[0] % p, x[1] * y[1] % p
        return ((t0 - t1) % p, ((x[0] + x[1]) * (y[0] + y[1]) - t0 - t1) % p)

    def sqr(self, x):
        self.m += 2; self.a += 3; self.calls += 1
        return ((x[0] + x[1]) * (x[0] - x[1]) % p, 2 * x[0] * x[1] % p)

    def add(self, x, y):
        self.a += 2
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

    def sub(self, x, y):
        self.a += 2
        return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)

    def inv(self, x):                  # 1/x = conj(x) / (x0^2 + x1^2), Fermat in Fp
        self.calls += 1
        n = (x[0] * x[0] + x[1] * x[1]) % p
        self.m += 2; self.a += 1
        e = p - 2
        self.m += (e.bit_length() - 1) + bin(e).count("1") - 1
        ni = pow(n, e, p)
        self.m += 2; self.a += 1
        return (x[0] * ni % p, (-x[1]) * ni % p)

    def instr(self):
        return self.m * MUL_I + self.a * ADD_I + self.calls * CALL_I


def f2inv_plain(x):
    n = (x[0] * x[0] + x[1] * x[1]) % p
    ni = pow(n, p - 2, p)
    return (x[0] * ni % p, (-x[1]) * ni % p)


def f2mul_plain(x, y):
    return ((x[0] * y[0] - x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)


B2 = f2mul_plain((3, 0), f2inv_plain((9, 1)))


def on_twist(P):
    x, y = P
    return f2mul_plain(y, y) == tuple((a + b) % p for a, b in zip(f2mul_plain(f2mul_plain(x, x), x), B2))


def aff_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1][0] + Q[1][0]) % p == 0 and (P[1][1] + Q[1][1]) % p == 0:
            return None
        lam = f2mul_plain(f2mul_plain((3, 0), f2mul_plain(P[0], P[0])), f2inv_plain(f2mul_plain((2, 0), P[1])))
    else:
        lam = f2mul_plain(tuple((a - b) % p for a, b in zip(Q[1], P[1])),
                          f2inv_plain(tuple((a - b) % p for a, b in zip(Q[0], P[0]))))
    x3 = tuple((a - b - c) % p for a, b, c in zip(f2mul_plain(lam, lam), P[0], Q[0]))
    y3 = tuple((a - b) % p for a, b in zip(f2mul_plain(lam, tuple((u - v) % p for u, v in zip(P[0], x3))), P[1]))
    return (x3, y3)


def aff_mul(k, P):
    R = None
    while k:
        if k & 1:
            R = aff_add(R, P)
        P = aff_add(P, P)
        k >>= 1
    return R


def jac_mul(k, P, F):
    X, Y, Z = P[0], P[1], (1, 0)
    for bit in bin(k)[3:]:
        F.calls += 1                                   # dbl-2009-l (a = 0): 2M + 5S
        A = F.sqr(X); B = F.sqr(Y); Cc = F.sqr(B)
        t = F.sub(F.sub(F.sqr(F.add(X, B)), A), Cc); D = F.add(t, t)
        E = F.add(F.add(A, A), A); Fv = F.sqr(E)
        X3 = F.sub(Fv, F.add(D, D))
        C8 = F.add(Cc, Cc); C8 = F.add(C8, C8); C8 = F.add(C8, C8)
        Y3 = F.sub(F.mul(E, F.sub(D, X3)), C8)
        Z3 = F.mul(Y, Z); Z3 = F.add(Z3, Z3)
        X, Y, Z = X3, Y3, Z3
        if bit == "1":
            F.calls += 1                               # madd-2007-bl: 7M + 4S
            Z1Z1 = F.sqr(Z); U2 = F.mul(P[0], Z1Z1); S2 = F.mul(P[1], F.mul(Z, Z1Z1))
            H = F.sub(U2, X)
            assert H != (0, 0), "exceptional addition"
            HH = F.sqr(H); I = F.add(HH, HH); I = F.add(I, I); J = F.mul(H, I)
            rr = F.sub(S2, Y); rr = F.add(rr, rr); V = F.mul(X, I)
            X3 = F.sub(F.sub(F.sqr(rr), J), F.add(V, V))
            t = F.mul(Y, J); t = F.add(t, t)
            Y3 = F.sub(F.mul(rr, F.sub(V, X3)), t)
            Z3 = F.sub(F.sub(F.sqr(F.add(Z, H)), Z1Z1), HH)
            X, Y, Z = X3, Y3, Z3
    zi = F.inv(Z); zi2 = F.sqr(zi); zi3 = F.mul(zi2, zi)
    return (F.mul(X, zi2), F.mul(Y, zi3))


def main():
    ok = on_twist(Q0) and aff_mul(r, Q0) is None
    edge = [1, 2, 3, r - 1, r - 2, (r - 1) // 2, 2 ** 253 - 1, 3 * 2 ** 252 - 1]
    xs = [x for x in edge if 1 <= x < r] + [rng.randrange(1, r) for _ in range(CASES)]
    fails, worst = 0, 0
    for i, x in enumerate(xs):
        F = F2()
        Q = jac_mul(x, Q0, F)
        fails += Q != aff_mul(x, Q0) or not on_twist(Q)
        if i < 3:
            fails += aff_mul(r, Q) is not None
        worst = max(worst, F.instr())
    # analytical worst case: the per-call counts are input independent
    Fd, Fa = F2(), F2()
    jac_mul(2, Q0, Fd)                  # one doubling + inversion + affine
    jac_mul(3, Q0, Fa)                  # one doubling + one addition + inversion + affine
    dbl_and_tail, add_only = Fd.instr(), Fa.instr() - Fd.instr()
    F1 = F2(); jac_mul(1, Q0, F1); tail = F1.instr()
    per_dbl = Fd.instr() - tail
    analytic = 253 * (per_dbl + add_only) + tail + 254 * 12
    # Three schoolbook Fp2 products: 12 Fp mul, 8 Fp add/sub incl. b',
    # three routine calls and 128 comparison/branch/outer-call instructions.
    # 12*159 + 8*25 + 3*32 + 128 = 2332 < the conservative 4096 charge.
    oncurve = 1 << 12
    serial = 6 * (MUL_I + 64 * 8)                             # P0 and Q: 6 coordinates out of Montgomery form, hex
    charged = 1 << 22
    out = {"seed": SEED, "Q0_on_twist_order_r": ok, "cases": len(xs), "failures": fails,
           "per_doubling_instr": per_dbl, "per_addition_instr": add_only, "inversion_affine_instr": tail,
           "analytic_worst_instr": analytic, "on_curve_check_instr": oncurve, "serialisation_instr": serial,
           "sample_max_instr": worst, "charged_instr": charged,
           "charge_ok": analytic + oncurve + serial <= charged and worst <= analytic,
           "identity_output": ["P0", None]}
    out["all_ok"] = ok and fails == 0 and out["charge_ok"]
    print(json.dumps(out, indent=1))
    sys.exit(0 if out["all_ok"] else 1)


main()
