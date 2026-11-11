"""Output step: P = [x]P0 in G1 by Jacobian double-and-add with counted field operations, Q = Q0, and the
instruction bounds of the output rows (imported by experiment.py). Standard library only."""

from progs import P, R

P0 = (1, 2)                                          # G1 generator; y^2 = x^3 + 3, cofactor 1
Q0 = ((0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed,
       0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2),
      (0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa,
       0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b))
ATE = 0x19d797039be763ba8                            # optimal ate loop count 6x + 2 (protected profile)
FEXP = (P ** 12 - 1) // R                            # final exponent
MUL_OP, ADD_OP, DBL_MUL = 11, 14, 7                  # products and sums per Jacobian step; products per doubling (checked)
CTRL_OP = 1 << 7                                     # bit test, branch and call per Jacobian step or exponent bit
AFF = 4                                              # products of the affine conversion besides the inversion
INV_X, F2M, F2M_P, F2M_S, F2A, F2A_S = 6, 12, 4, 2, 40, 2   # one affine G2 step: inversion extras, Fp2 products
G2_STEP = 2 * P.bit_length() + INV_X + F2M * (F2M_P + F2M_S) + F2A * F2A_S   # Fp operations of one G2 step with its line
FIX = 1 << 12                                        # encoding, comparisons and control of each output row


class Ops:
    """Field arithmetic modulo p that counts products (m) and sums or differences (a)."""
    def __init__(self):
        self.m = self.a = 0

    def mul(self, u, v):
        self.m += 1
        return u * v % P

    def add(self, u, v):
        self.a += 1
        return (u + v) % P

    def sub(self, u, v):
        self.a += 1
        return (u - v) % P


def dbl(o, X, Y, Z):
    """Jacobian doubling for a = 0 (2M + 5S)."""
    A, B = o.mul(X, X), o.mul(Y, Y)
    C = o.mul(B, B)
    t = o.add(X, B)
    D = o.sub(o.sub(o.mul(t, t), A), C)
    D = o.add(D, D)
    E = o.add(o.add(A, A), A)
    X3 = o.sub(o.mul(E, E), o.add(D, D))
    C8 = o.add(C, C); C8 = o.add(C8, C8); C8 = o.add(C8, C8)
    Y3 = o.sub(o.mul(E, o.sub(D, X3)), C8)
    Z3 = o.mul(Y, Z)
    return X3, Y3, o.add(Z3, Z3)


def madd(o, X1, Y1, Z1, x2, y2):
    """Mixed Jacobian-affine addition (7M + 4S); the operands are never equal or opposite (see proof.md)."""
    ZZ = o.mul(Z1, Z1)
    H = o.sub(o.mul(x2, ZZ), X1)
    S2 = o.mul(o.mul(y2, Z1), ZZ)
    HH = o.mul(H, H)
    I = o.add(HH, HH); I = o.add(I, I)
    J = o.mul(H, I)
    rr = o.sub(S2, Y1); rr = o.add(rr, rr)
    V = o.mul(X1, I)
    X3 = o.sub(o.sub(o.sub(o.mul(rr, rr), J), V), V)
    W = o.mul(Y1, J)
    Y3 = o.sub(o.mul(rr, o.sub(V, X3)), o.add(W, W))
    t = o.add(Z1, H)
    Z3 = o.sub(o.sub(o.mul(t, t), ZZ), HH)
    return X3, Y3, Z3


def point_out(x):
    """[x]P0 for 0 <= x < r: left-to-right double-and-add from the top bit, then one inversion by
    exponentiation to p - 2. Returns the affine point (None for x = 0), the Ops counter and the step count."""
    o, steps = Ops(), 0
    o.within = True
    if x == 0:
        return None, o, steps
    X, Y, Z = P0[0], P0[1], 1
    for b in bin(x)[3:]:
        for f, a in ((dbl, ()), (madd, P0)) if b == '1' else ((dbl, ()),):
            m0, a0 = o.m, o.a
            X, Y, Z = f(o, X, Y, Z, *a); steps += 1
            o.within &= o.m - m0 <= (DBL_MUL if f is dbl else MUL_OP) and o.a - a0 <= ADD_OP
    zi = 1
    for b in bin(P - 2)[2:]:
        zi = o.mul(zi, zi)
        if b == '1':
            zi = o.mul(zi, Z)
    z2 = o.mul(zi, zi)
    return (o.mul(X, z2), o.mul(o.mul(Y, z2), zi)), o, steps


def ref_mul(x):
    """Affine reference [x]P0 (uncounted)."""
    acc, add = None, P0
    while x:
        if x & 1:
            acc = add if acc is None else _aadd(acc, add)
        add, x = _aadd(add, add), x >> 1
    return acc


def _aadd(a, b):
    if a[0] == b[0]:
        if (a[1] + b[1]) % P == 0:
            return None
        lam = 3 * a[0] * a[0] * pow(2 * a[1], -1, P) % P
    else:
        lam = (b[1] - a[1]) * pow(b[0] - a[0], -1, P) % P
    x3 = (lam * lam - a[0] - b[0]) % P
    return x3, (lam * (a[0] - x3) - a[1]) % P


def check(xs, fp_mul):
    """Execute the output step for each x: result equals the reference, lies on the curve, field operations
    within the per-step bounds; the executed count stays within the scalar-multiplication bound."""
    rows, ok = [], True
    for x in xs:
        pt, o, steps = point_out(x)
        good = pt == ref_mul(x) and (pt is None or (pt[1] ** 2 - pt[0] ** 3 - 3) % P == 0)
        used = (o.m + o.a) * fp_mul + steps * CTRL_OP + FIX
        good &= o.within and steps <= 2 * (R.bit_length() - 1) and used <= bounds(fp_mul, 0)["point"]
        ok &= good
        rows.append({"x_bits": x.bit_length(), "jacobian_steps": steps, "products": o.m, "sums": o.a,
                     "matches_reference_on_curve": good})
    ok &= ref_mul(R) is None                         # P0 has order r
    return {"runs": rows, "identity_gives_null_and_Q0": point_out(0)[0] is None, "output_step_ok": ok}


def bounds(fp_mul, f12):
    """Every-run instruction bounds of the two output rows, from the executed unit counts."""
    pt = 2 * (R.bit_length() - 1) * ((MUL_OP + ADD_OP) * fp_mul + CTRL_OP) + (2 * P.bit_length() + AFF) * fp_mul + FIX
    na, nf = ATE.bit_length(), FEXP.bit_length()
    pc = (3 * na + 2 + 2 * nf) * f12 + nf * CTRL_OP + (2 * na + 2) * (G2_STEP * fp_mul + FIX) + FIX
    return {"point": pt, "pairing": pc, "ate_bits": na, "fexp_bits": nf}


def doc_lines(fp_mul, f12, pt_log, pc_log):
    """The lines of proof.md and machine.md that state output-step figures, build-filled from the constants
    above and the executed unit counts (experiment.py compares them in every mode)."""
    b, fm = bounds(fp_mul, f12), lambda n: f"{n:,}"
    nr, npb, na, nf = R.bit_length() - 1, P.bit_length(), b["ate_bits"], fm(b["fexp_bits"])
    proof = [
        "Jacobian coordinates (`outstep.py`: a doubling uses %d products, a mixed addition of P0 uses %d)," % (DBL_MUL, MUL_OP),
        "- *point output:* 2^%d for step 8's scalar multiplication: at most 2 x %d Jacobian steps for x < 2^%d," % (pt_log, nr, nr + 1),
        "each at most %d products and %d sums, every field operation charged at the executed `fp_mul` count, plus" % (MUL_OP, ADD_OP),
        "%d per step for the bit test, branch and call; the inversion (at most 2 x %d products) and %d products" % (CTRL_OP, npb, AFF),
        "of the affine conversion; %s for encoding and control:" % fm(FIX),
        "Scalar multiplication: 2 x %d x ((%d + %d) x %d + %d) + (2 x %d + %d) x %d + %s = %s <= 2^%d." % (
            nr, MUL_OP, ADD_OP, fp_mul, CTRL_OP, npb, AFF, fp_mul, fm(FIX), fm(b["point"]), pt_log),
        "- *pairing check:* 2^%d for the attack-side check of step 8. Fp12 products, each charged at the executed" % pc_log,
        "Fp12 product: the Miller loop over the %d bits of 6x + 2 (per bit a squaring and a line product, per" % na,
        "exponentiation by binary powering of (p^12 - 1)/r, %s bits, with %d per exponent bit for the bit" % (nf, CTRL_OP),
        "Fp inversion by exponentiation, two products, one negation), %d Fp2 products (%d products and %d sums" % (F2M, F2M_P, F2M_S),
        "each, enough for the Frobenius images of Q0 too) and %d Fp2 sums (%d sums each), every field operation" % (F2A, F2A_S),
        "charged at `fp_mul`, plus %s per step for embedding the line and control; %s for the comparison" % (fm(FIX), fm(FIX)),
        "One G2 step: 2 x %d + %d + %d x %d + %d x %d = %d field operations." % (
            npb, INV_X, F2M, F2M_P + F2M_S, F2A, F2A_S, G2_STEP),
        "Pairing check: (3 x %d + 2 + 2 x %s) x %s + %s x %d + (2 x %d + 2) x (%d x %d + %s) + %s = %s <= 2^%d." % (
            na, nf, fm(f12), nf, CTRL_OP, na, G2_STEP, fp_mul, fm(FIX), fm(FIX), fm(b["pairing"]), pc_log)]
    machine = [
        "on 1, r - 1 and a seeded x): Jacobian doubling for a = 0 (%d products, %d sums) and mixed addition of P0" % (DBL_MUL, ADD_OP),
        "(%d products, %d sums), counted operation by operation; the affine result must equal an independent" % (MUL_OP, ADD_OP)]
    return proof, machine
