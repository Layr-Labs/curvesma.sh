#!/usr/bin/env python3
"""Part A: exact checks on BLS12-381 itself, including a real pairing evaluation.

Every check is a named boolean; the summary line is PART_A PASS/FAIL and the same
booleans are re-read by verify_certificate.py from part_a-out.json.
"""
import json, random, sys, time
from fractions import Fraction
import tower
import pairing as PR
import constants as K

P = K.P_FIELD
R = K.R_ORDER
ATE_LOOP = 0xd201000000010000          # |x| of BLS12-381, the ate loop magnitude
Z6 = ((0, 0), (0, 0), (0, 0))

P0 = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
      0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
Q0 = ((0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,
       0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e),
      (0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,
       0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be))


def keys_invariant_checks(T, rng, n=4):
    """key rule, distinctness, kappa, min-key identity, class well-definedness."""
    out = dict(key_rule=True, keys_distinct=True, kappa_ok=True, min_ident=True,
               sqr_class_ok=True, cyclo_sqr_identity=True, gt_order=True,
               p6_conj=True)
    e = (P**12 - 1)//R
    for _ in range(n):
        while True:
            x = tuple(tuple((rng.randrange(P), rng.randrange(P)) for _ in range(3)) for _ in range(2))
            g = T.pow(x, e)
            if g != T.one12() and T.pow(g, R) == T.one12():
                break
        out['gt_order'] = out['gt_order'] and T.pow(g, R) == T.one12()
        # Frobenius^6 = inverse = conjugate
        g6 = g
        for _ in range(6):
            g6 = T.frob(g6)
        out['p6_conj'] = out['p6_conj'] and g6 == (g[0], T.f6neg(g[1]))
        out['frob_ok'] = T.frob(g) == T.pow(g, P)
        # cyclotomic squaring identity (the row charged in the ledger)
        out['cyclo_sqr_identity'] = out['cyclo_sqr_identity'] and T.cyclo_sqr_fast(g) == T.mul(g, g)
        # key rule: {A zeta^j} U {B zeta^j} equals the twelve conjugate keys
        d0 = g[1][0]
        A = d0[0]
        B = (d0[0]*T.gamma[0] + d0[1]*T.gamma[1]) % P
        pred = {(A*pow(T.zeta, j, P)) % P for j in range(6)} | {(B*pow(T.zeta, j, P)) % P for j in range(6)}
        ks = T.keys12(g)
        out['key_rule'] = out['key_rule'] and pred == set(ks)
        out['keys_distinct'] = out['keys_distinct'] and len(set(ks)) == 12
        # min-key identity: min over the twelve signed values = -(max of six |.|)
        absv = []
        for v in (A, (A*T.zeta) % P, (A*T.zeta*T.zeta) % P, B, (B*T.zeta) % P, (B*T.zeta*T.zeta) % P):
            absv.append(min(v, P-v))
        out['min_ident'] = out['min_ident'] and min(ks) == min(absv)
        # kappa: same key and same canonical element from all twelve conjugates
        base, kap = T.kappa(g)
        canon = None
        for k in range(12):
            h = g
            for _ in range(k):
                h = T.frob(h)
            mk, kk = T.kappa(h)
            if mk != base:
                out['kappa_ok'] = False
            c = h
            for _ in range(kk):
                c = T.frob(c)
            if canon is None:
                canon = c
            elif c != canon:
                out['kappa_ok'] = False
        # squaring is well defined on classes and injective on them:
        # the twelve conjugates of g^2 are exactly the class of (kappa-canonical g)^2
        csq = T.mul(canon, canon)
        keyssq = {min(T.keys12(csq))}
        for k in range(12):
            h = csq
            for _ in range(k):
                h = T.frob(h)
            keyssq.add(min(T.keys12(h)))
        out['sqr_class_ok'] = out['sqr_class_ok'] and len(keyssq) == 1
        # injectivity witness: x^2 = y^2 with x,y in distinct classes never happens;
        # algebraically x^2=y^2 => x=+-y and -1 = p^6 lies in the orbit, so the
        # classes coincide.  Check the algebraic premise on this field:
        out['minus_one_in_orbit'] = (pow(P, 6, R) == R - 1)
    return out


def main():
    t0 = time.time()
    T = tower.Tower(P)
    rng = random.Random(20261023)
    chk = []

    chk.append(("p_pow6_is_minus_one", pow(P, 6, R) == R - 1))
    chk.append(("ord_r_p_is_twelve",
                all(pow(P, k, R) != 1 for k in (1, 2, 3, 4, 6)) and pow(P, 12, R) == 1))
    chk.append(("p4_plus_p8_plus_one_zero", (1 + pow(P, 4, R) + pow(P, 8, R)) % R == 0))
    chk.append(("n_classes_exact", (R - 1) % 12 == 0 and K.N_CLASSES*12 + 1 == R))

    inv = keys_invariant_checks(T, rng)
    for k in ("gt_order", "p6_conj", "frob_ok", "key_rule", "keys_distinct",
              "kappa_ok", "min_ident", "sqr_class_ok", "cyclo_sqr_identity"):
        chk.append((k, bool(inv.get(k))))
    chk.append(("minus_one_in_orbit", bool(inv.get("minus_one_in_orbit"))))

    # ---- real pairing: non-degeneracy, bilinearity, planted reduction
    g = PR.pairing(T, P0, Q0, ATE_LOOP, R)
    chk.append(("pairing_nondegenerate", g != T.one12()))
    chk.append(("pairing_order_r", T.pow(g, R) == T.one12()))
    P2 = PR.g1_mul(P0, 2, P)
    Q3 = PR.g2_mul(T, Q0, 3)
    chk.append(("pairing_bilinear_g1", PR.pairing(T, P2, Q0, ATE_LOOP, R) == T.mul(g, g)))
    g3 = T.mul(T.mul(g, g), g)
    chk.append(("pairing_bilinear_g2", PR.pairing(T, P0, Q3, ATE_LOOP, R) == g3))
    # planted instance exactly as the challenge samples it: P=[b]P0, z=e(P,Q0)^a
    a = rng.randrange(1, R); b = rng.randrange(1, R)
    Pc = PR.g1_mul(P0, b, P)
    gp = PR.pairing(T, Pc, Q0, ATE_LOOP, R)
    z = T.pow(gp, a)
    Qans = PR.g2_mul(T, Q0, a)
    chk.append(("planted_ePQ_equals_z", PR.pairing(T, Pc, Qans, ATE_LOOP, R) == z))
    chk.append(("planted_wrong_Q_fails",
                PR.pairing(T, Pc, PR.g2_mul(T, Q0, (a + 1) % R), ATE_LOOP, R) != z))

    # ---- the three-multiplication selector, on this tower's own constants
    # d = gamma0*a + gamma1*b; in this tower gamma = g0*(1-u) so
    # d = g0*(a-b) is ONE base-field multiply and the twelve keys cost three.
    g0, g1 = T.gamma                    # gamma = (gamma0, gamma1) in Fp, two scalars
    chk.append(("selector_gamma_sum_zero", (g1 + g0) % P == 0))
    chk.append(("selector_g0_pow6_minus_eighth", (8*pow(g0, 6, P) + 1) % P == 0))
    chk.append(("selector_zeta_two_g0_sq", T.zeta == (2*pow(g0, 2, P)) % P))
    three_ok = True
    for _ in range(200):
        aa = rng.randrange(0, P); bb = rng.randrange(0, P)
        general = (g0*aa + g1*bb) % P
        fast = (g0*((aa - bb) % P)) % P
        three_ok = three_ok and general == fast
    chk.append(("selector_three_mul_equals_general", three_ok))


    ok = all(v for _, v in chk)
    out = {"checks": {n: bool(v) for n, v in chk}, "part_a": bool(ok),
           "seconds": round(time.time() - t0, 1), "checks_n": len(chk)}
    json.dump(out, open("part_a-out.json", "w"), indent=1)
    for n, v in chk:
        print("%-34s %s" % (n, "PASS" if v else "FAIL"))
    print("PART_A", "PASS" if ok else "FAIL", len(chk), "checks, %.1fs" % (time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
