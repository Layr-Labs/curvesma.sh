"""O_H-basis of a special-q module from the exact integer LLL rows R (DESCENT.md section 1):
(a1, b1) = R[0]; if the first O_H-independent later row gives |N(a1 b2 - a2 b1)| = q, that pair is an O_H-basis
(its determinant generates P); otherwise the basis is completed constructively by obasis.basis (generator of P by the
proven Fincke-Pohst bound, content division, exact Hermite complement, Babai reduction)."""
import audit_cm as _AC
import obasis as _OB


def ob_pair(q, ry, rt, R):
    Z = _AC
    a1, b1 = R[0][:8], R[0][8:]
    for v in R[1:]:
        d = Z.sub_phi15_z(Z.mul_phi15_z(a1, v[8:]), Z.mul_phi15_z(b1, v[:8]))
        if any(d):
            if abs(Z.norm_q15_z(d)) == q:
                return a1, b1, v[:8], v[8:]
            break
    pi = _OB.prime_generator(q, ry, Z)
    return _OB.basis(q, ry, rt, (a1, b1), Z, pi)
