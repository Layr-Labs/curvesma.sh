"""Exact-integer O_H-basis saturation selector for a special-q module L_q from the 99/100-LLL rows R[0..15]:
Sets (a1, b1) = R[0] and scans R[1..15] in pure exact integer arithmetic for the first row (a2, b2) satisfying
|N_{H/Q}(a1*b2 - b1*a2)| == q (i.e. [L_q : O_H v1 + O_H v2] == 1, so M == L_q identically).
Cost: <= 15 exact O_H multiplies + norms (< 2^18 Ir << LLL_NODE = 2^40 Ir; ZERO floating-point fallback!).
If no row in R[1..15] achieves |N(det)| == q, returns None so the node deterministically treats the attempt
as a certificate failure (fail_cert) and draws a fresh uniform boot exponent e."""
import audit_cm as _AC


def ob_pair(q, ry, rt, R):
    Z = _AC
    a1, b1 = R[0][:8], R[0][8:]
    for v in R[1:]:
        d = Z.sub_phi15_z(Z.mul_phi15_z(a1, v[8:]), Z.mul_phi15_z(b1, v[:8]))
        if any(d) and abs(Z.norm_q15_z(d)) == q:
            return a1, b1, v[:8], v[8:]
    return None
