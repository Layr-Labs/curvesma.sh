"""Toy witness: canonicalization leaves a negation 2-cycle; a mod-3 phase kills it.

Scope is the additive group Z/13Z with a GLV unit, not secp256k1, and not a
cost measurement. The production map and the integer work bound are in proof.md.

On y^2 = x^3 + 7 the order-6 automorphism orbit is represented here by
multiplication by the six units of Z[λ]/(λ^2+λ+1). The unphased adding map
with a single jump V=10 swaps the two canonical residues 1 and 2, which is
the structural identity can(can(C+V)+V)=C with both canonical units equal
to -1. The phased map uses the same jump but advances a counter modulo 3,
so no state satisfies F(F(s))=s.
"""

import json
import sys

Q = 13
LAM = 3  # 3^2 + 3 + 1 = 13
# 1, λ, λ^2, -1, -λ, -λ^2
UNITS = (1, 3, 9, 12, 10, 4)
JUMP = 10  # hand-checked: 1 <-> 2 under can(· + JUMP)


def canonicalize(v, q=Q):
    best_c = None
    best_u = 1
    best_i = 0
    for i, u in enumerate(UNITS):
        c = (u * v) % q
        if best_c is None or c < best_c or (c == best_c and i < best_i):
            best_c = c
            best_u = u
            best_i = i
    return best_c, best_u


def unphased(c, q=Q):
    return canonicalize((c + JUMP) % q, q)[0]


def phased(state, q=Q):
    c, i = state
    c2, _u = canonicalize((c + JUMP) % q, q)
    return (c2, (i + 1) % 3)


def solve(left, right, q):
    """left/right are (U, A, B). Returns k with den*k ≡ num, or None."""
    u, a, b = left
    u2, a2, b2 = right
    den = (u * a - u2 * a2) % q
    num = (u2 * b2 - u * b) % q
    if den == 0:
        return None
    return (num * pow(den, q - 2, q)) % q


def run(modulus):
    if modulus != Q:
        raise SystemExit("this witness is pinned to modulus 13")
    if (LAM * LAM + LAM + 1) % Q != 0:
        raise SystemExit("lambda")
    if len(set(UNITS)) != 6:
        raise SystemExit("units")

    # Structural 2-cycle of the unphased canonical walk.
    if unphased(1) != 2 or unphased(2) != 1:
        raise SystemExit("expected unphased negation cycle 1 <-> 2")
    c_at_11, u_at_11 = canonicalize((1 + JUMP) % Q)
    if c_at_11 != 2 or u_at_11 != Q - 1:
        raise SystemExit("first canonical unit is not negation")

    canonical = []
    for v in range(1, Q):
        c, _u = canonicalize(v)
        if c == v:
            canonical.append(c)
    if canonical != [1, 2]:
        raise SystemExit("unexpected canonical set")

    unphased_cycles = 0
    for c in canonical:
        img = unphased(c)
        if img != c and unphased(img) == c:
            unphased_cycles += 1

    phased_cycles = 0
    states = [(c, i) for c in canonical for i in range(3)]
    for s in states:
        if phased(phased(s)) == s:
            phased_cycles += 1

    # Planted collision: Q_pt = k*P with P=2, k=5, Q_pt=10.
    # Frame (U,A,B) represents U*(A*Q_pt + B*P).
    p_pt = 2
    k = 5
    q_pt = (k * p_pt) % Q
    frame_q = (1, 1, 0)  # 1*(1*Q + 0*P) = Q
    frame_p = (1, 0, k)  # 1*(0*Q + k*P) = Q
    got = solve(frame_q, frame_p, Q)
    if got != k:
        raise SystemExit("discrete log not recovered")
    # Same frames must describe the same group element.
    def val(frame):
        u, a, b = frame
        return (u * ((a * q_pt + b * p_pt) % Q)) % Q

    if val(frame_q) != val(frame_p):
        raise SystemExit("planted frames disagree")

    return {
        "group": "Z/13Z with GLV units",
        "modulus": Q,
        "unphased_proper_2cycles": unphased_cycles,
        "phased_proper_2cycles": phased_cycles,
        "recovered_k": got,
        "scope": "toy group only; does not validate production cost bounds",
    }


def main():
    spec = json.load(sys.stdin)
    result = run(spec["modulus"])
    if result["unphased_proper_2cycles"] < 1:
        raise SystemExit("canonicalization removed the negation cycle")
    if result["phased_proper_2cycles"] != 0:
        raise SystemExit("phase counter left a 2-cycle")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
