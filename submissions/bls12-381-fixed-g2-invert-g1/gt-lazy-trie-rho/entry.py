"""The child's 54 mul.s entry vectors (kernels/gen.py entry_ops): the per-op price charged in a trie node's
construction (proof.md section 5, experiment.py part D) and the replay mod p that parts K and C feed to mul.s."""

OP_ITEMS = {"3 + 3 ldp of the operands (a negation: 0 - b, one)": 6, "6-limb add or sub chain, carry word": 7,
            "3 ldp of p; subs/sbcs of p with borrow (add) or 6 ands of p by the borrow mask (sub)": 10,
            "6 csel (add) or 6-limb add of the masked p (sub)": 6, "3 stp, 2 address": 5, "margin": 11}
OP_INSTR = sum(OP_ITEMS.values())
OPS, COPIES = 113, 2
BUILD_INSTR = OPS * OP_INSTR + COPIES * 8
assert (OP_INSTR, BUILD_INSTR) == (45, 5101) and BUILD_INSTR <= 54 * 200


def mul_entry(lay, v, p):
    """a + s b (a = -1: zero) with one correction by p, as -p < a + s b < 2p is asserted; each of the 54 outputs is
    checked against its linear form in the 12 coordinates mod p."""
    v = list(v)
    for a, b, s in lay["entry_ops"]:
        t = (v[a] if a >= 0 else 0) + s * v[b]
        assert -p < t < 2 * p
        v.append(t - p if t >= p else t + p if t < 0 else t)
    out = [v[i] for i in lay["entry_index"]]
    env = {"m%d" % i: v[i] for i in range(12)}
    assert out == [sum(c * env[n] for n, c in L.items()) % p for L in lay["entry"]]
    return out
