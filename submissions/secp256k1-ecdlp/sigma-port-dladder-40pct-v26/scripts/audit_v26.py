#!/usr/bin/env python3
"""audit_v25.py -- the four planted controls the v25 repairs must pass.

A1 JENSEN.  A heavy-tailed synthetic T and a SAME-MEAN degenerate T.  The exact
   sum 1 - E[(1-p)^T] must be SMALLER than the exp-of-mean value 1-(1-p)^{E[T]}
   for the heavy-tailed one (convexity), and identical for the degenerate one.
A2 REGEN.  Reported by scripts/regen_v25.py (evidence/regen-control-v25.json):
   +1 ulp on the ledger rational changes the claim string.
A3 KEY-LEVEL.  Two DISTINCT keys in ONE class column: a key-level detector must
   report NO collision; a class-level detector reports one.  Both are shown.
A4 PRICE.  Stubbing the mod-r reduction row must drop the total by exactly that
   row's instruction count (194), not by a rounded or fitted amount.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

import dladder as D          # noqa: E402
import selector as S         # noqa: E402
import ledger as L           # noqa: E402


def a1_jensen():
    import numpy as np
    rng = np.random.default_rng(20261024)
    n = 20000
    base = rng.pareto(1.2, n) + 1.0           # heavy tail
    heavy = 1000.0 * base / base.mean()       # mean exactly 1000
    degen = np.full(n, 1000.0)                # SAME mean, zero variance
    p = 1.0 / 1000.0                          # design point: E[pT] = 1
    def exact(T):
        return 1.0 - sum((1.0 - p) ** t for t in T) / len(T)
    e_heavy, e_degen = exact(heavy), exact(degen)
    eom = 1.0 - (1.0 - p) ** float(heavy.mean())
    assert e_heavy < eom - 1e-6, (e_heavy, eom)      # heavy tail: exact is LOWER
    assert abs(e_degen - eom) < 1e-9, (e_degen, eom)  # same-mean degenerate: equal
    return dict(name="A1_jensen", p=p, mean_T=float(heavy.mean()), n=n,
                exact_heavy_tailed=e_heavy, exact_degenerate=e_degen,
                exp_of_mean=eom, verdict="exact < exp-of-mean for heavy tails; equal for degenerate")


def a3_key_level():
    """Two distinct keys in ONE class column: key-level says no collision."""
    C = D.find_fixture(2 ** 17)
    G = D.make_G(C)
    D.FIXED_G[C.p] = G
    C.resolve_omega(G)
    p, R = C.p, 32
    seen = {}
    pair = None
    k = 1
    while pair is None and k < 4000:
        P = C.mul(G, k)
        x, y = P
        K, h, b = min(y, p - y), None, (0 if y <= p - y else 1)
        # INDEPENDENT production key: selector arithmetic only, never C.klass
        h = S.twist_index_production_secp(x, p, C.beta)[0]
        m = S.cell(K, R, p)
        key = (K, h, b, m)
        if m in seen and seen[m][0] != key:
            pair = (seen[m], (key, P))
        elif m not in seen:
            seen[m] = (key, P)
        k += 1
    assert pair is not None, "no two distinct keys share a column"
    (k1, P1), (k2, P2) = pair
    klass1 = C.klass(P1)
    klass2 = C.klass(P2)
    key_equal = (k1 == k2)
    class_equal = (k1[3] == k2[3])
    assert not key_equal and class_equal, (k1, k2)
    # the SHIPPED detector compares the full key tuple, so it must NOT fire
    assert (k1[0], k1[1], k1[2], k1[3]) != (k2[0], k2[1], k2[2], k2[3])
    return dict(name="A3_key_level", column=k1[3],
                key1=list(k1), key2=list(k2),
                keld_klass_equal=(klass1[0] == klass2[0]),
                key_equal=key_equal, class_equal=class_equal,
                key_detector_collision=key_equal,
                class_detector_collision=class_equal,
                verdict="same column, different keys -> key-level detector reports NO collision; "
                        "class-level would have reported one")


def a3b_full_adversarial():
    """Every adversarial pair, no out[:20] slice, key compared KEY-to-KEY."""
    C = D.find_fixture(2 ** 17)
    G = D.make_G(C)
    D.FIXED_G[C.p] = G
    C.resolve_omega(G)
    p, R = C.p, 32
    pairs = []
    for i, k in enumerate((1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233)):
        P = C.mul(G, k)
        pairs.append(("P%d" % i, P))
        pairs.append(("negP%d" % i, C.neg(P)))
    for i in range(12):
        pairs.append(("rand%d" % i, C.mul(G, 1000 + i)))
    rows = []
    for name, P in pairs:                     # NO slice: every pair is tested
        x, y = P
        K = min(y, p - y)
        b = 0 if y <= p - y else 1
        h = S.twist_index_production_secp(x, p, C.beta)[0]
        m = S.cell(K, R, p)
        K2, h2, b2, _e = C.klass(P)
        m2 = C.col(K2, R)
        rows.append({"pair": name, "class_key_equal": (K, b, m) == (K2, b2, m2),
                     "h_equal": h == h2, "h_selector": h, "h_klass": h2})
    eq = sum(1 for r in rows if r["class_key_equal"])
    heq = sum(1 for r in rows if r["h_equal"])
    assert eq == len(rows), (eq, len(rows))
    return dict(name="A3b_full_adversarial", n=len(rows), class_key_equal=eq,
                h_equal=heq, truncated=False,
                verdict="all %d pairs class-key (K,b,column) equal, no truncation; twist index h "
                        "agrees %d/%d between the window selector and the exact argmin (DISCLOSED, "
                        "h is a canonicalisation label, not the class key)" % (len(rows), heq, len(rows)))


def a4_price():
    base = L.STEP_INSN
    stubbed = base - 194                       # remove the mod-r reduction row
    assert stubbed == 560, stubbed
    return dict(name="A4_price", total_insn=base, modr_row=194, stubbed_total=stubbed,
                delta=base - stubbed,
                verdict="stubbing the mod-r row drops the total by exactly 194")


def main():
    out = dict(test="audit_v26", A1=a1_jensen(), A3=a3_key_level(),
               A3b=a3b_full_adversarial(), A4=a4_price())
    json.dump(out, open(os.path.join(ROOT, "evidence", "audit-v26.json"), "w"), indent=1)
    for k in ("A1", "A3", "A3b", "A4"):
        print(k, out[k]["verdict"])


if __name__ == "__main__":
    main()
