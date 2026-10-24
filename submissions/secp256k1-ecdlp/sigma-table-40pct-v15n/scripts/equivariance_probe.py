#!/usr/bin/env python3
"""equivariance_probe.py — does the PRICED step realise the six-class quotient?

Three rules are compared on the package's own toy curve, importing only curve
arithmetic and class_of from scripts/toy_core.py (the walk itself is written here,
independently of toy_core's Trail.step):

  A'  v15l ALIGNED rule: table row = 3 alignment slots of the six-element orbit
      (tbl[j][u] = [u*e_j]P0), index j = (min x over the orbit) & (R-1), slot u =
      alignment of the current point; the state is never canonicalised, so the
      exponent pair stays attached to the stored raw point.
  B   the SHIPPED v15g rule: j = x[0] & (R-1), ONE point per index, even-y canon
      applied to the state (retro-audit TEST B).

A trial is a SPLIT if the orbit partners of one point do not land in the same
class after one step — that is what makes the walk a map on classes, which the
(r-1)/6 birthday model requires.  Also checks S vs -S.

Writes evidence/equivariance-v15l.json.
"""
import json, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import toy_core as tc

cur = tc.find_curve(24, 11)
rnd = random.Random(5)
G = tc.base_point(cur, rnd)
R = 32
tbl = tc.build_table(cur, R, G, rnd)
assert len(tbl) == R and len(tbl[0]) == 6


def step_A(Q):
    m, u = cur.class_of(Q)
    j = m & (R - 1)
    pt, e = tbl[j][u]
    return cur.add(Q, pt), e


def canon_even(Q):
    x, y = Q
    return (x, y if y % 2 == 0 else cur.p - y)


def step_B(Q):
    j = Q[0] & (R - 1)
    pt = tbl[j][0][0]
    return canon_even(cur.add(Q, pt)), tbl[j][0][1]


def probe(step, trials, depth):
    splits = seq = 0
    for _ in range(trials):
        k = rnd.randrange(1, cur.r)
        S = cur.mul(k, G)
        outs = []
        for u in cur.U:
            Q = cur.mul(u, S)
            e = 0
            for _ in range(depth):
                Q, de = step(Q)
                e += de
            outs.append((cur.class_of(Q)[0], e % cur.r))
        splits += len(set(o[0] for o in outs)) != 1
        seq += len(set(outs)) != 1
    return dict(trials=trials, depth=depth, class_splits=splits, label_splits=seq)


def pair(step, trials):
    bad = 0
    for _ in range(trials):
        k = rnd.randrange(1, cur.r)
        S = cur.mul(k, G)
        a = cur.class_of(step(S)[0])[0]
        b = cur.class_of(step(cur.mul(cur.r - 1, S))[0])[0]
        bad += a != b
    return bad


out = dict(kind="equivariance-v15l", curve=dict(p=cur.p, r=cur.r, beta=cur.beta),
           R=R, slots=6,
           A_aligned_depth1=probe(step_A, 600, 1),
           A_aligned_depth8=probe(step_A, 600, 8),
           B_shipped_depth1=probe(step_B, 600, 1),
           B_shipped_depth8=probe(step_B, 600, 8),
           A_aligned_S_vs_negS=pair(step_A, 600),
           B_shipped_S_vs_negS=pair(step_B, 600),
           trials_per_cell=600)
out["verdict"] = dict(
    aligned_rule_is_class_map=(out["A_aligned_depth1"]["class_splits"] == 0 and
                              out["A_aligned_depth8"]["class_splits"] == 0),
    aligned_rule_is_negation_equivariant=(out["A_aligned_S_vs_negS"] == 0),
    shipped_rule_is_not_a_class_map=(out["B_shipped_depth1"]["class_splits"] == 600 and
                                     out["B_shipped_depth8"]["class_splits"] == 600),
    shipped_rule_splits_negation=(out["B_shipped_S_vs_negS"] == 600))
json.dump(out, open(os.path.join(PKG, "evidence", "equivariance-v15l.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
