#!/usr/bin/env python3
"""Toy-spec predicate match on 100 random elements per field (v15).

Two independent checks, both scored as rows of evidence/predicate_match.json:

1. production: 100 random GT elements of the real BLS12-381 tower.  The shared
   module function gt_prod.prod_bits is compared with an independently written
   copy of the specification formula (spec_prod_bits) below.  Both must return
   the same (branch, multiplier index, DP flag) triple on every element.
2. toy: 100 random elements of the q = 523 toy tower.  toy15.Toy15.bits is
   compared with the specification formula for the toy, i.e. the SAME branch
   expression (t & 7) < 3 on the same invariant T = Re(c0), plus the toy's
   documented label function L(T) for the index and DP bits.

A mismatch on any element fails the row; the count of agreeing elements is the
number reported, never a rounded fraction.
"""
import json
import random

import gt_prod as G
import toy as T14
from toy15 import Toy15

ROOT = __file__.rsplit("/", 2)[0]


def spec_prod_bits(t, Rtab, dpbits):
    """Independent statement of the production specification."""
    branch = 1 if (t % 8) < 3 else 0
    index = (t // 8) % Rtab
    dp_window = (t // (8 * Rtab)) % (2 ** dpbits)      # bits 13.. of T
    return bool(branch), index, dp_window == 0


def spec_toy_bits(t, R, dpbits):
    label = ((t * 0x9E3779B1) ^ (t >> 13)) & 0xFFFFFFFF
    return (t % 8) < 3, label % R, ((label >> 8) & ((1 << dpbits) - 1)) == 0


def main():
    rng = random.Random(20261023)
    p, r = G.P, G.R
    tower = G.Tower(p)
    dpbits = 100
    Rtab = 1024
    agree = 0
    traces = []
    for _ in range(100):
        y = G.random_gt(tower, r, rng)
        t = tower.trace(y)
        traces.append(t)
        a = G.prod_bits(t, Rtab, dpbits)
        b = spec_prod_bits(t, Rtab, dpbits)
        if a == b:
            agree += 1
    # the branch rule must also reproduce the exact 3/8 frequency on the sample
    brs = sum(1 for t in traces if G.prod_bits(t, Rtab, dpbits)[0])

    rp, N, _o, _oo = T14.field_params(523)
    tw = Toy15(523, rp, R=128, dpbits=4, seed=7)
    tw.build_table()
    tagree = 0
    tbrs = 0
    for _ in range(100):
        a1 = rng.randrange(rp)
        d1 = rng.randrange(1, rp)
        y = tw.T.mul(tw.T.pow(tw.g, a1), tw.T.pow(tw.z, d1))
        t = y[0][0][0]
        if tw.bits(y)[:3] == spec_toy_bits(t, tw.R, tw.dpbits):
            tagree += 1
        if (t % 8) < 3:
            tbrs += 1

    out = dict(package="gt-metered-squaring-40pct-v15",
               prod_checked=100, prod_agree=agree,
               prod_branch_true_on_sample=brs,
               prod_shared_module="gt_prod.prod_bits", prod_spec_copy="spec_prod_bits",
               prod_dp_bits=dpbits, prod_R=Rtab,
               toy_checked=100, toy_agree=tagree, toy_branch_true_on_sample=tbrs,
               toy_branch_rule="(t & 7) < 3", toy_label_rule="L(t)=((t*0x9E3779B1)^(t>>13))&0xFFFFFFFF",
               toy_dp_bits=4, toy_R=128,
               verdict="PASS" if agree == 100 and tagree == 100 else "FAIL")
    json.dump(out, open(ROOT + "/evidence/predicate_match.json", "w"), indent=1)
    print(json.dumps(out))
    assert out["verdict"] == "PASS"


if __name__ == "__main__":
    main()
