#!/usr/bin/env python3
"""dp_sweep_v26.py -- the DIRECT-COUNT success row.

A trial counts as a recovery iff the walk decodes a candidate xx with
[xx]G == Q (dladder.run_trial's own gate).  The CP lower bound is computed from
that count alone, in the cell that produced it.  No transform, no ratio, no
pooling, no exp-of-mean, no stopping-time statistic.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
os.chdir(ROOT)

import dladder as D  # noqa: E402


def main():
    trials = 1200
    out = {"test": "dp_density_direct_count_v26",
           "law": "success = DIRECTLY COUNTED verified recoveries ([x]G==Q) / trials; "
                  "CP-lo from that count alone; no transform, no ratio, no pooling",
           "cells": []}
    for d in (0, 4, 8, 12, 16):
        r = D.cell(2 ** 17, trials, 12, dgate=d, L=16, R=32)
        n, k = r["trials"], r["successes"]
        out["cells"].append({"target_classes": 2 ** 17, "L": 16, "R": 32, "dgate": d,
                             "trials": n, "recovered": k, "rate": k / n,
                             "cp95_lower": D.cp_lower(k, n) if k else 0.0,
                             "visits": r["visits"], "N_CLASSES": r["N_CLASSES"],
                             "classes_log2": r["classes_log2"]})
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
