#!/usr/bin/env python3
"""Run the v24 censuses (family + unconditional recovery) and rewrite
evidence/census-v24.json.  Deterministic: per-worker seeds start at
BASE_SEED + 1000*worker, trial counts are the ones reported in the results."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v24_measure
v24_measure.main()
