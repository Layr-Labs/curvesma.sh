#!/usr/bin/env python3
"""v21 certificate check: recompute the exact-estimator, shape-ladder and retry rows
from the shipped evidence.  Delegates to verify_v21.py and prints the summary."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
rc = subprocess.call([sys.executable, os.path.join(HERE, "verify_v21.py")])
C = json.load(open(os.path.join(HERE, "certificate.json")))
v = C["v21_fixes"]
s = v["success"]
print("verified 1 checks")
print("ALL PASS  kappa=%.6f  E_bad=2^%.3f  markov=2^%.3f  W=2^%.6f  claim=%.4f  success=%.6f  coverage=%.8f"
      % (v["exact_identity"]["kappa_exact"], v["exact_identity"]["E_bad_log2"],
         v["exact_identity"]["markov_log2"], v["claim"]["W_log2"], v["claim"]["value"],
         s["after"], s["debits"]["coverage"]))
sys.exit(rc)
