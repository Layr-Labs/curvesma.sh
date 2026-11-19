"""Evaluator source is model16core.py inside lib.b126.

Local edits (this submission), each of which the package recomputes and prints
on every run: no measurement, kernel price, region parameter or heuristic is
changed.

1. The [5c]/[5d] row-cache debit is charged as the base's distribution-free
   per-slice worst case, with the two round caps the base charged above the
   values its own [5c] check reads: model16core.RLEN_STRONG = 600 and
   model16core.X_SHORT = 1.017e-4 become 584.1 and 9.93e-5.  [5c] prints
   N6/Nstrong <= 584.025935 and (N5 - Nstrong)/N6 <= 9.9274e-05 on every radius
   interval, so the charge
       d = 138.4/584.1 - 54.2875 x 9.93e-5 = 0.231555
   is the base's formula with the certified extremes.
2. The block-Wiedemann sizes are charged at the base's own linear-algebra model
   optimum: la_bn.BW_M, la_bn.BW_N become (768, 40) from (1024, 40); the ledger's
   linear algebra and generator rows move to the lower cost the model returns
   there.  H-X remains a declared premise at the charged lengths.
"""
import common12
from model16core import *
from model16core import _ctx
import model16core as _m16
import la_bn as _la

# 1. [5c]/[5d] caps: the certified extremes instead of the base's round numbers.
#    584.1 >= max over the 256 intervals of N6/Nstrong (printed 584.025935);
#    9.93e-5 >= max over the intervals of (N5 - Nstrong)/N6 (printed 9.9274e-05).
_m16.RLEN_STRONG = 584.1
_m16.X_SHORT = 9.93e-05
RLEN_STRONG = _m16.RLEN_STRONG
X_SHORT = _m16.X_SHORT

# 2. block-Wiedemann sizes at the shipped linear-algebra model's own optimum
_la.BW_M, _la.BW_N = 768, 40

# Reference points for the Price paragraph of proof.md (all recomputed in the
# package, none asserted):
#   flat (charged here):    138.4/584.1 - 54.2875 x 9.93e-5 = 0.2315549797
#   a-density alternative:  138.4 x 0.0019375846            = 0.2681617091
#   (the a-density is the reviewer-style weight of the certified ratios on the
#   six-dimensional a-population, density (lam - s)^2 s^3 per unit radius s;
#   the charged flat value is 13.2% conservative under it)
