"""Evaluator source is model16core.py inside lib.b126.

Local edits (this submission), each of which the package recomputes and prints
on every run: no measurement, price or region parameter is changed.

1. The [5c]/[5d] row-cache debit.  model16core.row_cache charges
   d = (155.4-64)/RLEN_STRONG + (64-17)/RLEN_STRONG - (64-155.4/16) X_SHORT per
   cell, with RLEN_STRONG = 600 and X_SHORT = 1.017e-4 as round numbers above
   the certified values.  [5c] recomputes both exactly on 256 radius intervals
   (lattice_count_bounds.strong_rows, exact rationals) and prints
   N6/Nstrong <= 584.025935 and (N5-Nstrong)/N6 <= 9.9274e-05, so the two round
   numbers can be replaced by 584.1 and 9.93e-05, both still above the
   certified values read by the same check.  Nothing else in [5c]/[5d] changes.
2. The block-Wiedemann sizes are charged at the base's own linear-algebra model
   optimum (768, 40) instead of (1024, 40): la_bn.linear_algebra is the base's
   cost model and returns a lower cost at (768, 40) for both the linear algebra
   and the generator row it drives.  H-X remains a declared premise at the
   charged lengths (RANK-BRIDGE.md states this).
"""
import common12
from model16core import *
from model16core import _ctx
import model16core as _m16
import la_bn as _la

# 1. [5c]/[5d]: the certified extremes instead of the base's round numbers
_m16.RLEN_STRONG = 584.1
_m16.X_SHORT = 9.93e-05
RLEN_STRONG = _m16.RLEN_STRONG
X_SHORT = _m16.X_SHORT

# 2. block Wiedemann sizes at the LA model's own optimum
_la.BW_M, _la.BW_N = 768, 40
