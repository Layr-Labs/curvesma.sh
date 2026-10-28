"""y3 fold: derive this revision's bloop_branchy.h from the paired-reduction revision's.

usage: gen13.py <paired-reduction measurement dir> <output dir>

A compare of the top limbs picks the order of y3 = y3' - y(W): y(W) - y3' = -y3 (no borrow) when y3' has the smaller
top limb, else the old order and its fix-up (site y3_tie). The edit must match exactly once (proof.md, "The y3 fold").
"""
import sys

LD = r'ldp x8,x9,[x24],#-64\n'
OLD = LD + r'subs x0,x0,x8\nsbcs x2,x2,x9\nsbcs x18,x18,x10\nsbcs x30,x30,x11\nb.cs Lq10" S "_%='
NEW = LD + (r'cmp x30,x11\nb.hs Lp10" S "_%=\nsubs x0,x8,x0\nsbcs x2,x9,x2\nsbcs x18,x10,x18\nsbc x30,x11,x30\n'
            r'b Lq10" S "_%=\nLp10" S "_%=:\n') + OLD[len(LD):]
t = open(sys.argv[1] + "/bloop_branchy.h").read()
assert t.count(OLD) == 1 and "Lp10" not in t
assert t.count("(paired reduction words)") == 1
t = t.replace(OLD, NEW).replace("(paired reduction words)", "(paired reduction words), then gen13.py (y3 fold)", 1)
open(sys.argv[2] + "/bloop_branchy.h", "w").write(t)
print("ok")
