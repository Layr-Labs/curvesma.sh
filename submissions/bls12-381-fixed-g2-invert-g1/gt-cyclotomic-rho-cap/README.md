# Pollard rho with cyclotomic squaring in GT

Given a nonidentity Q in G2 and z in GT, recover P in G1 with e(P, Q) = z
under the pinned pairing. The discrete logarithm is taken in GT by Pollard's
three-way distinguished-point walk. The squaring branch uses the
Granger–Scott cyclotomic formula, which applies because every walk state lies
in GT and the subgroup order divides p^4 - p^2 + 1. The other two branches are
ordinary dense multiplications by the base and by the target. The walk state
is the GT element itself: there is no Frobenius class and no representative
to select.

The work bound is a hard cap, not an expectation. The walk stops at
floor(1.2 sqrt(r)) steps or at ceil(11/16) of that many dense multiplications.
Each cyclotomic squaring is charged at 16384 scalar instructions plus 64
instructions of bookkeeping, against the reference price C = 20000. Each dense
branch is one reference operation plus the same bookkeeping. The arithmetic
is in [proof.md](proof.md). The declared bound is in [claim.json](claim.json).

experiment.py runs the same three branches on a small additive group.
Doubling there is only the image of the squaring branch. It does not measure
Fp12 costs.
