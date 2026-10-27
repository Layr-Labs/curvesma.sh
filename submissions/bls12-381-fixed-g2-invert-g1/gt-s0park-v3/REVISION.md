# S0 constant parking in the squaring kernel; re-rounded work cap 125.1169

This package reuses the complete Frobenius-class rho algorithm, sampler,
matching, caps, meter, empirical evidence and accounting of
`gt-rho-finecap-v2` (claim 125.11795). This is a new standalone claim, not a
prior-results acceptance reference.

The change is one measured kernel speedup plus its exact ledger
consequences. In `kernels/gen.py` (`sqr_kernel`, stored-sum loop S0) the six
single-limb const-keys each reloaded the 6-limb pool constant p (3 ldp) and
then chained one vector add and one vector subtract (6/6 instructions). The
Z-bank registers x8-x13 are untouched anywhere in S0, so this revision parks
p there once (3 ldp, pinned before the loop) and rebuilds each of the six
keys as A = V+ - V- + p: 3 ldp for V+, 6 chain for the subtract, 6 chain for
p from registers. Each key drops from 21 to 18 instructions; the parking
loads cost 3. Net: squaring kernel 6,767 -> **6,752** (-15), squaring step
6,993 -> **6,978** = 0.3489 units, average step 0.55309921875 ->
**0.552665625** units. `kernels/check_kernels.py` asserts the new emitted
count (mul unchanged at 14,398) and `check.py` re-verifies every output limb
(0 fails); the S0-only assembly diff is the park plus six 21->18
rebuilds, with all later phases byte-identical.

`experiment.py` re-derives every cap-dependent term at the new prices and
re-rounds the two-percentage-point drift-covering ledger upward on the
1/20,000-bit grid: nominal ledger 2^125.09186, drift ledger 2^125.1168526,
least grid point **2502338/20000 = 125.1169**, checked by the same exact
integer inequality (at most, and least). Drift costs 0.87395% per point;
the cap exceeds the nominal ledger by 1.75124% (0.02505 bits); walk-budget
slack is 0.009679 units per step over 2.0038 tolerated points; success after
every loss stays above 0.39 (model 0.4005, H1+ 0.01). The meter, caps,
matcher, sampler and heuristics are unchanged.

Material limitations are inherited unchanged: random-mapping, branch/cap
concentration and ideal-source heuristics; very small success margin;
toy/short production-field evidence rather than a production attack.
Historical acceptance does not certify this package. The manifest declares
the complete claim for fresh evaluation. All earlier kernel innovations and
recorded experiments are attributed to the source lineage, not to this
revision.
