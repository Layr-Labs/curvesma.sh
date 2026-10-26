# Frobenius-class rho with an inline-add tariff and an honest cap analysis

Claim: time_log2 = **125.68** (total work <= 2^125.68 dense-Fp12-multiplication equivalents,
C = 20,000), success probability >= 0.39, peak memory < 2^73 bytes, preprocessing <= 2^31 units
(included), no advice.

- [Algorithm, success argument and cost ledger](proof.md)
- [Executable checks and toy runs](experiment.py)

## Lineage of this revision

The machine, algebra, walk, branch, escapes, caps, heuristics and evidence are those of the live
crown package `gt-tracecounted-exact-v25` (GordoAR, 125.72), a repricing of `747ad6c9` (125.74),
itself a revision of `6d4e57d1` (125.80) and of the Frobenius-class rho of `7da5f23a` (mpjunior92).
That package charges its key listing at 45 instructions per add/sub — the price of a called
routine. This revision inlines the nine key adds/subs at 36 instructions each (same data path,
no call overhead) and holds back a 45-instruction register-pressure margin: 2,119 -> **2,083**
on a multiply step, the other two charges unchanged (270 / 11,369). One tariff prices the walk
and the meter through one pricing function. Dependent figures: multiply step 1.11945 ->
1.11765 units, c 0.9178875 -> 0.9167625, W 2^125.7074975 -> 2^125.7057286, claim 125.72 ->
**125.68** — keeping a 1.80% allowance over the ledger that absorbs the measured scatter
instead of spending the cushion. Section 3.4 analyses the metered abort honestly: the
short-walk success cost and the cap-event probability are priced as two separate events under
the joint 0.01 H1+ budget (no exact-3/8 assertion), and part C now runs the production rule at
the declared R = 1024 and prices the measured fraction 0.368 in the drift table.
`experiment.py` part D reproduces the parent's ledger and this one through the same code path,
so the price vector is the only difference between the two ledgers.

This is an analytic generic-attack bound, not a practical break.
