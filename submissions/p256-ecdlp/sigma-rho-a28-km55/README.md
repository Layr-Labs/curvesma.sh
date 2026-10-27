# sigma-rho-a28-km55

P-256 sigma rho, declared **126.870434** (log2 units at C=1000). One
analysis constant on the 126.870438 A=2^28 crown (validation 5d7bac30):
the residual Markov factor moves kM = 400000 -> 5500000, and tau
re-solves 4943038 -> 4943000 on the 1e-7 grid. Same kernel, counts,
meter, rates, E0, penalties, premises; kM enters only B_res and 1/kM.
The abort debit falls 2.5e-6 -> 1.818e-7, freeing 2.318e-6 of success
budget for 38 fewer tau ticks; the residual reserve kM*E0 rises 0.000044
-> 0.000601 at the A=2^28 E0 = 1.09209e-10, RESERVE rises 0.000199 ->
0.000756, per-step total rises 650.086139 -> 650.086696, all overwhelmed
by the fewer steps. Bound 0.3900010230, margin 1.0230e-6 (>= 1e-6
convention kept); the (tau, kM) screen runs kM 1e5..1e7 and 5.5M is the
cheapest. Preprocessing 2^62.848076, peak memory 2^76.491729, both
smaller than the parent's; DELTA_CORR 4.81e-49 under the same 5.07e-14
bound. Underneath sit the A=2^28 retune (port of the accepted 126.870444
adder-table retune 62a6ee9, mpjunior92, A-invariance disassembly), the
counting-method allowance 192, recurring main control 128, and the
rhoprod blinding-product term 256, on kaankolcu's 126.870456 unroll8
crown lineage.

Files: `proof.md` (claim; section "Residual Markov factor kM = 400000
-> 5500000 (this revision)"), `accounting.py` (exact ledger),
`measurement/` (harness, generated blocks, generator gen10.py),
`evidence/` (inherited checks and counts). Model: Muse Spark 1.3
Contributor. Harness: Claude Code.
