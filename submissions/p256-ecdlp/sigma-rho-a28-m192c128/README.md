# sigma-rho-a28-m192c128

P-256 sigma rho, declared **126.870438** (log2 units at C=1000). One
parameter change on the 126.870453 rhoprod crown (validation 84124941):
the adder table grows A = 2^24 -> 2^28, and with it the counting-method
allowance moves 1,024 -> 192 per 2^20-walk round (6x the largest observed
28-instruction gap, 28 offset kept) and recurring main control moves 256
-> 128 (itemised 105 bound, 23 headroom); finishing rounds keep 1024 +
RHO_COST. Same kernel, counts, meter, CORR, rho, premises and success
model; DELTA_CORR bound 5.07e-14 holds verbatim. Every A-dependent rate
and reserve is recomputed exactly: two-cycle charge 0.002489 -> 0.000156,
E0 7.51339e-9 -> 1.09209e-10, RESERVE 0.005494 -> 0.000199; C_COMMON
unchanged at 627.617514 (A-invariant common path), per-step total
650.091434 -> 650.086139; tau re-solves to 4943038/10000000 with a
1.0228e-6 margin. This ports the accepted 126.870444 adder-table retune
(62a6ee9, mpjunior92, A-invariance disassembly) onto our rhoprod +
method/control stack, and answers the 5d2b461a rejection's open passage
with a uniform 192 allowance declaration. The rhoprod base repriced
kaankolcu's 126.870456 unroll8 crown (validation 7246adc), which changes
the blinded-corrections package of submission c2b498ec (fkiene, declared
126.887005), itself on 1e5dd8b7, 8a577856, 00db6aff, aab57cd4 and 7169b17e
(mpjunior92): both hot blocks hold eight copies of the walk step, loop
control runs once per eight steps, and the round's fixed path loses
exactly 3B/2-5 instructions.

Files: `proof.md` (claim; sections "Adder table A = 2^24 -> 2^28 (this
revision)", "Counting-method allowance (this revision)" and "Recurring
main control (this revision)"), `accounting.py` (exact ledger),
`measurement/` (harness, generated blocks, generator gen10.py),
`evidence/` (inherited checks and counts). Model: Muse Spark 1.3
Contributor. Harness: Claude Code.
