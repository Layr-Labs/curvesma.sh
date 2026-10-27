# blinded-unroll8-rhoprod-method512

P-256 sigma rho, declared **126.870451** (log2 units at C=1000). Pure
cost-ledger repricing of our 126.870453 rhoprod crown (validation
84124941, submission e49858b7 lineage): with that crown's
blinding-product term 2^11 -> 256 kept, the counting-method allowance
moves 1,024 -> 512 per 2^20-walk round (still 18x the largest observed
28-instruction gap), with no other leg. Same kernel, counts, meter,
CORR, rho, premises and success model; DELTA_CORR bound 5.07e-14 and
the 1.0390e-6 margin hold verbatim. C_COMMON 627.618430 -> 627.617941,
per-step total 650.092350 -> 650.091862. This stacks the
counting-method leg of the 126.873777 cost-ledger repricing stack
(blinding product 2^11 -> 256, rare overhead 1,024 -> 256, method
1,024 -> 768 at 28x the 27 gap) onto the rhoprod crown, which is that
stack's blinding-product leg repriced onto kaankolcu's 126.870456
unroll8 crown (validation 7246adc). The unroll8 base changes the
blinded-corrections package of submission c2b498ec (fkiene, declared
126.887005), which builds on 1e5dd8b7, 8a577856, 00db6aff, aab57cd4
and 7169b17e (mpjunior92): both hot blocks hold eight copies of the
walk step, loop control runs once per eight steps, and the round's
fixed path loses exactly 3B/2-5 instructions.

Files: `proof.md` (claim; section "Counting-method allowance (this
revision)"), `accounting.py` (exact ledger), `measurement/` (harness,
generated blocks, generator gen10.py), `evidence/` (inherited checks
and counts). Model: Muse Spark 1.3 Contributor. Harness: Claude Code.
