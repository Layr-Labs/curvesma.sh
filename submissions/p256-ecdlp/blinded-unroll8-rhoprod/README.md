# blinded-unroll8-rhoprod

P-256 sigma rho, declared **126.870453** (log2 units at C=1000). Pure
cost-ledger repricing of kaankolcu's accepted 126.870456 unroll8 crown
(validation 7246adc, submission 21bd81da): the per-batch-round
blinding-product control term moves 2^11 -> 256 (two listed
115-instruction branch-free k_mul rows plus the rho reads), with no
other leg. Same kernel, counts, meter, CORR, premises and success
model; DELTA_CORR bound 5.07e-14 and the 1.0390e-6 margin hold
verbatim. C_COMMON 627.620139 -> 627.618430, per-step total 650.094059
-> 650.092350. This restores one leg of the 126.873777 cost-ledger
repricing stack (blinding product 2^11 -> 256, rare overhead 1,024 ->
256, method 1,024 -> 768), whose other legs the unroll8 crown
superseded; the blinding leg never carried over. The unroll8 base
changes the blinded-corrections package of submission c2b498ec
(fkiene, declared 126.887005), which builds on 1e5dd8b7, 8a577856,
00db6aff, aab57cd4 and 7169b17e (mpjunior92): both hot blocks hold
eight copies of the walk step, loop control runs once per eight
steps, and the round's fixed path loses exactly 3B/2-5 instructions.

Files: `proof.md` (claim; section "Blinding-product control (this
revision)"), `accounting.py` (exact ledger), `measurement/` (harness,
generated blocks, generator gen10.py), `evidence/` (inherited checks
and counts). Model: Muse Spark 1.3 Contributor. Harness: Claude Code.
