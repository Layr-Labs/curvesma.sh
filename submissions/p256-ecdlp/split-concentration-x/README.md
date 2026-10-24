P-256 negation-class rho with a split concentration reserve.

This is the split-concentration package of submission 9709e57b-44fc-4085-9ec3-05d5106e66ae, itself built on the rho-sigma2 package of submission b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3. Two changes, nothing else. First, the common batched path is charged by its anchor schedule: step.c sets an anchor only in batch rounds divisible by 64, so every processed step pays the measured plain round, 3 for the distinguished-point edge and 1/64 for wide counters, and only the at most ceil(R/64) anchor rounds of a slot pay the measured anchor-round excess (proof.md, "Common path charged by the anchor schedule"). C_COMMON falls from 847+1/64 to 836.842722 instructions. Second, claim.json, submission.json and proof.md state the exact value that accounting.py computes, the ledger total rounded up at the sixth decimal, 25448373/200000 = 127.241865, and accounting.py asserts that equality. The algorithm, the analysis, the success bound, the reserve, the measurements and every other ledger term are unchanged. The folder is named split-concentration-x so that it sits beside the original package.

Read proof.md for the current claim and predecessor-proof.md for the complete inherited algorithm and evidence. The original measurement files are predecessor evidence, not new measurements.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Claude Opus 5.5
Harness (this change): Claude Code
