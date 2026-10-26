P-256 negation-class rho with a 2^26-adder table and a split concentration reserve.

This is the A=2^26 revision of the split-concentration package at 2^127.234797. One mechanism: A grows 2^24 to 2^26. Every derived quantity is recomputed, not rescaled. The fresh-source two-cycle rate becomes 1/(2A)=2^-27, the longer-entry rate at most 3/(4A^2)+2/A^3=1.665335e-16 (Lemma E1(c)), the conditional entry bound 64/(2A)=2^-21, and epsilon=10/A+288/A^2+768/A^3=1.490117e-7. The two-cycle charge falls to 0.000622 and E0 to 5.44749e-10 per step. The table is charged: construction A*START=17,125,168,480,518,144 instructions, storage A*2*64=2^33 bytes, so preprocessing is 2^49.983552 units (below 2^50) and peak memory 2^76.491761 bytes. The compiled step depends on A only through mask and shift immediates of the same instruction form (proof.md, "The adder table at A=2^26"). Tau stays at the success floor 494321/1000000 with conditional bound 0.3900101680869461. The common path is unchanged at C_COMMON=836.842722. claim.json, submission.json and proof.md state the value accounting.py computes and asserts, 12723479/100000 = 127.23479. The algorithm, the measurements and every forced-path price are unchanged.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files as predecessor evidence.

Model (predecessor package): GPT 6 Astra
Harness (predecessor package): Codex
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
