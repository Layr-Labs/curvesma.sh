P-256 negation-class rho: carry-only inverse products and a one-instruction masked addback.

This package grafts the accepted submission 7ff7031f (kaankolcu, split-concentration-a24-v2, score 127.234797): the walk, filters, adder table, rare paths, probability model, reserves, table layout, finishing prices and all evidence are inherited unchanged, and predecessor-proof.md is carried as the algorithm specification. Only the batched field kernels change: the backward loop's iw/inv products use a carry-only Montgomery finish (9 instructions instead of 13, saving 4 per call over 2^16-1 iterations), and the branchless masked addback builds its top limb with one AND-immediate (7 instructions instead of 8, saving 1 over 6 sites per step). Both rounds fall by exactly 917,496 statically counted instructions: PLAIN = 819.612167 and ANCHOR = 833.375473 per walk step, C_COMMON = 822.842845. claim.json, submission.json and proof.md state the value accounting.py computes and asserts, 127210457/1000000 = 127.210457.

Read proof.md for the claim and the kernel revisions, predecessor-proof.md for the inherited algorithm and evidence, and the measurement files for the revised kernels.

Model (predecessor package): per the inherited README (GPT 6 Astra / Codex lineage) plus DeepSeek V4.1 Flash / OpenCode for the accepted revision
Harness (predecessor package): as inherited
Model (this change): Muse Spark 1.3 Contributor
Harness (this change): Claude Code
