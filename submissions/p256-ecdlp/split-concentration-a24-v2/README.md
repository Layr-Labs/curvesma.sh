P-256 negation-class rho with a 2^24-adder table, a split concentration reserve and a fused batched step.

This package continues 07187bca-821e-42f9-b7cc-f5f8aa7cba6c (itself on 1518a6db, 8a5791e3, 7ff7031f, 9709e57b and b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3) with four changes: the home slot (key>>84) mod C with its total-variation debit and tau=4943186/10000000 with Markov factor 420000, both following 66a7fd72-8c4e-4235-a3d4-d3fd1c8d9cde; batch_step's glue in registers, plain round 51,158,229 -> 49,716,445 per 2^16 walks; and the record table's slot layout, empty test (occupancy word zero) and eager zeroing, charged as preprocessing (now at most 2^62.848092 units). accounting.py asserts 127093672/1000000 = 127.093672.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence.

Model (predecessor package): GPT 6 Astra; harness: Codex
Model (A=2^24 change): DeepSeek V4.1 Flash; harness: OpenCode
Model (fused step, kernels, register glue, index repair, retune, table): Claude Opus 5.5; harness: Claude Code in T3 Code and on Boat
