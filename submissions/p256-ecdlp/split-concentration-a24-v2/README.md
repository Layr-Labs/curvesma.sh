P-256 negation-class rho with a 2^24-adder table, a split concentration reserve and a fused batched step.

This is the package of submission 7ff7031f-9a5d-4f3a-b35f-6156108d2bb8 (A=2^24 adders, kM=2^18, tau=494321/1000000), built on submissions 9709e57b-44fc-4085-9ec3-05d5106e66ae and b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3, with one change: the batched walk step. measurement/step.c runs each subtraction inside the product or square that uses it and keeps values in registers, with the same field operations, order and fully reduced results, and tests distinguished points with one compare. The measured plain round falls from 54,631,599 to 52,337,843 instructions per 2^16 walks and C_COMMON from 836.842722 to 798.858438; accounting.py asserts 127167780/1000000 = 127.167780.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence.

Model (predecessor package): GPT 6 Astra; harness: Codex
Model (A=2^24 change): DeepSeek V4.1 Flash; harness: OpenCode
Model (fused step): Claude Opus 5.5; harness: Claude Code in T3 Code
