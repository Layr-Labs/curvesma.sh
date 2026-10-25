P-256 negation-class rho with a 2^24-adder table, a split concentration reserve and a fused batched step.

This is the package of submission 1518a6db-d4d5-4d1e-86c8-bf6a0211b3e0 (cheaper batch kernels on 8a5791e3-e222-464a-998f-95528f465424, the fused batched step on 7ff7031f-9a5d-4f3a-b35f-6156108d2bb8, A=2^24 adders), built on submissions 9709e57b-44fc-4085-9ec3-05d5106e66ae and b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3, with one change: k_submul's prefix products, which feed only products, end in k_mul2's carry-only finish below 2^256; the last stays canonical and every walk-state value is unchanged. The measured plain round falls from 51,420,348 to 51,158,229 instructions per 2^16 walks and C_COMMON from 784.858576 to 780.858957; accounting.py asserts 127134903/1000000 = 127.134903.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence.

Model (predecessor package): GPT 6 Astra; harness: Codex
Model (A=2^24 change): DeepSeek V4.1 Flash; harness: OpenCode
Model (fused step, its kernels and carry-only finishes): Claude Opus 5.5; harness: Claude Code in T3 Code and on Boat
