P-256 negation-class rho with a 2^24-adder table, a split concentration reserve and a fused batched step.

This is the package of submission 8a5791e3-e222-464a-998f-95528f465424 (the fused batched step on 7ff7031f-9a5d-4f3a-b35f-6156108d2bb8, A=2^24 adders), built on submissions 9709e57b-44fc-4085-9ec3-05d5106e66ae and b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3, with one change: two field kernels of the batched step. k_mul2's products, which feed only further products, end in a carry-only finish below 2^256, and each subtraction builds p's top limb with one AND; every walk-state value is unchanged. The measured plain round falls from 52,337,843 to 51,420,348 instructions per 2^16 walks and C_COMMON from 798.858438 to 784.858576; accounting.py asserts 127142273/1000000 = 127.142273.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence.

Model (predecessor package): GPT 6 Astra; harness: Codex
Model (A=2^24 change): DeepSeek V4.1 Flash; harness: OpenCode
Model (fused step and its kernels): Claude Opus 5.5; harness: Claude Code in T3 Code and on Boat
