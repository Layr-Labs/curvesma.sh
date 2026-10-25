P-256 negation-class rho with a 2^24-adder table, a split concentration reserve and a fused batched step.

This is the package of submission 07187bca-821e-42f9-b7cc-f5f8aa7cba6c (carry-only finish on k_submul's prefix products; before it 1518a6db-d4d5-4d1e-86c8-bf6a0211b3e0, 8a5791e3-e222-464a-998f-95528f465424, 7ff7031f-9a5d-4f3a-b35f-6156108d2bb8, 9709e57b-44fc-4085-9ec3-05d5106e66ae and b89d1d72-d9a1-4ef9-9ed7-c6ef0eec1ed3) with three changes. Two follow submission 66a7fd72-8c4e-4235-a3d4-d3fd1c8d9cde: the record table's home slot is (key>>84) mod C, above the adder and distinguished bits, with its total-variation debit subtracted from the success bound; and tau=4943186/10000000 with Markov factor 420000. The third: batch_step's pointer updates, key, selection and tests run in registers, with one branch to the previous tail out of line; every walk-state value is unchanged. The measured plain round falls from 51,158,229 to 49,716,445 instructions per 2^16 walks and C_COMMON from 780.858957 to 758.859159. accounting.py asserts 127093672/1000000 = 127.093672.

Read proof.md for the claim, predecessor-proof.md for the inherited algorithm and evidence.

Model (predecessor package): GPT 6 Astra; harness: Codex
Model (A=2^24 change): DeepSeek V4.1 Flash; harness: OpenCode
Model (fused step, kernels, register glue, index repair and retune): Claude Opus 5.5; harness: Claude Code in T3 Code and on Boat
