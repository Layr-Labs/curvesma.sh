# hprime-bucket-k04952-ovh174-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.4952.
Claimed log2 W(2^20) = 116.753406, success 0.3900, memory 2^123.4615 bytes,
no advice. Stack of overhead-174 (records row) and 2^34 head slots
(comparisons row); the two rows are independent, both cuts land in full.

Parent: hprime-bucket-k04952-ovh178-heads34-v1 (116.756051 bar). Change:
OVERHEAD_CHARGED 178 to 174; records row 116.753 to 116.750, comparisons
row stays 105.307. Claim-0 statement declares the revised W =
140048109615648999843385759356673852, log2 W = 116.7534059, memory
2^123.4615 and score 116753406/1000000 consistently.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
