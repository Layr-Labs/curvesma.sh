# hprime-bucket-k04952-ovh176-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.4952.
Claimed log2 W(2^20) = 116.754729, success 0.3900, memory 2^123.4615 bytes,
no advice. Stack of overhead-176 (records row) and 2^34 head slots
(comparisons row); the two rows are independent, both cuts land in full.

Parent: hprime-bucket-k04952-ovh178-heads34-v1 (116.756051 bar). Change:
OVERHEAD_CHARGED 178 to 176; records row 116.753 to 116.751, comparisons
row stays 105.307. Claim-0 statement declares the revised W =
140176596557446780765046853246809801, log2 W = 116.7547289, memory
2^123.4615 and score 116754729/1000000 consistently.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
