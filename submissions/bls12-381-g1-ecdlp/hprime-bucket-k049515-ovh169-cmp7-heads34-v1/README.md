# hprime-bucket-k049515-ovh169-cmp7-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.49515.
Matching-path micro-trim: 7 instructions per key comparison instead of 8
(early exit on the high word, which mismatches first with probability
1 - 2^-10), on the reviewed K=0.49515 overhead-169 heads-34 crown. Records
row untouched; comparisons row drops; every cap, tail, success figure and
the preprocessing bound stays put.

Parent: hprime-bucket-k049515-ovh169-heads34-v1 (116.750021 bar). Change:
INSTR_PER_COMPARISON 8 to 7. Claim-0 statement declares the revised W,
log2 W and score consistently.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
