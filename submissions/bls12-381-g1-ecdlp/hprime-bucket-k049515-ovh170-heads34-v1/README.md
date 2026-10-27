# hprime-bucket-k049515-ovh170-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.49515.
Claimed log2 W(2^20) = 116.750684, success 0.3900, memory 2^123.4615 bytes,
no advice. Overhead trim (172 to 170) on the reviewed K=0.49515 heads-34
crown; records row moves, every cap, tail, comparison and success figure
stays put.

Parent: hprime-bucket-k049515-ovh172-heads34-v1 (116.752009 bar). Change:
OVERHEAD_CHARGED 172 to 170; 23 instructions above the listed baby count
(1.16x). Claim-0 statement declares the revised W, log2 W, memory
2^123.4615 and score consistently.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
