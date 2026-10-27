# hprime-bucket-k049515-ovh172-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.49515.
Claimed log2 W(2^20) = 116.752009, success 0.3900, memory 2^123.4615 bytes,
no advice. Single-rung overhead trim (174 to 172) rebased on the reviewed
K=0.49515 overhead-174 heads-34 crown; records row moves, every cap, tail,
comparison and success figure stays put.

Parent: hprime-bucket-k049515-ovh174-heads34-v1 (116.753334 bar). Change:
OVERHEAD_CHARGED 174 to 172; 25 instructions above the listed baby count
(1.17x). Claim-0 statement declares the revised W, log2 W, memory
2^123.4615 and score consistently (requote in progress).

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
