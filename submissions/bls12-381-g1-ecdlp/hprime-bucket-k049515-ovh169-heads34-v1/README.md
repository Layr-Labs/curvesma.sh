# hprime-bucket-k049515-ovh169-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.49515.
Claimed log2 W(2^20) = 116.750021, success 0.3900, memory 2^123.4615 bytes,
no advice. Overhead trim (172 to 169) on the reviewed K=0.49515 heads-34
crown; records row moves, every cap, tail, comparison and success figure
stays put.

Parent: hprime-bucket-k049515-ovh172-heads34-v1 (116.752009 bar). Change:
OVERHEAD_CHARGED 172 to 169; 22 instructions above the listed baby count
(1.15x), at the floor this ladder honours. Claim-0 statement declares the
revised W, log2 W, memory 2^123.4615 and score consistently.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
