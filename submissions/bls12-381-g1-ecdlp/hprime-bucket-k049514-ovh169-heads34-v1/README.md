# hprime-bucket-k049514-ovh169-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.49514.
Claimed log2 W(2^20) = 116.750007, success 0.3900, memory 2^123.4615 bytes,
no advice. Coverage sub-notch (49515 to 49514) on the reviewed K=0.49515
overhead-169 heads-34 crown; N shrinks with sqrt(K), records row moves,
every cap, tail, comparison and overhead figure stays put.

Parent: hprime-bucket-k049515-ovh169-heads34-v1 (116.750021 bar). Change:
K_NUM 49515 to 49514; success after H1 allowance 0.390013, margin 1.3e-05,
still positive and asserted. Claim-0 statement declares the revised K,
success, W, log2 W and score consistently.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
