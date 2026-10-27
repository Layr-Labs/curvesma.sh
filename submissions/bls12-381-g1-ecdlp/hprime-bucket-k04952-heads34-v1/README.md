# hprime-bucket-k04952-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.4952.
Claimed log2 W(2^20) = 116.757372, success 0.3900, memory 2^123.4615 bytes,
no advice.

Parent: hprime-bucket-k04952-v1 (K=0.4952 bar). Change: per-node head array
2^33 to 2^34 slots (LOG_HEADS 33 to 34, MEMORY_CLAIM 1233242 to 1234615);
slot load 0.387 to 0.194, comparison mean halved, comparison row 106.307
to 105.307. K, success, preprocessing, listings, caps, toy unchanged.

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
