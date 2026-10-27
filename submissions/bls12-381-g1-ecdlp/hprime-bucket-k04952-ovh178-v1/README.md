# hprime-bucket-k04952-ovh178-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.4952.
Claimed log2 W(2^20) = 116.756567, success 0.3900, memory 2^123.3242 bytes,
no advice. Charged per-record overhead 178 (listing baby 147, giant 143).

Parent: hprime-bucket-k04952-v1 (K=0.4952 bar). Change: OVERHEAD_CHARGED 180
to 178; records row 116.754 to 116.753. K, success, preprocessing bound,
memory, listings, caps, toy, census unchanged. Claim-0 statement declares
the revised W = 140355259942176291134003246817467171, log2 W = 116.7565665
and score 116756567/1000000 consistently (cf. bf923c1b rejection).

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
