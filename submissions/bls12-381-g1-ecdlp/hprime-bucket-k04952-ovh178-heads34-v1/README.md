# hprime-bucket-k04952-ovh178-heads34-v1

Cheon d=2^20 with H-prime semiprime walks and bucket matching at K=0.4952.
Claimed log2 W(2^20) = 116.756051, success 0.3900, memory 2^123.4615 bytes,
no advice. Stack of overhead-178 (records row) and 2^34 head slots
(comparisons row); the two rows are independent, both cuts land in full.

Parent: hprime-bucket-k04952-ovh178-v1 (116.756567 bar). Change: per-node
head array 2^33 to 2^34 slots (LOG_HEADS 33 to 34, MEMORY_CLAIM 1233242
to 1234615); comparison row 106.307 to 105.307. Claim-0 statement declares
the revised W = 140305083499244561686707947136945751, log2 W = 116.7560507,
memory 2^123.4615 and score 116756051/1000000 consistently (cf. dba1f42b).

Files: submission.json, claim.json, evidence.py, proof.md.
Recompute: echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
