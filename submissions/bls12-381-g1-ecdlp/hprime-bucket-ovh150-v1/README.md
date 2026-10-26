# hprime-bucket-ovh150-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Score lever: the charged per-record overhead drops from 180 to 150
instructions against the unchanged 143/139 itemized listing (7-instruction
margin over the baby side; byte transport stays itemised). The per-record
factor moves 1.090 to 1.075. Same proof shape as the reviewed 162/158/154
steps, one notch tighter; every listed item already rounds its own
amortisation up, so the remaining allowance still covers only unlisted slack.
The preprocessing figure is upward-rounded (ceil to 4 decimals) and
integer-asserted, answering the 934a4565 rejection, which faulted
round-to-nearest preprocessing accounting. Algorithm, walks, routing,
matching, keys, caps, charges otherwise, listings, and the toy procedure are
the accepted package.
Claimed log2 W(2^20) = 116.773150, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
