# hprime-bucket-ovh162b-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Score lever: the charged per-record overhead drops from 180 to 162
instructions against the unchanged 143/139 itemized listing (19-instruction,
1.13x margin; byte transport stays itemised). The per-record factor moves
1.090 to 1.081. The preprocessing figure is upward-rounded (ceil to 4
decimals) and integer-asserted, answering the 934a4565 rejection, whose sole
point faulted round-to-nearest preprocessing accounting while leaving the
total-work score unrefuted. Algorithm, walks, routing, matching, keys, caps,
and the toy procedure are the accepted package.

Claimed log2 W(2^20) = 116.781155, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
