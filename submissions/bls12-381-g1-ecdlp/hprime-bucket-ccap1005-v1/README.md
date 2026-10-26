# hprime-bucket-ccap1005-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole score lever: the per-node comparison cap drops from 1.01x to 1.005x the
per-node mean. The Bernstein tail moves from 2^-10732 to below 2^-1465, still
astronomically inside the 2^-40 tail budget, and the comparisons row shrinks
0.5%. Separately, the preprocessing figure is now upward-rounded (ceil to 4
decimals) and integer-asserted, answering the 934a4565 rejection, which faulted
round-to-nearest preprocessing accounting. Algorithm, walks, routing, matching,
keys, all other caps, and the toy procedure are the accepted package.

Claimed log2 W(2^20) = 116.793075, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
