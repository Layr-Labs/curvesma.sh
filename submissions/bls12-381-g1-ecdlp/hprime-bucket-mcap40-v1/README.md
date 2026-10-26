# hprime-bucket-mcap40-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole change: the global candidate-pair cap drops from 2^44 to 2^40. The
false-match mean is 2^39.33, so the Chebyshev cap-overflow tail moves from
below 2^-47.5 to 2^-36.8, still inside an explicit 2^-36 tail budget, and
every cap-proportional recovery row (candidate forwarding, both-chain
re-walks with table rebuilds, kept-value forwarding, per-candidate scalar
multiplications) shrinks 16x. Algorithm, walks, routing, matching, keys,
caps otherwise, and the toy procedure are byte-for-byte the accepted package.

Claimed log2 W(2^20) = 116.793079, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
