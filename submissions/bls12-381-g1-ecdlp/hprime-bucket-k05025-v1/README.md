# hprime-bucket-k05025-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole score lever: the coverage parameter drops from K = 0.52 to K = 0.5025,
shrinking both lists by 1.7% (N scales as sqrt(K)). Predicted success moves
from 0.4055 to 0.3950, claimed at 0.392 (0.0030 H1-error budget, 0.002 over
the required 0.39). K moves below that package's stated 0.51 search boundary,
which was its authors' stopping choice, not a proved barrier: H1 is
K-independent and the budget structure is unchanged. Separately, the
preprocessing figure is upward-rounded (ceil to 4 decimals) and
integer-asserted, answering the 934a4565 rejection, which faulted
round-to-nearest preprocessing accounting. Algorithm, walks, routing,
matching, keys, caps, charges, and the toy procedure are the accepted package.

Claimed log2 W(2^20) = 116.768426, success 0.392, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
