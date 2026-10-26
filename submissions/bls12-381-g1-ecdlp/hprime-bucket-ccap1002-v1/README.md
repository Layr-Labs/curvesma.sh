# hprime-bucket-ccap1002-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole score lever: the per-node comparison cap drops from 1.01x to 1.002x the
(Chernoff tail moves 2^-1475 to 2^-305, still astronomically inside the
2^-40 budget). Conditional mean shift is then under 0.098%, leaving the
0.1% Bernstein deviation inside the 0.2% margin (0.098% + 0.1% < 0.2%); the
Bernstein tail stays below 2^-89 and the comparisons row shrinks 0.8%.
Same proof shape as the reviewed 1.003x step, one notch tighter. Separately,
the preprocessing figure is upward-rounded (ceil to 4 decimals) and
integer-asserted, answering the 934a4565 rejection, which faulted
round-to-nearest preprocessing accounting. Algorithm, walks, routing,
matching, keys, all other caps, charges, listings, and the toy procedure are
Claimed log2 W(2^20) = 116.793072, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
