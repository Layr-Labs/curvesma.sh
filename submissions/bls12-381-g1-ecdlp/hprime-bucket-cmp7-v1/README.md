# hprime-bucket-cmp7-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole score lever: the per-comparison charge drops from 8 to 7 instructions.
The pointer advance folds into the key load's post-indexed addressing
(`ldp xA, xB, [ptr], #16`), so the comparison itemizes as post-indexed ldp,
ldr, cmp/ccmp chain, conditional branch, and counter decrement: every op of
the accepted listing stays covered with one fewer instruction. The
comparisons row shrinks 12.5%. Separately, the preprocessing figure is
upward-rounded (ceil to 4 decimals) and integer-asserted, answering the
934a4565 rejection, which faulted round-to-nearest preprocessing accounting.
Algorithm, walks, routing, matching, keys, caps, and the toy procedure are
the accepted package.

Claimed log2 W(2^20) = 116.792948, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
