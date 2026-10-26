# hprime-bucket-heads34-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole score lever: head slots per owner node double from 2^33 to 2^34, so total
head slots go 2^117 to 2^118, slot load halves (0.396 to 0.198), and bucket
comparisons halve (row 106.377 to 105.377). Keys stay 192 bits (84 node + 34
slot + 74 stored), so false matches are unchanged; the owner node grows to
2^39.46 B, still inside its 2^40 B budget (asserted). Head zeroing is counted
explicitly for the doubled array (2^33 stp + 6 x 2^29 loop) with setup/flushes
at the accepted 3 x 2^33; no hot path is re-itemized and no new instruction
type is introduced. 2^35 slots would quarter comparisons but the doubled
zeroing cost exceeds the gain, so 2^34 is the chosen point. Separately, the
preprocessing and memory figures are upward-rounded and integer-asserted,
accounting. Algorithm, walks, routing, matching, keys, caps, charges, and the
toy procedure are the accepted package.

Claimed log2 W(2^20) = 116.792864, success 0.40, memory 2^123.4615 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
