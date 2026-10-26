# hprime-bucket-hzero-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks
(one batched affine addition per step) and a bucket-table matching stage.
Base: the hprime-bucket-v1 package on HEAD (K = 0.52, accepted at 116.793080).
Sole score lever: the per-node fixed-work charge drops from a flat 2^35 to an
explicit count. Head zeroing (2^33 8-byte slots = 2^36 bytes) is 2^32 plain
X-register stp plus a 16x-unrolled loop (2^28 iterations at 6 instructions),
charged exactly as 2^32 + 6*2^28 with no new instruction types and no
re-itemized hot path; bucket-pointer setup and end-of-hop flushes stay at the
accepted 3x2^33 per node. The fixed-work row shrinks accordingly. Separately,
the preprocessing figure is upward-rounded (ceil to 4 decimals) and
integer-asserted, answering the 934a4565 rejection, which faulted
round-to-nearest preprocessing accounting. Algorithm, walks, routing,
matching, keys, caps, charges, and the toy procedure are the accepted package.

Claimed log2 W(2^20) = 116.792820, success 0.40, memory 2^123.3242 bytes,
no advice. Numbers recomputed exactly by evidence.py; proof.md quotes them.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger,
heuristics), `evidence.py` (exact ledger plus a toy end-to-end run),
`submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2 (see proof.md).
