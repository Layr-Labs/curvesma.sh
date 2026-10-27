# hprime-bucket-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks (one batched affine addition per step) and a bucket-table matching stage. Claimed log2 W(2^20) = 116.793080, success 0.40, memory 2^123.3242 bytes, no advice.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger, heuristics), `evidence.py` (exact ledger plus a toy end-to-end run), `submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8 (see proof.md).
