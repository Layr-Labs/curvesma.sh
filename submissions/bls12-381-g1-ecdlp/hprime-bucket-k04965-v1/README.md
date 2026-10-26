# hprime-bucket-k04965-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp with H-prime semiprime walks (one batched affine addition per step) and a bucket-table matching stage. Claimed log2 W(2^20) = 116.759775, success 0.3908 (K = 0.4965), memory 2^123.3242 bytes, no advice.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger, heuristics), `evidence.py` (exact ledger, a toy end-to-end run, and an exact coverage census on toy primes), `submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: cbadde77, 6ecec8e3, cf6d45a8, 0a45ddb2, 17ed8646 (see proof.md).
