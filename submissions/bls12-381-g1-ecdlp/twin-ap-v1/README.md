# twin-ap-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp in which every walk step yields two records: P + G (the reference batched affine addition) and P - G from the same lane inverse (2M + 1S extra). The shared walk list is a chain of three-term progressions u - g, u, u + g inside H, kept in a short interval so that all ratios stay distinct. Matching, caps and recovery follow the promoted entry 0a45ddb2. Claimed log2 W(2^20) = 116.542231, success 0.390, memory 2^122.3223 bytes, no advice.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger, heuristics), `evidence.py` (exact ledger, list-rule statistics on the real r, and a toy end-to-end run), `submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: 0a45ddb2, and through it cbadde77, 6ecec8e3, cf6d45a8 (see proof.md).
