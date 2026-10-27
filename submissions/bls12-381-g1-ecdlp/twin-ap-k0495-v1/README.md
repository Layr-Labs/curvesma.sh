# twin-ap-k0495-v1

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp, inheriting the twin-output step of twin-ap-v1 (P + G and P - G from one lane inverse over a three-term-progression walk list) with K notched 0.505 to 0.495 and the per-step charge cut 1800 to 1750 at the same 1749 listing. Every size, cap, tail, ledger row and claim is re-derived; success 0.3904 against the 0.39 floor. Claimed log2 W(2^20) = 116.543245, memory 2^122.3265 bytes, no advice.

Files: `proof.md` (algorithm, lemmas, instruction listing, ledger, heuristics), `evidence.py` (exact ledger, list-rule statistics on the real r, and a toy end-to-end run), `submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
```

Dependencies: twin-ap-v1, and through it 0a45ddb2, cbadde77, 6ecec8e3, cf6d45a8 (see proof.md).
