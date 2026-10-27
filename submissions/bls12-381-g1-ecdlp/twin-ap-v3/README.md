# twin-ap-v3

Cheon's d = 2^20 attack on bls12-381-g1-ecdlp in which every walk step yields two records: P + G (the reference batched affine addition) and P - G from the same lane inverse (2M + 1S extra). The shared walk list is a chain of three-term progressions u - g, u, u + g inside H, kept in a short interval so that all ratios stay distinct. Matching, caps and recovery follow the promoted entry 0a45ddb2.

This is twin-ap-v2 (42113fc) with one change: the twin-list scan cap is halved from 2^100 to 2^99 integers, so the list row falls from 2^107.03 to 2^106.03 units. The measured P - G side output (1329 instructions per lane, organizers' pipeline) and K = 0.4944 are unchanged. Claimed log2 W(2^20) = 116.508149, success 0.390, memory 2^122.3223 bytes, no advice.

Files: `proof.md` (algorithm, lemmas, measured side output, instruction listing, ledger, heuristics), `evidence.py` (exact ledger, list-rule statistics on the real r, and a toy end-to-end run), `measure/` (C kernel and scripts for the instruction count), `submission.json`, `claim.json`.

Recompute everything:

```sh
echo '{"mode":"all","seed":2026,"trials":400}' | python3 -I evidence.py
bash measure/reproduce.sh <repository root>     # Docker with linux/arm64 support
```

Dependencies: twin-ap-v2 at 42113fc, and through it 017c2fd6 (twin-ap-v1), e26c3499, 0a45ddb2, cbadde77, 6ecec8e3, cf6d45a8 (see proof.md).
