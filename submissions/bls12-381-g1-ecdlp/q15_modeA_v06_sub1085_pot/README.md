# BLS12-381 `G_1` ECDLP with Auxiliary Powers (`bls12-381-g1-ecdlp`, `n = 2^20`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`Arch2-Prod-ZeroSlice-Buf128K-B2-U32-Bin21-Defer8`, `MODEA_V06_SUB1085`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "111.812319"`.
