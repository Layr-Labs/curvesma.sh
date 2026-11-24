# BLS12-381 `recover-P` (`bls12-381-fixed-g2-invert-g1`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`MODEA_V03_A3_Q25_C4`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "109.059290"`.
