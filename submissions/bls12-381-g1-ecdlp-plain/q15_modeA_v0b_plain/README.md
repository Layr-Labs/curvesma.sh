# BLS12-381 Plain `G_1` ECDLP (`bls12-381-g1-ecdlp-plain`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`MODEA_V0_ORIG_RADII`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "113.847018"` (`C = 2,000` `G_1`-reference units = `10x` dense-`F_{p^12}` core units).
