# BLS12-381 Plain `G_1` ECDLP (`bls12-381-g1-ecdlp-plain`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`MODEA_V02_A2_GLOBAL_2X`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "112.536418"` (`C_G1 = 2,000` `G_1`-reference units = `10x` dense-`F_{p^12}` core units).
