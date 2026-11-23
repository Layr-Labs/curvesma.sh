# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`MODEA_V0_EA_TAU_MU700`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "109.659690"`.
