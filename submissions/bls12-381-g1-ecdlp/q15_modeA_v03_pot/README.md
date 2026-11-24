# BLS12-381 `G_1` ECDLP with Auxiliary Powers (`bls12-381-g1-ecdlp`, `n = 2^20`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`MODEA_V03_A3_Q25_C4`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "112.381218"` (`G_1` units at `C_G1 = 2,000` Ir/u = `10x` core dense-`F_{p^12}` units `109.059290`; auxiliary-power-free, constant for all `n >= 0`).
