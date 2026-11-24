# BLS12-381 Plain `G_1` ECDLP (`bls12-381-g1-ecdlp-plain`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`Arch2-Prod-ZeroSlice-Buf64K-B2-U32-Orb1b`, `MODEA_V03_A3_Q25_C4`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks and `"track_claim": "112.161118"` (`G_1` units at `C_G1 = 2,000` Ir/u = `10x` core dense-`F_{p^12}` units `108.839190`).
