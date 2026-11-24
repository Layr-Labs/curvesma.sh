# BLS12-381 `bls12-381-g1-ecdlp` (`G_1` ECDLP with Auxiliary Powers, `n = 2^20`, `W(n) = W_cold` for all `n >= 0`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`Arch2-Prod-ZeroSlice-Buf128K-B2-U32-Bin21-Defer8-LSE`, `MODEA_V08_PSURV_K1_Q25_C4`)

Single-command pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'
```

Expected output: JSON with `"checks_ok": true` across all 16 checks, `"manifest_subbounds_ok": true`, and `"track_claim": "111.562518"` (`G_1` units, `C_G1 = 2,000` Ir/u = `10x` dense-`F_{p^12}` core units `108.240590`).
