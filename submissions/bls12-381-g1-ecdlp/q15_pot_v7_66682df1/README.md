# BLS12-381 `bls12-381-g1-ecdlp` (`G_1` ECDLP with Auxiliary Powers, `n = 2^20`, `W(n) = W_cold` for all `n >= 0`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`Arch2-Prod-ZeroSlice-Buf128K-B2-U32-Bin21-Defer8-LSE`, `MODEA_K2_L1B_C5IV_K1_Q25`)

Self-contained pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'                     # default core audit (initial_inputs)
python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}' # on-demand extended diagnostics (incl. 6 full-ball replays)
```

Expected output: JSON with `"checks_ok": true` across all core and on-demand diagnostic checks, `"manifest_subbounds_ok": true`, and `"track_claim": "110.859618"` (`G_1` units, `C_G1 = 2,000` Ir/u = `10x` dense-`F_{p^12}` core units `107.537690`).
