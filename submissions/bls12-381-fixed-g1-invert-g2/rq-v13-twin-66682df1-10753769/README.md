# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `(4, 3)` `Q(zeta_15)` CM-STNFS (`Arch2-Prod-ZeroSlice-Buf128K-B2-U32-Bin21-Defer8-LSE`, `MODEA_K2_L1B_C5IV_K1_Q25`)

Self-contained pure-Python-stdlib verification (`python3 3.11+`, zero external packages):

```bash
python3 -I -B audit_cm.py <<< '{}'                     # default core audit (initial_inputs)
python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}' # on-demand extended diagnostics (incl. 6 full-ball replays)
```

Expected output: JSON with `"checks_ok": true` across all core and on-demand diagnostic checks and `"track_claim": "107.537690"`.
