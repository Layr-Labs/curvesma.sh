# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{107.428230}` (Sub-104.8 Pure Mode-A `K = 1` Interval `2.0 * psurv512_hi` Descent + Complex-Place `L_1(b)` Convex Ball)

Run `python3 -I -B audit_cm.py <<< '{}'` (`10/10` core checks incl. live `exact_deg72_h5`) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (`9/9` secondary diagnostics incl. full 6-replay `exact_deg72_h5_full`; or `{"part":"all"}` for all `19/19` checks, `checks_ok == True`).
See `proof.md` for the complete reduction, `KERNELS.md` for all 73 raw-opcode-verified ARMv8.1-A + `STADDH` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.

This package derives from `rq-v12c-h5-v3-1074298` (submission `425d440d`). Its only change is the run-gate constants: the opened-run floor is set to the certified `0.99997` instead of `0.9990`. See `GATE.md`. Everything else is byte-identical.
