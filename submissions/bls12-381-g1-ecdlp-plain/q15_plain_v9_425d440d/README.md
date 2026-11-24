# BLS12-381 Plain `G_1` ECDLP (`bls12-381-g1-ecdlp-plain`) at `2^{110.751618}` `G_1` Units (Sub-104.8 Pure Mode-A `K = 1` Interval `2.0 * psurv512_hi` Descent + Complex-Place `L_1(b)` Convex Ball)

Run `python3 -I -B audit_cm.py <<< '{}'` (`10/10` core checks incl. live `exact_deg72_h5`) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (`9/9` secondary diagnostics incl. full 6-replay `exact_deg72_h5_full`; or `{"part":"all"}` for all `19/19` checks, `checks_ok == True`, `track_claim == "110.751618"`, where `1` internal dense-`F_{p^12}` core unit at `C = 20,000` Ir equals `10` `G_1`-reference units at `C_G1 = 2,000` Ir).
See `proof.md` for the complete reduction, `KERNELS.md` for all 73 raw-opcode-verified ARMv8.1-A + `STADDH` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.
