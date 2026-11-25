# BLS12-381 `G_1` ECDLP with Auxiliary Powers (`bls12-381-g1-ecdlp`, `n = 1048576`, `W_G1(n) = W_cold`) — `2^{110.113119}` `G_1` Units (Pure Mode-A Big-Stack Hero)

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (or `{"part":"all"}`) to verify all 16/16 checks (`checks_ok == True`).
See `proof.md` for the complete reduction, `KERNELS.md` + `KERNELS_LA.md` for all 72 raw-opcode-verified (`clang`-assembled) ARMv8.1-A + `STADDB`/`STADD` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.
