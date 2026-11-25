# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.760500}` (Pure Mode-A Big-Stack Hero)

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (or `{"part":"all"}`) to verify all 16/16 checks (`checks_ok == True`).
See `proof.md` for the complete reduction, `KERNELS.md` + `KERNELS_LA.md` for all 72 raw-opcode-verified (`clang`-assembled) ARMv8.1-A + `STADDB`/`STADD` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details.
