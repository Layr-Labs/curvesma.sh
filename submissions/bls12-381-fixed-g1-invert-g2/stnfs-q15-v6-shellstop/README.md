# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.768400}` (Shell-Ordered Early Stop on the Mode-A Big-Stack package)

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (or `{"part":"all"}`) to verify all checks (`checks_ok == True`).
See `proof.md` section 3.0 for what changed in this revision (shell-ordered early-stop collection, priced Jacobian-NAF ECM curve generation, gate floor `0.99998`), the rest of `proof.md` for the complete reduction, `KERNELS.md` + `KERNELS_LA.md` for all raw-opcode-verified ARMv8.1-A + `STADDB`/`STADD` blocks,
and `DESCENT.md` / `ADAPTIVE.md` / `STORAGE.md` / `RANK-BRIDGE.md` for full details. Built on gbotrel's `stnfs-q15-sub107-v5b-bigstack` (`2^{106.7915}`).
