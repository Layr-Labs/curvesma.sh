# BLS12-381 `recover-Q` (`bls12-381-fixed-g1-invert-g2`) — `2^{106.790975}`

Run `python3 -I -B audit_cm.py <<< '{}'` (`part="core"`, default) and `python3 -I -B audit_cm.py <<< '{"part":"diagnostics"}'` (or `{"part":"all"}`) to verify all checks (`checks_ok == True`).
`proof.md` section 0 lists what changed relative to gbotrel's `stnfs-q15-sub107-v5b-bigstack` (`2^{106.7915}`): priced Jacobian-NAF ECM curve generation (`STEPPER_MULS = 4,900`, complete over `Z/NZ`, live checks including forced exceptional cases) and gate floor `0.99998`. The rest of `proof.md` gives the complete reduction. See also `KERNELS.md`, `KERNELS_LA.md`, `DESCENT.md`, `ADAPTIVE.md`, `STORAGE.md` and `RANK-BRIDGE.md`.
