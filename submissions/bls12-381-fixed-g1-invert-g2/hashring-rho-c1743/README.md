# hashring-rho-c1743

One-parameter tightening of the accepted `hashring-rho-sel1exact` (125.82):
the per-slot charge drops from `1,749` to `1,743`, the itemised row sum
itself — the explicit 9-instruction margin is deleted. The charge still
covers every slot type: no accepted slot executes the probe ladder, so the
maximum over slot types is `1,740`, three instructions below the charge.
Claim `time_log2 = 125.818` (125818/1000).

Walk, selector, hash-ring window (24), trap lemma, slot table rows, budget
`S_max`, store, replay, additives, loss bounds, and every experiment are
unchanged. Only `C_SLOT`, the unit it determines, and the claim precision
move. No statistic is threshold-dependent on the charge, so all sim/toy
figures are bit-identical to the parent's.

Margin: none explicit; the charge is the row sum, 3 instructions above the
maximal slot type. Ceiling room is stated exactly in `proof.md`.
