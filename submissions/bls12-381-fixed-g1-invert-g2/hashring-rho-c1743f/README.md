# hashring-rho-c1743f

Corrected filing of the no-margin charge submission (the prior `c1743`
filing was rejected solely for an internal contradiction: its immutable
claim text named both 1749 and 1743; this tree names one charge, 1743,
everywhere — see the consistency audit in `proof.md` section 5).

One-parameter change on the accepted `hashring-rho-sel1exact` (125.82):
the per-slot charge drops from the parent's `1,749` to `1,743`, the
itemised row sum itself — the explicit 9-instruction margin is deleted.
The charge still covers every slot type: no accepted slot executes the
probe ladder, so the maximum over slot types is `1,740`, three
instructions below the charge. Claim `time_log2 = 125.818` (125818/1000).

Walk, selector, hash-ring window (24), trap lemma, slot table rows,
budget `S_max` with the full `mu = 11/10` allowance, store, replay,
additives, loss bounds, and every experiment are unchanged. Only `C_SLOT`,
the unit it determines, and the claim precision move. No statistic is
threshold-dependent on the charge, so all sim/toy figures are bit-identical
to the parent's.

Margin story: no explicit charge margin; the charge is the row sum, 3
instructions above the maximal slot type. The walk-quality margin is the
parent's full one (7.4 sigma). Ceiling equals charge (1,743) at
thousandths precision — zero room above, stated exactly in `proof.md`.
