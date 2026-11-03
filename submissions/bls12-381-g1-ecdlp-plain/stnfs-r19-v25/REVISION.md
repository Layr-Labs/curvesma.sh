# Changes against our d605c287 (124.548670069062)

Every number is recomputed by `audit_g1.py`; no prior verdict is a premise.
A = our plain R18 line (sub-cell descent integration), B = the port of
gbotrel 42e80f80's collection levers, C = new here. Deltas are cumulative on
this package's model (G1 bits).

| # | change | delta |
|---|---|---:|
| B1 | 40,000-draw pool on the 19/50 region (42e80f80 seeds 33/34, regenerated: 0 mismatches); allowance 1.06905 -> 1.0312 | -0.1314 (B1-B5) |
| B3 | Line sieve to 2^48 (scores only); survivors: trial division below 2^20, two 99.9% band-strip ECM rounds priced per cell | |
| B4, A5 | Row tail over the 7 occupied owners (N + 21 < 7 x 2^56): 1,317 counted, 1,320 charged; C_ENTRY 20.7; blocks (22, 4) | |
| B5 | alpha_g 0.844 -> 0.85 (more conservative) | |
| B2 | Descent inputs at 2-SE lower bounds: 2^63.97 relations, means 738.0/497.9 | +0.0038 |
| A1, A3 | Descent P integrated over 24 measured cf sub-cells and 10 boot deciles, boost capped 4x/1.5x; complete-orbit relation-weighted UDUP certificate (R18b: certified >= 1.11 nodes, >= 0.67 boot vs 0.25) with the live `udup_relw` check | -0.1158 |
| C4 | Red line h: missing rates as 1/50-bit band upper sums (base sums were 5-8% low) times 1.002, Chernoff configuration bound under H5m, Bonferroni keep-probability at thin nodes, node counter 2^25.7 | +0.0003 |
| C1 | Missing logs descended as prime special-q nodes in their own lattices (no composite Q q'), charged first-bin work and failure per 5-bit sub-bin; shared-set repeat sums sum_I w_I(i) w_I(j) m_I in the repeat term (2^-7.78) | +0.0007 |
| A2, C2 | Plan and boot re-searched on this model (47 entries moved) | -0.0910 |
| A4, B6 | Meter sized for modelled success >= 0.50 at rho = 0 and >= 0.39 at rho = 0.05 (H5', exchangeable failures) | +0.0055 |
| C3 | 13 modules: indentation halved and long docstrings/comments shortened to fit the archive cap; code AST (docstrings removed) identical | 0 |

Ledger 121.2253 -> **<= 2^120.897292**; G1 124.548670069062 -> **124.220662069061269**.
Collection residual at +half-ulp constants 1.000015x, ledger 1.000014x, inside the
1001/1000 reserve. Slacks (297, 233), excess 1.3, merge weight 300, pad 1.10,
kernels and the wrapper are unchanged.
