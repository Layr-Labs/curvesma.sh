# Conditional next-relation evidence (H5)

H5 requires P(F_i | M>=i, F_1..F_(i-1), pre-entry history) <= psi and a conditional child-work bound. We measure both directly. All runs use the production deterministic qualification order: relations are tried in cell-enumeration order, and the node stops at the first relation whose children all succeed. Every recursive descent runs to completion and none is censored.

**Campaign.** The tower is the reduced one over the production BN base ring, with f_t=x^4+w*x^2+w^2 and g_t=x^2+4w. Each node enumerates its complete radius-2 ball (4996 cells). Norms are factored exactly in the finite windows, and every child's degree-one root is checked. The root window is (lb,n,k)=(12,root_window,2). A child with q of b bits uses n=max(12,b-d) and k=1, so windows strictly shrink, and ideals <=2^12 are known. There are five window settings: d=1 is the production-like one-bit shrink at roots f20, g17 and f22, and d=2 and d=3 are tighter, with more failures. In total there are 2142 roots (all succeed), 10,000 entered nodes, 4602 first attempts with children and **627 attempts after an earlier failed relation (139 failures)**.

**Envelopes.** psi for an attempt is the sum, over ALL of its children, of the child-row failure rate (row = side and bit size). The child-work envelope is the sum over all children of the row's mean subtree work, with full charge even when the first child fails. Both are calibrated ONLY on children of first attempts. The conditional test set never enters the calibration. O/E means observed divided by envelope. Intervals are 95% root-cluster bootstrap.

| setting | first n | first fail | after-fail n | fail | rate | 1-sided 95% upper | O/E fail (95%) | ratio to first (95%) | work O/E (95%) |
|---|---|---|---|---|---|---|---|---|---|
| f20, d=1 | 1989 | 0.120 | 205 | 26 | 0.127 | 0.170 | 0.830 (0.524-1.151) | 0.861 (0.544-1.193) | 0.987 (0.927-1.049) |
| g17, d=1 | 1131 | 0.085 | 75 | 8 | 0.107 | 0.180 | 1.005 (0.514-1.592) | 1.059 (0.541-1.677) | 0.881 (0.797-0.975) |
| f22, d=1 | 245 | 0.118 | 22 | 2 | 0.091 | 0.241 | 0.507 (0.000-1.211) | 0.599 (0.000-1.431) | 1.007 (0.810-1.212) |
| f20, d=2 | 971 | 0.207 | 232 | 72 | 0.310 | 0.363 | 0.997 (0.797-1.183) | 1.076 (0.860-1.277) | 0.941 (0.899-0.982) |
| f20, d=3 | 266 | 0.271 | 93 | 31 | 0.333 | 0.418 | 0.795 (0.588-0.968) | 0.904 (0.668-1.100) | 0.910 (0.853-0.982) |
| **pooled** | 4602 | 0.138 | **627** | **139** | 0.222 | | | **0.964 (0.839-1.095)** | **0.950 (0.918-0.983)** |

After an earlier failed relation, the next relation fails as often as its children's unconditional row rates predict. The pooled ratio to first attempts is 0.964. Its 95% interval is 0.839-1.095 and contains no detectable elevation. Raw rates are higher in the d=2/3 settings because nodes reached after a failure carry more and harder children. psi tracks that row mix, which is why the comparison is row-specific and not a pooled rate. Conditional child work is at or below the fresh same-class envelope in every setting (pooled 0.950, 95% up to 0.983).

**Shared descendants.** For 248 after-failure attempts we recorded every child ideal's key. None shares a child ideal with an earlier relation of the same node, and none reuses an ideal already entered elsewhere in the run. Consecutive relations of a node have distinct window primes. A child ideal is determined by the prime and the cell's root, so sharing needs the same prime to divide both relations' norms with the same root. Production window bins hold vastly more ideals than the toy bins (where sharing is already zero), and production aborts any repeat as a failure, charged in the repeat loss. Difficult descendants therefore cannot be shared silently between consecutive relations.

**Production envelope slack.** We also evaluated the package formula (threshold_bound.failure_bound at the row's model mu, r=1/2) as the child envelope. Observed/envelope is 0.11-0.35 in every setting, both on first attempts and after an earlier failure. The conditional rates stay below the formula envelope by a factor of 3 or more. The production plan caps every node at phi<=0.0397 with node psi<=0.4581 (boot expected psi 0.4311).

**Trend.** Across window shrink d=1,2,3, the first-attempt failure rate rises from 0.12 to 0.27. The conditional O/E stays at 0.795-0.997 and the ratio to first attempts at 0.861-1.076, with no trend in d. Root size 17-22 bits shows the same. Conditioning effects do not grow as windows widen or failures become more common. Production rows have far larger M (raw mu >= 333 against 5-25 here), so each next relation draws on a larger pool of independent children.

**Scope.** The tower uses degree 4 instead of 8, ideals of 13-22 bits and window primes of 12-21 bits. These are reduced windows, not the production factor base. Transfer to production degree and prime sizes is H5. The recurrence uses these data at their stated strength.

**Replay.** part=h5 (h5cond.py with h5data.b85, 3 KB) recomputes every number above from per-root tallies in under 0.1 s, including the bootstrap. The campaign generator h5camp.py ships in our accepted stnfs-bn10ez (declared dependency).
