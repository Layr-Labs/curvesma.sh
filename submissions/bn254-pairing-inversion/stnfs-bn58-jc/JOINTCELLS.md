# Joint-cell survivor bound (hypothesis H-JC)

Change. node_cost charges every cell passing the sieve threshold T with the full cofactor chain. Expected survivors per cell were bounded on six cf sextiles, each at the lowest offsets on record (QBINS_LO); now on 35 joint cells k (7 cf strips, 5 cg bands cut inside each strip): sv = sum_k W[k] surv(48 c_lo + cf_lo[k] - q_f, T) surv(12 c_lo + 6 log2 u + cg_lo[k] - q_g, T). surv, T, yield side, shifts, prices and the sieve/ECM choice are unchanged.

Bound. surv is non-increasing in the norm, so a cell above a corner survives at most as often as the corner; joint frequencies as weights need no independence of cf and cg.

Data. Edges: cf quantiles 0.001, 1/6..5/6 and, per strip, cg quantiles 0.001, 0.01, 0.1, 0.5 of campaign A, fixed before campaign B (160 fresh special-q ideals, 61,440 ball cells, sizes and radii pooled), which gives the counts. Corner = lower edge minus the package's shifts (10, 4); with no lower edge, the lowest value on record minus the shift (cf -92.7, cg -55.4): 0.14% and 0.19% of the mass.

Weights. Simultaneous over the 35 cells at level 0.05 (Bonferroni, z = 2.9827): W = share + z SE (SE over ideals, not below binomial), empty cell 6.551/n, rounded up. The top cell (measured 0.079) takes the remaining 0.0145.

H-JC. For every node type, the joint law of (cf, cg) of a uniform ball cell, lowered by the shifts, puts at most W[k] on each cell but the top one and nothing below (-92.7, -55.4).

Limits. (1) A change of the model, not a re-measurement: with the plan fixed the two-log row falls by 1.08 bits (yield, phi, psi, E[N], losses unchanged), on H-JC alone. (2) The tail corners remain record minima (cg minima -50.0, -51.3, -57.1 on three campaigns; package -49.4); tables with tails 5.8 bits apart give rows 0.005 bits apart. (3) Transfer from 160 pooled ideals to every node type rests on the shifts only, as for the package's bins. (4) Other edge and fill sets move the row by at most 0.014 bits; the simultaneous allowance costs 0.006. (5) The campaigns used an independent sampler (exact norms, LLL), not the package's; pooled statistics agree with the package's bins within about 1 SE. Not shipped (size); the rows below are the constant QSV.

Cells: i j count 1e7.SE W cf_lo cg_lo
0 0 0 0 0.0001067 -92.7 -55.4
0 1 2 230 0.0001013 -92.7 -41.9
0 2 19 779 0.0005418 -92.7 -33.7
0 3 35 1009 0.0008708 -92.7 -23.4
0 4 0 0 0.0001067 -92.7 -14.2
1 0 16 651 0.0004546 -53.0 -55.4
1 1 94 1577 0.0020003 -53.0 -41.9
1 2 961 5006 0.0171345 -53.0 -33.7
1 3 4175 13138 0.0718712 -53.0 -23.4
1 4 5241 18650 0.0908656 -53.0 -14.2
2 0 7 431 0.0002424 -6.7 -55.4
2 1 95 1597 0.0020227 -6.7 -40.8
2 2 926 4915 0.0165378 -6.7 -31.4
2 3 3985 12945 0.0687213 -6.7 -21.5
2 4 5261 19493 0.0914426 -6.7 -12.2
3 0 14 609 0.0004095 3.4 -55.4
3 1 109 1717 0.0022862 3.4 -40.4
3 2 972 5122 0.0173481 3.4 -30.7
3 3 4336 11083 0.0738787 3.4 -20.8
3 4 5061 15559 0.0870138 3.4 -11.5
4 0 9 488 0.0002922 11.1 -55.4
4 1 70 1361 0.0015453 11.1 -38.6
4 2 937 5891 0.0170078 11.1 -30.9
4 3 4027 10524 0.0686827 11.1 -20.8
4 4 5093 14677 0.0872716 11.1 -11.2
5 0 20 728 0.0005426 18.4 -55.4
5 1 98 1610 0.0020753 18.4 -37.5
5 2 867 6125 0.0159382 18.4 -30.4
5 3 4088 14100 0.0707421 18.4 -20.7
5 4 5095 15551 0.0875649 18.4 -11.2
6 0 5 364 0.00019 27.2 -55.4
6 1 97 1602 0.0020566 27.2 -41.9
6 2 959 7711 0.0179087 27.2 -32.0
6 3 3912 20361 0.0697449 27.2 -21.2
6 4 4854 24604 0.0144824 27.2 -11.4
