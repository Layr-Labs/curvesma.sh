# Joint-cell survivor bound (H-JC)

node_cost charges every cell passing the sieve threshold T with the full cofactor chain. The survivors per cell, so far bounded on six cells at the record minima (QBINS_LO), are bounded on 33 joint cells (cf strips x cg bands inside each strip): sv = sum_k W[k] surv(48 c_lo + cf_lo[k] - q_f, T) surv(12 c_lo + 6 log2 u + cg_lo[k] - q_g, T). surv is non-increasing in the norm, so a cell above a corner survives at most as often as the corner; the weights are joint frequencies (no independence of cf and cg).

Data. Three campaigns of 160 special-q ideals x 384 ball cells (61,440 cells; sizes 56-200 bits, both sides, radii 8-128). A fixes the edges (cf quantiles 0.001, 1/6..5/6; per strip cg quantiles 0.001, 0.01, 0.1, 0.5), B the counts, C is held out. Corner = lower edge - shifts (10, 4).

Weights. No empty cell (two cells empty in B are joined to a neighbour at the lower corner). W = share + 2.9827 SE (level 0.05 simultaneously, Bonferroni); SE = larger of binomial and between-ideal (160 per-ideal shares), so clustering is charged (design effect at most 1.2 below 100 counts, 5.7 above); below 100 counts W is at least the exact Poisson limit at that level. The top cell (share 0.079) takes the rest.

Tails. Without a lower edge the corner is the lowest value of ANY record minus the shift: cf -82.6 (the package's), cg -57.1 (C; package -49.4, A -50.0, B -51.3). Mass below it is not assumed zero (none seen in 184,320 cells): at most 4^-k of the cell's weight lies more than 16k bits below the corner (k <= 40), the last 4^-41 at probability 1. Cells below cg -38/-46/-54/-62: A 36/2/0/0, B 35/2/0/0, C 48/4/3/0; below cf -43/-59/-75/-93: A 61/18/2/0, B 55/6/1/0, C 64/5/0/0.

Every node type. For each of the 49 sieve rows, the survivor expectation computed from the sampled cells of the two nearest sizes on the same side (3,072 cells each), offsets lowered by (10, 4), is at most 0.32 of the bound on B and on C (0.18 unshifted), and stays below it with 13 more bits of shift on both sides.

H-JC. For every node type the law of (cf, cg) of a uniform ball cell, lowered by the shifts, puts at most W[k] on each cell but the top one and at most 4^-k of a tail cell's weight 16k bits below its corner.

Limits. A model change: with the plan fixed the two-log row falls by 1.05 bits. The ladder rate is assumed (no data in its range). Independent sampler (exact norms, LLL), not shipped (size).

i j count ideals 1e7.SE W cf_lo cg_lo tail(1 cf, 2 cg)
0 1 2 2 230 0.0001759 -92.7 -61.2 3
0 2 19 17 779 0.0007059 -92.7 -33.7 1
0 3 35 29 1024 0.0010417 -92.7 -23.4 1
1 0 16 14 700 0.000602 -53.0 -61.2 2
1 1 94 71 1589 0.0020958 -53.0 -41.9 0
1 2 961 160 5109 0.0171651 -53.0 -33.7 0
1 3 4175 160 13471 0.0719706 -53.0 -23.4 0
1 4 5241 160 19579 0.0911425 -53.0 -14.2 0
2 0 7 7 431 0.0003107 -6.7 -61.2 2
2 1 95 67 1654 0.0022653 -6.7 -40.8 0
2 2 926 159 5723 0.0167788 -6.7 -31.4 0
2 3 3985 160 12747 0.0686622 -6.7 -21.5 0
2 4 5261 160 20109 0.0916261 -6.7 -12.2 0
3 0 14 13 628 0.000505 3.4 -61.2 2
3 1 109 77 1809 0.0023137 3.4 -40.4 0
3 2 972 159 5649 0.0175053 3.4 -30.7 0
3 3 4336 160 12007 0.0741544 3.4 -20.8 0
3 4 5061 160 16986 0.0874395 3.4 -11.5 0
4 0 9 9 488 0.0003595 11.1 -61.2 2
4 1 70 56 1361 0.0016078 11.1 -38.6 0
4 2 937 158 5958 0.0170278 11.1 -30.9 0
4 3 4027 160 11642 0.0690162 11.1 -20.8 0
4 4 5093 160 15958 0.0876536 11.1 -11.2 0
5 0 20 17 791 0.0007185 18.4 -61.2 2
5 1 98 72 1649 0.0022417 18.4 -37.5 0
5 2 867 160 6586 0.0160757 18.4 -30.4 0
5 3 4088 160 14757 0.0709382 18.4 -20.7 0
5 4 5095 160 16258 0.0877757 18.4 -11.2 0
6 0 5 5 364 0.0002597 27.2 -61.2 2
6 1 97 73 1634 0.0022057 27.2 -41.9 0
6 2 959 160 8093 0.0180227 27.2 -32.0 0
6 3 3912 160 21534 0.070095 27.2 -21.2 0
6 4 4854 160 26082 0.0095418 27.2 -11.4 0
