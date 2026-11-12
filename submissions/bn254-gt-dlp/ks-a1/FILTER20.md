# Filtering scales 2^21 and 2^22

Same profile as [6] (B = 2^58, relation norms 417 / 471 bits, excess 1.095, purge to 0.0001 U', kmax 64, merge to weight 400, |a|, |b| <= 255),
seeds 1..8 fixed before any run above 2^19, none dropped. Runs: filt29, a C transcription of filtsim.py's row loop, purge and merge (same MT19937
stream, float order and set/dict/heap tie-breaks) fed by filtsim.build_system's own prelude; the package's Python filtsim needs ~1.4 GB and
~4 min at 2^19, too large for the hosted 20 s / 256 MB check. Identity: filt29 reproduces N, max coefficient and mean weight of all 18 stored
runs (2^17 seeds 1-2, 2^18 and 2^19 seeds 1-8) exactly. filtbn.b85 ships every run; [6] recomputes the table and the charge from it, [6r]
replays filtsim at 2^12.

2^21: N(400) per seed 1:29468 2:28642 3:28277 4:29244 5:28721 6:28437 7:28706 8:28780; mean 0.01373, SE 0.00007, mean + 2 SE 0.013858, max coefficient 183.
2^22: N(400) per seed 1:54741 2:53820 3:53908 4:56068 5:53952 6:55308 7:55416 8:54266; mean 0.01304, SE 0.00007, mean + 2 SE 0.013179, max coefficient 200.
H3 rule unchanged (larger mean + 2 SE of the two largest scales): 0.013858, charged NU = 0.01386 (was 0.01669). N <= 2^45.5350; linear algebra
2^88.807039, generator 2^81.728608.
