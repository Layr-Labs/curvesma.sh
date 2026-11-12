# Filtering scales 2^21, 2^22 and 2^23

Same profile as [6] (B = 2^58, norms 417 / 471 bits, excess 1.095, purge to 0.0001 U', kmax 64, merge to weight 400, |a|, |b| <= 255); seeds 1..8 fixed before any run above 2^19, none dropped. 2^21 and 2^22: gbotrel's filt29, 2^23: our filtx, both exact C transcriptions of filtsim.py's row loop, purge and merge (same MT19937 stream, float order and set/dict/heap tie-breaks) fed by filtsim.build_system's own prelude (Python filtsim needs ~1.4 GB and ~4 min at 2^19). filt29 reproduces all 18 stored runs of 2^17..2^19 exactly, filtx all 34 of 2^17..2^22 (N, max coefficient, mean weight). filtbn.b85 ships every run; [6] recomputes the table and the charge, [6r] replays filtsim at 2^12.

2^21: N(400) per seed 1:29468 2:28642 3:28277 4:29244 5:28721 6:28437 7:28706 8:28780; mean 0.01373, SE 0.00007, mean + 2 SE 0.013858, max coefficient 183.
2^22: N(400) per seed 1:54741 2:53820 3:53908 4:56068 5:53952 6:55308 7:55416 8:54266; mean 0.01304, SE 0.00007, mean + 2 SE 0.013179, max coefficient 200.
2^23: N(400) per seed 1:105082 2:104763 3:103219 4:104854 5:104009 6:102997 7:104724 8:103956; mean 0.01242, SE 0.00003, mean + 2 SE 0.012488, max coefficient 195.
H3 rule unchanged (larger mean + 2 SE of the two largest scales): 0.013179, charged NU = 0.01318; linear algebra 2^88.661885, generator 2^81.648859.
