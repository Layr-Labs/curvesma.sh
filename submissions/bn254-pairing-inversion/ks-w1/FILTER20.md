# Filtering scales 2^22 and 2^23, merge weight 200

Profile of [6] (B = 2^58, norms 417 / 471 bits, excess 1.095, purge to 0.0001 U', kmax 64, |a|, |b| <= 255), merge to weight 200 (W_ROW 201): the
argmin over 150..400 of the 2^16 LA cost (seeds 1-4), fixed before any result above 2^16 at w < 400 existed; C_ENTRY_K 8.42 and the unit share stay
the weight-400 measurement (our entry measure at 2^16 is not higher at 150..350). Seeds 1..8, none dropped. Every N(200) in filtbn.b85 is a real filtx run (none estimated); filtx is our exact C++
transcription of filtsim.py fed by build_system's prelude (same MT19937 stream, float order, set/dict/heap ties), all targets on one trajectory: it
reproduces the 34 stored N(400) runs of 2^17..2^22 and filtsim at 8 targets 200..1000 on 2^14..2^16 (2 seeds each). [6] recomputes table and
charge from filtbn.b85; [6r] replays filtsim at 2^12 to weight 200.

2^22: N(200) per seed 1:68553 2:67343 3:67521 4:70194 5:67548 6:69274 7:69346 8:67931; mean 0.01632, SE 0.00009, mean + 2 SE 0.016500, max coefficient 85.
2^23: N(200) per seed 1:131776 2:131357 3:129429 4:131458 5:130406 6:129143 7:131284 8:130271; mean 0.01557, SE 0.00004, mean + 2 SE 0.015657, max coefficient 83.
H3 unchanged (larger mean + 2 SE of the two largest scales): 0.016500, charged NU = 0.01650 (weight 400: 0.01318); linear algebra 2^88.498008, generator 2^82.004912.
