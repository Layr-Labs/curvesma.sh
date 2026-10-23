# bls12-381-g1-ecdlp: Cheon's reduction with an exact baby-step giant-step

Cheon's reduction with A_d, d = 2^20 dividing r-1, followed by a deterministic baby-step giant-step in the
order-(r-1)/d subgroup of F_r^*, with classes under the six units and a full-key radix sort for matching.

Claimed log2 W(2^20) <= 118.23 (11823/100, rounded up from 118.221865), success probability exactly
m*G/M = 0.400000000 by counting (no walk, distinguished-point or hashing heuristic), preprocessing 2^87.650917
units included in the total, peak memory below 2^123.48 bytes, no advice.

The argument and the itemised worst-case ledger are in [proof.md](proof.md); [claim.json](claim.json) lists the
claim fields. [experiment.py](experiment.py) recomputes every production number and runs the same algorithm on a
small curve for every secret, checking the recovered count against d*6*m*G exactly.
