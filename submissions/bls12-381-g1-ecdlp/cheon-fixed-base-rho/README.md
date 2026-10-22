# bls12-381-g1-ecdlp: Cheon with a tame-wild class walk

Cheon's algorithm using A_d, d=2^20 dividing r-1: a collision search over the order-(r-1)/d subgroup of F_r^*, with classes under the six roots of unity and fixed-base table multiplication per step.
Claimed log2 W(2^20) = 119.8 with success probability 0.39; heuristics and ledger in [proof.md](proof.md), claim in [claim.json](claim.json).
Analytic exploratory bound; the toy experiment checks bookkeeping only.
