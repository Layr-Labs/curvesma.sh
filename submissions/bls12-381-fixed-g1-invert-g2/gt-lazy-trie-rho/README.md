# Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **124.210** (total work <= 2^124.210 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 taken as an ideal source) under the
heuristics H1, H2 and H1+ with tolerance t = 2 (proof.md, evidence.md), peak memory
< 2^62 bytes, preprocessing at most 2^49.06 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0: the port and the fixes since the rejected revisions)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments, including the registered toy census](evidence.md) ([census driver](census.py)); [checks, production lazy walk, toy runs and the ledger](experiment.py)
- [Kernel generators, emulator and checks](kernels/) (`gen_lazy.py` emits `lazy.s` 4,964 and `root.s` 1,561 instructions)

Port to this track of mpjunior92's gt-lazy-trie-rho (a77744bb, 2^124.183 on the mirror track): the walk only
sees (g, z), g = e(P, Q0) has the same distribution as e(P0, Q) there, and the answer is Q = [x]Q0 (proof.md
section 1). Steps are lazy: an explicit base y plus a node of a prefix-product trie of the 288
twisted multipliers (R = 24, depth 6); the class key needs only the w^1 coefficient of y * Pi_node, an Fp2-linear
form in y the node stores, so a step costs 5,138 instructions and a full multiplication every sixth step: 0.3487
units per step. No parent figure enters this ledger. Analytic generic-attack bound, not a break.
