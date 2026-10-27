# Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **124.012** (total work <= 2^124.012 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 assumed indistinguishable from an ideal source), peak memory
< 2^62 bytes, preprocessing at most 2^49.06 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 lists every change from the parent)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and their results](evidence.md); [executable checks, production lazy walk, toy runs and the ledger](experiment.py)
- [Kernel generators, emulator and checks](kernels/) (`gen_lazy.py` emits `lazy.s` 3,527 and `root.s` 3,425 instructions)

Built on mpjunior92's gt-dprestart-replay (submission 0b8665ae, 2^125.111). One design change: the walk is
evaluated lazily (tag tracing). A walk keeps an explicit base y and a node of a prefix-product trie of the
288 twisted multipliers (R = 24, depth 6); the class key of the current element y * Pi_node needs only the
w^1 coefficient of the product, an Fp2-linear form in y whose coefficients the node stores, so a step costs
3,697 instructions instead of 16,578, and one full multiplication makes the base explicit every sixth step:
0.3042 units per step instead of 0.5503. This revision of gt-lazy-trie-rho (a77744bb, 2^124.183) changes one
thing: the key kernel pairs the six terms of that linear form (Winograd), 9 products instead of 18, with the
x-only term stored per base and the gamma-only term stored per trie node. Analytic generic-attack bound, not a break.

The success argument charges one explicit 0.01 allowance (heuristic H1+, proof.md section 3.4) for every cap
failure, an early threshold and any drift of the walk's frequencies beyond the 1.25 factor the ledger already
carries on its fruitless-cycle charge, and takes no independence between iterations, reports, restarts or trails:
the restart cap is charged at its own worst case as steps spent on restarts, and the pair budget by Markov on an
expected count. The drift table of section 3.4 gives the success at each drift size, and section 4 floors the
success after every loss at 0.39010.
