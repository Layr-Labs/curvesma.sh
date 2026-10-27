# Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **124.161** (total work <= 2^124.161 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 assumed indistinguishable from an ideal source), peak memory
< 2^62 bytes, preprocessing at most 2^49.06 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 lists every change from the parent)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and their results](evidence.md); [executable checks, production lazy walk, toy runs and the ledger](experiment.py)
- [Kernel generators, emulator and checks](kernels/) (`gen_lazy.py` emits `lazy.s` 4,964 and `root.s` 1,561 instructions)

Reuses mpjunior92's gt-lazy-trie-rho (submission a77744bb, 2^124.183) whole with one change: the multiplication
kernel `mul.s` uses asymmetric products with the fixed operand's combinations stored in its entry, 14,402 -> 13,761
instructions (block of six steps 41,821 -> 41,180); internal trie nodes keep their 12 coordinates beside the entry.
That package is built on gt-dprestart-replay (submission 0b8665ae, 2^125.111) with one design change: the walk is
evaluated lazily (tag tracing). A walk keeps an explicit base y and a node of a prefix-product trie of the
288 twisted multipliers (R = 24, depth 6); the class key of the current element y * Pi_node needs only the
w^1 coefficient of the product, an Fp2-linear form in y whose coefficients the node stores, so a step costs
5,134 instructions instead of 16,578, and one full multiplication makes the base explicit every sixth step:
0.3432 units per step here instead of 0.5503. Analytic generic-attack bound, not a break.
