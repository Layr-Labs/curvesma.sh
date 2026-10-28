# Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **123.947** (total work <= 2^123.947 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 assumed indistinguishable from an ideal source), peak memory
< 2^62 bytes, preprocessing at most 2^50.04 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 lists every change from the parent)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and their results](evidence.md); [executable checks, production lazy walk, toy runs and the ledger](experiment.py)
- [Kernel generators, emulator and checks](kernels/) (`gen_lazy.py` emits `lazy.s` 3,527 and `root.s` 3,425 instructions)

Built on mpjunior92's gt-dprestart-replay (submission 0b8665ae, 2^125.111). One design change: the walk is
evaluated lazily (tag tracing). A walk keeps an explicit base y and a node of a prefix-product trie of the
288 twisted multipliers (R = 24, depth 6); the class key of the current element y * Pi_node needs only the
w^1 coefficient of the product, an Fp2-linear form in y whose coefficients the node stores, so a step costs
3,697 instructions instead of 16,578, and one full multiplication makes the base explicit every sixth step:
0.2906 units per step instead of 0.5503. An earlier revision of gt-lazy-trie-rho (a77744bb, 2^124.183) changed
the key kernel to pair the six terms of that linear form (Winograd), 9 products instead of 18, and a later
revision (2^124.012) rebuilt its success argument. The multiplication-kernel revision (2^123.947)
puts the fixed operand on the table side at all three tower levels with its values stored reduced, 12,774
instructions instead of 14,402; an internal trie node keeps its 12 coordinates and is charged two units (22,040
instructions). This revision keeps that package and changes one thing, its success argument (below).
Analytic generic-attack bound, not a break.

The success argument charges no allowance: every loss is a computed bound that needs no independence between
iterations, trails or reports (proof.md sections 3.4 and 4, computed by `bounds.py`). Expectations come from
deferred decisions under H2 (a fresh class's labels are uniform given everything revealed), H1 for merges, and linearity; each cap
is charged its Markov share (the report cap 617 times its expectation, the restart cap 986 times, the pair budget
2,048 times); the restart cap's worst case is also charged as steps spent on restarts; and the closed-walk waste is
bounded by Freedman's martingale inequality (tail below 2^-1000). Section 4 floors the success after every loss at
0.39689.
