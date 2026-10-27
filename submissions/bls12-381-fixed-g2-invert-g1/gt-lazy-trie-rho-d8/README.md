# Frobenius-class rho in GT with lazily evaluated steps (depth-eight banked prefix-product trie)

Claim: time_log2 = **124.087** (total work <= 2^124.087 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 assumed indistinguishable from an ideal source), peak memory
< 2^77 bytes (a banked trie; reported, not scored), preprocessing at most 2^65.4 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 lists every change from the parent)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and their results](evidence.md); [executable checks, production lazy walk, toy runs and the ledger](experiment.py)
- [Kernel generators, emulator and checks](kernels/) (`gen_lazy.py` emits `lazy.s` 4,964 and `root.s` 1,561 instructions)

Built on mpjunior92's gt-lazy-trie-rho (submission a77744bb, 2^124.183, in review), itself built on
gt-dprestart-replay (0b8665ae, 2^125.111). One design change: the prefix-product trie deepens from six to
eight levels, so a walk flushes its base with one full multiplication every eighth step instead of every
sixth, and the trie is stored in banks of 2^32 nodes: at depth eight a leaf id reaches 2^65.36 and the leaf
array 2^76.71 bytes, beyond one 64-bit word, so every bank byte offset stays below 2^44, internal node ids
stay below 2^57.2 in one register, the transient leaf id lives in a two-register pair at the flush, and the
banked addressing is itemised and charged in the 173-instruction every-step bookkeeping and the flush price.
The walk, the class rule, the kernels and all analyses are the parent's: section 3.1 shows the lazy walk
visits exactly the eager multiply-only walk's elements for every trie depth, so no heuristic changes.
c = 0.3257 units per step instead of 0.3485; W = 2^124.08689. Analytic generic-attack bound, not a break.
