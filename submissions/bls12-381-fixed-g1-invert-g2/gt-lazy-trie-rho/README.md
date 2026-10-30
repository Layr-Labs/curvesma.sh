# Frobenius-class rho in GT with lazily evaluated steps (prefix-product trie)

Claim: time_log2 = **124.209** (total work <= 2^124.209 dense-Fp12-multiplication equivalents, C = 20,000),
success >= 0.39 in the ideal model (ChaCha20 assumed indistinguishable from an ideal source) under the
heuristics H1, H2 and H1+ with tolerance t = 2 (stated and tested in proof.md and evidence.md), peak memory
< 2^62 bytes, preprocessing at most 2^49.06 units (included), no advice.

- [Algorithm, success argument and cost ledger](proof.md) (section 0 states the port and the fixes since 680b01db)
- [Implementation details: meter, report log, matcher, restart cap, sampler](proof-impl.md)
- [Experiments and their results, including the registered toy census](evidence.md) ([census driver](census.py)); [executable checks, production lazy walk, toy runs and the ledger](experiment.py)
- `experiment.py` reads its inputs as JSON on stdin and prints one JSON object per part as that part finishes (JSON lines; the last line carries the seconds per part), so a partial run still shows every finished part
- [Kernel generators, emulator and checks](kernels/) (`gen_lazy.py` emits `lazy.s` 4,964 and `root.s` 1,561 instructions)

Port to this track of mpjunior92's gt-lazy-trie-rho (submission a77744bb, 2^124.183 on the mirror track
bls12-381-fixed-g2-invert-g1): the walk only ever sees (g, z), and g = e(P, Q0) for uniform nonidentity P has the
same distribution as g = e(P0, Q) there, so the walk, success analysis and ledger carry over; the answer is
Q = [x]Q0 in G2 (proof.md section 1). That package was built on mpjunior92's gt-dprestart-replay
(submission 0b8665ae, 2^125.111), with one design change: the walk is
evaluated lazily (tag tracing). A walk keeps an explicit base y and a node of a prefix-product trie of the
288 twisted multipliers (R = 24, depth 6); the class key of the current element y * Pi_node needs only the
w^1 coefficient of the product, an Fp2-linear form in y whose coefficients the node stores, so a step costs
5,134 instructions instead of the parent's 14,402 + 2,014 = 16,416 (its README quotes 16,578 for the same step), and one full multiplication makes the base explicit every sixth step:
0.3485 units per step against the 0.5503 the parent's README states for its own step mix. Neither parent figure
enters this ledger. Analytic generic-attack bound, not a break.
