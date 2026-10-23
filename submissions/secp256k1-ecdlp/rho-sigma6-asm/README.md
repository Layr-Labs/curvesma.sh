# secp256k1-ecdlp: parallel Pollard rho on order-6 classes (sigma walk, assembly field arithmetic)

Parallel distinguished-point Pollard rho over the order-6 classes given by negation and the GLV endomorphism of
secp256k1, with batched affine additions on 2^32 walks and 2^16 adders.

The walk stores actual points and adds the image sigma(R_i) of each adder, where sigma maps the canonical point to the
current one; the class sequence is that of the usual canonical walk. Walks run in trails from fresh coins, each ending at
its first distinguished point. The first key match is resolved by replaying the two trails, so a degenerate match has
probability exactly 1/r. A fruitless cycle escapes to twice the canonical point of its least-key member, with no table.
Escape depth is controlled by the exact certificate of section 2.4 (no k below 2^56 has 2^k in Gamma, computed with m = 2^28)
and the depth tail of section 2.3. Long closures are bounded in Lemma F, abandonment in Lemma T'' (section 5).

The common batched step is branch-free and charged at its measured instruction count, 821 per walk step (plain 808.67,
anchor 818.66, plus 2 for the canon_pick branch). Rare work is charged at forced-path maxima, measured with a bounding
build, and a Markov reserve covers its expected amount (7.89 per walk step). The meter aborts with probability at most
1/64, counted in delta.

The submitted bound is **time_log2 = 126.479214**, with success >= 0.39 over the uniform target and the algorithm's coins.
The premises are M1 (fresh classes receive uniform, independent keys and orientations, conditional on any history),
L-ori (orientation law, a measured premise: section 4 and evidence/sign-law.json), and the algorithm's own randomness,
which comes from the operating system.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Exact certificate for escape depths below 2^56](evidence/c-ord.json) and [formal-structure test](evidence/formal-structure.json)
- [Sign-law counts behind section 4](evidence/sign-law.json), rebuilt by [evidence/sign_law.py](evidence/sign_law.py)
- [Per-step instruction measurement](measurement/), results in [evidence/measurement/](evidence/measurement/)
- [Forced-path worst cases](measurement/forced/summary.json) of the rare-path operations
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/trails/](evidence/trails/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is not a practical break.
