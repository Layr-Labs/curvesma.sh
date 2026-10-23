# secp256k1-ecdlp: parallel Pollard rho on order-6 classes (sigma walk, assembly field arithmetic)

Parallel distinguished-point Pollard rho over the order-6 classes given by negation and the GLV endomorphism of secp256k1. Steps are batched affine additions on 2^32 walks
with 2^16 adders (2^30 for BLS12-381 G1, 16384 walks per batch).

The walk stores actual points and adds the image sigma(R_i) of each adder, where sigma maps the canonical point to
the current one. The class sequence is that of the usual canonical walk. Which rotation is canonical is decided from
64-bit fixed-point fractions, with an exact fallback.

Walks run in trails from fresh coins. Each trail ends at its first distinguished point and records it once with its
coins, including a distinguished point first reached by a look-ahead. The first key match is resolved by replaying the
two trails, so a degenerate match has probability exactly 1/r. A fruitless cycle is escaped to a point of a separate
escape table, chosen by the cycle's key, so the point after an escape is fresh given the cycle, and a trail cannot be
trapped by its own escapes (Lemma T). The toy experiment runs this whole procedure, including resolution, and every
trail is replayed.

The common batched step is branch-free and charged at its measured instruction count, 817 per walk step (plain steps
806.7, anchor steps 816.7). Rare work is charged at forced-path maxima, measured with a bounding build of
each unbatched operation, and a Markov reserve covers its expected amount (7.89 per walk step).
The meter aborts with probability at most 1/64, counted in delta.

The submitted bound is **time_log2 = 126.472235**, with success >= 0.39 under the disclosed random-mapping heuristic H1
and the fresh-class premise H2.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Per-step instruction measurement](measurement/), results in [evidence/measurement/](evidence/measurement/)
- [Forced-path worst cases](measurement/forced/summary.json) of the rare-path operations
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/trails/](evidence/trails/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is not a practical break.
