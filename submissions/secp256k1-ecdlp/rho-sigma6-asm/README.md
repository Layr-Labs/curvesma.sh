# secp256k1-ecdlp: parallel Pollard rho on order-6 classes (sigma walk, assembly field arithmetic)

Parallel distinguished-point Pollard rho over the order-6 classes given by negation and the GLV endomorphism of secp256k1. Steps are batched affine additions on 2^32 walks
with 2^16 adders (2^30 for BLS12-381 G1, 16384 walks per batch).

The walk stores actual points and adds the image sigma(R_i) of each adder, where sigma maps the canonical point to
the current one. The class sequence is that of the usual canonical walk. Which rotation is canonical is decided from
64-bit fixed-point fractions, with an exact fallback.

Walks run in trails from fresh coins. Each trail ends at its first distinguished point and records it once with its
coins, including a distinguished point first reached by a look-ahead. The first key match is resolved by replaying the
two trails, so a degenerate match has probability exactly 1/r. A fruitless cycle is escaped to twice the canonical point of its least-key member, a
doubling with no table. Escape depth and escape-point repeats are controlled by the formal-structure argument of proof.md section 2
(D1, from an exact certificate that no k below 2^50 has 2^k in Gamma, and a depth tail), and a trail cannot be trapped by its own escapes (Lemma T', section 5). Within-trail
coincidences, including re-entry into an escaped cycle, are bounded in section 3. The toy experiment runs this whole procedure, including resolution, and every
trail is replayed.

The common batched step is branch-free and charged at its measured instruction count, 821 per walk step (plain steps
808.67, anchor steps 818.66 measured, plus 2 for the canon_pick branch). Rare work is charged at forced-path maxima, measured with a bounding build of
each unbatched operation, and a Markov reserve covers its expected amount (7.89 per walk step), including the exact-k fallback of the batch path.
The meter aborts with probability at most 1/64, counted in delta.

The submitted bound is **time_log2 = 126.479214**, with success >= 0.39 over the uniform target and the algorithm's coins, under the disclosed random-mapping heuristic H1,
the class-data heuristic M1 and the algorithm's own randomness (the table and every trail's coins, drawn from the operating system), and the exact certificate of section 2.4 for escape depths below 2^50.

- [Algorithm, success argument and cost](proof.md)
- [Exact accounting script](accounting.py) and its output in [evidence/accounting.json](evidence/accounting.json)
- [Formal-structure test on small builds](evidence/formal-structure.json) and [the exact certificate for escape depths below 2^50](evidence/c-ord.json)
- [Per-step instruction measurement](measurement/), results in [evidence/measurement/](evidence/measurement/)
- [Forced-path worst cases](measurement/forced/summary.json) of the rare-path operations
- [Declared toy-curve experiment](experiment.py); offline runs in [evidence/trails/](evidence/trails/)
- [Claim summary](claim.json)

This applies a known generic square-root attack and accounts for it carefully. It is not a practical break.
