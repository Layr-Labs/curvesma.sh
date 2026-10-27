# Cheon collision search with local storage

This candidate recovers the full BLS12-381 scalar from A_1 and A_(2^20).
Its bounded two-kind walk uses fixed-base point tables and exact distinguished
orbit labels. All tables, records and scratch fit below 2^44 local address bytes.
The submitted total-work bound is 2^123 target-group-operation equivalents,
with heuristic success at least 0.40 under the explicit assumptions in proof.md.

The proof includes cold table construction, scalar-word and hash kernels,
all local memory operations, caps, failed trials, root recovery and output.
Experiments check curve recovery, scalar walk statistics, production curve
identities, Montgomery kernels, SHA-256 and address bounds at their stated scope.
They do not execute a production search or prove its heuristic probability.

Run repository preflight and declared experiments with the pinned verifier as
described in the root README. Existing submissions and protected parameters are
unchanged. The first manifest claim supplies the score; external review decides
whether the full claim qualifies.
