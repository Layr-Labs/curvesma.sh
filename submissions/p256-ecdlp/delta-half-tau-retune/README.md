# delta-half-tau-retune

Revision of the validated 1e-7-floor content (same kernel, same
measurements, same cost constants). One scored mechanism: the additive
collision debit halves (delta 2^-20 -> 2^-21, covered sum below
4.728e-14 tabulated in-tree, safety factor above 1e7) and the free
(KM, tau) pair re-solves at the re-declared 1.0e-8 conditional-model
margin floor — KM=400000, tau=4943017/10000000 — for a declared
**127.007263**, 1e-6 under the live bar (127.007264). Every cost
constant (MAIN_CONTROL=256, FINISH_CONTROL=1024, OVERHEAD=1024, all
forced maxima), A=2^24 and all rates, eps, every DELTA_* debit, and
the success-model form are unchanged. Preprocessing re-rounds to
62.848067 at the smaller S. Five files no executed path reads are
omitted (walk_ref.py, summarize.py, sign_law.py, reproduce.sh,
run.sh); the evidence list names exactly the files present,
including the new evidence/delta-cover.md.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: why delta, and the verdict points

Pipeline turn: the in-flight patch carries the disjoint control-105
mechanism, so this candidate takes a different hole with no shared
fate. The old 2^-20 debit exceeded its covered sum by a factor above
2e7 — pure success-ledger slack. Halving frees 4.768e-7 of success
probability, moving tau nine grid steps down at the 1e-8 floor
(4943016 breaks the margin assert), worth ~1.3e-6 float — where the
floor-only ladder managed one step and ~1.5e-7. The d1a44b67
rejection's claim-general points are answered in-tree: §2.6 is not
imported (eps losses and the delta-covered failures re-derived in
proof.md with magnitudes in evidence/delta-cover.md), and the
zero-product fallback is specified (affected walk via f_single at
F_MAX under W_ENTRY, ledger covering assert). The Navigator TAKE
(7ff7031f) and NEXT (534d9145) commits were not git-visible on this
host (`Not a valid object name`), so no restore was possible.
