# unroll16: sixteen-copy blocks on the three-arm anchor schedule

P-256 sigma rho, declared **126.731865** (25346373/200000 at C=1000) on the three-arm anchor-schedule package
(126.732475). One change: both per-pass blocks hold sixteen copies instead of eight, with a four-bit `tbz`/`tbnz`
entry tree on `(x23 >> 6) mod 16` and both end tests counting sixteen steps. Every copy body, tail, slow block and
rare path is byte-identical, so the walk, correction meter, distinguished predicate, caps and success argument are
unchanged. Saving: B/4 fewer loop-control instructions per round, less a measured +2 on the dispatch's common path.
PLAIN_FIXED 605,581,542 -> 605,319,400; anchor arm 605,320,427; skip arm 603,223,276; C_COMMON 576.327740; CORR
14.210693; reserve 0.0075889; total step 590.546021. Paired organizer-image counts, all three arms at 2^16 on both seeds and the 2^20 pair on seed 1:
every round falls by exactly 16,382 at B=A=2^16 and 262,142 at B=A=2^20, with identical meters and byte-identical
dumps, and the pair is also identical at production B=A=2^20 with 16,411 real distinguished points. See proof.md and
evidence/anchor-schedule.json.
