# anchor-h1-complete

Resubmission of the measured-anchor charge (rejected `625051c8`) with the
success-argument gap closed. Two parts, one shippable claim:

1. Scored delta (unchanged arithmetic): the anchor batch round is charged
   at its seed-1 measurement (47,242,860 per 2^16 walks) instead of the
   derived 47,407,420, dropping the 164,560 derivation margin. Declared
   **127.007187**, 77e-6 under the live bar (127.007264).
2. Verdict response (new support): `evidence/h1-success.md` reproduces the
   precise H1 statement for this adaptive cycle-escaping multi-trail
   process, derives the eps margins from in-tree rate lemmas, decomposes
   the 2^-20 allowance into named bounded events, and assigns the DP_MAX
   overflow failure its Chernoff bound. No success debit changes value;
   the ledger's eps/delta stand as conservative envelopes of derived
   bounds, with new asserts checking the envelope relations.

Base: the validated `joint-retune-rival` crown content (127.007267).
KM=400000, tau=4943041/10000000, 1.0e-6 floor, kernel, measurements, and
all other prices unchanged. Five files no executed path reads are
omitted (walk_ref.py, summarize.py, sign_law.py, reproduce.sh, run.sh);
the evidence list names exactly the files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: verdict response

Pipeline turn: the in-flight patch carries the 1e-7-floor retune, which
shares the same inherited success text and will likely draw the same
rejection. This candidate therefore answers the `625051c8` verdict
directly instead of stacking another cost delta: the H1 statement, the
eps-margin derivations, and the allowance decomposition are now in-tree
with ledger asserts, so the reviewer can assess named premises rather
than missing sections. The Navigator TAKE/NEXT commits were not
git-visible on this host, so no restore was possible.
