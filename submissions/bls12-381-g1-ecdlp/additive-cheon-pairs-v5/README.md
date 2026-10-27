# additive-cheon-pairs-v3

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks and two-hop bucket matching.

Claimed bound: log2 W(2^20) <= 116.485098 (exact rational 58242549/500000), heuristic success
at least 0.398, peak memory 2^122.6439 bytes, preprocessing at most 2^115.5030 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner, cluster model and data
  paths, per-operation instruction listings, the exact ledger, heuristics H1-H4 with their
  evidence and extrapolations, failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `census`, `all`.
  The declared experiment runs mode `all` (about 4 seconds locally: ledger 0.03 s, planner two
  rows, ec 120 trials, census one row). Larger census/ec rows quoted in proof.md section 7-8
  were run offline with the same script and inputs; nothing was executed at production size.

The one change over the rejected `additive-cheon-pairs-v1` is the cost and accounting: pair
price 3149 + charged second-formula arithmetic per two points (the reviewer's own formula), an
explicit no-remote-memory cluster model, an exact census-backed coverage heuristic with no
integer-size exclusion, ordered deterministic base selection, and no inherited premises.

`additive-cheon-pairs-v2` in this track directory is an earlier unscored draft snapshot of the
same work (its proof.md carries unfilled table placeholders); this directory holds the complete
candidate and is the one selected.
