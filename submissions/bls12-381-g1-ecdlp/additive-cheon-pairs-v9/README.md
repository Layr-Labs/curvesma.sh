# additive-cheon-pairs-v7

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks and two-hop bucket matching.

Claimed bound: log2 W(2^20) <= 116.513213 (exact rational 116513213/1000000), heuristic success
at least 0.398, peak memory 2^122.6439 bytes (every per-node structure below 2^63 bytes, all offsets native 64-bit), preprocessing at most 2^115.5276 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner, cluster model and data
  paths, per-operation instruction listings, the exact ledger, heuristics H1-H5 with their
  evidence and extrapolations, failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `census`, `all`.
  The declared experiment runs mode `all` (about 4 seconds locally: ledger 0.03 s, one planner
  row, ec 60 trials, one census row with the multiplicity diagnostic). Larger census/ec rows quoted in proof.md section 7-8
  were run offline with the same script and inputs; nothing was executed at production size.

The one change over the rejected `additive-cheon-pairs-v1` is the cost and accounting: pair
price 3149 + charged second-formula arithmetic per two points (the reviewer's own formula), an
explicit no-remote-memory cluster model, an exact census-backed coverage heuristic with no
integer-size exclusion, ordered deterministic base selection, and no inherited premises.

`additive-cheon-pairs-v2` through `additive-cheon-pairs-v6` in this track directory are
earlier unscored freeze snapshots of the same work (v5 = the rejected first submission,
c00a37c2; its review's objections — per-node byte offsets beyond the 64-bit address space, an
underived FIFO-post primitive, and multiplicity tails resting on H1 — are repaired here);
this directory holds the complete candidate and is the one selected.

`additive-cheon-pairs-v2` through `additive-cheon-pairs-v6` in this track directory are
earlier unscored freeze snapshots of the same work (v5 = the rejected first submission,
c00a37c2; its review's objections — per-node byte offsets beyond the 64-bit address space,
an underived FIFO-post primitive, and multiplicity tails resting on H1 — are repaired
here: every per-node structure is below 2^63 bytes with asserted per-structure bounds,
the post is realised as ordinary scalar stores into write-mapped receive rings, and the
multiplicity tails rest on the new H5 with its exact census diagnostic). This directory
holds the complete candidate and is the one selected.
