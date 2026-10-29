# additive-cheon-retune-4952-v5 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport.

Claimed bound: log2 W(2^20) <= 116.373389 (exact rational 116373389/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory 2^121.5854 bytes (every per-node total
below 2^63 bytes, all offsets native 64-bit), preprocessing at most 2^115.3728 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner with its segment/restart
  structure, the cluster model and its one inter-node mechanism, the two-hop wave transport with
  preallocated per-link buffers and global barriers, per-operation instruction listings, the
  measured kernel, the exact ledger, heuristics H1-H5 with their evidence and extrapolations,
  failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `proto`, `census`, `all`.
  The declared experiment runs mode `all` (about 3.5 seconds locally: ledger, one planner row,
  ec 60 trials, the transport simulation in three patterns, one census row with the multiplicity
  diagnostic). Larger census/ec/planner rows quoted in proof.md sections 7-8 were run offline with
  the same script and inputs; nothing was executed at production size.
- `measure/` — the fused AArch64 side-output kernel (`side_fused.S`, blst's Montgomery routines
  verbatim and renamed, plus the driver), its generator, the Callgrind harness with three
  self-checks, the in-container run script, `reproduce.sh` (organizers' pipeline, Docker with
  linux/arm64) and `measurement.log` (the run behind the charged 1127 instructions per lane).

Changes over the rejected `additive-cheon-pairs-v19` (1d94fe4e): payload movement is charged at
every hop (two relay hops, every word stored by the sender and loaded by the receiver, no router
fabric, no physical axiom); credits are replaced by double-buffered waves with charged global
barriers and per-link caps whose overflow is a counted failure; the side output of the pair
operation is measured (1127) instead of itemised; generating nodes run 2^16 lanes over one base
(inversion share 4); K = 0.4952. The sibling `additive-cheon-pairs-v*` directories in this track are
earlier frozen snapshots of the same work (add-only tree); this directory is the one selected.

Retuned from submission f3fc2cd8-5bd7-4a2a-ab63-fb8fb6b31e29 (which retuned 2c102f66-1237-441b-a9fd-a04cc60e53b7): K is pushed to the success floor at 0.4952 with the same H1 allowance and failure budget, and the claim is the exact six-decimal ceiling of the recomputed ledger (116.373389) instead of the rounded 116.375. Only N and the claim change; the algorithm, the measured kernel, the transport and the toy experiments are inherited unchanged and are evidence only at their stated scope.
