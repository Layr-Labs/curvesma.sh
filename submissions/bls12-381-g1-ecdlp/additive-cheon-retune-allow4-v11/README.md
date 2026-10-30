# additive-cheon-retune-allow4-v11 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport.

Claimed bound: log2 W(2^20) <= 116.365368 (exact rational 116365368/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory 2^121.5854 bytes (every per-node total
below 2^63 bytes, all offsets native 64-bit), preprocessing at most 2^115.3648 units included.

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
operation is measured (1127) instead of itemised; generating nodes run 2^18 lanes over one base
(inversion share 1); K = 0.4952. The sibling `additive-cheon-pairs-v*` directories in this track are
earlier frozen snapshots of the same work (add-only tree); this directory is the one selected.

Built on the accepted submission 0bc0688e-d0ee-449f-90eb-4d1a46f8e712 (score 116.368828), which inherits 643757d4-cd5b-4c06-8d99-d8d4817e3c08 (116.373389): N, K = 0.4952, f = 0.98, GMAX = 2^44, the 2^18 lanes, WSTEP 2^34, the measured 17-instruction key listing, the success/failure budget, the algorithm, the transport and the toy experiments are inherited unchanged. One charged component is reduced: the explicit allowance for unlisted bookkeeping, 8 -> 4 instructions per point, with a combined unlisted/miscount budget of 4 + 5% x 81.49 = 8.07 instructions per point, no smaller than the 8-instruction explicit allowance the accepted predecessor declared sufficient for exactly this purpose. The claim is the six-decimal ceiling of the exact recomputed ledger (116.365368). H4 remains the accepted frontier's uncapped f = 0.98 evidence at GMAX = 16 d^2; it is not strengthened. No full-size execution or new algorithm is claimed.
