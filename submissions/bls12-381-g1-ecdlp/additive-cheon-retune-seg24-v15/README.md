# additive-cheon-retune-seg24-v15 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport; restart segment SEG = 2^24.

Claimed bound: log2 W(2^20) <= 116.368701 (exact rational 116368701/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory
2^121.5854 bytes (every per-node total below 2^63 bytes, all offsets native 64-bit),
preprocessing at most 2^115.3681 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner with its segment/restart
  structure, the cluster model and its one inter-node mechanism, the two-hop wave transport with
  preallocated per-link buffers and global barriers, per-operation instruction listings, the
  measured kernel, the exact ledger, heuristics H1-H5 with their evidence and extrapolations,
  failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `proto`, `census`, `all`.
  The declared experiment runs mode `all` (about 6 seconds locally: ledger, one planner row,
  ec 60 trials, the transport simulation in three patterns, one census row with the multiplicity
  diagnostic). Larger census/ec/planner rows quoted in proof.md sections 7-8 were run offline with
  the same script and inputs; nothing was executed at production size.
- `measure/` — the fused AArch64 side-output kernel (`side_fused.S`, blst's Montgomery routines
  verbatim and renamed, plus the driver), its generator, the Callgrind harness with three
  self-checks, the in-container run script, `reproduce.sh` (organizers' pipeline, Docker with
  linux/arm64) and `measurement.log` (the run behind the charged 1127 instructions per lane).

This is the accepted package `additive-cheon-retune-key17-lanes18-v10` (submission
`0bc0688e-d0ee-449f-90eb-4d1a46f8e712`, official 116.368828) with exactly one changed quantity:
the restart segment is SEG = 2^24 entries instead of 2^25. The planner, padding lemma and segment
directory bound are stated for the new length in proof.md section 3; the padding row falls
(one segment per lane per base is cheaper), the restart row doubles (2N/SEG scalar
multiplications), and the per-wave, barrier and candidate re-walk rows move with SEGLEN/WAVES.
All of these are charged rows recomputed from exact integers by `evidence.py` mode `ledger`.
N, K = 0.4952, f = 0.98 (H4, inherited uncapped evidence at GMAX = 16 d^2, not strengthened),
GMAX = 2^44, ALLOW = 8, LANES = 2^18, WSTEP = 2^34, the 17-instruction key listing, the success
and failure budget, the algorithm, the transport and the toy experiments are inherited unchanged.
The claim is the six-decimal ceiling of the exact recomputed ledger (116.368701). No full-size
execution or new algorithm is claimed. A double-kernel pair pricing (our submission 2cf2cdc2) was
rejected because it omitted the 3(L-1) field multiplications of the batched inversion; this
package does not use it.
