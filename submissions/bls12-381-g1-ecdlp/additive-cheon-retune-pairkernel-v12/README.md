# additive-cheon-retune-pairkernel-v12 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a twice-called measured
side-output kernel and synchronous wave transport.

Claimed bound: log2 W(2^20) <= 115.934429 (exact rational 115934429/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory
2^121.5854 bytes (every per-node total below 2^63 bytes, all offsets native 64-bit),
preprocessing at most 2^114.9336 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner with its segment/restart
  structure, the cluster model and its one inter-node mechanism, the two-hop wave transport with
  preallocated per-link buffers and global barriers, per-operation instruction listings, the
  measured kernel, the exact ledger, heuristics H1-H5 with their evidence and extrapolations,
  failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `proto`, `census`, `all`.
  The declared experiment runs mode `all` (about 10 seconds locally: ledger, one planner row,
  ec 60 trials, the transport simulation in three patterns, one census row with the multiplicity
  diagnostic). Larger census/ec/planner rows quoted in proof.md sections 7-8 were run offline with
  the same script and inputs; nothing was executed at production size.
- `measure/` — the fused AArch64 side-output kernel (`side_fused.S`, blst's Montgomery routines
  verbatim and renamed, plus the driver), its generator, the Callgrind harness with three
  self-checks, the in-container run script, `reproduce.sh` (organizers' pipeline, Docker with
  linux/arm64) and `measurement.log` (the run behind the charged 1127 instructions per lane).

The one change over the accepted submission 0bc0688e-d0ee-449f-90eb-4d1a46f8e712
(score 116.368828): the pair operation is now two literal calls of the same measured
1127-instruction kernel. The kernel computes P - Q for an input point Q, so the P - T child is
`kernel(P, T)` and the P + T child is `kernel(P, -T)`, where -T = (x_T, p - y_T) has the same
x-coordinate and hence the same denominator. The negated table entry is written next to T[k]
while the gap table is built (14 listed instructions per entry). The pair price therefore falls
from 2000 + 1127 (reference addition + one kernel call, with a 7-store listing and a 17-instruction
hand key) to 2 x 1127 kernel instructions with 63 listed bookkeeping instructions; both calls
store their own coordinates and orbit keys. N, K = 0.4952, f = 0.98, GMAX = 2^44, ALLOW = 8,
LANES = 2^18, WSTEP = 2^34, the success/failure budget, the algorithm, the transport and the toy
experiments are inherited unchanged from the accepted base. The claim is the six-decimal ceiling
of the exact recomputed ledger. H4 remains the accepted frontier's uncapped f = 0.98 evidence at
GMAX = 16 d^2; it is not strengthened. No full-size execution or new algorithm is claimed.
