# additive-cheon-retune-999-capscale-v9 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport.

Claimed bound: log2 W(2^20) <= 116.369361 (exact rational 116369361/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory 2^121.5854 bytes (every per-node total
below 2^63 bytes, all offsets native 64-bit), preprocessing at most 2^115.3688 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner with its segment/restart
  structure, the cluster model and its one inter-node mechanism, the two-hop wave transport with
  preallocated per-link buffers and global barriers, per-operation instruction listings, the
  measured kernel, the exact ledger, heuristics H1-H5 with their evidence and extrapolations,
  failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `proto`, `census`,
  `all`. The declared experiment runs mode `all` in about 34 s: ledger; one no-cap planner spec
  at d=16, logX2=17 (span/2GMAX = 8, 4 classes, seed 0); three cap-active planner specs at
  d=8,16,32 with `max_tree = d` (window/tree-cap = 64 as at production; the size cap binds and
  rejects 15544-29676 merges); the same-list no-cap control at d=16; ec 60 trials; the transport
  simulation in three patterns; one census row with the multiplicity diagnostic. The d=64/128
  and multi-seed planner rows quoted in proof.md section 7 were run offline with the same script
  and the inputs listed there; nothing was executed at production size.
- `measure/` — the fused AArch64 side-output kernel (`side_fused.S`, blst's Montgomery routines
  verbatim and renamed, plus the driver), its generator, the Callgrind harness with three
  self-checks, the in-container run script, `reproduce.sh` (organizers' pipeline, Docker with
  linux/arm64) and `measurement.log` (the run behind the charged 1127 instructions per lane).

Changes over the pending retune 0bf85e4a-a77d-45f5-ba86-d29ca0c3e5be (116.369361), which
retuned the accepted frontier 643757d4-cd5b-4c06-8d99-d8d4817e3c08 (116.373389): the H4 evidence
only. The planner rows now include cap-active runs with the tree cap set to the class density d,
so the production window/tree-cap ratio 2*GMAX/d^2 = 64 is exercised while the size cap binds
(15544-29676 rejections); against the same-list no-cap control at d=16 the cap excess is
0.1631, 0.0828, 0.0414, 0.0227 and 0.0138 at d = 8, 16, 32, 64 and 128, i.e. B/d with B =
1.30-1.77, bounding the production cap excess by 2.1e-6. The no-cap span-8 minimum is 0.9901,
so f >= 0.9901 - 2.1e-6 > 0.9900. GMAX = 2^45 = 32 d^2, f = 0.9900, the 2^18 lanes, WSTEP =
2^34, N, K and the whole ledger are unchanged; the claimed score 116.369361 is the six-decimal
ceiling of the exact recomputed ledger, which is the v8 ledger.

Retuned from submission f3fc2cd8-5bd7-4a2a-ab63-fb8fb6b31e29 (which retuned
2c102f66-1237-441b-a9fd-a04cc60e53b7): K is at the success floor 0.4952 with the same H1 allowance
and failure budget. The algorithm, the measured kernel, the transport and the toy experiments are
inherited; the new planner evidence is evidence at toy scope, and H4 remains a disclosed
extrapolation to the production regime.
