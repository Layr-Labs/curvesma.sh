# additive-cheon-retune-990-lanes18-v8 (selected complete candidate)

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
- `evidence.py` — standard library only and unchanged in structure; modes `ledger`, `planner`,
  `ec`, `proto`, `census`, `all`. The declared experiment runs mode `all` in under 20 s: ledger,
  one planner spec at d=16, logX2=17 (span/2GMAX = 8, 8 classes, seed 0), ec 60 trials, the
  transport simulation in three patterns, one census row with the multiplicity diagnostic. The
  span-16/32 and multi-seed planner rows quoted in proof.md section 7 were run offline with the
  same script and the inputs listed there; nothing was executed at production size.
- `measure/` — the fused AArch64 side-output kernel (`side_fused.S`, blst's Montgomery routines
  verbatim and renamed, plus the driver), its generator, the Callgrind harness with three
  self-checks, the in-container run script, `reproduce.sh` (organizers' pipeline, Docker with
  linux/arm64) and `measurement.log` (the run behind the charged 1127 instructions per lane).

Changes over the accepted frontier submission 643757d4-cd5b-4c06-8d99-d8d4817e3c08
(score 116.373389): the pair-gap bound is raised from GMAX = 2^44 = 16 d^2 to GMAX = 2^45 = 32 d^2;
the claimed minimum pair fraction from f = 0.98 to f = 0.9900, supported by planner rows at
GMAX/d^2 = 32 (span/2GMAX = 8: 32 class fractions over seeds 0-3 between 0.9901 and 0.9927; span 16 seed 0: min 0.9921); and the
generating-node lane count from 2^16 to 2^18 with the wave length from 2^36 to 2^34 steps, which
divides the per-lane batched-inversion share by four while holding the records per wave at 2^53 (so
every cap, wave-load bound and transport argument is unchanged; the wave count moves from 2170982
to 2171080). The plan entry widens its gap field to 46 bits, and the records, gap-table, planner,
candidate re-walk, inversion-share, padding, strand-start and wave rows are recomputed for the new
GMAX, f and lane count; N and K are unchanged. The claimed score is the six-decimal ceiling of the
exact recomputed ledger.

Retuned from submission f3fc2cd8-5bd7-4a2a-ab63-fb8fb6b31e29 (which retuned
2c102f66-1237-441b-a9fd-a04cc60e53b7): K is at the success floor 0.4952 with the same H1 allowance
and failure budget. The algorithm, the measured kernel, the transport and the toy experiments are
inherited; the new planner evidence is evidence at toy scope, and H4 remains a disclosed
extrapolation to the production regime.
