# additive-cheon-retune-rel32-v20 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport; per-record transport listing re-counted (average 43.5 -> 40.5 listed
instructions per record), SEG = 2^24 as in the direct base, matched capped/uncapped planner
controls at the production window/tree-cap ratio as H4 evidence, an explicit per-element planner
state whose parent, partner, child-head and next fields are signed 32-bit offsets from the
element's own index (so the 32-byte record is representable), the repeated descendant traversals
charged in their own ledger row, and the pair midpoint load corrected from 8 to 7 instructions
(one `madd` address with a hoisted 96-byte stride).

Claimed bound: log2 W(2^20) <= 116.365531 (exact rational 116365531/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory
2^121.5854 bytes (every per-node total below 2^63 bytes; the planner pointer fields are
32-bit relative offsets and all other offsets are native 64-bit),
preprocessing at most 2^115.3650 units included.

- `proof.md` — the complete argument: reduction, algorithm, planner with its segment/restart
  structure, the cluster model and its one inter-node mechanism, the two-hop wave transport with
  preallocated per-link buffers and global barriers, per-operation instruction listings, the
  measured kernel, the exact ledger, heuristics H1-H5 with their evidence, the matched cap-active
  H4 controls, failure budget, limitations.
- `evidence.py` — standard library only; modes `ledger`, `planner`, `ec`, `proto`, `census`, `all`.
  The declared experiment runs mode `all` (about 11 seconds locally: ledger, five planner rows
  (the uncapped d=16 row plus the d=8 and d=16 capped/uncapped pairs with the `max_tree` spec),
  ec 60 trials, the transport simulation in three patterns, one census row with the multiplicity
  diagnostic). Larger census/ec/planner rows quoted in proof.md sections 7-8 were run offline with
  the same script and inputs; nothing was executed at production size.
- `measure/` — the fused AArch64 side-output kernel (`side_fused.S`, blst's Montgomery routines
  verbatim and renamed, plus the driver), its generator, the Callgrind harness with three
  self-checks, the in-container run script, `reproduce.sh` (organizers' pipeline, Docker with
  linux/arm64) and `measurement.log` (the run behind the charged 1127 instructions per lane).

This is the package of our accepted submission `97f53a10-2f83-409f-a654-040f283fb31f`
(official 116.365977), which is the accepted `676b4898-08d2-4b06-9abd-5c05f0b862b1` with the
per-record transport recount (13/16/16/7, average 43.5 -> 40.5) and the matched cap-active H4
controls; the H4 evidence and f = 0.98 are inherited unchanged. Relative to that accepted base
this revision changes three ledger rows:

1. **Pair midpoint load.** The midpoint P load is 6 ldp plus one `madd` address (slot x 96 + the
   hoisted lane-stack base) = 7 instructions, not the 8 charged in the base (a second address
   instruction); the T[g] load keeps its 12 (the gap table is built for 2 <= k <= GMAX, so a
   g = 1 access needs its base-point select). The pair listing is 74 (was 75) and the root
   listing 62 (unchanged), so per-point work falls 1654.6445 -> 1654.13 (exactly 165413/100) and
   the records row falls by 66106531554958283896266054357026477 instructions; the padding row
   (PADOPS x PAIR) falls by 17497954643939450004989755470643.
2. **Planner state and accounting.** The planner's per-element state is specified (parent, partner,
   depth, height, pair-slot, tree size = 24 bytes; child lists as adjacency arrays head/next =
   8 bytes; 32 bytes per element, charged in the planning-node memory row, 2^62.59 B per node,
   still below 2^63). The four pointer fields are signed 32-bit offsets from the element's own
   absolute class-list index, with -2^31 as the null sentinel: parent[i] = i_parent - i,
   partner[i] = i_partner - i, head[k] = i_first - k and next[a] = i_next - a. Every target lies
   in the 2^26-member planner window (the H2 window cap), so |offset| < 2^27 < 2^31; the
   executable `plan()` stores and decodes exactly these fields and its `rel32()` asserts the bound
   on every store. The repeated `set_depth`/height traversals over those child lists are charged by
   an explicit worst-case row: at most 4 l max_tree subtree visits per class at 16 instructions
   plus l DEPTH height increments at 8, i.e. 64 max_tree + 8 DEPTH = 2^26.03 instructions per
   element, adding 10161011842454401296372011630592 instructions, with the pointer decodes and
   encodes inside those per-visit and per-emission counts. The 8/point allowance is not relied
   upon for it.
3. Nothing else changes. The exact total is 214064784550171335983996421040764310916 instructions
   (log2 116.3655308008...); the claim is its six-decimal ceiling 116.365531. N, K = 0.4952,
   f = 0.98, GMAX = 2^44, ALLOW = 8, LANES = 2^18, WSTEP = 2^34, SEG = 2^24, the measured 1127
   kernel, the 17-instruction key, the success/failure budget, the algorithm, the transport and
   the toy experiments are otherwise inherited unchanged. Total memory stays 2^121.5854 bytes.
   No full-size execution or new algorithm is claimed. A double-kernel pair pricing (our
   submission 2cf2cdc2) was rejected because it omitted the 3(L-1) field multiplications of the
   batched inversion; this package does not use it.
