# additive-cheon-retune-capmatch-v17 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport; per-record transport listing re-counted (average 43.5 -> 40.5 listed
instructions per record), SEG = 2^24 as in the direct base, and matched capped/uncapped planner
controls at the production window/tree-cap ratio as H4 evidence.

Claimed bound: log2 W(2^20) <= 116.365977 (exact rational 116365977/1000000),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory
2^121.5854 bytes (every per-node total below 2^63 bytes, all offsets native 64-bit),
preprocessing at most 2^115.3654 units included.

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

This is the package of our submission `676b4898-08d2-4b06-9abd-5c05f0b862b1` (claimed 116.368701),
which is the accepted package `additive-cheon-retune-key17-lanes18-v10` (submission
`0bc0688e-d0ee-449f-90eb-4d1a46f8e712`, official 116.368828) with the restart segment
SEG = 2^24, with two changes and no ledger number:

1. **Per-record transport listing.** The generating-node, relay, owner-insert and owner-lookup
   counts become 13, 16, 16 and 7 (was 14, 17, 17, 8), average 40.5 instead of 43.5 listed per
   record. The savings are data-layout/loop-form only: (i) a link's window base and fill counter
   are adjacent 64-bit fields of the same 24-byte window-table entry, so they are one `ldp` (not
   two `ldr`) at the generating node and the relay; (ii) the owner receive loop already knows its
   record count from the per-link count word, and the five-word record load is post-indexed, so
   its loop control is `subs` + `b.ne` (2) instead of pointer advance + compare + branch (3). No
   word, cap check or counter update is removed; each hop still loads five words and stores five
   words. The records row falls by exactly 2N x 3.15 = 404733866663009901405710536879753935
   instructions; the candidate/reply row falls by 2^47 = 140737488355328 with REC_RELAY 17 -> 16;
   together the exact total falls by 404733866663009901405851274368109263 and every other row is
   unchanged and recomputed by `evidence.py` mode `ledger`.
2. **H4 evidence.** Matched capped/uncapped planner controls at GMAX = 16 d^2 on the *same* real
   chi-class lists for d = 8, 16, 32, 64, 128, with `max_tree = d` so that the window/cap ratio is
   32, exactly the production ratio (2 GMAX/d^2 = 32). The cap binds in every capped row (largest
   tree = d-1, 2135-4635 size rejections, zero in the uncapped controls), and the measured cap
   excess decays as B/d with B = 1.31-1.65. Inflating B by 33% to 2.2 bounds the production cap
   excess by 2.2/2^20 = 2.1e-6, while the same-list uncapped fraction rises with d to 0.9879 at
   d = 128 (production d = 2^20 and span 2^31 are far beyond), so f_prod >= 0.9879 - 2.1e-6 >
   0.9878 > 0.98. The declared experiment includes the d=8 and d=16 capped/uncapped pairs; the
   d=32..128 rows are offline with the same code and `max_tree` spec field. H4 remains a disclosed
   extrapolation, now with cap-sensitive evidence where the cap binds.

N, K = 0.4952, f = 0.98, GMAX = 2^44, ALLOW = 8, LANES = 2^18, WSTEP = 2^34, SEG = 2^24, the
17-instruction key listing, the success and failure budget, the algorithm, the transport and the
toy experiments are otherwise inherited unchanged. The claim is the six-decimal ceiling of the
exact recomputed ledger (116.365977). No full-size execution or new algorithm is claimed. A
double-kernel pair pricing (our submission 2cf2cdc2) was rejected because it omitted the 3(L-1)
field multiplications of the batched inversion; this package does not use it.
