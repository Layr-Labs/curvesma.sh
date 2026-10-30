# additive-cheon-retune-rel32-v23 (selected complete candidate)

Cheon's attack on bls12-381-g1-ecdlp with planned pair walks, a measured side-output kernel and
synchronous wave transport. This revision re-counts the charged kernel against the measured
Callgrind price and charges the itemised listings exactly.

Claimed bound: log2 W(2^20) <= 116.357427 (exact rational 116357427/1000000,
exact total 212865630936303925449710601903919952394 instructions),
heuristic success at least 0.3900 (K = 0.4952, computed floor 0.390050), peak memory
2^121.5854 bytes (every per-node total below 2^63 bytes; the planner pointer fields are
32-bit relative offsets and all other offsets are native 64-bit),
preprocessing at most 2^115.3608 units included (the baby half is charged at its own
45-instruction transport listing, not the 40.5 two-phase average; the combined total is unchanged).

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
  linux/arm64) and `measurement.log` (the run behind the measured 1127 = 1116 call-only + 11 C
  lane loop instructions per lane).

This is the package of our accepted submission `d421a23d-dcc9-4896-adea-06bbd1947c57`
(official 116.365049, candidate `additive-cheon-retune-rel32-v22`), which is the accepted
`97f53a10-2f83-409f-a654-040f283fb31f` with the planner rel32 repair, the pair lane-buffer
12 -> 11 recount, the root table-address 2 -> 0 recount and the baby-side preprocessing fix.
The algorithm, the transport, the planner, H1-H5, f = 0.98, K = 0.4952 and every unlisted row are
inherited unchanged from that accepted base. This revision changes two ledger items:

1. **Kernel charged call-only (1127 -> 1116) and the lane loop charged once.** The Callgrind
   price of the counted harness loop is 1127 per lane: 181 driver + 2 x 341 + 99 + 154 = 1116
   for the straight-line `side_fused` call, plus the 11-instruction C lane loop
   (`mov x3; mov x2; mov x1; mov x0; add; bl; add; add; add; cmp; b.ne`) that the harness
   executes once per lane. The previous revision charged KERNEL = 1127 (which already contains
   that loop) *and* a separate 10-instruction "lane control and loop" item, i.e. it charged the
   loop twice. The corrected ledger charges KERNEL = 1116 and lane control = 11: together 1127,
   exactly the measured price, with no double count and no undercount. The measured total
   (4616213 / 9232405 -> 1127) and every self-check are unchanged.
2. **Listed counts charged exactly (margin 21/20 removed).** The listings of section 5.2 are
   instruction-level counts of the implementation; the 21/20 factor was a 5% cushion on those
   listed counts, not a measured quantity. Removing it lowers the pair listing charge from
   1.05 x 74 to 74 and the root/record listings similarly. The separate explicit
   "8 instructions per point" allowance row is unchanged, still charged for both sides, and it
   exceeds the removed margin (the removed 5% is 3.9 instructions per point: 0.98/2 x 3.7 +
   0.02 x 3 + 0.05 x 40.5), so the claimed total keeps a strictly larger unallocated allowance
   than the margin it drops.

Consequences: per-point work falls 1653.5735 -> 1644.8 (exactly 8224/5); the exact total falls by
1127633132802301512254138618906798018 instructions (626603898686230592662282751351758416 from
the kernel/lane-control recount at the old margin, 501029234116070919591855867555039602 from the
margin removal); log2 changes from 116.3650487062493326 to 116.3574263536077..., and the claim is
its six-decimal ceiling 116.357427. N, K = 0.4952, f = 0.98, GMAX = 2^44, ALLOW = 8,
LANES = 2^18, WSTEP = 2^34, SEG = 2^24, the 17-instruction key, the success/failure budget, the
algorithm, the transport and the toy experiments are otherwise inherited unchanged. Total memory
stays 2^121.5854 bytes. No full-size execution or new algorithm is claimed.

Checks run on this package: `python3 evidence.py` mode `all` (ledger PASS, planner/ec/proto/census
toy evidence), the preflight, and the bounded Docker experiment (exit 0, ledger score 116.357427).
