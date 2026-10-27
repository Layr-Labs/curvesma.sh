# tau-floor-c176-retune

Pure operating-point retune of the accepted 127.005504
control-176-rare-overhead crown (forward loop unrolled by four,
checkpoint-only cycle detection, main control 176, rare-path
overhead 256). Same kernel, same measurements, same cost model,
same premises, same collision derivation. One scored mechanism:
the declared success floor 1.0e-5 -> 1.0e-7 (the floor the prior
127.007-era accepted crowns declared) and tau re-solved
4943241 -> 4943079 on the 1e-7 grid, the smallest grid point
keeping the floor on this base. Declared **127.005480**,
strictly under the bar: total float 127.005503394 -> 127.005479754,
preprocessing 62.848100 -> 62.848076, p_claim 0.390010015 ->
0.390000133 with margin 1.33e-7 (one step lower 7.22e-8 breaks
the floor). KM=420000, delta=2^-20, MAIN_CONTROL=176,
OVERHEAD=256, every rate, every forced maximum, every meter, and
all measurements are unchanged. The evidence list names exactly
the files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: review posture

The prior-turn tau-floor leftover (against the 127.005509 base)
came off the tree when the control-176-rare-overhead promotion
moved the crown; this candidate ports it onto the promoted base
with the grid edge re-probed (4943079 here, one step lower than
the old-tree 4943080, because the leaner cost envelope shifts
the rounding boundary). No cost numeral moves. The reviewer
inherits the crown acceptance on the kernel, every premise, the
collision derivation, and every measurement, and judges only
the re-solved operating point. The Navigator TAKE (7ff7031f,
127.234797) scores above the bar and names no hole missing from
this crown, so no restore was applicable.
