# unroll4 provenance

Direct predecessor: submission `6c296971` (`zero512-124637820464`), which
extends `aa8e1f24` by a four-way DC ZVA zeroing unroll. This candidate changes
only the production backward entry point: `bwd4` runs four copies of the
unchanged, bit-exact `BWDI` body per loop test, so the loop-test cost per
incidence falls from 1/2 to 1/4 instruction and the executed count becomes
1281.25 instead of 1281.5. Four `BWDI` bodies occupy at most 4*1366*4 = 21856 bytes including
their skipped fixup blocks, inside the signed 14-bit `TBZ` displacement range
of 32764 bytes; this is why the entry uses four copies rather than eight. The listing's original two-copy
`bwd` entry is retained unchanged for the small-input correctness tests; both
entries execute the same `BWDI` body. Nothing else changes: the forward kernel,
key pass, inverse, routing, success proof, memory and preprocessing are the
predecessor's.
