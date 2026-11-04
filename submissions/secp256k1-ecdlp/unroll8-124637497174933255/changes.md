# unroll8 provenance

Direct predecessor: submission `6c29697` (`zero512-124637820464`), which itself
extends `aa8e1f24` by a four-way DC ZVA zeroing unroll. This candidate changes
only the production backward entry point: `bwd8` runs eight copies of the
unchanged, bit-exact `BWDI` body per loop test, so the loop-test cost per
incidence falls from 1/2 to 1/8 instruction and the executed count becomes
1281.125 instead of 1281.5. The listing's original two-copy `bwd` entry is
retained unchanged for the small-input correctness tests; both entries execute
the same `BWDI` body. Nothing else changes: the forward kernel, key pass,
inverse, routing, success proof, memory and preprocessing are the predecessor's.
