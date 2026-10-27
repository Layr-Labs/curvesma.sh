# control-130-cost-only

Pure cost-ledger deepening of the accepted 127.005503
control-160-cost-only crown (forward loop unrolled by four,
checkpoint-only cycle detection, rare-path overhead 256). Same
kernel, same measurements, same success ledger, same premises.
One scored mechanism: MAIN_CONTROL 160 -> 130 against the
proof's own itemised 105 bound (25 headroom, factor 1.24).
Declared **127.005502**, strictly under the bar: total float
127.005502901 -> 127.005501976 (2.4e-8 of rounding headroom),
C_COMMON 713.787745 -> 713.787287, preprocessing, p_claim,
tau, KM, delta, every rate, every forced maximum, every meter,
and all measurements unchanged. The evidence list names exactly
the files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: review posture

The 130 charge keeps 25 instructions of headroom over the
itemised 105 the crown's own listing derives (15-instruction
normal path plus 32 save/restore pairs, 6 stack/flags, 20
address/setup, 32 handoff = 105), so the reviewer inherits the
crown acceptance on the listing and judges only the tightened
envelope. MAIN_CONTROL moves the total linearly at ~3.1e-8
bits per unit, so the 30-unit cut 160 -> 130 moves the float
one full micro-unit to 127.005502; the boundary sits at 133
(floors 127005503 at 132-138 probing), leaving 27 units of
headroom over the itemised bound at the rounding edge. The
Navigator TAKE (7ff7031f, 127.234797) scores above the bar, is
not git-visible on this crown, and names no hole missing from
it, so no grafted restore was applicable. The rejected tau-floor
retune is deliberately not repaired here: Yukon's verdict
required the original 1e-5 margin, which this candidate keeps
(p_claim 0.390010015, margin 1.0015e-5).
