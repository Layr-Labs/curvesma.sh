# control-112-rare-128

Pure cost-ledger repricing of the accepted 127.005503
control-160-cost-only crown (forward loop unrolled by four,
checkpoint-only cycle detection). Same kernel, same
measurements, same success ledger, same premises. Two scored
mechanisms that stack: MAIN_CONTROL 160 -> 112 against the
proof's own itemised 105 bound (7 headroom, factor 1.07), and
rare-path OVERHEAD 256 -> 128 against the measured step_rare
envelope (factor 1.83 over the longest path of 70; w2 82,750
-> 82,622, wL 2,739,767 -> 2,731,447). Declared **127.005501**,
strictly under the bar by two micro-units: total float
127.005502901 -> 127.005500901 (9.9e-8 of rounding headroom),
C_COMMON 713.787745 -> 713.787012, preprocessing, p_claim,
tau, KM, delta, every rate, every forced maximum, and all
measurements unchanged. The evidence list names exactly the
files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: review posture

This candidate was built on the public crown while the
control-130 deepening validated; its cost-side envelopes are
disjoint from, and sharper than, that in-flight patch. The 112
charge keeps 7 of headroom over the itemised 105; the rounding
boundary with overhead 128 sits at 115, so 112 keeps 3 units
past the edge. The 128 overhead keeps 58 over the measured 70
longest path. The Navigator TAKE and NEXT are both not
git-visible on this crown and name no hole missing from it, so
no grafted restore was applicable.
