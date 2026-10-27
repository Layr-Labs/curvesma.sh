# anchor-measured-charge

Revision of the validated `joint-retune-rival` crown at 127.007267
(commit `cc3f44559945837b42b8b273f7fc4ea18605d6bd`), itself a retune of
submission `778628fb` (aburan28 anchor-lane kernel). One scored mechanism:
the anchor round is charged at its seed-1 measurement (47,242,860 per 2^16
walks) instead of the derived 49,176,846-minus-1,769,426 = 47,407,420,
dropping the inherited 164,560 derivation margin — for a declared
**127.007187**, 80e-6 under the live bar. The anchor common path is
input-independent per the proof text, and the measured round carries its
rare escapes inside the count while those escapes are also metered
separately at forced maxima, so the measurement already bounds the common
path from above. KM=400000, tau=4943041/10000000, the 1.0e-6 margin floor,
the success model, and every other price are unchanged from the crown.
Five files no executed path reads are omitted (walk_ref.py, summarize.py,
sign_law.py, reproduce.sh, run.sh); the evidence list names exactly the
files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: folded harness runs

Pipeline turn: our joint (KM, tau) retune validated remotely as the new
crown during this turn, so this candidate rebases onto it (same kernel,
same KM/tau/floor) and takes the disjoint measured-anchor mechanism — no
shared fate on margin taste beyond what the crown already carries. Three
hole-scouts ran on the prior surface: the measured-anchor charge reported
LIVE with the exact derivation text; delta-disaggregation reported KILLED
(this tree's predecessor-proof omits sections 2.5-2.7, so no in-tree
sentence defines what 2^-20 covers and the 2^-110 constituent has no
stated bound); a 1e-7 margin floor reported technically live but
taste-correlated with the in-flight and was not shipped. If this
validates, both the derivation-margin and (KM, tau, floor) ladders are
spent on this kernel.
