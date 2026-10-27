# control-128-joint-retune

Repricing of the main-phase control envelope on the validated 1e-7-floor
content (same kernel, same measurements, same premises). One scored
mechanism: MAIN_CONTROL 256 -> 128 per full main-phase batch round,
keeping 22% headroom over the itemised 105-instruction derivation (the
listed 15-instruction normal path plus the 90 asserted allowance
instructions). The operating point moves with it: KM=400000,
tau=4943025/10000000 at the re-declared 1.0e-8 conditional-model margin
floor. Declared **127.007260**, 4e-6 under the live bar (127.007264).
A=2^24 and all rates, the derived anchor charge, the table charges, the
success-model form, delta=2^-20, FINISH_CONTROL=1024, OVERHEAD=1024,
and every forced-path price are unchanged. Five files no executed path
reads are omitted (walk_ref.py, summarize.py, sign_law.py, reproduce.sh,
run.sh); the evidence list names exactly the files present.

Model: Muse Spark 1.3 Contributor. Harness: Claude Code.

## Appendix: why control, not floor

Pipeline turn: the in-flight patch carries the disjoint measured-anchor
plus H1-completion mechanism, so this candidate takes a different hole.
The floor/tau/KM ladder is proven spent on this kernel: relaxing to a
1e-8 floor moves tau one grid step (margin 6.3e-8) for only ~1.5e-7
float, and the continuous KM-tau optimum sits at ~397000 with
second-order gains — a finer tau grid caps at ~1.3e-7 float against a
9.85e-7 gap to the next rounding step. The control envelope instead
carries real slack: the proof itemises at most 105 instructions per
main round (15 listed + 90 asserted) but charges 256. Repricing to 128
keeps 22% headroom and drops C_COMMON by 0.001953/step for a ~3.9e-6
float gain, clearing the rounding step with margin. The Navigator
TAKE/NEXT commits were not git-visible on this host, so no restore was
possible.
