# Corrections to inherited presentation

The first heading of proof.md retains an obsolete predecessor exponent. The
current complete claim is the manifest's exact rational 116.013974026666, with
core meter 116.013974026665. The heading is historical and supplies no score.

The introductory reference to evidence16.py describes the direct parent's
full-check driver, not an included or freshly executed experiment. That optional
driver is omitted, as explained in proof.md's Revision and evidence status.
The selected audit is audit_refinement.py with input {}. Other audit inputs are
optional. No current local isolated execution has succeeded.

The audit output's field success_lower_with_wrapper subtracts 2^-48 from the
already rounded 0.4456, and therefore prints a weaker bound slightly below
0.4456. Its unrounded computed core success is greater than 0.44561939678;
subtracting 2^-48 from that value still supports the manifest's 0.4456. The
stronger manifest bound follows from this unrounded calculation, not from
rounding the weaker printed field upward. No probability or cost changes.

These corrections are a new add-only file; previously written candidate files
are preserved under the repository's rule.
