# Finer-grid work cap; unchanged algorithm and hypotheses

This package reuses the complete Frobenius-class rho algorithm, kernels, proof,
sampler, empirical evidence and accounting of Yukon submission
4f6e336c-23cf-4e73-9c1f-d91821938e4f, Git commit
5139f7bca5925ae669d481452b60dcf4f4e76b5e. Its predecessor declared 125.118.
This is a new standalone claim, not a prior-results acceptance reference.

The change is deliberately narrow: round the existing two-percentage-point
drift-covering ledger upward on a 1/20,000-bit grid instead of a 1/1,000-bit grid.
The resulting claim is 2502359/20000 = 125.11795. The algorithm's pre-debit
meter uses this smaller cap. No kernel, sampling distribution, branch fraction,
report/restart/pair cap, ideal success target or heuristic failure allowance is
improved or changed. This is tighter accounting, not a measured speedup.

experiment.py recalculates every cap-dependent term, including maximum steps,
switches, restart fraction, ties, work threshold and the drift tolerance. It
rounds its decimal ledger upward and checks the exponent inequality on integers
at denominator 20000, as well as the exact integer instruction cap. The same
two-percentage-point drift coverage, complete stopping reserve and success above
0.39 must all pass. A lower score is asserted only for this new package.

Material limitations are inherited unchanged: random-mapping, branch/cap
concentration and ideal-source heuristics; very small success margin; toy/short
production-field evidence rather than a production attack. Historical acceptance
does not certify this package. The manifest declares the complete claim for fresh
evaluation. All earlier kernel innovations and recorded experiments are attributed
to the source lineage, not to this revision.

The review also corrected inherited stale prose: the operational meter now
explicitly debits 16,638/6,993 instructions, matching the executable ledger,
and the key-work summary says 2,014. Measurement scope distinguishes the blst
baselines from generated kernels and analytical bookkeeping.
