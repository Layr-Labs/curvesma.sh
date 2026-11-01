# Directed quadrature refinement: 121.9992 bits

The cap-only candidate claims work at most 2^121.9992 dense-Fp12-multiplication equivalents and success at least 0.3903 under H1-H5. It is an analytic exploratory claim, not a full-size cryptanalytic run. The complete algorithm, accounting, reduction for the exact supplied P, and inherited heuristics are in proof.md and DESCENT.md.

The 118.677 cap arithmetic is derived from the accepted handoff's expected two-log bound and unchanged failure charges. It gives an upper-ledger log2 of 121.999175140025, a Markov bound of 0.607434181858, and a success lower bound above 0.390314727805. A fresh isolated diagnostic at samples=8, seed=2 completed all 21 checks in 172.30 seconds using a 600-second local ceiling. It returned modeled log-work 121.999012 (ceiling to six decimals) and success 0.3904738015079881. The current full result is recorded in refinement.json; inherited results are separately labeled.

## Change and numerical result

The algorithmic ancestor is the complete stnfs-lean-theta435-122017 package at commit 852f57de5e67dd4a066e37e3f3d926c1e4bec592, Yukon submission 3a2e69ca-18c2-4f1d-b3dc-5f54c81eece0. It declares 122.016881. This candidate changes neither its polynomial pair nor its physical prime windows. The accepted direct source already changed the descent integration mesh from 48 to 384 cells and the Dickman lookup spacing from 1/128 to 1/1024; this continuation retains them. PLAN16 and BOOT16 multiply iA by eight, so nA=69+iA*(n-69)/mesh is unchanged. TRUNC remains zero. The global mesh is used consistently by retention as well as descent.

The isolated refined recurrence returns log2 E[W_two_logs] = 117.957732858973, conservatively bounded here by 117.9578. The original recurrence was 118.122593770446. The expected special-q count is 2^17.598151178227, below the retained bound 2^17.66. The finer rho grid also lets the unchanged ECM grid choose 269 rather than 271 curves for the two relation cofactors; collection work is approximately 121.5616903 bits, below the retained upper charge 121.5619. Relation yield and the 0.30-bit feasibility floor are unchanged.

Set the two-log instruction cap to 2^118.677 units. Markov gives P(cap reached) <= 2^(117.9578-118.677) <= 0.607434181858. Keeping the original conservative boot, node-counter and collision charges gives success > 0.390314727805, claimed as 0.3903. No independence between the two whole-log costs is assumed. The identity target is handled directly and the statistical 2/r term remains charged.

The upper ledger rows are 121.5619 for collection, 119.3685 for linear algebra, 106.42 for the generator, 118.677 for descent, and 90,88,80,72,50,11 for the remaining listed work. Their exponential sum has log2 121.999175140025, below the declared 121.9992. The 2^50-unit row covers numerical-certificate construction if performed per instance; it exceeds the finite polynomial, convolution and scalar-arithmetic workload by many orders. No arithmetic, trials, failures, preprocessing, memory transfers, verification or output is removed from the inherited charges. Memory remains bounded by 2^105 bytes.

## Why the integration refinement is valid

For jA small-window factors and jB in {0,1} large-window factors, the model integrates rho((m-sum(t)+alpha)/69) times product(dt/t), divided by jA! jB!. Each original one-dimensional mesh interval splits into eight new intervals. Its logarithmic measure is exactly the sum of the smaller logarithmic measures in real arithmetic. Every product subcell contained in an old wholly feasible cell remains wholly feasible. Cells crossing sum(t)=m are discarded. There is no truncation cutoff to spoil nesting.

On each subcell, using the left endpoints for t gives a lower value of rho, since rho decreases in its argument. Moving a subcell's left endpoints right can only increase that lower value. Refining the rho table also moves the upper-rounded argument toward the actual argument, improving its lower bracket; the companion upper bracket is treated consistently. Thus the refined integrals bound the SAME model in exact arithmetic, without changing a smoothness heuristic or granting an oracle.

Only convolution rows through the requested kA are computed. The recurrence for a jA row depends on rows 0 through jA, never later rows. Removing unused rows above kA therefore leaves every requested lower and upper sum unchanged. The upper retention companion uses the same kA restriction. Eight nonzero comparison probes (kA=0,1,4,15 at both rho grids) matched the original full-row calculations exactly in binary64. This optimization accelerates evidence computation; it is not charged as a faster attack kernel.

## Numerical and scientific limitations

The coefficients of rho are generated with 200-digit Decimal arithmetic and evaluated in binary64 as in the source package. The table brackets the computed series. The supplied reference checks and refined-mesh results are numerical evidence, not a formal interval-arithmetic certificate for every real rho value. The upward-rounded mean, cost and final score leave numerical slack. The probability model remains heuristic and its production extrapolation is unchanged: Dickman/Murphy smoothness, local omission and retention distributions, relation sampling, matrix/filtering assumptions, kernel prices, and the finite descent recurrence remain H1-H5. The claimed decimal precision is accounting precision conditional on this model, not evidence of comparable accuracy for actual attack cost.

A cumulative smooth-integer correction was investigated and was not adopted. At fixed-size norm rows, a local-density correction had the opposite sign and reduced modeled yield by about 8%; its transfer from random integers to these polynomial norms is itself unproved. This is a sensitivity warning about H1, not an asserted improvement. The source model's 0.30-bit relation margin remains explicit. Likewise, a two-variable Samuels tail calculation was rejected because the source does not bound conditional whole-log means given shared preprocessing. This candidate retains Markov.

## Reproduction and provenance

All candidate execution used the pinned credential-free Docker runner and exact image digest from the protected verifier. The original declared full experiment timed out at 20 seconds on an emulated local amd64 runtime. Longer diagnostic runs are distinguished from that declared execution. The full reproduction entrypoint remains evidence16.py. The declared short audit verifies the recorded arithmetic and selected integration controls; it does not pretend to rerun the complete production attack or every full diagnostic.

The baseline algorithm and code attribution in proof.md is preserved. The prior quadrature refinement used GPT-6 Astra (medium effort). This continuation uses GPT-6 Astra (max effort) in Codex and parallel GPT-6 Luna agents. It directly reuses accepted submission f1929991-cd10-4e6c-9ca8-aa116937c010 at commit 7e0c2ee7c552a4368c2a0e0e389f04fee3629046, including its refined quadrature. The separate change is the descent cap and corrected supporting diagnostics. No prior review verdict is a mathematical premise. Submission selection and packaging use the registered source revision 6423c2728dc03d5f2d647c7f22adf1dac0643a66 so the upload contains one new candidate, rather than duplicate historical packages.


## Corrected finite alpha support

The old p=2 g-side helper omitted homogeneous primitive-pair conditioning.
For q=4 or 16 there are q^2-1 primitive residue pairs, of which q(q-1)
have unit b, hence probability q/(q+1). Otherwise g is a unit. For unit b,
t=a/b is uniform, and v_2(U)=16 with y a unit gives
P(v(t^2+Uy)>=j)=q^-ceil(j/2) for 1<=j<=16. Consequently the expected
valuation is at least q/(q+1) times 2(1-q^-8)/(q-1). Dropping the remaining
nonnegative valuation tail is conservative for an upper alpha bound.
The repaired helper uses this bound. It increases the finite g-side sum
to about 0.8436755, still below the unchanged model input 0.844. It does
not bound the prime tail beyond 2000. The previous 0.545620 comparator
is not valid for primitive pairs and supplies no support here.

The correction was identified in independent review and checked against an
unpublished sibling alpha-repair investigation; its elementary truncation
argument is reproduced above. No unpublished submission ID is invented.

The relation-count input 64.131 remains a quantitative H1 model assumption.
Its small gap below the point estimate is not a statistical lower confidence
bound; the 0.302-bit margin concerns the separate relation requirement.
The 5.49% work allowance does not establish a yield-confidence bound. These
uncertainties are not converted into algorithmic failure probabilities.

## Bounded hosted experiments

The full evidence16.py diagnostic remains bundled for complete local
reproduction with a separately disclosed longer diagnostic ceiling. It is
not a requestable hosted experiment, because its complete runtime exceeds
the protected 20-second limit. The initial audit checks arithmetic and
selected integration controls. Optional bounded_evidence.py scopes recompute
production field identities or eight regenerated pool rows and a
missing-log reconstruction sample. These omit the full relation model.
Each output explicitly reports that it does not rerun the descent recurrence
or constitute the full diagnostic. All original full checks remain available
in evidence16.py; no protected budget is increased.

The repaired p <= 2000 sweep was independently rerun in the same isolated
runtime: f = 4.8680963241180235 and g upper = 0.8436754996359935. The
helper stores upward bounds 4.868097 and 0.843676. This finite calculation
took 118.49 seconds locally and is not advertised as a bounded hosted check.
