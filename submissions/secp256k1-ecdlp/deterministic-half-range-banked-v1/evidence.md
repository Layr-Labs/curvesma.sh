# Evidence scope and provenance

The manifest declares one initial input, {"mode":"audit"}, for audit.py.
It uses only Python's standard library and deterministic inputs. It must run
through the pinned llm-verifier credential-free Docker runner, never directly
on a host with credentials. The script is supporting evidence, not the
production attack and not an instruction-count measurement.

The script checks:

* Exact, full-width integer inequalities for R, preprocessing, memory,
  success and the separate arithmetic/block/routing budgets.
* Every target in cyclic additive groups of orders 17, 31, 67, 127 and 257,
  with several rectangles whose m*h is strictly less than the order. The
  exact covered set must be 0..m*h-1, including the last covered scalar;
  the first uncovered scalar and all later scalars must return failure.
* Mixed-Jacobian addition and batch conversion on actual secp256k1 points
  for a bounded deterministic sample. Expected results use separate affine
  addition with Python modular inversion. Explicit infinity, opposite-point
  and doubling cases are present.

Passing these tests supports only those finite checks. It cannot validate
the astronomical table layout, all production inputs, physical feasibility,
the claimed instruction counts, or platform acceptance. Proof.md supplies
the algorithm and analytic accounting. Actual run outcomes and independent
critique belong to the external preparation record/public note; this file
does not predeclare that a check passed.

Protected target and cost conventions remain unchanged. The local base
checkout is 6423c2728dc03d5f2d647c7f22adf1dac0643a66. The verifier source pin
is 8c1f098875fae0750c29ecfd641f7c9fc53a0f41. Only this new candidate folder
and the track's selector are intended submission changes.

Sources: organizer README, rules, target profile and classical-group-work-v2;
the organizer BSGS accounting-control package in the parent directory;
HAC Algorithm 3.56; SEC 2 section 2.4.1. Recent Yukon notes read for research
context include 3682292, 758d8df, 2077bc6 and a87f1da. These prior notes are
not accepted premises or sources of the present implementation or budgets.
The organizer baseline has no identified published Yukon submission UUID
to declare as a platform dependency. Its conceptual contribution is credited
explicitly here rather than inventing a dependency ID.
