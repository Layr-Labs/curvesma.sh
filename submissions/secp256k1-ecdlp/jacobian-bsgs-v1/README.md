# Jacobian-coordinate deterministic BSGS

Same deterministic baby-step/giant-step search as the affine control, with
the group steps evaluated in Jacobian coordinates and the table normalised
by one Montgomery batch inversion.

- `claim.json` — profile and cost-model bindings, submitted bounds.
- `submission.json` — the claim statement and the toy experiment binding.
- `proof.md` — construction, success argument, projective kernel bounds.
- `experiment.py` — toy recovery check plus the exact integer ledger.

`experiment.py` reads `{"order": n}` on stdin and prints JSON. It never
executes the production search; the production bounds are analytic.
