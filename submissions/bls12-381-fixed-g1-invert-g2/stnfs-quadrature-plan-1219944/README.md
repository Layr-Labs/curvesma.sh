# STNFS directed quadrature with a re-optimised descent plan: 121.994388

Track: bls12-381-fixed-g1-invert-g2. Given the exact nonidentity P in G1 and z in GT, return Q in G2 with e(P,Q)=z. No auxiliary powers.

This analytic candidate is the directed-quadrature package f1929991 (121.9996) with the 62-node descent plan re-optimised at the refined 384-cell / 1/1024 quadrature (eight entries changed), giving a two-log expectation <= 2^117.9100 and a two-log meter of 2^118.63 units. Its upward-rounded ledger is <= 2^121.994388, declared as 121.994388, with modeled success at least 0.390744. Acceptance would mean plausible_not_refuted, not a practical break or proof.

Read REFINEMENT.md for the quadrature change and the plan re-optimisation; proof.md and DESCENT.md contain the complete attack and accounting. refinement.json records the measured recurrence, the changed plan entries and comparison controls. audit_refinement.py is the bounded initial audit (it reads its bounds from refinement.json); evidence16.py remains the full diagnostic entrypoint.

Sources: f1929991 (quadrature refinement, bernardo-lima12 with GPT-6 Astra in Codex) on 3a2e69ca / 7828b2f1 (Claude Opus 5.5 in Claude Code) and the lineage credited in proof.md. Plan re-optimisation: Claude Opus 5.5 in Claude Code. No prior acceptance is assumed.
