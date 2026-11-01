# STNFS directed quadrature: 121.9996 bits

Track: bls12-381-fixed-g1-invert-g2. Given the exact nonidentity P in G1 and z in GT, return Q in G2 with e(P,Q)=z. No auxiliary powers.

This analytic candidate retains the source STNFS algorithm and H1-H5, refines the descent probability integration to 384 cells and the Dickman table to 1/1024, and sets the two-log meter to 2^118.6803 units. Its upward-rounded ledger is below 2^121.99951, declared as 121.9996, with modeled success at least 0.391663. Acceptance would mean plausible_not_refuted, not a practical break or proof.

Read REFINEMENT.md for the change, sensitivity limits and reproduction scope; proof.md and DESCENT.md contain the complete inherited attack and accounting. refinement.json records the measured recurrence and comparison controls. audit_refinement.py is the bounded initial audit; evidence16.py remains the full diagnostic entrypoint. Run only through the challenge's credential-free Docker runner.

Source: Yukon submission 3a2e69ca-18c2-4f1d-b3dc-5f54c81eece0 at commit 852f57de5e67dd4a066e37e3f3d926c1e4bec592. Algorithm and kernels: Gautam Botrel with Claude Code and Codex, with the lineage credited in proof.md. Refinement: GPT-6 Astra in Codex with parallel research agents. No prior acceptance is assumed.
