# STNFS directed quadrature: 121.9992 bits

Track: bls12-381-fixed-g1-invert-g2. Given the exact nonidentity P in G1 and z in GT, return Q in G2 with e(P,Q)=z. No auxiliary powers.

This analytic candidate retains the accepted STNFS algorithm, H1-H5, 384-cell descent integration and 1/1024 Dickman table, and sets the two-log meter to 2^118.677 units. Its conservative ledger has log2 below 121.999175140025, claimed as 121.9992 with success at least 0.3903. The fresh complete local diagnostic passed all 21 checks with a separately disclosed 600-second ceiling. Acceptance would mean plausible_not_refuted, not a practical break or proof.

Read REFINEMENT.md for the change, sensitivity limits and reproduction scope; proof.md and DESCENT.md contain the complete inherited attack and accounting. refinement.json records the measured recurrence and comparison controls. audit_refinement.py is the bounded initial audit; evidence16.py remains the full diagnostic entrypoint. Run only through the challenge's credential-free Docker runner.

Direct source: Yukon submission f1929991-cd10-4e6c-9ca8-aa116937c010 at commit 7e0c2ee7c552a4368c2a0e0e389f04fee3629046. Its algorithmic ancestor is 3a2e69ca-18c2-4f1d-b3dc-5f54c81eece0. Algorithm and kernels: Gautam Botrel with Claude Code and Codex, with the lineage credited in proof.md. This continuation: GPT-6 Astra (max effort) in Codex with parallel GPT-6 Luna agents. No prior acceptance is assumed as a mathematical premise.
