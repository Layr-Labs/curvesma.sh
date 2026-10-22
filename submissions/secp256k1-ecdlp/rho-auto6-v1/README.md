# secp256k1-ecdlp: truncated Pollard rho on automorphism classes

Classical Pollard rho on the order-6 automorphism classes of secp256k1
(`⟨−1, φ⟩`), truncated to success probability ≥ 0.39, accounted under
`classical-group-work-v2`.

- Claimed `time_log2`: **127.868282** (organiser BSGS control: 129.000802)
- Memory `log2` bytes: 85.868282
- Success probability: ≥ 0.42 under disclosed standard heuristics
- [Algorithm and cost argument](proof.md)
- [Submission claim](claim.json)
- Toy experiment: `experiment.py` (j=0 curves, 24–30 bit order)

This is not a practical break, a full-size execution or a claim of a new
cryptanalytic technique. Acceptance means plausible_not_refuted under the
track's exploratory review.
