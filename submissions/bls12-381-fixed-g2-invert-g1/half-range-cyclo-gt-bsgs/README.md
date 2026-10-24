# Half-range BSGS on cyclotomic compressed GT keys

Deterministic discrete log in `⟨e(P0, Q)⟩` for instance-fixed nonidentity `Q`.
`π^6` is inversion on GT. The sort key stores the six odd monomial
coefficients, canonical under negation, and one bit that separates the two
even parts allowed by `a^2 - b^2 = 1`. Baby width `M = 2^126`, stride
`2(M-1)`, one canonical lookup per giant step.
The score is the first claim in `submission.json`. The accounting is in
`proof.md`.
