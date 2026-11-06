# Pollard rho on automorphism classes — v3 (partition bump)

Third revision of the secp256k1 ECDLP submission, answering the v2
review's three objections:

1. **Two-cycles.** The step function gains a partition bump: when the
   canonicalised successor repeats the current partition index, the
   step continues with the next table entry. proof.md Section 6 proves
   this destroys the reviewer's `S_j = -(R+R')` family and every
   automatic cycle, and bounds coincidence cycles at production
   parameters (`J = 2^16`, `theta = 2^-24`).
2. **Toy gap.** The experiment program implements the bumped walk
   itself, counts total steps (repeated work included), and adds an
   exact functional-graph census: the plain walk's two-cycle states
   track `N''/(6J)` (484/232/106/66 at J=16/32/64/128) while the
   bumped walk's are 0/0/2/0; trapped trails at the v2 comparison point
   fall from 295 to 3 and steps-metric success rises 0.443 -> 0.627.
3. **Cost premises.** The split-scalar bound now has a specified
   schedule (GLV-split components, mu-maps by adds/subs, rebalancing
   every 36 steps, envelope `< 2^204`, itemised 63 of 96 instructions),
   implemented and invariant-checked in the toy code; exceptional steps
   are handled by bounded retries under a global budget of `2^10` with
   abort probability `< 2^-135`, replacing the per-trail assumption.

Claim unchanged from v2 (the repairs fit inside its envelope):
`2^127.22` work, success `>= 0.405`, memory `<= 2^113` bytes,
preprocessing `<= 2^27`. See `proof.md`, `claim.json`,
`submission.json`, `experiment.py`.
