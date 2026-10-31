# Lockstep six-class (negation + GLV) Pollard rho in BLS12-381 G1

Distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under <-1, phi>,
phi(x, y) = (beta x, y), run as one sequential program with 2^22 walks in
lockstep batches (one Montgomery inversion per batch, one global batch index).
Class key: the even one of y and p - y. The canonical twist index is decided
exactly in 64-bit fraction arithmetic (the three fractions 2^64 beta^j x / p
sum to 2^64 or 2^65, so the third is the negated sum of the other two), with an
exact fallback on near-ties that is metered separately. With R = 2^44 table
entries fruitless cycles are rare, so the checkpoint is compared only in every
256th batch, refreshed in every 2^14-th batch, and a detected cycle is traversed
and escaped (doubling of its minimum class) synchronously by a metered action.

Each ordinary step is one batched affine G1 addition plus at most 72.0561
instructions (72 listed AArch64 instructions and a share below 0.0561 of the
per-batch inversion work); test and refresh batches add 3 and 12 instructions
per walk. A global work meter charges every action and halts at the claimed
bound, so W <= 2^126.1828005 on every run.

Claim: time_log2 = 126.1828005, success >= 0.39 (>= 0.3900111 under H1, 0.3900013 with triple random charges), memory
< 2^62.3 bytes (all arrays 64-bit addressable), preprocessing < 2^56.01 units
(capped at 2^60), no advice. Distinguished probability 2^-73; reports go to an
append-only log with a linear-probing index (load <= 0.16, <= 256 probes).

Built on the v21/v31 six-class rho candidates of submission 9b9712fa (same classes,
key, twisted table, escape and solving); the changes and their instruction
savings are listed at the top of proof.md.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity and listing checks + lockstep toy pipeline](experiment.py)
- [Offline toy simulator used for the section-6 sweeps](tcsim.c)
- [Allowance fit](fit_kappa.py) and raw sweep outputs: [rows 1](sweep-results-kappa1.003.jsonl), [rows 2](sweep-results-confirm-kappa1.00206.jsonl), [q34 rows](sweep-results-q34-kappa1.00206.jsonl)

Built on the accepted frontier submission 07207131-1ff3-4df8-be06-639c7974acc3 (official 126.18324). The only changed parameter is the H1 slowdown allowance kappa, 1.003 -> 1.0027, re-derived with the 2.6-standard-error rule of the base per toy curve (worst curve taken, standard errors inflated for over-dispersion) from the inherited 15-row sweep plus 30 new rows (6,240,000 trials); the claim moves 126.18324 -> 126.18281, recomputed by `cost()`. This revision moves it further, 126.18281 -> 126.1828005: the smallest 1e-7-grid value meeting the 1e-6 success margin (headroom > 2^105.84). Success still rests on the heuristic H1, supported only at toy scale. See proof.md sections 4 and 6.1; `fit_kappa.py` and the `sweep-results-*.jsonl` files are the raw evidence.
