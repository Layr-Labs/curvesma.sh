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
bound, so W <= 2^126.18266 on every run.

Claim: time_log2 = 126.18266, success >= 0.39 (>= 0.390014 under H1), memory
< 2^62.3 bytes (all arrays 64-bit addressable), preprocessing < 2^56.01 units
(capped at 2^60), no advice. Distinguished probability 2^-73; reports go to an
append-only log with a linear-probing index (load <= 0.16, <= 256 probes).

Built on the v21/v31 six-class rho candidates of submission 9b9712fa (same classes,
key, twisted table, escape and solving); the changes and their instruction
savings are listed at the top of proof.md.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity and listing checks + lockstep toy pipeline](experiment.py)
- [Offline toy simulator used for the section-6 sweep](tcsim.c)

Retune of the accepted frontier submissions 07207131-1ff3-4df8-be06-639c7974acc3 (official 126.18324, the 0.39002 quantile retune of 1a3470dc-0997-4986-86ae-0bd81c6fa651 at official 126.18326, itself the kappa = 1.003 retune of dbc45587-0975-46ec-9599-126119972d31 at 126.18613): the only changed budget parameter is the H1 slowdown allowance kappa, 1.003 -> 1.0026, set at the one-sided 95% upper confidence limit (0.259%, rounded up to 0.26%) of the 15-row fitted toy intercept and below the accepted 0.30%; the quantile stays at the accepted 0.39002. S_ops, pi (S_ops + 2^22) = 2^126.182654, the headroom (> 2^107.92), the meter bound (< 0.0000012), success (> 0.390014) and the claimed cap (126.18266) are recomputed, everything else is unchanged. See proof.md sections 4-6.
