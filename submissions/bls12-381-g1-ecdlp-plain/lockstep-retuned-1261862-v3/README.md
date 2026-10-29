# Lockstep six-class (negation + GLV) Pollard rho in BLS12-381 G1

Distinguished-point Pollard rho on the (r-1)/6 orbits of G1 under <-1, phi>,
phi(x, y) = (beta x, y), run as one sequential program with 2^22 walks in
lockstep batches (one Montgomery inversion per batch, one global batch index).
Class key: the even one of y and p - y. The canonical twist index is decided
exactly in 64-bit fraction arithmetic (the three fractions 2^64 beta^j x / p
sum to 2^64 or 2^65, so the third is the negated sum of the other two), with an
exact fallback on near-ties that is metered separately. With R = 2^40 table
entries fruitless cycles are rare, so the checkpoint is compared only in every
256th batch, refreshed in every 2^14-th batch, and a detected cycle is traversed
and escaped (doubling of its minimum class) synchronously by a metered action.

Each ordinary step is one batched affine G1 addition plus at most 72.0561
instructions (72 listed AArch64 instructions and a share below 0.0561 of the
per-batch inversion work); test and refresh batches add 3 and 12 instructions
per walk. A global work meter charges every action and halts at the claimed
bound, so W <= 2^126.1862 on every run.

Claim: time_log2 = 126.1862, success >= 0.39 (>= 0.390024 under H1), memory
< 2^62.3 bytes (all arrays 64-bit addressable), preprocessing < 2^52.05 units
(capped at 2^56), no advice. Distinguished probability 2^-73; reports go to an
append-only log with a linear-probing index (load <= 0.16, <= 256 probes).

Built on the v21/v31 six-class rho candidates of submission 9b9712fa (same classes,
key, twisted table, escape and solving); the changes and their instruction
savings are listed at the top of proof.md.

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity and listing checks + lockstep toy pipeline](experiment.py)
- [Offline toy simulator used for the section-6 sweep](tcsim.c)

Retuned success budget based on submission ae35f5fe-c8ba-420e-88ab-c464bdafea8f; only the phase-one quantile and work cap are changed. See proof.md for updated accounting.
