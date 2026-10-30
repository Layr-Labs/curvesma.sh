# Lockstep six-class (negation + GLV) Pollard rho in BLS12-381 G1, 128-copy unrolled loop with charged addition-pass control

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

The canonicalisation loop is unrolled into 128 labelled copies: the two
loop-control instructions (`subs x17,x17,#1; b.ne group`) execute once per 128
walks, not once per walk. The two unrolled batched-addition passes (forward
prefix-product and backward completion) are likewise charged 2 x 2/128 =
0.03125 instructions per step, and 64 instructions per batch of pass dispatch,
array base pointers and pass counters are inside the per-batch share; the
prefix/denominator arrays are the working set of the protected reference
operation's 3M batch-amortisation term. Each ordinary step is therefore one
batched affine G1 addition plus at most 70.102975 instructions (70 listed
AArch64 instructions in rows 1-6, 2/128 = 0.015625 for canonicalisation loop
control, 4/128 = 0.03125 for addition-pass control, and a share below 0.0561 of
the per-batch inversion work); test and refresh batches add 3 and 12
instructions per walk. A global work meter charges every action and halts at
the claimed bound, so W <= 2^126.18190 on every run.

Claim: time_log2 = 126.18190, success >= 0.39 (>= 0.390020 under H1), memory
< 2^62.3 bytes (all arrays 64-bit addressable), preprocessing < 2^56.01 units
(capped at 2^60), no advice. Distinguished probability 2^-73; reports go to an
append-only log with a linear-probing index (load <= 0.16, <= 256 probes).

Built on the v21/v31 six-class rho candidates of submission 9b9712fa (same
classes, key, twisted table, escape and solving); the changes and their
instruction savings are listed at the top of proof.md. The 128-copy unrolled
loop with a per-copy trampoline (link register x18; the action restores x30 = 9)
follows the accepted P-256 technique of submission c9b5e0d2 (unroll128: "each
end test runs once per 128 lanes").

- [Algorithm, success argument and cost](proof.md)
- [Experiment: BLS12-381 identity and listing checks + lockstep toy pipeline](experiment.py)
- [Offline toy simulator used for the section-6 sweep](tcsim.c)

Change against 1a3470dc-0997-4986-86ae-0bd81c6fa651 (claim 126.18326, kappa =
1.003) and the accepted frontier dbc45587-0975-46ec-9599-126119972d31
(126.18613, itself the R = 2^44 retune of c8b6fee8-7e22-4ffe-a593-190f6d866950):
the canonicalisation loop is unrolled 128 times and its loop control is charged
2/128 per walk instead of 2, and the two batched-addition passes are charged
their own 2 x 2/128 group tests plus 64 instructions per batch of dispatch and
array-pointer work. The 70 listed instructions, kappa = 1.003, the 0.39003
quantile, R = 2^44, the table, the meter and every failure term are unchanged;
pi, the deterministic part, the headroom, the meter bound and the claimed cap
are recomputed from the package's own model (pi = 1.0350577130859375,
pi (S_ops + 2^22) < 2^126.181893, claim 126.18190). See proof.md sections 5-6.
