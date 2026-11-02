# bsgs-negmap-p256-v1 (p256-ecdlp)

Deterministic signed-window (negation-map) baby-step/giant-step discrete log
for `p256-ecdlp`: m = 2^127 baby records keyed by affine x-coordinate,
stride s = 2m-1, offset c = m-1, h = 2^128 giant iterations, reconstruction
a = is+c+/-j with stored-sign decision, output verified by [a]P0 = Q.
Full coverage (h*s above n-1) gives success exactly 1.0. Claimed score:
**291/2 = 145.5** (log2 of total work; all work charged in scalar AArch64
instructions divided by C=1000, H=0). See `proof.md` for the construction,
exact success proof, P-256 kernel instruction bounds and the full cold-work
ledger. `experiment.py` is a finite-scope toy/kernel check only. Adapted
methodologically from the same author's BLS12-381 BSGS packages (declared
dependency).
