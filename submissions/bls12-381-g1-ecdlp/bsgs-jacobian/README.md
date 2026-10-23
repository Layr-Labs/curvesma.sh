# bsgs-jacobian

Candidate for `bls12-381-g1-ecdlp`: deterministic baby-step/giant-step scalar
recovery on `T = [s]P0`, using only the public `A_1 = T`.

The construction walks the baby and giant lists in Jacobian coordinates with
inversion-free mixed additions, normalises to affine x in blocks of `2^16` with
Montgomery batch inversion, sorts both lists by the 381-bit x key with a 24-digit
radix sort whose bucket counters are 128-bit multiword values, and finds the match
with one linear merge. Success probability is 1; every instruction is counted once
in `R` and the submitted bound is `W = R/2000`.

`proof.md` carries the complete ledger, the exceptional-case handling, the memory
and preprocessing figures, and the covering argument. `experiment.py` exhaustively
checks the same construction on small Weierstrass groups for every `s` in `1..r-1`.
