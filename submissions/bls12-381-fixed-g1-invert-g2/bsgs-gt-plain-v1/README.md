# bsgs-gt-plain-v1 (bls12-381-fixed-g1-invert-g2)

Deterministic plain baby-step/giant-step discrete log in GT for
`bls12-381-fixed-g1-invert-g2`: D=e(P,Q0) generates GT, table D^j for
j in [0,2^128), sorted by canonical coefficient key, giant steps by
D^(-m), reconstruction a'=im+j, output Q=[a']Q0 verified by D^(a')=z.
Claimed score: **287/2 = 143.5** (log2 of total work; all work charged in
scalar AArch64 instructions divided by C=20000, H=0), success exactly 1.0,
no heuristics. See `proof.md` for the construction, exact success proof,
Fp12/Fp2/G2/Miller kernel instruction bounds and the full cold-work ledger.
`experiment.py` is a finite-scope toy/kernel check only. Adapted
methodologically from the same author's BLS12-381 G1 BSGS packages
(declared dependency).
