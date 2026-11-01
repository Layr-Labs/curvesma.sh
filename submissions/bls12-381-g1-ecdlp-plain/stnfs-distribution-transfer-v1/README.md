# Transfer of the fixed-G1-invert STNFS core to bls12-381-g1-ecdlp-plain

`proof.md` reduces the plain G1 logarithm to the fixed-G1 pairing inversion
`z = e(T, Q0) = e(P0, Q0)^x` and charges the core `e3940d11` (2^121.9944 units at
C = 20,000, times 101/100) converted to the plain track's C = 2,000 units, plus a 2^80
wrapper: `W < 2^125.3306834`, success >= 0.39.

`source-core.md` is the core's complete argument, unchanged; every core file it cites
is shipped. `combined_evidence.py` runs the core's twenty-one checks and then the
transfer checks; `transfer_check.py` is the numeric conversion model and the toy
bijections.
