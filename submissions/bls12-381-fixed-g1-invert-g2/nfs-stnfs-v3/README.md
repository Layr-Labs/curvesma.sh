# nfs-stnfs-v3: hard-capped STNFS discrete-log inversion in GT

Reduces `e(P,Q)=z` to a discrete logarithm in the prime-order group GT
(base `g=e(P,Q0)`), solves it with the published Special-Tower Number Field
Sieve in `F_{p^12}` using the Guillevic-Singh BLS12-381 polynomial and Table 10
parameters, returns `Q=[a]Q0` and verifies the pairing.

The algorithm is **hard-truncated**: it counts its own instruction equivalents
and stops with FAILURE at `R_cap = 2^135.146602479037` scalar instructions.
`R_cap = 2 (Markov) x 2 (charged retries) x 128 x 2^126.146602479037`, where
`2^126.146602479037` is the GS21 STNFS total recomputed from the published
parameters, and `128 = 16 cycle-equivalents per model operation x 8 scalar
instructions per cycle` (16 is 1.41x the worst measured model/hardware ratio
published for this model family, BD19 Table 6 `log2(c_sieve) <= 3.5`; 8 is the
widest current AArch64 scalar retirement width).

Claimed total `2^120.858891` target-group-operation equivalents after the fixed
`C = 20000` price and an explicit `H = 2^14` verification charge; success
probability at least `1/2` over uniform `P`, uniform `z` and fresh coins (Markov
applied to total work), above the required 0.39; peak memory `2^72.389801`
bytes; preprocessing `2^120.858891`, charged in full with no reuse credit;
nonuniform advice 0; no auxiliary powers.

See `proof.md` for the reduction, the published-parameter reproduction, the
instruction ledger with every factor sourced, and the truncation argument.
`experiment.py` recomputes the GS21 cost and the submitted hard-capped ledger
from the Table 10 parameters; it is arithmetic evidence, not an attack run.
