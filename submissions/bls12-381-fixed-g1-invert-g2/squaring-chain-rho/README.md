# Cyclotomic squaring-chain distinguished-point walk

Package for `bls12-381-fixed-g1-invert-g2`. Parent commit of the public tree:
`6423c2728dc03d5f2d647c7f22adf1dac0643a66`.

`proof.md` carries the reduction, the block walk, the collision algebra, the single
heuristic, the instruction ledger and the exact integers. `submission.json` declares
the one claim and its exact rational score. `experiment.py` is a toy check of the
block update only, at the scope of a small cyclic group.

The claimed bound is `time_log2 = 127.553439` in target-group-operation
equivalents, with peak storage `2^119.524868` bytes, preprocessing
`2^13.356590`, nonuniform advice 0 bytes and a claimed success probability of 0.5.

The mechanism is one dense Fp12 multiplication per block plus three Granger-Scott
squarings in the cyclotomic subgroup. The three extra elements per block are
credited only through the seven fixed ratio events they can create between block
pairs; nothing is divided by the chain length or by any automorphism orbit.
