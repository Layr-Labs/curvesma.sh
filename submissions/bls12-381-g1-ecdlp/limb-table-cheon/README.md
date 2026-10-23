# limb-table-cheon

Candidate for `bls12-381-g1-ecdlp`. The scored claim is in `submission.json`
(claim 0, score `148987/1250`, i.e. `time_log2 = 119.1896`) and the complete
argument, cost ledger, coverage count and exact integer checks are in
`proof.md`.

The attack is the Cheon `d`-th power reduction at `d = 2^20` (supplied as `A_n`
at the scoring length), with a baby-step giant-step rectangle inside the
order-`m = (r-1)/d` subgroup, quotiented by the `mu_6` endomorphism orbit. Every
point is a four-limb fixed-base evaluation against affine tables, and the two
sides use distinct tables: multiples of `P0` for the baby side, multiples of
`A_d` for the giant side. Matching is a deterministic radix sort plus merge on
the canonical orbit key, not a hash table.

Nothing outside this directory is changed except the track selector.
