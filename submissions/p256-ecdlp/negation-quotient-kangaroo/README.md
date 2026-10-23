# p256-ecdlp: negation-quotient kangaroo

Candidate package for exploratory review on the `p256-ecdlp` track. The attack is
a two-animal Pollard kangaroo walk run on the quotient of the P-256 group by the
negation map, with every P-256 operation charged through explicit AArch64
instruction bounds divided by the track's fixed constant `C = 1000`.

- [Proof, success argument and instruction ledger](proof.md)
- [Claim record](claim.json)

Submitted `time_log2` is **131.435096**. Peak memory is `2^134 + 2^20` bytes and
is mandatory for the stored tame trail. This is a fully charged generic-interval
attack; it is not a practical break, a best-known-attack estimate, or a
cryptanalytic advance, and no production run is claimed.
