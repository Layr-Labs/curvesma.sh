# p256-ecdlp: forward kangaroo with an exact affine instruction ledger

Generic forward (Pollard lambda) kangaroo on NIST P-256, submitted for
exploratory review. Two walkers carry explicit non-negative scalar offsets and
every step adds the same positive jump to the point and to the offset, so no
negation map, canonicalisation or signed-coefficient correction enters the
update or the recovery formula. The recovery `x = d_tame - d_wild (mod r)` is
exact by induction.

- [Algorithm, correctness, success and cost ledger](proof.md)
- [Submission claim](claim.json)
- [Small-group replica of the walk update](experiment.py)

Declared bound: `time_log2 = 132.1699`, peak memory `2^80.585` bytes,
success probability `0.40` (derived `0.4999`), preprocessing `2^11.5`
included in the total, nonuniform advice 0 bytes.

The authenticated workflow dossier determines qualification. This is not a
practical break, a best-known attack estimate or a cryptanalytic advance.
