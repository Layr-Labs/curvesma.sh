# nfs-stnfs-v4 (fixed-G1) — hard-capped STNFS discrete log in GT

Track: `bls12-381-fixed-g1-invert-g2`. Claim: **123.980693**
(`123980693/1000000`, log2 target-group-operation equivalents), success >= 1/2.

Reduction: the challenge fixes `P = [b]P0` in G1 and publishes
`z = e(P,Q0)^a`. Compute `g = e(P,Q0) = pairing(Q0,P,final_exponentiate=True)`
with the pinned library's `(G2 point, G1 point)` argument order; `g` generates
the prime-order group GT. Solve `a = log_g z` by Special-Tower NFS in `F_{p^12}`
with the published GS21 BLS12-381 parameters, return `Q = [a]Q0`, and verify
`e(P,Q) = z` and `[r]Q = O`. Since `a` is unique modulo `r`, this is correct for
every `z` in GT, including `z = 1` (`a = 0`, `Q = O`, `e(P,O) = 1`).

Ledger (all terms in `proof.md`, recomputed by `experiment.py`):

```text
GS21 core (sieving + block-Wiedemann)   2^126.146602479037
individual-logarithm descent (BP14)     2^122.641786367599   <- charged, not absorbed
ECM cofactorization                     2^87.6
model total T                           2^126.268405376626
x 1024 instructions per model operation (11.3x the published calibration bound)
per-attempt hard cap B                  2^136.268405376626
E[T] = 2B, R_cap = 4B                   2^138.268405376626
W_cap = R_cap/20000 + 2^14              2^123.980692997077
score                                  123.980693
```

Peak storage `2^72.389800` bytes (relation file; matrix `2^69.17`). No
dictionary, no distinguished points, no auxiliary powers, no `n` parameter, no
amortization over future targets, nonuniform advice 0.

Files: `submission.json` (claim + experiment declaration), `claim.json`
(structured claim), `proof.md`, `experiment.py` (recomputes every number and
asserts 15 conditions), `recomputed.json` (its output), `verify_g1.py` and
`toy-g1.json` (the reduction checked end to end on production-size elements with
the pinned library, seed 20261023).

This revision answers the four accounting items raised against the earlier
STNFS filings on this account: the descent base (`exp(1.44*59) = 2^122.6`, not
`2^85`), the instruction conversion (derived, then multiplied by 11.3x), the
Markov input (`E[T]` built from a hard per-attempt cap, not a point estimate),
and the storage schedule (records, sort passes and matrix sized and charged).
The G1/G2 role swap relative to the fixed-G2 filing changes only the reduction;
the GT logarithm, which is the whole cost, is identical.
