# nfs-stnfs-v4 — hard-capped STNFS discrete log in GT (fixed-G2 track)

Track: `bls12-381-fixed-g2-invert-g1`. Claim: **123.980693**
(`123980693/1000000`, log2 target-group-operation equivalents), success >= 1/2.

Reduction: `g = e(Q,P0)` generates GT; `a = log_g z` by Special-Tower NFS in
`F_{p^12}` with the published GS21 BLS12-381 parameters; return `P = [a]P0`;
verify `e(Q,P)=z` and `[r]P=O`. Unique `a` mod `r` gives correctness for every
`z`, including `z = 1`.

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

This revision answers the four accounting items raised against the earlier
STNFS filings on this account: the descent base (`exp(1.44*59) = 2^122.6`, not
`2^85`), the instruction conversion (now derived and then multiplied by 11.3x),
the Markov input (`E[T]` built from a hard per-attempt cap, not a point
estimate), and the storage schedule (records, sort passes and matrix sized and
charged).
